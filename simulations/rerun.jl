# rerun.jl -- resolved reruns of Figs 2-3, convergence tests, and robustness checks.
#
# Solves the same equations as two_dimensional_module_interaction.jl (frac = true, equal and
# constant D, K = 1, explicit Euler, growth then diffusion, 5-point Laplacian, zero-flux edges),
# but with a fused in-place kernel and without storing the full density field: the front
# height h(x,t) and front composition f(x,t) are recorded on the fly, and invasion velocities
# are fitted at the end of the run. Run `selftest` once to confirm that the kernel matches
# the original module step for step.
#
# usage:
#   julia rerun.jl count MODE           # number of tasks in MODE (used by submit_rerun.sh)
#   julia rerun.jl MODE TASK_ID         # run one task (0-based)
#   julia rerun.jl list MODE            # print the parameter table for MODE
#   julia rerun.jl selftest             # compare the kernel with the original module
# MODE: fig2a | fig2b | fig2c | fig3 | conv | ic | strip
# env:  OUT_DIR (default ./results), SAVE_SNAPSHOT=1 to also save the final density fields

using Printf
using LinearAlgebra

const MODULE_FILE = "two_dimensional_module_interaction.jl"

# --------------------------------------------------------------------------------------
# parameters
# --------------------------------------------------------------------------------------
Base.@kwdef struct Run
    mode::String = ""
    task::Int = 0
    label::String = ""
    rf::Float64 = 1.0
    rs::Float64 = 0.8
    afs::Float64 = 0.0       # interaction_12: benefit to the fast type from the slow type
    asf::Float64 = 0.0       # interaction_21: benefit to the slow type from the fast type
    D::Float64 = 1.0
    dx::Float64 = 1.0
    dt::Float64 = 0.025
    T::Float64 = 400.0
    Lx::Float64 = 1600.0
    Ly::Float64 = 950.0
    rec::Float64 = 1.0       # recording interval
    init::Symbol = :half     # :half (fast left, slow right) or :strip (fast strip in slow)
    h0::Float64 = 10.0       # initial height of the populated band
    offset::Float64 = 0.0    # extra initial height of the fast type (:half only)
    width::Float64 = 0.0     # strip width (:strip only)
end

function tasks(mode::AbstractString)
    R = Run[]
    add!(; kw...) = push!(R, Run(; mode=String(mode), task=length(R), kw...))
    if mode == "fig2a"            # a_fs = -a_sf = a,  Delta r = 0.2
        for a in range(-0.5, 0.5, length=51)
            add!(rs=0.8, afs=a, asf=-a)
        end
    elseif mode == "fig2b"        # a_fs fixed, a_sf varied,  Delta r = 0.2
        dr = 0.2
        for afs in (-dr, -4dr, -16dr), asf in range(-0.5, 0.5, length=41)
            add!(rs=1 - dr, afs=afs, asf=asf, label=@sprintf("afs%+.2f", afs))
        end
    elseif mode == "fig2c"        # a_sf fixed, a_fs varied,  Delta r = 0.2
        dr = 0.2
        for asf in (-dr, -4dr, -16dr), afs in range(-0.5, 0.5, length=41)
            add!(rs=1 - dr, afs=afs, asf=asf, label=@sprintf("asf%+.2f", asf))
        end
    elseif mode == "fig3"         # three trajectories of Fig 3C,  Delta r = 0.2
        for asf in range(0.0, 0.5, length=26)
            add!(afs=-0.1, asf=asf, label="traj1")     # circular bulge -> parasitic mixture
        end
        for asf in range(0.0, 0.5, length=26)
            add!(afs=0.2, asf=asf, label="traj2")      # composite bulge -> mutualistic mixture
        end
        for afs in range(-0.5, 0.5, length=41)
            add!(afs=afs, asf=0.4, label="traj3")      # dent -> parasitic -> mutualistic
        end
    elseif mode == "conv"         # one point per regime: time step, grid spacing, run length
        for (lab, afs, asf) in (("dent", -0.8, 0.4), ("circ", -0.8, 0.0), ("comp", 0.3, -0.2))
            for dt in (0.1, 0.05, 0.025, 0.0125)
                add!(afs=afs, asf=asf, dt=dt, label=lab)
            end
            add!(afs=afs, asf=asf, dx=0.5, dt=0.0125, label=lab * "_dx0.5")
            add!(afs=afs, asf=asf, T=1000.0, Lx=3200.0, Ly=2150.0, label=lab * "_long")
        end
    elseif mode == "ic"           # does the outcome near a_sf = Delta r depend on which type starts ahead?
        for asf in (0.15, 0.25, 0.3), off in (-20.0, 0.0, 20.0)
            add!(afs=-0.8, asf=asf, h0=30.0, offset=off, label=@sprintf("off%+g", off))
        end
    elseif mode == "strip"        # do narrow fast sectors survive in the bistable region?
        for asf in (0.0, 0.1), w in (2.0, 4.0, 8.0, 16.0, 32.0, 64.0)
            add!(afs=-0.8, asf=asf, init=:strip, width=w, Lx=1200.0, label=@sprintf("w%g", w))
        end
    else
        error("unknown mode $mode")
    end
    return R
end

# Eqs (8) and (11) of the main text. Returns (regime, u, u2): for coexistence u = u_ms
# (mixture into slow) and u2 = u_mf (mixture into fast); otherwise u2 = NaN.
function theory(p::Run)
    dr = p.rf - p.rs; A = p.afs + p.asf; D = p.D
    if p.afs > -dr && p.asf > dr
        ums = p.afs <= 0 ? 2sqrt(D * p.asf * (dr + p.afs) / A) : 2sqrt(D * (dr + p.afs))
        return ("coex", ums, 2sqrt(D * (p.asf - dr)))
    elseif p.asf > dr
        return ("dent", -2sqrt(D * (p.asf - dr)), NaN)
    else
        return (p.afs <= 0 ? "circ" : "comp", 2sqrt(D * max(dr, dr + p.afs)), NaN)
    end
end

# --------------------------------------------------------------------------------------
# kernel: one time step, identical to the module's update_densities! with frac = true
# --------------------------------------------------------------------------------------
function step!(n1, n2, t1, t2, p)
    Nx, Ny = size(n1)
    dt = p.dt; c = p.D * p.dt / p.dx^2
    rf, rs, afs, asf = p.rf, p.rs, p.afs, p.asf
    @inbounds for j in 1:Ny, i in 1:Nx                 # growth (in place)
        a = n1[i, j]; b = n2[i, j]
        tot = a + b
        tt = tot == 0.0 ? 1.0 : tot
        g = 1.0 - a - b
        da = a * (rf + afs * b / tt) * g
        db = b * (rs + asf * a / tt) * g
        n1[i, j] = a + da * dt
        n2[i, j] = b + db * dt
    end
    @inbounds for j in 1:Ny                            # diffusion into the scratch arrays
        jm = max(j - 1, 1); jp = min(j + 1, Ny)
        for i in 1:Nx
            il = max(i - 1, 1); ir = min(i + 1, Nx)
            v = n1[i, j]
            w = v + c * (n1[il, j] + n1[ir, j] + n1[i, jm] + n1[i, jp] - 4v)
            t1[i, j] = w < 0.0 ? 0.0 : w
            v = n2[i, j]
            w = v + c * (n2[il, j] + n2[ir, j] + n2[i, jm] + n2[i, jp] - 4v)
            t2[i, j] = w < 0.0 ? 0.0 : w
        end
    end
    return t1, t2, n1, n2                              # new state first, scratch second
end

# front height (largest y with n > thr) and composition there, for every column
function front!(H, F, k, n1, n2, thr, dx)
    Nx, Ny = size(n1)
    @inbounds for i in 1:Nx
        jj = 0
        for j in Ny:-1:1
            if n1[i, j] + n2[i, j] > thr
                jj = j
                break
            end
        end
        if jj == 0
            H[i, k] = NaN32; F[i, k] = NaN32
        else
            H[i, k] = (jj - 1) * dx
            F[i, k] = n1[i, jj] / (n1[i, jj] + n2[i, jj])
        end
    end
end

# first crossing from the left where f drops below c; last crossing from the right where f exceeds c
function first_cross(f, xs, c)
    @inbounds for i in 2:length(f)
        a, b = f[i-1], f[i]
        if isfinite(a) && isfinite(b) && a >= c && b < c
            return xs[i-1] + (c - a) * (xs[i] - xs[i-1]) / (b - a)
        end
    end
    return NaN
end
function last_cross(f, xs, c)
    @inbounds for i in length(f)-1:-1:1
        a, b = f[i], f[i+1]
        if isfinite(a) && isfinite(b) && a > c && b <= c
            return xs[i] + (a - c) * (xs[i+1] - xs[i]) / (a - b)
        end
    end
    return NaN
end

function slope(t, x)
    ok = isfinite.(x)
    t, x = t[ok], x[ok]
    length(t) < 3 && return NaN
    tm = sum(t) / length(t); xm = sum(x) / length(x)
    return sum((t .- tm) .* (x .- xm)) / sum((t .- tm) .^ 2)
end
function slope_sqrt(t, x)                       # x = x0 + u t + c sqrt(t); returns u
    ok = isfinite.(x)
    t, x = t[ok], x[ok]
    length(t) < 5 && return NaN
    A = hcat(ones(length(t)), t, sqrt.(t))
    return (A \ x)[2]
end

# --------------------------------------------------------------------------------------
# one run
# --------------------------------------------------------------------------------------
function simulate(p::Run)
    xs = collect(-p.Lx/2:p.dx:p.Lx/2)
    ys = collect(0.0:p.dx:p.Ly)
    Nx, Ny = length(xs), length(ys)
    n1 = zeros(Nx, Ny); n2 = zeros(Nx, Ny)
    t1 = zeros(Nx, Ny); t2 = zeros(Nx, Ny)
    nrow(h) = clamp(round(Int, h / p.dx), 1, Ny)
    if p.init == :half
        half = div(Nx, 2)
        n1[1:half, 1:nrow(p.h0 + p.offset)] .= 0.99
        n2[half+1:end, 1:nrow(p.h0)] .= 0.99
    else
        for i in 1:Nx
            if abs(xs[i]) <= p.width / 2
                n1[i, 1:nrow(p.h0)] .= 0.99
            else
                n2[i, 1:nrow(p.h0)] .= 0.99
            end
        end
    end

    nsteps = round(Int, p.T / p.dt)
    every = max(1, round(Int, p.rec / p.dt))
    nrec = nsteps ÷ every + 1
    ts = zeros(nrec)
    H99 = Matrix{Float32}(undef, Nx, nrec); F99 = similar(H99)
    H50 = similar(H99); F50 = similar(H99)
    front!(H99, F99, 1, n1, n2, 0.99, p.dx)
    front!(H50, F50, 1, n1, n2, 0.5, p.dx)

    k = 1
    for s in 1:nsteps
        n1, n2, t1, t2 = step!(n1, n2, t1, t2, p)
        if s % every == 0 && k < nrec
            k += 1
            ts[k] = s * p.dt
            front!(H99, F99, k, n1, n2, 0.99, p.dx)
            front!(H50, F50, k, n1, n2, 0.5, p.dx)
        end
    end
    return xs, ys, ts[1:k], H99[:, 1:k], F99[:, 1:k], H50[:, 1:k], F50[:, 1:k], n1, n2
end

function analyze(p::Run, xs, ts, H, F)
    regime, _, _ = theory(p)
    dr = p.rf - p.rs; A = p.afs + p.asf
    coex = regime == "coex"
    fm = coex ? (dr + p.afs) / A : 0.5
    cL = coex ? (1 + fm) / 2 : 0.5            # fast | (mixture or slow)
    cR = coex ? fm / 2 : 0.5                  # (mixture or fast) | slow
    nrec = length(ts)
    xL = fill(NaN, nrec); xR = fill(NaN, nrec); width = fill(NaN, nrec); hmax = fill(NaN, nrec)
    for k in 1:nrec
        f = @view F[:, k]
        xL[k] = first_cross(f, xs, cL)
        xR[k] = last_cross(f, xs, cR)
        width[k] = count(v -> isfinite(v) && v > 0.5, f) * p.dx
        hk = filter(isfinite, @view H[:, k])
        hmax[k] = isempty(hk) ? NaN : maximum(hk)
    end

    # records are usable while the front is away from the top wall and the boundaries are
    # away from the side walls; use the first contiguous stretch of usable records
    margin = 20.0
    ok_top = isfinite.(hmax) .& (hmax .< p.Ly - 30.0)
    if p.init == :half
        ok_side = isfinite.(xL) .& isfinite.(xR) .& (xL .> xs[1] + margin) .& (xR .< xs[end] - margin)
    else
        ok_side = (width .> 0) .& (width .< p.Lx - 2margin)
    end
    valid = ok_top .& ok_side
    k0 = findfirst(valid)
    if k0 === nothing
        return (xL=xL, xR=xR, width=width, hmax=hmax, valid=valid, te=NaN,
                u=(NaN, NaN, NaN), u2=(NaN, NaN, NaN), vf=NaN, vs=NaN)
    end
    k1 = k0
    while k1 < nrec && valid[k1+1]
        k1 += 1
    end
    for k in 1:nrec                       # mark only the contiguous stretch as valid
        valid[k] = k0 <= k <= k1
    end
    te = ts[k1]
    late = (ts .>= te / 2) .& valid
    mid = (ts .>= te / 4) .& (ts .<= te / 2) .& valid
    wide = (ts .>= te / 5) .& valid

    fits(x, sgn) = (sgn * slope(ts[late], x[late]), sgn * slope(ts[mid], x[mid]),
                    sgn * slope_sqrt(ts[wide], x[wide]))
    if p.init == :strip
        u = fits(width, 0.5)                        # speed of each boundary of the strip
        u2 = (NaN, NaN, NaN)
    elseif coex
        u = fits(xR, 1.0)                           # u_ms: mixture into slow (moves right)
        u2 = fits(xL, -1.0)                         # u_mf: mixture into fast (moves left)
    else
        u = fits(xL, 1.0)                           # > 0 when the fast type advances
        u2 = (NaN, NaN, NaN)
    end

    iw = max(1, round(Int, margin / p.dx))           # columns near the side walls
    iF = p.init == :strip ? div(length(xs), 2) + 1 : iw
    iS = length(xs) - iw
    vf = slope(ts[late], Float64.(H[iF, late]))
    vs = slope(ts[late], Float64.(H[iS, late]))
    return (xL=xL, xR=xR, width=width, hmax=hmax, valid=valid, te=te, u=u, u2=u2, vf=vf, vs=vs)
end

fmt(x) = isa(x, AbstractString) || isa(x, Symbol) ? string(x) : @sprintf("%.6g", x)

function run_task(mode, task_id)
    R = tasks(mode)
    0 <= task_id < length(R) || error("task $task_id out of range 0:$(length(R)-1) for $mode")
    p = R[task_id+1]
    regime, uth, u2th = theory(p)
    println("mode=$mode task=$task_id label=$(p.label) regime=$regime rf=$(p.rf) rs=$(p.rs) ",
            "afs=$(p.afs) asf=$(p.asf) D=$(p.D) dx=$(p.dx) dt=$(p.dt) T=$(p.T) Lx=$(p.Lx) Ly=$(p.Ly)")
    flush(stdout)

    t0 = time()
    xs, ys, ts, H99, F99, H50, F50, n1, n2 = simulate(p)
    runtime = time() - t0
    a99 = analyze(p, xs, ts, H99, F99)
    a50 = analyze(p, xs, ts, H50, F50)

    out = joinpath(get(ENV, "OUT_DIR", joinpath(pwd(), "results")), mode)
    mkpath(out)
    stem = joinpath(out, @sprintf("task_%04d", task_id))

    cols = ["mode", "task", "label", "regime", "rf", "rs", "afs", "asf", "D", "dx", "dt", "T",
            "Lx", "Ly", "init", "h0", "offset", "width",
            "u_theory", "u2_theory", "vf_theory", "vs_theory",
            "t_valid_end", "u_late", "u_mid", "u_sqrt", "u2_late", "u2_mid", "u2_sqrt",
            "vf_meas", "vs_meas", "u_late_thr50", "u2_late_thr50", "final_width", "runtime_s"]
    vals = Any[mode, task_id, p.label, regime, p.rf, p.rs, p.afs, p.asf, p.D, p.dx, p.dt, p.T,
               p.Lx, p.Ly, p.init, p.h0, p.offset, p.width,
               uth, u2th, 2sqrt(p.D * p.rf), 2sqrt(p.D * p.rs),
               a99.te, a99.u..., a99.u2..., a99.vf, a99.vs, a50.u[1], a50.u2[1],
               a99.width[end], runtime]
    open(stem * ".csv", "w") do io
        println(io, join(cols, ","))
        println(io, join(fmt.(vals), ","))
    end
    open(stem * "_traj.csv", "w") do io
        println(io, "t,xL_99,xR_99,width_99,hmax_99,valid_99,xL_50,xR_50")
        for k in eachindex(ts)
            println(io, join(fmt.(Any[ts[k], a99.xL[k], a99.xR[k], a99.width[k], a99.hmax[k],
                                      Int(a99.valid[k]), a50.xL[k], a50.xR[k]]), ","))
        end
    end
    if get(ENV, "SAVE_SNAPSHOT", "0") == "1"           # final fields as raw Float32, (Nx, Ny) column-major
        open(io -> write(io, Float32.(n1)), stem * "_n1.f32", "w")
        open(io -> write(io, Float32.(n2)), stem * "_n2.f32", "w")
        open(io -> println(io, "Nx=$(length(xs)) Ny=$(length(ys)) dx=$(p.dx)"), stem * "_shape.txt", "w")
    end
    @printf("done in %.1f s: t_valid_end=%.0f  u_late=%.4f  (theory %.4f)  u2_late=%.4f  (theory %.4f)  vf=%.4f vs=%.4f\n",
            runtime, a99.te, a99.u[1], uth, a99.u2[1], u2th, a99.vf, a99.vs)
end

# --------------------------------------------------------------------------------------
# self-test against the original module
# --------------------------------------------------------------------------------------
function selftest()
    path = joinpath(@__DIR__, MODULE_FILE)
    isfile(path) || error("put $MODULE_FILE next to rerun.jl for the self-test")
    include(path)                                   # needs Distributions installed
    M = getfield(Main, :two_type_two_dimensional_invasion)
    p = Run(rs=0.8, afs=0.3, asf=0.5, D=0.1, dx=1.0, dt=0.1, Lx=60.0, Ly=60.0)
    par = Base.invokelatest(M.Parameters; Lx=p.Lx, dx=p.dx, Ly=p.Ly, dy=p.dx, T=10.0, dt=p.dt,
        growth_rate_1=p.rf, growth_rate_2=p.rs, carrying_capacity_1=1.0, carrying_capacity_2=1.0,
        interaction_12=p.afs, interaction_21=p.asf, diffusion_constant_1=p.D, diffusion_constant_2=p.D,
        mutation_rate_12=0.0, mutation_rate_21=0.0, frac=true, thin_by=1, stochastic=false, sigma=0.1)
    Nx, Ny = length(par.xs), length(par.ys)
    d1 = zeros(Nx, Ny); d2 = zeros(Nx, Ny)
    d1[1:div(Nx, 2), 1:10] .= 0.99; d2[div(Nx, 2)+1:end, 1:10] .= 0.99
    d1[div(Nx, 2)-3:div(Nx, 2)+3, 11:14] .= 0.3     # an irregular patch, to exercise all terms
    n1 = copy(d1); n2 = copy(d2); t1 = zeros(Nx, Ny); t2 = zeros(Nx, Ny)
    tmp1 = zeros(Nx, Ny); tmp2 = zeros(Nx, Ny)
    for s in 1:300
        Base.invokelatest(M.update_densities!, d1, d2, tmp1, tmp2, par)
        n1, n2, t1, t2 = step!(n1, n2, t1, t2, p)
    end
    err = max(maximum(abs.(d1 .- n1)), maximum(abs.(d2 .- n2)))
    @printf("max |module - kernel| after 300 steps: %.3e  ->  %s\n", err, err < 1e-10 ? "PASS" : "FAIL")
end

# --------------------------------------------------------------------------------------
function main(args)
    isempty(args) && error("usage: julia rerun.jl (count MODE | list MODE | selftest | MODE TASK_ID)")
    if args[1] == "count"
        println(length(tasks(args[2])))
    elseif args[1] == "list"
        for p in tasks(args[2])
            r, u, u2 = theory(p)
            @printf("%4d  %-10s %-5s rs=%.2f afs=%+.3f asf=%+.3f dx=%.3g dt=%.4g T=%g  u_th=%+.4f u2_th=%.4f\n",
                    p.task, p.label, r, p.rs, p.afs, p.asf, p.dx, p.dt, p.T, u, u2)
        end
    elseif args[1] == "selftest"
        selftest()
    else
        run_task(args[1], parse(Int, args[2]))
    end
end

main(ARGS)
