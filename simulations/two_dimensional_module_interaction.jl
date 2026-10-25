module two_type_two_dimensional_invasion
export Parameters, run_simulation!, calculate_invasion_velocity

using Distributions

struct Parameters
    Lx::Float64 # Length of domain
    dx::Float64 # lattice spacing
    Ly::Float64 # lattice in y direciton length
    dy::Float64 # spacing in the y direction
    T::Float64 # time running simulation
    dt::Float64 # time step
    growth_rate_1::Float64 # growth rate for type 1
    growth_rate_2::Float64 # and for type 2
    carrying_capacity_1::Float64 # carrying capacity of type 1 
    carrying_capacity_2::Float64 # and for type 2
    interaction_12::Float64 # specifies density mediated interactions 
    interaction_21::Float64 # specifies density mediated interactions 
    diffusion_constant_1::Float64 # rate of diffusion in space for type 1
    diffusion_constant_2::Float64 # rate of diffusion in space for type 2
    mutation_rate_12::Float64 # rate of differentiation between the types
    mutation_rate_21::Float64 # likewise
    frac::Bool
    thin_by::Int64 # save every thin_by'th time step
    xs::Vector{Float64} # grid of spatial locations
    ys::Vector{Float64}
    stochastic::Bool
    sigma::Float64
    # Custom constructor to initialize xs based on L and dx
    function Parameters(;Lx::Float64, dx::Float64, Ly::Float64, dy::Float64,  T::Float64, dt::Float64, 
                        growth_rate_1::Float64, growth_rate_2::Float64, 
                        carrying_capacity_1::Float64, carrying_capacity_2::Float64, 
                        interaction_12::Float64, interaction_21::Float64, 
                        diffusion_constant_1::Float64, diffusion_constant_2::Float64, 
                        mutation_rate_12::Float64, mutation_rate_21::Float64, 
                        frac::Bool=true, 
                        thin_by::Int64=1, stochastic::Bool=false, sigma::Float64=1.0)
        xs = collect(-Lx/2:dx:Lx/2) # Create the grid
        ys = collect(-Ly/2:dy:Ly/2)
        new(Lx, dx, Ly, dy, T, dt, growth_rate_1, growth_rate_2, carrying_capacity_1, carrying_capacity_2,
            interaction_12, interaction_21, diffusion_constant_1, diffusion_constant_2,
            mutation_rate_12, mutation_rate_21, frac, thin_by, xs, ys, stochastic, sigma)
    end
end

function run_simulation!(density_1, density_2, par::Parameters)
    # these are dummy arrays to minimize memory allocation
    temp_data_1 = zeros(Float64, size(density_1))
    temp_data_2 = zeros(Float64, size(density_2))

    record_times = 0:par.dt*par.thin_by:par.T # store the times at which we will record the result
    result = fill(NaN,length(par.xs), length(par.ys), length(record_times),2) # we will store the result here

    result[:,:,1,1] .= density_1; result[:,:,1,2] .= density_2 # put the initial conditions in result

    # This loop does the simulation up to the record times and then records them
    for time_idx = 2:length(record_times)
        update_densities!(density_1, density_2, temp_data_1, temp_data_2, par) 
        result[:,:,time_idx,1] .= density_1; result[:,:,time_idx,2] .= density_2
    end 

    return result, record_times
end

function update_densities!(density_1, density_2, temp_data_1, temp_data_2, par::Parameters)
    
    for i = 1:par.thin_by
        if par.frac
            apply_growth_and_mutation_frac!(density_1, density_2, temp_data_1, temp_data_2, par)
            apply_diffusion_frac!(density_1, density_2, temp_data_1, temp_data_2, par)
        else
            apply_growth_and_mutation_den!(density_1, density_2, temp_data_1, temp_data_2, par)
            apply_diffusion_den!(density_1, density_2, temp_data_1, temp_data_2, par)
        end

        if par.stochastic
            apply_demographic_noise!(density_1, density_2, par)
        end
        # remove negative values
        map!(v -> v < 0 ? 0.0 : v, density_1, density_1)
        map!(v -> v < 0 ? 0.0 : v, density_2, density_2)
    end

end

function apply_diffusion!(density_1, density_2, temp_data_1, temp_data_2, par::Parameters)

    temp_data_1 .= 0
    temp_data_2 .= 0
    # pre compute important quantities
    dx, dy, dt = par.dx, par.dy, par.dt
    D1, D2 = par.diffusion_constant_1, par.diffusion_constant_2
    coeff1x = D1 * dt / dx^2
    coeff1y = D1 * dt / dy^2
    coeff2x = D2 * dt / dx^2
    coeff2y = D2 * dt / dy^2

    # Use temp_data for the diffusion calculation
    for j = 2:length(par.ys)-1
        for i = 2:length(par.xs)-1
            temp_data_1[i,j] = density_1[i,j] + coeff1x*(density_1[i-1,j] - 2 * density_1[i,j] + density_1[i+1,j]) + coeff1y*(density_1[i,j-1] - 2 * density_1[i,j] + density_1[i,j+1])
            temp_data_2[i,j] = density_2[i,j] + coeff2x*(density_2[i-1,j] - 2 * density_2[i,j] + density_2[i+1,j]) + coeff2y*(density_2[i,j-1] - 2 * density_2[i,j] + density_2[i,j+1])
        end 
    end
    # if any(isnan.(temp_data_1))
    #     println("NAN1")
    # end
    for i = 2:length(par.xs)-1
        temp_data_1[i,1] = density_1[i,1] + coeff1x*(density_1[i-1,1] - 2*density_1[i,1] + density_1[i+1,1]) + coeff1y*(density_1[i,1] - 2*density_1[i,1] + density_1[i,2])
        temp_data_1[i,end] = density_1[i,end] + coeff1x*(density_1[i-1,end] - 2*density_1[i,end] + density_1[i+1,end]) + coeff1y*(density_1[i,end-1] - 2*density_1[i,end] + density_1[i,end])
    
        temp_data_2[i,1] = density_2[i,1] + coeff2x*(density_2[i-1,1] - 2*density_2[i,1] + density_2[i+1,1]) + coeff2y*(density_2[i,1] - 2*density_2[i,1] + density_2[i,2])
        temp_data_2[i,end] = density_2[i,end] + coeff2x*(density_2[i-1,end] - 2*density_2[i,end] + density_2[i+1,end]) + coeff2y*(density_2[i,end-1] - 2*density_2[i,end] + density_2[i,end])
    end
    # if any(isnan.(temp_data_1))
    #     println("NAN2")
    # end
    # Left and right boundaries
    for j = 2:length(par.ys)-1
        temp_data_1[1,j] = density_1[1,j] + coeff1x*(density_1[1,j] - 2*density_1[1,j] + density_1[2,j]) + coeff1y*(density_1[1,j-1] - 2*density_1[1,j] + density_1[1,j+1])
        temp_data_1[end,j] = density_1[end,j] + coeff1x*(density_1[end-1,j] - 2*density_1[end,j] + density_1[end,j]) + coeff1y*(density_1[end,j-1] - 2*density_1[end,j] + density_1[end,j+1])
    
        temp_data_2[1,j] = density_2[1,j] + coeff2x*(density_2[1,j] - 2*density_2[1,j] + density_2[2,j]) + coeff2y*(density_2[1,j-1] - 2*density_2[1,j] + density_2[1,j+1])
        temp_data_2[end,j] = density_2[end,j] + coeff2x*(density_2[end-1,j] - 2*density_2[end,j] + density_2[end,j]) + coeff2y*(density_2[end,j-1] - 2*density_2[end,j] + density_2[end,j+1])
    end
    # if any(isnan.(temp_data_1))
    #     println("NAN3")
    # end
    # Top-left corner
    temp_data_1[1,1] = density_1[1,1] + coeff1x*(density_1[1,1] - 2*density_1[1,1] + density_1[2,1]) + coeff1y*(density_1[1,1] - 2*density_1[1,1] + density_1[1,2])
    temp_data_2[1,1] = density_2[1,1] + coeff2x*(density_2[1,1] - 2*density_2[1,1] + density_2[2,1]) + coeff2y*(density_2[1,1] - 2*density_2[1,1] + density_2[1,2])
    # if any(isnan.(temp_data_1))
    #     println("NAN4")
    # end
    # Top-right corner
    temp_data_1[end,1] = density_1[end,1] + coeff1x*(density_1[end-1,1] - 2*density_1[end,1] + density_1[end,1]) + coeff1y*(density_1[end,1] - 2*density_1[end,1] + density_1[end,2])
    temp_data_2[end,1] = density_2[end,1] + coeff2x*(density_2[end-1,1] - 2*density_2[end,1] + density_2[end,1]) + coeff2y*(density_2[end,1] - 2*density_2[end,1] + density_2[end,2])
    # if any(isnan.(temp_data_1))
    #     println("NAN5")
    # end
    # Bottom-left corner
    temp_data_1[1,end] = density_1[1,end] + coeff1x*(density_1[1,end] - 2*density_1[1,end] + density_1[2,end]) + coeff1y*(density_1[1,end-1] - 2*density_1[1,end] + density_1[1,end])
    temp_data_2[1,end] = density_2[1,end] + coeff2x*(density_2[1,end] - 2*density_2[1,end] + density_2[2,end]) + coeff2y*(density_2[1,end-1] - 2*density_2[1,end] + density_2[1,end])
    # if any(isnan.(temp_data_1))
    #     println("NAN6")
    # end
    # Bottom-right corner
    temp_data_1[end,end] = density_1[end,end] + coeff1x*(density_1[end-1,end] - 2*density_1[end,end] + density_1[end,end]) + coeff1y*(density_1[end,end-1] - 2*density_1[end,end] + density_1[end,end])
    temp_data_2[end,end] = density_2[end,end] + coeff2x*(density_2[end-1,end] - 2*density_2[end,end] + density_2[end,end]) + coeff2y*(density_2[end,end-1] - 2*density_2[end,end] + density_2[end,end])
    # if any(isnan.(temp_data_1))
    #     println("NAN7")
    # end
    # Copy the updated values from temp_data back to the original density arrays
    density_1 .= temp_data_1
    density_2 .= temp_data_2
end


function apply_diffusion_den!(density_1, density_2, temp_data_1, temp_data_2, par::Parameters)
    # density-dependent diffusion D(n) = D*(n1 + n2)

    temp_data_1 .= 0
    temp_data_2 .= 0
    # pre compute important quantities
    dx, dy, dt = par.dx, par.dy, par.dt
    D1, D2 = par.diffusion_constant_1, par.diffusion_constant_2
    coeff1x = D1 * (density_1 .+ density_2) * dt / dx^2 / 2 # divide by 2 
    coeff1y = D1 * (density_1 .+ density_2) * dt / dy^2 / 2
    coeff2x = D2 * (density_1 .+ density_2) * dt / dx^2 / 2
    coeff2y = D2 * (density_1 .+ density_2) * dt / dy^2 / 2

    K_factor = (1 .- density_1 / par.carrying_capacity_1 .- density_2 / par.carrying_capacity_2 ) 
    K_factor = K_factor .* (K_factor .> 0) + 0 .* (K_factor .< 0)

    ed_1x = coeff1x .* density_1 #.* K_factor # effective density
    ed_1y = coeff1y .* density_1 #.* K_factor 
    ed_2x = coeff2x .* density_2 #.* K_factor
    ed_2y = coeff2y .* density_2 #.* K_factor



    # Use temp_data for the diffusion calculation
    for j = 2:length(par.ys)-1
        for i = 2:length(par.xs)-1
            temp_data_1[i,j] = density_1[i,j] + (ed_1x[i-1,j] - 2 * ed_1x[i,j] + ed_1x[i+1,j]) + (ed_1y[i,j-1] - 2 * ed_1y[i,j] + ed_1y[i,j+1])
            temp_data_2[i,j] = density_2[i,j] + (ed_2x[i-1,j] - 2 * ed_2x[i,j] + ed_2x[i+1,j]) + (ed_2y[i,j-1] - 2 * ed_2y[i,j] + ed_2y[i,j+1])
        end 
    end
    # if any(isnan.(temp_data_1))
    #     println("NAN1")
    # end
    for i = 2:length(par.xs)-1
        temp_data_1[i,1] = density_1[i,1] + (ed_1x[i-1,1] - 2*ed_1x[i,1] + ed_1x[i+1,1]) + (ed_1y[i,1] - 2*ed_1y[i,1] + ed_1y[i,2])
        temp_data_1[i,end] = density_1[i,end] + (ed_1x[i-1,end] - 2*ed_1x[i,end] + ed_1x[i+1,end]) + (ed_1y[i,end-1] - 2*ed_1y[i,end] + ed_1y[i,end])
    
        temp_data_2[i,1] = density_2[i,1] + (ed_2x[i-1,1] - 2*ed_2x[i,1] + ed_2x[i+1,1]) + (ed_2y[i,1] - 2*ed_2y[i,1] + ed_2y[i,2])
        temp_data_2[i,end] = density_2[i,end] + (ed_2x[i-1,end] - 2*ed_2x[i,end] + ed_2x[i+1,end]) + (ed_2y[i,end-1] - 2*ed_2y[i,end] + ed_2y[i,end])
    end
    # if any(isnan.(temp_data_1))
    #     println("NAN2")
    # end
    # Left and right boundaries
    for j = 2:length(par.ys)-1
        temp_data_1[1,j] = density_1[1,j] + (ed_1x[1,j] - 2*ed_1x[1,j] + ed_1x[2,j]) + (ed_1y[1,j-1] - 2*ed_1y[1,j] + ed_1y[1,j+1])
        temp_data_1[end,j] = density_1[end,j] + (ed_1x[end-1,j] - 2*ed_1x[end,j] + ed_1x[end,j]) + (ed_1y[end,j-1] - 2*ed_1y[end,j] + ed_1y[end,j+1])
    
        temp_data_2[1,j] = density_2[1,j] + (ed_2x[1,j] - 2*ed_2x[1,j] + ed_2x[2,j]) + (ed_2y[1,j-1] - 2*ed_2y[1,j] + ed_2y[1,j+1])
        temp_data_2[end,j] = density_2[end,j] + (ed_2x[end-1,j] - 2*ed_2x[end,j] + ed_2x[end,j]) + (ed_2y[end,j-1] - 2*ed_2y[end,j] + ed_2y[end,j+1])
    end
    # if any(isnan.(temp_data_1))
    #     println("NAN3")
    # end
    # Top-left corner
    temp_data_1[1,1] = density_1[1,1] + (ed_1x[1,1] - 2*ed_1x[1,1] + ed_1x[2,1]) + (ed_1y[1,1] - 2*ed_1y[1,1] + ed_1y[1,2])
    temp_data_2[1,1] = density_2[1,1] + (ed_2x[1,1] - 2*ed_2x[1,1] + ed_2x[2,1]) + (ed_2y[1,1] - 2*ed_2y[1,1] + ed_2y[1,2])
    # if any(isnan.(temp_data_1))
    #     println("NAN4")
    # end
    # Top-right corner
    temp_data_1[end,1] = density_1[end,1] + (ed_1x[end-1,1] - 2*ed_1x[end,1] + ed_1x[end,1]) + (ed_1y[end,1] - 2*ed_1y[end,1] + ed_1y[end,2])
    temp_data_2[end,1] = density_2[end,1] + (ed_2x[end-1,1] - 2*ed_2x[end,1] + ed_2x[end,1]) + (ed_2y[end,1] - 2*ed_2y[end,1] + ed_2y[end,2])
    # if any(isnan.(temp_data_1))
    #     println("NAN5")
    # end
    # Bottom-left corner
    temp_data_1[1,end] = density_1[1,end] + (ed_1x[1,end] - 2*ed_1x[1,end] + ed_1x[2,end]) + (ed_1y[1,end-1] - 2*ed_1y[1,end] + ed_1y[1,end])
    temp_data_2[1,end] = density_2[1,end] + (ed_2x[1,end] - 2*ed_2x[1,end] + ed_2x[2,end]) + (ed_2y[1,end-1] - 2*ed_2y[1,end] + ed_2y[1,end])
    # if any(isnan.(temp_data_1))
    #     println("NAN6")
    # end
    # Bottom-right corner
    temp_data_1[end,end] = density_1[end,end] + (ed_1x[end-1,end] - 2*ed_1x[end,end] + ed_1x[end,end]) + (ed_1y[end,end-1] - 2*ed_1y[end,end] + ed_1y[end,end])
    temp_data_2[end,end] = density_2[end,end] + (ed_2x[end-1,end] - 2*ed_2x[end,end] + ed_2x[end,end]) + (ed_2y[end,end-1] - 2*ed_2y[end,end] + ed_2y[end,end])
    # if any(isnan.(temp_data_1))
    #     println("NAN7")
    # end
    # Copy the updated values from temp_data back to the original density arrays
    density_1 .= temp_data_1
    density_2 .= temp_data_2
end


function apply_diffusion_frac!(density_1, density_2, temp_data_1, temp_data_2, par::Parameters)
    # density-dependent diffusion D(n) = D*(n1 + n2)

    temp_data_1 .= 0
    temp_data_2 .= 0
    # pre compute important quantities
    dx, dy, dt = par.dx, par.dy, par.dt
    D1, D2 = par.diffusion_constant_1, par.diffusion_constant_2
    coeff1x = D1 * dt / dx^2 
    coeff1y = D1 * dt / dy^2 
    coeff2x = D2 * dt / dx^2 
    coeff2y = D2 * dt / dy^2 

    K_factor = (1 .- density_1 / par.carrying_capacity_1 .- density_2 / par.carrying_capacity_2 ) 
    K_factor = K_factor .* (K_factor .> 0) + 0 .* (K_factor .< 0)

    ed_1x = coeff1x .* density_1 #.* K_factor # effective density
    ed_1y = coeff1y .* density_1 #.* K_factor 
    ed_2x = coeff2x .* density_2 #.* K_factor
    ed_2y = coeff2y .* density_2 #.* K_factor



    # Use temp_data for the diffusion calculation
    for j = 2:length(par.ys)-1
        for i = 2:length(par.xs)-1
            temp_data_1[i,j] = density_1[i,j] + (ed_1x[i-1,j] - 2 * ed_1x[i,j] + ed_1x[i+1,j]) + (ed_1y[i,j-1] - 2 * ed_1y[i,j] + ed_1y[i,j+1])
            temp_data_2[i,j] = density_2[i,j] + (ed_2x[i-1,j] - 2 * ed_2x[i,j] + ed_2x[i+1,j]) + (ed_2y[i,j-1] - 2 * ed_2y[i,j] + ed_2y[i,j+1])
        end 
    end
    # if any(isnan.(temp_data_1))
    #     println("NAN1")
    # end
    for i = 2:length(par.xs)-1
        temp_data_1[i,1] = density_1[i,1] + (ed_1x[i-1,1] - 2*ed_1x[i,1] + ed_1x[i+1,1]) + (ed_1y[i,1] - 2*ed_1y[i,1] + ed_1y[i,2])
        temp_data_1[i,end] = density_1[i,end] + (ed_1x[i-1,end] - 2*ed_1x[i,end] + ed_1x[i+1,end]) + (ed_1y[i,end-1] - 2*ed_1y[i,end] + ed_1y[i,end])
    
        temp_data_2[i,1] = density_2[i,1] + (ed_2x[i-1,1] - 2*ed_2x[i,1] + ed_2x[i+1,1]) + (ed_2y[i,1] - 2*ed_2y[i,1] + ed_2y[i,2])
        temp_data_2[i,end] = density_2[i,end] + (ed_2x[i-1,end] - 2*ed_2x[i,end] + ed_2x[i+1,end]) + (ed_2y[i,end-1] - 2*ed_2y[i,end] + ed_2y[i,end])
    end
    # if any(isnan.(temp_data_1))
    #     println("NAN2")
    # end
    # Left and right boundaries
    for j = 2:length(par.ys)-1
        temp_data_1[1,j] = density_1[1,j] + (ed_1x[1,j] - 2*ed_1x[1,j] + ed_1x[2,j]) + (ed_1y[1,j-1] - 2*ed_1y[1,j] + ed_1y[1,j+1])
        temp_data_1[end,j] = density_1[end,j] + (ed_1x[end-1,j] - 2*ed_1x[end,j] + ed_1x[end,j]) + (ed_1y[end,j-1] - 2*ed_1y[end,j] + ed_1y[end,j+1])
    
        temp_data_2[1,j] = density_2[1,j] + (ed_2x[1,j] - 2*ed_2x[1,j] + ed_2x[2,j]) + (ed_2y[1,j-1] - 2*ed_2y[1,j] + ed_2y[1,j+1])
        temp_data_2[end,j] = density_2[end,j] + (ed_2x[end-1,j] - 2*ed_2x[end,j] + ed_2x[end,j]) + (ed_2y[end,j-1] - 2*ed_2y[end,j] + ed_2y[end,j+1])
    end
    # if any(isnan.(temp_data_1))
    #     println("NAN3")
    # end
    # Top-left corner
    temp_data_1[1,1] = density_1[1,1] + (ed_1x[1,1] - 2*ed_1x[1,1] + ed_1x[2,1]) + (ed_1y[1,1] - 2*ed_1y[1,1] + ed_1y[1,2])
    temp_data_2[1,1] = density_2[1,1] + (ed_2x[1,1] - 2*ed_2x[1,1] + ed_2x[2,1]) + (ed_2y[1,1] - 2*ed_2y[1,1] + ed_2y[1,2])
    # if any(isnan.(temp_data_1))
    #     println("NAN4")
    # end
    # Top-right corner
    temp_data_1[end,1] = density_1[end,1] + (ed_1x[end-1,1] - 2*ed_1x[end,1] + ed_1x[end,1]) + (ed_1y[end,1] - 2*ed_1y[end,1] + ed_1y[end,2])
    temp_data_2[end,1] = density_2[end,1] + (ed_2x[end-1,1] - 2*ed_2x[end,1] + ed_2x[end,1]) + (ed_2y[end,1] - 2*ed_2y[end,1] + ed_2y[end,2])
    # if any(isnan.(temp_data_1))
    #     println("NAN5")
    # end
    # Bottom-left corner
    temp_data_1[1,end] = density_1[1,end] + (ed_1x[1,end] - 2*ed_1x[1,end] + ed_1x[2,end]) + (ed_1y[1,end-1] - 2*ed_1y[1,end] + ed_1y[1,end])
    temp_data_2[1,end] = density_2[1,end] + (ed_2x[1,end] - 2*ed_2x[1,end] + ed_2x[2,end]) + (ed_2y[1,end-1] - 2*ed_2y[1,end] + ed_2y[1,end])
    # if any(isnan.(temp_data_1))
    #     println("NAN6")
    # end
    # Bottom-right corner
    temp_data_1[end,end] = density_1[end,end] + (ed_1x[end-1,end] - 2*ed_1x[end,end] + ed_1x[end,end]) + (ed_1y[end,end-1] - 2*ed_1y[end,end] + ed_1y[end,end])
    temp_data_2[end,end] = density_2[end,end] + (ed_2x[end-1,end] - 2*ed_2x[end,end] + ed_2x[end,end]) + (ed_2y[end,end-1] - 2*ed_2y[end,end] + ed_2y[end,end])
    # if any(isnan.(temp_data_1))
    #     println("NAN7")
    # end
    # Copy the updated values from temp_data back to the original density arrays
    density_1 .= temp_data_1
    density_2 .= temp_data_2
end

function apply_growth_and_mutation_den!(density_1, density_2, temp_data_1, temp_data_2, par::Parameters)
    #temp_data_1 .= par.growth_rate_1 .* density_1 .* (1 .- density_1 ./ par.carrying_capacity_1 .+ par.interaction_12 .* density_2 ./ par.carrying_capacity_2) .+ par.mutation_rate_12 .* density_2 .- par.mutation_rate_21 .* density_1
    #temp_data_2 .= par.growth_rate_2 .* density_2 .* (1 .- density_2 ./ par.carrying_capacity_2 .+ par.interaction_21 .* density_1 ./ par.carrying_capacity_1) .+ par.mutation_rate_21 .* density_1 .- par.mutation_rate_12 .* density_2

    temp_data_1 .= density_1 .* (par.growth_rate_1 .+ par.interaction_12 .* density_2 ./ par.carrying_capacity_2) .* (1 .- density_1 ./ par.carrying_capacity_1 .- density_2 ./ par.carrying_capacity_2)
    temp_data_2 .= density_2 .* (par.growth_rate_2 .+ par.interaction_21 .* density_1 ./ par.carrying_capacity_1) .* (1 .- density_1 ./ par.carrying_capacity_1 .- density_2 ./ par.carrying_capacity_2)

    density_1 .+= temp_data_1 .* par.dt
    density_2 .+= temp_data_2 .* par.dt
 end
 

 function apply_growth_and_mutation_frac!(density_1, density_2, temp_data_1, temp_data_2, par::Parameters)
    #temp_data_1 .= par.growth_rate_1 .* density_1 .* (1 .- density_1 ./ par.carrying_capacity_1 .+ par.interaction_12 .* density_2 ./ par.carrying_capacity_2) .+ par.mutation_rate_12 .* density_2 .- par.mutation_rate_21 .* density_1
    #temp_data_2 .= par.growth_rate_2 .* density_2 .* (1 .- density_2 ./ par.carrying_capacity_2 .+ par.interaction_21 .* density_1 ./ par.carrying_capacity_1) .+ par.mutation_rate_21 .* density_1 .- par.mutation_rate_12 .* density_2
    tot = density_1 .+ density_2
    tot[tot.==0] .= 1

    temp_data_1 .= density_1 .* (par.growth_rate_1 .+ par.interaction_12 .* density_2 ./ tot ./ par.carrying_capacity_2) .* (1 .- density_1 ./ par.carrying_capacity_1 .- density_2 ./ par.carrying_capacity_2)
    temp_data_2 .= density_2 .* (par.growth_rate_2 .+ par.interaction_21 .* density_1 ./ tot ./ par.carrying_capacity_1) .* (1 .- density_1 ./ par.carrying_capacity_1 .- density_2 ./ par.carrying_capacity_2)

    density_1 .+= temp_data_1 .* par.dt
    density_2 .+= temp_data_2 .* par.dt
 end
 

 function apply_demographic_noise!(density_1, density_2, par::Parameters)
     for i = eachindex(density_1)
        if(density_1[i] > 0)
            density_1[i] = par.dt * par.sigma^2 / par.dx * rand(Poisson(density_1[i] * par.dx / par.dt / par.sigma^2))
        end
        if(density_2[i] > 0)
            density_2[i] = par.dt * par.sigma^2 / par.dx * rand(Poisson(density_2[i] * par.dx / par.dt / par.sigma^2))
        end
     end
 end

 function h_xt(result, par::Parameters, thres=0.99)
    Lx, Ly, nt, dumi = size(result)
    Nt = result[:,:,:,1]/par.carrying_capacity_1 .+ result[:,:,:,2]/par.carrying_capacity_2
    # Find the indices of the maximal y value exceeding the threshold for each x
    h = fill(0, Lx, nt)
    for t in 1:nt
        h[:,t] .= [let idx = findlast(y -> y > thres, Nt[x, :, t])
                    idx === nothing ? 1 : idx
                end for x in 1:Lx]
    end
    return h
 end

 function f_xt(result, h)
    L, L, nt, dumi = size(result)
    f = fill(0.0, L, nt)
    f[1:div(L,2),1] .= 1.; f[div(L,2):end,1] .= 0.; 
    for t = 2:nt
        for x = 1:L
            if(h[x,t] <= 1)
                f[x,t] = f[x, t-1]
            else
                if( (result[x, h[x,t], t, 1] + result[x, h[x,t], t, 2]) > 0 )
                    f[x,t] = result[x, h[x,t], t, 1]/(result[x, h[x,t], t, 1] + result[x, h[x,t], t, 2])
                end
            end
        end
    end
    return f
end

end