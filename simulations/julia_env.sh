# Makes `julia` available, on the login node and on compute nodes.
# Edit the module name if `module avail julia` shows a specific version (e.g. julia/1.10.4).
module load julia 2>/dev/null || true
if ! command -v julia >/dev/null 2>&1 && [ -x "$HOME/.juliaup/bin/julia" ]; then
    export PATH="$HOME/.juliaup/bin:$PATH"       # fallback: user install via juliaup
fi
if ! command -v julia >/dev/null 2>&1; then
    echo "julia not found: run 'module avail julia', or install with" >&2
    echo "  curl -fsSL https://install.julialang.org | sh" >&2
    return 1 2>/dev/null || exit 1
fi
