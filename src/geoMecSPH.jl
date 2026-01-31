module geoMecSPH

# Define module location as constant
const SRC = @__DIR__

# Export main simulation functions
export sim, ic_slump
# Export plotting functions
export plot_ρ, plot_coh, plot_Δu, plot_P, plot_Δϵp
# Export utility functions
export get_Δt, get_g

# Include boot file
include(joinpath(SRC, "boot/boot.jl"))

@doc raw"""
    geoMecSPH

A geomechanical code within an explicit SPH/MPM numerical framework.

## Features
- Elastoplastic material models
- Updated Lagrangian formulation
- Multiple shape function bases (GIMP, B-spline)
- PIC-FLIP mapping schemes

## Quick Start
```julia
using geoMecSPH

# Run default simulation
sim()

# Or create custom initial conditions
setup = ic_slump(nel=200, lx=64.1584, lz=12.80)
```

## Reference
Based on the Material Point Method (MPM) and Smoothed Particle Hydrodynamics (SPH).
Refactored following ElastoPlasm.jl architecture.

## Usage
Run with: `julia -O3 --threads=auto --check-bounds=no --project=.`
"""
geoMecSPH

function __init__()
    welcome_log()
end

end # module geoMecSPH