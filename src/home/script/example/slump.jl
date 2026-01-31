export ic_slump

"""
    ic_slump(; kwargs...) -> String

Create initial conditions for a slump simulation problem.

# Keyword Arguments
- `nel`: Number of elements (default: 200)
- `lx`: Domain length in x-direction (default: 64.1584)
- `lz`: Domain length in z-direction (default: 12.80)
- Additional kwargs for simulation parameters

# Returns
- Path to simulation setup file

# Example
```julia
sim_file = ic_slump(nel=200, lx=64.1584, lz=12.80)
```
"""
function ic_slump(; nel=200, lx=64.1584, lz=12.80, kwargs...)
    @info "Setting up mesh & material point system for slump problem"
    
    # Simulation parameters
    typeD = Float64
    ν = 0.3          # Poisson's ratio
    ni = 3           # Material points per element direction
    nstr = 4         # Number of stress components
    
    # Physical constants
    g = 9.81         # Gravitational acceleration [m/s^2]
    E = 1.0e6        # Young's modulus [Pa]
    ρ0 = 2700.0      # Density [kg/m^3]
    c0 = 20.0e3      # Cohesion [Pa]
    ϕ0 = 20.0 * π / 180  # Friction angle [rad]
    cr = 4.0e3       # Residual cohesion [Pa]
    ϕr = 7.5 * π / 180   # Residual friction angle [rad]
    
    # Derived parameters
    K, G, Del = elasticity_matrix(E, ν)
    
    # Setup mesh and material points
    meD, bc = setup_mesh(nel, lx, lz, typeD)
    mpD = setup_mpts(meD, ni, lz, c0, cr, ϕ0, ϕr, ρ0, nstr, typeD)
    
    @info "Mesh & material points initialized" mpD.nmp meD.nno[end]
    
    return (mesh=meD, mpts=mpD, bc=bc, params=(;E=E, ν=ν, ρ0=ρ0, g=g, c0=c0, ϕ0=ϕ0))
end
