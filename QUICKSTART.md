# Quick Start Guide - geoMecSPH.jl

## Getting Started with the Refactored Structure

### 1. Installation

```bash
cd /path/to/geoMecSPH
julia --project=.
```

In Julia REPL:
```julia
julia> using Pkg
julia> Pkg.instantiate()
julia> using geoMecSPH
```

You should see:
```
┌ Welcome to geoMecSPH 👻
│ Geomechanical SPH/MPM solver
│ New comer ? Try:
│   using geoMecSPH
└   sim()
```

### 2. Run Your First Simulation

```julia
# Simple one-liner
sim()
```

### 3. Custom Simulation Setup

```julia
using geoMecSPH

# Create initial conditions for slump test
setup = ic_slump(
    nel = 200,          # Number of elements
    lx = 64.1584,      # Domain length [m]
    lz = 12.80         # Domain height [m]
)

# Access components
mesh = setup.mesh      # Computational mesh
mpts = setup.mpts      # Material points (particles)
bc = setup.bc          # Boundary conditions
params = setup.params  # Physical parameters
```

### 4. Understanding the Setup

```julia
# Mesh information
println("Number of nodes: ", mesh.nno[end])
println("Number of elements: ", mesh.nel[3])
println("Element size: ", mesh.h)

# Material point information
println("Number of particles: ", mpts.nmp)
println("Particle positions: ", size(mpts.x))
println("Particle velocities: ", size(mpts.v))

# Material parameters
println("Young's modulus: ", params.E, " Pa")
println("Density: ", params.ρ0, " kg/m³")
println("Cohesion: ", params.c0, " Pa")
println("Friction angle: ", params.ϕ0 * 180/π, " degrees")
```

### 5. Project Structure Overview

```
src/
├── boot/                    # Bootstrap (loads everything)
│   └── needs/types/        # Type definitions
│
└── home/
    ├── init/               # Setup functions
    │   ├── mesh/          # setup_mesh()
    │   └── mpts/          # setup_mpts(), RFS()
    │
    ├── api/               # Public functions
    │   ├── materials.jl   # elasticity_matrix()
    │   └── paths.jl       # setup_paths()
    │
    ├── core/              # Core algorithms (TO BE POPULATED)
    │   ├── workflow/     # Main solvers
    │   ├── constitutive/ # Material models
    │   └── mapping/      # Particle-grid transfer
    │
    └── script/example/   # Example scripts
        ├── sim.jl        # Main driver
        └── slump.jl      # Slump test
```

### 6. Next Steps

The refactoring provides a clean foundation. To complete the migration:

1. **Migrate core algorithms** from `fun_fs/` to `home/core/`
2. **Add workflow functions** for time integration
3. **Create more examples** in `home/script/example/`

### 7. Key Functions Reference

#### Initialization
- `ic_slump(; nel, lx, lz, kwargs...)` - Setup slump test
- `setup_mesh(nel, lx, lz, typeD)` - Initialize mesh
- `setup_mpts(...)` - Initialize material points
- `RFS(xp, zp, coh0, cohr, phi0, phir)` - Random fields

#### API
- `elasticity_matrix(E, ν)` - Compute elastic stiffness
- `setup_paths(path)` - Create output directories

#### Types
- `Mesh` - Computational grid
- `Point` - Material points/particles  
- `Boundary` - Boundary conditions

### 8. Comparison: Old vs New

#### Old Way
```julia
# Scattered across multiple files
include("_superORCH.jl")
include(joinpath("fun", "types.jl"))
include(joinpath("fun", "misc.jl"))
# Functions had unclear organization
```

#### New Way
```julia
# Clean, organized imports
using geoMecSPH

# Everything is automatically loaded
# Clear API boundaries
# Easy to find functionality
```

### 9. Development Workflow

```julia
# 1. Make changes to source files

# 2. Reload (if using Revise.jl)
julia> using Revise
julia> using geoMecSPH

# 3. Test changes
julia> sim()

# 4. Run with optimization for production
# From terminal:
# julia -O3 --threads=auto --check-bounds=no --project=.
```

### 10. Troubleshooting

**Issue: Package not loading**
```julia
# Solution: Reinstantiate
using Pkg
Pkg.instantiate()
```

**Issue: Can't find functions**
```julia
# Check what's exported
using geoMecSPH
names(geoMecSPH)
```

**Issue: Need to see structure**
```bash
# Run from terminal
bash STRUCTURE.sh
```

### 11. Documentation

- **README_NEW.md** - Overview and features
- **REFACTORING.md** - Detailed migration notes
- **MIGRATION_CHECKLIST.md** - Remaining tasks
- **STRUCTURE.sh** - Visual directory tree

### 12. Contributing

To add new features:

1. **New constitutive model**: Add to `src/home/core/constitutive/`
2. **New example**: Add to `src/home/script/example/`
3. **New type**: Add to `src/boot/needs/types/`

### 13. Performance Tips

```julia
# Run with multiple threads
julia> Threads.nthreads()  # Check available threads

# From terminal with optimization:
julia -O3 --threads=auto --check-bounds=no --project=.
```

### 14. Further Reading

- ElastoPlasm.jl: https://github.com/ewyser/ElastoPlasm.jl
- Material Point Method references in `refs/`
- Julia package development: https://pkgdocs.julialang.org/

---

**Happy simulating! 🎉**

For questions or issues, check:
- REFACTORING.md for migration details
- MIGRATION_CHECKLIST.md for remaining work
- Original papers in `refs/` directory
