# geoMecSPH.jl 👻

A geomechanical code within an explicit SPH/MPM numerical framework, refactored following [ElastoPlasm.jl](https://github.com/ewyser/ElastoPlasm.jl) architecture.

## ✨ Overview

**geoMecSPH.jl** provides a modern Julia implementation for elastoplastic simulations using:

- **Updated Lagrangian formulation** for large deformation problems
- **Multiple shape function bases**: Linear, GIMP, B-spline
- **PIC-FLIP mapping schemes** for particle-grid transfer
- **Elastoplastic constitutive models**: Drucker-Prager, Mohr-Coulomb
- **Gaussian Random Fields** for heterogeneous material properties

## 🏗️ Project Structure

The project follows a clean, modular architecture inspired by ElastoPlasm.jl:

```
geoMecSPH/
├── src/
│   ├── geoMecSPH.jl          # Main module file
│   ├── boot/                  # Bootstrap & initialization
│   │   ├── boot.jl            # Main boot orchestrator
│   │   ├── include.jl         # File inclusion utilities
│   │   └── needs/             # Core dependencies
│   │       ├── types/         # Type definitions
│   │       │   └── types.jl   # Mesh, Point, Boundary structs
│   │       └── utils.jl       # Utility functions
│   │
│   └── home/                  # Main functionality
│       ├── init/              # Initialization routines
│       │   ├── mesh/          # Mesh setup
│       │   │   └── setup_mesh.jl
│       │   └── mpts/          # Material points setup
│       │       ├── setup_mpts.jl
│       │       └── RFS.jl     # Random field generation
│       │
│       ├── api/               # Public API & configuration
│       │   ├── materials.jl   # Material parameters
│       │   └── paths.jl       # Path management
│       │
│       ├── core/              # Core algorithms
│       │   ├── workflow/      # Main simulation workflows
│       │   ├── constitutive/  # Constitutive models
│       │   └── mapping/       # Particle-grid mapping
│       │
│       └── script/            # Example scripts
│           └── example/
│               ├── sim.jl     # Main simulation driver
│               └── slump.jl   # Slump test setup
│
├── script/                    # Legacy scripts (to be migrated)
├── refs/                      # Reference materials
├── Project.toml              # Package dependencies
├── Manifest.toml             # Dependency lock file
└── README.md                 # This file
```

## 🛠️ Installation

1. **Install Julia**: [Download here](https://julialang.org/downloads/)

2. **Clone the repository**:
   ```bash
   git clone <repository-url>
   cd geoMecSPH
   ```

3. **Activate and instantiate**:
   ```bash
   julia --project=.
   ```
   ```julia
   julia> using Pkg
   julia> Pkg.instantiate()
   ```

4. **Load the package**:
   ```julia
   julia> using geoMecSPH
   ```

## 🚀 Quick Start

### Basic Usage

```julia
using geoMecSPH

# Run default slump simulation
sim()
```

### Custom Simulation

```julia
using geoMecSPH

# Create initial conditions
setup = ic_slump(
    nel = 200,           # Number of elements
    lx = 64.1584,       # Domain length [m]
    lz = 12.80          # Domain height [m]
)

# Access components
mesh = setup.mesh      # Mesh structure
mpts = setup.mpts      # Material points
bc = setup.bc          # Boundary conditions
params = setup.params  # Physical parameters
```

## 📊 Key Components

### Type Definitions

- **`Mesh`**: Computational mesh/grid structure
- **`Point`**: Material point (particle) structure
- **`Boundary`**: Boundary condition structure

### Initialization Functions

- `setup_mesh(nel, lx, lz, typeD)`: Initialize computational mesh
- `setup_mpts(...)`: Initialize material point system
- `RFS(...)`: Generate random fields for heterogeneous properties

### API Functions

- `elasticity_matrix(E, ν)`: Compute elastic stiffness matrix
- `setup_paths(path)`: Create output directories
- `ic_slump(...)`: Initialize slump test problem

## 🔬 Refactoring Benefits

This refactored structure provides:

1. **Modularity**: Clear separation of concerns
2. **Maintainability**: Easy to locate and update functionality
3. **Extensibility**: Simple to add new features
4. **Discoverability**: Organized file structure matches functionality
5. **Testing**: Easier to write unit tests for isolated components
6. **Documentation**: Clear API boundaries

## 📚 Architecture Comparison

### Old Structure
```
src/
├── geoMecSPH.jl
├── _superORCH.jl
├── fun/
│   ├── types.jl
│   ├── misc.jl
│   └── RFS.jl
└── fun_fs/
    ├── types.jl
    ├── elast.jl
    └── plast.jl
```

### New Structure (ElastoPlasm-inspired)
```
src/
├── geoMecSPH.jl
├── boot/               # Bootstrap system
│   ├── boot.jl
│   └── needs/
└── home/               # Main functionality
    ├── init/           # Initialization
    ├── api/            # Public API
    ├── core/           # Core algorithms
    └── script/         # Examples
```

## 🎯 Migration Notes

The refactoring organizes code by **functionality** rather than arbitrary groupings:

- **`fun/` → `home/init/` + `home/api/`**: Initialization and configuration
- **`fun_fs/` → `home/core/`**: Core computational algorithms
- **`script/` → `home/script/example/`**: Example scripts
- **Types consolidated** in `boot/needs/types/`

## 🔧 Development

Run with optimization flags:
```bash
julia -O3 --threads=auto --check-bounds=no --project=.
```

## 📝 References

- Inspired by [ElastoPlasm.jl](https://github.com/ewyser/ElastoPlasm.jl)
- Material Point Method (MPM)
- Smoothed Particle Hydrodynamics (SPH)

## 👥 Contributing

Contributions are welcome! The new modular structure makes it easier to:
- Add new constitutive models in `home/core/constitutive/`
- Implement new mapping schemes in `home/core/mapping/`
- Create example problems in `home/script/example/`

## 📄 License

See LICENSE file for details.
