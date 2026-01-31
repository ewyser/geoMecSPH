# Refactoring Summary: geoMecSPH → ElastoPlasm.jl Architecture

## Overview

The geoMecSPH project has been refactored to follow the modern, well-organized structure of ElastoPlasm.jl. This document summarizes the changes and migration path.

## New Directory Structure

```
src/
├── geoMecSPH.jl                      # Main module (refactored)
│
├── boot/                              # NEW: Bootstrap system
│   ├── boot.jl                       # Orchestrates package initialization
│   ├── include.jl                    # File inclusion utilities (superInc, superDir)
│   └── needs/                        # Core dependencies
│       ├── types/
│       │   └── types.jl              # Type definitions (Mesh, Point, Boundary)
│       └── utils.jl                  # Utility functions (welcome_log, exit_log)
│
└── home/                              # NEW: Main functionality
    ├── init/                          # NEW: Initialization routines
    │   ├── mesh/
    │   │   └── setup_mesh.jl         # Mesh initialization (from misc.jl)
    │   └── mpts/
    │       ├── setup_mpts.jl         # Material points (from misc.jl)
    │       └── RFS.jl                # Random fields (from RFS.jl)
    │
    ├── api/                           # NEW: Public API
    │   ├── materials.jl              # Material parameters (from misc.jl D())
    │   └── paths.jl                  # Path management
    │
    ├── core/                          # NEW: Core algorithms
    │   ├── workflow/                 # Simulation workflows (TO BE MIGRATED)
    │   ├── constitutive/             # Constitutive models (TO BE MIGRATED from fun_fs/)
    │   └── mapping/                  # Particle-grid mapping (TO BE MIGRATED from fun_fs/)
    │
    └── script/                        # NEW: Example scripts
        └── example/
            ├── sim.jl                # Main simulation driver (refactored)
            └── slump.jl              # Slump test (refactored)
```

## File Migration Map

### Completed Migrations

| Old Location | New Location | Status |
|-------------|--------------|---------|
| `src/fun/types.jl` | `src/boot/needs/types/types.jl` | ✅ Migrated |
| `src/fun/misc.jl` (meshSetup) | `src/home/init/mesh/setup_mesh.jl` | ✅ Migrated |
| `src/fun/misc.jl` (pointSetup) | `src/home/init/mpts/setup_mpts.jl` | ✅ Migrated |
| `src/fun/misc.jl` (D function) | `src/home/api/materials.jl` | ✅ Migrated |
| `src/fun/RFS.jl` | `src/home/init/mpts/RFS.jl` | ✅ Migrated |
| `script/sim.jl` | `src/home/script/example/sim.jl` | ✅ Refactored |
| `src/_superORCH.jl` | `src/boot/boot.jl` | ✅ Replaced |
| `src/geoMecSPH.jl` | Updated with new architecture | ✅ Updated |

### Pending Migrations

| Old Location | New Location | Status |
|-------------|--------------|---------|
| `src/fun_fs/elast.jl` | `src/home/core/constitutive/elast.jl` | 🔄 To be migrated |
| `src/fun_fs/plast.jl` | `src/home/core/constitutive/plast.jl` | 🔄 To be migrated |
| `src/fun_fs/solve.jl` | `src/home/core/workflow/solve.jl` | 🔄 To be migrated |
| `src/fun_fs/topol.jl` | `src/home/core/mapping/topology.jl` | 🔄 To be migrated |
| `src/fun_fs/flip.jl` | `src/home/core/mapping/flip.jl` | 🔄 To be migrated |
| `src/fun_fs/accum.jl` | `src/home/core/mapping/accumulate.jl` | 🔄 To be migrated |
| `src/fun_fs/BSpline.jl` | `src/home/core/mapping/bspline.jl` | 🔄 To be migrated |
| `src/fun_fs/GIMP.jl` | `src/home/core/mapping/gimp.jl` | 🔄 To be migrated |
| `src/fun_fs/DMBC.jl` | `src/home/core/workflow/boundary.jl` | 🔄 To be migrated |

## Key Changes

### 1. Module Structure

**Old:**
```julia
module geoMecSPH
include("_superORCH.jl")
global path_plot = "./out/"
sim()
include(joinpath("../script","sim.jl"))
end
```

**New:**
```julia
module geoMecSPH
const SRC = @__DIR__
export sim, ic_slump
include(joinpath(SRC, "boot/boot.jl"))
function __init__()
    welcome_log()
end
end
```

### 2. Type Definitions

**Old:** Split between `src/fun/types.jl` and `src/fun_fs/types.jl`

**New:** Consolidated in `src/boot/needs/types/types.jl`

### 3. Initialization

**Old:** Functions scattered in `src/fun/misc.jl`

**New:** Organized in dedicated files:
- `src/home/init/mesh/setup_mesh.jl`
- `src/home/init/mpts/setup_mpts.jl`
- `src/home/init/mpts/RFS.jl`

### 4. API Functions

**Old:** Mixed with implementation in `src/fun/misc.jl`

**New:** Clear API in `src/home/api/`:
- `materials.jl` - Material parameters
- `paths.jl` - File system management

### 5. Bootstrap System

**New Addition:** `src/boot/` directory with:
- Automatic file inclusion via `superInc()`
- Clean initialization sequence
- Welcome messages

## Benefits of Refactoring

1. **Clear Organization**: Files grouped by functionality, not arbitrary labels
2. **Discoverability**: Easy to find where things are
3. **Maintainability**: Changes localized to specific modules
4. **Extensibility**: Easy to add new features in appropriate locations
5. **Testing**: Can test isolated components
6. **Documentation**: Clear boundaries between public API and internals

## Migration Guide for Developers

### Using the New Structure

**Initialize simulation:**
```julia
using geoMecSPH

# Quick start
sim()

# Or custom setup
setup = ic_slump(nel=200, lx=64.1584, lz=12.80)
```

### Adding New Features

1. **New constitutive model**: Add to `src/home/core/constitutive/`
2. **New mapping scheme**: Add to `src/home/core/mapping/`
3. **New example**: Add to `src/home/script/example/`
4. **New type**: Add to `src/boot/needs/types/`

### Old Code Compatibility

The old structure is preserved in:
- `src/fun/` (old utilities)
- `src/fun_fs/` (old finite strain code)
- `script/sim.jl` (old simulation script)

These can be removed once remaining migrations are complete.

## Next Steps

1. **Migrate core algorithms** from `src/fun_fs/` to `src/home/core/`
2. **Create workflow functions** for elastodynamic and elastoplastic solvers
3. **Add comprehensive examples** in `src/home/script/example/`
4. **Write tests** for each module
5. **Update documentation** with new API
6. **Remove old directories** once migration is complete

## References

- [ElastoPlasm.jl](https://github.com/ewyser/ElastoPlasm.jl) - Inspiration for structure
- Material Point Method literature
- Julia package development best practices
