# Migration Checklist for geoMecSPH Refactoring

## ✅ Completed Tasks

- [x] Created new directory structure following ElastoPlasm.jl
- [x] Implemented bootstrap system (`boot/` directory)
- [x] Consolidated type definitions in `boot/needs/types/`
- [x] Migrated mesh initialization to `home/init/mesh/`
- [x] Migrated material point initialization to `home/init/mpts/`
- [x] Migrated random field generation to `home/init/mpts/RFS.jl`
- [x] Created API module for materials and paths
- [x] Refactored main module file (`geoMecSPH.jl`)
- [x] Created example scripts (`home/script/example/`)
- [x] Added comprehensive documentation (README, REFACTORING.md)
- [x] Added utility functions and logging

## 🔄 Remaining Tasks

### High Priority

- [ ] Migrate `fun_fs/elast.jl` to `home/core/constitutive/elast.jl`
  - Update function signatures to use new types
  - Add proper exports
  - Add documentation

- [ ] Migrate `fun_fs/plast.jl` to `home/core/constitutive/plast.jl`
  - Multiple plasticity implementations to organize
  - Document different return mapping algorithms

- [ ] Migrate `fun_fs/solve.jl` to `home/core/workflow/solve.jl`
  - Main time integration loop
  - Critical for simulation execution

### Medium Priority

- [ ] Migrate shape function implementations:
  - `fun_fs/BSpline.jl` → `home/core/mapping/bspline.jl`
  - `fun_fs/GIMP.jl` → `home/core/mapping/gimp.jl`

- [ ] Migrate mapping schemes:
  - `fun_fs/flip.jl` → `home/core/mapping/flip.jl`
  - `fun_fs/accum.jl` → `home/core/mapping/accumulate.jl`
  - `fun_fs/topol.jl` → `home/core/mapping/topology.jl`

- [ ] Migrate boundary conditions:
  - `fun_fs/DMBC.jl` → `home/core/workflow/boundary.jl`

### Low Priority

- [ ] Create workflow orchestrator:
  - `home/core/workflow/elastodynamic.jl`
  - `home/core/workflow/elastoplastic.jl`
  - Pattern: `elastoplasm!(setup; workflow=[elastodynamic!, elastoplastic!])`

- [ ] Add plotting utilities:
  - `home/api/plotting.jl`
  - Visualization functions

- [ ] Create additional examples:
  - Column collapse
  - Dam break
  - Custom geometries

### Testing & Documentation

- [ ] Write unit tests for:
  - Type constructors
  - Mesh initialization
  - Material point setup
  - Random field generation
  - Constitutive models

- [ ] Add inline documentation:
  - Docstrings for all exported functions
  - Usage examples

- [ ] Create tutorials:
  - Basic simulation setup
  - Custom material models
  - Post-processing

### Cleanup

- [ ] Remove old directories once migration complete:
  - `src/fun/` (after verifying all functions migrated)
  - `src/fun_fs/` (after verifying all functions migrated)
  - `script/sim.jl` (superseded by new structure)
  - `src/_superORCH.jl` (replaced by `boot/boot.jl`)

- [ ] Update `Project.toml`:
  - Add version number
  - Update description
  - Add keywords

- [ ] Code quality:
  - Run formatter on all new files
  - Check for unused imports
  - Verify all exports are intentional

## 📋 Testing Checklist

Before removing old code:

- [ ] Verify `sim()` function works with new structure
- [ ] Test mesh generation with different parameters
- [ ] Test material point initialization
- [ ] Verify random field generation produces expected results
- [ ] Check that all types are properly exported
- [ ] Ensure backward compatibility where needed

## 🎯 Success Criteria

The refactoring is complete when:

1. All functionality from `fun/` and `fun_fs/` is migrated
2. All tests pass with new structure
3. Documentation is comprehensive
4. Example scripts demonstrate key features
5. Old directories can be safely removed
6. Package can be loaded with `using geoMecSPH`
7. Main simulation runs without errors

## 📝 Notes

- Keep old code until fully tested
- Maintain git history of migrations
- Document any breaking changes
- Consider semantic versioning for release
- Update citations and references

## 🔗 References

- ElastoPlasm.jl: https://github.com/ewyser/ElastoPlasm.jl
- Julia package structure: https://pkgdocs.julialang.org/
- Material Point Method resources in `refs/`
