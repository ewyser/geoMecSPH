module geoMecSPH

export sim # host code


#include dependencies
include("_superORCH.jl")

# relative path for figs & data
global path_plot = "./out/"
if isdir(path_plot)==false
    mkdir(path_plot)    
end

#include main
# include geoflow routine in saintVenant module
@doc raw"""
geoflow(lx::Float64,ly::Float64,nx::Int64,T::Float64,tC::Float64,rheoType::String,solvType::String,isViz::Bool): solves a non-linear hyperbolic 2D Saint-Venant problem considering a Coulomb-type rheology within a finite volume framework on a Cartesian grid
# args:
- lx       : dimension along the x-direciton in [m].
- ly       : dimension along the y-direciton in [m].
- nx       : number of grid nodes along the x-direction  in [-].
- T        : total simulation time in [s].
- tC       : time interval for saving/plotting
- rheoType : select the rheology, i.e., "coulomb", "newtonian" or "plastic"
- solveType: select the numerical flux, i.e., "Rusanov", "HLL" or "HLLC"
- isViz    : plot or save, true or false
To run geoflow() on a GPU, add *_D, i.e., geoflow_D(lx::Float64,ly::Float64,nx::Int64,T::Float64,tC::Float64,rheoType::String,solvType::String,isViz::Bool)
"""
sim()
include(joinpath("../script","sim.jl"))

end # module geoMecSPH

# julia -O3 --threads=auto --check-bounds=no --project=.