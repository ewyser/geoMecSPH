# Boot script for geoMecSPH
# This file handles initial package setup and file inclusion

# Include dependencies
using LinearAlgebra
using Plots
using LaTeXStrings
using Base.Threads
using ProgressMeter

# Include utilities and type system
include(joinpath(SRC, "boot/include.jl"))
include(joinpath(SRC, "boot/needs/utils.jl"))

# Include type definitions
sucess = superInc(["boot/needs/types"]; root=SRC)

# Include main functionality
lists = ["home/init", "home/api", "home/core", "home/script"]
@info join(superInc(lists; root=SRC), "\n")
