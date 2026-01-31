# Export type definitions
export Mesh, Point, Boundary

"""
    Mesh

Mutable struct representing the computational mesh/grid.
"""
mutable struct Mesh
    nel::Array{Float64}   # number of elements
    nno::Array{Int64}     # number of nodes
    nn ::Int64            # nodes per element
    L  ::Array{Float64}   # domain dimensions
    h  ::Array{Float64}   # element sizes
    xn ::Array{Float64}   # node x-coordinates
    zn ::Array{Float64}   # node z-coordinates
    e2n::Array{Int64}     # element to node connectivity
    xB ::Array{Float64}   # boundary coordinates
end

"""
    Point

Mutable struct representing material points (particles).
"""
mutable struct Point
    # scalars & vectors
    nmp ::Int64              # number of material points
    h   ::Array{Float64}     # characteristic lengths
    v0  ::Array{Float64}     # initial volumes
    m   ::Array{Float64}     # masses
    x   ::Array{Float64}     # positions
    u   ::Array{Float64}     # displacements
    v   ::Array{Float64}     # velocities
    coh ::Array{Float64}     # cohesion
    cohr::Array{Float64}     # residual cohesion
    phi ::Array{Float64}     # friction angles
    epII::Array{Float64}     # equivalent plastic strain
    J   ::Array{Float64}     # Jacobian determinant
    # tensors
    ΔF  ::Array{Float64}     # incremental deformation gradient
    F   ::Array{Float64}     # deformation gradient
    b   ::Array{Float64}     # left Cauchy-Green tensor
    bT  ::Array{Float64}     # transpose of b
    ϵ   ::Array{Float64}     # strain tensor
    ω   ::Array{Float64}     # spin tensor
    σ   ::Array{Float64}     # stress tensor
    τ   ::Array{Float64}     # Kirchhoff stress
    dev ::Array{Float64}     # deviatoric stress
    ep  ::Array{Float64}     # plastic strain
    # additional quantities
    ϕ   ::Array{Float64}     # shape functions
    ∂ϕx ::Array{Float64}     # shape function x-derivatives
    ∂ϕz ::Array{Float64}     # shape function z-derivatives
    B   ::Array{Float64}     # B-matrix (strain-displacement)
    # connectivity
    p2e ::Array{Int64}       # particle to element
    p2n ::Array{Int64}       # particle to node
end

"""
    Boundary

Mutable struct for boundary conditions.
"""
mutable struct Boundary
    x  ::Array{Int64}        # x-direction boundary nodes
    z  ::Array{Int64}        # z-direction boundary nodes
end
