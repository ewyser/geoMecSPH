export setup_mpts

"""
    setup_mpts(meD, ni, lz, coh0, cohr, phi0, phir, rho0, nstr, typeD) -> Point

Initialize material point (particle) system.

# Arguments
- `meD`: Mesh structure
- `ni`: Number of material points per element direction
- `lz`: Domain height
- `coh0`: Initial cohesion
- `cohr`: Residual cohesion
- `phi0`: Initial friction angle
- `phir`: Residual friction angle
- `rho0`: Initial density
- `nstr`: Number of stress components
- `typeD`: Data type

# Returns
- `mpD`: Material point structure
"""
function setup_mpts(meD, ni, lz, coh0, cohr, phi0, phir, rho0, nstr, typeD)
    # Material point initialization
    xL = meD.xB[1]+(0.5*meD.h[1]/ni):meD.h[1]/ni:meD.xB[2]
    zL = meD.xB[3]+(0.5*meD.h[2]/ni):meD.h[2]/ni:lz-0.5*meD.h[2]/ni
    npx = length(xL)
    npz = length(zL)
    
    xp = vec(xL' .* ones(npz, 1))
    zp = vec(ones(npx, 1)' .* zL)
    
    # Water level
    wl = 0.15 * lz
    x = LinRange(minimum(xp), maximum(xp), 200)
    a = -1.25
    z = a .* x
    x = x .+ 0.5 .* meD.L[1]
    
    xlt = Float64[]
    zlt = Float64[]
    vxlt = Float64[]
    vzlt = Float64[]
    clt = Float64[]
    
    # Filter material points based on geometry
    for mp in 1:length(xp)
        xc0, zc0 = 0.4 * (meD.xB[2] - meD.xB[1]), 0.4 * lz
        xc1, zc1 = 0.6 * (meD.xB[2] - meD.xB[1]), 0.6 * lz
        dx, dz = xc1 - xc0, zc1 - zc0
        
        vx = dx / sqrt(dx * dx + dz * dz)
        vz = dz / sqrt(dx * dx + dz * dz)
        v = 20.5
        
        r = 0.1 * lz
        Δx, Δz = xp[mp] - xc0, zp[mp] - zc0
        d = sqrt(Δx * Δx + Δz * Δz)
        
        if d <= r
            push!(xlt, xp[mp])
            push!(zlt, zp[mp])
            push!(vxlt, v * vx)
            push!(vzlt, v * vz)
        end
        
        Δx, Δz = xp[mp] - xc1, zp[mp] - zc1
        d = sqrt(Δx * Δx + Δz * Δz)
        
        if d <= r
            push!(xlt, xp[mp])
            push!(zlt, zp[mp])
            push!(vxlt, -v * vx)
            push!(vzlt, -v * vz)
        end
    end
    
    # Material point quantities
    nmp = length(xlt)
    l0 = ones(typeD, nmp, 2) .* 0.5 .* (meD.h[1] ./ ni)
    v0 = prod(2 .* l0, dims=2)
    h0 = l0
    m = vec(rho0 .* v0)
    xp = hcat(xlt, zlt)
    up = zeros(typeD, nmp, 2)
    vp = hcat(vxlt, vzlt)
    
    # Material properties
    coh = ones(typeD, nmp, 1) .* coh0
    cohr = ones(typeD, nmp, 1) .* cohr
    phi = ones(typeD, nmp, 1) .* phi0
    epII = zeros(typeD, nmp, 1)
    J = ones(typeD, nmp, 1)
    
    # Tensors
    dF = zeros(typeD, nmp, 4)
    F = repeat([1 0 0 1], nmp, 1)
    b = zeros(typeD, nmp, 4)
    bT = zeros(typeD, nmp, 4)
    e = zeros(typeD, nstr, nmp)
    ome = zeros(typeD, 1, nmp)
    s = zeros(typeD, nstr, nmp)
    τ = zeros(typeD, nstr, nmp)
    dev = zeros(typeD, nstr, nmp)
    ep = zeros(typeD, nstr, nmp)
    
    # Additional quantities
    nn = convert(UInt64, meD.nn)
    ϕ = zeros(typeD, nmp, nmp)
    ∂ϕx = zeros(typeD, nmp, nmp)
    ∂ϕz = zeros(typeD, nmp, nmp)
    B = zeros(typeD, nstr, nn .* 2, nmp)
    
    # Connectivity
    p2e = zeros(UInt64, nmp, 1)
    p2n = zeros(UInt64, nmp, nn)
    
    # Create point structure
    mpD = Point(nmp, h0, v0, m, xp, up, vp, coh, cohr, phi, epII, J,
                dF, F, b, bT, e, ome, s, τ, dev, ep, ϕ, ∂ϕx, ∂ϕz, B, p2e, p2n)
    
    return mpD
end
