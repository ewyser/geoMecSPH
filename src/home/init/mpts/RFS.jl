export RFS, GRF

"""
    RFS(xp, zp, coh0, cohr, phi0, phir) -> (coh, phi)

Generate Random Field using Gaussian Random Field (GRF) method.

# Arguments
- `xp`: x-coordinates of material points
- `zp`: z-coordinates of material points
- `coh0`: Mean cohesion
- `cohr`: Residual cohesion (unused in current implementation)
- `phi0`: Mean friction angle
- `phir`: Residual friction angle (unused in current implementation)

# Returns
- `coh`: Cohesion field
- `phi`: Friction angle field
"""
function RFS(xp, zp, coh0, cohr, phi0, phir)
    # Parameters
    θx = 20.0
    θz = 2.0
    β = 45.0 * pi / 180
    μc = coh0
    σc = μc / 5.0
    μϕ = phi0
    σϕ = μϕ / 10.0
    
    # Vector format
    xp = vec(xp)
    zp = vec(zp)
    nmp = length(xp)
    
    # Relative distance
    Δx = (xp .- xp')
    Δz = zp .- zp'
    
    # Exponential covariance matrix
    if β != 0.0
        C = real.(exp.(-sqrt.(complex.(((Δx .* cos(β) .+ Δz .* sin(β)) ./ θx).^2 .+ 
                                       ((-Δx .* cos(β) .+ Δz .* sin(β)) ./ θz).^2))))
    else
        C = real.(exp.(-sqrt.(complex.((Δx ./ θx).^2 .+ (Δz ./ θz).^2))))
    end
    
    C[diagind(C)] .= 1.0
    cϕ = cholesky(C).L * randn(Float64, nmp, 2)
    p = 0.5
    R = [1.0 0.0; p sqrt(1.0 - p^2)]
    cϕ = R * cϕ'
    c = μc .+ σc .* cϕ[1, :]
    ϕ = μϕ .+ σϕ .* cϕ[2, :]
    
    return (c, ϕ)
end

# Alias for consistency
const GRF = RFS
