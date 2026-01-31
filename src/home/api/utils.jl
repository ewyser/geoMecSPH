export get_Δt, get_g

"""
    get_Δt(vp, h, yd) -> Float64

Compute stable timestep based on CFL condition.

# Arguments
- `vp`: Particle velocities
- `h`: Grid spacing [Δx, Δz]
- `yd`: Elastic wave speed

# Returns
- `Δt`: Stable timestep
"""
@views function get_Δt(vp, h, yd)
    Δx = h[1]
    Δz = h[2]
    vmax = [abs.(vp[:, 1]) abs.(vp[:, 2])]
    cmax = [maximum(vmax[:, 1]) maximum(vmax[:, 2])]
    cmax = [Δx/(cmax[1] + yd) Δz/(cmax[2] + yd)]
    Δt = 0.5 * maximum(cmax)
    return Δt
end

"""
    get_g(tw::Float64, tg::Float64) -> Float64

Compute gravity with gradual increase.

# Arguments
- `tw`: Current time
- `tg`: Time to reach full gravity

# Returns
- `g`: Gravitational acceleration
"""
function get_g(tw::Float64, tg::Float64)
    g = 0.0
    if tw <= tg
        g = 9.81 * tw / tg
    else
        g = 9.81
    end
    return g
end
