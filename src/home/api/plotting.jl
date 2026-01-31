export plot_ρ, plot_coh, plot_Δu, plot_P, plot_Δϵp

"""
    plot_ρ(x, ρ, clim, path_plot) -> Nothing

Plot density field.

# Arguments
- `x`: Material point positions
- `ρ`: Density values
- `clim`: Color limits (min, max)
- `path_plot`: Output path for saving plots
"""
@views function plot_ρ(x, ρ, clim, path_plot="./out/")
    gr(size=(2*250, 2*125), legend=true, markersize=1.0)
    scatter(x[:, 1], x[:, 2], zcolor=ρ,
        markershape=:circle,
        label="",
        show=true,
        aspect_ratio=1,
        c=:viridis,
        clims=(clim[1], clim[2]),
        ylim=(0.0, 10.0),
    )     
end

"""
    plot_coh(xp, coh, phi, ϕ0) -> Nothing

Plot cohesion and friction angle fields.

# Arguments
- `xp`: Material point positions
- `coh`: Cohesion values
- `phi`: Friction angle values
- `ϕ0`: Reference friction angle
"""
@views function plot_coh(xp, coh, phi, ϕ0)
    path_plot = setup_paths()
    
    gr(size=(2.0*250, 2*125), legend=true, markersize=2.25, markerstrokecolor=:auto)
    scatter(xp[:, 1], xp[:, 2], zcolor=coh./1e3,
        markershape=:circle,
        label="",
        show=true,
        aspect_ratio=1,
        c=:vik,
        clims=(10.0, 30.0),
        markerstrokecolor=:auto,
        markerstrokewidth=0,
        ylim=(-10, 20),
    )
    savefig(path_plot * "coh0.png")
    
    gr(size=(2.0*250, 2*125), legend=true, markersize=2.25, markerstrokecolor=:auto)
    scatter(xp[:, 1], xp[:, 2], zcolor=phi,
        markershape=:circle,
        label="",
        show=true,
        aspect_ratio=1,
        c=:vik,
        clims=(ϕ0 - ϕ0/5, ϕ0 + ϕ0/5),
        markerstrokecolor=:auto,
        markerstrokewidth=0,
        ylim=(-10, 20),
    )
    savefig(path_plot * "phi0.png")
end

"""
    plot_Δu(xp, up) -> Nothing

Plot displacement magnitude field.

# Arguments
- `xp`: Material point positions
- `up`: Displacement vectors
"""
@views function plot_Δu(xp, up)
    Δu = sqrt.(up[:, 1].^2 + up[:, 2].^2)
    gr(size=(2*250, 2*125), legend=true, markersize=2.5)
    scatter(xp[:, 1], xp[:, 2], zcolor=Δu,
        markershape=:circle,
        label="",
        show=true,
        aspect_ratio=1,
        c=:viridis,
        ylim=(-10.0, 20.0),
    )     
end

"""
    plot_P(xp, σ) -> Nothing

Plot pressure field.

# Arguments
- `xp`: Material point positions
- `σ`: Stress tensor
"""
@views function plot_P(xp, σ)
    p = -(σ[1, :] + σ[2, :] + σ[3, :]) / 3 / 1e3
    gr(size=(2*250, 2*125), legend=true, markersize=2.5)
    scatter(xp[:, 1], xp[:, 2], zcolor=p,
        markershape=:circle,
        label="",
        show=true,
        aspect_ratio=1,
        c=:viridis,
        ylim=(-10.0, 20.0),
    )     
end

"""
    plot_Δϵp(xp, epII) -> Nothing

Plot equivalent plastic strain field.

# Arguments
- `xp`: Material point positions
- `epII`: Equivalent plastic strain
"""
@views function plot_Δϵp(xp, epII)
    gr(size=(2*250, 2*125), legend=true, markersize=2.5)
    scatter(xp[:, 1], xp[:, 2], zcolor=epII,
        markershape=:circle,
        label="",
        show=true,
        aspect_ratio=1,
        c=:viridis,
        clims=(0.0, 2.0),
        ylim=(-10.0, 20.0),
    )     
end
