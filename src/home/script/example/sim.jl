export sim

"""
    sim() -> Nothing

Main SPH simulation function with Wendland kernel.

Implements a complete SPH simulation with:
- Wendland kernel function
- Density evolution
- Time integration
- Real-time visualization
"""
@views function sim()
    # Arithmetic precision (double=Float64 or single=Float32)
    typeD = Float64  
    
    # Parameters
    nel = 200
    ν = 0.3          # Poisson's ratio                                                        
    ni = 3           # Number of material points along 1d
    nstr = 4         # Number of stresses
    
    # Physical constants
    g = 9.81         # Gravitational acceleration [m/s^2]
    E = 1.0e6        # Young's modulus [Pa]
    K, G, Del = elasticity_matrix(E, ν)
    ρ0 = 2700.0      # Density [kg/m^3]
    yd = sqrt((K + 4.0/3.0*G) / ρ0)  # Elastic wave speed [m/s]
    c0 = 20.0e3      # Cohesion [Pa]
    ϕ0 = 20.0 * pi / 180    # Friction angle [Rad]
    ψ0 = 0.0         # Dilatancy angle
    H = -60.0e3      # Softening modulus [Pa]
    cr = 4.0e3       # Residual cohesion [Pa]
    ϕr = 7.5 * pi / 180     # Residual friction angle [Rad]
    t = 15.0         # Simulation time [s]
    te = 10.0        # Elastic loading [s]
    tg = te / 1.5    # Gravity increase 
    
    # Mesh & material point setup
    lx = 64.1584     # Domain length along x-direction
    lz = 12.80       # Domain length along z-direction
    meD, bc = setup_mesh(nel, lx, lz, typeD)
    mpD = setup_mpts(meD, ni, lz, c0, cr, ϕ0, ϕr, ρ0, nstr, typeD)
    Hp = H * meD.h[1]  # Softening modulus
    
    # Display parameters & runtime
    Δt = 0.5 * meD.h[1] / yd  # Unconditionally stable timestep
    nit = ceil(t / Δt)         # Maximum number of iterations
    nf = max(2, ceil(round(1/Δt) / 25))  # Number of frame intervals
    
    # Runtime parameters
    it = 1
    tw = 0.0
    rho = ρ0 * ones(typeD, mpD.nmp, 2)
    
    tw = 0.0
    t = 0.25
    it = 0
    nout = 5
    
    # Setup output path
    path_plot = setup_paths()
    
    println("o---------------------------------------------o")
    println("|             ** geoMecSPH v1.0 **            |")
    println("|      -- finite strain formulation --        |")
    println("o---------------------------------------------o")
    @info "initial geometry:" nel=Int64(meD.nel[3]) nno=meD.nno[3] nmp=mpD.nmp
    println("o---------------------------------------------o") 
    println("[=> action!")
    
    prog = ProgressUnknown("working hard:", spinner=true, showspeed=true)
    Q = zeros(size(mpD.ϕ))
    
    # Main time loop
    while tw < t
        # 1st SPH kernel using a Wendland type kernel
        α = (7.0 / (4.0 * π))
        for i in 1:mpD.nmp
            Δρ = 0.0
            for j in 1:mpD.nmp
                rx = (mpD.x[i, 1] - mpD.x[j, 1])
                rz = (mpD.x[i, 2] - mpD.x[j, 2])
                r = sqrt(rx^2 + rz^2)
                q = r / mpD.h[i]
                
                if 0.0 <= q <= 2.0
                    f = ((1.0 - 0.5*q)^4) * (1.0 + 2.0*q)
                    ∂f = -5.0 * q * (1.0 - 0.5*q)^3
                elseif 2.0 < q
                    f = 0.0
                    ∂f = 0.0
                end
                
                Q[i, j] = q
                mpD.ϕ[i, j] = α * f
                mpD.∂ϕx[i, j] = α * ∂f / mpD.h[i]
                mpD.∂ϕz[i, j] = α * ∂f / mpD.h[i]
                
                Δρ += mpD.m[j] * ((mpD.v[i, 1] - mpD.v[j, 1]) * mpD.∂ϕx[i, j] + 
                                  (mpD.v[i, 2] - mpD.v[j, 2]) * mpD.∂ϕz[i, j])
            end
            rho[i] += Δt * Δρ
        end
        
        mpD.x .= mpD.x .+ Δt * mpD.v
        
        tw += Δt
        it += 1
        
        if mod(it, nout) == 0
            plot_ρ(mpD.x, rho .- ρ0, (-10.0, 10.0), path_plot)    
        end 
        
        next!(prog; showvalues=[("[nel,np]", (round(Int64, meD.nel[1]*meD.nel[2]), mpD.nmp)),
                                 ("iteration(s)", it),
                                 ("(✗) t/T", round(tw/t, digits=2))])
    end
    
    ProgressMeter.finish!(prog, spinner='✓', 
                         showvalues=[("[nel,np]", (round(Int64, meD.nel[1]*meD.nel[2]), mpD.nmp)),
                                    ("iteration(s)", it),
                                    ("(✓) t/T", 1.0)])
    
    savefig(path_plot * "plot.png")
    @info "Figs saved in" path_plot
    
    # Plot kernel functions
    gr(size=(50, 100), legend=true, markersize=2.5)
    default(fontfamily="Computer Modern", linewidth=2, framestyle=:box, 
            label=nothing, grid=true)
    
    i = 10
    q = hcat(Q[i, :])
    w = hcat(mpD.ϕ[i, :])
    gr(size=(200, 200))
    plot(q, w, markershape=:circle, label="", show=true,
         xlim=(0.0, 3.0), ylim=(-0.8, 0.8),
         xlabel=L"r_{kl}=\dfrac{1}{h} \| \| \mathbf{r}_{k}-\mathbf{r}_{l} \| \| _2",
         ylabel=L"\omega_k(r_{kl})") 
    savefig(path_plot * "wendland_kernel.png")
    
    # Plot kernel derivative
    default(fontfamily="Computer Modern", linewidth=2, framestyle=:box,
            label=nothing, grid=true)
    w = hcat(mpD.∂ϕx[i, :])
    gr(size=(200, 200))
    plot(q, w, markershape=:circle, label="", show=true,
         xlim=(0.0, 3.0), ylim=(-0.8, 0.8),
         xlabel=L"r_{kl}=\dfrac{1}{h} \| \| \mathbf{r}_{k}-\mathbf{r}_{l} \| \| _2",
         ylabel=L"\partial\omega_k(r_{kl})") 
    savefig(path_plot * "wendland_kernel_derivative.png")
    
    println("[=> done! exiting...")
    return nothing
end
