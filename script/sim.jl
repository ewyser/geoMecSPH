@views function sim()
    # arithmetic precision (double=Float64 or single=Float32)
    typeD = Float64  
    # ---------------------------------------------------------------------------
    nel = 200
    # ---------------------------------------------------------------------------
    # non-dimensional constant 
    # ---------------------------------------------------------------------------
    ν       = 0.3                                                               # poisson's ratio                                                        
    ni      = 3                                                                 # number of material point along 1d
    nstr    = 4                                                                 # number of stresses
    # ---------------------------------------------------------------------------
    # independant physical constant
    # ---------------------------------------------------------------------------
    g       = 9.81                                                              # gravitationnal acceleration [m/s^2]
    E       = 1.0e6                                                             # young's modulus             [Pa]
    K,G,Del = D(E,ν)                                                            # elastic matrix    
    ρ0      = 2700.0                                                            # density                     [kg/m^3]
    yd      = sqrt((K+4.0/3.0*G)/ρ0)                                            # elastic wave speed          [m/s]
    c0      = 20.0e3                                                            # cohesion                    [Pa]
    ϕ0      = 20.0*pi/180                                                       # friction angle              [Rad]
    ψ0      = 0.0                                                               # dilantancy angle
    H       = -60.0e3                                                           # softening modulus           [Pa]
    cr      =   4.0e3                                                           # residual cohesion           [Pa]
    ϕr      = 7.5*pi/180                                                        # residual friction angle     [Rad]
    t       = 15.0                                                              # simulation time             [s]
    te      = 10.0                                                              # elastic loading             [s]
    tg      = te/1.5                                                            # gravity increase 
    # ---------------------------------------------------------------------------
    # mesh & mp setup
    # ---------------------------------------------------------------------------
    lx      = 64.1584                                                           # domain length along the x-direction
    lz      = 12.80                                                             # domain length along the z-direction
    meD,bc  = meshSetup(nel,lx,lz,typeD)                                        # mesh geometry setup
    mpD     = pointSetup(meD,ni,lz,c0,cr,ϕ0,ϕr,ρ0,nstr,typeD)                   # material point geometry setup
    Hp      = H*meD.h[1]                                                        # softening modulus
    # ---------------------------------------------------------------------------
    # display parameters & runtime
    # ---------------------------------------------------------------------------                                                            
    Δt      = 0.5*meD.h[1]/yd                                                   # unconditionally stable timestep
    nit     = ceil(t/Δt)                                                        # maximum number of interation
    nf      = max(2,ceil(round(1/Δt)/25))                                       # number of frame interval
    # runtime parameters
    it      = 1                                                                 # initialize iteration
    tw      = 0.0                                                               # initialize time while statement
    
    #char    = save2txt(meD,mpD,bc)
    #p       = [g;ρ0;ψ0;ν;E;Kc;Gc;cr;Hp;t;te;tg]
    #writedlm("/Users/manuwyser/Dropbox/PhD_Thesis/git_local/work_mpm/C_code_2D/scripts/setting_Exp2b/phys.txt" ,p)

    rho   = ρ0*ones(typeD,mpD.nmp,2) 

    tw = 0.0
    t  = 0.25
    it = 0
    itps = 0.0
    nout = 5
    wct  = 0.0
    flag = 0

    println("o---------------------------------------------o")
    println("|             ** geoMecSPH v1.0 **            |")
    println("|      -- finite strain formulation --        |")
    println("o---------------------------------------------o")
    @info "initial geometry:" nel=Int64(meD.nel[3]) nno=meD.nno[3] nmp=mpD.nmp
    println("o---------------------------------------------o") 
    println("[=> action!")
    prog  = ProgressUnknown("working hard:", spinner=true,showspeed=true)
    Q = zeros(size(mpD.ϕ))
    while tw<t
        ## 1st SPH kernel using a Wendland type kernel
        # ------------------------------------------------------------------------------------------------------------------------------------------------------------#
        α = (7.0/(4.0*π))
        for i in 1:mpD.nmp
            Δρ = 0.0
            for j in 1:mpD.nmp
                    rx = (mpD.x[i,1]-mpD.x[j,1])
                    rz = (mpD.x[i,2]-mpD.x[j,2])
                    r  = sqrt(rx^2+rz^2)
                    q  = r/mpD.h[i]
                    if 0.0<=q<=2.0
                        f = ((1.0-0.5*q)^4)*(1.0+2.0*q)
                        ∂f= -5.0*q*(1.0-0.5*q)^3
                    elseif 2.0 < q
                        f = 0.0
                        ∂f= 0.0
                    end
                    Q[i,j]       = q
                    mpD.ϕ[i,j]   = α*f
                    mpD.∂ϕx[i,j] = α*∂f/mpD.h[i]#*(rx/r)
                    mpD.∂ϕz[i,j] = α*∂f/mpD.h[i]#*(rz/r)
                    
                    Δρ += mpD.m[j]*((mpD.v[i,1]-mpD.v[j,1])*mpD.∂ϕx[i,j]+(mpD.v[i,2]-mpD.v[j,2])*mpD.∂ϕz[i,j])
                    
            end
            rho[i]+=Δt*Δρ
        end

        mpD.x .= mpD.x.+Δt*mpD.v

        tw += Δt
        it += 1
        if(mod(it,nout)==0)
            plot_ρ(mpD.x,rho.-ρ0,(-10.0,10.0))    
        end 
        next!(prog;showvalues = [("[nel,np]",(round(Int64,meD.nel[1]*meD.nel[2]),mpD.nmp)),("iteration(s)",it),("(✗) t/T",round(tw/t,digits=2))])
    end
    ProgressMeter.finish!(prog, spinner = '✓',showvalues = [("[nel,np]",(round(Int64,meD.nel[1]*meD.nel[2]),mpD.nmp)),("iteration(s)",it),("(✓) t/T",1.0)])
    savefig(path_plot*"plot.png")
    @info "Figs saved in" path_plot

    gr(size=(50,100),legend=true,markersize=2.5)



    default(
        fontfamily="Computer Modern",
        linewidth=2,
        framestyle=:box,
        label=nothing,
        grid=true
        )
        i = 10
        q = hcat(Q[i,:])
        w = hcat(mpD.ϕ[i,:])
        gr(size=(200,200))
        plot(q,w,
            markershape=:circle,
            label="",
            show=true,
            xlim=(0.0,3.0),
            ylim=(-0.8,0.8),
            xlabel=L"r_{kl}=\dfrac{1}{h} \| \| \mathbf{r}_{k}-\mathbf{r}_{l} \| \| _2",
            ylabel=L"\omega_k(r_{kl})"
        ) 
        savefig(path_plot*"wendland_kernel.png")
    default(
        fontfamily="Computer Modern",
        linewidth=2,
        framestyle=:box,
        label=nothing,
        grid=true
        )
        i = 10
        q = hcat(Q[i,:])
        w = hcat(mpD.∂ϕx[i,:])
        gr(size=(200,200))
        plot(q,w,
            markershape=:circle,
            label="",
            show=true,
            xlim=(0.0,3.0),
            ylim=(-0.8,0.8),
            xlabel=L"r_{kl}=\dfrac{1}{h} \| \| \mathbf{r}_{k}-\mathbf{r}_{l} \| \| _2",
            ylabel=L"\partial\omega_k(r_{kl})"
        ) 
        savefig(path_plot*"wendland_kernel_derivative.png")



    println("[=> done! exiting...")
end



# https://techytok.com/lesson-parallel-computing/