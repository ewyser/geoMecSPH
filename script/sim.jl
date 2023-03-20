@views function sim()
    # arithmetic precision (double=Float64 or single=Float32)
    typeD = Float64  
    # ---------------------------------------------------------------------------
    nel = 200
    # ---------------------------------------------------------------------------
    # non-dimensional constant 
    # ---------------------------------------------------------------------------
    ν       = 0.3                                                               # poisson's ratio                                                        
    ni      = 2                                                                 # number of material point along 1d
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

    

    tw = 0.0
    it = 0
    itps = 0.0
    nout = 50
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
    while tw<t

        tw += Δt
        it += 1
        if(mod(it,nout)==0)
            plot_Δϵp(mpD.xp,mpD.epII)       
        end 
        next!(prog;showvalues = [("[nel,np]",(round(Int64,meD.nel[1]*meD.nel[2]),mpD.nmp)),("iteration(s)",it),("(✗) t/T",round(tw/t,digits=2))])
    end
    ProgressMeter.finish!(prog, spinner = '✓',showvalues = [("[nel,np]",(round(Int64,meD.nel[1]*meD.nel[2]),mpD.nmp)),("iteration(s)",it),("(✓) t/T",1.0)])
    savefig(path_plot*"plot.png")
    @info "Figs saved in" path_plot
    println("[=> done! exiting...")
end



# https://techytok.com/lesson-parallel-computing/