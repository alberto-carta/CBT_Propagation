##
# ---------------------------------------------------------- 
using PyCall
# using Plots
using LaTeXStrings

using Distributions
using LinearAlgebra
using BlockArrays
using DelimitedFiles 
using BenchmarkTools
using LsqFit
using Expokit
using SparseArrays

using ProgressMeter


include("./modules/read_config.jl")
include("./modules/tripole.jl")
include("./modules/plotroutines.jl")
include("./modules/builder.jl")
include("./modules/propagation.jl")
include("./modules/types.jl")
include("./modules/spectral_calculator.jl")

print("Loaded everything \n")



##

# λ = 0.75 # Square root of the Huang-Rhys factor
# hw_vibr = 0.19375*1 #vibrational energy

# γ = 0.026 #0.0969*1.2 # HWHM of the lorentzian shape function OR standard deviation of a gaussian dynamic disorder
# dz = 3.4 # distance between monomers in the chain in Angstrom
# delz = 0.4 # taken from Saikin, accounts for an effective overlap of the wavefunctions
# trp_fac= 4.1 # 4.7 tripole interaction multiplicative factor (influences q^2, where q is the tripole charge)
# neighbors = 1 # number of neighbors in the chain
# θ = 28 #rotation angle in degrees
# Mers = 3000 #number of monomers in the chain

# Merschain = 1:Mers #iterable object useful for the chain indexing 

# D=1 #number of orthogonal degenerate states
# n_vib= 2 #how many basis function in the target displaced oscillator you want to approximate, 5 give convergence to actual spectrum

# coupling = true
# #randomization parameters
# # σ = hw_vibr*0.6 #static noise  0.08075
# σ = 0.075
# σ_conv = hw_vibr*0.58 #convolution for quick estimation of disorder effect, equivalent to no correlation and usable only for small sigma
# D_shift = 0.36 # gas-to-crystal shift
# hw_00 = 2.681*1.005-D_shift
# L0 =  20 #70  correlation length


# kbT = 0.026 #thermal energy in eV
# wc = 100 #cutoff freq for ohmic bath
# Wo = 1 #coupling to bath constant

# startpos = Mers ÷ 2  #starting position

# sigmastart = 340/2.355/0.34 #broadening of the gaussian at the beginning
# diff_lim = 206/2.355/0.34 #PSF taken from optical FWHM = sigma*2.355 = 206nm, transformed in the position basis

# sigmastart = 200
# diff_lim = 120

function main_loop()
    blas_n_threads = BLAS.get_num_threads()
    println("Multithreading on $blas_n_threads threads")
    
    global Dt = 100000 # timestep measured in hbar/eV, about 0.6 femtoseconds 

    global calc_params, prop_params, full_dict = read_config_file() # this also sets the calcdir folder
    global Mers =  calc_params.Mers
    global Merschain =  1:Mers
    global xchain = (Merschain.-Mers/2).*(calc_params.dz*10^(-4)) #calculate the position along the axis in μm

    #name of the computation
    mkpath(calcdir)

    
    open("$calcdir/config.toml", "w") do io
        TOML.print(io, full_dict)
    end

    #choose whether to save or load an existing propagation file 
    savefile = true
    savepath = string("propagations/En_pos_L0 = ", calc_params.L0, "3rd prop_trapping" )
    #match(r^"(L0 )+\d+", savepath)
    # loadfile = false
    #loadpath = string("propagations/En_pos_L0 = 70/En_pos.txt")
    loadpath = string("./calculations/Trial_05/En_pos.dat")
    ;
    ##
    # ---------------------------------------------------------- 

    # μ_E = ones(calc_params.Mers)*hw_00
    if loadfile == false
        En_pos = build_En_chain(calc_params) 
        println("The correlation length is " , calc_params.L0)


    else
        println("Overriding specifications, loading file")
        En_pos = readdlm(loadpath)#[end:-1:1]
    end

    if savefile == true
        mkpath(calcdir)
        filename = string(calcdir, "/En_pos.dat"  )
        writedlm(filename,  En_pos)
    end

    MerschainPlot =  plot_Merschain(xchain, En_pos, calc_params)

    ## Hamiltonian 

    # build the coupling part of the Hamiltonian
    Coupl_M, Coupl_Arr = build_neigh_block(calc_params)
    # build the Hamiltonian
    Ham_M = Hermitian(build_hamiltonian(En_pos, Coupl_Arr, calc_params))
    println("Printing to line the Hamiltonian")
    display(Ham_M)

    ## Diagonalization
    println("Starting diagonalization")
    eig_decomp = @time eigen(Ham_M)
    global eig_En = eig_decomp.values
    global eig_V = eig_decomp.vectors
    println("Finished diagonalization")

    ## Absorption Spectrum
    println("Starting computing compute_absorption")
    x_intens, spectrum = @time compute_absorption(eig_V, eig_En, calc_params)
    println("Finished computing compute_absorption")

    AbsPlot =  plot_absorption(x_intens, spectrum, calc_params)

    ## Propagation

    if compute_propagation

        println("Building R matrix")
        Renorm_M= @time build_R_matrix(eig_En, eig_V, calc_params, prop_params)
        # Renorm_sparse = sparse(Renorm_M)
        # droptol!(Renorm_sparse, 10e-6)
        display(Renorm_M)

        pstart_e, pstart_n = build_start_vectors(Mers ÷ 2, eig_V, calc_params, prop_params)
        ####################################################
        # println("Exponentiating R matrix")
        # Mark_M  =  @time exp(Wo*Dt*Renorm_M)
        # time_vec, parr_n = @time propagate_power(pstart_e, Mark_M, 100, calc_params, prop_params) 
        ####################################################
        println("Exponentiating with expokit")
        time_vec, coarse_n, parr_n,  parr_e, MSD_vec = @time propagate_expokit(pstart_e, Renorm_M, 100, Dt, calc_params, prop_params) 



        plot_propagation_rspace(xchain, time_vec, parr_n, "unconvoluted_heatmap")
        plot_propagation_rspace(xchain, time_vec, coarse_n, "convoluted_heatmap")
        plot_propagation_rspace(eig_En, time_vec, parr_e, "energy_heatmap", "Blues")
        plot_slices_rspace(xchain, parr_n, time_vec)
        plot_slices_espace(eig_En, parr_e, time_vec, prop_params)

        println("Interpolating MSD")
        coeffs = fit_msd(time_vec, MSD_vec, cut=40)
        plot_MSD(time_vec, MSD_vec, coeffs)
    end

end

