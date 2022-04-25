Base.@kwdef struct input_parameters
    ## Chain parameters
    Mers      :: Int64   # number of monomers in the chain
    dz        :: Float64 # distance between monomers in the chain in Angstrom
    coupling  :: Bool    # if coupling is present
    J_coupl   :: Float64 # nearestneighbour coupling
    neighbors :: Int64   # number of nearest neighbors that are affected by the interaction
    
    ## Hamiltonian parameters
    λ         :: Float64 # Square root of the Huang-Rhys factor
    hw_vibr   :: Float64 # vibrational energy
    n_vib     :: Int64   # how many basis function in the target displaced oscillator you want to approximate, 5 give convergence to actual spectrum
    
    ##Spectral parameters
    γ         :: Float64 #0.0969*1.2 # HWHM of the lorentzian shape function OR standard deviation of a gaussian dynamic disorder
    hw_00     :: Float64 # excited-ground state shift
    D_shift   :: Float64 # gas-to-crystal shift
    
    ##Tripole model specific parameters
    delz      :: Float64 # taken from Saikin, accounts for an effective overlap of the wavefunctions
    trp_fac   :: Float64 # 4.7 tripole interaction multiplicative factor (influences q^2, where q is the tripole charge)
    D_orth    :: Int64   # number of orthogonal degenerate states
    θ         :: Float64 #rotation angle in degrees
    
    ##Randomization parameters
    σ         :: Float64 #static noise  08075
    # σ_conv    :: Float64 #convolution for quick estimation of disorder effect, equivalent to no correlation and usable only for small sigma
    L0        :: Float64 #70  correlation length in monomers
end

Base.@kwdef struct propagation_parameters
    ## bath parameters
    kbT        :: Float64   #thermal energy in eV
    wc         :: Float64   #cutoff freq for ohmic bath
    Wo         :: Float64   #coupling to bath constant 
    ## positioning
    startpos   :: Int64     #starting position
    sigmastart :: Float64   #broadening of the gaussian at the beginning
    diff_lim   :: Float64   #PSF taken from optical FWHM = sigma*2.355 = 206nm, transformed in the position basis

end