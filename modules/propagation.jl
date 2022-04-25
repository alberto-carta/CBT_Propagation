using Expokit
using SparseArrays
using ProgressMeter

include("./tripole.jl")
include("./types.jl")

function propagate_power(pstart_e, Mark_M, iterations, params::input_parameters, prop_params::propagation_parameters)
    """
    Propagate the dynamics by repeated  matrix multiplication
    """
    
    #extract quantities from parameters
    Mers    = params.Mers
    n_vib   = params.n_vib
    D       = params.D_orth
    dz = params.dz
    diff_lim = prop_params.diff_lim
    Wo = prop_params.Wo


    
    block_size = D * n_vib #gives the size of the monomer part of the hamiltonian
    mat_size = block_size * Mers
    Merschain = 1:Mers 
    diff_lim_eff = diff_lim#*0.69

    parr_e = zeros(mat_size, iterations)
    parr_n = zeros(Mers, iterations)
    coarse_arr = zeros(Mers, iterations)

    time_vec = LinRange(1 : iterations)*Dt*6.582*(10^(-7) ) #convert to nanoseconds in the timescale
    MSD_vec = zeros(iterations)
    
    Dens_M = abs.(eig_V).^2 #computes the exciton probability density of the chain for each eigenstate

    p_e = pstart_e

    Mark_M = Mark_M^Wo

    for i in 1:iterations
        # println("currently in iteration $i")
    
        p_e = Mark_M * p_e
    
        p_n = Dens_M * p_e
    
        psingle_n = reshape(p_n, (block_size, Mers))'#reshape vector to couple states on the same site
        psingle_n = sum( psingle_n, dims=2) #add together the occupancies in the first and second vibrational states
    
        #convolute with a gausssian to account for diffraction limit
        coarse_n = convolute_difflim(psingle_n, Merschain, diff_lim_eff) 
        
        #even if you shine the laser on in the center of the fiber, the initial exciton population is determined
        #by the projection of the eigenstates on the position basis, this does not guarantee that your startpos
        #coincides with the position of the actual first excitation, calculate the mean position of the population
        #after the first measurement you can perform (which is what is done experimentally) to correct for that.
        # if i == 1
        #     first_meanpos = sum(coarse_n .*(Merschain))
        #     println("Printing first mean position")
        #     println( (first_meanpos-Mers/2)*dz*10^(-4))
        # end
        
        # MSD_eff =sum(coarse_n.*((Merschain .- first_meanpos).^2) .* (dz*10^(-4) )  )   #calculate the msd 
        coarse_n = coarse_n/maximum(coarse_n) #set maximum to 1 to compare with experiments
        
        # MSD_vec[i] = MSD_eff
        parr_e[:, i] = p_e 
        parr_n[:, i] = psingle_n
        coarse_arr[:, i] = coarse_n
    end

  return time_vec, coarse_arr 
end

function propagate_expokit(pstart_e, Renorm_M, iterations,Dt, params::input_parameters, prop_params::propagation_parameters)
    #extract quantities from parameters
    Mers    = params.Mers
    n_vib   = params.n_vib
    D       = params.D_orth
    dz = params.dz
    diff_lim = prop_params.diff_lim
    Wo = prop_params.Wo


    
    block_size = D * n_vib #gives the size of the monomer part of the hamiltonian
    mat_size = block_size * Mers
    Merschain = 1:Mers 
    diff_lim_eff = diff_lim*1#*0.69

    parr_e = zeros(mat_size, iterations)
    parr_n = zeros(Mers, iterations)
    coarse_arr = zeros(Mers, iterations)

    time_vec = LinRange(1 : iterations)*Dt*6.582*(10^(-7) ) #convert to nanoseconds in the timescale
    MSD_vec = zeros(iterations)
    
    Dens_M = abs.(eig_V).^2 #computes the exciton probability density of the chain for each eigenstate

    p_e = copy(pstart_e)

    Rsparse = Renorm_M
    println("Exponentiating R matrix with pade approximation")
    Mark_M = padm(Dt*Renorm_M, p=6)

    # Rsparse = sparse(Renorm_M)
    # droptol!(Rsparse, 10e-15)
    # display(Rsparse) 

    @showprogress for i in 1:iterations
        
        # expmv!(Dt, Rsparse, p_e)
        p_e = Mark_M * p_e
    
        p_n = Dens_M * p_e
        
        psingle_n = reshape(p_n, (block_size, Mers))'#reshape vector to couple states on the same site
        psingle_n = sum( psingle_n, dims=2) #add together the occupancies in the first and second vibrational states
    
        #convolute with a gausssian to account for diffraction limit
        coarse_n = convolute_difflim(psingle_n, Merschain, diff_lim_eff) 
        
        # compute initial timestep mean positon
        if i == 1
            global first_meanpos = sum(coarse_n .*(Merschain))
        end
        
        
        # compute MSD
        MSD_eff = sum( coarse_n.*((Merschain .- first_meanpos).^2) .* (dz*10^(-4) )^2 )   #calculate the msd 
        
        # normalize 
        # coarse_n = coarse_n/maximum(coarse_n) #set maximum to 1 to compare with experiments
        
        
        
        # save the results to an array
        coarse_arr[:, i]    = coarse_n
        parr_e[:, i]        = p_e
        parr_n[:, i]        = psingle_n
        MSD_vec[i]       = MSD_eff
  
    end
    
    return time_vec, coarse_arr, parr_n, parr_e, MSD_vec


end

function convolute_difflim(y, x, sigma)
    @assert length(y)==length(x) "Convolute error: array and chain of different dimensions"
    arrlen = length(y)
    convoluted_arr = zeros(Mers)
        
    for m in Merschain    
            convoluted_arr = convoluted_arr + gaussian.(x, m, sigma )*y[m] 
    end

    convoluted_arr = convoluted_arr/sum(convoluted_arr)

    return convoluted_arr
end