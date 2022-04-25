
using BlockArrays

include("./tripole.jl")
include("./types.jl")

##
# ---------------------------------------------------------- 
"""
    build_En_chain(Mers, σ, hw_00, L0)


Build the random energy fluctuation on the chain
inputs:

    * Mers (integer): number of monomers in the chain 
    * σ (real): standard deviation of the energy distribution [eV]
    * hw_00 (real): excitation energy (mean of the distribution) [eV]
    * LO (integer)= correlation length [in Monomers] 

returns:
    * En_pos(real,array): vector of the local energies
"""
function build_En_chain(params::input_parameters)
    Mers = params.Mers
    σ = params.σ
    hw_00 = params.hw_00
    L0 = params.L0

    Merschain = 1:Mers


    μ_E = ones(Mers)*hw_00
    if L0==0
        Distro_static = MvNormal(μ_E, σ)

    else
        Pos_M = zeros(Mers, Mers ) # initialize a Mers x Mers matrix that will contain the relative positions btw the monomers

        for i ∈ Merschain
            Pos_M[i, :] = abs.( Merschain .- i) # every monomer is distant 0 monomers with itself 
        end
        Var_M = exp.(-Pos_M/L0)*(σ^2)
        Distro_static = MvNormal(μ_E, Var_M)
    end
    En_pos = rand(Distro_static)
    
    return En_pos
end

function build_vib_vect(params::input_parameters)
    n_vib  = params.n_vib
    hw_vibr = params.hw_vibr
    D = params.D_orth

    vib_vect = ( 0 : ( n_vib-1) ) .* hw_vibr #vector containing the vibrations
    vib_vect =  repeat(vib_vect, D) # repeat the vector if you want to consider degenerate states 
    
    return vib_vect
end

function build_neigh_block(params::input_parameters)

    Mers       = params.Mers
    neighbors  = params.neighbors
    coupling   = params.coupling
    n_vib      = params.n_vib
    D          = params.D_orth
    λ          = params.λ
    trp_fac    = params.trp_fac          
    θ          = params.θ          
    dz         = params.dz           
    delz       = params.delz         
    J_coupl          = params.J_coupl          
    
    block_size = D * n_vib #gives the size of the monomer part of the hamiltonian
    subblock_size = n_vib # allows you to work with the vibrational states
    mat_size = block_size * Mers
    mat_div = ones(Int8, Mers)*block_size #divide the matrix in 'Mers' blocks of size 'block_size'
    block_div = ones(Int8, D)*subblock_size #divide the block into sub-blocks





    Coupl_M = BlockArray(zeros(block_size, block_size), block_div, block_div) #initialize the coupling matrix
    Coupl_Arr = zeros(neighbors, block_size, block_size)

    #create the coupling matrix for the m-th neighbor, we make a block matrix
    if coupling == true 
        FC_factor = FC_overlap_0K.(Array(0:(n_vib-1)), λ) #calculate the vector of Frank_Condon factors for each n
        FC_M = kron(FC_factor, FC_factor') # calculate the Frank-Condon Matrix
        for m in 1:neighbors
            if tripole
                # println("using tripole")
                Coupl_Block_same = trp_fac.*FC_M.*same_int(m*(dz-delz), m, (mod(m*θ+60,120)-60) ) #modulo is for symmetrization after a rotation of 120 degrees
                Coupl_Block_cross = trp_fac.*FC_M.*cross_int(m*(dz-delz), m,  (mod(m*θ+60,120)-60)) #the mod expression keeps θ btw -60 and 60
            else
                # println("using dipole")
                Coupl_Block_same =  J_coupl * (1/m)^3 .*FC_M
                Coupl_Block_cross = J_coupl * (1/m)^3 .*FC_M
            end

                for i in 1:D
                    for j in 1:D
                        if i == j
                            setblock!(Coupl_M, Coupl_Block_same, i, i )
                        else
                            setblock!(Coupl_M, Coupl_Block_cross, i, j )
                        end
                    end
                end
            Coupl_Arr[m,:,:] = Coupl_M #create a vector of coupling matrices for each neighbor
            
        end
    end
    return  Coupl_M, Coupl_Arr
end

function build_hamiltonian(En_pos, Coupl_Arr, params::input_parameters)
    #iterate now over the chain to write the complete hamiltonian, we split the cycle in two parts:
    #the first up to Mers-neighbors monomers and then for the last monomers since they do not have m neighbors to the right
    neighbors = params.neighbors
    D         = params.D_orth
    n_vib      = params.n_vib
    D          = params.D_orth
    
    Mers = length(En_pos)    
    
    block_size = D * n_vib #gives the size of the monomer part of the hamiltonian
    subblock_size = n_vib # allows you to work with the vibrational states
    mat_size = block_size * Mers
    mat_div = ones(Int8, Mers)*block_size #divide the matrix in 'Mers' blocks of size 'block_size'
    block_div = ones(Int8, D)*subblock_size #divide the block into sub-blocks
    
    Ham_M = BlockArray(zeros(mat_size, mat_size), mat_div, mat_div) #initialize the Hamiltonian
    
    vib_vect = build_vib_vect(params)
    
    
    for i in 1:(Mers-neighbors)
        stat_M = Diagonal(ones(block_size)) .* En_pos[i]  .+ Diagonal(vib_vect) #write the energy of the excited monomers with the vibration present
        setblock!(Ham_M, stat_M, i, i) # set the diagonal blocks of the hamiltonian with the static energies
        
        for m in 1:neighbors
            setblock!(Ham_M, Coupl_Arr[m,:,:], i+m, i) #set the coupling to the m-th neighbor
            setblock!(Ham_M, Coupl_Arr[m,:,:], i, i+m)
        end

    end

    #here we consider the last monomers of the chain that lack at least one of the right neighbors
    for i in (Mers-neighbors+1):Mers 
        stat_M = Diagonal(ones(block_size)) .* En_pos[i]  .+ Diagonal(vib_vect)
        setblock!(Ham_M, stat_M, i, i)
        
        for m in 1:(Mers-i)
            setblock!(Ham_M, Coupl_Arr[m,:,:], i+m, i)
            setblock!(Ham_M, Coupl_Arr[m,:,:], i, i+m)
        end
        
        
        
    end


    Ham_M = Array(Ham_M);

    return Ham_M
    
end

function build_R_matrix(eig_En, eig_V, params, prop_params)
    
    kbT = prop_params.kbT 
    wc = prop_params.wc
    Wo = prop_params.Wo
    J_coupl = params.J_coupl

    
    Dens_M = abs.(eig_V).^2 #computes the exciton probability density of the chain for each eigenstate
    Overlap_M =  Dens_M' * Dens_M #this computes the overlap of the exciton densities of the eigenstates

    #julia broadcasting takes care of the outer product
    Ediff_M =  eig_En .- eig_En' .+ 0.000001

    #wc = 0.109
    Wscatt_M = Wo .* Bose_weight.(Ediff_M, kbT) .* (OhmSD.(Ediff_M, wc))./J_coupl .* Overlap_M

    sum_scatt = sum(Wscatt_M, dims=1) #sum over all scattering rates
    sum_scatt = dropdims(sum_scatt, dims=1) #necessary for vectorization, sum preserves the dimesion of the array in Julia

    Renorm_M = - Diagonal(sum_scatt)  + Wscatt_M ;
    return Renorm_M
end

function build_start_vectors(startpos, eig_V, params::input_parameters, prop_params::propagation_parameters)
    """
    Create an excitation via a gaussian shaped illumination:
    provides real space and energy space representation of a gaussian shaped 
    illumination to start the propagation
    """
    sigmastart = prop_params.sigmastart
    D = params.D_orth
    n_vib = params.n_vib
    Mers = params.Mers
    
    
    Dens_M = abs.(eig_V).^2 #computes the exciton probability density of the chain for each eigenstate

    # illuminate a gaussian shape within the chain
    # println("$Mers, $startpos, $sigmastart")
    pstart_n = gaussian.(1:Mers, startpos, sigmastart)
    pstart_n = repeat(pstart_n, inner = (D*n_vib) )
    pstart_n = pstart_n/sum(pstart_n) # renormalize

    #morph now to energy space
    pstart_e = Dens_M' * pstart_n

    #now delete the eigenstates that cannot be excited by a laser of energy 2.7589 ± 0.006 eV
    # unexcitable_En = findall(x-> (x <= (2.7589-0.006) || x >= (2.7589+0.006)), eig_En)
    # pstart_e[unexcitable_En] .= 0
    # pstart_e = pstart_e/sum(pstart_e);

    return pstart_e, pstart_n

    
end