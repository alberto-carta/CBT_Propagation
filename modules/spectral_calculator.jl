
include("./tripole.jl")
include("./types.jl")

function compute_absorption(eig_V, eig_En, params, x_size=1000)
    
    Mers = params.Mers
    D = params.D_orth
    n_vib = params.n_vib
    λ =  params.λ
    γ =  params.γ

    mat_size = D * n_vib * Mers 
    
    FC_factor = FC_overlap_0K.(Array(0:(n_vib-1)), λ) #calculate the vector of Frank_Condon factors for each n

    Peak_intens = (eig_V' * repeat(FC_factor, Mers*D)).^2  #calculate the intensity of each peak according to the Frank-Condon principle


    x_intens = LinRange(0,5, x_size)
    envelope = zeros(x_size)

    for i in 1:mat_size
        envelope = envelope + gaussian.(x_intens, eig_En[i], γ )*Peak_intens[i]
    end

    envelope = envelope/maximum(envelope)

    return x_intens, envelope

end

function compute_IPR(eig_V, params)
    
    mat_size = params.D_orth * params.n_vib * params.Mers 
    
    IPR = zeros(mat_size) #calculate the inverse participation ratio
    for i in 1:mat_size
        IPR[i] = 1/(sum( eig_V[:,i].^(4) ))
    end

    L_D = mean(IPR)# evaluate the mean delocalization length
    return L_D, IPR
end