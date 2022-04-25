using PyCall
using PyPlot
# using PlotlyJS
using LaTeXStrings

include("./tripole.jl")

# suppress matplotlib output 
ioff() 
# color palette from matplotlib
colorcycle = ["#1f77b4", "#ff7f0e","#2ca02c", "#d62728", "#9467bd","#8c564b", "#e377c2","#7f7f7f","#bcbd22", "#17becf"]

function plot_Merschain(Merschain, En_pos, calc_params, savefile=true)

    """
    Plots the Monomer chain
    """
    println("Plotting the monomers chain local energy")
    L0 = calc_params.L0
    
    plt.plot(Merschain,
                En_pos,
                label = L"L_0 ="*"$L0",
                # legendfontsize = 10,
                )
    plt.xlabel("Monomer number"),
    plt.ylabel("Local energy [eV]"),
    plt.rc("axes", labelsize=14) 
    plt.savefig("$calcdir/chain_energy.png", dpi=300)    
    plt.close()
    
end

function plot_absorption(x_intens, envelope, calc_params,  savefile=true)

    """
    Plots the absorption spectrum
    """
    println("Plotting the absorption spectrum")

    if calc_params.coupling == false
        exp_data = readdlm("./spectral_data/AbsTHF.dat", skipstart = 1)
        exp_label = "Experimental monomer data in THF"
    else
        exp_data = readdlm("./spectral_data/AbsDod.dat", skipstart = 1)
        exp_label = "Experimental polymer data in "
    end

    intensity = exp_data[:,2]
    omega_cm1 = exp_data[:,1]
    energy_exp= (c.*ħ./(omega_cm1.*10^(-9)))./eV 
    intensity = intensity./energy_exp

    intensity =  intensity ./ maximum(intensity)

#load the experimental data

    plt.plot(energy_exp, intensity, label = exp_label)
    plt.plot(x_intens,
                envelope,
                label = "Computed spectrum",
                linewidth=2,
                )
    plt.xlabel("Energy [eV]"),
    plt.ylabel("Intensity [a.u.]"),
    plt.rc("axes", labelsize=14) 
    plt.xlim(left = 2.0, right=3.2)
    plt.savefig("$calcdir/abs.png", dpi=300)    
    plt.close()
    
end

function plot_propagation_rspace(chain, time, parr_n, savename, cmap="inferno", savefile=true)
    """
    Plots the exciton propagation in real space
    """
    println("Plotting the propagation: heatmap")

    plt.pcolor(chain, time,   parr_n' , cmap = cmap, shading="auto")
    plt.ylim(6,0.3)
    plt.xlabel("Position along x [μm]", size = 16)
    plt.ylabel("Time [ns]", size = 16)
    plt.colorbar().set_label("Normalized intensity", size = 16)
    plt.savefig("$calcdir/$savename.png", dpi=300)
    plt.close()

end

function plot_slices_rspace(chain, parr_n, time_vec,  savefile=true)
    """
    Plots the exciton propagation in real space, slincing at timestep
    """
    
    println("Plotting real space propagation: slices")
    
    for t in [1, 10, 30, 100]
        time = round(time_vec[t], digits=3)
        p_n = parr_n[:,t]
        plt.plot(chain,
                    p_n,
                    label = "Time = $time ns",
                    linewidth=2,
                    )
    end
    plt.xlabel("Position along x [μm]")
    plt.ylabel("Normalized intensity")
    plt.legend()
    plt.rc("axes", labelsize=12) 
    plt.savefig("$calcdir/slices_r.png", dpi=300)    
    plt.close()

end

function plot_slices_espace(energies, parr_e, time_vec, prop_params::propagation_parameters,   savefile=true)
    """
    Plots the exciton propagation in energy space
    """
    kbT = prop_params.kbT
    
    println("Plotting energy space propagation: slices")
    
    for t in [1, 10, 30, 100]
        time = round(time_vec[t], digits=3)
        p_e = parr_e[:,t]
        plt.plot(energies,
                    p_e,
                    label = "Time = $time ns",
                    linewidth=2,
                    )
    end
    # plot boltzmann distro 
    Z_partition = sum(Boltz_weight.(energies[1:end], kbT))
    plot(eig_En[1:end], Boltz_weight.(energies[1:end], kbT)/Z_partition, label="Boltzmann Distribution", linewidth = 2, color = colorcycle[10])




    plt.xlabel("Energy [eV]")
    plt.ylabel("Normalized intensity")
    plt.legend()
    plt.rc("axes", labelsize=12) 
    plt.savefig("$calcdir/slices_e.png", dpi=300)    
    plt.close()

end

function plot_MSD(time_vec, MSD_vec, coeffs_fit, savefile=true)
    """
    Plots the MSD evolution
    
    Inputs:
    * time_vec (array of floats): timesteps in nanoseconds 
    * MSD_vec (array of floats): Value of the MSD at the timesteps 
    * coeffs_fit (array of floats) : array containing [α , A]

    """
    norm_time = time_vec[1:end].-time_vec[1]
    norm_MSD = MSD_vec[1:end].-MSD_vec[1]
    
    println("Plotting MSD")

    coeffs = coeffs_fit

    plt.plot(norm_time,  norm_MSD, label = "Simulation data")
    plt.plot(norm_time[1:end],  subdif_msd.(norm_time[1:end], coeffs[1], coeffs[2]), label = string("Subdiffusion under 1 ns: α = ", round(coeffs[1]; digits = 3) ) )
    println("The subdiffusion coefficient is: ", coeffs[1])
    plt.xlabel("Time [ns]", size = 14)
    plt.ylabel("MSD [μm^2]", size = 14)
    plt.legend()
    plt.savefig("$calcdir/MSD.png", dpi=300)    
    plt.close()
    
end


