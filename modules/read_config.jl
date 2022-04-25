using TOML
include("./types.jl")
include("./tripole.jl")

function read_config_file(filename="./config.toml")

    println("#######################################################################  ")
    println("Setting globals")
    println("")
    
    input_dict = TOML.parsefile(filename)
    global loadfile = input_dict["general"]["loadfile"]
    global calcdir = input_dict["general"]["calcdir"]
    global tripole = input_dict["general"]["chain"]["tripole"]
    global compute_propagation = input_dict["propagation"]["compute_propagation"]
    
    println("Loading from previous file: $loadfile")
    println("The output directory is set to $calcdir")
    
    if tripole
        println("Using tripole coupling")
        
        trp_fac = input_dict["general"]["tripole"]["trp_fac"]
        θ = input_dict["general"]["tripole"]["θ"]
        dz        = input_dict["general"]["chain"]["dz"] 
        delz        = input_dict["general"]["tripole"]["delz"] 
        J = trp_fac*same_int(1*(dz-delz), 1, (mod(1*θ+60,120)-60) ) 

        calc_params =  input_parameters(
           Mers      =  input_dict["general"]["chain"]["Mers"],
           dz        =  input_dict["general"]["chain"]["dz"],
           coupling  =  input_dict["general"]["chain"]["coupling"],
           J_coupl   =  J,
           neighbors =  input_dict["general"]["chain"]["neighbors"],
           
           λ         = input_dict["general"]["chain"]["λ"],
           hw_vibr   = input_dict["general"]["chain"]["hw_vibr"],
           n_vib     = input_dict["general"]["chain"]["n_vib"],

           
           γ         = input_dict["general"]["spectrum"]["γ"],
           D_shift   = input_dict["general"]["spectrum"]["D_shift"],
           hw_00     = input_dict["general"]["spectrum"]["hw_00"]-input_dict["general"]["spectrum"]["D_shift"],
           
           delz      = input_dict["general"]["tripole"]["delz"],
           trp_fac   = input_dict["general"]["tripole"]["trp_fac"],
           D_orth    = input_dict["general"]["tripole"]["D_orth"],
           θ         = input_dict["general"]["tripole"]["θ"],
           
           σ         = input_dict["general"]["chain"]["σ"],
           L0        = input_dict["general"]["chain"]["L0"]
        )
    else
        println("Using dipole coupling")
        
        calc_params =  input_parameters(
           Mers      =  input_dict["general"]["chain"]["Mers"],
           dz        =  input_dict["general"]["chain"]["dz"],
           coupling  =  input_dict["general"]["chain"]["coupling"],
           J_coupl   =  input_dict["general"]["chain"]["J_coupl"],
           neighbors =  input_dict["general"]["chain"]["neighbors"],
           
           λ         = input_dict["general"]["chain"]["λ"],
           hw_vibr   = input_dict["general"]["chain"]["hw_vibr"],
           n_vib     = input_dict["general"]["chain"]["n_vib"],

           
           γ         = input_dict["general"]["spectrum"]["γ"],
           D_shift   = input_dict["general"]["spectrum"]["D_shift"],
           hw_00     = input_dict["general"]["spectrum"]["hw_00"]-input_dict["general"]["spectrum"]["D_shift"],
           
           delz      = 0.0, 
           trp_fac   = 0.0, 
           D_orth    = 1, 
           θ         = 0.0, 
           
           σ         = input_dict["general"]["chain"]["σ"],
           L0        = input_dict["general"]["chain"]["L0"]
        )
    
    end

    if compute_propagation
        println("Computing propagation set to true")
        startpos = calc_params.Mers ÷ 2
        prop_params =  propagation_parameters(
       
            kbT        =  input_dict["propagation"]["bath"]["kbT"],
            wc         =  input_dict["propagation"]["bath"]["wc"],
            Wo         =  input_dict["propagation"]["bath"]["Wo"],
            
            startpos   =  startpos,
            sigmastart =  input_dict["propagation"]["chain"]["sigmastart"],
            diff_lim   =  input_dict["propagation"]["chain"]["diff_lim"]
        )
    else
        println("Not computing propagation")
        
        prop_params =  propagation_parameters(
       
            kbT        =  0.0,
            wc         =  0.0,
            Wo         =  0.0,
            
            startpos   =  0,
            sigmastart =  0.0,
            diff_lim   =  0.0
        )

    end
    
    println("")
    println("#######################################################################  ")
    println("")
        

    return calc_params, prop_params, input_dict
end