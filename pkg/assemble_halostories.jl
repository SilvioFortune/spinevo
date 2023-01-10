
@doc """
DESCRIPTION:\n
INPUT:\n 
OUTPUT:\n
""" -> 
function assemble_halostories(; mass_ST_thr=2e10, outfile=replace("assembly_Mstar_$(mass_ST_thr)_$(today()).jld", "-" => ""), central_only=false,
    outdir="./OUT_assemble_halostories", indir=current_dir_stories,
    simbox=current_dir_simbox, verbose=true
    )
    println("Initiating assemble_halostories with:\n 
                mass_ST_thr = $(mass_ST_thr)\n 
                outfile = $(outfile)\n 
                central_only = $(central_only)\n 
                outdir = $(outdir)\n 
                indir = $(indir)\n 
                simbox = $(simbox)\n 
                verbose = $(verbose)
                ")

    # Setup
    assembly_STARS = Dict(
        "snapNR"          => missings(Int64  , 0),
        "switch"        => missings(Int64  , 0),
        #"border"        => missings(Int64  , 0),
        #"jump"          => missings(Int64  , 0),
        #"exceed"        => missings(Int64  , 0),
        "FAKEFLIP"      => missings(Int64  , 0),
        "subID"         => missings(Int64  , 0),
        "N_MERGERS"     => missings(Int64  , 0),
        "ID_ISUB"       => missings(Int64  , 0),
        "ID_M"     => missings(Float64, 0),
        "ID_M2"         => missings(Float64, 0),
        "ID_Mpeak"         => missings(Float64, 0),
        "redshift"      => missings(Float64, 0),
        "lookbacktime"  => missings(Float64, 0),
        "M_MM"          => missings(Float64, 0),
        "M2_MM"         => missings(Float64, 0),
        "Mpeak_MM"         => missings(Float64, 0),
        "ΔM"      => missings(Float64, 0), 
        "ΔM2"     => missings(Float64, 0), 
        "ΔMpeak"     => missings(Float64, 0), 
        "ΔM_fromJ"      => missings(Float64, 0), 
        "M"       => missings(Float64, 0), 
        "M2"      => missings(Float64, 0), 
        "Mpeak"      => missings(Float64, 0), 
        "M_fromJ"       => missings(Float64, 0), 
        "ϕ_flip"        => missings(Float64, 0), 
        "M_MERGERS"     => missings(Float64, 0),
        "M_MISSED"      => missings(Float64, 0), 
        "M_CONSIDERED"  => missings(Float64, 0),  
        "M2_MERGERS"    => missings(Float64, 0),
        "M2_MISSED"     => missings(Float64, 0), 
        "M2_CONSIDERED" => missings(Float64, 0),  
        "Mpeak_MERGERS"    => missings(Float64, 0),
        "Mpeak_MISSED"     => missings(Float64, 0), 
        "Mpeak_CONSIDERED" => missings(Float64, 0),  
        "BVAL"          => missings(Float64, 0), 
        "BVAL_0"          => missings(Float64, 0), 
        "ΔBVAL"         => missings(Float64, 0), 
        "ΔBVAL_0"         => missings(Float64, 0), 
        
        "Δlookbacktime"         => missings(Float64, 0), 

        "SFR"           => missings(Float64 , 0), 
        "M200"          => missings(Float64 , 0), 
        "MCRI"          => missings(Float64 , 0), 
        "M500"          => missings(Float64 , 0), 
        "MVIR"          => missings(Float64 , 0), 
        "SAGE"          => missings(Float64 , 0), 
        "DSUB"          => missings(Float64 , 0), 
        "VMAX"          => missings(Float64 , 0), 
        "RMAX"          => missings(Float64 , 0), 
        "RVIR"          => missings(Float64 , 0), 
        "SPOS"          => missings(Float64, 3, 0),
        "SVEL"          => missings(Float64, 3, 0),
        "SPIN"          => missings(Float64, 3, 0),

        "J_MM"   => missings(Float64, 12, 0), 
        "J_SUM"  => missings(Float64, 12, 0), 
        "ΔJ_main"       => missings(Float64, 3, 0), 
        "Δj_main"       => missings(Float64, 3, 0),
        "J_main"        => missings(Float64, 3, 0), 
        "j_main"        => missings(Float64, 3, 0), 
        "J_vir"        => missings(Float64, 3, 0), 
        "j_vir"        => missings(Float64, 3, 0), 
        "ΔJ_vir"       => missings(Float64, 3, 0), 
        "Δj_vir"       => missings(Float64, 3, 0),
        "merger map"    => missings(Float64, 29, 0),

        "J_MM_0"   => missings(Float64, 12, 0), 
        "J_SUM_0"  => missings(Float64, 12, 0), 
        "ΔJ_main_0"       => missings(Float64, 3, 0), 
        "Δj_main_0"       => missings(Float64, 3, 0),
        "J_main_0"        => missings(Float64, 3, 0), 
        "j_main_0"        => missings(Float64, 3, 0), 
        "J_vir_0"        => missings(Float64, 3, 0), 
        "j_vir_0"        => missings(Float64, 3, 0), 
        "ΔJ_vir_0"       => missings(Float64, 3, 0), 
        "Δj_vir_0"       => missings(Float64, 3, 0),
        "merger map 0"    => missings(Float64, 29, 0))
    
    assembly_GAS = Dict(
        "snapNR"          => missings(Int64  , 0),
        "switch"        => missings(Int64  , 0),
        #"border"        => missings(Int64  , 0),
        #"jump"          => missings(Int64  , 0),
        #"exceed"        => missings(Int64  , 0),
        "FAKEFLIP"      => missings(Int64  , 0),
        "subID"         => missings(Int64  , 0),
        "N_MERGERS"     => missings(Int64  , 0),
        "ID_ISUB"       => missings(Int64  , 0),
        "ID_M"     => missings(Float64, 0),
        "ID_M2"         => missings(Float64, 0),
        "ID_Mpeak"         => missings(Float64, 0),
        "redshift"      => missings(Float64, 0),
        "lookbacktime"  => missings(Float64, 0),
        "M_MM"          => missings(Float64, 0),
        "M2_MM"         => missings(Float64, 0),
        "Mpeak_MM"         => missings(Float64, 0),
        "ΔM"      => missings(Float64, 0), 
        "ΔM2"     => missings(Float64, 0), 
        "ΔMpeak"     => missings(Float64, 0), 
        "ΔM_fromJ"      => missings(Float64, 0), 
        "M"       => missings(Float64, 0), 
        "M2"      => missings(Float64, 0), 
        "Mpeak"      => missings(Float64, 0), 
        "M_fromJ"       => missings(Float64, 0), 
        "ϕ_flip"        => missings(Float64, 0), 
        "M_MERGERS"     => missings(Float64, 0),
        "M_MISSED"      => missings(Float64, 0), 
        "M_CONSIDERED"  => missings(Float64, 0),  
        "M2_MERGERS"    => missings(Float64, 0),
        "M2_MISSED"     => missings(Float64, 0), 
        "M2_CONSIDERED" => missings(Float64, 0),  
        "Mpeak_MERGERS"    => missings(Float64, 0),
        "Mpeak_MISSED"     => missings(Float64, 0), 
        "Mpeak_CONSIDERED" => missings(Float64, 0),  
        "J_MM"   => missings(Float64, 12, 0), 
        "J_SUM"  => missings(Float64, 12, 0), 
        "ΔJ_main"       => missings(Float64, 3, 0), 
        "J_main"        => missings(Float64, 3, 0), 
        "j_main"        => missings(Float64, 3, 0), 
        "Δj_main"       => missings(Float64, 3, 0),
        "J_vir"        => missings(Float64, 3, 0), 
        "j_vir"        => missings(Float64, 3, 0), 
        "ΔJ_vir"       => missings(Float64, 3, 0), 
        "Δj_vir"       => missings(Float64, 3, 0),
        "merger map"    => missings(Float64, 29, 0),

        "J_MM_0"   => missings(Float64, 12, 0), 
        "J_SUM_0"  => missings(Float64, 12, 0), 
        "ΔJ_main_0"       => missings(Float64, 3, 0), 
        "Δj_main_0"       => missings(Float64, 3, 0),
        "J_main_0"        => missings(Float64, 3, 0), 
        "j_main_0"        => missings(Float64, 3, 0), 
        "J_vir_0"        => missings(Float64, 3, 0), 
        "j_vir_0"        => missings(Float64, 3, 0), 
        "ΔJ_vir_0"       => missings(Float64, 3, 0), 
        "Δj_vir_0"       => missings(Float64, 3, 0),
        "merger map 0"    => missings(Float64, 29, 0))
    
    assembly_DM = Dict(
        "snapNR"          => missings(Int64  , 0),
        "switch"        => missings(Int64  , 0),
        #"border"        => missings(Int64  , 0),
        #"jump"          => missings(Int64  , 0),
        #"exceed"        => missings(Int64  , 0),
        "FAKEFLIP"      => missings(Int64  , 0),
        "subID"         => missings(Int64  , 0),
        "N_MERGERS"     => missings(Int64  , 0),
        "ID_ISUB"       => missings(Int64  , 0),
        "ID_M"     => missings(Float64, 0),
        "ID_M2"         => missings(Float64, 0),
        "ID_Mpeak"         => missings(Float64, 0),
        "redshift"      => missings(Float64, 0),
        "lookbacktime"  => missings(Float64, 0),
        "M_MM"          => missings(Float64, 0),
        "M2_MM"         => missings(Float64, 0), 
        "Mpeak_MM"         => missings(Float64, 0), 
        "ΔM"      => missings(Float64, 0), 
        "ΔM2"     => missings(Float64, 0), 
        "ΔMpeak"     => missings(Float64, 0), 
        "ΔM_fromJ"      => missings(Float64, 0), 
        "M"       => missings(Float64, 0), 
        "M2"      => missings(Float64, 0), 
        "Mpeak"      => missings(Float64, 0), 
        "M_fromJ"       => missings(Float64, 0), 
        "ϕ_flip"        => missings(Float64, 0), 
        "M_MERGERS"     => missings(Float64, 0), 
        "M_MISSED"      => missings(Float64, 0), 
        "M_CONSIDERED"  => missings(Float64, 0), 
        "M2_MERGERS"    => missings(Float64, 0), 
        "M2_MISSED"     => missings(Float64, 0), 
        "M2_CONSIDERED" => missings(Float64, 0),  
        "Mpeak_MERGERS"    => missings(Float64, 0), 
        "Mpeak_MISSED"     => missings(Float64, 0), 
        "Mpeak_CONSIDERED" => missings(Float64, 0),  
        "J_MM"   => missings(Float64, 12, 0), 
        "J_SUM"  => missings(Float64, 12, 0), 
        "ΔJ_main"       => missings(Float64, 3, 0), 
        "J_main"        => missings(Float64, 3, 0), 
        "j_main"        => missings(Float64, 3, 0), 
        "Δj_main"       => missings(Float64, 3, 0),
        "merger map"    => missings(Float64, 29, 0),

        "J_MM_0"   => missings(Float64, 12, 0), 
        "J_SUM_0"  => missings(Float64, 12, 0), 
        "ΔJ_main_0"       => missings(Float64, 3, 0), 
        "Δj_main_0"       => missings(Float64, 3, 0),
        "J_main_0"        => missings(Float64, 3, 0), 
        "j_main_0"        => missings(Float64, 3, 0), 
        "merger map 0"    => missings(Float64, 29, 0))
    
    storyfilelist   = readdir(indir)
    noncentralsSUM  = 0
    centralsSUM     = 0
    if verbose
        println("   Scanning through $(length(storyfilelist)) stories.")
        flush(stdout)
    end
    for ii in 1:length(storyfilelist)
        if verbose
            print("$ii", ifelse(ii % 30 == 0, "\n", " "))
            flush(stdout)
        end
        merger_collection_STARS = load(joinpath(indir, storyfilelist[ii]), "merger_collection_STARS")
        merger_collection_GAS   = load(joinpath(indir, storyfilelist[ii]), "merger_collection_GAS")
        merger_collection_DM    = load(joinpath(indir, storyfilelist[ii]), "merger_collection_DM")

        # Condition of ending up as central halo in group (central_only)
        if central_only && get_first_subhalo(get_group(Galaxy(Snapshot(simbox, merger_collection_STARS["snapNR"][end]), merger_collection_STARS["subID"][end]))).isub != merger_collection_STARS["subID"][end]
            noncentralsSUM += 1
        else
            centralsSUM    += 1
            for i in 1:length(merger_collection_STARS["snapNR"])
                if merger_collection_STARS["M"][i] > mass_ST_thr
                    if length( findcs( assembly_STARS["snapNR"], eq=merger_collection_STARS["snapNR"][i], comparewith=findcs( assembly_STARS["subID"], eq=merger_collection_STARS["subID"][i] )) ) == 0
                        assembly_STARS["FAKEFLIP"     ]  = vcat( assembly_STARS["FAKEFLIP"   ], fake_flip_finder(merger_collection_STARS, i) )
                        assembly_STARS["snapNR"         ] = vcat( assembly_STARS["snapNR"         ], merger_collection_STARS["snapNR"         ][i] )
                        assembly_STARS["subID"        ] = vcat( assembly_STARS["subID"        ], merger_collection_STARS["subID"        ][i] )
                        assembly_STARS["switch"       ] = vcat( assembly_STARS["switch"       ], merger_collection_STARS["switch"       ][i] )
                        #assembly_STARS["border"       ] = vcat( assembly_STARS["border"       ], merger_collection_STARS["border"       ][i] )
                        #assembly_STARS["jump"         ] = vcat( assembly_STARS["jump"         ], merger_collection_STARS["jump"         ][i] )
                        #assembly_STARS["exceed"       ] = vcat( assembly_STARS["exceed"       ], merger_collection_STARS["exceed"       ][i] )
                        assembly_STARS["ID_ISUB"      ] = vcat( assembly_STARS["ID_ISUB"      ], merger_collection_STARS["ID_ISUB"      ][i] )
                        assembly_STARS["ID_M"    ] = vcat( assembly_STARS["ID_M"    ], merger_collection_STARS["ID_M"    ][i] )
                        assembly_STARS["ID_M2"        ] = vcat( assembly_STARS["ID_M2"        ], merger_collection_STARS["ID_M2"        ][i] )
                        assembly_STARS["ID_Mpeak"        ] = vcat( assembly_STARS["ID_Mpeak"        ], merger_collection_STARS["ID_Mpeak"        ][i] )
                        assembly_STARS["redshift"     ] = vcat( assembly_STARS["redshift"     ], merger_collection_STARS["redshift"     ][i] )
                        assembly_STARS["lookbacktime" ] = vcat( assembly_STARS["lookbacktime" ], merger_collection_STARS["lookbacktime" ][i] )
                        assembly_STARS["M_MM"         ] = vcat( assembly_STARS["M_MM"         ], merger_collection_STARS["M_MM"         ][i] )
                        assembly_STARS["M2_MM"        ] = vcat( assembly_STARS["M2_MM"        ], merger_collection_STARS["M2_MM"        ][i] )
                        assembly_STARS["Mpeak_MM"        ] = vcat( assembly_STARS["Mpeak_MM"        ], merger_collection_STARS["Mpeak_MM"        ][i] )
                        assembly_STARS["ΔM"     ] = vcat( assembly_STARS["ΔM"     ], merger_collection_STARS["ΔM"     ][i] )
                        assembly_STARS["ΔM2"    ] = vcat( assembly_STARS["ΔM2"    ], merger_collection_STARS["ΔM2"    ][i] )
                        assembly_STARS["ΔMpeak"    ] = vcat( assembly_STARS["ΔMpeak"    ], merger_collection_STARS["ΔMpeak"    ][i] )
                        assembly_STARS["ΔM_fromJ"     ] = vcat( assembly_STARS["ΔM_fromJ"     ], merger_collection_STARS["ΔM_fromJ"     ][i] )
                        assembly_STARS["M"      ] = vcat( assembly_STARS["M"      ], merger_collection_STARS["M"      ][i] )
                        assembly_STARS["M2"     ] = vcat( assembly_STARS["M2"     ], merger_collection_STARS["M2"     ][i] )
                        assembly_STARS["Mpeak"     ] = vcat( assembly_STARS["Mpeak"     ], merger_collection_STARS["Mpeak"     ][i] )
                        assembly_STARS["M_fromJ"      ] = vcat( assembly_STARS["M_fromJ"      ], merger_collection_STARS["M_fromJ"      ][i] )
                        assembly_STARS["ϕ_flip"       ] = vcat( assembly_STARS["ϕ_flip"       ], merger_collection_STARS["ϕ_flip"       ][i] )
                        assembly_STARS["N_MERGERS"    ] = vcat( assembly_STARS["N_MERGERS"    ], merger_collection_STARS["N_MERGERS"    ][i] )
                        assembly_STARS["M_MERGERS"    ] = vcat( assembly_STARS["M_MERGERS"    ], merger_collection_STARS["M_MERGERS"    ][i] )
                        assembly_STARS["M_MISSED"     ] = vcat( assembly_STARS["M_MISSED"     ], merger_collection_STARS["M_MISSED"     ][i] )
                        assembly_STARS["M_CONSIDERED" ] = vcat( assembly_STARS["M_CONSIDERED" ], merger_collection_STARS["M_CONSIDERED" ][i] )
                        assembly_STARS["M2_MERGERS"   ] = vcat( assembly_STARS["M2_MERGERS"   ], merger_collection_STARS["M2_MERGERS"   ][i] )
                        assembly_STARS["M2_MISSED"    ] = vcat( assembly_STARS["M2_MISSED"    ], merger_collection_STARS["M2_MISSED"    ][i] )
                        assembly_STARS["M2_CONSIDERED"] = vcat( assembly_STARS["M2_CONSIDERED"], merger_collection_STARS["M2_CONSIDERED"][i] )
                        assembly_STARS["Mpeak_MERGERS"   ] = vcat( assembly_STARS["Mpeak_MERGERS"   ], merger_collection_STARS["Mpeak_MERGERS"   ][i] )
                        assembly_STARS["Mpeak_MISSED"    ] = vcat( assembly_STARS["Mpeak_MISSED"    ], merger_collection_STARS["Mpeak_MISSED"    ][i] )
                        assembly_STARS["Mpeak_CONSIDERED"] = vcat( assembly_STARS["Mpeak_CONSIDERED"], merger_collection_STARS["Mpeak_CONSIDERED"][i] )
                        assembly_STARS["BVAL"         ] = vcat( assembly_STARS["BVAL"         ], merger_collection_STARS["BVAL"         ][i] )
                        assembly_STARS["ΔBVAL"        ] = vcat( assembly_STARS["ΔBVAL"        ], merger_collection_STARS["ΔBVAL"        ][i] )
                        assembly_STARS["Δlookbacktime"] = vcat( assembly_STARS["Δlookbacktime"], merger_collection_STARS["Δlookbacktime"][i] )
                        assembly_STARS["BVAL_0"       ] = vcat( assembly_STARS["BVAL_0"       ], merger_collection_STARS["BVAL_0"       ][i] )
                        assembly_STARS["ΔBVAL_0"      ] = vcat( assembly_STARS["ΔBVAL_0"      ], merger_collection_STARS["ΔBVAL_0"      ][i] )
                        assembly_STARS["SFR"          ] = vcat( assembly_STARS["SFR"          ], merger_collection_STARS["SFR"          ][i] )
                        assembly_STARS["M200"         ] = vcat( assembly_STARS["M200"         ], merger_collection_STARS["M200"         ][i] )
                        assembly_STARS["MCRI"         ] = vcat( assembly_STARS["MCRI"         ], merger_collection_STARS["MCRI"         ][i] )
                        assembly_STARS["M500"         ] = vcat( assembly_STARS["M500"         ], merger_collection_STARS["M500"         ][i] )
                        assembly_STARS["MVIR"         ] = vcat( assembly_STARS["MVIR"         ], merger_collection_STARS["MVIR"         ][i] )
                        assembly_STARS["SAGE"         ] = vcat( assembly_STARS["SAGE"         ], merger_collection_STARS["SAGE"         ][i] )
                        assembly_STARS["DSUB"         ] = vcat( assembly_STARS["DSUB"         ], merger_collection_STARS["DSUB"         ][i] )
                        assembly_STARS["VMAX"         ] = vcat( assembly_STARS["VMAX"         ], merger_collection_STARS["VMAX"         ][i] )
                        assembly_STARS["RMAX"         ] = vcat( assembly_STARS["RMAX"         ], merger_collection_STARS["RMAX"         ][i] )
                        assembly_STARS["RVIR"         ] = vcat( assembly_STARS["RVIR"         ], merger_collection_STARS["RVIR"         ][i] )
                        assembly_STARS["SPOS"         ] = hcat( assembly_STARS["SPOS"         ], merger_collection_STARS["SPOS"         ][:,i] )
                        assembly_STARS["SVEL"         ] = hcat( assembly_STARS["SVEL"         ], merger_collection_STARS["SVEL"         ][:,i] )
                        assembly_STARS["SPIN"         ] = hcat( assembly_STARS["SPIN"         ], merger_collection_STARS["SPIN"         ][:,i] )
                        assembly_STARS["J_MM"  ] = hcat( assembly_STARS["J_MM"  ], merger_collection_STARS["J_MM"  ][:,i] )
                        assembly_STARS["J_SUM" ] = hcat( assembly_STARS["J_SUM" ], merger_collection_STARS["J_SUM" ][:,i] )
                        assembly_STARS["ΔJ_main"      ] = hcat( assembly_STARS["ΔJ_main"      ], merger_collection_STARS["ΔJ_main"      ][:,i] )
                        assembly_STARS["Δj_main"      ] = hcat( assembly_STARS["Δj_main"      ], merger_collection_STARS["Δj_main"      ][:,i] )
                        assembly_STARS["J_main"       ] = hcat( assembly_STARS["J_main"       ], merger_collection_STARS["J_main"       ][:,i] )
                        assembly_STARS["j_main"       ] = hcat( assembly_STARS["j_main"       ], merger_collection_STARS["j_main"       ][:,i] )
                        assembly_STARS["J_vir"        ] = hcat( assembly_STARS["J_vir"        ], merger_collection_STARS["J_vir"        ][:,i] )
                        assembly_STARS["j_vir"        ] = hcat( assembly_STARS["j_vir"        ], merger_collection_STARS["j_vir"        ][:,i] )
                        assembly_STARS["merger map"   ] = hcat( assembly_STARS["merger map"   ], merger_collection_STARS["merger map"   ][:,1+sum(merger_collection_STARS["N_MERGERS"][1:i-1]):sum(merger_collection_STARS["N_MERGERS"][1:i])] )
                        #println("$(storyfilelist[ii])   $i   $(merger_collection_STARS["merger map"][:,i])")
                        assembly_STARS["ΔJ_vir"       ] = hcat( assembly_STARS["ΔJ_vir"      ], merger_collection_STARS["ΔJ_vir"        ][:,i] )
                        assembly_STARS["Δj_vir"       ] = hcat( assembly_STARS["Δj_vir"      ], merger_collection_STARS["Δj_vir"        ][:,i] )
                        assembly_STARS["J_MM_0"  ] = hcat( assembly_STARS["J_MM_0"  ], merger_collection_STARS["J_MM_0"  ][:,i] )
                        assembly_STARS["J_SUM_0" ] = hcat( assembly_STARS["J_SUM_0" ], merger_collection_STARS["J_SUM_0" ][:,i] )
                        assembly_STARS["ΔJ_main_0"      ] = hcat( assembly_STARS["ΔJ_main_0"      ], merger_collection_STARS["ΔJ_main_0"      ][:,i] )
                        assembly_STARS["Δj_main_0"      ] = hcat( assembly_STARS["Δj_main_0"      ], merger_collection_STARS["Δj_main_0"      ][:,i] )
                        assembly_STARS["J_main_0"       ] = hcat( assembly_STARS["J_main_0"       ], merger_collection_STARS["J_main_0"       ][:,i] )
                        assembly_STARS["j_main_0"       ] = hcat( assembly_STARS["j_main_0"       ], merger_collection_STARS["j_main_0"       ][:,i] )
                        assembly_STARS["J_vir_0"        ] = hcat( assembly_STARS["J_vir_0"        ], merger_collection_STARS["J_vir_0"        ][:,i] )
                        assembly_STARS["j_vir_0"        ] = hcat( assembly_STARS["j_vir_0"        ], merger_collection_STARS["j_vir_0"        ][:,i] )
                        assembly_STARS["merger map 0"   ] = hcat( assembly_STARS["merger map 0"   ], merger_collection_STARS["merger map 0"   ][:,1+sum(merger_collection_STARS["N_MERGERS"][1:i-1]):sum(merger_collection_STARS["N_MERGERS"][1:i])] )
                        assembly_STARS["ΔJ_vir_0"       ] = hcat( assembly_STARS["ΔJ_vir_0"      ], merger_collection_STARS["ΔJ_vir_0"        ][:,i] )
                        assembly_STARS["Δj_vir_0"       ] = hcat( assembly_STARS["Δj_vir_0"      ], merger_collection_STARS["Δj_vir_0"        ][:,i] )
                        assembly_GAS["FAKEFLIP"     ]   = vcat( assembly_GAS["FAKEFLIP"   ], fake_flip_finder(merger_collection_GAS, i) )
                        assembly_GAS["snapNR"         ]   = vcat( assembly_GAS["snapNR"         ], merger_collection_GAS["snapNR"         ][i] )
                        assembly_GAS["switch"       ]   = vcat( assembly_GAS["switch"       ], merger_collection_GAS["switch"       ][i] )
                        #assembly_GAS["border"       ]   = vcat( assembly_GAS["border"       ], merger_collection_GAS["border"       ][i] )
                        #assembly_GAS["jump"         ]   = vcat( assembly_GAS["jump"         ], merger_collection_GAS["jump"         ][i] )
                        #assembly_GAS["exceed"       ]   = vcat( assembly_GAS["exceed"       ], merger_collection_GAS["exceed"       ][i] )
                        assembly_GAS["subID"        ]   = vcat( assembly_GAS["subID"        ], merger_collection_GAS["subID"        ][i] )
                        assembly_GAS["ID_ISUB"      ]   = vcat( assembly_GAS["ID_ISUB"      ], merger_collection_GAS["ID_ISUB"      ][i] )
                        assembly_GAS["ID_M"    ]   = vcat( assembly_GAS["ID_M"    ], merger_collection_GAS["ID_M"    ][i] )
                        assembly_GAS["ID_M2"        ]   = vcat( assembly_GAS["ID_M2"        ], merger_collection_GAS["ID_M2"        ][i] )
                        assembly_GAS["ID_Mpeak"        ]   = vcat( assembly_GAS["ID_Mpeak"        ], merger_collection_GAS["ID_Mpeak"        ][i] )
                        assembly_GAS["redshift"     ]   = vcat( assembly_GAS["redshift"     ], merger_collection_GAS["redshift"     ][i] )
                        assembly_GAS["lookbacktime" ]   = vcat( assembly_GAS["lookbacktime" ], merger_collection_GAS["lookbacktime" ][i] )
                        assembly_GAS["M_MM"         ]   = vcat( assembly_GAS["M_MM"         ], merger_collection_GAS["M_MM"         ][i] )
                        assembly_GAS["M2_MM"        ]   = vcat( assembly_GAS["M2_MM"        ], merger_collection_GAS["M2_MM"        ][i] )
                        assembly_GAS["Mpeak_MM"        ]   = vcat( assembly_GAS["Mpeak_MM"        ], merger_collection_GAS["Mpeak_MM"        ][i] )
                        assembly_GAS["ΔM"     ]   = vcat( assembly_GAS["ΔM"     ], merger_collection_GAS["ΔM"     ][i] )
                        assembly_GAS["ΔM2"    ]   = vcat( assembly_GAS["ΔM2"    ], merger_collection_GAS["ΔM2"    ][i] )
                        assembly_GAS["ΔMpeak"    ]   = vcat( assembly_GAS["ΔMpeak"    ], merger_collection_GAS["ΔMpeak"    ][i] )
                        assembly_GAS["ΔM_fromJ"     ]   = vcat( assembly_GAS["ΔM_fromJ"     ], merger_collection_GAS["ΔM_fromJ"     ][i] )
                        assembly_GAS["M"      ]   = vcat( assembly_GAS["M"      ], merger_collection_GAS["M"      ][i] )
                        assembly_GAS["M2"     ]   = vcat( assembly_GAS["M2"     ], merger_collection_GAS["M2"     ][i] )
                        assembly_GAS["Mpeak"     ]   = vcat( assembly_GAS["Mpeak"     ], merger_collection_GAS["Mpeak"     ][i] )
                        assembly_GAS["M_fromJ"      ]   = vcat( assembly_GAS["M_fromJ"      ], merger_collection_GAS["M_fromJ"      ][i] )
                        assembly_GAS["ϕ_flip"       ]   = vcat( assembly_GAS["ϕ_flip"       ], merger_collection_GAS["ϕ_flip"       ][i] )
                        assembly_GAS["N_MERGERS"    ]   = vcat( assembly_GAS["N_MERGERS"    ], merger_collection_GAS["N_MERGERS"    ][i] )
                        assembly_GAS["M_MERGERS"    ]   = vcat( assembly_GAS["M_MERGERS"    ], merger_collection_GAS["M_MERGERS"    ][i] )
                        assembly_GAS["M_MISSED"     ]   = vcat( assembly_GAS["M_MISSED"     ], merger_collection_GAS["M_MISSED"     ][i] )
                        assembly_GAS["M_CONSIDERED" ]   = vcat( assembly_GAS["M_CONSIDERED" ], merger_collection_GAS["M_CONSIDERED" ][i] )
                        assembly_GAS["M2_MERGERS"   ]   = vcat( assembly_GAS["M2_MERGERS"   ], merger_collection_GAS["M2_MERGERS"   ][i] )
                        assembly_GAS["M2_MISSED"    ]   = vcat( assembly_GAS["M2_MISSED"    ], merger_collection_GAS["M2_MISSED"    ][i] )
                        assembly_GAS["M2_CONSIDERED"]   = vcat( assembly_GAS["M2_CONSIDERED"], merger_collection_GAS["M2_CONSIDERED"][i] )
                        assembly_GAS["Mpeak_MERGERS"   ]   = vcat( assembly_GAS["Mpeak_MERGERS"   ], merger_collection_GAS["Mpeak_MERGERS"   ][i] )
                        assembly_GAS["Mpeak_MISSED"    ]   = vcat( assembly_GAS["Mpeak_MISSED"    ], merger_collection_GAS["Mpeak_MISSED"    ][i] )
                        assembly_GAS["Mpeak_CONSIDERED"]   = vcat( assembly_GAS["Mpeak_CONSIDERED"], merger_collection_GAS["Mpeak_CONSIDERED"][i] )
                        assembly_GAS["J_MM"  ]   = hcat( assembly_GAS["J_MM"  ], merger_collection_GAS["J_MM"  ][:,i] )
                        assembly_GAS["J_SUM" ]   = hcat( assembly_GAS["J_SUM" ], merger_collection_GAS["J_SUM" ][:,i] )
                        assembly_GAS["ΔJ_main"      ]   = hcat( assembly_GAS["ΔJ_main"      ], merger_collection_GAS["ΔJ_main"      ][:,i] )
                        assembly_GAS["J_main"       ]   = hcat( assembly_GAS["J_main"       ], merger_collection_GAS["J_main"       ][:,i] )
                        assembly_GAS["j_main"       ]   = hcat( assembly_GAS["j_main"       ], merger_collection_GAS["j_main"       ][:,i] )
                        assembly_GAS["Δj_main"      ]   = hcat( assembly_GAS["Δj_main"      ], merger_collection_GAS["Δj_main"      ][:,i] )
                        assembly_GAS["J_vir"        ]   = hcat( assembly_GAS["J_vir"        ], merger_collection_GAS["J_vir"        ][:,i] )
                        assembly_GAS["j_vir"        ]   = hcat( assembly_GAS["j_vir"        ], merger_collection_GAS["j_vir"        ][:,i] )
                        assembly_GAS["ΔJ_vir"       ]   = hcat( assembly_GAS["ΔJ_vir"      ], merger_collection_GAS["ΔJ_vir"        ][:,i] )
                        assembly_GAS["Δj_vir"       ]   = hcat( assembly_GAS["Δj_vir"      ], merger_collection_GAS["Δj_vir"        ][:,i] )
                        assembly_GAS["merger map"   ]   = hcat( assembly_GAS["merger map"   ], merger_collection_GAS["merger map"   ][:,1+sum(merger_collection_GAS["N_MERGERS"][1:i-1]):sum(merger_collection_GAS["N_MERGERS"][1:i])] )
                        assembly_GAS["J_MM_0"  ] = hcat( assembly_GAS["J_MM_0"  ], merger_collection_GAS["J_MM_0"  ][:,i] )
                        assembly_GAS["J_SUM_0" ] = hcat( assembly_GAS["J_SUM_0" ], merger_collection_GAS["J_SUM_0" ][:,i] )
                        assembly_GAS["ΔJ_main_0"      ] = hcat( assembly_GAS["ΔJ_main_0"      ], merger_collection_GAS["ΔJ_main_0"      ][:,i] )
                        assembly_GAS["Δj_main_0"      ] = hcat( assembly_GAS["Δj_main_0"      ], merger_collection_GAS["Δj_main_0"      ][:,i] )
                        assembly_GAS["J_main_0"       ] = hcat( assembly_GAS["J_main_0"       ], merger_collection_GAS["J_main_0"       ][:,i] )
                        assembly_GAS["j_main_0"       ] = hcat( assembly_GAS["j_main_0"       ], merger_collection_GAS["j_main_0"       ][:,i] )
                        assembly_GAS["J_vir_0"        ] = hcat( assembly_GAS["J_vir_0"        ], merger_collection_GAS["J_vir_0"        ][:,i] )
                        assembly_GAS["j_vir_0"        ] = hcat( assembly_GAS["j_vir_0"        ], merger_collection_GAS["j_vir_0"        ][:,i] )
                        assembly_GAS["merger map 0"   ] = hcat( assembly_GAS["merger map 0"   ], merger_collection_GAS["merger map 0"   ][:,1+sum(merger_collection_GAS["N_MERGERS"][1:i-1]):sum(merger_collection_GAS["N_MERGERS"][1:i])] )
                        assembly_GAS["ΔJ_vir_0"       ] = hcat( assembly_GAS["ΔJ_vir_0"      ], merger_collection_GAS["ΔJ_vir_0"        ][:,i] )
                        assembly_GAS["Δj_vir_0"       ] = hcat( assembly_GAS["Δj_vir_0"      ], merger_collection_GAS["Δj_vir_0"        ][:,i] )
                        assembly_DM["FAKEFLIP"     ]    = vcat( assembly_DM["FAKEFLIP"   ], fake_flip_finder(merger_collection_DM, i) )
                        assembly_DM["snapNR"         ]    = vcat( assembly_DM["snapNR"         ], merger_collection_DM["snapNR"         ][i] )
                        assembly_DM["switch"       ]   = vcat( assembly_DM["switch"       ], merger_collection_DM["switch"       ][i] )
                        #assembly_DM["border"       ]   = vcat( assembly_DM["border"       ], merger_collection_GAS["border"       ][i] )
                        #assembly_DM["jump"         ]   = vcat( assembly_DM["jump"         ], merger_collection_GAS["jump"         ][i] )
                        #assembly_DM["exceed"       ]   = vcat( assembly_DM["exceed"       ], merger_collection_GAS["exceed"       ][i] )
                        assembly_DM["subID"        ]    = vcat( assembly_DM["subID"        ], merger_collection_DM["subID"        ][i] )
                        assembly_DM["ID_ISUB"      ]    = vcat( assembly_DM["ID_ISUB"      ], merger_collection_DM["ID_ISUB"      ][i] )
                        assembly_DM["ID_M"    ]    = vcat( assembly_DM["ID_M"    ], merger_collection_DM["ID_M"    ][i] )
                        assembly_DM["ID_M2"        ]    = vcat( assembly_DM["ID_M2"        ], merger_collection_DM["ID_M2"        ][i] )
                        assembly_DM["ID_Mpeak"        ]    = vcat( assembly_DM["ID_Mpeak"        ], merger_collection_DM["ID_Mpeak"        ][i] )
                        assembly_DM["redshift"     ]    = vcat( assembly_DM["redshift"     ], merger_collection_DM["redshift"     ][i] )
                        assembly_DM["lookbacktime" ]    = vcat( assembly_DM["lookbacktime" ], merger_collection_DM["lookbacktime" ][i] )
                        assembly_DM["M_MM"         ]    = vcat( assembly_DM["M_MM"         ], merger_collection_DM["M_MM"         ][i] )
                        assembly_DM["M2_MM"        ]    = vcat( assembly_DM["M2_MM"        ], merger_collection_DM["M2_MM"        ][i] )
                        assembly_DM["Mpeak_MM"        ]    = vcat( assembly_DM["Mpeak_MM"        ], merger_collection_DM["Mpeak_MM"        ][i] )
                        assembly_DM["ΔM"     ]    = vcat( assembly_DM["ΔM"     ], merger_collection_DM["ΔM"     ][i] )
                        assembly_DM["ΔM2"    ]    = vcat( assembly_DM["ΔM2"    ], merger_collection_DM["ΔM2"    ][i] )
                        assembly_DM["ΔMpeak"    ]    = vcat( assembly_DM["ΔMpeak"    ], merger_collection_DM["ΔMpeak"    ][i] )
                        assembly_DM["ΔM_fromJ"     ]    = vcat( assembly_DM["ΔM_fromJ"     ], merger_collection_DM["ΔM_fromJ"     ][i] )
                        assembly_DM["M"      ]    = vcat( assembly_DM["M"      ], merger_collection_DM["M"      ][i] )
                        assembly_DM["M2"     ]    = vcat( assembly_DM["M2"     ], merger_collection_DM["M2"     ][i] )
                        assembly_DM["Mpeak"     ]    = vcat( assembly_DM["Mpeak"     ], merger_collection_DM["Mpeak"     ][i] )
                        assembly_DM["M_fromJ"      ]    = vcat( assembly_DM["M_fromJ"      ], merger_collection_DM["M_fromJ"      ][i] )
                        assembly_DM["ϕ_flip"       ]    = vcat( assembly_DM["ϕ_flip"       ], merger_collection_DM["ϕ_flip"       ][i] )
                        assembly_DM["N_MERGERS"    ]    = vcat( assembly_DM["N_MERGERS"    ], merger_collection_DM["N_MERGERS"    ][i] )
                        assembly_DM["M_MERGERS"    ]    = vcat( assembly_DM["M_MERGERS"    ], merger_collection_DM["M_MERGERS"    ][i] )
                        assembly_DM["M_MISSED"     ]    = vcat( assembly_DM["M_MISSED"     ], merger_collection_DM["M_MISSED"     ][i] )
                        assembly_DM["M_CONSIDERED" ]    = vcat( assembly_DM["M_CONSIDERED" ], merger_collection_DM["M_CONSIDERED" ][i] )
                        assembly_DM["M2_MERGERS"   ]    = vcat( assembly_DM["M2_MERGERS"   ], merger_collection_DM["M2_MERGERS"   ][i] )
                        assembly_DM["M2_MISSED"    ]    = vcat( assembly_DM["M2_MISSED"    ], merger_collection_DM["M2_MISSED"    ][i] )
                        assembly_DM["M2_CONSIDERED"]    = vcat( assembly_DM["M2_CONSIDERED"], merger_collection_DM["M2_CONSIDERED"][i] )
                        assembly_DM["Mpeak_MERGERS"   ]    = vcat( assembly_DM["Mpeak_MERGERS"   ], merger_collection_DM["Mpeak_MERGERS"   ][i] )
                        assembly_DM["Mpeak_MISSED"    ]    = vcat( assembly_DM["Mpeak_MISSED"    ], merger_collection_DM["Mpeak_MISSED"    ][i] )
                        assembly_DM["Mpeak_CONSIDERED"]    = vcat( assembly_DM["Mpeak_CONSIDERED"], merger_collection_DM["Mpeak_CONSIDERED"][i] )
                        assembly_DM["J_MM"  ]    = hcat( assembly_DM["J_MM"  ], merger_collection_DM["J_MM"  ][:,i] )
                        assembly_DM["J_SUM" ]    = hcat( assembly_DM["J_SUM" ], merger_collection_DM["J_SUM" ][:,i] )
                        assembly_DM["ΔJ_main"      ]    = hcat( assembly_DM["ΔJ_main"      ], merger_collection_DM["ΔJ_main"      ][:,i] )
                        assembly_DM["J_main"       ]    = hcat( assembly_DM["J_main"       ], merger_collection_DM["J_main"       ][:,i] )
                        assembly_DM["j_main"       ]    = hcat( assembly_DM["j_main"       ], merger_collection_DM["j_main"       ][:,i] )
                        assembly_DM["Δj_main"      ]    = hcat( assembly_DM["Δj_main"      ], merger_collection_DM["Δj_main"      ][:,i] )
                        assembly_DM["merger map"   ]    = hcat( assembly_DM["merger map"   ], merger_collection_DM["merger map"   ][:,1+sum(merger_collection_DM["N_MERGERS"][1:i-1]):sum(merger_collection_DM["N_MERGERS"][1:i])] )
                        assembly_DM["J_MM_0"  ] = hcat( assembly_DM["J_MM_0"  ], merger_collection_DM["J_MM_0"  ][:,i] )
                        assembly_DM["J_SUM_0" ] = hcat( assembly_DM["J_SUM_0" ], merger_collection_DM["J_SUM_0" ][:,i] )
                        assembly_DM["ΔJ_main_0"      ] = hcat( assembly_DM["ΔJ_main_0"      ], merger_collection_DM["ΔJ_main_0"      ][:,i] )
                        assembly_DM["Δj_main_0"      ] = hcat( assembly_DM["Δj_main_0"      ], merger_collection_DM["Δj_main_0"      ][:,i] )
                        assembly_DM["J_main_0"       ] = hcat( assembly_DM["J_main_0"       ], merger_collection_DM["J_main_0"       ][:,i] )
                        assembly_DM["j_main_0"       ] = hcat( assembly_DM["j_main_0"       ], merger_collection_DM["j_main_0"       ][:,i] )
                        assembly_DM["merger map 0"   ] = hcat( assembly_DM["merger map 0"   ], merger_collection_DM["merger map 0"   ][:,1+sum(merger_collection_DM["N_MERGERS"][1:i-1]):sum(merger_collection_DM["N_MERGERS"][1:i])] )
                    elseif verbose
                        print("d", ifelse(i == length(merger_collection_STARS["snapNR"]), " ", ""))
                        flush(stdout)
                    end
                end
            end
        end
    end
    
    save(joinpath(outdir, outfile), 
        "assembly_STARS",   assembly_STARS,
        "assembly_DM",      assembly_DM,
        "assembly_GAS",     assembly_GAS)

    if verbose
        ifelse( central_only, println("\n\n$centralsSUM Centrals and $noncentralsSUM Non-Centrals"), println("\n") )
        println("J_main size for stars$(size(assembly_STARS["J_main"]))")
        println("J_main size for gas  $(size(assembly_GAS["J_main"]))")
        println("J_main size for dm   $(size(assembly_DM["J_main"]))")
        println("Merger map & N_MERGERS size for stars$(size(assembly_STARS["merger map"]))   =   $(sum(assembly_STARS["N_MERGERS"]))   ?")
        println("Merger map & N_MERGERS size for gas  $(size(assembly_GAS["merger map"]))   =   $(sum(assembly_GAS["N_MERGERS"]))   ?")
        println("Merger map & N_MERGERS size for dm   $(size(assembly_DM["merger map"]))   =   $(sum(assembly_DM["N_MERGERS"]))   ?")
    end
    
    println("\n\n\n---------------------------\n\nNow witness the firepower of this fully armed and operational battle station!\n\n")
    return nothing
end

print("'assemble_halostories'   ")
