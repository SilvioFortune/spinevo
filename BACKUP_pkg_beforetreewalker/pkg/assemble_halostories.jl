
@doc """
DESCRIPTION:\n
INPUT:\n 
OUTPUT:\n
""" -> 
function assemble_halostories(; mass_ST_thr=1e10, outfile=replace("assembly_Mstar_$(mass_ST_thr)_$(today()).jld", "-" => ""), central_only=false,
    outdir="./OUT_assemble_halostories", indir=current_dir_jld,
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
        "SNAP"          => missings(Int64  , 0),
        "SWITCH"        => missings(Int64  , 0),
        "BORDER"        => missings(Int64  , 0),
        "JUMP"          => missings(Int64  , 0),
        "EXCEED"        => missings(Int64  , 0),
        "FAKEFLIP"      => missings(Int64  , 0),
        "I_SUB"         => missings(Int64  , 0),
        "N_MERGERS"     => missings(Int64  , 0),
        "ID_ISUB"       => missings(Int64  , 0),
        "ID_Mfelix"     => missings(Float64, 0),
        "ID_M2"         => missings(Float64, 0),
        "REDSHIFT"      => missings(Float64, 0),
        "LOOKBACKTIME"  => missings(Float64, 0),
        "M_MM"          => missings(Float64, 0),
        "M2_MM"         => missings(Float64, 0),
        "ΔM_felix"      => missings(Float64, 0), 
        "ΔM2_felix"     => missings(Float64, 0), 
        "ΔM_fromJ"      => missings(Float64, 0), 
        "M_felix"       => missings(Float64, 0), 
        "M2_felix"      => missings(Float64, 0), 
        "M_fromJ"       => missings(Float64, 0), 
        "ϕ_flip"        => missings(Float64, 0), 
        "M_MERGERS"     => missings(Float64, 0),
        "M_MISSED"      => missings(Float64, 0), 
        "M_CONSIDERED"  => missings(Float64, 0),  
        "M2_MERGERS"    => missings(Float64, 0),
        "M2_MISSED"     => missings(Float64, 0), 
        "M2_CONSIDERED" => missings(Float64, 0),  
        "BVAL"          => missings(Float64, 0), 
        "BVAL_0"          => missings(Float64, 0), 
        "ΔBVAL"         => missings(Float64, 0), 
        "ΔBVAL_0"         => missings(Float64, 0), 
        "SFR"           => missings(Float64 , 0), 
        "sSFR"          => missings(Float64 , 0), 
        "J_MMorbital"   => missings(Float64, 3, 0), 
        "J_SUMorbital"  => missings(Float64, 3, 0), 
        "ΔJ_main"       => missings(Float64, 3, 0), 
        "Δj_main"       => missings(Float64, 3, 0),
        "J_main"        => missings(Float64, 3, 0), 
        "j_main"        => missings(Float64, 3, 0), 
        "J_vir"        => missings(Float64, 3, 0), 
        "j_vir"        => missings(Float64, 3, 0), 
        "ΔJ_vir"       => missings(Float64, 3, 0), 
        "Δj_vir"       => missings(Float64, 3, 0),
        "Merger_Map"    => missings(Float64, 8, 0))
    
    assembly_GAS = Dict(
        "SNAP"          => missings(Int64  , 0),
        "SWITCH"        => missings(Int64  , 0),
        "BORDER"        => missings(Int64  , 0),
        "JUMP"          => missings(Int64  , 0),
        "EXCEED"        => missings(Int64  , 0),
        "FAKEFLIP"      => missings(Int64  , 0),
        "I_SUB"         => missings(Int64  , 0),
        "N_MERGERS"     => missings(Int64  , 0),
        "ID_ISUB"       => missings(Int64  , 0),
        "ID_Mfelix"     => missings(Float64, 0),
        "ID_M2"         => missings(Float64, 0),
        "REDSHIFT"      => missings(Float64, 0),
        "LOOKBACKTIME"  => missings(Float64, 0),
        "M_MM"          => missings(Float64, 0),
        "M2_MM"         => missings(Float64, 0),
        "ΔM_felix"      => missings(Float64, 0), 
        "ΔM2_felix"     => missings(Float64, 0), 
        "ΔM_fromJ"      => missings(Float64, 0), 
        "M_felix"       => missings(Float64, 0), 
        "M2_felix"      => missings(Float64, 0), 
        "M_fromJ"       => missings(Float64, 0), 
        "ϕ_flip"        => missings(Float64, 0), 
        "M_MERGERS"     => missings(Float64, 0),
        "M_MISSED"      => missings(Float64, 0), 
        "M_CONSIDERED"  => missings(Float64, 0),  
        "M2_MERGERS"    => missings(Float64, 0),
        "M2_MISSED"     => missings(Float64, 0), 
        "M2_CONSIDERED" => missings(Float64, 0),  
        "J_MMorbital"   => missings(Float64, 3, 0), 
        "J_SUMorbital"  => missings(Float64, 3, 0), 
        "ΔJ_main"       => missings(Float64, 3, 0), 
        "J_main"        => missings(Float64, 3, 0), 
        "j_main"        => missings(Float64, 3, 0), 
        "Δj_main"       => missings(Float64, 3, 0),
        "J_vir"        => missings(Float64, 3, 0), 
        "j_vir"        => missings(Float64, 3, 0), 
        "ΔJ_vir"       => missings(Float64, 3, 0), 
        "Δj_vir"       => missings(Float64, 3, 0),
        "Merger_Map"    => missings(Float64, 8, 0))
    
    assembly_DM = Dict(
        "SNAP"          => missings(Int64  , 0),
        "SWITCH"        => missings(Int64  , 0),
        "BORDER"        => missings(Int64  , 0),
        "JUMP"          => missings(Int64  , 0),
        "EXCEED"        => missings(Int64  , 0),
        "FAKEFLIP"      => missings(Int64  , 0),
        "I_SUB"         => missings(Int64  , 0),
        "N_MERGERS"     => missings(Int64  , 0),
        "ID_ISUB"       => missings(Int64  , 0),
        "ID_Mfelix"     => missings(Float64, 0),
        "ID_M2"         => missings(Float64, 0),
        "REDSHIFT"      => missings(Float64, 0),
        "LOOKBACKTIME"  => missings(Float64, 0),
        "M_MM"          => missings(Float64, 0),
        "M2_MM"         => missings(Float64, 0), 
        "ΔM_felix"      => missings(Float64, 0), 
        "ΔM2_felix"     => missings(Float64, 0), 
        "ΔM_fromJ"      => missings(Float64, 0), 
        "M_felix"       => missings(Float64, 0), 
        "M2_felix"      => missings(Float64, 0), 
        "M_fromJ"       => missings(Float64, 0), 
        "ϕ_flip"        => missings(Float64, 0), 
        "M_MERGERS"     => missings(Float64, 0), 
        "M_MISSED"      => missings(Float64, 0), 
        "M_CONSIDERED"  => missings(Float64, 0), 
        "M2_MERGERS"    => missings(Float64, 0), 
        "M2_MISSED"     => missings(Float64, 0), 
        "M2_CONSIDERED" => missings(Float64, 0),  
        "J_MMorbital"   => missings(Float64, 3, 0), 
        "J_SUMorbital"  => missings(Float64, 3, 0), 
        "ΔJ_main"       => missings(Float64, 3, 0), 
        "J_main"        => missings(Float64, 3, 0), 
        "j_main"        => missings(Float64, 3, 0), 
        "Δj_main"       => missings(Float64, 3, 0),
        "Merger_Map"    => missings(Float64, 8, 0))
    
    storyfilelist   = readdir(indir)
    noncentralsSUM  = 0
    centralsSUM     = 0
    for ii in 1:length(storyfilelist)
        if verbose
            print("$ii", ifelse(ii % 20 == 0, "\n", " "))
            flush(stdout)
        end
        merger_collection_STARS = load(joinpath(indir, storyfilelist[ii]), "merger_collection_STARS")
        merger_collection_GAS   = load(joinpath(indir, storyfilelist[ii]), "merger_collection_GAS")
        merger_collection_DM    = load(joinpath(indir, storyfilelist[ii]), "merger_collection_DM")

        # Condition of ending up as central halo in group (central_only)
        if central_only && get_first_subhalo(get_group(Galaxy(Snapshot(simbox, merger_collection_STARS["SNAP"][end]), merger_collection_STARS["I_SUB"][end]))).isub != merger_collection_STARS["I_SUB"][end]
            noncentralsSUM += 1
        else
            centralsSUM    += 1
            for i in 1:length(merger_collection_STARS["SNAP"])
                #if merger_collection_STARS["M2_felix"][i] > mass_ST_thr
                if merger_collection_STARS["M_fromJ"][i] > mass_ST_thr
                    assembly_STARS["FAKEFLIP"     ]  = vcat( assembly_STARS["FAKEFLIP"   ], fake_flip_finder(merger_collection_STARS, i) )
                    assembly_STARS["SNAP"         ] = vcat( assembly_STARS["SNAP"         ], merger_collection_STARS["SNAP"         ][i] )
                    assembly_STARS["I_SUB"        ] = vcat( assembly_STARS["I_SUB"        ], merger_collection_STARS["I_SUB"        ][i] )
                    assembly_STARS["SWITCH"       ] = vcat( assembly_STARS["SWITCH"       ], merger_collection_STARS["SWITCH"       ][i] )
                    assembly_STARS["BORDER"       ] = vcat( assembly_STARS["BORDER"       ], merger_collection_STARS["BORDER"       ][i] )
                    assembly_STARS["JUMP"         ] = vcat( assembly_STARS["JUMP"         ], merger_collection_STARS["JUMP"         ][i] )
                    assembly_STARS["EXCEED"       ] = vcat( assembly_STARS["EXCEED"       ], merger_collection_STARS["EXCEED"       ][i] )
                    assembly_STARS["ID_ISUB"      ] = vcat( assembly_STARS["ID_ISUB"      ], merger_collection_STARS["ID_ISUB"      ][i] )
                    assembly_STARS["ID_Mfelix"    ] = vcat( assembly_STARS["ID_Mfelix"    ], merger_collection_STARS["ID_Mfelix"    ][i] )
                    assembly_STARS["ID_M2"        ] = vcat( assembly_STARS["ID_M2"        ], merger_collection_STARS["ID_M2"        ][i] )
                    assembly_STARS["REDSHIFT"     ] = vcat( assembly_STARS["REDSHIFT"     ], merger_collection_STARS["REDSHIFT"     ][i] )
                    assembly_STARS["LOOKBACKTIME" ] = vcat( assembly_STARS["LOOKBACKTIME" ], merger_collection_STARS["LOOKBACKTIME" ][i] )
                    assembly_STARS["M_MM"         ] = vcat( assembly_STARS["M_MM"         ], merger_collection_STARS["M_MM"         ][i] )
                    assembly_STARS["M2_MM"        ] = vcat( assembly_STARS["M2_MM"        ], merger_collection_STARS["M2_MM"        ][i] )
                    assembly_STARS["ΔM_felix"     ] = vcat( assembly_STARS["ΔM_felix"     ], merger_collection_STARS["ΔM_felix"     ][i] )
                    assembly_STARS["ΔM2_felix"    ] = vcat( assembly_STARS["ΔM2_felix"    ], merger_collection_STARS["ΔM2_felix"    ][i] )
                    assembly_STARS["ΔM_fromJ"     ] = vcat( assembly_STARS["ΔM_fromJ"     ], merger_collection_STARS["ΔM_fromJ"     ][i] )
                    assembly_STARS["M_felix"      ] = vcat( assembly_STARS["M_felix"      ], merger_collection_STARS["M_felix"      ][i] )
                    assembly_STARS["M2_felix"     ] = vcat( assembly_STARS["M2_felix"     ], merger_collection_STARS["M2_felix"     ][i] )
                    assembly_STARS["M_fromJ"      ] = vcat( assembly_STARS["M_fromJ"      ], merger_collection_STARS["M_fromJ"      ][i] )
                    assembly_STARS["ϕ_flip"       ] = vcat( assembly_STARS["ϕ_flip"       ], merger_collection_STARS["ϕ_flip"       ][i] )
                    assembly_STARS["N_MERGERS"    ] = vcat( assembly_STARS["N_MERGERS"    ], merger_collection_STARS["N_MERGERS"    ][i] )
                    assembly_STARS["M_MERGERS"    ] = vcat( assembly_STARS["M_MERGERS"    ], merger_collection_STARS["M_MERGERS"    ][i] )
                    assembly_STARS["M_MISSED"     ] = vcat( assembly_STARS["M_MISSED"     ], merger_collection_STARS["M_MISSED"     ][i] )
                    assembly_STARS["M_CONSIDERED" ] = vcat( assembly_STARS["M_CONSIDERED" ], merger_collection_STARS["M_CONSIDERED" ][i] )
                    assembly_STARS["M2_MERGERS"   ] = vcat( assembly_STARS["M2_MERGERS"   ], merger_collection_STARS["M2_MERGERS"   ][i] )
                    assembly_STARS["M2_MISSED"    ] = vcat( assembly_STARS["M2_MISSED"    ], merger_collection_STARS["M2_MISSED"    ][i] )
                    assembly_STARS["M2_CONSIDERED"] = vcat( assembly_STARS["M2_CONSIDERED"], merger_collection_STARS["M2_CONSIDERED"][i] )
                    assembly_STARS["BVAL"         ] = vcat( assembly_STARS["BVAL"         ], merger_collection_STARS["BVAL"         ][i] )
                    assembly_STARS["ΔBVAL"        ] = vcat( assembly_STARS["ΔBVAL"        ], merger_collection_STARS["ΔBVAL"        ][i] )
                    assembly_STARS["BVAL_0"       ] = vcat( assembly_STARS["BVAL_0"       ], merger_collection_STARS["BVAL_0"       ][i] )
                    assembly_STARS["ΔBVAL_0"      ] = vcat( assembly_STARS["ΔBVAL_0"      ], merger_collection_STARS["ΔBVAL_0"      ][i] )
                    assembly_STARS["SFR"          ] = vcat( assembly_STARS["SFR"          ], merger_collection_STARS["SFR"          ][i] )
                    assembly_STARS["sSFR"         ] = vcat( assembly_STARS["sSFR"         ], merger_collection_STARS["sSFR"         ][i] )
                    assembly_STARS["BVAL_0"       ] = vcat( assembly_STARS["BVAL_0"       ], merger_collection_STARS["BVAL_0"       ][i] )
                    assembly_STARS["J_MMorbital"  ] = hcat( assembly_STARS["J_MMorbital"  ], merger_collection_STARS["J_MMorbital"  ][:,i] )
                    assembly_STARS["J_SUMorbital" ] = hcat( assembly_STARS["J_SUMorbital" ], merger_collection_STARS["J_SUMorbital" ][:,i] )
                    assembly_STARS["ΔJ_main"      ] = hcat( assembly_STARS["ΔJ_main"      ], merger_collection_STARS["ΔJ_main"      ][:,i] )
                    assembly_STARS["Δj_main"      ] = hcat( assembly_STARS["Δj_main"      ], merger_collection_STARS["Δj_main"      ][:,i] )
                    assembly_STARS["J_main"       ] = hcat( assembly_STARS["J_main"       ], merger_collection_STARS["J_main"       ][:,i] )
                    assembly_STARS["j_main"       ] = hcat( assembly_STARS["j_main"       ], merger_collection_STARS["j_main"       ][:,i] )
                    assembly_STARS["J_vir"        ] = hcat( assembly_STARS["J_vir"        ], merger_collection_STARS["J_vir"        ][:,i] )
                    assembly_STARS["j_vir"        ] = hcat( assembly_STARS["j_vir"        ], merger_collection_STARS["j_vir"        ][:,i] )
                    assembly_STARS["Merger_Map"   ] = hcat( assembly_STARS["Merger_Map"   ], merger_collection_STARS["Merger_Map"   ][:,1+sum(merger_collection_STARS["N_MERGERS"][1:i-1]):sum(merger_collection_STARS["N_MERGERS"][1:i])] )
                    #println("$(storyfilelist[ii])   $i   $(merger_collection_STARS["Merger_Map"][:,i])")
                    assembly_STARS["ΔJ_vir"       ] = hcat( assembly_STARS["ΔJ_vir"      ], merger_collection_STARS["ΔJ_vir"        ][:,i] )
                    assembly_STARS["Δj_vir"       ] = hcat( assembly_STARS["Δj_vir"      ], merger_collection_STARS["Δj_vir"        ][:,i] )
                    assembly_GAS["ΔJ_vir"       ]   = hcat( assembly_GAS["ΔJ_vir"      ], merger_collection_GAS["ΔJ_vir"        ][:,i] )
                    assembly_GAS["Δj_vir"       ]   = hcat( assembly_GAS["Δj_vir"      ], merger_collection_GAS["Δj_vir"        ][:,i] )
                    assembly_GAS["FAKEFLIP"     ]   = vcat( assembly_GAS["FAKEFLIP"   ], fake_flip_finder(merger_collection_GAS, i) )
                    assembly_GAS["SNAP"         ]   = vcat( assembly_GAS["SNAP"         ], merger_collection_GAS["SNAP"         ][i] )
                    assembly_GAS["SWITCH"       ]   = vcat( assembly_GAS["SWITCH"       ], merger_collection_GAS["SWITCH"       ][i] )
                    assembly_GAS["BORDER"       ]   = vcat( assembly_GAS["BORDER"       ], merger_collection_GAS["BORDER"       ][i] )
                    assembly_GAS["JUMP"         ]   = vcat( assembly_GAS["JUMP"         ], merger_collection_GAS["JUMP"         ][i] )
                    assembly_GAS["EXCEED"       ]   = vcat( assembly_GAS["EXCEED"       ], merger_collection_GAS["EXCEED"       ][i] )
                    assembly_GAS["I_SUB"        ]   = vcat( assembly_GAS["I_SUB"        ], merger_collection_GAS["I_SUB"        ][i] )
                    assembly_GAS["ID_ISUB"      ]   = vcat( assembly_GAS["ID_ISUB"      ], merger_collection_GAS["ID_ISUB"      ][i] )
                    assembly_GAS["ID_Mfelix"    ]   = vcat( assembly_GAS["ID_Mfelix"    ], merger_collection_GAS["ID_Mfelix"    ][i] )
                    assembly_GAS["ID_M2"        ]   = vcat( assembly_GAS["ID_M2"        ], merger_collection_GAS["ID_M2"        ][i] )
                    assembly_GAS["REDSHIFT"     ]   = vcat( assembly_GAS["REDSHIFT"     ], merger_collection_GAS["REDSHIFT"     ][i] )
                    assembly_GAS["LOOKBACKTIME" ]   = vcat( assembly_GAS["LOOKBACKTIME" ], merger_collection_GAS["LOOKBACKTIME" ][i] )
                    assembly_GAS["M_MM"         ]   = vcat( assembly_GAS["M_MM"         ], merger_collection_GAS["M_MM"         ][i] )
                    assembly_GAS["M2_MM"        ]   = vcat( assembly_GAS["M2_MM"        ], merger_collection_GAS["M2_MM"        ][i] )
                    assembly_GAS["ΔM_felix"     ]   = vcat( assembly_GAS["ΔM_felix"     ], merger_collection_GAS["ΔM_felix"     ][i] )
                    assembly_GAS["ΔM2_felix"    ]   = vcat( assembly_GAS["ΔM2_felix"    ], merger_collection_GAS["ΔM2_felix"    ][i] )
                    assembly_GAS["ΔM_fromJ"     ]   = vcat( assembly_GAS["ΔM_fromJ"     ], merger_collection_GAS["ΔM_fromJ"     ][i] )
                    assembly_GAS["M_felix"      ]   = vcat( assembly_GAS["M_felix"      ], merger_collection_GAS["M_felix"      ][i] )
                    assembly_GAS["M2_felix"     ]   = vcat( assembly_GAS["M2_felix"     ], merger_collection_GAS["M2_felix"     ][i] )
                    assembly_GAS["M_fromJ"      ]   = vcat( assembly_GAS["M_fromJ"      ], merger_collection_GAS["M_fromJ"      ][i] )
                    assembly_GAS["ϕ_flip"       ]   = vcat( assembly_GAS["ϕ_flip"       ], merger_collection_GAS["ϕ_flip"       ][i] )
                    assembly_GAS["N_MERGERS"    ]   = vcat( assembly_GAS["N_MERGERS"    ], merger_collection_GAS["N_MERGERS"    ][i] )
                    assembly_GAS["M_MERGERS"    ]   = vcat( assembly_GAS["M_MERGERS"    ], merger_collection_GAS["M_MERGERS"    ][i] )
                    assembly_GAS["M_MISSED"     ]   = vcat( assembly_GAS["M_MISSED"     ], merger_collection_GAS["M_MISSED"     ][i] )
                    assembly_GAS["M_CONSIDERED" ]   = vcat( assembly_GAS["M_CONSIDERED" ], merger_collection_GAS["M_CONSIDERED" ][i] )
                    assembly_GAS["M2_MERGERS"   ]   = vcat( assembly_GAS["M2_MERGERS"   ], merger_collection_GAS["M2_MERGERS"   ][i] )
                    assembly_GAS["M2_MISSED"    ]   = vcat( assembly_GAS["M2_MISSED"    ], merger_collection_GAS["M2_MISSED"    ][i] )
                    assembly_GAS["M2_CONSIDERED"]   = vcat( assembly_GAS["M2_CONSIDERED"], merger_collection_GAS["M2_CONSIDERED"][i] )
                    assembly_GAS["J_MMorbital"  ]   = hcat( assembly_GAS["J_MMorbital"  ], merger_collection_GAS["J_MMorbital"  ][:,i] )
                    assembly_GAS["J_SUMorbital" ]   = hcat( assembly_GAS["J_SUMorbital" ], merger_collection_GAS["J_SUMorbital" ][:,i] )
                    assembly_GAS["ΔJ_main"      ]   = hcat( assembly_GAS["ΔJ_main"      ], merger_collection_GAS["ΔJ_main"      ][:,i] )
                    assembly_GAS["J_main"       ]   = hcat( assembly_GAS["J_main"       ], merger_collection_GAS["J_main"       ][:,i] )
                    assembly_GAS["j_main"       ]   = hcat( assembly_GAS["j_main"       ], merger_collection_GAS["j_main"       ][:,i] )
                    assembly_GAS["Δj_main"      ]   = hcat( assembly_GAS["Δj_main"      ], merger_collection_GAS["Δj_main"      ][:,i] )
                    assembly_GAS["J_vir"        ]   = hcat( assembly_GAS["J_vir"        ], merger_collection_GAS["J_vir"        ][:,i] )
                    assembly_GAS["j_vir"        ]   = hcat( assembly_GAS["j_vir"        ], merger_collection_GAS["j_vir"        ][:,i] )
                    assembly_GAS["Merger_Map"   ]   = hcat( assembly_GAS["Merger_Map"   ], merger_collection_GAS["Merger_Map"   ][:,1+sum(merger_collection_GAS["N_MERGERS"][1:i-1]):sum(merger_collection_GAS["N_MERGERS"][1:i])] )
                    assembly_DM["FAKEFLIP"     ]    = vcat( assembly_DM["FAKEFLIP"   ], fake_flip_finder(merger_collection_DM, i) )
                    assembly_DM["SNAP"         ]    = vcat( assembly_DM["SNAP"         ], merger_collection_DM["SNAP"         ][i] )
                    assembly_DM["I_SUB"        ]    = vcat( assembly_DM["I_SUB"        ], merger_collection_DM["I_SUB"        ][i] )
                    assembly_DM["ID_ISUB"      ]    = vcat( assembly_DM["ID_ISUB"      ], merger_collection_DM["ID_ISUB"      ][i] )
                    assembly_DM["ID_Mfelix"    ]    = vcat( assembly_DM["ID_Mfelix"    ], merger_collection_DM["ID_Mfelix"    ][i] )
                    assembly_DM["ID_M2"        ]    = vcat( assembly_DM["ID_M2"        ], merger_collection_DM["ID_M2"        ][i] )
                    assembly_DM["REDSHIFT"     ]    = vcat( assembly_DM["REDSHIFT"     ], merger_collection_DM["REDSHIFT"     ][i] )
                    assembly_DM["LOOKBACKTIME" ]    = vcat( assembly_DM["LOOKBACKTIME" ], merger_collection_DM["LOOKBACKTIME" ][i] )
                    assembly_DM["M_MM"         ]    = vcat( assembly_DM["M_MM"         ], merger_collection_DM["M_MM"         ][i] )
                    assembly_DM["M2_MM"        ]    = vcat( assembly_DM["M2_MM"        ], merger_collection_DM["M2_MM"        ][i] )
                    assembly_DM["ΔM_felix"     ]    = vcat( assembly_DM["ΔM_felix"     ], merger_collection_DM["ΔM_felix"     ][i] )
                    assembly_DM["ΔM2_felix"    ]    = vcat( assembly_DM["ΔM2_felix"    ], merger_collection_DM["ΔM2_felix"    ][i] )
                    assembly_DM["ΔM_fromJ"     ]    = vcat( assembly_DM["ΔM_fromJ"     ], merger_collection_DM["ΔM_fromJ"     ][i] )
                    assembly_DM["M_felix"      ]    = vcat( assembly_DM["M_felix"      ], merger_collection_DM["M_felix"      ][i] )
                    assembly_DM["M2_felix"     ]    = vcat( assembly_DM["M2_felix"     ], merger_collection_DM["M2_felix"     ][i] )
                    assembly_DM["M_fromJ"      ]    = vcat( assembly_DM["M_fromJ"      ], merger_collection_DM["M_fromJ"      ][i] )
                    assembly_DM["ϕ_flip"       ]    = vcat( assembly_DM["ϕ_flip"       ], merger_collection_DM["ϕ_flip"       ][i] )
                    assembly_DM["N_MERGERS"    ]    = vcat( assembly_DM["N_MERGERS"    ], merger_collection_DM["N_MERGERS"    ][i] )
                    assembly_DM["M_MERGERS"    ]    = vcat( assembly_DM["M_MERGERS"    ], merger_collection_DM["M_MERGERS"    ][i] )
                    assembly_DM["M_MISSED"     ]    = vcat( assembly_DM["M_MISSED"     ], merger_collection_DM["M_MISSED"     ][i] )
                    assembly_DM["M_CONSIDERED" ]    = vcat( assembly_DM["M_CONSIDERED" ], merger_collection_DM["M_CONSIDERED" ][i] )
                    assembly_DM["M2_MERGERS"   ]    = vcat( assembly_DM["M2_MERGERS"   ], merger_collection_DM["M2_MERGERS"   ][i] )
                    assembly_DM["M2_MISSED"    ]    = vcat( assembly_DM["M2_MISSED"    ], merger_collection_DM["M2_MISSED"    ][i] )
                    assembly_DM["M2_CONSIDERED"]    = vcat( assembly_DM["M2_CONSIDERED"], merger_collection_DM["M2_CONSIDERED"][i] )
                    assembly_DM["J_MMorbital"  ]    = hcat( assembly_DM["J_MMorbital"  ], merger_collection_DM["J_MMorbital"  ][:,i] )
                    assembly_DM["J_SUMorbital" ]    = hcat( assembly_DM["J_SUMorbital" ], merger_collection_DM["J_SUMorbital" ][:,i] )
                    assembly_DM["ΔJ_main"      ]    = hcat( assembly_DM["ΔJ_main"      ], merger_collection_DM["ΔJ_main"      ][:,i] )
                    assembly_DM["J_main"       ]    = hcat( assembly_DM["J_main"       ], merger_collection_DM["J_main"       ][:,i] )
                    assembly_DM["j_main"       ]    = hcat( assembly_DM["j_main"       ], merger_collection_DM["j_main"       ][:,i] )
                    assembly_DM["Δj_main"      ]    = hcat( assembly_DM["Δj_main"      ], merger_collection_DM["Δj_main"      ][:,i] )
                    assembly_DM["Merger_Map"   ]    = hcat( assembly_DM["Merger_Map"   ], merger_collection_DM["Merger_Map"   ][:,1+sum(merger_collection_DM["N_MERGERS"][1:i-1]):sum(merger_collection_DM["N_MERGERS"][1:i])] )
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
        println("Merger map & N_MERGERS size for stars$(size(assembly_STARS["Merger_Map"]))   =   $(sum(assembly_STARS["N_MERGERS"]))   ?")
        println("Merger map & N_MERGERS size for gas  $(size(assembly_GAS["Merger_Map"]))   =   $(sum(assembly_GAS["N_MERGERS"]))   ?")
        println("Merger map & N_MERGERS size for dm   $(size(assembly_DM["Merger_Map"]))   =   $(sum(assembly_DM["N_MERGERS"]))   ?")
    end
    
    println("\n\n\n---------------------------\n\nNow witness the firepower of this fully armed and operational battle station!")
    return nothing
end

print("'assemble_halostories'   ")
