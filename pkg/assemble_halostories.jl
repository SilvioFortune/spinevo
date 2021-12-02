
@doc """
DESCRIPTION:\n
INPUT:\n
OUTPUT:\n
""" ->
function assemble_halostories(; mass_ST_thr=1e10, central_only=false,
    outdir="./OUT_assemble_halostories", indir="/home/moon/sfortune/spinevo/halostories_v20211127_min0.0Gyr",
    simbox="/HydroSims/Magneticum/Box4/uhr_test", verbose=true
    )

    # Setup
    assembly_STARS = Dict(
        "SNAP"          => missings(Int64  , 0),
        "I_SUB"         => missings(Int64  , 0),
        "N_MERGERS"     => missings(Int64  , 0),
        "ID_ISUB"       => missings(Int64  , 0),
        "ID_Mfelix"     => missings(Float64, 0),
        "ID_M2"         => missings(Float64, 0),
        "REDSHIFT"      => missings(Float64, 0),
        "LOOKBACKTIME"  => missings(Float64, 0),
        "M_MM"          => missings(Float64, 0),
        "M2_MM"         => missings(Float64, 0),
        "δM_felix"      => missings(Float64, 0), 
        "δM2_felix"     => missings(Float64, 0), 
        "δM_fromJ"      => missings(Float64, 0), 
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
        "δBVAL"         => missings(Float64, 0), 
        "J_MMorbital"   => missings(Float64, 3, 0), 
        "J_SUMorbital"  => missings(Float64, 3, 0), 
        "δJ_main"       => missings(Float64, 3, 0), 
        "J_main"        => missings(Float64, 3, 0), 
        "j_main"        => missings(Float64, 3, 0), 
        "δj_main"       => missings(Float64, 3, 0),
        "Merger_Map"    => missings(Float64, 7, 0))
    
    assembly_GAS = Dict(
        "SNAP"          => missings(Int64  , 0),
        "I_SUB"         => missings(Int64  , 0),
        "N_MERGERS"     => missings(Int64  , 0),
        "ID_ISUB"       => missings(Int64  , 0),
        "ID_Mfelix"     => missings(Float64, 0),
        "ID_M2"         => missings(Float64, 0),
        "REDSHIFT"      => missings(Float64, 0),
        "LOOKBACKTIME"  => missings(Float64, 0),
        "M_MM"          => missings(Float64, 0),
        "M2_MM"         => missings(Float64, 0),
        "δM_felix"      => missings(Float64, 0), 
        "δM2_felix"     => missings(Float64, 0), 
        "δM_fromJ"      => missings(Float64, 0), 
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
        "δBVAL"         => missings(Float64, 0), 
        "J_MMorbital"   => missings(Float64, 3, 0), 
        "J_SUMorbital"  => missings(Float64, 3, 0), 
        "δJ_main"       => missings(Float64, 3, 0), 
        "J_main"        => missings(Float64, 3, 0), 
        "j_main"        => missings(Float64, 3, 0), 
        "δj_main"       => missings(Float64, 3, 0),
        "Merger_Map"    => missings(Float64, 7, 0))
    
    assembly_DM = Dict(
        "SNAP"          => missings(Int64  , 0),
        "I_SUB"         => missings(Int64  , 0),
        "N_MERGERS"     => missings(Int64  , 0),
        "ID_ISUB"       => missings(Int64  , 0),
        "ID_Mfelix"     => missings(Float64, 0),
        "ID_M2"         => missings(Float64, 0),
        "REDSHIFT"      => missings(Float64, 0),
        "LOOKBACKTIME"  => missings(Float64, 0),
        "M_MM"          => missings(Float64, 0),
        "M2_MM"         => missings(Float64, 0),
        "δM_felix"      => missings(Float64, 0), 
        "δM2_felix"     => missings(Float64, 0), 
        "δM_fromJ"      => missings(Float64, 0), 
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
        "δBVAL"         => missings(Float64, 0), 
        "J_MMorbital"   => missings(Float64, 3, 0), 
        "J_SUMorbital"  => missings(Float64, 3, 0), 
        "δJ_main"       => missings(Float64, 3, 0), 
        "J_main"        => missings(Float64, 3, 0), 
        "j_main"        => missings(Float64, 3, 0), 
        "δj_main"       => missings(Float64, 3, 0),
        "Merger_Map"    => missings(Float64, 7, 0))
    
    storyfilelist   = readdir(indir)
    noncentralsSUM  = 0
    centralsSUM     = 0
    for ii in 1:length(storyfilelist)
        if verbose
            println("$ii")
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
                if merger_collection_STARS["M2_felix"][i] > mass_ST_thr
                    assembly_STARS["SNAP"         ] = vcat( assembly_STARS["SNAP"         ], merger_collection_STARS["SNAP"         ][i] )
                    assembly_STARS["I_SUB"        ] = vcat( assembly_STARS["I_SUB"        ], merger_collection_STARS["I_SUB"        ][i] )
                    assembly_STARS["ID_ISUB"      ] = vcat( assembly_STARS["ID_ISUB"      ], merger_collection_STARS["ID_ISUB"      ][i] )
                    assembly_STARS["ID_Mfelix"    ] = vcat( assembly_STARS["ID_Mfelix"    ], merger_collection_STARS["ID_Mfelix"    ][i] )
                    assembly_STARS["ID_M2"        ] = vcat( assembly_STARS["ID_M2"        ], merger_collection_STARS["ID_M2"        ][i] )
                    assembly_STARS["REDSHIFT"     ] = vcat( assembly_STARS["REDSHIFT"     ], merger_collection_STARS["REDSHIFT"     ][i] )
                    assembly_STARS["LOOKBACKTIME" ] = vcat( assembly_STARS["LOOKBACKTIME" ], merger_collection_STARS["LOOKBACKTIME" ][i] )
                    assembly_STARS["M_MM"         ] = vcat( assembly_STARS["M_MM"         ], merger_collection_STARS["M_MM"         ][i] )
                    assembly_STARS["M2_MM"        ] = vcat( assembly_STARS["M2_MM"        ], merger_collection_STARS["M2_MM"        ][i] )
                    assembly_STARS["δM_felix"     ] = vcat( assembly_STARS["δM_felix"     ], merger_collection_STARS["δM_felix"     ][i] )
                    assembly_STARS["δM2_felix"    ] = vcat( assembly_STARS["δM2_felix"    ], merger_collection_STARS["δM2_felix"    ][i] )
                    assembly_STARS["δM_fromJ"     ] = vcat( assembly_STARS["δM_fromJ"     ], merger_collection_STARS["δM_fromJ"     ][i] )
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
                    assembly_STARS["δBVAL"        ] = vcat( assembly_STARS["δBVAL"        ], merger_collection_STARS["δBVAL"        ][i] )
                    assembly_STARS["J_MMorbital"  ] = hcat( assembly_STARS["J_MMorbital"  ], merger_collection_STARS["J_MMorbital"  ][:,i] )
                    assembly_STARS["J_SUMorbital" ] = hcat( assembly_STARS["J_SUMorbital" ], merger_collection_STARS["J_SUMorbital" ][:,i] )
                    assembly_STARS["δJ_main"      ] = hcat( assembly_STARS["δJ_main"      ], merger_collection_STARS["δJ_main"      ][:,i] )
                    assembly_STARS["J_main"       ] = hcat( assembly_STARS["J_main"       ], merger_collection_STARS["J_main"       ][:,i] )
                    assembly_STARS["j_main"       ] = hcat( assembly_STARS["j_main"       ], merger_collection_STARS["j_main"       ][:,i] )
                    assembly_STARS["δj_main"      ] = hcat( assembly_STARS["δj_main"      ], merger_collection_STARS["δj_main"      ][:,i] )
                    #println("$(storyfilelist[ii])   $i   $(merger_collection_STARS["Merger_Map"][:,i])")
                    assembly_STARS["Merger_Map"   ] = hcat( assembly_STARS["Merger_Map"   ], merger_collection_STARS["Merger_Map"   ][:,1+sum(merger_collection_STARS["N_MERGERS"][1:i-1]):sum(merger_collection_STARS["N_MERGERS"][1:i])] )
                    assembly_GAS["SNAP"         ]   = vcat( assembly_GAS["SNAP"         ], merger_collection_GAS["SNAP"         ][i] )
                    assembly_GAS["I_SUB"        ]   = vcat( assembly_GAS["I_SUB"        ], merger_collection_GAS["I_SUB"        ][i] )
                    assembly_GAS["ID_ISUB"      ]   = vcat( assembly_GAS["ID_ISUB"      ], merger_collection_GAS["ID_ISUB"      ][i] )
                    assembly_GAS["ID_Mfelix"    ]   = vcat( assembly_GAS["ID_Mfelix"    ], merger_collection_GAS["ID_Mfelix"    ][i] )
                    assembly_GAS["ID_M2"        ]   = vcat( assembly_GAS["ID_M2"        ], merger_collection_GAS["ID_M2"        ][i] )
                    assembly_GAS["REDSHIFT"     ]   = vcat( assembly_GAS["REDSHIFT"     ], merger_collection_GAS["REDSHIFT"     ][i] )
                    assembly_GAS["LOOKBACKTIME" ]   = vcat( assembly_GAS["LOOKBACKTIME" ], merger_collection_GAS["LOOKBACKTIME" ][i] )
                    assembly_GAS["M_MM"         ]   = vcat( assembly_GAS["M_MM"         ], merger_collection_GAS["M_MM"         ][i] )
                    assembly_GAS["M2_MM"        ]   = vcat( assembly_GAS["M2_MM"        ], merger_collection_GAS["M2_MM"        ][i] )
                    assembly_GAS["δM_felix"     ]   = vcat( assembly_GAS["δM_felix"     ], merger_collection_GAS["δM_felix"     ][i] )
                    assembly_GAS["δM2_felix"    ]   = vcat( assembly_GAS["δM2_felix"    ], merger_collection_GAS["δM2_felix"    ][i] )
                    assembly_GAS["δM_fromJ"     ]   = vcat( assembly_GAS["δM_fromJ"     ], merger_collection_GAS["δM_fromJ"     ][i] )
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
                    assembly_GAS["BVAL"         ]   = vcat( assembly_GAS["BVAL"         ], merger_collection_GAS["BVAL"         ][i] )
                    assembly_GAS["δBVAL"        ]   = vcat( assembly_GAS["δBVAL"        ], merger_collection_GAS["δBVAL"        ][i] )
                    assembly_GAS["J_MMorbital"  ]   = hcat( assembly_GAS["J_MMorbital"  ], merger_collection_GAS["J_MMorbital"  ][:,i] )
                    assembly_GAS["J_SUMorbital" ]   = hcat( assembly_GAS["J_SUMorbital" ], merger_collection_GAS["J_SUMorbital" ][:,i] )
                    assembly_GAS["δJ_main"      ]   = hcat( assembly_GAS["δJ_main"      ], merger_collection_GAS["δJ_main"      ][:,i] )
                    assembly_GAS["J_main"       ]   = hcat( assembly_GAS["J_main"       ], merger_collection_GAS["J_main"       ][:,i] )
                    assembly_GAS["j_main"       ]   = hcat( assembly_GAS["j_main"       ], merger_collection_GAS["j_main"       ][:,i] )
                    assembly_GAS["δj_main"      ]   = hcat( assembly_GAS["δj_main"      ], merger_collection_GAS["δj_main"      ][:,i] )
                    assembly_GAS["Merger_Map"   ]   = hcat( assembly_GAS["Merger_Map"   ], merger_collection_GAS["Merger_Map"   ][:,1+sum(merger_collection_GAS["N_MERGERS"][1:i-1]):sum(merger_collection_GAS["N_MERGERS"][1:i])] )
                    assembly_DM["SNAP"         ]    = vcat( assembly_DM["SNAP"         ], merger_collection_DM["SNAP"         ][i] )
                    assembly_DM["I_SUB"        ]    = vcat( assembly_DM["I_SUB"        ], merger_collection_DM["I_SUB"        ][i] )
                    assembly_DM["ID_ISUB"      ]    = vcat( assembly_DM["ID_ISUB"      ], merger_collection_DM["ID_ISUB"      ][i] )
                    assembly_DM["ID_Mfelix"    ]    = vcat( assembly_DM["ID_Mfelix"    ], merger_collection_DM["ID_Mfelix"    ][i] )
                    assembly_DM["ID_M2"        ]    = vcat( assembly_DM["ID_M2"        ], merger_collection_DM["ID_M2"        ][i] )
                    assembly_DM["REDSHIFT"     ]    = vcat( assembly_DM["REDSHIFT"     ], merger_collection_DM["REDSHIFT"     ][i] )
                    assembly_DM["LOOKBACKTIME" ]    = vcat( assembly_DM["LOOKBACKTIME" ], merger_collection_DM["LOOKBACKTIME" ][i] )
                    assembly_DM["M_MM"         ]    = vcat( assembly_DM["M_MM"         ], merger_collection_DM["M_MM"         ][i] )
                    assembly_DM["M2_MM"        ]    = vcat( assembly_DM["M2_MM"        ], merger_collection_DM["M2_MM"        ][i] )
                    assembly_DM["δM_felix"     ]    = vcat( assembly_DM["δM_felix"     ], merger_collection_DM["δM_felix"     ][i] )
                    assembly_DM["δM2_felix"    ]    = vcat( assembly_DM["δM2_felix"    ], merger_collection_DM["δM2_felix"    ][i] )
                    assembly_DM["δM_fromJ"     ]    = vcat( assembly_DM["δM_fromJ"     ], merger_collection_DM["δM_fromJ"     ][i] )
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
                    assembly_DM["BVAL"         ]    = vcat( assembly_DM["BVAL"         ], merger_collection_DM["BVAL"         ][i] )
                    assembly_DM["δBVAL"        ]    = vcat( assembly_DM["δBVAL"        ], merger_collection_DM["δBVAL"        ][i] )
                    assembly_DM["J_MMorbital"  ]    = hcat( assembly_DM["J_MMorbital"  ], merger_collection_DM["J_MMorbital"  ][:,i] )
                    assembly_DM["J_SUMorbital" ]    = hcat( assembly_DM["J_SUMorbital" ], merger_collection_DM["J_SUMorbital" ][:,i] )
                    assembly_DM["δJ_main"      ]    = hcat( assembly_DM["δJ_main"      ], merger_collection_DM["δJ_main"      ][:,i] )
                    assembly_DM["J_main"       ]    = hcat( assembly_DM["J_main"       ], merger_collection_DM["J_main"       ][:,i] )
                    assembly_DM["j_main"       ]    = hcat( assembly_DM["j_main"       ], merger_collection_DM["j_main"       ][:,i] )
                    assembly_DM["δj_main"      ]    = hcat( assembly_DM["δj_main"      ], merger_collection_DM["δj_main"      ][:,i] )
                    assembly_DM["Merger_Map"   ]    = hcat( assembly_DM["Merger_Map"   ], merger_collection_DM["Merger_Map"   ][:,1+sum(merger_collection_DM["N_MERGERS"][1:i-1]):sum(merger_collection_DM["N_MERGERS"][1:i])] )
                end
            end
        end
    end
    
    save(joinpath(outdir, "assembly_Mstar_$(mass_ST_thr).jld"), 
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
