
###################
#  FUNCTIONS from Lucas
###################
function read_nsubhalos(dir, snapnum)
    snapnumpad = lpad(snapnum, 3, "0")
    nsubhalos_cum = [read_subfind_header(joinpath(dir, "groups_$snapnumpad/sub_$snapnumpad.$filenr")).nsubhalos for filenr in 0:15] |> cumsum
    [0; nsubhalos_cum[1:end-1]]
end

function read_all_nsubhalos(dir)
    d = Dict{Int, Vector{Int}}()
    for snapnum in 0:136
        d[snapnum] = read_nsubhalos(dir, snapnum)
    end
    return d
end

get_global_index(d::Dict, snapnum, filenr, subindex) = d[snapnum][filenr + 1] + subindex
get_global_index(d::Dict, halonode) = get_global_index(d, halonode.snapnum, halonode.filenr, halonode.subhalo_index)
##############################################################################################################################




@doc """
DESCRIPTION:\n
INPUT:\n
OUTPUT:\n
""" ->
function initialize_halostory(treeID, trees, subID_map;
    simbox = current_dir_simbox
    )
    return Dict(
        "rootID"        => get_global_index(subID_map, trees.halos[treeID]), 
        "box"           => simbox, 
        "snapNR"        => Vector{Int64}(undef, 0), 
        "snapFP"        => Vector{Int64}(undef, 0), 
        "subID"         => Vector{Int64}(undef, 0), 
        "treeID"        => Vector{Int64}(undef, 0), 
        "fileNR"        => Vector{Int64}(undef, 0), 
        "mmp"           => Vector{Int64}(undef, 0), 
        "switch"        => Vector{Int64}(undef, 0), 
        "border"        => Vector{Int64}(undef, 0), 
        "jump"          => Vector{Int64}(undef, 0), 
        "exceed"        => Vector{Int64}(undef, 0), 
        "redshift"      => Vector{Float64}(undef, 0), 
        "lookbacktime"  => Vector{Float64}(undef, 0), 
        "M_STARS"       => Vector{Float64}(undef, 0), 
        "M_GAS"         => Vector{Float64}(undef, 0), 
        "M_DM"          => Vector{Float64}(undef, 0), 
        "M_STARS_2"     => missings(Float64, 0),
        "M_GAS_2"       => missings(Float64, 0),
        "M_DM_2"        => missings(Float64, 0),
        "M_STARS_peak"  => Vector{Float64}(undef, 0), 
        "M_GAS_peak"    => Vector{Float64}(undef, 0), 
        "M_DM_peak"     => Vector{Float64}(undef, 0), 
        "SFR"           => Vector{Float64}(undef, 0), 
        "sSFR"          => Vector{Float64}(undef, 0), 
        "M200"          => missings(Float64, 0),
        "MCRI"          => missings(Float64, 0),
        "M500"          => missings(Float64, 0),
        "MVIR"          => missings(Float64, 0),
        "SAGE"          => Vector{Float64}(undef, 0), 
        "DSUB"          => Vector{Float64}(undef, 0), 
        "VMAX"          => missings(Float64, 0),
        "RMAX"          => missings(Float64, 0),
        "SPOS"          => Array{Float64}(undef, 3, 0), 
        "SVEL"          => Array{Float64}(undef, 3, 0), 
        "SPIN"          => Array{Float64}(undef, 3, 0), 
        "RVIR"          => missings(Float64, 0),
        "RHMS_DM"       => missings(Float64, 0),
        "J_DM"          => missings(Float64, 3, 0),
        "j_DM"          => missings(Float64, 3, 0),
        "RHMS_GAS"      => missings(Float64, 0),
        "J_GAS"         => missings(Float64, 3, 0),
        "j_GAS"         => missings(Float64, 3, 0),
        "J_GASvir"      => missings(Float64, 3, 0),
        "j_GASvir"      => missings(Float64, 3, 0),
        "RHMS_STARS"    => missings(Float64, 0),
        "J_STARS"       => missings(Float64, 3, 0),
        "j_STARS"       => missings(Float64, 3, 0),
        "J_STARSvir"    => missings(Float64, 3, 0),
        "j_STARSvir"    => missings(Float64, 3, 0),
        "merger_spin_map_STARS"     => missings(Float64, 18, 0),
        "merger_spin_map_GAS"     => missings(Float64, 18, 0),
        "merger_spin_map_DM"     => missings(Float64, 18, 0),
        "BVAL"          => missings(Float64, 0),
        "BVAL_0"        => missings(Float64, 0),
        )
end

print("'initialize_halostory'   ")





#########################################################################




@doc """
DESCRIPTION:\n
INPUT:\n
OUTPUT:\n
""" ->
function append_to_halostory!(halo_story, treeID, trees, subID_map;
    simbox=current_dir_simbox, snap_diff=3, mmp=" ", refID=" ", treeIDcentral=" ",
    snapFP=0,
    )
    if mmp==" "
        error("mmp must be set.")
    end
    #haloID      = get_haloid(trees.halos[treeID])
    #println("Root ID: $haloID")
    snapNR      = trees.halos[treeID].snapnum
    #println(snapNR)
    #subnode     = get_halonode(trees, snapNR, haloID)
    head        = read_header("$simbox/groups_$(lpad(snapNR, 3, "0"))/sub_$(lpad(snapNR, 3, "0"))")
    snapshot    = Snapshot(simbox, snapNR)
    g           = Galaxy(snapshot, get_global_index(subID_map, trees.halos[treeID]))

    # Checkout starting point
    halo_story["snapNR"     ] = vcat(   halo_story["snapNR"     ],   snapNR                 )
    halo_story["snapFP"     ] = vcat(   halo_story["snapFP"     ],   ifelse(mmp==1,snapNR,snapFP)                 )
    halo_story["subID"      ] = vcat(   halo_story["subID"      ],   get_global_index(subID_map, trees.halos[treeID])                 )
    halo_story["treeID"     ] = vcat(   halo_story["treeID"     ],   treeID                 )
    halo_story["fileNR"     ] = vcat(   halo_story["fileNR"     ],   trees.halos[treeID].filenr )
    halo_story["mmp"        ] = vcat(   halo_story["mmp"        ],   mmp                      )
    halo_story["redshift"   ] = vcat(   halo_story["redshift"   ],   head.z                 )
    halo_story["lookbacktime"   ] = vcat(   halo_story["lookbacktime"   ],   ustrip(lookback_time(cosmology(h=head.h0, OmegaM=head.omega_0), head.z))                 )
    halo_story["M_STARS"    ] = vcat(   halo_story["M_STARS"    ],   convert_units_physical_mass(trees.halos[treeID].submasstab[5], head)   )
    halo_story["M_GAS"      ] = vcat(   halo_story["M_GAS"      ],   convert_units_physical_mass(trees.halos[treeID].submasstab[1], head)   )
    halo_story["M_DM"       ] = vcat(   halo_story["M_DM"       ],   convert_units_physical_mass(trees.halos[treeID].submasstab[2], head)   )
    peak_idx, diff_idx = trace_back(trees, treeID, snap_diff=snap_diff)
    halo_story["M_STARS_2"  ] = vcat(   halo_story["M_STARS_2"  ],   convert_units_physical_mass(trees.halos[diff_idx].submasstab[5], 
                                                                        read_header("$simbox/groups_$(lpad(trees.halos[diff_idx].snapnum,3,"0"))/sub_$(lpad(trees.halos[diff_idx].snapnum,3,"0"))"))   )
    halo_story["M_GAS_2"    ] = vcat(   halo_story["M_GAS_2"    ],   convert_units_physical_mass(trees.halos[diff_idx].submasstab[1], 
                                                                        read_header("$simbox/groups_$(lpad(trees.halos[diff_idx].snapnum,3,"0"))/sub_$(lpad(trees.halos[diff_idx].snapnum,3,"0"))"))   )
    halo_story["M_DM_2"     ] = vcat(   halo_story["M_DM_2"     ],   convert_units_physical_mass(trees.halos[diff_idx].submasstab[2], 
                                                                        read_header("$simbox/groups_$(lpad(trees.halos[diff_idx].snapnum,3,"0"))/sub_$(lpad(trees.halos[diff_idx].snapnum,3,"0"))"))   )
    halo_story["M_STARS_peak"  ] = vcat(   halo_story["M_STARS_peak"  ],   convert_units_physical_mass(trees.halos[peak_idx].submasstab[5], 
                                                                            read_header("$simbox/groups_$(lpad(trees.halos[peak_idx].snapnum,3,"0"))/sub_$(lpad(trees.halos[peak_idx].snapnum,3,"0"))"))   )
    halo_story["M_GAS_peak"    ] = vcat(   halo_story["M_GAS_peak"    ],   convert_units_physical_mass(trees.halos[peak_idx].submasstab[1], 
                                                                            read_header("$simbox/groups_$(lpad(trees.halos[peak_idx].snapnum,3,"0"))/sub_$(lpad(trees.halos[peak_idx].snapnum,3,"0"))"))   )
    halo_story["M_DM_peak"     ] = vcat(   halo_story["M_DM_peak"     ],   convert_units_physical_mass(trees.halos[peak_idx].submasstab[2], 
                                                                            read_header("$simbox/groups_$(lpad(trees.halos[peak_idx].snapnum,3,"0"))/sub_$(lpad(trees.halos[peak_idx].snapnum,3,"0"))"))   )
    halo_story["SFR"        ] = vcat(   halo_story["SFR"        ],   read_galaxy_prop(g, "SSFR", :physical) )
    halo_story["M200"       ] = vcat(   halo_story["M200"       ],   read_galaxy_prop(get_group(g), "M200", :physical)  )
    halo_story["MCRI"       ] = vcat(   halo_story["MCRI"       ],   read_galaxy_prop(get_group(g), "MCRI", :physical)  )
    halo_story["M500"       ] = vcat(   halo_story["M500"       ],   read_galaxy_prop(get_group(g), "M500", :physical)  )
    halo_story["MVIR"       ] = vcat(   halo_story["MVIR"       ],   read_galaxy_prop(get_group(g), "MVIR", :physical)  )
    halo_story["VMAX"       ] = vcat(   halo_story["VMAX"       ],   read_galaxy_prop(g, "VMAX", :physical)  )
    halo_story["RMAX"       ] = vcat(   halo_story["RMAX"       ],   read_galaxy_prop(g, "RMAX", :physical)  )
    halo_story["RVIR"       ] = vcat(   halo_story["RVIR"       ],   read_galaxy_prop(get_group(g), "RVIR", :physical)  )
    sph_small   = GadgetGalaxies.Sphere(0.1*halo_story["RVIR"][end])
    sph_large   = GadgetGalaxies.Sphere(halo_story["RVIR"][end])
    halo_story["SAGE"       ] = vcat(   halo_story["SAGE"       ],   read_galaxy_prop(g, "SAGE", :physical) )
    halo_story["DSUB"       ] = vcat(   halo_story["DSUB"       ],   read_galaxy_prop(g, "DSUB", :physical) )
    halo_story["SPOS"       ] = hcat(   halo_story["SPOS"       ],   read_galaxy_prop(g, "SPOS", :physical) )
    halo_story["SVEL"       ] = hcat(   halo_story["SVEL"       ],   read_galaxy_prop(g, "SVEL", :physical) )
    halo_story["SPIN"       ] = hcat(   halo_story["SPIN"       ],   read_galaxy_prop(g, "SPIN", :physical) )
    # Particle Dependent
    if in(snapNR, box4_snaplist)
        try
            @suppress begin
                read_halo!(g, units=:physical, props=((:gas, ["POS", "VEL", "MASS"]), (:stars, ["POS", "VEL", "MASS"]), (:dm, ["POS", "VEL"])), radius_units=:physical, radius=halo_story["RVIR"][end]) 
            end
        catch
            g = " "
        end
        if g != " " && length(g.dm.id) > 2# && GadgetGalaxies.total_mass(g.dm.pos, g.dm.mass) > 0.0
            halo_story["RHMS_DM"    ] = vcat(   halo_story["RHMS_DM"    ],   half_mass_radius(g.dm)                     )
            halo_story["J_DM"       ] = hcat(   halo_story["J_DM"       ],   angular_momentum(g.dm, sph_large)          )
            halo_story["j_DM"       ] = hcat(   halo_story["j_DM"       ],   specific_angular_momentum(g.dm, sph_large) )
        else
            halo_story["RHMS_DM"    ] = vcat(   halo_story["RHMS_DM"    ],   missings(Float64, 1) )
            halo_story["J_DM"       ] = hcat(   halo_story["J_DM"       ],   missings(Float64, 3) )
            halo_story["j_DM"       ] = hcat(   halo_story["j_DM"       ],   missings(Float64, 3) )
        end
        if g != " " && length(g.gas.id) > 2# && GadgetGalaxies.total_mass(g.gas.pos, g.gas.mass) > 0.0
            halo_story["RHMS_GAS"   ] = vcat(   halo_story["RHMS_GAS"   ],   half_mass_radius(g.gas)                        )
            try
                halo_story["J_GAS"      ] = hcat(   halo_story["J_GAS"      ],   angular_momentum(g.gas, sph_small)             )
            catch
                halo_story["J_GAS"      ] = hcat(   halo_story["J_GAS"      ],   missings(Float64, 3)    )
            end
            try
                halo_story["j_GAS"      ] = hcat(   halo_story["j_GAS"      ],   specific_angular_momentum(g.gas, sph_small)    )
            catch
                halo_story["j_GAS"      ] = hcat(   halo_story["j_GAS"      ],   missings(Float64, 3)    )
            end
            halo_story["J_GASvir"   ] = hcat(   halo_story["J_GASvir"   ],   angular_momentum(g.gas, sph_large)             )
            halo_story["j_GASvir"   ] = hcat(   halo_story["j_GASvir"   ],   specific_angular_momentum(g.gas, sph_large)    )
        else
            halo_story["RHMS_GAS"   ] = vcat(   halo_story["RHMS_GAS"   ],   missings(Float64, 1)    )
            halo_story["J_GAS"      ] = hcat(   halo_story["J_GAS"      ],   missings(Float64, 3)    )
            halo_story["j_GAS"      ] = hcat(   halo_story["j_GAS"      ],   missings(Float64, 3)    )
            halo_story["J_GASvir"   ] = hcat(   halo_story["J_GASvir"   ],   missings(Float64, 3)    )
            halo_story["j_GASvir"   ] = hcat(   halo_story["j_GASvir"   ],   missings(Float64, 3)    )
        end
        if g != " " && length(g.stars.id) > 2# && GadgetGalaxies.total_mass(g.stars.pos, g.stars.mass) > 0.0
            try
                halo_story["BVAL"       ] = vcat(   halo_story["BVAL"       ],   b_value(g.stars, sph_small)                      )
            catch
                halo_story["BVAL"       ] = vcat(   halo_story["BVAL"       ],   missings(Float64, 1)    )
            end
            try
                halo_story["J_STARS"    ] = hcat(   halo_story["J_STARS"    ],   angular_momentum(g.stars, sph_small)             )
            catch
                halo_story["J_STARS"    ] = hcat(   halo_story["J_STARS"    ],   missings(Float64, 3)    )
            end
            try
                halo_story["j_STARS"    ] = hcat(   halo_story["j_STARS"    ],   specific_angular_momentum(g.stars, sph_small)    )
            catch
                halo_story["j_STARS"    ] = hcat(   halo_story["j_STARS"    ],   missings(Float64, 3)    )
            end
            halo_story["RHMS_STARS" ] = vcat(   halo_story["RHMS_STARS" ],   half_mass_radius(g.stars)                        )
            halo_story["J_STARSvir" ] = hcat(   halo_story["J_STARSvir" ],   angular_momentum(g.stars, sph_large)             )
            halo_story["j_STARSvir" ] = hcat(   halo_story["j_STARSvir" ],   specific_angular_momentum(g.stars, sph_large)    )
            halo_story["BVAL_0"     ] = vcat(   halo_story["BVAL_0"     ],   halo_story["BVAL"][end] + 0.5*log10(1+head.z)    )
        else
            halo_story["RHMS_STARS" ] = vcat(   halo_story["RHMS_STARS" ],   missings(Float64, 1)    )
            halo_story["BVAL"       ] = vcat(   halo_story["BVAL"       ],   missings(Float64, 1)    )
            halo_story["BVAL_0"     ] = vcat(   halo_story["BVAL_0"     ],   missings(Float64, 1)    )
            halo_story["J_STARS"    ] = hcat(   halo_story["J_STARS"    ],   missings(Float64, 3)    )
            halo_story["j_STARS"    ] = hcat(   halo_story["j_STARS"    ],   missings(Float64, 3)    )
            halo_story["J_STARSvir" ] = hcat(   halo_story["J_STARSvir" ],   missings(Float64, 3)    )
            halo_story["j_STARSvir" ] = hcat(   halo_story["j_STARSvir" ],   missings(Float64, 3)    )
        end
    else
        halo_story["RHMS_DM"    ] = vcat(   halo_story["RHMS_DM"    ],   missings(Float64, 1) )
        halo_story["J_DM"       ] = hcat(   halo_story["J_DM"       ],   missings(Float64, 3) )
        halo_story["j_DM"       ] = hcat(   halo_story["j_DM"       ],   missings(Float64, 3) )
        halo_story["RHMS_GAS"   ] = vcat(   halo_story["RHMS_GAS"   ],   missings(Float64, 1)    )
        halo_story["J_GAS"      ] = hcat(   halo_story["J_GAS"      ],   missings(Float64, 3)    )
        halo_story["j_GAS"      ] = hcat(   halo_story["j_GAS"      ],   missings(Float64, 3)    )
        halo_story["J_GASvir"   ] = hcat(   halo_story["J_GASvir"   ],   missings(Float64, 3)    )
        halo_story["j_GASvir"   ] = hcat(   halo_story["j_GASvir"   ],   missings(Float64, 3)    )
        halo_story["RHMS_STARS" ] = vcat(   halo_story["RHMS_STARS" ],   missings(Float64, 1)    )
        halo_story["BVAL"       ] = vcat(   halo_story["BVAL"       ],   missings(Float64, 1)    )
        halo_story["BVAL_0"     ] = vcat(   halo_story["BVAL_0"     ],   missings(Float64, 1)    )
        halo_story["J_STARS"    ] = hcat(   halo_story["J_STARS"    ],   missings(Float64, 3)    )
        halo_story["j_STARS"    ] = hcat(   halo_story["j_STARS"    ],   missings(Float64, 3)    )
        halo_story["J_STARSvir" ] = hcat(   halo_story["J_STARSvir" ],   missings(Float64, 3)    )
        halo_story["j_STARSvir" ] = hcat(   halo_story["j_STARSvir" ],   missings(Float64, 3)    )
    end

    if mmp == 0
        halo_story["merger_spin_map_STARS"  ] = hcat(   halo_story["merger_spin_map_STARS"  ],   get_merger_spins(trees, subID_map, treeID, treeIDcentral, parttype=:stars, simbox=simbox, verbose=false)    )
        halo_story["merger_spin_map_GAS"  ] = hcat(   halo_story["merger_spin_map_GAS"  ],   get_merger_spins(trees, subID_map, treeID, treeIDcentral, parttype=:gas, simbox=simbox, verbose=false)    )
        halo_story["merger_spin_map_DM"  ] = hcat(   halo_story["merger_spin_map_DM"  ],   get_merger_spins(trees, subID_map, treeID, treeIDcentral, parttype=:dm, simbox=simbox, verbose=false)    )
    else
        halo_story["merger_spin_map_STARS"  ] = hcat(   halo_story["merger_spin_map_STARS"  ],   missings(Float64, 18)    )
        halo_story["merger_spin_map_GAS"  ] = hcat(   halo_story["merger_spin_map_GAS"  ],   missings(Float64, 18)    )
        halo_story["merger_spin_map_DM"  ] = hcat(   halo_story["merger_spin_map_DM"  ],   missings(Float64, 18)    )
    end

    return halo_story
end

print("'append_to_halostory!'   ")





#########################################################################






@doc """
DESCRIPTION:\n
INPUT:\n
OUTPUT:\n
""" ->
function trace_back(trees, treeID; 
    part_type = "stars", snap_diff = 3,
    )
    # Setup
    peak_idx    = treeID
    peak_mass   = trees.halos[peak_idx].submasstab[convert_parttype_to_idx(str=part_type)]
    diff_idx    = treeID
    #haloID      = get_haloid(trees.halos[treeID])
    #snapNR      = trees.halos[treeID].snapnum
    #subnode     = get_halonode(trees, snapNR, haloID)
    #fp_idx      = subnode.first_progenitor
    fp_idx      = trees.halos[treeID].first_progenitor
    while fp_idx > -1
        #print("$fp_idx ")
        #println(trees.halos[fp_idx].submasstab[convert_parttype_to_idx(str=part_type)])
        # Mass Check
        if trees.halos[fp_idx].submasstab[convert_parttype_to_idx(str=part_type)] > trees.halos[peak_idx].submasstab[convert_parttype_to_idx(str=part_type)]
            peak_idx    = fp_idx
            peak_mass   = trees.halos[peak_idx].submasstab[convert_parttype_to_idx(str=part_type)]
        end
        # snap diff check
        if trees.halos[fp_idx].snapnum == trees.halos[treeID].snapnum - snap_diff
            #println(trees.halos[fp_idx].snapnum)
            diff_idx    = fp_idx
        end
        # move to to next Snap
        #haloID      = get_haloid(trees.halos[fp_idx])
        #snapNR      = trees.halos[fp_idx].snapnum
        #subnode     = get_halonode(trees, snapNR, haloID)
        #fp_idx      = subnode.first_progenitor
        fp_idx      = trees.halos[fp_idx].first_progenitor
    end

    # return
    return peak_idx, diff_idx
end

print("'trace_back'   ")



#######################################################################################################
######################################################################################################
#######################################################################################################



@doc """
DESCRIPTION:\n
INPUT:\n
OUTPUT:\n
""" ->
function select_root_halos(; start=" ", stop=" ", root_snap=136,
    root_mass_thr=0, part_type="stars", id_list=" ",
    outdir=current_dir_data, simbox=current_dir_simbox,
    iterstep=1, 
    verbose=true, output_type="file", filename_appendix="", 
    )
    # Setup
    if verbose
        println("Initiating select_root_halos() with:\n 
                    start           = $(start)\n 
                    stop            = $(stop)\n 
                    root_mass_thr   = $(root_mass_thr)\n 
                    outdir          = $(outdir)\n 
                    simbox          = $(simbox)\n 
                    part_type       = $(part_type)\n 
                    output_type     = $(output_type)\n 
                    iterstep        = $(iterstep)\n
                    verbose         = $(verbose)\n
                    filename_appendix    = $(filename_appendix)\n
                    ")
        flush(stdout)
    end
    trees       = read_tree(simbox, 136)
    subID_map   = read_all_nsubhalos(simbox)
    if start == " "
        start   = 1
    end
    if stop == " "
        stop    = length(trees.halos)
    end

    root_list = Dict(
        "simbox"        => simbox, 
        "snapNR"        => root_snap, 
        "part_type"     => part_type, 
        "root_mass_thr" => root_mass_thr, 
        "all_idx"       => Vector{Int64}(undef, 0), 
        "rootID"        => Vector{Int64}(undef, 0), 
        "central_idx"   => Vector{Int64}(undef, 0), 
        "smst"          => Array{Float64}(undef, 6, 0), 
        )
    root_head   = read_header("$simbox/groups_$(lpad(root_snap, 3, "0"))/sub_$(lpad(root_snap, 3, "0"))")
    snapshot    = Snapshot(simbox, root_snap)

    if verbose
        println("   Initiating loop from $start to $stop. ")
        flush(stdout)
    end
    if typeof(id_list) == String
        for i in start:stop
            #rootID  = get_haloid(trees.halos[i]).id-1
            rootID  = get_global_index(subID_map, trees.halos[i])
            smst    = convert_units_physical_mass(trees.halos[i].submasstab, root_head)
            if trees.halos[i].snapnum == root_snap && smst[convert_parttype_to_idx(str=part_type)] > root_mass_thr
                root_list["all_idx"     ] = vcat( root_list["all_idx"       ], i )
                root_list["rootID"      ] = vcat( root_list["rootID"        ], rootID )
                root_list["smst"        ] = hcat( root_list["smst"          ], smst )
                root_list["central_idx" ] = vcat( root_list["central_idx"   ], ifelse(rootID == get_first_subhalo(get_group(Galaxy(snapshot, rootID))).isub, 1, 0) )
            end
        end
    else
        for i in start:stop
            #rootID  = get_haloid(trees.halos[i]).id-1
            rootID  = get_global_index(subID_map, trees.halos[i])
            smst    = convert_units_physical_mass(trees.halos[i].submasstab, root_head)
            if trees.halos[i].snapnum == root_snap && smst[convert_parttype_to_idx(str=part_type)] > root_mass_thr && rootID in id_list
                root_list["all_idx"     ] = vcat( root_list["all_idx"       ], i )
                root_list["rootID"      ] = vcat( root_list["rootID"        ], rootID )
                root_list["smst"        ] = hcat( root_list["smst"          ], smst )
                root_list["central_idx" ] = vcat( root_list["central_idx"   ], ifelse(rootID == get_first_subhalo(get_group(Galaxy(snapshot, rootID))).isub, 1, 0) )
            end
        end
    end

    
    if verbose
        println("\n\n\n---------------------------\n\nNow witness the firepower of this fully armed and operational battle station!\n---------------------------\n\n")
        flush(stdout)
    end
    
    if output_type == "file"
        save(
            joinpath(outdir, "root_list_"*replace("$(today())", "-" => "" )*filename_appendix*".jld"), 
            "root_list",   root_list
            )
        return nothing
    elseif output_type == "return"
        return root_list
    elseif output_type == "both"
        save(
            joinpath(outdir, "root_list_"*replace("$(today())", "-" => "" )*filename_appendix*".jld"), 
            "root_list",   root_list
            )
        return root_list
    else
        error("Unknown output type: $output_type")
    end
end

print("'select_root_halos'   ")



#######################################################################################################
######################################################################################################
#######################################################################################################





#########################################################################




include("/home/moon/sfortune/spinevo/pkg/treewalker/collect_mergers.jl")



###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################






@doc """
DESCRIPTION:\n
INPUT:\n
OUTPUT:\n
""" ->
function write_halo_stories(; subtype = "all", root_list = " ", 
    start=" ", stop=" ", iterstep=1,
    outdir=current_dir_stories,
    root_mass_thr=0, part_type="stars", merger_mass_thr=2e8, min_time=0.0,
    simbox=current_dir_simbox, 
    snap_diff = 3,
    verbose=true, skipexisting=false, central_switch=false,
    # root list settings
    root_snap=136, 
    )
    # Setup

    ### Global Info
    trees       = read_tree(simbox, root_snap)
    subID_map    = read_all_nsubhalos(simbox)
    if root_list == " "
        root_list = select_root_halos(output_type="return",root_snap=root_snap, root_mass_thr=root_mass_thr, part_type=part_type, simbox=simbox, verbose=false)
    elseif typeof(root_list) == Dict{String, Any}
    elseif typeof(root_list) == String
        root_list = load(root_list, "root_list")
    else
        error("Unknown root_list: $root_list")
        flush(stdout)
    end

    ### Selection
    if occursin("central", lowercase(subtype))
        selection = findcs( root_list["central_idx"], eq=1, comparewith=findcs(root_list["smst"][convert_parttype_to_idx(str=part_type), :], geq=root_mass_thr) )
    elseif occursin("all", lowercase(subtype))
        selection = findcs(root_list["smst"][convert_parttype_to_idx(str=part_type), :], geq=root_mass_thr)
    elseif occursin("non-central", lowercase(subtype))
        selection = findcs(root_list["central_idx"], eq=0, comparewith=findcs(root_list["smst"][convert_parttype_to_idx(str=part_type), :], geq=root_mass_thr) )
    end
    if start == " "
        start   = 1
    end
    if stop == " "
        stop    = length(selection)
    end

    if verbose
        println("Initiating select_root_halos() with:\n 
                    start           = $(start)\n 
                    stop            = $(stop)\n 
                    root_snap       = $(root_snap)\n 
                    snap_diff       = $(snap_diff)\n 
                    root_mass_thr   = $(root_mass_thr)\n 
                    merger_mass_thr = $(merger_mass_thr)\n 
                    outdir          = $(outdir)\n 
                    simbox          = $(simbox)\n 
                    part_type       = $(part_type)\n 
                    subtype         = $(subtype)\n
                    verbose         = $(verbose)\n
                    skipexisting    = $(skipexisting)\n
                    central_switch  = $(central_switch)\n
                    ")
                    #root_list       = $(root_list)\n 
        flush(stdout)
    end

    if verbose
        println("Initiating main loop for $(iterstep*(stop-start)+1) halos from entries $start to $stop.")
        flush(stdout)
    end
    # Loop over root halos

    hcount = 0
    for i in root_list["all_idx"][selection][start:iterstep:stop]
        if verbose
            println("$(start+hcount*iterstep)   ---   $(get_global_index(subID_map, trees.halos[i]))")
            flush(stdout)
        end
        hcount += 1
        switch = 0
        fp_idx = i
        halo_story = initialize_halostory(fp_idx, trees, subID_map, simbox=simbox)
        append_to_halostory!(halo_story, fp_idx, trees, subID_map, mmp=1, simbox=simbox, snap_diff=snap_diff)
        #fp_idx      = subnode.first_progenitor
        if central_switch && trees.halos[fp_idx].first_progenitor != -1
            if trees.halos[fp_idx].first_progenitor == trees.halos[trees.halos[fp_idx].first_progenitor].first_halo_in_fof_group
                fp_idx  = trees.halos[fp_idx].first_progenitor
                switch  = 0
                halo_story["switch" ] = vcat(   halo_story["switch" ],   switch  )
            else
                fp_idx  = trees.halos[trees.halos[fp_idx].first_progenitor].first_halo_in_fof_group
                switch  = 1
                halo_story["switch" ] = vcat(   halo_story["switch" ],   switch  )
            end
        else
            fp_idx      = trees.halos[fp_idx].first_progenitor
            switch  = 0
            halo_story["switch" ] = vcat(   halo_story["switch" ],   switch  )
        end

        snapcount = 0
        while fp_idx > -1
            snapcount += 1
            if verbose
                print("$(trees.halos[fp_idx].snapnum)", ifelse(snapcount % 30 == 0, "\n", " "))
                flush(stdout)
            end
            # Check out halo info
            #print("FP ")
            append_to_halostory!(halo_story, fp_idx, trees, subID_map, mmp=1, simbox=simbox, snap_diff=snap_diff)
            
            # Snap update
            #fp_snap     = trees.halos[fp_idx].snapnum
            snapNR      = trees.halos[fp_idx].snapnum
            #println("$snapNR")

            # move to to next Snap
            #haloID      = get_haloid(trees.halos[fp_idx])
            #subnode     = get_halonode(trees, fp_snap, haloID)
            #fp_idx      = subnode.first_progenitor
            pre_fp_idx = fp_idx
            if central_switch && trees.halos[fp_idx].first_progenitor != -1
                if trees.halos[fp_idx].first_progenitor == trees.halos[trees.halos[fp_idx].first_progenitor].first_halo_in_fof_group
                    fp_idx  = trees.halos[fp_idx].first_progenitor
                    switch  = 0
                    halo_story["switch" ] = vcat(   halo_story["switch" ],   switch  )
                else
                    fp_idx  = trees.halos[trees.halos[fp_idx].first_progenitor].first_halo_in_fof_group
                    switch  = 1
                    halo_story["switch" ] = vcat(   halo_story["switch" ],   switch  )
                end
            elseif !central_switch && trees.halos[fp_idx].first_progenitor != -1
                if trees.halos[fp_idx].first_progenitor == trees.halos[trees.halos[fp_idx].first_progenitor].first_halo_in_fof_group
                    switch  = 0
                else
                    switch  = 2
                end
                fp_idx      = trees.halos[fp_idx].first_progenitor
                halo_story["switch" ] = vcat(   halo_story["switch" ],   switch  )
            else
                fp_idx      = trees.halos[fp_idx].first_progenitor
                switch  = 0
                halo_story["switch" ] = vcat(   halo_story["switch" ],   switch  )
            end
        
            #println("FP = $(trees.halos[fp_idx].subhalo_index) @ snap $(fp_snap)   ---   $(get_global_index(subID_map, trees.halos[fp_idx]))")
            #print(" NP: ")
            np_idx  = trees.halos[pre_fp_idx].next_progenitor
            while np_idx > -1
                #np_node = get_halonode(trees, trees.halos[np_idx].snapnum, get_haloid(trees.halos[np_idx]))
                if convert_units_physical_mass(trees.halos[np_idx].submasstab[5], read_header("$simbox/groups_$(lpad(snapNR, 3, "0"))/sub_$(lpad(snapNR, 3, "0"))")) ≥ merger_mass_thr
                    #print("$(trees.halos[np_idx].subhalo_index) / $(get_global_index(subID_map, trees.halos[np_idx]))   ---   ")
                    #print("   NP readout: $(trees.halos[pre_fp_idx].snapnum)   ---   $(trees.halos[np_idx].snapnum)      ------      ")
                    append_to_halostory!(halo_story, np_idx, trees, subID_map, mmp=0, snapFP=trees.halos[pre_fp_idx].snapnum, treeIDcentral=pre_fp_idx, simbox=simbox)
                    halo_story["switch" ] = vcat(   halo_story["switch" ],   switch  )
                end
                np_idx  = trees.halos[np_idx].next_progenitor
            end
            #println()
        
    
        end
        if verbose
            println()
            flush(stdout)
        end


        ### Transitional
        halo_story["switch" ] = vcat(   halo_story["switch" ],   0  )
        #halo_story["border" ] = vcat(   halo_story["border" ],   0  )
        #halo_story["jump"   ] = vcat(   halo_story["jump"   ],   0  )
        #halo_story["exceed" ] = vcat(   halo_story["exceed" ],   0  )
        

        merger_collection_STARS, merger_collection_GAS, merger_collection_DM = collect_mergers(halo_story, min_time=min_time)
        save(joinpath(outdir, "halo_$(halo_story["rootID"]).jld"), 
            "halo_story",   halo_story,
            "merger_collection_DM",     merger_collection_DM,
            "merger_collection_GAS",    merger_collection_GAS,
            "merger_collection_STARS",  merger_collection_STARS)
    end
    

    
    println("\n\n\n---------------------------\n\nNow witness the firepower of this fully armed and operational battle station!")
end

print("'write_halo_stories'   ")



##############################################################################################################################




@doc """
DESCRIPTION:\n
INPUT:\n
OUTPUT:\n
""" ->
function update_mergercollections(; verbose=true,
    story_dir=current_dir_stories, min_time=0.0,
    start=" ", stop=" ",
    )

    storyfilelist    = readdir(story_dir)
    if start==" "
        start = 1
    end
    if stop==" "
        stop = length(storyfilelist)
    end
    if verbose
        println("Updating $(length(storyfilelist)) files.")
        flush(stdout)
    end
    for i in start:stop
        if verbose
            print("$(i)-$(chop(storyfilelist[i],head=5,tail=4))", ifelse(i % 10 == 0, "\n", "   "))
            flush(stdout)
        end
        halo_story = load(joinpath(story_dir, storyfilelist[i]), "halo_story")
        merger_collection_STARS, merger_collection_GAS, merger_collection_DM = collect_mergers(halo_story, min_time=min_time)
        save(joinpath(story_dir, "halo_$(halo_story["rootID"]).jld"), 
            "halo_story",   halo_story,
            "merger_collection_DM",     merger_collection_DM,
            "merger_collection_GAS",    merger_collection_GAS,
            "merger_collection_STARS",  merger_collection_STARS)

    end

    return nothing
end

print("'update_mergercollections'   ")



##############################################################################################################################




@doc """
DESCRIPTION:\n
INPUT:\n
OUTPUT:\n
""" ->
function repair_halostories(; verbose=true,
    story_dir=current_dir_stories,
    subtype="central",central_switch=true, root_mass_thr=1e10, root_list=current_root_list, # @write_halo_stories()
    )
    storyfilelist    = readdir(story_dir)
    if verbose
        println("Scanning and repairing $(length(storyfilelist)) files.")
        flush(stdout)
    end
    for i in 1:length(storyfilelist)
        try
            halo_story = load(joinpath(story_dir, storyfilelist[i]), "halo_story")
            if verbose
                print("$(i)-$(chop(storyfilelist[i],head=5,tail=4))", ifelse(i % 10 == 0, "\n", "   "))
                flush(stdout)
            end
        catch
            if verbose
                print("\n$(i)-$(chop(storyfilelist[i],head=5,tail=4)) needs repairs...\n")
                flush(stdout)
            end
            write_halo_stories(start=i, stop=i, 
                central_switch=central_switch, root_mass_thr=root_mass_thr, subtype = subtype, 
                root_list=root_list, 
                outdir=story_dir)
        end
    end
    
    println("\n\n\n---------------------------\n\nNow witness the firepower of this fully armed and operational battle station!")
    return nothing
end

print("'repair_halostories'   ")



##############################################################################################################################




@doc """
DESCRIPTION:\n
INPUT:\n
OUTPUT:\n
""" ->
function get_merger_spins(trees, subID_map, merger_itree, central_itree; 
    simbox=current_dir_simbox, parttype=:stars,
    verbose=true,
    min_t_diff=0.9, # Gyr
    )

     # 1 = snap immediate, 2 = lbt immediate, 3 = mass immediate, 4-6 = orbital immediate, 7-9 intrinsic immediate, 
     # 10 = snap later, 11 = lbt later, 12 = mass later, 13-15 orbital later, 16-18, intrinsic later
    j_coll = missings(Float64, 18)
    if trees.halos[central_itree].snapnum ≥ minimum(box4_snaplist) # check for available snapshots
        # aux from input
        in_head = read_header("$simbox/groups_$(lpad(trees.halos[central_itree].snapnum, 3, "0"))/sub_$(lpad(trees.halos[central_itree].snapnum, 3, "0"))")
        in_lbt  = ustrip(lookback_time(cosmology(h=in_head.h0, OmegaM=in_head.omega_0), in_head.z))
        
        # immediate values
        snap_imm    = maximum(box4_snaplist[findcs(box4_snaplist, leq=trees.halos[central_itree].snapnum)]) # defined by FP
        try
            m_imm_tID   = trace_back(trees, merger_itree, snap_diff=trees.halos[central_itree].snapnum-snap_imm)[2]
            c_imm_tID   = trace_back(trees, central_itree, snap_diff=trees.halos[central_itree].snapnum-trees.halos[m_imm_tID].snapnum)[2] # match merger snap
            j_coll[1]   = trees.halos[m_imm_tID].snapnum
            j_coll[2]   = box4_fulllbt[findcs(box4_fullsnaplist, eq = Int64(j_coll[1]))][end]
            
            head        = read_header("$simbox/groups_$(lpad(trees.halos[m_imm_tID].snapnum, 3, "0"))/sub_$(lpad(trees.halos[m_imm_tID].snapnum, 3, "0"))")
            j_coll[3]   = convert_units_physical_mass(trees.halos[m_imm_tID].submasstab[convert_parttype_to_idx(sym=parttype)], head)

            g_m         = Galaxy(Snapshot(simbox, trees.halos[m_imm_tID].snapnum), get_global_index(subID_map, trees.halos[m_imm_tID]))
            g_c         = Galaxy(Snapshot(simbox, trees.halos[c_imm_tID].snapnum), get_global_index(subID_map, trees.halos[c_imm_tID]))
            j_coll[4:6] = relative_j(g_m, g_c)
            try
                @suppress begin
                    read_halo!(g_m, units=:physical, props=((parttype, ["POS", "VEL", "MASS"]),)) 
                end
            catch
                g_m = " "
            end
            if g_m != " " && length(g_m[parttype].id) > 2
                j_coll[7:9] = angular_momentum(g_m[parttype])
            end
        catch 
            if verbose
                println("Cannot trace back halos to snap $(j_coll[1])")
            end
        end
        

        if maximum(box4_lbt) - in_lbt ≥ min_t_diff
            lbt_min     = minimum(box4_lbt[findcs(box4_lbt .- in_lbt, geq = min_t_diff)]) # defined by FP (-> in_lbt)
            snap_min    = box4_snaplist[findcs(box4_lbt, eq = lbt_min)][end]
            try
                m_imm_tID   = trace_back(trees, merger_itree, snap_diff=trees.halos[central_itree].snapnum-snap_min)[2]
                c_imm_tID   = trace_back(trees, central_itree, snap_diff=trees.halos[central_itree].snapnum-trees.halos[m_imm_tID].snapnum)[2] # match merger snap
                j_coll[10]  = trees.halos[m_imm_tID].snapnum
                j_coll[11]  = box4_fulllbt[findcs(box4_fullsnaplist, eq = Int64(j_coll[10]))][end]
                head        = read_header("$simbox/groups_$(lpad(trees.halos[m_imm_tID].snapnum, 3, "0"))/sub_$(lpad(trees.halos[m_imm_tID].snapnum, 3, "0"))")
                j_coll[12]  = convert_units_physical_mass(trees.halos[m_imm_tID].submasstab[convert_parttype_to_idx(sym=parttype)], head)
    
                g_m             = Galaxy(Snapshot(simbox, trees.halos[m_imm_tID].snapnum), get_global_index(subID_map, trees.halos[m_imm_tID]))
                g_c             = Galaxy(Snapshot(simbox, trees.halos[c_imm_tID].snapnum), get_global_index(subID_map, trees.halos[c_imm_tID]))
                j_coll[13:15]   = relative_j(g_m, g_c)
                try
                    @suppress begin
                        read_halo!(g_m, units=:physical, props=((parttype, ["POS", "VEL", "MASS"]),)) 
                    end
                catch
                    g_m = " "
                end
                if g_m != " " && length(g_m[parttype].id) > 2
                    j_coll[16:18]   = angular_momentum(g_m[parttype])
                end
            catch 
                if verbose
                    println("Cannot trace back halos to snap $(j_coll[10])")
                end
            end
        end


        #try
        #    for i_snap in sort(box4_snaplist[findcs(box4_snaplist, leq=j_coll[1])], rev=true) # from j_coll[1] to 12
        #        late_head = read_header("$simbox/groups_$(lpad(trees.halos[merger_itree].snapnum, 3, "0"))/sub_$(lpad(trees.halos[merger_itree].snapnum, 3, "0"))")
        #        late_lbt  = ustrip(lookback_time(cosmology(h=late_head.h0, OmegaM=late_head.omega_0), late_head.z))
        #        if late_lbt - in_lbt ≥ min_t_diff
        #            late_snap   = i_snap
        #            j_coll[10]  = late_snap
        #            j_coll[11]  = late_lbt
        #            break
        #        end
        #    end
        #catch
        #end
    
        return j_coll
    end




    return j_coll # done: 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18
end

print("'get_merger_spins'   ")