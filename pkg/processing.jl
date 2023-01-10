println("   Processing...")


@doc """
DESCRIPTION:\n
    - reduce particle number by distribution to a 3D-grid
INPUT:\n
    - n_cells:          approximate number of cells for the final grid
    - positions:        (3 x n)D array for the positions of n reduce_n_particles
    - masses:           nD array for the particle masses
    - focus_type:       "center", "outerrim" or "linear"; determines where the cell resolution should be highest
    - focus_strength:   [0. , 1.0] determines degree of focus
    - m_th:             minimum cell mass fraction of total mass to accept the cell in the final output, used to get rid of empty cells
    - verbose:          true or false; progress output
OUTPUT:\n
    - pos_new:  (3 x m)D array of mass-weighted cell positions
    - mnew:     mD array of cell masses
    - n_part:   mD array of number of particles in each cell
""" ->
function reduce_n_particles(n_cells, positions, masses; focus_type="center", focus_strength=0.7, m_th=1e-6, verbose=true)
    
    pos_new = Array{Float64}(undef, 3, 0)
    mnew    = Array{Float64}(undef, 0)
    n_part      = Array{Int64}(undef, 0)
    mtot    = sum(masses)


    x_ext   = maximum(positions[1,:])-minimum(positions[1,:])
    y_ext   = maximum(positions[2,:])-minimum(positions[2,:])
    z_ext   = maximum(positions[3,:])-minimum(positions[3,:])
    cell_length = maximum([z_ext,y_ext,x_ext]) / cbrt(n_cells)
    if focus_type == "center"
        limits  = 0.7 + focus_strength*0.8    # x-range for tan function: 0.7 to 1.5
        xborders    = borders_tan(positions[1,:], cell_length, limits)
        yborders    = borders_tan(positions[2,:], cell_length, limits)
        zborders    = borders_tan(positions[3,:], cell_length, limits)
    elseif focus_type == "outerrim"
        limits  = 0.5 + focus_strength*4.5    # x-range for arctan function: 0.5 to 5
        xborders    = borders_arctan(positions[1,:], cell_length, limits)
        yborders    = borders_arctan(positions[2,:], cell_length, limits)
        zborders    = borders_arctan(positions[3,:], cell_length, limits)
    elseif focus_type == "linear"
        xborders    = LinRange(minimum(positions[1,:]), maximum(positions[1,:]), Int(round(x_ext/cell_length))+1)
        yborders    = LinRange(minimum(positions[2,:]), maximum(positions[2,:]), Int(round(y_ext/cell_length))+1)
        zborders    = LinRange(minimum(positions[3,:]), maximum(positions[3,:]), Int(round(z_ext/cell_length))+1)
    else
        error("Unknown Focus type")
    end
    #cell_IDs    = Array{Int64}(undef, 3, length(masses)) # assign each particle to a cell

    # looping over upper cell limits to assign to cells
    if verbose == true
        println("$(length(xborders))   ---   $(length(yborders))   ---   $(length(zborders))")
    end
    #Threads.@threads for i in 2:length(xborders)
    for i in 2:length(xborders)
        if verbose == true
            println("$i ")
        end
        for ii in 2:length(yborders)
            for iii in 2:length(zborders)
                m_temp      = Array{Float64}(undef, 0)
                pos_temp    = Array{Float64}(undef, 3, 0)
                n_temp      = 0
                for iiii in 1:length(masses)
                    if xborders[i-1] <= positions[1,iiii] <= xborders[i] && yborders[ii-1] <= positions[2,iiii] <= yborders[ii] && zborders[iii-1] <= positions[3,iiii] <= zborders[iii]
                        #cell_IDs[:,iiii] = [i,ii,iii]
                        m_temp           = vcat(m_temp, masses[iiii])
                        pos_temp         = hcat(pos_temp, positions[:,iiii])
                        n_temp          += 1
                    end
                end
                if sum(m_temp) > mtot*m_th
                    n_part  = vcat(n_part, n_temp)
                    mnew    = vcat(mnew, sum(m_temp))
                    #pos_new = hcat(pos_new, [ mean(pos_temp[1,:], weights(m_temp)), mean(pos_temp[2,:], weights(m_temp)), mean(pos_temp[3,:], weights(m_temp)) ])
                    pos_new = hcat(pos_new, mean(pos_temp, weights(m_temp); dims=2))
                end
            end
        end
    end
    println("$(size(pos_new))   $(size(mnew))   $(size(n_part))")
    return pos_new, mnew, n_part
end

print("'reduce_n_particles'   ")



########################################################################################################################################
########################################################################################################################################
########################################################################################################################################
########################################################################################################################################
########################################################################################################################################
########################################################################################################################################
########################################################################################################################################
########################################################################################################################################
########################################################################################################################################



@doc """
DESCRIPTION:\n
INPUT:\n
OUTPUT:\n
""" ->
function fake_flip_finder(indict, idx; verbose=false)
    flip1 = 3
    if 1 < idx < length(indict["M"]) && 
    count(ismissing, indict["j_main"][:,idx-1]) + count(ismissing, indict["j_main"][:,idx]) + 
    count(ismissing, indict["j_main"][:,idx+1]) == 0 # can be checked
        flip1 = ifelse( norm(indict["j_main"][:,idx+1]-indict["j_main"][:,idx-1]) < norm(indict["j_main"][:,idx]-indict["j_main"][:,idx-1]) && # && indict["M"][idx+1]-indict["M"][idx-1] < indict["M"][idx]-indict["M"][idx-1], 
                        norm(indict["j_main"][:,idx+1]-indict["j_main"][:,idx-1]) < norm(indict["j_main"][:,idx+1]-indict["j_main"][:,idx]),
                        1, 0 
                        )
    end

    flip2 = 3
    if 2 < idx && 
    count(ismissing, indict["j_main"][:,idx-2]) + count(ismissing, indict["j_main"][:,idx-1]) + 
    count(ismissing, indict["j_main"][:,idx]) == 0 # can be checked
        flip2 = ifelse( norm(indict["j_main"][:,idx]-indict["j_main"][:,idx-2]) < norm(indict["j_main"][:,idx-1]-indict["j_main"][:,idx-2]) && # && indict["M"][idx]-indict["M"][idx-2] < indict["M"][idx-1]-indict["M"][idx-2], 
                        norm(indict["j_main"][:,idx]-indict["j_main"][:,idx-2]) < norm(indict["j_main"][:,idx]-indict["j_main"][:,idx-1]),
                        1, 0 
                        )
    end
    if verbose
        println("Snap $(indict["snapNR"][idx])   ---   in = $flip1, out = $flip2")
    end

    if flip1 + flip2 == 0 # both not fake
        return 0
    elseif flip1 + flip2 < 3 # any is fake
        return 1
    else
        return missing
    end
end

print("'fake_flip_finder'   ")



########################################################################################################################################
########################################################################################################################################
########################################################################################################################################
########################################################################################################################################
########################################################################################################################################
########################################################################################################################################
########################################################################################################################################
########################################################################################################################################
########################################################################################################################################



# Functions that are too long and have own file
include("/home/moon/sfortune/spinevo/pkg/tree2story.jl")
include("/home/moon/sfortune/spinevo/pkg/assemble_halostories.jl")
include("/home/moon/sfortune/spinevo/pkg/treewalker.jl")



########################################################################################################################################
########################################################################################################################################
########################################################################################################################################
########################################################################################################################################
########################################################################################################################################
########################################################################################################################################
########################################################################################################################################
########################################################################################################################################
########################################################################################################################################



@doc """
DESCRIPTION:\n
INPUT:\n
OUTPUT:\n
""" ->
function quickadd_feature(type; target=" ", 
                                simbox=current_dir_simbox,
                                verbose=true,
                                min_time=0.0,
                                )
    if type == "hs"
        storyfilelist   = readdir(target)
        count = 0
        for i in storyfilelist
            count += 1
            halo_story = load(joinpath(target, i), "halo_story")
            if verbose
                print("$count-$(halo_story["rootID"])", ifelse(count % 30 == 0, "\n", " "))
                flush(stdout)
            end
            merger_collection_STARS = load(joinpath(target, i), "merger_collection_STARS")
            merger_collection_GAS   = load(joinpath(target, i), "merger_collection_GAS")
            merger_collection_DM    = load(joinpath(target, i), "merger_collection_DM")
            #halo_story["J_DM_0"                     ] = missings(Float64, 3, length(halo_story["subID"]))
            #halo_story["j_DM_0"                     ] = missings(Float64, 3, length(halo_story["subID"]))
            #halo_story["J_GAS_0"                    ] = missings(Float64, 3, length(halo_story["subID"]))
            #halo_story["j_GAS_0"                    ] = missings(Float64, 3, length(halo_story["subID"]))
            #halo_story["J_GASvir_0"                 ] = missings(Float64, 3, length(halo_story["subID"]))
            #halo_story["j_GASvir_0"                 ] = missings(Float64, 3, length(halo_story["subID"]))
            #halo_story["J_STARS_0"                  ] = missings(Float64, 3, length(halo_story["subID"]))
            #halo_story["j_STARS_0"                  ] = missings(Float64, 3, length(halo_story["subID"]))
            #halo_story["J_STARSvir_0"               ] = missings(Float64, 3, length(halo_story["subID"]))
            #halo_story["j_STARSvir_0"               ] = missings(Float64, 3, length(halo_story["subID"]))
            #halo_story["merger_spin_map_STARS_0"    ] = missings(Float64, 18, length(halo_story["subID"]))
            #halo_story["merger_spin_map_GAS_0"      ] = missings(Float64, 18, length(halo_story["subID"]))
            #halo_story["merger_spin_map_DM_0"       ] = missings(Float64, 18, length(halo_story["subID"]))
            halo_story["J_DM_0"                     ] = deepcopy(halo_story["J_DM"                     ])
            halo_story["j_DM_0"                     ] = deepcopy(halo_story["j_DM"                     ])
            halo_story["J_GAS_0"                    ] = deepcopy(halo_story["J_GAS"                    ])
            halo_story["j_GAS_0"                    ] = deepcopy(halo_story["j_GAS"                    ])
            halo_story["J_GASvir_0"                 ] = deepcopy(halo_story["J_GASvir"                 ])
            halo_story["j_GASvir_0"                 ] = deepcopy(halo_story["j_GASvir"                 ])
            halo_story["J_STARS_0"                  ] = deepcopy(halo_story["J_STARS"                  ])
            halo_story["j_STARS_0"                  ] = deepcopy(halo_story["j_STARS"                  ])
            halo_story["J_STARSvir_0"               ] = deepcopy(halo_story["J_STARSvir"               ])
            halo_story["j_STARSvir_0"               ] = deepcopy(halo_story["j_STARSvir"               ])
            halo_story["merger_spin_map_STARS_0"    ] = deepcopy(halo_story["merger_spin_map_STARS"    ])
            halo_story["merger_spin_map_GAS_0"      ] = deepcopy(halo_story["merger_spin_map_GAS"      ])
            halo_story["merger_spin_map_DM_0"       ] = deepcopy(halo_story["merger_spin_map_DM"       ])

            for ii in 1:length(halo_story["subID"])
                halo_story["J_DM_0"                     ][:,ii] = halo_story["J_DM_0"                     ][:,ii] .* (sqrt(1 + halo_story["redshift"][ii]))
                halo_story["j_DM_0"                     ][:,ii] = halo_story["j_DM_0"                     ][:,ii] .* (sqrt(1 + halo_story["redshift"][ii]))
                halo_story["J_GAS_0"                    ][:,ii] = halo_story["J_GAS_0"                    ][:,ii] .* (sqrt(1 + halo_story["redshift"][ii]))
                halo_story["j_GAS_0"                    ][:,ii] = halo_story["j_GAS_0"                    ][:,ii] .* (sqrt(1 + halo_story["redshift"][ii]))
                halo_story["J_GASvir_0"                 ][:,ii] = halo_story["J_GASvir_0"                 ][:,ii] .* (sqrt(1 + halo_story["redshift"][ii]))
                halo_story["j_GASvir_0"                 ][:,ii] = halo_story["j_GASvir_0"                 ][:,ii] .* (sqrt(1 + halo_story["redshift"][ii]))
                halo_story["J_STARS_0"                  ][:,ii] = halo_story["J_STARS_0"                  ][:,ii] .* (sqrt(1 + halo_story["redshift"][ii]))
                halo_story["j_STARS_0"                  ][:,ii] = halo_story["j_STARS_0"                  ][:,ii] .* (sqrt(1 + halo_story["redshift"][ii]))
                halo_story["J_STARSvir_0"               ][:,ii] = halo_story["J_STARSvir_0"               ][:,ii] .* (sqrt(1 + halo_story["redshift"][ii]))
                halo_story["j_STARSvir_0"               ][:,ii] = halo_story["j_STARSvir_0"               ][:,ii] .* (sqrt(1 + halo_story["redshift"][ii]))
                # immediate
                if !ismissing(halo_story["merger_spin_map_STARS_0"    ][1,ii])
                    halo_story["merger_spin_map_STARS_0"    ][4:9,ii] = halo_story["merger_spin_map_STARS_0"    ][4:9,ii] .* (sqrt(1 + read_header("$simbox/groups_$(lpad(Int(halo_story["merger_spin_map_STARS_0"    ][1,ii]), 3, "0"))/sub_$(lpad(Int(halo_story["merger_spin_map_STARS_0"    ][1,ii]), 3, "0"))").z))
                    halo_story["merger_spin_map_GAS_0"      ][4:9,ii] = halo_story["merger_spin_map_GAS_0"      ][4:9,ii] .* (sqrt(1 + read_header("$simbox/groups_$(lpad(Int(halo_story["merger_spin_map_STARS_0"    ][1,ii]), 3, "0"))/sub_$(lpad(Int(halo_story["merger_spin_map_STARS_0"    ][1,ii]), 3, "0"))").z))
                    halo_story["merger_spin_map_DM_0"       ][4:9,ii] = halo_story["merger_spin_map_DM_0"       ][4:9,ii] .* (sqrt(1 + read_header("$simbox/groups_$(lpad(Int(halo_story["merger_spin_map_STARS_0"    ][1,ii]), 3, "0"))/sub_$(lpad(Int(halo_story["merger_spin_map_STARS_0"    ][1,ii]), 3, "0"))").z))
                end
                # earlier
                if !ismissing(halo_story["merger_spin_map_STARS_0"    ][10,ii])
                    halo_story["merger_spin_map_STARS_0"    ][13:18,ii] = halo_story["merger_spin_map_STARS_0"    ][13:18,ii] .* (sqrt(1 + read_header("$simbox/groups_$(lpad(Int(halo_story["merger_spin_map_STARS_0"    ][10,ii]), 3, "0"))/sub_$(lpad(Int(halo_story["merger_spin_map_STARS_0"    ][10,ii]), 3, "0"))").z))
                    halo_story["merger_spin_map_GAS_0"      ][13:18,ii] = halo_story["merger_spin_map_GAS_0"      ][13:18,ii] .* (sqrt(1 + read_header("$simbox/groups_$(lpad(Int(halo_story["merger_spin_map_STARS_0"    ][10,ii]), 3, "0"))/sub_$(lpad(Int(halo_story["merger_spin_map_STARS_0"    ][10,ii]), 3, "0"))").z))
                    halo_story["merger_spin_map_DM_0"       ][13:18,ii] = halo_story["merger_spin_map_DM_0"       ][13:18,ii] .* (sqrt(1 + read_header("$simbox/groups_$(lpad(Int(halo_story["merger_spin_map_STARS_0"    ][10,ii]), 3, "0"))/sub_$(lpad(Int(halo_story["merger_spin_map_STARS_0"    ][10,ii]), 3, "0"))").z))
                end
            end

            merger_collection_STARS, merger_collection_GAS, merger_collection_DM = collect_mergers(halo_story, min_time=min_time)
            save(joinpath(target, i), 
                "halo_story",   halo_story,
                "merger_collection_DM",     merger_collection_DM,
                "merger_collection_GAS",    merger_collection_GAS,
                "merger_collection_STARS",  merger_collection_STARS)
        end
    else
        error("   Unknown target type:$type \n   Accepted types: hs\n")
    end
    if verbose
        println("\n\n\n---------------------------\n\nNow witness the firepower of this fully armed and operational battle station!\n---------------------------\n\n")
        flush(stdout)
    end
    return nothing
end

print("'quickadd_feature'   ")



########################################################################################################################################
########################################################################################################################################
########################################################################################################################################
########################################################################################################################################
########################################################################################################################################
########################################################################################################################################
########################################################################################################################################
########################################################################################################################################
########################################################################################################################################



println()