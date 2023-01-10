
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



@doc """
DESCRIPTION:\n
INPUT:\n
OUTPUT:\n
""" ->
function select_root_halos(; start=" ", stop=" ", root_snap=136,
    root_mass_thr=0, part_type="stars", 
    outdir=current_dir_data, simbox=current_dir_simbox,
    iterstep=1, 
    verbose=true, skipexisting=true, output_type="file", filename_appendix="",
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
                    skipexisting    = $(skipexisting)\n
                    filename_appendix    = $(filename_appendix)\n
                    ")
        println("   Reading trees. ")
        flush(stdout)
    end
    trees = read_tree(simbox, root_snap)
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
    @showprogress for i in start:stop
        rootID  = get_haloid(trees.halos[i]).id-1
        smst    = convert_units_physical_mass(trees.halos[i].submasstab, root_head)
        if trees.halos[i].snapnum == root_snap && smst[convert_parttype_to_idx(part_type)] > root_mass_thr
            if verbose
                flush(stdout)
            end
            root_list["all_idx"     ] = vcat( root_list["all_idx"       ], i )
            root_list["rootID"      ] = vcat( root_list["rootID"        ], rootID )
            root_list["smst"        ] = hcat( root_list["smst"          ], smst )
            root_list["central_idx" ] = vcat( root_list["central_idx"   ], ifelse(rootID == get_first_subhalo(get_group(Galaxy(snapshot, rootID))).isub, 1, 0) )
        end
    end

    
    println("\n\n\n---------------------------\n\nNow witness the firepower of this fully armed and operational battle station!\n---------------------------\n\n")
    
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



@doc """
DESCRIPTION:\n
INPUT:\n
OUTPUT:\n
""" ->
function halo_stories(; subtype = "all", root_list = " ", 
    start=" ", stop=" ", 
    outdir=current_dir_data,
    root_mass_thr=0, part_type="stars", 
    simbox=current_dir_simbox, 
    verbose=true, skipexisting=false, 
    # root list settings
    root_snap=136, 
    )
    # Setup
    if verbose
        println("Initiating select_root_halos() with:\n 
                    start           = $(start)\n 
                    stop            = $(stop)\n 
                    root_snap       = $(root_snap)\n 
                    root_mass_thr   = $(root_mass_thr)\n 
                    outdir          = $(outdir)\n 
                    simbox          = $(simbox)\n 
                    root_list       = $(root_list)\n 
                    part_type       = $(part_type)\n 
                    subtype         = $(subtype)\n
                    verbose         = $(verbose)\n
                    skipexisting    = $(skipexisting)\n
                    ")
        flush(stdout)
    end

    ### Global Info
    trees       = read_tree(simbox, root_snap)
    subID_map    = read_all_nsubhalos(simbox)
    if root_list == " "
        root_list = select_root_halos(output_type="return",root_snap=root_snap, root_mass_thr=root_mass_thr, part_type=part_type, simbox=simbox, verbose=verbose)
    elseif typeof(root_list) == Dict{String, Any}
    elseif typeof(root_list) == String
        root_list == load(root_list, "root_list")
    else
        error("Unknown root_list: $root_list")
        flush(stdout)
    end

    ### Selection
    if occursin("central", lowercase(subtype))
        selection = find_indices( root_list["central_idx"], eq=1, comparewith=find_indices(root_list["smst"][:, convert_parttype_to_idx(part_type)], geq=root_mass_thr) )
    elseif occursin("all", lowercase(subtype))
        selection = find_indices(root_list["smst"][:, convert_parttype_to_idx(part_type)], geq=root_mass_thr)
    elseif occursin("non-central", lowercase(subtype))
        selection = find_indices(root_list["central_idx"], eq=0, comparewith=find_indices(root_list["smst"][:, convert_parttype_to_idx(part_type)], geq=root_mass_thr) )
    end


    error("Initiating main loop.")
    flush(stdout)
    # Loop over root halos
    for i in selection
        snapNR      = trees.halos[i].snapnum
        rootID      = get_global_index( subID_map, get_halonode(trees, snapNR, get_haloid(trees.halos[i])) )
        halo_story = Dict(
            "rootID"    => rootID, 
            "felixID"   => find_felixID(rootID), 
            "box"       => simbox, 
            "snapNR"      => Vector{Int64}(undef, 0), 
            "subID"     => Vector{Int64}(undef, 0), 
            "treeID"    => Vector{Int64}(undef, 0), 
            "FILE_NR"   => Vector{Int64}(undef, 0), 
            "REDSHIFT"  => Vector{Float64}(undef, 0), 
            "M_STARS"   => Vector{Float64}(undef, 0), 
            "M_STARS_2" => Vector{Float64}(undef, 0), 
            "M_GAS"     => Vector{Float64}(undef, 0), 
            "M_GAS_2"   => Vector{Float64}(undef, 0), 
            "M_DM"      => Vector{Float64}(undef, 0), 
            "M_DM_2"    => Vector{Float64}(undef, 0), 
            "POS"       => Array{Float64}(undef, 3, 0), 
            "MMP"       => Vector{Int64}(undef, 0), 
            "SWITCH"    => Vector{Int64}(undef, 0), 
            "BORDER"    => Vector{Int64}(undef, 0), 
            "JUMP"      => Vector{Int64}(undef, 0), 
            "EXCEED"    => Vector{Int64}(undef, 0), 

            "RVIR"          => missings(Float64, length(treefile_df[:, :I_SUB])),
            "RHMS_DM"       => missings(Float64, length(treefile_df[:, :I_SUB])),
            "J_DM"          => missings(Float64, 3, length(treefile_df[:, :I_SUB])),
            "j_DM"          => missings(Float64, 3, length(treefile_df[:, :I_SUB])),
            "RHMS_GAS"      => missings(Float64, length(treefile_df[:, :I_SUB])),
            "J_GAS"         => missings(Float64, 3, length(treefile_df[:, :I_SUB])),
            "j_GAS"         => missings(Float64, 3, length(treefile_df[:, :I_SUB])),
            "J_GASvir"      => missings(Float64, 3, length(treefile_df[:, :I_SUB])),
            "j_GASvir"      => missings(Float64, 3, length(treefile_df[:, :I_SUB])),
            "RHMS_STARS"    => missings(Float64, length(treefile_df[:, :I_SUB])),
            "J_STARS"       => missings(Float64, 3, length(treefile_df[:, :I_SUB])),
            "j_STARS"       => missings(Float64, 3, length(treefile_df[:, :I_SUB])),
            "J_STARSvir"    => missings(Float64, 3, length(treefile_df[:, :I_SUB])),
            "j_STARSvir"    => missings(Float64, 3, length(treefile_df[:, :I_SUB])),
            "j_orbital"     => missings(Float32, 3, length(treefile_df[:, :I_SUB])),
            "SFR"           => missings(Float64, length(treefile_df[:, :I_SUB])),
            "sSFR"          => missings(Float64, length(treefile_df[:, :I_SUB])),
            "BVAL"          => missings(Float64, length(treefile_df[:, :I_SUB])),
            "BVAL_0"        => missings(Float64, length(treefile_df[:, :I_SUB])))



            i = 1
            haloID      = get_haloid(trees.halos[i])
            println("Root ID: $haloID")
            snapNR      = trees.halos[i].snapnum
            subnode     = get_halonode(trees, snapNR, haloID)
            fp_idx      = subnode.first_progenitor
            while fp_idx > -1
                # fp update
                fp_snap     = trees.halos[fp_idx].snapnum
                snapNR  = fp_snap
        
                println("FP = $(trees.halos[fp_idx].subhalo_index) @ snap $(fp_snap)   ---   $(haloID)   ---   $(get_global_index(subID_map, subnode))")
                print(" NP: ")
                np_idx  = trees.halos[fp_idx].next_progenitor
                while np_idx > -1
                    np_node = get_halonode(trees, trees.halos[np_idx].snapnum, get_haloid(trees.halos[np_idx]))
                    if convert_units_physical_mass(trees.halos[np_idx].submasstab[5], read_header("$simbox/groups_$(lpad(snapNR, 3, "0"))/sub_$(lpad(snapNR, 3, "0"))")) ≥ 2e8
                        print("$(trees.halos[np_idx].subhalo_index) / $(get_global_index(subID_map, np_node))   ---   ")
                    end
                    np_idx  = np_node.next_progenitor
                end
                println()
        
                # move to to next Snap
                haloID      = get_haloid(trees.halos[fp_idx])
                subnode     = get_halonode(trees, fp_snap, haloID)
                fp_idx      = subnode.first_progenitor
            end
    end
    

    
    println("\n\n\n---------------------------\n\nNow witness the firepower of this fully armed and operational battle station!")
end

print("'halo_stories'   ")
