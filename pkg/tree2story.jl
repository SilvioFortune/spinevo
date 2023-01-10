
@doc """
DESCRIPTION:\n
INPUT:\n
OUTPUT:\n
""" ->
function tree2story(; start=" ", stop=" ", min_time=0.0, 
    outdir=current_dir_stories, indir=current_dir_trees, 
    simbox=current_dir_simbox, iterstep=1
    , verbose=true, skipexisting=true
    )
    # Setup
    halofiles   = readdir(indir)
    if start == " "
        start   = 1
    end
    if stop == " "
        stop    = length(halofiles)
    end

    println("Initiating tree2story with:\n 
                start = $(start)\n 
                stop = $(stop)\n 
                min_time = $(min_time)\n 
                outdir = $(outdir)\n 
                indir = $(indir)\n 
                simbox = $(simbox)\n 
                iterstep = $(iterstep)\n
                verbose = $(verbose)\n
                skipexisting = $(skipexisting)
                ")
        
    if skipexisting
        outfiles    = readdir(outdir)
    end

    if verbose
        println("\n\nInitiating main loop.\n")
        flush(stdout)
    end
    # Main
    for iii in start:iterstep:stop
        treefile_df = CSV.read(joinpath(indir, halofiles[iii]), DataFrame; delim=' ', ignorerepeated=true, header=2)
        halo_story = Dict(
            "rootID"    => parse(Int64, chop(halofiles[iii], head=5, tail=4)), 
            "box"       => simbox, 
            "snapNR"      => treefile_df[:, :SNAP], 
            "subID"     => treefile_df[:, :I_SUB], 
            "treeID"    => treefile_df[:, :I_TREE], 
            "fileNR"   => treefile_df[:, :FILE_NR], 
            "redshift"  => treefile_df[:, :REDSHIFT], 
            "M_STARS"   => treefile_df[:, :M_STARS], 
            "M_STARS_2" => treefile_df[:, :M_STAR_2], 
            "M_GAS"     => treefile_df[:, :M_GAS], 
            "M_GAS_2"   => treefile_df[:, :M_GAS_2], 
            "M_DM"      => treefile_df[:, :M_DM], 
            "M_DM_2"    => treefile_df[:, :M_DM_2], 
            "mmp"       => treefile_df[:, :mmp], 
            "switch"       => treefile_df[:, :SWITCH], 
            "border"       => treefile_df[:, :BORDER], 
            "jump"       => treefile_df[:, :JUMP], 
            "exceed"       => treefile_df[:, :EXCEED], 
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
            "BVAL"          => missings(Float64, length(treefile_df[:, :I_SUB])),
            "BVAL_0"        => missings(Float64, length(treefile_df[:, :I_SUB])))
        halo_story["mmp"] = Int.(halo_story["mmp"])
        
        if skipexisting
            if in(outfiles).("halo_$(halo_story["rootID"])_$(halo_story["subID"][1]).jld")
                if verbose
                    println("Existing file skipped: halo_$(halo_story["rootID"])_$(halo_story["subID"][1]).jld")
                end
                continue
            end
        end 
    
        # Final Subhalo First
        snapshot                = Snapshot(simbox, halo_story["snapNR"][1])
        g                       = Galaxy(snapshot, halo_story["subID"][1])
        halo_story["RVIR"][1]   = read_galaxy_prop(get_group(g), "RVIR", :physical)
        halo_story["SFR"][1]    = read_galaxy_prop(g, "SSFR", :physical)
        sph_small               = GadgetGalaxies.Sphere(0.1*halo_story["RVIR"][1])
        sph_large               = GadgetGalaxies.Sphere(halo_story["RVIR"][1])
        try
            @suppress begin
                read_halo!(g, units=:physical, props=((:gas, ["POS", "VEL", "MASS"]), (:stars, ["POS", "VEL", "MASS"])), radius_units=:physical, radius=halo_story["RVIR"][1]) 
                read_halo!(g, units=:physical, props=((:dm, ["POS", "VEL"]),), radius_units=:physical, radius=halo_story["RVIR"][1]) 
            end
            head    = read_header("$simbox/groups_$(@sprintf("%03i", halo_story["snapNR"][1]))/sub_$(@sprintf("%03i", halo_story["snapNR"][1]))")
            if GadgetGalaxies.total_mass(g.dm.pos, g.dm.mass) > 0.0
                halo_story["RHMS_DM"][1]    = half_mass_radius(g.dm)
                halo_story["J_DM"][:,1]     = angular_momentum(g.dm, sph_large)
                halo_story["j_DM"][:,1]     = specific_angular_momentum(g.dm, sph_large)
            end
            if GadgetGalaxies.total_mass(g.gas.pos, g.gas.mass) > 0.0
                halo_story["RHMS_GAS"][1]   = half_mass_radius(g.gas)
                halo_story["J_GAS"][:,1]    = angular_momentum(g.gas, sph_small)
                halo_story["j_GAS"][:,1]    = specific_angular_momentum(g.gas, sph_small)
                halo_story["J_GASvir"][:,1]    = angular_momentum(g.gas, sph_large)
                halo_story["j_GASvir"][:,1]    = specific_angular_momentum(g.gas, sph_large)
            end
            if GadgetGalaxies.total_mass(g.stars.pos, g.stars.mass) > 0.0
                halo_story["RHMS_STARS"][1] = half_mass_radius(g.stars)
                halo_story["J_STARS"][:,1]  = angular_momentum(g.stars, sph_small)
                halo_story["j_STARS"][:,1]  = specific_angular_momentum(g.stars, sph_small)
                halo_story["J_STARSvir"][:,1]    = angular_momentum(g.gas, sph_large)
                halo_story["j_STARSvir"][:,1]    = specific_angular_momentum(g.gas, sph_large)
                halo_story["BVAL"][1]       = b_value(g.stars, sph_small)
                halo_story["BVAL_0"][1]     = halo_story["BVAL"][1] + 0.5*log10(1+head.z)
            end
        catch
        end
    
        # Loop over Rest
        cent_id = 0
        for i in 2:length(halo_story["subID"])
            if halo_story["snapNR"][i] < halo_story["snapNR"][i-1]   # This corresponds to an earlier snap and therefore the most massive progenitor
                if verbose
                    println("$(halofiles[iii]) ($iii/$(stop))   ---   $i / $(length(treefile_df[:, :I_SUB]))\r")
                    flush(stdout)
                end
                cent_id                 = halo_story["subID"][i]
                snapshot                = Snapshot(simbox, halo_story["snapNR"][i])
                g                       = Galaxy(snapshot, halo_story["subID"][i])
                halo_story["RVIR"][i]   = read_galaxy_prop(get_group(g), "RVIR", :physical)
                sph_small               = GadgetGalaxies.Sphere(0.1*halo_story["RVIR"][i])
                sph_large               = GadgetGalaxies.Sphere(halo_story["RVIR"][i])
                try
                    @suppress begin
                        read_halo!(g, units=:physical, props=((:gas, ["POS", "VEL", "MASS"]), (:stars, ["POS", "VEL", "MASS"]), (:dm, ["POS", "VEL"])), radius_units=:physical, radius=halo_story["RVIR"][i]) 
                        #read_halo!(g, units=:physical, props=((:dm, ["POS", "VEL"]),), radius_units=:physical, radius=halo_story["RVIR"][i])  
                    end
                    head        = read_header("$simbox/groups_$(@sprintf("%03i", halo_story["snapNR"][i]))/sub_$(@sprintf("%03i", halo_story["snapNR"][i]))")
                    if GadgetGalaxies.total_mass(g.dm.pos, g.dm.mass) > 0.0
                        halo_story["RHMS_DM"][i]    = half_mass_radius(g.dm)
                        halo_story["J_DM"][:,i]     = angular_momentum(g.dm, sph_large)
                        halo_story["j_DM"][:,i]     = specific_angular_momentum(g.dm, sph_large)
                    end
                    if GadgetGalaxies.total_mass(g.gas.pos, g.gas.mass) > 0.0
                        halo_story["RHMS_GAS"][i]   = half_mass_radius(g.gas)
                        halo_story["J_GAS"][:,i]    = angular_momentum(g.gas, sph_small)
                        halo_story["j_GAS"][:,i]    = specific_angular_momentum(g.gas, sph_small)
                        halo_story["J_GASvir"][:,i]    = angular_momentum(g.gas, sph_large)
                        halo_story["j_GASvir"][:,i]    = specific_angular_momentum(g.gas, sph_large)
                    end
                    if GadgetGalaxies.total_mass(g.stars.pos, g.stars.mass) > 0.0
                        halo_story["RHMS_STARS"][i] = half_mass_radius(g.stars)
                        halo_story["J_STARS"][:,i]  = angular_momentum(g.stars, sph_small)
                        halo_story["j_STARS"][:,i]  = specific_angular_momentum(g.stars, sph_small)
                        halo_story["J_STARSvir"][:,i]  = angular_momentum(g.stars, sph_large)
                        halo_story["j_STARSvir"][:,i]  = specific_angular_momentum(g.stars, sph_large)
                        halo_story["BVAL"][i]       = b_value(g.stars, sph_small)
                        halo_story["BVAL_0"][i]     = halo_story["BVAL"][i] + 0.5*log10(1+head.z)
                    end
                    halo_story["SFR"][i]    = read_galaxy_prop(g, "SSFR", :physical)
                catch e
                    #halo_story["ERRMSG"][i] = string(e)
                    #println(e,"\r")
                end
            else   # Smaller progenitors
                halo_story["j_orbital"][:,i]    = orbit_j(halo_story["subID"][i], cent_id, halo_story["snapNR"][i])
                try
                    snapshot                = Snapshot(simbox, halo_story["snapNR"][i])
                    g                       = Galaxy(snapshot, halo_story["subID"][i])
                    halo_story["SFR"][i]    = read_galaxy_prop(g, "SSFR", :physical)
                catch e
                end
            end
        end
        


























        ##################################################################
        # Merger Data
        ##################################################################

        head        = read_header("$simbox/groups_$(@sprintf("%03i", halo_story["snapNR"][1]))/sub_$(@sprintf("%03i", halo_story["snapNR"][1]))")
        head2        = read_header("$simbox/groups_$(@sprintf("%03i", halo_story["snapNR"][1]-3))/sub_$(@sprintf("%03i", halo_story["snapNR"][1]-3))")
        id_mfelix   = convert_units_physical_mass(halo_story["M_STARS"][1], head)
        id_m2       = convert_units_physical_mass(halo_story["M_STARS_2"][1], head2)
        
        
        # Identify First Progenitors and mergers
    
        # Make lookbacktime array
        lbt = Array{Float64}(undef, 0)
        for i in 1:length(halo_story["snapNR"])
            head    = read_header("$simbox/groups_$(@sprintf("%03i", halo_story["snapNR"][i]))/sub_$(@sprintf("%03i", halo_story["snapNR"][i]))")
            lbt     = vcat( lbt, ustrip(lookback_time(cosmology(h=head.h0, OmegaM=head.omega_0), head.z)) )
        end
    
        # First Progenitors, backward in time
        fp_indices      = Array{Int64}(undef, 0)
        fp_snaps        = Array{Int64}(undef, 0)
        loop_length     = length(halo_story["snapNR"])
        for i in 1:loop_length
            #if i==1 || halo_story["snapNR"][i-1] > halo_story["snapNR"][i] && count(ismissing, halo_story["J_STARS"][:,i]) == 0 # First Progenitor
            if halo_story["mmp"][i] == 1 && count(ismissing, halo_story["J_STARS"][:,i]) == 0 # First Progenitor
                #print("FP@SNAP$(halo_story["snapNR"][i]) ")
                #flush(stdout)
                if length(fp_indices) == 0# && lbt[i]-lbt[1] > min_time
                    fp_indices  = vcat( fp_indices, i )
                    fp_snaps    = vcat( fp_snaps, halo_story["snapNR"][i] )
                elseif length(fp_indices) > 0 && lbt[i]-lbt[fp_indices[end]] > min_time
                    fp_indices  = vcat( fp_indices, i )
                    fp_snaps    = vcat( fp_snaps, halo_story["snapNR"][i] )
                end
            end
        end
        if length(fp_indices) ≤ 1 # skip if there is only one entry in the FP list
            continue
        end
        # Reverse to forward in time
        fp_indices      = fp_indices[end:-1:1]
        fp_snaps        = fp_snaps[end:-1:1]
        #println(fp_indices)
        #println(fp_snaps)
        
        merger_count        = 0
        n_mergers           = Array{Int64}(undef, 0)
        merger_index_map    = Array{Int64}(undef, 2, 0)    # Merger indices, same-snap FP indices
        fp_index    = 2
        n_mergers   = vcat( n_mergers, merger_count )
        # check if first snap already contains a first progenitor
        #println(fp_indices)
        #flush(stdout)
        #if fp_snaps[1] == halo_story["snapNR"][end]
            #fp_index    = 2
            #n_mergers   = vcat( n_mergers, merger_count )
        #else
            #fp_index    = 1
        #end

        # Merger Map, forward in time by starting at the bottom
        for i in loop_length:-1:1
            #print("$i ")
            #println("$i     $(fp_snaps[end]) > $(halo_story["snapNR"][i]) > $(fp_snaps[1])")
            if fp_snaps[end] > halo_story["snapNR"][i] ≥ fp_snaps[1]  # only accept entries within snaps with available stellar spins
                #println("   $(halo_story["snapNR"][i]) == $(fp_snaps[fp_index])")
                if halo_story["snapNR"][i] == fp_snaps[fp_index] # First Progenitor
                    # Check out Merger Info and assign to next
                    n_mergers       = vcat( n_mergers, merger_count )
                    merger_count    = 0
                    fp_index       += 1
                end
            
                # Identify Mergers
                #println("       $(halo_story["mmp"][i])     $(halo_story["snapNR"][i]) == $(halo_story["snapNR"][i-1])")
                if halo_story["mmp"][i] == 1 # FP
                elseif halo_story["mmp"][i] == 0 # secondary Progenitor
                    merger_count     += 1
                    merger_index_map  = hcat( merger_index_map, [i, ssFPfinder(i, halo_story)] )
                else
                    println("Error for i = $i, $(halo_story["snapNR"][i]), $(fp_snaps[fp_index]), $(halo_story["snapNR"][i-1]), $(halo_story["subID"][i-1])")
                end
            elseif fp_snaps[end] == halo_story["snapNR"][i] # final checkout
                n_mergers       = vcat( n_mergers, merger_count )
                break
            end
        end
        
        #println("$(length(fp_indices)) $(length(fp_snaps)) $(length(n_mergers))")
        #println("$(length(merger_index_map)) $(sum(n_mergers)) ")
        
        
        # Fill the STARS dictionary
    
        merger_collection_STARS = Dict(
                "N_MERGERS"     => n_mergers,
                "snapNR"          => halo_story["snapNR"][fp_indices],
                "ID_ISUB"       => Int.(ones(length(fp_indices)) .* halo_story["subID"][1]),
                "ID_Mfelix"     => ones(length(fp_indices)) .* id_mfelix,
                "ID_M2"         => ones(length(fp_indices)) .* id_m2,
                "subID"         => halo_story["subID"][fp_indices],
                "switch"        => halo_story["switch"][fp_indices],
                "border"        => halo_story["border"][fp_indices],
                "jump"          => halo_story["jump"][fp_indices],
                "exceed"        => halo_story["exceed"][fp_indices],
                "redshift"      => halo_story["redshift"][fp_indices],
                "BVAL"          => halo_story["BVAL"][fp_indices], 
                "BVAL_0"        => halo_story["BVAL_0"][fp_indices], 
                "J_main"        => halo_story["J_STARS"][:,fp_indices], 
                "j_main"        => halo_story["j_STARS"][:,fp_indices], 
                "J_vir"        => halo_story["J_STARSvir"][:,fp_indices], 
                "j_vir"        => halo_story["j_STARSvir"][:,fp_indices], 
                "LOOKBACKTIME"  => missings(Float64 , 0),
                "ΔM_felix"      => missings(Float64 , 0), 
                "ΔM2_felix"     => missings(Float64 , 0), 
                "ΔM_fromJ"      => missings(Float64 , 0), 
                "M_felix"       => missings(Float64 , 0), 
                "M2_felix"      => missings(Float64 , 0), 
                "M_fromJ"       => missings(Float64 , 0), 
                "ϕ_flip"        => missings(Float64 , 0), 
                "M_MERGERS"     => missings(Float64 , 0),
                "M_MISSED"      => missings(Float64 , 0), 
                "M_CONSIDERED"  => missings(Float64 , 0),  
                "M2_MERGERS"    => missings(Float64 , 0),
                "M2_MISSED"     => missings(Float64 , 0), 
                "M2_CONSIDERED" => missings(Float64 , 0),  
                "ΔBVAL"         => missings(Float64 , 0), 
                "ΔBVAL_0"       => missings(Float64 , 0), 
                "SFR"           => missings(Float64 , 0), 
                "M_MM"          => missings(Float64 , 2, 0), # 1 = mass, 2 = FPmass
                "M2_MM"         => missings(Float64 , 2, 0), # 1 = mass, 2 = FPmass
                "J_MMorbital"   => missings(Float64 , 3, 0), 
                "J_SUMorbital"  => missings(Float64 , 3, 0), 
                "ΔJ_main"       => missings(Float64 , 3, 0), 
                "Δj_main"       => missings(Float64 , 3, 0),
                "ΔJ_vir"       => missings(Float64 , 3, 0), 
                "Δj_vir"       => missings(Float64 , 3, 0),
                "Merger_Map"    => missings(Float64 , 8, 0)) # 1=mass, 2=mass2, 3=FPmass, 4=FPmass2, 5=SNAP, 6=z, 7=lbt, 8=subID
    
        for i in 1:length(fp_indices)
            #print("$i ")
            #flush(stdout)/home/moon/sfortune/spinevo/halostories_v20211204_min0.0Gyr/halo_1_1414.jld
    
            # Basic Info
            head    = read_header("$simbox/groups_$(@sprintf("%03i", halo_story["snapNR"][fp_indices[i]]))/sub_$(@sprintf("%03i", halo_story["snapNR"][fp_indices[i]]))")
            head2    = read_header("$simbox/groups_$(@sprintf("%03i", halo_story["snapNR"][fp_indices[i]]-3))/sub_$(@sprintf("%03i", halo_story["snapNR"][fp_indices[i]]-3))")
            merger_collection_STARS["LOOKBACKTIME"] = vcat( merger_collection_STARS["LOOKBACKTIME"], ustrip(lookback_time(cosmology(h=head.h0, OmegaM=head.omega_0), head.z)) )
            merger_collection_STARS["M_felix"]      = vcat( merger_collection_STARS["M_felix"], convert_units_physical_mass(halo_story["M_STARS"][fp_indices[i]], head) )
            merger_collection_STARS["M2_felix"]     = vcat( merger_collection_STARS["M2_felix"], convert_units_physical_mass(halo_story["M_STARS_2"][fp_indices[i]], head2) )
            merger_collection_STARS["SFR"]          = vcat( merger_collection_STARS["SFR"], halo_story["SFR"][fp_indices[i]] )
            
            if count(ismissing, halo_story["j_STARS"][:,fp_indices[i]]) == 0 
                merger_collection_STARS["M_fromJ"]  = vcat( merger_collection_STARS["M_fromJ"], norm(halo_story["J_STARS"][:,fp_indices[i]]) / norm(halo_story["j_STARS"][:,fp_indices[i]]) )
            else
                merger_collection_STARS["M_fromJ"]  = vcat( merger_collection_STARS["M_fromJ"], missing )
            end
            
            # Transitional Data
            if i == 1
                merger_collection_STARS["ΔJ_main"]          = hcat( merger_collection_STARS["ΔJ_main"], missings(Float64, 3) )
                merger_collection_STARS["Δj_main"]          = hcat( merger_collection_STARS["Δj_main"], missings(Float64, 3) )
                merger_collection_STARS["ΔJ_vir"]          = hcat( merger_collection_STARS["ΔJ_vir"], missings(Float64, 3) )
                merger_collection_STARS["Δj_vir"]          = hcat( merger_collection_STARS["Δj_vir"], missings(Float64, 3) )
                merger_collection_STARS["ΔBVAL"]            = vcat( merger_collection_STARS["ΔBVAL"], missing )
                merger_collection_STARS["ΔBVAL_0"]            = vcat( merger_collection_STARS["ΔBVAL_0"], missing )
                merger_collection_STARS["ϕ_flip"]           = vcat( merger_collection_STARS["ϕ_flip"], missing )
                merger_collection_STARS["ΔM_felix"]         = vcat( merger_collection_STARS["ΔM_felix"], missing )
                merger_collection_STARS["ΔM2_felix"]        = vcat( merger_collection_STARS["ΔM2_felix"], missing )
                merger_collection_STARS["ΔM_fromJ"]         = vcat( merger_collection_STARS["ΔM_fromJ"], missing )
                merger_collection_STARS["M_MM"]             = hcat( merger_collection_STARS["M_MM"], [missing, missing] )
                merger_collection_STARS["M2_MM"]            = hcat( merger_collection_STARS["M2_MM"], [missing, missing] )
                merger_collection_STARS["J_MMorbital"]      = hcat( merger_collection_STARS["J_MMorbital"], missings(Float64, 3) )
                merger_collection_STARS["J_SUMorbital"]     = hcat( merger_collection_STARS["J_SUMorbital"], missings(Float64, 3) )
                merger_collection_STARS["M_MERGERS"]        = vcat( merger_collection_STARS["M_MERGERS"], missing )
                merger_collection_STARS["M_MISSED"]         = vcat( merger_collection_STARS["M_MISSED"], missing )
                merger_collection_STARS["M_CONSIDERED"]     = vcat( merger_collection_STARS["M_CONSIDERED"], missing )
                merger_collection_STARS["M2_MERGERS"]       = vcat( merger_collection_STARS["M2_MERGERS"], missing )
                merger_collection_STARS["M2_MISSED"]        = vcat( merger_collection_STARS["M2_MISSED"], missing )
                merger_collection_STARS["M2_CONSIDERED"]    = vcat( merger_collection_STARS["M2_CONSIDERED"], missing )
            else
                merger_collection_STARS["ΔJ_main"]  = hcat( merger_collection_STARS["ΔJ_main"], halo_story["J_STARS"][:,fp_indices[i]] .- halo_story["J_STARS"][:,fp_indices[i-1]] )
                merger_collection_STARS["Δj_main"]  = hcat( merger_collection_STARS["Δj_main"], halo_story["j_STARS"][:,fp_indices[i]] .- halo_story["j_STARS"][:,fp_indices[i-1]] )
                merger_collection_STARS["ΔJ_vir"]  = hcat( merger_collection_STARS["ΔJ_vir"], halo_story["J_STARSvir"][:,fp_indices[i]] .- halo_story["J_STARSvir"][:,fp_indices[i-1]] )
                merger_collection_STARS["Δj_vir"]  = hcat( merger_collection_STARS["Δj_vir"], halo_story["j_STARSvir"][:,fp_indices[i]] .- halo_story["j_STARSvir"][:,fp_indices[i-1]] )
                merger_collection_STARS["ΔBVAL"]    = vcat( merger_collection_STARS["ΔBVAL"], halo_story["BVAL"][fp_indices[i]] - halo_story["BVAL"][fp_indices[i-1]] )
                merger_collection_STARS["ΔBVAL_0"]    = vcat( merger_collection_STARS["ΔBVAL_0"], halo_story["BVAL_0"][fp_indices[i]] - halo_story["BVAL_0"][fp_indices[i-1]] )
                if count(ismissing, halo_story["j_STARS"][:,fp_indices[i]]) == 0 && count(ismissing, halo_story["j_STARS"][:,fp_indices[i-1]]) == 0 
                    merger_collection_STARS["ϕ_flip"]   = vcat( merger_collection_STARS["ϕ_flip"], 180/π * angle(halo_story["J_STARS"][:,fp_indices[i]], halo_story["J_STARS"][:,fp_indices[i-1]] ) )
                    merger_collection_STARS["ΔM_fromJ"] = vcat( merger_collection_STARS["ΔM_fromJ"], (norm(halo_story["J_STARS"][:,fp_indices[i]]) / norm(halo_story["j_STARS"][:,fp_indices[i]])) - (norm(halo_story["J_STARS"][:,fp_indices[i-1]]) / norm(halo_story["j_STARS"][:,fp_indices[i-1]])) )
                else
                    #println(halo_story["j_STARS"][:,fp_indices[i]])
                    #println("missing j_STARS for $(storyfilelist[iii]) @ index $(fp_indices[i])")
                    merger_collection_STARS["ϕ_flip"]   = vcat( merger_collection_STARS["ϕ_flip"], missing )
                    merger_collection_STARS["ΔM_fromJ"] = vcat( merger_collection_STARS["ΔM_fromJ"], missing )
                end
                #println("\n$(halo_story["J_STARS"][:,fp_indices[i]]) $(halo_story["J_STARS"][:,fp_indices[i-1]])")
                merger_collection_STARS["ΔM_felix"] = vcat( merger_collection_STARS["ΔM_felix"], convert_units_physical_mass(halo_story["M_STARS"][fp_indices[i]], head) - convert_units_physical_mass(halo_story["M_STARS"][fp_indices[i-1]], head) )
                merger_collection_STARS["ΔM2_felix"] = vcat( merger_collection_STARS["ΔM2_felix"], convert_units_physical_mass(halo_story["M_STARS_2"][fp_indices[i]], head2) - convert_units_physical_mass(halo_story["M_STARS_2"][fp_indices[i-1]], head2) )
                
                # Merger Data
                # Setup
                merger_collection_STARS["M_MM"]             = vcat( merger_collection_STARS["M_MM"], 0 )
                merger_collection_STARS["M2_MM"]            = vcat( merger_collection_STARS["M2_MM"], 0 )
                merger_collection_STARS["J_MMorbital"]      = hcat( merger_collection_STARS["J_MMorbital"], zeros(3) )
                merger_collection_STARS["J_SUMorbital"]     = hcat( merger_collection_STARS["J_SUMorbital"], zeros(3) )
                merger_collection_STARS["M_MERGERS"]        = vcat( merger_collection_STARS["M_MERGERS"], 0 )
                merger_collection_STARS["M_MISSED"]         = vcat( merger_collection_STARS["M_MISSED"], 0 )
                merger_collection_STARS["M_CONSIDERED"]     = vcat( merger_collection_STARS["M_CONSIDERED"], 0 )
                merger_collection_STARS["M2_MERGERS"]       = vcat( merger_collection_STARS["M2_MERGERS"], 0 )
                merger_collection_STARS["M2_MISSED"]        = vcat( merger_collection_STARS["M2_MISSED"], 0 )
                merger_collection_STARS["M2_CONSIDERED"]    = vcat( merger_collection_STARS["M2_CONSIDERED"], 0 )
                # Actual check
                pos_in_mim  = 1+sum(n_mergers[1:i-1]) # position / index in merger index map
                for ii in merger_index_map[1,1+sum(n_mergers[1:i-1]):sum(n_mergers[1:i])]
                    head    = read_header("$simbox/groups_$(@sprintf("%03i", halo_story["snapNR"][ii]))/sub_$(@sprintf("%03i", halo_story["snapNR"][ii]))")
                    head2   = read_header("$simbox/groups_$(@sprintf("%03i", halo_story["snapNR"][ii]-3))/sub_$(@sprintf("%03i", halo_story["snapNR"][ii]-3))")
                    # Merger Map
                    #println("$(halo_story["snapNR"][ii])   i = $i")
                    #@show merger_index_map
                    #@show merger_index_map[1,1+sum(n_mergers[1:i-1]):sum(n_mergers[1:i])]
                    #flush(stdout)
                    merger_collection_STARS["Merger_Map"]   = hcat( merger_collection_STARS["Merger_Map"], 
                                                                        [   convert_units_physical_mass(halo_story["M_STARS"][ii], head),
                                                                            convert_units_physical_mass(halo_story["M_STARS_2"][ii], head2), 
                                                                            convert_units_physical_mass(halo_story["M_STARS"][merger_index_map[2,pos_in_mim]], head), 
                                                                            convert_units_physical_mass(halo_story["M_STARS_2"][merger_index_map[2,pos_in_mim]], head2), 
                                                                            halo_story["snapNR"][ii], 
                                                                            halo_story["redshift"][ii], 
                                                                            ustrip(lookback_time(cosmology(h=head.h0, OmegaM=head.omega_0), head.z)) , 
                                                                            halo_story["subID"][ii]
                                                                            ] )
                    # Most Massive condition
                    if convert_units_physical_mass(halo_story["M_STARS_2"][ii], head2) > merger_collection_STARS["M2_MM"][end]
                        merger_collection_STARS["M_MM"][end]            = convert_units_physical_mass(halo_story["M_STARS"][ii], head)
                        merger_collection_STARS["M2_MM"][end]           = convert_units_physical_mass(halo_story["M_STARS_2"][ii], head2)
                        merger_collection_STARS["J_MMorbital"][:,end]   = halo_story["j_orbital"][:,ii] .* convert_units_physical_mass(halo_story["M_STARS_2"][ii], head2)
                    end
                    # Orbital Data
                    if count(ismissing, halo_story["j_orbital"][:,ii]) == 0
                        merger_collection_STARS["J_SUMorbital"][:,end] .+= ( halo_story["j_orbital"][:,ii] .* convert_units_physical_mass(halo_story["M_STARS_2"][ii], head2) )
                        merger_collection_STARS["M_MERGERS"][end]       += convert_units_physical_mass(halo_story["M_STARS"][ii], head)
                        merger_collection_STARS["M_CONSIDERED"][end]    += convert_units_physical_mass(halo_story["M_STARS"][ii], head)
                        merger_collection_STARS["M2_MERGERS"][end]      += convert_units_physical_mass(halo_story["M_STARS_2"][ii], head2)
                        merger_collection_STARS["M2_CONSIDERED"][end]   += convert_units_physical_mass(halo_story["M_STARS_2"][ii], head2)
                        #println(halo_story["j_orbital"][:,ii])
                    else
                        #println("$ii")
                        merger_collection_STARS["M_MERGERS"][end]   += convert_units_physical_mass(halo_story["M_STARS"][ii], head)
                        merger_collection_STARS["M_MISSED"][end]    += convert_units_physical_mass(halo_story["M_STARS"][ii], head)
                        merger_collection_STARS["M2_MERGERS"][end]  += convert_units_physical_mass(halo_story["M_STARS_2"][ii], head2)
                        merger_collection_STARS["M2_MISSED"][end]   += convert_units_physical_mass(halo_story["M_STARS_2"][ii], head2)
                    end
                    pos_in_mim += 1
                end
            end
        end
        # Replace zeros with missings
        replace!(merger_collection_STARS["M_MM"]         , 0. => missing)
        replace!(merger_collection_STARS["M2_MM"]        , 0. => missing)
        replace!(merger_collection_STARS["J_MMorbital"]  , 0. => missing)
        replace!(merger_collection_STARS["J_SUMorbital"] , 0. => missing)
        replace!(merger_collection_STARS["M_MERGERS"]    , 0. => missing)
        replace!(merger_collection_STARS["M_MISSED"]     , 0. => missing)
        replace!(merger_collection_STARS["M_CONSIDERED"] , 0. => missing)
        replace!(merger_collection_STARS["M2_MERGERS"]   , 0. => missing)
        replace!(merger_collection_STARS["M2_MISSED"]    , 0. => missing)
        replace!(merger_collection_STARS["M2_CONSIDERED"], 0. => missing)
        #replace!(merger_collection_STARS["Merger_Map"][1:4,:], 0. => missing)
    
    
        # Fill the DM dictionary
    
        merger_collection_DM = Dict(
                "N_MERGERS"     => n_mergers,
                "snapNR"          => halo_story["snapNR"][fp_indices],
                "ID_ISUB"       => Int.(ones(length(fp_indices)) .* halo_story["subID"][1]),
                "ID_Mfelix"     => ones(length(fp_indices)) .* id_mfelix,
                "ID_M2"         => ones(length(fp_indices)) .* id_m2,
                "subID"         => halo_story["subID"][fp_indices],
                "switch"        => halo_story["switch"][fp_indices],
                "border"        => halo_story["border"][fp_indices],
                "jump"          => halo_story["jump"][fp_indices],
                "exceed"        => halo_story["exceed"][fp_indices],
                "redshift"      => halo_story["redshift"][fp_indices],
                "J_main"        => halo_story["J_DM"][:,fp_indices], 
                "j_main"        => halo_story["j_DM"][:,fp_indices], 
                "LOOKBACKTIME"  => missings(Float64 , 0),
                "ΔM_felix"      => missings(Float64 , 0), 
                "ΔM2_felix"     => missings(Float64 , 0), 
                "ΔM_fromJ"      => missings(Float64 , 0), 
                "M_felix"       => missings(Float64 , 0), 
                "M2_felix"      => missings(Float64 , 0), 
                "M_fromJ"       => missings(Float64 , 0), 
                "ϕ_flip"        => missings(Float64 , 0), 
                "M_MERGERS"     => missings(Float64 , 0),
                "M_MISSED"      => missings(Float64 , 0), 
                "M_CONSIDERED"  => missings(Float64 , 0),  
                "M2_MERGERS"    => missings(Float64 , 0),
                "M2_MISSED"     => missings(Float64 , 0), 
                "M2_CONSIDERED" => missings(Float64 , 0),  
                "M_MM"          => missings(Float64 , 2, 0), # 1 = mass, 2 = FPmass
                "M2_MM"         => missings(Float64 , 2, 0), # 1 = mass, 2 = FPmass
                "J_MMorbital"   => missings(Float64 , 3, 0), 
                "J_SUMorbital"  => missings(Float64 , 3, 0), 
                "ΔJ_main"       => missings(Float64 , 3, 0), 
                "Δj_main"       => missings(Float64 , 3, 0),
                "Merger_Map"    => missings(Float64 , 8, 0)) # 1=mass, 2=mass2, 3=FPmass, 4=FPmass2, 5=SNAP, 6=z, 7=lbt, 8=subID
    
        for i in 1:length(fp_indices)
            #print("$i ")
            #flush(stdout)/home/moon/sfortune/spinevo/halostories_v20211204_min0.0Gyr/halo_1_1414.jld
    
            # Basic Info
            head    = read_header("$simbox/groups_$(@sprintf("%03i", halo_story["snapNR"][fp_indices[i]]))/sub_$(@sprintf("%03i", halo_story["snapNR"][fp_indices[i]]))")
            head2    = read_header("$simbox/groups_$(@sprintf("%03i", halo_story["snapNR"][fp_indices[i]]-3))/sub_$(@sprintf("%03i", halo_story["snapNR"][fp_indices[i]]-3))")
            merger_collection_DM["LOOKBACKTIME"] = vcat( merger_collection_DM["LOOKBACKTIME"], ustrip(lookback_time(cosmology(h=head.h0, OmegaM=head.omega_0), head.z)) )
            merger_collection_DM["M_felix"]      = vcat( merger_collection_DM["M_felix"], convert_units_physical_mass(halo_story["M_DM"][fp_indices[i]], head) )
            merger_collection_DM["M2_felix"]     = vcat( merger_collection_DM["M2_felix"], convert_units_physical_mass(halo_story["M_DM_2"][fp_indices[i]], head2) )
            if count(ismissing, halo_story["j_DM"][:,fp_indices[i]]) == 0 
                merger_collection_DM["M_fromJ"]  = vcat( merger_collection_DM["M_fromJ"], norm(halo_story["J_DM"][:,fp_indices[i]]) / norm(halo_story["j_DM"][:,fp_indices[i]]) )
            else
                merger_collection_DM["M_fromJ"]  = vcat( merger_collection_DM["M_fromJ"], missing )
            end
            
            # Transitional Data
            if i == 1
                merger_collection_DM["ΔJ_main"]          = hcat( merger_collection_DM["ΔJ_main"], missings(Float64, 3) )
                merger_collection_DM["Δj_main"]          = hcat( merger_collection_DM["Δj_main"], missings(Float64, 3) )
                merger_collection_DM["ϕ_flip"]           = vcat( merger_collection_DM["ϕ_flip"], missing )
                merger_collection_DM["ΔM_felix"]         = vcat( merger_collection_DM["ΔM_felix"], missing )
                merger_collection_DM["ΔM2_felix"]        = vcat( merger_collection_DM["ΔM2_felix"], missing )
                merger_collection_DM["ΔM_fromJ"]         = vcat( merger_collection_DM["ΔM_fromJ"], missing )
                merger_collection_DM["M_MM"]             = hcat( merger_collection_DM["M_MM"], [missing, missing] )
                merger_collection_DM["M2_MM"]            = hcat( merger_collection_DM["M2_MM"], [missing, missing] )
                merger_collection_DM["J_MMorbital"]      = hcat( merger_collection_DM["J_MMorbital"], missings(Float64, 3) )
                merger_collection_DM["J_SUMorbital"]     = hcat( merger_collection_DM["J_SUMorbital"], missings(Float64, 3) )
                merger_collection_DM["M_MERGERS"]        = vcat( merger_collection_DM["M_MERGERS"], missing )
                merger_collection_DM["M_MISSED"]         = vcat( merger_collection_DM["M_MISSED"], missing )
                merger_collection_DM["M_CONSIDERED"]     = vcat( merger_collection_DM["M_CONSIDERED"], missing )
                merger_collection_DM["M2_MERGERS"]       = vcat( merger_collection_DM["M2_MERGERS"], missing )
                merger_collection_DM["M2_MISSED"]        = vcat( merger_collection_DM["M2_MISSED"], missing )
                merger_collection_DM["M2_CONSIDERED"]    = vcat( merger_collection_DM["M2_CONSIDERED"], missing )
            else
                merger_collection_DM["ΔJ_main"]  = hcat( merger_collection_DM["ΔJ_main"], halo_story["J_DM"][:,fp_indices[i]] .- halo_story["J_DM"][:,fp_indices[i-1]] )
                merger_collection_DM["Δj_main"]  = hcat( merger_collection_DM["Δj_main"], halo_story["j_DM"][:,fp_indices[i]] .- halo_story["j_DM"][:,fp_indices[i-1]] )
                if count(ismissing, halo_story["j_DM"][:,fp_indices[i]]) == 0 && count(ismissing, halo_story["j_DM"][:,fp_indices[i-1]]) == 0 
                    merger_collection_DM["ϕ_flip"]   = vcat( merger_collection_DM["ϕ_flip"], 180/π * angle( halo_story["J_DM"][:,fp_indices[i]], halo_story["J_DM"][:,fp_indices[i-1]] ) )
                    merger_collection_DM["ΔM_fromJ"] = vcat( merger_collection_DM["ΔM_fromJ"], (norm(halo_story["J_DM"][:,fp_indices[i]]) / norm(halo_story["j_DM"][:,fp_indices[i]])) - (norm(halo_story["J_DM"][:,fp_indices[i-1]]) / norm(halo_story["j_DM"][:,fp_indices[i-1]])) )
                else
                    #println(halo_story["j_DM"][:,fp_indices[i]])
                    #println("missing j_DM for $(storyfilelist[iii]) @ index $(fp_indices[i])")
                    merger_collection_DM["ϕ_flip"]   = vcat( merger_collection_DM["ϕ_flip"], missing )
                    merger_collection_DM["ΔM_fromJ"] = vcat( merger_collection_DM["ΔM_fromJ"], missing )
                end
                #println("\n$(halo_story["J_DM"][:,fp_indices[i]]) $(halo_story["J_DM"][:,fp_indices[i-1]])")
                merger_collection_DM["ΔM_felix"] = vcat( merger_collection_DM["ΔM_felix"], convert_units_physical_mass(halo_story["M_DM"][fp_indices[i]], head) - convert_units_physical_mass(halo_story["M_DM"][fp_indices[i-1]], head) )
                merger_collection_DM["ΔM2_felix"] = vcat( merger_collection_DM["ΔM2_felix"], convert_units_physical_mass(halo_story["M_DM_2"][fp_indices[i]], head2) - convert_units_physical_mass(halo_story["M_DM_2"][fp_indices[i-1]], head2) )
                
                # Merger Data
                # Setup
                merger_collection_DM["M_MM"]             = vcat( merger_collection_DM["M_MM"], 0 )
                merger_collection_DM["M2_MM"]            = vcat( merger_collection_DM["M2_MM"], 0 )
                merger_collection_DM["J_MMorbital"]      = hcat( merger_collection_DM["J_MMorbital"], zeros(3) )
                merger_collection_DM["J_SUMorbital"]     = hcat( merger_collection_DM["J_SUMorbital"], zeros(3) )
                merger_collection_DM["M_MERGERS"]        = vcat( merger_collection_DM["M_MERGERS"], 0 )
                merger_collection_DM["M_MISSED"]         = vcat( merger_collection_DM["M_MISSED"], 0 )
                merger_collection_DM["M_CONSIDERED"]     = vcat( merger_collection_DM["M_CONSIDERED"], 0 )
                merger_collection_DM["M2_MERGERS"]       = vcat( merger_collection_DM["M2_MERGERS"], 0 )
                merger_collection_DM["M2_MISSED"]        = vcat( merger_collection_DM["M2_MISSED"], 0 )
                merger_collection_DM["M2_CONSIDERED"]    = vcat( merger_collection_DM["M2_CONSIDERED"], 0 )
                # Actual check
                pos_in_mim  = 1+sum(n_mergers[1:i-1]) # position / index in merger index map
                for ii in merger_index_map[1,1+sum(n_mergers[1:i-1]):sum(n_mergers[1:i])]
                    head    = read_header("$simbox/groups_$(@sprintf("%03i", halo_story["snapNR"][ii]))/sub_$(@sprintf("%03i", halo_story["snapNR"][ii]))")
                    head2   = read_header("$simbox/groups_$(@sprintf("%03i", halo_story["snapNR"][ii]-3))/sub_$(@sprintf("%03i", halo_story["snapNR"][ii]-3))")
                    # Merger Map
                    #println("$(halo_story["snapNR"][ii])   i = $i")
                    #@show merger_index_map
                    #@show merger_index_map[1,1+sum(n_mergers[1:i-1]):sum(n_mergers[1:i])]
                    #flush(stdout)
                    merger_collection_DM["Merger_Map"]   = hcat( merger_collection_DM["Merger_Map"], 
                                                                        [   convert_units_physical_mass(halo_story["M_DM"][ii], head),
                                                                            convert_units_physical_mass(halo_story["M_DM_2"][ii], head2), 
                                                                            convert_units_physical_mass(halo_story["M_DM"][merger_index_map[2,pos_in_mim]], head), 
                                                                            convert_units_physical_mass(halo_story["M_DM_2"][merger_index_map[2,pos_in_mim]], head2), 
                                                                            halo_story["snapNR"][ii], 
                                                                            halo_story["redshift"][ii], 
                                                                            ustrip(lookback_time(cosmology(h=head.h0, OmegaM=head.omega_0), head.z)) , 
                                                                            halo_story["subID"][ii]
                                                                            ] )
                    # Most Massive condition
                    if convert_units_physical_mass(halo_story["M_DM_2"][ii], head2) > merger_collection_DM["M2_MM"][end]
                        merger_collection_DM["M_MM"][end]            = convert_units_physical_mass(halo_story["M_DM"][ii], head)
                        merger_collection_DM["M2_MM"][end]           = convert_units_physical_mass(halo_story["M_DM_2"][ii], head2)
                        merger_collection_DM["J_MMorbital"][:,end]   = halo_story["j_orbital"][:,ii] .* convert_units_physical_mass(halo_story["M_DM_2"][ii], head2)
                    end
                    # Orbital Data
                    if count(ismissing, halo_story["j_orbital"][:,ii]) == 0
                        merger_collection_DM["J_SUMorbital"][:,end] .+= ( halo_story["j_orbital"][:,ii] .* convert_units_physical_mass(halo_story["M_DM_2"][ii], head2) )
                        merger_collection_DM["M_MERGERS"][end]       += convert_units_physical_mass(halo_story["M_DM"][ii], head)
                        merger_collection_DM["M_CONSIDERED"][end]    += convert_units_physical_mass(halo_story["M_DM"][ii], head)
                        merger_collection_DM["M2_MERGERS"][end]      += convert_units_physical_mass(halo_story["M_DM_2"][ii], head2)
                        merger_collection_DM["M2_CONSIDERED"][end]   += convert_units_physical_mass(halo_story["M_DM_2"][ii], head2)
                        #println(halo_story["j_orbital"][:,ii])
                    else
                        #println("$ii")
                        merger_collection_DM["M_MERGERS"][end]   += convert_units_physical_mass(halo_story["M_DM"][ii], head)
                        merger_collection_DM["M_MISSED"][end]    += convert_units_physical_mass(halo_story["M_DM"][ii], head)
                        merger_collection_DM["M2_MERGERS"][end]  += convert_units_physical_mass(halo_story["M_DM_2"][ii], head2)
                        merger_collection_DM["M2_MISSED"][end]   += convert_units_physical_mass(halo_story["M_DM_2"][ii], head2)
                    end
                    pos_in_mim += 1
                end
            end
        end
        # Replace zeros with missings
        replace!(merger_collection_DM["M_MM"]         , 0. => missing)
        replace!(merger_collection_DM["M2_MM"]        , 0. => missing)
        replace!(merger_collection_DM["J_MMorbital"]  , 0. => missing)
        replace!(merger_collection_DM["J_SUMorbital"] , 0. => missing)
        replace!(merger_collection_DM["M_MERGERS"]    , 0. => missing)
        replace!(merger_collection_DM["M_MISSED"]     , 0. => missing)
        replace!(merger_collection_DM["M_CONSIDERED"] , 0. => missing)
        replace!(merger_collection_DM["M2_MERGERS"]   , 0. => missing)
        replace!(merger_collection_DM["M2_MISSED"]    , 0. => missing)
        replace!(merger_collection_DM["M2_CONSIDERED"], 0. => missing)
        #replace!(merger_collection_DM["Merger_Map"][1:4,:], 0. => missing)
    
    
        # Fill the GAS dictionary
    
        merger_collection_GAS = Dict(
                "N_MERGERS"     => n_mergers,
                "snapNR"          => halo_story["snapNR"][fp_indices],
                "ID_ISUB"       => Int.(ones(length(fp_indices)) .* halo_story["subID"][1]),
                "ID_Mfelix"     => ones(length(fp_indices)) .* id_mfelix,
                "ID_M2"         => ones(length(fp_indices)) .* id_m2,
                "subID"         => halo_story["subID"][fp_indices],
                "switch"        => halo_story["switch"][fp_indices],
                "border"        => halo_story["border"][fp_indices],
                "jump"          => halo_story["jump"][fp_indices],
                "exceed"        => halo_story["exceed"][fp_indices],
                "redshift"      => halo_story["redshift"][fp_indices],
                "BVAL"          => halo_story["BVAL"][fp_indices], 
                "BVAL_0"        => halo_story["BVAL_0"][fp_indices], 
                "J_main"        => halo_story["J_GAS"][:,fp_indices], 
                "j_main"        => halo_story["j_GAS"][:,fp_indices], 
                "J_vir"        => halo_story["J_GASvir"][:,fp_indices], 
                "j_vir"        => halo_story["j_GASvir"][:,fp_indices], 
                "LOOKBACKTIME"  => missings(Float64 , 0),
                "ΔM_felix"      => missings(Float64 , 0), 
                "ΔM2_felix"     => missings(Float64 , 0), 
                "ΔM_fromJ"      => missings(Float64 , 0), 
                "M_felix"       => missings(Float64 , 0), 
                "M2_felix"      => missings(Float64 , 0), 
                "M_fromJ"       => missings(Float64 , 0), 
                "ϕ_flip"        => missings(Float64 , 0), 
                "M_MERGERS"     => missings(Float64 , 0),
                "M_MISSED"      => missings(Float64 , 0), 
                "M_CONSIDERED"  => missings(Float64 , 0),  
                "M2_MERGERS"    => missings(Float64 , 0),
                "M2_MISSED"     => missings(Float64 , 0), 
                "M2_CONSIDERED" => missings(Float64 , 0),  
                "ΔBVAL"         => missings(Float64 , 0), 
                "ΔBVAL_0"         => missings(Float64 , 0), 
                "M_MM"          => missings(Float64 , 2, 0), # 1 = mass, 2 = FPmass
                "M2_MM"         => missings(Float64 , 2, 0), # 1 = mass, 2 = FPmass
                "J_MMorbital"   => missings(Float64 , 3, 0), 
                "J_SUMorbital"  => missings(Float64 , 3, 0), 
                "ΔJ_main"       => missings(Float64 , 3, 0), 
                "Δj_main"       => missings(Float64 , 3, 0),
                "ΔJ_vir"       => missings(Float64 , 3, 0), 
                "Δj_vir"       => missings(Float64 , 3, 0),
                "Merger_Map"    => missings(Float64 , 8, 0)) # 1=mass, 2=mass2, 3=FPmass, 4=FPmass2, 5=SNAP, 6=z, 7=lbt, 8=subID
    
        for i in 1:length(fp_indices)
            #print("$i ")
            #flush(stdout)/home/moon/sfortune/spinevo/halostories_v20211204_min0.0Gyr/halo_1_1414.jld
    
            # Basic Info
            head    = read_header("$simbox/groups_$(@sprintf("%03i", halo_story["snapNR"][fp_indices[i]]))/sub_$(@sprintf("%03i", halo_story["snapNR"][fp_indices[i]]))")
            head2    = read_header("$simbox/groups_$(@sprintf("%03i", halo_story["snapNR"][fp_indices[i]]-3))/sub_$(@sprintf("%03i", halo_story["snapNR"][fp_indices[i]]-3))")
            merger_collection_GAS["LOOKBACKTIME"] = vcat( merger_collection_GAS["LOOKBACKTIME"], ustrip(lookback_time(cosmology(h=head.h0, OmegaM=head.omega_0), head.z)) )
            merger_collection_GAS["M_felix"]      = vcat( merger_collection_GAS["M_felix"], convert_units_physical_mass(halo_story["M_GAS"][fp_indices[i]], head) )
            merger_collection_GAS["M2_felix"]     = vcat( merger_collection_GAS["M2_felix"], convert_units_physical_mass(halo_story["M_GAS_2"][fp_indices[i]], head2) )
            if count(ismissing, halo_story["j_GAS"][:,fp_indices[i]]) == 0 
                merger_collection_GAS["M_fromJ"]  = vcat( merger_collection_GAS["M_fromJ"], norm(halo_story["J_GAS"][:,fp_indices[i]]) / norm(halo_story["j_GAS"][:,fp_indices[i]]) )
            else
                merger_collection_GAS["M_fromJ"]  = vcat( merger_collection_GAS["M_fromJ"], missing )
            end
            
            # Transitional Data
            if i == 1
                merger_collection_GAS["ΔJ_main"]          = hcat( merger_collection_GAS["ΔJ_main"], missings(Float64, 3) )
                merger_collection_GAS["Δj_main"]          = hcat( merger_collection_GAS["Δj_main"], missings(Float64, 3) )
                merger_collection_GAS["ΔJ_vir"]          = hcat( merger_collection_GAS["ΔJ_vir"], missings(Float64, 3) )
                merger_collection_GAS["Δj_vir"]          = hcat( merger_collection_GAS["Δj_vir"], missings(Float64, 3) )
                merger_collection_GAS["ΔBVAL"]            = vcat( merger_collection_GAS["ΔBVAL"], missing )
                merger_collection_GAS["ΔBVAL_0"]            = vcat( merger_collection_GAS["ΔBVAL_0"], missing )
                merger_collection_GAS["ϕ_flip"]           = vcat( merger_collection_GAS["ϕ_flip"], missing )
                merger_collection_GAS["ΔM_felix"]         = vcat( merger_collection_GAS["ΔM_felix"], missing )
                merger_collection_GAS["ΔM2_felix"]        = vcat( merger_collection_GAS["ΔM2_felix"], missing )
                merger_collection_GAS["ΔM_fromJ"]         = vcat( merger_collection_GAS["ΔM_fromJ"], missing )
                merger_collection_GAS["M_MM"]             = hcat( merger_collection_GAS["M_MM"], [missing, missing] )
                merger_collection_GAS["M2_MM"]            = hcat( merger_collection_GAS["M2_MM"], [missing, missing] )
                merger_collection_GAS["J_MMorbital"]      = hcat( merger_collection_GAS["J_MMorbital"], missings(Float64, 3) )
                merger_collection_GAS["J_SUMorbital"]     = hcat( merger_collection_GAS["J_SUMorbital"], missings(Float64, 3) )
                merger_collection_GAS["M_MERGERS"]        = vcat( merger_collection_GAS["M_MERGERS"], missing )
                merger_collection_GAS["M_MISSED"]         = vcat( merger_collection_GAS["M_MISSED"], missing )
                merger_collection_GAS["M_CONSIDERED"]     = vcat( merger_collection_GAS["M_CONSIDERED"], missing )
                merger_collection_GAS["M2_MERGERS"]       = vcat( merger_collection_GAS["M2_MERGERS"], missing )
                merger_collection_GAS["M2_MISSED"]        = vcat( merger_collection_GAS["M2_MISSED"], missing )
                merger_collection_GAS["M2_CONSIDERED"]    = vcat( merger_collection_GAS["M2_CONSIDERED"], missing )
            else
                merger_collection_GAS["ΔJ_main"]  = hcat( merger_collection_GAS["ΔJ_main"], halo_story["J_GAS"][:,fp_indices[i]] .- halo_story["J_GAS"][:,fp_indices[i-1]] )
                merger_collection_GAS["Δj_main"]  = hcat( merger_collection_GAS["Δj_main"], halo_story["j_GAS"][:,fp_indices[i]] .- halo_story["j_GAS"][:,fp_indices[i-1]] )
                merger_collection_GAS["ΔJ_vir"]  = hcat( merger_collection_GAS["ΔJ_vir"], halo_story["J_GASvir"][:,fp_indices[i]] .- halo_story["J_GASvir"][:,fp_indices[i-1]] )
                merger_collection_GAS["Δj_vir"]  = hcat( merger_collection_GAS["Δj_vir"], halo_story["j_GASvir"][:,fp_indices[i]] .- halo_story["j_GASvir"][:,fp_indices[i-1]] )
                merger_collection_GAS["ΔBVAL"]    = vcat( merger_collection_GAS["ΔBVAL"], halo_story["BVAL"][fp_indices[i]] - halo_story["BVAL"][fp_indices[i-1]] )
                merger_collection_GAS["ΔBVAL_0"]    = vcat( merger_collection_GAS["ΔBVAL_0"], halo_story["BVAL_0"][fp_indices[i]] - halo_story["BVAL_0"][fp_indices[i-1]] )
                if count(ismissing, halo_story["j_GAS"][:,fp_indices[i]]) == 0 && count(ismissing, halo_story["j_GAS"][:,fp_indices[i-1]]) == 0 
                    merger_collection_GAS["ϕ_flip"]   = vcat( merger_collection_GAS["ϕ_flip"], 180/π * angle( halo_story["J_GAS"][:,fp_indices[i]], halo_story["J_GAS"][:,fp_indices[i-1]] ) )
                    merger_collection_GAS["ΔM_fromJ"] = vcat( merger_collection_GAS["ΔM_fromJ"], (norm(halo_story["J_GAS"][:,fp_indices[i]]) / norm(halo_story["j_GAS"][:,fp_indices[i]])) - (norm(halo_story["J_GAS"][:,fp_indices[i-1]]) / norm(halo_story["j_GAS"][:,fp_indices[i-1]])) )
                else
                    #println(halo_story["j_GAS"][:,fp_indices[i]])
                    #println("missing j_GAS for $(storyfilelist[iii]) @ index $(fp_indices[i])")
                    merger_collection_GAS["ϕ_flip"]   = vcat( merger_collection_GAS["ϕ_flip"], missing )
                    merger_collection_GAS["ΔM_fromJ"] = vcat( merger_collection_GAS["ΔM_fromJ"], missing )
                end
                #println("\n$(halo_story["J_GAS"][:,fp_indices[i]]) $(halo_story["J_GAS"][:,fp_indices[i-1]])")
                merger_collection_GAS["ΔM_felix"] = vcat( merger_collection_GAS["ΔM_felix"], convert_units_physical_mass(halo_story["M_GAS"][fp_indices[i]], head) - convert_units_physical_mass(halo_story["M_GAS"][fp_indices[i-1]], head) )
                merger_collection_GAS["ΔM2_felix"] = vcat( merger_collection_GAS["ΔM2_felix"], convert_units_physical_mass(halo_story["M_GAS_2"][fp_indices[i]], head2) - convert_units_physical_mass(halo_story["M_GAS_2"][fp_indices[i-1]], head2) )
                
                # Merger Data
                # Setup
                merger_collection_GAS["M_MM"]             = vcat( merger_collection_GAS["M_MM"], 0 )
                merger_collection_GAS["M2_MM"]            = vcat( merger_collection_GAS["M2_MM"], 0 )
                merger_collection_GAS["J_MMorbital"]      = hcat( merger_collection_GAS["J_MMorbital"], zeros(3) )
                merger_collection_GAS["J_SUMorbital"]     = hcat( merger_collection_GAS["J_SUMorbital"], zeros(3) )
                merger_collection_GAS["M_MERGERS"]        = vcat( merger_collection_GAS["M_MERGERS"], 0 )
                merger_collection_GAS["M_MISSED"]         = vcat( merger_collection_GAS["M_MISSED"], 0 )
                merger_collection_GAS["M_CONSIDERED"]     = vcat( merger_collection_GAS["M_CONSIDERED"], 0 )
                merger_collection_GAS["M2_MERGERS"]       = vcat( merger_collection_GAS["M2_MERGERS"], 0 )
                merger_collection_GAS["M2_MISSED"]        = vcat( merger_collection_GAS["M2_MISSED"], 0 )
                merger_collection_GAS["M2_CONSIDERED"]    = vcat( merger_collection_GAS["M2_CONSIDERED"], 0 )
                # Actual check
                pos_in_mim  = 1+sum(n_mergers[1:i-1]) # position / index in merger index map
                for ii in merger_index_map[1,1+sum(n_mergers[1:i-1]):sum(n_mergers[1:i])]
                    head    = read_header("$simbox/groups_$(@sprintf("%03i", halo_story["snapNR"][ii]))/sub_$(@sprintf("%03i", halo_story["snapNR"][ii]))")
                    head2   = read_header("$simbox/groups_$(@sprintf("%03i", halo_story["snapNR"][ii]-3))/sub_$(@sprintf("%03i", halo_story["snapNR"][ii]-3))")
                    # Merger Map
                    #println("$(halo_story["snapNR"][ii])   i = $i")
                    #@show merger_index_map
                    #@show merger_index_map[1,1+sum(n_mergers[1:i-1]):sum(n_mergers[1:i])]
                    #flush(stdout)
                    merger_collection_GAS["Merger_Map"]   = hcat( merger_collection_GAS["Merger_Map"], 
                                                                        [   convert_units_physical_mass(halo_story["M_GAS"][ii], head),
                                                                            convert_units_physical_mass(halo_story["M_GAS_2"][ii], head2), 
                                                                            convert_units_physical_mass(halo_story["M_GAS"][merger_index_map[2,pos_in_mim]], head), 
                                                                            convert_units_physical_mass(halo_story["M_GAS_2"][merger_index_map[2,pos_in_mim]], head2), 
                                                                            halo_story["snapNR"][ii], 
                                                                            halo_story["redshift"][ii], 
                                                                            ustrip(lookback_time(cosmology(h=head.h0, OmegaM=head.omega_0), head.z)) , 
                                                                            halo_story["subID"][ii]
                                                                            ] )
                    # Most Massive condition
                    if convert_units_physical_mass(halo_story["M_GAS_2"][ii], head2) > merger_collection_GAS["M2_MM"][end]
                        merger_collection_GAS["M_MM"][end]            = convert_units_physical_mass(halo_story["M_GAS"][ii], head)
                        merger_collection_GAS["M2_MM"][end]           = convert_units_physical_mass(halo_story["M_GAS_2"][ii], head2)
                        merger_collection_GAS["J_MMorbital"][:,end]   = halo_story["j_orbital"][:,ii] .* convert_units_physical_mass(halo_story["M_GAS_2"][ii], head2)
                    end
                    # Orbital Data
                    if count(ismissing, halo_story["j_orbital"][:,ii]) == 0
                        merger_collection_GAS["J_SUMorbital"][:,end] .+= ( halo_story["j_orbital"][:,ii] .* convert_units_physical_mass(halo_story["M_GAS_2"][ii], head2) )
                        merger_collection_GAS["M_MERGERS"][end]       += convert_units_physical_mass(halo_story["M_GAS"][ii], head)
                        merger_collection_GAS["M_CONSIDERED"][end]    += convert_units_physical_mass(halo_story["M_GAS"][ii], head)
                        merger_collection_GAS["M2_MERGERS"][end]      += convert_units_physical_mass(halo_story["M_GAS_2"][ii], head2)
                        merger_collection_GAS["M2_CONSIDERED"][end]   += convert_units_physical_mass(halo_story["M_GAS_2"][ii], head2)
                        #println(halo_story["j_orbital"][:,ii])
                    else
                        #println("$ii")
                        merger_collection_GAS["M_MERGERS"][end]   += convert_units_physical_mass(halo_story["M_GAS"][ii], head)
                        merger_collection_GAS["M_MISSED"][end]    += convert_units_physical_mass(halo_story["M_GAS"][ii], head)
                        merger_collection_GAS["M2_MERGERS"][end]  += convert_units_physical_mass(halo_story["M_GAS_2"][ii], head2)
                        merger_collection_GAS["M2_MISSED"][end]   += convert_units_physical_mass(halo_story["M_GAS_2"][ii], head2)
                    end
                    pos_in_mim += 1
                end
            end
        end
        # Replace zeros with missings
        replace!(merger_collection_GAS["M_MM"]         , 0. => missing)
        replace!(merger_collection_GAS["M2_MM"]        , 0. => missing)
        replace!(merger_collection_GAS["J_MMorbital"]  , 0. => missing)
        replace!(merger_collection_GAS["J_SUMorbital"] , 0. => missing)
        replace!(merger_collection_GAS["M_MERGERS"]    , 0. => missing)
        replace!(merger_collection_GAS["M_MISSED"]     , 0. => missing)
        replace!(merger_collection_GAS["M_CONSIDERED"] , 0. => missing)
        replace!(merger_collection_GAS["M2_MERGERS"]   , 0. => missing)
        replace!(merger_collection_GAS["M2_MISSED"]    , 0. => missing)
        replace!(merger_collection_GAS["M2_CONSIDERED"], 0. => missing)
        #replace!(merger_collection_GAS["Merger_Map"][1:4,:], 0. => missing)



        save(joinpath(outdir, "halo_$(halo_story["rootID"])_$(halo_story["subID"][1]).jld"), 
            "halo_story",   halo_story,
            "merger_collection_DM",     merger_collection_DM,
            "merger_collection_GAS",    merger_collection_GAS,
            "merger_collection_STARS",  merger_collection_STARS)
    end



    println("\n\n\n---------------------------\n\nNow witness the firepower of this fully armed and operational battle station!")
    return nothing
end

print("'felix2jld'   ")
