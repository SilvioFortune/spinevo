
@doc """
DESCRIPTION:\n
INPUT:\n
OUTPUT:\n
""" ->
function collect_mergers(halo_story; min_time=0.0
    )

    id_m   = halo_story["M_STARS"][1]
    id_m2       = halo_story["M_STARS_2"][1]
    id_mpeak       = halo_story["M_STARS_peak"][1]
    
    
    # Identify First Progenitors and mergers

    # Make lookbacktime array
    #lbt = halo_story["lookbacktime"]

    # First Progenitors, backward in time
    fp_indices      = Array{Int64}(undef, 0)
    fp_snaps        = Array{Int64}(undef, 0)
    loop_length     = length(halo_story["snapFP"])
    for i in 1:loop_length
        #if i==1 || halo_story["snapFP"][i-1] > halo_story["snapFP"][i] && count(ismissing, halo_story["J_STARS"][:,i]) == 0 # First Progenitor
        if halo_story["mmp"][i] == 1 && count(ismissing, halo_story["J_STARS"][:,i]) == 0 # First Progenitor
            #print("FP@SNAP$(halo_story["snapFP"][i]) ")
            #flush(stdout)
            #if length(fp_indices) == 0# && lbt[i]-lbt[1] > min_time
            #    fp_indices  = vcat( fp_indices, i )
            #    fp_snaps    = vcat( fp_snaps, halo_story["snapFP"][i] )
            #elseif length(fp_indices) > 0# && lbt[i]-lbt[fp_indices[end]] > min_time
            #    fp_indices  = vcat( fp_indices, i )
            #    fp_snaps    = vcat( fp_snaps, halo_story["snapFP"][i] )
            #end
            fp_indices  = vcat( fp_indices, i )
            fp_snaps    = vcat( fp_snaps, halo_story["snapFP"][i] )
        end
    end

    if length(fp_indices) ≤ 1 # skip if there is only one entry in the FP list
        return merger_collection_STARS, merger_collection_GAS, merger_collection_DM
    else
        # Reverse to forward in time
        fp_indices      = fp_indices[end:-1:1]
        fp_snaps        = fp_snaps[end:-1:1]
        #println(fp_indices)
        #println(fp_snaps)
        
        merger_count        = 0
        switch_count        = 0
        n_mergers           = Array{Int64}(undef, 0)
        switches            = Array{Int64}(undef, 0)
        merger_index_map    = Array{Int64}(undef, 2, 0)    # Merger indices, same-snap FP indices
        fp_index    = 2
        n_mergers   = vcat( n_mergers, merger_count )
        # check if first snap already contains a first progenitor
        #println(fp_indices)
        #flush(stdout)
        #if fp_snaps[1] == halo_story["snapFP"][end]
            #fp_index    = 2
            #n_mergers   = vcat( n_mergers, merger_count )
        #else
            #fp_index    = 1
        #end
    
        # Merger Map, forward in time by starting at the bottom
        for i in loop_length:-1:1
            #print("$i ")
            #println("$i     $(fp_snaps[end]) > $(halo_story["snapFP"][i]) > $(fp_snaps[1])")
            if fp_snaps[end] > halo_story["snapFP"][i] ≥ fp_snaps[1]  # only accept entries within snaps with available stellar spins
                #println("   $(halo_story["snapFP"][i]) == $(fp_snaps[fp_index])")
                if halo_story["snapFP"][i] == fp_snaps[fp_index] # First Progenitor snap
                    # Check out Merger Info and assign to next
                    n_mergers       = vcat( n_mergers, merger_count )
                    merger_count    = 0
                    switches        = vcat( switches, switch_count )
                    switch_count    = 0
                    fp_index       += 1
                end
            
                # Identify Mergers
                #println("       $(halo_story["mmp"][i])     $(halo_story["snapFP"][i]) == $(halo_story["snapFP"][i-1])")
                if halo_story["mmp"][i] == 1 # FP
                    switch_count   += ifelse(halo_story["switch"][i] > 0, 1, 0)
                elseif halo_story["mmp"][i] == 0 # secondary Progenitor
                    merger_count     += 1
                    merger_index_map  = hcat( merger_index_map, [i, ssFPfinder(i, halo_story)] )
                else
                    println("Error for i = $i, $(halo_story["snapFP"][i]), $(fp_snaps[fp_index]), $(halo_story["snapFP"][i-1]), $(halo_story["subID"][i-1])")
                end
            elseif fp_snaps[end] == halo_story["snapFP"][i] # final checkout
                n_mergers       = vcat( n_mergers, merger_count )
                break
            end
        end
        
        #println("$(length(fp_indices)) $(length(fp_snaps)) $(length(n_mergers))")
        #println("$(length(merger_index_map)) $(sum(n_mergers)) ")

        merger_collection_STARS = Dict(
                "N_MERGERS"     => n_mergers,
                "snapNR"        => halo_story["snapNR"][fp_indices],
                "ID_ISUB"       => Int.(ones(length(fp_indices)) .* halo_story["subID"][1]),
                "ID_M"     => ones(length(fp_indices)) .* id_m,
                "ID_M2"         => ones(length(fp_indices)) .* id_m2,
                "ID_Mpeak"         => ones(length(fp_indices)) .* id_mpeak,
                "subID"         => halo_story["subID"][fp_indices],
                #"switch"        => halo_story["switch"][fp_indices],
                "redshift"      => halo_story["redshift"][fp_indices],
                "BVAL"          => halo_story["BVAL"][fp_indices], 
                "BVAL_0"        => halo_story["BVAL_0"][fp_indices], 
                "J_main"        => halo_story["J_STARS"][:,fp_indices], 
                "j_main"        => halo_story["j_STARS"][:,fp_indices], 
                "J_vir"        => halo_story["J_STARSvir"][:,fp_indices], 
                "j_vir"        => halo_story["j_STARSvir"][:,fp_indices], 
                "J_main_0"        => halo_story["J_STARS_0"][:,fp_indices], 
                "j_main_0"        => halo_story["j_STARS_0"][:,fp_indices], 
                "J_vir_0"        => halo_story["J_STARSvir_0"][:,fp_indices], 
                "j_vir_0"        => halo_story["j_STARSvir_0"][:,fp_indices], 

                "lookbacktime"  => missings(Float64 , 0),
                "ΔM"      => missings(Float64 , 0), 
                "ΔM2"     => missings(Float64 , 0), 
                "ΔMpeak"     => missings(Float64 , 0), 
                "ΔM_fromJ"      => missings(Float64 , 0), 
                "M"       => missings(Float64 , 0), 
                "M2"      => missings(Float64 , 0), 
                "Mpeak"      => missings(Float64 , 0), 
                "M_fromJ"       => missings(Float64 , 0), 
                "ϕ_flip"        => missings(Float64 , 0), 
                "M_MERGERS"     => missings(Float64 , 0),
                "M_MISSED"      => missings(Float64 , 0), 
                "M_CONSIDERED"  => missings(Float64 , 0),  
                "M2_MERGERS"    => missings(Float64 , 0),
                "M2_MISSED"     => missings(Float64 , 0), 
                "M2_CONSIDERED" => missings(Float64 , 0),   
                "Mpeak_MERGERS"    => missings(Float64 , 0),
                "Mpeak_MISSED"     => missings(Float64 , 0), 
                "Mpeak_CONSIDERED" => missings(Float64 , 0),   
                "ΔBVAL"         => missings(Float64 , 0), 
                "ΔBVAL_0"       => missings(Float64 , 0), 

                "switch" => missings(Int64 , 0),
                "Δlookbacktime" => missings(Float64 , 0),

                
                "SFR"           => halo_story["SFR"][fp_indices], 
                #"sSFR"          => halo_story["sSFR"][fp_indices], 
                "M200"          => halo_story["M200"][fp_indices],
                "MCRI"          => halo_story["MCRI"][fp_indices],
                "M500"          => halo_story["M500"][fp_indices],
                "MVIR"          => halo_story["MVIR"][fp_indices],
                "SAGE"          => halo_story["SAGE"][fp_indices], 
                "DSUB"          => halo_story["DSUB"][fp_indices], 
                "VMAX"          => halo_story["VMAX"][fp_indices], 
                "RMAX"          => halo_story["RMAX"][fp_indices], 
                "RVIR"          => halo_story["RVIR"][fp_indices], 
                "SPOS"          => halo_story["SPOS"][:,fp_indices], 
                "SVEL"          => halo_story["SVEL"][:,fp_indices], 
                "SPIN"          => halo_story["SPIN"][:,fp_indices], 

                "M_MM"          => missings(Float64, 0), # WTF 1 = mass, 2 = FPmass
                "M2_MM"         => missings(Float64, 0), # WTF 1 = mass, 2 = FPmass
                "Mpeak_MM"         => missings(Float64, 0), # WTF 1 = mass, 2 = FPmass
                "J_MM"   => missings(Float64 , 12, 0), 
                "J_SUM"  => missings(Float64 , 12, 0), 
                "ΔJ_main"       => missings(Float64 , 3, 0), 
                "Δj_main"       => missings(Float64 , 3, 0),
                "ΔJ_vir"       => missings(Float64 , 3, 0), 
                "Δj_vir"       => missings(Float64 , 3, 0),
                "J_MM_0"   => missings(Float64 , 12, 0), 
                "J_SUM_0"  => missings(Float64 , 12, 0), 
                "ΔJ_main_0"       => missings(Float64 , 3, 0), 
                "Δj_main_0"       => missings(Float64 , 3, 0),
                "ΔJ_vir_0"       => missings(Float64 , 3, 0), 
                "Δj_vir_0"       => missings(Float64 , 3, 0),
                "merger map 0"    => missings(Float64 , 29, 0), 
                "merger map"    => missings(Float64 , 29, 0)) # 1=mass, 2=mass2, 3=FPmass, 4=FPmass2, 5=FP SNAP, 6=FP z, 7= FP lbt, 8=subID, 9=mass_peak, 10=FPmass_peak, 11=snapNR, 12-29=merger spin map
                
    
        merger_collection_DM = Dict(
                "N_MERGERS"     => n_mergers,
                "snapNR"          => halo_story["snapNR"][fp_indices],
                "ID_ISUB"       => Int.(ones(length(fp_indices)) .* halo_story["subID"][1]),
                "ID_M"     => ones(length(fp_indices)) .* id_m,
                "ID_M2"         => ones(length(fp_indices)) .* id_m2,
                "ID_Mpeak"         => ones(length(fp_indices)) .* id_mpeak,
                "subID"         => halo_story["subID"][fp_indices],
                #"switch"        => halo_story["switch"][fp_indices],
                "redshift"      => halo_story["redshift"][fp_indices],
                "J_main"        => halo_story["J_DM"][:,fp_indices], 
                "j_main"        => halo_story["j_DM"][:,fp_indices], 

                "J_main_0"        => halo_story["J_DM_0"][:,fp_indices], 
                "j_main_0"        => halo_story["j_DM_0"][:,fp_indices], 

                "lookbacktime"  => missings(Float64 , 0),
                "ΔM"      => missings(Float64 , 0), 
                "ΔM2"     => missings(Float64 , 0), 
                "ΔMpeak"     => missings(Float64 , 0), 
                "ΔM_fromJ"      => missings(Float64 , 0), 
                "M"       => missings(Float64 , 0), 
                "M2"      => missings(Float64 , 0), 
                "Mpeak"      => missings(Float64 , 0), 
                "M_fromJ"       => missings(Float64 , 0), 
                "ϕ_flip"        => missings(Float64 , 0), 

                "switch" => missings(Int64 , 0),
                "Δlookbacktime" => missings(Float64 , 0),

                "M_MERGERS"     => missings(Float64 , 0),
                "M_MISSED"      => missings(Float64 , 0), 
                "M_CONSIDERED"  => missings(Float64 , 0),  
                "M2_MERGERS"    => missings(Float64 , 0),
                "M2_MISSED"     => missings(Float64 , 0), 
                "M2_CONSIDERED" => missings(Float64 , 0),  
                "Mpeak_MERGERS"    => missings(Float64 , 0),
                "Mpeak_MISSED"     => missings(Float64 , 0), 
                "Mpeak_CONSIDERED" => missings(Float64 , 0),  
                "M_MM"          => missings(Float64, 0), # WTF 1 = mass, 2 = FPmass
                "M2_MM"         => missings(Float64, 0), # WTF 1 = mass, 2 = FPmass
                "Mpeak_MM"         => missings(Float64, 0), # WTF 1 = mass, 2 = FPmass
                "J_MM"   => missings(Float64 , 12, 0), 
                "J_SUM"  => missings(Float64 , 12, 0), 
                "ΔJ_main"       => missings(Float64 , 3, 0), 
                "Δj_main"       => missings(Float64 , 3, 0),
                "merger map"    => missings(Float64 , 29, 0), # 1=mass, 2=mass2, 3=FPmass, 4=FPmass2, 5=FP SNAP, 6=FP z, 7=FP lbt, 8=subID, 9=mass_peak, 10=FPmass_peak, 11=snapNR
                
                "J_MM_0"   => missings(Float64 , 12, 0), 
                "J_SUM_0"  => missings(Float64 , 12, 0), 
                "ΔJ_main_0"       => missings(Float64 , 3, 0), 
                "Δj_main_0"       => missings(Float64 , 3, 0),
                "merger map 0"    => missings(Float64 , 29, 0))
    
        merger_collection_GAS = Dict(
                "N_MERGERS"     => n_mergers,
                "snapNR"          => halo_story["snapNR"][fp_indices],
                "ID_ISUB"       => Int.(ones(length(fp_indices)) .* halo_story["subID"][1]),
                "ID_M"     => ones(length(fp_indices)) .* id_m,
                "ID_M2"         => ones(length(fp_indices)) .* id_m2,
                "ID_Mpeak"         => ones(length(fp_indices)) .* id_mpeak,
                "subID"         => halo_story["subID"][fp_indices],
                #"switch"        => halo_story["switch"][fp_indices],
                "redshift"      => halo_story["redshift"][fp_indices],
                "J_main"        => halo_story["J_GAS"][:,fp_indices], 
                "j_main"        => halo_story["j_GAS"][:,fp_indices], 
                "J_vir"        => halo_story["J_GASvir"][:,fp_indices], 
                "j_vir"        => halo_story["j_GASvir"][:,fp_indices], 

                "J_main_0"        => halo_story["J_GAS_0"][:,fp_indices], 
                "j_main_0"        => halo_story["j_GAS_0"][:,fp_indices], 
                "J_vir_0"        => halo_story["J_GASvir_0"][:,fp_indices], 
                "j_vir_0"        => halo_story["j_GASvir_0"][:,fp_indices], 

                "lookbacktime"  => missings(Float64 , 0),
                "ΔM"      => missings(Float64 , 0), 
                "ΔM2"     => missings(Float64 , 0), 
                "ΔMpeak"     => missings(Float64 , 0), 
                "ΔM_fromJ"      => missings(Float64 , 0), 
                "M"       => missings(Float64 , 0), 
                "M2"      => missings(Float64 , 0), 
                "Mpeak"      => missings(Float64 , 0), 
                "M_fromJ"       => missings(Float64 , 0), 
                "ϕ_flip"        => missings(Float64 , 0), 

                "switch" => missings(Int64 , 0),
                "Δlookbacktime" => missings(Float64 , 0),

                "M_MERGERS"     => missings(Float64 , 0),
                "M_MISSED"      => missings(Float64 , 0), 
                "M_CONSIDERED"  => missings(Float64 , 0),  
                "M2_MERGERS"    => missings(Float64 , 0),
                "M2_MISSED"     => missings(Float64 , 0), 
                "M2_CONSIDERED" => missings(Float64 , 0),  
                "Mpeak_MERGERS"    => missings(Float64 , 0),
                "Mpeak_MISSED"     => missings(Float64 , 0), 
                "Mpeak_CONSIDERED" => missings(Float64 , 0),  
                #"ΔBVAL"         => missings(Float64 , 0), 
                #"ΔBVAL_0"         => missings(Float64 , 0), 
                "M_MM"          => missings(Float64, 0), # WTF 1 = mass, 2 = FPmass
                "M2_MM"         => missings(Float64, 0), # WTF 1 = mass, 2 = FPmass
                "Mpeak_MM"         => missings(Float64, 0), # WTF 1 = mass, 2 = FPmass
                "J_MM"   => missings(Float64 , 12, 0), 
                "J_SUM"  => missings(Float64 , 12, 0), 
                "ΔJ_main"       => missings(Float64 , 3, 0), 
                "Δj_main"       => missings(Float64 , 3, 0),
                "ΔJ_vir"       => missings(Float64 , 3, 0), 
                "Δj_vir"       => missings(Float64 , 3, 0),
                "merger map"    => missings(Float64 , 29, 0), 

                "J_MM_0"   => missings(Float64 , 12, 0), 
                "J_SUM_0"  => missings(Float64 , 12, 0), 
                "ΔJ_main_0"       => missings(Float64 , 3, 0), 
                "Δj_main_0"       => missings(Float64 , 3, 0),
                "ΔJ_vir_0"       => missings(Float64 , 3, 0), 
                "Δj_vir_0"       => missings(Float64 , 3, 0),
                "merger map 0"    => missings(Float64 , 29, 0)) # 1=mass, 2=mass2, 3=FPmass, 4=FPmass2, 5=FP SNAP, 6=FP z, 7=FP lbt, 8=subID, 9=mass_peak, 10=FPmass_peak, 11=snapNR
        
        
        # Fill the STARS dictionary
    
        for i in 1:length(fp_indices)
            #print("$i ")
            #flush(stdout)/home/moon/sfortune/spinevo/halostories_v20211204_min0.0Gyr/halo_1_1414.jld
    
            # Basic Info
            merger_collection_STARS["lookbacktime"] = vcat( merger_collection_STARS["lookbacktime"], halo_story["lookbacktime"][fp_indices[i]] )
            merger_collection_STARS["M"]      = vcat( merger_collection_STARS["M"], halo_story["M_STARS"][fp_indices[i]] )
            merger_collection_STARS["M2"]     = vcat( merger_collection_STARS["M2"], halo_story["M_STARS_2"][fp_indices[i]] )
            merger_collection_STARS["Mpeak"]     = vcat( merger_collection_STARS["Mpeak"], halo_story["M_STARS_peak"][fp_indices[i]] )
            #merger_collection_STARS["SFR"]          = vcat( merger_collection_STARS["SFR"], halo_story["SFR"][fp_indices[i]] )
            
            if count(ismissing, halo_story["j_STARS"][:,fp_indices[i]]) == 0 
                merger_collection_STARS["M_fromJ"]  = vcat( merger_collection_STARS["M_fromJ"], norm(halo_story["J_STARS"][:,fp_indices[i]]) / norm(halo_story["j_STARS"][:,fp_indices[i]]) )
            else
                merger_collection_STARS["M_fromJ"]  = vcat( merger_collection_STARS["M_fromJ"], missing )
            end
            
            # Transitional Data
            if i == 1
                merger_collection_STARS["ΔJ_main"]          = hcat( merger_collection_STARS["ΔJ_main"], missings(Float64, 3) )
                merger_collection_STARS["Δj_main"]          = hcat( merger_collection_STARS["Δj_main"], missings(Float64, 3) )
                merger_collection_STARS["ΔJ_main_0"]          = hcat( merger_collection_STARS["ΔJ_main_0"], missings(Float64, 3) )
                merger_collection_STARS["Δj_main_0"]          = hcat( merger_collection_STARS["Δj_main_0"], missings(Float64, 3) )
                merger_collection_STARS["ΔJ_vir"]          = hcat( merger_collection_STARS["ΔJ_vir"], missings(Float64, 3) )
                merger_collection_STARS["Δj_vir"]          = hcat( merger_collection_STARS["Δj_vir"], missings(Float64, 3) )
                merger_collection_STARS["ΔJ_vir_0"]          = hcat( merger_collection_STARS["ΔJ_vir_0"], missings(Float64, 3) )
                merger_collection_STARS["Δj_vir_0"]          = hcat( merger_collection_STARS["Δj_vir_0"], missings(Float64, 3) )
                merger_collection_STARS["ΔBVAL"]            = vcat( merger_collection_STARS["ΔBVAL"], missing )
                merger_collection_STARS["ΔBVAL_0"]            = vcat( merger_collection_STARS["ΔBVAL_0"], missing )
                merger_collection_STARS["ϕ_flip"]           = vcat( merger_collection_STARS["ϕ_flip"], missing )
                merger_collection_STARS["Δlookbacktime"]           = vcat( merger_collection_STARS["Δlookbacktime"], missing )
                merger_collection_STARS["switch"]           = vcat( merger_collection_STARS["switch"], missing )
                merger_collection_STARS["ΔM"]         = vcat( merger_collection_STARS["ΔM"], missing )
                merger_collection_STARS["ΔM2"]        = vcat( merger_collection_STARS["ΔM2"], missing )
                merger_collection_STARS["ΔMpeak"]        = vcat( merger_collection_STARS["ΔMpeak"], missing )
                merger_collection_STARS["ΔM_fromJ"]         = vcat( merger_collection_STARS["ΔM_fromJ"], missing )
                merger_collection_STARS["M_MM"]             = vcat( merger_collection_STARS["M_MM"], missing )
                merger_collection_STARS["M2_MM"]            = vcat( merger_collection_STARS["M2_MM"], missing )
                merger_collection_STARS["Mpeak_MM"]            = vcat( merger_collection_STARS["Mpeak_MM"], missing )
                merger_collection_STARS["J_MM"]      = hcat( merger_collection_STARS["J_MM"], missings(Float64, 12) )
                merger_collection_STARS["J_MM_0"]      = hcat( merger_collection_STARS["J_MM_0"], missings(Float64, 12) )
                merger_collection_STARS["J_SUM"]     = hcat( merger_collection_STARS["J_SUM"], zeros(12) )
                merger_collection_STARS["J_SUM_0"]     = hcat( merger_collection_STARS["J_SUM_0"], zeros(12) )
                merger_collection_STARS["M_MERGERS"]        = vcat( merger_collection_STARS["M_MERGERS"], missing )
                merger_collection_STARS["M_MISSED"]         = vcat( merger_collection_STARS["M_MISSED"], missing )
                merger_collection_STARS["M_CONSIDERED"]     = vcat( merger_collection_STARS["M_CONSIDERED"], missing )
                merger_collection_STARS["M2_MERGERS"]       = vcat( merger_collection_STARS["M2_MERGERS"], missing )
                merger_collection_STARS["M2_MISSED"]        = vcat( merger_collection_STARS["M2_MISSED"], missing )
                merger_collection_STARS["M2_CONSIDERED"]    = vcat( merger_collection_STARS["M2_CONSIDERED"], missing )
                merger_collection_STARS["Mpeak_MERGERS"]       = vcat( merger_collection_STARS["Mpeak_MERGERS"], missing )
                merger_collection_STARS["Mpeak_MISSED"]        = vcat( merger_collection_STARS["Mpeak_MISSED"], missing )
                merger_collection_STARS["Mpeak_CONSIDERED"]    = vcat( merger_collection_STARS["Mpeak_CONSIDERED"], missing )
            else 
                mintidx = find_mintime_idx( halo_story["lookbacktime"], fp_indices, i, min_time )
                merger_collection_STARS["ΔJ_main"]  = hcat( merger_collection_STARS["ΔJ_main"], halo_story["J_STARS"][:,fp_indices[i]] .- halo_story["J_STARS"][:,mintidx] )
                merger_collection_STARS["Δj_main"]  = hcat( merger_collection_STARS["Δj_main"], halo_story["j_STARS"][:,fp_indices[i]] .- halo_story["j_STARS"][:,mintidx] )
                merger_collection_STARS["ΔJ_main_0"]  = hcat( merger_collection_STARS["ΔJ_main_0"], halo_story["J_STARS_0"][:,fp_indices[i]] .- halo_story["J_STARS_0"][:,mintidx] )
                merger_collection_STARS["Δj_main_0"]  = hcat( merger_collection_STARS["Δj_main_0"], halo_story["j_STARS_0"][:,fp_indices[i]] .- halo_story["j_STARS_0"][:,mintidx] )
                merger_collection_STARS["ΔJ_vir"]  = hcat( merger_collection_STARS["ΔJ_vir"], halo_story["J_STARSvir"][:,fp_indices[i]] .- halo_story["J_STARSvir"][:,mintidx] )
                merger_collection_STARS["Δj_vir"]  = hcat( merger_collection_STARS["Δj_vir"], halo_story["j_STARSvir"][:,fp_indices[i]] .- halo_story["j_STARSvir"][:,mintidx] )
                merger_collection_STARS["ΔJ_vir_0"]  = hcat( merger_collection_STARS["ΔJ_vir_0"], halo_story["J_STARSvir_0"][:,fp_indices[i]] .- halo_story["J_STARSvir_0"][:,mintidx] )
                merger_collection_STARS["Δj_vir_0"]  = hcat( merger_collection_STARS["Δj_vir_0"], halo_story["j_STARSvir_0"][:,fp_indices[i]] .- halo_story["j_STARSvir_0"][:,mintidx] )
                merger_collection_STARS["ΔBVAL"]    = vcat( merger_collection_STARS["ΔBVAL"], halo_story["BVAL"][fp_indices[i]] - halo_story["BVAL"][mintidx] )
                merger_collection_STARS["ΔBVAL_0"]    = vcat( merger_collection_STARS["ΔBVAL_0"], halo_story["BVAL_0"][fp_indices[i]] - halo_story["BVAL_0"][mintidx] )
                merger_collection_STARS["Δlookbacktime"]           = vcat( merger_collection_STARS["Δlookbacktime"], halo_story["lookbacktime"][fp_indices[i]] - halo_story["lookbacktime"][mintidx] )
                merger_collection_STARS["switch"]           = vcat( merger_collection_STARS["switch"], sum( halo_story["switch"][fp_indices[findcs( fp_indices, geq=fp_indices[i], lt=mintidx)]] ) )
                if count(ismissing, halo_story["j_STARS"][:,fp_indices[i]]) == 0 && count(ismissing, halo_story["j_STARS"][:,mintidx]) == 0 
                    merger_collection_STARS["ϕ_flip"]   = vcat( merger_collection_STARS["ϕ_flip"], 180/π * angle(halo_story["J_STARS"][:,fp_indices[i]], halo_story["J_STARS"][:,mintidx] ) )
                    merger_collection_STARS["ΔM_fromJ"] = vcat( merger_collection_STARS["ΔM_fromJ"], (norm(halo_story["J_STARS"][:,fp_indices[i]]) / norm(halo_story["j_STARS"][:,fp_indices[i]])) - (norm(halo_story["J_STARS"][:,mintidx]) / norm(halo_story["j_STARS"][:,mintidx])) )
                else
                    #println(halo_story["j_STARS"][:,fp_indices[i]])
                    #println("missing j_STARS for $(storyfilelist[iii]) @ index $(fp_indices[i])")
                    merger_collection_STARS["ϕ_flip"]   = vcat( merger_collection_STARS["ϕ_flip"], missing )
                    merger_collection_STARS["ΔM_fromJ"] = vcat( merger_collection_STARS["ΔM_fromJ"], missing )
                end
                #println("\n$(halo_story["J_STARS"][:,fp_indices[i]]) $(halo_story["J_STARS"][:,mintidx])")
                merger_collection_STARS["ΔM"] = vcat( merger_collection_STARS["ΔM"], halo_story["M_STARS"][fp_indices[i]] - halo_story["M_STARS"][mintidx] )
                merger_collection_STARS["ΔM2"] = vcat( merger_collection_STARS["ΔM2"], halo_story["M_STARS_2"][fp_indices[i]] - halo_story["M_STARS_2"][mintidx] )
                merger_collection_STARS["ΔMpeak"] = vcat( merger_collection_STARS["ΔMpeak"], halo_story["M_STARS_peak"][fp_indices[i]] - halo_story["M_STARS_peak"][mintidx] )
                
                # Merger Data
                # Setup
                merger_collection_STARS["M_MM"]             = vcat( merger_collection_STARS["M_MM"], 0 )
                merger_collection_STARS["M2_MM"]            = vcat( merger_collection_STARS["M2_MM"], 0 )
                merger_collection_STARS["Mpeak_MM"]            = vcat( merger_collection_STARS["Mpeak_MM"], 0 )
                merger_collection_STARS["J_MM"]      = hcat( merger_collection_STARS["J_MM"], missings(Float64, 12) )
                merger_collection_STARS["J_MM_0"]      = hcat( merger_collection_STARS["J_MM_0"], missings(Float64, 12) )
                merger_collection_STARS["J_SUM"]     = hcat( merger_collection_STARS["J_SUM"], zeros(12) )
                merger_collection_STARS["J_SUM_0"]     = hcat( merger_collection_STARS["J_SUM_0"], zeros(12) )
                merger_collection_STARS["M_MERGERS"]        = vcat( merger_collection_STARS["M_MERGERS"], 0 )
                merger_collection_STARS["M_MISSED"]         = vcat( merger_collection_STARS["M_MISSED"], 0 )
                merger_collection_STARS["M_CONSIDERED"]     = vcat( merger_collection_STARS["M_CONSIDERED"], 0 )
                merger_collection_STARS["M2_MERGERS"]       = vcat( merger_collection_STARS["M2_MERGERS"], 0 )
                merger_collection_STARS["M2_MISSED"]        = vcat( merger_collection_STARS["M2_MISSED"], 0 )
                merger_collection_STARS["M2_CONSIDERED"]    = vcat( merger_collection_STARS["M2_CONSIDERED"], 0 )
                merger_collection_STARS["Mpeak_MERGERS"]       = vcat( merger_collection_STARS["Mpeak_MERGERS"], 0 )
                merger_collection_STARS["Mpeak_MISSED"]        = vcat( merger_collection_STARS["Mpeak_MISSED"], 0 )
                merger_collection_STARS["Mpeak_CONSIDERED"]    = vcat( merger_collection_STARS["Mpeak_CONSIDERED"], 0 )
                # Actual check
                pos_in_mim  = 1+sum(n_mergers[1:i-1]) # position / index in merger index map
                for ii in merger_index_map[1,1+sum(n_mergers[1:i-1]):sum(n_mergers[1:i])]
                    # Merger Map
                    #println("$(halo_story["snapNR"][ii])   i = $i")
                    #@show merger_index_map
                    #@show merger_index_map[1,1+sum(n_mergers[1:i-1]):sum(n_mergers[1:i])]
                    #flush(stdout)
                    merger_collection_STARS["merger map"]   = hcat( merger_collection_STARS["merger map"], 
                                                                        vcat( [   halo_story["M_STARS"][ii],
                                                                            halo_story["M_STARS_2"][ii], 
                                                                            halo_story["M_STARS"][merger_index_map[2,pos_in_mim]], 
                                                                            halo_story["M_STARS_2"][merger_index_map[2,pos_in_mim]], 
                                                                            halo_story["snapFP"][ii], 
                                                                            halo_story["redshift"][fp_indices[i]], 
                                                                            halo_story["lookbacktime"][fp_indices[i]] , 
                                                                            halo_story["subID"][ii],
                                                                            halo_story["M_STARS_peak"][ii], 
                                                                            halo_story["M_STARS_peak"][merger_index_map[2,pos_in_mim]], 
                                                                            halo_story["snapNR"][ii]
                                                                            ], halo_story["merger_spin_map_STARS"][:,ii] )
                                                                            )
                    merger_collection_STARS["merger map 0"]   = hcat( merger_collection_STARS["merger map 0"], 
                                                                        vcat( [   halo_story["M_STARS"][ii],
                                                                            halo_story["M_STARS_2"][ii], 
                                                                            halo_story["M_STARS"][merger_index_map[2,pos_in_mim]], 
                                                                            halo_story["M_STARS_2"][merger_index_map[2,pos_in_mim]], 
                                                                            halo_story["snapFP"][ii], 
                                                                            halo_story["redshift"][fp_indices[i]], 
                                                                            halo_story["lookbacktime"][fp_indices[i]] , 
                                                                            halo_story["subID"][ii],
                                                                            halo_story["M_STARS_peak"][ii], 
                                                                            halo_story["M_STARS_peak"][merger_index_map[2,pos_in_mim]], 
                                                                            halo_story["snapNR"][ii]
                                                                            ], halo_story["merger_spin_map_STARS_0"][:,ii] )
                                                                            )
                    # Most Massive condition
                    if halo_story["M_STARS_peak"][ii] > merger_collection_STARS["Mpeak_MM"][end]
                        merger_collection_STARS["M_MM"][end]            = halo_story["M_STARS"][ii]
                        merger_collection_STARS["M2_MM"][end]           = halo_story["M_STARS_2"][ii]
                        merger_collection_STARS["Mpeak_MM"][end]           = halo_story["M_STARS_peak"][ii]
                        merger_collection_STARS["J_MM"][:,end]   = vcat(halo_story["merger_spin_map_STARS"][4:6,ii] .*halo_story["merger_spin_map_STARS"][3,ii] , 
                                                                        halo_story["merger_spin_map_STARS"][7:9,ii], 
                                                                        halo_story["merger_spin_map_STARS"][13:15,ii] .*halo_story["merger_spin_map_STARS"][12,ii] , 
                                                                        halo_story["merger_spin_map_STARS"][16:18,ii]
                                                                        ) 
                        merger_collection_STARS["J_MM_0"][:,end]   = vcat(halo_story["merger_spin_map_STARS_0"][4:6,ii] .*halo_story["merger_spin_map_STARS_0"][3,ii] , 
                                                                        halo_story["merger_spin_map_STARS_0"][7:9,ii], 
                                                                        halo_story["merger_spin_map_STARS_0"][13:15,ii] .*halo_story["merger_spin_map_STARS_0"][12,ii] , 
                                                                        halo_story["merger_spin_map_STARS_0"][16:18,ii]
                                                                        ) 
                    end
                    # Orbital Data
                    merger_collection_STARS["M_MERGERS"][end]   += halo_story["M_STARS"][ii]
                    merger_collection_STARS["M2_MERGERS"][end]  += halo_story["M_STARS_2"][ii]
                    merger_collection_STARS["Mpeak_MERGERS"][end]  += halo_story["M_STARS_peak"][ii]
                    merger_collection_STARS["J_SUM"][:,end]       .+= replace(vcat(halo_story["merger_spin_map_STARS"][4:6,ii] .*halo_story["merger_spin_map_STARS"][3,ii] , 
                                                                                halo_story["merger_spin_map_STARS"][7:9,ii], 
                                                                                halo_story["merger_spin_map_STARS"][13:15,ii] .*halo_story["merger_spin_map_STARS"][12,ii] , 
                                                                                halo_story["merger_spin_map_STARS"][16:18,ii]
                                                                                ),
                                                                            missing => 0.0)
                    merger_collection_STARS["J_SUM_0"][:,end]       .+= replace(vcat(halo_story["merger_spin_map_STARS_0"][4:6,ii] .*halo_story["merger_spin_map_STARS_0"][3,ii] , 
                                                                                halo_story["merger_spin_map_STARS_0"][7:9,ii], 
                                                                                halo_story["merger_spin_map_STARS_0"][13:15,ii] .*halo_story["merger_spin_map_STARS_0"][12,ii] , 
                                                                                halo_story["merger_spin_map_STARS_0"][16:18,ii]
                                                                                ),
                                                                            missing => 0.0)
                    if count(ismissing, halo_story["merger_spin_map_STARS"][:,ii]) == 0
                        merger_collection_STARS["M_CONSIDERED"][end]    += halo_story["M_STARS"][ii]
                        merger_collection_STARS["M2_CONSIDERED"][end]   += halo_story["M_STARS_2"][ii]
                        merger_collection_STARS["Mpeak_CONSIDERED"][end]   += halo_story["M_STARS_peak"][ii]
                    else
                        #println("$ii")
                        merger_collection_STARS["M_MISSED"][end]    += halo_story["M_STARS"][ii]
                        merger_collection_STARS["M2_MISSED"][end]   += halo_story["M_STARS_2"][ii]
                        merger_collection_STARS["Mpeak_MISSED"][end]   += halo_story["M_STARS_peak"][ii]
                    end
                    pos_in_mim += 1
                end
            end
        end
        # Replace zeros with missings
        replace!(merger_collection_STARS["M_MM"]         , 0. => missing)
        replace!(merger_collection_STARS["M2_MM"]        , 0. => missing)
        replace!(merger_collection_STARS["Mpeak_MM"]        , 0. => missing)
        replace!(merger_collection_STARS["J_MM"]  , 0. => missing)
        replace!(merger_collection_STARS["J_MM_0"]  , 0. => missing)
        #replace!(merger_collection_STARS["J_SUM"] , 0. => missing)
        #replace!(merger_collection_STARS["M_MERGERS"]    , 0. => missing)
        #replace!(merger_collection_STARS["M_MISSED"]     , 0. => missing)
        #replace!(merger_collection_STARS["M_CONSIDERED"] , 0. => missing)
        #replace!(merger_collection_STARS["M2_MERGERS"]   , 0. => missing)
        #replace!(merger_collection_STARS["M2_MISSED"]    , 0. => missing)
        #replace!(merger_collection_STARS["M2_CONSIDERED"], 0. => missing)
        #replace!(merger_collection_STARS["Mpeak_MERGERS"]   , 0. => missing)
        #replace!(merger_collection_STARS["Mpeak_MISSED"]    , 0. => missing)
        #replace!(merger_collection_STARS["Mpeak_CONSIDERED"], 0. => missing)
        #replace!(merger_collection_STARS["merger map"][1:4,:], 0. => missing)
    
    
        # Fill the DM dictionary
    
        for i in 1:length(fp_indices)
            #print("$i ")
            #flush(stdout)/home/moon/sfortune/spinevo/halostories_v20211204_min0.0Gyr/halo_1_1414.jld
    
            # Basic Info
            merger_collection_DM["lookbacktime"] = vcat( merger_collection_DM["lookbacktime"], halo_story["lookbacktime"][fp_indices[i]] )
            merger_collection_DM["M"]      = vcat( merger_collection_DM["M"],  halo_story["M_DM"][fp_indices[i]] )
            merger_collection_DM["M2"]     = vcat( merger_collection_DM["M2"], halo_story["M_DM_2"][fp_indices[i]] )
            merger_collection_DM["Mpeak"]     = vcat( merger_collection_DM["Mpeak"], halo_story["M_DM_peak"][fp_indices[i]] )
            if count(ismissing, halo_story["j_DM"][:,fp_indices[i]]) == 0 
                merger_collection_DM["M_fromJ"]  = vcat( merger_collection_DM["M_fromJ"], norm(halo_story["J_DM"][:,fp_indices[i]]) / norm(halo_story["j_DM"][:,fp_indices[i]]) )
            else
                merger_collection_DM["M_fromJ"]  = vcat( merger_collection_DM["M_fromJ"], missing )
            end
            
            # Transitional Data
            if i == 1
                merger_collection_DM["ΔJ_main"]          = hcat( merger_collection_DM["ΔJ_main"], missings(Float64, 3) )
                merger_collection_DM["Δj_main"]          = hcat( merger_collection_DM["Δj_main"], missings(Float64, 3) )
                merger_collection_DM["ΔJ_main_0"]          = hcat( merger_collection_DM["ΔJ_main_0"], missings(Float64, 3) )
                merger_collection_DM["Δj_main_0"]          = hcat( merger_collection_DM["Δj_main_0"], missings(Float64, 3) )
                merger_collection_DM["ϕ_flip"]           = vcat( merger_collection_DM["ϕ_flip"], missing )
                merger_collection_DM["Δlookbacktime"]    = vcat( merger_collection_DM["Δlookbacktime"], missing )
                merger_collection_DM["switch"]           = vcat( merger_collection_DM["switch"], missing )
                merger_collection_DM["ΔM"]         = vcat( merger_collection_DM["ΔM"], missing )
                merger_collection_DM["ΔM2"]        = vcat( merger_collection_DM["ΔM2"], missing )
                merger_collection_DM["ΔMpeak"]        = vcat( merger_collection_DM["ΔMpeak"], missing )
                merger_collection_DM["ΔM_fromJ"]         = vcat( merger_collection_DM["ΔM_fromJ"], missing )
                merger_collection_DM["M_MM"]             = vcat( merger_collection_DM["M_MM"], missing )
                merger_collection_DM["M2_MM"]            = vcat( merger_collection_DM["M2_MM"], missing )
                merger_collection_DM["Mpeak_MM"]            = vcat( merger_collection_DM["Mpeak_MM"], missing )
                merger_collection_DM["J_MM"]      = hcat( merger_collection_DM["J_MM"], missings(Float64, 12) )
                merger_collection_DM["J_MM_0"]      = hcat( merger_collection_DM["J_MM_0"], missings(Float64, 12) )
                merger_collection_DM["J_SUM"]     = hcat( merger_collection_DM["J_SUM"], zeros(12) )
                merger_collection_DM["J_SUM_0"]     = hcat( merger_collection_DM["J_SUM_0"], zeros(12) )
                merger_collection_DM["M_MERGERS"]        = vcat( merger_collection_DM["M_MERGERS"], missing )
                merger_collection_DM["M_MISSED"]         = vcat( merger_collection_DM["M_MISSED"], missing )
                merger_collection_DM["M_CONSIDERED"]     = vcat( merger_collection_DM["M_CONSIDERED"], missing )
                merger_collection_DM["M2_MERGERS"]       = vcat( merger_collection_DM["M2_MERGERS"], missing )
                merger_collection_DM["M2_MISSED"]        = vcat( merger_collection_DM["M2_MISSED"], missing )
                merger_collection_DM["M2_CONSIDERED"]    = vcat( merger_collection_DM["M2_CONSIDERED"], missing )
                merger_collection_DM["Mpeak_MERGERS"]       = vcat( merger_collection_DM["Mpeak_MERGERS"], missing )
                merger_collection_DM["Mpeak_MISSED"]        = vcat( merger_collection_DM["Mpeak_MISSED"], missing )
                merger_collection_DM["Mpeak_CONSIDERED"]    = vcat( merger_collection_DM["Mpeak_CONSIDERED"], missing )
            else
                mintidx = find_mintime_idx( halo_story["lookbacktime"], fp_indices, i, min_time )
                merger_collection_DM["ΔJ_main_0"]  = hcat( merger_collection_DM["ΔJ_main_0"], halo_story["J_DM_0"][:,fp_indices[i]] .- halo_story["J_DM_0"][:,mintidx] )
                merger_collection_DM["Δj_main_0"]  = hcat( merger_collection_DM["Δj_main_0"], halo_story["j_DM_0"][:,fp_indices[i]] .- halo_story["j_DM_0"][:,mintidx] )
                merger_collection_DM["ΔJ_main"]  = hcat( merger_collection_DM["ΔJ_main"], halo_story["J_DM"][:,fp_indices[i]] .- halo_story["J_DM"][:,mintidx] )
                merger_collection_DM["Δj_main"]  = hcat( merger_collection_DM["Δj_main"], halo_story["j_DM"][:,fp_indices[i]] .- halo_story["j_DM"][:,mintidx] )
                merger_collection_DM["Δlookbacktime"]           = vcat( merger_collection_DM["Δlookbacktime"], halo_story["lookbacktime"][fp_indices[i]] - halo_story["lookbacktime"][mintidx] )
                merger_collection_DM["switch"]           = vcat( merger_collection_DM["switch"], sum( halo_story["switch"][fp_indices[findcs( fp_indices, geq=fp_indices[i], lt=mintidx)]] ) )
                if count(ismissing, halo_story["j_DM"][:,fp_indices[i]]) == 0 && count(ismissing, halo_story["j_DM"][:,mintidx]) == 0 
                    merger_collection_DM["ϕ_flip"]   = vcat( merger_collection_DM["ϕ_flip"], 180/π * angle( halo_story["J_DM"][:,fp_indices[i]], halo_story["J_DM"][:,mintidx] ) )
                    merger_collection_DM["ΔM_fromJ"] = vcat( merger_collection_DM["ΔM_fromJ"], (norm(halo_story["J_DM"][:,fp_indices[i]]) / norm(halo_story["j_DM"][:,fp_indices[i]])) - (norm(halo_story["J_DM"][:,mintidx]) / norm(halo_story["j_DM"][:,mintidx])) )
                else
                    #println(halo_story["j_DM"][:,fp_indices[i]])
                    #println("missing j_DM for $(storyfilelist[iii]) @ index $(fp_indices[i])")
                    merger_collection_DM["ϕ_flip"]   = vcat( merger_collection_DM["ϕ_flip"], missing )
                    merger_collection_DM["ΔM_fromJ"] = vcat( merger_collection_DM["ΔM_fromJ"], missing )
                end
                #println("\n$(halo_story["J_DM"][:,fp_indices[i]]) $(halo_story["J_DM"][:,mintidx])")
                merger_collection_DM["ΔM"] = vcat( merger_collection_DM["ΔM"],  halo_story["M_DM"][fp_indices[i]] - halo_story["M_DM"][mintidx] )
                merger_collection_DM["ΔM2"] = vcat( merger_collection_DM["ΔM2"], halo_story["M_DM_2"][fp_indices[i]] - halo_story["M_DM_2"][mintidx])
                merger_collection_DM["ΔMpeak"] = vcat( merger_collection_DM["ΔMpeak"], halo_story["M_DM_peak"][fp_indices[i]] - halo_story["M_DM_peak"][mintidx])
                
                # Merger Data
                # Setup
                merger_collection_DM["M_MM"]             = vcat( merger_collection_DM["M_MM"], 0 )
                merger_collection_DM["M2_MM"]            = vcat( merger_collection_DM["M2_MM"], 0 )
                merger_collection_DM["Mpeak_MM"]            = vcat( merger_collection_DM["Mpeak_MM"], 0 )
                merger_collection_DM["J_MM"]      = hcat( merger_collection_DM["J_MM"], missings(Float64, 12) )
                merger_collection_DM["J_MM_0"]      = hcat( merger_collection_DM["J_MM_0"], missings(Float64, 12) )
                merger_collection_DM["J_SUM"]     = hcat( merger_collection_DM["J_SUM"], zeros(12) )
                merger_collection_DM["J_SUM_0"]     = hcat( merger_collection_DM["J_SUM_0"], zeros(12) )
                merger_collection_DM["M_MERGERS"]        = vcat( merger_collection_DM["M_MERGERS"], 0 )
                merger_collection_DM["M_MISSED"]         = vcat( merger_collection_DM["M_MISSED"], 0 )
                merger_collection_DM["M_CONSIDERED"]     = vcat( merger_collection_DM["M_CONSIDERED"], 0 )
                merger_collection_DM["M2_MERGERS"]       = vcat( merger_collection_DM["M2_MERGERS"], 0 )
                merger_collection_DM["M2_MISSED"]        = vcat( merger_collection_DM["M2_MISSED"], 0 )
                merger_collection_DM["M2_CONSIDERED"]    = vcat( merger_collection_DM["M2_CONSIDERED"], 0 )
                merger_collection_DM["Mpeak_MERGERS"]       = vcat( merger_collection_DM["Mpeak_MERGERS"], 0 )
                merger_collection_DM["Mpeak_MISSED"]        = vcat( merger_collection_DM["Mpeak_MISSED"], 0 )
                merger_collection_DM["Mpeak_CONSIDERED"]    = vcat( merger_collection_DM["Mpeak_CONSIDERED"], 0 )
                # Actual check
                pos_in_mim  = 1+sum(n_mergers[1:i-1]) # position / index in merger index map
                for ii in merger_index_map[1,1+sum(n_mergers[1:i-1]):sum(n_mergers[1:i])]
                    # Merger Map
                    #println("$(halo_story["snapNR"][ii])   i = $i")
                    #@show merger_index_map
                    #@show merger_index_map[1,1+sum(n_mergers[1:i-1]):sum(n_mergers[1:i])]
                    #flush(stdout)
                    merger_collection_DM["merger map"]   = hcat( merger_collection_DM["merger map"], 
                                                                        vcat( [   halo_story["M_DM"][ii],
                                                                            halo_story["M_DM_2"][ii],
                                                                            halo_story["M_DM"][merger_index_map[2,pos_in_mim]], 
                                                                            halo_story["M_DM_2"][merger_index_map[2,pos_in_mim]], 
                                                                            halo_story["snapFP"][ii], 
                                                                            halo_story["redshift"][fp_indices[i]], 
                                                                            halo_story["lookbacktime"][fp_indices[i]], 
                                                                            halo_story["subID"][ii],
                                                                            halo_story["M_DM_peak"][ii],
                                                                            halo_story["M_DM_peak"][merger_index_map[2,pos_in_mim]], 
                                                                            halo_story["snapNR"][ii]
                                                                            ], halo_story["merger_spin_map_DM"][:,ii] )
                                                                            )
                    merger_collection_DM["merger map 0"]   = hcat( merger_collection_DM["merger map 0"], 
                                                                        vcat( [   halo_story["M_DM"][ii],
                                                                            halo_story["M_DM_2"][ii],
                                                                            halo_story["M_DM"][merger_index_map[2,pos_in_mim]], 
                                                                            halo_story["M_DM_2"][merger_index_map[2,pos_in_mim]], 
                                                                            halo_story["snapFP"][ii], 
                                                                            halo_story["redshift"][fp_indices[i]], 
                                                                            halo_story["lookbacktime"][fp_indices[i]], 
                                                                            halo_story["subID"][ii],
                                                                            halo_story["M_DM_peak"][ii],
                                                                            halo_story["M_DM_peak"][merger_index_map[2,pos_in_mim]], 
                                                                            halo_story["snapNR"][ii]
                                                                            ], halo_story["merger_spin_map_DM_0"][:,ii] )
                                                                            )
                    # Most Massive condition
                    if halo_story["M_DM_peak"][ii] > merger_collection_DM["Mpeak_MM"][end]
                        merger_collection_DM["M_MM"][end]            = halo_story["M_DM"][ii]
                        merger_collection_DM["M2_MM"][end]           = halo_story["M_DM_2"][ii]
                        merger_collection_DM["Mpeak_MM"][end]           = halo_story["M_DM_peak"][ii]
                        merger_collection_DM["J_MM"][:,end]   = vcat(halo_story["merger_spin_map_DM"][4:6,ii] .*halo_story["merger_spin_map_DM"][3,ii] , 
                                                                        halo_story["merger_spin_map_DM"][7:9,ii], 
                                                                        halo_story["merger_spin_map_DM"][13:15,ii] .*halo_story["merger_spin_map_DM"][12,ii] , 
                                                                        halo_story["merger_spin_map_DM"][16:18,ii]
                                                                        ) 
                        merger_collection_DM["J_MM"][:,end]   = vcat(halo_story["merger_spin_map_DM_0"][4:6,ii] .*halo_story["merger_spin_map_DM_0"][3,ii] , 
                                                                        halo_story["merger_spin_map_DM_0"][7:9,ii], 
                                                                        halo_story["merger_spin_map_DM_0"][13:15,ii] .*halo_story["merger_spin_map_DM_0"][12,ii] , 
                                                                        halo_story["merger_spin_map_DM_0"][16:18,ii]
                                                                        ) 
                    end
                    # Orbital Data
                    merger_collection_DM["M_MERGERS"][end]   += halo_story["M_DM"][ii]
                    merger_collection_DM["M2_MERGERS"][end]  += halo_story["M_DM_2"][ii]
                    merger_collection_DM["Mpeak_MERGERS"][end]  += halo_story["M_DM_peak"][ii]
                    merger_collection_DM["J_SUM"][:,end]       .+= replace(vcat(halo_story["merger_spin_map_DM"][4:6,ii] .*halo_story["merger_spin_map_DM"][3,ii] , 
                                                                                halo_story["merger_spin_map_DM"][7:9,ii], 
                                                                                halo_story["merger_spin_map_DM"][13:15,ii] .*halo_story["merger_spin_map_DM"][12,ii] , 
                                                                                halo_story["merger_spin_map_DM"][16:18,ii]
                                                                                ),
                                                                            missing => 0.0)
                    merger_collection_DM["J_SUM_0"][:,end]       .+= replace(vcat(halo_story["merger_spin_map_DM_0"][4:6,ii] .*halo_story["merger_spin_map_DM_0"][3,ii] , 
                                                                                halo_story["merger_spin_map_DM_0"][7:9,ii], 
                                                                                halo_story["merger_spin_map_DM_0"][13:15,ii] .*halo_story["merger_spin_map_DM_0"][12,ii] , 
                                                                                halo_story["merger_spin_map_DM_0"][16:18,ii]
                                                                                ),
                                                                            missing => 0.0)
                    if count(ismissing, halo_story["merger_spin_map_DM"][:,ii]) == 0
                        merger_collection_DM["M_CONSIDERED"][end]    += halo_story["M_DM"][ii]
                        merger_collection_DM["M2_CONSIDERED"][end]   += halo_story["M_DM_2"][ii]
                        merger_collection_DM["Mpeak_CONSIDERED"][end]   += halo_story["M_DM_peak"][ii]
                    else
                        #println("$ii")
                        merger_collection_DM["M_MISSED"][end]    += halo_story["M_DM"][ii]
                        merger_collection_DM["M2_MISSED"][end]   += halo_story["M_DM_2"][ii]
                        merger_collection_DM["Mpeak_MISSED"][end]   += halo_story["M_DM_peak"][ii]
                    end
                    pos_in_mim += 1
                end
            end
        end
        # Replace zeros with missings
        replace!(merger_collection_DM["M_MM"]         , 0. => missing)
        replace!(merger_collection_DM["M2_MM"]        , 0. => missing)
        replace!(merger_collection_DM["Mpeak_MM"]        , 0. => missing)
        replace!(merger_collection_DM["J_MM"]  , 0. => missing)
        #replace!(merger_collection_DM["J_SUM"] , 0. => missing)
        #replace!(merger_collection_DM["M_MERGERS"]    , 0. => missing)
        #replace!(merger_collection_DM["M_MISSED"]     , 0. => missing)
        #replace!(merger_collection_DM["M_CONSIDERED"] , 0. => missing)
        #replace!(merger_collection_DM["M2_MERGERS"]   , 0. => missing)
        #replace!(merger_collection_DM["M2_MISSED"]    , 0. => missing)
        #replace!(merger_collection_DM["M2_CONSIDERED"], 0. => missing)
        #replace!(merger_collection_DM["Mpeak_MERGERS"]   , 0. => missing)
        #replace!(merger_collection_DM["Mpeak_MISSED"]    , 0. => missing)
        #replace!(merger_collection_DM["Mpeak_CONSIDERED"], 0. => missing)
        #replace!(merger_collection_DM["merger map"][1:4,:], 0. => missing)
    
    
        # Fill the GAS dictionary
    
        for i in 1:length(fp_indices)
            #print("$i ")
            #flush(stdout)/home/moon/sfortune/spinevo/halostories_v20211204_min0.0Gyr/halo_1_1414.jld
    
            # Basic Info
            merger_collection_GAS["lookbacktime"] = vcat( merger_collection_GAS["lookbacktime"], halo_story["lookbacktime"][fp_indices[i]] )
            merger_collection_GAS["M"]      = vcat( merger_collection_GAS["M"], halo_story["M_GAS"][fp_indices[i]] )
            merger_collection_GAS["M2"]     = vcat( merger_collection_GAS["M2"],halo_story["M_GAS_2"][fp_indices[i]])
            merger_collection_GAS["Mpeak"]     = vcat( merger_collection_GAS["Mpeak"],halo_story["M_GAS_peak"][fp_indices[i]])
            if count(ismissing, halo_story["j_GAS"][:,fp_indices[i]]) == 0 
                merger_collection_GAS["M_fromJ"]  = vcat( merger_collection_GAS["M_fromJ"], norm(halo_story["J_GAS"][:,fp_indices[i]]) / norm(halo_story["j_GAS"][:,fp_indices[i]]) )
            else
                merger_collection_GAS["M_fromJ"]  = vcat( merger_collection_GAS["M_fromJ"], missing )
            end
            
            # Transitional Data
            if i == 1
                merger_collection_GAS["ΔJ_main"]          = hcat( merger_collection_GAS["ΔJ_main"], missings(Float64, 3) )
                merger_collection_GAS["Δj_main"]          = hcat( merger_collection_GAS["Δj_main"], missings(Float64, 3) )
                merger_collection_GAS["ΔJ_main_0"]          = hcat( merger_collection_GAS["ΔJ_main_0"], missings(Float64, 3) )
                merger_collection_GAS["Δj_main_0"]          = hcat( merger_collection_GAS["Δj_main_0"], missings(Float64, 3) )
                merger_collection_GAS["ΔJ_vir"]          = hcat( merger_collection_GAS["ΔJ_vir"], missings(Float64, 3) )
                merger_collection_GAS["Δj_vir"]          = hcat( merger_collection_GAS["Δj_vir"], missings(Float64, 3) )
                merger_collection_GAS["ΔJ_vir_0"]          = hcat( merger_collection_GAS["ΔJ_vir_0"], missings(Float64, 3) )
                merger_collection_GAS["Δj_vir_0"]          = hcat( merger_collection_GAS["Δj_vir_0"], missings(Float64, 3) )
                #merger_collection_GAS["ΔBVAL"]            = vcat( merger_collection_GAS["ΔBVAL"], missing )
                #merger_collection_GAS["ΔBVAL_0"]            = vcat( merger_collection_GAS["ΔBVAL_0"], missing )
                merger_collection_GAS["ϕ_flip"]           = vcat( merger_collection_GAS["ϕ_flip"], missing )
                merger_collection_GAS["Δlookbacktime"]    = vcat( merger_collection_GAS["Δlookbacktime"], missing )
                merger_collection_GAS["switch"]           = vcat( merger_collection_GAS["switch"], missing )
                merger_collection_GAS["ΔM"]         = vcat( merger_collection_GAS["ΔM"], missing )
                merger_collection_GAS["ΔM2"]        = vcat( merger_collection_GAS["ΔM2"], missing )
                merger_collection_GAS["ΔMpeak"]        = vcat( merger_collection_GAS["ΔMpeak"], missing )
                merger_collection_GAS["ΔM_fromJ"]         = vcat( merger_collection_GAS["ΔM_fromJ"], missing )
                merger_collection_GAS["M_MM"]             = vcat( merger_collection_GAS["M_MM"], missing )
                merger_collection_GAS["M2_MM"]            = vcat( merger_collection_GAS["M2_MM"], missing )
                merger_collection_GAS["Mpeak_MM"]            = vcat( merger_collection_GAS["Mpeak_MM"], missing )
                merger_collection_GAS["J_MM"]      = hcat( merger_collection_GAS["J_MM"], missings(Float64, 12) )
                merger_collection_GAS["J_MM_0"]      = hcat( merger_collection_GAS["J_MM_0"], missings(Float64, 12) )
                merger_collection_GAS["J_SUM"]     = hcat( merger_collection_GAS["J_SUM"], zeros(12) )
                merger_collection_GAS["J_SUM_0"]     = hcat( merger_collection_GAS["J_SUM_0"], zeros(12) )
                merger_collection_GAS["M_MERGERS"]        = vcat( merger_collection_GAS["M_MERGERS"], missing )
                merger_collection_GAS["M_MISSED"]         = vcat( merger_collection_GAS["M_MISSED"], missing )
                merger_collection_GAS["M_CONSIDERED"]     = vcat( merger_collection_GAS["M_CONSIDERED"], missing )
                merger_collection_GAS["M2_MERGERS"]       = vcat( merger_collection_GAS["M2_MERGERS"], missing )
                merger_collection_GAS["M2_MISSED"]        = vcat( merger_collection_GAS["M2_MISSED"], missing )
                merger_collection_GAS["M2_CONSIDERED"]    = vcat( merger_collection_GAS["M2_CONSIDERED"], missing )
                merger_collection_GAS["Mpeak_MERGERS"]       = vcat( merger_collection_GAS["Mpeak_MERGERS"], missing )
                merger_collection_GAS["Mpeak_MISSED"]        = vcat( merger_collection_GAS["Mpeak_MISSED"], missing )
                merger_collection_GAS["Mpeak_CONSIDERED"]    = vcat( merger_collection_GAS["Mpeak_CONSIDERED"], missing )
            else
                mintidx = find_mintime_idx( halo_story["lookbacktime"], fp_indices, i, min_time )
                merger_collection_GAS["ΔJ_main"]  = hcat( merger_collection_GAS["ΔJ_main"], halo_story["J_GAS"][:,fp_indices[i]] .- halo_story["J_GAS"][:,mintidx] )
                merger_collection_GAS["Δj_main"]  = hcat( merger_collection_GAS["Δj_main"], halo_story["j_GAS"][:,fp_indices[i]] .- halo_story["j_GAS"][:,mintidx] )
                merger_collection_GAS["ΔJ_main_0"]  = hcat( merger_collection_GAS["ΔJ_main_0"], halo_story["J_GAS_0"][:,fp_indices[i]] .- halo_story["J_GAS_0"][:,mintidx] )
                merger_collection_GAS["Δj_main_0"]  = hcat( merger_collection_GAS["Δj_main_0"], halo_story["j_GAS_0"][:,fp_indices[i]] .- halo_story["j_GAS_0"][:,mintidx] )
                merger_collection_GAS["ΔJ_vir"]  = hcat( merger_collection_GAS["ΔJ_vir"], halo_story["J_GASvir"][:,fp_indices[i]] .- halo_story["J_GASvir"][:,mintidx] )
                merger_collection_GAS["Δj_vir"]  = hcat( merger_collection_GAS["Δj_vir"], halo_story["j_GASvir"][:,fp_indices[i]] .- halo_story["j_GASvir"][:,mintidx] )
                merger_collection_GAS["ΔJ_vir_0"]  = hcat( merger_collection_GAS["ΔJ_vir_0"], halo_story["J_GASvir_0"][:,fp_indices[i]] .- halo_story["J_GASvir_0"][:,mintidx] )
                merger_collection_GAS["Δj_vir_0"]  = hcat( merger_collection_GAS["Δj_vir_0"], halo_story["j_GASvir_0"][:,fp_indices[i]] .- halo_story["j_GASvir_0"][:,mintidx] )
                merger_collection_GAS["Δlookbacktime"]  = vcat( merger_collection_GAS["Δlookbacktime"], halo_story["lookbacktime"][fp_indices[i]] - halo_story["lookbacktime"][mintidx] )
                merger_collection_GAS["switch"]           = vcat( merger_collection_GAS["switch"], sum( halo_story["switch"][fp_indices[findcs( fp_indices, geq=fp_indices[i], lt=mintidx)]] ) )
                #merger_collection_GAS["ΔBVAL"]    = vcat( merger_collection_GAS["ΔBVAL"], halo_story["BVAL"][fp_indices[i]] - halo_story["BVAL"][mintidx] )
                #merger_collection_GAS["ΔBVAL_0"]    = vcat( merger_collection_GAS["ΔBVAL_0"], halo_story["BVAL_0"][fp_indices[i]] - halo_story["BVAL_0"][mintidx] )
                if count(ismissing, halo_story["j_GAS"][:,fp_indices[i]]) == 0 && count(ismissing, halo_story["j_GAS"][:,mintidx]) == 0 
                    merger_collection_GAS["ϕ_flip"]   = vcat( merger_collection_GAS["ϕ_flip"], 180/π * angle( halo_story["J_GAS"][:,fp_indices[i]], halo_story["J_GAS"][:,mintidx] ) )
                    merger_collection_GAS["ΔM_fromJ"] = vcat( merger_collection_GAS["ΔM_fromJ"], (norm(halo_story["J_GAS"][:,fp_indices[i]]) / norm(halo_story["j_GAS"][:,fp_indices[i]])) - (norm(halo_story["J_GAS"][:,mintidx]) / norm(halo_story["j_GAS"][:,mintidx])) )
                else
                    #println(halo_story["j_GAS"][:,fp_indices[i]])
                    #println("missing j_GAS for $(storyfilelist[iii]) @ index $(fp_indices[i])")
                    merger_collection_GAS["ϕ_flip"]   = vcat( merger_collection_GAS["ϕ_flip"], missing )
                    merger_collection_GAS["ΔM_fromJ"] = vcat( merger_collection_GAS["ΔM_fromJ"], missing )
                end
                #println("\n$(halo_story["J_GAS"][:,fp_indices[i]]) $(halo_story["J_GAS"][:,mintidx])")
                merger_collection_GAS["ΔM"] = vcat( merger_collection_GAS["ΔM"], halo_story["M_GAS"][fp_indices[i]] - halo_story["M_GAS"][mintidx] )
                merger_collection_GAS["ΔM2"] = vcat( merger_collection_GAS["ΔM2"],halo_story["M_GAS_2"][fp_indices[i]] - halo_story["M_GAS_2"][mintidx] )
                merger_collection_GAS["ΔMpeak"] = vcat( merger_collection_GAS["ΔMpeak"],halo_story["M_GAS_peak"][fp_indices[i]] - halo_story["M_GAS_peak"][mintidx] )
                
                # Merger Data
                # Setup
                merger_collection_GAS["M_MM"]             = vcat( merger_collection_GAS["M_MM"], 0 )
                merger_collection_GAS["M2_MM"]            = vcat( merger_collection_GAS["M2_MM"], 0 )
                merger_collection_GAS["Mpeak_MM"]            = vcat( merger_collection_GAS["Mpeak_MM"], 0 )
                merger_collection_GAS["J_MM"]      = hcat( merger_collection_GAS["J_MM"], missings(Float64, 12) )
                merger_collection_GAS["J_MM_0"]      = hcat( merger_collection_GAS["J_MM_0"], missings(Float64, 12) )
                merger_collection_GAS["J_SUM"]     = hcat( merger_collection_GAS["J_SUM"], zeros(12) )
                merger_collection_GAS["J_SUM_0"]     = hcat( merger_collection_GAS["J_SUM_0"], zeros(12) )
                merger_collection_GAS["M_MERGERS"]        = vcat( merger_collection_GAS["M_MERGERS"], 0 )
                merger_collection_GAS["M_MISSED"]         = vcat( merger_collection_GAS["M_MISSED"], 0 )
                merger_collection_GAS["M_CONSIDERED"]     = vcat( merger_collection_GAS["M_CONSIDERED"], 0 )
                merger_collection_GAS["M2_MERGERS"]       = vcat( merger_collection_GAS["M2_MERGERS"], 0 )
                merger_collection_GAS["M2_MISSED"]        = vcat( merger_collection_GAS["M2_MISSED"], 0 )
                merger_collection_GAS["M2_CONSIDERED"]    = vcat( merger_collection_GAS["M2_CONSIDERED"], 0 )
                merger_collection_GAS["Mpeak_MERGERS"]       = vcat( merger_collection_GAS["Mpeak_MERGERS"], 0 )
                merger_collection_GAS["Mpeak_MISSED"]        = vcat( merger_collection_GAS["Mpeak_MISSED"], 0 )
                merger_collection_GAS["Mpeak_CONSIDERED"]    = vcat( merger_collection_GAS["Mpeak_CONSIDERED"], 0 )
                # Actual check
                pos_in_mim  = 1+sum(n_mergers[1:i-1]) # position / index in merger index map
                for ii in merger_index_map[1,1+sum(n_mergers[1:i-1]):sum(n_mergers[1:i])]
                    # Merger Map
                    #println("$(halo_story["snapNR"][ii])   i = $i")
                    #@show merger_index_map
                    #@show merger_index_map[1,1+sum(n_mergers[1:i-1]):sum(n_mergers[1:i])]
                    #flush(stdout)
                    merger_collection_GAS["merger map 0"]   = hcat( merger_collection_GAS["merger map 0"], 
                                                                        vcat( [   halo_story["M_GAS"][ii],
                                                                            halo_story["M_GAS_2"][ii], 
                                                                            halo_story["M_GAS"][merger_index_map[2,pos_in_mim]], 
                                                                            halo_story["M_GAS_2"][merger_index_map[2,pos_in_mim]], 
                                                                            halo_story["snapFP"][ii], 
                                                                            halo_story["redshift"][fp_indices[i]], 
                                                                            halo_story["lookbacktime"][fp_indices[i]] , 
                                                                            halo_story["subID"][ii],
                                                                            halo_story["M_GAS_peak"][ii], 
                                                                            halo_story["M_GAS_peak"][merger_index_map[2,pos_in_mim]], 
                                                                            halo_story["snapNR"][ii]
                                                                            ], halo_story["merger_spin_map_GAS_0"][:,ii] )
                                                                            )
                    merger_collection_GAS["merger map"]   = hcat( merger_collection_GAS["merger map"], 
                                                                        vcat( [   halo_story["M_GAS"][ii],
                                                                            halo_story["M_GAS_2"][ii], 
                                                                            halo_story["M_GAS"][merger_index_map[2,pos_in_mim]], 
                                                                            halo_story["M_GAS_2"][merger_index_map[2,pos_in_mim]], 
                                                                            halo_story["snapFP"][ii], 
                                                                            halo_story["redshift"][fp_indices[i]], 
                                                                            halo_story["lookbacktime"][fp_indices[i]] , 
                                                                            halo_story["subID"][ii],
                                                                            halo_story["M_GAS_peak"][ii], 
                                                                            halo_story["M_GAS_peak"][merger_index_map[2,pos_in_mim]], 
                                                                            halo_story["snapNR"][ii]
                                                                            ], halo_story["merger_spin_map_GAS"][:,ii] )
                                                                            )
                    # Most Massive condition
                    if halo_story["M_GAS_peak"][ii] > merger_collection_GAS["Mpeak_MM"][end]
                        merger_collection_GAS["M_MM"][end]            = halo_story["M_GAS"][ii]
                        merger_collection_GAS["M2_MM"][end]           = halo_story["M_GAS_2"][ii]
                        merger_collection_GAS["Mpeak_MM"][end]           = halo_story["M_GAS_peak"][ii]
                        merger_collection_GAS["J_MM"][:,end]   = vcat(halo_story["merger_spin_map_GAS"][4:6,ii] .*halo_story["merger_spin_map_GAS"][3,ii] , 
                                                                        halo_story["merger_spin_map_GAS"][7:9,ii], 
                                                                        halo_story["merger_spin_map_GAS"][13:15,ii] .*halo_story["merger_spin_map_GAS"][12,ii] , 
                                                                        halo_story["merger_spin_map_GAS"][16:18,ii]
                                                                        ) 
                        merger_collection_GAS["J_MM"][:,end]   = vcat(halo_story["merger_spin_map_GAS_0"][4:6,ii] .*halo_story["merger_spin_map_GAS_0"][3,ii] , 
                                                                        halo_story["merger_spin_map_GAS_0"][7:9,ii], 
                                                                        halo_story["merger_spin_map_GAS_0"][13:15,ii] .*halo_story["merger_spin_map_GAS_0"][12,ii] , 
                                                                        halo_story["merger_spin_map_GAS_0"][16:18,ii]
                                                                        ) 
                    end
                    # Orbital Data
                    merger_collection_GAS["M_MERGERS"][end]   += halo_story["M_GAS"][ii]
                    merger_collection_GAS["M2_MERGERS"][end]  += halo_story["M_GAS_2"][ii]
                    merger_collection_GAS["Mpeak_MERGERS"][end]  += halo_story["M_GAS_peak"][ii]
                    merger_collection_GAS["J_SUM"][:,end]       .+= replace(vcat(halo_story["merger_spin_map_GAS"][4:6,ii] .*halo_story["merger_spin_map_GAS"][3,ii] , 
                                                                                halo_story["merger_spin_map_GAS"][7:9,ii], 
                                                                                halo_story["merger_spin_map_GAS"][13:15,ii] .*halo_story["merger_spin_map_GAS"][12,ii] , 
                                                                                halo_story["merger_spin_map_GAS"][16:18,ii]
                                                                                ),
                                                                            missing => 0.0)
                    merger_collection_GAS["J_SUM_0"][:,end]       .+= replace(vcat(halo_story["merger_spin_map_GAS_0"][4:6,ii] .*halo_story["merger_spin_map_GAS_0"][3,ii] , 
                                                                                halo_story["merger_spin_map_GAS_0"][7:9,ii], 
                                                                                halo_story["merger_spin_map_GAS_0"][13:15,ii] .*halo_story["merger_spin_map_GAS_0"][12,ii] , 
                                                                                halo_story["merger_spin_map_GAS_0"][16:18,ii]
                                                                                ),
                                                                            missing => 0.0)
                    if count(ismissing, halo_story["merger_spin_map_GAS"][:,ii]) == 0
                        merger_collection_GAS["M_CONSIDERED"][end]    += halo_story["M_GAS"][ii]
                        merger_collection_GAS["M2_CONSIDERED"][end]   += halo_story["M_GAS_2"][ii]
                        merger_collection_GAS["Mpeak_CONSIDERED"][end]   += halo_story["M_GAS_peak"][ii]
                    else
                        #println("$ii")
                        merger_collection_GAS["M_MISSED"][end]    += halo_story["M_GAS"][ii]
                        merger_collection_GAS["M2_MISSED"][end]   += halo_story["M_GAS_2"][ii]
                        merger_collection_GAS["Mpeak_MISSED"][end]   += halo_story["M_GAS_peak"][ii]
                    end
                    pos_in_mim += 1
                end
            end
        end
        # Replace zeros with missings
        replace!(merger_collection_GAS["M_MM"]         , 0. => missing)
        replace!(merger_collection_GAS["M2_MM"]        , 0. => missing)
        replace!(merger_collection_GAS["Mpeak_MM"]        , 0. => missing)
        replace!(merger_collection_GAS["J_MM"]  , 0. => missing)
        #replace!(merger_collection_GAS["J_SUM"] , 0. => missing)
        #replace!(merger_collection_GAS["M_MERGERS"]    , 0. => missing)
        #replace!(merger_collection_GAS["M_MISSED"]     , 0. => missing)
        #replace!(merger_collection_GAS["M_CONSIDERED"] , 0. => missing)
        #replace!(merger_collection_GAS["M2_MERGERS"]   , 0. => missing)
        #replace!(merger_collection_GAS["M2_MISSED"]    , 0. => missing)
        #replace!(merger_collection_GAS["M2_CONSIDERED"], 0. => missing)
        #replace!(merger_collection_GAS["Mpeak_MERGERS"]   , 0. => missing)
        #replace!(merger_collection_GAS["Mpeak_MISSED"]    , 0. => missing)
        #replace!(merger_collection_GAS["Mpeak_CONSIDERED"], 0. => missing)
        #replace!(merger_collection_GAS["merger map"][1:4,:], 0. => missing)
    
        # return
        return merger_collection_STARS, merger_collection_GAS, merger_collection_DM

    end
end

print("'collect_mergers'   ")







