### Settings

box             = "/HydroSims/Magneticum/Box4/uhr_test"
input_dir       = "/home/moon/sfortune/spinevo/halostories"
output_dir      = "/home/moon/sfortune/spinevo/halostories_update_stars"


### Packages

using Printf
using LinearAlgebra
using JLD
using QuadGK
using Statistics
using PyCall
using PyPlot
using LaTeXStrings
using GadgetIO
using GadgetUnits
using GadgetGalaxies
using Unitful
using UnitfulAstro
using Missings
using HypothesisTests
using Distributions
using CSV
using DataFrames
using Cosmology
using Suppressor



storyfilelist   = readdir(input_dir)




limit_filelist  = length(storyfilelist)
println("$(limit_filelist) Runs in total.")
for ii in 1:limit_filelist
    halo_story  = load(joinpath(input_dir, "$(storyfilelist[ii])"), "halo_story")
    head        = read_header("$box/groups_$(@sprintf("%03i", halo_story["SNAP"][1]))/sub_$(@sprintf("%03i", halo_story["SNAP"][1]))")
    snapshot    = Snapshot(box, halo_story["SNAP"][1])
    snapshot.snapbase
    snapshot.subbase
    g           = Galaxy(snapshot, halo_story["I_SUB"][1])
    if get_first_subhalo(get_group(g)).isub == halo_story["I_SUB"][1]
        println("Starting RUN $ii (Halo $(halo_story["I_SUB"][1]))")
        flush(stdout)
        
        # Read and correct orbital spins for stellar part only
        for i in 1:length(halo_story["SNAP"])
            if !ismissing(halo_story["J_orbital"][1,i])
                path_to_file    = "$box/groups_$(@sprintf("%03i", halo_story["SNAP"][i]))/sub_$(@sprintf("%03i", halo_story["SNAP"][i]))"
                head            = read_header("$box/groups_$(@sprintf("%03i", halo_story["SNAP"][i]))/sub_$(@sprintf("%03i", halo_story["SNAP"][i]))")
                loop_smst       = convert_units_physical(read_subfind(path_to_file, "SMST"), :mass, head)
                halo_story["J_orbital"][:,i]    ./= ( (loop_smst[1,halo_story["I_SUB"][i]+1] + loop_smst[2,halo_story["I_SUB"][i]+1] + loop_smst[5,halo_story["I_SUB"][i]+1]) )
            end
        end
        merger_data = Dict(
            "J_main"        => Array{Float64}(undef, 3, 0), 
            "J_main_abs"    => Array{Float64}(undef, 0), 
            "j_main_abs"    => Array{Float64}(undef, 0), 
            "j_main"        => Array{Float64}(undef, 3, 0), 
            "J_orbits"      => Array{Float64}(undef, 3, 0), 
            "J_sum"         => Array{Float64}(undef, 3, 0), 
            "J_sum_abs"     => Array{Float64}(undef, 0), 
            "J_int_sum"     => Array{Float64}(undef, 3, 0), 
            "J_int_sum_abs" => Array{Float64}(undef, 0), 
            "ANG_neighbors" => Array{Float64}(undef, 0), 
            "ANG_main"      => Array{Float64}(undef, 0), 
            "δJ_main"       => Array{Float64}(undef, 3, 0), 
            "δJ_main_abs"   => Array{Float64}(undef, 0), 
            "δj_main"       => Array{Float64}(undef, 3, 0), 
            "N_MERGERS"     => Array{Int64}(undef, 0),
            "M_MERGERS"     => Array{Float64}(undef, 0),
            "M_felix"       => Array{Float64}(undef, 0), 
            "M2_felix"      => Array{Float64}(undef, 0), 
            "M_fromJ"       => Array{Float64}(undef, 0), 
            "REDSHIFT"      => Array{Float64}(undef, 0),
            "SNAP"          => Array{Int64}(undef, 0),
            "ADDED_MERGERS" => Array{Int64}(undef, 0),
            "MISSING_MERGERS" => Array{Int64}(undef, 0),
            "LOOKBACKTIME"  => Array{Float64}(undef, 0),
            "I_FILE"        => parse(Int64, chop(storyfilelist[ii], head=5, tail=4)), 
            "M_MISSED"      => 0., 
            "M_CONSIDERED"  => 0., 
            "I_SUB"         => Array{Int64}(undef, 0) )
        
        
        
        limit_loop      = length(halo_story["I_SUB"])
        loop_nmergers   = 0
        loop_mmergers   = 0.
        loop_J_sum      = zeros(3)
        loop_J_int_sum  = zeros(3)
        mass_missed     = 0.
        mass_added      = 0.
        for i in 1:limit_loop
            #print("$i ")
            #flush(stdout)
            if !ismissing(halo_story["J_STARS"][1,i]) && !ismissing(halo_story["J_STARS"][2,i]) && !ismissing(halo_story["J_STARS"][3,i])   # most massive
                # Check out the previous merger history
                merger_data["N_MERGERS"]    = vcat( merger_data["N_MERGERS"],       loop_nmergers )
                merger_data["M_MERGERS"]    = vcat( merger_data["M_MERGERS"],       loop_mmergers )
                merger_data["J_sum"]        = hcat( merger_data["J_sum"],           loop_J_sum )
                merger_data["J_sum_abs"]    = vcat( merger_data["J_sum_abs"],       norm(loop_J_sum) )
                merger_data["J_int_sum"]    = hcat( merger_data["J_int_sum"],       loop_J_int_sum )
                merger_data["J_int_sum_abs"]= vcat( merger_data["J_int_sum_abs"],   norm(loop_J_int_sum) )
                # Append Central halo Info
                head    = read_header("$box/groups_$(@sprintf("%03i", halo_story["SNAP"][i]))/sub_$(@sprintf("%03i", halo_story["SNAP"][i]))")
                merger_data["J_main"]       = hcat( merger_data["J_main"],          halo_story["J_STARS"][:,i] )
                merger_data["j_main"]       = hcat( merger_data["j_main"],          halo_story["j_STARS"][:,i] )
                merger_data["J_main_abs"]   = vcat( merger_data["J_main_abs"],      norm(halo_story["J_STARS"][:,i]) )
                merger_data["j_main_abs"]   = vcat( merger_data["j_main_abs"],      norm(halo_story["j_STARS"][:,i]) )
                merger_data["M_felix"]      = vcat( merger_data["M_felix"],         convert_units_physical_mass(halo_story["M_STARS"][i], head) )
                merger_data["M2_felix"]     = vcat( merger_data["M2_felix"],        convert_units_physical_mass(halo_story["M_STAR_2"][i], head) )
                merger_data["M_fromJ"]      = vcat( merger_data["M_fromJ"],         norm(halo_story["J_STARS"][:,i]) / norm(halo_story["j_STARS"][:,i]) )
                merger_data["LOOKBACKTIME"] = vcat( merger_data["LOOKBACKTIME"],    ustrip(lookback_time(cosmology(h=head.h0, OmegaM=head.omega_0), head.z)) )
                merger_data["REDSHIFT"]     = vcat( merger_data["REDSHIFT"],        halo_story["REDSHIFT"][i] )
                merger_data["I_SUB"]        = vcat( merger_data["I_SUB"],           halo_story["I_SUB"][i] )
                merger_data["SNAP"]         = vcat( merger_data["SNAP"],            halo_story["SNAP"][i] )
                # Reset merger values
                loop_J_sum      = zeros(3)
                loop_J_int_sum  = zeros(3)
                loop_nmergers   = 0
                loop_mmergers   = 0.
            elseif !ismissing(halo_story["J_orbital"][1,i]) && !ismissing(halo_story["J_orbital"][2,i]) && !ismissing(halo_story["J_orbital"][3,i]) # existing neighbor
                loop_nmergers  += 1
                loop_mmergers  += convert_units_physical_mass(halo_story["M_STAR_2"][i], head)
                loop_J_sum     += halo_story["J_orbital"][:,i]
        
                head    = read_header("$box/groups_$(@sprintf("%03i", halo_story["SNAP"][i]))/sub_$(@sprintf("%03i", halo_story["SNAP"][i]))")
                try
                    snapshot    = Snapshot(box, halo_story["SNAP"][i])
                    snapshot.snapbase
                    snapshot.subbase
                    g           = Galaxy(snapshot, halo_story["I_SUB"][i])
                    @suppress begin
                        read_halo!(g, units=:physical, props=((:stars, ["POS", "VEL", "MASS"]),))
                    end
                    halo_story["J_STARS"][:,i]  = angular_momentum(g.stars)#, sph_small)
                    halo_story["j_STARS"][:,i]  = specific_angular_momentum(g.stars)#, sph_small)
                    loop_J_int_sum += halo_story["J_STARS"][:,i]  
                    merger_data["ADDED_MERGERS"]    = vcat( merger_data["ADDED_MERGERS"], i)
                    mass_added  += convert_units_physical_mass(halo_story["M_STAR_2"][i], head)
                catch e
                    #println(e)
                    merger_data["MISSING_MERGERS"]  = vcat( merger_data["MISSING_MERGERS"], i)
                    mass_missed += convert_units_physical_mass(halo_story["M_STAR_2"][i], head)
                end
            elseif i > 1
                if halo_story["SNAP"][i-1] < halo_story["SNAP"][i] || !ismissing(halo_story["J_STARS"][1,i]) || !ismissing(halo_story["J_STARS"][2,i]) || !ismissing(halo_story["J_STARS"][3,i])
                    println("Problem with line $i")
                end
            end
        end
        # Final Check out of previous merger history
        merger_data["N_MERGERS"]    = vcat( merger_data["N_MERGERS"],       loop_nmergers )
        merger_data["M_MERGERS"]    = vcat( merger_data["M_MERGERS"],       loop_mmergers )
        merger_data["J_sum"]        = hcat( merger_data["J_sum"],           loop_J_sum )
        merger_data["J_sum_abs"]    = vcat( merger_data["J_sum_abs"],       norm(loop_J_sum) )
        merger_data["J_int_sum"]    = hcat( merger_data["J_int_sum"],       loop_J_int_sum )
        merger_data["J_int_sum_abs"]= vcat( merger_data["J_int_sum_abs"],   norm(loop_J_int_sum) )
        
        # Calculate changes
        for i in 1:length(merger_data["SNAP"])-1
            if norm(merger_data["J_sum"][:,i+1]) > 0
                merger_data["ANG_neighbors"]    = vcat( merger_data["ANG_neighbors"],
                    acosd((transpose(merger_data["J_main"][:,i+1])*merger_data["J_sum"][:,i+1]) / (norm(merger_data["J_main"][:,i+1])*norm(merger_data["J_sum"][:,i+1]))) )
            else
                merger_data["ANG_neighbors"]    = vcat( merger_data["ANG_neighbors"], 0. )
            end
            merger_data["ANG_main"] = vcat( merger_data["ANG_main"],
                acosd((transpose(merger_data["J_main"][:,i])*merger_data["J_main"][:,i+1]) / (norm(merger_data["J_main"][:,i])*norm(merger_data["J_main"][:,i+1]))) )
            merger_data["δJ_main"]  = hcat( merger_data["δJ_main"],
                merger_data["J_main"][:,i] .- merger_data["J_main"][:,i+1] )
            merger_data["δj_main"]  = hcat( merger_data["δj_main"],
                merger_data["j_main"][:,i] .- merger_data["j_main"][:,i+1] )
            merger_data["δJ_main_abs"]  = vcat( merger_data["δJ_main_abs"],
                norm(merger_data["J_main"][:,i] .- merger_data["J_main"][:,i+1]) )
        end
        
        merger_data["M_MISSED"]     = mass_missed
        merger_data["M_CONSIDERED"] = mass_added
        println("RUN $ii (Halo $(halo_story["I_SUB"][1])): Missed Mass = $mass_missed   ---  Considered Mass = $mass_added")
        flush(stdout)


        # Rename J_orbital to j_orbital since no mass defined
        halo_story["j_orbital"] = halo_story["J_orbital"]
        delete!(halo_story, "J_orbital")

        
        save(joinpath(output_dir, "halo_$(merger_data["I_FILE"])_$(halo_story["I_SUB"][1]).jld"), 
            "halo_story",   halo_story,
            "merger_data",  merger_data)
    end
end



println("\n\nNow witness the firepower of this fully armed and operational battle station!")