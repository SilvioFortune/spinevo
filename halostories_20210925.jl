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


### Functions

# Tranform Cartesian vector into [r,θ,ϕ] Spherical vector
function cartesian_to_spherical(x)
    s = zeros(3)
    s[1]    = norm(x) # Radius
    s[2]    = acosd(x[3] / s[1])# * 180 / π # θ[°]
    s[3]    = acosd(x[1] / sqrt(x[1]*x[1] + x[2]*x[2]))# * 180 / π # ϕ[°]
    return s
end

# Rotation Matrix to align x with reference vector
function aligner(x, ref=[0,0,1])
    x       = x ./ norm(x)
    ref     = ref ./ norm(ref)
    vx      = x × ref
    cx      = transpose(x) * ref
    Vx      = [0 -vx[3] vx[2]; vx[3] 0 -vx[1]; -vx[2] vx[1] 0]
    aligner = I + Vx + (Vx * Vx ./ (1+cx))
    return aligner
end

function orbit_J(subID, centID, pathtofile) # Find central sub maybe using SUBFIND -> centID = FSUB[GRNR[subID+1]+1]
    # Read subhalo features
    head        = read_header(pathtofile)
    smst        = convert_units_physical(read_subfind(pathtofile, "SMST"), :mass, head)
    spos        = convert_units_physical(read_subfind(pathtofile, "SPOS"), :pos, head)
    svel        = read_subfind(pathtofile, "SVEL") # no conversion since already physical units??
    # Return with masses added
    return ((spos[:,subID+1] .- spos[:,centID+1]) × (svel[:,subID+1] .- svel[:,centID+1])) .* (smst[1,subID+1]+smst[2,subID+1]+smst[5,subID+1])
end



### Settings
box         = "/HydroSims/Magneticum/Box4/uhr_test"
input_dir   = "/home/moon/sfortune/spinevo/halostories_update_stars"
output_dir  = "/home/moon/sfortune/spinevo/halostories_20210925"
min_time    = 600



### Main

storyfilelist   = readdir(input_dir)


limit_filelist  = length(storyfilelist)
println("$(limit_filelist) Runs in total.")
flush(stdout)
for iii in 1:limit_filelist
    println("$iii")
    flush(stdout)
    merger_data = load(joinpath(input_dir, storyfilelist[iii]), "merger_data")
    halo_story  = load(joinpath(input_dir, storyfilelist[iii]), "halo_story")

    head        = read_header("$box/groups_$(@sprintf("%03i", halo_story["SNAP"][1]))/sub_$(@sprintf("%03i", halo_story["SNAP"][1]))")
    id_mfelix   = convert_units_physical_mass(halo_story["M_STARS"][1], head)
    id_m2       = convert_units_physical_mass(halo_story["M_STAR_2"][1], head)
    
    
    # Identify First Progenitors and mergers

    loop_nmergers   = 0
    loop_mmergers   = 0.
    loop_J_sum      = zeros(3)
    mass_missed     = 0.
    mass_added      = 0.
    
    fp_indices      = Array{Int64}(undef, 0)
    fp_snaps        = Array{Int64}(undef, 0)
    n_mergers       = Array{Int64}(undef, 0)
    merger_indices  = Array{Int64}(undef, 0)
    
    loop_length     = length(halo_story["SNAP"])
    # First Progenitors, forward in time by starting at the bottom
    for i in loop_length:-1:1
        if i==1 || halo_story["SNAP"][i-1] > halo_story["SNAP"][i] && !ismissing(halo_story["J_STARS"][1,i]) # First Progenitor
            #print("FP@SNAP$(halo_story["SNAP"][i]) ")
            #flush(stdout)
    
            fp_indices  = vcat( fp_indices, i )
            fp_snaps    = vcat( fp_snaps, halo_story["SNAP"][i] )
        end
    end
    
    
    merger_count    = 0
    # check if first snap already contains a first progenitor
    if fp_snaps[1] == halo_story["SNAP"][end]
        fp_index    = 2
        n_mergers   = vcat( n_mergers, merger_count )
    else
        fp_index    = 1
    end
    # Merger Map, forward in time by starting at the bottom
     for i in loop_length:-1:1
        # Check out Merger Info and assign to next
        if halo_story["SNAP"][i] == fp_snaps[fp_index] # First Progenitor
            n_mergers   = vcat( n_mergers, merger_count )
            merger_count = 0
            fp_index   += 1
        end
    
        # Identify Mergers
        if i == 1 || halo_story["SNAP"][i] < halo_story["SNAP"][i-1]
            print("")
        elseif halo_story["SNAP"][i] == halo_story["SNAP"][i-1] # Merger
            merger_count   += 1
            merger_indices  = vcat( merger_indices, i )
        else
            println("Error for i = $i, $(halo_story["SNAP"][i]), $(fp_snaps[fp_index]), $(halo_story["SNAP"][i-1]), $(halo_story["I_SUB"][i-1])")
        end
    end
    
    #println("$(length(fp_indices)) $(length(fp_snaps)) $(length(n_mergers))")
    #println("$(length(merger_indices)) $(sum(n_mergers)) ")
    
    
    # Fill the dictionary

    merger_collection_STARS = Dict(
            "SNAP"          => missings(Int64   , 0),
            "ID_ISUB"       => missings(Int64   , 0),
            "I_SUB"         => missings(Int64   , 0),
            "N_MERGERS"     => missings(Int64   , 0),
            "ID_Mfelix"     => missings(Float64 , 0),
            "ID_M2"         => missings(Float64 , 0),
            "REDSHIFT"      => missings(Float64 , 0),
            "LOOKBACKTIME"  => missings(Float64 , 0),
            "M_MM"          => missings(Float64 , 0),
            "M2_MM"         => missings(Float64 , 0),
            "δM_felix"      => missings(Float64 , 0), 
            "δM2_felix"     => missings(Float64 , 0), 
            "δM_fromJ"      => missings(Float64 , 0), 
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
            "BVAL"          => missings(Float64 , 0), 
            "δBVAL"         => missings(Float64 , 0), 
            "J_MMorbital"   => missings(Float64 , 3, 0), 
            "J_SUMorbital"  => missings(Float64 , 3, 0), 
            "δJ_main"       => missings(Float64 , 3, 0), 
            "J_main"        => missings(Float64 , 3, 0), 
            "j_main"        => missings(Float64 , 3, 0), 
            "δj_main"       => missings(Float64 , 3, 0))

    for i in 1:length(fp_indices)
        #print("$i ")
        #flush(stdout)

        # Basic Info
        head    = read_header("$box/groups_$(@sprintf("%03i", halo_story["SNAP"][fp_indices[i]]))/sub_$(@sprintf("%03i", halo_story["SNAP"][fp_indices[i]]))")
        merger_collection_STARS["LOOKBACKTIME"] = vcat( merger_collection_STARS["LOOKBACKTIME"], ustrip(lookback_time(cosmology(h=head.h0, OmegaM=head.omega_0), head.z)) )
        merger_collection_STARS["SNAP"]         = vcat( merger_collection_STARS["SNAP"], halo_story["SNAP"][fp_indices[i]] )
        merger_collection_STARS["ID_ISUB"]      = vcat( merger_collection_STARS["ID_ISUB"], halo_story["I_SUB"][1] )
        merger_collection_STARS["ID_Mfelix"]    = vcat( merger_collection_STARS["ID_Mfelix"], id_mfelix )
        merger_collection_STARS["ID_M2"]        = vcat( merger_collection_STARS["ID_M2"], id_m2 )
        merger_collection_STARS["REDSHIFT"]     = vcat( merger_collection_STARS["REDSHIFT"], halo_story["REDSHIFT"][fp_indices[i]] )
        merger_collection_STARS["I_SUB"]        = vcat( merger_collection_STARS["I_SUB"], halo_story["I_SUB"][fp_indices[i]] )
        merger_collection_STARS["M_felix"]      = vcat( merger_collection_STARS["M_felix"], convert_units_physical_mass(halo_story["M_STARS"][fp_indices[i]], head) )
        merger_collection_STARS["M2_felix"]     = vcat( merger_collection_STARS["M2_felix"], convert_units_physical_mass(halo_story["M_STAR_2"][fp_indices[i]], head) )
        if !ismissing(halo_story["j_STARS"][1,fp_indices[i]]) && !ismissing(halo_story["j_STARS"][2,fp_indices[i]]) && !ismissing(halo_story["j_STARS"][3,fp_indices[i]])
            merger_collection_STARS["M_fromJ"]  = vcat( merger_collection_STARS["M_fromJ"], norm(halo_story["J_STARS"][:,fp_indices[i]]) / norm(halo_story["j_STARS"][:,fp_indices[i]]) )
        else
            merger_collection_STARS["M_fromJ"]  = vcat( merger_collection_STARS["M_fromJ"], missing )
        end
        merger_collection_STARS["BVAL"]         = vcat( merger_collection_STARS["BVAL"], halo_story["BVAL"][fp_indices[i]] )
        merger_collection_STARS["J_main"]       = hcat( merger_collection_STARS["J_main"], halo_story["J_STARS"][:,fp_indices[i]] )
        merger_collection_STARS["j_main"]       = hcat( merger_collection_STARS["j_main"], halo_story["j_STARS"][:,fp_indices[i]] )
        merger_collection_STARS["N_MERGERS"]    = vcat( merger_collection_STARS["N_MERGERS"], length(merger_indices[1+sum(n_mergers[1:i-1]):sum(n_mergers[1:i])]) )

        # Transitional Data
        if i == 1
            merger_collection_STARS["δJ_main"]          = hcat( merger_collection_STARS["δJ_main"], missings(Float64, 3) )
            merger_collection_STARS["δj_main"]          = hcat( merger_collection_STARS["δj_main"], missings(Float64, 3) )
            merger_collection_STARS["δBVAL"]            = vcat( merger_collection_STARS["δBVAL"], missing )
            merger_collection_STARS["ϕ_flip"]           = vcat( merger_collection_STARS["ϕ_flip"], missing )
            merger_collection_STARS["δM_felix"]         = vcat( merger_collection_STARS["δM_felix"], missing )
            merger_collection_STARS["δM2_felix"]        = vcat( merger_collection_STARS["δM2_felix"], missing )
            merger_collection_STARS["δM_fromJ"]         = vcat( merger_collection_STARS["δM_fromJ"], missing )
            merger_collection_STARS["M_MM"]             = vcat( merger_collection_STARS["M_MM"], missing )
            merger_collection_STARS["M2_MM"]            = vcat( merger_collection_STARS["M2_MM"], missing )
            merger_collection_STARS["J_MMorbital"]      = hcat( merger_collection_STARS["J_MMorbital"], missings(Float64, 3) )
            merger_collection_STARS["J_SUMorbital"]     = hcat( merger_collection_STARS["J_SUMorbital"], missings(Float64, 3) )
            merger_collection_STARS["M_MERGERS"]        = vcat( merger_collection_STARS["M_MERGERS"], missing )
            merger_collection_STARS["M_MISSED"]         = vcat( merger_collection_STARS["M_MISSED"], missing )
            merger_collection_STARS["M_CONSIDERED"]     = vcat( merger_collection_STARS["M_CONSIDERED"], missing )
            merger_collection_STARS["M2_MERGERS"]       = vcat( merger_collection_STARS["M2_MERGERS"], missing )
            merger_collection_STARS["M2_MISSED"]        = vcat( merger_collection_STARS["M2_MISSED"], missing )
            merger_collection_STARS["M2_CONSIDERED"]    = vcat( merger_collection_STARS["M2_CONSIDERED"], missing )
        else
            merger_collection_STARS["δJ_main"]  = hcat( merger_collection_STARS["δJ_main"], halo_story["J_STARS"][:,fp_indices[i]] .- halo_story["J_STARS"][:,fp_indices[i-1]] )
            merger_collection_STARS["δj_main"]  = hcat( merger_collection_STARS["δj_main"], halo_story["j_STARS"][:,fp_indices[i]] .- halo_story["j_STARS"][:,fp_indices[i-1]] )
            merger_collection_STARS["δBVAL"]    = vcat( merger_collection_STARS["δBVAL"], halo_story["BVAL"][fp_indices[i]] - halo_story["BVAL"][fp_indices[i-1]] )
            if !ismissing(halo_story["j_STARS"][1,fp_indices[i]]) && !ismissing(halo_story["j_STARS"][2,fp_indices[i]]) && !ismissing(halo_story["j_STARS"][3,fp_indices[i]])
                merger_collection_STARS["ϕ_flip"]   = vcat( merger_collection_STARS["ϕ_flip"], acosd( transpose(halo_story["J_STARS"][:,fp_indices[i]]) * halo_story["J_STARS"][:,fp_indices[i-1]] / norm(halo_story["J_STARS"][:,fp_indices[i]]) / norm(halo_story["J_STARS"][:,fp_indices[i-1]]) ) )
                merger_collection_STARS["δM_fromJ"] = vcat( merger_collection_STARS["δM_fromJ"], (norm(halo_story["J_STARS"][:,fp_indices[i]]) / norm(halo_story["j_STARS"][:,fp_indices[i]])) - (norm(halo_story["J_STARS"][:,fp_indices[i-1]]) / norm(halo_story["j_STARS"][:,fp_indices[i-1]])) )
            else
                #println(halo_story["j_STARS"][:,fp_indices[i]])
                println("missing j_STARS for $(storyfilelist[iii]) @ index $(fp_indices[i])")
                merger_collection_STARS["ϕ_flip"]   = vcat( merger_collection_STARS["ϕ_flip"], missing )
                merger_collection_STARS["δM_fromJ"] = vcat( merger_collection_STARS["δM_fromJ"], missing )
            end
            #println("\n$(halo_story["J_STARS"][:,fp_indices[i]]) $(halo_story["J_STARS"][:,fp_indices[i-1]])")
            merger_collection_STARS["δM_felix"] = vcat( merger_collection_STARS["δM_felix"], convert_units_physical_mass(halo_story["M_STARS"][fp_indices[i]], head) - convert_units_physical_mass(halo_story["M_STARS"][fp_indices[i-1]], head) )
            merger_collection_STARS["δM2_felix"] = vcat( merger_collection_STARS["δM2_felix"], convert_units_physical_mass(halo_story["M_STAR_2"][fp_indices[i]], head) - convert_units_physical_mass(halo_story["M_STAR_2"][fp_indices[i-1]], head) )
            
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
            for ii in merger_indices[1+sum(n_mergers[1:i-1]):sum(n_mergers[1:i])]
                # Most Massive condition
                head    = read_header("$box/groups_$(@sprintf("%03i", halo_story["SNAP"][ii]))/sub_$(@sprintf("%03i", halo_story["SNAP"][ii]))")
                if convert_units_physical_mass(halo_story["M_STAR_2"][ii], head) > merger_collection_STARS["M2_MM"][end]
                    merger_collection_STARS["M_MM"][end]            = convert_units_physical_mass(halo_story["M_STARS"][ii], head)
                    merger_collection_STARS["M2_MM"][end]           = convert_units_physical_mass(halo_story["M_STAR_2"][ii], head)
                    merger_collection_STARS["J_MMorbital"][:,end]   = halo_story["J_orbital"][:,ii]
                end
                # Orbital Data
                if !ismissing(halo_story["J_orbital"][1,ii]) || !ismissing(halo_story["J_orbital"][2,ii]) || !ismissing(halo_story["J_orbital"][3,ii])
                    merger_collection_STARS["J_SUMorbital"][:,end] .+= halo_story["J_orbital"][:,ii]
                    merger_collection_STARS["M_MERGERS"][end]       += convert_units_physical_mass(halo_story["M_STARS"][ii], head)
                    merger_collection_STARS["M_CONSIDERED"][end]    += convert_units_physical_mass(halo_story["M_STARS"][ii], head)
                    merger_collection_STARS["M2_MERGERS"][end]      += convert_units_physical_mass(halo_story["M_STAR_2"][ii], head)
                    merger_collection_STARS["M2_CONSIDERED"][end]   += convert_units_physical_mass(halo_story["M_STAR_2"][ii], head)
                    #println(halo_story["J_orbital"][:,ii])
                else
                    #println("$ii")
                    merger_collection_STARS["M_MERGERS"][end]   += convert_units_physical_mass(halo_story["M_STARS"][ii], head)
                    merger_collection_STARS["M_MISSED"][end]    += convert_units_physical_mass(halo_story["M_STARS"][ii], head)
                    merger_collection_STARS["M2_MERGERS"][end]  += convert_units_physical_mass(halo_story["M_STAR_2"][ii], head)
                    merger_collection_STARS["M2_MISSED"][end]   += convert_units_physical_mass(halo_story["M_STAR_2"][ii], head)
                end
            end
        end
    end
    # Sanity check
    #println("$(size(merger_collection_STARS["SNAP"]))")
    #println("$(size(merger_collection_STARS["δJ_main"]))")
    #println("$(size(merger_collection_STARS["J_SUMorbital"]))")
    #println("$(size(merger_collection_STARS["M2_MERGERS"]))")
    #@show merger_collection_STARS["M2_MISSED"]

    save(joinpath(output_dir, "halo_$(merger_data["I_FILE"])_$(halo_story["I_SUB"][1]).jld"), 
        "halo_story",               halo_story,
        "merger_data",              merger_data,
        "merger_collection_STARS",  merger_collection_STARS)
    
end


println("\n\nNow witness the firepower of this fully armed and operational battle station!")
flush(stdout)