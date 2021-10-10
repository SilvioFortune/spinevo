# This version attributes mergers to the shared snap, so this includes snaps that cannot be read out

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

function ssFPfinder(mID, dict_fct) # Find same-snap First Progenitor depending on merger ID
    result_ID = 0
    for finder_ID in mID:-1:2    # walking the array towards the first progenitor
        if dict_fct["SNAP"][finder_ID] < dict_fct["SNAP"][finder_ID-1]
            result_ID   = finder_ID
            break
        end
    end
    return result_ID
end



### Settings
box         = "/HydroSims/Magneticum/Box4/uhr_test"
input_dir   = "/home/moon/sfortune/spinevo/halostories_update_stars"
min_time    = 0.9   # Gyr
max_time    = 1.5   # Gyr, not used yet
output_dir  = "/home/moon/sfortune/spinevo/halostories_v20211007_min$(min_time)Gyr"



### Main

storyfilelist   = readdir(input_dir)


limit_filelist  = length(storyfilelist)
println("$(limit_filelist) Runs in total.")
flush(stdout)
for iii in 1:limit_filelist
    println("$iii $(storyfilelist[iii])")
    flush(stdout)
    merger_data = load(joinpath(input_dir, storyfilelist[iii]), "merger_data")
    halo_story  = load(joinpath(input_dir, storyfilelist[iii]), "halo_story")

    head        = read_header("$box/groups_$(@sprintf("%03i", halo_story["SNAP"][1]))/sub_$(@sprintf("%03i", halo_story["SNAP"][1]))")
    head2        = read_header("$box/groups_$(@sprintf("%03i", halo_story["SNAP"][1]-3))/sub_$(@sprintf("%03i", halo_story["SNAP"][1]-3))")
    id_mfelix   = convert_units_physical_mass(halo_story["M_STARS"][1], head)
    id_m2       = convert_units_physical_mass(halo_story["M_STAR_2"][1], head2)
    
    
    # Identify First Progenitors and mergers

    loop_nmergers   = 0
    loop_mmergers   = 0.
    loop_J_sum      = zeros(3)
    mass_missed     = 0.
    mass_added      = 0.

    # Make lookbacktime array
    lbt = Array{Float64}(undef, 0)
    for i in 1:length(halo_story["SNAP"])
        head    = read_header("$box/groups_$(@sprintf("%03i", halo_story["SNAP"][i]))/sub_$(@sprintf("%03i", halo_story["SNAP"][i]))")
        lbt     = vcat( lbt, ustrip(lookback_time(cosmology(h=head.h0, OmegaM=head.omega_0), head.z)) )
    end

    # First Progenitors, backward in time
    fp_indices      = Array{Int64}(undef, 0)
    fp_snaps        = Array{Int64}(undef, 0)
    loop_length     = length(halo_story["SNAP"])
    for i in 1:loop_length
        if i==1 || halo_story["SNAP"][i-1] > halo_story["SNAP"][i] && count(ismissing, halo_story["J_STARS"][:,i]) == 0 # First Progenitor
            #print("FP@SNAP$(halo_story["SNAP"][i]) ")
            #flush(stdout)
            if length(fp_indices) == 0# && lbt[i]-lbt[1] > min_time
                fp_indices  = vcat( fp_indices, i )
                fp_snaps    = vcat( fp_snaps, halo_story["SNAP"][i] )
            elseif length(fp_indices) > 0 && lbt[i]-lbt[fp_indices[end]] > min_time
                fp_indices  = vcat( fp_indices, i )
                fp_snaps    = vcat( fp_snaps, halo_story["SNAP"][i] )
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
    #if fp_snaps[1] == halo_story["SNAP"][end]
        #fp_index    = 2
        #n_mergers   = vcat( n_mergers, merger_count )
    #else
        #fp_index    = 1
    #end
    # Merger Map, forward in time by starting at the bottom
    for i in loop_length:-1:1
        #print("$i ")
        if halo_story["SNAP"][i] ≥ fp_snaps[1]  # only accept entries starting with the same SNAP as the first FP
            # Check out Merger Info and assign to next
            if halo_story["SNAP"][i] == fp_snaps[fp_index] # First Progenitor
                n_mergers       = vcat( n_mergers, merger_count )
                merger_count    = 0
                fp_index       += 1
            end
        
            # Identify Mergers
            if i == 1 || halo_story["SNAP"][i] < halo_story["SNAP"][i-1] # FP
                print("")
            elseif halo_story["SNAP"][i] == halo_story["SNAP"][i-1] # Merger
                merger_count     += 1
                merger_index_map  = hcat( merger_index_map, [i, ssFPfinder(i, halo_story)] )
            else
                println("Error for i = $i, $(halo_story["SNAP"][i]), $(fp_snaps[fp_index]), $(halo_story["SNAP"][i-1]), $(halo_story["I_SUB"][i-1])")
            end
        end
    end
    
    #println("$(length(fp_indices)) $(length(fp_snaps)) $(length(n_mergers))")
    #println("$(length(merger_index_map)) $(sum(n_mergers)) ")
    
    
    # Fill the STARS dictionary

    merger_collection_STARS = Dict(
            "N_MERGERS"     => n_mergers,
            "SNAP"          => missings(Int64   , 0),
            "ID_ISUB"       => missings(Int64   , 0),
            "I_SUB"         => missings(Int64   , 0),
            "ID_Mfelix"     => missings(Float64 , 0),
            "ID_M2"         => missings(Float64 , 0),
            "REDSHIFT"      => missings(Float64 , 0),
            "LOOKBACKTIME"  => missings(Float64 , 0),
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
            "M_MM"          => missings(Float64 , 2, 0), # 1 = mass, 2 = FPmass
            "M2_MM"         => missings(Float64 , 2, 0), # 1 = mass, 2 = FPmass
            "J_MMorbital"   => missings(Float64 , 3, 0), 
            "J_SUMorbital"  => missings(Float64 , 3, 0), 
            "δJ_main"       => missings(Float64 , 3, 0), 
            "J_main"        => missings(Float64 , 3, 0), 
            "j_main"        => missings(Float64 , 3, 0), 
            "δj_main"       => missings(Float64 , 3, 0),
            "Merger_Map"    => missings(Float64 , 7, 0)) # 1=mass, 2=mass2, 3=FPmass, 4=FPmass2, 5=SNAP, 6=z, 7=lbt

    for i in 1:length(fp_indices)
        #print("$i ")
        #flush(stdout)

        # Basic Info
        head    = read_header("$box/groups_$(@sprintf("%03i", halo_story["SNAP"][fp_indices[i]]))/sub_$(@sprintf("%03i", halo_story["SNAP"][fp_indices[i]]))")
        head2    = read_header("$box/groups_$(@sprintf("%03i", halo_story["SNAP"][fp_indices[i]]-3))/sub_$(@sprintf("%03i", halo_story["SNAP"][fp_indices[i]]-3))")
        merger_collection_STARS["LOOKBACKTIME"] = vcat( merger_collection_STARS["LOOKBACKTIME"], ustrip(lookback_time(cosmology(h=head.h0, OmegaM=head.omega_0), head.z)) )
        merger_collection_STARS["SNAP"]         = vcat( merger_collection_STARS["SNAP"], halo_story["SNAP"][fp_indices[i]] )
        merger_collection_STARS["ID_ISUB"]      = vcat( merger_collection_STARS["ID_ISUB"], halo_story["I_SUB"][1] )
        merger_collection_STARS["ID_Mfelix"]    = vcat( merger_collection_STARS["ID_Mfelix"], id_mfelix )
        merger_collection_STARS["ID_M2"]        = vcat( merger_collection_STARS["ID_M2"], id_m2 )
        merger_collection_STARS["REDSHIFT"]     = vcat( merger_collection_STARS["REDSHIFT"], halo_story["REDSHIFT"][fp_indices[i]] )
        merger_collection_STARS["I_SUB"]        = vcat( merger_collection_STARS["I_SUB"], halo_story["I_SUB"][fp_indices[i]] )
        merger_collection_STARS["M_felix"]      = vcat( merger_collection_STARS["M_felix"], convert_units_physical_mass(halo_story["M_STARS"][fp_indices[i]], head) )
        merger_collection_STARS["M2_felix"]     = vcat( merger_collection_STARS["M2_felix"], convert_units_physical_mass(halo_story["M_STAR_2"][fp_indices[i]], head2) )
        if count(ismissing, halo_story["j_STARS"][:,fp_indices[i]]) == 0 
            merger_collection_STARS["M_fromJ"]  = vcat( merger_collection_STARS["M_fromJ"], norm(halo_story["J_STARS"][:,fp_indices[i]]) / norm(halo_story["j_STARS"][:,fp_indices[i]]) )
        else
            merger_collection_STARS["M_fromJ"]  = vcat( merger_collection_STARS["M_fromJ"], missing )
        end
        merger_collection_STARS["BVAL"]         = vcat( merger_collection_STARS["BVAL"], halo_story["BVAL"][fp_indices[i]] )
        merger_collection_STARS["J_main"]       = hcat( merger_collection_STARS["J_main"], halo_story["J_STARS"][:,fp_indices[i]] )
        merger_collection_STARS["j_main"]       = hcat( merger_collection_STARS["j_main"], halo_story["j_STARS"][:,fp_indices[i]] )
        
        # Transitional Data
        if i == 1
            merger_collection_STARS["δJ_main"]          = hcat( merger_collection_STARS["δJ_main"], missings(Float64, 3) )
            merger_collection_STARS["δj_main"]          = hcat( merger_collection_STARS["δj_main"], missings(Float64, 3) )
            merger_collection_STARS["δBVAL"]            = vcat( merger_collection_STARS["δBVAL"], missing )
            merger_collection_STARS["ϕ_flip"]           = vcat( merger_collection_STARS["ϕ_flip"], missing )
            merger_collection_STARS["δM_felix"]         = vcat( merger_collection_STARS["δM_felix"], missing )
            merger_collection_STARS["δM2_felix"]        = vcat( merger_collection_STARS["δM2_felix"], missing )
            merger_collection_STARS["δM_fromJ"]         = vcat( merger_collection_STARS["δM_fromJ"], missing )
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
            merger_collection_STARS["δJ_main"]  = hcat( merger_collection_STARS["δJ_main"], halo_story["J_STARS"][:,fp_indices[i]] .- halo_story["J_STARS"][:,fp_indices[i-1]] )
            merger_collection_STARS["δj_main"]  = hcat( merger_collection_STARS["δj_main"], halo_story["j_STARS"][:,fp_indices[i]] .- halo_story["j_STARS"][:,fp_indices[i-1]] )
            merger_collection_STARS["δBVAL"]    = vcat( merger_collection_STARS["δBVAL"], halo_story["BVAL"][fp_indices[i]] - halo_story["BVAL"][fp_indices[i-1]] )
            if count(ismissing, halo_story["j_STARS"][:,fp_indices[i]]) == 0 && count(ismissing, halo_story["j_STARS"][:,fp_indices[i-1]]) == 0 
                merger_collection_STARS["ϕ_flip"]   = vcat( merger_collection_STARS["ϕ_flip"], acosd( transpose(halo_story["J_STARS"][:,fp_indices[i]]) * halo_story["J_STARS"][:,fp_indices[i-1]] / norm(halo_story["J_STARS"][:,fp_indices[i]]) / norm(halo_story["J_STARS"][:,fp_indices[i-1]]) ) )
                merger_collection_STARS["δM_fromJ"] = vcat( merger_collection_STARS["δM_fromJ"], (norm(halo_story["J_STARS"][:,fp_indices[i]]) / norm(halo_story["j_STARS"][:,fp_indices[i]])) - (norm(halo_story["J_STARS"][:,fp_indices[i-1]]) / norm(halo_story["j_STARS"][:,fp_indices[i-1]])) )
            else
                #println(halo_story["j_STARS"][:,fp_indices[i]])
                #println("missing j_STARS for $(storyfilelist[iii]) @ index $(fp_indices[i])")
                merger_collection_STARS["ϕ_flip"]   = vcat( merger_collection_STARS["ϕ_flip"], missing )
                merger_collection_STARS["δM_fromJ"] = vcat( merger_collection_STARS["δM_fromJ"], missing )
            end
            #println("\n$(halo_story["J_STARS"][:,fp_indices[i]]) $(halo_story["J_STARS"][:,fp_indices[i-1]])")
            merger_collection_STARS["δM_felix"] = vcat( merger_collection_STARS["δM_felix"], convert_units_physical_mass(halo_story["M_STARS"][fp_indices[i]], head) - convert_units_physical_mass(halo_story["M_STARS"][fp_indices[i-1]], head) )
            merger_collection_STARS["δM2_felix"] = vcat( merger_collection_STARS["δM2_felix"], convert_units_physical_mass(halo_story["M_STAR_2"][fp_indices[i]], head2) - convert_units_physical_mass(halo_story["M_STAR_2"][fp_indices[i-1]], head2) )
            
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
                head    = read_header("$box/groups_$(@sprintf("%03i", halo_story["SNAP"][ii]))/sub_$(@sprintf("%03i", halo_story["SNAP"][ii]))")
                head2   = read_header("$box/groups_$(@sprintf("%03i", halo_story["SNAP"][ii]-3))/sub_$(@sprintf("%03i", halo_story["SNAP"][ii]-3))")
                # Merger Map
                #println("$(halo_story["SNAP"][ii])   i = $i")
                #@show merger_index_map
                #@show merger_index_map[1,1+sum(n_mergers[1:i-1]):sum(n_mergers[1:i])]
                #flush(stdout)
                merger_collection_STARS["Merger_Map"]   = hcat( merger_collection_STARS["Merger_Map"], 
                                                                    [   convert_units_physical_mass(halo_story["M_STARS"][ii], head),
                                                                        convert_units_physical_mass(halo_story["M_STAR_2"][ii], head2), 
                                                                        convert_units_physical_mass(halo_story["M_STARS"][merger_index_map[2,pos_in_mim]], head), 
                                                                        convert_units_physical_mass(halo_story["M_STAR_2"][merger_index_map[2,pos_in_mim]], head2), 
                                                                        halo_story["SNAP"][ii], 
                                                                        halo_story["REDSHIFT"][ii], 
                                                                        ustrip(lookback_time(cosmology(h=head.h0, OmegaM=head.omega_0), head.z)) , ] )
                # Most Massive condition
                if convert_units_physical_mass(halo_story["M_STAR_2"][ii], head2) > merger_collection_STARS["M2_MM"][end]
                    merger_collection_STARS["M_MM"][end]            = convert_units_physical_mass(halo_story["M_STARS"][ii], head)
                    merger_collection_STARS["M2_MM"][end]           = convert_units_physical_mass(halo_story["M_STAR_2"][ii], head2)
                    merger_collection_STARS["J_MMorbital"][:,end]   = halo_story["j_orbital"][:,ii] .* convert_units_physical_mass(halo_story["M_STAR_2"][ii], head2)
                end
                # Orbital Data
                if count(ismissing, halo_story["j_orbital"][:,ii]) == 0
                    merger_collection_STARS["J_SUMorbital"][:,end] .+= ( halo_story["j_orbital"][:,ii] .* convert_units_physical_mass(halo_story["M_STAR_2"][ii], head2) )
                    merger_collection_STARS["M_MERGERS"][end]       += convert_units_physical_mass(halo_story["M_STARS"][ii], head)
                    merger_collection_STARS["M_CONSIDERED"][end]    += convert_units_physical_mass(halo_story["M_STARS"][ii], head)
                    merger_collection_STARS["M2_MERGERS"][end]      += convert_units_physical_mass(halo_story["M_STAR_2"][ii], head2)
                    merger_collection_STARS["M2_CONSIDERED"][end]   += convert_units_physical_mass(halo_story["M_STAR_2"][ii], head2)
                    #println(halo_story["j_orbital"][:,ii])
                else
                    #println("$ii")
                    merger_collection_STARS["M_MERGERS"][end]   += convert_units_physical_mass(halo_story["M_STARS"][ii], head)
                    merger_collection_STARS["M_MISSED"][end]    += convert_units_physical_mass(halo_story["M_STARS"][ii], head)
                    merger_collection_STARS["M2_MERGERS"][end]  += convert_units_physical_mass(halo_story["M_STAR_2"][ii], head2)
                    merger_collection_STARS["M2_MISSED"][end]   += convert_units_physical_mass(halo_story["M_STAR_2"][ii], head2)
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
            "SNAP"          => missings(Int64   , 0),
            "ID_ISUB"       => missings(Int64   , 0),
            "I_SUB"         => missings(Int64   , 0),
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
            "δj_main"       => missings(Float64 , 3, 0),
            "Merger_Map"    => missings(Float64 , 7, 0)) # 1=mass, 2=mass2, 3=FPmass, 4=FPmass2, 5=SNAP, 6=z, 7=lbt

    for i in 1:length(fp_indices)
        #print("$i ")
        #flush(stdout)

        # Basic Info
        head    = read_header("$box/groups_$(@sprintf("%03i", halo_story["SNAP"][fp_indices[i]]))/sub_$(@sprintf("%03i", halo_story["SNAP"][fp_indices[i]]))")
        head2    = read_header("$box/groups_$(@sprintf("%03i", halo_story["SNAP"][fp_indices[i]]-3))/sub_$(@sprintf("%03i", halo_story["SNAP"][fp_indices[i]]-3))")
        merger_collection_DM["LOOKBACKTIME"] = vcat( merger_collection_DM["LOOKBACKTIME"], ustrip(lookback_time(cosmology(h=head.h0, OmegaM=head.omega_0), head.z)) )
        merger_collection_DM["SNAP"]         = vcat( merger_collection_DM["SNAP"], halo_story["SNAP"][fp_indices[i]] )
        merger_collection_DM["ID_ISUB"]      = vcat( merger_collection_DM["ID_ISUB"], halo_story["I_SUB"][1] )
        merger_collection_DM["ID_Mfelix"]    = vcat( merger_collection_DM["ID_Mfelix"], id_mfelix )
        merger_collection_DM["ID_M2"]        = vcat( merger_collection_DM["ID_M2"], id_m2 )
        merger_collection_DM["REDSHIFT"]     = vcat( merger_collection_DM["REDSHIFT"], halo_story["REDSHIFT"][fp_indices[i]] )
        merger_collection_DM["I_SUB"]        = vcat( merger_collection_DM["I_SUB"], halo_story["I_SUB"][fp_indices[i]] )
        merger_collection_DM["M_felix"]      = vcat( merger_collection_DM["M_felix"], convert_units_physical_mass(halo_story["M_DM"][fp_indices[i]], head) )
        merger_collection_DM["M2_felix"]     = vcat( merger_collection_DM["M2_felix"], convert_units_physical_mass(halo_story["M_DM_2"][fp_indices[i]], head2) )
        if count(ismissing, halo_story["j_DM"][:,fp_indices[i]]) == 0
            merger_collection_DM["M_fromJ"]  = vcat( merger_collection_DM["M_fromJ"], norm(halo_story["J_DM"][:,fp_indices[i]]) / norm(halo_story["j_DM"][:,fp_indices[i]]) )
        else
            merger_collection_DM["M_fromJ"]  = vcat( merger_collection_DM["M_fromJ"], missing )
        end
        merger_collection_DM["BVAL"]         = vcat( merger_collection_DM["BVAL"], halo_story["BVAL"][fp_indices[i]] )
        merger_collection_DM["J_main"]       = hcat( merger_collection_DM["J_main"], halo_story["J_DM"][:,fp_indices[i]] )
        merger_collection_DM["j_main"]       = hcat( merger_collection_DM["j_main"], halo_story["j_DM"][:,fp_indices[i]] )
        
        # Transitional Data
        if i == 1
            merger_collection_DM["δJ_main"]          = hcat( merger_collection_DM["δJ_main"], missings(Float64, 3) )
            merger_collection_DM["δj_main"]          = hcat( merger_collection_DM["δj_main"], missings(Float64, 3) )
            merger_collection_DM["δBVAL"]            = vcat( merger_collection_DM["δBVAL"], missing )
            merger_collection_DM["ϕ_flip"]           = vcat( merger_collection_DM["ϕ_flip"], missing )
            merger_collection_DM["δM_felix"]         = vcat( merger_collection_DM["δM_felix"], missing )
            merger_collection_DM["δM2_felix"]        = vcat( merger_collection_DM["δM2_felix"], missing )
            merger_collection_DM["δM_fromJ"]         = vcat( merger_collection_DM["δM_fromJ"], missing )
            merger_collection_DM["M_MM"]             = vcat( merger_collection_DM["M_MM"], missing )
            merger_collection_DM["M2_MM"]            = vcat( merger_collection_DM["M2_MM"], missing )
            merger_collection_DM["J_MMorbital"]      = hcat( merger_collection_DM["J_MMorbital"], missings(Float64, 3) )
            merger_collection_DM["J_SUMorbital"]     = hcat( merger_collection_DM["J_SUMorbital"], missings(Float64, 3) )
            merger_collection_DM["M_MERGERS"]        = vcat( merger_collection_DM["M_MERGERS"], missing )
            merger_collection_DM["M_MISSED"]         = vcat( merger_collection_DM["M_MISSED"], missing )
            merger_collection_DM["M_CONSIDERED"]     = vcat( merger_collection_DM["M_CONSIDERED"], missing )
            merger_collection_DM["M2_MERGERS"]       = vcat( merger_collection_DM["M2_MERGERS"], missing )
            merger_collection_DM["M2_MISSED"]        = vcat( merger_collection_DM["M2_MISSED"], missing )
            merger_collection_DM["M2_CONSIDERED"]    = vcat( merger_collection_DM["M2_CONSIDERED"], missing )
        else
            merger_collection_DM["δJ_main"]  = hcat( merger_collection_DM["δJ_main"], halo_story["J_DM"][:,fp_indices[i]] .- halo_story["J_DM"][:,fp_indices[i-1]] )
            merger_collection_DM["δj_main"]  = hcat( merger_collection_DM["δj_main"], halo_story["j_DM"][:,fp_indices[i]] .- halo_story["j_DM"][:,fp_indices[i-1]] )
            merger_collection_DM["δBVAL"]    = vcat( merger_collection_DM["δBVAL"], halo_story["BVAL"][fp_indices[i]] - halo_story["BVAL"][fp_indices[i-1]] )
            if count(ismissing, halo_story["j_DM"][:,fp_indices[i]]) == 0 && count(ismissing, halo_story["j_DM"][:,fp_indices[i-1]]) == 0
                #println("checking error for index $(fp_indices[i])   ---   $(halo_story["J_DM"][:,fp_indices[i]])   ---   $(halo_story["J_DM"][:,fp_indices[i-1]])   ---   $(halo_story["j_DM"][:,fp_indices[i]])   ---   $(halo_story["j_DM"][:,fp_indices[i-1]])")
                merger_collection_DM["ϕ_flip"]   = vcat( merger_collection_DM["ϕ_flip"], acosd( transpose(halo_story["J_DM"][:,fp_indices[i]]) * halo_story["J_DM"][:,fp_indices[i-1]] / norm(halo_story["J_DM"][:,fp_indices[i]]) / norm(halo_story["J_DM"][:,fp_indices[i-1]]) ) )
                merger_collection_DM["δM_fromJ"] = vcat( merger_collection_DM["δM_fromJ"], (norm(halo_story["J_DM"][:,fp_indices[i]]) / norm(halo_story["j_DM"][:,fp_indices[i]])) - (norm(halo_story["J_DM"][:,fp_indices[i-1]]) / norm(halo_story["j_DM"][:,fp_indices[i-1]])) )
            else
                #println(halo_story["j_DM"][:,fp_indices[i]])
                #println("missing j_DM for $(storyfilelist[iii]) @ index $(fp_indices[i])")
                merger_collection_DM["ϕ_flip"]   = vcat( merger_collection_DM["ϕ_flip"], missing )
                merger_collection_DM["δM_fromJ"] = vcat( merger_collection_DM["δM_fromJ"], missing )
            end
            #println("\n$(halo_story["J_DM"][:,fp_indices[i]]) $(halo_story["J_DM"][:,fp_indices[i-1]])")
            merger_collection_DM["δM_felix"] = vcat( merger_collection_DM["δM_felix"], convert_units_physical_mass(halo_story["M_DM"][fp_indices[i]], head) - convert_units_physical_mass(halo_story["M_DM"][fp_indices[i-1]], head) )
            merger_collection_DM["δM2_felix"] = vcat( merger_collection_DM["δM2_felix"], convert_units_physical_mass(halo_story["M_DM_2"][fp_indices[i]], head2) - convert_units_physical_mass(halo_story["M_DM_2"][fp_indices[i-1]], head2) )
            
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
            pos_in_mim = 1+sum(n_mergers[1:i-1])
            for ii in merger_index_map[1,1+sum(n_mergers[1:i-1]):sum(n_mergers[1:i])]
                head    = read_header("$box/groups_$(@sprintf("%03i", halo_story["SNAP"][ii]))/sub_$(@sprintf("%03i", halo_story["SNAP"][ii]))")
                head2   = read_header("$box/groups_$(@sprintf("%03i", halo_story["SNAP"][ii]-3))/sub_$(@sprintf("%03i", halo_story["SNAP"][ii]-3))")
                # Merger Map
                merger_collection_DM["Merger_Map"]   = hcat( merger_collection_DM["Merger_Map"], 
                                                                    [   convert_units_physical_mass(halo_story["M_DM"][ii], head),
                                                                        convert_units_physical_mass(halo_story["M_DM_2"][ii], head2), 
                                                                        convert_units_physical_mass(halo_story["M_DM"][merger_index_map[2,pos_in_mim]], head), 
                                                                        convert_units_physical_mass(halo_story["M_DM_2"][merger_index_map[2,pos_in_mim]], head2), 
                                                                        halo_story["SNAP"][ii], 
                                                                        halo_story["REDSHIFT"][ii], 
                                                                        ustrip(lookback_time(cosmology(h=head.h0, OmegaM=head.omega_0), head.z)) , ] )
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
            "SNAP"          => missings(Int64   , 0),
            "ID_ISUB"       => missings(Int64   , 0),
            "I_SUB"         => missings(Int64   , 0),
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
            "δj_main"       => missings(Float64 , 3, 0),
            "Merger_Map"    => missings(Float64 , 7, 0)) # 1=mass, 2=mass2, 3=FPmass, 4=FPmass2, 5=SNAP, 6=z, 7=lbt

    for i in 1:length(fp_indices)
        #print("$i ")
        #flush(stdout)

        # Basic Info
        head    = read_header("$box/groups_$(@sprintf("%03i", halo_story["SNAP"][fp_indices[i]]))/sub_$(@sprintf("%03i", halo_story["SNAP"][fp_indices[i]]))")
        head2    = read_header("$box/groups_$(@sprintf("%03i", halo_story["SNAP"][fp_indices[i]]-3))/sub_$(@sprintf("%03i", halo_story["SNAP"][fp_indices[i]]-3))")
        merger_collection_GAS["LOOKBACKTIME"] = vcat( merger_collection_GAS["LOOKBACKTIME"], ustrip(lookback_time(cosmology(h=head.h0, OmegaM=head.omega_0), head.z)) )
        merger_collection_GAS["SNAP"]         = vcat( merger_collection_GAS["SNAP"], halo_story["SNAP"][fp_indices[i]] )
        merger_collection_GAS["ID_ISUB"]      = vcat( merger_collection_GAS["ID_ISUB"], halo_story["I_SUB"][1] )
        merger_collection_GAS["ID_Mfelix"]    = vcat( merger_collection_GAS["ID_Mfelix"], id_mfelix )
        merger_collection_GAS["ID_M2"]        = vcat( merger_collection_GAS["ID_M2"], id_m2 )
        merger_collection_GAS["REDSHIFT"]     = vcat( merger_collection_GAS["REDSHIFT"], halo_story["REDSHIFT"][fp_indices[i]] )
        merger_collection_GAS["I_SUB"]        = vcat( merger_collection_GAS["I_SUB"], halo_story["I_SUB"][fp_indices[i]] )
        merger_collection_GAS["M_felix"]      = vcat( merger_collection_GAS["M_felix"], convert_units_physical_mass(halo_story["M_GAS"][fp_indices[i]], head) )
        merger_collection_GAS["M2_felix"]     = vcat( merger_collection_GAS["M2_felix"], convert_units_physical_mass(halo_story["M_GAS_2"][fp_indices[i]], head2) )
        if count(ismissing, halo_story["j_GAS"][:,fp_indices[i]]) == 0
            merger_collection_GAS["M_fromJ"]  = vcat( merger_collection_GAS["M_fromJ"], norm(halo_story["J_GAS"][:,fp_indices[i]]) / norm(halo_story["j_GAS"][:,fp_indices[i]]) )
        else
            merger_collection_GAS["M_fromJ"]  = vcat( merger_collection_GAS["M_fromJ"], missing )
        end
        merger_collection_GAS["BVAL"]         = vcat( merger_collection_GAS["BVAL"], halo_story["BVAL"][fp_indices[i]] )
        merger_collection_GAS["J_main"]       = hcat( merger_collection_GAS["J_main"], halo_story["J_GAS"][:,fp_indices[i]] )
        merger_collection_GAS["j_main"]       = hcat( merger_collection_GAS["j_main"], halo_story["j_GAS"][:,fp_indices[i]] )
        
        # Transitional Data
        if i == 1
            merger_collection_GAS["δJ_main"]          = hcat( merger_collection_GAS["δJ_main"], missings(Float64, 3) )
            merger_collection_GAS["δj_main"]          = hcat( merger_collection_GAS["δj_main"], missings(Float64, 3) )
            merger_collection_GAS["δBVAL"]            = vcat( merger_collection_GAS["δBVAL"], missing )
            merger_collection_GAS["ϕ_flip"]           = vcat( merger_collection_GAS["ϕ_flip"], missing )
            merger_collection_GAS["δM_felix"]         = vcat( merger_collection_GAS["δM_felix"], missing )
            merger_collection_GAS["δM2_felix"]        = vcat( merger_collection_GAS["δM2_felix"], missing )
            merger_collection_GAS["δM_fromJ"]         = vcat( merger_collection_GAS["δM_fromJ"], missing )
            merger_collection_GAS["M_MM"]             = vcat( merger_collection_GAS["M_MM"], missing )
            merger_collection_GAS["M2_MM"]            = vcat( merger_collection_GAS["M2_MM"], missing )
            merger_collection_GAS["J_MMorbital"]      = hcat( merger_collection_GAS["J_MMorbital"], missings(Float64, 3) )
            merger_collection_GAS["J_SUMorbital"]     = hcat( merger_collection_GAS["J_SUMorbital"], missings(Float64, 3) )
            merger_collection_GAS["M_MERGERS"]        = vcat( merger_collection_GAS["M_MERGERS"], missing )
            merger_collection_GAS["M_MISSED"]         = vcat( merger_collection_GAS["M_MISSED"], missing )
            merger_collection_GAS["M_CONSIDERED"]     = vcat( merger_collection_GAS["M_CONSIDERED"], missing )
            merger_collection_GAS["M2_MERGERS"]       = vcat( merger_collection_GAS["M2_MERGERS"], missing )
            merger_collection_GAS["M2_MISSED"]        = vcat( merger_collection_GAS["M2_MISSED"], missing )
            merger_collection_GAS["M2_CONSIDERED"]    = vcat( merger_collection_GAS["M2_CONSIDERED"], missing )
        else
            merger_collection_GAS["δJ_main"]  = hcat( merger_collection_GAS["δJ_main"], halo_story["J_GAS"][:,fp_indices[i]] .- halo_story["J_GAS"][:,fp_indices[i-1]] )
            merger_collection_GAS["δj_main"]  = hcat( merger_collection_GAS["δj_main"], halo_story["j_GAS"][:,fp_indices[i]] .- halo_story["j_GAS"][:,fp_indices[i-1]] )
            merger_collection_GAS["δBVAL"]    = vcat( merger_collection_GAS["δBVAL"], halo_story["BVAL"][fp_indices[i]] - halo_story["BVAL"][fp_indices[i-1]] )
            if count(ismissing, halo_story["j_GAS"][:,fp_indices[i]]) == 0 && count(ismissing, halo_story["j_GAS"][:,fp_indices[i-1]]) == 0
                merger_collection_GAS["ϕ_flip"]   = vcat( merger_collection_GAS["ϕ_flip"], acosd( transpose(halo_story["J_GAS"][:,fp_indices[i]]) * halo_story["J_GAS"][:,fp_indices[i-1]] / norm(halo_story["J_GAS"][:,fp_indices[i]]) / norm(halo_story["J_GAS"][:,fp_indices[i-1]]) ) )
                merger_collection_GAS["δM_fromJ"] = vcat( merger_collection_GAS["δM_fromJ"], (norm(halo_story["J_GAS"][:,fp_indices[i]]) / norm(halo_story["j_GAS"][:,fp_indices[i]])) - (norm(halo_story["J_GAS"][:,fp_indices[i-1]]) / norm(halo_story["j_GAS"][:,fp_indices[i-1]])) )
            else
                #println(halo_story["j_GAS"][:,fp_indices[i]])
                #println("missing j_GAS for $(storyfilelist[iii]) @ index $(fp_indices[i])")
                merger_collection_GAS["ϕ_flip"]   = vcat( merger_collection_GAS["ϕ_flip"], missing )
                merger_collection_GAS["δM_fromJ"] = vcat( merger_collection_GAS["δM_fromJ"], missing )
            end
            #println("\n$(halo_story["J_GAS"][:,fp_indices[i]]) $(halo_story["J_GAS"][:,fp_indices[i-1]])")
            merger_collection_GAS["δM_felix"] = vcat( merger_collection_GAS["δM_felix"], convert_units_physical_mass(halo_story["M_GAS"][fp_indices[i]], head) - convert_units_physical_mass(halo_story["M_GAS"][fp_indices[i-1]], head) )
            merger_collection_GAS["δM2_felix"] = vcat( merger_collection_GAS["δM2_felix"], convert_units_physical_mass(halo_story["M_GAS_2"][fp_indices[i]], head2) - convert_units_physical_mass(halo_story["M_GAS_2"][fp_indices[i-1]], head2) )
            
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
            pos_in_mim = 1+sum(n_mergers[1:i-1])
            for ii in merger_index_map[1,1+sum(n_mergers[1:i-1]):sum(n_mergers[1:i])]
                head    = read_header("$box/groups_$(@sprintf("%03i", halo_story["SNAP"][ii]))/sub_$(@sprintf("%03i", halo_story["SNAP"][ii]))")
                head2    = read_header("$box/groups_$(@sprintf("%03i", halo_story["SNAP"][ii]-3))/sub_$(@sprintf("%03i", halo_story["SNAP"][ii]-3))")
                # Merger Map
                merger_collection_GAS["Merger_Map"]   = hcat( merger_collection_GAS["Merger_Map"], 
                                                                    [   convert_units_physical_mass(halo_story["M_GAS"][ii], head),
                                                                        convert_units_physical_mass(halo_story["M_GAS_2"][ii], head2), 
                                                                        convert_units_physical_mass(halo_story["M_GAS"][merger_index_map[2,pos_in_mim]], head), 
                                                                        convert_units_physical_mass(halo_story["M_GAS_2"][merger_index_map[2,pos_in_mim]], head2), 
                                                                        halo_story["SNAP"][ii], 
                                                                        halo_story["REDSHIFT"][ii], 
                                                                        ustrip(lookback_time(cosmology(h=head.h0, OmegaM=head.omega_0), head.z)) , ] )
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



    save(joinpath(output_dir, "halo_$(merger_data["I_FILE"])_$(halo_story["I_SUB"][1]).jld"), 
        "halo_story",               halo_story,
        "merger_data",              merger_data,
        "merger_collection_DM",     merger_collection_DM,
        "merger_collection_GAS",    merger_collection_GAS,
        "merger_collection_STARS",  merger_collection_STARS)
    
end


println("\n\nNow witness the firepower of this fully armed and operational battle station!")
flush(stdout)