### Basics

halofiles   = readdir("/home/moon/fschulze/MA/Data/thesis/Halo_trace_20200605_new/")
box         = "/HydroSims/Magneticum/Box4/uhr_test/"



### Packages

using Printf
using LinearAlgebra
using JLD
using QuadGK
using Statistics
using PyCall
using LaTeXStrings
using GadgetIO
using GadgetUnits
using GadgetGalaxies
using UnitfulAstro
using Missings
using Distributions
using CSV
using DataFrames



### Functions

# Calculate the orbital angular momentum with respect to the central subhalo
function orbit_J(subID, centID, pathtofile) # Find central sub maybe using SUBFIND -> centID = FSUB[GRNR[subID+1]+1]
    # Read subhalo features
    head        = read_header(pathtofile)
    smst        = convert_units_physical(read_subfind(pathtofile, "SMST"), :mass, head)
    spos        = convert_units_physical(read_subfind(pathtofile, "SPOS"), :pos, head)
    svel        = read_subfind(pathtofile, "SVEL") # no conversion since already physical units??
    # Return with masses added
    return ((spos[:,subID+1] .- spos[:,centID+1]) × (svel[:,subID+1] .- svel[:,centID+1])) .* (smst[1,subID+1]+smst[2,subID+1]+smst[5,subID+1])
end


### Main

println("\n\nInitiating main loop.\n")
for ii in length(halofiles):-1:1200
    filepath = string("/home/moon/fschulze/MA/Data/thesis/Halo_trace_20200605_new/", halofiles[ii])
    treefile_df = CSV.read(filepath, DataFrame; delim=' ', ignorerepeated=true, header=2)
    halo_story = Dict(
        "ID"        => parse(Int64, chop(halofiles[ii], head=5, tail=4)), 
        "BOX"       => box, 
        "SNAP"      => treefile_df[:, :SNAP], 
        "I_SUB"     => treefile_df[:, :I_SUB], 
        "I_TREE"    => treefile_df[:, :I_TREE], 
        "FILE_NR"   => treefile_df[:, :FILE_NR], 
        "REDSHIFT"  => treefile_df[:, :REDSHIFT], 
        "M_STARS"   => treefile_df[:, :M_STARS], 
        "M_STAR_2"  => treefile_df[:, :M_STAR_2], 
        "M_GAS"     => treefile_df[:, :M_GAS], 
        "M_GAS_2"   => treefile_df[:, :M_GAS_2], 
        "M_DM"      => treefile_df[:, :M_DM], 
        "M_DM_2"    => treefile_df[:, :M_DM_2], 
        "MMP"       => treefile_df[:, :MMP], 
        "RVIR"          => missings(Float64, length(treefile_df[:, :I_SUB])),
        "RHMS_DM"       => missings(Float64, length(treefile_df[:, :I_SUB])),
        "J_DM"          => missings(Float64, 3, length(treefile_df[:, :I_SUB])),
        "j_DM"          => missings(Float64, 3, length(treefile_df[:, :I_SUB])),
        "RHMS_GAS"      => missings(Float64, length(treefile_df[:, :I_SUB])),
        "J_GAS"         => missings(Float64, 3, length(treefile_df[:, :I_SUB])),
        "j_GAS"         => missings(Float64, 3, length(treefile_df[:, :I_SUB])),
        "RHMS_STARS"    => missings(Float64, length(treefile_df[:, :I_SUB])),
        "J_STARS"       => missings(Float64, 3, length(treefile_df[:, :I_SUB])),
        "j_STARS"       => missings(Float64, 3, length(treefile_df[:, :I_SUB])),
        "J_orbital"     => missings(Float32, 3, length(treefile_df[:, :I_SUB])),
        "BVAL"          => missings(Float64, length(treefile_df[:, :I_SUB])))#,
        #"ERRMSG"                => Vector{String}(undef, length(treefile_df[:, :I_SUB])))

    # Final Subhalo First
    snapshot                = Snapshot(box, halo_story["SNAP"][1])
    snapshot.snapbase
    snapshot.subbase
    g                       = Galaxy(snapshot, halo_story["I_SUB"][1])
    halo_story["RVIR"][1]   = read_galaxy_prop(get_group(g), "RVIR", :physical)
    sph_small               = Sphere(0.1*halo_story["RVIR"][1])
    sph_large               = Sphere(halo_story["RVIR"][1])
    try
        read_halo!(g, units=:physical, props=((:gas, ["POS", "VEL", "MASS"]), (:stars, ["POS", "VEL", "MASS"])), radius_units=:physical, radius=0.1*halo_story["RVIR"][1]) 
        read_halo!(g, units=:physical, props=((:dm, ["POS", "VEL"]),), radius_units=:physical, radius=halo_story["RVIR"][1])       
        if GadgetGalaxies.total_mass(g.dm.pos, g.dm.mass) > 0.0
            halo_story["RHMS_DM"][1]    = half_mass_radius(g.dm)
            halo_story["J_DM"][:,1]     = angular_momentum(g.dm, sph_large)
            halo_story["j_DM"][:,1]     = specific_angular_momentum(g.dm, sph_large)
        end
        if GadgetGalaxies.total_mass(g.gas.pos, g.gas.mass) > 0.0
            halo_story["RHMS_GAS"][1]   = half_mass_radius(g.gas)
            halo_story["J_GAS"][:,1]    = angular_momentum(g.gas, sph_small)
            halo_story["j_GAS"][:,1]    = specific_angular_momentum(g.gas, sph_small)
        end
        if GadgetGalaxies.total_mass(g.stars.pos, g.stars.mass) > 0.0
            halo_story["RHMS_STARS"][1] = half_mass_radius(g.stars)
            halo_story["J_STARS"][:,1]  = angular_momentum(g.stars, sph_small)
            halo_story["j_STARS"][:,1]  = specific_angular_momentum(g.stars, sph_small)
            halo_story["BVAL"][1]       = b_value(g.stars, sph_small)
        end
    catch e
        #halo_story["ERRMSG"][1] = string(e)
        #println(e,"\r")
    end

    # Loop over Rest
    cent_id = 0
    for i in 2:length(halo_story["I_SUB"])
        if halo_story["SNAP"][i] < halo_story["SNAP"][i-1]   # This corresponds to an earlier snap and therefore the most massive progenitor
            println("$(halofiles[ii]) ($ii/$(length(halofiles)))   ---   $i / $(length(treefile_df[:, :I_SUB]))\r")
            flush(stdout)
            cent_id                 = halo_story["I_SUB"][i]
            snapshot                = Snapshot(box, halo_story["SNAP"][i])
            snapshot.snapbase
            snapshot.subbase
            g                       = Galaxy(snapshot, halo_story["I_SUB"][i])
            halo_story["RVIR"][i]   = read_galaxy_prop(get_group(g), "RVIR", :physical)
            sph_small               = Sphere(0.1*halo_story["RVIR"][i])
            sph_large               = Sphere(halo_story["RVIR"][i])
            try
                read_halo!(g, units=:physical, props=((:gas, ["POS", "VEL", "MASS"]), (:stars, ["POS", "VEL", "MASS"])), radius_units=:physical, radius=0.1*halo_story["RVIR"][i]) 
                read_halo!(g, units=:physical, props=((:dm, ["POS", "VEL"]),), radius_units=:physical, radius=halo_story["RVIR"][i])    
                if GadgetGalaxies.total_mass(g.dm.pos, g.dm.mass) > 0.0
                    halo_story["RHMS_DM"][i]    = half_mass_radius(g.dm)
                    halo_story["J_DM"][:,i]     = angular_momentum(g.dm, sph_large)
                    halo_story["j_DM"][:,i]     = specific_angular_momentum(g.dm, sph_large)
                end
                if GadgetGalaxies.total_mass(g.gas.pos, g.gas.mass) > 0.0
                    halo_story["RHMS_GAS"][i]   = half_mass_radius(g.gas)
                    halo_story["J_GAS"][:,i]    = angular_momentum(g.gas, sph_small)
                    halo_story["j_GAS"][:,i]    = specific_angular_momentum(g.gas, sph_small)
                end
                if GadgetGalaxies.total_mass(g.stars.pos, g.stars.mass) > 0.0
                    halo_story["RHMS_STARS"][i] = half_mass_radius(g.stars)
                    halo_story["J_STARS"][:,i]  = angular_momentum(g.stars, sph_small)
                    halo_story["j_STARS"][:,i]  = specific_angular_momentum(g.stars, sph_small)
                    halo_story["BVAL"][i]       = b_value(g.stars, sph_small)
                end
            catch e
                #halo_story["ERRMSG"][i] = string(e)
                #println(e,"\r")
            end
        else   # Smaller progenitors
            halo_story["J_orbital"][:,i]    = orbit_J(halo_story["I_SUB"][i], cent_id, "$box/groups_$(@sprintf("%03i", halo_story["SNAP"][i]))/sub_$(@sprintf("%03i", halo_story["SNAP"][i]))")
        end
    end
    save(joinpath(@__DIR__, "halostories/$(chop(halofiles[ii], tail=4)).jld"), 
        "halo_story",   halo_story)
end



### Output

println("\n\n\n---------------------------\n\nNow witness the firepower of this fully armed and operational battle station!")
