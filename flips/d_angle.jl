### Packages

using Printf
using LinearAlgebra
using JLD
using Statistics
using PyCall
using PyPlot
using LaTeXStrings
using Missings
using GadgetIO
using GadgetUnits
using GadgetGalaxies
using Unitful
using UnitfulAstro
using Cosmology


### Setup

flipangle       = 30
th_MSOLAR       = 1e10
dir_halostories = "/home/moon/sfortune/spinevo/halostories_update_stars"
box             = "/HydroSims/Magneticum/Box4/uhr_test/"

storyfilelist   = readdir(dir_halostories)

### Spin change histograms

angle_data = Dict(
    "CENTRALS"      => Array{Int64}(undef, 2, 0), 
    "ANGLES"        => Array{Float64}(undef, 0), 
    "M_felix"       => Array{Float64}(undef, 0), 
    "δM_felix"      => Array{Float64}(undef, 0), 
    "M_fromJ"       => Array{Float64}(undef, 0), 
    "M2_mm"          => Array{Float64}(undef, 0), 
    "J_orbit_mm"    => Array{Float64}(undef, 3, 0), 
    "δM_fromJ"      => Array{Float64}(undef, 0),
    "REDSHIFT"      => Array{Float64}(undef, 0),
    "LOOKBACKTIME"  => Array{Float64}(undef, 0),
    "BVAL"          => Array{Float64}(undef, 0),
    "STORY_FLIPS"   => Array{Int64}(undef, 0), 
    "STORY_IFILE"   => Array{Int64}(undef, 0), 
    "STORY_BVAL"    => Array{Float64}(undef, 0),
    "STORY_BVAL_START"    => Array{Float64}(undef, 0),
    "STORY_M_fromJ" => Array{Float64}(undef, 0),
    "I_FILE"        => Array{Int64}(undef, 0), 
    "I_SUB"         => Array{Int64}(undef, 0) )


limit_filelist  = 10#length(storyfilelist)
for i in 1:limit_filelist
    halo_story  = load(joinpath(dir_halostories, storyfilelist[i]), "halo_story")
    head        = read_header("$box/groups_$(@sprintf("%03i", halo_story["SNAP"][1]))/sub_$(@sprintf("%03i", halo_story["SNAP"][1]))")
    snapshot    = Snapshot(box, halo_story["SNAP"][1])
    snapshot.snapbase
    snapshot.subbase
    g           = Galaxy(snapshot, halo_story["I_SUB"][1])
    println("$i")
    flush(stdout)
    #print("$(halo_story["I_SUB"][1]) ? $(get_first_subhalo(get_group(g)).isub)   ")
    if get_first_subhalo(get_group(g)).isub == halo_story["I_SUB"][1]
        angle_data["CENTRALS"] = hcat(angle_data["CENTRALS"], [halo_story["I_SUB"][1], halo_story["ID"]])
        if convert_units_physical_mass(halo_story["M_STARS"][1], head) > th_MSOLAR
            #print("$(convert_units_physical_mass(halo_story["M_STARS"][1], head))   yyy   ")
            spinTOT = Array{Float64}(undef, 3, 0)
            M_comp  = Array{Float64}(undef, 3, 0)
            n_flips = 0
            for ii in length(halo_story["SNAP"]):-1:1
                #if !ismissing(halo_story["J_DM"][1,ii])
                if !ismissing(halo_story["J_STARS"][1,ii])   &&   norm(halo_story["J_STARS"][:,ii]) > 0. && halo_story["SNAP"][ii]
                    #J_tot   = halo_story["J_DM"][:,ii]
                    #M_comp  = [norm(halo_story["J_DM"][:,ii]) / norm(halo_story["j_DM"][:,ii]), 0., 0.]
                    #if !ismissing(halo_story["J_GAS"][1,ii])
                        #J_tot         .+= halo_story["J_GAS"][:,ii]
                        #M_comp[2,end]   = norm(halo_story["J_GAS"][:,ii]) / norm(halo_story["j_GAS"][:,ii])
                    #end
                    #if !ismissing(halo_story["J_STARS"][1,ii])
                    J_tot           = halo_story["J_STARS"][:,ii]
                    spinTOT         = hcat(spinTOT, J_tot ./ norm(J_tot))
                    #@show spinTOT
                        #M_comp[3,end]   = norm(halo_story["J_STARS"][:,ii]) / norm(halo_story["j_STARS"][:,ii])
                    M_comp          = hcat(M_comp, [0., 0., norm(J_tot) / norm(halo_story["j_STARS"][:,ii])])
                    #else
                        #J_tot           = [0.,0.,0.]
                        #M_comp  = hcat(M_comp, [0., 0., 0.])
                    #end
                    head                = read_header("$box/groups_$(@sprintf("%03i", halo_story["SNAP"][ii]))/sub_$(@sprintf("%03i", halo_story["SNAP"][ii]))")
                    #println("$(length(M_comp[3,:])) | $(M_comp[3,end]) | $(convert_units_physical_mass(halo_story["M_STARS"][ii], head)) | ")
                    if length(M_comp[3,:]) > 1   &&   M_comp[3,end] > th_MSOLAR   &&   convert_units_physical_mass(halo_story["M_STARS"][ii], head) > th_MSOLAR
                        #println("$i/$limit_filelist   $ii/$(length(halo_story["SNAP"]))   --- Pass!")
                        angle_data["LOOKBACKTIME"]  = vcat(angle_data["LOOKBACKTIME"], ustrip(lookback_time(cosmology(h=head.h0, OmegaM=head.omega_0), head.z)))
                        angle_data["REDSHIFT"]      = vcat(angle_data["REDSHIFT"], halo_story["REDSHIFT"][ii])
                        angle_data["I_SUB"]         = vcat(angle_data["I_SUB"], halo_story["I_SUB"][ii])
                        angle_data["ANGLES"]        = vcat(angle_data["ANGLES"], acosd((transpose(spinTOT[:,end])*spinTOT[:,end-1])))
                        if angle_data["ANGLES"][end] > flipangle
                            n_flips              += 1
                        end
                        angle_data["M_fromJ"]       = vcat(angle_data["M_fromJ"], sum(M_comp[3,end]))
                        angle_data["δM_fromJ"]      = vcat(angle_data["δM_fromJ"], sum(M_comp[3,end]) - sum(M_comp[3,end-1]))
                        angle_data["M_felix"]       = vcat(angle_data["M_felix"], convert_units_physical_mass(halo_story["M_STARS"][ii], head))
                        angle_data["δM_felix"]      = vcat(angle_data["δM_felix"], convert_units_physical_mass(halo_story["M_STARS"][ii], head)-convert_units_physical_mass(halo_story["M_STARS"][ii+1], head))
                        angle_data["I_FILE"]        = vcat(angle_data["I_FILE"], halo_story["ID"])
                        angle_data["BVAL"]          = vcat(angle_data["BVAL"], halo_story["BVAL"][ii])
                    end
                end
            end
            angle_data["STORY_FLIPS"]   = vcat(angle_data["STORY_FLIPS"],   n_flips)
            angle_data["STORY_IFILE"]   = vcat(angle_data["STORY_IFILE"],   halo_story["ID"])
            angle_data["STORY_BVAL"]    = vcat(angle_data["STORY_BVAL"],    angle_data["BVAL"][end])
            angle_data["STORY_BVAL_START"]    = vcat(angle_data["STORY_BVAL_START"],    angle_data["BVAL"][1])
            angle_data["STORY_M_fromJ"] = vcat(angle_data["STORY_M_fromJ"], angle_data["M_fromJ"][end])
        end
    #else
        #print("$(get_first_subhalo(get_group(g)).isub)!=$(halo_story["I_SUB"][1])  $(convert_units_physical_mass(halo_story["M_STARS"][1], head))<$(th_MSOLAR)   ")
    end
end

println("\n\nNow witness the firepower of this fully armed and operational battle station!\n$(size(angle_data["ANGLES"]))")
println("$(size(angle_data["STORY_FLIPS"]))   ---   $(size(angle_data["STORY_IFILE"]))   ---   $(size(angle_data["STORY_BVAL"]))   ---   $(size(angle_data["STORY_M_fromJ"]))")

save(joinpath(@__DIR__, "angle_$(flipangle)_data_CENTRALS_MSTARSgt$(th_MSOLAR).jld"), 
    "angle_data",   angle_data)