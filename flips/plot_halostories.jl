
flipangle       = 180
th_MSOLAR       = 1e10
storyfilelist   = readdir(joinpath(@__DIR__, "halostories/"))
box             = "/HydroSims/Magneticum/Box4/uhr_test/"


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

# Functions

# Tranform Cartesian vector into [r,θ,ϕ] Spherical vector
function cartesian_to_spherical(x)
    s = zeros(3)
    s[1]    = norm(x) # Radius
    s[2]    = acosd(x[3] / s[1])# * 180 / π # θ[°]
    s[3]    = atand(x[2] / x[1])# * 180 / π # ϕ[°]
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


angle_data  = load(joinpath(@__DIR__, "angle_data_CENTRALS_MSTARgt10.jld"), "angle_data")

limit_filelist  = 10#length(angle_data["FLIPS_IFILE"])
println(limit_filelist)
for ii in 1:limit_filelist
    try
        haloID      = angle_data["FLIPS_IFILE"][ii]
        halo_story  = load(joinpath(@__DIR__, "halostories/halo_$haloID.jld"), "halo_story")
        #head        = read_header("$box/groups_$(@sprintf("%03i", halo_story["SNAP"][1]))/sub_$(@sprintf("%03i", halo_story["SNAP"][1]))")
        #snapshot    = Snapshot(box, halo_story["SNAP"][1])
        #snapshot.snapbase
        #snapshot.subbase
        #g           = Galaxy(snapshot, halo_story["I_SUB"][1])
        println("$ii $haloID")
        flush(stdout)
        #print("$(halo_story["I_SUB"][1]) ? $(get_first_subhalo(get_group(g)).isub)   ")
        snap                = Array{Int64}(undef, 0)
        z                   = Array{Float64}(undef, 0)
        spinTOT             = Array{Float64}(undef, 3, 0)
        j_dm                = Array{Float64}(undef, 0)
        j_gas               = Array{Float64}(undef, 0)
        j_stars             = Array{Float64}(undef, 0)
        M_tot               = Array{Float64}(undef, 0)
        spherical_spinTOT   = Array{Float64}(undef, 3, 0)
        δ_spinTOT           = zeros(1)
        flipsnaps           = Array{Int64}(undef, 0)
        refvec              = zeros(3)
        R                   = zeros(3,3)
        M_comp              = Array{Float64}(undef, 3, 0)
        progenitors         = Array{Int64}(undef, 0)
        lookbacktime        = Array{Float64}(undef, 0)
        for i in length(halo_story["SNAP"]):-1:1
            #if !ismissing(halo_story["J_DM"][1,i]) && !ismissing(halo_story["J_DM"][2,i]) && !ismissing(halo_story["J_DM"][3,i])
            if !ismissing(halo_story["J_STARS"][1,i]) && !ismissing(halo_story["J_STARS"][2,i]) && !ismissing(halo_story["J_STARS"][3,i])
                head                = read_header("$box/groups_$(@sprintf("%03i", halo_story["SNAP"][i]))/sub_$(@sprintf("%03i", halo_story["SNAP"][i]))")
                #lookbacktime        = vcat(lookbacktime, convert_units_age(1/(1+halo_story["REDSHIFT"][i]), read_header("$box/groups_$(@sprintf("%03i", halo_story["SNAP"][i]))/sub_$(@sprintf("%03i", halo_story["SNAP"][i]))"), :physical))
                lookbacktime        = vcat(lookbacktime, ustrip(lookback_time(cosmology(h=head.h0, OmegaM=head.omega_0), head.z)))
                progenitors         = vcat(progenitors, count(p->(p==halo_story["SNAP"][i]), halo_story["SNAP"]))
                #println("$(progenitors[end]) --- $(halo_story["SNAP"][i])")
                snap                = vcat(snap, halo_story["SNAP"][i])
                z                   = vcat(z, halo_story["REDSHIFT"][i])
                M_tot               = vcat(M_tot, halo_story["M_DM"][i]+halo_story["M_GAS"][i]+halo_story["M_STARS"][i])
                #J_tot               = halo_story["J_DM"][:,i]
                J_tot               = halo_story["J_STARS"][:,i]
                j_stars                = vcat(j_stars, norm(halo_story["j_STARS"][:,i]))
                M_comp              = hcat(M_comp, [0., 0., norm(halo_story["J_STARS"][:,i]) / norm(halo_story["j_STARS"][:,i])])
                if !ismissing(halo_story["J_GAS"][1,i]) && !ismissing(halo_story["J_GAS"][2,i]) && !ismissing(halo_story["J_GAS"][3,i])
                    #J_tot         .+= halo_story["J_GAS"][:,i]
                    j_gas           = vcat(j_gas, norm(halo_story["j_GAS"][:,i]))
                    M_comp[2,end]   = norm(halo_story["J_GAS"][:,i]) / norm(halo_story["j_GAS"][:,i])
                else
                    j_gas           = vcat(j_gas, 0)
                end
                if !ismissing(halo_story["J_DM"][1,i]) && !ismissing(halo_story["J_DM"][2,i]) && !ismissing(halo_story["J_DM"][3,i])
                    #J_tot         .+= halo_story["J_STARS"][:,i]
                    j_dm         = vcat(j_dm, norm(halo_story["j_DM"][:,i]))
                    M_comp[1,end]   = norm(halo_story["J_DM"][:,i]) / norm(halo_story["j_DM"][:,i])
                else
                    j_stars         = vcat(j_stars, 0)
                end
                if length(z) == 1
                    R               = aligner(J_tot)
                    spinTOT             = hcat(spinTOT, 
                        R*(J_tot ./ norm(J_tot)))
                    spherical_spinTOT   = hcat(spherical_spinTOT, 
                        cartesian_to_spherical(R*J_tot))
                    refvec          = spinTOT[:,end]
                else
                    spinTOT             = hcat(spinTOT, 
                        R*(J_tot ./ norm(J_tot)))
                    spherical_spinTOT   = hcat(spherical_spinTOT, 
                        cartesian_to_spherical(R*J_tot))
                    #δ_spinTOT       = vcat(δ_spinTOT, acosd((transpose(spinTOT[:,end])*refvec)))
                    δ_spinTOT       = vcat(δ_spinTOT, acosd((transpose(spinTOT[:,end])*spinTOT[:,end-1])))
                    if δ_spinTOT[end] > flipangle
                        snap        = vcat(snap, halo_story["SNAP"][i])
                        z           = vcat(z, halo_story["REDSHIFT"][i])
                        M_tot       = vcat(M_tot, halo_story["M_DM"][i]+halo_story["M_GAS"][i]+halo_story["M_STARS"][i])
                        spinTOT     = hcat(spinTOT,
                            R*(J_tot ./ norm(J_tot)))
                        spherical_spinTOT   = hcat(spherical_spinTOT,
                            cartesian_to_spherical(R*J_tot))
                        δ_spinTOT   = vcat(δ_spinTOT, 0)
                        flipsnaps   = vcat(flipsnaps, halo_story["SNAP"][i])
                        refvec      = spinTOT[:,end]
                        j_dm        = vcat(j_dm, norm(halo_story["j_DM"][:,i]))
                        if !ismissing(halo_story["j_GAS"][1,i])
                            j_gas           = vcat(j_gas, norm(halo_story["j_GAS"][:,i]))
                        end
                        if !ismissing(halo_story["j_STARS"][1,i])
                            j_stars         = vcat(j_stars, norm(halo_story["j_STARS"][:,i]))
                        end
                    end
                end
            end
        end
    
        ### Figure arrows
        scale           = 1
        proj_vec        = spinTOT[1,:]
        x_vec           = spinTOT[2,:]
        plotheight      = 3
        plotwidth       = 16
        z_2dec  = round.(z,digits=2) 
        
        #pltcolors = pyimport("matplotlib.colors")
        pltcm           = pyimport("matplotlib.cm")
        pltcolors       = pyimport("matplotlib.colors")
        colormap        = pltcm.get_cmap(name="coolwarm")
        
        
        
        
        #fig = figure(figsize=(16,9))
        #ax1 = fig.add_subplot()
        
        y = zeros(length(lookbacktime))
        x = zeros(length(lookbacktime))
        
        fig, ax = subplots()
        
        ax.set_title("Evolution of J_STARS Orientation of Halo $(haloID)")
        
        #ax.set_ylim(-(maximum(lookbacktime)-minimum(lookbacktime))/2*(plotheight/plotwidth), (maximum(lookbacktime)-minimum(lookbacktime))/2*(plotheight/plotwidth))
        ax.set_ylim(-1, 1)
        
        ax.quiver(lookbacktime, y, 
                x_vec, spinTOT[3,:],
                color=colormap(proj_vec.*0.5.+0.5), scale_units="y", scale=scale, alpha=0.7, width=0.005, headwidth=2, headlength=3, edgecolor="black", lw=0.5)
        
        ax.set_xlabel("Lookback Time [Gyr]")
        ax.invert_xaxis()
        ax.set_ylabel("Y-Axis Component")
        ax.grid()
        ax2 = ax.twiny()
        ax2.set_xlim(ax.get_xlim())
        ax2.set_xticks(lookbacktime[1:3:end])
        ax2.set_xticklabels(z_2dec[1:3:end])
        ax2.set_xlabel("z")
        clb = fig.colorbar(pltcm.ScalarMappable(norm=pltcolors.Normalize(vmin=-1, vmax=1), cmap="coolwarm"))
        clb.ax.set_ylabel("Z-Axis Component")
        clb.ax.yaxis.set_label_position("left")
        fig.gca().set_aspect("equal", adjustable="datalim")
        fig.set_size_inches(plotwidth, plotheight)
        
        fig.savefig(joinpath(@__DIR__, "plots/halo_$(haloID)_vectors.png"), bbox_inches="tight", pad_inches=.1)
        
        
        ### Figure plots
        N_plots = 4
        z_2dec  = round.(z,digits=2) 
        fig, ax = subplots(N_plots)
        
        
        ax[1].plot(lookbacktime, δ_spinTOT, 
                "b-", lw=5, label="J_STARS --- δ", alpha=1, zorder=3)
        ax[1].plot(lookbacktime, spherical_spinTOT[2,:], 
                "g-", lw=2, label="J_STARS --- θ", alpha=1, zorder=1)
        ax[1].plot(lookbacktime[2:end], spherical_spinTOT[3,2:end], 
                "c-", lw=2, label="J_STARS --- ϕ", alpha=1, zorder=2)
        ax[1].bar(flipsnaps, 150, color="black")
        ax[1].hlines(flipangle,minimum(lookbacktime), maximum(lookbacktime), label="$flipangle °", color="grey")
        #ax[1].set_title("J_STARS Orientation")
        ax[1].set_xlabel("Lookback Time [Gyr]")
        ax[1].invert_xaxis()
        ax[1].set_ylabel("Angle [°]")
        ax[1].grid()
        ax[1].axhline(color="black")
        ax[1].legend(loc="upper right", frameon=true, borderpad=1, handlelength=1.8)
        ax2 = ax[1].twiny()
        ax2.set_xlim(ax[1].get_xlim())
        ax2.set_xticks(lookbacktime[1:3:end])
        ax2.set_xticklabels(z_2dec[1:3:end])
        ax2.set_xlabel("z")
        #ax[1].invert_xaxis()
        
        
        #ax[2].plot(lookbacktime, progenitors, 
        #        "-", lw=2, label="Halo $haloID", color="navy", alpha=1, zorder=2)
        #ax[2].set_xlabel("Lookback Time [Gyr]")
        #ax[2].set_title("Number of Future Progenitors")
        #ax[2].invert_xaxis()
        #ax[2].set_ylabel("N_progenitors")
        ##ax[2].set_yscale("log")
        #ax[2].grid()
        #ax[2].legend(loc="upper right", frameon=true, borderpad=1, handlelength=1.8)
        #ax2 = ax[2].twiny()
        #ax2.set_xlim(ax[2].get_xlim())
        #ax2.set_xticks(lookbacktime[1:3:end])
        #ax2.set_xticklabels(z_2dec[1:3:end])
        #ax2.set_xlabel("z")
        
        #ax[2].plot(lookbacktime, abs_diff_spinTOT, 
                #"c-", lw=2, label="Halo_$haloID diffvec", alpha=1, zorder=1)
        ax[2].plot(lookbacktime, spherical_spinTOT[1,:], 
                "b-", lw=2, label="J_STARS Magnitude", alpha=1, zorder=2)
        #ax[2].set_title("J_STARS Magnitude")
        ax[2].set_xlabel("Lookback Time [Gyr]")
        ax[2].invert_xaxis()
        ax[2].set_yscale("log")
        ax[2].set_ylabel("J [M⊙ kpc km/s]")
        ax[2].grid()
        ax[2].legend(loc="upper right", frameon=true, borderpad=1, handlelength=1.8)
        ax2 = ax[2].twiny()
        ax2.set_xlim(ax[2].get_xlim())
        ax2.set_xticks(lookbacktime[1:3:end])
        ax2.set_xticklabels(z_2dec[1:3:end])
        ax2.set_xlabel("z")
        
        ax[3].plot(lookbacktime, j_dm, 
                "-", lw=2, label="j_DM Magnitude", color="darkred", alpha=1, zorder=4)
        ax[3].plot(lookbacktime, j_gas, 
                "-", lw=2, label="j_GAS Magnitude", color="orangered", alpha=1, zorder=3)
        ax[3].plot(lookbacktime, j_stars, 
                "-", lw=2, label="j_STARS Magnitude", color="gold", alpha=1, zorder=2)
        #ax[3].set_title("Specific Angular Momentum Magnitudes")
        ax[3].set_xlabel("Lookback Time [Gyr]")
        ax[3].invert_xaxis()
        ax[3].set_yscale("log")
        ax[3].set_ylabel("j [kpc km/s]")
        ax[3].grid()
        ax[3].legend(loc="upper right", frameon=true, borderpad=1, handlelength=1.8)
        ax2 = ax[3].twiny()
        ax2.set_xlim(ax[4].get_xlim())
        ax2.set_xticks(lookbacktime[1:3:end])
        ax2.set_xticklabels(z_2dec[1:3:end])
        ax2.set_xlabel("z")
        
        #ax[4].plot(lookbacktime, M_tot, 
                #"b-", lw=2, label="M_tot", alpha=1, zorder=2)
        ax[4].plot(lookbacktime, M_comp[1,:], 
                "-", lw=2, label="M_DM", color="darkred", alpha=1, zorder=2)
        ax[4].plot(lookbacktime, M_comp[2,:], 
                "-", lw=2, label="M_GAS", color="orangered", alpha=1, zorder=2)
        ax[4].plot(lookbacktime, M_comp[3,:], 
                "-", lw=2, label="M_STARS", color="gold", alpha=1, zorder=2)
        ax[4].set_xlabel("Lookback Time [Gyr]")
        #ax[4].set_title("Halo Mass")
        ax[4].invert_xaxis()
        ax[4].set_ylabel("Mass [M⊙]")
        ax[4].set_yscale("log")
        ax[4].grid()
        ax[4].legend(loc="upper right", frameon=true, borderpad=1, handlelength=1.8)
        ax2 = ax[4].twiny()
        ax2.set_xlim(ax[4].get_xlim())
        ax2.set_xticks(lookbacktime[1:3:end])
        ax2.set_xticklabels(z_2dec[1:3:end])
        ax2.set_xlabel("z")
        
        
        scale=0.7
        fig.set_size_inches(16scale, 5*N_plots*scale)
        fig.tight_layout()
        
        
        fig.savefig(joinpath(@__DIR__, "plots/halo_$(haloID)_graphs.png"), bbox_inches="tight", pad_inches=.1)
        
    catch e
    end
end


println("\n\nNow witness the firepower of this fully armed and operational battle station!\n")