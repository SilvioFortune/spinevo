
@doc """
DESCRIPTION:\n
INPUT:\n
OUTPUT:\n
""" ->
function plot_phasespace(; snapNR = "all", ymin=" ", ymax=" ",
    felixID=" ", rootID=" ", subID=" ", pointsize=2, opact=0.9, 
    parttype="STARS", r_by_rvir=0.1, weights=" ", 
    lognorm=lognorm, scale=15, gridres=(25, 15), colorlimits=(0.8,nothing), 
    indir=current_dir_stories, simbox=current_dir_simbox, 
    outdir="./OUT_plot_phasespace", outfile="PhSp", plottype="hexbin", 
    )
    #if colorlimits == ""

    
    if snapNR == "all"
        if typeof(felixID)==Int  # simple case since tree is provided
            mc = load(joinpath(current_dir_stories, find_rootID(felixID)[2]), "merger_collection_$parttype")
            rootID = find_rootID(felixID)[1]
        elseif typeof(rootID)==Int  # simple case since tree is provided
            mc = load(joinpath(current_dir_stories, "halo_$rootID.jld"), "merger_collection_$parttype")
        else
            error("Either felixID or rootID has to be provided.")
        end

        println("$(length(mc["snapNR"])) steps:")
        for i in 1:length(mc["snapNR"])
            print("$(mc["snapNR"][i]) ")
            flush(stdout)
            snapshot    = Snapshot(simbox, mc["snapNR"][i])
            galaxyID_list   = find_merging_progenitors(mc["snapNR"][i]; rootID=rootID, path_to_halostories=indir)
            g           = Galaxy(snapshot, galaxyID_list[1])
            rvir_group  = read_galaxy_prop(get_group(g), "RVIR", :physical)
            if parttype == "STARS"
                read_halo!(g, units=:physical, props=((:stars, ["POS", "VEL", "MASS", "ID"]),), radius_units=:physical, radius=r_by_rvir*rvir_group)
                #velocities  = norm.(eachcol(g.stars.vel))
                radii       = norm.(eachcol(g.stars.pos))
                inside      = findcs(radii, leq=r_by_rvir*rvir_group)
                velocities  = get_radial_velocity.(eachcol(g.stars.pos),eachcol(g.stars.vel))
            elseif parttype == "GAS"
                read_halo!(g, units=:physical, props=((:gas, ["POS", "VEL", "MASS", "ID"]),), radius_units=:physical, radius=r_by_rvir*rvir_group)
                #velocities  = norm.(eachcol(g.gas.vel))
                radii       = norm.(eachcol(g.gas.pos))
                inside      = findcs(radii, leq=r_by_rvir*rvir_group)
                velocities  = get_radial_velocity.(eachcol(g.gas.pos),eachcol(g.gas.vel))
            elseif parttype == "DM"
                read_halo!(g, units=:physical, props=((:dm, ["POS", "VEL", "MASS", "ID"]),), radius_units=:physical, radius=r_by_rvir*rvir_group)
                #velocities  = norm.(eachcol(g.dm.vel))
                radii       = norm.(eachcol(g.dm.pos))
                inside      = findcs(radii, leq=r_by_rvir*rvir_group)
                velocities  = get_radial_velocity.(eachcol(g.dm.pos),eachcol(g.dm.vel))
            else
                error("You have failed me for the last time, parttype = $parttype")
            end

            if weights == "MASS"
                weights = g.stars.mass[inside]
            end
            if opact == "mass"
                opact = g.stars.mass[inside] ./ maximum(g.stars.mass[inside])
            end
            if plottype == "scatter"
                fig, ax = subplots()
                p = ax.scatter(pyplottable(radii[inside]), 
                        pyplottable(velocities[inside]) , 
                        zorder=2   , s=pointsize .* opact    , alpha=opact, 
                        cmap="rainbow_r", c=g.stars.mass[inside]) # brg_r  
                colbar  = fig.colorbar(p, ax=ax)#, ticks=[minimum(as["BVAL_0"    ]), b_ell, b_disk, maximum(as["BVAL_0"    ])])
                #colbar.ax.set_yticklabels(["$(round(minimum(as["BVAL_0"    ]), digits=1))", "ETG", "LTG", "$(round(maximum(as["BVAL_0"    ]), digits=1))"])
                colbar.ax.tick_params(axis="y")
                #ax.legend(loc="lower right", frameon=true, borderpad=1, handlelength=1.8)
                ax.set_title("Halo $(rootID), Snap $(mc["snapNR"][i])")
                ax.set_xlabel("r [kpc]")
                ax.set_ylabel("v [km/s]")
                ax.grid()
                #ax.autoscale()
                fig.set_size_inches(16scale, 9scale)
                fig.tight_layout()
                fig.savefig(joinpath(outdir,outfile*"_halo$(rootID)_snap$(mc["snapNR"][i])"*".png"), bbox_inches="tight", pad_inches=.1)
            elseif plottype == "hexbin"
                plot_hexbins(radii[inside], 
                            velocities[inside], ymin=ymin, ymax=ymax, 
                            weights=weights, scale=scale, colorlimits=colorlimits, lognorm=lognorm, gridres=gridres, 
                            outfile=joinpath(outdir,outfile*"_halo$(rootID)_snap$(mc["snapNR"][i])"*".png"), 
                            xlabel="r [kpc]", 
                            ylabel="v [km/s]", 
                            title="Halo $(rootID), Snap $(mc["snapNR"][i])")#, R = $(round(r_by_rvir*rvir_group, digits=2))")
            else
                error("You have failed me for the last time, plottype = $(plottype)")
            end
        end
    else
        if typeof(felixID)==Int  # simple case since tree is provided
            galaxyID_list   = find_merging_progenitors(snapNR; felixID=felixID, path_to_halostories=indir)
            println(galaxyID_list)
            rootID = find_rootID(felixID)[1]
        elseif typeof(rootID)==Int  # simple case since tree is provided
            galaxyID_list   = find_merging_progenitors(snapNR; rootID=rootID, path_to_halostories=indir)
            println(galaxyID_list)
        elseif typeof(subID)==Int   # here we need to find the right tree first
            galaxyID_list   = find_merging_progenitors(snapNR; subID=subID, path_to_halostories=indir)
            println(galaxyID_list)
        else
            error("Either felixID, subID or rootID has to be provided.")
        end
    
        snapshot    = Snapshot(simbox, snapNR)
        g           = Galaxy(snapshot, galaxyID_list[1])
        rvir_group  = read_galaxy_prop(get_group(g), "RVIR", :physical)
        if parttype == "STARS"
            read_halo!(g, units=:physical, props=((:stars, ["POS", "VEL", "MASS", "ID"]),), radius_units=:physical, radius=r_by_rvir*rvir_group)
            velocities  = get_radial_velocity.(eachcol(g.stars.pos),eachcol(g.stars.vel))
            radii       = norm.(eachcol(g.stars.pos))
            inside      = findcs(radii, leq=r_by_rvir*rvir_group)
        elseif parttype == "GAS"
            read_halo!(g, units=:physical, props=((:gas, ["POS", "VEL", "MASS", "ID"]),), radius_units=:physical, radius=r_by_rvir*rvir_group)
            velocities  = get_radial_velocity.(eachcol(g.gas.pos),eachcol(g.gas.vel))
            radii       = norm.(eachcol(g.gas.pos))
            inside      = findcs(radii, leq=r_by_rvir*rvir_group)
        elseif parttype == "DM"
            read_halo!(g, units=:physical, props=((:dm, ["POS", "VEL", "MASS", "ID"]),), radius_units=:physical, radius=r_by_rvir*rvir_group)
            velocities  = get_radial_velocity.(eachcol(g.dm.pos),eachcol(g.dm.vel))
            radii       = norm.(eachcol(g.dm.pos))
            inside      = findcs(radii, leq=r_by_rvir*rvir_group)
        else
            error("You have failed me for the last time, parttype = $parttype")
        end

        if weights == "MASS"
            weights = g.stars.mass[inside]
        end
        if opact == "mass"
            opact = g.stars.mass[inside] ./ maximum(g.stars.mass[inside])
        end
    
        if plottype == "scatter"
            fig, ax = subplots()
            p = ax.scatter(pyplottable(radii[inside]), 
                    pyplottable(velocities[inside]) , 
                    zorder=2   , s=pointsize .* opact   , alpha=opact, 
                    cmap="rainbow_r", c=g.stars.mass[inside]) # brg_r  
            colbar  = fig.colorbar(p, ax=ax)#, ticks=[minimum(as["BVAL_0"    ]), b_ell, b_disk, maximum(as["BVAL_0"    ])])
            #colbar.ax.set_yticklabels(["$(round(minimum(as["BVAL_0"    ]), digits=1))", "ETG", "LTG", "$(round(maximum(as["BVAL_0"    ]), digits=1))"])
            colbar.ax.tick_params(axis="y")
            #ax.legend(loc="lower right", frameon=true, borderpad=1, handlelength=1.8)
            ax.set_title("Halo $(rootID), Snap $(snapNR)")
            ax.set_xlabel("r [kpc]")
            ax.set_ylabel("v [km/s]")
            ax.grid()
            #ax.autoscale()
            fig.set_size_inches(16scale, 9scale)
            fig.tight_layout()
            fig.savefig(joinpath(outdir,outfile*"_halo$(rootID)_snap$(snapNR)"*".png"), bbox_inches="tight", pad_inches=.1)
        elseif plottype == "hexbin"
            plot_hexbins(radii[inside], 
                        velocities[inside], ymin=ymin, ymax=ymax, 
                        weights=weights, scale=scale, colorlimits=colorlimits, lognorm=lognorm, gridres=gridres, 
                        outfile=joinpath(outdir,outfile*"_halo$(rootID)_snap$(snapNR)"*".png"), 
                        xlabel="r [kpc]", 
                        ylabel="v [km/s]", 
                        title="Halo $(rootID), Snap $(snapNR)")#, R = $(round(r_by_rvir*rvir_group, digits=2))")
        else
            error("You have failed me for the last time, plottype = $(plottype)")
        end
    end
    
    return nothing
end

print("'plot_phasespace'   ")

println()















@doc """
DESCRIPTION:\n
INPUT:\n
OUTPUT:\n
""" ->
function plot_scatter_nxm(xdata, ydata, selections, labels, rows, cols, colorbar; 
    outfile="./OUT_plot_hexbins.png", 
    xmin=" ", xmax=" ", ymin=" ", ymax=" ",
    scale=15, gridres=" ", xscale="linear", yscale="linear", xlabel="", ylabel="", grid=false, fontsize=20, 
    lognorm=true, 
    cmap="rainbow_r", cbarshift=2.0, pointsize=5,
    binned_median = " ", plot_bval=" ",
    calc_norm=true,    # pyplottable
    )

    PyPlot.matplotlib[:rc]("text", usetex=true)
    PyPlot.matplotlib[:rc]("font", family="serif", serif=["Computer Modern"])
    #   Check for dimensions
    if rows*cols != length(keys(selections))
        error("Mismatch of dimensions $n x $m and datasets $(length(keys(selections))).")
    end

    #   Create sanity selection
    all_sanes   = Vector{Int64}(undef, 0)
    sane        = Dict{Int64, Vector{Int64}}()
    for i in 1:rows*cols
        sane[i]     = findall(x-> !isnan.(x), pyplottable(xdata[selections[i]], calc_norm=calc_norm))[findall(x-> !isnan.(x), pyplottable(xdata[selections[i]], calc_norm=calc_norm)) .∈ Ref(Set(findall(x-> !isnan.(x), pyplottable(ydata[selections[i]], calc_norm=calc_norm))))]
        all_sanes   = vcat( all_sanes, selections[i][sane[i]] )
    end


    if xmin == " "
        xmin = minimum(xdata[all_sanes])
    end
    if xmax == " "
        xmax = maximum(xdata[all_sanes])
    end
    if ymin == " "
        ymin = minimum(ydata[all_sanes])
    end
    if ymax == " "
        ymax = maximum(ydata[all_sanes])
    end
    if gridres == " "
        gridres = (50,25)
    elseif typeof(gridres) == Int
        gridres = ( gridres, Int(cld(gridres,2)) )
    end

    ### Figure
    #rc          = pyimport("matplotlib.rc")
    pltcolors   = pyimport("matplotlib.colors")
    fig, axs    = subplots(rows, cols; sharex=true, sharey=true)

    colr    = pyplottable(colorbar)#./as09["Δlookbacktime"] ./-1000)
    maxcolr = maximum(colr[findall(x->!isnan.(x), colr)])
    slct    = findall(x->!isnan.(x), colr)
    # now the real thing
    scatters  = Dict{Int64, Any}()
    n   = 0
    m   = 1
    for i in 1:rows*cols
        n += 1
        if n == rows+1
            n   = 1
            m  += 1
        end
        xd = xdata[selections[i]][sane[i]]
        yd = ydata[selections[i]][sane[i]]
        scatters[i]   = axs[n, m].scatter( xd, yd, 
                                        s= (pointsize/maxcolr) .* colr[selections[i]][sane[i]]    , alpha=0.6,
                                        cmap="rainbow_r", c=colr[selections[i]][sane[i]])
        if plot_bval == "y"
            axs[n, m].hlines(b_disk, xmin=xmin, xmax=xmax, colors="blue")
            axs[n, m].hlines(b_ell, xmin=xmin, xmax=xmax, colors="red")
        elseif plot_bval == "x"
            axs[n, m].vlines(b_disk, ymin=ymin, ymax=ymax, colors="blue")
            axs[n, m].vlines(b_ell, ymin=ymin, ymax=ymax, colors="red")
        end
        axs[n, m].legend(loc="upper left", frameon=true, borderpad=1, handlelength=1.8, fontsize=floor(0.5*fontsize))

        if binned_median == "y"
            if yscale == "log"
                medians, edges, bin_number     = stats.binned_statistic(log10.(yd), xd, statistic="median", bins=gridres[2])
                bin_width   = edges[2] - edges[1]
                centers     = edges[2:end] .- (bin_width/2)
                pctl_68, edges_up, bin_number_up     = stats.binned_statistic(log10.(yd), xd, statistic=pctl68upper, bins=gridres[2])
                pctl_32, edges_lo, bin_number_lo     = stats.binned_statistic(log10.(yd), xd, statistic=pctl68lower, bins=gridres[2])
                axs[n, m].fill_betweenx(log10.(centers), pctl_32, pctl_68, color="blue", alpha=0.2)#, transform=ax.get_xaxis_transform())
                axs[n, m].plot(medians, log10.(centers), "b-", title="Median")
            elseif yscale == "linear" # tested
                medians, edges, bin_number     = stats.binned_statistic(yd, xd, statistic="median", bins=gridres[2])
                bin_width   = edges[2] - edges[1]
                centers     = edges[2:end] .- (bin_width/2)
                pctl_68, edges_up, bin_number_up     = stats.binned_statistic(yd, xd, statistic=pctl68upper, bins=gridres[2])
                pctl_32, edges_lo, bin_number_lo     = stats.binned_statistic(yd, xd, statistic=pctl68lower, bins=gridres[2])
                axs[n, m].fill_betweenx(centers, pctl_32, pctl_68, color="blue", alpha=0.2)#, transform=ax.get_xaxis_transform())
                axs[n, m].plot(medians, centers, "b-", label="Median")
            else
                error("Unknown yscale: $yscale")
            end
        elseif binned_median == "x"
            if xscale == "log"
                medians, edges, bin_number     = stats.binned_statistic(log10.(xd), yd, statistic="median", bins=gridres[1])
                bin_width   = edges[2] - edges[1]
                centers     = edges[2:end] .- (bin_width/2)
                pctl_68, edges_up, bin_number_up     = stats.binned_statistic(log10.(xd), yd, statistic=pctl68upper, bins=gridres[1])
                pctl_32, edges_lo, bin_number_lo     = stats.binned_statistic(log10.(xd), yd, statistic=pctl68lower, bins=gridres[1])
                axs[n, m].fill_between(log10.(centers), pctl_32, pctl_68, color="blue", alpha=0.2)
                axs[n, m].plot(log10.(centers), medians, "b-", title="Median", linewidth=5)
            elseif xscale == "linear" #tested
                medians, edges, bin_number     = stats.binned_statistic(xd, yd, statistic="median", bins=gridres[1])
                bin_width   = edges[2] - edges[1]
                centers     = edges[2:end] .- (bin_width/2)
                pctl_68, edges_up, bin_number_up     = stats.binned_statistic(xd, yd, statistic=pctl68upper, bins=gridres[1])
                pctl_32, edges_lo, bin_number_lo     = stats.binned_statistic(xd, yd, statistic=pctl68lower, bins=gridres[1])
                axs[n, m].fill_between(centers, pctl_32, pctl_68, color="blue", alpha=0.2)
                axs[n, m].plot(centers, medians, "b-", label="Median", linewidth=5)
            else
                error("Unknown xscale: $xscale")
            end
        end
    end
    

    #ax.set_yticks([1e-2,1e-1,1e0,1e1,1e2,1e3,1e4,1e5,1e6])
    #rc("axes", ticksize="40")
    fig.supxlabel(xlabel, fontsize="$(fontsize)")
    fig.supylabel(ylabel, fontsize="$(fontsize)")
    #axs[3,1].set_xscale(xscale)
    #axs[3,2].set_xscale(xscale)
    #axs[3,1].set_xlim([xmin, xmax])
    #axs[3,2].set_xlim([xmin, xmax])
    #axs[1,1].set_yscale(yscale)
    #axs[1,1].set_ylabel(ylabel, fontsize="$(2*scale)")
    #axs[1,1].set_ylim([ymin, ymax])

    #ticker   = pyimport("matplotlib.ticker")
    #nbins = length(axs[3,1].get_xticklabels()) # added 
    #axs[1,1].tick_params(axis="y", labelsize="$(scale)")
    #axs[2,1].tick_params(axis="y", labelsize="$(scale)")
    ##axs[2,1].yaxis.set_major_locator(ticker.MaxNLocator(prune="upper")) # added 
    #axs[3,1].tick_params(axis="both", labelsize="$(scale)")
    #axs[3,1].yaxis.set_major_locator(ticker.MaxNLocator(prune="lower")) # added 
    #axs[3,1].xaxis.set_major_locator(ticker.MaxNLocator(prune="upper")) # added 
    for i in 1:cols
        axs[rows,i].tick_params(axis="x", labelsize="$(scale)")
    end
    for i in 1:rows
        axs[i,1].tick_params(axis="y", labelsize="$(scale)")
    end
    if grid
        axs.grid()
    end
    colbar  = fig.colorbar(scatters[rows*cols], ax=axs[:,cols], anchor=(cbarshift, 0.5))
    colbar.ax.tick_params(axis="y", labelsize="$(1*scale)")
    
    #fig.tight_layout()
    fig.subplots_adjust(wspace=0, hspace=0)
    fig.set_size_inches(scale, 0.7scale)
    fig.savefig(outfile, bbox_inches="tight", pad_inches=.1)
    close("all")

    PyPlot.matplotlib[:rc]("text", usetex=false)
    
    return scatters
end

print("'plot_scatter_nxm'   ")