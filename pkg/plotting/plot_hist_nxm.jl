
@doc """
DESCRIPTION:\n
INPUT:\n
OUTPUT:\n
""" ->
function plot_hist_nxm(xdata, selections, labels, rows, cols, ; 
    outfile="./OUT_plot_histograms.png", 
    xmin=" ", xmax=" ", ymin=" ", ymax=" ", density=true,
    scale=15, n_bins=" ", xscale="linear", yscale="linear", xlabel="", ylabel="", grid=false, fontsize=20, 
    lognorm=true, 
    cmap="mediumblue", cbarshift=2.0,
    binned_median = " ", plot_bval=" ", plot_mergers=" ",
    calc_norm=true,    # pyplottable()
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
        ymin = nothing
    end
    if ymax == " "
        ymax = nothing
    end
    if n_bins == " "
        n_bins = 50
    end

    ### Figure
    #rc          = pyimport("matplotlib.rc")
    pltcolors   = pyimport("matplotlib.colors")
    fig, axs    = subplots(rows, cols; sharex=true, sharey=true)

    # now the real thing
    histograms  = Dict{Int64, Any}()
    n   = 0
    m   = 1
    for i in 1:rows*cols
        n += 1
        if n == rows+1
            n   = 1
            m  += 1
        end
        xd = xdata[selections[i]][sane[i]]
        histograms[i]   = PyPlot.hist( 
                            xd, bins=n_bins    , rwidth=0.9    , alpha=0.7, color="darkred", label="1 Gyr"
                            )
        bcenters    = 0.5 .* (histograms[i][2][2:end]  .+ histograms[i][2][1:end-1])
        axs[n,m].plot(bcenters , histograms[i][1] , ".-", ms=10, label=labels[i] , color=cmap   )
        axs[n,m].set_xlim([xmin,xmax])
        axs[n,m].set_ylim([ymin,ymax])
    end
    

    fig.supxlabel(xlabel, fontsize="$(fontsize)")
    fig.supylabel(ylabel, fontsize="$(fontsize)")
    for i in 1:cols
        axs[rows,i].tick_params(axis="x", labelsize="$(scale)")
    end
    for i in 1:rows
        axs[i,1].tick_params(axis="y", labelsize="$(scale)")
    end
    if grid
        axs.grid()
    end
    colbar  = fig.colorbar(histograms[rows*cols], ax=axs[:,cols], anchor=(cbarshift, 0.5))
    colbar.ax.tick_params(axis="y", labelsize="$(1*scale)")
    
    #fig.tight_layout()
    fig.subplots_adjust(wspace=0, hspace=0)
    fig.set_size_inches(scale, 0.7scale)
    fig.savefig(outfile, bbox_inches="tight", pad_inches=.1)
    close("all")

    PyPlot.matplotlib[:rc]("text", usetex=false)
    
    return histograms
end

print("'plot_hist_nxm'   ")