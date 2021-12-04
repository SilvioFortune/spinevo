println("   Plotting...")
using PyPlot
using ColorSchemes
using JSServe
using WGLMakie
WGLMakie.activate!()

include("/home/moon/sfortune/spinevo/pkg/plot_group_collection.jl")


@doc """
DESCRIPTION:\n
INPUT:\n
OUTPUT:\n
""" ->
function plot_halostory_arrows(; outfile="./OUT_plot_halostory_arrows.png", felixID=" ", lastID=" "
    , parttype="STARS", snap=false
    , scale=1, plotwidth=16, plotheight=3, gridspace=0.6 #Gyr
    , indir   = "/home/moon/sfortune/spinevo/halostories_v20211204_min0.0Gyr"
    )

    halo_filestring = find_halo_file(indir=indir,felixID=felixID,lastID=lastID)

    merger_collection = load(joinpath(indir, halo_filestring), string("merger_collection_", parttype))
    
    lbt     = Array{Float64}(undef, 1)
    lbt[1]  = merger_collection["LOOKBACKTIME"][1]
    rs      = Array{Float64}(undef, 1)
    rs[1]   = merger_collection["REDSHIFT"][1]
    snaps    = Array{Int64}(undef, 1)
    snaps[1] = merger_collection["SNAP"][1]
    for i in 2:length(merger_collection["SNAP"])-1
        if lbt[end] - merger_collection["LOOKBACKTIME"][i] > gridspace
            snaps    = vcat( snaps, merger_collection["SNAP"        ][i] )
            lbt     = vcat( lbt , merger_collection["LOOKBACKTIME"][i] )
            rs      = vcat( rs  , merger_collection["REDSHIFT"    ][i] )
        end
    end
    snaps    = vcat( snaps, merger_collection["SNAP"        ][end] )
    lbt     = vcat( lbt , merger_collection["LOOKBACKTIME"][end] )
    rs      = vcat( rs  , merger_collection["REDSHIFT"    ][end] )
    
    
    J_norm      = vcat( transpose(norm.(eachcol(merger_collection["J_main"]))), transpose(norm.(eachcol(merger_collection["J_main"]))), transpose(norm.(eachcol(merger_collection["J_main"]))) )
    
    proj_vec    = (merger_collection["J_main"]./J_norm)[1,:]
    x_vec       = (merger_collection["J_main"]./J_norm)[2,:]
    y_vec       = (merger_collection["J_main"]./J_norm)[3,:]
    
    pltcm       = pyimport("matplotlib.cm")
    pltcolors   = pyimport("matplotlib.colors")
    colormap    = pltcm.get_cmap(name="coolwarm")
    
    fig, ax = subplots()
    
    ax.set_ylim(-1, 1)
    ax.set_yticks([-1,0,1])
    ax.set_xlim( minimum(merger_collection["LOOKBACKTIME"])-0.2-2*maximum([0, x_vec[end]]), maximum(merger_collection["LOOKBACKTIME"])+0.2-2*minimum([0, x_vec[1]]) )
    ax.set_xticks( lbt )
    ax.set_xlabel("Lookback Time [Gyr]")
    
    ax.quiver(merger_collection["LOOKBACKTIME"], zeros(length(merger_collection["LOOKBACKTIME"])), 
            x_vec, y_vec,
            color=colormap(proj_vec.*0.5.+0.5), scale_units="y", scale=scale, alpha=0.7, width=0.005, headwidth=2, headlength=3, edgecolor="black", lw=0.5)
    
    ax.invert_xaxis()
    ax.grid()
    clb = fig.colorbar(pltcm.ScalarMappable(norm=pltcolors.Normalize(vmin=-1, vmax=1), cmap="coolwarm"))
    clb.ax.set_ylabel("Z-Axis Component")
    clb.ax.set_yticks([-1,0,1])
    clb.ax.yaxis.set_label_position("left")
    
    ax2 = ax.twiny()
    ax2.set_xlim(ax.get_xlim())
    ax2.set_xticks(ax.get_xticks())
    ax2.set_xticklabels(round.(rs, digits=2)[1:1:end])
    ax2.set_xlabel("z")
    
    if snap
        ax3 = ax.twiny()
        ax3.set_xlim(ax.get_xlim())
        ax3.set_xticks(ax.get_xticks())
        ax3.set_xticklabels(snaps[1:1:end])
        ax3.set_xlabel("Snap")
        ax3.xaxis.set_ticks_position("bottom")
        ax3.xaxis.set_label_position("bottom")
        ax3.spines["bottom"].set_position(("outward", 40))
    end
    
    fig.set_size_inches(plotwidth, plotheight)
    
    fig.savefig(outfile, bbox_inches="tight", pad_inches=.1)
    
    return nothing
end

print("'plot_halostory_arrows'   ")


println()