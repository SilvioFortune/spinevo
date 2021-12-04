
include("/home/moon/sfortune/spinevo/pkg/meta.jl")

fID     = 1
simbox  = "/HydroSims/Magneticum/Box4/uhr_test"
indir   = "/home/moon/sfortune/spinevo/halostories_v20211204_min0.0Gyr"


storyfilelist   = readdir(indir)
halo_filestring = " "
for i in 1:length(storyfilelist)
    if occursin("halo_$(fID)_", storyfilelist[i])
        halo_filestring = storyfilelist[i]
        println(i)
    end
end
merger_collection_STARS = load(joinpath(indir, halo_filestring), "merger_collection_STARS")

gridspace   = 0.6 #Gyr
lbt     = Array{Float64}(undef, 1)
lbt[1]  = merger_collection_STARS["LOOKBACKTIME"][1]
rs      = Array{Float64}(undef, 1)
rs[1]   = merger_collection_STARS["REDSHIFT"][1]
snap    = Array{Int64}(undef, 1)
snap[1] = merger_collection_STARS["SNAP"][1]
for i in 2:length(merger_collection_STARS["SNAP"])-1
    if lbt[end] - merger_collection_STARS["LOOKBACKTIME"][i] > gridspace
        snap    = vcat( snap, merger_collection_STARS["SNAP"        ][i] )
        lbt     = vcat( lbt , merger_collection_STARS["LOOKBACKTIME"][i] )
        rs      = vcat( rs  , merger_collection_STARS["REDSHIFT"    ][i] )
    end
end
snap    = vcat( snap, merger_collection_STARS["SNAP"        ][end] )
lbt     = vcat( lbt , merger_collection_STARS["LOOKBACKTIME"][end] )
rs      = vcat( rs  , merger_collection_STARS["REDSHIFT"    ][end] )


J_norm          = vcat( transpose(norm.(eachcol(merger_collection_STARS["J_main"]))), transpose(norm.(eachcol(merger_collection_STARS["J_main"]))), transpose(norm.(eachcol(merger_collection_STARS["J_main"]))) )

scale           = 1
proj_vec        = (merger_collection_STARS["J_main"]./J_norm)[1,:]
x_vec           = (merger_collection_STARS["J_main"]./J_norm)[2,:]
plotheight      = 3
plotwidth       = 16

pltcm           = pyimport("matplotlib.cm")
pltcolors       = pyimport("matplotlib.colors")
colormap        = pltcm.get_cmap(name="coolwarm")

y = zeros(length(merger_collection_STARS["LOOKBACKTIME"]))
x = zeros(length(merger_collection_STARS["LOOKBACKTIME"]))

fig, ax = subplots()

ax.set_ylim(-1, 1)
ax.set_yticks([-1,0,1])
ax.set_xlim( minimum(merger_collection_STARS["LOOKBACKTIME"])-1, maximum(merger_collection_STARS["LOOKBACKTIME"])+1 )
ax.set_xticks( lbt )
ax.set_xlabel("Lookback Time [Gyr]")


ax.quiver(merger_collection_STARS["LOOKBACKTIME"], zeros(length(merger_collection_STARS["LOOKBACKTIME"])), 
        x_vec, (merger_collection_STARS["J_main"]./J_norm)[3,:],
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

ax3 = ax.twiny()
ax3.set_xlim(ax.get_xlim())
ax3.set_xticks(ax.get_xticks())
ax3.set_xticklabels(snap[1:1:end])
ax3.set_xlabel("Snap")
ax3.xaxis.set_ticks_position("bottom")
ax3.xaxis.set_label_position("bottom")
ax3.spines["bottom"].set_position(("outward", 40))

fig.set_size_inches(plotwidth, plotheight)

fig.savefig(joinpath(@__DIR__, "test.png"), bbox_inches="tight", pad_inches=.1)
