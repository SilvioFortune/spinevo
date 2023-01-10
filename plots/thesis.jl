# pdf-crop-margins -mo 

include("/home/moon/sfortune/spinevo/pkg/meta.jl")

outputdir = "/e/ldata/users/sfortune/Dropbox/Apps/Overleaf/SpinEvolution/figures"
#outputdir = "/home/moon/sfortune/spinevo/plots/thesis"

include("/home/moon/sfortune/spinevo/plots/load_data.jl")

csw03[:start_thr]
csw09[:start_thr]
idxmatch(csw03[:start_thr], csw03[:switches])
idxmatch(csw09[:start_thr], csw09[:switches])
idxclude(csw03[:fakes], idxmatch(csw03[:start_thr], csw03[:switches]))
idxclude(csw09[:fakes], idxmatch(csw09[:start_thr], csw09[:switches]))

#sns = pyimport("seaborn")

# SFR / boxsize
#, fontsize=floor(0.5*scale)
#./ sfr_coll["boxsize"] ./ 1e9
outnm   = joinpath(outputdir,"SFRDcomoving_vs_redshift")
outfile = outnm*".pdf"
PyPlot.matplotlib[:rc]("text", usetex=true)
PyPlot.matplotlib[:rc]("font", family="serif", serif=["Computer Modern"])
pointsize = 5
scale = 0.5
fig, ax = subplots()
ax.axvline(z_madau[idxmax_madau], color="black", label="Madau, z$(latexstring("\$_\\textrm{peak}\$")) = $(round(z_madau[idxmax_madau], digits=1))")
ax.axvline(z_hr[idxmax_hr], color="darkred", label="Box4hr, z$(latexstring("\$_\\textrm{peak}\$")) = $(round(z_hr[idxmax_hr], digits=1))")
ax.axvline(z_uhr[idxmax_uhr], color="red", label="Box4uhr, z$(latexstring("\$_\\textrm{peak}\$")) = $(round(z_uhr[idxmax_uhr], digits=1))")
ax.axvline(global_sfr_coll["redshift"][idxmax_global_sfr_coll], color="mediumblue", label="All subhalos, z$(latexstring("\$_\\textrm{peak}\$")) = $(round(global_sfr_coll["redshift"][idxmax_global_sfr_coll], digits=1))")
ax.axvline(sub1e10_sfr_coll["redshift"][idxmax_sub1e10_sfr_coll], color="dodgerblue", label="All 1e10-cut, z$(latexstring("\$_\\textrm{peak}\$")) = $(round(sub1e10_sfr_coll["redshift"][idxmax_sub1e10_sfr_coll], digits=1))")
ax.axvline(sfr_coll["redshift"][idxmax_sfr_coll], color="cyan", label="Centrals 2e10-cut, z$(latexstring("\$_\\textrm{peak}\$")) = $(round(sfr_coll["redshift"][idxmax_sfr_coll], digits=1))")
#ax.plot( global_sfr_coll["redshift"], log10.(global_sfr_coll["totSFR"] ./ global_sfr_coll["boxsize"] .* 1e9), "-", color="darkred", markersize=5, label="Global")
ax.plot( z_madau, log10.(sfr_madau), "-", color="black")
ax.plot( z_hr, log10.(sfr_hr ./ (48/little_h0)^3), "-", color="darkred")
ax.plot( z_uhr, log10.(sfr_uhr ./ (48/little_h0)^3), "-", color="red")
ax.plot( global_sfr_coll["redshift"], log10.(global_sfr_coll["totSFR"] ./ (48/little_h0)^3), ".-", markersize=2, color="mediumblue")
ax.plot( sub1e10_sfr_coll["redshift"], log10.(sub1e10_sfr_coll["totSFR"] ./ (48/little_h0)^3), ".-", markersize=2, color="dodgerblue")
ax.plot( sfr_coll["redshift"], log10.(sfr_coll["totSFR"] ./ (48/little_h0)^3), ".-", markersize=2, color="cyan")
#ax.set_title("$(latexstring("\$\\Delta\$"))t $(latexstring("\$\\sim\$")) 0.3 Gyr")
ax.set_xlim([0,8.5])
ax.set_ylim([-2,-0.7])
ax.set_xlabel("z")
ax.set_xscale("linear")
ax.set_xticks([1,2,3,4,5,6,7,8])
ax.set_xticklabels([1,2,3,4,5,6,7,8])
ax.set_ylabel("log( SFRD [$(latexstring("\\frac{M_\\odot}{yr\\ (Mpc/h)^3}"))] )")
ax.legend(loc="upper right", frameon=false, borderpad=1, handlelength=1.8)
ax.grid(false)
fig.set_size_inches(16scale, 9scale)
fig.tight_layout()
fig.savefig(outfile, bbox_inches="tight", pad_inches=.1)
close("all")
PyPlot.matplotlib[:rc]("text", usetex=false)




# b-value vs SFR at different redshift bins
outnm   = joinpath(outputdir,"sSFR_vs_bval_zsplit")
outfile = outnm*".pdf"
x   = as03["BVAL_0"]
y   = log10.(as03["SFR"] ./ as03["M"] )
miny= minimum(pyplottable(y)[findall(x->x .!== NaN, pyplottable(y))])
maxy= maximum(pyplottable(y)[findall(x->x .!== NaN, pyplottable(y))])
minx= minimum(pyplottable(x)[findall(x->x .!== NaN, pyplottable(x))])
maxx= maximum(pyplottable(x)[findall(x->x .!== NaN, pyplottable(x))])
slcts       = Dict{Int64, Vector{Int64}}()
slcts[1]    = 1:length(as03["BVAL_0"])
slcts[2]    = csw03[:z_md]
slcts[3]    = csw03[:z_hi]
slcts[4]    = csw03[:z_lo]
labels      = Dict{Int64, String}()
labels[1]   = "Full Set, N = $(length(slcts[1]))"
labels[2]   = "$(round(minimum(as03["redshift"][slcts[3]]), digits=1))  $(latexstring("\$>\$")) z  $(latexstring("\$>\$")) $(round(maximum(as03["redshift"][slcts[4]]), digits=1)), N = $(length(slcts[2]))"
labels[3]   = "$(round(maximum(as03["redshift"][slcts[3]]), digits=1)) $(latexstring("\$\\geq\$")) z $(latexstring("\$\\geq\$")) $(round(minimum(as03["redshift"][slcts[3]]), digits=1)), N = $(length(slcts[3]))"
labels[4]   = "$(round(maximum(as03["redshift"][slcts[4]]), digits=1)) $(latexstring("\$\\geq\$")) z, N = $(length(slcts[4]))"
labpos      = Array{String}(undef, 4)
labpos     .= "upper left"
plot_hexbins_nxm( x, 
            y, 
            slcts, labels, 2, 2,
            #scale=0.0,
            cbarshift=2.5,
            label_pos=labpos, 
            outfile=outfile, 
            lognorm=true, gridres=(15,5), 
            ymin=miny, ymax=-7.5,
            xmin=minx, xmax=maxx, 
            ylabel="log( sSFR [$(latexstring("\\frac{1}{yr}"))] )",
            xlabel="b-value", fontsize=30,
            binned_median="y", plot_bval="x",
            )
#





gridspace   = 2 # Gyr
scale       = 0.5
outfile     = joinpath(outputdir,"Nsub_vs_lbt_z.pdf")
lbt     = Array{Float64}(undef, 1)
lbt[1]  = global_sfr_coll["lookbacktime"][1]
rs      = Array{Float64}(undef, 1)
rs[1]   = global_sfr_coll["redshift"][1]
for i in 2:length(global_sfr_coll["lookbacktime"])-1
        println(lbt[end] - global_sfr_coll["lookbacktime"][i], "   ", global_sfr_coll["lookbacktime"][i] - global_sfr_coll["lookbacktime"][end])
    if lbt[end] - global_sfr_coll["lookbacktime"][i] > gridspace && global_sfr_coll["lookbacktime"][i] - global_sfr_coll["lookbacktime"][end] > gridspace
        global lbt     = vcat( lbt , global_sfr_coll["lookbacktime"][i] )
        global rs      = vcat( rs  , global_sfr_coll["redshift"    ][i] )
    end
end
lbt     = vcat( lbt , global_sfr_coll["lookbacktime"][end] )
rs      = vcat( rs  , global_sfr_coll["redshift"    ][end] )

fig, ax = subplots()
PyPlot.matplotlib[:rc]("text", usetex=true)
PyPlot.matplotlib[:rc]("font", family="serif", serif=["Computer Modern"])
ax.plot( global_sfr_coll["lookbacktime"], global_sfr_coll["N_halos"], ".-", linewidth=1.5, color="darkblue", label="All Subhalos")
ax.plot( sfr_coll["lookbacktime"], sfr_coll["N_halos"], ".-", linewidth=1.5, color="mediumblue", label="This Sample")
ax.plot( sthr_coll["lookbacktime"], sthr_coll["N_halos"], ".-", linewidth=1.5, color="cyan", label="Transition Sample")
ax.invert_xaxis()
#ax.hist( pyplottable(as03["lookbacktime"])           , bins=100, linewidth=1.5, rwidth=1, alpha=0.6, histtype="step", color="mediumblue")
#ax.hist( pyplottable(as03["lookbacktime"][start_thr]), bins=200, linewidth=1.5, rwidth=1, alpha=0.6, histtype="step", color="darkred", label="Current \\& Past" )
ax.set_xlabel("Lookback Time [Gyr]")
ax.set_ylabel("Number of Subhalos")
ax.set_yscale("log")
ax.set_xscale("linear")
#ax.set_xlim([10, nothing])
ax.legend(loc="lower right", frameon=true, borderpad=1, handlelength=1.8)
ax.grid()
ax.set_xticks( round.(lbt,digits=1) )

ax2 = ax.twiny()
ax2.set_xlim(ax.get_xlim())
ax2.set_xticks(ax.get_xticks())
ax2.set_xticklabels(round.(rs, digits=2))
ax2.set_xlabel("z")

fig.set_size_inches(9scale, 9scale)
fig.tight_layout()
fig.savefig(outfile, bbox_inches="tight", pad_inches=.1)
close("all")
PyPlot.matplotlib[:rc]("text", usetex=false)




gridspace   = 2 # Gyr
scale       = 0.5
outfile     = joinpath(outputdir,"Mstars_vs_lbt_z.pdf")
lbt     = Array{Float64}(undef, 1)
lbt[1]  = global_sfr_coll["lookbacktime"][1]
rs      = Array{Float64}(undef, 1)
rs[1]   = global_sfr_coll["redshift"][1]
for i in 2:length(global_sfr_coll["lookbacktime"])-1
        println(lbt[end] - global_sfr_coll["lookbacktime"][i], "   ", global_sfr_coll["lookbacktime"][i] - global_sfr_coll["lookbacktime"][end])
    if lbt[end] - global_sfr_coll["lookbacktime"][i] > gridspace && global_sfr_coll["lookbacktime"][i] - global_sfr_coll["lookbacktime"][end] > gridspace
        global lbt     = vcat( lbt , global_sfr_coll["lookbacktime"][i] )
        global rs      = vcat( rs  , global_sfr_coll["redshift"    ][i] )
    end
end
lbt     = vcat( lbt , global_sfr_coll["lookbacktime"][end] )
rs      = vcat( rs  , global_sfr_coll["redshift"    ][end] )

fig, ax = subplots()
PyPlot.matplotlib[:rc]("text", usetex=true)
PyPlot.matplotlib[:rc]("font", family="serif", serif=["Computer Modern"])
ax.plot( global_sfr_coll["lookbacktime"], global_sfr_coll["totM"], ".-", linewidth=1.5, color="darkblue", label="All Subhalos")
ax.plot( sfr_coll["lookbacktime"], sfr_coll["totM"], ".-", linewidth=1.5, color="mediumblue", label="This Sample")
ax.plot( sthr_coll["lookbacktime"], sthr_coll["totM"], ".-", linewidth=1.5, color="cyan", label="Transition Sample")
ax.invert_xaxis()
#ax.hist( pyplottable(as03["lookbacktime"])           , bins=100, linewidth=1.5, rwidth=1, alpha=0.6, histtype="step", color="mediumblue")
#ax.hist( pyplottable(as03["lookbacktime"][start_thr]), bins=200, linewidth=1.5, rwidth=1, alpha=0.6, histtype="step", color="darkred", label="Current \\& Past" )
ax.set_xlabel("Lookback Time [Gyr]")
ax.set_ylabel("Sum of Subhalo M$(latexstring("\$_*\$")) [M$(latexstring("\$_\\odot\$"))]")
ax.set_yscale("log")
ax.set_xscale("linear")
ax.set_ylim([1e11, nothing])
ax.legend(loc="lower right", frameon=true, borderpad=1, handlelength=1.8)
ax.grid()
ax.set_xticks( round.(lbt,digits=1) )

ax2 = ax.twiny()
ax2.set_xlim(ax.get_xlim())
ax2.set_xticks(ax.get_xticks())
ax2.set_xticklabels(round.(rs, digits=2))
ax2.set_xlabel("z")

fig.set_size_inches(9scale, 9scale)
fig.tight_layout()
fig.savefig(outfile, bbox_inches="tight", pad_inches=.1)
close("all")
PyPlot.matplotlib[:rc]("text", usetex=false)




# mass diff stars vs merger mass vs flipangle color
PyPlot.matplotlib[:rc]("text", usetex=true)
PyPlot.matplotlib[:rc]("font", family="serif", serif=["Computer Modern"])
pointsize = 10
scale = 0.5
fig, ax = subplots()
colr    = pyplottable(as03["ϕ_flip"][csw03[:transition]])
order   = sortperm(colr)
colr    = colr[order]
y       = replace( symlog10.(as03["M_MERGERS"][csw03[:transition]]), missing => 0. )[order]
x       = replace( symlog10.(as03["ΔM"][csw03[:transition]]), missing => 0. )[order]
#slct    = findcs(x, gt=0., comparewith=findall(dummy->!isnan.(dummy), colr))
#slct    = findcs(y, gt=0., comparewith=slct)
maxcolr = maximum(colr)
#slct    = idxmatch(slct2, n_switches)
p = ax.scatter( x, y,
        zorder=2   , s= pointsize    , alpha=0.7,#colr[slct] ./ maxcolr, (pointsize/maxcolr) .* colr[slct]
        cmap="plasma_r", c=colr, rasterized=true) # brg_r  
ax.plot([8,12], [8,12], "-",color="black", linewidth=0.5, rasterized=true)
colbar  = fig.colorbar(p, ax=ax)#, ticks=[minimum(as03["BVAL_0"    ]), b_ell, b_disk, maximum(as03["BVAL_0"    ])])
#colbar.ax.set_yticklabels(["$(round(minimum(as03["BVAL_0"    ]), digits=1))", "ETG", "LTG", "$(round(maximum(as03["BVAL_0"    ]), digits=1))"])
colbar.ax.tick_params(axis="y")
colbar.ax.set_ylabel("Flip Angle [°]")
colbar.ax.yaxis.set_label_position("left")
#ax.legend(loc="lower right", frameon=true, borderpad=1, handlelength=1.8)
ax.set_xlim([5.5,11])
ax.set_ylim([8,11])
#ax.set_title("N = $(length(slct))")
#ax.set_aspect("equal")
ax.set_xlabel("symlog( $(latexstring("\$\\Delta\$"))M$(latexstring("\$_*\$")))")
ax.set_ylabel("log( M$(latexstring("\$_\\textrm{All Mergers}\$")) )")
ax.grid()
#ax.autoscale()
fig.set_size_inches(11*scale, 9*scale)
fig.tight_layout()
fig.savefig(joinpath(outputdir, "M_Mergers_vs_delM_g0_lores.pdf"), bbox_inches="tight", pad_inches=.1)
close("all")
PyPlot.matplotlib[:rc]("text", usetex=false)

# negatives
PyPlot.matplotlib[:rc]("text", usetex=true)
PyPlot.matplotlib[:rc]("font", family="serif", serif=["Computer Modern"])
pointsize = 10
scale = 0.5
fig, ax = subplots()
colr    = pyplottable(as03["ϕ_flip"][csw03[:transition]])
order   = sortperm(colr)
colr    = colr[order]
y       = replace( symlog10.(as03["M_MERGERS"][csw03[:transition]]), missing => 0. )[order]
x       = replace( symlog10.(as03["ΔM"][csw03[:transition]]), missing => 0. )[order]
#slct    = findcs(x, gt=0., comparewith=findall(dummy->!isnan.(dummy), colr))
#slct    = findcs(y, gt=0., comparewith=slct)
maxcolr = maximum(colr)
#slct    = idxmatch(slct2, n_switches)
p = ax.scatter( x, y,
        zorder=2   , s= pointsize    , alpha=0.7,#colr[slct] ./ maxcolr, (pointsize/maxcolr) .* colr[slct]
        cmap="plasma_r", c=colr, rasterized=true) # brg_r  
ax.plot([-8,-12], [8,12], "-",color="black", linewidth=0.5, rasterized=true)
colbar  = fig.colorbar(p, ax=ax)#, ticks=[minimum(as03["BVAL_0"    ]), b_ell, b_disk, maximum(as03["BVAL_0"    ])])
#colbar.ax.set_yticklabels(["$(round(minimum(as03["BVAL_0"    ]), digits=1))", "ETG", "LTG", "$(round(maximum(as03["BVAL_0"    ]), digits=1))"])
colbar.ax.tick_params(axis="y")
colbar.ax.set_ylabel("Flip Angle [°]")
colbar.ax.yaxis.set_label_position("left")
#ax.legend(loc="lower right", frameon=true, borderpad=1, handlelength=1.8)
ax.set_xlim([-11,-5.5])
ax.set_ylim([8,11])
#ax.set_title("N = $(length(slct))")
#ax.set_aspect("equal")
ax.set_xlabel("symlog( $(latexstring("\$\\Delta\$"))M$(latexstring("\$_*\$")))")
ax.set_ylabel("log( M$(latexstring("\$_\\textrm{All Mergers}\$")) )")
ax.grid()
#ax.autoscale()
fig.set_size_inches(11*scale, 9*scale)
fig.tight_layout()
fig.savefig(joinpath(outputdir, "M_Mergers_vs_delM_l0_lores.pdf"), bbox_inches="tight", pad_inches=.1)
close("all")
PyPlot.matplotlib[:rc]("text", usetex=false)









############################################################################################
# Hirschmann 2014

outfile     = joinpath(outputdir,"stellarmassfct_Hirschmann2014like.pdf")
scale=0.5
N_bins=10
pltcm       = pyimport("matplotlib.cm")
pltcolors   = pyimport("matplotlib.colors")
colormap    = pltcm.get_cmap(name="gnuplot")
fig, ax = subplots()
n_0, x_0      = PyPlot.hist( log10.(pyplottable(as03["M"][csw03[:z_hirsch_0] ])), bins=N_bins, density=true, log=false)
n_05,x_05     = PyPlot.hist( log10.(pyplottable(as03["M"][csw03[:z_hirsch_05]])), bins=N_bins, density=true, log=false)
n_1, x_1      = PyPlot.hist( log10.(pyplottable(as03["M"][csw03[:z_hirsch_1 ]])), bins=N_bins, density=true, log=false)
n_2, x_2      = PyPlot.hist( log10.(pyplottable(as03["M"][csw03[:z_hirsch_2 ]])), bins=N_bins, density=true, log=false)
n_3, x_3      = PyPlot.hist( log10.(pyplottable(as03["M"][csw03[:z_hirsch_3 ]])), bins=N_bins, density=true, log=false)
n_4, x_4      = PyPlot.hist( log10.(pyplottable(as03["M"][csw03[:z_hirsch_4 ]])), bins=N_bins, density=true, log=false)
bcenters_0     = 0.5 .* (x_0[2:end]  .+ x_0[1:end-1])
bcenters_05    = 0.5 .* (x_05[2:end] .+ x_05[1:end-1])
bcenters_1     = 0.5 .* (x_1[2:end]  .+ x_1[1:end-1])
bcenters_2     = 0.5 .* (x_2[2:end]  .+ x_2[1:end-1])
bcenters_3     = 0.5 .* (x_3[2:end]  .+ x_3[1:end-1])
bcenters_4     = 0.5 .* (x_4[2:end]  .+ x_4[1:end-1])
clf()
fig, ax = subplots()
PyPlot.matplotlib[:rc]("text", usetex=true)
PyPlot.matplotlib[:rc]("font", family="serif", serif=["Computer Modern"])
ax.grid(true)
ax.plot(bcenters_0 , log10.(n_0 .* (length(csw03[:z_hirsch_0])/(48/little_h0)^3) ./ length(unique(as03["redshift"][csw03[:z_hirsch_0]]))) , ".-", ms=10, label="z $(latexstring("\$\\in\$"))[$(round(minimum(as03["redshift"][csw03[:z_hirsch_0]]),digits=2)), $(round(maximum(as03["redshift"][csw03[:z_hirsch_0]]),digits=2 ))]         ", color=colormap(0/4)   )
ax.plot(bcenters_05 ,log10.(n_05.*(length(csw03[:z_hirsch_05])/(48/little_h0)^3) ./ length(unique(as03["redshift"][csw03[:z_hirsch_05]]))) ,".-", ms=10, label="z $(latexstring("\$\\in\$"))[$(round(minimum(as03["redshift"][csw03[:z_hirsch_05]]),digits=2)),$(round(maximum(as03["redshift"][csw03[:z_hirsch_05]]),digits=2))]   "      , color=colormap(0.5/4)   )
ax.plot(bcenters_1 , log10.(n_1 .* (length(csw03[:z_hirsch_1])/(48/little_h0)^3) ./ length(unique(as03["redshift"][csw03[:z_hirsch_1]]))) , ".-", ms=10, label="z $(latexstring("\$\\in\$"))[$(round(minimum(as03["redshift"][csw03[:z_hirsch_1]]),digits=2)), $(round(maximum(as03["redshift"][csw03[:z_hirsch_1]]),digits=2 ))]         ", color=colormap(1/4)   )
ax.plot(bcenters_2 , log10.(n_2 .* (length(csw03[:z_hirsch_2])/(48/little_h0)^3) ./ length(unique(as03["redshift"][csw03[:z_hirsch_2]]))) , ".-", ms=10, label="z $(latexstring("\$\\in\$"))[$(round(minimum(as03["redshift"][csw03[:z_hirsch_2]]),digits=2)), $(round(maximum(as03["redshift"][csw03[:z_hirsch_2]]),digits=2 ))]         ", color=colormap(2/4)   )
ax.plot(bcenters_3 , log10.(n_3 .* (length(csw03[:z_hirsch_3])/(48/little_h0)^3) ./ length(unique(as03["redshift"][csw03[:z_hirsch_3]]))) , ".-", ms=10, label="z $(latexstring("\$\\in\$"))[$(round(minimum(as03["redshift"][csw03[:z_hirsch_3]]),digits=2)), $(round(maximum(as03["redshift"][csw03[:z_hirsch_3]]),digits=2 ))]         ", color=colormap(3/4)   )
ax.plot(bcenters_4 , log10.(n_4 .* (length(csw03[:z_hirsch_4])/(48/little_h0)^3) ./ length(unique(as03["redshift"][csw03[:z_hirsch_4]]))) , ".-", ms=10, label="z $(latexstring("\$\\in\$"))[$(round(minimum(as03["redshift"][csw03[:z_hirsch_4]]),digits=2)), $(round(maximum(as03["redshift"][csw03[:z_hirsch_4]]),digits=2 ))]         ", color=colormap(4/4)   )
ax.set_xlabel("log( M$(latexstring("\$_*\$")) [ M$(latexstring("\$_\\odot\$")) ] )")
ax.set_ylabel("log( $(latexstring("\$\\Phi\$")) [ $(latexstring("\$\\frac{h^3}{(Mpc^3)}\$")) ] )")
ax.set_yscale("linear")
ax.set_ylim(-5.5,-1)
ax.set_xlim(8.5,12.5)
#ax.set_title("$(latexstring("\$\\Delta\$"))t $(latexstring("\$\\sim\$")) 0.3 Gyr within $(round(maximum(as["redshift"]),digits=2)) $(latexstring("\$\\geq\$")) z $(latexstring("\$\\geq\$")) $(round(minimum(as["redshift"]),digits=2))")
ax.legend(loc="lower left", frameon=true, borderpad=1, handlelength=1.8)
fig.set_size_inches(10scale, 9scale)
fig.tight_layout()
fig.savefig(outfile, bbox_inches="tight", pad_inches=.1)
close("all")
PyPlot.matplotlib[:rc]("text", usetex=false)
println("z $(latexstring("\$\\in\$"))[$(round(minimum(as03["redshift"][csw03[:z_hirsch_0]]),digits=2)), $(round(maximum(as03["redshift"][csw03[:z_hirsch_0]]),digits=2 ))], N = $(length(csw03[:z_hirsch_0] ))         ")
println("z $(latexstring("\$\\in\$"))[$(round(minimum(as03["redshift"][csw03[:z_hirsch_05]]),digits=2)),$(round(maximum(as03["redshift"][csw03[:z_hirsch_05]]),digits=2))], N = $(length(csw03[:z_hirsch_05]))   "      )
println("z $(latexstring("\$\\in\$"))[$(round(minimum(as03["redshift"][csw03[:z_hirsch_1]]),digits=2)), $(round(maximum(as03["redshift"][csw03[:z_hirsch_1]]),digits=2 ))], N = $(length(csw03[:z_hirsch_1] ))         ")
println("z $(latexstring("\$\\in\$"))[$(round(minimum(as03["redshift"][csw03[:z_hirsch_2]]),digits=2)), $(round(maximum(as03["redshift"][csw03[:z_hirsch_2]]),digits=2 ))], N = $(length(csw03[:z_hirsch_2] ))         ")
println("z $(latexstring("\$\\in\$"))[$(round(minimum(as03["redshift"][csw03[:z_hirsch_3]]),digits=2)), $(round(maximum(as03["redshift"][csw03[:z_hirsch_3]]),digits=2 ))], N = $(length(csw03[:z_hirsch_3] ))         ")
println("z $(latexstring("\$\\in\$"))[$(round(minimum(as03["redshift"][csw03[:z_hirsch_4]]),digits=2)), $(round(maximum(as03["redshift"][csw03[:z_hirsch_4]]),digits=2 ))], N = $(length(csw03[:z_hirsch_4] ))         ")




###########################################
# MfromJ vs M subfind
outnm = joinpath(outputdir,"MfromJ_vs_Msubfind_vs_bval_lores")
outfile     = outnm*".pdf"
scale = 0.5
PyPlot.matplotlib[:rc]("text", usetex=true)
PyPlot.matplotlib[:rc]("font", family="serif", serif=["Computer Modern"])
pointsize = 2
fig, ax = subplots()
#colr    = as03["SFR"] ./ as03["M"] .* 1e9
colr    = as03["BVAL_0"]
nomissings = idxclude(1:length(colr), findcs(colr, eq=missing))
colr    = colr[nomissings]
order   = sortperm(colr)#[length(colr):-1:1]
colr    = colr[order]
#y       = pyplottable(delJ[:,order] .- Jom[:,order])
y       = pyplottable(as03["M_fromJ"])
x       = pyplottable(as03["M"])
x    = x[nomissings]
y    = y[nomissings]
maxcolr = maximum(colr)
#slct    = idxmatch(slct2, n_switches)
ax.grid()
p = ax.scatter( log10.(x[order]), log10.(y[order]),
        zorder=2   , s= pointsize    , alpha=0.5,#colr ./ maxcolr, (pointsize/maxcolr) .* colr
        cmap=ColorMap(get_bvalue_cmap().colors), c=colr, vmin=-6.1,vmax=-3.9, rasterized=true) # plasma_r  
#
ax.plot([10,13], [10,13], "-",color="black", linewidth=0.5, rasterized=true)
colbar  = fig.colorbar(p, ax=ax)#, ticks=[minimum(as03["BVAL_0"    ]), b_ell, b_disk, maximum(as03["BVAL_0"    ])])
#colbar.ax.set_yticklabels(["$(round(minimum(as03["BVAL_0"    ]), digits=1))", "ETG", "LTG", "$(round(maximum(as03["BVAL_0"    ]), digits=1))"])
colbar.ax.tick_params(axis="y")
#colbar.ax.set_ylabel("sSFR [Gyr$(latexstring("\$^{-1}\$"))]")
colbar.ax.set_ylabel("b-value")
colbar.ax.yaxis.set_label_position("left")
#ax.legend(loc="lower right", frameon=true, borderpad=1, handlelength=1.8)
ax.set_ylim(extrema(log10.(y[order])))
ax.set_xlim(extrema(log10.(x[order])))
#ax.set_title("immediate read, all mergers, 03 Gyr time step")
#ax.set_aspect("equal")\\textrm{orbital,Mergers}
ax.set_ylabel("log( $(latexstring("\$M_\\textrm{read}\$")) [ $(latexstring("\$M_\\odot\$")) ] )")
ax.set_xlabel("log( $(latexstring("\$M_\\textrm{subfind}\$")) [ $(latexstring("\$M_\\odot\$")) ] )")
ax.grid(true)
#ax.autoscale()
fig.set_size_inches(11*scale, 9*scale)
fig.tight_layout()
fig.savefig(outfile, bbox_inches="tight", pad_inches=.1)
close("all")
PyPlot.matplotlib[:rc]("text", usetex=false)

# hexbins
outfile     = outnm*"_hexbins.pdf"
y       = log10.(pyplottable(as03["M_fromJ"]))
x       = log10.(pyplottable(as03["M"]))
plot_hexbins( x, y, 
    outfile=outfile,
    scale=15, gridres=(50,23), xscale="linear", yscale="linear", 
    xlabel="log( $(latexstring("\$M_\\textrm{subfind}\$")) [ $(latexstring("\$M_\\odot\$")) ] )", 
    ylabel="log( $(latexstring("\$M_\\textrm{read}\$")) [ $(latexstring("\$M_\\odot\$")) ] )", 
    grid=true,
    #xmin=0, xmax=180,# ymin=ymin, ymax=ymax,
    lognorm=true,
    binned_median = "x", plot_bval=" ", plot_line=hcat(LinRange(0, 20, 1000), LinRange(0, 20, 1000)),
    calc_norm=true,    # @pyplottable
    )
#



##########################################################################
# M hist
outfile     = joinpath(outputdir,"hist_M.pdf")
scale = 0.5
lw=2
PyPlot.matplotlib[:rc]("text", usetex=true)
PyPlot.matplotlib[:rc]("font", family="serif", serif=["Computer Modern"])
fig, ax = subplots()
ax.hist( pyplottable(log10.(pyplottable(ad03["M"])))    , bins=20, label="$(latexstring("\$M_\\textrm{DM,centrals}\$"))", histtype="step", linewidth=lw    , rwidth=1    , alpha=1, color="darkblue" )
ax.hist( pyplottable(log10.(pyplottable(as03["M"])))    , bins=20, label="$(latexstring("\$M_{*\\textrm{,centrals}}\$"))", histtype="step", linewidth=lw    , rwidth=1    , alpha=1, color="mediumblue" )
ax.hist( pyplottable(log10.(pyplottable(ad03["merger map"][23,:])))    , bins=20, label="$(latexstring("\$M_\\textrm{DM,mergers}\$"))", histtype="step", linewidth=lw    , rwidth=1    , alpha=1, color="seagreen" )
ax.hist( pyplottable(log10.(pyplottable(as03["merger map"][23,:])))    , bins=20, label="$(latexstring("\$M_{*\\textrm{,mergers}}\$"))", histtype="step", linewidth=lw    , rwidth=1    , alpha=1, color="mediumseagreen" )
ax.set_xlabel("log( Mass [M$(latexstring("\$M_\\odot\$")) ] )")
ax.set_yscale("log")
#ax.set_xlim([10, nothing])
ax.legend(loc="upper left", frameon=false, borderpad=1, handlelength=1.8)
ax.grid()
scale=0.5
fig.set_size_inches(11scale, 9scale)
fig.tight_layout()
fig.savefig(outfile, bbox_inches="tight", pad_inches=.1)
close("all")
PyPlot.matplotlib[:rc]("text", usetex=false)

outfile     = joinpath(outputdir,"hist_M_pre.pdf")
scale = 0.5
lw=2
PyPlot.matplotlib[:rc]("text", usetex=true)
PyPlot.matplotlib[:rc]("font", family="serif", serif=["Computer Modern"])
fig, ax = subplots()
ax.hist( pyplottable(log10.(pyplottable(as03["M"][csw03[:transition]] .-as03["ΔM"][csw03[:transition]])))    , bins=50, label="M$(latexstring("\$_\\textrm{*,pre}\$"))", histtype="step", linewidth=lw, alpha=1, color="mediumblue" )
ax.axvline(pctl68up_logMstarpre, color="black",linewidth=lw)
ax.axvline(pctl68lo_logMstarpre, color="black",linewidth=lw)
ax.set_xlabel("log( $(latexstring("\$M_{*\\textrm{,pre}}\$")) [ $(latexstring("\$M_\\odot\$")) ] )")
ax.set_yscale("log")
#ax.set_xlim([10, nothing])
ax.legend(loc="upper right", frameon=false, borderpad=1, handlelength=1.8)
ax.grid()
scale=0.5
fig.set_size_inches(11scale, 9scale)
fig.tight_layout()
fig.savefig(outfile, bbox_inches="tight", pad_inches=.1)
close("all")
PyPlot.matplotlib[:rc]("text", usetex=false)





##################################################
# J hist / j hist
outfile     = joinpath(outputdir,"hist_j.pdf")
scale = 0.5
lw=2
PyPlot.matplotlib[:rc]("text", usetex=true)
PyPlot.matplotlib[:rc]("font", family="serif", serif=["Computer Modern"])
x1 = pyplottable(log10.(pyplottable(as03["j_main"][:,csw03[:result_ell]])))
x2 = pyplottable(log10.(pyplottable(as03["j_main"][:,csw03[:result_int]])))
x3 = pyplottable(log10.(pyplottable(as03["j_main"][:,csw03[:result_disk]])))
fig, ax = subplots()
ax.axvline(jlim1, color="black")
ax.axvline(jlim2, color="black")
ax.hist( x1    , bins=50, label="Ellipticals, N = $(length(x1)-count(isnan,x1))", histtype="step", linewidth=lw    , rwidth=1    , alpha=1, color="darkred" )
ax.hist( x2    , bins=50, label="Intermediates, N = $(length(x2)-count(isnan,x2))", histtype="step", linewidth=lw    , rwidth=1    , alpha=1, color="orange" )
ax.hist( x3    , bins=50, label="Disks, N = $(length(x3)-count(isnan,x3))", histtype="step", linewidth=lw    , rwidth=1    , alpha=1, color="mediumblue" )
ax.set_xlabel("log( j$(latexstring("\$_*\$")) [ kpc km/s ] )")
#ax.set_xscale("log")
#ax.set_xlim([10, nothing])
ax.legend(loc="upper right", frameon=true, borderpad=1, handlelength=1.8)
ax.grid()
scale=0.5
fig.set_size_inches(10scale, 9scale)
fig.tight_layout()
fig.savefig(outfile, bbox_inches="tight", pad_inches=.1)
close("all")
PyPlot.matplotlib[:rc]("text", usetex=false)


outfile     = joinpath(outputdir,"hist_J_uppercase.pdf")
scale = 0.5
lw=2
PyPlot.matplotlib[:rc]("text", usetex=true)
PyPlot.matplotlib[:rc]("font", family="serif", serif=["Computer Modern"])
x1 = pyplottable(log10.(pyplottable(as03["J_main"][:,csw03[:result_ell]])))
x2 = pyplottable(log10.(pyplottable(as03["J_main"][:,csw03[:result_int]])))
x3 = pyplottable(log10.(pyplottable(as03["J_main"][:,csw03[:result_disk]])))
fig, ax = subplots()
ax.axvline(Jlim1, color="black")
ax.axvline(Jlim2, color="black")
ax.hist( x1    , bins=50, label="Ellipticals, N = $(length(x1)-count(isnan,x1))", histtype="step", linewidth=lw    , rwidth=1    , alpha=1, color="darkred" )
ax.hist( x2    , bins=50, label="Intermediates, N = $(length(x2)-count(isnan,x2))", histtype="step", linewidth=lw    , rwidth=1    , alpha=1, color="orange" )
ax.hist( x3    , bins=50, label="Disks, N = $(length(x3)-count(isnan,x3))", histtype="step", linewidth=lw    , rwidth=1    , alpha=1, color="mediumblue" )
ax.set_xlabel("log( $(latexstring("\$J_*\$")) [ $(latexstring("\$M_\\odot\$")) kpc km/s ] )")
#ax.set_xscale("log")
#ax.set_xlim([10, nothing])
ax.legend(loc="upper right", frameon=true, borderpad=1, handlelength=1.8)
ax.grid()
scale=0.5
fig.set_size_inches(10scale, 9scale)
fig.tight_layout()
fig.savefig(outfile, bbox_inches="tight", pad_inches=.1)
close("all")
PyPlot.matplotlib[:rc]("text", usetex=false)





####################################################################################
# m star - j star 
z2 = findcs(as03["redshift"], leq=2.5, geq=1.5)
z0 = findcs(as03["redshift"], eq=minimum(as03["redshift"]))

outfile     = joinpath(outputdir,"j_vs_M_z02_lores.pdf")
scale = 0.5
bins  = 10
PyPlot.matplotlib[:rc]("text", usetex=true)
PyPlot.matplotlib[:rc]("font", family="serif", serif=["Computer Modern"])
pointsize = 2
lw        = 2
fig, ax = subplots()
#y       = pyplottable(delJ[:,order] .- Jom[:,order])
y0       = pyplottable(as03["j_main"][:,z0])
x0       = pyplottable(as03["M"][z0])
y2       = pyplottable(as03["j_main"][:,z2])
x2       = pyplottable(as03["M"][z2])
#slct    = idxmatch(slct2, n_switches)
ax.grid(false)
ax.scatter( log10.(x0), log10.(y0),
        s= pointsize    , alpha=0.5,#colr ./ maxcolr, (pointsize/maxcolr) .* colr
        color="mediumblue", label="z $(latexstring("\$\\approx\$")) 0, N = $(length(x0))", rasterized=true) # brg_r  
ax.scatter( log10.(x2), log10.(y2),
        s= pointsize    , alpha=0.5,#colr ./ maxcolr, (pointsize/maxcolr) .* colr
        color="darkred", label="z $(latexstring("\$\\approx\$")) 2, N = $(length(x2))", rasterized=true) # brg_r  
ax.legend(loc="lower right", frameon=true, borderpad=1, handlelength=1.8)

medians0, edges0, bin_number0     = stats.binned_statistic(log10.(x0), log10.(y0), statistic="median", bins=bins)
bin_width0   = edges0[2] - edges0[1]
centers0     = edges0[2:end] .- (bin_width0/2)
pctl_680, edges_up0, bin_number_up0     = stats.binned_statistic(log10.(x0), log10.(y0), statistic=pctl68upper, bins=bins)
pctl_320, edges_lo0, bin_number_lo0     = stats.binned_statistic(log10.(x0), log10.(y0), statistic=pctl68lower, bins=bins)
ax.fill_between(centers0, pctl_320, pctl_680, color="mediumblue", alpha=0.2, rasterized=true)
ax.plot(centers0, medians0, "-", color="mediumblue", linewidth=lw, rasterized=true)

medians2, edges2, bin_number2     = stats.binned_statistic(log10.(x2), log10.(y2), statistic="median", bins=bins)
bin_width2   = edges2[2] - edges2[1]
centers2     = edges2[2:end] .- (bin_width2/2)
pctl_682, edges_up2, bin_number_up2     = stats.binned_statistic(log10.(x2), log10.(y2), statistic=pctl68upper, bins=bins)
pctl_322, edges_lo2, bin_number_lo2     = stats.binned_statistic(log10.(x2), log10.(y2), statistic=pctl68lower, bins=bins)
ax.fill_between(centers2, pctl_322, pctl_682, color="darkred", alpha=0.2, rasterized=true)
ax.plot(centers2, medians2, "-", color="darkred", linewidth=lw, rasterized=true)

#ax.set_ylim(extrema(log10.(y[order])))
#ax.set_xlim(extrema(log10.(x[order])))
#ax.set_title("immediate read, all mergers, 03 Gyr time step")
#ax.set_aspect("equal")\\textrm{orbital,Mergers}
ax.set_ylabel("log( j$(latexstring("\$_*\$")) [ kpc km/s ] )")
ax.set_xlabel("log( M$(latexstring("\$_*\$")) [ M$(latexstring("\$_\\odot\$")) ] )")
#ax.grid(true)
#ax.autoscale()
fig.set_size_inches(11*scale, 9*scale)
fig.tight_layout()
fig.savefig(outfile, bbox_inches="tight", pad_inches=.1)
close("all")
PyPlot.matplotlib[:rc]("text", usetex=false)



# sfc:
# scatter j end vs j start vs flip angle or something with mass

outfile     = joinpath(outputdir,"jpost_v_jpre_vs_phi.pdf")
scale = 0.5
PyPlot.matplotlib[:rc]("text", usetex=true)
PyPlot.matplotlib[:rc]("font", family="serif", serif=["Computer Modern"])
pointsize = 2
fig, ax = subplots()
slct    = csw03[:transition]#idxclude(csw03[:start_thr], csw03[:switches])
y       = pyplottable(as03["j_main"][:,slct])
x       = pyplottable(as03["j_main"][:,slct] .- as03["Δj_main"][:,slct])
colr    = as03["ϕ_flip"][slct]
nomissings = idxclude(1:length(colr), findcs(colr, eq=missing))
colr    = colr[nomissings]
x    = x[nomissings]
y    = y[nomissings]
order   = sortperm(colr)#[length(colr):-1:1]
colr    = colr[order]
maxcolr = maximum(colr)
ax.grid()
#ax.plot([10,13], [10,13], "-",color="black", linewidth=0.5)
p = ax.scatter( log10.(x[order]), log10.(y[order]),
        zorder=2   , s= pointsize    , alpha=0.5,#colr ./ maxcolr, (pointsize/maxcolr) .* colr
        cmap="plasma_r", c=colr) # brg_r  
colbar  = fig.colorbar(p, ax=ax)#, ticks=[minimum(as03["BVAL_0"    ]), b_ell, b_disk, maximum(as03["BVAL_0"    ])])
#colbar.ax.set_yticklabels(["$(round(minimum(as03["BVAL_0"    ]), digits=1))", "ETG", "LTG", "$(round(maximum(as03["BVAL_0"    ]), digits=1))"])
colbar.ax.tick_params(axis="y")
colbar.ax.set_ylabel("Flip Angle [°]")
colbar.ax.yaxis.set_label_position("left")
#ax.legend(loc="lower right", frameon=true, borderpad=1, handlelength=1.8)
#ax.set_ylim(extrema(log10.(y[order])))
#ax.set_xlim(extrema(log10.(x[order])))
#ax.set_title("immediate read, all mergers, 03 Gyr time step")
#ax.set_aspect("equal")
ax.set_ylabel("log( j$(latexstring("\$_\\textrm{*,post}\$")) [ kpc km/s ] )")
ax.set_xlabel("log( j$(latexstring("\$_\\textrm{*,pre}\$")) [ kpc km/s ] )")
ax.grid(true)
#ax.autoscale()
fig.set_size_inches(11*scale, 9*scale)
fig.tight_layout()
fig.savefig(outfile, bbox_inches="tight", pad_inches=.1)
close("all")
PyPlot.matplotlib[:rc]("text", usetex=false)



outfile     = joinpath(outputdir,"J_post_v_J_pre_vs_phi.pdf")
scale = 0.5
PyPlot.matplotlib[:rc]("text", usetex=true)
PyPlot.matplotlib[:rc]("font", family="serif", serif=["Computer Modern"])
pointsize = 2
fig, ax = subplots()
slct    = csw03[:transition]#idxclude(csw03[:start_thr], csw03[:switches])
y       = pyplottable(as03["J_main"][:,slct])
x       = pyplottable(as03["J_main"][:,slct] .- as03["ΔJ_main"][:,slct])
colr    = as03["ϕ_flip"][slct]
nomissings = idxclude(1:length(colr), findcs(colr, eq=missing))
colr    = colr[nomissings]
x    = x[nomissings]
y    = y[nomissings]
order   = sortperm(colr)#[length(colr):-1:1]
colr    = colr[order]
maxcolr = maximum(colr)
ax.grid()
#ax.plot([10,13], [10,13], "-",color="black", linewidth=0.5)
p = ax.scatter( log10.(x[order]), log10.(y[order]),
        zorder=2   , s= pointsize    , alpha=0.5,#colr ./ maxcolr, (pointsize/maxcolr) .* colr
        cmap="plasma_r", c=colr) # brg_r  
colbar  = fig.colorbar(p, ax=ax)#, ticks=[minimum(as03["BVAL_0"    ]), b_ell, b_disk, maximum(as03["BVAL_0"    ])])
#colbar.ax.set_yticklabels(["$(round(minimum(as03["BVAL_0"    ]), digits=1))", "ETG", "LTG", "$(round(maximum(as03["BVAL_0"    ]), digits=1))"])
colbar.ax.tick_params(axis="y")
colbar.ax.set_ylabel("Flip Angle [°]")
colbar.ax.yaxis.set_label_position("left")
#ax.legend(loc="lower right", frameon=true, borderpad=1, handlelength=1.8)
#ax.set_ylim(extrema(log10.(y[order])))
#ax.set_xlim(extrema(log10.(x[order])))
#ax.set_title("immediate read, all mergers, 03 Gyr time step")
#ax.set_aspect("equal")
ax.set_ylabel("log( J$(latexstring("\$_\\textrm{*,post}\$")) [ M$(latexstring("\$_\\odot\$")) kpc km/s ] )")
ax.set_xlabel("log( J$(latexstring("\$_\\textrm{*,pre}\$")) [ M$(latexstring("\$_\\odot\$")) kpc km/s ] )")
ax.grid(true)
#ax.autoscale()
fig.set_size_inches(11*scale, 9*scale)
fig.tight_layout()
fig.savefig(outfile, bbox_inches="tight", pad_inches=.1)
close("all")
PyPlot.matplotlib[:rc]("text", usetex=false)


outfile     = joinpath(outputdir,"bpost_v_bpre_vs_phi.pdf")
scale = 0.5
PyPlot.matplotlib[:rc]("text", usetex=true)
PyPlot.matplotlib[:rc]("font", family="serif", serif=["Computer Modern"])
pointsize = 2
fig, ax = subplots()
slct    = csw03[:transition]#idxclude(csw03[:start_thr], csw03[:switches])
y       = pyplottable(as03["BVAL_0"][slct])
x       = pyplottable(as03["BVAL_0"][slct] .- as03["ΔBVAL_0"][slct])
colr    = as03["ϕ_flip"][slct]
nomissings = idxclude(1:length(colr), findcs(colr, eq=missing))
colr    = colr[nomissings]
x    = x[nomissings]
y    = y[nomissings]
order   = sortperm(colr)#[length(colr):-1:1]
colr    = colr[order]
maxcolr = maximum(colr)
ax.grid()
#ax.plot([10,13], [10,13], "-",color="black", linewidth=0.5)
p = ax.scatter( x[order], y[order],
        zorder=2   , s= pointsize    , alpha=0.5,#colr ./ maxcolr, (pointsize/maxcolr) .* colr
        cmap="plasma_r", c=colr) # brg_r  
colbar  = fig.colorbar(p, ax=ax)#, ticks=[minimum(as03["BVAL_0"    ]), b_ell, b_disk, maximum(as03["BVAL_0"    ])])
#colbar.ax.set_yticklabels(["$(round(minimum(as03["BVAL_0"    ]), digits=1))", "ETG", "LTG", "$(round(maximum(as03["BVAL_0"    ]), digits=1))"])
colbar.ax.tick_params(axis="y")
colbar.ax.set_ylabel("Flip Angle [°]")
colbar.ax.yaxis.set_label_position("left")
#ax.legend(loc="lower right", frameon=true, borderpad=1, handlelength=1.8)
#ax.set_ylim(extrema(log10.(y[order])))
#ax.set_xlim(extrema(log10.(x[order])))
#ax.set_title("immediate read, all mergers, 03 Gyr time step")
#ax.set_aspect("equal")
ax.set_ylabel("b-value$(latexstring("\$_\\textrm{post}\$"))")
ax.set_xlabel("b-value$(latexstring("\$_\\textrm{pre}\$"))")
ax.grid(true)
#ax.autoscale()
fig.set_size_inches(11*scale, 9*scale)
fig.tight_layout()
fig.savefig(outfile, bbox_inches="tight", pad_inches=.1)
close("all")
PyPlot.matplotlib[:rc]("text", usetex=false)





###########################
# dj / j vs flipangle

outnm   = joinpath(outputdir,"djrel_v_flipangle_vs_bval_03Gyr")
outfile = outnm*".pdf"
scale = 0.5
PyPlot.matplotlib[:rc]("text", usetex=true)
PyPlot.matplotlib[:rc]("font", family="serif", serif=["Computer Modern"])
pointsize = 2
fig, ax = subplots()
slct    = csw03[:transition]#idxclude(csw03[:start_thr], csw03[:switches])
y       = pyplottable(as03["Δj_main"][:,slct]) ./ pyplottable(as03["j_main"][:,slct] .- as03["Δj_main"][:,slct])
x       = pyplottable(as03["ϕ_flip"][slct])
colr    = as03["BVAL_0"][slct] .- as03["ΔBVAL_0"][slct]
nomissings = idxclude(1:length(colr), findcs(colr, eq=missing))
colr    = colr[nomissings]
x    = x[nomissings]
y    = y[nomissings]
order   = sortperm(colr)#[length(colr):-1:1]
colr    = colr[order]
maxcolr = maximum(colr)
ax.grid()
p = ax.scatter( x[order], log10.(norm.(y[order])),
        zorder=2   , s= pointsize    , alpha=0.5,
        cmap=ColorMap(get_bvalue_cmap().colors), c=colr, vmin=-6.1,vmax=-3.9, label="Sample, $(latexstring("\$\\Delta\$"))t = 0.3 Gyr") 
ax.plot(LinRange(0, 90, 1000), log10.(sind.(LinRange(0, 90, 1000))), "-",color="black", linewidth=2, label="Limit")
colbar  = fig.colorbar(p, ax=ax)
colbar.ax.tick_params(axis="y")
colbar.ax.set_ylabel("b-value at Flip Start")
colbar.ax.yaxis.set_label_position("left")
#ax.set_ylim([ymin,ymax])
ax.set_xlim([xmin,xmax])
ax.legend(loc="lower right", frameon=true, borderpad=1, handlelength=1.8)
ax.set_ylabel("log( $(latexstring("\$\\Delta\$"))j$(latexstring("\$_*\$")) / j$(latexstring("\$_\\textrm{*,pre}\$")) )")
ax.set_xlabel("Instant Flip Angle [°]")
ax.grid(true)
fig.set_size_inches(11*scale, 9*scale)
fig.tight_layout()
fig.savefig(outfile, bbox_inches="tight", pad_inches=.1)
close("all")
PyPlot.matplotlib[:rc]("text", usetex=false)


plot_hexbins( x, log10.(norm.(y)), 
    outfile=outnm*"_hx.pdf",
    scale=15, gridres=(50,23), xscale="linear", yscale="linear", 
    xlabel="Instant Flip Angle [°]", 
    ylabel="log( $(latexstring("\$_\\Delta\$"))j$(latexstring("\$_*\$")) / j$(latexstring("\$_\\textrm{*,pre}\$")) )", 
    grid=true,
    xmin=0, xmax=180,# ymin=ymin, ymax=ymax,
    lognorm=true,
    binned_median = "x", plot_bval=" ", plot_line=hcat(LinRange(0, 180, 1000), log10.(min_dj_v_flipangled.(LinRange(0, 180, 1000)))),
    calc_norm=true,    # @pyplottable
    )
#



##########################################################################
# dt hist
outnm   = joinpath(outputdir,"hist_dt")
outfile = outnm*".pdf"
scale = 0.5
lw=2
PyPlot.matplotlib[:rc]("text", usetex=true)
PyPlot.matplotlib[:rc]("font", family="serif", serif=["Computer Modern"])
fig, ax = subplots()
inputslct03 = csw03[:transition]#idxclude(csw03[:transition], findcs(norm.(as03["Δlookbacktime"]), eq=0.0))# 
inputslct09 = csw09[:transition]#idxclude(csw09[:transition], findcs(norm.(as09["Δlookbacktime"]), eq=0.0))# 
x03 = norm.(pyplottable(as03["Δlookbacktime"][inputslct03]))
x09 = norm.(pyplottable(as09["Δlookbacktime"][inputslct09]))
ax.hist( x03   , bins=10, label="$(latexstring("\$\\Delta\$"))t = 0.3 Gyr" , rwidth=0.9    , alpha=1, color="seagreen" )
ax.hist( x09   , bins=10, label="$(latexstring("\$\\Delta\$"))t = 1 Gyr"   , rwidth=0.9    , alpha=1, color="darkblue" )
ax.axvline(median(x03), color="black")
ax.axvline(median(x09), color="black")
ax.set_xlabel("Time Step Size [Gyr]")
#ax.set_xscale("log")
#ax.set_xlim([10, nothing])
ax.legend(loc="upper center", frameon=true, borderpad=1, handlelength=1.8)
#ax.grid()
scale=0.5
fig.set_size_inches(11scale, 9scale)
fig.tight_layout()
fig.savefig(outfile, bbox_inches="tight", pad_inches=.1)
close("all")
PyPlot.matplotlib[:rc]("text", usetex=false)

#for i in sdt
#        println("$(as09["snapNR"][ifelse(i-1 > 0, i-1, 1)]) - $(as09["ID_ISUB"][ifelse(i-1 > 0, i-1, 1)]) - $(as09["lookbacktime"][ifelse(i-1 > 0, i-1, 1)])   ---   $(as09["snapNR"][i]) - $(as09["ID_ISUB"][i]) - $(as09["lookbacktime"][i])")
#end












##########################################
# 10 % rvir spin vs full rvir spin vs alignment

outnm   = joinpath(outputdir,"J10_v_J100_vs_align")
outfile = outnm*".pdf"
scale = 0.5
PyPlot.matplotlib[:rc]("text", usetex=true)
PyPlot.matplotlib[:rc]("font", family="serif", serif=["Computer Modern"])
pointsize = 2
fig, ax = subplots()
slct    = 1:length(as03["M"])
colr    = Vector{Float64}(undef, 0)
for i in slct
        colr    = vcat( colr, (180/π) .* angle(convert(Array{Float64,1},as03["J_vir"][:,i]), convert(Array{Float64,1},as03["J_main"][:,i]) ) )
end
y       = pyplottable(as03["J_main"][:,slct])
x       = pyplottable(as03["J_vir"][:,slct])
nomissings = idxclude(1:length(colr), findcs(colr, eq=missing))
nomissings = idxclude(nomissings, findcs(x, eq=NaN))
nomissings = idxclude(nomissings, findcs(y, eq=NaN))
colr    = colr[nomissings]
x       = x[nomissings]
y       = y[nomissings]
order   = sortperm(colr)#[length(colr):-1:1]
colr    = colr[order]
x       = x[order]
y       = y[order]
maxcolr = maximum(colr)
ax.grid()
p = ax.scatter( log10.(x), log10.(y),
        zorder=2   , s= pointsize    , alpha=0.5,
        cmap="plasma_r", c=colr, vmin=0,vmax=180) 
colbar  = fig.colorbar(p, ax=ax)
colbar.ax.tick_params(axis="y")
colbar.ax.set_ylabel("Alignment Angle [°]")
colbar.ax.yaxis.set_label_position("left")
#ax.set_ylim([ymin,ymax])
#ax.set_xlim([xmin,xmax])
#ax.legend(loc="lower right", frameon=true, borderpad=1, handlelength=1.8)
ax.set_ylabel("log( J$(latexstring("\$_\\textrm{10}\$")) )")
ax.set_xlabel("log( J$(latexstring("\$_\\textrm{100}\$")) )")
ax.grid(true)
fig.set_size_inches(11*scale, 9*scale)
fig.tight_layout()
fig.savefig(outfile, bbox_inches="tight", pad_inches=.1)
close("all")
PyPlot.matplotlib[:rc]("text", usetex=false)

outnm   = joinpath(outputdir,"J10_v_J100_vs_align")
outfile = outnm*"_hxcount.pdf"
y       = log10.(pyplottable(as03["J_main"]))
x       = log10.(pyplottable(as03["J_vir"]))
miny= minimum(y)
maxy= maximum(y)
minx= minimum(x)
maxx= maximum(x)
slcts       = Dict{Int64, Vector{Int64}}()
slcts[1]    = 1:length(x)
slcts[3]    = csw03[:result_disk]
slcts[2]    = csw03[:result_int]
slcts[4]    = csw03[:result_ell]
labels      = Dict{Int64, String}()
labels[1]   = "Full Set"#Full Set, N = $(length(slcts[1]))"
labels[3]   = "Disks"#$(round(minimum(as["redshift"][slcts[3]]), digits=1))  $(latexstring("\$>\$")) z  $(latexstring("\$>\$")) $(round(maximum(as["redshift"][slcts[4]]), digits=1)), N = $(length(slcts[2]))"
labels[2]   = "Intermediates"#$(round(maximum(as["redshift"][slcts[3]]), digits=1)) $(latexstring("\$\\geq\$")) z $(latexstring("\$\\geq\$")) $(round(minimum(as["redshift"][slcts[3]]), digits=1)), N = $(length(slcts[3]))"
labels[4]   = "Ellipticals"#$(round(maximum(as["redshift"][slcts[4]]), digits=1)) $(latexstring("\$\\geq\$")) z, N = $(length(slcts[4]))"
labpos      = Array{String}(undef, 4)
labpos     .= "upper left"
plot_hexbins_nxm( x, y, 
            slcts, labels, 2, 2,
            #scale=0.0,
            cbarshift=2.5,
            label_pos=labpos, 
            outfile=outfile, 
            lognorm=true, gridres=(20,7), 
            ymin=miny, ymax=maxy,
            xmin=minx, xmax=maxx, 
            grid=true,
            ylabel="log( J$(latexstring("\$_{10}\$")) [ M$(latexstring("\$_\\odot\$")) kpc km/s ] )",
            xlabel="log( J$(latexstring("\$_{100}\$")) [ M$(latexstring("\$_\\odot\$")) kpc km/s ] )", fontsize=30,
            binned_median="y",
            plot_line=hcat(LinRange(11, 16, 10), LinRange(11, 16, 10))
            )
#

outnm   = joinpath(outputdir,"J10_v_J100_vs_align")
outfile = outnm*"_hxmean.pdf"
plot_hexbins_nxm( x, y, 
            slcts, labels, 2, 2,
            #scale=0.0,
            cbarshift=2.5,
            label_pos=labpos, 
            outfile=outfile, 
            lognorm=false, gridres=(20,7), 
            ymin=miny, ymax=maxy,
            xmin=minx, xmax=maxx, 
            grid=true,
            ylabel="log( J$(latexstring("\$_{10}\$")) [ M$(latexstring("\$_\\odot\$")) kpc km/s ] )",
            xlabel="log( J$(latexstring("\$_{100}\$")) [ M$(latexstring("\$_\\odot\$")) kpc km/s ] )", fontsize=30,
            binned_median="y",
            cfunc=mean, colormod=colr, cmap="plasma_r", colorlimits=(0,180),
            plot_line=hcat(LinRange(11, 16, 10), LinRange(11, 16, 10))
            )
#


#outputdir = "/home/moon/sfortune/spinevo/plots/thesis"
outnm   = joinpath(outputdir,"J10_v_J100_vs_bval")
outfile = outnm*"_hxmean.pdf"
colr    = pyplottable(as03["BVAL_0"])
y       = log10.(pyplottable(as03["J_main"]))
x       = log10.(pyplottable(as03["J_vir"]))
miny= minimum(y)
maxy= maximum(y)
minx= minimum(x)
maxx= maximum(x)
slcts       = Dict{Int64, Vector{Int64}}()
slcts[1]    = 1:length(x)
slcts[3]    = csw03[:result_disk]
slcts[2]    = csw03[:result_int]
slcts[4]    = csw03[:result_ell]
labels      = Dict{Int64, String}()
labels[1]   = "Full Set"#Full Set, N = $(length(slcts[1]))"
labels[3]   = "Disks"#$(round(minimum(as["redshift"][slcts[3]]), digits=1))  $(latexstring("\$>\$")) z  $(latexstring("\$>\$")) $(round(maximum(as["redshift"][slcts[4]]), digits=1)), N = $(length(slcts[2]))"
labels[2]   = "Intermediates"#$(round(maximum(as["redshift"][slcts[3]]), digits=1)) $(latexstring("\$\\geq\$")) z $(latexstring("\$\\geq\$")) $(round(minimum(as["redshift"][slcts[3]]), digits=1)), N = $(length(slcts[3]))"
labels[4]   = "Ellipticals"#$(round(maximum(as["redshift"][slcts[4]]), digits=1)) $(latexstring("\$\\geq\$")) z, N = $(length(slcts[4]))"
labpos      = Array{String}(undef, 4)
labpos     .= "upper left"

plot_hexbins_nxm( x, y, 
            slcts, labels, 2, 2,
            #scale=0.0,
            cbarshift=2.5,
            label_pos=labpos, 
            outfile=outfile, 
            lognorm=false, gridres=(20,7), 
            ymin=miny, ymax=maxy,
            xmin=minx, xmax=maxx, 
            grid=true,
            ylabel="log( J$(latexstring("\$_{10}\$")) [ M$(latexstring("\$_\\odot\$")) kpc km/s ] )",
            xlabel="log( J$(latexstring("\$_{100}\$")) [ M$(latexstring("\$_\\odot\$")) kpc km/s ] )", fontsize=30,
            binned_median="y",
            cfunc=mean, colormod=colr, cmap=ColorMap(get_bvalue_cmap().colors), colorlimits=(-6.1,-3.9),
            plot_line=hcat(LinRange(11, 16, 10), LinRange(11, 16, 10))
            )
#







###########################
# j_STARS v M_STARS vs bvalue

outnm   = joinpath(outputdir,"jSTAR_v_MSTAR_vs_bval")
outfile = outnm*"_hxcount.pdf"
colr    = pyplottable(as03["BVAL_0"])
y       = pyplottable(log10.(pyplottable(as03["j_main"])))
x       = pyplottable(log10.(pyplottable(as03["M"])))
println("NaN values = $(length(findcs(y, eq=NaN))+length(findcs(x, eq=NaN))length(findcs(colr, eq=NaN)))")
miny= minimum(y)
maxy= maximum(y)
minx= minimum(x)
maxx= maximum(x)
slcts       = Dict{Int64, Vector{Int64}}()
slcts[1]    = 1:length(x)
slcts[3]    = csw03[:z_0]
slcts[2]    = findcs(as03["redshift"], gt=minimum(as03["redshift"]), lt=1.5)
slcts[4]    = csw03[:z_2]
labels      = Dict{Int64, String}()
labels[1]   = "Full Set, N = $(length(slcts[1]))"#Full Set, N = $(length(slcts[1]))"
labels[3]   = "z = 0, N = $(length(slcts[3]))"#$(round(minimum(as["redshift"][slcts[3]]), digits=1))  $(latexstring("\$>\$")) z  $(latexstring("\$>\$")) $(round(maximum(as["redshift"][slcts[4]]), digits=1)), N = $(length(slcts[2]))"
labels[2]   = "0 $(latexstring("<")) z $(latexstring("<")) 1.5, N = $(length(slcts[2]))"#Full Set, N = $(length(slcts[1]))"
labels[4]   = "z = 2, N = $(length(slcts[4]))"#$(round(maximum(as["redshift"][slcts[3]]), digits=1)) $(latexstring("\$\\geq\$")) z $(latexstring("\$\\geq\$")) $(round(minimum(as["redshift"][slcts[3]]), digits=1)), N = $(length(slcts[3]))"
labpos      = Array{String}(undef, 4)
labpos     .= "lower right"
plot_hexbins_nxm( x, y, 
            slcts, labels, 2, 2,
            #scale=0.0,
            cbarshift=2.5,
            label_pos=labpos, 
            outfile=outfile, 
            lognorm=true, gridres=(20,7), 
            ymin=miny, ymax=maxy,
            xmin=minx, xmax=maxx, 
            grid=true,
            ylabel="log( j$(latexstring("\$_*\$")) [  kpc km/s ] )",
            xlabel="log( M$(latexstring("\$_*\$")) [ M$(latexstring("\$_\\odot\$")) ] )", fontsize=30,
            binned_median="y",
            plot_line=hcat(LinRange(minx, maxx, 10), b_disk .+ (2//3 .* LinRange(minx, maxx, 10)))
            )
#
outfile = outnm*"_hxmean.pdf"
plot_hexbins_nxm( x, y, 
            slcts, labels, 2, 2,
            #scale=0.0,
            cbarshift=2.5,
            label_pos=labpos, 
            outfile=outfile, 
            lognorm=false, gridres=(20,7), 
            ymin=miny, ymax=maxy,
            xmin=minx, xmax=maxx, 
            grid=true,
            ylabel="log( j$(latexstring("\$_*\$")) [  kpc km/s ] )",
            xlabel="log( M$(latexstring("\$_*\$")) [ M$(latexstring("\$_\\odot\$")) ] )", fontsize=30,
            binned_median="y",
            cfunc=mean, colormod=colr, cmap=ColorMap(get_bvalue_cmap().colors), colorlimits=(-6.1,-3.9),
            plot_line=hcat(LinRange(minx, maxx, 10), b_disk .+ (2//3 .* LinRange(minx, maxx, 10)))
            )
#







##########################################################################
# SFR hist
outfile     = joinpath(outputdir,"hist_SFR.pdf")
scale = 0.5
lw=2
PyPlot.matplotlib[:rc]("text", usetex=true)
PyPlot.matplotlib[:rc]("font", family="serif", serif=["Computer Modern"])
fig, ax = subplots()
ax.hist( pyplottable(log10.(pyplottable(1e9 .* as03["SFR"])))    , bins=20    , rwidth=0.9    , alpha=1, color="mediumblue" )
ax.axvline(median_logSFR, color="black")
ax.axvline(pctl68up_logSFR, color="dimgrey")
ax.axvline(pctl68lo_logSFR, color="dimgrey")
ax.set_xlabel("log( SFR [ M$(latexstring("\$_\\odot\$")) / Gyr ] )")
ax.grid()
fig.set_size_inches(11scale, 9scale)
fig.tight_layout()
fig.savefig(outfile, bbox_inches="tight", pad_inches=.1)
close("all")
PyPlot.matplotlib[:rc]("text", usetex=false)


outfile     = joinpath(outputdir,"hist_sSFR.pdf")
scale = 0.5
lw=2
PyPlot.matplotlib[:rc]("text", usetex=true)
PyPlot.matplotlib[:rc]("font", family="serif", serif=["Computer Modern"])
fig, ax = subplots()
ax.hist( pyplottable(log10.(pyplottable(1e9 .* as03["SFR"] ./ as03["M"])))    , bins=20    , rwidth=0.9    , alpha=1, color="dodgerblue" )
ax.axvline(median_logsSFR, color="black")
ax.axvline(pctl68up_logsSFR, color="dimgrey")
ax.axvline(pctl68lo_logsSFR, color="dimgrey")
ax.set_xlabel("log( sSFR [ 1 / Gyr ] )")
ax.grid()
fig.set_size_inches(11scale, 9scale)
fig.tight_layout()
fig.savefig(outfile, bbox_inches="tight", pad_inches=.1)
close("all")
PyPlot.matplotlib[:rc]("text", usetex=false)



for i in csw03[:transition]#inputslct
        if -1.01*ifelse(ismissing(symlog10(as03["ΔM"][i])), 0, symlog10(as03["ΔM"][i])) > log10(abs(as03["M"][i] - as03["ΔM"][i]))
                println(i, "   ", as03["ΔM"][i], "   ", as03["M"][i] - as03["ΔM"][i], "   ", as03["M"][i], "   ", as03["ID_ISUB"][i], "   ", as03["snapNR"][i])
        end
end

####################################################################################
# del M stars vs M stars vs flip angle
PyPlot.matplotlib[:rc]("text", usetex=true)
PyPlot.matplotlib[:rc]("font", family="serif", serif=["Computer Modern"])
outnm   = joinpath(outputdir,"delM_v_Mpre_vs_flipangle_l0_lores")
outfile = outnm*".pdf"
pointsize = 10
scale = 0.5
fig, ax = subplots()
inputslct = csw03[:transition]#idxclude(csw03[:start_thr], csw03[:switches])
colr    = pyplottable(as03["ϕ_flip"][inputslct])
y       = replace(symlog10.(as03["ΔM"][inputslct]), missing => 0)
x       = log10.(abs.(as03["M"][inputslct] .- as03["ΔM"][inputslct]))
order   = sortperm(colr)
colr    = colr[order]
y       = y[order]
x       = x[order]
slct    = findcs(x, gt=0., comparewith=findall(dummy->!isnan.(dummy), colr))
slct    = findcs(y, lt=0., comparewith=slct)
maxcolr = maximum(colr[slct])
#slct    = idxmatch(slct2, n_switches)
p = ax.scatter( x[slct], y[slct],
        zorder=2   , s= pointsize    , alpha=0.7,#colr[slct] ./ maxcolr, (pointsize/maxcolr) .* colr[slct]
        cmap="plasma_r", c=colr[slct], rasterized=true) # brg_r  
ax.plot([8,12], [-8,-12], "-",color="black", linewidth=0.5, rasterized=true)
colbar  = fig.colorbar(p, ax=ax)#, ticks=[minimum(as03["BVAL_0"    ]), b_ell, b_disk, maximum(as03["BVAL_0"    ])])
#colbar.ax.set_yticklabels(["$(round(minimum(as03["BVAL_0"    ]), digits=1))", "ETG", "LTG", "$(round(maximum(as03["BVAL_0"    ]), digits=1))"])
colbar.ax.tick_params(axis="y")
colbar.ax.set_ylabel("Flip Angle [°]")
colbar.ax.yaxis.set_label_position("left")
#ax.legend(loc="lower right", frameon=true, borderpad=1, handlelength=1.8)
ax.set_xlim([10,13])
ax.set_ylim([-12,-4.5])
#ax.set_title("N = $(length(slct))")
#ax.set_aspect("equal")
ax.set_ylabel("symlog( $(latexstring("\$\\Delta\$"))$(latexstring("\$M_*\$")) [$(latexstring("\$M_\\odot\$"))] )")
ax.set_xlabel("log( $(latexstring("\$M_{*,\\textrm{pre}}\$")) [$(latexstring("\$M_\\odot\$"))] )")
ax.grid()
#ax.autoscale()
fig.set_size_inches(10*scale, 9*scale)
fig.tight_layout()
fig.savefig(outfile, bbox_inches="tight", pad_inches=.1)
close("all")
PyPlot.matplotlib[:rc]("text", usetex=false)



PyPlot.matplotlib[:rc]("text", usetex=true)
PyPlot.matplotlib[:rc]("font", family="serif", serif=["Computer Modern"])
outnm   = joinpath(outputdir,"delM_v_Mpre_vs_flipangle_g0_lores")
outfile = outnm*".pdf"
pointsize = 10
scale = 0.5
fig, ax = subplots()
inputslct = csw03[:transition]#idxclude(csw03[:start_thr], csw03[:switches])
colr    = pyplottable(as03["ϕ_flip"][inputslct])
y       = replace( symlog10.(as03["ΔM"][inputslct]), missing => 0 )
x       = log10.(abs.(as03["M"][inputslct] .- as03["ΔM"][inputslct]))
order   = sortperm(colr)
colr    = colr[order]
y       = y[order]
x       = x[order]
slct    = findcs(x, gt=0., comparewith=findall(dummy->!isnan.(dummy), colr))
slct    = findcs(y, gt=0., comparewith=slct)
maxcolr = maximum(colr[slct])
#slct    = idxmatch(slct2, n_switches)
p = ax.scatter( x[slct], y[slct],
        zorder=2   , s= pointsize    , alpha=0.7,#colr[slct] ./ maxcolr, (pointsize/maxcolr) .* colr[slct]
        cmap="plasma_r", c=colr[slct], rasterized=true) # brg_r  
ax.plot([8,12], [8,12], "-",color="black", linewidth=0.5, rasterized=true)
colbar  = fig.colorbar(p, ax=ax)#, ticks=[minimum(as03["BVAL_0"    ]), b_ell, b_disk, maximum(as03["BVAL_0"    ])])
#colbar.ax.set_yticklabels(["$(round(minimum(as03["BVAL_0"    ]), digits=1))", "ETG", "LTG", "$(round(maximum(as03["BVAL_0"    ]), digits=1))"])
colbar.ax.tick_params(axis="y")
colbar.ax.set_ylabel("Flip Angle [°]")
colbar.ax.yaxis.set_label_position("left")
#ax.legend(loc="lower right", frameon=true, borderpad=1, handlelength=1.8)
ax.set_xlim([10,13])
ax.set_ylim([4.5,12])
#ax.set_title("N = $(length(slct))")
#ax.set_aspect("equal")
ax.set_ylabel("symlog( $(latexstring("\$|\\Delta\$"))$(latexstring("\$M_*|\$")) [$(latexstring("\$M_\\odot\$"))] )")
ax.set_xlabel("log( $(latexstring("\$M_{*,\\textrm{pre}}\$")) [$(latexstring("\$M_\\odot\$"))] )")
ax.grid()
#ax.autoscale()
fig.set_size_inches(10*scale, 9*scale)
fig.tight_layout()
fig.savefig(outfile, bbox_inches="tight", pad_inches=.1)
close("all")
PyPlot.matplotlib[:rc]("text", usetex=false)



# del M stars vs del M dm vs flip angle
PyPlot.matplotlib[:rc]("text", usetex=true)
PyPlot.matplotlib[:rc]("font", family="serif", serif=["Computer Modern"])
outnm   = joinpath(outputdir,"delM_v_delMdm_vs_flipangle_lores")
outfile = outnm*".pdf"
pointsize = 10
scale = 0.5
fig, ax = subplots()
inputslct = csw03[:transition]#idxclude(csw03[:start_thr], csw03[:switches])
colr    = pyplottable(as03["ϕ_flip"][inputslct])
y       = log10.(abs.(as03["ΔM"][inputslct]))
x       = log10.(abs.(ad03["ΔM"][inputslct]))
order   = sortperm(colr)
colr    = colr[order]
y       = y[order]
x       = x[order]
slct    = findcs(x, gt=0., comparewith=findall(dummy->!isnan.(dummy), colr))
slct    = findcs(y, gt=0., comparewith=slct)
maxcolr = maximum(colr[slct])
#slct    = idxmatch(slct2, n_switches)
p = ax.scatter( x[slct], y[slct],
        zorder=2   , s= pointsize    , alpha=0.7,#colr[slct] ./ maxcolr, (pointsize/maxcolr) .* colr[slct]
        cmap="plasma_r", c=colr[slct], rasterized=true) # brg_r  
ax.plot([7.5,14], [7.5,14], "-",color="black", linewidth=0.5, rasterized=true)
colbar  = fig.colorbar(p, ax=ax)#, ticks=[minimum(as03["BVAL_0"    ]), b_ell, b_disk, maximum(as03["BVAL_0"    ])])
#colbar.ax.set_yticklabels(["$(round(minimum(as03["BVAL_0"    ]), digits=1))", "ETG", "LTG", "$(round(maximum(as03["BVAL_0"    ]), digits=1))"])
colbar.ax.tick_params(axis="y")
colbar.ax.set_ylabel("Flip Angle [°]")
colbar.ax.yaxis.set_label_position("left")
#ax.legend(loc="lower right", frameon=true, borderpad=1, handlelength=1.8)
ax.set_xlim([7.5,14])
ax.set_ylim([4.5,12])
#ax.set_title("N = $(length(slct))")
#ax.set_aspect("equal")
ax.set_ylabel("log( $(latexstring("\$|\\Delta\$"))$(latexstring("\$M_*|\$")) [$(latexstring("\$M_\\odot\$"))] )")
ax.set_xlabel("log( $(latexstring("\$|\\Delta\$"))$(latexstring("\$M_\\textrm{DM}|\$")) [$(latexstring("\$M_\\odot\$"))] )")
ax.grid()
#ax.autoscale()
fig.set_size_inches(10*scale, 9*scale)
fig.tight_layout()
fig.savefig(outfile, bbox_inches="tight", pad_inches=.1)
close("all")
PyPlot.matplotlib[:rc]("text", usetex=false)



# del M stars vs del M dm vs flip angle
PyPlot.matplotlib[:rc]("text", usetex=true)
PyPlot.matplotlib[:rc]("font", family="serif", serif=["Computer Modern"])
outnm   = joinpath(outputdir,"delM_v_delMdm_vs_flipangle_REL_lores")
outfile = outnm*".pdf"
pointsize = 10
scale = 0.5
fig, ax = subplots()
inputslct = csw03[:transition]#idxclude(csw03[:start_thr], csw03[:switches])
colr    = pyplottable(as03["ϕ_flip"][inputslct])
y       = log10.(abs.(as03["ΔM"][inputslct] ./ (as03["M"][inputslct] .- as03["ΔM"][inputslct])))
x       = log10.(abs.(ad03["ΔM"][inputslct] ./ (ad03["M"][inputslct] .- ad03["ΔM"][inputslct])))
order   = sortperm(colr)
colr    = colr[order]
y       = y[order]
x       = x[order]
maxcolr = maximum(colr)
#slct    = idxmatch(slct2, n_switches)
p = ax.scatter( x, y,
        zorder=2   , s= pointsize    , alpha=0.7,#colr[slct] ./ maxcolr, (pointsize/maxcolr) .* colr[slct]
        cmap="plasma_r", c=colr, rasterized=true) # brg_r  
ax.plot([-5,1], [-5,1], "-",color="black", linewidth=0.5, rasterized=true)#, rasterized=true)
colbar  = fig.colorbar(p, ax=ax)#, ticks=[minimum(as03["BVAL_0"    ]), b_ell, b_disk, maximum(as03["BVAL_0"    ])])
#colbar.ax.set_yticklabels(["$(round(minimum(as03["BVAL_0"    ]), digits=1))", "ETG", "LTG", "$(round(maximum(as03["BVAL_0"    ]), digits=1))"])
colbar.ax.tick_params(axis="y")
colbar.ax.set_ylabel("Flip Angle [°]")
colbar.ax.yaxis.set_label_position("left")
#ax.legend(loc="lower right", frameon=true, borderpad=1, handlelength=1.8)
#ax.set_xlim([7.5,14])
#ax.set_ylim([4.5,12])
#ax.set_title("N = $(length(slct))")
#ax.set_aspect("equal")
ax.set_ylabel("log( $(latexstring("\$|\\Delta\$"))$(latexstring("\$M_*|\$")) / $(latexstring("\$M_\\textrm{*,pre}\$")) )")
ax.set_xlabel("log( $(latexstring("\$|\\Delta\$"))$(latexstring("\$M_\\textrm{DM}|\$")) / $(latexstring("\$M_\\textrm{DM,pre}\$")) )")
#ax.set_ylabel("log( $(latexstring("\$|\\Delta M_*|\$")) / $(latexstring("\$M_\\textrm{*,pre}\$")) )")
#ax.set_xlabel("log( $(latexstring("\$|\\Delta M_\\textrm{DM}|\$")) / $(latexstring("\$M_\\textrm{DM,pre}\$")) )")
ax.grid()
#ax.autoscale()
fig.set_size_inches(10*scale, 9*scale)
fig.tight_layout()
fig.savefig(outfile, bbox_inches="tight", pad_inches=.1)
close("all")
PyPlot.matplotlib[:rc]("text", usetex=false)


##############################################################################################################
# d j vs d M_stars
##############################################################################################################
# combined
outnm   = joinpath(outputdir,"dj_vs_dM_dt-Mtype-split")
outfile = outnm*".pdf"
#y   = vcat(log10.(pyplottable(as03["Δj_main"][:,csw03[:transition]])), log10.(pyplottable(as03["Δj_main"][:,csw03[:transition]])), log10.(pyplottable(as03["Δj_main"][:,csw03[:transition]])), log10.(pyplottable(as09["Δj_main"][:,csw09[:transition]])), log10.(pyplottable(as09["Δj_main"][:,csw09[:transition]])), log10.(pyplottable(as09["Δj_main"][:,csw09[:transition]])))
#x   = vcat(log10.(abs.(pyplottable(as03["ΔM"][csw03[:transition]]))), log10.(abs.(pyplottable(as03["Mpeak_MERGERS"][csw03[:transition]]))), log10.(abs.(pyplottable(as03["Mpeak_MM"][csw03[:transition]]))), log10.(abs.(pyplottable(as09["ΔM"][csw09[:transition]]))), log10.(abs.(pyplottable(as09["Mpeak_MERGERS"][csw09[:transition]]))), log10.(abs.(pyplottable(as09["Mpeak_MM"][csw09[:transition]]))))
y   = vcat(log10.(pyplottable(as03["Δj_main"][:,csw03[:transition]])), log10.(pyplottable(as09["Δj_main"][:,csw09[:transition]])), log10.(pyplottable(as03["Δj_main"][:,csw03[:transition]])), log10.(pyplottable(as09["Δj_main"][:,csw09[:transition]])), log10.(pyplottable(as03["Δj_main"][:,csw03[:transition]])), log10.(pyplottable(as09["Δj_main"][:,csw09[:transition]])))
x   = vcat(log10.(abs.(pyplottable(as03["ΔM"][csw03[:transition]]))), log10.(abs.(pyplottable(as09["ΔM"][csw09[:transition]]))), log10.(abs.(pyplottable(as03["Mpeak_MERGERS"][csw03[:transition]]))), log10.(abs.(pyplottable(as09["Mpeak_MERGERS"][csw09[:transition]]))), log10.(abs.(pyplottable(as03["Mpeak_MM"][csw03[:transition]]))), log10.(abs.(pyplottable(as09["Mpeak_MM"][csw09[:transition]]))))
miny= minimum(pyplottable(y)[findall(x->x .!== NaN, pyplottable(y))])
maxy= maximum(pyplottable(y)[findall(x->x .!== NaN, pyplottable(y))])
minx= minimum(pyplottable(x)[findall(x->x .!== NaN, pyplottable(x))])
maxx= maximum(pyplottable(x)[findall(x->x .!== NaN, pyplottable(x))])
slcts       = Dict{Int64, Vector{Int64}}()
#slcts[1]    = 1:length(csw03[:transition])
#slcts[3]    = length(csw03[:transition])+1:2*length(csw03[:transition])
#slcts[5]    = 2*length(csw03[:transition])+1:3*length(csw03[:transition])
#slcts[2]    = 3*length(csw03[:transition])+1:3*length(csw03[:transition])+length(csw09[:transition])
#slcts[4]    = 3*length(csw09[:transition])+length(csw09[:transition])+1:3*length(csw09[:transition])+2*length(csw09[:transition])
#slcts[6]    = 3*length(csw09[:transition])+2*length(csw09[:transition])+1:3*length(csw09[:transition])+3*length(csw09[:transition])
slcts[1]    = 1:length(csw03[:transition])
slcts[2]    = length(slcts[1])+1:length(slcts[1])+length(csw09[:transition])
slcts[3]    = length(slcts[1])+length(slcts[2])+1:length(slcts[1])+length(slcts[2])+length(csw03[:transition])
slcts[4]    = length(slcts[1])+length(slcts[2])+length(slcts[3])+1:length(slcts[1])+length(slcts[2])+length(slcts[3])+length(csw09[:transition])
slcts[5]    = length(slcts[1])+length(slcts[2])+length(slcts[3])+length(slcts[4])+1:length(slcts[1])+length(slcts[2])+length(slcts[3])+length(slcts[4])+length(csw03[:transition])
slcts[6]    = length(slcts[1])+length(slcts[2])+length(slcts[3])+length(slcts[4])+length(slcts[5])+1:length(y)
labels      = Dict{Int64, String}()
labels[1]   = "$(latexstring("\$\\Delta\$"))t $(latexstring("\$\\sim\$")) 0.3 Gyr, $(latexstring("\$\\Delta\$"))M"
labels[3]   = "$(latexstring("\$\\Delta\$"))t $(latexstring("\$\\sim\$")) 0.3 Gyr, M$(latexstring("\$_{mergers}\$"))"
labels[5]   = "$(latexstring("\$\\Delta\$"))t $(latexstring("\$\\sim\$")) 0.3 Gyr, M$(latexstring("\$_{MM}\$"))"
labels[2]   = "$(latexstring("\$\\Delta\$"))t $(latexstring("\$\\sim\$")) 1 Gyr, $(latexstring("\$\\Delta\$"))M"
labels[4]   = "$(latexstring("\$\\Delta\$"))t $(latexstring("\$\\sim\$")) 1 Gyr, M$(latexstring("\$_{mergers}\$"))"
labels[6]   = "$(latexstring("\$\\Delta\$"))t $(latexstring("\$\\sim\$")) 1 Gyr, M$(latexstring("\$_{MM}\$"))"
labpos      = Array{String}(undef, 6)
labpos     .= "lower left"
labpos[1]   = "upper left"
plot_hexbins_nxm( x, 
            y, 
            slcts, labels, 2, 3,
            #scale=0.0,
            cbarshift=2.5,
            outfile=outfile, 
            lognorm=true, gridres=(20,10), label_pos=labpos,
            ymin=miny, ymax=maxy,
            xmin=minx, xmax=maxx, 
            ylabel="log( $(latexstring("\$\\Delta\$"))j [kpc km/s] )",
            xlabel="log( [$(latexstring("\$M_\\odot\$"))] )", fontsize=30,
            binned_median="y", plot_bval=" ",
            )
#

for i in csw09[:transition]
        if ismissing(as09["Mpeak_MERGERS"][i]) || isnan(as09["Mpeak_MERGERS"][i]) || ismissing(as09["Mpeak_MM"][i]) || isnan(as09["Mpeak_MM"][i])
        elseif as09["Mpeak_MERGERS"][i] < as09["Mpeak_MM"][i]
                println(i, "   ", as09["Mpeak_MERGERS"][i] ,"   ", as09["Mpeak_MM"][i])
        end
end


########################################################################
outnm   = "sSFR_vs_Mstar"
y       = pyplottable(log10.(as03["SFR"] ./ as03["M"] ))
x       = pyplottable(log10.(as03["M"]))

plot_hexbins( x, y, 
    outfile=joinpath(outputdir,outnm*"_hx.pdf"),
    scale=15, gridres=(25,11), xscale="linear", yscale="linear", 
    xlabel="log( $(latexstring("\$M_{*}\$")) [$(latexstring("\$M_\\odot\$"))] )", 
    ylabel="log( sSFR [$(latexstring("\\frac{1}{yr}"))] )", 
    grid=true,
    #xmin=0, xmax=180,# ymin=ymin, ymax=ymax,
    lognorm=true,
    binned_median = "y", plot_bval=" ", #plot_line=hcat(LinRange(0, 180, 1000), log10.(min_dj_v_flipangled.(LinRange(0, 180, 1000)))),
    calc_norm=true,    # @pyplottable
    )
#




########################################################################
outnm   = "dJ_vs_bval_03Gyr"
y       = log10.(pyplottable( as03["ΔJ_main"][:,csw03[:transition]] ))
x       = pyplottable( as03["ΔBVAL_0"][csw03[:transition]] )
plot_hexbins( x, y, 
    outfile=joinpath(outputdir,outnm*"_hx.pdf"),
    scale=15, gridres=(25,12), xscale="linear", yscale="linear", 
    xlabel="$(latexstring("\$\\Delta \\textrm{b-value}\$"))", 
    ylabel="log( $(latexstring("\$|\\Delta J_*|\$")) [ $(latexstring("\$M_\\odot\$")) kpc km/s ] )", 
    grid=true,
    xmin=-1.5, xmax=1.5, ymin=10.5, ymax=16,
    lognorm=true,
    binned_median = "y", plot_bval=" ", #plot_line=hcat(LinRange(0, 180, 1000), log10.(min_dj_v_flipangled.(LinRange(0, 180, 1000)))),
    calc_norm=true,    # @pyplottable
    )
#

outnm   = "dJrelpre_vs_bval_03Gyr"
y       = log10.(pyplottable( as03["ΔJ_main"][:,csw03[:transition]]) ./ pyplottable(as03["J_main"][:,csw03[:transition]] .- as03["ΔJ_main"][:,csw03[:transition]]) )
x       = pyplottable( as03["ΔBVAL_0"][csw03[:transition]] )
plot_hexbins( x, y, 
    outfile=joinpath(outputdir,outnm*"_hx.pdf"),
    scale=15, gridres=(25,12), xscale="linear", yscale="linear", 
    xlabel="$(latexstring("\$\\Delta \\textrm{b-value}\$"))", 
    ylabel="log( $(latexstring("\$|\\Delta J_*|\$")) / $(latexstring("\$|J_{*,\\textrm{pre}}|\$")) )", 
    grid=true,
    xmin=-1.5, xmax=1.5, ymin=-2, ymax=1.5,
    lognorm=true,
    binned_median = "y", plot_bval=" ", #plot_line=hcat(LinRange(0, 180, 1000), log10.(min_dj_v_flipangled.(LinRange(0, 180, 1000)))),
    calc_norm=true,    # @pyplottable
    )
#


outnm   = "dJ_vs_bval_09Gyr"
y       = log10.(pyplottable( as09["ΔJ_main"][:,csw09[:transition]] ))
x       = pyplottable( as09["ΔBVAL_0"][csw09[:transition]] )
plot_hexbins( x, y, 
    outfile=joinpath(outputdir,outnm*"_hx.pdf"),
    scale=15, gridres=(25,12), xscale="linear", yscale="linear", 
    xlabel="$(latexstring("\$\\Delta \\textrm{b-value}\$"))", 
    ylabel="log( $(latexstring("\$|\\Delta J_*|\$")) [ $(latexstring("\$M_\\odot\$")) kpc km/s ] )", 
    grid=true,
    xmin=-1.5, xmax=1.5, ymin=10.5, ymax=16,
    lognorm=true,
    binned_median = "y", plot_bval=" ", #plot_line=hcat(LinRange(0, 180, 1000), log10.(min_dj_v_flipangled.(LinRange(0, 180, 1000)))),
    calc_norm=true,    # @pyplottable
    )
#

outnm   = "dJrelpre_vs_bval_09Gyr"
y       = log10.(pyplottable( as09["ΔJ_main"][:,csw09[:transition]]) ./ pyplottable(as09["J_main"][:,csw09[:transition]] .- as09["ΔJ_main"][:,csw09[:transition]]) )
x       = pyplottable( as09["ΔBVAL_0"][csw09[:transition]] )
plot_hexbins( x, y, 
    outfile=joinpath(outputdir,outnm*"_hx.pdf"),
    scale=15, gridres=(25,12), xscale="linear", yscale="linear", 
    xlabel="$(latexstring("\$\\Delta \\textrm{b-value}\$"))", 
    ylabel="log( $(latexstring("\$|\\Delta J_*|\$")) / $(latexstring("\$|J_{*,\\textrm{pre}}|\$")) )", 
    grid=true,
    xmin=-1.5, xmax=1.5, ymin=-2, ymax=1.5,
    lognorm=true,
    binned_median = "y", plot_bval=" ", #plot_line=hcat(LinRange(0, 180, 1000), log10.(min_dj_v_flipangled.(LinRange(0, 180, 1000)))),
    calc_norm=true,    # @pyplottable
    )
#





##############################################################################################################
# Flipangle distribution
##############################################################################################################
outnm   = "hist_flipangle"
rand_res       = 10000
scale=0.5
N_bins=50
fig, ax = subplots()
n_03, x_03     = PyPlot.hist( pyplottable(as03["ϕ_flip" ][csw03[:transition]]), bins=N_bins, density=true, log=false)
n_09, x_09     = PyPlot.hist( pyplottable(as09["ϕ_flip" ][csw09[:transition]]), bins=N_bins, density=true, log=false)
bcenters_03    = 0.5 .* (x_03[2:end]    .+ x_03[1:end-1])
bcenters_09    = 0.5 .* (x_09[2:end]    .+ x_09[1:end-1])
bcenters_rd    = LinRange(0,180,rand_res)
n_rd      = sind.(bcenters_rd)
n_rd    .*= 0.5 / sum(sind.(bcenters_rd)) / N_bins * rand_res
clf()
fig, ax = subplots()
PyPlot.matplotlib[:rc]("text", usetex=true)
PyPlot.matplotlib[:rc]("font", family="serif", serif=["Computer Modern"])
ax.plot(bcenters_rd , n_rd , "-", linewidth=2, label= "Random Flips"   , color="black")
ax.plot(bcenters_03 , n_03 , "-", linewidth=2, label= "Instant Flips"   , color="mediumblue")
ax.plot(bcenters_09 , n_09 , "-", linewidth=2, label= "1 Gyr Flips"   , color="darkred")
ax.set_xlabel("$(latexstring("\$\\phi_{flip}\$")) [°]")
ax.set_ylabel("P")
ax.set_yscale("log")
ax.set_xlim(0,180)
ax.set_ylim(nothing,nothing)
ax.legend(loc="lower left", frameon=true, borderpad=1, handlelength=1.8)
ax.grid(false)
fig.set_size_inches(10scale, 9scale)
fig.tight_layout()
fig.savefig(joinpath(outputdir,outnm*".pdf"), bbox_inches="tight", pad_inches=.1)
close("all")
PyPlot.matplotlib[:rc]("text", usetex=false)




##################################################
# bvalue pre vs delta bvalue
outnm   = "bpre_vs_delb_03Gyr"
y       = pyplottable( as03["BVAL_0"][csw03[:transition]] .- as03["ΔBVAL_0"][csw03[:transition]]) 
x       = pyplottable(as03["ΔBVAL_0"][csw03[:transition]] )
plot_hexbins( x, y, 
    outfile=joinpath(outputdir,outnm*"_hxcount.pdf"),
    scale=15, gridres=(25,12), xscale="linear", yscale="linear", 
    xlabel="$(latexstring("\$\\Delta \\textrm{b-value}\$")) (0.3 Gyr)", 
    ylabel="$(latexstring("\$b_\\textrm{pre}\$"))", 
    grid=true,
    xmin=-1., xmax=1., ymin=-6.1, ymax=-3.9, 
    lognorm=true,
    binned_median = "y", plot_bval="y", #plot_line=hcat(LinRange(0, 180, 1000), log10.(min_dj_v_flipangled.(LinRange(0, 180, 1000)))),
    calc_norm=true,    # @pyplottable
    )
#

outnm   = "bpre_vs_delb_09Gyr"
y       = pyplottable( as09["BVAL_0"][csw09[:transition]] .- as09["ΔBVAL_0"][csw09[:transition]]) 
x       = pyplottable( as09["ΔBVAL_0"][csw09[:transition]] )
plot_hexbins( x, y, 
    outfile=joinpath(outputdir,outnm*"_hxcount.pdf"),
    scale=15, gridres=(25,12), xscale="linear", yscale="linear", 
    xlabel="$(latexstring("\$\\Delta \\textrm{b-value}\$")) (1 Gyr)", 
    ylabel="$(latexstring("\$b_\\textrm{pre}\$"))", 
    grid=true,
    xmin=-1., xmax=1., ymin=-6.1, ymax=-3.9, 
    lognorm=true,
    binned_median = "y", plot_bval="y", #plot_line=hcat(LinRange(0, 180, 1000), log10.(min_dj_v_flipangled.(LinRange(0, 180, 1000)))),
    calc_norm=true,    # @pyplottable
    )
#


outnm   = "bpre_vs_delb_03Gyr"
y       = pyplottable( as03["BVAL_0"][csw03[:transition]] .- as03["ΔBVAL_0"][csw03[:transition]]) 
x       = pyplottable(as03["ΔBVAL_0"][csw03[:transition]] )
colr    = pyplottable(as03["ϕ_flip"][csw03[:transition]] )
plot_hexbins( x, y, 
    outfile=joinpath(outputdir,outnm*"_hxmean.pdf"),
    scale=15, gridres=(25,12), xscale="linear", yscale="linear", 
    xlabel="$(latexstring("\$\\Delta \\textrm{b-value}\$")) (0.3 Gyr)", 
    ylabel="$(latexstring("\$b_\\textrm{pre}\$"))", 
    grid=true,
    xmin=-1., xmax=1., ymin=-6.1, ymax=-3.9, 
    lognorm=false,
    cfunc=mean, colormod=colr, cmap="plasma_r", colorlimits=(0,180),
    binned_median = "y", plot_bval="y", #plot_line=hcat(LinRange(0, 180, 1000), log10.(min_dj_v_flipangled.(LinRange(0, 180, 1000)))),
    calc_norm=true,    # @pyplottable
    )
#

outnm   = "bpre_vs_delb_09Gyr"
y       = pyplottable( as09["BVAL_0"][csw09[:transition]] .- as09["ΔBVAL_0"][csw09[:transition]]) 
x       = pyplottable( as09["ΔBVAL_0"][csw09[:transition]] )
colr    = pyplottable(as09["ϕ_flip"][csw09[:transition]] )
plot_hexbins( x, y, 
    outfile=joinpath(outputdir,outnm*"_hxmean.pdf"),
    scale=15, gridres=(25,12), xscale="linear", yscale="linear", 
    xlabel="$(latexstring("\$\\Delta \\textrm{b-value}\$")) (1 Gyr)", 
    ylabel="$(latexstring("\$b_\\textrm{pre}\$"))", 
    grid=true,
    xmin=-1., xmax=1., ymin=-6.1, ymax=-3.9, 
    lognorm=false,
    cfunc=mean, colormod=colr, cmap="plasma_r", colorlimits=(0,180),
    binned_median = "y", plot_bval="y", #plot_line=hcat(LinRange(0, 180, 1000), log10.(min_dj_v_flipangled.(LinRange(0, 180, 1000)))),
    calc_norm=true,    # @pyplottable
    )
#




########################################################################
outnm   = "MM_v_mergers_Mpeak_03"
y       = pyplottable(log10.(as03["Mpeak_MM"][csw03[:transition]]))
x       = pyplottable(log10.(as03["Mpeak_MERGERS"][csw03[:transition]]))

plot_hexbins( x, y, 
    outfile=joinpath(outputdir,outnm*"_hxcount.pdf"),
    scale=15, gridres=(50,22), xscale="linear", yscale="linear", 
    xlabel="log( $(latexstring("\$M_\\textrm{*,mergers}\$")) [$(latexstring("\$M_\\odot\$"))] )", 
    ylabel="log( $(latexstring("\$M_\\textrm{*,MM}\$")) [$(latexstring("\$M_\\odot\$"))] )", 
    grid=true,
    #xmin=0, xmax=180,# ymin=ymin, ymax=ymax,
    lognorm=true,
    binned_median = "y", plot_bval=" ", plot_line=hcat(LinRange(8, 12, 1000), LinRange(8, 12, 1000)),
    calc_norm=true,    # @pyplottable
    )
#


outnm   = "MM_v_mergers_Mpeak_09"
y       = pyplottable(log10.(as09["Mpeak_MM"][csw09[:transition]]))
x       = pyplottable(log10.(as09["Mpeak_MERGERS"][csw09[:transition]]))

plot_hexbins( x, y, 
    outfile=joinpath(outputdir,outnm*"_hxcount.pdf"),
    scale=15, gridres=(50,22), xscale="linear", yscale="linear", 
    xlabel="log( $(latexstring("\$M_\\textrm{*,mergers}\$")) [$(latexstring("\$M_\\odot\$"))] )", 
    ylabel="log( $(latexstring("\$M_\\textrm{*,MM}\$")) [$(latexstring("\$M_\\odot\$"))] )", 
    grid=true,
    #xmin=0, xmax=180,# ymin=ymin, ymax=ymax,
    lognorm=true,
    binned_median = "y", plot_bval=" ", plot_line=hcat(LinRange(8, 12, 1000), LinRange(8, 12, 1000)),
    calc_norm=true,    # @pyplottable
    )
#

















############################################################################################
# Hist flipangles split by mergernumber
for (ias, ia1, ia2, ia3, ic) in zip(as, as_label1, as_label2, as_label3, csw)
        println(ia3)
        outfile     = joinpath(outputdir,"hist_flip_Nmergersplit_$(ia1).pdf")
        scale=0.5
        N_bins=5
        pltcm       = pyimport("matplotlib.cm")
        pltcolors   = pyimport("matplotlib.colors")
        colormap    = pltcm.get_cmap(name="gnuplot")
        PyPlot.matplotlib[:rc]("text", usetex=true)
        PyPlot.matplotlib[:rc]("font", family="serif", serif=["Computer Modern"])
        fig, ax = subplots()
        PyPlot.matplotlib[:rc]("text", usetex=true)
        PyPlot.matplotlib[:rc]("font", family="serif", serif=["Computer Modern"])
        n_0, x_0        = PyPlot.hist( cosd.(pyplottable(ias["ϕ_flip"][idxmatch(ic[:Nmerger_0] , ic[:transition])])), bins=N_bins, range=(0,1), density=true, log=false)
        n_1,x_1         = PyPlot.hist( cosd.(pyplottable(ias["ϕ_flip"][idxmatch(ic[:Nmerger_1_2], ic[:transition])])), bins=N_bins, range=(0,1), density=true, log=false)
        n_3, x_3        = PyPlot.hist( cosd.(pyplottable(ias["ϕ_flip"][idxmatch(ic[:Nmerger_3_10 ], ic[:transition])])), bins=N_bins, range=(0,1), density=true, log=false)
        n_11, x_11      = PyPlot.hist( cosd.(pyplottable(ias["ϕ_flip"][idxmatch(ic[:Nmerger_11_ ], ic[:transition])])), bins=N_bins, range=(0,1), density=true, log=false)
        bcenters_0     = 0.5 .* (x_0[2:end]  .+ x_0[1:end-1])
        bcenters_1     = 0.5 .* (x_1[2:end]  .+ x_1[1:end-1])
        bcenters_3     = 0.5 .* (x_3[2:end]  .+ x_3[1:end-1])
        bcenters_11    = 0.5 .* (x_11[2:end]  .+ x_11[1:end-1])
        clf()
        PyPlot.matplotlib[:rc]("text", usetex=true)
        PyPlot.matplotlib[:rc]("font", family="serif", serif=["Computer Modern"])
        fig, ax = subplots()
        PyPlot.matplotlib[:rc]("text", usetex=true)
        PyPlot.matplotlib[:rc]("font", family="serif", serif=["Computer Modern"])
        ax.grid(true)
        ax.plot(bcenters_0 , log10.(n_0 ), ".-", ms=10, label="0 Mergers, N = $(length(idxmatch(ic[:Nmerger_0] , ic[:transition])))", color=colormap(0/3)   )
        ax.plot(bcenters_1 , log10.(n_1 ), ".-", ms=10, label="1-2 Mergers, N = $(length(idxmatch(ic[:Nmerger_1_2] , ic[:transition])))", color=colormap(1/3)   )
        ax.plot(bcenters_3 , log10.(n_3 ), ".-", ms=10, label="3-10 Mergers, N = $(length(idxmatch(ic[:Nmerger_3_10] , ic[:transition])))", color=colormap(2/3)   )
        ax.plot(bcenters_11, log10.(n_11), ".-", ms=10, label="$(latexstring(">"))11 Mergers, N = $(length(idxmatch(ic[:Nmerger_11_] , ic[:transition])))", color=colormap(3/3)   )
        
        ax.set_xlabel("cos( $(latexstring("\$\\phi_\\textrm{flip, $(ia2)}\$")) )")# $(latexstring("\$\\phi _\\textrm{flip}\$"))
        ax.set_ylabel("log( P )")
        ax.set_yscale("linear")
        
        ax.set_xlim(0,1)
        ax.set_ylim(-2.8,1)
        ax.legend(loc="upper left", frameon=true, borderpad=1, handlelength=1.8)
        fig.set_size_inches(10scale, 9scale)
        fig.tight_layout()
        fig.savefig(outfile, bbox_inches="tight", pad_inches=.1)
        PyPlot.matplotlib[:rc]("text", usetex=false)
        close("all")
end
#







##########################################################################################
# flipangle vs endmass bval-split
outnm   = joinpath(outputdir,"flipangle_vs_masspost_bvalstartsplit_03")
outfile = outnm*".pdf"
x   = log10.(as03["M"])
y   = as03["ϕ_flip"]
miny= 0#minimum(pyplottable(y)[findall(x->x .!== NaN, pyplottable(y))])
maxy= 180#maximum(pyplottable(y)[findall(x->x .!== NaN, pyplottable(y))])
minx= minimum(pyplottable(x)[findall(x->x .!== NaN, pyplottable(x))])
maxx= maximum(pyplottable(x)[findall(x->x .!== NaN, pyplottable(x))])
slcts       = Dict{Int64, Vector{Int64}}()
slcts[1]    = csw03[:transition]
slcts[2]    = idxmatch(csw03[:transition], csw03[:start_int ])
slcts[3]    = idxmatch(csw03[:transition], csw03[:start_disk])
slcts[4]    = idxmatch(csw03[:transition], csw03[:start_ell ])
labels      = Dict{Int64, String}()
labels[1]   = "Full Set, N = $(length(slcts[1]))"
labels[2]   = "Intermediates, N = $(length(slcts[2]))"
labels[3]   = "Disks, N = $(length(slcts[3]))"
labels[4]   = "Ellipticals, N = $(length(slcts[4]))"
labpos      = Array{String}(undef, 4)
labpos     .= "upper left"
plot_hexbins_nxm( x, 
            y, 
            slcts, labels, 2, 2,
            #scale=0.0,
            cbarshift=2.5,
            label_pos=labpos, 
            outfile=outfile, 
            lognorm=true, gridres=(15,5), 
            ymin=miny, ymax=maxy,
            xmin=minx, xmax=maxx, 
            ylabel="Flip Angle [°]",
            xlabel="log( $(latexstring("\$M_\\textrm{*,post}\$")) )", fontsize=30,
            binned_median="y", plot_bval=" ",
            )
#


outnm   = joinpath(outputdir,"flipangle_vs_masspost_bvalstartsplit_09")
outfile = outnm*".pdf"
x   = log10.(as09["M"])
y   = as09["ϕ_flip"]
miny= 0#minimum(pyplottable(y)[findall(x->x .!== NaN, pyplottable(y))])
maxy= 180#maximum(pyplottable(y)[findall(x->x .!== NaN, pyplottable(y))])
minx= minimum(pyplottable(x)[findall(x->x .!== NaN, pyplottable(x))])
maxx= maximum(pyplottable(x)[findall(x->x .!== NaN, pyplottable(x))])
slcts       = Dict{Int64, Vector{Int64}}()
slcts[1]    = csw09[:transition]
slcts[2]    = idxmatch(csw09[:transition], csw09[:start_int ])
slcts[3]    = idxmatch(csw09[:transition], csw09[:start_disk])
slcts[4]    = idxmatch(csw09[:transition], csw09[:start_ell ])
labels      = Dict{Int64, String}()
labels[1]   = "Full Set, N = $(length(slcts[1]))"
labels[2]   = "Intermediates, N = $(length(slcts[2]))"
labels[3]   = "Disks, N = $(length(slcts[3]))"
labels[4]   = "Ellipticals, N = $(length(slcts[4]))"
labpos      = Array{String}(undef, 4)
labpos     .= "upper left"
plot_hexbins_nxm( x, 
            y, 
            slcts, labels, 2, 2,
            #scale=0.0,
            cbarshift=2.5,
            label_pos=labpos, 
            outfile=outfile, 
            lognorm=true, gridres=(15,5), 
            ymin=miny, ymax=maxy,
            xmin=minx, xmax=maxx, 
            ylabel="Flip Angle [°]",
            xlabel="log( $(latexstring("\$M_\\textrm{*,post}\$")) )", fontsize=30,
            binned_median="y", plot_bval=" ",
            )
#


outnm   = joinpath(outputdir,"flipangle_vs_masspost_bvalresultsplit_03")
outfile = outnm*".pdf"
x   = log10.(as03["M"])
y   = as03["ϕ_flip"]
miny= 0#minimum(pyplottable(y)[findall(x->x .!== NaN, pyplottable(y))])
maxy= 180#maximum(pyplottable(y)[findall(x->x .!== NaN, pyplottable(y))])
minx= minimum(pyplottable(x)[findall(x->x .!== NaN, pyplottable(x))])
maxx= maximum(pyplottable(x)[findall(x->x .!== NaN, pyplottable(x))])
slcts       = Dict{Int64, Vector{Int64}}()
slcts[1]    = csw03[:transition]
slcts[2]    = idxmatch(csw03[:transition], csw03[:result_int ])
slcts[3]    = idxmatch(csw03[:transition], csw03[:result_disk])
slcts[4]    = idxmatch(csw03[:transition], csw03[:result_ell ])
labels      = Dict{Int64, String}()
labels[1]   = "Full Set, N = $(length(slcts[1]))"
labels[2]   = "Intermediates, N = $(length(slcts[2]))"
labels[3]   = "Disks, N = $(length(slcts[3]))"
labels[4]   = "Ellipticals, N = $(length(slcts[4]))"
labpos      = Array{String}(undef, 4)
labpos     .= "upper left"
plot_hexbins_nxm( x, 
            y, 
            slcts, labels, 2, 2,
            #scale=0.0,
            cbarshift=2.5,
            label_pos=labpos, 
            outfile=outfile, 
            lognorm=true, gridres=(15,5), 
            ymin=miny, ymax=maxy,
            xmin=minx, xmax=maxx, 
            ylabel="Flip Angle [°]",
            xlabel="log( $(latexstring("\$M_\\textrm{*,post}\$")) )", fontsize=30,
            binned_median="y", plot_bval=" ",
            )
#


outnm   = joinpath(outputdir,"flipangle_vs_masspost_bvalresultsplit_09")
outfile = outnm*".pdf"
x   = log10.(as09["M"])
y   = as09["ϕ_flip"]
miny= 0#minimum(pyplottable(y)[findall(x->x .!== NaN, pyplottable(y))])
maxy= 180#maximum(pyplottable(y)[findall(x->x .!== NaN, pyplottable(y))])
minx= minimum(pyplottable(x)[findall(x->x .!== NaN, pyplottable(x))])
maxx= maximum(pyplottable(x)[findall(x->x .!== NaN, pyplottable(x))])
slcts       = Dict{Int64, Vector{Int64}}()
slcts[1]    = csw09[:transition]
slcts[2]    = idxmatch(csw09[:transition], csw09[:result_int ])
slcts[3]    = idxmatch(csw09[:transition], csw09[:result_disk])
slcts[4]    = idxmatch(csw09[:transition], csw09[:result_ell ])
labels      = Dict{Int64, String}()
labels[1]   = "Full Set, N = $(length(slcts[1]))"
labels[2]   = "Intermediates, N = $(length(slcts[2]))"
labels[3]   = "Disks, N = $(length(slcts[3]))"
labels[4]   = "Ellipticals, N = $(length(slcts[4]))"
labpos      = Array{String}(undef, 4)
labpos     .= "upper left"
plot_hexbins_nxm( x, 
            y, 
            slcts, labels, 2, 2,
            #scale=0.0,
            cbarshift=2.5,
            label_pos=labpos, 
            outfile=outfile, 
            lognorm=true, gridres=(15,5), 
            ymin=miny, ymax=maxy,
            xmin=minx, xmax=maxx, 
            ylabel="Flip Angle [°]",
            xlabel="log( $(latexstring("\$M_\\textrm{*,post}\$")) )", fontsize=30,
            binned_median="y", plot_bval=" ",
            )
#









##########################################################################################
# flip - merger - splits multimultiplot A-C * 1-3
# A1
outnm   = joinpath(outputdir,"flip_vs_M_MMratio_jsplit_A1")
outfile = outnm*".pdf"
x   = log10.((as03["M"] .- as03["ΔM"]) ./ as03["M_MM"])
y   = as03["ϕ_flip"]
miny= 0#minimum(pyplottable(y)[findall(x->x .!== NaN, pyplottable(y))])
maxy= 180#maximum(pyplottable(y)[findall(x->x .!== NaN, pyplottable(y))])
minx= log10(0.8)#minimum(pyplottable(x)[findall(x->x .!== NaN, pyplottable(x))])
maxx= 3#maximum(pyplottable(x)[findall(x->x .!== NaN, pyplottable(x))])
slcts       = Dict{Int64, Vector{Int64}}()
slcts[1]    = csw03[:transition]
slcts[2]    = idxmatch(csw03[:transition], csw03[:j_md ])
slcts[3]    = idxmatch(csw03[:transition], csw03[:j_lo])
slcts[4]    = idxmatch(csw03[:transition], csw03[:j_hi ])
labels      = Dict{Int64, String}()
labels[1]   = "Full Set, N = $(length(slcts[1])-count(isnan,pyplottable(x[slcts[1]])))"
labels[2]   = "Medium j, N = $(length(slcts[2])-count(isnan,pyplottable(x[slcts[2]])))"
labels[3]   = "Low j, N = $(length(slcts[3])-count(isnan,pyplottable(x[slcts[3]])))"
labels[4]   = "High j, N = $(length(slcts[4])-count(isnan,pyplottable(x[slcts[4]])))"
labpos      = Array{String}(undef, 4)
labpos     .= "upper right"
plot_hexbins_nxm( x, 
            y, 
            slcts, labels, 2, 2,
            #scale=0.0,
            cbarshift=2.5,
            label_pos=labpos, 
            outfile=outfile, 
            lognorm=true, gridres=(30,10), 
            ymin=miny, ymax=maxy,
            xmin=minx, xmax=maxx, 
            ylabel="Flip Angle [°]",
            xlabel="log( $(latexstring("\$M_{*\\textrm{,pre}}\\ /\\ M_{*\\textrm{,MM}}\$")) )", fontsize=30,
            binned_median="y", plot_mergers="x",
            )
#

# B1
outnm   = joinpath(outputdir,"flip_vs_M_MERGERSratio_jsplit_B1")
outfile = outnm*".pdf"
x   = log10.((as03["M"] .- as03["ΔM"]) ./ as03["M_MERGERS"])
y   = as03["ϕ_flip"]
slcts       = Dict{Int64, Vector{Int64}}()
slcts[1]    = csw03[:transition]
slcts[2]    = idxmatch(csw03[:transition], csw03[:j_md ])
slcts[3]    = idxmatch(csw03[:transition], csw03[:j_lo])
slcts[4]    = idxmatch(csw03[:transition], csw03[:j_hi ])
labels      = Dict{Int64, String}()
labels[1]   = "Full Set, N = $(length(slcts[1])-count(isnan,pyplottable(x[slcts[1]])))"
labels[2]   = "Medium j, N = $(length(slcts[2])-count(isnan,pyplottable(x[slcts[2]])))"
labels[3]   = "Low j, N = $(length(slcts[3])-count(isnan,pyplottable(x[slcts[3]])))"
labels[4]   = "High j, N = $(length(slcts[4])-count(isnan,pyplottable(x[slcts[4]])))"
labpos      = Array{String}(undef, 4)
labpos     .= "upper right"
plot_hexbins_nxm( x, 
            y, 
            slcts, labels, 2, 2,
            #scale=0.0,
            cbarshift=2.5,
            label_pos=labpos, 
            outfile=outfile, 
            lognorm=true, gridres=(30,10), 
            ymin=miny, ymax=maxy,
            xmin=minx, xmax=maxx, 
            ylabel="Flip Angle [°]",
            xlabel="log( $(latexstring("\$M_{*\\textrm{,pre}}\\ /\\ M_{*\\textrm{,Mergers}}\$")) )", fontsize=30,
            binned_median="y", plot_mergers="x",
            )
#

# C1
outnm   = joinpath(outputdir,"flip_vs_delMratio_jsplit_C1")
outfile = outnm*".pdf"
x   = log10.((as03["M"] .- as03["ΔM"]) ./ abs.(as03["ΔM"]))
y   = as03["ϕ_flip"]
slcts       = Dict{Int64, Vector{Int64}}()
slcts[1]    = csw03[:transition]
slcts[2]    = idxmatch(csw03[:transition], csw03[:j_md ])
slcts[3]    = idxmatch(csw03[:transition], csw03[:j_lo])
slcts[4]    = idxmatch(csw03[:transition], csw03[:j_hi ])
labels      = Dict{Int64, String}()
labels[1]   = "Full Set, N = $(length(slcts[1])-count(isnan,pyplottable(x[slcts[1]])))"
labels[2]   = "Medium j, N = $(length(slcts[2])-count(isnan,pyplottable(x[slcts[2]])))"
labels[3]   = "Low j, N = $(length(slcts[3])-count(isnan,pyplottable(x[slcts[3]])))"
labels[4]   = "High j, N = $(length(slcts[4])-count(isnan,pyplottable(x[slcts[4]])))"
labpos      = Array{String}(undef, 4)
labpos     .= "upper right"
plot_hexbins_nxm( x, 
            y, 
            slcts, labels, 2, 2,
            #scale=0.0,
            cbarshift=2.5,
            label_pos=labpos, 
            outfile=outfile, 
            lognorm=true, gridres=(30,10), 
            ymin=miny, ymax=maxy,
            xmin=minx, xmax=maxx, 
            ylabel="Flip Angle [°]",
            xlabel="log( $(latexstring("\$M_{*\\textrm{,pre}}\\ /\\ \\Delta M_{*}\$")) )", fontsize=30,
            binned_median="y", plot_mergers="x",
            )
#

# A2
outnm   = joinpath(outputdir,"flip_vs_M_MMratio_Jsplit_A2")
outfile = outnm*".pdf"
x   = log10.((as03["M"] .- as03["ΔM"]) ./ as03["M_MM"])
y   = as03["ϕ_flip"]
slcts       = Dict{Int64, Vector{Int64}}()
slcts[1]    = csw03[:transition]
slcts[2]    = idxmatch(csw03[:transition], csw03[:J_md ])
slcts[3]    = idxmatch(csw03[:transition], csw03[:J_lo])
slcts[4]    = idxmatch(csw03[:transition], csw03[:J_hi ])
labels      = Dict{Int64, String}()
labels[1]   = "Full Set, N = $(length(slcts[1])-count(isnan,pyplottable(x[slcts[1]])))"
labels[2]   = "Medium J, N = $(length(slcts[2])-count(isnan,pyplottable(x[slcts[2]])))"
labels[3]   = "Low J, N = $(length(slcts[3])-count(isnan,pyplottable(x[slcts[3]])))"
labels[4]   = "High J, N = $(length(slcts[4])-count(isnan,pyplottable(x[slcts[4]])))"
labpos      = Array{String}(undef, 4)
labpos     .= "upper right"
plot_hexbins_nxm( x, 
            y, 
            slcts, labels, 2, 2,
            #scale=0.0,
            cbarshift=2.5,
            label_pos=labpos, 
            outfile=outfile, 
            lognorm=true, gridres=(30,10), 
            ymin=miny, ymax=maxy,
            xmin=minx, xmax=maxx, 
            ylabel="Flip Angle [°]",
            xlabel="log( $(latexstring("\$M_{*\\textrm{,pre}}\\ /\\ M_{*\\textrm{,MM}}\$")) )", fontsize=30,
            binned_median="y", plot_mergers="x",
            )
#

# B2
outnm   = joinpath(outputdir,"flip_vs_M_MERGERSratio_Jsplit_B2")
outfile = outnm*".pdf"
x   = log10.((as03["M"] .- as03["ΔM"]) ./ as03["M_MERGERS"])
y   = as03["ϕ_flip"]
slcts       = Dict{Int64, Vector{Int64}}()
slcts[1]    = csw03[:transition]
slcts[2]    = idxmatch(csw03[:transition], csw03[:j_md ])
slcts[3]    = idxmatch(csw03[:transition], csw03[:j_lo])
slcts[4]    = idxmatch(csw03[:transition], csw03[:j_hi ])
labels      = Dict{Int64, String}()
labels[1]   = "Full Set, N = $(length(slcts[1])-count(isnan,pyplottable(x[slcts[1]])))"
labels[2]   = "Medium J, N = $(length(slcts[2])-count(isnan,pyplottable(x[slcts[2]])))"
labels[3]   = "Low J, N = $(length(slcts[3])-count(isnan,pyplottable(x[slcts[3]])))"
labels[4]   = "High J, N = $(length(slcts[4])-count(isnan,pyplottable(x[slcts[4]])))"
labpos      = Array{String}(undef, 4)
labpos     .= "upper right"
plot_hexbins_nxm( x, 
            y, 
            slcts, labels, 2, 2,
            #scale=0.0,
            cbarshift=2.5,
            label_pos=labpos, 
            outfile=outfile, 
            lognorm=true, gridres=(30,10), 
            ymin=miny, ymax=maxy,
            xmin=minx, xmax=maxx, 
            ylabel="Flip Angle [°]",
            xlabel="log( $(latexstring("\$M_{*\\textrm{,pre}}\\ /\\ M_{*\\textrm{,Mergers}}\$")) )", fontsize=30,
            binned_median="y", plot_mergers="x",
            )
#

# C2
outnm   = joinpath(outputdir,"flip_vs_delMratio_Jsplit_C2")
outfile = outnm*".pdf"
x   = log10.((as03["M"] .- as03["ΔM"]) ./ abs.(as03["ΔM"]))
y   = as03["ϕ_flip"]
slcts       = Dict{Int64, Vector{Int64}}()
slcts[1]    = csw03[:transition]
slcts[2]    = idxmatch(csw03[:transition], csw03[:J_md ])
slcts[3]    = idxmatch(csw03[:transition], csw03[:J_lo])
slcts[4]    = idxmatch(csw03[:transition], csw03[:J_hi ])
labels      = Dict{Int64, String}()
labels[1]   = "Full Set, N = $(length(slcts[1])-count(isnan,pyplottable(x[slcts[1]])))"
labels[2]   = "Medium J, N = $(length(slcts[2])-count(isnan,pyplottable(x[slcts[2]])))"
labels[3]   = "Low J, N = $(length(slcts[3])-count(isnan,pyplottable(x[slcts[3]])))"
labels[4]   = "High J, N = $(length(slcts[4])-count(isnan,pyplottable(x[slcts[4]])))"
labpos      = Array{String}(undef, 4)
labpos     .= "upper right"
plot_hexbins_nxm( x, 
            y, 
            slcts, labels, 2, 2,
            #scale=0.0,
            cbarshift=2.5,
            label_pos=labpos, 
            outfile=outfile, 
            lognorm=true, gridres=(30,10), 
            ymin=miny, ymax=maxy,
            xmin=minx, xmax=maxx, 
            ylabel="Flip Angle [°]",
            xlabel="log( $(latexstring("\$M_{*\\textrm{,pre}}\\ /\\ \\Delta M_{*}\$")) )", fontsize=30,
            binned_median="y", plot_mergers="x",
            )
#


# A3
outnm   = joinpath(outputdir,"flip_vs_M_MMratio_bsplit_A3")
outfile = outnm*".pdf"
x   = log10.((as03["M"] .- as03["ΔM"]) ./ as03["M_MM"])
y   = as03["ϕ_flip"]
slcts       = Dict{Int64, Vector{Int64}}()
slcts[1]    = csw03[:transition]
slcts[2]    = idxmatch(csw03[:transition], csw03[:start_int ])
slcts[3]    = idxmatch(csw03[:transition], csw03[:start_ell])
slcts[4]    = idxmatch(csw03[:transition], csw03[:start_disk ])
labels      = Dict{Int64, String}()
labels[1]   = "Full Set, N = $(length(slcts[1])-count(isnan,pyplottable(x[slcts[1]])))"
labels[2]   = "Intermediates, N = $(length(slcts[2])-count(isnan,pyplottable(x[slcts[2]])))"
labels[3]   = "Ellipticals, N = $(length(slcts[3])-count(isnan,pyplottable(x[slcts[3]])))"
labels[4]   = "Disks, N = $(length(slcts[4])-count(isnan,pyplottable(x[slcts[4]])))"
labpos      = Array{String}(undef, 4)
labpos     .= "upper right"
plot_hexbins_nxm( x, 
            y, 
            slcts, labels, 2, 2,
            #scale=0.0,
            cbarshift=2.5,
            label_pos=labpos, 
            outfile=outfile, 
            lognorm=true, gridres=(30,10), 
            ymin=miny, ymax=maxy,
            xmin=minx, xmax=maxx, 
            ylabel="Flip Angle [°]",
            xlabel="log( $(latexstring("\$M_{*\\textrm{,pre}}\\ /\\ M_{*\\textrm{,MM}}\$")) )", fontsize=30,
            binned_median="y", plot_mergers="x",
            )
#

# B3
outnm   = joinpath(outputdir,"flip_vs_M_MERGERSratio_bsplit_B3")
outfile = outnm*".pdf"
x   = log10.((as03["M"] .- as03["ΔM"]) ./ as03["M_MERGERS"])
y   = as03["ϕ_flip"]
slcts       = Dict{Int64, Vector{Int64}}()
slcts[1]    = csw03[:transition]
slcts[2]    = idxmatch(csw03[:transition], csw03[:start_int ])
slcts[3]    = idxmatch(csw03[:transition], csw03[:start_ell])
slcts[4]    = idxmatch(csw03[:transition], csw03[:start_disk ])
labels      = Dict{Int64, String}()
labels[1]   = "Full Set, N = $(length(slcts[1])-count(isnan,pyplottable(x[slcts[1]])))"
labels[2]   = "Intermediates, N = $(length(slcts[2])-count(isnan,pyplottable(x[slcts[2]])))"
labels[3]   = "Ellipticals, N = $(length(slcts[3])-count(isnan,pyplottable(x[slcts[3]])))"
labels[4]   = "Disks, N = $(length(slcts[4])-count(isnan,pyplottable(x[slcts[4]])))"
labpos      = Array{String}(undef, 4)
labpos     .= "upper right"
plot_hexbins_nxm( x, 
            y, 
            slcts, labels, 2, 2,
            #scale=0.0,
            cbarshift=2.5,
            label_pos=labpos, 
            outfile=outfile, 
            lognorm=true, gridres=(30,10), 
            ymin=miny, ymax=maxy,
            xmin=minx, xmax=maxx, 
            ylabel="Flip Angle [°]",
            xlabel="log( $(latexstring("\$M_{*\\textrm{,pre}}\\ /\\ M_{*\\textrm{,Mergers}}\$")) )", fontsize=30,
            binned_median="y", plot_mergers="x",
            )
#

# C3
outnm   = joinpath(outputdir,"flip_vs_delMratio_bsplit_C3")
outfile = outnm*".pdf"
x   = log10.((as03["M"] .- as03["ΔM"]) ./ abs.(as03["ΔM"]))
y   = as03["ϕ_flip"]
slcts       = Dict{Int64, Vector{Int64}}()
slcts[1]    = csw03[:transition]
slcts[2]    = idxmatch(csw03[:transition], csw03[:start_int ])
slcts[3]    = idxmatch(csw03[:transition], csw03[:start_ell])
slcts[4]    = idxmatch(csw03[:transition], csw03[:start_disk ])
labels      = Dict{Int64, String}()
labels[1]   = "Full Set, N = $(length(slcts[1])-count(isnan,pyplottable(x[slcts[1]])))"
labels[2]   = "Intermediates, N = $(length(slcts[2])-count(isnan,pyplottable(x[slcts[2]])))"
labels[3]   = "Ellipticals, N = $(length(slcts[3])-count(isnan,pyplottable(x[slcts[3]])))"
labels[4]   = "Disks, N = $(length(slcts[4])-count(isnan,pyplottable(x[slcts[4]])))"
labpos      = Array{String}(undef, 4)
labpos     .= "upper right"
plot_hexbins_nxm( x, 
            y, 
            slcts, labels, 2, 2,
            #scale=0.0,
            cbarshift=2.5,
            label_pos=labpos, 
            outfile=outfile, 
            lognorm=true, gridres=(30,10), 
            ymin=miny, ymax=maxy,
            xmin=minx, xmax=maxx, 
            ylabel="Flip Angle [°]",
            xlabel="log( $(latexstring("\$M_{*\\textrm{,pre}}\\ /\\ \\Delta M_{*}\$")) )", fontsize=30,
            binned_median="y", plot_mergers="x",
            )
#












# repeat for 09 timesteps

# A1
outnm   = joinpath(outputdir,"flip_vs_M_MMratio_jsplit_A1_09")
outfile = outnm*".pdf"
x   = log10.((as09["M"] .- as09["ΔM"]) ./ as09["M_MM"])
y   = as09["ϕ_flip"]
slcts       = Dict{Int64, Vector{Int64}}()
slcts[1]    = csw09[:transition]
slcts[2]    = idxmatch(csw09[:transition], csw09[:j_md ])
slcts[3]    = idxmatch(csw09[:transition], csw09[:j_lo])
slcts[4]    = idxmatch(csw09[:transition], csw09[:j_hi ])
labels      = Dict{Int64, String}()
labels[1]   = "Full Set, N = $(length(slcts[1])-count(isnan,pyplottable(x[slcts[1]])))"
labels[2]   = "Medium j, N = $(length(slcts[2])-count(isnan,pyplottable(x[slcts[2]])))"
labels[3]   = "Low j, N = $(length(slcts[3])-count(isnan,pyplottable(x[slcts[3]])))"
labels[4]   = "High j, N = $(length(slcts[4])-count(isnan,pyplottable(x[slcts[4]])))"
labpos      = Array{String}(undef, 4)
labpos     .= "upper right"
plot_hexbins_nxm( x, 
            y, 
            slcts, labels, 2, 2,
            #scale=0.0,
            cbarshift=2.5,
            label_pos=labpos, 
            outfile=outfile, 
            lognorm=true, gridres=(30,10), 
            ymin=miny, ymax=maxy,
            xmin=minx, xmax=maxx, 
            ylabel="Flip Angle [°]",
            xlabel="log( $(latexstring("\$M_{*\\textrm{,pre}}\\ /\\ M_{*\\textrm{,MM}}\$")) )", fontsize=30,
            binned_median="y", plot_mergers="x",
            )
#

# B1
outnm   = joinpath(outputdir,"flip_vs_M_MERGERSratio_jsplit_B1_09")
outfile = outnm*".pdf"
x   = log10.((as09["M"] .- as09["ΔM"]) ./ as09["M_MERGERS"])
y   = as09["ϕ_flip"]
slcts       = Dict{Int64, Vector{Int64}}()
slcts[1]    = csw09[:transition]
slcts[2]    = idxmatch(csw09[:transition], csw09[:j_md ])
slcts[3]    = idxmatch(csw09[:transition], csw09[:j_lo])
slcts[4]    = idxmatch(csw09[:transition], csw09[:j_hi ])
labels      = Dict{Int64, String}()
labels[1]   = "Full Set, N = $(length(slcts[1])-count(isnan,pyplottable(x[slcts[1]])))"
labels[2]   = "Medium j, N = $(length(slcts[2])-count(isnan,pyplottable(x[slcts[2]])))"
labels[3]   = "Low j, N = $(length(slcts[3])-count(isnan,pyplottable(x[slcts[3]])))"
labels[4]   = "High j, N = $(length(slcts[4])-count(isnan,pyplottable(x[slcts[4]])))"
labpos      = Array{String}(undef, 4)
labpos     .= "upper right"
plot_hexbins_nxm( x, 
            y, 
            slcts, labels, 2, 2,
            #scale=0.0,
            cbarshift=2.5,
            label_pos=labpos, 
            outfile=outfile, 
            lognorm=true, gridres=(30,10), 
            ymin=miny, ymax=maxy,
            xmin=minx, xmax=maxx, 
            ylabel="Flip Angle [°]",
            xlabel="log( $(latexstring("\$M_{*\\textrm{,pre}}\\ /\\ M_{*\\textrm{,Mergers}}\$")) )", fontsize=30,
            binned_median="y", plot_mergers="x",
            )
#

# C1
outnm   = joinpath(outputdir,"flip_vs_delMratio_jsplit_C1_09")
outfile = outnm*".pdf"
x   = log10.((as09["M"] .- as09["ΔM"]) ./ abs.(as09["ΔM"]))
y   = as09["ϕ_flip"]
slcts       = Dict{Int64, Vector{Int64}}()
slcts[1]    = csw09[:transition]
slcts[2]    = idxmatch(csw09[:transition], csw09[:j_md ])
slcts[3]    = idxmatch(csw09[:transition], csw09[:j_lo])
slcts[4]    = idxmatch(csw09[:transition], csw09[:j_hi ])
labels      = Dict{Int64, String}()
labels[1]   = "Full Set, N = $(length(slcts[1])-count(isnan,pyplottable(x[slcts[1]])))"
labels[2]   = "Medium j, N = $(length(slcts[2])-count(isnan,pyplottable(x[slcts[2]])))"
labels[3]   = "Low j, N = $(length(slcts[3])-count(isnan,pyplottable(x[slcts[3]])))"
labels[4]   = "High j, N = $(length(slcts[4])-count(isnan,pyplottable(x[slcts[4]])))"
labpos      = Array{String}(undef, 4)
labpos     .= "upper right"
plot_hexbins_nxm( x, 
            y, 
            slcts, labels, 2, 2,
            #scale=0.0,
            cbarshift=2.5,
            label_pos=labpos, 
            outfile=outfile, 
            lognorm=true, gridres=(30,10), 
            ymin=miny, ymax=maxy,
            xmin=minx, xmax=maxx, 
            ylabel="Flip Angle [°]",
            xlabel="log( $(latexstring("\$M_{*\\textrm{,pre}}\\ /\\ \\Delta M_{*}\$")) )", fontsize=30,
            binned_median="y", plot_mergers="x",
            )
#

# A2
outnm   = joinpath(outputdir,"flip_vs_M_MMratio_Jsplit_A2_09")
outfile = outnm*".pdf"
x   = log10.((as09["M"] .- as09["ΔM"]) ./ as09["M_MM"])
y   = as09["ϕ_flip"]
slcts       = Dict{Int64, Vector{Int64}}()
slcts[1]    = csw09[:transition]
slcts[2]    = idxmatch(csw09[:transition], csw09[:J_md ])
slcts[3]    = idxmatch(csw09[:transition], csw09[:J_lo])
slcts[4]    = idxmatch(csw09[:transition], csw09[:J_hi ])
labels      = Dict{Int64, String}()
labels[1]   = "Full Set, N = $(length(slcts[1])-count(isnan,pyplottable(x[slcts[1]])))"
labels[2]   = "Medium J, N = $(length(slcts[2])-count(isnan,pyplottable(x[slcts[2]])))"
labels[3]   = "Low J, N = $(length(slcts[3])-count(isnan,pyplottable(x[slcts[3]])))"
labels[4]   = "High J, N = $(length(slcts[4])-count(isnan,pyplottable(x[slcts[4]])))"
labpos      = Array{String}(undef, 4)
labpos     .= "upper right"
plot_hexbins_nxm( x, 
            y, 
            slcts, labels, 2, 2,
            #scale=0.0,
            cbarshift=2.5,
            label_pos=labpos, 
            outfile=outfile, 
            lognorm=true, gridres=(30,10), 
            ymin=miny, ymax=maxy,
            xmin=minx, xmax=maxx, 
            ylabel="Flip Angle [°]",
            xlabel="log( $(latexstring("\$M_{*\\textrm{,pre}}\\ /\\ M_{*\\textrm{,MM}}\$")) )", fontsize=30,
            binned_median="y", plot_mergers="x",
            )
#

# B2
outnm   = joinpath(outputdir,"flip_vs_M_MERGERSratio_Jsplit_B2_09")
outfile = outnm*".pdf"
x   = log10.((as09["M"] .- as09["ΔM"]) ./ as09["M_MERGERS"])
y   = as09["ϕ_flip"]
slcts       = Dict{Int64, Vector{Int64}}()
slcts[1]    = csw09[:transition]
slcts[2]    = idxmatch(csw09[:transition], csw09[:j_md ])
slcts[3]    = idxmatch(csw09[:transition], csw09[:j_lo])
slcts[4]    = idxmatch(csw09[:transition], csw09[:j_hi ])
labels      = Dict{Int64, String}()
labels[1]   = "Full Set, N = $(length(slcts[1])-count(isnan,pyplottable(x[slcts[1]])))"
labels[2]   = "Medium J, N = $(length(slcts[2])-count(isnan,pyplottable(x[slcts[2]])))"
labels[3]   = "Low J, N = $(length(slcts[3])-count(isnan,pyplottable(x[slcts[3]])))"
labels[4]   = "High J, N = $(length(slcts[4])-count(isnan,pyplottable(x[slcts[4]])))"
labpos      = Array{String}(undef, 4)
labpos     .= "upper right"
plot_hexbins_nxm( x, 
            y, 
            slcts, labels, 2, 2,
            #scale=0.0,
            cbarshift=2.5,
            label_pos=labpos, 
            outfile=outfile, 
            lognorm=true, gridres=(30,10), 
            ymin=miny, ymax=maxy,
            xmin=minx, xmax=maxx, 
            ylabel="Flip Angle [°]",
            xlabel="log( $(latexstring("\$M_{*\\textrm{,pre}}\\ /\\ M_{*\\textrm{,Mergers}}\$")) )", fontsize=30,
            binned_median="y", plot_mergers="x",
            )
#

# C2
outnm   = joinpath(outputdir,"flip_vs_delMratio_Jsplit_C2_09")
outfile = outnm*".pdf"
x   = log10.((as09["M"] .- as09["ΔM"]) ./ abs.(as09["ΔM"]))
y   = as09["ϕ_flip"]
slcts       = Dict{Int64, Vector{Int64}}()
slcts[1]    = csw09[:transition]
slcts[2]    = idxmatch(csw09[:transition], csw09[:J_md ])
slcts[3]    = idxmatch(csw09[:transition], csw09[:J_lo])
slcts[4]    = idxmatch(csw09[:transition], csw09[:J_hi ])
labels      = Dict{Int64, String}()
labels[1]   = "Full Set, N = $(length(slcts[1])-count(isnan,pyplottable(x[slcts[1]])))"
labels[2]   = "Medium J, N = $(length(slcts[2])-count(isnan,pyplottable(x[slcts[2]])))"
labels[3]   = "Low J, N = $(length(slcts[3])-count(isnan,pyplottable(x[slcts[3]])))"
labels[4]   = "High J, N = $(length(slcts[4])-count(isnan,pyplottable(x[slcts[4]])))"
labpos      = Array{String}(undef, 4)
labpos     .= "upper right"
plot_hexbins_nxm( x, 
            y, 
            slcts, labels, 2, 2,
            #scale=0.0,
            cbarshift=2.5,
            label_pos=labpos, 
            outfile=outfile, 
            lognorm=true, gridres=(30,10), 
            ymin=miny, ymax=maxy,
            xmin=minx, xmax=maxx, 
            ylabel="Flip Angle [°]",
            xlabel="log( $(latexstring("\$M_{*\\textrm{,pre}}\\ /\\ \\Delta M_{*}\$")) )", fontsize=30,
            binned_median="y", plot_mergers="x",
            )
#


# A3
outnm   = joinpath(outputdir,"flip_vs_M_MMratio_bsplit_A3_09")
outfile = outnm*".pdf"
x   = log10.((as09["M"] .- as09["ΔM"]) ./ as09["M_MM"])
y   = as09["ϕ_flip"]
slcts       = Dict{Int64, Vector{Int64}}()
slcts[1]    = csw09[:transition]
slcts[2]    = idxmatch(csw09[:transition], csw09[:start_int ])
slcts[3]    = idxmatch(csw09[:transition], csw09[:start_ell])
slcts[4]    = idxmatch(csw09[:transition], csw09[:start_disk ])
labels      = Dict{Int64, String}()
labels[1]   = "Full Set, N = $(length(slcts[1])-count(isnan,pyplottable(x[slcts[1]])))"
labels[2]   = "Intermediates, N = $(length(slcts[2])-count(isnan,pyplottable(x[slcts[2]])))"
labels[3]   = "Ellipticals, N = $(length(slcts[3])-count(isnan,pyplottable(x[slcts[3]])))"
labels[4]   = "Disks, N = $(length(slcts[4])-count(isnan,pyplottable(x[slcts[4]])))"
labpos      = Array{String}(undef, 4)
labpos     .= "upper right"
plot_hexbins_nxm( x, 
            y, 
            slcts, labels, 2, 2,
            #scale=0.0,
            cbarshift=2.5,
            label_pos=labpos, 
            outfile=outfile, 
            lognorm=true, gridres=(30,10), 
            ymin=miny, ymax=maxy,
            xmin=minx, xmax=maxx, 
            ylabel="Flip Angle [°]",
            xlabel="log( $(latexstring("\$M_{*\\textrm{,pre}}\\ /\\ M_{*\\textrm{,MM}}\$")) )", fontsize=30,
            binned_median="y", plot_mergers="x",
            )
#

# B3
outnm   = joinpath(outputdir,"flip_vs_M_MERGERSratio_bsplit_B3_09")
outfile = outnm*".pdf"
x   = log10.((as09["M"] .- as09["ΔM"]) ./ as09["M_MERGERS"])
y   = as09["ϕ_flip"]
slcts       = Dict{Int64, Vector{Int64}}()
slcts[1]    = csw09[:transition]
slcts[2]    = idxmatch(csw09[:transition], csw09[:start_int ])
slcts[3]    = idxmatch(csw09[:transition], csw09[:start_ell])
slcts[4]    = idxmatch(csw09[:transition], csw09[:start_disk ])
labels      = Dict{Int64, String}()
labels[1]   = "Full Set, N = $(length(slcts[1])-count(isnan,pyplottable(x[slcts[1]])))"
labels[2]   = "Intermediates, N = $(length(slcts[2])-count(isnan,pyplottable(x[slcts[2]])))"
labels[3]   = "Ellipticals, N = $(length(slcts[3])-count(isnan,pyplottable(x[slcts[3]])))"
labels[4]   = "Disks, N = $(length(slcts[4])-count(isnan,pyplottable(x[slcts[4]])))"
labpos      = Array{String}(undef, 4)
labpos     .= "upper right"
plot_hexbins_nxm( x, 
            y, 
            slcts, labels, 2, 2,
            #scale=0.0,
            cbarshift=2.5,
            label_pos=labpos, 
            outfile=outfile, 
            lognorm=true, gridres=(30,10), 
            ymin=miny, ymax=maxy,
            xmin=minx, xmax=maxx, 
            ylabel="Flip Angle [°]",
            xlabel="log( $(latexstring("\$M_{*\\textrm{,pre}}\\ /\\ M_{*\\textrm{,Mergers}}\$")) )", fontsize=30,
            binned_median="y", plot_mergers="x",
            )
#

# C3
outnm   = joinpath(outputdir,"flip_vs_delMratio_bsplit_C3_09")
outfile = outnm*".pdf"
x   = log10.((as09["M"] .- as09["ΔM"]) ./ abs.(as09["ΔM"]))
y   = as09["ϕ_flip"]
slcts       = Dict{Int64, Vector{Int64}}()
slcts[1]    = csw09[:transition]
slcts[2]    = idxmatch(csw09[:transition], csw09[:start_int ])
slcts[3]    = idxmatch(csw09[:transition], csw09[:start_ell])
slcts[4]    = idxmatch(csw09[:transition], csw09[:start_disk ])
labels      = Dict{Int64, String}()
labels[1]   = "Full Set, N = $(length(slcts[1])-count(isnan,pyplottable(x[slcts[1]])))"
labels[2]   = "Intermediates, N = $(length(slcts[2])-count(isnan,pyplottable(x[slcts[2]])))"
labels[3]   = "Ellipticals, N = $(length(slcts[3])-count(isnan,pyplottable(x[slcts[3]])))"
labels[4]   = "Disks, N = $(length(slcts[4])-count(isnan,pyplottable(x[slcts[4]])))"
labpos      = Array{String}(undef, 4)
labpos     .= "upper right"
plot_hexbins_nxm( x, 
            y, 
            slcts, labels, 2, 2,
            #scale=0.0,
            cbarshift=2.5,
            label_pos=labpos, 
            outfile=outfile, 
            lognorm=true, gridres=(30,10), 
            ymin=miny, ymax=maxy,
            xmin=minx, xmax=maxx, 
            ylabel="Flip Angle [°]",
            xlabel="log( $(latexstring("\$M_{*\\textrm{,pre}}\\ /\\ \\Delta M_{*}\$")) )", fontsize=30,
            binned_median="y", plot_mergers="x",
            )
#








####################################################################################################
# 
outnm   = "flipangle_vs_lookbacktime_03"
y       = pyplottable(as03["ϕ_flip"][csw03[:transition]])
x       = pyplottable(as03["lookbacktime"][csw03[:transition]])
plot_hexbins( x, y, 
    outfile=joinpath(outputdir,outnm*"_hxcount.pdf"),
    scale=15, gridres=(15,7), xscale="linear", yscale="linear", 
    xlabel="Lookback Time [Gyr]", 
    ylabel="Flip Angle [°]", 
    grid=false,
    #xmin=0, xmax=180,# ymin=ymin, ymax=ymax,
    lognorm=true,
    binned_median = "y", plot_bval=" ", 
    calc_norm=true,    # @pyplottable
    )
#

outnm   = "flipangle_vs_lookbacktime_09"
y       = pyplottable(as09["ϕ_flip"][csw09[:transition]])
x       = pyplottable(as09["lookbacktime"][csw09[:transition]])
plot_hexbins( x, y, 
    outfile=joinpath(outputdir,outnm*"_hxcount.pdf"),
    scale=15, gridres=(15,7), xscale="linear", yscale="linear", 
    xlabel="Lookback Time [Gyr]", 
    ylabel="Flip Angle [°]", 
    grid=false,
    xmax=12,#xmin=0, xmax=180,# ymin=ymin, ymax=ymax,
    lognorm=true,
    binned_median = "y", plot_bval=" ", 
    calc_norm=true,    # @pyplottable
    )
#







#####################################################################################################################
# flip vs mm ratio bval & z split
# 03
outnm   = joinpath(outputdir,"flip_vs_Mmergers_bzsplit_03")
outfile = outnm*".pdf"
x   = log10.((as03["M"] .- as03["ΔM"]) ./ abs.(as03["M_MERGERS"]))
y   = as03["ϕ_flip"]
miny= 0#minimum(pyplottable(y)[findall(x->x .!== NaN, pyplottable(y))])
maxy= 180#maximum(pyplottable(y)[findall(x->x .!== NaN, pyplottable(y))])
minx= log10(0.8)#minimum(pyplottable(x)[findall(x->x .!== NaN, pyplottable(x))])
maxx= 3#maximum(pyplottable(x)[findall(x->x .!== NaN, pyplottable(x))])
slcts       = Dict{Int64, Vector{Int64}}()
slcts[1]    = idxmatch(idxmatch(csw03[:transition], csw03[:start_ell])  , csw03[:z_lo])#split15_l])
slcts[2]    = idxmatch(idxmatch(csw03[:transition], csw03[:start_int ]) , csw03[:z_lo])#split15_l])
slcts[3]    = idxmatch(idxmatch(csw03[:transition], csw03[:start_disk ]), csw03[:z_lo])#split15_l])
slcts[4]    = idxmatch(idxmatch(csw03[:transition], csw03[:start_ell])  , csw03[:z_hi])#split15_g])
slcts[5]    = idxmatch(idxmatch(csw03[:transition], csw03[:start_int ]) , csw03[:z_hi])#split15_g])
slcts[6]    = idxmatch(idxmatch(csw03[:transition], csw03[:start_disk ]), csw03[:z_hi])#split15_g])
labels      = Dict{Int64, String}()
labels[1]   = "Ellipticals, $(latexstring("\$z<0.45\$")), N = $(length(slcts[1]))" # 1.5
labels[2]   = "Intermediates, $(latexstring("\$z<0.45\$")), N = $(length(slcts[2]))"
labels[3]   = "Disks, $(latexstring("\$z<0.45\$")), N = $(length(slcts[3]))"
labels[4]   = "Ellipticals, $(latexstring("\$z>1\$")), N = $(length(slcts[4]))"
labels[5]   = "Intermediates, $(latexstring("\$z>1\$")), N = $(length(slcts[5]))"
labels[6]   = "Disks, $(latexstring("\$z>1\$")), N = $(length(slcts[6]))"
labpos      = Array{String}(undef, 6)
labpos     .= "upper right"
plot_hexbins_nxm( x, 
            y, 
            slcts, labels, 3, 2,
            #scale=0.0,
            cbarshift=2.5,
            label_pos=labpos, 
            outfile=outfile, 
            lognorm=true, gridres=(15,4), 
            ymin=miny, ymax=maxy,
            xmin=minx, xmax=maxx, 
            ylabel="Flip Angle [°]",
            xlabel="log( $(latexstring("\$M_{*\\textrm{,pre}}\\ /\\ M_{*\\textrm{,Mergers}}\$")) )", fontsize=30,
            binned_median="y", plot_mergers="x",
            )
#


# 09
outnm   = joinpath(outputdir,"flip_vs_Mmergers_bzsplit_09")
outfile = outnm*".pdf"
x   = log10.((as09["M"] .- as09["ΔM"]) ./ abs.(as09["M_MERGERS"]))
y   = as09["ϕ_flip"]
miny= 0#minimum(pyplottable(y)[findall(x->x .!== NaN, pyplottable(y))])
maxy= 180#maximum(pyplottable(y)[findall(x->x .!== NaN, pyplottable(y))])
minx= log10(0.8)#minimum(pyplottable(x)[findall(x->x .!== NaN, pyplottable(x))])
maxx= 3#maximum(pyplottable(x)[findall(x->x .!== NaN, pyplottable(x))])
slcts       = Dict{Int64, Vector{Int64}}()
slcts[1]    = idxmatch(idxmatch(csw09[:transition], csw09[:start_ell])  , csw09[:z_lo])#split15_l])
slcts[2]    = idxmatch(idxmatch(csw09[:transition], csw09[:start_int ]) , csw09[:z_lo])#split15_l])
slcts[3]    = idxmatch(idxmatch(csw09[:transition], csw09[:start_disk ]), csw09[:z_lo])#split15_l])
slcts[4]    = idxmatch(idxmatch(csw09[:transition], csw09[:start_ell])  , csw09[:z_hi])#split15_g])
slcts[5]    = idxmatch(idxmatch(csw09[:transition], csw09[:start_int ]) , csw09[:z_hi])#split15_g])
slcts[6]    = idxmatch(idxmatch(csw09[:transition], csw09[:start_disk ]), csw09[:z_hi])#split15_g])
labels      = Dict{Int64, String}()
labels[1]   = "Ellipticals, $(latexstring("\$z<0.45\$")), N = $(length(slcts[1]))" # 1.5
labels[2]   = "Intermediates, $(latexstring("\$z<0.45\$")), N = $(length(slcts[2]))"
labels[3]   = "Disks, $(latexstring("\$z<0.45\$")), N = $(length(slcts[3]))"
labels[4]   = "Ellipticals, $(latexstring("\$z>1\$")), N = $(length(slcts[4]))"
labels[5]   = "Intermediates, $(latexstring("\$z>1\$")), N = $(length(slcts[5]))"
labels[6]   = "Disks, $(latexstring("\$z>1\$")), N = $(length(slcts[6]))"
labpos      = Array{String}(undef, 6)
labpos     .= "upper right"
plot_hexbins_nxm( x, 
            y, 
            slcts, labels, 3, 2,
            #scale=0.0,
            cbarshift=2.5,
            label_pos=labpos, 
            outfile=outfile, 
            lognorm=true, gridres=(15,4), 
            ymin=miny, ymax=maxy,
            xmin=minx, xmax=maxx, 
            ylabel="Flip Angle [°]",
            xlabel="log( $(latexstring("\$M_{*\\textrm{,pre}}\\ /\\ M_{*\\textrm{,Mergers}}\$")) )", fontsize=30,
            binned_median="y", plot_mergers="x",
            )
#







############################################################################################
# Welker 2014, welker

# z_welker, 03, stars
outnm   = joinpath(outputdir,"welker2014_Mstars_zwelker_03")
outfile = outnm*".pdf"
scale=0.5
N_bins=5
pltcm       = pyimport("matplotlib.cm")
pltcolors   = pyimport("matplotlib.colors")
colormap    = pltcm.get_cmap(name="gnuplot")
fig, ax = subplots()
n_1, x_1        = PyPlot.hist( cosd.(pyplottable(as03["ϕ_flip"][idxmatch(idxmatch(csw03[:mergerSTARS_acc_eq0] , csw03[:transition]), csw03[:z_welker])])), range=(0,1), bins=N_bins, density=true, log=false)
n_2, x_2        = PyPlot.hist( cosd.(pyplottable(as03["ϕ_flip"][idxmatch(idxmatch(csw03[:mergerSTARS_acc_g0_l5 ], csw03[:transition]), csw03[:z_welker])])), range=(0,1), bins=N_bins, density=true, log=false)
n_3, x_3        = PyPlot.hist( cosd.(pyplottable(as03["ϕ_flip"][idxmatch(idxmatch(csw03[:mergerSTARS_acc_g5_l10 ], csw03[:transition]), csw03[:z_welker])])), range=(0,1), bins=N_bins, density=true, log=false)
n_4, x_4        = PyPlot.hist( cosd.(pyplottable(as03["ϕ_flip"][idxmatch(idxmatch(csw03[:mergerSTARS_acc_g10], csw03[:transition]), csw03[:z_welker])])), range=(0,1), bins=N_bins, density=true, log=false)
bcenters_1     = 0.5 .* (x_1[2:end]  .+ x_1[1:end-1])
bcenters_2     = 0.5 .* (x_2[2:end]  .+ x_2[1:end-1])
bcenters_3     = 0.5 .* (x_3[2:end]  .+ x_3[1:end-1])
bcenters_4     = 0.5 .* (x_4[2:end]  .+ x_4[1:end-1])
clf()
PyPlot.matplotlib[:rc]("text", usetex=true)
PyPlot.matplotlib[:rc]("font", family="serif", serif=["Computer Modern"])
fig, ax = subplots()
ax.grid(true)

ax.plot(bcenters_1 , log10.(n_1), ".-", ms=10, label="$(latexstring("\$\\delta m = 0\$\\%")), N = $(length(idxmatch(idxmatch(csw03[:mergerSTARS_acc_eq0] , csw03[:transition]), csw03[:z_welker])))", color="red")
ax.plot(bcenters_2 , log10.(n_2), ".-", ms=10, label="$(latexstring("\$0 < \\delta m < 5\$\\%")), N = $(length(idxmatch(idxmatch(csw03[:mergerSTARS_acc_g0_l5 ], csw03[:transition]), csw03[:z_welker])))", color="orange"   )
ax.plot(bcenters_3 , log10.(n_3), ".-", ms=10, label="$(latexstring("\$5 < \\delta m < 10\$\\%")), N = $(length(idxmatch(idxmatch(csw03[:mergerSTARS_acc_g5_l10 ], csw03[:transition]), csw03[:z_welker])))", color="purple"   )
ax.plot(bcenters_4 , log10.(n_4), ".-", ms=10, label="$(latexstring("\$\\delta m > 10\$\\%")), N = $(length(idxmatch(idxmatch(csw03[:mergerSTARS_acc_g10], csw03[:transition]), csw03[:z_welker])))", color="mediumblue"   )

ax.set_xlabel("cos($(latexstring("\$\\phi\$")))")# $(latexstring("\$\\phi _\\textrm{flip}\$"))
ax.set_ylabel("log( P )")
ax.set_yscale("linear")
ax.set_xlim(0,1)
ax.set_ylim(-2.8,1)
ax.legend(loc="upper left", frameon=true, borderpad=1, handlelength=1.8)
fig.set_size_inches(10scale, 9scale)
fig.tight_layout()
fig.savefig(outfile, bbox_inches="tight", pad_inches=.1)
PyPlot.matplotlib[:rc]("text", usetex=false)
close("all")
#

# z_welker, 03, DM
outnm   = joinpath(outputdir,"welker2014_Mdm_zwelker_03")
outfile = outnm*".pdf"
scale=0.5
N_bins=5
pltcm       = pyimport("matplotlib.cm")
pltcolors   = pyimport("matplotlib.colors")
colormap    = pltcm.get_cmap(name="gnuplot")
fig, ax = subplots()
n_1, x_1        = PyPlot.hist( cosd.(pyplottable(ad03["ϕ_flip"][idxmatch(idxmatch(csw03[:mergerDM_acc_eq0] , csw03[:transition]), csw03[:z_welker])])), range=(0,1), bins=N_bins, density=true, log=false)
n_2, x_2        = PyPlot.hist( cosd.(pyplottable(ad03["ϕ_flip"][idxmatch(idxmatch(csw03[:mergerDM_acc_g0_l5 ], csw03[:transition]), csw03[:z_welker])])), range=(0,1), bins=N_bins, density=true, log=false)
n_3, x_3        = PyPlot.hist( cosd.(pyplottable(ad03["ϕ_flip"][idxmatch(idxmatch(csw03[:mergerDM_acc_g5_l10 ], csw03[:transition]), csw03[:z_welker])])), range=(0,1), bins=N_bins, density=true, log=false)
n_4, x_4        = PyPlot.hist( cosd.(pyplottable(ad03["ϕ_flip"][idxmatch(idxmatch(csw03[:mergerDM_acc_g10], csw03[:transition]), csw03[:z_welker])])), range=(0,1), bins=N_bins, density=true, log=false)
bcenters_1     = 0.5 .* (x_1[2:end]  .+ x_1[1:end-1])
bcenters_2     = 0.5 .* (x_2[2:end]  .+ x_2[1:end-1])
bcenters_3     = 0.5 .* (x_3[2:end]  .+ x_3[1:end-1])
bcenters_4     = 0.5 .* (x_4[2:end]  .+ x_4[1:end-1])
clf()
PyPlot.matplotlib[:rc]("text", usetex=true)
PyPlot.matplotlib[:rc]("font", family="serif", serif=["Computer Modern"])
fig, ax = subplots()
ax.grid(true)

ax.plot(bcenters_1 , log10.(n_1), ".-", ms=10, label="$(latexstring("\$\\delta m = 0\$\\%")), N = $(length(idxmatch(idxmatch(csw03[:mergerDM_acc_eq0] , csw03[:transition]), csw03[:z_welker])))", color="red")
ax.plot(bcenters_2 , log10.(n_2), ".-", ms=10, label="$(latexstring("\$0 < \\delta m < 5\$\\%")), N = $(length(idxmatch(idxmatch(csw03[:mergerDM_acc_g0_l5 ], csw03[:transition]), csw03[:z_welker])))", color="orange"   )
ax.plot(bcenters_3 , log10.(n_3), ".-", ms=10, label="$(latexstring("\$5 < \\delta m < 10\$\\%")), N = $(length(idxmatch(idxmatch(csw03[:mergerDM_acc_g5_l10 ], csw03[:transition]), csw03[:z_welker])))", color="purple"   )
ax.plot(bcenters_4 , log10.(n_4), ".-", ms=10, label="$(latexstring("\$\\delta m > 10\$\\%")), N = $(length(idxmatch(idxmatch(csw03[:mergerDM_acc_g10], csw03[:transition]), csw03[:z_welker])))", color="mediumblue"   )

ax.set_xlabel("cos($(latexstring("\$\\phi\$")))")# $(latexstring("\$\\phi _\\textrm{flip}\$"))
ax.set_ylabel("log( P )")
ax.set_yscale("linear")
ax.set_xlim(0,1)
ax.set_ylim(-2.8,1)
ax.legend(loc="upper left", frameon=true, borderpad=1, handlelength=1.8)
fig.set_size_inches(10scale, 9scale)
fig.tight_layout()
fig.savefig(outfile, bbox_inches="tight", pad_inches=.1)
PyPlot.matplotlib[:rc]("text", usetex=false)
close("all")
#



# z_all, 03, stars
outnm   = joinpath(outputdir,"welker2014_Mstars_zall_03")
outfile = outnm*".pdf"
scale=0.5
N_bins=5
pltcm       = pyimport("matplotlib.cm")
pltcolors   = pyimport("matplotlib.colors")
colormap    = pltcm.get_cmap(name="gnuplot")
fig, ax = subplots()
n_1, x_1        = PyPlot.hist( cosd.(pyplottable(as03["ϕ_flip"][idxmatch(csw03[:mergerSTARS_acc_eq0] , csw03[:transition])])), range=(0,1), bins=N_bins, density=true, log=false)
n_2, x_2        = PyPlot.hist( cosd.(pyplottable(as03["ϕ_flip"][idxmatch(csw03[:mergerSTARS_acc_g0_l5 ], csw03[:transition])])), range=(0,1), bins=N_bins, density=true, log=false)
n_3, x_3        = PyPlot.hist( cosd.(pyplottable(as03["ϕ_flip"][idxmatch(csw03[:mergerSTARS_acc_g5_l10 ], csw03[:transition])])), range=(0,1), bins=N_bins, density=true, log=false)
n_4, x_4        = PyPlot.hist( cosd.(pyplottable(as03["ϕ_flip"][idxmatch(csw03[:mergerSTARS_acc_g10], csw03[:transition])])), range=(0,1), bins=N_bins, density=true, log=false)
bcenters_1     = 0.5 .* (x_1[2:end]  .+ x_1[1:end-1])
bcenters_2     = 0.5 .* (x_2[2:end]  .+ x_2[1:end-1])
bcenters_3     = 0.5 .* (x_3[2:end]  .+ x_3[1:end-1])
bcenters_4     = 0.5 .* (x_4[2:end]  .+ x_4[1:end-1])
clf()
PyPlot.matplotlib[:rc]("text", usetex=true)
PyPlot.matplotlib[:rc]("font", family="serif", serif=["Computer Modern"])
fig, ax = subplots()
ax.grid(true)

ax.plot(bcenters_1 , log10.(n_1), ".-", ms=10, label="$(latexstring("\$\\delta m = 0\$\\%")), N = $(length(idxmatch(csw03[:mergerSTARS_acc_eq0] , csw03[:transition])))", color="red")
ax.plot(bcenters_2 , log10.(n_2), ".-", ms=10, label="$(latexstring("\$0 < \\delta m < 5\$\\%")), N = $(length(idxmatch(csw03[:mergerSTARS_acc_g0_l5 ], csw03[:transition])))", color="orange"   )
ax.plot(bcenters_3 , log10.(n_3), ".-", ms=10, label="$(latexstring("\$5 < \\delta m < 10\$\\%")), N = $(length(idxmatch(csw03[:mergerSTARS_acc_g5_l10 ], csw03[:transition])))", color="purple"   )
ax.plot(bcenters_4 , log10.(n_4), ".-", ms=10, label="$(latexstring("\$\\delta m > 10\$\\%")), N = $(length(idxmatch(csw03[:mergerSTARS_acc_g10], csw03[:transition])))", color="mediumblue"   )

ax.set_xlabel("cos($(latexstring("\$\\phi\$")))")# $(latexstring("\$\\phi _\\textrm{flip}\$"))
ax.set_ylabel("log( P )")
ax.set_yscale("linear")
ax.set_xlim(0,1)
ax.set_ylim(-2.8,1)
ax.legend(loc="upper left", frameon=true, borderpad=1, handlelength=1.8)
fig.set_size_inches(10scale, 9scale)
fig.tight_layout()
fig.savefig(outfile, bbox_inches="tight", pad_inches=.1)
PyPlot.matplotlib[:rc]("text", usetex=false)
close("all")
#

# z_all, 03, DM
outnm   = joinpath(outputdir,"welker2014_Mdm_zall_03")
outfile = outnm*".pdf"
scale=0.5
N_bins=5
pltcm       = pyimport("matplotlib.cm")
pltcolors   = pyimport("matplotlib.colors")
colormap    = pltcm.get_cmap(name="gnuplot")
fig, ax = subplots()
n_1, x_1        = PyPlot.hist( cosd.(pyplottable(ad03["ϕ_flip"][idxmatch(csw03[:mergerDM_acc_eq0] , csw03[:transition])])), range=(0,1), bins=N_bins, density=true, log=false)
n_2, x_2        = PyPlot.hist( cosd.(pyplottable(ad03["ϕ_flip"][idxmatch(csw03[:mergerDM_acc_g0_l5 ], csw03[:transition])])), range=(0,1), bins=N_bins, density=true, log=false)
n_3, x_3        = PyPlot.hist( cosd.(pyplottable(ad03["ϕ_flip"][idxmatch(csw03[:mergerDM_acc_g5_l10 ], csw03[:transition])])), range=(0,1), bins=N_bins, density=true, log=false)
n_4, x_4        = PyPlot.hist( cosd.(pyplottable(ad03["ϕ_flip"][idxmatch(csw03[:mergerDM_acc_g10], csw03[:transition])])), range=(0,1), bins=N_bins, density=true, log=false)
bcenters_1     = 0.5 .* (x_1[2:end]  .+ x_1[1:end-1])
bcenters_2     = 0.5 .* (x_2[2:end]  .+ x_2[1:end-1])
bcenters_3     = 0.5 .* (x_3[2:end]  .+ x_3[1:end-1])
bcenters_4     = 0.5 .* (x_4[2:end]  .+ x_4[1:end-1])
clf()
PyPlot.matplotlib[:rc]("text", usetex=true)
PyPlot.matplotlib[:rc]("font", family="serif", serif=["Computer Modern"])
fig, ax = subplots()
ax.grid(true)

ax.plot(bcenters_1 , log10.(n_1), ".-", ms=10, label="$(latexstring("\$\\delta m = 0\$\\%")), N = $(length(idxmatch(csw03[:mergerDM_acc_eq0] , csw03[:transition])))", color="red")
ax.plot(bcenters_2 , log10.(n_2), ".-", ms=10, label="$(latexstring("\$0 < \\delta m < 5\$\\%")), N = $(length(idxmatch(csw03[:mergerDM_acc_g0_l5 ], csw03[:transition])))", color="orange"   )
ax.plot(bcenters_3 , log10.(n_3), ".-", ms=10, label="$(latexstring("\$5 < \\delta m < 10\$\\%")), N = $(length(idxmatch(csw03[:mergerDM_acc_g5_l10 ], csw03[:transition])))", color="purple"   )
ax.plot(bcenters_4 , log10.(n_4), ".-", ms=10, label="$(latexstring("\$\\delta m > 10\$\\%")), N = $(length(idxmatch(csw03[:mergerDM_acc_g10], csw03[:transition])))", color="mediumblue"   )

ax.set_xlabel("cos($(latexstring("\$\\phi\$")))")# $(latexstring("\$\\phi _\\textrm{flip}\$"))
ax.set_ylabel("log( P )")
ax.set_yscale("linear")
ax.set_xlim(0,1)
ax.set_ylim(-2.8,1)
ax.legend(loc="upper left", frameon=true, borderpad=1, handlelength=1.8)
fig.set_size_inches(10scale, 9scale)
fig.tight_layout()
fig.savefig(outfile, bbox_inches="tight", pad_inches=.1)
PyPlot.matplotlib[:rc]("text", usetex=false)
close("all")
#






#####################################################################################
# sSFR vs delta J bval split
outnm   = joinpath(outputdir,"sSFR_vs_delJup_bpostsplit_03")
outfile = outnm*"_hxcount.pdf"
x   = pyplottable(log10.(pyplottable(as03["ΔJ_main"])))
y   = pyplottable(log10.(as03["SFR"] .*1e9 ./ as03["M"]))
miny= minimum(pyplottable(y)[findall(x->x .!== NaN, pyplottable(y))])
maxy= maximum(pyplottable(y)[findall(x->x .!== NaN, pyplottable(y))])
minx= minimum(pyplottable(x)[findall(x->x .!== NaN, pyplottable(x))])
maxx= maximum(pyplottable(x)[findall(x->x .!== NaN, pyplottable(x))])
slcts       = Dict{Int64, Vector{Int64}}()
slcts[1]    = csw03[:transition]
slcts[2]    = idxmatch(csw03[:transition], csw03[:result_int ])
slcts[3]    = idxmatch(csw03[:transition], csw03[:result_ell])
slcts[4]    = idxmatch(csw03[:transition], csw03[:result_disk ])
labels      = Dict{Int64, String}()
labels[1]   = "Full Set, N = $(length(slcts[1])-count(isnan,y[slcts[1]]))"
labels[2]   = "Intermediates, N = $(length(slcts[2])-count(isnan,y[slcts[2]]))"
labels[3]   = "Ellipticals, N = $(length(slcts[3])-count(isnan,y[slcts[3]]))"
labels[4]   = "Disks, N = $(length(slcts[4])-count(isnan,y[slcts[4]]))"
labpos      = Array{String}(undef, 4)
labpos     .= "lower right"
plot_hexbins_nxm( x, 
            y, 
            slcts, labels, 2, 2,
            #scale=0.0,
            cbarshift=2.5,
            label_pos=labpos, 
            outfile=outfile, 
            lognorm=true, gridres=(15,5), 
            ymin=miny, ymax=maxy,
            xmin=minx, xmax=maxx, 
            ylabel="log( sSFR [1/Gyr] )",
            xlabel="log( $(latexstring("\$\\Delta J_*\$")) [ $(latexstring("\$M_\\odot\$")) kpc km/s ] )", fontsize=30,
            binned_median="y", plot_bval=" ",
            )
#
# sSFR vs delta j bval split
outnm   = joinpath(outputdir,"sSFR_vs_delj_bpostsplit_03")
outfile = outnm*"_hxcount.pdf"
x   = pyplottable(log10.(pyplottable(as03["Δj_main"])))
y   = pyplottable(log10.(as03["SFR"] .*1e9 ./ as03["M"]))
miny= minimum(pyplottable(y)[findall(x->x .!== NaN, pyplottable(y))])
maxy= maximum(pyplottable(y)[findall(x->x .!== NaN, pyplottable(y))])
minx= minimum(pyplottable(x)[findall(x->x .!== NaN, pyplottable(x))])
maxx= maximum(pyplottable(x)[findall(x->x .!== NaN, pyplottable(x))])
slcts       = Dict{Int64, Vector{Int64}}()
slcts[1]    = csw03[:transition]
slcts[2]    = idxmatch(csw03[:transition], csw03[:result_int ])
slcts[3]    = idxmatch(csw03[:transition], csw03[:result_ell])
slcts[4]    = idxmatch(csw03[:transition], csw03[:result_disk ])
labels      = Dict{Int64, String}()
labels[1]   = "Full Set, N = $(length(slcts[1])-count(isnan,y[slcts[1]]))"
labels[2]   = "Intermediates, N = $(length(slcts[2])-count(isnan,y[slcts[2]]))"
labels[3]   = "Ellipticals, N = $(length(slcts[3])-count(isnan,y[slcts[3]]))"
labels[4]   = "Disks, N = $(length(slcts[4])-count(isnan,y[slcts[4]]))"
labpos      = Array{String}(undef, 4)
labpos     .= "lower right"
plot_hexbins_nxm( x, 
            y, 
            slcts, labels, 2, 2,
            #scale=0.0,
            cbarshift=2.5,
            label_pos=labpos, 
            outfile=outfile, 
            lognorm=true, gridres=(15,5), 
            ymin=miny, ymax=maxy,
            xmin=minx, xmax=maxx, 
            ylabel="log( sSFR [1/Gyr] )",
            xlabel="log( $(latexstring("\$\\Delta j_*\$")) [ kpc km/s ] )", fontsize=30,
            binned_median="y", plot_bval=" ",
            )
#
outnm   = joinpath(outputdir,"sSFR_vs_delJup_bpostsplit_09")
outfile = outnm*"_hxcount.pdf"
x   = pyplottable(log10.(pyplottable(as09["ΔJ_main"])))
y   = pyplottable(log10.(as09["SFR"] .*1e9 ./ as09["M"]))
miny= minimum(pyplottable(y)[findall(x->x .!== NaN, pyplottable(y))])
maxy= maximum(pyplottable(y)[findall(x->x .!== NaN, pyplottable(y))])
minx= minimum(pyplottable(x)[findall(x->x .!== NaN, pyplottable(x))])
maxx= maximum(pyplottable(x)[findall(x->x .!== NaN, pyplottable(x))])
slcts       = Dict{Int64, Vector{Int64}}()
slcts[1]    = csw09[:transition]
slcts[2]    = idxmatch(csw09[:transition], csw09[:result_int ])
slcts[3]    = idxmatch(csw09[:transition], csw09[:result_ell])
slcts[4]    = idxmatch(csw09[:transition], csw09[:result_disk ])
labels      = Dict{Int64, String}()
labels[1]   = "Full Set, N = $(length(slcts[1])-count(isnan,y[slcts[1]]))"
labels[2]   = "Intermediates, N = $(length(slcts[2])-count(isnan,y[slcts[2]]))"
labels[3]   = "Ellipticals, N = $(length(slcts[3])-count(isnan,y[slcts[3]]))"
labels[4]   = "Disks, N = $(length(slcts[4])-count(isnan,y[slcts[4]]))"
labpos      = Array{String}(undef, 4)
labpos     .= "lower right"
plot_hexbins_nxm( x, 
            y, 
            slcts, labels, 2, 2,
            #scale=0.0,
            cbarshift=2.5,
            label_pos=labpos, 
            outfile=outfile, 
            lognorm=true, gridres=(15,5), 
            ymin=miny, ymax=maxy,
            xmin=minx, xmax=maxx, 
            ylabel="log( sSFR [1/Gyr] )",
            xlabel="log( $(latexstring("\$\\Delta J_*\$")) [ $(latexstring("\$M_\\odot\$")) kpc km/s ] )", fontsize=30,
            binned_median="y", plot_bval=" ",
            )
#
# sSFR vs delta j bval split
outnm   = joinpath(outputdir,"sSFR_vs_delj_bpostsplit_09")
outfile = outnm*"_hxcount.pdf"
x   = pyplottable(log10.(pyplottable(as09["Δj_main"])))
y   = pyplottable(log10.(as09["SFR"] .*1e9 ./ as09["M"]))
miny= minimum(pyplottable(y)[findall(x->x .!== NaN, pyplottable(y))])
maxy= maximum(pyplottable(y)[findall(x->x .!== NaN, pyplottable(y))])
minx= minimum(pyplottable(x)[findall(x->x .!== NaN, pyplottable(x))])
maxx= maximum(pyplottable(x)[findall(x->x .!== NaN, pyplottable(x))])
slcts       = Dict{Int64, Vector{Int64}}()
slcts[1]    = csw09[:transition]
slcts[2]    = idxmatch(csw09[:transition], csw09[:result_int ])
slcts[3]    = idxmatch(csw09[:transition], csw09[:result_ell])
slcts[4]    = idxmatch(csw09[:transition], csw09[:result_disk ])
labels      = Dict{Int64, String}()
labels[1]   = "Full Set, N = $(length(slcts[1])-count(isnan,y[slcts[1]]))"
labels[2]   = "Intermediates, N = $(length(slcts[2])-count(isnan,y[slcts[2]]))"
labels[3]   = "Ellipticals, N = $(length(slcts[3])-count(isnan,y[slcts[3]]))"
labels[4]   = "Disks, N = $(length(slcts[4])-count(isnan,y[slcts[4]]))"
labpos      = Array{String}(undef, 4)
labpos     .= "lower right"
plot_hexbins_nxm( x, 
            y, 
            slcts, labels, 2, 2,
            #scale=0.0,
            cbarshift=2.5,
            label_pos=labpos, 
            outfile=outfile, 
            lognorm=true, gridres=(15,5), 
            ymin=miny, ymax=maxy,
            xmin=minx, xmax=maxx, 
            ylabel="log( sSFR [1/Gyr] )",
            xlabel="log( $(latexstring("\$\\Delta j_*\$")) [ kpc km/s ] )", fontsize=30,
            binned_median="y", plot_bval=" ",
            )
#
# sSFR vs delta jrel bval split
outnm   = joinpath(outputdir,"sSFR_vs_deljrel_bpostsplit_03")
outfile = outnm*"_hxcount.pdf"
x   = pyplottable(log10.(pyplottable(as03["Δj_main"] ./ (as03["j_main"] .- as03["Δj_main"]))))
y   = pyplottable(log10.(as03["SFR"] .*1e9 ./ as03["M"]))
miny= minimum(pyplottable(y)[findall(x->x .!== NaN, pyplottable(y))])
maxy= maximum(pyplottable(y)[findall(x->x .!== NaN, pyplottable(y))])
minx= minimum(pyplottable(x)[findall(x->x .!== NaN, pyplottable(x))])
maxx= maximum(pyplottable(x)[findall(x->x .!== NaN, pyplottable(x))])
slcts       = Dict{Int64, Vector{Int64}}()
slcts[1]    = csw03[:transition]
slcts[2]    = idxmatch(csw03[:transition], csw03[:result_int ])
slcts[3]    = idxmatch(csw03[:transition], csw03[:result_ell])
slcts[4]    = idxmatch(csw03[:transition], csw03[:result_disk ])
labels      = Dict{Int64, String}()
labels[1]   = "Full Set, N = $(length(slcts[1])-count(isnan,y[slcts[1]]))"
labels[2]   = "Intermediates, N = $(length(slcts[2])-count(isnan,y[slcts[2]]))"
labels[3]   = "Ellipticals, N = $(length(slcts[3])-count(isnan,y[slcts[3]]))"
labels[4]   = "Disks, N = $(length(slcts[4])-count(isnan,y[slcts[4]]))"
labpos      = Array{String}(undef, 4)
labpos     .= "lower right"
plot_hexbins_nxm( x, 
            y, 
            slcts, labels, 2, 2,
            #scale=0.0,
            cbarshift=2.5,
            label_pos=labpos, 
            outfile=outfile, 
            lognorm=true, gridres=(15,5), 
            ymin=miny, ymax=maxy,
            xmin=minx, xmax=maxx, 
            ylabel="log( sSFR [1/Gyr] )",
            xlabel="log( $(latexstring("\$\\Delta j_\\textrm{rel}\$"))", fontsize=30,
            binned_median="y", plot_bval=" ",
            )
#

insane_y = findcs(pyplottable(log10.(as03["SFR"] .*1e9 ./ as03["M"])), eq=NaN)
insane_x = findcs(pyplottable(log10.(pyplottable(as03["Δj_main"] ./ (as03["j_main"] .- as03["Δj_main"])))), eq=NaN)
sane     = idxclude(idxclude(1:length(as03["M"]), insane_x), insane_y)
println("end as ELL: ", length(idxmatch(idxmatch(sane, idxmatch(csw03[:transition], csw03[:result_ell  ])), csw03[:start_ell ])))
println("end as ELL: ", length(idxmatch(idxmatch(sane, idxmatch(csw03[:transition], csw03[:result_ell  ])), csw03[:start_int ])))
println("end as ELL: ", length(idxmatch(idxmatch(sane, idxmatch(csw03[:transition], csw03[:result_ell  ])), csw03[:start_disk])))
println("end as INT: ", length(idxmatch(idxmatch(sane, idxmatch(csw03[:transition], csw03[:result_int  ])), csw03[:start_ell ])))
println("end as INT: ", length(idxmatch(idxmatch(sane, idxmatch(csw03[:transition], csw03[:result_int  ])), csw03[:start_int ])))
println("end as INT: ", length(idxmatch(idxmatch(sane, idxmatch(csw03[:transition], csw03[:result_int  ])), csw03[:start_disk])))
println("end as DSK: ", length(idxmatch(idxmatch(sane, idxmatch(csw03[:transition], csw03[:result_disk ])), csw03[:start_ell ])))
println("end as DSK: ", length(idxmatch(idxmatch(sane, idxmatch(csw03[:transition], csw03[:result_disk ])), csw03[:start_int ])))
println("end as DSK: ", length(idxmatch(idxmatch(sane, idxmatch(csw03[:transition], csw03[:result_disk ])), csw03[:start_disk])))


insane_y = findcs(pyplottable(log10.(as09["SFR"] .*1e9 ./ as09["M"])), eq=NaN)
insane_x = findcs(pyplottable(log10.(pyplottable(as09["Δj_main"] ./ (as09["j_main"] .- as09["Δj_main"])))), eq=NaN)
sane     = idxclude(idxclude(1:length(as09["M"]), insane_x), insane_y)
println("end as ELL: ", length(idxmatch(idxmatch(sane, idxmatch(csw09[:transition], csw09[:result_ell  ])), csw09[:start_ell ])))
println("end as ELL: ", length(idxmatch(idxmatch(sane, idxmatch(csw09[:transition], csw09[:result_ell  ])), csw09[:start_int ])))
println("end as ELL: ", length(idxmatch(idxmatch(sane, idxmatch(csw09[:transition], csw09[:result_ell  ])), csw09[:start_disk])))
println("end as INT: ", length(idxmatch(idxmatch(sane, idxmatch(csw09[:transition], csw09[:result_int  ])), csw09[:start_ell ])))
println("end as INT: ", length(idxmatch(idxmatch(sane, idxmatch(csw09[:transition], csw09[:result_int  ])), csw09[:start_int ])))
println("end as INT: ", length(idxmatch(idxmatch(sane, idxmatch(csw09[:transition], csw09[:result_int  ])), csw09[:start_disk])))
println("end as DSK: ", length(idxmatch(idxmatch(sane, idxmatch(csw09[:transition], csw09[:result_disk ])), csw09[:start_ell ])))
println("end as DSK: ", length(idxmatch(idxmatch(sane, idxmatch(csw09[:transition], csw09[:result_disk ])), csw09[:start_int ])))
println("end as DSK: ", length(idxmatch(idxmatch(sane, idxmatch(csw09[:transition], csw09[:result_disk ])), csw09[:start_disk])))


# sSFR vs delta jrel bval split
outnm   = joinpath(outputdir,"sSFR_vs_deljrel_bpostsplit_09")
outfile = outnm*"_hxcount.pdf"
x   = pyplottable(log10.(pyplottable(as09["Δj_main"] ./ (as09["j_main"] .- as09["Δj_main"]))))
y   = pyplottable(log10.(as09["SFR"] .*1e9 ./ as09["M"]))
miny= minimum(pyplottable(y)[findall(x->x .!== NaN, pyplottable(y))])
maxy= maximum(pyplottable(y)[findall(x->x .!== NaN, pyplottable(y))])
minx= minimum(pyplottable(x)[findall(x->x .!== NaN, pyplottable(x))])
maxx= maximum(pyplottable(x)[findall(x->x .!== NaN, pyplottable(x))])
slcts       = Dict{Int64, Vector{Int64}}()
slcts[1]    = csw09[:transition]
slcts[2]    = idxmatch(csw09[:transition], csw09[:result_int ])
slcts[3]    = idxmatch(csw09[:transition], csw09[:result_ell])
slcts[4]    = idxmatch(csw09[:transition], csw09[:result_disk ])
labels      = Dict{Int64, String}()
labels[1]   = "Full Set, N = $(length(slcts[1])-count(isnan,y[slcts[1]]))"
labels[2]   = "Intermediates, N = $(length(slcts[2])-count(isnan,y[slcts[2]]))"
labels[3]   = "Ellipticals, N = $(length(slcts[3])-count(isnan,y[slcts[3]]))"
labels[4]   = "Disks, N = $(length(slcts[4])-count(isnan,y[slcts[4]]))"
labpos      = Array{String}(undef, 4)
labpos     .= "lower right"
plot_hexbins_nxm( x, 
            y, 
            slcts, labels, 2, 2,
            #scale=0.0,
            cbarshift=2.5,
            label_pos=labpos, 
            outfile=outfile, 
            lognorm=true, gridres=(15,5), 
            ymin=miny, ymax=maxy,
            xmin=minx, xmax=maxx, 
            ylabel="log( sSFR [1/Gyr] )",
            xlabel="log( $(latexstring("\$\\Delta j_\\textrm{rel}\$"))", fontsize=30,
            binned_median="y", plot_bval=" ",
            )
#


# z-low
outnm   = joinpath(outputdir,"sSFR_vs_delJup_bpostsplit_zlow_03")
outfile = outnm*"_hxcount.pdf"
x   = pyplottable(log10.(pyplottable(as03["ΔJ_main"])))
y   = pyplottable(log10.(pyplottable(as03["SFR"] .*1e9 ./ as03["M"])))
miny= minimum(pyplottable(y)[findall(x->x .!== NaN, pyplottable(y))])
maxy= maximum(pyplottable(y)[findall(x->x .!== NaN, pyplottable(y))])
minx= minimum(pyplottable(x)[findall(x->x .!== NaN, pyplottable(x))])
maxx= maximum(pyplottable(x)[findall(x->x .!== NaN, pyplottable(x))])
slcts       = Dict{Int64, Vector{Int64}}()
slcts[1]    = idxmatch(csw03[:transition], csw03[:z_lo]) # z_split15_l
slcts[2]    = idxmatch(idxmatch(csw03[:transition], csw03[:result_int ]) , csw03[:z_lo])
slcts[3]    = idxmatch(idxmatch(csw03[:transition], csw03[:result_ell])  , csw03[:z_lo])
slcts[4]    = idxmatch(idxmatch(csw03[:transition], csw03[:result_disk ]), csw03[:z_lo])
labels      = Dict{Int64, String}()
labels[1]   = "Full Set, N = $(length(slcts[1])-count(isnan,y[slcts[1]]))"
labels[2]   = "Intermediates, N = $(length(slcts[2])-count(isnan,y[slcts[2]]))"
labels[3]   = "Ellipticals, N = $(length(slcts[3])-count(isnan,y[slcts[3]]))"
labels[4]   = "Disks, N = $(length(slcts[4])-count(isnan,y[slcts[4]]))"
labpos      = Array{String}(undef, 4)
labpos     .= "lower right"
plot_hexbins_nxm( x, 
            y, 
            slcts, labels, 2, 2,
            #scale=0.0,
            cbarshift=2.5,
            label_pos=labpos, 
            outfile=outfile, 
            lognorm=true, gridres=(15,5), 
            ymin=miny, ymax=maxy,
            xmin=minx, xmax=maxx, 
            ylabel="log( sSFR [1/Gyr] )",
            xlabel="log( $(latexstring("\$\\Delta J_*\$")) [$(latexstring("\$M_\\odot\$")) kpc km/s ] )", fontsize=30,
            binned_median="y", plot_bval=" ",
            )
#




# z-high
outnm   = joinpath(outputdir,"sSFR_vs_delJup_bpostsplit_zhigh_03")
outfile = outnm*"_hxcount.pdf"
x   = pyplottable(log10.(pyplottable(as03["ΔJ_main"])))
y   = pyplottable(log10.(pyplottable(as03["SFR"] .*1e9 ./ as03["M"])))
miny= minimum(pyplottable(y)[findall(x->x .!== NaN, pyplottable(y))])
maxy= maximum(pyplottable(y)[findall(x->x .!== NaN, pyplottable(y))])
minx= minimum(pyplottable(x)[findall(x->x .!== NaN, pyplottable(x))])
maxx= maximum(pyplottable(x)[findall(x->x .!== NaN, pyplottable(x))])
slcts       = Dict{Int64, Vector{Int64}}()
slcts[1]    = idxmatch(csw03[:transition], csw03[:z_hi]) # z_split15_g
slcts[2]    = idxmatch(idxmatch(csw03[:transition], csw03[:result_int ]) , csw03[:z_hi])
slcts[3]    = idxmatch(idxmatch(csw03[:transition], csw03[:result_ell])  , csw03[:z_hi])
slcts[4]    = idxmatch(idxmatch(csw03[:transition], csw03[:result_disk ]), csw03[:z_hi])
labels      = Dict{Int64, String}()
labels[1]   = "Full Set, N = $(length(slcts[1])-count(isnan,y[slcts[1]]))"
labels[2]   = "Intermediates, N = $(length(slcts[2])-count(isnan,y[slcts[2]]))"
labels[3]   = "Ellipticals, N = $(length(slcts[3])-count(isnan,y[slcts[3]]))"
labels[4]   = "Disks, N = $(length(slcts[4])-count(isnan,y[slcts[4]]))"
labpos      = Array{String}(undef, 4)
labpos     .= "lower right"
plot_hexbins_nxm( x, 
            y, 
            slcts, labels, 2, 2,
            #scale=0.0,
            cbarshift=2.5,
            label_pos=labpos, 
            outfile=outfile, 
            lognorm=true, gridres=(15,5), 
            ymin=miny, ymax=maxy,
            xmin=minx, xmax=maxx, 
            ylabel="log( sSFR [1/Gyr] )",
            xlabel="log( $(latexstring("\$\\Delta J_*\$")) [$(latexstring("\$M_\\odot\$")) kpc km/s ] )", fontsize=30,
            binned_median="y", plot_bval=" ",
            )
#
# z-low
outnm   = joinpath(outputdir,"sSFR_vs_delj_bpostsplit_zlow_03")
outfile = outnm*"_hxcount.pdf"
x   = pyplottable(log10.(pyplottable(as03["Δj_main"])))
y   = pyplottable(log10.(pyplottable(as03["SFR"] .*1e9 ./ as03["M"])))
miny= minimum(pyplottable(y)[findall(x->x .!== NaN, pyplottable(y))])
maxy= maximum(pyplottable(y)[findall(x->x .!== NaN, pyplottable(y))])
minx= minimum(pyplottable(x)[findall(x->x .!== NaN, pyplottable(x))])
maxx= maximum(pyplottable(x)[findall(x->x .!== NaN, pyplottable(x))])
slcts       = Dict{Int64, Vector{Int64}}()
slcts[1]    = idxmatch(csw03[:transition], csw03[:z_lo]) # z_split15_l
slcts[2]    = idxmatch(idxmatch(csw03[:transition], csw03[:result_int ]) , csw03[:z_lo])
slcts[3]    = idxmatch(idxmatch(csw03[:transition], csw03[:result_ell])  , csw03[:z_lo])
slcts[4]    = idxmatch(idxmatch(csw03[:transition], csw03[:result_disk ]), csw03[:z_lo])
labels      = Dict{Int64, String}()
labels[1]   = "Full Set, N = $(length(slcts[1])-count(isnan,y[slcts[1]]))"
labels[2]   = "Intermediates, N = $(length(slcts[2])-count(isnan,y[slcts[2]]))"
labels[3]   = "Ellipticals, N = $(length(slcts[3])-count(isnan,y[slcts[3]]))"
labels[4]   = "Disks, N = $(length(slcts[4])-count(isnan,y[slcts[4]]))"
labpos      = Array{String}(undef, 4)
labpos     .= "lower right"
plot_hexbins_nxm( x, 
            y, 
            slcts, labels, 2, 2,
            #scale=0.0,
            cbarshift=2.5,
            label_pos=labpos, 
            outfile=outfile, 
            lognorm=true, gridres=(15,5), 
            ymin=miny, ymax=maxy,
            xmin=minx, xmax=maxx, 
            ylabel="log( sSFR [1/Gyr] )",
            xlabel="log( $(latexstring("\$\\Delta j_*\$")) [ kpc km/s ] )", fontsize=30,
            binned_median="y", plot_bval=" ",
            )
#




# z-high
outnm   = joinpath(outputdir,"sSFR_vs_delj_bpostsplit_zhigh_03")
outfile = outnm*"_hxcount.pdf"
x   = pyplottable(log10.(pyplottable(as03["Δj_main"])))
y   = pyplottable(log10.(pyplottable(as03["SFR"] .*1e9 ./ as03["M"])))
miny= minimum(pyplottable(y)[findall(x->x .!== NaN, pyplottable(y))])
maxy= maximum(pyplottable(y)[findall(x->x .!== NaN, pyplottable(y))])
minx= minimum(pyplottable(x)[findall(x->x .!== NaN, pyplottable(x))])
maxx= maximum(pyplottable(x)[findall(x->x .!== NaN, pyplottable(x))])
slcts       = Dict{Int64, Vector{Int64}}()
slcts[1]    = idxmatch(csw03[:transition], csw03[:z_hi]) # z_split15_g
slcts[2]    = idxmatch(idxmatch(csw03[:transition], csw03[:result_int ]) , csw03[:z_hi])
slcts[3]    = idxmatch(idxmatch(csw03[:transition], csw03[:result_ell])  , csw03[:z_hi])
slcts[4]    = idxmatch(idxmatch(csw03[:transition], csw03[:result_disk ]), csw03[:z_hi])
labels      = Dict{Int64, String}()
labels[1]   = "Full Set, N = $(length(slcts[1])-count(isnan,y[slcts[1]]))"
labels[2]   = "Intermediates, N = $(length(slcts[2])-count(isnan,y[slcts[2]]))"
labels[3]   = "Ellipticals, N = $(length(slcts[3])-count(isnan,y[slcts[3]]))"
labels[4]   = "Disks, N = $(length(slcts[4])-count(isnan,y[slcts[4]]))"
labpos      = Array{String}(undef, 4)
labpos     .= "lower right"
plot_hexbins_nxm( x, 
            y, 
            slcts, labels, 2, 2,
            #scale=0.0,
            cbarshift=2.5,
            label_pos=labpos, 
            outfile=outfile, 
            lognorm=true, gridres=(15,5), 
            ymin=miny, ymax=maxy,
            xmin=minx, xmax=maxx, 
            ylabel="log( sSFR [1/Gyr] )",
            xlabel="log( $(latexstring("\$\\Delta j_*\$")) [ kpc km/s ] )", fontsize=30,
            binned_median="y", plot_bval=" ",
            )
#










#####################################################################################
# sSFR vs flipangle bval split
outnm   = joinpath(outputdir,"sSFR_vs_flip_bpostsplit_03")
outfile = outnm*"_hxcount.pdf"
x   = pyplottable(as03["ϕ_flip"])
y   = pyplottable(log10.(pyplottable(as03["SFR"] .*1e9 ./ as03["M"])))
miny= minimum(pyplottable(y)[findall(x->x .!== NaN, pyplottable(y))])
maxy= maximum(pyplottable(y)[findall(x->x .!== NaN, pyplottable(y))])
minx= 0#minimum(pyplottable(x)[findall(x->x .!== NaN, pyplottable(x))])
maxx= 180#maximum(pyplottable(x)[findall(x->x .!== NaN, pyplottable(x))])
slcts       = Dict{Int64, Vector{Int64}}()
slcts[1]    = csw03[:transition]
slcts[2]    = idxmatch(csw03[:transition], csw03[:result_int ])
slcts[3]    = idxmatch(csw03[:transition], csw03[:result_ell])
slcts[4]    = idxmatch(csw03[:transition], csw03[:result_disk ])
labels      = Dict{Int64, String}()
labels[1]   = "Full Set, N = $(length(slcts[1])-count(isnan,y[slcts[1]]))"
labels[2]   = "Intermediates, N = $(length(slcts[2])-count(isnan,y[slcts[2]]))"
labels[3]   = "Ellipticals, N = $(length(slcts[3])-count(isnan,y[slcts[3]]))"
labels[4]   = "Disks, N = $(length(slcts[4])-count(isnan,y[slcts[4]]))"
labpos      = Array{String}(undef, 4)
labpos     .= "lower right"
plot_hexbins_nxm( x, 
            y, 
            slcts, labels, 2, 2,
            #scale=0.0,
            cbarshift=2.5,
            label_pos=labpos, 
            outfile=outfile, 
            grid=true,
            lognorm=true, gridres=(15,5), 
            ymin=miny, ymax=maxy,
            xmin=minx, xmax=maxx, 
            ylabel="log( sSFR [1/Gyr] )",
            xlabel="Flip Angle [°]", fontsize=30,
            binned_median="y", plot_bval=" ",
            )
#


# z-low
outnm   = joinpath(outputdir,"sSFR_vs_flip_bpostsplit_zlow_03")
outfile = outnm*"_hxcount.pdf"
slcts       = Dict{Int64, Vector{Int64}}()
slcts[1]    = idxmatch(csw03[:transition], csw03[:z_lo]) # z_split15_l
slcts[2]    = idxmatch(idxmatch(csw03[:transition], csw03[:result_int ]) , csw03[:z_lo])
slcts[3]    = idxmatch(idxmatch(csw03[:transition], csw03[:result_ell])  , csw03[:z_lo])
slcts[4]    = idxmatch(idxmatch(csw03[:transition], csw03[:result_disk ]), csw03[:z_lo])
labels      = Dict{Int64, String}()
labels[1]   = "Full Set, N = $(length(slcts[1])-count(isnan,y[slcts[1]]))"
labels[2]   = "Intermediates, N = $(length(slcts[2])-count(isnan,y[slcts[2]]))"
labels[3]   = "Ellipticals, N = $(length(slcts[3])-count(isnan,y[slcts[3]]))"
labels[4]   = "Disks, N = $(length(slcts[4])-count(isnan,y[slcts[4]]))"
labpos      = Array{String}(undef, 4)
labpos     .= "lower right"
plot_hexbins_nxm( x, 
            y, 
            slcts, labels, 2, 2,
            #scale=0.0,
            cbarshift=2.5,
            label_pos=labpos, 
            outfile=outfile, 
            grid=true,
            lognorm=true, gridres=(15,5), 
            ymin=miny, ymax=maxy,
            xmin=minx, xmax=maxx, 
            ylabel="log( sSFR [1/Gyr] )",
            xlabel="Flip Angle [°]", fontsize=30,
            binned_median="y", plot_bval=" ",
            )
#




# z-low
outnm   = joinpath(outputdir,"sSFR_vs_flip_bpostsplit_zhigh_03")
outfile = outnm*"_hxcount.pdf"
slcts       = Dict{Int64, Vector{Int64}}()
slcts[1]    = idxmatch(csw03[:transition], csw03[:z_hi]) # z_split15_g
slcts[2]    = idxmatch(idxmatch(csw03[:transition], csw03[:result_int ]) , csw03[:z_hi])
slcts[3]    = idxmatch(idxmatch(csw03[:transition], csw03[:result_ell])  , csw03[:z_hi])
slcts[4]    = idxmatch(idxmatch(csw03[:transition], csw03[:result_disk ]), csw03[:z_hi])
labels      = Dict{Int64, String}()
labels[1]   = "Full Set, N = $(length(slcts[1])-count(isnan,y[slcts[1]]))"
labels[2]   = "Intermediates, N = $(length(slcts[2])-count(isnan,y[slcts[2]]))"
labels[3]   = "Ellipticals, N = $(length(slcts[3])-count(isnan,y[slcts[3]]))"
labels[4]   = "Disks, N = $(length(slcts[4])-count(isnan,y[slcts[4]]))"
labpos      = Array{String}(undef, 4)
labpos     .= "lower right"
plot_hexbins_nxm( x, 
            y, 
            slcts, labels, 2, 2,
            #scale=0.0,
            cbarshift=2.5,
            label_pos=labpos, 
            outfile=outfile, 
            grid=true,
            lognorm=true, gridres=(15,5), 
            ymin=miny, ymax=maxy,
            xmin=minx, xmax=maxx, 
            ylabel="log( sSFR [1/Gyr] )",
            xlabel="Flip Angle [°]", fontsize=30,
            binned_median="y", plot_bval=" ",
            )
#









#####################################################################################
# sSFR vs flipangle bval z split
# bpost
outnm   = joinpath(outputdir,"sSFR_vs_flip_zbpostsplit_03")
outfile = outnm*"_hxcount.pdf"
x   = pyplottable(as03["ϕ_flip"])
y   = pyplottable(log10.(pyplottable(as03["SFR"] .*1e9 ./ as03["M"])))
miny= minimum(pyplottable(y)[findall(x->x .!== NaN, pyplottable(y))])
maxy= maximum(pyplottable(y)[findall(x->x .!== NaN, pyplottable(y))])
minx= 0#minimum(pyplottable(x)[findall(x->x .!== NaN, pyplottable(x))])
maxx= 180#maximum(pyplottable(x)[findall(x->x .!== NaN, pyplottable(x))])
slcts       = Dict{Int64, Vector{Int64}}()
slcts[1]    = idxmatch(idxmatch(csw03[:transition], csw03[:result_ell ])  , csw03[:z_lo])
slcts[2]    = idxmatch(idxmatch(csw03[:transition], csw03[:result_int ])  , csw03[:z_lo])
slcts[3]    = idxmatch(idxmatch(csw03[:transition], csw03[:result_disk])  , csw03[:z_lo])
slcts[4]    = idxmatch(idxmatch(csw03[:transition], csw03[:result_ell ])  , csw03[:z_md])
slcts[5]    = idxmatch(idxmatch(csw03[:transition], csw03[:result_int ])  , csw03[:z_md])
slcts[6]    = idxmatch(idxmatch(csw03[:transition], csw03[:result_disk])  , csw03[:z_md])
slcts[7]    = idxmatch(idxmatch(csw03[:transition], csw03[:result_ell ])  , csw03[:z_hi])
slcts[8]    = idxmatch(idxmatch(csw03[:transition], csw03[:result_int ])  , csw03[:z_hi])
slcts[9]    = idxmatch(idxmatch(csw03[:transition], csw03[:result_disk])  , csw03[:z_hi])
labels      = Dict{Int64, String}()
labels[1]   = "Ell., $(latexstring("z < 0.45")), N = $(length(slcts[1])-count(isnan,y[slcts[1]]))"
labels[2]   = "Int., $(latexstring("z < 0.45")), N = $(length(slcts[2])-count(isnan,y[slcts[2]]))"
labels[3]   = "Disks, $(latexstring("z < 0.45")), N = $(length(slcts[3])-count(isnan,y[slcts[3]]))"
labels[4]   = "Ell., $(latexstring("0.45 < z < 1")), N = $(length(slcts[4])-count(isnan,y[slcts[4]]))"
labels[5]   = "Int., $(latexstring("0.45 < z < 1")), N = $(length(slcts[5])-count(isnan,y[slcts[5]]))"
labels[6]   = "Disks, $(latexstring("0.45 < z < 1")), N = $(length(slcts[6])-count(isnan,y[slcts[6]]))"
labels[7]   = "Ell., $(latexstring("z > 1")), N = $(length(slcts[7])-count(isnan,y[slcts[7]]))"
labels[8]   = "Int., $(latexstring("z > 1")), N = $(length(slcts[8])-count(isnan,y[slcts[8]]))"
labels[9]   = "Disks, $(latexstring("z > 1")), N = $(length(slcts[9])-count(isnan,y[slcts[9]]))"
labpos      = Array{String}(undef, 9)
labpos     .= "lower right"
plot_hexbins_nxm( x, 
            y, 
            slcts, labels, 3, 3,
            #scale=0.0,
            cbarshift=8,
            label_pos=labpos, 
            outfile=outfile, 
            grid=true,
            lognorm=true, gridres=(15,5), 
            ymin=miny, ymax=maxy,
            xmin=minx, xmax=maxx, 
            ylabel="log( sSFR [1/Gyr] )",
            xlabel="Flip Angle [°]", fontsize=30,
            binned_median="y", plot_bval=" ",
            )
#


# bpre
outnm   = joinpath(outputdir,"sSFR_vs_flip_zbpresplit_03")
outfile = outnm*"_hxcount.pdf"
x   = pyplottable(as03["ϕ_flip"])
y   = pyplottable(log10.(pyplottable(as03["SFR"] .*1e9 ./ as03["M"])))
miny= minimum(pyplottable(y)[findall(x->x .!== NaN, pyplottable(y))])
maxy= maximum(pyplottable(y)[findall(x->x .!== NaN, pyplottable(y))])
minx= 0#minimum(pyplottable(x)[findall(x->x .!== NaN, pyplottable(x))])
maxx= 180#maximum(pyplottable(x)[findall(x->x .!== NaN, pyplottable(x))])
slcts       = Dict{Int64, Vector{Int64}}()
slcts[1]    = idxmatch(idxmatch(csw03[:transition], csw03[:start_ell ])  , csw03[:z_lo])
slcts[2]    = idxmatch(idxmatch(csw03[:transition], csw03[:start_int ])  , csw03[:z_lo])
slcts[3]    = idxmatch(idxmatch(csw03[:transition], csw03[:start_disk])  , csw03[:z_lo])
slcts[4]    = idxmatch(idxmatch(csw03[:transition], csw03[:start_ell ])  , csw03[:z_md])
slcts[5]    = idxmatch(idxmatch(csw03[:transition], csw03[:start_int ])  , csw03[:z_md])
slcts[6]    = idxmatch(idxmatch(csw03[:transition], csw03[:start_disk])  , csw03[:z_md])
slcts[7]    = idxmatch(idxmatch(csw03[:transition], csw03[:start_ell ])  , csw03[:z_hi])
slcts[8]    = idxmatch(idxmatch(csw03[:transition], csw03[:start_int ])  , csw03[:z_hi])
slcts[9]    = idxmatch(idxmatch(csw03[:transition], csw03[:start_disk])  , csw03[:z_hi])
labels      = Dict{Int64, String}()
labels[1]   = "Ell., $(latexstring("z < 0.45")), N = $(length(slcts[1])-count(isnan,y[slcts[1]]))"
labels[2]   = "Int., $(latexstring("z < 0.45")), N = $(length(slcts[2])-count(isnan,y[slcts[2]]))"
labels[3]   = "Disks, $(latexstring("z < 0.45")), N = $(length(slcts[3])-count(isnan,y[slcts[3]]))"
labels[4]   = "Ell., $(latexstring("0.45 < z < 1")), N = $(length(slcts[4])-count(isnan,y[slcts[4]]))"
labels[5]   = "Int., $(latexstring("0.45 < z < 1")), N = $(length(slcts[5])-count(isnan,y[slcts[5]]))"
labels[6]   = "Disks, $(latexstring("0.45 < z < 1")), N = $(length(slcts[6])-count(isnan,y[slcts[6]]))"
labels[7]   = "Ell., $(latexstring("z > 1")), N = $(length(slcts[7])-count(isnan,y[slcts[7]]))"
labels[8]   = "Int., $(latexstring("z > 1")), N = $(length(slcts[8])-count(isnan,y[slcts[8]]))"
labels[9]   = "Disks, $(latexstring("z > 1")), N = $(length(slcts[9])-count(isnan,y[slcts[9]]))"
labpos      = Array{String}(undef, 9)
labpos     .= "lower right"
plot_hexbins_nxm( x, 
            y, 
            slcts, labels, 3, 3,
            #scale=0.0,
            cbarshift=8,
            label_pos=labpos, 
            outfile=outfile, 
            grid=true,
            lognorm=true, gridres=(15,5), 
            ymin=miny, ymax=maxy,
            xmin=minx, xmax=maxx, 
            ylabel="log( sSFR [1/Gyr] )",
            xlabel="Flip Angle [°]", fontsize=30,
            binned_median="y", plot_bval=" ",
            )
#



###########################################################################
# merger predictions
include("/home/moon/sfortune/spinevo/plots/plot_merger_predictions.jl")

extrema(log10.(as03["M"] .- as03["ΔM"])[csw03[:start_thr]])







#####################################################################################
# J & j vs flip

for (ias, ia1, ia2, ia3, ic) in zip(as, as_label1, as_label2, as_label3, csw)
        outnm   = joinpath(outputdir,"Juppre_vs_flip_$(ia1)")
        outfile = outnm*"_hxcount.pdf"
        x   = pyplottable(ias["ϕ_flip"])
        y   = log10.(pyplottable(ias["J_main"] .- ias["ΔJ_main"]))
        miny= minimum(pyplottable(y)[findall(x->x .!== NaN, pyplottable(y))])
        maxy= maximum(pyplottable(y)[findall(x->x .!== NaN, pyplottable(y))])
        minx= 0#minimum(pyplottable(x)[findall(x->x .!== NaN, pyplottable(x))])
        maxx= 180#maximum(pyplottable(x)[findall(x->x .!== NaN, pyplottable(x))])
        slcts       = Dict{Int64, Vector{Int64}}()
        slcts[1]    = ic[:transition]
        slcts[2]    = idxmatch(ic[:transition], ic[:start_int ])
        slcts[3]    = idxmatch(ic[:transition], ic[:start_ell])
        slcts[4]    = idxmatch(ic[:transition], ic[:start_disk ])
        labels      = Dict{Int64, String}()
        labels[1]   = "Full Set, N = $(length(slcts[1])-count(isnan,y[slcts[1]]))"
        labels[2]   = "Intermediates, N = $(length(slcts[2])-count(isnan,y[slcts[2]]))"
        labels[3]   = "Ellipticals, N = $(length(slcts[3])-count(isnan,y[slcts[3]]))"
        labels[4]   = "Disks, N = $(length(slcts[4])-count(isnan,y[slcts[4]]))"
        labpos      = Array{String}(undef, 4)
        labpos     .= "lower right"
        plot_hexbins_nxm( x, 
                    y, 
                    slcts, labels, 2, 2,
                    #scale=0.0,
                    cbarshift=2.5,
                    label_pos=labpos, 
                    outfile=outfile, 
                    grid=true,
                    lognorm=true, gridres=(15,5), 
                    ymin=miny, ymax=maxy,
                    xmin=minx, xmax=maxx, 
                    ylabel="log( $(latexstring("\$J_{*\\textrm{,pre}}\$")) [ $(latexstring("\$M_\\odot\$")) kpc km/s ] )",
                    xlabel="$(ia3) Flip Angle [°]", fontsize=30,
                    binned_median="y", plot_bval=" ",
                    )
        #

        outnm   = joinpath(outputdir,"jpre_vs_flip_$(ia1)")
        outfile = outnm*"_hxcount.pdf"
        x   = pyplottable(ias["ϕ_flip"])
        y   = log10.(pyplottable(ias["j_main"] .- ias["Δj_main"]))
        miny= minimum(pyplottable(y)[findall(x->x .!== NaN, pyplottable(y))])
        maxy= maximum(pyplottable(y)[findall(x->x .!== NaN, pyplottable(y))])
        minx= 0#minimum(pyplottable(x)[findall(x->x .!== NaN, pyplottable(x))])
        maxx= 180#maximum(pyplottable(x)[findall(x->x .!== NaN, pyplottable(x))])
        labpos      = Array{String}(undef, 4)
        labpos     .= "lower right"
        plot_hexbins_nxm( x, 
                    y, 
                    slcts, labels, 2, 2,
                    #scale=0.0,
                    cbarshift=2.5,
                    label_pos=labpos, 
                    outfile=outfile, 
                    grid=true,
                    lognorm=true, gridres=(15,5), 
                    ymin=miny, ymax=maxy,
                    xmin=minx, xmax=maxx, 
                    ylabel="log( $(latexstring("\$j_{*\\textrm{,pre}}\$")) [ kpc km/s ] )",
                    xlabel="$(ia3) Flip Angle [°]", fontsize=30,
                    binned_median="y", plot_bval=" ",
                    )
        #
end

# after flip
for (ias, ia1, ia2, ia3, ic) in zip(as, as_label1, as_label2, as_label3, csw)
        outnm   = joinpath(outputdir,"Juppost_vs_flip_$(ia1)")
        outfile = outnm*"_hxcount.pdf"
        x   = pyplottable(ias["ϕ_flip"])
        y   = log10.(pyplottable(ias["J_main"]))
        miny= minimum(pyplottable(y)[findall(x->x .!== NaN, pyplottable(y))])
        maxy= maximum(pyplottable(y)[findall(x->x .!== NaN, pyplottable(y))])
        minx= 0#minimum(pyplottable(x)[findall(x->x .!== NaN, pyplottable(x))])
        maxx= 180#maximum(pyplottable(x)[findall(x->x .!== NaN, pyplottable(x))])
        slcts       = Dict{Int64, Vector{Int64}}()
        slcts[1]    = ic[:transition]
        slcts[2]    = idxmatch(ic[:transition], ic[:result_int ])
        slcts[3]    = idxmatch(ic[:transition], ic[:result_ell])
        slcts[4]    = idxmatch(ic[:transition], ic[:result_disk ])
        labels      = Dict{Int64, String}()
        labels[1]   = "Full Set, N = $(length(slcts[1])-count(isnan,y[slcts[1]]))"
        labels[2]   = "Intermediates, N = $(length(slcts[2])-count(isnan,y[slcts[2]]))"
        labels[3]   = "Ellipticals, N = $(length(slcts[3])-count(isnan,y[slcts[3]]))"
        labels[4]   = "Disks, N = $(length(slcts[4])-count(isnan,y[slcts[4]]))"
        labpos      = Array{String}(undef, 4)
        labpos     .= "lower right"
        labpos[1]   = "upper right"
        labpos[2]   = "upper right"
        labpos[3]   = "upper right"
        plot_hexbins_nxm( x, 
                    y, 
                    slcts, labels, 2, 2,
                    #scale=0.0,
                    cbarshift=2.5,
                    label_pos=labpos, 
                    outfile=outfile, 
                    grid=true,
                    lognorm=true, gridres=(15,5), 
                    ymin=miny, ymax=maxy,
                    xmin=minx, xmax=maxx, 
                    ylabel="log( $(latexstring("\$J_{*\\textrm{,post}}\$")) [ $(latexstring("\$M_\\odot\$")) kpc km/s ] )",
                    xlabel="$(ia3) Flip Angle [°]", fontsize=30,
                    binned_median="y", plot_bval=" ",
                    )
        #

        outnm   = joinpath(outputdir,"jpost_vs_flip_$(ia1)")
        outfile = outnm*"_hxcount.pdf"
        x   = pyplottable(ias["ϕ_flip"])
        y   = log10.(pyplottable(ias["j_main"]))
        miny= minimum(pyplottable(y)[findall(x->x .!== NaN, pyplottable(y))])
        maxy= maximum(pyplottable(y)[findall(x->x .!== NaN, pyplottable(y))])
        minx= 0#minimum(pyplottable(x)[findall(x->x .!== NaN, pyplottable(x))])
        maxx= 180#maximum(pyplottable(x)[findall(x->x .!== NaN, pyplottable(x))])
        labpos      = Array{String}(undef, 4)
        labpos     .= "lower right"
        labpos[1]   = "upper right"
        labpos[2]   = "upper right"
        labpos[3]   = "upper right"
        plot_hexbins_nxm( x, 
                    y, 
                    slcts, labels, 2, 2,
                    #scale=0.0,
                    cbarshift=2.5,
                    label_pos=labpos, 
                    outfile=outfile, 
                    grid=true,
                    lognorm=true, gridres=(15,5), 
                    ymin=miny, ymax=maxy,
                    xmin=minx, xmax=maxx, 
                    ylabel="log( $(latexstring("\$j_{*\\textrm{,post}}\$")) [ kpc km/s ] )",
                    xlabel="$(ia3) Flip Angle [°]", fontsize=30,
                    binned_median="y", plot_bval=" ",
                    )
        #
end


# pre flip but different redshifts
for (ias, ia1, ia2, ia3, ic) in zip(as, as_label1, as_label2, as_label3, csw)
        outnm   = joinpath(outputdir,"Juppre_vs_flip_zsplit_$(ia1)")
        outfile = outnm*"_hxcount.pdf"
        x   = pyplottable(ias["ϕ_flip"])
        y   = log10.(pyplottable(ias["J_main"] .- ias["ΔJ_main"]))
        miny= minimum(pyplottable(y)[findall(x->x .!== NaN, pyplottable(y))])
        maxy= maximum(pyplottable(y)[findall(x->x .!== NaN, pyplottable(y))])
        minx= 0#minimum(pyplottable(x)[findall(x->x .!== NaN, pyplottable(x))])
        maxx= 180#maximum(pyplottable(x)[findall(x->x .!== NaN, pyplottable(x))])
        slcts       = Dict{Int64, Vector{Int64}}()
        slcts[1]    = ic[:transition]
        slcts[2]    = idxmatch(ic[:transition], ic[:z_md ])
        slcts[3]    = idxmatch(ic[:transition], ic[:z_hi ])
        slcts[4]    = idxmatch(ic[:transition], ic[:z_lo])
        labels      = Dict{Int64, String}()
        labels[1]   = "Full Set, N = $(length(slcts[1])-count(isnan,y[slcts[1]]))"
        labels[2]   = "z = $(round(median(ias["redshift"][slcts[2]]), digits=2)), N = $(length(slcts[2])-count(isnan,y[slcts[2]]))"
        labels[3]   = "z = $(round(median(ias["redshift"][slcts[3]]), digits=2)), N = $(length(slcts[3])-count(isnan,y[slcts[3]]))"
        labels[4]   = "z = $(round(median(ias["redshift"][slcts[4]]), digits=2)), N = $(length(slcts[4])-count(isnan,y[slcts[4]]))"
        labpos      = Array{String}(undef, 4)
        labpos     .= "lower right"
        plot_hexbins_nxm( x, 
                    y, 
                    slcts, labels, 2, 2,
                    #scale=0.0,
                    cbarshift=2.5,
                    label_pos=labpos, 
                    outfile=outfile, 
                    grid=true,
                    lognorm=true, gridres=(15,5), 
                    ymin=miny, ymax=maxy,
                    xmin=minx, xmax=maxx, 
                    ylabel="log( $(latexstring("\$J_{*\\textrm{,pre}}\$")) [ $(latexstring("\$M_\\odot\$")) kpc km/s ] )",
                    xlabel="$(ia3) Flip Angle [°]", fontsize=30,
                    binned_median="y", plot_bval=" ",
                    )
        #

        outnm   = joinpath(outputdir,"jpre_vs_flip_zsplit_$(ia1)")
        outfile = outnm*"_hxcount.pdf"
        x   = pyplottable(ias["ϕ_flip"])
        y   = log10.(pyplottable(ias["j_main"] .- ias["Δj_main"]))
        miny= minimum(pyplottable(y)[findall(x->x .!== NaN, pyplottable(y))])
        maxy= maximum(pyplottable(y)[findall(x->x .!== NaN, pyplottable(y))])
        minx= 0#minimum(pyplottable(x)[findall(x->x .!== NaN, pyplottable(x))])
        maxx= 180#maximum(pyplottable(x)[findall(x->x .!== NaN, pyplottable(x))])
        labpos      = Array{String}(undef, 4)
        labpos     .= "lower right"
        plot_hexbins_nxm( x, 
                    y, 
                    slcts, labels, 2, 2,
                    #scale=0.0,
                    cbarshift=2.5,
                    label_pos=labpos, 
                    outfile=outfile, 
                    grid=true,
                    lognorm=true, gridres=(15,5), 
                    ymin=miny, ymax=maxy,
                    xmin=minx, xmax=maxx, 
                    ylabel="log( $(latexstring("\$j_{*\\textrm{,pre}}\$")) [ kpc km/s ] )",
                    xlabel="$(ia3) Flip Angle [°]", fontsize=30,
                    binned_median="y", plot_bval=" ",
                    )
        #
end









#####################################################################################
# bval pre flip

for (ias, ia1, ia2, ia3, ic) in zip(as, as_label1, as_label2, as_label3, csw)

        outnm   = joinpath(outputdir,"bvalpre_vs_flip_$(ia1)")
        outfile = outnm*"_hxcount.pdf"
        x   = pyplottable(ias["ϕ_flip"][ic[:transition]])
        y   = pyplottable(ias["BVAL_0"][ic[:transition]] .- ias["ΔBVAL_0"][ic[:transition]])
        miny= -6.1#minimum(pyplottable(y)[findall(x->x .!== NaN, pyplottable(y))])#-6.1#
        maxy= -3.9#maximum(pyplottable(y)[findall(x->x .!== NaN, pyplottable(y))])#-3.9#

        plot_hexbins( x, y, 
            outfile=joinpath(outputdir,outnm*"_hxcount.pdf"),
            scale=15, gridres=(15,7), xscale="linear", yscale="linear", 
            ylabel="$(latexstring("\$b_\\textrm{pre}\$"))", 
            xlabel="$(ia3) Flip Angle [°]", 
            grid=false,
            xmin=0, xmax=180, ymin=miny, ymax=maxy,
            lognorm=true,
            binned_median = "y", plot_bval="y", 
            calc_norm=true,    # @pyplottable
            )
        #

        colr = pyplottable(ias["BVAL_0"][ic[:transition]])
        plot_hexbins( x, y, 
            outfile=joinpath(outputdir,outnm*"_hxbpost.pdf"),
            scale=15, gridres=(15,7), xscale="linear", yscale="linear", 
            ylabel="$(latexstring("\$b_\\textrm{pre}\$"))", 
            xlabel="$(ia3) Flip Angle [°]", 
            grid=false,
            xmin=0, xmax=180, ymin=miny, ymax=maxy,
            lognorm=false,
            binned_median = "y", plot_bval="y", 
            cfunc=mean, colormod=colr, cmap=ColorMap(get_bvalue_cmap().colors), colorlimits=(-6.1,-3.9),
            calc_norm=true,    # @pyplottable
            )
        #
end









#####################################################################################
# flip vs delta b
# as03 csw03
for (ias, ia1, ia2, ia3, ic) in zip(as, as_label1, as_label2, as_label3, csw)

        println("  $(ia3)")
        outnm   = joinpath(outputdir,"flip_vs_delb_$(ia1)")
        outfile = outnm*"_hxcount.pdf"
        y   = pyplottable(ias["ϕ_flip"][ic[:transition]])
        x   = pyplottable(ias["ΔBVAL_0"][ic[:transition]])
        miny= 0#minimum(pyplottable(y)[findall(x->x .!== NaN, pyplottable(y))])#-6.1#
        maxy= 180#maximum(pyplottable(y)[findall(x->x .!== NaN, pyplottable(y))])#-3.9#
        minx= -1#minimum(pyplottable(x)[findall(x->x .!== NaN, pyplottable(x))])#-6.1#
        maxx= 1#maximum(pyplottable(x)[findall(x->x .!== NaN, pyplottable(x))])#-3.9#

        plot_hexbins( x, y, 
            outfile=joinpath(outputdir,outnm*"_hxcount.pdf"),
            scale=15, gridres=(15,7), xscale="linear", yscale="linear", 
            xlabel="$(ia3) $(latexstring("\$\\Delta b\$"))", 
            ylabel="$(ia3) Flip Angle [°]", 
            grid=false,
            xmin=minx, xmax=maxx, ymin=miny, ymax=maxy, 
            lognorm=true,
            binned_median = "y", plot_bval=" ", 
            calc_norm=true,    # @pyplottable
            )
        #

        colr = pyplottable(ias["BVAL_0"][ic[:transition]] .- ias["ΔBVAL_0"][ic[:transition]])
        plot_hexbins( x, y, 
            outfile=joinpath(outputdir,outnm*"_hxbpre.pdf"),
            scale=15, gridres=(15,7), xscale="linear", yscale="linear", 
            xlabel="$(ia3) $(latexstring("\$\\Delta b\$"))", 
            ylabel="$(ia3) Flip Angle [°]", 
            clabel="$(latexstring("\$b_\\textrm{pre}\$"))",
            grid=false,
            xmin=minx, xmax=maxx, ymin=miny, ymax=maxy,
            lognorm=false,
            binned_median = "y", plot_bval=" ", 
            cfunc=mean, colormod=colr, cmap=ColorMap(get_bvalue_cmap().colors), colorlimits=(-6.1,-3.9),
            calc_norm=true,    # @pyplottable
            )
        #
end









#####################################################################################
# endstate n flips
outnm   = joinpath(outputdir,"endstate_bval_Nflips")
outfile = outnm*".pdf"
y = vcat(
        pyplottable(endstate["nflips30"]), pyplottable(endstate["nflips45"]), pyplottable(endstate["nflips90"]), pyplottable(endstate["nflips135"]), 
        pyplottable(endstate["G09_nflips30"]), pyplottable(endstate["G09_nflips45"]), pyplottable(endstate["G09_nflips90"]), pyplottable(endstate["G09_nflips135"])
)
x = vcat(
        pyplottable(endstate["bval_end"]), pyplottable(endstate["bval_end"]), pyplottable(endstate["bval_end"]), pyplottable(endstate["bval_end"]), 
        pyplottable(endstate["bval_end"]), pyplottable(endstate["bval_end"]), pyplottable(endstate["bval_end"]), pyplottable(endstate["bval_end"])
)
miny= minimum(pyplottable(y)[findall(x->x .!== NaN, pyplottable(y))])
maxy= maximum(pyplottable(y)[findall(x->x .!== NaN, pyplottable(y))])
minx= -6.1#minimum(pyplottable(x)[findall(x->x .!== NaN, pyplottable(x))])
maxx= -3.9#maximum(pyplottable(x)[findall(x->x .!== NaN, pyplottable(x))])
slcts       = Dict{Int64, Vector{Int64}}()
slcts[1]    = 1:length(endstate["bval_end"])
slcts[2]    = length(slcts[1])+1:length(slcts[1])+length(endstate["bval_end"])
slcts[3]    = length(slcts[1])+length(slcts[2])+1:length(slcts[1])+length(slcts[2])+length(endstate["bval_end"])
slcts[4]    = length(slcts[1])+length(slcts[2])+length(slcts[3])+1:length(slcts[1])+length(slcts[2])+length(slcts[3])+length(endstate["bval_end"])
slcts[5]    = length(slcts[1])+length(slcts[2])+length(slcts[3])+length(slcts[4])+1:length(slcts[1])+length(slcts[2])+length(slcts[3])+length(slcts[4])+length(endstate["bval_end"])
slcts[6]    = length(endstate["bval_end"])*5+1:length(endstate["bval_end"])*6
slcts[7]    = length(endstate["bval_end"])*6+1:length(endstate["bval_end"])*7
slcts[8]    = length(endstate["bval_end"])*7+1:length(y)
labels      = Dict{Int64, String}()
labels[1]   = "$(latexstring("\$\\phi_\\textrm{flip} >\$ 30°, Instant"))"
labels[2]   = "$(latexstring("\$\\phi_\\textrm{flip} >\$ 45°, Instant"))"
labels[3]   = "$(latexstring("\$\\phi_\\textrm{flip} >\$ 90°, Instant"))"
labels[4]   = "$(latexstring("\$\\phi_\\textrm{flip} >\$ 135°, Instant"))"
labels[5]   = "$(latexstring("\$\\phi_\\textrm{flip} >\$ 30°, Long-term"))"
labels[6]   = "$(latexstring("\$\\phi_\\textrm{flip} >\$ 45°, Long-term"))"
labels[7]   = "$(latexstring("\$\\phi_\\textrm{flip} >\$ 90°, Long-term"))"
labels[8]   = "$(latexstring("\$\\phi_\\textrm{flip} >\$ 135°, Long-term"))"
labpos      = Array{String}(undef, 8)
labpos     .= "upper right"
plot_hexbins_nxm( x, 
            y, 
            slcts, labels, 4, 2,
            #scale=0.0, cbarshift=2.65
            cbarshift=23, cpad=0.05,
            outfile=outfile, 
            lognorm=true, gridres=(20,8), label_pos=labpos,
            ymin=miny, ymax=maxy,
            xmin=minx, xmax=maxx, 
            ylabel="N$(latexstring("\$_\\textrm{flip}\$"))",
            xlabel="$(latexstring("\$b(z=0)\$"))", fontsize=30,
            binned_median="y", plot_bval="x",
            dim_ratio=1.31, scale=15, tick_size=1.2,
            )
#


# endstate n mergers
outnm   = joinpath(outputdir,"endstate_bval_Nmergers")
outfile = outnm*".pdf"
y = vcat(
        pyplottable(endstate["nmergers2e8"]), pyplottable(endstate["nmergers2e9"]), pyplottable(endstate["nmergers2e10"]), 
        pyplottable(endstate["G09_nflips30"]), pyplottable(endstate["G09_nflips45"]), pyplottable(endstate["G09_nflips90"])
)
x = vcat(
        pyplottable(endstate["bval_end"]), pyplottable(endstate["bval_end"]), pyplottable(endstate["bval_end"]), 
        pyplottable(endstate["bval_end"]), pyplottable(endstate["bval_end"]), pyplottable(endstate["bval_end"])
)
miny= 0#minimum(pyplottable(y)[findall(x->x .!== NaN, pyplottable(y))])
maxy= 20#maximum(pyplottable(y)[findall(x->x .!== NaN, pyplottable(y))])
minx= -6.1#minimum(pyplottable(x)[findall(x->x .!== NaN, pyplottable(x))])
maxx= -3.9#maximum(pyplottable(x)[findall(x->x .!== NaN, pyplottable(x))])
slcts       = Dict{Int64, Vector{Int64}}()
slcts[1]    = 1:length(endstate["bval_end"])
slcts[2]    = length(slcts[1])+1:length(slcts[1])+length(endstate["bval_end"])
slcts[3]    = length(slcts[1])+length(slcts[2])+1:length(slcts[1])+length(slcts[2])+length(endstate["bval_end"])
slcts[4]    = length(slcts[1])+length(slcts[2])+length(slcts[3])+1:length(slcts[1])+length(slcts[2])+length(slcts[3])+length(endstate["bval_end"])
slcts[5]    = length(slcts[1])+length(slcts[2])+length(slcts[3])+length(slcts[4])+1:length(slcts[1])+length(slcts[2])+length(slcts[3])+length(slcts[4])+length(endstate["bval_end"])
slcts[6]    = length(endstate["bval_end"])*5+1:length(y)
labels      = Dict{Int64, String}()
labels[1]   = "$(latexstring("\$M_{*\\textrm{,mergers}} >\$ 2e8, Instant"))"
labels[2]   = "$(latexstring("\$M_{*\\textrm{,mergers}} >\$ 2e9, Instant"))"
labels[3]   = "$(latexstring("\$M_{*\\textrm{,mergers}} >\$ 2e10, Instant"))"
labels[4]   = "$(latexstring("\$M_{*\\textrm{,mergers}} >\$ 2e8, Long-term"))"
labels[5]   = "$(latexstring("\$M_{*\\textrm{,mergers}} >\$ 2e9, Long-term"))"
labels[6]   = "$(latexstring("\$M_{*\\textrm{,mergers}} >\$ 2e10, Long-term"))"
labpos      = Array{String}(undef, 6)
labpos     .= "upper right"
plot_hexbins_nxm( x, 
            y, 
            slcts, labels, 3, 2,
            #scale=0.0,
            cbarshift=23, cpad=0.05,
            outfile=outfile, 
            lognorm=true, gridres=(20,10), label_pos=labpos,
            ymin=miny, ymax=maxy,
            xmin=minx, xmax=maxx, 
            ylabel="N$(latexstring("\$_\\textrm{mergers}\$"))",
            xlabel="$(latexstring("\$b(z=0)\$"))", fontsize=30,
            binned_median="y", plot_bval="x",
            dim_ratio=1.31, scale=15, tick_size=1.2,
            )
#









#####################################################################################
# spin change vs mass pre
#
for (ias, ia1, ia2, ia3, ic) in zip(as, as_label1, as_label2, as_label3, csw)
        for (sp, lbl, ylbl,ylbl2) in zip(["J", "j"], ["Jup", "j"], ["log( $(latexstring("\$|\\Delta J_*|\$")) [ $(latexstring("\$M_\\odot\$")) kpc km/s ] )", "log( $(latexstring("\$|\\Delta j_*|\$")) [ kpc km/s ] )"], ["log( $(latexstring("\$\\Delta J_*\$")) / $(latexstring("\$J_\\textrm{*,pre}\$")) )", "log( $(latexstring("\$\\Delta j_*\$")) / $(latexstring("\$j_\\textrm{*,pre}\$")) )"])
                println("  $(ia3)")
                outnm   = joinpath(outputdir,"del$(lbl)_vs_Mpre_$(ia1)")
                outfile = outnm*"_hxcount.pdf"
                y   = log10.(pyplottable(ias["Δ$(sp)_main"][:,ic[:transition]]))
                x   = log10.(pyplottable(ias["M"][ic[:transition]] .- ias["ΔM"][ic[:transition]]))
                miny= minimum(pyplottable(y[findall(x->x .!== NaN, pyplottable(y))]))#-6.1#
                maxy= maximum(pyplottable(y[findall(x->x .!== NaN, pyplottable(y))]))#-3.9#
                minx= minimum(pyplottable(x[findall(x->x .!== NaN, pyplottable(x))]))#-6.1#
                maxx= maximum(pyplottable(x[findall(x->x .!== NaN, pyplottable(x))]))#-3.9#
                plot_hexbins( x, y, 
                    outfile=outfile,
                    scale=15, gridres=(15,7), xscale="linear", yscale="linear", 
                    xlabel="log( $(latexstring("\$M_{*,\\textrm{pre}}\$")) [$(latexstring("\$M_\\odot\$"))] )", 
                    ylabel = ylbl, 
                    clabel="N",
                    grid=true,
                    xmin=minx, xmax=maxx, ymin=miny, ymax=maxy, 
                    lognorm=true,
                    binned_median = "y", plot_bval=" ", 
                    calc_norm=true,    # @pyplottable
                    )
                #
        
                colr = pyplottable(ias["ϕ_flip"][ic[:transition]])
                outfile = outnm*"_hxflip.pdf"
                plot_hexbins( x, y, 
                    outfile=outfile,
                    scale=15, gridres=(15,7), xscale="linear", yscale="linear", 
                    xlabel="log( $(latexstring("\$M_{*,\\textrm{pre}}\$")) [$(latexstring("\$M_\\odot\$"))] )", 
                    ylabel = ylbl, 
                    clabel="Flip Angle [°]",
                    grid=true,
                    xmin=minx, xmax=maxx, ymin=miny, ymax=maxy,
                    lognorm=false,
                    binned_median = "y", plot_bval=" ", 
                    cfunc=mean, colormod=colr, cmap="plasma_r", colorlimits=(0,180),
                    calc_norm=true,    # @pyplottable
                    )
                #


                # relative
                
                outnm   = joinpath(outputdir,"del$(lbl)rel_vs_Mpre_$(ia1)")
                outfile = outnm*"_hxcount.pdf"
                y   = log10.( pyplottable(ias["Δ$(sp)_main"][:,ic[:transition]]) ./ pyplottable(ias["$(sp)_main"][:,ic[:transition]] .- ias["Δ$(sp)_main"][:,ic[:transition]])  )
                x   = log10.(pyplottable(ias["M"][ic[:transition]] .- ias["ΔM"][ic[:transition]]))
                miny= minimum(pyplottable(y[findall(x->x .!== NaN, pyplottable(y))]))#-6.1#
                maxy= maximum(pyplottable(y[findall(x->x .!== NaN, pyplottable(y))]))#-3.9#
                minx= minimum(pyplottable(x[findall(x->x .!== NaN, pyplottable(x))]))#-6.1#
                maxx= maximum(pyplottable(x[findall(x->x .!== NaN, pyplottable(x))]))#-3.9#
                plot_hexbins( x, y, 
                    outfile=outfile,
                    scale=15, gridres=(15,7), xscale="linear", yscale="linear", 
                    xlabel="log( $(latexstring("\$M_{*,\\textrm{pre}}\$")) [$(latexstring("\$M_\\odot\$"))] )", 
                    ylabel = ylbl2, 
                    clabel="N",
                    grid=true,
                    xmin=minx, xmax=maxx, ymin=miny, ymax=maxy, 
                    lognorm=true,
                    binned_median = "y", plot_bval=" ", 
                    calc_norm=true,    # @pyplottable
                    )
                #
        
                colr = pyplottable(ias["ϕ_flip"][ic[:transition]])
                outfile = outnm*"_hxflip.pdf"
                plot_hexbins( x, y, 
                    outfile=outfile,
                    scale=15, gridres=(15,7), xscale="linear", yscale="linear", 
                    xlabel="log( $(latexstring("\$M_{*,\\textrm{pre}}\$")) [$(latexstring("\$M_\\odot\$"))] )", 
                    ylabel = ylbl2, 
                    clabel="Flip Angle [°]",
                    grid=true,
                    xmin=minx, xmax=maxx, ymin=miny, ymax=maxy,
                    lognorm=false,
                    binned_median = "y", plot_bval=" ", 
                    cfunc=mean, colormod=colr, cmap="plasma_r", colorlimits=(0,180),
                    calc_norm=true,    # @pyplottable
                    )
                #
        end
end









#####################################################################################
# delta j vs flip ssfr split

for (ias, ia1, ia2, ia3, ic) in zip(as, as_label1, as_label2, as_label3, csw)
        println(ia3)
        outnm   = joinpath(outputdir,"delj_vs_flip_sSFRsplit_$(ia1)")
        outfile = outnm*"_hxcount.pdf"
        x   = pyplottable(ias["ϕ_flip"])
        y   = pyplottable(log10.(pyplottable(ias["Δj_main"])))
        insane_y = findcs(y, eq=NaN)
        insane_x = findcs(x, eq=NaN)
        sane     = idxclude(idxclude(1:length(x), insane_x), insane_y)
        miny= minimum(pyplottable(y)[findall(x->x .!== NaN, pyplottable(y))])
        maxy= maximum(pyplottable(y)[findall(x->x .!== NaN, pyplottable(y))])
        minx= 0#minimum(pyplottable(x)[findall(x->x .!== NaN, pyplottable(x))])
        maxx= 180#maximum(pyplottable(x)[findall(x->x .!== NaN, pyplottable(x))])
        slcts       = Dict{Int64, Vector{Int64}}()
        slcts[1]    = idxmatch(ic[:transition], vcat(ic[:sSFR_lo ], ic[:sSFR_md ], ic[:sSFR_hi ]))
        slcts[3]    = idxmatch(ic[:transition], ic[:sSFR_lo ])
        slcts[2]    = idxmatch(ic[:transition], ic[:sSFR_md ])
        slcts[4]    = idxmatch(ic[:transition], ic[:sSFR_hi ])
        labels      = Dict{Int64, String}()
        labels[1]   = "Full Set, N = $(length(slcts[1])-count(isnan,y[slcts[1]]))"
        labels[3]   = "Low sSFR, N = $(length(slcts[2])-count(isnan,y[slcts[2]]))"
        labels[2]   = "Medium sSFR, N = $(length(slcts[3])-count(isnan,y[slcts[3]]))"
        labels[4]   = "High sSFR, N = $(length(slcts[4])-count(isnan,y[slcts[4]]))"
        labpos      = Array{String}(undef, 4)
        labpos     .= "lower right"
        plot_hexbins_nxm( x, 
                    y, 
                    slcts, labels, 2, 2,
                    #scale=0.0,
                    cbarshift=2.5,
                    label_pos=labpos, 
                    outfile=outfile, 
                    grid=true,
                    lognorm=true, gridres=(15,5), 
                    ymin=miny, ymax=maxy,
                    xmin=minx, xmax=maxx, 
                    ylabel="log( $(latexstring("\$\\Delta j_*\$")) [ kpc km/s ] )",
                    xlabel="$(ia3) Flip Angle [°]", fontsize=30,
                    binned_median="y", plot_bval=" ",
                    )
        #

        # Jup
        outnm   = joinpath(outputdir,"delJup_vs_flip_sSFRsplit_$(ia1)")
        outfile = outnm*"_hxcount.pdf"
        x   = pyplottable(ias["ϕ_flip"])
        y   = pyplottable(log10.(pyplottable(ias["ΔJ_main"])))
        miny= minimum(pyplottable(y)[findall(x->x .!== NaN, pyplottable(y))])
        maxy= maximum(pyplottable(y)[findall(x->x .!== NaN, pyplottable(y))])
        minx= 0#minimum(pyplottable(x)[findall(x->x .!== NaN, pyplottable(x))])
        maxx= 180#maximum(pyplottable(x)[findall(x->x .!== NaN, pyplottable(x))])
        slcts       = Dict{Int64, Vector{Int64}}()
        slcts[1]    = idxmatch(ic[:transition], vcat(ic[:sSFR_lo ], ic[:sSFR_md ], ic[:sSFR_hi ]))
        slcts[3]    = idxmatch(ic[:transition], ic[:sSFR_lo ])
        slcts[2]    = idxmatch(ic[:transition], ic[:sSFR_md ])
        slcts[4]    = idxmatch(ic[:transition], ic[:sSFR_hi ])
        labels      = Dict{Int64, String}()
        labels[1]   = "Full Set, N = $(length(slcts[1])-count(isnan,y[slcts[1]]))"
        labels[3]   = "Low sSFR, N = $(length(slcts[2])-count(isnan,y[slcts[2]]))"
        labels[2]   = "Medium sSFR, N = $(length(slcts[3])-count(isnan,y[slcts[3]]))"
        labels[4]   = "High sSFR, N = $(length(slcts[4])-count(isnan,y[slcts[4]]))"
        labpos      = Array{String}(undef, 4)
        labpos     .= "lower right"
        plot_hexbins_nxm( x, 
                    y, 
                    slcts, labels, 2, 2,
                    #scale=0.0,
                    cbarshift=2.5,
                    label_pos=labpos, 
                    outfile=outfile, 
                    grid=true,
                    lognorm=true, gridres=(15,5), 
                    ymin=miny, ymax=maxy,
                    xmin=minx, xmax=maxx, 
                    ylabel="log( $(latexstring("\$\\Delta J_*\$")) [ $(latexstring("\$M_\\odot\$")) kpc km/s ] )",
                    xlabel="$(ia3) Flip Angle [°]", fontsize=30,
                    binned_median="y", plot_bval=" ",
                    )
        #
end

flip_lim = 30
allssfrratio = length(idxmatch( idxmatch(csw09[:transition],vcat(csw09[:sSFR_lo ], csw09[:sSFR_md ], csw09[:sSFR_hi ])), findcs(as09["ϕ_flip"],geq=flip_lim) )) / length(idxmatch( idxmatch(csw09[:transition],vcat(csw09[:sSFR_lo ], csw09[:sSFR_md ], csw09[:sSFR_hi ])), findcs(as09["ϕ_flip"],lt=flip_lim) ))
lowssfrratio = length(idxmatch( idxmatch(csw09[:transition],csw09[:sSFR_lo]), findcs(as09["ϕ_flip"],geq=flip_lim) )) / length(idxmatch( idxmatch(csw09[:transition],csw09[:sSFR_lo]), findcs(as09["ϕ_flip"],lt=flip_lim) ))
midssfrratio = length(idxmatch( idxmatch(csw09[:transition],csw09[:sSFR_md]), findcs(as09["ϕ_flip"],geq=flip_lim) )) / length(idxmatch( idxmatch(csw09[:transition],csw09[:sSFR_md]), findcs(as09["ϕ_flip"],lt=flip_lim) ))
highssfrratio = length(idxmatch( idxmatch(csw09[:transition],csw09[:sSFR_hi]), findcs(as09["ϕ_flip"],geq=flip_lim) )) / length(idxmatch( idxmatch(csw09[:transition],csw09[:sSFR_hi]), findcs(as09["ϕ_flip"],lt=flip_lim) ))





##########################################################################################
#  flip vs M
for (ias, ia1, ia2, ia3, ic) in zip(as, as_label1, as_label2, as_label3, csw)
        println("  $(ia3)")
        for (mlbl, x_src) in zip(["pre", "post"], [ ias["M"][ic[:transition]] .- ias["ΔM"][ic[:transition]], ias["M"][ic[:transition]] ])
                outnm   = joinpath(outputdir,"flip_vs_M$(mlbl)_$(ia1)")
                outfile = outnm*"_hxcount.pdf"
                y   = pyplottable(ias["ϕ_flip"][ic[:transition]])
                x   = log10.(pyplottable(x_src))
                miny= 0#minimum(pyplottable(y[findall(x->x .!== NaN, pyplottable(y))]))#-6.1#
                maxy= 180#maximum(pyplottable(y[findall(x->x .!== NaN, pyplottable(y))]))#-3.9#
                minx= minimum(pyplottable(x[findall(x->x .!== NaN, pyplottable(x))]))#-6.1#
                maxx= maximum(pyplottable(x[findall(x->x .!== NaN, pyplottable(x))]))#-3.9#
                plot_hexbins( x, y, 
                    outfile=outfile,
                    scale=15, gridres=(30,13), xscale="linear", yscale="linear", 
                    xlabel="log( $(latexstring("\$M_{*,\\textrm{$(mlbl)}}\$")) [$(latexstring("\$M_\\odot\$"))] )", 
                    ylabel = "$(ia3) Flip Angle [°]", 
                    clabel="N",
                    grid=false,
                    xmin=minx, xmax=maxx, ymin=miny, ymax=maxy, 
                    lognorm=true,
                    binned_median = "y", plot_bval=" ", 
                    calc_norm=true,    # @pyplottable
                    )
                #
        
                #colr = pyplottable(abs.(ias["ΔM"][ic[:transition]] ./ x_src ))
                colr = pyplottable(symlog10.(ias["ΔM"][ic[:transition]] ))
                outfile = outnm*"_hxflip.pdf"
                plot_hexbins( x, y, 
                    outfile=outfile,
                    scale=15, gridres=(30,13), xscale="linear", yscale="linear", 
                    xlabel="log( $(latexstring("\$M_{*,\\textrm{$(mlbl)}}\$")) [$(latexstring("\$M_\\odot\$"))] )", 
                    ylabel = "$(ia3) Flip Angle [°]", 
                    clabel="log( $(latexstring("\$\\Delta M_*\$")) [$(latexstring("\$M_\\odot\$"))] )",
                    #clabel="log( $(latexstring("\$\\Delta M_*\$")) / $(latexstring("\$M_{*,\\textrm{$(mlbl)}}\$")) )",
                    grid=false,
                    xmin=minx, xmax=maxx, ymin=miny, ymax=maxy,
                    lognorm=false,
                    binned_median = "y", plot_bval=" ", 
                    cfunc=mean, colormod=colr, cmap="plasma_r", #colorlimits=(0,180),
                    calc_norm=true,    # @pyplottable
                    )
                #
        end
end








#####################################################################################
# b-value vmin=-6.1,vmax=-3.9
outfile     = joinpath(outputdir,"hist_bvalues.pdf")
scale = 0.5
lw=2
N_bins = 50
PyPlot.matplotlib[:rc]("text", usetex=true)
PyPlot.matplotlib[:rc]("font", family="serif", serif=["Computer Modern"])
fig, ax = subplots()
ax.hist( as03["BVAL"]    , bins=N_bins    , rwidth=0.9          ,linewidth=lw+2, alpha=1, color="darkred", histtype="step", label="Uncorrected", density=false )
ax.hist( as03["BVAL_0"]    , bins=N_bins    , rwidth=0.9        ,linewidth=lw+1, alpha=1, color="navy", histtype="step", label="Full Sample", density=false )
ax.hist( as03["BVAL_0"][csw03[:transition]]    , bins=N_bins    ,linewidth=lw-1, rwidth=0.9    , alpha=1, color="cyan", histtype="step", label="Transition Sample", density=false  )
ax.hist( as03["BVAL_0"][csw03[:z_0]]    , bins=N_bins    ,linewidth=lw, rwidth=0.9    , alpha=1, color="black", histtype="step", label="Sample at $(latexstring("\$z = 0\$"))", density=false  )
ax.axvline(b_disk, color="mediumblue",linewidth=lw+2)
ax.axvline(b_ell, color="red",linewidth=lw+1)
ax.axvline(-6.1, color="black",linewidth=lw+1)
ax.axvline(-3.9, color="black",linewidth=lw+1)
ax.set_xlabel("$(latexstring("\$b\$"))")
ax.grid(false)
ax.legend(loc="upper left", frameon=true, borderpad=1, handlelength=1.8)
fig.set_size_inches(11scale, 9scale)
fig.tight_layout()
fig.savefig(outfile, bbox_inches="tight", pad_inches=.1)
close("all")
PyPlot.matplotlib[:rc]("text", usetex=false)





##############################################################################
# bval vs time

peak_idcs = findcs(as03["redshift"], lt=2.1, gt=1.7)
extrema(as03["redshift"][peak_idcs])
extrema(as03["lookbacktime"][peak_idcs])

outnm   = joinpath(outputdir,"bval_vs_redshift")
outfile = outnm*"_hxcount.pdf"
y   = as03["BVAL_0"]
x   = as03["redshift"]
miny= minimum(pyplottable(y[findall(x->x .!== NaN, pyplottable(y))]))#-6.1#
maxy= maximum(pyplottable(y[findall(x->x .!== NaN, pyplottable(y))]))#-3.9#
minx= minimum(pyplottable(x[findall(x->x .!== NaN, pyplottable(x))]))#-6.1#
maxx= maximum(pyplottable(x[findall(x->x .!== NaN, pyplottable(x))]))#-3.9#
plot_hexbins( x, y, 
    outfile=outfile,
    scale=15, gridres=(15,7), xscale="linear", yscale="linear", 
    xlabel="$(latexstring("\$z\$"))", 
    ylabel = "$(latexstring("\$b\$"))", 
    clabel="N",
    grid=false,
    xmin=minx, xmax=maxx, ymin=miny, ymax=maxy, 
    lognorm=true,
    binned_median = "y", plot_bval="y", plot_vline=1.9,
    calc_norm=true,    # @pyplottable
    )
#


outnm   = joinpath(outputdir,"bval_vs_lbt")
outfile = outnm*"_hxcount.pdf"
y   = as03["BVAL_0"]
x   = as03["lookbacktime"]
miny= minimum(pyplottable(y[findall(x->x .!== NaN, pyplottable(y))]))#-6.1#
maxy= maximum(pyplottable(y[findall(x->x .!== NaN, pyplottable(y))]))#-3.9#
minx= minimum(pyplottable(x[findall(x->x .!== NaN, pyplottable(x))]))#-6.1#
maxx= maximum(pyplottable(x[findall(x->x .!== NaN, pyplottable(x))]))#-3.9#
plot_hexbins( x, y, 
    outfile=outfile,
    scale=15, gridres=(15,7), xscale="linear", yscale="linear", 
    xlabel="Lookback Time [Gyr]", 
    ylabel = "$(latexstring("\$b\$"))", 
    clabel="N",
    grid=false,
    xmin=minx, xmax=maxx, ymin=miny, ymax=maxy, 
    lognorm=true,
    binned_median = "y", plot_bval="y", plot_vline=10.0,
    calc_norm=true,    # @pyplottable
    )
#





# b-value vs SFR at different redshift bins
outnm   = joinpath(outputdir,"sSFR_vs_bval_zsplit")
outfile = outnm*"_hxcount.pdf"
x   = pyplottable(as03["BVAL_0"])
y   = pyplottable(log10.(1e9 .* as03["SFR"] ./ as03["M"] ))
insane_y = findcs(y, eq=NaN)
insane_x = findcs(x, eq=NaN)
sane     = idxclude(idxclude(1:length(x), insane_x), insane_y)
#vmin=-6.1,vmax=-3.9
miny= minimum(pyplottable(y)[findall(x->x .!== NaN, pyplottable(y))])
maxy= maximum(pyplottable(y)[findall(x->x .!== NaN, pyplottable(y))])
minx= -6.1#minimum(pyplottable(x)[findall(x->x .!== NaN, pyplottable(x))])
maxx= -3.9#maximum(pyplottable(x)[findall(x->x .!== NaN, pyplottable(x))])
slcts       = Dict{Int64, Vector{Int64}}()
slcts[1]    = sane
slcts[2]    = idxmatch(sane, idxmatch(csw03[:z_postpeak], findcs( as03["redshift"         ], gt=0.5 )))
slcts[4]    = idxmatch(sane, csw03[:z_prepeak] )
slcts[3]    = idxmatch(sane, idxmatch(csw03[:z_postpeak], findcs( as03["redshift"         ], leq=0.5 )))
labels      = Dict{Int64, String}()
labels[1]   = "Full Set, N = $(length(slcts[1]))"
labels[2]   = "$(z_peak) $(latexstring("\$>\$")) z $(latexstring("\$>\$")) 0.5, N = $(length(slcts[2]))"
labels[4]   = "z $(latexstring("\$\\geq\$")) $(z_peak), N = $(length(slcts[4]))"
labels[3]   = "z $(latexstring("\$\\leq\$")) 0.5, N = $(length(slcts[3]))"
labpos      = Array{String}(undef, 4)
labpos     .= "lower right"
plot_hexbins_nxm( x, 
            y, 
            slcts, labels, 2, 2,
            label_pos=labpos, 
            #scale=0.0,
            cbarshift=2.5,
            outfile=outfile, 
            lognorm=true, gridres=(30,13), 
            ymin=miny, ymax=maxy,
            xmin=minx, xmax=maxx, 
            ylabel="log( sSFR [ 1 / Gyr] )",
            xlabel="$(latexstring("\$b\$"))", fontsize=30,
            binned_median="y", plot_bval="x",
            )
#

# fitting to the above
linfct(x,p) = p[1] .+ p[2] .* (x .+ 0)

xdata_al =      x[sane]
ydata_al =      y[sane]
xdata_hi =      x[idxmatch(sane,slcts[4])]
ydata_hi =      y[idxmatch(sane,slcts[4])]
xdata_md =      x[idxmatch(sane,slcts[2])]
ydata_md =      y[idxmatch(sane,slcts[2])]
xdata_lo =      x[idxmatch(sane,slcts[3])]
ydata_lo =      y[idxmatch(sane,slcts[3])]
p0 = [1.5, 0.5]
fit_al  = curve_fit(linfct, xdata_al, ydata_al, p0)
fit_hi  = curve_fit(linfct, xdata_hi, ydata_hi, p0)
fit_md  = curve_fit(linfct, xdata_md, ydata_md, p0)
fit_lo  = curve_fit(linfct, xdata_lo, ydata_lo, p0)

# Plot it
outfile = outnm*"_linfits.pdf"
pointsize = 1
scale = 0.5
PyPlot.matplotlib[:rc]("text", usetex=true)
PyPlot.matplotlib[:rc]("font", family="serif", serif=["Computer Modern"])
fig, ax = subplots()
PyPlot.matplotlib[:rc]("text", usetex=true)
PyPlot.matplotlib[:rc]("font", family="serif", serif=["Computer Modern"])
ax.axvline(b_disk, color="blue")
ax.axvline(b_ell, color="red")
#ax.plot( xdata_al, linfct(xdata_al, fit_al.param), "-", color="black", label="$(labels[1])")
xplot = LinRange(-7, 1, 1000)
#ax.scatter( x[idxmatch(sane,slcts[3])], y[idxmatch(sane,slcts[3])]  , s= pointsize    , alpha=0.1, c="mediumblue") # brg_r  
#ax.scatter( x[idxmatch(sane,slcts[2])], y[idxmatch(sane,slcts[2])]  , s= pointsize    , alpha=0.1, c="darkorchid") # brg_r  
#ax.scatter( x[idxmatch(sane,slcts[4])], y[idxmatch(sane,slcts[4])]  , s= pointsize    , alpha=0.1, c="darkred") # brg_r  
ax.plot( xplot, linfct(xplot, fit_lo.param), "-", color="mediumblue", label="z $(latexstring("\$\\leq\$")) 0.5")
ax.plot( xplot, linfct(xplot, fit_md.param), "-", color="darkorchid", label="$(z_peak) $(latexstring("\$>\$")) z $(latexstring("\$>\$")) 0.5")
ax.plot( xplot, linfct(xplot, fit_hi.param), "-", color="darkred"   , label="z $(latexstring("\$\\geq\$")) $(z_peak)")
em1 = 1
em2 = 0
ax.fill_between(xplot, linfct(xplot, fit_lo.param .+ (standard_errors(fit_lo).*[-em1,em2])), y2=linfct(xplot, fit_lo.param .+ (standard_errors(fit_lo).*[em1,-em2])), color="mediumblue", alpha=0.2)
ax.fill_between(xplot, linfct(xplot, fit_md.param .+ (standard_errors(fit_md).*[-em1,em2])), y2=linfct(xplot, fit_md.param .+ (standard_errors(fit_md).*[em1,-em2])), color="darkorchid", alpha=0.2)
ax.fill_between(xplot, linfct(xplot, fit_hi.param .+ (standard_errors(fit_hi).*[-em1,em2])), y2=linfct(xplot, fit_hi.param .+ (standard_errors(fit_hi).*[em1,-em2])), color="darkred"   , alpha=0.2)
ax.set_xlabel("$(latexstring("\$b\$"))")
ax.set_ylabel("log( sSFR [ 1 / Gyr] )")
ax.set_xlim([-6.1,-3.9])
ax.set_ylim([-3,0])
ax.legend(loc="lower right", frameon=true, borderpad=1, handlelength=1.8)
ax.grid(false)
fig.set_size_inches(9scale, 9scale)
fig.tight_layout()
fig.savefig(outfile, bbox_inches="tight", pad_inches=.1)
close("all")
PyPlot.matplotlib[:rc]("text", usetex=false)







###########################################################################################
# halo stories
id_list = [1414, 3212, 18369, 21029, 22253, 14587, 14635, 14679, 14755, 16508]

for i in id_list
        plot_halostory_full(outfile=joinpath(outputdir, "halo_$(i)_halostory_full.pdf"), snap=false, rootID=i)
end

hf_list = readdir("/home/moon/sfortune/spinevo/data/silvio_stories_spinmap")
outputdir = "/home/moon/sfortune/spinevo/plots/hs_full"
@showprogress for i in 1:length(hf_list)
        isub            = parse( Int64, chop( hf_list[i], head=5,tail=4 ) )
        plot_halostory_full(outfile=joinpath(outputdir, "halo_$(isub)_halostory_full.pdf"), snap=true, rootID=isub, title=true)
end
outputdir = "/e/ldata/users/sfortune/Dropbox/Apps/Overleaf/SpinEvolution/figures"

#####################################################################################
# sSFR vs flipangle bval z split












# central-switch-Erratic subs:
# 25271, 26278