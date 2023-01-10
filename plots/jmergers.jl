
#include("/home/moon/sfortune/spinevo/pkg/meta.jl")
#
#outputdir = "/home/moon/sfortune/spinevo/plots/jmergers"
#
#include("/home/moon/sfortune/spinevo/data/load_data.jl")




# merger J contribution
# components: 
#       J_orb, J_int, 
#       immediate, before, 
#       sum, most massive, 
#       0.3 Gyr, 0.9 Gyr
# 1=mass, 2=mass2, 3=FPmass, 4=FPmass2, 5=FP SNAP, 6=FP z, 7= FP lbt, 8=subID, 9=mass_peak, 10=FPmass_peak, 11=snapNR
     # 12 = snap immediate, 13 = lbt immediate, 14 = mass immediate, 15:17 = orbital immediate, 18:20 intrinsic immediate, 
     # 21 = snap later, 22 = lbt later, 23 = mass later, 24:26 orbital later, 27:29, intrinsic later
mmi_as03 = mergermap_indices(as03)
mmi_as09 = mergermap_indices(as09)
mmi_asn03 = mergermap_indices(asn03)
mmi_asn09 = mergermap_indices(asn09)

# orb, immediate, sum, 0.3
idx_base = idxclude(idxclude(csw03[:start_thr], csw03[:switches]), csw03[:fakes])
delJ = convert( Array{Float64,2}, as03["ΔJ_main"][:,idxmatch(mmi_as03["main"], idx_base )] ) # all entries with mergers and matching idx_base
J_0  = convert( Array{Float64,2}, as03["J_main"][:,idxmatch(mmi_as03["main"], idx_base )] .- as03["ΔJ_main"][:,idxmatch(mmi_as03["main"], idx_base )] )
J_1  = convert( Array{Float64,2}, as03["J_main"][:,idxmatch(mmi_as03["main"], idx_base )])
idx = 1:length(mmi_as03["main"])
idx = idx[mmi_as03["main"] .∈ Ref(Set(idx_base))] # selection for merger map indices
Jom  = missings(Float64, 3, 0)
JoMM = missings(Float64, 3, 0)
Jim  = missings(Float64, 3, 0)
JiMM = missings(Float64, 3, 0)
Jsm  = missings(Float64, 3, 0)
JsMM = missings(Float64, 3, 0)
idx_missed = Vector{Int64}(undef, 0) # in merger map
idx_caught = Vector{Int64}(undef, 0) # in merger map
for i in idx # iterate over y entries
        JoMM = hcat( JoMM,  as03["merger map"][15:17,mmi_as03["most_massive"][i]] .* as03["merger map"][14,mmi_as03["most_massive"][i]])
        JiMM = hcat( JiMM,  as03["merger map"][18:20,mmi_as03["most_massive"][i]] )
        if count(ismissing, JiMM[:,end]) == 0
                JsMM = hcat( JsMM, JiMM[:,end] .+ JoMM[:,end] )
        else
                JsMM = hcat( JsMM, JoMM[:,end] )
        end
        Jo_sum = zeros(3)
        Ji_sum = zeros(3)
        for ii in mmi_as03["first"][i]:mmi_as03["last"][i] # iterate over merger entries in merger map
                if count(ismissing, as03["merger map"][15:17,ii]) == 0
                        Jo_sum .+= as03["merger map"][15:17,ii] .* as03["merger map"][14,ii]
                end
                if count(ismissing, as03["merger map"][18:20,ii]) == 0
                        Ji_sum .+= as03["merger map"][18:20,ii]
                        idx_caught = vcat( idx_caught, ii )
                else 
                        idx_missed = vcat( idx_missed, ii )
                end
        end
        Jom = hcat( Jom,  ifelse(count(iszero, Jo_sum) == 3, missings(3), Jo_sum) )
        Jim = hcat( Jim,  ifelse(count(iszero, Ji_sum) == 3, missings(3), Ji_sum) ) # missings if zero vector
        Jsm = hcat( Jsm,  ifelse(count(iszero, Jo_sum) == 3, missings(3), Jo_sum) .+ Ji_sum )
end

angle_J0_Jom  = missings(Float64, length(J_0[1,:]))
angle_J0_Jim  = missings(Float64, length(J_0[1,:]))
angle_J0_Jsm  = missings(Float64, length(J_0[1,:]))
angle_J0_JoMM = missings(Float64, length(J_0[1,:]))
angle_J0_JiMM = missings(Float64, length(J_0[1,:]))
angle_J0_JsMM = missings(Float64, length(J_0[1,:]))
angle_J1_Jom  = missings(Float64, length(J_0[1,:]))
angle_J1_Jim  = missings(Float64, length(J_0[1,:]))
angle_J1_Jsm  = missings(Float64, length(J_0[1,:]))
angle_J1_JoMM = missings(Float64, length(J_0[1,:]))
angle_J1_JiMM = missings(Float64, length(J_0[1,:]))
angle_J1_JsMM = missings(Float64, length(J_0[1,:]))
angle_Jim_Jom = missings(Float64, length(J_0[1,:]))
angle_JiMM_JoMM= missings(Float64, length(J_0[1,:]))
for i in 1:length(J_0[1,:])
        if count(ismissing, Jom[:,i]) == 0
                angle_J0_Jom[i] = 180/π * angle(J_0[:,i], J_0[:,i] .+ convert(Array{Float64,1},Jom[:,i]) ) 
                angle_J1_Jom[i] = 180/π * angle(J_1[:,i], J_0[:,i] .+ convert(Array{Float64,1},Jom[:,i]) ) 
                if count(ismissing, Jim[:,i]) == 0
                        angle_Jim_Jom[i] = 180/π * angle(convert(Array{Float64,1},Jim[:,i]), convert(Array{Float64,1},Jom[:,i]) ) 
                end
        end
        if count(ismissing, JoMM[:,i]) == 0
                angle_J0_JoMM[i] = 180/π * angle(J_0[:,i], J_0[:,i] .+ convert(Array{Float64,1},JoMM[:,i]) ) 
                angle_J1_JoMM[i] = 180/π * angle(J_1[:,i], J_0[:,i] .+ convert(Array{Float64,1},JoMM[:,i]) ) 
                if count(ismissing, JiMM[:,i]) == 0
                        angle_JiMM_JoMM[i] = 180/π * angle(convert(Array{Float64,1},JiMM[:,i]), convert(Array{Float64,1},JoMM[:,i]) ) 
                end
        end
        if count(ismissing, Jsm[:,i]) == 0
                angle_J0_Jsm[i] = 180/π * angle(J_0[:,i], J_0[:,i] .+ convert(Array{Float64,1},Jsm[:,i]) ) 
                angle_J1_Jsm[i] = 180/π * angle(J_1[:,i], J_0[:,i] .+ convert(Array{Float64,1},Jsm[:,i]) ) 
        end
        if count(ismissing, JsMM[:,i]) == 0
                angle_J0_JsMM[i] = 180/π * angle(J_0[:,i], J_0[:,i] .+ convert(Array{Float64,1},JsMM[:,i]) ) 
                angle_J1_JsMM[i] = 180/π * angle(J_1[:,i], J_0[:,i] .+ convert(Array{Float64,1},JsMM[:,i]) ) 
        end
        if count(ismissing, Jim[:,i]) == 0
                angle_J0_Jim[i] = 180/π * angle(J_0[:,i], J_0[:,i] .+ convert(Array{Float64,1},Jim[:,i]) ) 
                angle_J1_Jim[i] = 180/π * angle(J_1[:,i], J_0[:,i] .+ convert(Array{Float64,1},Jim[:,i]) ) 
        end
        if count(ismissing, JiMM[:,i]) == 0
                angle_J0_JiMM[i] = 180/π * angle(J_0[:,i], J_0[:,i] .+ convert(Array{Float64,1},JiMM[:,i]) ) 
                angle_J1_JiMM[i] = 180/π * angle(J_1[:,i], J_0[:,i] .+ convert(Array{Float64,1},JiMM[:,i]) ) 
        end
end

# hexbin: mean() statt count
outfile     = joinpath(outputdir,"Jim_vs_Jom-imm-03_hexmean.png")
pointsize = 5
fig, ax = subplots()
colr    = angle_Jim_Jom
nomissings = idxclude(1:length(colr), findcs(colr, eq=missing))
colr    = colr[nomissings]
y       = pyplottable(Jim[:,nomissings])
x       = pyplottable(Jom[:,nomissings])
plot_hexbins( log10.(x), log10.(y), 
    outfile=outfile,
    #ymin=0, ymax=180, xmin=-1, xmax=1, 
    title="$(latexstring("\$\\Delta\$"))t $(latexstring("\$\\sim\$")) 0.3 Gyr",
    scale=15, gridres=(30,13), xscale="linear", yscale="linear", 
    xlabel="log( J$(latexstring("\$_\\textrm{orbital,Mergers}\$")) [ M$(latexstring("\$_\\odot\$")) kpc km/s ] )", 
    ylabel="log( J$(latexstring("\$_\\textrm{intrinsic,Mergers}\$")) [ M$(latexstring("\$_\\odot\$")) kpc km/s ] )", 
    clabel= "Mean Alignment Angle [°]",
    grid=true,
    xmin=10, xmax=16, ymin=7, ymax=15,
    lognorm=false, colorlimits=(0,180),
    cmap="plasma_r", cfunc=mean, colormod=colr,
    binned_median = "x", plot_bval=" ", plot_line=[10,15,10,15],
    calc_norm=true,    # @pyplottable
    )
#
# hexbin: mean() statt count
outfile     = joinpath(outputdir,"JiMM_vs_JoMM-imm-03_hexmean.png")
pointsize = 5
fig, ax = subplots()
colr    = angle_JiMM_JoMM
nomissings = idxclude(1:length(colr), findcs(colr, eq=missing))
colr    = colr[nomissings]
y       = pyplottable(JiMM[:,nomissings])
x       = pyplottable(JoMM[:,nomissings])
plot_hexbins( log10.(x), log10.(y), 
    outfile=outfile,
    #ymin=0, ymax=180, xmin=-1, xmax=1, 
    title="$(latexstring("\$\\Delta\$"))t $(latexstring("\$\\sim\$")) 0.3 Gyr",
    scale=15, gridres=(30,13), xscale="linear", yscale="linear", 
    xlabel="log( J$(latexstring("\$_\\textrm{orbital,MM}\$")) [ M$(latexstring("\$_\\odot\$")) kpc km/s ] )", 
    ylabel="log( J$(latexstring("\$_\\textrm{intrinsic,MM}\$")) [ M$(latexstring("\$_\\odot\$")) kpc km/s ] )", 
    clabel= "Mean Alignment Angle [°]",
    grid=true,
    xmin=10, xmax=16, ymin=7, ymax=15,
    lognorm=false, colorlimits=(0,180),
    cmap="plasma_r", cfunc=mean, colormod=colr,
    binned_median = "x", plot_bval=" ", plot_line=[10,15,10,15],
    calc_norm=true,    # @pyplottable
    )
#


outfile     = joinpath(outputdir,"Jim_vs_Jom-imm-03.png")
scale = 0.5
PyPlot.matplotlib[:rc]("text", usetex=true)
PyPlot.matplotlib[:rc]("font", family="serif", serif=["Computer Modern"])
pointsize = 5
fig, ax = subplots()
colr    = angle_Jim_Jom
nomissings = idxclude(1:length(colr), findcs(colr, eq=missing))
colr    = colr[nomissings]
order   = sortperm(colr)[length(colr):-1:1]
colr    = colr[order]
#y       = pyplottable(delJ[:,order] .- Jom[:,order])
y       = pyplottable(Jim[:,nomissings])
x       = pyplottable(Jom[:,nomissings])
maxcolr = maximum(colr)
#slct    = idxmatch(slct2, n_switches)
p = ax.scatter( log10.(x[order]), log10.(y[order]),
        zorder=2   , s= pointsize    , alpha=0.5,#colr ./ maxcolr, (pointsize/maxcolr) .* colr
        cmap="plasma_r", c=colr) # brg_r  
ax.plot([7,16], [7,16], "-",color="black", linewidth=0.5)
colbar  = fig.colorbar(p, ax=ax)#, ticks=[minimum(as03["BVAL_0"    ]), b_ell, b_disk, maximum(as03["BVAL_0"    ])])
#colbar.ax.set_yticklabels(["$(round(minimum(as03["BVAL_0"    ]), digits=1))", "ETG", "LTG", "$(round(maximum(as03["BVAL_0"    ]), digits=1))"])
colbar.ax.tick_params(axis="y")
colbar.ax.set_ylabel("Alignment Angle [°]")
colbar.ax.yaxis.set_label_position("left")
#ax.legend(loc="lower right", frameon=true, borderpad=1, handlelength=1.8)
ax.set_ylim(7,15)
ax.set_xlim(10,16)
ax.set_title("immediate read, all mergers, 03 Gyr time step")
#ax.set_aspect("equal")
ax.set_ylabel("log( J$(latexstring("\$_\\textrm{intrinsic,Mergers}\$")) [ M$(latexstring("\$_\\odot\$")) kpc km/s ] )")
ax.set_xlabel("log( J$(latexstring("\$_\\textrm{orbital,Mergers}\$")) [ M$(latexstring("\$_\\odot\$")) kpc km/s ] )")
ax.grid()
#ax.autoscale()
fig.set_size_inches(11*scale, 9*scale)
fig.tight_layout()
fig.savefig(outfile, bbox_inches="tight", pad_inches=.1)
close("all")
PyPlot.matplotlib[:rc]("text", usetex=false)


outfile     = joinpath(outputdir,"JiMM_vs_JoMM-imm-03.png")
scale = 0.5
PyPlot.matplotlib[:rc]("text", usetex=true)
PyPlot.matplotlib[:rc]("font", family="serif", serif=["Computer Modern"])
pointsize = 5
fig, ax = subplots()
colr    = angle_JiMM_JoMM
nomissings = idxclude(1:length(colr), findcs(colr, eq=missing))
colr    = colr[nomissings]
order   = sortperm(colr)[length(colr):-1:1]
colr    = colr[order]
#y       = pyplottable(delJ[:,order] .- Jom[:,order])
y       = pyplottable(JiMM[:,nomissings])
x       = pyplottable(JoMM[:,nomissings])
maxcolr = maximum(colr)
#slct    = idxmatch(slct2, n_switches)
p = ax.scatter( log10.(x[order]), log10.(y[order]),
        zorder=2   , s= pointsize    , alpha=0.5,#colr ./ maxcolr, (pointsize/maxcolr) .* colr
        cmap="plasma_r", c=colr) # brg_r  
ax.plot([7,16], [7,16], "-",color="black", linewidth=0.5)
colbar  = fig.colorbar(p, ax=ax)#, ticks=[minimum(as03["BVAL_0"    ]), b_ell, b_disk, maximum(as03["BVAL_0"    ])])
#colbar.ax.set_yticklabels(["$(round(minimum(as03["BVAL_0"    ]), digits=1))", "ETG", "LTG", "$(round(maximum(as03["BVAL_0"    ]), digits=1))"])
colbar.ax.tick_params(axis="y")
colbar.ax.set_ylabel("Alignment Angle [°]")
colbar.ax.yaxis.set_label_position("left")
#ax.legend(loc="lower right", frameon=true, borderpad=1, handlelength=1.8)
ax.set_ylim(7,15)
ax.set_xlim(10,16)
ax.set_title("immediate read, Most Massive, 03 Gyr time step")
#ax.set_aspect("equal")
ax.set_ylabel("log( J$(latexstring("\$_\\textrm{intrinsic,MM}\$")) [ M$(latexstring("\$_\\odot\$")) kpc km/s ] )")
ax.set_xlabel("log( J$(latexstring("\$_\\textrm{orbital,MM}\$")) [ M$(latexstring("\$_\\odot\$")) kpc km/s ] )")
ax.grid()
#ax.autoscale()
fig.set_size_inches(11*scale, 9*scale)
fig.tight_layout()
fig.savefig(outfile, bbox_inches="tight", pad_inches=.1)
close("all")
PyPlot.matplotlib[:rc]("text", usetex=false)


 


plt_scale = "log"
xmin = 1e-1
ymin = xmin
xmax = 180
ymax = xmax
outfile     = joinpath(outputdir,"angJ0Jsm_v_angJdelta_$(plt_scale)_vs_bval0-imm-03.png")
scale = 0.5
PyPlot.matplotlib[:rc]("text", usetex=true)
PyPlot.matplotlib[:rc]("font", family="serif", serif=["Computer Modern"])
pointsize = 2
fig, ax = subplots()
#y       = pyplottable(delJ[:,order] .- Jom[:,order])
colr    = as03["BVAL_0"][idxmatch(mmi_as03["main"], idx_base)] .- as03["ΔBVAL_0"][idxmatch(mmi_as03["main"], idx_base)] # log10.(angle_J1_Jsm ./ x)
nomissings = idxclude(1:length(colr), findcs(colr, eq=missing))
colr    = colr[nomissings]
order   = sortperm(colr)
colr    = colr[order]
x       = pyplottable(as03["ϕ_flip"][idxmatch(mmi_as03["main"], idx_base)])[nomissings]
y       = pyplottable(angle_J0_Jsm)[nomissings]
maxcolr = maximum(colr)
#slct    = idxmatch(slct2, n_switches)
ax.plot([xmin,xmax], [ymin,ymax], "-",color="black", linewidth=0.5)
p = ax.scatter( x[order], y[order],
        zorder=2   , s= pointsize    , alpha=0.5,#colr ./ maxcolr, (pointsize/maxcolr) .* colr
        cmap=ColorMap(get_bvalue_cmap().colors), c=colr, vmin=-6.1,vmax=-3.9) # brg_r  , vmin=-6.1,vmax=-3.9
#ax.plot([11,16], [11,16], "-",color="black", linewidth=0.5)
colbar  = fig.colorbar(p, ax=ax)#, ticks=[minimum(as03["BVAL_0"    ]), b_ell, b_disk, maximum(as03["BVAL_0"    ])])
#colbar.ax.set_yticklabels(["$(round(minimum(as03["BVAL_0"    ]), digits=1))", "ETG", "LTG", "$(round(maximum(as03["BVAL_0"    ]), digits=1))"])
colbar.ax.tick_params(axis="y")
colbar.ax.set_ylabel("b-value at Flip Start")
colbar.ax.yaxis.set_label_position("left")
#ax.legend(loc="lower right", frameon=true, borderpad=1, handlelength=1.8)
#ax.set_ylim(extrema(y))
ax.set_title("total, immediate, all mergers, 03 Gyr")
#ax.set_aspect("equal")
ax.set_ylabel("Flip by Merger Contribution [°]")
ax.set_xlabel("Flip Angle [°]")
ax.set_ylim([ymin,ymax])
ax.set_xlim([xmin,xmax])
ax.set_xscale(plt_scale)
ax.set_yscale(plt_scale)
ax.grid()
#ax.autoscale()
fig.set_size_inches(11*scale, 9*scale)
fig.tight_layout()
fig.savefig(outfile, bbox_inches="tight", pad_inches=.1)
close("all")
PyPlot.matplotlib[:rc]("text", usetex=false)


 


plt_scale = "linear"
xmin = 0
ymin = xmin
xmax = 30
ymax = ymax
outfile     = joinpath(outputdir,"angJ0Jsm_v_angJdelta_$(plt_scale)_vs_bval0-imm-03.png")
scale = 0.5
PyPlot.matplotlib[:rc]("text", usetex=true)
PyPlot.matplotlib[:rc]("font", family="serif", serif=["Computer Modern"])
pointsize = 2
fig, ax = subplots()
#y       = pyplottable(delJ[:,order] .- Jom[:,order])
colr    = as03["BVAL_0"][idxmatch(mmi_as03["main"], idx_base)] .- as03["ΔBVAL_0"][idxmatch(mmi_as03["main"], idx_base)] # log10.(angle_J1_Jsm ./ x)
nomissings = idxclude(1:length(colr), findcs(colr, eq=missing))
colr    = colr[nomissings]
order   = sortperm(colr)
colr    = colr[order]
x       = pyplottable(as03["ϕ_flip"][idxmatch(mmi_as03["main"], idx_base)])[nomissings]
y       = pyplottable(angle_J0_Jsm)[nomissings]
maxcolr = maximum(colr)
#slct    = idxmatch(slct2, n_switches)
ax.plot([xmin,xmax], [ymin,ymax], "-",color="black", linewidth=0.5)
p = ax.scatter( x[order], y[order],
        zorder=2   , s= pointsize    , alpha=0.5,#colr ./ maxcolr, (pointsize/maxcolr) .* colr
        cmap=ColorMap(get_bvalue_cmap().colors), c=colr, vmin=-6.1,vmax=-3.9) # brg_r  , vmin=-6.1,vmax=-3.9
#ax.plot([11,16], [11,16], "-",color="black", linewidth=0.5)
colbar  = fig.colorbar(p, ax=ax)#, ticks=[minimum(as03["BVAL_0"    ]), b_ell, b_disk, maximum(as03["BVAL_0"    ])])
#colbar.ax.set_yticklabels(["$(round(minimum(as03["BVAL_0"    ]), digits=1))", "ETG", "LTG", "$(round(maximum(as03["BVAL_0"    ]), digits=1))"])
colbar.ax.tick_params(axis="y")
colbar.ax.set_ylabel("b-value at Flip Start")
colbar.ax.yaxis.set_label_position("left")
#ax.legend(loc="lower right", frameon=true, borderpad=1, handlelength=1.8)
#ax.set_ylim(extrema(y))
ax.set_title("total, immediate, all mergers, 03 Gyr")
#ax.set_aspect("equal")
ax.set_ylabel("Flip by Merger Contribution [°]")
ax.set_xlabel("Flip Angle [°]")
ax.set_ylim([ymin,ymax])
ax.set_xlim([xmin,xmax])
ax.set_xscale(plt_scale)
ax.set_yscale(plt_scale)
ax.grid()
#ax.autoscale()
fig.set_size_inches(11*scale, 9*scale)
fig.tight_layout()
fig.savefig(outfile, bbox_inches="tight", pad_inches=.1)
close("all")
PyPlot.matplotlib[:rc]("text", usetex=false)


plt_scale = "linear"
xmin = 0
ymin = xmin
xmax = 30
ymax = ymax
outfile     = joinpath(outputdir,"angJ0JsMM_v_angJdelta_$(plt_scale)_vs_bval0-imm-03.png")
scale = 0.5
PyPlot.matplotlib[:rc]("text", usetex=true)
PyPlot.matplotlib[:rc]("font", family="serif", serif=["Computer Modern"])
pointsize = 2
fig, ax = subplots()
#y       = pyplottable(delJ[:,order] .- Jom[:,order])
colr    = as03["BVAL_0"][idxmatch(mmi_as03["main"], idx_base)] .- as03["ΔBVAL_0"][idxmatch(mmi_as03["main"], idx_base)] # log10.(angle_J1_Jsm ./ x)
nomissings = idxclude(1:length(colr), findcs(colr, eq=missing))
colr    = colr[nomissings]
order   = sortperm(colr)
colr    = colr[order]
x       = pyplottable(as03["ϕ_flip"][idxmatch(mmi_as03["main"], idx_base)])[nomissings]
y       = pyplottable(angle_J0_JsMM)[nomissings]
maxcolr = maximum(colr)
#slct    = idxmatch(slct2, n_switches)
ax.plot([xmin,xmax], [ymin,ymax], "-",color="black", linewidth=0.5)
p = ax.scatter( x[order], y[order],
        zorder=2   , s= pointsize    , alpha=0.5,#colr ./ maxcolr, (pointsize/maxcolr) .* colr
        cmap=ColorMap(get_bvalue_cmap().colors), c=colr, vmin=-6.1,vmax=-3.9) # brg_r  , vmin=-6.1,vmax=-3.9
#ax.plot([11,16], [11,16], "-",color="black", linewidth=0.5)
colbar  = fig.colorbar(p, ax=ax)#, ticks=[minimum(as03["BVAL_0"    ]), b_ell, b_disk, maximum(as03["BVAL_0"    ])])
#colbar.ax.set_yticklabels(["$(round(minimum(as03["BVAL_0"    ]), digits=1))", "ETG", "LTG", "$(round(maximum(as03["BVAL_0"    ]), digits=1))"])
colbar.ax.tick_params(axis="y")
colbar.ax.set_ylabel("b-value at Flip Start")
colbar.ax.yaxis.set_label_position("left")
#ax.legend(loc="lower right", frameon=true, borderpad=1, handlelength=1.8)
ax.set_ylim([ymin,ymax])
ax.set_xlim([xmin,xmax])
ax.set_xscale(plt_scale)
ax.set_yscale(plt_scale)
ax.set_title("total, immediate, most massive, 03 Gyr")
#ax.set_aspect("equal")
ax.set_ylabel("Flip by Merger Contribution [°]")
ax.set_xlabel("Flip Angle [°]")
ax.grid()
#ax.autoscale()
fig.set_size_inches(11*scale, 9*scale)
fig.tight_layout()
fig.savefig(outfile, bbox_inches="tight", pad_inches=.1)
close("all")
PyPlot.matplotlib[:rc]("text", usetex=false)


 

##################################################################
# now with color by J
plt_scale = "linear"
xmin = 0
ymin = xmin
xmax = 30
ymax = xmax
outfile     = joinpath(outputdir,"angJ0Jsm_v_angJdelta_$(plt_scale)_vs_J_pre-imm-03.png")
scale = 0.5
PyPlot.matplotlib[:rc]("text", usetex=true)
PyPlot.matplotlib[:rc]("font", family="serif", serif=["Computer Modern"])
pointsize = 2
fig, ax = subplots()
#y       = pyplottable(delJ[:,order] .- Jom[:,order])
colr    = pyplottable(as03["J_main"][:,idxmatch(mmi_as03["main"], idx_base)] .- as03["ΔJ_main"][:,idxmatch(mmi_as03["main"], idx_base)]) # log10.(angle_J1_Jsm ./ x)
nomissings = idxclude(1:length(colr), findcs(colr, eq=missing))
colr    = colr[nomissings]
order   = sortperm(colr)
colr    = colr[order]
x       = pyplottable(as03["ϕ_flip"][idxmatch(mmi_as03["main"], idx_base)])[nomissings]
y       = pyplottable(angle_J0_Jsm)[nomissings]
maxcolr = maximum(colr)
#slct    = idxmatch(slct2, n_switches)
ax.plot([xmin,xmax], [ymin,ymax], "-",color="black", linewidth=0.5)
p = ax.scatter( x[order], y[order],
        zorder=2   , s= pointsize    , alpha=0.5,#colr ./ maxcolr, (pointsize/maxcolr) .* colr
        cmap="plasma_r", c=colr) # brg_r  , vmin=-6.1,vmax=-3.9
#ax.plot([11,16], [11,16], "-",color="black", linewidth=0.5)
colbar  = fig.colorbar(p, ax=ax)#, ticks=[minimum(as03["BVAL_0"    ]), b_ell, b_disk, maximum(as03["BVAL_0"    ])])
#colbar.ax.set_yticklabels(["$(round(minimum(as03["BVAL_0"    ]), digits=1))", "ETG", "LTG", "$(round(maximum(as03["BVAL_0"    ]), digits=1))"])
colbar.ax.tick_params(axis="y")
colbar.ax.set_ylabel("log( J$(latexstring("\$_*\$")) [ M$(latexstring("\$_\\odot\$")) kpc km/s ] ) at Flip Start")
colbar.ax.yaxis.set_label_position("left")
#ax.legend(loc="lower right", frameon=true, borderpad=1, handlelength=1.8)
#ax.set_ylim(extrema(y))
ax.set_title("total, immediate, all mergers, 03 Gyr")
#ax.set_aspect("equal")
ax.set_ylabel("Flip by Merger Contribution [°]")
ax.set_xlabel("Flip Angle [°]")
ax.set_ylim([ymin,ymax])
ax.set_xlim([xmin,xmax])
ax.set_xscale(plt_scale)
ax.set_yscale(plt_scale)
ax.grid()
#ax.autoscale()
fig.set_size_inches(11*scale, 9*scale)
fig.tight_layout()
fig.savefig(outfile, bbox_inches="tight", pad_inches=.1)
close("all")
PyPlot.matplotlib[:rc]("text", usetex=false)


plt_scale = "linear"
xmin = 0
ymin = xmin
xmax = 30
ymax = ymax
outfile     = joinpath(outputdir,"angJ0JsMM_v_angJdelta_$(plt_scale)_vs_J_pre-imm-03.png")
scale = 0.5
PyPlot.matplotlib[:rc]("text", usetex=true)
PyPlot.matplotlib[:rc]("font", family="serif", serif=["Computer Modern"])
pointsize = 2
fig, ax = subplots()
#y       = pyplottable(delJ[:,order] .- Jom[:,order])
colr    = pyplottable(as03["J_main"][:,idxmatch(mmi_as03["main"], idx_base)] .- as03["ΔJ_main"][:,idxmatch(mmi_as03["main"], idx_base)]) # log10.(angle_J1_Jsm ./ x)
nomissings = idxclude(1:length(colr), findcs(colr, eq=missing))
colr    = colr[nomissings]
order   = sortperm(colr)
colr    = colr[order]
x       = pyplottable(as03["ϕ_flip"][idxmatch(mmi_as03["main"], idx_base)])[nomissings]
y       = pyplottable(angle_J0_JsMM)[nomissings]
maxcolr = maximum(colr)
#slct    = idxmatch(slct2, n_switches)
ax.plot([xmin,xmax], [ymin,ymax], "-",color="black", linewidth=0.5)
p = ax.scatter( x[order], y[order],
        zorder=2   , s= pointsize    , alpha=0.5,#colr ./ maxcolr, (pointsize/maxcolr) .* colr
        cmap="plasma_r", c=colr) # brg_r  , vmin=-6.1,vmax=-3.9
#ax.plot([11,16], [11,16], "-",color="black", linewidth=0.5)
colbar  = fig.colorbar(p, ax=ax)#, ticks=[minimum(as03["BVAL_0"    ]), b_ell, b_disk, maximum(as03["BVAL_0"    ])])
#colbar.ax.set_yticklabels(["$(round(minimum(as03["BVAL_0"    ]), digits=1))", "ETG", "LTG", "$(round(maximum(as03["BVAL_0"    ]), digits=1))"])
colbar.ax.tick_params(axis="y")
colbar.ax.set_ylabel("log( J$(latexstring("\$_*\$")) [ M$(latexstring("\$_\\odot\$")) kpc km/s ] ) at Flip Start")
colbar.ax.yaxis.set_label_position("left")
#ax.legend(loc="lower right", frameon=true, borderpad=1, handlelength=1.8)
ax.set_ylim([ymin,ymax])
ax.set_xlim([xmin,xmax])
ax.set_xscale(plt_scale)
ax.set_yscale(plt_scale)
ax.set_title("total, immediate, most massive, 03 Gyr")
#ax.set_aspect("equal")
ax.set_ylabel("Flip by Merger Contribution [°]")
ax.set_xlabel("Flip Angle [°]")
ax.grid()
#ax.autoscale()
fig.set_size_inches(11*scale, 9*scale)
fig.tight_layout()
fig.savefig(outfile, bbox_inches="tight", pad_inches=.1)
close("all")
PyPlot.matplotlib[:rc]("text", usetex=false)


 

##################################################################
# now with color by j
plt_scale = "linear"
xmin = 0
ymin = xmin
xmax = 30
ymax = xmax
outfile     = joinpath(outputdir,"angJ0Jsm_v_angJdelta_$(plt_scale)_vs_jpre-imm-03.png")
scale = 0.5
PyPlot.matplotlib[:rc]("text", usetex=true)
PyPlot.matplotlib[:rc]("font", family="serif", serif=["Computer Modern"])
pointsize = 2
fig, ax = subplots()
#y       = pyplottable(delJ[:,order] .- Jom[:,order])
colr    = pyplottable(as03["j_main"][:,idxmatch(mmi_as03["main"], idx_base)] .- as03["Δj_main"][:,idxmatch(mmi_as03["main"], idx_base)]) # log10.(angle_J1_Jsm ./ x)
nomissings = idxclude(1:length(colr), findcs(colr, eq=missing))
colr    = colr[nomissings]
order   = sortperm(colr)
colr    = colr[order]
x       = pyplottable(as03["ϕ_flip"][idxmatch(mmi_as03["main"], idx_base)])[nomissings]
y       = pyplottable(angle_J0_Jsm)[nomissings]
maxcolr = maximum(colr)
#slct    = idxmatch(slct2, n_switches)
ax.plot([xmin,xmax], [ymin,ymax], "-",color="black", linewidth=0.5)
p = ax.scatter( x[order], y[order],
        zorder=2   , s= pointsize    , alpha=0.5,#colr ./ maxcolr, (pointsize/maxcolr) .* colr
        cmap="plasma_r", c=colr) # brg_r  , vmin=-6.1,vmax=-3.9
#ax.plot([11,16], [11,16], "-",color="black", linewidth=0.5)
colbar  = fig.colorbar(p, ax=ax)#, ticks=[minimum(as03["BVAL_0"    ]), b_ell, b_disk, maximum(as03["BVAL_0"    ])])
#colbar.ax.set_yticklabels(["$(round(minimum(as03["BVAL_0"    ]), digits=1))", "ETG", "LTG", "$(round(maximum(as03["BVAL_0"    ]), digits=1))"])
colbar.ax.tick_params(axis="y")
colbar.ax.set_ylabel("log( j$(latexstring("\$_*\$")) [ kpc km/s ] ) at Flip Start")
colbar.ax.yaxis.set_label_position("left")
#ax.legend(loc="lower right", frameon=true, borderpad=1, handlelength=1.8)
#ax.set_ylim(extrema(y))
ax.set_title("total, immediate, all mergers, 03 Gyr")
#ax.set_aspect("equal")
ax.set_ylabel("Flip by Merger Contribution [°]")
ax.set_xlabel("Flip Angle [°]")
ax.set_ylim([ymin,ymax])
ax.set_xlim([xmin,xmax])
ax.set_xscale(plt_scale)
ax.set_yscale(plt_scale)
ax.grid()
#ax.autoscale()
fig.set_size_inches(11*scale, 9*scale)
fig.tight_layout()
fig.savefig(outfile, bbox_inches="tight", pad_inches=.1)
close("all")
PyPlot.matplotlib[:rc]("text", usetex=false)


plt_scale = "linear"
xmin = 0
ymin = xmin
xmax = 30
ymax = ymax
outfile     = joinpath(outputdir,"angJ0JsMM_v_angJdelta_$(plt_scale)_vs_jpre-imm-03.png")
scale = 0.5
PyPlot.matplotlib[:rc]("text", usetex=true)
PyPlot.matplotlib[:rc]("font", family="serif", serif=["Computer Modern"])
pointsize = 2
fig, ax = subplots()
#y       = pyplottable(delJ[:,order] .- Jom[:,order])
colr    = pyplottable(as03["j_main"][:,idxmatch(mmi_as03["main"], idx_base)] .- as03["Δj_main"][:,idxmatch(mmi_as03["main"], idx_base)]) # log10.(angle_J1_Jsm ./ x)
nomissings = idxclude(1:length(colr), findcs(colr, eq=missing))
colr    = colr[nomissings]
order   = sortperm(colr)
colr    = colr[order]
x       = pyplottable(as03["ϕ_flip"][idxmatch(mmi_as03["main"], idx_base)])[nomissings]
y       = pyplottable(angle_J0_JsMM)[nomissings]
maxcolr = maximum(colr)
#slct    = idxmatch(slct2, n_switches)
ax.plot([xmin,xmax], [ymin,ymax], "-",color="black", linewidth=0.5)
p = ax.scatter( x[order], y[order],
        zorder=2   , s= pointsize    , alpha=0.5,#colr ./ maxcolr, (pointsize/maxcolr) .* colr
        cmap="plasma_r", c=colr) # brg_r  , vmin=-6.1,vmax=-3.9
#ax.plot([11,16], [11,16], "-",color="black", linewidth=0.5)
colbar  = fig.colorbar(p, ax=ax)#, ticks=[minimum(as03["BVAL_0"    ]), b_ell, b_disk, maximum(as03["BVAL_0"    ])])
#colbar.ax.set_yticklabels(["$(round(minimum(as03["BVAL_0"    ]), digits=1))", "ETG", "LTG", "$(round(maximum(as03["BVAL_0"    ]), digits=1))"])
colbar.ax.tick_params(axis="y")
colbar.ax.set_ylabel("log( j$(latexstring("\$_*\$")) [ kpc km/s ] ) at Flip Start")
colbar.ax.yaxis.set_label_position("left")
#ax.legend(loc="lower right", frameon=true, borderpad=1, handlelength=1.8)
ax.set_ylim([ymin,ymax])
ax.set_xlim([xmin,xmax])
ax.set_xscale(plt_scale)
ax.set_yscale(plt_scale)
ax.set_title("total, immediate, most massive, 03 Gyr")
#ax.set_aspect("equal")
ax.set_ylabel("Flip by Merger Contribution [°]")
ax.set_xlabel("Flip Angle [°]")
ax.grid()
#ax.autoscale()
fig.set_size_inches(11*scale, 9*scale)
fig.tight_layout()
fig.savefig(outfile, bbox_inches="tight", pad_inches=.1)
close("all")
PyPlot.matplotlib[:rc]("text", usetex=false)






###################################################################
# predict delta J
xmin = " "#1e0
ymin = " "#xmin
xmax = " "#180
ymax = " "#xmax
outfile     = joinpath(outputdir,"J0Jsm_v_Jdelta-imm-03.png")
x       = pyplottable(as03["ΔJ_main"][:,idxmatch(mmi_as03["main"], idx_base)])
y       = pyplottable(Jsm)
#x       = pyplottable(as03["J_main"][:,idxmatch(mmi_as03["main"], idx_base)])
#y       = pyplottable(as03["J_main"][:,idxmatch(mmi_as03["main"], idx_base)] .- as03["ΔJ_main"][:,idxmatch(mmi_as03["main"], idx_base)] .+ Jsm)
plot_hexbins( log10.(x), log10.(y), 
    outfile=outfile,
    title="$(latexstring("\$\\Delta\$"))t $(latexstring("\$\\sim\$")) 0.3 Gyr",
    scale=15, gridres=(30,13), xscale="linear", yscale="linear", 
    xlabel="log( $(latexstring("\$\\Delta\$"))J [ M$(latexstring("\$_\\odot\$")) kpc km/s ] )", 
    ylabel="log( J$(latexstring("\$_\\textrm{Mergers}\$")) [ M$(latexstring("\$_\\odot\$")) kpc km/s ] )", 
    grid=true,
    xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax,
    lognorm=false,
    binned_median = "x", plot_bval=" ", plot_line=[9,18,9,18],
    calc_norm=true,    # @pyplottable
    )
#
plt_scale = "linear"
xmin = nothing
ymin = nothing
xmax = nothing
ymax = nothing
outfile     = joinpath(outputdir,"J0Jsm_v_Jdelta_vs_angle-imm-03.png")
scale = 0.5
PyPlot.matplotlib[:rc]("text", usetex=true)
PyPlot.matplotlib[:rc]("font", family="serif", serif=["Computer Modern"])
pointsize = 2
fig, ax = subplots()
#y       = pyplottable(delJ[:,order] .- Jom[:,order])
colr    = Vector{Float64}(undef, 0)
for i in 1:length(idxmatch(mmi_as03["main"], idx_base))
        colr    = vcat( colr, (180/π) .* angle(convert(Array{Float64,1},as03["ΔJ_main"][:,idxmatch(mmi_as03["main"], idx_base)][:,i]), convert(Array{Float64,1},Jsm[:,i]) ) )
end
order   = sortperm(colr)
colr    = colr[order]
maxcolr = maximum(colr)
#slct    = idxmatch(slct2, n_switches)
ax.plot([11,15], [11,15], "-",color="black", linewidth=0.5)
p = ax.scatter( log10.(x[order]), log10.(y[order]),
        zorder=2   , s= pointsize    , alpha=0.5,#colr ./ maxcolr, (pointsize/maxcolr) .* colr
        cmap="plasma_r", c=colr) # brg_r  , vmin=-6.1,vmax=-3.9
#ax.plot([11,16], [11,16], "-",color="black", linewidth=0.5)
colbar  = fig.colorbar(p, ax=ax)#, ticks=[minimum(as03["BVAL_0"    ]), b_ell, b_disk, maximum(as03["BVAL_0"    ])])
#colbar.ax.set_yticklabels(["$(round(minimum(as03["BVAL_0"    ]), digits=1))", "ETG", "LTG", "$(round(maximum(as03["BVAL_0"    ]), digits=1))"])
colbar.ax.tick_params(axis="y")
colbar.ax.set_ylabel("Alignment Angle [°]")
colbar.ax.yaxis.set_label_position("left")
#ax.legend(loc="lower right", frameon=true, borderpad=1, handlelength=1.8)
ax.set_ylim([ymin,ymax])
ax.set_xlim([xmin,xmax])
#ax.set_aspect("equal")
ax.set_ylabel("log( J$(latexstring("\$_\\textrm{Mergers}\$")) [ M$(latexstring("\$_\\odot\$")) kpc km/s ] )")
ax.set_xlabel("log( $(latexstring("\$\\Delta\$"))J [ M$(latexstring("\$_\\odot\$")) kpc km/s ] )")
ax.grid()
#ax.autoscale()
fig.set_size_inches(11*scale, 9*scale)
fig.tight_layout()
fig.savefig(outfile, bbox_inches="tight", pad_inches=.1)
close("all")
PyPlot.matplotlib[:rc]("text", usetex=false)



































#outfile     = joinpath(outputdir,"Jim_vs_Jom-imm-03.png")
#scale = 0.5
#PyPlot.matplotlib[:rc]("text", usetex=true)
#PyPlot.matplotlib[:rc]("font", family="serif", serif=["Computer Modern"])
#pointsize = 10
#fig, ax = subplots()
#colr    = as03["ϕ_flip"][idxmatch(mmi_as03["main"], idx_base)]
#order   = sortperm(colr)
#colr    = colr[order]
##y       = pyplottable(delJ[:,order] .- Jom[:,order])
#y       = pyplottable(Jim[:,order])
#x       = pyplottable(Jom[:,order])
#maxcolr = maximum(colr)
##slct    = idxmatch(slct2, n_switches)
#p = ax.scatter( log10.(x), log10.(y),
#        zorder=2   , s= pointsize    , alpha=0.7,#colr ./ maxcolr, (pointsize/maxcolr) .* colr
#        cmap="plasma_r", c=colr) # brg_r  
#ax.plot([10,16], [10,16], "-",color="black", linewidth=0.5)
#colbar  = fig.colorbar(p, ax=ax)#, ticks=[minimum(as03["BVAL_0"    ]), b_ell, b_disk, maximum(as03["BVAL_0"    ])])
##colbar.ax.set_yticklabels(["$(round(minimum(as03["BVAL_0"    ]), digits=1))", "ETG", "LTG", "$(round(maximum(as03["BVAL_0"    ]), digits=1))"])
#colbar.ax.tick_params(axis="y")
#colbar.ax.set_ylabel("Flip Angle [°]")
#colbar.ax.yaxis.set_label_position("left")
##ax.legend(loc="lower right", frameon=true, borderpad=1, handlelength=1.8)
##ax.set_ylim(extrema(y))
#ax.set_title("orb, immediate, sum, 0.3")
##ax.set_aspect("equal")
#ax.set_xlabel("log( $(latexstring("\$\\Delta\$"))J$(latexstring("\$_*\$")) [ M$(latexstring("\$_\\odot\$")) kpc km/s ] )")
#ax.set_ylabel("log( $(latexstring("\$\\Delta\$"))J$(latexstring("\$_*\$")) - $(latexstring("\$\\Delta\$"))J$(latexstring("\$_\\textrm{orbital,Mergers}\$")))")
#ax.grid()
##ax.autoscale()
#fig.set_size_inches(25*scale, 9*scale)
#fig.tight_layout()
#fig.savefig(outfile, bbox_inches="tight", pad_inches=.1)
#close("all")
#PyPlot.matplotlib[:rc]("text", usetex=false)



#outfile     = joinpath(outputdir,"JiMM_vs_JoMM-imm-03.png")
#plot_hexbins( log10.(pyplottable(JoMM)), log10.(pyplottable(JiMM)), 
#    outfile=outfile,
#    #ymin=0, ymax=180, xmin=-1, xmax=1, 
#    title="$(latexstring("\$\\Delta\$"))t $(latexstring("\$\\sim\$")) 0.3 Gyr, immediate",
#    scale=15, gridres=(30,13), xscale="linear", yscale="linear", 
#    xlabel="J$(latexstring("\$_\\textrm{orbital,MM}\$")) [ M$(latexstring("\$_\\odot\$")) kpc km/s ]", 
#    ylabel="J$(latexstring("\$_\\textrm{intrinsic,MM}\$")) [ M$(latexstring("\$_\\odot\$")) kpc km/s ]", grid=true,
#    lognorm=true, colorlimits=(0.8,nothing),
#    cmap="BuPu",
#    binned_median = "x", plot_bval=" ",
#    weights=" ",
#    calc_norm=true,    # @pyplottable
#    )
##



# sfc: to do
## angles:
### phi_flip vs ang J-delJ & J-delJ+Jm



#PyPlot.matplotlib[:rc]("text", usetex=true)
#PyPlot.matplotlib[:rc]("font", family="serif", serif=["Computer Modern"])
#fig, ax = subplots()
##ax.hist( angle_J0_delJ_Jom    , bins=50, label="Disks"    , rwidth=0.9    , alpha=0.9, color="mediumblue" )
#ax.scatter( as03["ϕ_flip"][idxmatch(mmi_as03["main"], idx_base)], angle_J0_delJ_Jom        , alpha=0.9, color="mediumblue" )
##ax.set_xlabel("Angle [] between $(latexstring("\$\\Delta\$"))J & J$(latexstring("\$_\\textrm{orbital,Mergers}\$"))")
#ax.set_xlabel("Angle")
##ax.set_xscale("linear")
##ax.set_xlim([0, 180])
#ax.legend(loc="upper right", frameon=true, borderpad=1, handlelength=1.8)
#ax.grid()
#fig.set_size_inches(16scale, 9scale)
#fig.tight_layout()
#fig.savefig(outfile, bbox_inches="tight", pad_inches=.1)
#close("all")
#PyPlot.matplotlib[:rc]("text", usetex=false)



#plot_hexbins( log10.(x), log10.(y), 
#    outfile=outfile,
#    #ymin=0, ymax=180, xmin=-1, xmax=1, 
#    title="$(latexstring("\$\\Delta\$"))t $(latexstring("\$\\sim\$")) 0.3 Gyr",
#    scale=15, gridres=(30,13), xscale="linear", yscale="linear", 
#    xlabel="J$(latexstring("\$_\\textrm{orbital,Mergers}\$")) [ M$(latexstring("\$_\\odot\$")) kpc km/s ]", 
#    ylabel="$(latexstring("\$\\Delta\$"))J [ M$(latexstring("\$_\\odot\$")) kpc km/s ]", grid=true,
#    lognorm=true, colorlimits=(0.8,nothing),
#    cmap="BuPu",
#    binned_median = "x", plot_bval=" ",
#    weights=" ",
#    calc_norm=true,    # @pyplottable
#    )



## histogram lookbacktimes
#outfile     = joinpath(outputdir,"hist_lbt.png")
#PyPlot.matplotlib[:rc]("text", usetex=true)
#PyPlot.matplotlib[:rc]("font", family="serif", serif=["Computer Modern"])
#scale = 0.5
#fig, ax = subplots()
#ax.hist( as03["merger map"][22,mmi_as03["most_massive"][i]]    , bins=50, label="Ellipticals"    , rwidth=0.9    , alpha=0.5, color="darkred" )
#ax.hist( as03["merger map"][13,mmi_as03["most_massive"][i]]    , bins=50, label="Disks"    , rwidth=0.9    , alpha=0.5, color="mediumblue" )
#ax.set_xlabel("Time Difference [Gyr]")you
##ax.set_xscale("linear")
##ax.set_xlim([10, nothing])
#ax.legend(loc="upper right", frameon=true, borderpad=1, handlelength=1.8)
#ax.grid()
#scale=0.5
#fig.set_size_inches(16scale, 9scale)
#fig.tight_layout()
#fig.savefig(outfile, bbox_inches="tight", pad_inches=.1)
#close("all")
#PyPlot.matplotlib[:rc]("text", usetex=false)


######################################
# endstate vmin=-6.1,vmax=-3.9
xmin = -6.1
ymin = " "
xmax = -3.9
ymax = " "
outnm     = joinpath(outputdir,"endstate_flip30-imm-03")
outfile = outnm*".png"
x       = pyplottable(endstate["bval_end"])
y       = pyplottable(endstate["nflips30"])
#x       = pyplottable(as03["J_main"][:,idxmatch(mmi_as03["main"], idx_base)])
#y       = pyplottable(as03["J_main"][:,idxmatch(mmi_as03["main"], idx_base)] .- as03["ΔJ_main"][:,idxmatch(mmi_as03["main"], idx_base)] .+ Jsm)
plot_hexbins( x, y, 
    outfile=outfile,
    title="$(latexstring("\$\\Delta\$"))t $(latexstring("\$\\sim\$")) 0.3 Gyr",
    scale=15, gridres=(20,8), xscale="linear", yscale="linear", 
    xlabel="b-value at z = 0", 
    ylabel="N$(latexstring("\$_\\textrm{flips, 30°}\$"))", 
    grid=true,
    xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax,
    lognorm=true,
    binned_median = "x", plot_bval="x",
    calc_norm=true,    # @pyplottable
    )
#

outnm     = joinpath(outputdir,"endstate_flip90-imm-03")
outfile = outnm*".png"
x       = pyplottable(endstate["bval_end"])
y       = pyplottable(endstate["nflips90"])
#x       = pyplottable(as03["J_main"][:,idxmatch(mmi_as03["main"], idx_base)])
#y       = pyplottable(as03["J_main"][:,idxmatch(mmi_as03["main"], idx_base)] .- as03["ΔJ_main"][:,idxmatch(mmi_as03["main"], idx_base)] .+ Jsm)
plot_hexbins( x, y, 
    outfile=outfile,
    title="$(latexstring("\$\\Delta\$"))t $(latexstring("\$\\sim\$")) 0.3 Gyr",
    scale=15, gridres=(20,8), xscale="linear", yscale="linear", 
    xlabel="b-value at z = 0", 
    ylabel="N$(latexstring("\$_\\textrm{flips, 90°}\$"))", 
    grid=true,
    xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax,
    lognorm=true,
    binned_median = "x", plot_bval="x",
    calc_norm=true,    # @pyplottable
    )
#


