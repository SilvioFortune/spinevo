
#include("/home/moon/sfortune/spinevo/pkg/meta.jl")
#
#outputdir = "/home/moon/sfortune/spinevo/plots/friday"
#
#include("/home/moon/sfortune/spinevo/data/load_data.jl")


# 14, 15:17, 18:20, Immediate, imm
# 23, 24:26, 27:29, Earlier, earl
# 03, 0.3 Gyr, Instant
# 09, 1 Gyr, Long-term

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




























###################################################################
# predict delta J
xmin = 11.5#" "#1e0
ymin = 10#" "#xmin
xmax = 15#" "#180
ymax = 16#" "#xmax
outnm     = joinpath(outputdir,"J0Jsm_v_Jdelta-imm-03")
outfile = outnm*"hxcount.png"
x       = pyplottable(as03["ΔJ_main"][:,idxmatch(mmi_as03["main"], idx_base)])
y       = pyplottable(Jsm)
#x       = pyplottable(as03["J_main"][:,idxmatch(mmi_as03["main"], idx_base)])
#y       = pyplottable(as03["J_main"][:,idxmatch(mmi_as03["main"], idx_base)] .- as03["ΔJ_main"][:,idxmatch(mmi_as03["main"], idx_base)] .+ Jsm)
plot_hexbins( log10.(x), log10.(y), 
    outfile=outfile,
    title="$(latexstring("\$\\Delta\$"))t $(latexstring("\$\\sim\$")) 0.3 Gyr, Immediate Read",
    scale=15, gridres=(20,8), xscale="linear", yscale="linear", 
    xlabel="log( $(latexstring("\$\\Delta\$"))J [ M$(latexstring("\$_\\odot\$")) kpc km/s ] )", 
    ylabel="log( J$(latexstring("\$_\\textrm{Mergers}\$")) [ M$(latexstring("\$_\\odot\$")) kpc km/s ] )", 
    grid=true,
    xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax,
    lognorm=true,
    binned_median = "x", plot_bval=" ", plot_line=hcat([9,18],[9,18]),
    calc_norm=true,    # @pyplottable
    )
#


colr    = Vector{Float64}(undef, 0)
for i in 1:length(idxmatch(mmi_as03["main"], idx_base))
        colr    = vcat( colr, (180/π) .* angle(convert(Array{Float64,1},as03["ΔJ_main"][:,idxmatch(mmi_as03["main"], idx_base)][:,i]), convert(Array{Float64,1},Jsm[:,i]) ) )
end

outfile = outnm*"hxalign.png"
plot_hexbins( log10.(x), log10.(y), 
outfile=outfile,
title="$(latexstring("\$\\Delta\$"))t $(latexstring("\$\\sim\$")) 0.3 Gyr, Immediate Read",
scale=15, gridres=(30,13), xscale="linear", yscale="linear", 
xlabel="log( $(latexstring("\$\\Delta\$"))J [ M$(latexstring("\$_\\odot\$")) kpc km/s ] )", 
ylabel="log( J$(latexstring("\$_\\textrm{Mergers}\$")) [ M$(latexstring("\$_\\odot\$")) kpc km/s ] )", 
grid=true,
xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax,
lognorm=false,
binned_median = "x", plot_bval=" ", plot_line=hcat([9,18],[9,18]),
calc_norm=true,    # @pyplottable
cmap = "plasma_r", colormod = colr, colorlimits=(0,180), cfunc=mean,
)
#


order   = sortperm(colr)
colr    = colr[order]
maxcolr = maximum(colr)
#slct    = idxmatch(slct2, n_switches)
scale = 0.5
PyPlot.matplotlib[:rc]("text", usetex=true)
PyPlot.matplotlib[:rc]("font", family="serif", serif=["Computer Modern"])
pointsize = 2
fig, ax = subplots()
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
fig.savefig(outnm*"sctalign.png", bbox_inches="tight", pad_inches=.1)
close("all")
PyPlot.matplotlib[:rc]("text", usetex=false)


#colr = as03["BVAL_0"][idxmatch(mmi_as03["main"], idx_base)]
outfile = outnm*"hxbvalue.png"
plot_hexbins( log10.(x), log10.(y), 
outfile=outfile,
title="$(latexstring("\$\\Delta\$"))t $(latexstring("\$\\sim\$")) 0.3 Gyr, Immediate Read",
scale=15, gridres=(30,13), xscale="linear", yscale="linear", 
xlabel="log( $(latexstring("\$\\Delta\$"))J [ M$(latexstring("\$_\\odot\$")) kpc km/s ] )", 
ylabel="log( J$(latexstring("\$_\\textrm{Mergers}\$")) [ M$(latexstring("\$_\\odot\$")) kpc km/s ] )", 
grid=true,
xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax,
lognorm=false,
binned_median = "x", plot_bval=" ", plot_line=hcat([9,18],[9,18]),
calc_norm=true,    # @pyplottable
cmap = ColorMap(get_bvalue_cmap().colors), colormod = as03["BVAL_0"][idxmatch(mmi_as03["main"], idx_base)], colorlimits=(-6.1,-3.9), cfunc=mean,
)





###################################################################
# predict angle
xmin = 1e-2#" "#1e0
ymin = 1e-2#" "#xmin
xmax = 180#" "#180
ymax = 180#" "#xmax
#y       = pyplottable(delJ[:,order] .- Jom[:,order])
colr    = as03["BVAL_0"][idxmatch(mmi_as03["main"], idx_base)] .- as03["ΔBVAL_0"][idxmatch(mmi_as03["main"], idx_base)] # log10.(angle_J1_Jsm ./ x)
nomissings = idxclude(1:length(colr), findcs(colr, eq=missing))
colr    = colr[nomissings]
order   = sortperm(colr)
colr    = colr[order]
x       = pyplottable(as03["ϕ_flip"][idxmatch(mmi_as03["main"], idx_base)])[nomissings][order]
y       = pyplottable(angle_J0_Jsm)[nomissings][order]
outnm     = joinpath(outputdir,"angmerger_v_ang-imm-03")
outfile = outnm*"hxcount.png"
#x       = pyplottable(as03["J_main"][:,idxmatch(mmi_as03["main"], idx_base)])
#y       = pyplottable(as03["J_main"][:,idxmatch(mmi_as03["main"], idx_base)] .- as03["ΔJ_main"][:,idxmatch(mmi_as03["main"], idx_base)] .+ Jsm)
plot_hexbins( log10.(x), log10.(y), 
    outfile=outfile,
    title="$(latexstring("\$\\Delta\$"))t $(latexstring("\$\\sim\$")) 0.3 Gyr, Immediate Read",
    scale=15, gridres=(20,8), xscale="linear", yscale="linear", 
    ylabel="Flip by Merger Contribution [°]",
    xlabel="Flip Angle [°]",
    grid=true,
    xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax,
    lognorm=true,
    binned_median = "x", plot_bval=" ", plot_line=hcat([xmin,xmax],[ymin,ymax]),
    calc_norm=true,    # @pyplottable
    )
#


outfile = outnm*"hxbvalue.png"
plot_hexbins( log10.(x), log10.(y), 
    outfile=outfile,
    title="$(latexstring("\$\\Delta\$"))t $(latexstring("\$\\sim\$")) 0.3 Gyr, Immediate Read",
    scale=15, gridres=(20,8), xscale="linear", yscale="linear", 
    ylabel="Flip by Merger Contribution [°]",
    xlabel="Flip Angle [°]",
    grid=true,
    xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax,
    lognorm=false,
    binned_median = "x", plot_bval=" ", plot_line=hcat([xmin,xmax],[ymin,ymax]),
    calc_norm=true,    # @pyplottable
    cmap=ColorMap(get_bvalue_cmap().colors), colormod = colr, cfunc=mean, colorlimits=(-6.1,-3.9),
)
#


xmin = 0#" "#1e0
ymin = 0#" "#xmin
xmax = 10#" "#180
ymax = 10#" "#xmax
#y       = pyplottable(delJ[:,order] .- Jom[:,order])
colr    = as03["BVAL_0"][idxmatch(mmi_as03["main"], idx_base)] .- as03["ΔBVAL_0"][idxmatch(mmi_as03["main"], idx_base)] # log10.(angle_J1_Jsm ./ x)
nomissings = idxclude(1:length(colr), findcs(colr, eq=missing))
colr    = colr[nomissings]
order   = sortperm(colr)
colr    = colr[order]
x       = pyplottable(as03["ϕ_flip"][idxmatch(mmi_as03["main"], idx_base)])[nomissings][order]
y       = pyplottable(angle_J0_Jsm)[nomissings][order]
outnm     = joinpath(outputdir,"angmerger_v_ang-imm-03")
outfile = outnm*"hxcount.png"
#x       = pyplottable(as03["J_main"][:,idxmatch(mmi_as03["main"], idx_base)])
#y       = pyplottable(as03["J_main"][:,idxmatch(mmi_as03["main"], idx_base)] .- as03["ΔJ_main"][:,idxmatch(mmi_as03["main"], idx_base)] .+ Jsm)
plot_hexbins( x, y, 
    outfile=outfile,
    title="$(latexstring("\$\\Delta\$"))t $(latexstring("\$\\sim\$")) 0.3 Gyr, Immediate Read",
    scale=15, gridres=(20,8), xscale="linear", yscale="linear", 
    ylabel="Flip by Merger Contribution [°]",
    xlabel="Flip Angle [°]",
    grid=true,
    xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax,
    lognorm=true,
    binned_median = "x", plot_bval=" ", plot_line=hcat([xmin,xmax],[ymin,ymax]),
    calc_norm=true,    # @pyplottable
    )
#


outfile = outnm*"hxbvalue.png"
plot_hexbins( x, y, 
    outfile=outfile,
    title="$(latexstring("\$\\Delta\$"))t $(latexstring("\$\\sim\$")) 0.3 Gyr, Immediate Read",
    scale=15, gridres=(20,8), xscale="linear", yscale="linear", 
    ylabel="Flip by Merger Contribution [°]",
    xlabel="Flip Angle [°]",
    grid=true,
    xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax,
    lognorm=false,
    binned_median = "x", plot_bval=" ", plot_line=hcat([xmin,xmax],[ymin,ymax]),
    calc_norm=true,    # @pyplottable
    cmap=ColorMap(get_bvalue_cmap().colors), colormod = colr, colorlimits=(-6.1,-3.9), cfunc=mean,
)




###########################
# dj / j vs flipangle

outnm   = joinpath(outputdir,"djrel_v_flipangle_vs_bval_03Gyr")
outfile = outnm*".png"
xmin    = 0
xmax    = 180
ymin    = -2
ymax    = 1
scale = 0.5
PyPlot.matplotlib[:rc]("text", usetex=true)
PyPlot.matplotlib[:rc]("font", family="serif", serif=["Computer Modern"])
pointsize = 2
fig, ax = subplots()
slct    = idxclude(idxclude(csw03[:start_thr], csw03[:switches]), csw03[:fakes])
y       = pyplottable(as03["Δj_main"][:,slct]) ./ pyplottable(as03["j_main"][:,slct] .- as03["Δj_main"][:,slct])
x       = pyplottable(as03["ϕ_flip"][slct])
colr    = as03["BVAL_0"][slct] .- as03["ΔBVAL_0"][slct]
nomissings = idxclude(1:length(colr), findcs(colr, eq=missing))
colr    = colr[nomissings]
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
ax.set_ylim([ymin,ymax])
ax.set_xlim([xmin,xmax])
ax.legend(loc="lower right", frameon=true, borderpad=1, handlelength=1.8)
ax.set_ylabel("log( $(latexstring("\$_\\Delta\$"))j$(latexstring("\$_*\$")) / j$(latexstring("\$_\\textrm{*,pre}\$")) )")
ax.set_xlabel("Instant Flip Angle [°]")
ax.grid(true)
fig.set_size_inches(11*scale, 9*scale)
fig.tight_layout()
fig.savefig(outfile, bbox_inches="tight", pad_inches=.1)
close("all")
PyPlot.matplotlib[:rc]("text", usetex=false)


plot_hexbins( x[order], log10.(norm.(y[order])), 
    outfile=outnm*"_hxcount.png",
    scale=15, gridres=(20,8), xscale="linear", yscale="linear", 
    xlabel="Instant Flip Angle [°]", 
    ylabel="log( $(latexstring("\$_\\Delta\$"))j$(latexstring("\$_*\$")) / j$(latexstring("\$_\\textrm{*,pre}\$")) )", 
    grid=true,
    xmin=0, xmax=180, ymin=ymin, ymax=ymax,
    lognorm=true,
    binned_median = "x", plot_bval=" ", plot_line=hcat(LinRange(0, 180, 1000), log10.(min_dj_v_flipangled.(LinRange(0, 180, 1000)))),
    calc_norm=true,    # @pyplottable
    )
#


plot_hexbins( x[order], log10.(norm.(y[order])), 
outfile=outnm*"_hxbvalpost.png",
scale=15, gridres=(20,8), xscale="linear", yscale="linear", 
xlabel="Instant Flip Angle [°]", 
ylabel="log( $(latexstring("\$_\\Delta\$"))j$(latexstring("\$_*\$")) / j$(latexstring("\$_\\textrm{*,post}\$")) )", 
grid=true,
xmin=0, xmax=180, ymin=ymin, ymax=ymax,
lognorm=false,
binned_median = "x", plot_bval=" ", plot_line=hcat(LinRange(0, 180, 1000), log10.(min_dj_v_flipangled.(LinRange(0, 180, 1000)))),
calc_norm=true,    # @pyplottable
cmap = ColorMap(get_bvalue_cmap().colors), colormod = colr, colorlimits=(-6.1,-3.9), cfunc=mean,
)
#

# bval at end
colr    = as03["BVAL_0"][slct]
nomissings = idxclude(1:length(colr), findcs(colr, eq=missing))
colr    = colr[nomissings]
order   = sortperm(colr)#[length(colr):-1:1]
colr    = colr[order]
maxcolr = maximum(colr)
plot_hexbins( x[order], log10.(norm.(y[order])), 
outfile=outnm*"_hxbval.png",
scale=15, gridres=(20,8), xscale="linear", yscale="linear", 
xlabel="Instant Flip Angle [°]", 
ylabel="log( $(latexstring("\$_\\Delta\$"))j$(latexstring("\$_*\$")) / j$(latexstring("\$_\\textrm{*,pre}\$")) )", 
grid=true,
xmin=0, xmax=180, ymin=ymin, ymax=ymax,
lognorm=false,
binned_median = "x", plot_bval=" ", plot_line=hcat(LinRange(0, 180, 1000), log10.(min_dj_v_flipangled.(LinRange(0, 180, 1000)))),
calc_norm=true,    # @pyplottable
cmap = ColorMap(get_bvalue_cmap().colors), colormod = colr, colorlimits=(-6.1,-3.9), cfunc=mean,
)
#


outnm   = joinpath(outputdir,"flipangle_v_djrel_vs_bval_03Gyr")
outfile = outnm*".png"
ymin    = 0
ymax    = 180
xmin    = -2
xmax    = 1
scale = 0.5
PyPlot.matplotlib[:rc]("text", usetex=true)
PyPlot.matplotlib[:rc]("font", family="serif", serif=["Computer Modern"])
pointsize = 2
fig, ax = subplots()
slct    = idxclude(idxclude(csw03[:start_thr], csw03[:switches]), csw03[:fakes])
x       = pyplottable(as03["Δj_main"][:,slct]) ./ pyplottable(as03["j_main"][:,slct] .- as03["Δj_main"][:,slct])
y       = pyplottable(as03["ϕ_flip"][slct])
colr    = as03["BVAL_0"][slct] .- as03["ΔBVAL_0"][slct]
nomissings = idxclude(1:length(colr), findcs(colr, eq=missing))
colr    = colr[nomissings]
order   = sortperm(colr)#[length(colr):-1:1]
colr    = colr[order]
maxcolr = maximum(colr)
ax.grid()
p = ax.scatter(  log10.(norm.(x[order])),y[order],
        zorder=2   , s= pointsize    , alpha=0.5,
        cmap=ColorMap(get_bvalue_cmap().colors), c=colr, vmin=-6.1,vmax=-3.9, label="Sample, $(latexstring("\$\\Delta\$"))t = 0.3 Gyr") 
ax.plot(log10.(LinRange(0.01, 10, 1000)), max_flipangled_v_dj.(LinRange(0.01, 10, 1000)), "-",color="black", linewidth=2, label="Limit")
colbar  = fig.colorbar(p, ax=ax)
colbar.ax.tick_params(axis="y")
colbar.ax.set_ylabel("b-value at Flip Start")
colbar.ax.yaxis.set_label_position("left")
ax.set_ylim([ymin,ymax])
ax.set_xlim([xmin,xmax])
ax.legend(loc="upper left", frameon=true, borderpad=1, handlelength=1.8)
ax.set_ylabel("log( $(latexstring("\$_\\Delta\$"))j$(latexstring("\$_*\$")) / j$(latexstring("\$_\\textrm{*,pre}\$")) )")
ax.set_xlabel("Instant Flip Angle [°]")
ax.grid(true)
fig.set_size_inches(11*scale, 9*scale)
fig.tight_layout()
fig.savefig(outfile, bbox_inches="tight", pad_inches=.1)
close("all")
PyPlot.matplotlib[:rc]("text", usetex=false)


plot_hexbins( log10.(norm.(x[order])), y[order],
    outfile=outnm*"_hxcount.png",
    scale=15, gridres=(20,8), xscale="linear", yscale="linear", 
    ylabel="Instant Flip Angle [°]", 
    xlabel="log( $(latexstring("\$_\\Delta\$"))j$(latexstring("\$_*\$")) / j$(latexstring("\$_\\textrm{*,pre}\$")) )", 
    grid=true,
    xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax,
    lognorm=true,
    binned_median = "x", plot_bval=" ", plot_line=hcat(log10.(LinRange(0.01, 10, 1000)), max_flipangled_v_dj.(LinRange(0.01, 10, 1000))),
    calc_norm=true,    # @pyplottable
    )
#



plot_hexbins( log10.(norm.(x[order])), y[order],
    outfile=outnm*"_hxbval.png",
    scale=15, gridres=(20,8), xscale="linear", yscale="linear", 
    ylabel="Instant Flip Angle [°]", 
    xlabel="log( $(latexstring("\$_\\Delta\$"))j$(latexstring("\$_*\$")) / j$(latexstring("\$_\\textrm{*,pre}\$")) )", 
    grid=true,
    xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax,
    lognorm=false,
    binned_median = "x", plot_bval=" ", plot_line=hcat(log10.(LinRange(0.01, 10, 1000)), max_flipangled_v_dj.(LinRange(0.01, 10, 1000))),
    calc_norm=true,    # @pyplottable
    cmap = ColorMap(get_bvalue_cmap().colors), colormod = colr, colorlimits=(-6.1,-3.9), cfunc=mean,
    )
#

# bval at end
colr    = as03["BVAL_0"][slct]
nomissings = idxclude(1:length(colr), findcs(colr, eq=missing))
colr    = colr[nomissings]
order   = sortperm(colr)#[length(colr):-1:1]
colr    = colr[order]
maxcolr = maximum(colr)
plot_hexbins( log10.(norm.(x[order])), y[order],
    outfile=outnm*"_hxbvalpost.png",
    scale=15, gridres=(20,8), xscale="linear", yscale="linear", 
    ylabel="Instant Flip Angle [°]", 
    xlabel="log( $(latexstring("\$_\\Delta\$"))j$(latexstring("\$_*\$")) / j$(latexstring("\$_\\textrm{*,post}\$")) )", 
    grid=true,
    xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax,
    lognorm=false,
    binned_median = "x", plot_bval=" ", plot_line=hcat(log10.(LinRange(0.01, 10, 1000)), max_flipangled_v_dj.(LinRange(0.01, 10, 1000))),
    calc_norm=true,    # @pyplottable
    cmap = ColorMap(get_bvalue_cmap().colors), colormod = colr, cfunc=mean, colorlimits=(-6.1,-3.9), 
    )
#






























######################
# Aftermath

# flip limit multiplots
outnm   = joinpath(outputdir,"flipangle_v_djrelpre_03Gyr")
outfile = outnm*"_hxcount.png"
y       = pyplottable(as03["ϕ_flip"])
x       = pyplottable(log10.(pyplottable(as03["Δj_main"]) ./ pyplottable(as03["j_main"] .- as03["Δj_main"])))
colr    = as03["BVAL_0"] .- as03["ΔBVAL_0"]
slct    = idxclude(idxclude(csw03[:start_thr], csw03[:switches]), csw03[:fakes])
miny= 0#minimum(y[slct])
maxy= 180#maximum(y[slct])
minx= -2.5#minimum(x[slct])
maxx= 1#maximum(x[slct])
slcts       = Dict{Int64, Vector{Int64}}()
slcts[1]    = slct
slcts[3]    = idxmatch(slct, csw03[:start_disk])
slcts[2]    = idxmatch(slct, csw03[:start_int])
slcts[4]    = idxmatch(slct, csw03[:start_ell])
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
            xlabel="log( $(latexstring("\$\\Delta\$"))j$(latexstring("\$_*\$")) / j$(latexstring("\$_\\textrm{*,pre}\$")) )",
            ylabel="Instant Flip Angle [°]",
            fontsize=30,
            binned_median="y",
            plot_line=hcat(log10.(LinRange(0.01, 10, 1000)), max_flipangled_v_dj.(LinRange(0.01, 10, 1000)))
            )
#

outfile = outnm*"_hxbval.png"
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
            xlabel="log( $(latexstring("\$\\Delta\$"))j$(latexstring("\$_*\$")) / j$(latexstring("\$_\\textrm{*,pre}\$")) )",
            ylabel="Instant Flip Angle [°]",
            fontsize=30,
            binned_median="y",
            cfunc=mean, colormod=colr, cmap=ColorMap(get_bvalue_cmap().colors), colorlimits=(-6.1,-3.9),
            plot_line=hcat(log10.(LinRange(0.01, 10, 1000)), max_flipangled_v_dj.(LinRange(0.01, 10, 1000)))
            )
#



outnm   = joinpath(outputdir,"flipangle_v_djrelpost_03Gyr")
outfile = outnm*"_hxcount.png"
y       = pyplottable(as03["ϕ_flip"])
x       = pyplottable(log10.(pyplottable(as03["Δj_main"]) ./ pyplottable(as03["j_main"])))
colr    = as03["BVAL_0"]
slct    = idxclude(idxclude(csw03[:start_thr], csw03[:switches]), csw03[:fakes])
miny= 0#minimum(y[slct])
maxy= 180#maximum(y[slct])
minx= -2.5#minimum(x[slct])
maxx= 1#maximum(x[slct])
slcts       = Dict{Int64, Vector{Int64}}()
slcts[1]    = slct
slcts[3]    = idxmatch(slct, csw03[:start_disk])
slcts[2]    = idxmatch(slct, csw03[:start_int])
slcts[4]    = idxmatch(slct, csw03[:start_ell])
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
            xlabel="log( $(latexstring("\$\\Delta\$"))j$(latexstring("\$_*\$")) / j$(latexstring("\$_\\textrm{*,post}\$")) )",
            ylabel="Instant Flip Angle [°]",
            fontsize=30,
            binned_median="y",
            plot_line=hcat(log10.(LinRange(0.01, 10, 1000)), max_flipangled_v_dj.(LinRange(0.01, 10, 1000)))
            )
#

outfile = outnm*"_hxbval.png"
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
            xlabel="log( $(latexstring("\$\\Delta\$"))j$(latexstring("\$_*\$")) / j$(latexstring("\$_\\textrm{*,post}\$")) )",
            ylabel="Instant Flip Angle [°]",
            fontsize=30,
            binned_median="y",
            cfunc=mean, colormod=colr, cmap=ColorMap(get_bvalue_cmap().colors), colorlimits=(-6.1,-3.9),
            plot_line=hcat(log10.(LinRange(0.01, 10, 1000)), max_flipangled_v_dj.(LinRange(0.01, 10, 1000)))
            )
#










###################################################################
# predict delta J
xmin = 11.5#" "#1e0
ymin = 10#" "#xmin
xmax = 15#" "#180
ymax = 16#" "#xmax
outnm     = joinpath(outputdir,"J0Jsm_v_dJ10-imm-03")
outfile = outnm*"_hxcount.png"
x       = log10.(pyplottable(as03["ΔJ_main"][:,idxmatch(mmi_as03["main"], idx_base)]))
y       = log10.(pyplottable(Jsm))
#x       = pyplottable(as03["J_main"][:,idxmatch(mmi_as03["main"], idx_base)])
#y       = pyplottable(as03["J_main"][:,idxmatch(mmi_as03["main"], idx_base)] .- as03["ΔJ_main"][:,idxmatch(mmi_as03["main"], idx_base)] .+ Jsm)
plot_hexbins( x, y, 
    outfile=outfile,
    title="$(latexstring("\$\\Delta\$"))t $(latexstring("\$\\sim\$")) 0.3 Gyr, Immediate Read",
    scale=15, gridres=(20,8), xscale="linear", yscale="linear", 
    xlabel="log( $(latexstring("\$\\Delta\$"))J$(latexstring("\$_\\textrm{10}\$")) [ M$(latexstring("\$_\\odot\$")) kpc km/s ] )", 
    ylabel="log( J$(latexstring("\$_\\textrm{Mergers}\$")) [ M$(latexstring("\$_\\odot\$")) kpc km/s ] )", 
    grid=true,
    xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax,
    lognorm=true,
    binned_median = "y", plot_bval=" ", plot_line=hcat([9,18],[9,18]),
    calc_norm=true,    # @pyplottable
    )
#

###################################################################
# predict delta J
xmin = 11.5#" "#1e0
ymin = 10#" "#xmin
xmax = 15#" "#180
ymax = 16#" "#xmax
outnm     = joinpath(outputdir,"J0Jsm_v_dJ100-imm-03")
outfile = outnm*"_hxcount.png"
x       = log10.(pyplottable(as03["ΔJ_vir"][:,idxmatch(mmi_as03["main"], idx_base)]))
y       = log10.(pyplottable(Jsm))
#x       = pyplottable(as03["J_main"][:,idxmatch(mmi_as03["main"], idx_base)])
#y       = pyplottable(as03["J_main"][:,idxmatch(mmi_as03["main"], idx_base)] .- as03["ΔJ_main"][:,idxmatch(mmi_as03["main"], idx_base)] .+ Jsm)
plot_hexbins( x, y, 
    outfile=outfile,
    title="$(latexstring("\$\\Delta\$"))t $(latexstring("\$\\sim\$")) 0.3 Gyr, Immediate Read",
    scale=15, gridres=(20,8), xscale="linear", yscale="linear", 
    xlabel="log( $(latexstring("\$\\Delta\$"))J$(latexstring("\$_\\textrm{100}\$")) [ M$(latexstring("\$_\\odot\$")) kpc km/s ] )", 
    ylabel="log( J$(latexstring("\$_\\textrm{Mergers}\$")) [ M$(latexstring("\$_\\odot\$")) kpc km/s ] )", 
    grid=true,
    xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax,
    lognorm=true,
    binned_median = "y", plot_bval=" ", plot_line=hcat([9,18],[9,18]),
    calc_norm=true,    # @pyplottable
    )
#

# bval mean
xmin = 11.5#" "#1e0
ymin = 10#" "#xmin
xmax = 15#" "#180
ymax = 16#" "#xmax
outnm     = joinpath(outputdir,"J0Jsm_v_dJ10-imm-03")
outfile = outnm*"_hxbval.png"
x       = log10.(pyplottable(as03["ΔJ_main"][:,idxmatch(mmi_as03["main"], idx_base)]))
y       = log10.(pyplottable(Jsm))
#x       = pyplottable(as03["J_main"][:,idxmatch(mmi_as03["main"], idx_base)])
#y       = pyplottable(as03["J_main"][:,idxmatch(mmi_as03["main"], idx_base)] .- as03["ΔJ_main"][:,idxmatch(mmi_as03["main"], idx_base)] .+ Jsm)
plot_hexbins( x, y, 
    outfile=outfile,
    title="$(latexstring("\$\\Delta\$"))t $(latexstring("\$\\sim\$")) 0.3 Gyr, Immediate Read",
    scale=15, gridres=(20,8), xscale="linear", yscale="linear", 
    xlabel="log( $(latexstring("\$\\Delta\$"))J$(latexstring("\$_\\textrm{10}\$")) [ M$(latexstring("\$_\\odot\$")) kpc km/s ] )", 
    ylabel="log( J$(latexstring("\$_\\textrm{Mergers}\$")) [ M$(latexstring("\$_\\odot\$")) kpc km/s ] )", 
    grid=true,
    xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax,
    lognorm=false,
    binned_median = "y", plot_bval=" ", plot_line=hcat([9,18],[9,18]),
    cmap = ColorMap(get_bvalue_cmap().colors), colormod = as03["BVAL_0"][idxmatch(mmi_as03["main"], idx_base)] .-as03["ΔBVAL_0"][idxmatch(mmi_as03["main"], idx_base)], colorlimits=(-6.1,-3.9), cfunc=mean,
    calc_norm=true,    # @pyplottable
    )
#

###################################################################
# predict delta J
xmin = 11.5#" "#1e0
ymin = 10#" "#xmin
xmax = 15#" "#180
ymax = 16#" "#xmax
outnm     = joinpath(outputdir,"J0Jsm_v_dJ100-imm-03")
outfile = outnm*"_hxbval.png"
x       = log10.(pyplottable(as03["ΔJ_vir"][:,idxmatch(mmi_as03["main"], idx_base)]))
y       = log10.(pyplottable(Jsm))
#x       = pyplottable(as03["J_main"][:,idxmatch(mmi_as03["main"], idx_base)])
#y       = pyplottable(as03["J_main"][:,idxmatch(mmi_as03["main"], idx_base)] .- as03["ΔJ_main"][:,idxmatch(mmi_as03["main"], idx_base)] .+ Jsm)
plot_hexbins( x, y, 
    outfile=outfile,
    title="$(latexstring("\$\\Delta\$"))t $(latexstring("\$\\sim\$")) 0.3 Gyr, Immediate Read",
    scale=15, gridres=(20,8), xscale="linear", yscale="linear", 
    xlabel="log( $(latexstring("\$\\Delta\$"))J$(latexstring("\$_\\textrm{100}\$")) [ M$(latexstring("\$_\\odot\$")) kpc km/s ] )", 
    ylabel="log( J$(latexstring("\$_\\textrm{Mergers}\$")) [ M$(latexstring("\$_\\odot\$")) kpc km/s ] )", 
    grid=true,
    xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax,
    lognorm=false,
    binned_median = "y", plot_bval=" ", plot_line=hcat([9,18],[9,18]),
    calc_norm=true,    # @pyplottable
    cmap = ColorMap(get_bvalue_cmap().colors), colormod = as03["BVAL_0"][idxmatch(mmi_as03["main"], idx_base)] .-as03["ΔBVAL_0"][idxmatch(mmi_as03["main"], idx_base)], colorlimits=(-6.1,-3.9), cfunc=mean,
    )
#

outnm   = joinpath(outputdir,"mergerflip_v_flipangle_imm-03Gyr")
outfile = outnm*"_hxcount.png"
slct    = idxmatch(mmi_as03["main"], idx_base)
x       = pyplottable(as03["ϕ_flip"][slct])
y       = pyplottable(angle_J0_Jsm)
colr    = as03["BVAL_0"] .- as03["ΔBVAL_0"]
disks   = findall(dummy -> dummy .≥ b_disk, as03["BVAL_0"][slct])
ints   = findall(dummy -> b_ell .≤ dummy .≤ b_disk, as03["BVAL_0"][slct])
ells   = findall(dummy -> dummy .< b_ell, as03["BVAL_0"][slct])
miny= 0#minimum(y[slct])
maxy= 10#maximum(y[slct])
minx= 0#minimum(x[slct])
maxx= 10#maximum(x[slct])
slcts       = Dict{Int64, Vector{Int64}}()
slcts[1]    = 1:length(x)
slcts[3]    = disks
slcts[2]    = ints
slcts[4]    = ells
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
            ylabel="Predicted Flip Angle [°]",
            xlabel="Instant Flip Angle [°]",
            fontsize=30,
            binned_median="y",
            plot_line=hcat(LinRange(0, 180, 10), LinRange(0, 180, 10) )
            )
#
