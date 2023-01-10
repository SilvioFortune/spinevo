
old_as  = load("/home/moon/sfortune/spinevo/data/assembly_centrals2co_20220112.jld", "assembly_STARS")

old_ad  = load("/home/moon/sfortune/spinevo/data/assembly_centrals2co_20220112.jld", "assembly_DM")

include("/home/moon/sfortune/spinevo/pkg/meta.jl")

as      = load("/home/moon/sfortune/spinevo/data/updated_assembly_M2e10.jld", "assembly_STARS")

ad      = load("/home/moon/sfortune/spinevo/data/updated_assembly_M2e10.jld", "assembly_DM")



# Selections

# borders: snap52=z1.18, snap100=z0.42
z_hi        = findcs(as["snapNR"      ], leq=52)
z_md        = findcs(as["snapNR"      ], gt=52, lt=100)
z_lo        = findcs(as["snapNR"      ], geq=100)
snap136     = findcs(as["snapNR"      ], eq=136)
snap36      = findcs(as["snapNR"      ], eq=36)

halo3212    = findcs(as["ID_ISUB"      ], eq=3212)
halo1414    = findcs(as["ID_ISUB"      ], eq=1414)

switches    = findcs(as["switch"      ], eq=1)
n_switches  = findcs(as["switch"      ], eq=0)
fakes       = findcs(as["FAKEFLIP"    ], eq=1)
n_fakes     = findcs(as["FAKEFLIP"    ], eq=0)
old_n_fakes = findcs(old_as["FAKEFLIP"], eq=0)

result_disks    = findcs(as["BVAL"  ], geq=b_disk)
result_int      = findcs(as["BVAL"  ], gt=b_ell, lt=b_disk)
result_ell      = findcs(as["BVAL"  ], leq=b_ell)

start_disks     = findcs(as["BVAL_0"  ] .- as["ΔBVAL_0"    ], geq=b_disk)
start_int       = findcs(as["BVAL_0"  ] .- as["ΔBVAL_0"    ], gt=b_ell, lt=b_disk)
start_ell       = findcs(as["BVAL_0"  ] .- as["ΔBVAL_0"    ], leq=b_ell)
stay_disk       = start_disks[  start_disks .∈  Ref(Set(result_disks    ))]
stay_int        = start_int[    start_int   .∈  Ref(Set(result_int      ))]
stay_ell        = start_ell[    start_ell   .∈  Ref(Set(result_ell      ))]

result_disk_fakes    = fakes[   fakes   .∈  Ref(Set(result_disks    ))]
result_int_fakes     = fakes[   fakes   .∈  Ref(Set(result_int      ))]
result_ell_fakes     = fakes[   fakes   .∈  Ref(Set(result_ell      ))]
result_disk_nfakes   = n_fakes[ n_fakes .∈  Ref(Set(result_disks    ))]
result_int_nfakes    = n_fakes[ n_fakes .∈  Ref(Set(result_int      ))]
result_ell_nfakes    = n_fakes[ n_fakes .∈  Ref(Set(result_ell      ))]

start_disk_fakes    = fakes[   fakes   .∈  Ref(Set(start_disks    ))]
start_int_fakes     = fakes[   fakes   .∈  Ref(Set(start_int      ))]
start_ell_fakes     = fakes[   fakes   .∈  Ref(Set(start_ell      ))]
start_disk_nfakes   = n_fakes[ n_fakes .∈  Ref(Set(start_disks    ))]
start_int_nfakes    = n_fakes[ n_fakes .∈  Ref(Set(start_int      ))]
start_ell_nfakes    = n_fakes[ n_fakes .∈  Ref(Set(start_ell      ))]

stay_disk_fakes     = fakes[   fakes   .∈  Ref(Set(stay_disk     ))]
stay_int_fakes      = fakes[   fakes   .∈  Ref(Set(stay_int      ))]
stay_ell_fakes      = fakes[   fakes   .∈  Ref(Set(stay_ell      ))]
stay_disk_nfakes    = n_fakes[ n_fakes .∈  Ref(Set(stay_disk     ))]
stay_int_nfakes     = n_fakes[ n_fakes .∈  Ref(Set(stay_int      ))]
stay_ell_nfakes     = n_fakes[ n_fakes .∈  Ref(Set(stay_ell      ))]

ell_to_disk_nfakes  = start_ell_nfakes[start_ell_nfakes     .∈ Ref(Set(result_disk_nfakes))]
ell_to_int_nfakes   = start_ell_nfakes[start_ell_nfakes     .∈ Ref(Set(result_int_nfakes))]
int_to_disk_nfakes  = start_int_nfakes[start_int_nfakes     .∈ Ref(Set(result_disk_nfakes))]
int_to_ell_nfakes   = start_int_nfakes[start_int_nfakes     .∈ Ref(Set(result_ell_nfakes))]
disk_to_ell_nfakes  = start_disk_nfakes[start_disk_nfakes   .∈ Ref(Set(result_ell_nfakes))]
disk_to_int_nfakes  = start_disk_nfakes[start_disk_nfakes   .∈ Ref(Set(result_int_nfakes))]

mergersDM_1to50 = findcs(ad["M2"]./ad["Mpeak_MM"], geq=1, leq=50)

# b-value vs j_STARS - Data Sanity
y   = as["BVAL_0"]
x   = log10.(pyplottable(as["j_main"]))
plot_hexbins(   x, y,
            #selection = findcs(pyplottable(testas2["snapNR"]), eq=snap),
            outfile="/home/moon/sfortune/spinevo/plots/sanity/bval0_vs_jSTARS.png",
            lognorm=true, gridres=(20,10), grid=true,
            #xmin=0, xmax=3.5, ymin=0, ymax=180,
            ylabel="b-value",
            xlabel="log( j_STARS )", 
            #title="Uncorrected & Using M_tree",
            binned_median="x", plot_bval="y",
            colorlimits=(0.8,nothing)
            )

# b-value + M vs j_STARS - Data Sanity
y   = as["BVAL_0"] .+ ((2/3) .* log10.(as["M"]))
x   = log10.(pyplottable(as["j_main"]))
plot_hexbins(   x, y,
            #selection = findcs(pyplottable(testas2["snapNR"]), eq=snap),
            outfile="/home/moon/sfortune/spinevo/plots/sanity/bval0_M_vs_jSTARS.png",
            lognorm=true, gridres=(20,10), grid=true,
            #xmin=0, xmax=3.5, ymin=0, ymax=180,
            ylabel="b-value + 2/3*log(M)",
            xlabel="log( j_STARS )", 
            #title="Uncorrected & Using M_tree",
            binned_median="x", #plot_bval="y",
            )

# b-value + M2 vs j_STARS - Data Sanity
y   = as["BVAL_0"] .+ ((2/3) .* log10.(as["M2"]))
x   = log10.(pyplottable(as["j_main"]))
plot_hexbins(   x, y,
            #selection = findcs(pyplottable(testas2["snapNR"]), eq=snap),
            outfile="/home/moon/sfortune/spinevo/plots/sanity/bval0_M2_vs_jSTARS.png",
            lognorm=true, gridres=(20,10), grid=true,
            #xmin=0, xmax=3.5, ymin=0, ymax=180,
            ylabel="b-value + 2/3*log(M2)",
            xlabel="log( j_STARS )", 
            #title="Uncorrected & Using M_tree",
            binned_median="x", #plot_bval="y",
            )


# b-value + Mpeak vs j_STARS - Data Sanity
y   = as["BVAL_0"] .+ ((2/3) .* log10.(as["Mpeak"]))
x   = log10.(pyplottable(as["j_main"]))
plot_hexbins(   x, y,
            #selection = findcs(pyplottable(testas2["snapNR"]), eq=snap),
            outfile="/home/moon/sfortune/spinevo/plots/sanity/bval0_Mpeak_vs_jSTARS.png",
            lognorm=true, gridres=(20,10), grid=true,
            #xmin=0, xmax=3.5, ymin=0, ymax=180,
            ylabel="b-value + 2/3*log(Mpeak)",
            xlabel="log( j_STARS )", 
            #title="Uncorrected & Using M_tree",
            binned_median="x", #plot_bval="y",
            )


# M v M2 v Mpeak
fig, ax = subplots()
#ax.hist( log10.(pyplottable(as["M"          ]))                             , bins=50, label="M"     , rwidth=0.9    , alpha=0.5, color="navy" )
#ax.hist( log10.(pyplottable(as["Mpeak"      ]))                             , bins=50, label="Mpeak" , rwidth=0.7    , alpha=0.5, color="green" )
#ax.hist( log10.(pyplottable(as["M2" ][findcs(as["M2"], geq=1e10)]))   , bins=50, label="M2"    , rwidth=0.5    , alpha=0.5, color="red" )
ax.hist( pyplottable(log10.(abs.(as["M"          ] .- as["M2" ]) ./ as["M2"      ] ))    , bins=50, label="(M - M2) / M2"       , rwidth=0.7    , alpha=0.5, color="green" )
ax.hist( pyplottable(log10.(abs.(as["M"          ] .- as["Mpeak" ]) ./ as["Mpeak"   ] ))    , bins=50, label="(M - Mpeak) / Mpeak"    , rwidth=0.9    , alpha=0.5, color="navy" )
ax.hist( pyplottable(log10.(abs.(as["Mpeak"      ] .- as["M2" ]) ./ as["M2"  ] ))    , bins=50, label="(Mpeak - M2) / M2"   , rwidth=0.5    , alpha=0.5, color="red" )
ax.set_xlabel("log( ΔM / M )")
#ax.set_xscale("log")
ax.set_yscale("log")
#ax.set_xlim([10, nothing])
ax.legend(loc="upper right", frameon=true, borderpad=1, handlelength=1.8)
ax.grid()
scale=0.5
fig.set_size_inches(16scale, 9scale)
fig.tight_layout()
fig.savefig("/home/moon/sfortune/spinevo/plots/sanity/Mmain_comparison2.png", bbox_inches="tight", pad_inches=.1)




# Mergers vs Flip
y   = as["ϕ_flip"]
x   = ad["M2"]./ad["Mpeak_MM"]
slct= n_fakes[ n_fakes .∈  Ref(Set(mergersDM_1to50     ))]
plot_hexbins(   x, y,
            selection = slct,
            outfile="/home/moon/sfortune/spinevo/plots/sanity/mergerRatio2DM_vs_flipangle_zoom.png",
            lognorm=true, gridres=(20,10), grid=true,
            ymin=0, ymax=180, 
            #xmin=1, xmax=50, 
            xlabel="Merger DM Ratio",
            ylabel="ϕ_flip [°]", 
            #title="Uncorrected & Using M_tree",
            binned_median="x",
            )





# j_STARS vs M_STARS
slctn = 1:length(as["BVAL_0"])
#slctn = snap136[ snap136 .∈  Ref(Set(result_disk_nfakes      ))]
fig, ax = subplots()
ax.scatter(pyplottable(log10.(replace(as["M"    ][result_disks], 0.0 => NaN))), 
    pyplottable(log10.(pyplottable(as["j_main"   ][:,result_disks])))    ,  label="Disks"           , 
    zorder=3   , s=0.6   , alpha=0.6, color="mediumblue" )
ax.scatter(pyplottable(log10.(replace(as["M"    ][result_int]  , 0.0 => NaN))), 
    pyplottable(log10.(pyplottable(as["j_main"   ][:,result_int]  )))    ,  label="Intermediates"   , 
    zorder=1   , s=0.6   , alpha=0.6, color="darkorchid" )
ax.scatter(pyplottable(log10.(replace(as["M"    ][result_ell]  , 0.0 => NaN))), 
    pyplottable(log10.(pyplottable(as["j_main"   ][:,result_ell]  )))    ,  label="Ellipticals"     , 
    zorder=2   , s=0.6   , alpha=0.6, color="darkred" )
#p = ax.scatter(pyplottable(log10.(replace(as["M"    ][slctn]   , 0.0 => NaN))), 
#    pyplottable(log10.(pyplottable(as["j_main"   ][:,slctn]  ))) , 
#    zorder=2   , s=2   , alpha=0.9, 
#    cmap="brg_r", vmin=minimum(as["BVAL"    ]), vmax=maximum(as["BVAL"    ]), c=as["BVAL"    ][slctn]) # brg_r   
x_min = ax.get_xlim()[1]
x_max = ax.get_xlim()[2]
ax.plot([x_min, x_max], [b_disk+2/3*x_min, b_disk+2/3*x_max],"-" , color="mediumblue")
ax.plot([x_min, x_max], [b_ell+2/3*x_min, b_ell+2/3*x_max],"-" , color="darkred")
#colbar  = fig.colorbar(p, ax=ax, ticks=[minimum(as["BVAL_0"    ]), b_ell, b_disk, maximum(as["BVAL_0"    ])])
#colbar.ax.set_yticklabels(["$(round(minimum(as["BVAL_0"    ]), digits=1))", "ETG", "LTG", "$(round(maximum(as["BVAL_0"    ]), digits=1))"])
scale=0.5
#colbar.ax.tick_params(axis="y")
###colbar.ax.set_yticks("ticks")
ax.legend(loc="lower right", frameon=true, borderpad=1, handlelength=1.8)

#ax.set_title("Halo 3212")
ax.set_xlabel("log( M_⊙ )")
ax.set_ylabel("log( j_STARS )")
ax.grid()
ax.autoscale()
fig.set_size_inches(16scale, 9scale)
fig.tight_layout()
fig.savefig("/home/moon/sfortune/spinevo/plots/sanity/jSTARS_vs_M_BVALdiscrete.png", bbox_inches="tight", pad_inches=.1)








# Testing
slcts       = Dict{Int64, Vector{Int64}}()
slcts[1]    = result_disks 
slcts[2]    = result_int 
slcts[3]    = result_ell
slcts[4]    = result_disk_nfakes
slcts[5]    = result_int_nfakes 
slcts[6]    = result_ell_nfakes 
labels      = Dict{Int64, String}()
labels[1]   = "d f"
labels[2]   = "i f"
labels[3]   = "e f"
labels[4]   = "d n"
labels[5]   = "i n"
labels[6]   = "e n"
plot_hexbins_nxm( log10.(as["M"] ./ as["M_MM"]), 
            as["ϕ_flip"], 
            slcts, labels, 3, 2,
            #scale=0.0,
            cbarshift=2.5,
            outfile="/home/moon/sfortune/spinevo/plots/correlations/flipangle_vs_MergerSTARS.png", 
            lognorm=true, gridres=(30,10), 
            xmin=0, xmax=3.5, ymin=0, ymax=180,
            ylabel="Flip [°]",
            xlabel="log( Stellar Merger Ratio )", 
            binned_median="x",
            )






# Phase Space plot
res = 5
plot_phasespace(; #snapNR = 128, 
    felixID=" ", rootID=3212, subID=" " , ymin=-1000, ymax=1000,
    parttype="STARS", r_by_rvir=0.1, #weights="MASS" ,
    lognorm=true, scale=15, gridres=(25*res, 15*res), colorlimits=(0.7,300), 
    outdir="/home/moon/sfortune/spinevo/plots/phase_space/halo_3212", outfile="PhSp", 
    )





# Plot halo 3212
halo = 3212
snap = 36
prop = "vel"
plot_group(snap; rootID=halo, 
property=prop, stepsize=10, arsize=5, ptsize=5000, rad=200, res=(1920,1080), spin=false, 
bgcolor=:black, indir=current_dir_stories, boxfix=false, 
file="/home/moon/sfortune/spinevo/plots/halo3212/halo$(halo)_snap$(snap)_$(prop).html")
