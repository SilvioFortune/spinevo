
include("/home/moon/sfortune/spinevo/pkg/meta.jl")

as      = load("/home/moon/sfortune/spinevo/data/assembly_softfix.jld", "assembly_STARS")
ad      = load("/home/moon/sfortune/spinevo/data/assembly_softfix.jld", "assembly_DM")
as09    = load("/home/moon/sfortune/spinevo/data/assembly_softfix_09Gyr.jld", "assembly_STARS")
ad09    = load("/home/moon/sfortune/spinevo/data/assembly_softfix_09Gyr.jld", "assembly_DM")

start_thr   = findcs(log10.(abs.(as["M"] .- as["ΔM"])), geq=log10.(2e10))
# borders: snap52=z1.18, snap100=z0.42
z_hi        = findcs(as["snapNR"      ], leq=52, comparewith=start_thr)
z_md        = findcs(as["snapNR"      ], gt=52, lt=100, comparewith=start_thr)
z_lo        = findcs(as["snapNR"      ], geq=100, comparewith=start_thr)
snap136     = findcs(as["snapNR"      ], eq=136, comparewith=start_thr)
snap36      = findcs(as["snapNR"      ], eq=36, comparewith=start_thr)

halo3212    = findcs(as["ID_ISUB"      ], eq=3212, comparewith=start_thr)
halo1414    = findcs(as["ID_ISUB"      ], eq=1414, comparewith=start_thr)

switches    = findcs(as["switch"      ], eq=1, comparewith=start_thr)
n_switches  = findcs(as["switch"      ], eq=0, comparewith=start_thr)
fakes       = findcs(as["FAKEFLIP"    ], eq=1, comparewith=start_thr)
n_fakes     = findcs(as["FAKEFLIP"    ], eq=0, comparewith=start_thr)

result_disk     = findcs(as["BVAL"  ], geq=b_disk, comparewith=start_thr)
result_int      = findcs(as["BVAL"  ], gt=b_ell, lt=b_disk, comparewith=start_thr)
result_ell      = findcs(as["BVAL"  ], leq=b_ell, comparewith=start_thr)

start_disk      = findcs(as["BVAL_0"  ] .- as["ΔBVAL_0"    ], geq=b_disk, comparewith=start_thr)
start_int       = findcs(as["BVAL_0"  ] .- as["ΔBVAL_0"    ], gt=b_ell, lt=b_disk, comparewith=start_thr)
start_ell       = findcs(as["BVAL_0"  ] .- as["ΔBVAL_0"    ], leq=b_ell, comparewith=start_thr)
stay_disk       = start_disk[  start_disk  .∈  Ref(Set(result_disk     ))]
stay_int        = start_int[    start_int   .∈  Ref(Set(result_int      ))]
stay_ell        = start_ell[    start_ell   .∈  Ref(Set(result_ell      ))]
result_disk_fakes    = fakes[   fakes   .∈  Ref(Set(result_disk     ))]
result_int_fakes     = fakes[   fakes   .∈  Ref(Set(result_int      ))]
result_ell_fakes     = fakes[   fakes   .∈  Ref(Set(result_ell      ))]
result_disk_nfakes   = n_fakes[ n_fakes .∈  Ref(Set(result_disk     ))]
result_int_nfakes    = n_fakes[ n_fakes .∈  Ref(Set(result_int      ))]
result_ell_nfakes    = n_fakes[ n_fakes .∈  Ref(Set(result_ell      ))]
start_disk_fakes    = fakes[   fakes   .∈  Ref(Set(start_disk     ))]
start_int_fakes     = fakes[   fakes   .∈  Ref(Set(start_int      ))]
start_ell_fakes     = fakes[   fakes   .∈  Ref(Set(start_ell      ))]
start_disk_nfakes   = n_fakes[ n_fakes .∈  Ref(Set(start_disk     ))]
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

mergersDM_1to50         = findcs(ad["M2"]./ad["Mpeak_MM"], geq=1, leq=50, comparewith=start_thr)

mergerSTARS_acc_g10     = findcs(         as["Mpeak_MERGERS"] ./ as["M2"],       geq=0.1, comparewith=start_thr)
mergerSTARS_acc_g5_l10  = findcs(         as["Mpeak_MERGERS"] ./ as["M2"],       geq=0.05, lt=0.1, comparewith=start_thr)
mergerSTARS_acc_g0_l5   = findcs(         as["Mpeak_MERGERS"] ./ as["M2"],       gt=0.0, lt=0.05, comparewith=start_thr)
mergerSTARS_acc_eq0     = findcs(replace( as["Mpeak_MERGERS"], missing => 0.0),  eq=0.0, comparewith=start_thr)











M9_start_thr   = findcs(log10.(abs.(as09["M"] .- as09["ΔM"])), geq=log10.(2e10))
# borders: snap52=z1.18, snap100=z0.42
M9_z_hi        = findcs(as09["snapNR"      ], leq=52, comparewith=M9_start_thr)
M9_z_md        = findcs(as09["snapNR"      ], gt=52, lt=100, comparewith=M9_start_thr)
M9_z_lo        = findcs(as09["snapNR"      ], geq=100, comparewith=M9_start_thr)
M9_snap136     = findcs(as09["snapNR"      ], eq=136, comparewith=M9_start_thr)
M9_snap36      = findcs(as09["snapNR"      ], eq=36, comparewith=M9_start_thr)


M9_halo3212    = findcs(as09["ID_ISUB"      ], eq=3212, comparewith=M9_start_thr)
M9_halo1414    = findcs(as09["ID_ISUB"      ], eq=1414, comparewith=M9_start_thr)

M9_switches    = findcs(as09["switch"      ], eq=1, comparewith=M9_start_thr)
M9_n_switches  = findcs(as09["switch"      ], eq=0, comparewith=M9_start_thr)
M9_fakes       = findcs(as09["FAKEFLIP"    ], eq=1, comparewith=M9_start_thr)
M9_n_fakes     = findcs(as09["FAKEFLIP"    ], eq=0, comparewith=M9_start_thr)

M9_result_disk     = findcs(as09["BVAL"  ], geq=b_disk, comparewith=M9_start_thr)
M9_result_int      = findcs(as09["BVAL"  ], gt=b_ell, lt=b_disk, comparewith=M9_start_thr)
M9_result_ell      = findcs(as09["BVAL"  ], leq=b_ell, comparewith=M9_start_thr)

M9_start_disk      = findcs(as09["BVAL_0"  ] .- as09["ΔBVAL_0"    ], geq=b_disk, comparewith=M9_start_thr)
M9_start_int       = findcs(as09["BVAL_0"  ] .- as09["ΔBVAL_0"    ], gt=b_ell, lt=b_disk, comparewith=M9_start_thr)
M9_start_ell       = findcs(as09["BVAL_0"  ] .- as09["ΔBVAL_0"    ], leq=b_ell, comparewith=M9_start_thr)
M9_stay_disk       = M9_start_disk[  M9_start_disk  .∈  Ref(Set(M9_result_disk     ))]
M9_stay_int        = M9_start_int[    M9_start_int   .∈  Ref(Set(M9_result_int      ))]
M9_stay_ell        = M9_start_ell[    M9_start_ell   .∈  Ref(Set(M9_result_ell      ))]
M9_result_disk_fakes    = M9_fakes[   M9_fakes   .∈  Ref(Set(M9_result_disk     ))]
M9_result_int_fakes     = M9_fakes[   M9_fakes   .∈  Ref(Set(M9_result_int      ))]
M9_result_ell_fakes     = M9_fakes[   M9_fakes   .∈  Ref(Set(M9_result_ell      ))]
M9_result_disk_nfakes   = M9_n_fakes[ M9_n_fakes .∈  Ref(Set(M9_result_disk     ))]
M9_result_int_nfakes    = M9_n_fakes[ M9_n_fakes .∈  Ref(Set(M9_result_int      ))]
M9_result_ell_nfakes    = M9_n_fakes[ M9_n_fakes .∈  Ref(Set(M9_result_ell      ))]
M9_start_disk_fakes    = M9_fakes[   M9_fakes   .∈  Ref(Set(M9_start_disk     ))]
M9_start_int_fakes     = M9_fakes[   M9_fakes   .∈  Ref(Set(M9_start_int      ))]
M9_start_ell_fakes     = M9_fakes[   M9_fakes   .∈  Ref(Set(M9_start_ell      ))]
M9_start_disk_nfakes   = M9_n_fakes[ M9_n_fakes .∈  Ref(Set(M9_start_disk     ))]
M9_start_int_nfakes    = M9_n_fakes[ M9_n_fakes .∈  Ref(Set(M9_start_int      ))]
M9_start_ell_nfakes    = M9_n_fakes[ M9_n_fakes .∈  Ref(Set(M9_start_ell      ))]
M9_stay_disk_fakes     = M9_fakes[   M9_fakes   .∈  Ref(Set(M9_stay_disk     ))]
M9_stay_int_fakes      = M9_fakes[   M9_fakes   .∈  Ref(Set(M9_stay_int      ))]
M9_stay_ell_fakes      = M9_fakes[   M9_fakes   .∈  Ref(Set(M9_stay_ell      ))]
M9_stay_disk_nfakes    = M9_n_fakes[ M9_n_fakes .∈  Ref(Set(M9_stay_disk     ))]
M9_stay_int_nfakes     = M9_n_fakes[ M9_n_fakes .∈  Ref(Set(M9_stay_int      ))]
M9_stay_ell_nfakes     = M9_n_fakes[ M9_n_fakes .∈  Ref(Set(M9_stay_ell      ))]

M9_ell_to_disk_nfakes  = M9_start_ell_nfakes[M9_start_ell_nfakes     .∈ Ref(Set(M9_result_disk_nfakes))]
M9_ell_to_int_nfakes   = M9_start_ell_nfakes[M9_start_ell_nfakes     .∈ Ref(Set(M9_result_int_nfakes))]
M9_int_to_disk_nfakes  = M9_start_int_nfakes[M9_start_int_nfakes     .∈ Ref(Set(M9_result_disk_nfakes))]
M9_int_to_ell_nfakes   = M9_start_int_nfakes[M9_start_int_nfakes     .∈ Ref(Set(M9_result_ell_nfakes))]
M9_disk_to_ell_nfakes  = M9_start_disk_nfakes[M9_start_disk_nfakes   .∈ Ref(Set(M9_result_ell_nfakes))]
M9_disk_to_int_nfakes  = M9_start_disk_nfakes[M9_start_disk_nfakes   .∈ Ref(Set(M9_result_int_nfakes))]

M9_mergersDM_1to50         = findcs(ad09["M2"]./ad09["Mpeak_MM"], geq=1, leq=50, comparewith=M9_start_thr)

M9_mergerSTARS_acc_g10     = findcs(         as09["Mpeak_MERGERS"] ./ as09["M2"],       geq=0.1, comparewith=M9_start_thr)
M9_mergerSTARS_acc_g5_l10  = findcs(         as09["Mpeak_MERGERS"] ./ as09["M2"],       geq=0.05, lt=0.1, comparewith=M9_start_thr)
M9_mergerSTARS_acc_g0_l5   = findcs(         as09["Mpeak_MERGERS"] ./ as09["M2"],       gt=0.0, lt=0.05, comparewith=M9_start_thr)
M9_mergerSTARS_acc_eq0     = findcs(replace( as09["Mpeak_MERGERS"], missing => 0.0),  eq=0.0, comparewith=M9_start_thr)



mergerDM_acc_g10     = findcs(         ad["Mpeak_MERGERS"] ./ ad["M2"],       geq=0.1, comparewith=start_thr)
mergerDM_acc_g5_l10  = findcs(         ad["Mpeak_MERGERS"] ./ ad["M2"],       geq=0.05, lt=0.1, comparewith=start_thr)
mergerDM_acc_g0_l5   = findcs(         ad["Mpeak_MERGERS"] ./ ad["M2"],       gt=0.0, lt=0.05, comparewith=start_thr)
mergerDM_acc_eq0     = findcs(replace( ad["Mpeak_MERGERS"], missing => 0.0),  eq=0.0, comparewith=start_thr)
M9_mergerDM_acc_g10     = findcs(         ad09["Mpeak_MERGERS"] ./ ad09["M2"],       geq=0.1, comparewith=M9_start_thr)
M9_mergerDM_acc_g5_l10  = findcs(         ad09["Mpeak_MERGERS"] ./ ad09["M2"],       geq=0.05, lt=0.1, comparewith=M9_start_thr)
M9_mergerDM_acc_g0_l5   = findcs(         ad09["Mpeak_MERGERS"] ./ ad09["M2"],       gt=0.0, lt=0.05, comparewith=M9_start_thr)
M9_mergerDM_acc_eq0     = findcs(replace( ad09["Mpeak_MERGERS"], missing => 0.0),  eq=0.0, comparewith=M9_start_thr)

M9_sfr_hi     = findcs(         log10.(as09["SFR"]),       geq=0.5, comparewith=M9_start_thr)
M9_sfr_lo     = findcs(         log10.(as09["SFR"]),       lt=0.5, comparewith=M9_start_thr)


z_welker    = findcs(as["redshift"    ], geq=1.2, leq=3.8, comparewith=start_thr)
z_lower     = findcs(as["redshift"    ], lt=1.2, comparewith=start_thr)
z_higher    = findcs(as["redshift"    ], gt=3.8, comparewith=start_thr)
M9_z_welker    = findcs(as09["redshift"    ], geq=1.2, leq=3.8, comparewith=M9_start_thr)
M9_z_lower     = findcs(as09["redshift"    ], lt=1.2, comparewith=M9_start_thr)
M9_z_higher    = findcs(as09["redshift"    ], gt=3.8, comparewith=M9_start_thr)



z_prepeak   = findcs(as["redshift"    ], geq=1.9, comparewith=start_thr)
z_postpeak  = findcs(as["redshift"    ], lt=1.9, comparewith=start_thr)
M9_z_prepeak   = findcs(as09["redshift"    ], geq=1.9, comparewith=M9_start_thr)
M9_z_postpeak  = findcs(as09["redshift"    ], lt=1.9, comparewith=M9_start_thr)
z_pre047   = findcs(as["redshift"    ], geq=0.47, comparewith=start_thr)
z_post047  = findcs(as["redshift"    ], lt=0.47, comparewith=start_thr)
M9_z_pre047   = findcs(as09["redshift"    ], geq=0.47, comparewith=M9_start_thr)
M9_z_post047  = findcs(as09["redshift"    ], lt=0.47, comparewith=M9_start_thr)
z_pre064   = findcs(as["redshift"    ], geq=0.64, comparewith=start_thr)
z_post064  = findcs(as["redshift"    ], lt=0.64, comparewith=start_thr)
M9_z_pre064   = findcs(as09["redshift"    ], geq=0.64, comparewith=M9_start_thr)
M9_z_post064  = findcs(as09["redshift"    ], lt=0.64, comparewith=M9_start_thr)



# new
j_hi    = findcs(log10.(pyplottable(as["j_main"]  .- as["Δj_main"] )), gt=2.55, comparewith=start_thr)
j_md    = findcs(log10.(pyplottable(as["j_main"]  .- as["Δj_main"] )), gt=2.2 ,leq=2.55, comparewith=start_thr)
j_lo    = findcs(log10.(pyplottable(as["j_main"]  .- as["Δj_main"] )), leq=2.2, comparewith=start_thr)

J_hi    = findcs(log10.(pyplottable(as["J_main"]  .- as["ΔJ_main"] )), gt=12.9, comparewith=start_thr)
J_md    = findcs(log10.(pyplottable(as["J_main"]  .- as["ΔJ_main"] )), gt=12.6 ,leq=12.9, comparewith=start_thr)
J_lo    = findcs(log10.(pyplottable(as["J_main"]  .- as["ΔJ_main"] )), leq=12.6, comparewith=start_thr)

z_hirsch_0  = findcs(as["redshift"], lt=0.25)
z_hirsch_05 = findcs(as["redshift"], geq=0.25, lt=0.75)
z_hirsch_1  = findcs(as["redshift"], geq=0.75, lt=1.5)
z_hirsch_2  = findcs(as["redshift"], geq=1.5, lt=2.5)
z_hirsch_3  = findcs(as["redshift"], geq=2.5, lt=3.5)
z_hirsch_4  = findcs(as["redshift"], geq=3.5)








# SFR vs Flipangle
x   = as["ϕ_flip"]
y   = log10.(as["SFR"])# ./ as["M"]
miny= minimum(pyplottable(log10.(as["SFR"]))[findall(x->x .!== NaN, pyplottable(log10.(as["SFR"])))])
maxy= maximum(pyplottable(log10.(as["SFR"]))[findall(x->x .!== NaN, pyplottable(log10.(as["SFR"])))])
#slct= n_fakes[ n_fakes .∈  Ref(Set(mergersDM_1to50     ))]
slct = stay_ell
plot_hexbins(   x, y,
            selection = slct,
            outfile="/home/moon/sfortune/spinevo/plots/correlations/SFR_vs_flipangle_ETGs.png",
            lognorm=true, gridres=(20,10), grid=true,
            ymin=miny, ymax=maxy, 
            xmin=0, xmax=180, 
            xlabel="Flip Angle [°]",
            ylabel="log( SFR [M_{\odot}/Gyr] )", 
            title="Ellipticals",
            binned_median="x", 
            )




# b-value vs SFR at different redshift bins
x   = as["BVAL_0"]
y   = log10.(as["SFR"] ./ as["M"] )
miny= minimum(pyplottable(y)[findall(x->x .!== NaN, pyplottable(y))])
maxy= maximum(pyplottable(y)[findall(x->x .!== NaN, pyplottable(y))])
minx= minimum(pyplottable(x)[findall(x->x .!== NaN, pyplottable(x))])
maxx= maximum(pyplottable(x)[findall(x->x .!== NaN, pyplottable(x))])
slcts       = Dict{Int64, Vector{Int64}}()
slcts[1]    = start_thr
slcts[2]    = z_md
slcts[3]    = z_hi 
slcts[4]    = z_lo
labels      = Dict{Int64, String}()
labels[1]   = "Full Set, N = $(length(slcts[1]))"
labels[2]   = "$(round(maximum(as["redshift"][slcts[2]]), digits=2)) $(latexstring("\$\\geq\$")) z $(latexstring("\$\\geq\$")) $(round(minimum(as["redshift"][slcts[2]]), digits=2)), N = $(length(slcts[2]))"
labels[3]   = "$(round(maximum(as["redshift"][slcts[3]]), digits=2)) $(latexstring("\$\\geq\$")) z $(latexstring("\$\\geq\$")) $(round(minimum(as["redshift"][slcts[3]]), digits=2)), N = $(length(slcts[3]))"
labels[4]   = "$(round(maximum(as["redshift"][slcts[4]]), digits=2)) $(latexstring("\$\\geq\$")) z $(latexstring("\$\\geq\$")) $(round(minimum(as["redshift"][slcts[4]]), digits=2)), N = $(length(slcts[4]))"
plot_hexbins_nxm( x, 
            y, 
            slcts, labels, 2, 2,
            #scale=0.0,
            cbarshift=2.5,
            outfile="/home/moon/sfortune/spinevo/plots/correlations/sSFR_vs_bval_zsplit.png", 
            lognorm=true, gridres=(30,10), 
            ymin=miny, ymax=maxy,
            xmin=minx, xmax=maxx, 
            ylabel="log( sSFR [$(latexstring("\\frac{1}{yr}"))] )",
            xlabel="b-value", fontsize=30,
            binned_median="x", plot_bval="x",
            )

# fitting to the above
linfct(x,p) = p[1] .+ p[2] .* (x .+ 4.5)
sane_hi = findall(x->x .!== NaN, pyplottable(as["BVAL_0"][slcts[3]]))[findall(x->x .!== NaN, pyplottable(as["BVAL_0"][slcts[3]])) .∈ Ref(Set(findall(x->x .!== NaN, pyplottable(log10.(as["SFR"] ./ as["M"] )[slcts[3]]))))]
sane_md = findall(x->x .!== NaN, pyplottable(as["BVAL_0"][slcts[2]]))[findall(x->x .!== NaN, pyplottable(as["BVAL_0"][slcts[2]])) .∈ Ref(Set(findall(x->x .!== NaN, pyplottable(log10.(as["SFR"] ./ as["M"] )[slcts[2]]))))]
sane_lo = findall(x->x .!== NaN, pyplottable(as["BVAL_0"][slcts[4]]))[findall(x->x .!== NaN, pyplottable(as["BVAL_0"][slcts[4]])) .∈ Ref(Set(findall(x->x .!== NaN, pyplottable(log10.(as["SFR"] ./ as["M"] )[slcts[4]]))))]
xdata_hi =                  as["BVAL_0"][slcts[3]][sane_hi]
ydata_hi = log10.(as["SFR"] ./ as["M"] )[slcts[3]][sane_hi]
xdata_md =                  as["BVAL_0"][slcts[2]][sane_md]
ydata_md = log10.(as["SFR"] ./ as["M"] )[slcts[2]][sane_md]
xdata_lo =                  as["BVAL_0"][slcts[4]][sane_lo]
ydata_lo = log10.(as["SFR"] ./ as["M"] )[slcts[4]][sane_lo]
p0 = [-10.5, 1.0]
fit_hi  = curve_fit(linfct, xdata_hi, ydata_hi, p0)
fit_md  = curve_fit(linfct, xdata_md, ydata_md, p0)
fit_lo  = curve_fit(linfct, xdata_lo, ydata_lo, p0)

# Plot it
PyPlot.matplotlib[:rc]("text", usetex=true)
PyPlot.matplotlib[:rc]("font", family="serif", serif=["Computer Modern"])
pointsize = 5
scale = 0.7
fig, ax = subplots()
ax.axvline(b_disk, color="blue")
ax.axvline(b_ell, color="red")
ax.fill_between(xdata_lo, linfct(xdata_lo, fit_lo.param .+ diag(estimate_covar(fit_lo)) .* [1,-1]), linfct(xdata_lo, fit_lo.param .+ diag(estimate_covar(fit_lo)) .* [-1,1]), color="mediumblue", alpha=0.2)
ax.fill_between(xdata_md, linfct(xdata_md, fit_md.param .+ diag(estimate_covar(fit_md)) .* [1,-1]), linfct(xdata_md, fit_md.param .+ diag(estimate_covar(fit_md)) .* [-1,1]), color="darkorchid", alpha=0.2)
ax.fill_between(xdata_hi, linfct(xdata_hi, fit_hi.param .+ diag(estimate_covar(fit_hi)) .* [1,-1]), linfct(xdata_hi, fit_hi.param .+ diag(estimate_covar(fit_hi)) .* [-1,1]), color="darkred"   , alpha=0.2)
ax.plot( xdata_lo, linfct(xdata_lo, fit_lo.param), "-", color="mediumblue", label="$(labels[4])")
ax.plot( xdata_md, linfct(xdata_md, fit_md.param), "-", color="darkorchid", label="$(labels[2])")
ax.plot( xdata_hi, linfct(xdata_hi, fit_hi.param), "-", color="darkred"   , label="$(labels[3])")
ax.fill_between(xdata_lo, linfct(xdata_lo, fit_lo.param .+ diag(estimate_covar(fit_lo)) .* [1,-1]), linfct(xdata_lo, fit_lo.param .+ diag(estimate_covar(fit_lo)) .* [-1,1]), color="mediumblue", alpha=0.2)
ax.fill_between(xdata_md, linfct(xdata_md, fit_md.param .+ diag(estimate_covar(fit_md)) .* [1,-1]), linfct(xdata_md, fit_md.param .+ diag(estimate_covar(fit_md)) .* [-1,1]), color="darkorchid", alpha=0.2)
ax.fill_between(xdata_hi, linfct(xdata_hi, fit_hi.param .+ diag(estimate_covar(fit_hi)) .* [1,-1]), linfct(xdata_hi, fit_hi.param .+ diag(estimate_covar(fit_hi)) .* [-1,1]), color="darkred"   , alpha=0.2)
ax.set_title("$(latexstring("\$\\Delta\$"))t $(latexstring("\$\\sim\$")) 0.3 Gyr")
ax.set_xlabel("b-value")
ax.set_ylabel("log( sSFR [$(latexstring("\\frac{1}{yr}"))] )")
ax.legend(loc="upper left", frameon=true, borderpad=1, handlelength=1.8)
ax.grid()
fig.set_size_inches(9scale, 9scale)
fig.tight_layout()
fig.savefig("/home/moon/sfortune/spinevo/plots/correlations/sSFR_vs_bval_zsplit_linfits.png", bbox_inches="tight", pad_inches=.1)
close("all")
PyPlot.matplotlib[:rc]("text", usetex=false)







# Welker 2014

e0 = mergerSTARS_acc_eq0#[   mergerSTARS_acc_eq0   .∈  Ref(Set(z_welker     ))]
g0 = mergerSTARS_acc_g0_l5#[   mergerSTARS_acc_g0_l5   .∈  Ref(Set(z_welker     ))]
g5 = mergerSTARS_acc_g5_l10#[   mergerSTARS_acc_g5_l10   .∈  Ref(Set(z_welker     ))]
g10= mergerSTARS_acc_g10#[   mergerSTARS_acc_g10   .∈  Ref(Set(z_welker     ))]
scale=5
N_bins= 5
n_e0, x_e0      = PyPlot.hist( cosd.(pyplottable(as["ϕ_flip"][e0])), 
    range=(0,1),#maximum(cosd.(mergers_e0["ϕ_flip"]))), 
    bins=N_bins, density=true, log=true, rwidth=0.9, histtype="step", alpha=1, zorder=1)#, edgecolor="black")
n_g0, x_g0      = PyPlot.hist( cosd.(pyplottable(as["ϕ_flip"][g0])), 
    range=(0,1),#maximum(cosd.(mergers_g0["ϕ_flip"]))), 
    bins=N_bins, density=true, log=true, rwidth=0.9, histtype="step", alpha=1, zorder=1)#, edgecolor="black")
n_g5, x_g5      = PyPlot.hist( cosd.(pyplottable(as["ϕ_flip"][g5])), 
    range=(0,1),#maximum(cosd.(mergers_g5["ϕ_flip"]))), 
    bins=N_bins, density=true, log=true, rwidth=0.9, histtype="step", alpha=1, zorder=1)#, edgecolor="black")
n_g10, x_g10    = PyPlot.hist( cosd.(pyplottable(as["ϕ_flip"][g10])), 
    range=(0,1),#maximum(cosd.(mergers_g10["ϕ_flip"]))), 
    bins=N_bins, density=true,log=true,rwidth=0.9, histtype="step", alpha=1, zorder=1)#, edgecolor="black")
bcenters_e0     = 0.5 .* (x_e0[2:end]  .+ x_e0[1:end-1])
bcenters_g0     = 0.5 .* (x_g0[2:end]  .+ x_g0[1:end-1])
bcenters_g5     = 0.5 .* (x_g5[2:end]  .+ x_g5[1:end-1])
bcenters_g10    = 0.5 .* (x_g10[2:end] .+ x_g10[1:end-1])
clf()
PyPlot.matplotlib[:rc]("text", usetex=true)
PyPlot.matplotlib[:rc]("font", family="serif", serif=["Computer Modern"])
fig, ax = subplots()
ax.plot(bcenters_e0 , n_e0 , ".-", ms=10, label="$(latexstring("\$\\delta\$"))m = 0\\%,       N = $(length(e0))" , color="red"   )
ax.plot(bcenters_g0 , n_g0 , ".-", ms=10, label="5\\% $(latexstring("\$\\geq\$")) $(latexstring("\$\\delta\$"))m $(latexstring("\$\\geq\$")) 0\\%,  N = $(length(g0))" , color="orange")
ax.plot(bcenters_g5 , n_g5 , ".-", ms=10, label="10\\% $(latexstring("\$\\geq\$")) $(latexstring("\$\\delta\$"))m $(latexstring("\$\\geq\$")) 5\\%, N = $(length(g5))" , color="purple")
ax.plot(bcenters_g10, n_g10, ".-", ms=10, label="$(latexstring("\$\\delta\$"))m $(latexstring("\$\\geq\$")) 10\\%,      N = $(length(g10))" , color="blue"  )
ax.set_xlabel("cos($(latexstring("\$\\phi\$")))")
ax.set_ylabel("P")
ax.set_yscale("log")
ax.set_ylim(1e-3,)
ax.set_title("$(latexstring("\$\\Delta\$"))t $(latexstring("\$\\sim\$")) 0.3 Gyr within $(round(maximum(as["redshift"]),digits=2)) $(latexstring("\$\\geq\$")) z $(latexstring("\$\\geq\$")) $(round(minimum(as["redshift"]),digits=2))")
ax.legend(loc="upper left", frameon=true, borderpad=1, handlelength=1.8)
ax.grid()
fig.set_size_inches(scale, scale)
fig.tight_layout()
fig.savefig("/home/moon/sfortune/spinevo/plots/correlations/flipPDF_by_mergerRatioSTARS_Welker2014_03_bins5_zall.png", bbox_inches="tight", pad_inches=.1)
close("all")
PyPlot.matplotlib[:rc]("text", usetex=false)




# Welker 2014 0.9 Gyr


scale=5
N_bins= 10
n_e0, x_e0      = PyPlot.hist( cosd.(pyplottable(as09["ϕ_flip"][M9_mergerDM_acc_eq0])), 
    range=(0,1),#maximum(cosd.(mergers_e0["ϕ_flip"]))), 
    bins=N_bins, density=true, log=true, rwidth=0.9, histtype="step", alpha=1, zorder=1)#, edgecolor="black")
n_g0, x_g0      = PyPlot.hist( cosd.(pyplottable(as09["ϕ_flip"][M9_mergerDM_acc_g0_l5])), 
    range=(0,1),#maximum(cosd.(mergers_g0["ϕ_flip"]))), 
    bins=N_bins, density=true, log=true, rwidth=0.9, histtype="step", alpha=1, zorder=1)#, edgecolor="black")
n_g5, x_g5      = PyPlot.hist( cosd.(pyplottable(as09["ϕ_flip"][M9_mergerDM_acc_g5_l10])), 
    range=(0,1),#maximum(cosd.(mergers_g5["ϕ_flip"]))), 
    bins=N_bins, density=true, log=true, rwidth=0.9, histtype="step", alpha=1, zorder=1)#, edgecolor="black")
n_g10, x_g10    = PyPlot.hist( cosd.(pyplottable(as09["ϕ_flip"][M9_mergerDM_acc_g10])), 
    range=(0,1),#maximum(cosd.(mergers_g10["ϕ_flip"]))), 
    bins=N_bins, density=true,log=true,rwidth=0.9, histtype="step", alpha=1, zorder=1)#, edgecolor="black")
bcenters_e0     = 0.5*(x_e0[2:end]  + x_e0[1:end-1])
bcenters_g0     = 0.5*(x_g0[2:end]  + x_g0[1:end-1])
bcenters_g5     = 0.5*(x_g5[2:end]  + x_g5[1:end-1])
bcenters_g10    = 0.5*(x_g10[2:end] + x_g10[1:end-1])
clf()
fig, ax = subplots()
ax.plot(bcenters_e0 , n_e0 , ".-", ms=10, label="δm = 0%,       N = $(length(M9_mergerDM_acc_eq0))" , color="red"   )
ax.plot(bcenters_g0 , n_g0 , ".-", ms=10, label="5% > δm > 0%,  N = $(length(M9_mergerDM_acc_g0_l5))" , color="orange")
ax.plot(bcenters_g5 , n_g5 , ".-", ms=10, label="10% > δm > 5%, N = $(length(M9_mergerDM_acc_g5_l10))" , color="purple")
ax.plot(bcenters_g10, n_g10, ".-", ms=10, label="δm > 10%,      N = $(length(M9_mergerDM_acc_g10))" , color="blue"  )
ax.set_xlabel("cos(ϕ)")
ax.set_ylabel("P")
ax.set_yscale("log")
ax.set_ylim(1e-3,)
ax.set_title("1 Gyr Time Steps")
ax.legend(loc="upper left", frameon=true, borderpad=1, handlelength=1.8)
ax.grid()
fig.set_size_inches(scale, scale)
fig.tight_layout()
fig.savefig("/home/moon/sfortune/spinevo/plots/correlations/flipPDF_by_mergerRatioDM_Welker2014_09.png", bbox_inches="tight", pad_inches=.1)





# mass diff DM vs STARS with color in flip 0.9 Gyr
PyPlot.matplotlib[:rc]("text", usetex=true)
PyPlot.matplotlib[:rc]("font", family="serif", serif=["Computer Modern"])
pointsize = 10
scale = 0.5
fig, ax = subplots()
#y       = pyplottable(log10.(abs.(ad09["ΔM"])))
#x       = pyplottable(log10.(abs.(as09["ΔM"])))
y       = pyplottable(log10.(abs.(ad09["ΔM"]) ./ (ad09["M"] .- ad09["ΔM"])))
x       = pyplottable(log10.(abs.(as09["ΔM"]) ./ (as09["M"] .- as09["ΔM"])))
#colr    = log10.(pyplottable(as09["ϕ_flip"]./as09["Δlookbacktime"] ./-1000)) .+ abs(minimum(log10.(pyplottable(as09["ϕ_flip"]./as09["Δlookbacktime"] ./-1000))[findall(x->x .!== NaN, log10.(pyplottable(as["ϕ_flip"]./as["Δlookbacktime"] ./-1000)))]))
colr    = pyplottable(as09["ϕ_flip"]./as09["Δlookbacktime"] ./-1000)
maxcolr = maximum(colr[findall(x->!isnan.(x), colr)])
slct    = findall(x->!isnan.(x), colr)
p = ax.scatter( x[slct], y[slct],
        zorder=2   , s= (pointsize/maxcolr) .* colr[slct]    , alpha=0.6,#colr[slct] ./ maxcolr, 
        cmap="plasma_r", c=colr[slct]) # brg_r  
colbar  = fig.colorbar(p, ax=ax)#, ticks=[minimum(as["BVAL_0"    ]), b_ell, b_disk, maximum(as["BVAL_0"    ])])
#colbar.ax.set_yticklabels(["$(round(minimum(as["BVAL_0"    ]), digits=1))", "ETG", "LTG", "$(round(maximum(as["BVAL_0"    ]), digits=1))"])
colbar.ax.tick_params(axis="y")
#ax.legend(loc="lower right", frameon=true, borderpad=1, handlelength=1.8)
ax.set_title("1 Gyr Flip")
ax.set_xlabel("log( ΔM_STARS / M_DM )")
ax.set_ylabel("log( ΔM_DM / M_DM )")
ax.grid()
#ax.autoscale()
fig.set_size_inches(16scale, 9scale)
fig.tight_layout()
fig.savefig("/home/moon/sfortune/spinevo/plots/correlations/RELATIVE_delMdm_vs_delMstars_vs_flipbyMyr_09.png", bbox_inches="tight", pad_inches=.1)
close("all")
PyPlot.matplotlib[:rc]("text", usetex=false)




# delta M STARS vs M STARS with color in flip
PyPlot.matplotlib[:rc]("text", usetex=true)
PyPlot.matplotlib[:rc]("font", family="serif", serif=["Computer Modern"])
pointsize = 5
scale = 0.6
fig, ax = subplots()
y       = pyplottable(log10.(abs.(ad["ΔM"])))[findcs(log10.(abs.(as["M"] .- as["ΔM"])), geq=log10.(2e10))]
x       = pyplottable(log10.(abs.(as["ΔM"])))[findcs(log10.(abs.(as["M"] .- as["ΔM"])), geq=log10.(2e10))]
#y       = pyplottable(log10.(abs.(ad["ΔM"]) ./ (ad["M"] .- ad["ΔM"])))
#x       = pyplottable(log10.(abs.(as["ΔM"]) ./ (as["M"] .- as["ΔM"])))
colr    = pyplottable(as["ϕ_flip"])[findcs(log10.(abs.(as["M"] .- as["ΔM"])), geq=log10.(2e10))]#./as["Δlookbacktime"] ./-1000)
maxcolr = maximum(colr[findall(x->!isnan.(x), colr)])
slct    = findall(x->!isnan.(x), colr)
order   = sortperm(colr[slct])
p = ax.scatter( x[slct][order], y[slct][order],
        zorder=2   , s= (pointsize),#/maxcolr) .* colr[slct][order]    , 
        alpha=0.6,#colr[slct][order] ./ maxcolr, 
        cmap="plasma_r", c=colr[slct][order]) # brg_r  
colbar  = fig.colorbar(p, ax=ax, label="Flip [°]")#, ticks=[minimum(as["BVAL_0"    ]), b_ell, b_disk, maximum(as["BVAL_0"    ])])
#colbar.ax.set_yticklabels(["$(round(minimum(as["BVAL_0"    ]), digits=1))", "ETG", "LTG", "$(round(maximum(as["BVAL_0"    ]), digits=1))"])
colbar.ax.tick_params(axis="y")
#ax.legend(loc="lower right", frameon=true, borderpad=1, handlelength=1.8)
ax.set_title("$(latexstring("\$\\Delta\$"))t $(latexstring("\$\\sim\$")) 0.3 Gyr")
ax.set_xlabel("log( $(latexstring("\$\\Delta\$"))M_STARS )")#/ M_STARS )")
ax.set_ylabel("log( $(latexstring("\$\\Delta\$"))M_DM )")#/ M_DM )")
ax.grid()
#ax.autoscale()
fig.set_size_inches(13scale, 9scale)
fig.tight_layout()
fig.savefig("/home/moon/sfortune/spinevo/plots/correlations/delMdm_vs_delMstars_ABS_vs_flipangle_03.png", bbox_inches="tight", pad_inches=.1)
close("all")
PyPlot.matplotlib[:rc]("text", usetex=false)





# mass diff DM vs STARS with color in flip
PyPlot.matplotlib[:rc]("text", usetex=true)
PyPlot.matplotlib[:rc]("font", family="serif", serif=["Computer Modern"])
pointsize = 5
scale = 0.5
fig, ax = subplots()
y       = pyplottable(log10.(abs.(ad09["ΔM"])))
x       = pyplottable(log10.(abs.(as09["ΔM"])))
#y       = pyplottable(log10.(abs.(ad09["ΔM"]) ./ (ad09["M"] .- ad09["ΔM"])))
#x       = pyplottable(log10.(abs.(as09["ΔM"]) ./ (as09["M"] .- as09["ΔM"])))
colr    = pyplottable(as09["ϕ_flip"])#./as["Δlookbacktime"] ./-1000)
maxcolr = maximum(colr[findall(x->!isnan.(x), colr)])
slct    = findall(x->!isnan.(x), colr)
order   = sortperm(colr[slct])
p = ax.scatter( x[slct][order], y[slct][order],
        zorder=2   , s= (pointsize),#/maxcolr) .* colr[slct][order]    , 
        alpha=0.6,#colr[slct][order] ./ maxcolr, 
        cmap="plasma_r", c=colr[slct][order]) # brg_r  
colbar  = fig.colorbar(p, ax=ax, label="Flip [°]")#, ticks=[minimum(as["BVAL_0"    ]), b_ell, b_disk, maximum(as["BVAL_0"    ])])
#colbar.ax.set_yticklabels(["$(round(minimum(as["BVAL_0"    ]), digits=1))", "ETG", "LTG", "$(round(maximum(as["BVAL_0"    ]), digits=1))"])
colbar.ax.tick_params(axis="y")
#ax.legend(loc="lower right", frameon=true, borderpad=1, handlelength=1.8)
ax.set_title("$(latexstring("\$\\Delta\$"))t $(latexstring("\$\\sim\$")) 1 Gyr")
ax.set_xlabel("log( $(latexstring("\$\\Delta\$"))M_STARS )")#/ M_STARS )")
ax.set_ylabel("log( $(latexstring("\$\\Delta\$"))M_DM )")#/ M_DM )")
ax.set_ylabel("log( $(latexstring("\$\\Delta\$"))M_DM )")#/ M_DM )")
ax.grid()
#ax.autoscale()
fig.set_size_inches(16scale, 9scale)
fig.tight_layout()
fig.savefig("/home/moon/sfortune/spinevo/plots/correlations/delMdm_vs_delMstars_ABS_vs_flipangle_09.png", bbox_inches="tight", pad_inches=.1)
close("all")
PyPlot.matplotlib[:rc]("text", usetex=false)



# SFR
fig, ax = subplots()
#ax.hist( pyplottable(log10.(as09["SFR"]))    , bins=50, label="1 Gyr"    , rwidth=0.9    , alpha=0.5, color="red" )
ax.hist( pyplottable(log10.(as["SFR"]))    , bins=50, rwidth=0.9    , alpha=0.5, color="navy" )
ax.set_xlabel("log( SFR )")
#ax.set_xscale("log")
#ax.set_xlim([10, nothing])
#ax.legend(loc="upper right", frameon=true, borderpad=1, handlelength=1.8)
ax.grid()
scale=0.5
fig.set_size_inches(16scale, 9scale)
fig.tight_layout()
fig.savefig("/home/moon/sfortune/spinevo/plots/correlations/histogram_sfr.png", bbox_inches="tight", pad_inches=.1)




# delta j stars vs flip
x   = as["ϕ_flip"]
y   = log10.(as["Δj_main"])
miny= minimum(pyplottable(y)[findall(x->x .!== NaN, pyplottable(y))])
maxy= maximum(pyplottable(y)[findall(x->x .!== NaN, pyplottable(y))])
minx= minimum(pyplottable(x)[findall(x->x .!== NaN, pyplottable(x))])
maxx= maximum(pyplottable(x)[findall(x->x .!== NaN, pyplottable(x))])
slcts       = Dict{Int64, Vector{Int64}}()
slcts[1]    = 1:length(as["ϕ_flip"]) 
slcts[2]    = stay_disk
slcts[3]    = stay_ell 
slcts[4]    = z_lo
labels      = Dict{Int64, String}()
labels[1]   = "Full Set"
labels[2]   = "$(round(maximum(as["redshift"][slcts[2]]), digits=2)) $(latexstring("\$\\geq\$")) z $(latexstring("\$\\geq\$")) $(round(minimum(as["redshift"][slcts[2]]), digits=2)), N = $(length(slcts[2]))"
labels[3]   = "$(round(maximum(as["redshift"][slcts[3]]), digits=2)) $(latexstring("\$\\geq\$")) z $(latexstring("\$\\geq\$")) $(round(minimum(as["redshift"][slcts[3]]), digits=2)), N = $(length(slcts[3]))"
labels[4]   = "$(round(maximum(as["redshift"][slcts[4]]), digits=2)) $(latexstring("\$\\geq\$")) z $(latexstring("\$\\geq\$")) $(round(minimum(as["redshift"][slcts[4]]), digits=2)), N = $(length(slcts[4]))"
plot_hexbins_nxm( x, 
            y, 
            slcts, labels, 3, 2,
            #scale=0.0,
            cbarshift=2.5,
            outfile="/home/moon/sfortune/spinevo/plots/correlations/sSFR_vs_bval_zsplit.png", 
            lognorm=true, gridres=(30,10), 
            ymin=miny, ymax=maxy,
            xmin=minx, xmax=maxx, 
            ylabel="|Δj_STARS|",
            xlabel="Flip [°]", fontsize=30,
            binned_median="x", plot_bval="x",
            )




# j_hist
PyPlot.matplotlib[:rc]("text", usetex=true)
PyPlot.matplotlib[:rc]("font", family="serif", serif=["Computer Modern"])
scale = 0.5
fig, ax = subplots()
#ax.hist( pyplottable(log10.(as09["SFR"]))    , bins=50, label="1 Gyr"    , rwidth=0.9    , alpha=0.5, color="red" )
ax.hist( log10.(pyplottable(as["j_main"][:,start_ell]  .- as["Δj_main"][:,start_ell] ))    , bins=50, label="Ellipticals"    , rwidth=0.9    , alpha=0.5, color="darkred" )
ax.hist( log10.(pyplottable(as["j_main"][:,start_int]  .- as["Δj_main"][:,start_int] ))    , bins=50, label="Intermediates"    , rwidth=0.9    , alpha=0.5, color="darkorchid" )
ax.hist( log10.(pyplottable(as["j_main"][:,start_disk] .- as["Δj_main"][:,start_disk]))    , bins=50, label="Disks"    , rwidth=0.9    , alpha=0.5, color="mediumblue" )
ax.axvline( 2.2, color ="black")
ax.axvline( 2.55, color ="black")
ax.set_xlabel("log( j )")
#ax.set_xscale("log")
#ax.set_xlim([10, nothing])
ax.legend(loc="upper right", frameon=true, borderpad=1, handlelength=1.8)
ax.grid()
scale=0.5
fig.set_size_inches(16scale, 9scale)
fig.tight_layout()
fig.savefig("/home/moon/sfortune/spinevo/plots/correlations/histogram_j_start.png", bbox_inches="tight", pad_inches=.1)
close("all")
PyPlot.matplotlib[:rc]("text", usetex=false)












# delta J vs sSFR at different redshift bins
x   = log10.(pyplottable(as["ΔJ_main"]))
y   = log10.(as["SFR"] ./ as["M"] )
miny= minimum(pyplottable(y)[findall(x->x .!== NaN, pyplottable(y))])
maxy= maximum(pyplottable(y)[findall(x->x .!== NaN, pyplottable(y))])
minx= minimum(pyplottable(x)[findall(x->x .!== NaN, pyplottable(x))])
maxx= maximum(pyplottable(x)[findall(x->x .!== NaN, pyplottable(x))])
slcts       = Dict{Int64, Vector{Int64}}()
slcts[1]    = z_lo
slcts[2]    = z_lo[ z_lo .∈  Ref(Set(start_int_nfakes      ))]
slcts[3]    = z_lo[ z_lo .∈  Ref(Set(start_disk_nfakes      ))] 
slcts[4]    = z_lo[ z_lo .∈  Ref(Set(start_ell_nfakes      ))]
labels      = Dict{Int64, String}()
labels[1]   = "Full Set, N = $(length(slcts[1]))"
labels[2]   = "Intermediates, N = $(length(slcts[2]))"
labels[3]   = "Disks, N = $(length(slcts[3]))"
labels[4]   = "Ellipticals, N = $(length(slcts[4]))"
plot_hexbins_nxm( x, 
            y, 
            slcts, labels, 2, 2,
            #scale=0.0,
            cbarshift=2.5,
            outfile="/home/moon/sfortune/spinevo/plots/correlations/sSFR_vs_delJ_BVALsplit_zlow.png", 
            lognorm=true, gridres=(30,10), 
            ymin=miny, ymax=maxy,
            xmin=minx, xmax=maxx, 
            ylabel="log( sSFR [$(latexstring("\\frac{1}{yr}"))] )",
            xlabel="log( $(latexstring("\$\\Delta\$"))J_STARS )", fontsize=30,
            binned_median="x", #plot_bval="x",
            )
            





# The Proof flip by MM depends on j, binned by j
y   = as["ϕ_flip"]
x   = log10.( abs.(as["M2"] ./ as["Mpeak_MM"]) )
miny= minimum(pyplottable(y)[findall(x->x .!== NaN, pyplottable(y))])
maxy= maximum(pyplottable(y)[findall(x->x .!== NaN, pyplottable(y))])
minx= minimum(pyplottable(x)[findall(x->x .!== NaN, pyplottable(x))])
maxx= maximum(pyplottable(x)[findall(x->x .!== NaN, pyplottable(x))])
slcts       = Dict{Int64, Vector{Int64}}()
slcts[1]    = start_thr
slcts[2]    = j_md
slcts[3]    = j_lo
slcts[4]    = j_hi
labels      = Dict{Int64, String}()
labels[1]   = "Full Set N = $(length(slcts[1]))"
labels[2]   = "Medium j, N = $(length(slcts[2]))"
labels[3]   = "Low j, N = $(length(slcts[3]))"
labels[4]   = "High j, N = $(length(slcts[4]))"
plot_hexbins_nxm( x, y, 
            slcts, labels, 2, 2,
            #scale=0.0,
            cbarshift=2.5,
            outfile="/home/moon/sfortune/spinevo/plots/correlations/flip_vs_M_MMRatio_jsplit.png", 
            lognorm=true, gridres=(30,10), 
            ymin=miny, ymax=maxy,
            xmin=minx, xmax=maxx, 
            ylabel="Flip [°]",
            xlabel="log( Most Massive Merger Ratio )", fontsize=30,
            binned_median="x", plot_mergers="x",
            )





# J_hist
PyPlot.matplotlib[:rc]("text", usetex=true)
PyPlot.matplotlib[:rc]("font", family="serif", serif=["Computer Modern"])
scale = 0.5
fig, ax = subplots()
#ax.hist( pyplottable(log10.(as09["SFR"]))    , bins=50, label="1 Gyr"    , rwidth=0.9    , alpha=0.5, color="red" )
ax.hist( log10.(pyplottable(as["J_main"][:,start_ell]  .- as["ΔJ_main"][:,start_ell] ))    , bins=50, label="Ellipticals"    , rwidth=0.9    , alpha=0.5, color="darkred" )
ax.hist( log10.(pyplottable(as["J_main"][:,start_int]  .- as["ΔJ_main"][:,start_int] ))    , bins=50, label="Intermediates"    , rwidth=0.9    , alpha=0.5, color="darkorchid" )
ax.hist( log10.(pyplottable(as["J_main"][:,start_disk] .- as["ΔJ_main"][:,start_disk]))    , bins=50, label="Disks"    , rwidth=0.9    , alpha=0.5, color="mediumblue" )
ax.axvline( 12.6, color ="black")
ax.axvline( 12.9, color ="black")
ax.set_xlabel("log( J )")
#ax.set_xscale("log")
#ax.set_xlim([10, nothing])
ax.legend(loc="upper right", frameon=true, borderpad=1, handlelength=1.8)
ax.grid()
scale=0.5
fig.set_size_inches(16scale, 9scale)
fig.tight_layout()
fig.savefig("/home/moon/sfortune/spinevo/plots/correlations/histogram_J_start.png", bbox_inches="tight", pad_inches=.1)
close("all")
PyPlot.matplotlib[:rc]("text", usetex=false)












# The Proof flip by MM depends on J, binned by J
y   = as["ϕ_flip"]
x   = log10.( abs.((as["M"] .- as["ΔM"]) ./ as["ΔM"]) )
miny= minimum(pyplottable(y)[findall(x->x .!== NaN, pyplottable(y))])
maxy= maximum(pyplottable(y)[findall(x->x .!== NaN, pyplottable(y))])
minx= minimum(pyplottable(x)[findall(x->x .!== NaN, pyplottable(x))])
maxx= maximum(pyplottable(x)[findall(x->x .!== NaN, pyplottable(x))])
slcts       = Dict{Int64, Vector{Int64}}()
slcts[1]    = start_thr
slcts[2]    = J_md
slcts[3]    = J_lo
slcts[4]    = J_hi
labels      = Dict{Int64, String}()
labels[1]   = "Full Set N = $(length(slcts[1]))"
labels[2]   = "Medium J, N = $(length(slcts[2]))"
labels[3]   = "Low J, N = $(length(slcts[3]))"
labels[4]   = "High J, N = $(length(slcts[4]))"
plot_hexbins_nxm( x, y, 
            slcts, labels, 2, 2,
            #scale=0.0,
            cbarshift=2.5,
            outfile="/home/moon/sfortune/spinevo/plots/correlations/flip_vs_delMRatio_Jsplit.png", 
            lognorm=true, gridres=(20,7), 
            ymin=miny, ymax=maxy,
            xmin=minx, xmax=maxx, 
            ylabel="Flip [°]",
            xlabel="log( M / $(latexstring("\$\\Delta\$"))M )", fontsize=30,
            binned_median="x", plot_mergers="x",
            )







# The Proof flip by MM depends on b-value, binned by b
y   = as09["ϕ_flip"]
x   = log10.( abs.(as09["M2"] ./ as09["Mpeak_MM"]) )
miny= minimum(pyplottable(y)[findall(x->x .!== NaN, pyplottable(y))])
maxy= maximum(pyplottable(y)[findall(x->x .!== NaN, pyplottable(y))])
minx= -0.1#minimum(pyplottable(x)[findall(x->x .!== NaN, pyplottable(x))])
maxx= 3.6#maximum(pyplottable(x)[findall(x->x .!== NaN, pyplottable(x))])
slcts       = Dict{Int64, Vector{Int64}}()
slcts[1]    = M9_start_thr
slcts[2]    = M9_start_int_nfakes
slcts[3]    = M9_start_ell_nfakes
slcts[4]    = M9_start_disk_nfakes
labels      = Dict{Int64, String}()
labels[1]   = "Full Set N = $(length(slcts[1]))"
labels[2]   = "Intermediates, N = $(length(slcts[2]))"
labels[3]   = "Ellipticals, N = $(length(slcts[3]))"
labels[4]   = "Disks, N = $(length(slcts[4]))"
plot_hexbins_nxm( x, y, 
            slcts, labels, 2, 2,
            #scale=0.0,
            cbarshift=2.5,
            outfile="/home/moon/sfortune/spinevo/plots/correlations/flip_vs_MMRatio_BVALsplit_0.9Gyr.png", 
            lognorm=true, gridres=(30,13), 
            ymin=miny, ymax=maxy,
            xmin=minx, xmax=maxx, 
            ylabel="Flip [°]",
            xlabel="log( Most Massive Merger Ratio )", fontsize=30,
            binned_median="x", plot_mergers="x",
            )






# flip by MM split by b-value and redshift
y   = replace(as09["ϕ_flip"], 0.0 => NaN)
x   = log10.( abs.(as09["M2"] ./ as09["Mpeak_MM"]) )
miny= minimum(pyplottable(y)[findall(x->x .!== NaN, pyplottable(y))])
maxy= maximum(pyplottable(y)[findall(x->x .!== NaN, pyplottable(y))])
minx= minimum(pyplottable(x)[findall(x->x .!== NaN, pyplottable(x))])
maxx= maximum(pyplottable(x)[findall(x->x .!== NaN, pyplottable(x))])
slcts       = Dict{Int64, Vector{Int64}}()
slcts[1]    = M9_start_ell_nfakes[  M9_start_ell_nfakes     .∈  Ref(Set(M9_z_prepeak      ))]
slcts[2]    = M9_start_int_nfakes[  M9_start_int_nfakes     .∈  Ref(Set(M9_z_prepeak      ))]
slcts[3]    = M9_start_disk_nfakes[ M9_start_disk_nfakes    .∈  Ref(Set(M9_z_prepeak      ))]
slcts[4]    = M9_start_ell_nfakes[  M9_start_ell_nfakes     .∈  Ref(Set(M9_z_postpeak      ))]
slcts[5]    = M9_start_int_nfakes[  M9_start_int_nfakes     .∈  Ref(Set(M9_z_postpeak      ))]
slcts[6]    = M9_start_disk_nfakes[ M9_start_disk_nfakes    .∈  Ref(Set(M9_z_postpeak      ))]
labels      = Dict{Int64, String}()
labels[1]   = "ETG, z $(latexstring("\$\\geq\$")) 1.9, N = $(length(slcts[1]))"
labels[2]   = "Int, z $(latexstring("\$\\geq\$")) 1.9, N = $(length(slcts[2]))"
labels[3]   = "LTG, z $(latexstring("\$\\geq\$")) 1.9, N = $(length(slcts[3]))"
labels[4]   = "ETG, z $(latexstring("\$\\leq\$")) 1.9, N = $(length(slcts[4]))"
labels[5]   = "Int, z $(latexstring("\$\\leq\$")) 1.9, N = $(length(slcts[5]))"
labels[6]   = "LTG, z $(latexstring("\$\\leq\$")) 1.9, N = $(length(slcts[6]))"
plot_hexbins_nxm( x, y, 
            slcts, labels, 3, 2,
            #scale=0.0,
            cbarshift=2.5,
            outfile="/home/moon/sfortune/spinevo/plots/correlations/flip_vs_MMRatio_BVAL+z-split1.9_0.9Gyr.png", 
            lognorm=true, gridres=(18,7), 
            ymin=miny, ymax=maxy,
            xmin=minx, xmax=maxx, 
            ylabel="Flip [°]",
            xlabel="log( Most Massive Merger Ratio )", fontsize=30,
            binned_median="x", plot_mergers="x",
            )






# flip by MM split by b-value and redshift
y   = replace(as09["ϕ_flip"], 0.0 => NaN)
x   = log10.( abs.(as09["M2"] ./ as09["Mpeak_MM"]) )
miny= minimum(pyplottable(y)[findall(x->x .!== NaN, pyplottable(y))])
maxy= maximum(pyplottable(y)[findall(x->x .!== NaN, pyplottable(y))])
minx= minimum(pyplottable(x)[findall(x->x .!== NaN, pyplottable(x))])
maxx= maximum(pyplottable(x)[findall(x->x .!== NaN, pyplottable(x))])
slcts       = Dict{Int64, Vector{Int64}}()
slcts[1]    = M9_start_ell_nfakes[  M9_start_ell_nfakes     .∈  Ref(Set(M9_z_pre047      ))]
slcts[2]    = M9_start_int_nfakes[  M9_start_int_nfakes     .∈  Ref(Set(M9_z_pre047      ))]
slcts[3]    = M9_start_disk_nfakes[ M9_start_disk_nfakes    .∈  Ref(Set(M9_z_pre047      ))]
slcts[4]    = M9_start_ell_nfakes[  M9_start_ell_nfakes     .∈  Ref(Set(M9_z_post047      ))]
slcts[5]    = M9_start_int_nfakes[  M9_start_int_nfakes     .∈  Ref(Set(M9_z_post047      ))]
slcts[6]    = M9_start_disk_nfakes[ M9_start_disk_nfakes    .∈  Ref(Set(M9_z_post047      ))]
labels      = Dict{Int64, String}()
labels[1]   = "ETG, z $(latexstring("\$\\geq\$")) 0.47, N = $(length(slcts[1]))"
labels[2]   = "Int, z $(latexstring("\$\\geq\$")) 0.47, N = $(length(slcts[2]))"
labels[3]   = "LTG, z $(latexstring("\$\\geq\$")) 0.47, N = $(length(slcts[3]))"
labels[4]   = "ETG, z $(latexstring("\$\\leq\$")) 0.47, N = $(length(slcts[4]))"
labels[5]   = "Int, z $(latexstring("\$\\leq\$")) 0.47, N = $(length(slcts[5]))"
labels[6]   = "LTG, z $(latexstring("\$\\leq\$")) 0.47, N = $(length(slcts[6]))"
plot_hexbins_nxm( x, y, 
            slcts, labels, 3, 2,
            #scale=0.0,
            cbarshift=2.5,
            outfile="/home/moon/sfortune/spinevo/plots/correlations/flip_vs_MMRatio_BVAL+z-split047_0.9Gyr.png", 
            lognorm=true, gridres=(18,7), 
            ymin=miny, ymax=maxy,
            xmin=minx, xmax=maxx, 
            ylabel="Flip [°]",
            xlabel="log( Most Massive Merger Ratio )", fontsize=30,
            binned_median="x", plot_mergers="x",
            )







# flip by MM split by b-value and redshift
y   = replace(as09["ϕ_flip"], 0.0 => NaN)
x   = log10.( abs.(as09["M2"] ./ as09["Mpeak_MM"]) )
miny= minimum(pyplottable(y)[findall(x->x .!== NaN, pyplottable(y))])
maxy= maximum(pyplottable(y)[findall(x->x .!== NaN, pyplottable(y))])
minx= minimum(pyplottable(x)[findall(x->x .!== NaN, pyplottable(x))])
maxx= maximum(pyplottable(x)[findall(x->x .!== NaN, pyplottable(x))])
slcts       = Dict{Int64, Vector{Int64}}()
slcts[1]    = M9_start_ell_nfakes[  M9_start_ell_nfakes     .∈  Ref(Set(M9_z_pre064      ))]
slcts[2]    = M9_start_int_nfakes[  M9_start_int_nfakes     .∈  Ref(Set(M9_z_pre064      ))]
slcts[3]    = M9_start_disk_nfakes[ M9_start_disk_nfakes    .∈  Ref(Set(M9_z_pre064      ))]
slcts[4]    = M9_start_ell_nfakes[  M9_start_ell_nfakes     .∈  Ref(Set(M9_z_post064      ))]
slcts[5]    = M9_start_int_nfakes[  M9_start_int_nfakes     .∈  Ref(Set(M9_z_post064      ))]
slcts[6]    = M9_start_disk_nfakes[ M9_start_disk_nfakes    .∈  Ref(Set(M9_z_post064      ))]
labels      = Dict{Int64, String}()
labels[1]   = "ETG, z $(latexstring("\$\\geq\$")) 0.64, N = $(length(slcts[1]))"
labels[2]   = "Int, z $(latexstring("\$\\geq\$")) 0.64, N = $(length(slcts[2]))"
labels[3]   = "LTG, z $(latexstring("\$\\geq\$")) 0.64, N = $(length(slcts[3]))"
labels[4]   = "ETG, z $(latexstring("\$\\leq\$")) 0.64, N = $(length(slcts[4]))"
labels[5]   = "Int, z $(latexstring("\$\\leq\$")) 0.64, N = $(length(slcts[5]))"
labels[6]   = "LTG, z $(latexstring("\$\\leq\$")) 0.64, N = $(length(slcts[6]))"
plot_hexbins_nxm( x, y, 
            slcts, labels, 3, 2,
            #scale=0.0,
            cbarshift=2.5,
            outfile="/home/moon/sfortune/spinevo/plots/correlations/flip_vs_MMRatio_BVAL+z-split064_0.9Gyr.png", 
            lognorm=true, gridres=(18,7), 
            ymin=miny, ymax=maxy,
            xmin=minx, xmax=maxx, 
            ylabel="Flip [°]",
            xlabel="log( Most Massive Merger Ratio )", fontsize=30,
            binned_median="x", plot_mergers="x",
            )







# flip by MM split by b-value and redshift
y   = replace(as["ϕ_flip"], 0.0 => NaN)
x   = log10.( abs.(as["M2"] ./ as["Mpeak_MM"]) )
miny= minimum(pyplottable(y)[findall(x->x .!== NaN, pyplottable(y))])
maxy= maximum(pyplottable(y)[findall(x->x .!== NaN, pyplottable(y))])
minx= minimum(pyplottable(x)[findall(x->x .!== NaN, pyplottable(x))])
maxx= maximum(pyplottable(x)[findall(x->x .!== NaN, pyplottable(x))])
slcts       = Dict{Int64, Vector{Int64}}()
slcts[1]    = start_ell_nfakes[  start_ell_nfakes     .∈  Ref(Set(z_prepeak      ))]
slcts[2]    = start_int_nfakes[  start_int_nfakes     .∈  Ref(Set(z_prepeak      ))]
slcts[3]    = start_disk_nfakes[ start_disk_nfakes    .∈  Ref(Set(z_prepeak      ))]
slcts[4]    = start_ell_nfakes[  start_ell_nfakes     .∈  Ref(Set(z_postpeak      ))]
slcts[5]    = start_int_nfakes[  start_int_nfakes     .∈  Ref(Set(z_postpeak      ))]
slcts[6]    = start_disk_nfakes[ start_disk_nfakes    .∈  Ref(Set(z_postpeak      ))]
labels      = Dict{Int64, String}()
labels[1]   = "ETG, z $(latexstring("\$\\geq\$")) 1.9, N = $(length(slcts[1]))"
labels[2]   = "Int, z $(latexstring("\$\\geq\$")) 1.9, N = $(length(slcts[2]))"
labels[3]   = "LTG, z $(latexstring("\$\\geq\$")) 1.9, N = $(length(slcts[3]))"
labels[4]   = "ETG, z $(latexstring("\$\\leq\$")) 1.9, N = $(length(slcts[4]))"
labels[5]   = "Int, z $(latexstring("\$\\leq\$")) 1.9, N = $(length(slcts[5]))"
labels[6]   = "LTG, z $(latexstring("\$\\leq\$")) 1.9, N = $(length(slcts[6]))"
plot_hexbins_nxm( x, y, 
            slcts, labels, 3, 2,
            #scale=0.0,
            cbarshift=2.5,
            outfile="/home/moon/sfortune/spinevo/plots/correlations/flip_vs_MMRatio_BVAL+z-split1.9_0.3Gyr.png", 
            lognorm=true, gridres=(18,7), 
            ymin=miny, ymax=maxy,
            xmin=minx, xmax=maxx, 
            ylabel="Flip [°]",
            xlabel="log( Most Massive Merger Ratio )", fontsize=30,
            binned_median="x", plot_mergers="x",
            )






# flip by MM split by b-value and redshift
y   = replace(as["ϕ_flip"], 0.0 => NaN)
x   = log10.( abs.(as["M2"] ./ as["Mpeak_MM"]) )
miny= minimum(pyplottable(y)[findall(x->x .!== NaN, pyplottable(y))])
maxy= maximum(pyplottable(y)[findall(x->x .!== NaN, pyplottable(y))])
minx= minimum(pyplottable(x)[findall(x->x .!== NaN, pyplottable(x))])
maxx= maximum(pyplottable(x)[findall(x->x .!== NaN, pyplottable(x))])
slcts       = Dict{Int64, Vector{Int64}}()
slcts[1]    = start_ell_nfakes[  start_ell_nfakes     .∈  Ref(Set(z_pre047      ))]
slcts[2]    = start_int_nfakes[  start_int_nfakes     .∈  Ref(Set(z_pre047      ))]
slcts[3]    = start_disk_nfakes[ start_disk_nfakes    .∈  Ref(Set(z_pre047      ))]
slcts[4]    = start_ell_nfakes[  start_ell_nfakes     .∈  Ref(Set(z_post047      ))]
slcts[5]    = start_int_nfakes[  start_int_nfakes     .∈  Ref(Set(z_post047      ))]
slcts[6]    = start_disk_nfakes[ start_disk_nfakes    .∈  Ref(Set(z_post047      ))]
labels      = Dict{Int64, String}()
labels[1]   = "ETG, z $(latexstring("\$\\geq\$")) 0.47, N = $(length(slcts[1]))"
labels[2]   = "Int, z $(latexstring("\$\\geq\$")) 0.47, N = $(length(slcts[2]))"
labels[3]   = "LTG, z $(latexstring("\$\\geq\$")) 0.47, N = $(length(slcts[3]))"
labels[4]   = "ETG, z $(latexstring("\$\\leq\$")) 0.47, N = $(length(slcts[4]))"
labels[5]   = "Int, z $(latexstring("\$\\leq\$")) 0.47, N = $(length(slcts[5]))"
labels[6]   = "LTG, z $(latexstring("\$\\leq\$")) 0.47, N = $(length(slcts[6]))"
plot_hexbins_nxm( x, y, 
            slcts, labels, 3, 2,
            #scale=0.0,
            cbarshift=2.5,
            outfile="/home/moon/sfortune/spinevo/plots/correlations/flip_vs_MMRatio_BVAL+z-split047_0.3Gyr.png", 
            lognorm=true, gridres=(18,7), 
            ymin=miny, ymax=maxy,
            xmin=minx, xmax=maxx, 
            ylabel="Flip [°]",
            xlabel="log( Most Massive Merger Ratio )", fontsize=30,
            binned_median="x", plot_mergers="x",
            )







# flip by MM split by b-value and redshift
y   = replace(as["ϕ_flip"], 0.0 => NaN)
x   = log10.( abs.(as["M2"] ./ as["Mpeak_MM"]) )
miny= minimum(pyplottable(y)[findall(x->x .!== NaN, pyplottable(y))])
maxy= maximum(pyplottable(y)[findall(x->x .!== NaN, pyplottable(y))])
minx= minimum(pyplottable(x)[findall(x->x .!== NaN, pyplottable(x))])
maxx= maximum(pyplottable(x)[findall(x->x .!== NaN, pyplottable(x))])
slcts       = Dict{Int64, Vector{Int64}}()
slcts[1]    = start_ell_nfakes[  start_ell_nfakes     .∈  Ref(Set(z_pre064      ))]
slcts[2]    = start_int_nfakes[  start_int_nfakes     .∈  Ref(Set(z_pre064      ))]
slcts[3]    = start_disk_nfakes[ start_disk_nfakes    .∈  Ref(Set(z_pre064      ))]
slcts[4]    = start_ell_nfakes[  start_ell_nfakes     .∈  Ref(Set(z_post064      ))]
slcts[5]    = start_int_nfakes[  start_int_nfakes     .∈  Ref(Set(z_post064      ))]
slcts[6]    = start_disk_nfakes[ start_disk_nfakes    .∈  Ref(Set(z_post064      ))]
labels      = Dict{Int64, String}()
labels[1]   = "ETG, z $(latexstring("\$\\geq\$")) 0.64, N = $(length(slcts[1]))"
labels[2]   = "Int, z $(latexstring("\$\\geq\$")) 0.64, N = $(length(slcts[2]))"
labels[3]   = "LTG, z $(latexstring("\$\\geq\$")) 0.64, N = $(length(slcts[3]))"
labels[4]   = "ETG, z $(latexstring("\$\\leq\$")) 0.64, N = $(length(slcts[4]))"
labels[5]   = "Int, z $(latexstring("\$\\leq\$")) 0.64, N = $(length(slcts[5]))"
labels[6]   = "LTG, z $(latexstring("\$\\leq\$")) 0.64, N = $(length(slcts[6]))"
plot_hexbins_nxm( x, y, 
            slcts, labels, 3, 2,
            #scale=0.0,
            cbarshift=2.5,
            outfile="/home/moon/sfortune/spinevo/plots/correlations/flip_vs_MMRatio_BVAL+z-split064_0.3Gyr.png", 
            lognorm=true, gridres=(18,7), 
            ymin=miny, ymax=maxy,
            xmin=minx, xmax=maxx, 
            ylabel="Flip [°]",
            xlabel="log( Most Massive Merger Ratio )", fontsize=30,
            binned_median="x", plot_mergers="x",
            )















# Flips vs sSFR at different redshift bins
x   = as["ϕ_flip"]
y   = log10.(as["SFR"] ./ as["M"] )
miny= minimum(pyplottable(y)[findall(x->x .!== NaN, pyplottable(y))])
maxy= maximum(pyplottable(y)[findall(x->x .!== NaN, pyplottable(y))])
minx= minimum(pyplottable(x)[findall(x->x .!== NaN, pyplottable(x))])
maxx= maximum(pyplottable(x)[findall(x->x .!== NaN, pyplottable(x))])
slcts       = Dict{Int64, Vector{Int64}}()
slcts[1]    = start_thr
slcts[2]    = z_md
slcts[3]    = z_hi 
slcts[4]    = z_lo
labels      = Dict{Int64, String}()
labels[1]   = "Full Set, N = $(length(slcts[1]))"
labels[2]   = "$(round(maximum(as["redshift"][slcts[2]]), digits=2)) $(latexstring("\$\\geq\$")) z $(latexstring("\$\\geq\$")) $(round(minimum(as["redshift"][slcts[2]]), digits=2)), N = $(length(slcts[2]))"
labels[3]   = "$(round(maximum(as["redshift"][slcts[3]]), digits=2)) $(latexstring("\$\\geq\$")) z $(latexstring("\$\\geq\$")) $(round(minimum(as["redshift"][slcts[3]]), digits=2)), N = $(length(slcts[3]))"
labels[4]   = "$(round(maximum(as["redshift"][slcts[4]]), digits=2)) $(latexstring("\$\\geq\$")) z $(latexstring("\$\\geq\$")) $(round(minimum(as["redshift"][slcts[4]]), digits=2)), N = $(length(slcts[4]))"
plot_hexbins_nxm( x, 
            y, 
            slcts, labels, 2, 2,
            #scale=0.0,
            cbarshift=2.5,
            outfile="/home/moon/sfortune/spinevo/plots/correlations/sSFR_vs_flips_zsplit.png", 
            lognorm=true, gridres=(30,10), 
            ymin=miny, ymax=maxy,
            xmin=minx, xmax=maxx, 
            ylabel="log( sSFR [$(latexstring("\\frac{1}{yr}"))] )",
            xlabel="Flip [°]", fontsize=30,
            binned_median="x", #plot_bval="x",
            )





# delta J vs sSFR at different redshift bins
x   = log10.(pyplottable(as["ΔJ_main"]))
y   = log10.(as["SFR"] ./ as["M"] )
miny= minimum(pyplottable(y)[findall(x->x .!== NaN, pyplottable(y))])
maxy= maximum(pyplottable(y)[findall(x->x .!== NaN, pyplottable(y))])
minx= minimum(pyplottable(x)[findall(x->x .!== NaN, pyplottable(x))])
maxx= maximum(pyplottable(x)[findall(x->x .!== NaN, pyplottable(x))])
slcts       = Dict{Int64, Vector{Int64}}()
slcts[1]    = z_lo
slcts[2]    = z_lo[ z_lo .∈  Ref(Set(start_int_nfakes      ))]
slcts[3]    = z_lo[ z_lo .∈  Ref(Set(start_disk_nfakes      ))] 
slcts[4]    = z_lo[ z_lo .∈  Ref(Set(start_ell_nfakes      ))]
labels      = Dict{Int64, String}()
labels[1]   = "Full Set, N = $(length(slcts[1]))"
labels[2]   = "Intermediates, N = $(length(slcts[2]))"
labels[3]   = "Disks, N = $(length(slcts[3]))"
labels[4]   = "Ellipticals, N = $(length(slcts[4]))"
plot_hexbins_nxm( x, 
            y, 
            slcts, labels, 2, 2,
            #scale=0.0,
            cbarshift=2.5,
            outfile="/home/moon/sfortune/spinevo/plots/correlations/sSFR_vs_delJ_BVALsplit_zlow.png", 
            lognorm=true, gridres=(30,10), 
            ymin=miny, ymax=maxy,
            xmin=minx, xmax=maxx, 
            ylabel="log( sSFR [$(latexstring("\\frac{1}{yr}"))] )",
            xlabel="log( $(latexstring("\$\\Delta\$"))J_STARS )", fontsize=30,
            binned_median="x", #plot_bval="x",
            )







# Hirschmann 2014

pltcm       = pyimport("matplotlib.cm")
pltcolors   = pyimport("matplotlib.colors")
colormap    = pltcm.get_cmap(name="gnuplot")
scale=0.5
N_bins=10
fig, ax = subplots()
n_0, x_0      = PyPlot.hist( log10.(pyplottable(as["M"][z_hirsch_0])), bins=N_bins, density=true, log=false)
n_05,x_05     = PyPlot.hist( log10.(pyplottable(as["M"][z_hirsch_05])), bins=N_bins, density=true, log=false)
n_1, x_1      = PyPlot.hist( log10.(pyplottable(as["M"][z_hirsch_1])), bins=N_bins, density=true, log=false)
n_2, x_2      = PyPlot.hist( log10.(pyplottable(as["M"][z_hirsch_2])), bins=N_bins, density=true, log=false)
n_3, x_3      = PyPlot.hist( log10.(pyplottable(as["M"][z_hirsch_3])), bins=N_bins, density=true, log=false)
n_4, x_4      = PyPlot.hist( log10.(pyplottable(as["M"][z_hirsch_4])), bins=N_bins, density=true, log=false)
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
ax.plot(bcenters_0 , log10.(n_0 .* (length(z_hirsch_0)/48^3) ./ length(unique(as["redshift"][z_hirsch_0]))) , ".-", ms=10, label="z $(latexstring("\$\\in\$"))[$(round(minimum(as["redshift"][z_hirsch_0]),digits=2)), $(round(maximum(as["redshift"][z_hirsch_0]),digits=2))], N = $(length(z_hirsch_0))         ", color=colormap(0/4)   )
ax.plot(bcenters_05 ,log10.(n_05.*(length(z_hirsch_05)/48^3) ./ length(unique(as["redshift"][z_hirsch_05]))) ,".-", ms=10, label="z $(latexstring("\$\\in\$"))[$(round(minimum(as["redshift"][z_hirsch_05]),digits=2)), $(round(maximum(as["redshift"][z_hirsch_05]),digits=2))], N = $(length(z_hirsch_05))   ", color=colormap(0.5/4)   )
ax.plot(bcenters_1 , log10.(n_1 .* (length(z_hirsch_1)/48^3) ./ length(unique(as["redshift"][z_hirsch_1]))) , ".-", ms=10, label="z $(latexstring("\$\\in\$"))[$(round(minimum(as["redshift"][z_hirsch_1]),digits=2)), $(round(maximum(as["redshift"][z_hirsch_1]),digits=2))], N = $(length(z_hirsch_1))         ", color=colormap(1/4)   )
ax.plot(bcenters_2 , log10.(n_2 .* (length(z_hirsch_2)/48^3) ./ length(unique(as["redshift"][z_hirsch_2]))) , ".-", ms=10, label="z $(latexstring("\$\\in\$"))[$(round(minimum(as["redshift"][z_hirsch_2]),digits=2)), $(round(maximum(as["redshift"][z_hirsch_2]),digits=2))], N = $(length(z_hirsch_2))         ", color=colormap(2/4)   )
ax.plot(bcenters_3 , log10.(n_3 .* (length(z_hirsch_3)/48^3) ./ length(unique(as["redshift"][z_hirsch_3]))) , ".-", ms=10, label="z $(latexstring("\$\\in\$"))[$(round(minimum(as["redshift"][z_hirsch_3]),digits=2)), $(round(maximum(as["redshift"][z_hirsch_3]),digits=2))], N = $(length(z_hirsch_3))         ", color=colormap(3/4)   )
ax.plot(bcenters_4 , log10.(n_4 .* (length(z_hirsch_4)/48^3) ./ length(unique(as["redshift"][z_hirsch_4]))) , ".-", ms=10, label="z $(latexstring("\$\\in\$"))[$(round(minimum(as["redshift"][z_hirsch_4]),digits=2)), $(round(maximum(as["redshift"][z_hirsch_4]),digits=2))], N = $(length(z_hirsch_4))         ", color=colormap(4/4)   )
ax.set_xlabel("log( M_STARS [ $(latexstring("\$M_\\odot\$"))] )")
ax.set_ylabel("log( $(latexstring("\$\\Phi\$")) [ $(latexstring("\$\\frac{h^3}{(Mpc^3)}\$")) ] )")
ax.set_yscale("linear")
ax.set_ylim(nothing,nothing)
#ax.set_title("$(latexstring("\$\\Delta\$"))t $(latexstring("\$\\sim\$")) 0.3 Gyr within $(round(maximum(as["redshift"]),digits=2)) $(latexstring("\$\\geq\$")) z $(latexstring("\$\\geq\$")) $(round(minimum(as["redshift"]),digits=2))")
ax.legend(loc="upper right", frameon=true, borderpad=1, handlelength=1.8)
ax.grid(false)
fig.set_size_inches(16scale, 9scale)
fig.tight_layout()
fig.savefig("/home/moon/sfortune/spinevo/plots/correlations/stellarmassfct_Hirschmann2014like.png", bbox_inches="tight", pad_inches=.1)
close("all")
PyPlot.matplotlib[:rc]("text", usetex=false)














# Robustness breaking
# flip by mass vs del j
# flip by mass vs del J
# flip by del j











# Flips impact 
bval_evo      = Dict{Int64, Vector{Int64}(undef, 2)}()
for i in as["ID_ISUB"][snap136]
    max = findcs( as["BVAL_0"] )
end