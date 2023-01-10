
include("/home/moon/sfortune/spinevo/pkg/meta.jl")
assembly_STARS  = load("/home/moon/sfortune/spinevo/mergerimpact/assembly_Mstar_1.0e10.jld", "assembly_STARS")
centrals_136 = filter(assembly_STARS, condition="snap")


fig, ax = subplots()

ax.hist(assembly_STARS["M2_felix"][centrals_136["main"]], bins=100, label="M2_felix  at z ≈ $(@sprintf("%.2f", assembly_STARS["redshift"][centrals_136["main"]][1]))", rwidth=0.9, color="navy")


ax.set_xlabel("M_⊙")
ax.set_xscale("log")
ax.set_yscale("log")
ax.legend(loc="upper right", frameon=true, borderpad=1, handlelength=1.8)
ax.grid()

scale=0.5

fig.set_size_inches(16scale, 9scale)
fig.tight_layout()

println(minimum(assembly_STARS["M2_felix"][centrals_136["main"]]))
println(maximum(assembly_STARS["M2_felix"][centrals_136["main"]]))
println(Statistics.median(assembly_STARS["M2_felix"][centrals_136["main"]]))
#fig.savefig(joinpath(@__DIR__, "masses_z$(@sprintf("%.2f", head.z)).png"), bbox_inches="tight", pad_inches=.1)

m_smalls    = 10^10.5
m_intermediates = 10^11.5
smalls          = missings(Int64,0)
intermediates   = missings(Int64,0)
bigs            = missings(Int64,0)

for i in centrals_136["main"]
    if count(ismissing, assembly_STARS["J_SUMorbital"][:,i]) > 0 || count(ismissing, assembly_STARS["J_main"][:,i]) > 0
    elseif assembly_STARS["M2_felix"][i] < m_smalls
        smalls  = vcat(smalls, i)
    elseif assembly_STARS["M2_felix"][i] < m_intermediates
        intermediates  = vcat(intermediates, i)
    else
        bigs    = vcat(bigs,i)
    end
end



asm = load("/home/moon/sfortune/spinevo/newtrees_test/assembly_Mstar_1.0e10.jld", "assembly_STARS")

#search_index(; table=load("/home/moon/sfortune/spinevo/newtrees_test/assembly_Mstar_1.0e10.jld", "assembly_STARS")
#    , condition=
#    )

findall(x->x.==maximum(norm.(eachcol(replace(asm["j_main"], missing => 0.)))), norm.(eachcol(replace(asm["j_main"], missing => 0.))))

println("$(asm["ID_ISUB"][6230])   ---   $(asm["snapNR"][6230])")

find_felixID(3212)

plot_hexbins(ams["ϕ_flip"], ams["LOOKBACKTIME"], outfile="./time_vs_flipangle.png", 
    lognorm=true, gridres=10, xlabel="Flip [°]", ylabel="Lookback Time [Gyr]")


plot_hexbins(ams["BVAL"], log10.(pyplottable(ams["δj_main"])./pyplottable(ams["j_main"])), outfile="./djrel_vs_flipangle.png", 
    lognorm=true, gridres=20, xlabel="Flip [°]", ylabel="log10( dj / j )")


# Relative Mass change
plot_hexbins(ams["ϕ_flip"], log10.(abs.(ams["δM_felix"])./ams["M_felix"]), outfile="./stellar_relmassdiff_vs_flipangle.png", 
    lognorm=true, gridres=20, xlabel="Flip [°]", ylabel="log10( dM / M )")

plot_hexbins(ams["ϕ_flip"], log10.(ams["M2_felix"]./ams["M2_MM"]), outfile="./stellar_mergerratio_vs_flipangle.png", 
    lognorm=true, gridres=20, xlabel="Flip [°]", ylabel="log10( Merger Ratio )", title="Stellar")

plot_hexbins(ams["ϕ_flip"], log10.(amd["M2_felix"]./amd["M2_MM"]), outfile="./dm_mergerratio_vs_flipangle.png", 
    lognorm=true, gridres=20, xlabel="Flip [°]", ylabel="log10( Merger Ratio )", title="Dark Matter")

min_tdiff = 100.
lbt = 0.
max_tdiff = 0.
for i in 1:length(hs["redshift"])
    head    = read_header("$current_dir_simbox/groups_$(@sprintf("%03i", hs["snapNR"][i]))/sub_$(@sprintf("%03i", hs["snapNR"][i]))")
    lbt_old = lbt
    lbt     = ustrip(lookback_time(cosmology(h=head.h0, OmegaM=head.omega_0), head.z))
    println("SNAP = $(hs["snapNR"][i])   ---   TIMEDIFF = $(lbt-lbt_old)")
    max_tdiff = ifelse(lbt-lbt_old > max_tdiff && hs["snapNR"][i] < 136, lbt-lbt_old, max_tdiff)
    min_tdiff = ifelse(0. < lbt-lbt_old < min_tdiff && hs["snapNR"][i] < 136, lbt-lbt_old, min_tdiff)
    #if lbt-lbt_old > max_tdiff
    #    max_tdiff = lbt-lbt_old
end
println(min_tdiff)
println(max_tdiff)

###################################################################################
###################################################################################
###################################################################################


for i in 1:100
    println("$(old_halo3212["snapNR"][i])      $(old_halo3212["M_STARS"][i])          $(old_halo3212["j_STARS"][:,i])    ---   $(halo3212["snapNR"][i])      $(halo3212["M_STARS"][i])          $(halo3212["j_STARS"][:,i])")
end


mcs3212 = load(joinpath(current_dir_stories, find_felixID(3212)[2]), "merger_collection_STARS")
testas = load("/home/moon/sfortune/spinevo/data/assembly_centrals2_20220112.jld", "assembly_STARS")
include("/home/moon/sfortune/spinevo/pkg/meta.jl")
testas2 = load("/home/moon/sfortune/spinevo/data/assembly_centrals2co_20220112.jld", "assembly_STARS")
testas2 = load("/home/moon/sfortune/spinevo/data/assembly_20220112_centrals2co_mfromj.jld", "assembly_STARS")

for i in 1:length(testas2["j_main"][1,:])
    if testas2["ID_ISUB"][i] == 3212.0
        println("$(fake_flip_finder(testas2, i, verbose = true))")
    end
end

jumps       = findcs(testas["jump"])
n_jumps     = findcs(testas["jump"], val=0)
exceeds     = findcs(testas["exceed"])
n_exceeds   = findcs(testas["exceed"], val=0)
switches    = findcs(testas["switch"])
n_switches  = findcs(testas["switch"], val=0)
borders     = findcs(testas["border"])
n_borders   = findcs(testas["border"], val=0)
fakes       = findcs(testas["FAKEFLIP"])
n_fakes     = findcs(testas["FAKEFLIP"], val=0)

jumps2       = findcs(testas2["jump"])
n_jumps2     = findcs(testas2["jump"], eq=0)
exceeds2     = findcs(testas2["exceed"])
n_exceeds2   = findcs(testas2["exceed"], eq=0)
switches2    = findcs(testas2["switch"])
n_switches2  = findcs(testas2["switch"], eq=0)
borders2     = findcs(testas2["border"])
n_borders2   = findcs(testas2["border"], eq=0)
fakes2       = findcs(testas2["FAKEFLIP"])
n_fakes2     = findcs(testas2["FAKEFLIP"], eq=0)


#all_notin   = ALL_STARS["ID"] .∉ Ref(Set(particleID_list))
nf_nj   = n_fakes2 .∈ Ref(Set(n_jumps2))
nf_ne   = n_fakes2 .∈ Ref(Set(n_exceeds2))

plot_hexbins(testas2["ϕ_flip"][:], 
            log10.(testas2["M2_felix"][:]./testas2["M2_MM"][:]), 
            outfile="./monday/flip_vs_mmratio_co_all.png", 
            lognorm=true, gridres=20, 
            xlabel="Flip [°]", 
            ylabel="log10( Merger Ratio )", 
            title="All Flips ($(length(testas2["ϕ_flip"])))")

plot_hexbins(testas2["ϕ_flip"][n_fakes2], 
            log10.(testas2["M2_felix"][n_fakes2]./testas2["M2_MM"][n_fakes2]), 
            outfile="./monday/flip_vs_mmratio_co_nfakes.png", 
            lognorm=true, gridres=20, 
            xlabel="Flip [°]", 
            ylabel="log10( Merger Ratio )", 
            title="Real Flips ($(length(n_fakes2)) / $(length(testas2["ϕ_flip"])))")

plot_hexbins(testas2["ϕ_flip"][fakes2], 
            log10.(testas2["M2_felix"][fakes2]./testas2["M2_MM"][fakes2]), 
            outfile="./monday/flip_vs_mmratio_co_fakes.png", 
            lognorm=true, gridres=20, 
            xlabel="Flip [°]", 
            ylabel="log10( Merger Ratio )", 
            title="Fake Flips ($(length(fakes2)) / $(length(testas2["ϕ_flip"])))")

plot_hexbins(testas2["ϕ_flip"][borders2], 
            log10.(testas2["M2_felix"][borders2]./testas2["M2_MM"][borders2]), 
            outfile="./monday/flip_vs_mmratio_co_borders.png", 
            lognorm=true, gridres=20, 
            xlabel="Flip [°]", 
            ylabel="log10( Merger Ratio )", 
            title="Border Crosser ($(length(borders2)) / $(length(testas2["ϕ_flip"])))")

plot_hexbins(testas2["ϕ_flip"][jumps2], 
            log10.(testas2["M2_felix"][jumps2]./testas2["M2_MM"][jumps2]), 
            outfile="./monday/flip_vs_mmratio_co_jump.png", 
            lognorm=true, gridres=20, 
            xlabel="Flip [°]", 
            ylabel="log10( Merger Ratio )", 
            title="Δx > 200 kpc ($(length(jumps2)) / $(length(testas2["ϕ_flip"])))")

plot_hexbins(testas2["ϕ_flip"][exceeds2], 
            log10.(testas2["M2_felix"][exceeds2]./testas2["M2_MM"][exceeds2]), 
            outfile="./monday/flip_vs_mmratio_co_exceed.png", 
            lognorm=true, gridres=20, 
            xlabel="Flip [°]", 
            ylabel="log10( Merger Ratio )", 
            title="Δx > t*v ($(length(exceeds2)) / $(length(testas2["ϕ_flip"])))")

plot_hexbins(testas2["ϕ_flip"][n_borders2], 
            log10.(testas2["M2_felix"][n_borders2]./testas2["M2_MM"][n_borders2]), 
            outfile="./monday/flip_vs_mmratio_co_nborders.png", 
            lognorm=true, gridres=20, 
            xlabel="Flip [°]", 
            ylabel="log10( Merger Ratio )", 
            title="No Border Crosser ($(length(n_borders2)) / $(length(testas2["ϕ_flip"])))")

plot_hexbins(testas2["ϕ_flip"][n_jumps2], 
            log10.(testas2["M2_felix"][n_jumps2]./testas2["M2_MM"][n_jumps2]), 
            outfile="./monday/flip_vs_mmratio_co_njump.png", 
            lognorm=true, gridres=20, 
            xlabel="Flip [°]", 
            ylabel="log10( Merger Ratio )", 
            title="Δx < 200 kpc ($(length(n_jumps2)) / $(length(testas2["ϕ_flip"])))")

plot_hexbins(testas2["ϕ_flip"][n_exceeds2], 
            log10.(testas2["M2_felix"][n_exceeds2]./testas2["M2_MM"][n_exceeds2]), 
            outfile="./monday/flip_vs_mmratio_co_nexceed.png", 
            lognorm=true, gridres=20, 
            xlabel="Flip [°]", 
            ylabel="log10( Merger Ratio )", 
            title="Δx < t*v ($(length(n_exceeds2)) / $(length(testas2["ϕ_flip"])))")

plot_hexbins(testas2["ϕ_flip"][n_fakes2[nf_nj]], 
            log10.(testas2["M2_felix"][n_fakes2[nf_nj]]./testas2["M2_MM"][n_fakes2[nf_nj]]), 
            outfile="./monday/flip_vs_mmratio_co_nfakes_njumps.png", 
            lognorm=true, gridres=20, 
            xlabel="Flip [°]", 
            ylabel="log10( Merger Ratio )", 
            title="Real Flips ($(length(n_fakes2[nf_nj])) / $(length(testas2["ϕ_flip"])))")

plot_hexbins(testas2["ϕ_flip"][n_fakes2[nf_ne]], 
            log10.(testas2["M2_felix"][n_fakes2[nf_ne]]./testas2["M2_MM"][n_fakes2[nf_ne]]), 
            outfile="./monday/flip_vs_mmratio_co_nfakes_nexceed.png", 
            lognorm=true, gridres=20, 
            xlabel="Flip [°]", 
            ylabel="log10( Merger Ratio )", 
            title="Real Flips ($(length(n_fakes2[nf_ne])) / $(length(testas2["ϕ_flip"])))")

plot_hexbins(testas2["ϕ_flip"][n_fakes2], 
            log10.(testas2["M_felix"][n_fakes2]./abs.(testas2["δM_felix"][n_fakes2])), 
            outfile="./monday/flip_vs_delMratio_co_nfakes.png", 
            lognorm=true, gridres=20, 
            xlabel="Flip [°]", 
            ylabel="log10( M / ΔM )", 
            title="Real Flips ($(length(n_fakes2)) / $(length(testas2["ϕ_flip"])))")

plot_hexbins(testas2["ϕ_flip"][:], 
            log10.(testas2["M_felix"][:]./abs.(testas2["δM_felix"][:])), 
            outfile="./monday/flip_vs_delMratio_co_all.png", 
            lognorm=true, gridres=20, 
            xlabel="Flip [°]", 
            ylabel="log10( M / ΔM )", 
            title="All Flips ($(length(testas2["ϕ_flip"])))")

plot_hexbins(testas2["ϕ_flip"][fakes2], 
            log10.(testas2["M_felix"][fakes2]./abs.(testas2["δM_felix"][fakes2])), 
            outfile="./monday/flip_vs_delMratio_co_fakes.png", 
            lognorm=true, gridres=20, 
            xlabel="Flip [°]", 
            ylabel="log10( M / ΔM )", 
            title="Fake Flips ($(length(fakes2)) / $(length(testas2["ϕ_flip"])))")


for i in 1:length(testas2["j_main"][1,n_fakes2])
    if testas2["ϕ_flip"][i] > 100 && (log10.(replace(testas2["M2_felix"]./testas2["M2_MM"], missing => 1)))[i] > 3
        println("$(testas2["ID_ISUB"][i])   ---   $(testas2["snapNR"][i])   ---   $(testas2["ϕ_flip"][i])   ---   $((log10.(replace(testas2["M2_felix"]./testas2["M2_MM"], missing => 1)))[i])")
    end
end

plot_group(80; rootID=2323, property="pos", stepsize=50, arsize=10, ptsize=5000, rad=200, res=(1600,900), spin=true)


# sfc: distinguish between b-values



result_disks    = findcs(testas2["BVAL"] .+ 0.5 .*log10.(1 .+testas2["redshift"]), geq=b_disk)
disk_fakes2     = findcs(testas2["FAKEFLIP"], comparewith = result_disks)
disk_nfakes2    = findcs(testas2["FAKEFLIP"], eq=0, comparewith = result_disks)
result_int      = findcs(testas2["BVAL"] .+ 0.5 .*log10.(1 .+testas2["redshift"]), gt=b_ell, lt=b_disk)
int_fakes2     = findcs(testas2["FAKEFLIP"], comparewith = result_int)
int_nfakes2    = findcs(testas2["FAKEFLIP"], eq=0, comparewith = result_int)
result_ell      = findcs(testas2["BVAL"] .+ 0.5 .*log10.(1 .+testas2["redshift"]), leq=b_ell)
ell_fakes2     = findcs(testas2["FAKEFLIP"], comparewith = result_ell)
ell_nfakes2    = findcs(testas2["FAKEFLIP"], eq=0, comparewith = result_ell)

start_disk      = findcs(testas2["BVAL"] .+ 0.5 .*log10.(1 .+testas2["redshift"]) .- testas2["δBVAL_0"], geq=b_disk)
start_int       = findcs(testas2["BVAL"] .+ 0.5 .*log10.(1 .+testas2["redshift"]) .- testas2["δBVAL_0"], gt=b_ell, lt=b_disk)
start_ell       = findcs(testas2["BVAL"] .+ 0.5 .*log10.(1 .+testas2["redshift"]) .- testas2["δBVAL_0"], leq=b_ell)
sdisk_nfakes2   = findcs(testas2["FAKEFLIP"], eq=0, comparewith = start_disk)
sint_nfakes2    = findcs(testas2["FAKEFLIP"], eq=0, comparewith = start_int)
sell_nfakes2    = findcs(testas2["FAKEFLIP"], eq=0, comparewith = start_ell)
stay_disk_nfakes= sdisk_nfakes2[sdisk_nfakes2 .∈ Ref(Set(result_disks))]
stay_int_nfakes = sint_nfakes2[sint_nfakes2 .∈ Ref(Set(result_int))]
stay_ell_nfakes = sell_nfakes2[sell_nfakes2 .∈ Ref(Set(result_ell))]
disk_to_ell     = sdisk_nfakes2[sdisk_nfakes2 .∈ Ref(Set(result_ell))]
ell_to_disk     = sell_nfakes2[sell_nfakes2 .∈ Ref(Set(result_disks))]
int_to_disk     = sint_nfakes2[sint_nfakes2 .∈ Ref(Set(result_disks))]
int_to_ell      = sint_nfakes2[sint_nfakes2 .∈ Ref(Set(result_ell))]
sdisk_fakes2    = findcs(testas2["FAKEFLIP"], eq=1, comparewith = start_disk)
sint_fakes2     = findcs(testas2["FAKEFLIP"], eq=1, comparewith = start_int)
sell_fakes2     = findcs(testas2["FAKEFLIP"], eq=1, comparewith = start_ell)
stay_disk_fakes = sdisk_fakes2[sdisk_fakes2 .∈ Ref(Set(result_disks))]
stay_int_fakes  = sint_fakes2[sint_fakes2 .∈ Ref(Set(result_int))]
stay_ell_fakes  = sell_fakes2[sell_fakes2 .∈ Ref(Set(result_ell))]

selct = disk_to_ell
plot_hexbins(log10.(testas2["M2_felix"]./testas2["M2_MM"]), 
            testas2["ϕ_flip"], 
            outfile="./monday/next/flips_LTGtoETG.png", 
            selection=selct,
            lognorm=true, gridres=(10,5), 
            xmin=0, xmax=3.5, ymin=0, ymax=180,
            ylabel="Flip [°]",
            xlabel="log( Merger Ratio )", 
            title="LTG-to-ETG Flips ($(length(testas2["ϕ_flip"][selct])) / $(length(testas2["ϕ_flip"])))",
            binned_median="x",
            )

flip_gt20    = findcs(testas2["ϕ_flip"], geq=20)
flip_gt45    = findcs(testas2["ϕ_flip"], geq=45)

disk_nfakes2_f20    = disk_nfakes2[disk_nfakes2 .∈ Ref(Set(flip_gt20))]
disk_nfakes2_f45    = disk_nfakes2[disk_nfakes2 .∈ Ref(Set(flip_gt45))]
ell_nfakes2_f20    = ell_nfakes2[ell_nfakes2 .∈ Ref(Set(flip_gt20))]
ell_nfakes2_f45    = ell_nfakes2[ell_nfakes2 .∈ Ref(Set(flip_gt45))]



plot_hexbins_nxm(log10.(testas2["M2_felix"] ./ testas2["M2_MM"]), 
            testas2["ϕ_flip"], 
            
            title_col1= "LTGs", title_col2= "ETGs",
            outfile="./monday/next/flips_ltg_etg.png", 
            selection11=disk_nfakes2, 
            selection12=ell_nfakes2, 
            selection21=disk_nfakes2_f20, 
            selection22=ell_nfakes2_f20, 
            selection31=disk_nfakes2_f45,
            selection32=ell_nfakes2_f45,
            lognorm=true, gridres=(20,10), 
            xmin=0, xmax=3.5, ymin=0, ymax=180,
            #ylabel="Flip [°]",
            #xlabel="log10( Merger Ratio )", 
            #title="Ellipticals after Flip ($(length(testas2["ϕ_flip"][ell_nfakes2])) / $(length(testas2["ϕ_flip"])))",
            #binned_median="x",
            )



#minimum(norm.(eachcol(testas2["j_main"][])))



fig, ax = subplots()
ax.hist(pyplottable(testas2["redshift"]), bins=50, rwidth=0.9, color="navy", edgecolor="black", alpha=1, zorder=3)
#ax.set_title("Short-term Spin Orientation Change")
ax.set_xlabel("z")
ax.set_ylabel("N")
#ax.set_yscale("log")
#ax.set_xscale("log")
ax.legend(loc="upper right", frameon=true, borderpad=1, handlelength=1.8)
ax.grid()
scale=0.5
fig.set_size_inches(16scale, 9scale)
fig.tight_layout()
fig.savefig("/home/moon/sfortune/spinevo/monday/update/hist_z.png", bbox_inches="tight", pad_inches=.1)
##########
# move to 2d
##########
snap = 92
plot_hexbins(log10.(pyplottable(testas2["j_main"])), 
            testas2["BVAL"] .+ 0.5 .*log10.(1 .+testas2["redshift"]) .+ ((2/3) .* log10.(testas2["M_felix"])), 
            #selection = findcs(pyplottable(testas2["snapNR"]), eq=snap),
            outfile="./monday/next/old_bval0_Mfelix_vs_j.png", 
            lognorm=true, gridres=(20,10), grid=true,
            #xmin=0, xmax=3.5, ymin=0, ymax=180,
            ylabel="b-value + 2/3*log(M)",
            xlabel="log( j_STARS )", 
            title="Uncorrected & Using M_tree",
            binned_median="x", #plot_bval="y"
            )

j_small = findcs(log10.(pyplottable(testas2["j_main"])), leq=2)
j_big   = findcs(log10.(pyplottable(testas2["j_main"])), gt=2)

z_lt05  = findcs(pyplottable(testas2["redshift"]), lt=0.5)
z_geq05_lt1 = findcs(pyplottable(testas2["redshift"]), geq=0.5, lt=1)
z_geq1  = findcs(pyplottable(testas2["redshift"]), geq=1)

rhea = findcs(testas2["BVAL"] .+ 0.5 .* log10.(1 .+testas2["redshift"]), geq=-3, comparewith=findcs(pyplottable(testas2["snapNR"]), eq=92))

# sfc: plot phase space and find relaxed systems

plot_group(92; rootID=testas2["ID_ISUB"][rhea][1], property="vel", stepsize=10, arsize=10, ptsize=5000, rad=100, res=(1920,1080), spin=false)
testas2["M_fromJ"][rhea]


##############
# Plot phase space # hexbin
##############

rootID = 3212
simbox=current_dir_simbox
r_by_rvir=0.1
scale=0.7
indir="/home/moon/sfortune/spinevo/data/halostories_v20211127_min0.0Gyr"
old_mcs = load(joinpath("/home/moon/sfortune/spinevo/data/halostories_v20211127_min0.0Gyr", find_felixID(rootID)[2]), "merger_collection_STARS")
for snapNR in old_mcs["snapNR"][1:end]
    println(snapNR)
    #snapNR = 128
    galaxyID_list   = find_merging_progenitors(snapNR; rootID=rootID, path_to_halostories=indir)
    snapshot    = Snapshot(simbox, snapNR)
    g           = Galaxy(snapshot, galaxyID_list[1])
    rvir_group  = read_galaxy_prop(get_group(g), "RVIR", :physical)
    read_halo!(g, units=:physical, props=((:stars, ["POS", "VEL", "MASS", "ID"]),), radius_units=:physical, radius=r_by_rvir*rvir_group)
    velocities  = norm.(eachcol(g.stars.vel))
    radii       = norm.(eachcol(g.stars.pos))
    println("Radius = $(r_by_rvir*rvir_group)")
    plot_hexbins(radii, 
                velocities, 
                weights=log10.(g.stars.mass),
                outfile="./test_weights/phase_space_hexbin_snap$(snapNR).png", 
                lognorm=false, gridres=(25, 15), 
                xlabel="r [kpc]", 
                ylabel="v [km/s]", 
                title="Snap $(snapNR), 0.1*rvir = $(round(r_by_rvir*rvir_group, digits=0))")
    
end


res = 5
plot_phasespace(; #snapNR = 128, 
    felixID=" ", rootID=3212, subID=" " ,
    parttype="STARS", r_by_rvir=0.1, #weights="MASS" ,
    lognorm=true, scale=15, gridres=(25*res, 15*res), colorlimits=(0.8,200), 
    outdir="/home/moon/sfortune/spinevo/plots/phase_space/halo_3212", outfile="PhSp", 
    )





#####################################
# plot galaxy arrows in box 
####################################


#pos1    = 
#pos2    =
#pos3    =
#vel1    =
#vel2    =
#vel3    =

stepsize = 1
JSServe.configure_server!(listen_port=1688, forwarded_port=1688)
set_theme!(resolution=(1600,900), backgroundcolor = :black)#, xgridcolor = :black)
scene = WGLMakie.arrows(
    as["SPOS"][1,snap136][1:stepsize:end], 
    as["SPOS"][2,snap136][1:stepsize:end], 
    as["SPOS"][3,snap136][1:stepsize:end],
    as["j_main"][1,snap136][1:stepsize:end], 
    as["j_main"][2,snap136][1:stepsize:end], 
    as["j_main"][3,snap136][1:stepsize:end],
    # convert kms/s to 10 kpc/Myr by 0.010227047347441423
    lengthscale = 3e0,
    arrowsize   = 5e1 .* log10.(as["M"][snap136][1:stepsize:end]),#4e-10*as["M"][snap136][1:stepsize:end],
    linewidth   = 2e2,
    linecolor   = (1e0/maximum(log10.(as["M"][snap136][1:stepsize:end]))) .* log10.(as["M"][snap136][1:stepsize:end]),#pyplottable(as["J_main"][:,snap136])[1:stepsize:end],
    arrowcolor  = (1e0/maximum(log10.(as["M"][snap136][1:stepsize:end]))) .* log10.(as["M"][snap136][1:stepsize:end]),#pyplottable(as["J_main"][:,snap136])[1:stepsize:end],
    colormap    = :autumn1,
    align       = :head,
    quality     = 3,
    shading     = false,
    ssao        = false,
    normalize   = false,
    transparency = true,
)
scene



# plot all the halo story datas
outdirfiles = readdir("/home/moon/sfortune/spinevo/plots/hs_data")
hs_dir = readdir(current_dir_stories)
for i in hs_dir
    rID = parse(Int, chop(i, head=5,tail=4))
    if sum(occursin.("$(rID)_halostory_data.png", outdirfiles)) == 0
        plot_halostory_data(outfile=joinpath("/home/moon/sfortune/spinevo/plots/hs_data", "halo_$(rID)_halostory_data.png"), snap=true, rootID=rID)
    else
        println("   Skipping $i")
    end
end












include("/home/moon/sfortune/spinevo/pkg/meta.jl")
testconnection()

h = 14587
s = 84
plot_group(s; rootID=h, 
       property="vel", stepsize=2, arsize=5, ptsize=5000, rad=300, res=(1920,1080), spin=true, 
       bgcolor=:black, indir=current_dir_stories, boxfix=false)




for i in 2:length(hs["lookbacktime"])
    if hs["snapNR"][i] > hs["snapNR"][i-1]
        println("       $(hs["switch"][i])   $(hs["mmp"][i])   $(hs["snapNR"][i])   $(hs["subID"][i])   $(hs["fileNR"][i])")#   $(hs["M_STARS"][i])")#   $(hs["lookbacktime"][i])")#$(hs["lookbacktime"][i]-hs["lookbacktime"][i-1])")
    else
        println("$(hs["switch"][i])   $(hs["mmp"][i])   $(hs["snapNR"][i])   $(hs["subID"][i])   $(hs["treeID"][i])   $(hs["fileNR"][i])")#   $(hs["M_STARS"][i])")#   $(hs["lookbacktime"][i])")#$(hs["lookbacktime"][i]-hs["lookbacktime"][i-1])")
    end
end

write_halo_stories(start=1, stop=1, central_switch=false, root_mass_thr=1e10, subtype = "central", root_list="/home/moon/sfortune/spinevo/data/root_list_20220306.jld", outdir="/home/moon/sfortune/spinevo/data/test")


for i in a
    println("$i   ---    $(length(findcs(as["redshift"], eq=i)))")
end

mergermap_indices(as; 
    condition="mergers",
    min=-1, max=1e16,
    mtype=convert_parttype_to_idx(str="stars"), snap=136,
    )
    
function quickndirty()
    stories         = readdir(current_dir_stories)
    bval_end        = Array{Float64}(undef, 0)
    nflips30        = Array{Int64}(undef, 0)
    nflips45        = Array{Int64}(undef, 0)
    nflips90        = Array{Int64}(undef, 0)
    nflips135       = Array{Int64}(undef, 0)
    G09_nflips30    = Array{Int64}(undef, 0)
    G09_nflips45    = Array{Int64}(undef, 0)
    G09_nflips90    = Array{Int64}(undef, 0)
    G09_nflips135   = Array{Int64}(undef, 0)
    nmergers2e8     = Array{Int64}(undef, 0)
    nmergers2e9     = Array{Int64}(undef, 0)
    nmergers2e10    = Array{Int64}(undef, 0)
    G09_nmergers2e8 = Array{Int64}(undef, 0)
    G09_nmergers2e9 = Array{Int64}(undef, 0)
    G09_nmergers2e10= Array{Int64}(undef, 0)

    sstacked03    = Array{Int64}(undef, 0)
    istacked03    = Array{Int64}(undef, 0)
    istacked09    = Array{Int64}(undef, 0)
    sstacked09    = Array{Int64}(undef, 0)
    for i in 1:length(stories)
        print("$(i)", ifelse(i % 30 == 0, "\n", " "))
        flush(stdout)
        mcs03 = load(joinpath("/home/moon/sfortune/spinevo/data/silvio_stories_softfix", stories[i]), "merger_collection_STARS")
        mcs09 = load(joinpath("/home/moon/sfortune/spinevo/data/silvio_stories_softfix_mintime09Gyr", stories[i]), "merger_collection_STARS")
        noswitch03  = findcs(mcs03["switch"         ], eq=0)
        noswitch09  = findcs(mcs09["switch"         ], eq=0)
        start_thr03 = findcs(mcs03["M"] .- mcs03["ΔM"], geq=2e10)
        start_thr03 = findcs(mcs09["M"] .- mcs09["ΔM"], geq=2e10)
        end_thr03   = findcs(mcs03["M"], geq=2e10)
        end_thr09   = findcs(mcs09["M"], geq=2e10)
        nofake03    = Array{Int64}(undef, 0)
        for ii in 1:length(mcs03["snapNR"])
            if fake_flip_finder(mcs03, ii; verbose=false) === missing
            elseif fake_flip_finder(mcs03, ii; verbose=false) == 0 && mcs03["snapNR"][ii] ∉ sstacked03[findcs(istacked03, eq = mcs03["subID"][ii])]
                nofake03 = vcat( nofake03, ii)
            end
        end
        nofake09    = Array{Int64}(undef, 0)
        for ii in 1:length(mcs09["snapNR"])
            if fake_flip_finder(mcs09, ii; verbose=false) === missing
            elseif fake_flip_finder(mcs09, ii; verbose=false) == 0 && mcs09["snapNR"][ii] ∉ sstacked09[findcs(istacked09, eq = mcs09["subID"][ii])]
                nofake09 = vcat( nofake09, ii)
            end
        end
        slct03 = idxmatch(idxmatch(idxmatch(noswitch03,nofake03),end_thr03), start_thr03)
        slct09 = idxmatch(idxmatch(idxmatch(noswitch09,nofake09),end_thr09), start_thr09)
        sstacked03   = vcat( sstacked03, mcs03["snapNR"][slct03] )
        sstacked09   = vcat( sstacked09, mcs09["snapNR"][slct09] )
        istacked09   = vcat( istacked09, mcs09["subID" ][slct09] )
        istacked03   = vcat( istacked03, mcs03["subID" ][slct03] )
        bval_end        = vcat( bval_end        , mcs03["BVAL_0"][end])
        nflips30        = vcat( nflips30        , length(findcs(mcs03["ϕ_flip"         ][slct03],geq=30  )) )
        nflips45        = vcat( nflips45        , length(findcs(mcs03["ϕ_flip"         ][slct03],geq=45  )) )
        nflips90        = vcat( nflips90        , length(findcs(mcs03["ϕ_flip"         ][slct03],geq=90  )) )
        nflips135       = vcat( nflips135       , length(findcs(mcs03["ϕ_flip"         ][slct03],geq=135 )) )
        G09_nflips30    = vcat( G09_nflips30    , length(findcs(mcs09["ϕ_flip"         ][slct09],geq=30  )) )
        G09_nflips45    = vcat( G09_nflips45    , length(findcs(mcs09["ϕ_flip"         ][slct09],geq=45  )) )
        G09_nflips90    = vcat( G09_nflips90    , length(findcs(mcs09["ϕ_flip"         ][slct09],geq=90  )) )
        G09_nflips135   = vcat( G09_nflips135   , length(findcs(mcs09["ϕ_flip"         ][slct09],geq=135 )) )
        nmergers2e8     = vcat( nmergers2e8     , length(findcs(mcs03["merger map"][9,slct03],geq=2e8 )) )
        nmergers2e9     = vcat( nmergers2e9     , length(findcs(mcs03["merger map"][9,slct03],geq=2e9 )) )
        nmergers2e10    = vcat( nmergers2e10    , length(findcs(mcs03["merger map"][9,slct03],geq=2e10)) )
        G09_nmergers2e8 = vcat( G09_nmergers2e8 , length(findcs(mcs09["merger map"][9,slct09],geq=2e8 )) )
        G09_nmergers2e9 = vcat( G09_nmergers2e9 , length(findcs(mcs09["merger map"][9,slct09],geq=2e9 )) )
        G09_nmergers2e10= vcat( G09_nmergers2e10, length(findcs(mcs09["merger map"][9,slct09],geq=2e10)) )
    end
    save("/home/moon/sfortune/spinevo/data/endstate.jld", 
        "bval_end",         bval_end,
        "nflips30",         nflips30,
        "nflips45",         nflips45,
        "nflips90",         nflips90,
        "nflips135",        nflips135,
        "G09_nflips30",     G09_nflips30,
        "G09_nflips45",     G09_nflips45,
        "G09_nflips90",     G09_nflips90,
        "G09_nflips135",    G09_nflips135,
        "nmergers2e8",      nmergers2e8,
        "nmergers2e9",      nmergers2e9,
        "nmergers2e10",     nmergers2e10,
        "G09_nmergers2e8",  G09_nmergers2e8,
        "G09_nmergers2e9",  G09_nmergers2e9,
        "G09_nmergers2e10", G09_nmergers2e10
        )
    println("done")
    return nothing
end

hs3212 = load("/home/moon/sfortune/spinevo/data/silvio_stories_spinmap/halo_3212.jld", "halo_story")
for i in 1:length(hs3212["snapFP"])
    #println("$(hs3212["snapFP"][i])   $(hs3212["snapNR"][i])   $(hs3212["treeID"][i])   $(hs3212["subID"][i])")
    println("$(hs3212["merger_spin_map_STARS"][1,i])")
end



j_test = get_merger_spins(trees, subID_map, 729075, 428458)



for i in 1:length(root_list["all_idx"][selection])
    if 1414 == get_global_index(subID_map, trees.halos[root_list["all_idx"][selection][i]])
        println("$i    $(get_global_index(subID_map, trees.halos[root_list["all_idx"][selection][i]]))")
    end
end



#box4_fulllbt    = Array{Float64}(undef, length(box4_fullsnaplist))
#box4_lbt        = Array{Float64}(undef, 0)
#for i in 1:length(box4_fullsnaplist)
#    head = read_header("$current_dir_simbox/groups_$(lpad(box4_fullsnaplist[i], 3, "0"))/sub_$(lpad(box4_fullsnaplist[i], 3, "0"))")
#    box4_fulllbt[i] = ustrip(lookback_time(cosmology(h=head.h0, OmegaM=head.omega_0), head.z))
#    if in(box4_fullsnaplist[i], box4_snaplist)
#        box4_lbt = vcat(box4_lbt, box4_fulllbt[i])
#    end
#end
#save(
#    "/home/moon/sfortune/spinevo/pkg/processed_variables.jld", 
#    "box4_fulllbt"  ,   box4_fulllbt,
#    "box4_lbt"      ,   box4_lbt
#    )
