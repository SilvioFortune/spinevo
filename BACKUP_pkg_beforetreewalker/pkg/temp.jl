
include("/home/moon/sfortune/spinevo/pkg/meta.jl")
assembly_STARS  = load("/home/moon/sfortune/spinevo/mergerimpact/assembly_Mstar_1.0e10.jld", "assembly_STARS")
centrals_136 = filter(assembly_STARS, condition="snap")


fig, ax = subplots()

ax.hist(assembly_STARS["M2_felix"][centrals_136["main"]], bins=100, label="M2_felix  at z ≈ $(@sprintf("%.2f", assembly_STARS["REDSHIFT"][centrals_136["main"]][1]))", rwidth=0.9, color="navy")


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

println("$(asm["ID_ISUB"][6230])   ---   $(asm["SNAP"][6230])")

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
for i in 1:length(hs["REDSHIFT"])
    head    = read_header("$current_dir_simbox/groups_$(@sprintf("%03i", hs["SNAP"][i]))/sub_$(@sprintf("%03i", hs["SNAP"][i]))")
    lbt_old = lbt
    lbt     = ustrip(lookback_time(cosmology(h=head.h0, OmegaM=head.omega_0), head.z))
    println("SNAP = $(hs["SNAP"][i])   ---   TIMEDIFF = $(lbt-lbt_old)")
    max_tdiff = ifelse(lbt-lbt_old > max_tdiff && hs["SNAP"][i] < 136, lbt-lbt_old, max_tdiff)
    min_tdiff = ifelse(0. < lbt-lbt_old < min_tdiff && hs["SNAP"][i] < 136, lbt-lbt_old, min_tdiff)
    #if lbt-lbt_old > max_tdiff
    #    max_tdiff = lbt-lbt_old
end
println(min_tdiff)
println(max_tdiff)

###################################################################################
###################################################################################
###################################################################################


for i in 1:100
    println("$(old_halo3212["SNAP"][i])      $(old_halo3212["M_STARS"][i])          $(old_halo3212["j_STARS"][:,i])    ---   $(halo3212["SNAP"][i])      $(halo3212["M_STARS"][i])          $(halo3212["j_STARS"][:,i])")
end


mcs3212 = load(joinpath(current_dir_jld, find_felixID(3212)[2]), "merger_collection_STARS")
testas = load("/home/moon/sfortune/spinevo/data/assembly_centrals2_20220112.jld", "assembly_STARS")
include("/home/moon/sfortune/spinevo/pkg/meta.jl")
testas2 = load("/home/moon/sfortune/spinevo/data/assembly_centrals2co_20220112.jld", "assembly_STARS")
testas2 = load("/home/moon/sfortune/spinevo/data/assembly_20220112_centrals2co_mfromj.jld", "assembly_STARS")

for i in 1:length(testas2["j_main"][1,:])
    if testas2["ID_ISUB"][i] == 3212.0
        println("$(fake_flip_finder(testas2, i, verbose = true))")
    end
end

jumps       = find_indices(testas["JUMP"])
n_jumps     = find_indices(testas["JUMP"], val=0)
exceeds     = find_indices(testas["EXCEED"])
n_exceeds   = find_indices(testas["EXCEED"], val=0)
switches    = find_indices(testas["SWITCH"])
n_switches  = find_indices(testas["SWITCH"], val=0)
borders     = find_indices(testas["BORDER"])
n_borders   = find_indices(testas["BORDER"], val=0)
fakes       = find_indices(testas["FAKEFLIP"])
n_fakes     = find_indices(testas["FAKEFLIP"], val=0)

jumps2       = find_indices(testas2["JUMP"])
n_jumps2     = find_indices(testas2["JUMP"], eq=0)
exceeds2     = find_indices(testas2["EXCEED"])
n_exceeds2   = find_indices(testas2["EXCEED"], eq=0)
switches2    = find_indices(testas2["SWITCH"])
n_switches2  = find_indices(testas2["SWITCH"], eq=0)
borders2     = find_indices(testas2["BORDER"])
n_borders2   = find_indices(testas2["BORDER"], eq=0)
fakes2       = find_indices(testas2["FAKEFLIP"])
n_fakes2     = find_indices(testas2["FAKEFLIP"], eq=0)


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
        println("$(testas2["ID_ISUB"][i])   ---   $(testas2["SNAP"][i])   ---   $(testas2["ϕ_flip"][i])   ---   $((log10.(replace(testas2["M2_felix"]./testas2["M2_MM"], missing => 1)))[i])")
    end
end

plot_group(80; lastID=2323, property="pos", stepsize=50, arsize=10, ptsize=5000, rad=200, res=(1600,900), spin=true)


# sfc: distinguish between b-values

b_disk  = -4.357
b_ell   = -4.732



result_disks    = find_indices(testas2["BVAL"] .+ 0.5 .*log10.(1 .+testas2["REDSHIFT"]), geq=b_disk)
disk_fakes2     = find_indices(testas2["FAKEFLIP"], comparewith = result_disks)
disk_nfakes2    = find_indices(testas2["FAKEFLIP"], eq=0, comparewith = result_disks)
result_int      = find_indices(testas2["BVAL"] .+ 0.5 .*log10.(1 .+testas2["REDSHIFT"]), gt=b_ell, lt=b_disk)
int_fakes2     = find_indices(testas2["FAKEFLIP"], comparewith = result_int)
int_nfakes2    = find_indices(testas2["FAKEFLIP"], eq=0, comparewith = result_int)
result_ell      = find_indices(testas2["BVAL"] .+ 0.5 .*log10.(1 .+testas2["REDSHIFT"]), leq=b_ell)
ell_fakes2     = find_indices(testas2["FAKEFLIP"], comparewith = result_ell)
ell_nfakes2    = find_indices(testas2["FAKEFLIP"], eq=0, comparewith = result_ell)

start_disk      = find_indices(testas2["BVAL"] .+ 0.5 .*log10.(1 .+testas2["REDSHIFT"]) .- testas2["δBVAL_0"], geq=b_disk)
start_int       = find_indices(testas2["BVAL"] .+ 0.5 .*log10.(1 .+testas2["REDSHIFT"]) .- testas2["δBVAL_0"], gt=b_ell, lt=b_disk)
start_ell       = find_indices(testas2["BVAL"] .+ 0.5 .*log10.(1 .+testas2["REDSHIFT"]) .- testas2["δBVAL_0"], leq=b_ell)
sdisk_nfakes2   = find_indices(testas2["FAKEFLIP"], eq=0, comparewith = start_disk)
sint_nfakes2    = find_indices(testas2["FAKEFLIP"], eq=0, comparewith = start_int)
sell_nfakes2    = find_indices(testas2["FAKEFLIP"], eq=0, comparewith = start_ell)
stay_disk_nfakes= sdisk_nfakes2[sdisk_nfakes2 .∈ Ref(Set(result_disks))]
stay_int_nfakes = sint_nfakes2[sint_nfakes2 .∈ Ref(Set(result_int))]
stay_ell_nfakes = sell_nfakes2[sell_nfakes2 .∈ Ref(Set(result_ell))]
disk_to_ell     = sdisk_nfakes2[sdisk_nfakes2 .∈ Ref(Set(result_ell))]
ell_to_disk     = sell_nfakes2[sell_nfakes2 .∈ Ref(Set(result_disks))]
int_to_disk     = sint_nfakes2[sint_nfakes2 .∈ Ref(Set(result_disks))]
int_to_ell      = sint_nfakes2[sint_nfakes2 .∈ Ref(Set(result_ell))]
sdisk_fakes2    = find_indices(testas2["FAKEFLIP"], eq=1, comparewith = start_disk)
sint_fakes2     = find_indices(testas2["FAKEFLIP"], eq=1, comparewith = start_int)
sell_fakes2     = find_indices(testas2["FAKEFLIP"], eq=1, comparewith = start_ell)
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

flip_gt20    = find_indices(testas2["ϕ_flip"], geq=20)
flip_gt45    = find_indices(testas2["ϕ_flip"], geq=45)

disk_nfakes2_f20    = disk_nfakes2[disk_nfakes2 .∈ Ref(Set(flip_gt20))]
disk_nfakes2_f45    = disk_nfakes2[disk_nfakes2 .∈ Ref(Set(flip_gt45))]
ell_nfakes2_f20    = ell_nfakes2[ell_nfakes2 .∈ Ref(Set(flip_gt20))]
ell_nfakes2_f45    = ell_nfakes2[ell_nfakes2 .∈ Ref(Set(flip_gt45))]



plot_hexbins_3x2(log10.(testas2["M2_felix"] ./ testas2["M2_MM"]), 
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
ax.hist(pyplottable(testas2["REDSHIFT"]), bins=50, rwidth=0.9, color="navy", edgecolor="black", alpha=1, zorder=3)
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
            testas2["BVAL"] .+ 0.5 .*log10.(1 .+testas2["REDSHIFT"]) .+ ((2/3) .* log10.(testas2["M_felix"])), 
            #selection = find_indices(pyplottable(testas2["SNAP"]), eq=snap),
            outfile="./monday/next/old_bval0_Mfelix_vs_j.png", 
            lognorm=true, gridres=(20,10), grid=true,
            #xmin=0, xmax=3.5, ymin=0, ymax=180,
            ylabel="b-value + 2/3*log(M)",
            xlabel="log( j_STARS )", 
            title="Uncorrected & Using M_tree",
            binned_median="x", #plot_bval="y"
            )

j_small = find_indices(log10.(pyplottable(testas2["j_main"])), leq=2)
j_big   = find_indices(log10.(pyplottable(testas2["j_main"])), gt=2)

z_lt05  = find_indices(pyplottable(testas2["REDSHIFT"]), lt=0.5)
z_geq05_lt1 = find_indices(pyplottable(testas2["REDSHIFT"]), geq=0.5, lt=1)
z_geq1  = find_indices(pyplottable(testas2["REDSHIFT"]), geq=1)

rhea = find_indices(testas2["BVAL"] .+ 0.5 .* log10.(1 .+testas2["REDSHIFT"]), geq=-3, comparewith=find_indices(pyplottable(testas2["SNAP"]), eq=92))

# sfc: plot phase space and find relaxed systems

plot_group(92; lastID=testas2["ID_ISUB"][rhea][1], property="vel", stepsize=10, arsize=10, ptsize=5000, rad=100, res=(1920,1080), spin=false)
testas2["M_fromJ"][rhea]


##############
# Plot phase space # hexbin
##############

lastID = 3212
simbox=current_dir_simbox
r_by_rvir=0.1
scale=0.7
indir="/home/moon/sfortune/spinevo/data/halostories_v20211127_min0.0Gyr"
old_mcs = load(joinpath("/home/moon/sfortune/spinevo/data/halostories_v20211127_min0.0Gyr", find_felixID(lastID)[2]), "merger_collection_STARS")
for snapNR in old_mcs["SNAP"][1:end]
    println(snapNR)
    #snapNR = 128
    galaxyID_list   = find_merging_progenitors(snapNR; lastID=lastID, path_to_halostories=indir)
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
                title="Snap $(snapNR), 0.1*rvir = $(round(r_by_rvir*rvir_group, digits=2))")
    
end



plot_phasespace(; snapNR = 128
    , felixID=" ", lastID=3212, subID=" "
    , parttype="STARS", r_by_rvir=0.1, weights="MASS"
    , lognorm=true, scale=15, gridres=(25, 15), colorlimits=(0.8,nothing)
    , outdir="./test_phasespace", outfile="OUT"
    )





#####################################
# plot galaxy arrows in box 
####################################

set_theme!(resolution=res, backgroundcolor = white)

pos1    =
pos2    =
pos3    =
vel1    =
vel2    =
vel3    =

scene = WGLMakie.arrows(
#arrows!(ax, 
    @view(ALL_STARS["POS"][1,all_notin][1:stepsize:end]), 
    @view(ALL_STARS["POS"][2,all_notin][1:stepsize:end]), 
    @view(ALL_STARS["POS"][3,all_notin][1:stepsize:end]),
    @view(ALL_STARS["VEL"][1,all_notin][1:stepsize:end]), 
    @view(ALL_STARS["VEL"][2,all_notin][1:stepsize:end]), 
    @view(ALL_STARS["VEL"][3,all_notin][1:stepsize:end]),
    # convert kms/s to 10 kpc/Myr by 0.010227047347441423
    lengthscale = 0.010227047347441423,#arsize*1e-2,#arsize * (maximum(ALL_STARS["POS"][1,all_notin][1,1:stepsize:end])-minimum(ALL_STARS["POS"][1,all_notin][1,1:stepsize:end])) / maximum(ALL_STARS["VEL"][:,all_notin][:,1:stepsize:end]),
    arrowsize   = (arsize/maximum(galaxies[1].stars.mass[1:stepsize:end])) .* ALL_STARS["MASS"][all_notin][1:stepsize:end],
    linewidth   = 0.1,
    linecolor   = @view(ALL_STARS["MASS"][all_notin][1:stepsize:end]) .* norm.(eachcol(@view(ALL_STARS["VEL"][:,all_notin][:,1:stepsize:end]))) ./ maximum(@view(galaxies[1].stars.mass[1:stepsize:end]) .* norm.(eachcol(@view(galaxies[1].stars.vel[:,1:stepsize:end])))),
    arrowcolor  = @view(ALL_STARS["MASS"][all_notin][1:stepsize:end]) .* norm.(eachcol(@view(ALL_STARS["VEL"][:,all_notin][:,1:stepsize:end]))) ./ maximum(@view(galaxies[1].stars.mass[1:stepsize:end]) .* norm.(eachcol(@view(galaxies[1].stars.vel[:,1:stepsize:end])))),
    colormap    = :autumn1,
    align       = :head,
    quality     = 3,
    shading     = false,
    ssao        = false,
    normalize   = false,
    transparency = true,
)
scene