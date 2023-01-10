
include("/home/moon/sfortune/spinevo/pkg/meta.jl")

as      = load("/home/moon/sfortune/spinevo/data/assembly_spinmap_03Gyr_patch.jld", "assembly_STARS")
ad      = load("/home/moon/sfortune/spinevo/data/assembly_spinmap_03Gyr_patch.jld", "assembly_DM")
as09    = load("/home/moon/sfortune/spinevo/data/assembly_spinmap_09Gyr_patch.jld", "assembly_STARS")
ad09    = load("/home/moon/sfortune/spinevo/data/assembly_spinmap_09Gyr_patch.jld", "assembly_DM")

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











G09_start_thr   = findcs(log10.(abs.(as09["M"] .- as09["ΔM"])), geq=log10.(2e10))
# borders: snap52=z1.18, snap100=z0.42
G09_z_hi        = findcs(as09["snapNR"      ], leq=52, comparewith=G09_start_thr)
G09_z_md        = findcs(as09["snapNR"      ], gt=52, lt=100, comparewith=G09_start_thr)
G09_z_lo        = findcs(as09["snapNR"      ], geq=100, comparewith=G09_start_thr)
G09_snap136     = findcs(as09["snapNR"      ], eq=136, comparewith=G09_start_thr)
G09_snap36      = findcs(as09["snapNR"      ], eq=36, comparewith=G09_start_thr)


G09_halo3212    = findcs(as09["ID_ISUB"      ], eq=3212, comparewith=G09_start_thr)
G09_halo1414    = findcs(as09["ID_ISUB"      ], eq=1414, comparewith=G09_start_thr)

G09_switches    = findcs(as09["switch"      ], eq=1, comparewith=G09_start_thr)
G09_n_switches  = findcs(as09["switch"      ], eq=0, comparewith=G09_start_thr)
G09_fakes       = findcs(as09["FAKEFLIP"    ], eq=1, comparewith=G09_start_thr)
G09_n_fakes     = findcs(as09["FAKEFLIP"    ], eq=0, comparewith=G09_start_thr)

G09_result_disk     = findcs(as09["BVAL"  ], geq=b_disk, comparewith=G09_start_thr)
G09_result_int      = findcs(as09["BVAL"  ], gt=b_ell, lt=b_disk, comparewith=G09_start_thr)
G09_result_ell      = findcs(as09["BVAL"  ], leq=b_ell, comparewith=G09_start_thr)

G09_start_disk      = findcs(as09["BVAL_0"  ] .- as09["ΔBVAL_0"    ], geq=b_disk, comparewith=G09_start_thr)
G09_start_int       = findcs(as09["BVAL_0"  ] .- as09["ΔBVAL_0"    ], gt=b_ell, lt=b_disk, comparewith=G09_start_thr)
G09_start_ell       = findcs(as09["BVAL_0"  ] .- as09["ΔBVAL_0"    ], leq=b_ell, comparewith=G09_start_thr)
G09_stay_disk       = G09_start_disk[  G09_start_disk  .∈  Ref(Set(G09_result_disk     ))]
G09_stay_int        = G09_start_int[    G09_start_int   .∈  Ref(Set(G09_result_int      ))]
G09_stay_ell        = G09_start_ell[    G09_start_ell   .∈  Ref(Set(G09_result_ell      ))]
G09_result_disk_fakes    = G09_fakes[   G09_fakes   .∈  Ref(Set(G09_result_disk     ))]
G09_result_int_fakes     = G09_fakes[   G09_fakes   .∈  Ref(Set(G09_result_int      ))]
G09_result_ell_fakes     = G09_fakes[   G09_fakes   .∈  Ref(Set(G09_result_ell      ))]
G09_result_disk_nfakes   = G09_n_fakes[ G09_n_fakes .∈  Ref(Set(G09_result_disk     ))]
G09_result_int_nfakes    = G09_n_fakes[ G09_n_fakes .∈  Ref(Set(G09_result_int      ))]
G09_result_ell_nfakes    = G09_n_fakes[ G09_n_fakes .∈  Ref(Set(G09_result_ell      ))]
G09_start_disk_fakes    = G09_fakes[   G09_fakes   .∈  Ref(Set(G09_start_disk     ))]
G09_start_int_fakes     = G09_fakes[   G09_fakes   .∈  Ref(Set(G09_start_int      ))]
G09_start_ell_fakes     = G09_fakes[   G09_fakes   .∈  Ref(Set(G09_start_ell      ))]
G09_start_disk_nfakes   = G09_n_fakes[ G09_n_fakes .∈  Ref(Set(G09_start_disk     ))]
G09_start_int_nfakes    = G09_n_fakes[ G09_n_fakes .∈  Ref(Set(G09_start_int      ))]
G09_start_ell_nfakes    = G09_n_fakes[ G09_n_fakes .∈  Ref(Set(G09_start_ell      ))]
G09_stay_disk_fakes     = G09_fakes[   G09_fakes   .∈  Ref(Set(G09_stay_disk     ))]
G09_stay_int_fakes      = G09_fakes[   G09_fakes   .∈  Ref(Set(G09_stay_int      ))]
G09_stay_ell_fakes      = G09_fakes[   G09_fakes   .∈  Ref(Set(G09_stay_ell      ))]
G09_stay_disk_nfakes    = G09_n_fakes[ G09_n_fakes .∈  Ref(Set(G09_stay_disk     ))]
G09_stay_int_nfakes     = G09_n_fakes[ G09_n_fakes .∈  Ref(Set(G09_stay_int      ))]
G09_stay_ell_nfakes     = G09_n_fakes[ G09_n_fakes .∈  Ref(Set(G09_stay_ell      ))]

G09_ell_to_disk_nfakes  = G09_start_ell_nfakes[G09_start_ell_nfakes     .∈ Ref(Set(G09_result_disk_nfakes))]
G09_ell_to_int_nfakes   = G09_start_ell_nfakes[G09_start_ell_nfakes     .∈ Ref(Set(G09_result_int_nfakes))]
G09_int_to_disk_nfakes  = G09_start_int_nfakes[G09_start_int_nfakes     .∈ Ref(Set(G09_result_disk_nfakes))]
G09_int_to_ell_nfakes   = G09_start_int_nfakes[G09_start_int_nfakes     .∈ Ref(Set(G09_result_ell_nfakes))]
G09_disk_to_ell_nfakes  = G09_start_disk_nfakes[G09_start_disk_nfakes   .∈ Ref(Set(G09_result_ell_nfakes))]
G09_disk_to_int_nfakes  = G09_start_disk_nfakes[G09_start_disk_nfakes   .∈ Ref(Set(G09_result_int_nfakes))]

G09_mergersDM_1to50         = findcs(ad09["M2"]./ad09["Mpeak_MM"], geq=1, leq=50, comparewith=G09_start_thr)

G09_mergerSTARS_acc_g10     = findcs(         as09["Mpeak_MERGERS"] ./ as09["M2"],       geq=0.1, comparewith=G09_start_thr)
G09_mergerSTARS_acc_g5_l10  = findcs(         as09["Mpeak_MERGERS"] ./ as09["M2"],       geq=0.05, lt=0.1, comparewith=G09_start_thr)
G09_mergerSTARS_acc_g0_l5   = findcs(         as09["Mpeak_MERGERS"] ./ as09["M2"],       gt=0.0, lt=0.05, comparewith=G09_start_thr)
G09_mergerSTARS_acc_eq0     = findcs(replace( as09["Mpeak_MERGERS"], missing => 0.0),  eq=0.0, comparewith=G09_start_thr)



mergerDM_acc_g10     = findcs(         ad["Mpeak_MERGERS"] ./ ad["M2"],       geq=0.1, comparewith=start_thr)
mergerDM_acc_g5_l10  = findcs(         ad["Mpeak_MERGERS"] ./ ad["M2"],       geq=0.05, lt=0.1, comparewith=start_thr)
mergerDM_acc_g0_l5   = findcs(         ad["Mpeak_MERGERS"] ./ ad["M2"],       gt=0.0, lt=0.05, comparewith=start_thr)
mergerDM_acc_eq0     = findcs(replace( ad["Mpeak_MERGERS"], missing => 0.0),  eq=0.0, comparewith=start_thr)
G09_mergerDM_acc_g10     = findcs(         ad09["Mpeak_MERGERS"] ./ ad09["M2"],       geq=0.1, comparewith=G09_start_thr)
G09_mergerDM_acc_g5_l10  = findcs(         ad09["Mpeak_MERGERS"] ./ ad09["M2"],       geq=0.05, lt=0.1, comparewith=G09_start_thr)
G09_mergerDM_acc_g0_l5   = findcs(         ad09["Mpeak_MERGERS"] ./ ad09["M2"],       gt=0.0, lt=0.05, comparewith=G09_start_thr)
G09_mergerDM_acc_eq0     = findcs(replace( ad09["Mpeak_MERGERS"], missing => 0.0),  eq=0.0, comparewith=G09_start_thr)


z_welker    = findcs(as["redshift"    ], geq=1.2, leq=3.8, comparewith=start_thr)
z_lower     = findcs(as["redshift"    ], lt=1.2, comparewith=start_thr)
z_higher    = findcs(as["redshift"    ], gt=3.8, comparewith=start_thr)
G09_z_welker    = findcs(as09["redshift"    ], geq=1.2, leq=3.8, comparewith=G09_start_thr)
G09_z_lower     = findcs(as09["redshift"    ], lt=1.2, comparewith=G09_start_thr)
G09_z_higher    = findcs(as09["redshift"    ], gt=3.8, comparewith=G09_start_thr)



z_prepeak   = findcs(as["redshift"    ], geq=1.9, comparewith=start_thr)
z_postpeak  = findcs(as["redshift"    ], lt=1.9, comparewith=start_thr)
G09_z_prepeak   = findcs(as09["redshift"    ], geq=1.9, comparewith=G09_start_thr)
G09_z_postpeak  = findcs(as09["redshift"    ], lt=1.9, comparewith=G09_start_thr)
z_pre047   = findcs(as["redshift"    ], geq=0.47, comparewith=start_thr)
z_post047  = findcs(as["redshift"    ], lt=0.47, comparewith=start_thr)
G09_z_pre047   = findcs(as09["redshift"    ], geq=0.47, comparewith=G09_start_thr)
G09_z_post047  = findcs(as09["redshift"    ], lt=0.47, comparewith=G09_start_thr)
z_pre064   = findcs(as["redshift"    ], geq=0.64, comparewith=start_thr)
z_post064  = findcs(as["redshift"    ], lt=0.64, comparewith=start_thr)
G09_z_pre064   = findcs(as09["redshift"    ], geq=0.64, comparewith=G09_start_thr)
G09_z_post064  = findcs(as09["redshift"    ], lt=0.64, comparewith=G09_start_thr)



# new
j_hi    = findcs(log10.(pyplottable(as["j_main"]  .- as["Δj_main"] )), gt=2.55, comparewith=start_thr)
j_md    = findcs(log10.(pyplottable(as["j_main"]  .- as["Δj_main"] )), gt=2.2 ,leq=2.55, comparewith=start_thr)
j_lo    = findcs(log10.(pyplottable(as["j_main"]  .- as["Δj_main"] )), leq=2.2, comparewith=start_thr)

J_hi    = findcs(log10.(pyplottable(as["J_main"]  .- as["ΔJ_main"] )), gt=12.9, comparewith=start_thr)
J_md    = findcs(log10.(pyplottable(as["J_main"]  .- as["ΔJ_main"] )), gt=12.6 ,leq=12.9, comparewith=start_thr)
J_lo    = findcs(log10.(pyplottable(as["J_main"]  .- as["ΔJ_main"] )), leq=12.6, comparewith=start_thr)

z_hirsch_0  = findcs(as["redshift"], leq=0.25)
z_hirsch_05 = findcs(as["redshift"], gt=0.25, lt=0.75)
z_hirsch_1  = findcs(as["redshift"], geq=0.75, lt=1.5)
z_hirsch_2  = findcs(as["redshift"], geq=1.5, lt=2.5)
z_hirsch_3  = findcs(as["redshift"], geq=2.5, lt=3.5)
z_hirsch_4  = findcs(as["redshift"], geq=3.5)

# SFR collection
sfr_coll    = Dict{String, Any}()
sfr_coll["subIDs"]          = Dict{Int64, Any}()
sfr_coll["redshift"]        = Array{Float64}(undef, 0)
sfr_coll["lookbacktime"]    = Array{Float64}(undef, 0)
sfr_coll["totSFR"]          = Array{Float64}(undef, 0)
sfr_coll["mean_sSFR"]       = Array{Float64}(undef, 0)
sfr_coll["N_halos"]         = Array{Int64}(undef, 0)
sfr_coll["totM"]            = Array{Float64}(undef, 0)
sfr_coll["boxsize"]         = Array{Float64}(undef, 0)
for i in sort(unique(as["snapNR"])) # snaps
    head    = read_header("$current_dir_simbox/groups_$(lpad(i, 3, "0"))/sub_$(lpad(i, 3, "0"))")
    idcs    = findcs(as["snapNR"], eq=i)#, comparewith=start_thr)
    #idcs    = Array{Int64}(undef, 0)
    #for ii in idcs_raw
    #    if !in(as["subID"][ii], as["subID"][idcs])
    #        idcs    = vcat(idcs, ii)
    #    end
    #end
    sfr_coll["subIDs"][i]   = as["subID"][idcs]
    sfr_coll["N_halos"]     = vcat(sfr_coll["N_halos"],     length(idcs) )
    sfr_coll["redshift"]    = vcat(sfr_coll["redshift"],    mean(as["redshift"  ][idcs]))
    sfr_coll["lookbacktime"]= vcat(sfr_coll["lookbacktime"],mean(as["lookbacktime"][idcs]))
    sfr_coll["boxsize"]     = vcat(sfr_coll["boxsize"],     prod(convert_units_physical([48000.0, 48000.0, 48000.0], :pos, head)) ) #kpc^3
    sfr_coll["totSFR"]      = vcat(sfr_coll["totSFR"],      sum( as["SFR"       ][idcs]))
    sfr_coll["mean_sSFR"]   = vcat(sfr_coll["mean_sSFR"],   mean(as["SFR"       ][idcs] ./ as["M"][idcs]))
    sfr_coll["totM"]        = vcat(sfr_coll["totM"],        sum( as["M"         ][idcs]))
end

for i in sort(unique(as["snapNR"]))
    if length(unique(sfr_coll["subIDs"][i])) != length(sfr_coll["subIDs"][i])
        println("\n$i   ---   $(length(unique(sfr_coll["subIDs"][i])))   ---   $(length(sfr_coll["subIDs"][i]))")
    else
        print("$i ")
    end
end


# global SFR collection
global_sfr_coll    = Dict{String, Any}()
global_sfr_coll["subIDs"]          = Dict{Int64, Any}()
global_sfr_coll["N_halos"]         = Array{Int64}(undef, 0)
global_sfr_coll["redshift"]        = Array{Float64}(undef, 0)
global_sfr_coll["lookbacktime"]    = Array{Float64}(undef, 0)
global_sfr_coll["totSFR"]          = Array{Float64}(undef, 0)
global_sfr_coll["mean_sSFR"]       = Array{Float64}(undef, 0)
global_sfr_coll["N_halos"]         = Array{Int64}(undef, 0)
global_sfr_coll["totM"]            = Array{Float64}(undef, 0)
global_sfr_coll["boxsize"]         = Array{Float64}(undef, 0)
sub1e10_sfr_coll    = Dict{String, Any}()
sub1e10_sfr_coll["subIDs"]          = Dict{Int64, Any}()
sub1e10_sfr_coll["N_halos"]         = Array{Int64}(undef, 0)
sub1e10_sfr_coll["redshift"]        = Array{Float64}(undef, 0)
sub1e10_sfr_coll["lookbacktime"]    = Array{Float64}(undef, 0)
sub1e10_sfr_coll["totSFR"]          = Array{Float64}(undef, 0)
sub1e10_sfr_coll["mean_sSFR"]       = Array{Float64}(undef, 0)
sub1e10_sfr_coll["N_halos"]         = Array{Int64}(undef, 0)
sub1e10_sfr_coll["totM"]            = Array{Float64}(undef, 0)
sub1e10_sfr_coll["boxsize"]         = Array{Float64}(undef, 0)
for i in 4:136 #sort(box4_snaplist) # snaps
    print("$i ")
    head    = read_header("$current_dir_simbox/groups_$(lpad(i, 3, "0"))/sub_$(lpad(i, 3, "0"))")
    ssfr    = read_subfind("$current_dir_simbox/groups_$(lpad(i,3,"0"))/sub_$(lpad(i,3,"0"))", "SSFR")
    mass    = convert_units_physical_mass(read_subfind("$current_dir_simbox/groups_$(lpad(i,3,"0"))/sub_$(lpad(i,3,"0"))", "SMST")[5,:], head)
    idx1e10  = findcs(mass, geq=1e10)
    global_sfr_coll["subIDs"][i]   = 0:length(mass)-1
    global_sfr_coll["N_halos"]     = vcat( global_sfr_coll["N_halos"],     length(ssfr) )
    global_sfr_coll["redshift"]    = vcat( global_sfr_coll["redshift"],    head.z )
    global_sfr_coll["lookbacktime"]= vcat( global_sfr_coll["lookbacktime"],ustrip(lookback_time(cosmology(h=head.h0, OmegaM=head.omega_0), head.z)) )
    global_sfr_coll["boxsize"]     = vcat( global_sfr_coll["boxsize"],     prod(convert_units_physical([48000.0, 48000.0, 48000.0], :pos, read_header("$current_dir_simbox/groups_$(lpad(i, 3, "0"))/sub_$(lpad(i, 3, "0"))"))) ) #kpc^3
    global_sfr_coll["totSFR"]      = vcat( global_sfr_coll["totSFR"],      sum(Float64.(ssfr)) )
    global_sfr_coll["totM"]        = vcat( global_sfr_coll["totM"],        Float64(sum(mass)) )
    global_sfr_coll["mean_sSFR"]   = vcat( global_sfr_coll["mean_sSFR"],   global_sfr_coll["totSFR"][end] / global_sfr_coll["totM"][end] )
    
    sub1e10_sfr_coll["subIDs"][i]    = idx1e10 .- 1
    sub1e10_sfr_coll["N_halos"]     = vcat( sub1e10_sfr_coll["N_halos"],     length(ssfr[idx1e10]) )
    sub1e10_sfr_coll["redshift"]    = vcat( sub1e10_sfr_coll["redshift"],    head.z )
    sub1e10_sfr_coll["lookbacktime"]= vcat( sub1e10_sfr_coll["lookbacktime"],ustrip(lookback_time(cosmology(h=head.h0, OmegaM=head.omega_0), head.z)) )
    sub1e10_sfr_coll["boxsize"]     = vcat( sub1e10_sfr_coll["boxsize"],     prod(convert_units_physical([48000.0, 48000.0, 48000.0], :pos, read_header("$current_dir_simbox/groups_$(lpad(i, 3, "0"))/sub_$(lpad(i, 3, "0"))"))) ) #kpc^3
    sub1e10_sfr_coll["totSFR"]      = vcat( sub1e10_sfr_coll["totSFR"],      sum(Float64.(ssfr[idx1e10])) )
    sub1e10_sfr_coll["totM"]        = vcat( sub1e10_sfr_coll["totM"],        Float64(sum(convert_units_physical_mass(read_subfind("$current_dir_simbox/groups_$(lpad(i,3,"0"))/sub_$(lpad(i,3,"0"))", "SMST")[5,:], head)[idx1e10])) )
    sub1e10_sfr_coll["mean_sSFR"]   = vcat( sub1e10_sfr_coll["mean_sSFR"],   sub1e10_sfr_coll["totSFR"][end] / sub1e10_sfr_coll["totM"][end] )
end

i = 136
println("$(sort(sfr_coll["subIDs"][i]))")

function sfrd_madau(z)
    # M_⊙ * yr-1 * Mpc-3 * h3
    return 0.015*((1+z)^2.7)/(1+((1+z)/2.9)^5.6)
end

uhr = CSV.read("/HydroSims/Magneticum/Box4/uhr_test/sfr.txt", DataFrame; delim=' ', ignorerepeated=true, header=0)
hr  = CSV.read("/HydroSims/Magneticum/Box4/hr_kobayashi/sfr.txt", DataFrame; delim=' ', ignorerepeated=true, header=0)

z_uhr       = (1 ./ uhr[:,:Column1]) .- 1
sfr_uhr     = uhr[:,:Column3]
z_hr        = (1 ./ hr[:,:Column1]) .- 1
sfr_hr      = hr[:,:Column3]

save("/home/moon/sfortune/spinevo/plots/magneticum.jld", 
    "global_sfr_coll",         global_sfr_coll,
    "sub1e10_sfr_coll",         sub1e10_sfr_coll,
    "z_uhr",         z_uhr,
    "sfr_uhr",         sfr_uhr,
    "z_hr",         z_hr,
    "sfr_hr",         sfr_hr
    )

z_madau     = LinRange(0, 12, 1000)
sfr_madau   = sfrd_madau.(z_madau)
idxmax_sfr_coll         = findcs(sfr_coll["totSFR"],  eq=maximum(sfr_coll["totSFR"]))[end]
idxmax_global_sfr_coll  = findcs(global_sfr_coll["totSFR"],   eq=maximum(global_sfr_coll["totSFR"]))[end]
idxmax_sub1e10_sfr_coll  = findcs(sub1e10_sfr_coll["totSFR"],   eq=maximum(sub1e10_sfr_coll["totSFR"]))[end]
idxmax_madau            = findcs(sfr_madau,  eq=maximum(sfr_madau))[end]
idxmax_uhr              = findcs(sfr_uhr,  eq=maximum(sfr_uhr))[end]
idxmax_hr               = findcs(sfr_hr,  eq=maximum(sfr_hr))[end]



mmi_as = mergermap_indices(as)



##############################################################################################################
##############################################################################################################
##############################################################################################################
##############################################################################################################
##############################################################################################################
##############################################################################################################
##############################################################################################################
##############################################################################################################




##############################################################################################################
# SFR / boxsize
#, fontsize=floor(0.5*scale)
#./ sfr_coll["boxsize"] ./ 1e9
PyPlot.matplotlib[:rc]("text", usetex=true)
PyPlot.matplotlib[:rc]("font", family="serif", serif=["Computer Modern"])
pointsize = 5
scale = 0.6
fig, ax = subplots()
ax.axvline(z_madau[idxmax_madau], color="black", label="Madau, z$(latexstring("\$_{peak}\$")) = $(round(z_madau[idxmax_madau], digits=1))")
ax.axvline(z_hr[idxmax_hr], color="darkred", label="Box4hr, z$(latexstring("\$_{peak}\$")) = $(round(z_hr[idxmax_hr], digits=1))")
ax.axvline(z_uhr[idxmax_uhr], color="red", label="Box4uhr, z$(latexstring("\$_{peak}\$")) = $(round(z_uhr[idxmax_uhr], digits=1))")
ax.axvline(global_sfr_coll["redshift"][idxmax_global_sfr_coll], color="mediumblue", label="All subhalos, z$(latexstring("\$_{peak}\$")) = $(round(global_sfr_coll["redshift"][idxmax_global_sfr_coll], digits=1))")
ax.axvline(sub1e10_sfr_coll["redshift"][idxmax_sub1e10_sfr_coll], color="dodgerblue", label="All 1e10-cut, z$(latexstring("\$_{peak}\$")) = $(round(sub1e10_sfr_coll["redshift"][idxmax_sub1e10_sfr_coll], digits=1))")
ax.axvline(sfr_coll["redshift"][idxmax_sfr_coll], color="cyan", label="Centrals 2e10-cut, z$(latexstring("\$_{peak}\$")) = $(round(sfr_coll["redshift"][idxmax_sfr_coll], digits=1))")
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
ax.legend(loc="upper right", frameon=true, borderpad=1, handlelength=1.8)
ax.grid(false)
fig.set_size_inches(16scale, 9scale)
fig.tight_layout()
fig.savefig("/home/moon/sfortune/spinevo/plots/meeting_20220405/SFRDcomoving_vs_redshift.png", bbox_inches="tight", pad_inches=.1)
close("all")
PyPlot.matplotlib[:rc]("text", usetex=false)




##############################################################################################################
# SFR histogram
##############################################################################################################
sfr_g0      = findcs(         as["SFR"],       gt=0.0, comparewith=start_thr)
G09_sfr_g0  = findcs(         as09["SFR"],       gt=0.0, comparewith=G09_start_thr)

sfr_hi      = findcs(         log10.(as["SFR"]),       geq=median(log10.(as["SFR"][sfr_g0])), comparewith=start_thr)
sfr_lo      = findcs(         log10.(as["SFR"]),        lt=median(log10.(as["SFR"][sfr_g0])), comparewith=start_thr)
G09_sfr_hi  = findcs(         log10.(as09["SFR"]),       geq=median(log10.(as09["SFR"][G09_sfr_g0])), comparewith=G09_start_thr)
G09_sfr_lo  = findcs(         log10.(as09["SFR"]),       lt=median(log10.(as09["SFR"][G09_sfr_g0])), comparewith=G09_start_thr)

fig, ax = subplots()
PyPlot.matplotlib[:rc]("text", usetex=true)
PyPlot.matplotlib[:rc]("font", family="serif", serif=["Computer Modern"])
ax.hist( pyplottable(log10.(as09["SFR"][G09_sfr_g0]))    , bins=50, rwidth=0.9    , alpha=0.5, color="darkred" )
ax.hist( pyplottable(log10.(as["SFR"][sfr_g0]))    , bins=50, rwidth=0.5    , alpha=0.5, color="mediumblue" )
ax.axvline( median(log10.(as09["SFR"][G09_sfr_g0])), linewidth=2, color ="darkred", label="$(latexstring("\$\\Delta\$"))t $(latexstring("\$\\sim\$")) 1 Gyr, "    )
ax.axvline( median(log10.(as["SFR"][sfr_g0])), linewidth=2, color ="mediumblue", label="$(latexstring("\$\\Delta\$"))t $(latexstring("\$\\sim\$")) 0.3 Gyr"    )
ax.set_xlabel("log( SFR [ $(latexstring("\$M_\\odot\\ /\\ yr\$"))] )")
#ax.set_xscale("log")
#ax.set_xlim([10, nothing])
ax.legend(loc="upper left", frameon=true, borderpad=1, handlelength=1.8)
ax.grid()
scale=0.5
fig.set_size_inches(16scale, 9scale)
fig.tight_layout()
fig.savefig("/home/moon/sfortune/spinevo/plots/meeting_20220405/histogram_sfr.png", bbox_inches="tight", pad_inches=.1)
close("all")
PyPlot.matplotlib[:rc]("text", usetex=false)




##############################################################################################################
# Mass histogram
##############################################################################################################
mS_hi      = findcs(         log10.(as["M"]),       geq=median(log10.(as["M"])), comparewith=start_thr)
mS_lo      = findcs(         log10.(as["M"]),        lt=median(log10.(as["M"])), comparewith=start_thr)
G09_mS_hi  = findcs(         log10.(as09["M"]),       geq=median(log10.(as09["M"])), comparewith=G09_start_thr)
G09_mS_lo  = findcs(         log10.(as09["M"]),       lt=median(log10.(as09["M"])), comparewith=G09_start_thr)
mD_hi      = findcs(         log10.(ad["M"]),       geq=median(log10.(ad["M"])), comparewith=start_thr)
mD_lo      = findcs(         log10.(ad["M"]),        lt=median(log10.(ad["M"])), comparewith=start_thr)
G09_mD_hi  = findcs(         log10.(ad09["M"]),       geq=median(log10.(ad09["M"])), comparewith=G09_start_thr)
G09_mD_lo  = findcs(         log10.(ad09["M"]),       lt=median(log10.(ad09["M"])), comparewith=G09_start_thr)

fig, ax = subplots()
PyPlot.matplotlib[:rc]("text", usetex=true)
PyPlot.matplotlib[:rc]("font", family="serif", serif=["Computer Modern"])
ax.hist( log10.(ad["M"])    , bins=50, rwidth=1    , alpha=0.8, color="navy" )
ax.hist( log10.(as["M"])    , bins=50, rwidth=1    , alpha=0.8, color="lightseagreen" )
ax.axvline( median(log10.(ad["M"])), linewidth=2, color ="navy", label="M$(latexstring("\$_{DM}\$"))")
ax.axvline( median(log10.(as["M"])), linewidth=2, color ="lightseagreen", label="M$(latexstring("\$_{*}\$"))")
ax.set_xlabel("log( M [ $(latexstring("M\$_\\odot\$"))] )")
#ax.set_xscale("log")
#ax.set_xlim([10, nothing])
ax.legend(loc="upper left", frameon=true, borderpad=1, handlelength=1.8)
ax.grid(false)
scale=0.5
fig.set_size_inches(16scale, 9scale)
fig.tight_layout()
fig.savefig("/home/moon/sfortune/spinevo/plots/meeting_20220405/histogram_M.png", bbox_inches="tight", pad_inches=.1)
close("all")
PyPlot.matplotlib[:rc]("text", usetex=false)




##############################################################################################################
# Merger Mass histogram
##############################################################################################################
mergerS_g0      = findcs(         log10.(as["merger map"][9,:]),       gt=0.0, comparewith=start_thr)
G09_mergerS_g0  = findcs(         log10.(as09["merger map"][9,:]),       gt=0.0, comparewith=start_thr)
mergerD_g0      = findcs(         log10.(ad["merger map"][9,:]),       gt=0.0, comparewith=start_thr)
G09_mergerD_g0  = findcs(         log10.(ad09["merger map"][9,:]),       gt=0.0, comparewith=start_thr)
mergerS_hi      = findcs(         log10.(as["merger map"][9,:]),       geq=median(log10.(as["merger map"][9,:])), comparewith=mergerS_g0)
mergerS_lo      = findcs(         log10.(as["merger map"][9,:]),        lt=median(log10.(as["merger map"][9,:])), comparewith=mergerS_g0)
G09_mergerS_hi  = findcs(         log10.(as09["merger map"][9,:]),       geq=median(log10.(as09["merger map"][9,:])), comparewith=G09_mergerS_g0)
G09_mergerS_lo  = findcs(         log10.(as09["merger map"][9,:]),       lt=median(log10.(as09["merger map"][9,:])), comparewith=G09_mergerS_g0)
mergerD_hi      = findcs(         log10.(ad["merger map"][9,:]),       geq=median(log10.(ad["merger map"][9,:])), comparewith=mergerD_g0)
mergerD_lo      = findcs(         log10.(ad["merger map"][9,:]),        lt=median(log10.(ad["merger map"][9,:])), comparewith=mergerD_g0)
G09_mergerD_hi  = findcs(         log10.(ad09["merger map"][9,:]),       geq=median(log10.(ad09["merger map"][9,:])), comparewith=G09_mergerD_g0)
G09_mergerD_lo  = findcs(         log10.(ad09["merger map"][9,:]),       lt=median(log10.(ad09["merger map"][9,:])), comparewith=G09_mergerD_g0)

fig, ax = subplots()
PyPlot.matplotlib[:rc]("text", usetex=true)
PyPlot.matplotlib[:rc]("font", family="serif", serif=["Computer Modern"])
ax.hist( pyplottable(log10.(ad["merger map"][9,mergerD_g0]))    , bins=50, rwidth=1    , alpha=0.8, color="navy" )
ax.hist( pyplottable(log10.(as["merger map"][9,mergerS_g0]))    , bins=50, rwidth=1    , alpha=0.8, color="lightseagreen" )
ax.axvline( median(log10.(ad["merger map"][9,mergerD_g0])), linewidth=2, color ="navy", label="M$(latexstring("\$_{DM}\$"))")
ax.axvline( median(log10.(as["merger map"][9,mergerS_g0])), linewidth=2, color ="lightseagreen", label="M$(latexstring("\$_{*}\$"))")
ax.set_xlabel("log( M$(latexstring("\$_{merger}\$")) [M$(latexstring("\$_\\odot\$"))] )")
#ax.set_xscale("log")
#ax.set_xlim([10, nothing])
ax.legend(loc="upper left", frameon=true, borderpad=1, handlelength=1.8)
ax.grid(false)
scale=0.5
fig.set_size_inches(16scale, 9scale)
fig.tight_layout()
fig.savefig("/home/moon/sfortune/spinevo/plots/meeting_20220405/histogram_M_MERGER.png", bbox_inches="tight", pad_inches=.1)
close("all")
PyPlot.matplotlib[:rc]("text", usetex=false)




##############################################################################################################
# Hirschmann 2014
##############################################################################################################
pltcm       = pyimport("matplotlib.cm")
pltcolors   = pyimport("matplotlib.colors")
colormap    = pltcm.get_cmap(name="gnuplot")
scale=0.6
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
ax.plot(bcenters_0 , log10.(n_0 .* (length(z_hirsch_0)/(48/little_h0)^3) ./ length(unique(as["redshift"][z_hirsch_0]))) , "-", linewidth=2, ms=10, label= "$(round(minimum(as["redshift"][z_hirsch_0]),digits=2)) $(latexstring("\$\\leq\$")) z $(latexstring("\$\\leq\$")) $(round(maximum(as["redshift"][z_hirsch_0]),digits=2)), N = $(length(z_hirsch_0))", color=colormap(0/maximum(as["redshift"]))   )
ax.plot(bcenters_05 ,log10.(n_05.*(length(z_hirsch_05)/(48/little_h0)^3) ./ length(unique(as["redshift"][z_hirsch_05]))) ,"-", linewidth=2, ms=10, label="$(round(minimum(as["redshift"][z_hirsch_05]),digits=2)) $(latexstring("\$\\leq\$")) z $(latexstring("\$\\leq\$")) $(round(maximum(as["redshift"][z_hirsch_05]),digits=2)), N = $(length(z_hirsch_05))", color=colormap(0.5/maximum(as["redshift"]))   )
ax.plot(bcenters_1 , log10.(n_1 .* (length(z_hirsch_1)/(48/little_h0)^3) ./ length(unique(as["redshift"][z_hirsch_1]))) , "-", linewidth=2, ms=10, label= "$(round(minimum(as["redshift"][z_hirsch_1]),digits=2)) $(latexstring("\$\\leq\$")) z $(latexstring("\$\\leq\$")) $(round(maximum(as["redshift"][z_hirsch_1]),digits=2)), N = $(length(z_hirsch_1))", color=colormap(1/maximum(as["redshift"]))   )
ax.plot(bcenters_2 , log10.(n_2 .* (length(z_hirsch_2)/(48/little_h0)^3) ./ length(unique(as["redshift"][z_hirsch_2]))) , "-", linewidth=2, ms=10, label= "$(round(minimum(as["redshift"][z_hirsch_2]),digits=2)) $(latexstring("\$\\leq\$")) z $(latexstring("\$\\leq\$")) $(round(maximum(as["redshift"][z_hirsch_2]),digits=2)), N = $(length(z_hirsch_2))", color=colormap(2/maximum(as["redshift"]))   )
ax.plot(bcenters_3 , log10.(n_3 .* (length(z_hirsch_3)/(48/little_h0)^3) ./ length(unique(as["redshift"][z_hirsch_3]))) , "-", linewidth=2, ms=10, label= "$(round(minimum(as["redshift"][z_hirsch_3]),digits=2)) $(latexstring("\$\\leq\$")) z $(latexstring("\$\\leq\$")) $(round(maximum(as["redshift"][z_hirsch_3]),digits=2)), N = $(length(z_hirsch_3))", color=colormap(3/maximum(as["redshift"]))   )
ax.plot(bcenters_4 , log10.(n_4 .* (length(z_hirsch_4)/(48/little_h0)^3) ./ length(unique(as["redshift"][z_hirsch_4]))) , "-", linewidth=2, ms=10, label= "$(round(minimum(as["redshift"][z_hirsch_4]),digits=2)) $(latexstring("\$\\leq\$")) z $(latexstring("\$\\leq\$")) $(round(maximum(as["redshift"][z_hirsch_4]),digits=2)), N = $(length(z_hirsch_4))", color=colormap(4/maximum(as["redshift"]))   )
ax.set_xlabel("log( M$(latexstring("\$_*\$")) [ M$(latexstring("\$_\\odot\$"))] )")
ax.set_ylabel("log( $(latexstring("\$\\Phi\$")) [ $(latexstring("\$(Mpc\\ /\\ h)^{-3}\$")) ] )")
ax.set_yscale("linear")
ax.set_xlim(8.5,12.5)
ax.set_ylim(-5.5,-1)
#ax.set_title("$(latexstring("\$\\Delta\$"))t $(latexstring("\$\\sim\$")) 0.3 Gyr within $(round(maximum(as["redshift"]),digits=2)) $(latexstring("\$\\geq\$")) z $(latexstring("\$\\geq\$")) $(round(minimum(as["redshift"]),digits=2))")
ax.legend(loc="lower left", frameon=true, borderpad=1, handlelength=1.8)
ax.grid(false)
fig.set_size_inches(10scale, 9scale)
fig.tight_layout()
fig.savefig("/home/moon/sfortune/spinevo/plots/meeting_20220405/stellarmassfct_Hirschmann2014_comparison.png", bbox_inches="tight", pad_inches=.1)
close("all")
PyPlot.matplotlib[:rc]("text", usetex=false)




##############################################################################################################
# Flipangle distribution
##############################################################################################################
#rand_res       = 100000
#rand_ang=ones(rand_res)
#for i in 1:rand_res
#    rand_ang[i] = 180/π * angle([rand(Uniform(-1,1)), rand(Uniform(-1,1)), rand(Uniform(-1,1))],[rand(Uniform(-1,1)), rand(Uniform(-1,1)), rand(Uniform(-1,1))])
#end
rand_res       = 10000
scale=0.5
N_bins=50
fig, ax = subplots()
n_00, x_00     = PyPlot.hist( pyplottable(as["ϕ_flip"   ][n_fakes     ]), bins=N_bins, density=true, log=false)
n_09, x_09     = PyPlot.hist( pyplottable(as09["ϕ_flip" ][G09_n_fakes ]), bins=N_bins, density=true, log=false)
#n_rd, x_rd     = PyPlot.hist( rand_ang,                                     bins=N_bins, density=true, log=false)
#n_00         .*= count( !isnan, pyplottable(as["ϕ_flip"   ][n_fakes     ]) )
#n_09         .*= count( !isnan, pyplottable(as09["ϕ_flip" ][G09_n_fakes ]) )
bcenters_00    = 0.5 .* (x_00[2:end]    .+ x_00[1:end-1])
bcenters_09    = 0.5 .* (x_09[2:end]    .+ x_09[1:end-1])
#bcenters_rd    = 0.5 .* (x_rd[2:end]    .+ x_rd[1:end-1])
bcenters_rd    = LinRange(0,180,rand_res)
n_rd      = sind.(bcenters_rd)
n_rd    .*= 0.5 / sum(sind.(bcenters_rd)) / N_bins * rand_res
clf()
fig, ax = subplots()
PyPlot.matplotlib[:rc]("text", usetex=true)
PyPlot.matplotlib[:rc]("font", family="serif", serif=["Computer Modern"])
#ax.plot(bcenters_rd , n_rd , "-", linewidth=3, label= "Random Flips"   , color="black")
ax.plot(bcenters_rd , n_rd , "-", linewidth=2, label= "Random Flips"   , color="black")
ax.plot(bcenters_00 , n_00 , "-", linewidth=2, label= "Instant Flips"   , color="mediumblue")
ax.plot(bcenters_09 , n_09 , "-", linewidth=2, label= "1 Gyr Flips"   , color="darkred")
ax.set_xlabel("$(latexstring("\$\\phi_{flip}\$")) [°]")
ax.set_ylabel("P")
ax.set_yscale("log")
ax.set_xlim(0,180)
ax.set_ylim(nothing,nothing)
#ax.set_title("$(latexstring("\$\\Delta\$"))t $(latexstring("\$\\sim\$")) 0.3 Gyr within $(round(maximum(as["redshift"]),digits=2)) $(latexstring("\$\\geq\$")) z $(latexstring("\$\\geq\$")) $(round(minimum(as["redshift"]),digits=2))")
ax.legend(loc="lower left", frameon=true, borderpad=1, handlelength=1.8)
ax.grid(false)
fig.set_size_inches(10scale, 9scale)
fig.tight_layout()
fig.savefig("/home/moon/sfortune/spinevo/plots/meeting_20220405/flipangle_dist.png", bbox_inches="tight", pad_inches=.1)
close("all")
PyPlot.matplotlib[:rc]("text", usetex=false)




##############################################################################################################
# Welker 2014 comparison 
##############################################################################################################
# welker z
e0 = idxmatch(mergerSTARS_acc_eq0       , z_welker)
g0 = idxmatch(mergerSTARS_acc_g0_l5     , z_welker)
g5 = idxmatch(mergerSTARS_acc_g5_l10    , z_welker)
g10= idxmatch(mergerSTARS_acc_g10       , z_welker)
scale=0.5
N_bins= 5
n_e0, x_e0      = PyPlot.hist( cosd.(pyplottable(as["ϕ_flip"][e0])), 
    range=(0,1),#maximum(cosd.(mergers_e0["ϕ_flip"]))), 
    bins=N_bins, density=true, log=false, rwidth=0.9, histtype="step", alpha=1, zorder=1)#, edgecolor="black")
n_g0, x_g0      = PyPlot.hist( cosd.(pyplottable(as["ϕ_flip"][g0])), 
    range=(0,1),#maximum(cosd.(mergers_g0["ϕ_flip"]))), 
    bins=N_bins, density=true, log=false, rwidth=0.9, histtype="step", alpha=1, zorder=1)#, edgecolor="black")
n_g5, x_g5      = PyPlot.hist( cosd.(pyplottable(as["ϕ_flip"][g5])), 
    range=(0,1),#maximum(cosd.(mergers_g5["ϕ_flip"]))), 
    bins=N_bins, density=true, log=false, rwidth=0.9, histtype="step", alpha=1, zorder=1)#, edgecolor="black")
n_g10, x_g10    = PyPlot.hist( cosd.(pyplottable(as["ϕ_flip"][g10])), 
    range=(0,1),#maximum(cosd.(mergers_g10["ϕ_flip"]))), 
    bins=N_bins, density=true,log=false,rwidth=0.9, histtype="step", alpha=1, zorder=1)#, edgecolor="black")
bcenters_e0     = 0.5 .* (x_e0[2:end]  .+ x_e0[1:end-1])
bcenters_g0     = 0.5 .* (x_g0[2:end]  .+ x_g0[1:end-1])
bcenters_g5     = 0.5 .* (x_g5[2:end]  .+ x_g5[1:end-1])
bcenters_g10    = 0.5 .* (x_g10[2:end] .+ x_g10[1:end-1])
clf()
PyPlot.matplotlib[:rc]("text", usetex=true)
PyPlot.matplotlib[:rc]("font", family="serif", serif=["Computer Modern"])
fig, ax = subplots()
ax.plot(bcenters_e0 , log10.(n_e0 ), ".-", ms=10, label="$(latexstring("\$\\delta\$"))m = 0\\%,       N = $(length(e0))" , color="red"   )
ax.plot(bcenters_g0 , log10.(n_g0 ), ".-", ms=10, label="5\\% $(latexstring("\$\\geq\$")) $(latexstring("\$\\delta\$"))m $(latexstring("\$\\geq\$")) 0\\%,  N = $(length(g0))" , color="orange")
ax.plot(bcenters_g5 , log10.(n_g5 ), ".-", ms=10, label="10\\% $(latexstring("\$\\geq\$")) $(latexstring("\$\\delta\$"))m $(latexstring("\$\\geq\$")) 5\\%, N = $(length(g5))" , color="purple")
ax.plot(bcenters_g10, log10.(n_g10), ".-", ms=10, label="$(latexstring("\$\\delta\$"))m $(latexstring("\$\\geq\$")) 10\\%,      N = $(length(g10))" , color="blue"  )
ax.set_xlabel("cos($(latexstring("\$\\phi_{flip}\$")))")
ax.set_ylabel("log( P )")
ax.set_yscale("linear")
ax.set_ylim(-3,nothing)
ax.set_title("$(latexstring("\$\\Delta\$"))t $(latexstring("\$\\sim\$")) 0.3 Gyr within $(round(maximum(as["redshift"][z_welker]),digits=2)) $(latexstring("\$\\geq\$")) z $(latexstring("\$\\geq\$")) $(round(minimum(as["redshift"][z_welker]),digits=2))")
ax.legend(loc="upper left", frameon=true, borderpad=1, handlelength=1.8)
ax.grid()
fig.set_size_inches(9scale, 9scale)
fig.tight_layout()
fig.savefig("/home/moon/sfortune/spinevo/plots/meeting_20220405/flipPDF_by_mergerRatioSTARS_Welker2014_zwelker_03.png", bbox_inches="tight", pad_inches=.1)
close("all")
PyPlot.matplotlib[:rc]("text", usetex=false)

# all z
e0 = idxmatch(mergerSTARS_acc_eq0       , start_thr)
g0 = idxmatch(mergerSTARS_acc_g0_l5     , start_thr)
g5 = idxmatch(mergerSTARS_acc_g5_l10    , start_thr)
g10= idxmatch(mergerSTARS_acc_g10       , start_thr)
scale=0.5
N_bins= 5
n_e0, x_e0      = PyPlot.hist( cosd.(pyplottable(as["ϕ_flip"][e0])), 
    range=(0,1),#maximum(cosd.(mergers_e0["ϕ_flip"]))), 
    bins=N_bins, density=true, log=false, rwidth=0.9, histtype="step", alpha=1, zorder=1)#, edgecolor="black")
n_g0, x_g0      = PyPlot.hist( cosd.(pyplottable(as["ϕ_flip"][g0])), 
    range=(0,1),#maximum(cosd.(mergers_g0["ϕ_flip"]))), 
    bins=N_bins, density=true, log=false, rwidth=0.9, histtype="step", alpha=1, zorder=1)#, edgecolor="black")
n_g5, x_g5      = PyPlot.hist( cosd.(pyplottable(as["ϕ_flip"][g5])), 
    range=(0,1),#maximum(cosd.(mergers_g5["ϕ_flip"]))), 
    bins=N_bins, density=true, log=false, rwidth=0.9, histtype="step", alpha=1, zorder=1)#, edgecolor="black")
n_g10, x_g10    = PyPlot.hist( cosd.(pyplottable(as["ϕ_flip"][g10])), 
    range=(0,1),#maximum(cosd.(mergers_g10["ϕ_flip"]))), 
    bins=N_bins, density=true,log=false,rwidth=0.9, histtype="step", alpha=1, zorder=1)#, edgecolor="black")
bcenters_e0     = 0.5 .* (x_e0[2:end]  .+ x_e0[1:end-1])
bcenters_g0     = 0.5 .* (x_g0[2:end]  .+ x_g0[1:end-1])
bcenters_g5     = 0.5 .* (x_g5[2:end]  .+ x_g5[1:end-1])
bcenters_g10    = 0.5 .* (x_g10[2:end] .+ x_g10[1:end-1])
clf()
PyPlot.matplotlib[:rc]("text", usetex=true)
PyPlot.matplotlib[:rc]("font", family="serif", serif=["Computer Modern"])
fig, ax = subplots()
ax.plot(bcenters_e0 , log10.(n_e0 ), ".-", ms=10, label="$(latexstring("\$\\delta\$"))m = 0\\%,       N = $(length(e0))" , color="red"   )
ax.plot(bcenters_g0 , log10.(n_g0 ), ".-", ms=10, label="5\\% $(latexstring("\$\\geq\$")) $(latexstring("\$\\delta\$"))m $(latexstring("\$\\geq\$")) 0\\%,  N = $(length(g0))" , color="orange")
ax.plot(bcenters_g5 , log10.(n_g5 ), ".-", ms=10, label="10\\% $(latexstring("\$\\geq\$")) $(latexstring("\$\\delta\$"))m $(latexstring("\$\\geq\$")) 5\\%, N = $(length(g5))" , color="purple")
ax.plot(bcenters_g10, log10.(n_g10), ".-", ms=10, label="$(latexstring("\$\\delta\$"))m $(latexstring("\$\\geq\$")) 10\\%,      N = $(length(g10))" , color="blue"  )
ax.set_xlabel("cos($(latexstring("\$\\phi_{flip}\$")))")
ax.set_ylabel("log( P )")
ax.set_yscale("linear")
ax.set_ylim(-3,nothing)
ax.set_title("$(latexstring("\$\\Delta\$"))t $(latexstring("\$\\sim\$")) 0.3 Gyr within $(round(maximum(as["redshift"][start_thr]),digits=2)) $(latexstring("\$\\geq\$")) z $(latexstring("\$\\geq\$")) $(round(minimum(as["redshift"][start_thr]),digits=2))")
ax.legend(loc="upper left", frameon=true, borderpad=1, handlelength=1.8)
ax.grid()
fig.set_size_inches(9scale, 9scale)
fig.tight_layout()
fig.savefig("/home/moon/sfortune/spinevo/plots/meeting_20220405/flipPDF_by_mergerRatioSTARS_Welker2014_zall_03.png", bbox_inches="tight", pad_inches=.1)
close("all")
PyPlot.matplotlib[:rc]("text", usetex=false)

# 09 Gyr flips
# welker z
e0 = idxmatch(G09_mergerSTARS_acc_eq0       , G09_z_welker)
g0 = idxmatch(G09_mergerSTARS_acc_g0_l5     , G09_z_welker)
g5 = idxmatch(G09_mergerSTARS_acc_g5_l10    , G09_z_welker)
g10= idxmatch(G09_mergerSTARS_acc_g10       , G09_z_welker)
scale=0.5
N_bins= 5
n_e0, x_e0      = PyPlot.hist( cosd.(pyplottable(as09["ϕ_flip"][e0])), 
    range=(0,1),#maximum(cosd.(mergers_e0["ϕ_flip"]))), 
    bins=N_bins, density=true, log=false, rwidth=0.9, histtype="step", alpha=1, zorder=1)#, edgecolor="black")
n_g0, x_g0      = PyPlot.hist( cosd.(pyplottable(as09["ϕ_flip"][g0])), 
    range=(0,1),#maximum(cosd.(mergers_g0["ϕ_flip"]))), 
    bins=N_bins, density=true, log=false, rwidth=0.9, histtype="step", alpha=1, zorder=1)#, edgecolor="black")
n_g5, x_g5      = PyPlot.hist( cosd.(pyplottable(as09["ϕ_flip"][g5])), 
    range=(0,1),#maximum(cosd.(mergers_g5["ϕ_flip"]))), 
    bins=N_bins, density=true, log=false, rwidth=0.9, histtype="step", alpha=1, zorder=1)#, edgecolor="black")
n_g10, x_g10    = PyPlot.hist( cosd.(pyplottable(as09["ϕ_flip"][g10])), 
    range=(0,1),#maximum(cosd.(mergers_g10["ϕ_flip"]))), 
    bins=N_bins, density=true,log=false,rwidth=0.9, histtype="step", alpha=1, zorder=1)#, edgecolor="black")
bcenters_e0     = 0.5 .* (x_e0[2:end]  .+ x_e0[1:end-1])
bcenters_g0     = 0.5 .* (x_g0[2:end]  .+ x_g0[1:end-1])
bcenters_g5     = 0.5 .* (x_g5[2:end]  .+ x_g5[1:end-1])
bcenters_g10    = 0.5 .* (x_g10[2:end] .+ x_g10[1:end-1])
clf()
PyPlot.matplotlib[:rc]("text", usetex=true)
PyPlot.matplotlib[:rc]("font", family="serif", serif=["Computer Modern"])
fig, ax = subplots()
ax.plot(bcenters_e0 , log10.(n_e0 ), ".-", ms=10, label="$(latexstring("\$\\delta\$"))m = 0\\%,       N = $(length(e0))" , color="red"   )
ax.plot(bcenters_g0 , log10.(n_g0 ), ".-", ms=10, label="5\\% $(latexstring("\$\\geq\$")) $(latexstring("\$\\delta\$"))m $(latexstring("\$\\geq\$")) 0\\%,  N = $(length(g0))" , color="orange")
ax.plot(bcenters_g5 , log10.(n_g5 ), ".-", ms=10, label="10\\% $(latexstring("\$\\geq\$")) $(latexstring("\$\\delta\$"))m $(latexstring("\$\\geq\$")) 5\\%, N = $(length(g5))" , color="purple")
ax.plot(bcenters_g10, log10.(n_g10), ".-", ms=10, label="$(latexstring("\$\\delta\$"))m $(latexstring("\$\\geq\$")) 10\\%,      N = $(length(g10))" , color="blue"  )
ax.set_xlabel("cos($(latexstring("\$\\phi_{flip}\$")))")
ax.set_ylabel("log( P )")
ax.set_yscale("linear")
ax.set_ylim(-3,nothing)
ax.set_title("$(latexstring("\$\\Delta\$"))t $(latexstring("\$\\sim\$")) 1 Gyr within $(round(maximum(as09["redshift"][G09_z_welker]),digits=2)) $(latexstring("\$\\geq\$")) z $(latexstring("\$\\geq\$")) $(round(minimum(as09["redshift"][G09_z_welker]),digits=2))")
ax.legend(loc="upper left", frameon=true, borderpad=1, handlelength=1.8)
ax.grid()
fig.set_size_inches(9scale, 9scale)
fig.tight_layout()
fig.savefig("/home/moon/sfortune/spinevo/plots/meeting_20220405/flipPDF_by_mergerRatioSTARS_Welker2014_zwelker_09.png", bbox_inches="tight", pad_inches=.1)
close("all")
PyPlot.matplotlib[:rc]("text", usetex=false)

# all z
e0 = idxmatch(G09_mergerSTARS_acc_eq0       , G09_start_thr)
g0 = idxmatch(G09_mergerSTARS_acc_g0_l5     , G09_start_thr)
g5 = idxmatch(G09_mergerSTARS_acc_g5_l10    , G09_start_thr)
g10= idxmatch(G09_mergerSTARS_acc_g10       , G09_start_thr)
scale=0.5
N_bins= 5
n_e0, x_e0      = PyPlot.hist( cosd.(pyplottable(as09["ϕ_flip"][e0])), 
    range=(0,1),#maximum(cosd.(mergers_e0["ϕ_flip"]))), 
    bins=N_bins, density=true, log=false, rwidth=0.9, histtype="step", alpha=1, zorder=1)#, edgecolor="black")
n_g0, x_g0      = PyPlot.hist( cosd.(pyplottable(as09["ϕ_flip"][g0])), 
    range=(0,1),#maximum(cosd.(mergers_g0["ϕ_flip"]))), 
    bins=N_bins, density=true, log=false, rwidth=0.9, histtype="step", alpha=1, zorder=1)#, edgecolor="black")
n_g5, x_g5      = PyPlot.hist( cosd.(pyplottable(as09["ϕ_flip"][g5])), 
    range=(0,1),#maximum(cosd.(mergers_g5["ϕ_flip"]))), 
    bins=N_bins, density=true, log=false, rwidth=0.9, histtype="step", alpha=1, zorder=1)#, edgecolor="black")
n_g10, x_g10    = PyPlot.hist( cosd.(pyplottable(as09["ϕ_flip"][g10])), 
    range=(0,1),#maximum(cosd.(mergers_g10["ϕ_flip"]))), 
    bins=N_bins, density=true,log=false,rwidth=0.9, histtype="step", alpha=1, zorder=1)#, edgecolor="black")
bcenters_e0     = 0.5 .* (x_e0[2:end]  .+ x_e0[1:end-1])
bcenters_g0     = 0.5 .* (x_g0[2:end]  .+ x_g0[1:end-1])
bcenters_g5     = 0.5 .* (x_g5[2:end]  .+ x_g5[1:end-1])
bcenters_g10    = 0.5 .* (x_g10[2:end] .+ x_g10[1:end-1])
clf()
PyPlot.matplotlib[:rc]("text", usetex=true)
PyPlot.matplotlib[:rc]("font", family="serif", serif=["Computer Modern"])
fig, ax = subplots()
ax.plot(bcenters_e0 , log10.(n_e0 ), ".-", ms=10, label="$(latexstring("\$\\delta\$"))m = 0\\%,       N = $(length(e0))" , color="red"   )
ax.plot(bcenters_g0 , log10.(n_g0 ), ".-", ms=10, label="5\\% $(latexstring("\$\\geq\$")) $(latexstring("\$\\delta\$"))m $(latexstring("\$\\geq\$")) 0\\%,  N = $(length(g0))" , color="orange")
ax.plot(bcenters_g5 , log10.(n_g5 ), ".-", ms=10, label="10\\% $(latexstring("\$\\geq\$")) $(latexstring("\$\\delta\$"))m $(latexstring("\$\\geq\$")) 5\\%, N = $(length(g5))" , color="purple")
ax.plot(bcenters_g10, log10.(n_g10), ".-", ms=10, label="$(latexstring("\$\\delta\$"))m $(latexstring("\$\\geq\$")) 10\\%,      N = $(length(g10))" , color="blue"  )
ax.set_xlabel("cos($(latexstring("\$\\phi_{flip}\$")))")
ax.set_ylabel("log( P )")
ax.set_yscale("linear")
ax.set_ylim(-3,nothing)
ax.set_title("$(latexstring("\$\\Delta\$"))t $(latexstring("\$\\sim\$")) 1 Gyr within $(round(maximum(as09["redshift"][G09_start_thr]),digits=2)) $(latexstring("\$\\geq\$")) z $(latexstring("\$\\geq\$")) $(round(minimum(as09["redshift"][G09_start_thr]),digits=2))")
ax.legend(loc="upper left", frameon=true, borderpad=1, handlelength=1.8)
ax.grid()
fig.set_size_inches(9scale, 9scale)
fig.tight_layout()
fig.savefig("/home/moon/sfortune/spinevo/plots/meeting_20220405/flipPDF_by_mergerRatioSTARS_Welker2014_zall_09.png", bbox_inches="tight", pad_inches=.1)
close("all")
PyPlot.matplotlib[:rc]("text", usetex=false)




##############################################################################################################
# flipangle delta b-value
##############################################################################################################
plot_hexbins(as["ΔBVAL_0"], as["ϕ_flip"]; selection=n_fakes, 
    outfile="/home/moon/sfortune/spinevo/plots/meeting_20220405/flipangle_vs_dbval_03.png",
    ymin=0, ymax=180, xmin=-1, xmax=1, title="$(latexstring("\$\\Delta\$"))t $(latexstring("\$\\sim\$")) 0.3 Gyr",
    scale=15, gridres=(30,13), xscale="linear", yscale="linear", 
    ylabel="$(latexstring("\$\\phi_{flip}\$")) [°]", xlabel="$(latexstring("\$\\Delta\$"))b-value", grid=false,
    lognorm=true, colorlimits=(0.8,nothing),
    cmap="BuPu",
    binned_median = "x", plot_bval=" ",
    weights=" ",
    calc_norm=true,    # @pyplottable
    )

# 0.9 Gyr
plot_hexbins(as09["ΔBVAL_0"], as09["ϕ_flip"]; selection=G09_n_fakes, 
    outfile="/home/moon/sfortune/spinevo/plots/meeting_20220405/flipangle_vs_dbval_09.png",
    ymin=0, ymax=180, xmin=-1, xmax=1, title="$(latexstring("\$\\Delta\$"))t $(latexstring("\$\\sim\$")) 1 Gyr",
    scale=15, gridres=(30,13), xscale="linear", yscale="linear", 
    ylabel="$(latexstring("\$\\phi_{flip}\$")) [°]", xlabel="$(latexstring("\$\\Delta\$"))b-value", grid=false,
    lognorm=true, colorlimits=(0.8,nothing),
    cmap="BuPu",
    binned_median = "x", plot_bval=" ",
    weights=" ",
    calc_norm=true,    # @pyplottable
    )




##############################################################################################################
# flipangle delta j
##############################################################################################################
plot_hexbins(as["ϕ_flip"], log10.(pyplottable(as["Δj_main"])); selection=n_fakes, 
    outfile="/home/moon/sfortune/spinevo/plots/meeting_20220405/djSTARS_vs_flipangle_03.png",
    xmin=0, xmax=180, ymin=" ", ymax=" ", title="$(latexstring("\$\\Delta\$"))t $(latexstring("\$\\sim\$")) 0.3 Gyr",
    scale=15, gridres=(30,13), xscale="linear", yscale="linear", 
    xlabel="$(latexstring("\$\\phi_{flip}\$")) [°]", ylabel="log( $(latexstring("\$\\Delta\$"))j [kpc km/s] )", grid=false,
    lognorm=true, colorlimits=(0.8,nothing),
    cmap="BuPu",
    binned_median = "x", plot_bval=" ",
    weights=" ",
    calc_norm=true,    # @pyplottable
    )

# 0.9 Gyr
plot_hexbins(as09["ϕ_flip"], log10.(pyplottable(as09["Δj_main"])); selection=G09_n_fakes, 
    outfile="/home/moon/sfortune/spinevo/plots/meeting_20220405/djSTARS_vs_flipangle_09.png",
    xmin=0, xmax=180, ymin=" ", ymax=" ", title="$(latexstring("\$\\Delta\$"))t $(latexstring("\$\\sim\$")) 1 Gyr",
    scale=15, gridres=(30,13), xscale="linear", yscale="linear", 
    xlabel="$(latexstring("\$\\phi_{flip}\$")) [°]", ylabel="log( $(latexstring("\$\\Delta\$"))j [kpc km/s] )", grid=false,
    lognorm=true, colorlimits=(0.8,nothing),
    cmap="BuPu",
    binned_median = "x", plot_bval=" ",
    weights=" ",
    calc_norm=true,    # @pyplottable
    )




##############################################################################################################
# flipangle delta J
##############################################################################################################
plot_hexbins(as["ϕ_flip"], log10.(pyplottable(as["ΔJ_main"])); selection=n_fakes, 
    outfile="/home/moon/sfortune/spinevo/plots/meeting_20220405/dJ_STARS_vs_flipangle_03.png",
    xmin=0, xmax=180, ymin=" ", ymax=" ", title="$(latexstring("\$\\Delta\$"))t $(latexstring("\$\\sim\$")) 0.3 Gyr",
    scale=15, gridres=(30,13), xscale="linear", yscale="linear", 
    xlabel="$(latexstring("\$\\phi_{flip}\$")) [°]", ylabel="log( $(latexstring("\$\\Delta\$"))J [$(latexstring("\$M_\\odot\\ kpc\\ km/s\$"))] )", grid=false,
    lognorm=true, colorlimits=(0.8,nothing),
    cmap="BuPu",
    binned_median = "x", plot_bval=" ",
    weights=" ",
    calc_norm=true,    # @pyplottable
    )

# 0.9 Gyr
plot_hexbins(as09["ϕ_flip"], log10.(pyplottable(as09["ΔJ_main"])); selection=G09_n_fakes, 
    outfile="/home/moon/sfortune/spinevo/plots/meeting_20220405/dJ_STARS_vs_flipangle_09.png",
    xmin=0, xmax=180, ymin=" ", ymax=" ", title="$(latexstring("\$\\Delta\$"))t $(latexstring("\$\\sim\$")) 1 Gyr",
    scale=15, gridres=(30,13), xscale="linear", yscale="linear", 
    xlabel="$(latexstring("\$\\phi_{flip}\$")) [°]", ylabel="log( $(latexstring("\$\\Delta\$"))J [$(latexstring("\$M_\\odot\\ kpc\\ km/s\$"))] )", grid=false,
    lognorm=true, colorlimits=(0.8,nothing),
    cmap="BuPu",
    binned_median = "x", plot_bval=" ",
    weights=" ",
    calc_norm=true,    # @pyplottable
    )




##############################################################################################################
# flipangle delta j rel
##############################################################################################################
plot_hexbins(as["ϕ_flip"], log10.(pyplottable(as["Δj_main"] ./ (as["j_main"].-as["Δj_main"]))); selection=n_fakes, 
    outfile="/home/moon/sfortune/spinevo/plots/meeting_20220405/djrel1_vs_flipangle_03.png",
    xmin=0, xmax=180, ymin=-2, ymax=4, title="$(latexstring("\$\\Delta\$"))t $(latexstring("\$\\sim\$")) 0.3 Gyr",
    scale=15, gridres=(30,13), xscale="linear", yscale="linear", 
    xlabel="$(latexstring("\$\\phi_{flip}\$")) [°]", ylabel="log( $(latexstring("\$\\Delta\$"))j / j1 )", 
    grid=true,
    lognorm=true, colorlimits=(0.8,nothing),
    cmap="BuPu",
    binned_median = "x", plot_bval=" ",
    weights=" ",
    calc_norm=true,    # @pyplottable
    )

# 0.9 Gyr
plot_hexbins(as09["ϕ_flip"], log10.(pyplottable(as09["Δj_main"] ./ (as09["j_main"].-as09["Δj_main"]))); selection=G09_n_fakes, 
    outfile="/home/moon/sfortune/spinevo/plots/meeting_20220405/djrel1_vs_flipangle_09.png",
    xmin=0, xmax=180, ymin=-2, ymax=4, title="$(latexstring("\$\\Delta\$"))t $(latexstring("\$\\sim\$")) 1 Gyr",
    scale=15, gridres=(30,13), xscale="linear", yscale="linear", 
    xlabel="$(latexstring("\$\\phi_{flip}\$")) [°]", ylabel="log( $(latexstring("\$\\Delta\$"))j / j1 )", 
    grid=true,
    lognorm=true, colorlimits=(0.8,nothing),
    cmap="BuPu",
    binned_median = "x", plot_bval=" ",
    weights=" ",
    calc_norm=true,    # @pyplottable
    )




##############################################################################################################
# d bval delta j rel
##############################################################################################################
plot_hexbins(as["ΔBVAL_0"], log10.(pyplottable(as["Δj_main"] ./ (as["j_main"].-as["Δj_main"]))); selection=n_fakes, 
    outfile="/home/moon/sfortune/spinevo/plots/meeting_20220405/djrel1_vs_dbval_03.png",
    xmin=" ", xmax=" ", ymin=-2, ymax=4, title="$(latexstring("\$\\Delta\$"))t $(latexstring("\$\\sim\$")) 0.3 Gyr",
    scale=15, gridres=(30,13), xscale="linear", yscale="linear", 
    xlabel="$(latexstring("\$\\Delta\$"))b-value", ylabel="log( $(latexstring("\$\\Delta\$"))j / j1 )", 
    grid=true,
    lognorm=true, colorlimits=(0.8,nothing),
    cmap="BuPu",
    binned_median = "x", plot_bval=" ",
    weights=" ",
    calc_norm=true,    # @pyplottable
    )

# 0.9 Gyr
plot_hexbins(as09["ΔBVAL_0"], log10.(pyplottable(as09["Δj_main"] ./ (as09["j_main"].-as09["Δj_main"]))); selection=G09_n_fakes, 
    outfile="/home/moon/sfortune/spinevo/plots/meeting_20220405/djrel1_vs_dbval_09.png",
    xmin=" ", xmax=" ", ymin=-2, ymax=4, title="$(latexstring("\$\\Delta\$"))t $(latexstring("\$\\sim\$")) 1 Gyr",
    scale=15, gridres=(30,13), xscale="linear", yscale="linear", 
    xlabel="$(latexstring("\$\\Delta\$"))b-value", ylabel="log( $(latexstring("\$\\Delta\$"))j / j1 )", 
    grid=true,
    lognorm=true, colorlimits=(0.8,nothing),
    cmap="BuPu",
    binned_median = "x", plot_bval=" ",
    weights=" ",
    calc_norm=true,    # @pyplottable
    )




##############################################################################################################
# d bval delta j abs
##############################################################################################################
plot_hexbins(as["ΔBVAL_0"], log10.(pyplottable(as["Δj_main"])); selection=n_fakes, 
    outfile="/home/moon/sfortune/spinevo/plots/meeting_20220405/djSTARS_vs_dbval_03.png",
    xmin=-2, xmax=2, ymin=-0.5, ymax=4, title="$(latexstring("\$\\Delta\$"))t $(latexstring("\$\\sim\$")) 0.3 Gyr",
    scale=15, gridres=(30,13), xscale="linear", yscale="linear", 
    xlabel="$(latexstring("\$\\Delta\$"))b-value", ylabel="log( $(latexstring("\$\\Delta\$"))j [kpc km/s] )", 
    grid=true,
    lognorm=true, colorlimits=(0.8,nothing),
    cmap="BuPu",
    binned_median = "x", plot_bval=" ",
    weights=" ",
    calc_norm=true,    # @pyplottable
    )

# 0.9 Gyr
plot_hexbins(as09["ΔBVAL_0"], log10.(pyplottable(as09["Δj_main"])); selection=G09_n_fakes, 
    outfile="/home/moon/sfortune/spinevo/plots/meeting_20220405/djSTARS_vs_dbval_09.png",
    xmin=-2, xmax=2, ymin=-0.5, ymax=4, title="$(latexstring("\$\\Delta\$"))t $(latexstring("\$\\sim\$")) 1 Gyr",
    scale=15, gridres=(30,13), xscale="linear", yscale="linear", 
    xlabel="$(latexstring("\$\\Delta\$"))b-value", ylabel="log( $(latexstring("\$\\Delta\$"))j [kpc km/s] )", 
    grid=true,
    lognorm=true, colorlimits=(0.8,nothing),
    cmap="BuPu",
    binned_median = "x", plot_bval=" ",
    weights=" ",
    calc_norm=true,    # @pyplottable
    )






##############################################################################################################
# flipangle N mergers
##############################################################################################################
mmi_as      = mergermap_indices(as)
mmi_as09    = mergermap_indices(as09)
n_all       = zeros(Int64,length(as["M"]))
n_g2e9      = zeros(Int64,length(as["M"]))
n_g2e10     = zeros(Int64,length(as["M"]))
G09_n_all   = zeros(Int64,length(as09["M"]))
G09_n_g2e9  = zeros(Int64,length(as09["M"]))
G09_n_g2e10 = zeros(Int64,length(as09["M"]))
for i in 1:length(mmi_as["main"])
    n_all[mmi_as["main"][i]]    = length(findcs(as["merger map"][9,mmi_as["first"][i]:mmi_as["last"][i]], geq=1e8 ))
    n_g2e9[mmi_as["main"][i]]   = length(findcs(as["merger map"][9,mmi_as["first"][i]:mmi_as["last"][i]], geq=2e9 ))
    n_g2e10[mmi_as["main"][i]]  = length(findcs(as["merger map"][9,mmi_as["first"][i]:mmi_as["last"][i]], geq=2e10))
end
for i in 1:length(mmi_as09["main"])
    G09_n_all[mmi_as09["main"][i]]    = length(findcs(as09["merger map"][9,mmi_as09["first"][i]:mmi_as09["last"][i]], geq=1e8 ))
    G09_n_g2e9[mmi_as09["main"][i]]   = length(findcs(as09["merger map"][9,mmi_as09["first"][i]:mmi_as09["last"][i]], geq=2e9 ))
    G09_n_g2e10[mmi_as09["main"][i]]  = length(findcs(as09["merger map"][9,mmi_as09["first"][i]:mmi_as09["last"][i]], geq=2e10))
end

x   = vcat(as["ϕ_flip"], as["ϕ_flip"], as["ϕ_flip"], as09["ϕ_flip"], as09["ϕ_flip"], as09["ϕ_flip"])
y   = vcat(n_all, n_g2e9, n_g2e10, G09_n_all, G09_n_g2e9, G09_n_g2e10)
miny= minimum(pyplottable(y)[findall(x->x .!== NaN, pyplottable(y))])
maxy= maximum(pyplottable(y)[findall(x->x .!== NaN, pyplottable(y))])
minx= 0
maxx= 180
slcts       = Dict{Int64, Vector{Int64}}()
slcts[1]    = 1:length(n_all)
slcts[2]    = length(n_all)+1:2*length(n_all)
slcts[3]    = 2*length(n_all)+1:3*length(n_all)
slcts[4]    = 3*length(n_all)+1:3*length(n_all)+length(G09_n_all)
slcts[5]    = 3*length(n_all)+length(G09_n_all)+1:3*length(n_all)+2*length(G09_n_all)
slcts[6]    = 3*length(n_all)+2*length(G09_n_all)+1:3*length(n_all)+3*length(G09_n_all)
labels      = Dict{Int64, String}()
labels[1]   = "$(latexstring("\$\\Delta\$"))t $(latexstring("\$\\sim\$")) 0.3 Gyr, All Mergers"
labels[2]   = "log( M$(latexstring("\$_*\$")) [M$(latexstring("\$_\\odot\$"))] ) $(latexstring("\$\\geq\$")) $(round(log10(2e9), digits=2))"
labels[3]   = "log( M$(latexstring("\$_*\$")) [M$(latexstring("\$_\\odot\$"))] ) $(latexstring("\$\\geq\$")) $(round(log10(2e10), digits=2))"
labels[4]   = "$(latexstring("\$\\Delta\$"))t $(latexstring("\$\\sim\$")) 1 Gyr, All Mergers"
labels[5]   = "log( M$(latexstring("\$_*\$")) [M$(latexstring("\$_\\odot\$"))] ) $(latexstring("\$\\geq\$")) $(round(log10(2e9), digits=2))"
labels[6]   = "log( M$(latexstring("\$_*\$")) [M$(latexstring("\$_\\odot\$"))] ) $(latexstring("\$\\geq\$")) $(round(log10(2e10), digits=2))"
labpos      = Array{String}(undef, 6)
labpos     .= "upper right"
plot_hexbins_nxm( x, 
            y, 
            slcts, labels, 3, 2,
            #scale=0.0,
            cbarshift=2.5,
            outfile="/home/moon/sfortune/spinevo/plots/meeting_20220405/N_mergers_vs_flipangle_dt-M-split.png", 
            lognorm=true, gridres=(30,10), label_pos=labpos,
            ymin=miny, ymax=maxy,
            xmin=minx, xmax=maxx, 
            ylabel="N$(latexstring("\$_{mergers}\$"))",
            xlabel="$(latexstring("\$\\phi_{flip}\$"))", fontsize=30,
            binned_median="y", plot_bval=" ",
            )




##############################################################################################################
# delta j by flipangle
##############################################################################################################
x   = as["ϕ_flip"]
y   = log10.(pyplottable(as["Δj_main"]))
miny= minimum(pyplottable(y)[findall(x->x .!== NaN, pyplottable(y))])
maxy= maximum(pyplottable(y)[findall(x->x .!== NaN, pyplottable(y))])
minx= 0
maxx= 180
slcts       = Dict{Int64, Vector{Int64}}()
slcts[1]    = result_ell_nfakes
slcts[2]    = idxmatch(mS_lo, n_fakes)
slcts[3]    = idxmatch(sfr_lo , n_fakes)
slcts[4]    = result_disk_nfakes
slcts[5]    = idxmatch(mS_hi, n_fakes) 
slcts[6]    = idxmatch(sfr_hi , n_fakes)
labels      = Dict{Int64, String}()
labels[1]   = "Ellipticals, N = $(length(slcts[1]))"
labels[2]   = "log( M$(latexstring("\$_*\$")) [M$(latexstring("\$_\\odot\$"))] ) $(latexstring("\$\\leq\$")) $(round(median(log10.(as["M"])), digits=2)), N = $(length(slcts[2]))"
labels[3]   = "log( SFR [M$(latexstring("\$_\\odot\\ /\\ yr\$"))] ) $(latexstring("\$\\leq\$")) $(round(median(log10.(as["SFR"])), digits=2)), N = $(length(slcts[3]))"
labels[4]   = "Disks, N = $(length(slcts[4]))"
labels[5]   = "log( M$(latexstring("\$_*\$")) [M$(latexstring("\$_\\odot\$"))] ) $(latexstring("\$\\geq\$")) $(round(median(log10.(as["M"])), digits=2)), N = $(length(slcts[5]))"
labels[6]   = "log( SFR [M$(latexstring("\$_\\odot\\ /\\ yr\$"))] ) $(latexstring("\$\\geq\$")) $(round(median(log10.(as["SFR"])), digits=2)), N = $(length(slcts[6]))"
plot_hexbins_nxm( x, 
            y, 
            slcts, labels, 3, 2,
            #scale=0.0,
            cbarshift=2.5,
            outfile="/home/moon/sfortune/spinevo/plots/meeting_20220405/dj_STARS_vs_flipangle_b-M-sfr-split_03Gyr.png", 
            lognorm=true, gridres=(30,10), 
            ymin=miny, ymax=maxy,
            xmin=minx, xmax=maxx, 
            ylabel="log( $(latexstring("\$\\Delta\$"))j [kpc km/s] )",
            xlabel="$(latexstring("\$\\phi_{flip}\$"))", fontsize=30,
            binned_median="x", plot_bval=" ",
            )

# 0.9 Gyr
x   = as09["ϕ_flip"]
y   = log10.(pyplottable(as09["Δj_main"]))
miny= minimum(pyplottable(y)[findall(x->x .!== NaN, pyplottable(y))])
maxy= maximum(pyplottable(y)[findall(x->x .!== NaN, pyplottable(y))])
minx= 0
maxx= 180
slcts       = Dict{Int64, Vector{Int64}}()
slcts[1]    = G09_result_ell_nfakes
slcts[2]    = idxmatch(G09_mS_lo, G09_n_fakes)
slcts[3]    = idxmatch(G09_sfr_lo , G09_n_fakes)
slcts[4]    = G09_result_disk_nfakes
slcts[5]    = idxmatch(G09_mS_hi, G09_n_fakes) 
slcts[6]    = idxmatch(G09_sfr_hi , G09_n_fakes)
labels      = Dict{Int64, String}()
labels[1]   = "Ellipticals, N = $(length(slcts[1]))"
labels[2]   = "log( M$(latexstring("\$_*\$")) [M$(latexstring("\$_\\odot\$"))] ) $(latexstring("\$\\leq\$")) $(round(median(log10.(as09["M"])), digits=2)), N = $(length(slcts[2]))"
labels[3]   = "log( SFR [M$(latexstring("\$_\\odot\\ /\\ yr\$"))] ) $(latexstring("\$\\leq\$")) $(round(median(log10.(as09["SFR"])), digits=2)), N = $(length(slcts[3]))"
labels[4]   = "Disks, N = $(length(slcts[4]))"
labels[5]   = "log( M$(latexstring("\$_*\$")) [M$(latexstring("\$_\\odot\$"))] ) $(latexstring("\$\\geq\$")) $(round(median(log10.(as09["M"])), digits=2)), N = $(length(slcts[5]))"
labels[6]   = "log( SFR [M$(latexstring("\$_\\odot\\ /\\ yr\$"))] ) $(latexstring("\$\\geq\$")) $(round(median(log10.(as09["SFR"])), digits=2)), N = $(length(slcts[6]))"
plot_hexbins_nxm( x, 
            y, 
            slcts, labels, 3, 2,
            #scale=0.0,
            cbarshift=2.5,
            outfile="/home/moon/sfortune/spinevo/plots/meeting_20220405/dj_STARS_vs_flipangle_b-M-sfr-split_09Gyr.png", 
            lognorm=true, gridres=(30,10), 
            ymin=miny, ymax=maxy,
            xmin=minx, xmax=maxx, 
            ylabel="log( $(latexstring("\$\\Delta\$"))j [kpc km/s] )",
            xlabel="$(latexstring("\$\\phi_{flip}\$"))", fontsize=30,
            binned_median="x", plot_bval=" ",
            )






##############################################################################################################
# M_star vs flipangle
##############################################################################################################
plot_hexbins(as["ϕ_flip"], log10.(pyplottable(as["M"])); selection=n_fakes, 
    outfile="/home/moon/sfortune/spinevo/plots/meeting_20220405/Mstars_vs_flipangle_03.png",
    xmin=0, xmax=180, ymin=" ", ymax=" ", title="$(latexstring("\$\\Delta\$"))t $(latexstring("\$\\sim\$")) 0.3 Gyr",
    scale=15, gridres=(30,13), xscale="linear", yscale="linear", 
    xlabel="$(latexstring("\$\\phi_{flip}\$")) [°]", ylabel="log( M$(latexstring("\$_*\$")) [M$(latexstring("\$_\\odot\$"))] )", 
    grid=false,
    lognorm=true, colorlimits=(0.8,nothing),
    cmap="BuPu",
    binned_median = "x", plot_bval=" ",
    weights=" ",
    calc_norm=true,    # @pyplottable
    )

# 0.9 Gyr
plot_hexbins(as09["ϕ_flip"], log10.(pyplottable(as09["M"])); selection=G09_n_fakes, 
    outfile="/home/moon/sfortune/spinevo/plots/meeting_20220405/Mstars_vs_flipangle_09.png",
    xmin=0, xmax=180, ymin=" ", ymax=" ", title="$(latexstring("\$\\Delta\$"))t $(latexstring("\$\\sim\$")) 1 Gyr",
    scale=15, gridres=(30,13), xscale="linear", yscale="linear", 
    xlabel="$(latexstring("\$\\phi_{flip}\$")) [°]", ylabel="log( M$(latexstring("\$_*\$")) [M$(latexstring("\$_\\odot\$"))] )", 
    grid=false,
    lognorm=true, colorlimits=(0.8,nothing),
    cmap="BuPu",
    binned_median = "x", plot_bval=" ",
    weights=" ",
    calc_norm=true,    # @pyplottable
    )




##############################################################################################################
# d j vs d M_stars
##############################################################################################################
# combined
y   = vcat(log10.(pyplottable(as["Δj_main"])), log10.(pyplottable(as["Δj_main"])), log10.(pyplottable(as["Δj_main"])), log10.(pyplottable(as09["Δj_main"])), log10.(pyplottable(as09["Δj_main"])), log10.(pyplottable(as09["Δj_main"])))
x   = vcat(log10.(abs.(pyplottable(as["ΔM"]))), log10.(abs.(pyplottable(as["Mpeak_MERGERS"]))), log10.(abs.(pyplottable(as["Mpeak_MM"]))), log10.(abs.(pyplottable(as09["ΔM"]))), log10.(abs.(pyplottable(as09["Mpeak_MERGERS"]))), log10.(abs.(pyplottable(as09["Mpeak_MM"]))))
miny= minimum(pyplottable(y)[findall(x->x .!== NaN, pyplottable(y))])
maxy= maximum(pyplottable(y)[findall(x->x .!== NaN, pyplottable(y))])
minx= minimum(pyplottable(x)[findall(x->x .!== NaN, pyplottable(x))])
maxx= maximum(pyplottable(x)[findall(x->x .!== NaN, pyplottable(x))])
slcts       = Dict{Int64, Vector{Int64}}()
slcts[1]    = 1:length(as["ΔM"])
slcts[2]    = length(as["ΔM"])+1:2*length(as["ΔM"])
slcts[3]    = 2*length(as["ΔM"])+1:3*length(as["ΔM"])
slcts[4]    = 3*length(as["ΔM"])+1:3*length(as["ΔM"])+length(as09["ΔM"])
slcts[5]    = 3*length(as["ΔM"])+length(as09["ΔM"])+1:3*length(as["ΔM"])+2*length(as09["ΔM"])
slcts[6]    = 3*length(as["ΔM"])+2*length(as09["ΔM"])+1:3*length(as["ΔM"])+3*length(as09["ΔM"])
labels      = Dict{Int64, String}()
labels[1]   = "$(latexstring("\$\\Delta\$"))t $(latexstring("\$\\sim\$")) 0.3 Gyr, $(latexstring("\$\\Delta\$"))M"
labels[2]   = "$(latexstring("\$\\Delta\$"))t $(latexstring("\$\\sim\$")) 0.3 Gyr, M$(latexstring("\$_{mergers}\$"))"
labels[3]   = "$(latexstring("\$\\Delta\$"))t $(latexstring("\$\\sim\$")) 0.3 Gyr, M$(latexstring("\$_{MM}\$"))"
labels[4]   = "$(latexstring("\$\\Delta\$"))t $(latexstring("\$\\sim\$")) 1 Gyr, $(latexstring("\$\\Delta\$"))M"
labels[5]   = "$(latexstring("\$\\Delta\$"))t $(latexstring("\$\\sim\$")) 1 Gyr, M$(latexstring("\$_{mergers}\$"))"
labels[6]   = "$(latexstring("\$\\Delta\$"))t $(latexstring("\$\\sim\$")) 1 Gyr, M$(latexstring("\$_{MM}\$"))"
labpos      = Array{String}(undef, 6)
labpos     .= "upper left"
plot_hexbins_nxm( x, 
            y, 
            slcts, labels, 3, 2,
            #scale=0.0,
            cbarshift=2.5,
            outfile="/home/moon/sfortune/spinevo/plots/meeting_20220405/dj_vs_dM_dt-Mtype-split.png", 
            lognorm=true, gridres=(30,10), label_pos=labpos,
            ymin=miny, ymax=maxy,
            xmin=minx, xmax=maxx, 
            ylabel="log( $(latexstring("\$\\Delta\$"))j [kpc km/s] )",
            xlabel="log( M$(latexstring("\$_*\$")) [M$(latexstring("\$_\\odot\$"))] )", fontsize=30,
            binned_median="x", plot_bval=" ",
            )




##############################################################################################################
# b-value vs N_neighbors
##############################################################################################################
neigh_l2g8      = Array{Int64}(undef, 0)
neigh_g2g8l2e9  = Array{Int64}(undef, 0)
neigh_g2g9l2e10 = Array{Int64}(undef, 0)
neigh_g2g10     = Array{Int64}(undef, 0)
println(length(as["snapNR"]))
for i in 1:length(as["snapNR"])
    println("$i")
    groupID = Int(get_group(Galaxy(Snapshot(current_dir_simbox, as["snapNR"][i]), as["subID"][i])).igroup)
    head    = read_header("$current_dir_simbox/groups_$(lpad(as["snapNR"][i], 3, "0"))/sub_$(lpad(as["snapNR"][i], 3, "0"))")
    nsub    = read_subfind("$current_dir_simbox/groups_$(lpad(as["snapNR"][i],3,"0"))/sub_$(lpad(as["snapNR"][i],3,"0"))", "NSUB")
    smst    = convert_units_physical_mass(read_subfind("$current_dir_simbox/groups_$(lpad(as["snapNR"][i],3,"0"))/sub_$(lpad(as["snapNR"][i],3,"0"))", "SMST")[5,:], head)
    neigh_l2g8      = vcat( neigh_l2g8     , length(findcs(smst[2+sum(nsub[1:groupID]):sum(nsub[1:groupID+1])], leq=2e8)) )
    neigh_g2g8l2e9  = vcat( neigh_g2g8l2e9 , length(findcs(smst[2+sum(nsub[1:groupID]):sum(nsub[1:groupID+1])], gt=2e8, lt=2e9)) )
    neigh_g2g9l2e10 = vcat( neigh_g2g9l2e10, length(findcs(smst[2+sum(nsub[1:groupID]):sum(nsub[1:groupID+1])], geq=2e9, lt=2e10)) )
    neigh_g2g10     = vcat( neigh_g2g10    , length(findcs(smst[2+sum(nsub[1:groupID]):sum(nsub[1:groupID+1])], geq=2e10)) )
end

#sfc: error in following plot
x   = vcat(pyplottable(log.(neigh_l2g8)), pyplottable(log.(neigh_g2g8l2e9)), pyplottable(log.(neigh_g2g9l2e10)), pyplottable(log.(neigh_g2g10)))
y   = vcat(as["BVAL_0"], as["BVAL_0"], as["BVAL_0"], as["BVAL_0"])
miny= minimum(pyplottable(y)[findall(x->x .!== NaN, pyplottable(y))])
maxy= maximum(pyplottable(y)[findall(x->x .!== NaN, pyplottable(y))])
minx= minimum(pyplottable(x)[findall(x->x .!== NaN, pyplottable(x))])
maxx= maximum(pyplottable(x)[findall(x->x .!== NaN, pyplottable(x))])
slcts       = Dict{Int64, Vector{Int64}}()
slcts[1]    = 1:length(as["BVAL_0"])
slcts[2]    = 1+length(as["BVAL_0"]):2*length(as["BVAL_0"])
slcts[3]    = 1+2*length(as["BVAL_0"]):3*length(as["BVAL_0"]) 
slcts[4]    = 1+3*length(as["BVAL_0"]):4*length(as["BVAL_0"])
labels      = Dict{Int64, String}()
labels[1]   = "N( M$(latexstring("\$_*\\leq\\ 2e8\\ M_\\odot\$")) )"
labels[2]   = "N( $(latexstring("\$2e8\\ \\lt\$")) M$(latexstring("\$_*\\ \\lt\\ 2e9\\ M_\\odot\$")) )"
labels[3]   = "N( $(latexstring("\$2e9\\ \\leq\$")) M$(latexstring("\$_*\\ \\lt\\ 2e10\\ M_\\odot\$")) )"
labels[4]   = "N( $(latexstring("\$2e10\\ M_\\odot\\ \\leq\$")) M$(latexstring("\$_*\$")) )"
plot_hexbins_nxm( x, 
            y, 
            slcts, labels, 2, 2,
            #scale=0.0,
            cbarshift=2.5,
            outfile="/home/moon/sfortune/spinevo/plots/meeting_20220405/bvalue_vs_Nneighbors.png", 
            lognorm=true, gridres=(30,10), 
            ymin=miny, ymax=maxy,
            xmin=minx, xmax=maxx, 
            ylabel="b-value",
            xlabel="log( N$(latexstring("\$_{neighbors}\$")) )", fontsize=30,
            binned_median="x", plot_bval="x",
            )
























##############################################################################################################
# bval at z=0 vs N_flips (different angles) N_mergers (different masses)
##############################################################################################################
# sfc: endstate data outdated
bval_end            = load("/home/moon/sfortune/spinevo/data/endstate.jld", "bval_end")
nflips30            = load("/home/moon/sfortune/spinevo/data/endstate.jld", "nflips30")
nflips45            = load("/home/moon/sfortune/spinevo/data/endstate.jld", "nflips45")
nflips90            = load("/home/moon/sfortune/spinevo/data/endstate.jld", "nflips90")
nflips135           = load("/home/moon/sfortune/spinevo/data/endstate.jld", "nflips135")
G09_nflips30        = load("/home/moon/sfortune/spinevo/data/endstate.jld", "G09_nflips30")
G09_nflips45        = load("/home/moon/sfortune/spinevo/data/endstate.jld", "G09_nflips45")
G09_nflips90        = load("/home/moon/sfortune/spinevo/data/endstate.jld", "G09_nflips90")
G09_nflips135       = load("/home/moon/sfortune/spinevo/data/endstate.jld", "G09_nflips135")
nmergers2e8         = load("/home/moon/sfortune/spinevo/data/endstate.jld", "nmergers2e8")
nmergers2e9         = load("/home/moon/sfortune/spinevo/data/endstate.jld", "nmergers2e9")
nmergers2e10        = load("/home/moon/sfortune/spinevo/data/endstate.jld", "nmergers2e10")
G09_nmergers2e8     = load("/home/moon/sfortune/spinevo/data/endstate.jld", "G09_nmergers2e8 ")
G09_nmergers2e9     = load("/home/moon/sfortune/spinevo/data/endstate.jld", "G09_nmergers2e9 ")
G09_nmergers2e10    = load("/home/moon/sfortune/spinevo/data/endstate.jld", "G09_nmergers2e10")
# flips
y   = vcat( nflips30, nflips45, nflips90, nflips135, G09_nflips30, G09_nflips45, G09_nflips90, G09_nflips135 )
x   = vcat( bval_end, bval_end, bval_end, bval_end, bval_end, bval_end, bval_end, bval_end )
miny= minimum(pyplottable(y)[findall(x->x .!== NaN, pyplottable(y))])
maxy= maximum(pyplottable(y)[findall(x->x .!== NaN, pyplottable(y))])
minx= minimum(pyplottable(x)[findall(x->x .!== NaN, pyplottable(x))])
maxx= maximum(pyplottable(x)[findall(x->x .!== NaN, pyplottable(x))])
slcts       = Dict{Int64, Vector{Int64}}()
slcts[1]    = 1:length(bval_end)
slcts[2]    = length(bval_end)+1:2*length(bval_end)
slcts[3]    = 2*length(bval_end)+1:3*length(bval_end)
slcts[4]    = 3*length(bval_end)+1:4*length(bval_end)
slcts[5]    = 4*length(bval_end)+1:5*length(bval_end)
slcts[6]    = 5*length(bval_end)+1:6*length(bval_end)
slcts[7]    = 6*length(bval_end)+1:7*length(bval_end)
slcts[8]    = 7*length(bval_end)+1:8*length(bval_end)
labels      = Dict{Int64, String}()
labels[1]   = "$(latexstring("\$\\Delta\$"))t $(latexstring("\$\\sim\$")) 0.3 Gyr, N( $(latexstring("\$\\phi_{flip}\\ \\geq\\ 30° \$")) )"
labels[2]   = "N( $(latexstring("\$\\phi_{flip}\\ \\geq\\ 45° \$")) )"
labels[3]   = "N( $(latexstring("\$\\phi_{flip}\\ \\geq\\ 90° \$")) )"
labels[4]   = "N( $(latexstring("\$\\phi_{flip}\\ \\geq\\ 135°\$")) )"
labels[5]   = "$(latexstring("\$\\Delta\$"))t $(latexstring("\$\\sim  \$")) 1 Gyr, N( $(latexstring("\$\\phi_{flip}\\ \\geq\\ 30° \$")) )"
labels[6]   = "N( $(latexstring("\$\\phi_{flip}\\ \\geq\\ 45° \$")) )"
labels[7]   = "N( $(latexstring("\$\\phi_{flip}\\ \\geq\\ 90° \$")) )"
labels[8]   = "N( $(latexstring("\$\\phi_{flip}\\ \\geq\\ 135°\$")) )"
labpos      = Array{String}(undef, 8)
labpos     .= "upper right"
plot_hexbins_nxm( x, 
            y, 
            slcts, labels, 4, 2,
            #scale=0.0,
            cbarshift=2.5,
            outfile="/home/moon/sfortune/spinevo/plots/meeting_20220405/endstate_bval_vs_Nflips.png", 
            lognorm=true, gridres=(40,10), label_pos=labpos,
            ymin=miny, ymax=maxy,
            xmin=minx, xmax=maxx, 
            ylabel="Number of Flips",
            xlabel="b-value", fontsize=30,
            binned_median="x", plot_bval="x",
            )

# n_mergers
# sfc:fix
y   = vcat( nmergers2e8, nmergers2e9, nmergers2e10, G09_nmergers2e8, G09_nmergers2e9, G09_nmergers2e10 )
x   = vcat( bval_end, bval_end, bval_end, bval_end, bval_end, bval_end)
miny= minimum(pyplottable(y)[findall(x->x .!== NaN, pyplottable(y))])
maxy= maximum(pyplottable(y)[findall(x->x .!== NaN, pyplottable(y))])
minx= minimum(pyplottable(x)[findall(x->x .!== NaN, pyplottable(x))])
maxx= maximum(pyplottable(x)[findall(x->x .!== NaN, pyplottable(x))])
slcts       = Dict{Int64, Vector{Int64}}()
slcts[1]    = 1:length(bval_end)
slcts[2]    = length(bval_end)+1:2*length(bval_end)
slcts[3]    = 2*length(bval_end)+1:3*length(bval_end)
slcts[4]    = 3*length(bval_end)+1:4*length(bval_end)
slcts[5]    = 4*length(bval_end)+1:5*length(bval_end)
slcts[6]    = 5*length(bval_end)+1:6*length(bval_end)
labels      = Dict{Int64, String}()
labels[1]   = "$(latexstring("\$\\Delta\$"))t $(latexstring("\$\\sim\$")) 0.3 Gyr, N( M$(latexstring("\$_*\\geq\\ 2e8\\ M_\\odot\$")) )"
labels[2]   = "N( M$(latexstring("\$_*\\geq\\ 2e9\\ M_\\odot\$")) )"
labels[3]   = "N( M$(latexstring("\$_*\\geq\\ 2e10\\ M_\\odot\$")) )"
labels[4]   = "$(latexstring("\$\\Delta\$"))t $(latexstring("\$\\sim  \$")) 1 Gyr, N( M$(latexstring("\$_*\\geq\\ 2e8\\ M_\\odot\$")) )"
labels[5]   = " N( M$(latexstring("\$_*\\geq\\ 2e9\\ M_\\odot\$")) )"
labels[6]   = " N( M$(latexstring("\$_*\\geq\\ 2e10\\ M_\\odot\$")) )"
labpos      = Array{String}(undef, 6)
labpos     .= "upper right"
plot_hexbins_nxm( x, 
            y, 
            slcts, labels, 3, 2,
            #scale=0.0,
            cbarshift=2.5,
            outfile="/home/moon/sfortune/spinevo/plots/meeting_20220405/endstate_bval_vs_Nmergers.png", 
            lognorm=true, gridres=(30,10), label_pos=labpos,
            ymin=miny, ymax=maxy,
            xmin=minx, xmax=maxx, 
            ylabel="Number of Mergers",
            xlabel="b-value", fontsize=30,
            binned_median="x", plot_bval="x",
            )






##############################################################################################################
# M_star vs MfromJ
##############################################################################################################
plot_hexbins(log10.(pyplottable(as["M"])), log10.(pyplottable(as["M_fromJ"])); selection=start_thr, 
    outfile="/home/moon/sfortune/spinevo/plots/meeting_20220405/MfromJ_vs_Msubfind.png",
    xmin=" ", xmax=" ", ymin=" ", ymax=" ",
    scale=15, gridres=(30,13), xscale="linear", yscale="linear", 
    xlabel="log( M$(latexstring("\$_*\$")) [M$(latexstring("\$_\\odot\$"))] )", ylabel="log( M$(latexstring("\$_*\$")) [M$(latexstring("\$_\\odot\$"))] )", 
    grid=true,
    lognorm=true, colorlimits=(0.8,nothing),
    cmap="BuPu",
    binned_median = "x", plot_bval=" ",
    weights=" ",
    calc_norm=true,    # @pyplottable
    )













##############################################################################################################
##############################################################################################################
##############################################################################################################
##############################################################################################################
# from correlations


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
fig.savefig("/home/moon/sfortune/spinevo/plots/meeting_20220405/histogram_j_start.png", bbox_inches="tight", pad_inches=.1)
close("all")
PyPlot.matplotlib[:rc]("text", usetex=false)


##############################################################################################################
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
fig.savefig("/home/moon/sfortune/spinevo/plots/meeting_20220405/histogram_Jstart.png", bbox_inches="tight", pad_inches=.1)
close("all")
PyPlot.matplotlib[:rc]("text", usetex=false)




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
fig.savefig("/home/moon/sfortune/spinevo/plots/meeting_20220405/RELATIVE_delMdm_vs_delMstars_vs_flipbyMyr_09.png", bbox_inches="tight", pad_inches=.1)
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
fig.savefig("/home/moon/sfortune/spinevo/plots/meeting_20220405/delMdm_vs_delMstars_ABS_vs_flipangle_03.png", bbox_inches="tight", pad_inches=.1)
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
fig.savefig("/home/moon/sfortune/spinevo/plots/meeting_20220405/delMdm_vs_delMstars_ABS_vs_flipangle_03.png", bbox_inches="tight", pad_inches=.1)
close("all")
PyPlot.matplotlib[:rc]("text", usetex=false)


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
labels[2]   = "$(round(minimum(as["redshift"][slcts[3]]), digits=1))  $(latexstring("\$>\$")) z  $(latexstring("\$>\$")) $(round(maximum(as["redshift"][slcts[4]]), digits=1)), N = $(length(slcts[2]))"
labels[3]   = "$(round(maximum(as["redshift"][slcts[3]]), digits=1)) $(latexstring("\$\\geq\$")) z $(latexstring("\$\\geq\$")) $(round(minimum(as["redshift"][slcts[3]]), digits=1)), N = $(length(slcts[3]))"
labels[4]   = "$(round(maximum(as["redshift"][slcts[4]]), digits=1)) $(latexstring("\$\\geq\$")) z, N = $(length(slcts[4]))"
labpos      = Array{String}(undef, 4)
labpos     .= "upper left"
plot_hexbins_nxm( x, 
            y, 
            slcts, labels, 2, 2,
            #scale=0.0,
            cbarshift=2.5,
            label_pos=labpos, 
            outfile="/home/moon/sfortune/spinevo/plots/meeting_20220405/sSFR_vs_bval_zsplit.png", 
            lognorm=true, gridres=(15,5), 
            ymin=miny, ymax=-7.5,
            xmin=minx, xmax=maxx, 
            ylabel="log( sSFR [$(latexstring("\\frac{1}{yr}"))] )",
            xlabel="b-value", fontsize=30,
            binned_median="x", plot_bval="x",
            )
