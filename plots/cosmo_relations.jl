
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


# new

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


j_hi    = findcs(log10.(pyplottable(as["j_main"]  .- as["Δj_main"] )), gt=2.55, comparewith=start_thr)
j_md    = findcs(log10.(pyplottable(as["j_main"]  .- as["Δj_main"] )), gt=2.2 ,leq=2.55, comparewith=start_thr)
j_lo    = findcs(log10.(pyplottable(as["j_main"]  .- as["Δj_main"] )), leq=2.2, comparewith=start_thr)

J_hi    = findcs(log10.(pyplottable(as["J_main"]  .- as["ΔJ_main"] )), gt=12.9, comparewith=start_thr)
J_md    = findcs(log10.(pyplottable(as["J_main"]  .- as["ΔJ_main"] )), gt=12.6 ,leq=12.9, comparewith=start_thr)
J_lo    = findcs(log10.(pyplottable(as["J_main"]  .- as["ΔJ_main"] )), leq=12.6, comparewith=start_thr)


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






# SFR collection
sfr_coll    = Dict{String, Any}()
sfr_coll["redshift"]        = Array{Float64}(undef, 0)
sfr_coll["lookbacktime"]    = Array{Float64}(undef, 0)
sfr_coll["totSFR"]          = Array{Float64}(undef, 0)
sfr_coll["mean_sSFR"]       = Array{Float64}(undef, 0)
sfr_coll["N_halos"]         = Array{Int64}(undef, 0)
sfr_coll["totM"]            = Array{Float64}(undef, 0)
sfr_coll["boxsize"]         = Array{Float64}(undef, 0)
for i in sort(unique(as["snapNR"])) # snaps
    head    = read_header("$current_dir_simbox/groups_$(lpad(i, 3, "0"))/sub_$(lpad(i, 3, "0"))")
    idcs    = findcs(as["snapNR"], eq=i, comparewith=start_thr)
    sfr_coll["N_halos"]     = vcat(sfr_coll["N_halos"],     length(idcs) )
    sfr_coll["redshift"]    = vcat(sfr_coll["redshift"],    mean(as["redshift"  ][idcs]))
    sfr_coll["lookbacktime"]= vcat(sfr_coll["lookbacktime"],mean(as["lookbacktime"][idcs]))
    sfr_coll["boxsize"]     = vcat(sfr_coll["boxsize"],     prod(convert_units_physical([48000.0, 48000.0, 48000.0], :pos, head)) ) #kpc^3
    sfr_coll["totSFR"]      = vcat(sfr_coll["totSFR"],      sum( as["SFR"       ][idcs]))
    sfr_coll["mean_sSFR"]   = vcat(sfr_coll["mean_sSFR"],   mean(as["SFR"       ][idcs] ./ as["M"][idcs]))
    sfr_coll["totM"]        = vcat(sfr_coll["totM"],        sum( as["M"         ][idcs]))
end

# global SFR collection
global_sfr_coll    = Dict{String, Any}()
global_sfr_coll["N_halos"]         = Array{Int64}(undef, 0)
global_sfr_coll["redshift"]        = Array{Float64}(undef, 0)
global_sfr_coll["lookbacktime"]    = Array{Float64}(undef, 0)
global_sfr_coll["totSFR"]          = Array{Float64}(undef, 0)
global_sfr_coll["mean_sSFR"]       = Array{Float64}(undef, 0)
global_sfr_coll["N_halos"]         = Array{Int64}(undef, 0)
global_sfr_coll["totM"]            = Array{Float64}(undef, 0)
global_sfr_coll["boxsize"]         = Array{Float64}(undef, 0)
sub1e10_sfr_coll    = Dict{String, Any}()
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
    global_sfr_coll["N_halos"]     = vcat( global_sfr_coll["N_halos"],     length(ssfr) )
    global_sfr_coll["redshift"]    = vcat( global_sfr_coll["redshift"],    head.z )
    global_sfr_coll["lookbacktime"]= vcat( global_sfr_coll["lookbacktime"],ustrip(lookback_time(cosmology(h=head.h0, OmegaM=head.omega_0), head.z)) )
    global_sfr_coll["boxsize"]     = vcat( global_sfr_coll["boxsize"],     prod(convert_units_physical([48000.0, 48000.0, 48000.0], :pos, read_header("$current_dir_simbox/groups_$(lpad(i, 3, "0"))/sub_$(lpad(i, 3, "0"))"))) ) #kpc^3
    global_sfr_coll["totSFR"]      = vcat( global_sfr_coll["totSFR"],      sum(Float64.(ssfr)) )
    global_sfr_coll["totM"]        = vcat( global_sfr_coll["totM"],        Float64(sum(mass)) )
    global_sfr_coll["mean_sSFR"]   = vcat( global_sfr_coll["mean_sSFR"],   global_sfr_coll["totSFR"][end] / global_sfr_coll["totM"][end] )
    
    sub1e10_sfr_coll["N_halos"]     = vcat( sub1e10_sfr_coll["N_halos"],     length(ssfr[idx1e10]) )
    sub1e10_sfr_coll["redshift"]    = vcat( sub1e10_sfr_coll["redshift"],    head.z )
    sub1e10_sfr_coll["lookbacktime"]= vcat( sub1e10_sfr_coll["lookbacktime"],ustrip(lookback_time(cosmology(h=head.h0, OmegaM=head.omega_0), head.z)) )
    sub1e10_sfr_coll["boxsize"]     = vcat( sub1e10_sfr_coll["boxsize"],     prod(convert_units_physical([48000.0, 48000.0, 48000.0], :pos, read_header("$current_dir_simbox/groups_$(lpad(i, 3, "0"))/sub_$(lpad(i, 3, "0"))"))) ) #kpc^3
    sub1e10_sfr_coll["totSFR"]      = vcat( sub1e10_sfr_coll["totSFR"],      sum(Float64.(ssfr[idx1e10])) )
    sub1e10_sfr_coll["totM"]        = vcat( sub1e10_sfr_coll["totM"],        Float64(sum(convert_units_physical_mass(read_subfind("$current_dir_simbox/groups_$(lpad(i,3,"0"))/sub_$(lpad(i,3,"0"))", "SMST")[5,:], head)[idx1e10])) )
    sub1e10_sfr_coll["mean_sSFR"]   = vcat( sub1e10_sfr_coll["mean_sSFR"],   sub1e10_sfr_coll["totSFR"][end] / sub1e10_sfr_coll["totM"][end] )
end


function sfrd_madau(z)
    # M_⊙ * yr-1 * Mpc-3 * h3
    return 0.015*((1+z)^2.7)/(1+((1+z)/2.9)^5.6)
end

uhr = CSV.read("/HydroSims/Magneticum/Box4/uhr_test/sfr.txt", DataFrame; delim=' ', ignorerepeated=true, header=0)
hr  = CSV.read("/HydroSims/Magneticum/Box4/hr_kobayashi/sfr.txt", DataFrame; delim=' ', ignorerepeated=true, header=0)

z_madau     = LinRange(0, 12, 1000)
sfr_madau   = sfrd_madau.(z_madau)
z_uhr       = (1 ./ uhr[:,:Column1]) .- 1
sfr_uhr     = uhr[:,:Column3]
z_hr        = (1 ./ hr[:,:Column1]) .- 1
sfr_hr      = hr[:,:Column3]
idxmax_sfr_coll         = findcs(sfr_coll["totSFR"],  eq=maximum(sfr_coll["totSFR"]))[end]
idxmax_global_sfr_coll  = findcs(global_sfr_coll["totSFR"],   eq=maximum(global_sfr_coll["totSFR"]))[end]
idxmax_sub1e10_sfr_coll  = findcs(sub1e10_sfr_coll["totSFR"],   eq=maximum(sub1e10_sfr_coll["totSFR"]))[end]
idxmax_madau            = findcs(sfr_madau,  eq=maximum(sfr_madau))[end]
idxmax_uhr              = findcs(sfr_uhr,  eq=maximum(sfr_uhr))[end]
idxmax_hr               = findcs(sfr_hr,  eq=maximum(sfr_hr))[end]











# SFR / boxsize
#, fontsize=floor(0.5*scale)
#./ sfr_coll["boxsize"] ./ 1e9
PyPlot.matplotlib[:rc]("text", usetex=true)
PyPlot.matplotlib[:rc]("font", family="serif", serif=["Computer Modern"])
pointsize = 5
scale = 0.6
fig, ax = subplots()
ax.axvline(z_madau[idxmax_madau], color="black", label="Madau, z = $(round(z_madau[idxmax_madau], digits=1))")
ax.axvline(z_hr[idxmax_hr], color="darkred", label="Box4hr, z = $(round(z_hr[idxmax_hr], digits=1))")
ax.axvline(z_uhr[idxmax_uhr], color="red", label="Box4uhr, z = $(round(z_uhr[idxmax_uhr], digits=1))")
ax.axvline(global_sfr_coll["redshift"][idxmax_global_sfr_coll], color="mediumblue", label="All subhalos, z = $(round(global_sfr_coll["redshift"][idxmax_global_sfr_coll], digits=1))")
ax.axvline(sub1e10_sfr_coll["redshift"][idxmax_sub1e10_sfr_coll], color="dodgerblue", label="All 1e10-cut, z = $(round(sub1e10_sfr_coll["redshift"][idxmax_sub1e10_sfr_coll], digits=1))")
ax.axvline(sfr_coll["redshift"][idxmax_sfr_coll], color="cyan", label="Centrals 1e10-cut, z = $(round(sfr_coll["redshift"][idxmax_sfr_coll], digits=1))")
#ax.plot( global_sfr_coll["redshift"], log10.(global_sfr_coll["totSFR"] ./ global_sfr_coll["boxsize"] .* 1e9), "-", color="darkred", markersize=5, label="Global")
ax.plot( z_madau, log10.(sfr_madau), "-", color="black")
ax.plot( z_hr, log10.(sfr_hr ./ 48^3), "-", color="darkred")
ax.plot( z_uhr, log10.(sfr_uhr ./ 48^3), "-", color="red")
ax.plot( global_sfr_coll["redshift"], log10.(global_sfr_coll["totSFR"] ./ 48^3), ".-", markersize=2, color="mediumblue")
ax.plot( sub1e10_sfr_coll["redshift"], log10.(sub1e10_sfr_coll["totSFR"] ./ 48^3), ".-", markersize=2, color="dodgerblue")
ax.plot( sfr_coll["redshift"], log10.(sfr_coll["totSFR"] ./ 48^3), ".-", markersize=2, color="cyan")
#ax.set_title("$(latexstring("\$\\Delta\$"))t $(latexstring("\$\\sim\$")) 0.3 Gyr")
ax.set_xlim([0,8.5])
ax.set_ylim([-2,nothing])
ax.set_xlabel("redshift")
ax.set_xscale("linear")
ax.set_xticks([1,2,3,4,5,6,7,8])
ax.set_xticklabels([1,2,3,4,5,6,7,8])
ax.set_ylabel("log( SFRD [$(latexstring("\\frac{M_\\odot}{yr\\ (Mpc/h)^3}"))] )")
ax.legend(loc="upper right", frameon=true, borderpad=1, handlelength=1.8)
ax.grid()
fig.set_size_inches(16scale, 9scale)
fig.tight_layout()
fig.savefig("/home/moon/sfortune/spinevo/plots/cosmo_relations/SFRDcomoving_vs_redshift.png", bbox_inches="tight", pad_inches=.1)
close("all")
PyPlot.matplotlib[:rc]("text", usetex=false)


# specific Star formation rate vs redshift
PyPlot.matplotlib[:rc]("text", usetex=true)
PyPlot.matplotlib[:rc]("font", family="serif", serif=["Computer Modern"])
pointsize = 10
scale = 0.7
fig, ax = subplots()
#ax.axvline(b_disk, color="blue")
ax.axvline(1.9, color="black", label="z = 1.9")
ax.plot( global_sfr_coll["redshift"], global_sfr_coll["mean_sSFR"] .* 1e9, "-", color="darkred", ms=pointsize, label="Global")
ax.plot( sfr_coll["redshift"], sfr_coll["mean_sSFR"] .* 1e9, "-", color="mediumblue", ms=pointsize, label="Centrals")
#ax.set_title("$(latexstring("\$\\Delta\$"))t $(latexstring("\$\\sim\$")) 0.3 Gyr")
#ax.set_xlim([maximum])
ax.set_yscale("log")
ax.set_xlabel("z")
ax.set_ylabel("sSFR [$(latexstring("\\frac{1}{Gyr}"))]")
ax.legend(loc="upper left", frameon=true, borderpad=1, handlelength=1.8)
ax.grid()
fig.set_size_inches(9scale, 9scale)
fig.tight_layout()
fig.savefig("/home/moon/sfortune/spinevo/plots/cosmo_relations/sSFR_vs_redshift.png", bbox_inches="tight", pad_inches=.1)
close("all")
PyPlot.matplotlib[:rc]("text", usetex=false)


# SFR / n_halos / boxsize
PyPlot.matplotlib[:rc]("text", usetex=true)
PyPlot.matplotlib[:rc]("font", family="serif", serif=["Computer Modern"])
pointsize = 5
scale = 0.7
fig, ax = subplots()
#ax.axvline(b_disk, color="blue")
ax.axvline(1.9, color="black", label="z = 1.9")
ax.plot( global_sfr_coll["redshift"], log10.(global_sfr_coll["totSFR"] ./ global_sfr_coll["N_halos"] ./ global_sfr_coll["boxsize"] .* 1e9), "-", color="darkred", markersize=5, label="Global")
ax.plot( sfr_coll["redshift"], log10.(sfr_coll["totSFR"] ./ sfr_coll["N_halos"] ./ sfr_coll["boxsize"] .* 1e9), "-", color="mediumblue", markersize=5, label="Centrals")
#ax.set_title("$(latexstring("\$\\Delta\$"))t $(latexstring("\$\\sim\$")) 0.3 Gyr")
#ax.set_xlim([maximum])
ax.set_xlabel("z")
ax.set_xscale("linear")
ax.set_ylabel("log( SFRD [$(latexstring("\\frac{M_\\odot}{yr\\ Mpc^3\\ N_{halos}}"))] )")
ax.legend(loc="upper left", frameon=true, borderpad=1, handlelength=1.8)
ax.grid()
fig.set_size_inches(9scale, 9scale)
fig.tight_layout()
fig.savefig("/home/moon/sfortune/spinevo/plots/cosmo_relations/SFRDNHphysical_vs_redshift.png", bbox_inches="tight", pad_inches=.1)
close("all")
PyPlot.matplotlib[:rc]("text", usetex=false)







# z_hist
PyPlot.matplotlib[:rc]("text", usetex=true)
PyPlot.matplotlib[:rc]("font", family="serif", serif=["Computer Modern"])
scale = 0.5
fig, ax = subplots()
ax.hist( as["redshift"][start_thr], bins=length(unique(as["redshift"][start_thr]))    , rwidth=0.9    , alpha=1, color="mediumblue" )
ax.axvline( 1.9, color ="black", label= "z = 1.9")
ax.axvline( median(as["redshift"][start_thr]), color ="darkred", label= "z = $(round(median(as["redshift"][start_thr]), digits=2))")
ax.axvline( mean(as["redshift"][start_thr]), color ="darkorchid", label= "z = $(round(mean(as["redshift"][start_thr]), digits=2))")
ax.set_xlabel("z")
ax.set_xscale("linear")
#ax.set_xlim([10, nothing])
ax.legend(loc="upper left", frameon=true, borderpad=1, handlelength=1.8)
ax.grid()
scale=0.5
fig.set_size_inches(16scale, 9scale)
fig.tight_layout()
fig.savefig("/home/moon/sfortune/spinevo/plots/cosmo_relations/histogram_redshift.png", bbox_inches="tight", pad_inches=.1)
close("all")
PyPlot.matplotlib[:rc]("text", usetex=false)















# ϕ_flip / ΔJ
PyPlot.matplotlib[:rc]("text", usetex=true)
PyPlot.matplotlib[:rc]("font", family="serif", serif=["Computer Modern"])
x1  = log10.(replace(pyplottable(  as[   "ϕ_flip"][start_thr] ./ log10.(pyplottable(as[     "ΔJ_main"][:,start_thr]))), 0.0 => NaN))
x2  = log10.(replace(pyplottable(  as09[ "ϕ_flip"][start_thr] ./ log10.(pyplottable(as09[   "ΔJ_main"][:,start_thr]))), 0.0 => NaN))
scale = 0.5
fig, ax = subplots()
ax.hist( x2, bins=10    , rwidth=0.9    , alpha=0.7, color="darkred", label="1 Gyr" )
ax.hist( x1, bins=10    , rwidth=0.9    , alpha=0.7, color="mediumblue", label="0.3 Gyr" )
#ax.axvline( 1.9, color ="black", label= "z = 1.9")
#ax.axvline( median(as["redshift"][start_thr]), color ="darkred", label= "z = $(round(median(as["redshift"][start_thr]), digits=2))")
#ax.axvline( mean(as["redshift"][start_thr]), color ="darkorchid", label= "z = $(round(mean(as["redshift"][start_thr]), digits=2))")
ax.set_xlabel("flip by delta J")
ax.set_xscale("linear")
ax.set_yscale("linear")
#ax.set_xlim([10, nothing])
ax.legend(loc="upper left", frameon=true, borderpad=1, handlelength=1.8)
ax.grid()
scale=0.5
fig.set_size_inches(16scale, 9scale)
fig.tight_layout()
fig.savefig("/home/moon/sfortune/spinevo/plots/cosmo_relations/histogram_flipangleOVERdelJ.png", bbox_inches="tight", pad_inches=.1)
close("all")
PyPlot.matplotlib[:rc]("text", usetex=false)
