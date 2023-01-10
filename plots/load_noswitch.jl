z_peak  = 1.9
jlim1   = 2.2
jlim2   = 2.55
Jlim1   = 12.6
Jlim2   = 12.9

asn03    = load("/home/moon/sfortune/spinevo/data/assembly_spinmap_noswitch_03Gyr.jld", "assembly_STARS")
println("Loaded asn03")
nsw03   = Dict{Symbol, Any}()
nsw03[:halo3212                 ] = findcs( asn03["ID_ISUB"          ], eq=3212 )
nsw03[:halo1414                 ] = findcs( asn03["ID_ISUB"          ], eq=1414 )
nsw03[:z_hi                     ] = findcs( asn03["snapNR"           ], leq=52 ) # borders: snap52=z1.18, snap100=z0.42
nsw03[:z_md                     ] = findcs( asn03["snapNR"           ], gt=52, lt=100 ) # borders: snap52=z1.18, snap100=z0.42
nsw03[:z_lo                     ] = findcs( asn03["snapNR"           ], geq=100 ) # borders: snap52=z1.18, snap100=z0.42
nsw03[:snap136                  ] = findcs( asn03["snapNR"           ], eq=136 )
nsw03[:snap36                   ] = findcs( asn03["snapNR"           ], eq=36 )
nsw03[:switches                 ] = findcs( asn03["switch"           ], eq=1 )
nsw03[:fakes                    ] = findcs( asn03["FAKEFLIP"         ], eq=1 )
nsw03[:result_disk              ] = findcs( asn03["BVAL"             ], geq=b_disk )
nsw03[:result_int               ] = findcs( asn03["BVAL"             ], gt=b_ell, lt=b_disk )
nsw03[:result_ell               ] = findcs( asn03["BVAL"             ], leq=b_ell )
nsw03[:z_welker                 ] = findcs( asn03["redshift"         ], geq=1.2, leq=3.8 )
nsw03[:z_lower                  ] = findcs( asn03["redshift"         ], lt=1.2 )
nsw03[:z_higher                 ] = findcs( asn03["redshift"         ], gt=3.8 )
nsw03[:z_prepeak                ] = findcs( asn03["redshift"         ], geq=z_peak )
nsw03[:z_postpeak               ] = findcs( asn03["redshift"         ], lt=z_peak )
nsw03[:z_hirsch_0               ] = findcs( asn03["redshift"         ], leq=0.25 )
nsw03[:z_hirsch_05              ] = findcs( asn03["redshift"         ], gt=0.25, lt=0.75 )
nsw03[:z_hirsch_1               ] = findcs( asn03["redshift"         ], geq=0.75, lt=1.5 )
nsw03[:z_hirsch_2               ] = findcs( asn03["redshift"         ], geq=1.5, lt=2.5 )
nsw03[:z_hirsch_3               ] = findcs( asn03["redshift"         ], geq=2.5, lt=3.5 )
nsw03[:z_hirsch_4               ] = findcs( asn03["redshift"         ], geq=3.5 )
nsw03[:mergerSTARS_acc_eq0      ] = findcs( asn03["Mpeak_MERGERS"    ], eq=0.0 )
nsw03[:mergerSTARS_acc_g10      ] = findcs( asn03["Mpeak_MERGERS"    ] ./ asn03["M2"      ], geq=0.1)
nsw03[:mergerSTARS_acc_g5_l10   ] = findcs( asn03["Mpeak_MERGERS"    ] ./ asn03["M2"      ], geq=0.05, lt=0.1)
nsw03[:mergerSTARS_acc_g0_l5    ] = findcs( asn03["Mpeak_MERGERS"    ] ./ asn03["M2"      ], gt=0.0, lt=0.05)
nsw03[:start_thr                ] = findcs( asn03["M"                ] .- asn03["ΔM"      ], geq=2e10)
nsw03[:start_disk               ] = findcs( asn03["BVAL_0"           ] .- asn03["ΔBVAL_0" ], geq=b_disk)
nsw03[:start_int                ] = findcs( asn03["BVAL_0"           ] .- asn03["ΔBVAL_0" ], gt=b_ell, lt=b_disk)
nsw03[:start_ell                ] = findcs( asn03["BVAL_0"           ] .- asn03["ΔBVAL_0" ], leq=b_ell)
nsw03[:j_hi                     ] = findcs( log10.(pyplottable( asn03["j_main"]  .- asn03["Δj_main"] )), gt=jlim2 )
nsw03[:j_md                     ] = findcs( log10.(pyplottable( asn03["j_main"]  .- asn03["Δj_main"] )), gt=jlim1 ,leq=jlim2 )
nsw03[:j_lo                     ] = findcs( log10.(pyplottable( asn03["j_main"]  .- asn03["Δj_main"] )), leq=jlim1 )
nsw03[:J_hi                     ] = findcs( log10.(pyplottable( asn03["J_main"]  .- asn03["ΔJ_main"] )), gt=Jlim2 )
nsw03[:J_md                     ] = findcs( log10.(pyplottable( asn03["J_main"]  .- asn03["ΔJ_main"] )), gt=Jlim1 ,leq=Jlim2 )
nsw03[:J_lo                     ] = findcs( log10.(pyplottable( asn03["J_main"]  .- asn03["ΔJ_main"] )), leq=Jlim1 )


# Dark matter 0.3 Gyr
#ad      = load("/home/moon/sfortune/spinevo/data/assembly_spinmap_noswitch_03Gyr.jld", "assembly_DM")



# 0.9Gyr data set
asn09    = load("/home/moon/sfortune/spinevo/data/assembly_spinmap_noswitch_09Gyr.jld", "assembly_STARS")
println("Loaded asn09")
nsw09 = Dict{Symbol, Any}()
nsw09[:halo3212                 ] = findcs( asn09["ID_ISUB"          ], eq=3212 )
nsw09[:halo1414                 ] = findcs( asn09["ID_ISUB"          ], eq=1414 )
nsw09[:z_hi                     ] = findcs( asn09["snapNR"           ], leq=52 ) # borders: snap52=z1.18, snap100=z0.42
nsw09[:z_md                     ] = findcs( asn09["snapNR"           ], gt=52, lt=100 ) # borders: snap52=z1.18, snap100=z0.42
nsw09[:z_lo                     ] = findcs( asn09["snapNR"           ], geq=100 ) # borders: snap52=z1.18, snap100=z0.42
nsw09[:snap136                  ] = findcs( asn09["snapNR"           ], eq=136 )
nsw09[:snap36                   ] = findcs( asn09["snapNR"           ], eq=36 )
nsw09[:switches                 ] = findcs( asn09["switch"           ], eq=1 )
nsw09[:fakes                    ] = findcs( asn09["FAKEFLIP"         ], eq=1 )
nsw09[:result_disk              ] = findcs( asn09["BVAL"             ], geq=b_disk )
nsw09[:result_int               ] = findcs( asn09["BVAL"             ], gt=b_ell, lt=b_disk )
nsw09[:result_ell               ] = findcs( asn09["BVAL"             ], leq=b_ell )
nsw09[:z_welker                 ] = findcs( asn09["redshift"         ], geq=1.2, leq=3.8 )
nsw09[:z_lower                  ] = findcs( asn09["redshift"         ], lt=1.2 )
nsw09[:z_higher                 ] = findcs( asn09["redshift"         ], gt=3.8 )
nsw09[:z_prepeak                ] = findcs( asn09["redshift"         ], geq=z_peak )
nsw09[:z_postpeak               ] = findcs( asn09["redshift"         ], lt=z_peak )
nsw09[:z_hirsch_0               ] = findcs( asn09["redshift"         ], leq=0.25 )
nsw09[:z_hirsch_05              ] = findcs( asn09["redshift"         ], gt=0.25, lt=0.75 )
nsw09[:z_hirsch_1               ] = findcs( asn09["redshift"         ], geq=0.75, lt=1.5 )
nsw09[:z_hirsch_2               ] = findcs( asn09["redshift"         ], geq=1.5, lt=2.5 )
nsw09[:z_hirsch_3               ] = findcs( asn09["redshift"         ], geq=2.5, lt=3.5 )
nsw09[:z_hirsch_4               ] = findcs( asn09["redshift"         ], geq=3.5 )
nsw09[:mergerSTARS_acc_eq0      ] = findcs( asn09["Mpeak_MERGERS"    ], eq=0.0 )
nsw09[:mergerSTARS_acc_g10      ] = findcs( asn09["Mpeak_MERGERS"    ] ./ asn09["M2"      ], geq=0.1)
nsw09[:mergerSTARS_acc_g5_l10   ] = findcs( asn09["Mpeak_MERGERS"    ] ./ asn09["M2"      ], geq=0.05, lt=0.1)
nsw09[:mergerSTARS_acc_g0_l5    ] = findcs( asn09["Mpeak_MERGERS"    ] ./ asn09["M2"      ], gt=0.0, lt=0.05)
nsw09[:start_thr                ] = findcs( asn09["M"                ] .- asn09["ΔM"      ], geq=2e10)
nsw09[:start_disk               ] = findcs( asn09["BVAL_0"           ] .- asn09["ΔBVAL_0" ], geq=b_disk)
nsw09[:start_int                ] = findcs( asn09["BVAL_0"           ] .- asn09["ΔBVAL_0" ], gt=b_ell, lt=b_disk)
nsw09[:start_ell                ] = findcs( asn09["BVAL_0"           ] .- asn09["ΔBVAL_0" ], leq=b_ell)
nsw09[:j_hi                     ] = findcs( log10.(pyplottable( asn09["j_main"]  .- asn09["Δj_main"] )), gt=jlim2 )
nsw09[:j_md                     ] = findcs( log10.(pyplottable( asn09["j_main"]  .- asn09["Δj_main"] )), gt=jlim1 ,leq=jlim2 )
nsw09[:j_lo                     ] = findcs( log10.(pyplottable( asn09["j_main"]  .- asn09["Δj_main"] )), leq=jlim1 )
nsw09[:J_hi                     ] = findcs( log10.(pyplottable( asn09["J_main"]  .- asn09["ΔJ_main"] )), gt=Jlim2 )
nsw09[:J_md                     ] = findcs( log10.(pyplottable( asn09["J_main"]  .- asn09["ΔJ_main"] )), gt=Jlim1 ,leq=Jlim2 )
nsw09[:J_lo                     ] = findcs( log10.(pyplottable( asn09["J_main"]  .- asn09["ΔJ_main"] )), leq=Jlim1 )

# Dark matter 0.3 Gyr
#ad09    = load("/home/moon/sfortune/spinevo/data/assembly_spinmap_noswitch_09Gyr.jld", "assembly_DM")








#nsw03[:stay_disk       = start_disk[  start_disk  .∈  Ref(Set(result_disk     ))]
#nsw03[:stay_int        = start_int[    start_int   .∈  Ref(Set(result_int      ))]
#nsw03[:stay_ell        = start_ell[    start_ell   .∈  Ref(Set(result_ell      ))]
#nsw03[:result_disk_fakes    = fakes[   fakes   .∈  Ref(Set(result_disk     ))]
#nsw03[:result_int_fakes     = fakes[   fakes   .∈  Ref(Set(result_int      ))]
#nsw03[:result_ell_fakes     = fakes[   fakes   .∈  Ref(Set(result_ell      ))]
#nsw03[:result_disk_nfakes   ]   = n_fakes[ n_fakes .∈  Ref(Set(result_disk     ))]
#nsw03[:result_int_nfakes    = n_fakes[ n_fakes .∈  Ref(Set(result_int      ))]
#nsw03[:result_ell_nfakes    = n_fakes[ n_fakes .∈  Ref(Set(result_ell      ))]
#nsw03[:start_disk_fakes    = fakes[   fakes   .∈  Ref(Set(start_disk     ))]
#nsw03[:start_int_fakes     = fakes[   fakes   .∈  Ref(Set(start_int      ))]
#nsw03[:start_ell_fakes     = fakes[   fakes   .∈  Ref(Set(start_ell      ))]
#nsw03[:start_disk_nfakes   = n_fakes[ n_fakes .∈  Ref(Set(start_disk     ))]
#nsw03[:start_int_nfakes    = n_fakes[ n_fakes .∈  Ref(Set(start_int      ))]
#nsw03[:start_ell_nfakes    = n_fakes[ n_fakes .∈  Ref(Set(start_ell      ))]
#nsw03[:stay_disk_fakes     = fakes[   fakes   .∈  Ref(Set(stay_disk     ))]
#nsw03[:stay_int_fakes      = fakes[   fakes   .∈  Ref(Set(stay_int      ))]
#nsw03[:stay_ell_fakes      = fakes[   fakes   .∈  Ref(Set(stay_ell      ))]
#nsw03[:stay_disk_nfakes    = n_fakes[ n_fakes .∈  Ref(Set(stay_disk     ))]
#nsw03[:stay_int_nfakes     = n_fakes[ n_fakes .∈  Ref(Set(stay_int      ))]
#nsw03[:stay_ell_nfakes     = n_fakes[ n_fakes .∈  Ref(Set(stay_ell      ))]
#nsw03[:ell_to_disk_nfakes  = start_ell_nfakes[start_ell_nfakes     .∈ Ref(Set(result_disk_nfakes))]
#nsw03[:ell_to_int_nfakes   = start_ell_nfakes[start_ell_nfakes     .∈ Ref(Set(result_int_nfakes))]
#nsw03[:int_to_disk_nfakes  = start_int_nfakes[start_int_nfakes     .∈ Ref(Set(result_disk_nfakes))]
#nsw03[:int_to_ell_nfakes   = start_int_nfakes[start_int_nfakes     .∈ Ref(Set(result_ell_nfakes))]
#nsw03[:disk_to_ell_nfakes  = start_disk_nfakes[start_disk_nfakes   .∈ Ref(Set(result_ell_nfakes))]
#nsw03[:disk_to_int_nfakes  = start_disk_nfakes[start_disk_nfakes   .∈ Ref(Set(result_int_nfakes))]
#nsw03[:mergersDM_1to50         = findcs(ad03["M2"]./ad03["Mpeak_MM"], geq=1, leq=50, comparewith=start_thr)











#G09_start_thr   = findcs(log10.(abs.(asn09["M"] .- asn09["ΔM"])), geq=log10.(2e10))
## borders: snap52=z1.18, snap100=z0.42
#G09_z_hi        = findcs(asn09["snapNR"      ], leq=52, comparewith=G09_start_thr)
#G09_z_md        = findcs(asn09["snapNR"      ], gt=52, lt=100, comparewith=G09_start_thr)
#G09_z_lo        = findcs(asn09["snapNR"      ], geq=100, comparewith=G09_start_thr)
#G09_snap136     = findcs(asn09["snapNR"      ], eq=136, comparewith=G09_start_thr)
#G09_snap36      = findcs(asn09["snapNR"      ], eq=36, comparewith=G09_start_thr)
#
#
#G09_halo3212    = findcs(asn09["ID_ISUB"      ], eq=3212, comparewith=G09_start_thr)
#G09_halo1414    = findcs(asn09["ID_ISUB"      ], eq=1414, comparewith=G09_start_thr)
#
#G09_switches    = findcs(asn09["switch"      ], eq=1, comparewith=G09_start_thr)
#G09_n_switches  = findcs(asn09["switch"      ], eq=0, comparewith=G09_start_thr)
#G09_fakes       = findcs(asn09["FAKEFLIP"    ], eq=1, comparewith=G09_start_thr)
#G09_n_fakes     = findcs(asn09["FAKEFLIP"    ], eq=0, comparewith=G09_start_thr)
#
#G09_result_disk     = findcs(asn09["BVAL"  ], geq=b_disk, comparewith=G09_start_thr)
#G09_result_int      = findcs(asn09["BVAL"  ], gt=b_ell, lt=b_disk, comparewith=G09_start_thr)
#G09_result_ell      = findcs(asn09["BVAL"  ], leq=b_ell, comparewith=G09_start_thr)
#
#G09_start_disk      = findcs(asn09["BVAL_0"  ] .- asn09["ΔBVAL_0"    ], geq=b_disk, comparewith=G09_start_thr)
#G09_start_int       = findcs(asn09["BVAL_0"  ] .- asn09["ΔBVAL_0"    ], gt=b_ell, lt=b_disk, comparewith=G09_start_thr)
#G09_start_ell       = findcs(asn09["BVAL_0"  ] .- asn09["ΔBVAL_0"    ], leq=b_ell, comparewith=G09_start_thr)
#G09_stay_disk       = G09_start_disk[  G09_start_disk  .∈  Ref(Set(G09_result_disk     ))]
#G09_stay_int        = G09_start_int[    G09_start_int   .∈  Ref(Set(G09_result_int      ))]
#G09_stay_ell        = G09_start_ell[    G09_start_ell   .∈  Ref(Set(G09_result_ell      ))]
#G09_result_disk_fakes    = G09_fakes[   G09_fakes   .∈  Ref(Set(G09_result_disk     ))]
#G09_result_int_fakes     = G09_fakes[   G09_fakes   .∈  Ref(Set(G09_result_int      ))]
#G09_result_ell_fakes     = G09_fakes[   G09_fakes   .∈  Ref(Set(G09_result_ell      ))]
#G09_result_disk_nfakes   = G09_n_fakes[ G09_n_fakes .∈  Ref(Set(G09_result_disk     ))]
#G09_result_int_nfakes    = G09_n_fakes[ G09_n_fakes .∈  Ref(Set(G09_result_int      ))]
#G09_result_ell_nfakes    = G09_n_fakes[ G09_n_fakes .∈  Ref(Set(G09_result_ell      ))]
#G09_start_disk_fakes    = G09_fakes[   G09_fakes   .∈  Ref(Set(G09_start_disk     ))]
#G09_start_int_fakes     = G09_fakes[   G09_fakes   .∈  Ref(Set(G09_start_int      ))]
#G09_start_ell_fakes     = G09_fakes[   G09_fakes   .∈  Ref(Set(G09_start_ell      ))]
#G09_start_disk_nfakes   = G09_n_fakes[ G09_n_fakes .∈  Ref(Set(G09_start_disk     ))]
#G09_start_int_nfakes    = G09_n_fakes[ G09_n_fakes .∈  Ref(Set(G09_start_int      ))]
#G09_start_ell_nfakes    = G09_n_fakes[ G09_n_fakes .∈  Ref(Set(G09_start_ell      ))]
#G09_stay_disk_fakes     = G09_fakes[   G09_fakes   .∈  Ref(Set(G09_stay_disk     ))]
#G09_stay_int_fakes      = G09_fakes[   G09_fakes   .∈  Ref(Set(G09_stay_int      ))]
#G09_stay_ell_fakes      = G09_fakes[   G09_fakes   .∈  Ref(Set(G09_stay_ell      ))]
#G09_stay_disk_nfakes    = G09_n_fakes[ G09_n_fakes .∈  Ref(Set(G09_stay_disk     ))]
#G09_stay_int_nfakes     = G09_n_fakes[ G09_n_fakes .∈  Ref(Set(G09_stay_int      ))]
#G09_stay_ell_nfakes     = G09_n_fakes[ G09_n_fakes .∈  Ref(Set(G09_stay_ell      ))]
#
#G09_ell_to_disk_nfakes  = G09_start_ell_nfakes[G09_start_ell_nfakes     .∈ Ref(Set(G09_result_disk_nfakes))]
#G09_ell_to_int_nfakes   = G09_start_ell_nfakes[G09_start_ell_nfakes     .∈ Ref(Set(G09_result_int_nfakes))]
#G09_int_to_disk_nfakes  = G09_start_int_nfakes[G09_start_int_nfakes     .∈ Ref(Set(G09_result_disk_nfakes))]
#G09_int_to_ell_nfakes   = G09_start_int_nfakes[G09_start_int_nfakes     .∈ Ref(Set(G09_result_ell_nfakes))]
#G09_disk_to_ell_nfakes  = G09_start_disk_nfakes[G09_start_disk_nfakes   .∈ Ref(Set(G09_result_ell_nfakes))]
#G09_disk_to_int_nfakes  = G09_start_disk_nfakes[G09_start_disk_nfakes   .∈ Ref(Set(G09_result_int_nfakes))]
#
#G09_mergersDM_1to50         = findcs(ad09["M2"]./ad09["Mpeak_MM"], geq=1, leq=50, comparewith=G09_start_thr)
#
#G09_mergerSTARS_acc_g10     = findcs(         asn09["Mpeak_MERGERS"] ./ asn09["M2"],       geq=0.1, comparewith=G09_start_thr)
#G09_mergerSTARS_acc_g5_l10  = findcs(         asn09["Mpeak_MERGERS"] ./ asn09["M2"],       geq=0.05, lt=0.1, comparewith=G09_start_thr)
#G09_mergerSTARS_acc_g0_l5   = findcs(         asn09["Mpeak_MERGERS"] ./ asn09["M2"],       gt=0.0, lt=0.05, comparewith=G09_start_thr)
#G09_mergerSTARS_acc_eq0     = findcs(replace( asn09["Mpeak_MERGERS"], missing => 0.0),  eq=0.0, comparewith=G09_start_thr)



#mergerDM_acc_g10     = findcs( ad03["Mpeak_MERGERS"] ./ ad03["M2"],       geq=0.1, comparewith=start_thr)
#mergerDM_acc_g5_l10  = findcs( ad03["Mpeak_MERGERS"] ./ ad03["M2"],       geq=0.05, lt=0.1, comparewith=start_thr)
#mergerDM_acc_g0_l5   = findcs( ad03["Mpeak_MERGERS"] ./ ad03["M2"],       gt=0.0, lt=0.05, comparewith=start_thr)
#mergerDM_acc_eq0     = findcs( ad03["Mpeak_MERGERS"], missing => 0.0),  eq=0.0, comparewith=start_thr)
#G09_mergerDM_acc_g10     = findcs(         ad09["Mpeak_MERGERS"] ./ ad09["M2"],       geq=0.1, comparewith=G09_start_thr)
#G09_mergerDM_acc_g5_l10  = findcs(         ad09["Mpeak_MERGERS"] ./ ad09["M2"],       geq=0.05, lt=0.1, comparewith=G09_start_thr)
#G09_mergerDM_acc_g0_l5   = findcs(         ad09["Mpeak_MERGERS"] ./ ad09["M2"],       gt=0.0, lt=0.05, comparewith=G09_start_thr)
#G09_mergerDM_acc_eq0     = findcs(replace( ad09["Mpeak_MERGERS"], missing => 0.0),  eq=0.0, comparewith=G09_start_thr)



#z_pre047   = findcs(asn03["redshift"    ], geq=0.47)
#z_post047  = findcs(asn03["redshift"    ], lt=0.47)
#z_pre064   = findcs(asn03["redshift"    ], geq=0.64)
#z_post064  = findcs(asn03["redshift"    ], lt=0.64)