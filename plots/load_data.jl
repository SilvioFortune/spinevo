z_peak  = 1.9
jlim1   = 2.2
jlim2   = 2.55
Jlim1   = 12.6
Jlim2   = 12.9
z_split15 = 1.5
as03    = load("/home/moon/sfortune/spinevo/data/assembly_spinmap_03Gyr.jld", "assembly_STARS")
median_logSFR  = median(pyplottable(log10.(pyplottable(1e9 .* as03["SFR"])))[findall(x->x .!== NaN, pyplottable(log10.(pyplottable(1e9 .* as03["SFR"]))))])
median_logsSFR = median(pyplottable(pyplottable(log10.(pyplottable(1e9 .* as03["SFR"] ./ as03["M"]))))[findall(x->x .!== NaN, pyplottable(pyplottable(log10.(pyplottable(1e9 .* as03["SFR"] ./ as03["M"])))))])
pctl68up_logSFR  = pctl68upper(pyplottable(log10.(pyplottable(1e9 .* as03["SFR"])))[findall(x->x .!== NaN, pyplottable(log10.(pyplottable(1e9 .* as03["SFR"]))))])
pctl68up_logsSFR = pctl68upper(pyplottable(pyplottable(log10.(pyplottable(1e9 .* as03["SFR"] ./ as03["M"]))))[findall(x->x .!== NaN, pyplottable(pyplottable(log10.(pyplottable(1e9 .* as03["SFR"] ./ as03["M"])))))])
pctl68lo_logSFR  = pctl68lower(pyplottable(log10.(pyplottable(1e9 .* as03["SFR"])))[findall(x->x .!== NaN, pyplottable(log10.(pyplottable(1e9 .* as03["SFR"]))))])
pctl68lo_logsSFR = pctl68lower(pyplottable(pyplottable(log10.(pyplottable(1e9 .* as03["SFR"] ./ as03["M"]))))[findall(x->x .!== NaN, pyplottable(pyplottable(log10.(pyplottable(1e9 .* as03["SFR"] ./ as03["M"])))))])


println("Loaded as03")

# Dark matter 0.3 Gyr
ad03    = load("/home/moon/sfortune/spinevo/data/assembly_spinmap_03Gyr.jld", "assembly_DM")
println("Loaded ad03")

csw03   = Dict{Symbol, Any}()
csw03[:SFR_lo                   ] = findcs( log10.(1e9 .* as03["SFR"]), leq=pctl68lo_logSFR )
csw03[:SFR_md                   ] = findcs( log10.(1e9 .* as03["SFR"]), gt=pctl68lo_logSFR, lt=pctl68up_logSFR )
csw03[:SFR_hi                   ] = findcs( log10.(1e9 .* as03["SFR"]), geq=pctl68up_logSFR )
csw03[:sSFR_lo                  ] = findcs( log10.(1e9 .* as03["SFR"] ./ as03["M"]), leq=pctl68lo_logsSFR )
csw03[:sSFR_md                  ] = findcs( log10.(1e9 .* as03["SFR"] ./ as03["M"]), gt=pctl68lo_logsSFR, lt=pctl68up_logsSFR )
csw03[:sSFR_hi                  ] = findcs( log10.(1e9 .* as03["SFR"] ./ as03["M"]), geq=pctl68up_logsSFR )
csw03[:z_split15_l              ] = findcs( as03["redshift"         ], lt=z_split15 )
csw03[:z_split15_g              ] = findcs( as03["redshift"         ], geq=z_split15 )
csw03[:halo3212                 ] = findcs( as03["ID_ISUB"          ], eq=3212 )
csw03[:halo1414                 ] = findcs( as03["ID_ISUB"          ], eq=1414 )
csw03[:z_0                      ] = findcs( as03["redshift"         ], eq=minimum(as03["redshift"]) )
csw03[:z_2                      ] = findcs( as03["redshift"         ], geq=1.5, leq=2.5 )
csw03[:z_hi                     ] = findcs( as03["snapNR"           ], leq=52 ) # borders: snap52=z1.18, snap100=z0.42
csw03[:z_md                     ] = findcs( as03["snapNR"           ], gt=52, lt=100 ) # borders: snap52=z1.18, snap100=z0.42
csw03[:z_lo                     ] = findcs( as03["snapNR"           ], geq=100 ) # borders: snap52=z1.18, snap100=z0.42
csw03[:z_hi2                    ] = findcs( as03["redshift"         ], geq=1.5 )
csw03[:z_md2                    ] = findcs( as03["redshift"         ], lt=1.5, gt=0.5 )
csw03[:z_lo2                    ] = findcs( as03["redshift"         ], leq=0.5 )
csw03[:snap136                  ] = findcs( as03["snapNR"           ], eq=136 )
csw03[:snap36                   ] = findcs( as03["snapNR"           ], eq=36 )
csw03[:switches                 ] = findcs( as03["switch"           ], gt=0 )
csw03[:fakes                    ] = findcs( as03["FAKEFLIP"         ], eq=1 )
csw03[:result_disk              ] = findcs( as03["BVAL"             ], geq=b_disk )
csw03[:result_int               ] = findcs( as03["BVAL"             ], gt=b_ell, lt=b_disk )
csw03[:result_ell               ] = findcs( as03["BVAL"             ], leq=b_ell )
csw03[:z_welker                 ] = findcs( as03["redshift"         ], geq=1.2, leq=3.8 )
csw03[:z_lower                  ] = findcs( as03["redshift"         ], lt=1.2 )
csw03[:z_higher                 ] = findcs( as03["redshift"         ], gt=3.8 )
csw03[:z_prepeak                ] = findcs( as03["redshift"         ], geq=z_peak )
csw03[:z_postpeak               ] = findcs( as03["redshift"         ], lt=z_peak )
csw03[:z_hirsch_0               ] = findcs( as03["redshift"         ], leq=0.25 )
csw03[:z_hirsch_05              ] = findcs( as03["redshift"         ], gt=0.25, lt=0.75 )
csw03[:z_hirsch_1               ] = findcs( as03["redshift"         ], geq=0.75, lt=1.5 )
csw03[:z_hirsch_2               ] = findcs( as03["redshift"         ], geq=1.5, lt=2.5 )
csw03[:z_hirsch_3               ] = findcs( as03["redshift"         ], geq=2.5, lt=3.5 )
csw03[:z_hirsch_4               ] = findcs( as03["redshift"         ], geq=3.5 )
csw03[:mergerSTARS_acc_eq0      ] = findcs( as03["Mpeak_MERGERS"    ], eq=0.0 )
csw03[:mergerSTARS_acc_g10      ] = findcs( as03["Mpeak_MERGERS"    ] ./ as03["M2"      ], geq=0.1)
csw03[:mergerSTARS_acc_g5_l10   ] = findcs( as03["Mpeak_MERGERS"    ] ./ as03["M2"      ], geq=0.05, lt=0.1)
csw03[:mergerSTARS_acc_g0_l5    ] = findcs( as03["Mpeak_MERGERS"    ] ./ as03["M2"      ], gt=0.0, lt=0.05)
csw03[:mergerDM_acc_eq0         ] = findcs( ad03["Mpeak_MERGERS"    ], eq=0.0 )
csw03[:mergerDM_acc_g10         ] = findcs( ad03["Mpeak_MERGERS"    ] ./ ad03["M2"      ], geq=0.1)
csw03[:mergerDM_acc_g5_l10      ] = findcs( ad03["Mpeak_MERGERS"    ] ./ ad03["M2"      ], geq=0.05, lt=0.1)
csw03[:mergerDM_acc_g0_l5       ] = findcs( ad03["Mpeak_MERGERS"    ] ./ ad03["M2"      ], gt=0.0, lt=0.05)
csw03[:start_thr                ] = findcs( as03["M"                ] .- as03["ΔM"      ], geq=2e10)
csw03[:start_disk               ] = findcs( as03["BVAL_0"           ] .- as03["ΔBVAL_0" ], geq=b_disk)
csw03[:start_int                ] = findcs( as03["BVAL_0"           ] .- as03["ΔBVAL_0" ], gt=b_ell, lt=b_disk)
csw03[:start_ell                ] = findcs( as03["BVAL_0"           ] .- as03["ΔBVAL_0" ], leq=b_ell)
csw03[:j_hi                     ] = findcs( log10.(pyplottable( as03["j_main"]  .- as03["Δj_main"] )), gt=jlim2 )
csw03[:j_md                     ] = findcs( log10.(pyplottable( as03["j_main"]  .- as03["Δj_main"] )), gt=jlim1 ,leq=jlim2 )
csw03[:j_lo                     ] = findcs( log10.(pyplottable( as03["j_main"]  .- as03["Δj_main"] )), leq=jlim1 )
csw03[:J_hi                     ] = findcs( log10.(pyplottable( as03["J_main"]  .- as03["ΔJ_main"] )), gt=Jlim2 )
csw03[:J_md                     ] = findcs( log10.(pyplottable( as03["J_main"]  .- as03["ΔJ_main"] )), gt=Jlim1 ,leq=Jlim2 )
csw03[:J_lo                     ] = findcs( log10.(pyplottable( as03["J_main"]  .- as03["ΔJ_main"] )), leq=Jlim1 )

csw03[:Nmerger_0   ] = findcs( as03["N_MERGERS"], eq=0)
csw03[:Nmerger_1_2 ] = findcs( as03["N_MERGERS"], gt=0,leq=2)
csw03[:Nmerger_3_10] = findcs( as03["N_MERGERS"], gt=2,leq=10)
csw03[:Nmerger_11_ ] = findcs( as03["N_MERGERS"], gt=10)

csw03[:transition               ] = idxclude(idxclude(csw03[:start_thr], csw03[:switches]), findcs(norm.(as03["Δlookbacktime"]), eq=0.0) )

a = as03["M"][csw03[:transition]] .- as03["ΔM"][csw03[:transition]]
median_logMstarpre  = median(pyplottable(log10.(pyplottable(a)))[findall(x->x .!== NaN, pyplottable(log10.(pyplottable(a))))])
pctl68up_logMstarpre  = pctl68upper(pyplottable(log10.(pyplottable(a)))[findall(x->x .!== NaN, pyplottable(log10.(pyplottable(a))))])
pctl68lo_logMstarpre  = pctl68lower(pyplottable(log10.(pyplottable(a)))[findall(x->x .!== NaN, pyplottable(log10.(pyplottable(a))))])
csw03[:Mpre_lo                  ] = idxmatch(csw03[:transition],findcs( log10.(as03["M"] .- as03["ΔM"]), leq=pctl68lo_logMstarpre ))
csw03[:Mpre_md                  ] = idxmatch(csw03[:transition],findcs( log10.(as03["M"] .- as03["ΔM"]), gt=pctl68lo_logMstarpre, lt=pctl68up_logMstarpre ))
csw03[:Mpre_hi                  ] = idxmatch(csw03[:transition],findcs( log10.(as03["M"] .- as03["ΔM"]), geq=pctl68up_logMstarpre ))
# sfc: change j limits and decide on naming before and after flip




# 0.9Gyr data set
as09    = load("/home/moon/sfortune/spinevo/data/assembly_spinmap_09Gyr.jld", "assembly_STARS")
println("Loaded as09")

# Dark matter 1 Gyr
ad09    = load("/home/moon/sfortune/spinevo/data/assembly_spinmap_09Gyr_patch.jld", "assembly_DM")
println("Loaded ad09")


csw09 = Dict{Symbol, Any}()
csw09[:SFR_lo                   ] = findcs( log10.(1e9 .* as09["SFR"]), leq=pctl68lo_logSFR )
csw09[:SFR_md                   ] = findcs( log10.(1e9 .* as09["SFR"]), gt=pctl68lo_logSFR, lt=pctl68up_logSFR )
csw09[:SFR_hi                   ] = findcs( log10.(1e9 .* as09["SFR"]), geq=pctl68up_logSFR )
csw09[:sSFR_lo                  ] = findcs( log10.(1e9 .* as09["SFR"] ./ as09["M"]), leq=pctl68lo_logsSFR )
csw09[:sSFR_md                  ] = findcs( log10.(1e9 .* as09["SFR"] ./ as09["M"]), gt=pctl68lo_logsSFR, lt=pctl68up_logsSFR )
csw09[:sSFR_hi                  ] = findcs( log10.(1e9 .* as09["SFR"] ./ as09["M"]), geq=pctl68up_logsSFR )
csw09[:z_split15_l              ] = findcs( as09["redshift"         ], lt=z_split15 )
csw09[:z_split15_g              ] = findcs( as09["redshift"         ], geq=z_split15 )
csw09[:halo3212                 ] = findcs( as09["ID_ISUB"          ], eq=3212 )
csw09[:halo1414                 ] = findcs( as09["ID_ISUB"          ], eq=1414 )
csw09[:z_0                     ]  = findcs( as09["redshift"         ], eq=minimum(as09["redshift"]) )
csw03[:z_2                     ]  = findcs( as09["redshift"         ], geq=1.5, leq=2.5 )
csw09[:z_hi                     ] = findcs( as09["snapNR"           ], leq=52 ) # borders: snap52=z1.18, snap100=z0.42
csw09[:z_md                     ] = findcs( as09["snapNR"           ], gt=52, lt=100 ) # borders: snap52=z1.18, snap100=z0.42
csw09[:z_lo                     ] = findcs( as09["snapNR"           ], geq=100 ) # borders: snap52=z1.18, snap100=z0.42
csw09[:snap136                  ] = findcs( as09["snapNR"           ], eq=136 )
csw09[:snap36                   ] = findcs( as09["snapNR"           ], eq=36 )
csw09[:switches                 ] = findcs( as09["switch"           ], gt=0 )
csw09[:fakes                    ] = findcs( as09["FAKEFLIP"         ], eq=1 )
csw09[:result_disk              ] = findcs( as09["BVAL"             ], geq=b_disk )
csw09[:result_int               ] = findcs( as09["BVAL"             ], gt=b_ell, lt=b_disk )
csw09[:result_ell               ] = findcs( as09["BVAL"             ], leq=b_ell )
csw09[:z_welker                 ] = findcs( as09["redshift"         ], geq=1.2, leq=3.8 )
csw09[:z_lower                  ] = findcs( as09["redshift"         ], lt=1.2 )
csw09[:z_higher                 ] = findcs( as09["redshift"         ], gt=3.8 )
csw09[:z_prepeak                ] = findcs( as09["redshift"         ], geq=z_peak )
csw09[:z_postpeak               ] = findcs( as09["redshift"         ], lt=z_peak )
csw09[:z_hirsch_0               ] = findcs( as09["redshift"         ], leq=0.25 )
csw09[:z_hirsch_05              ] = findcs( as09["redshift"         ], gt=0.25, lt=0.75 )
csw09[:z_hirsch_1               ] = findcs( as09["redshift"         ], geq=0.75, lt=1.5 )
csw09[:z_hirsch_2               ] = findcs( as09["redshift"         ], geq=1.5, lt=2.5 )
csw09[:z_hirsch_3               ] = findcs( as09["redshift"         ], geq=2.5, lt=3.5 )
csw09[:z_hirsch_4               ] = findcs( as09["redshift"         ], geq=3.5 )
csw09[:mergerSTARS_acc_eq0      ] = findcs( as09["Mpeak_MERGERS"    ], eq=0.0 )
csw09[:mergerSTARS_acc_g10      ] = findcs( as09["Mpeak_MERGERS"    ] ./ as09["M2"      ], geq=0.1)
csw09[:mergerSTARS_acc_g5_l10   ] = findcs( as09["Mpeak_MERGERS"    ] ./ as09["M2"      ], geq=0.05, lt=0.1)
csw09[:mergerSTARS_acc_g0_l5    ] = findcs( as09["Mpeak_MERGERS"    ] ./ as09["M2"      ], gt=0.0, lt=0.05)
csw09[:mergerDM_acc_eq0         ] = findcs( ad09["Mpeak_MERGERS"    ], eq=0.0 )
csw09[:mergerDM_acc_g10         ] = findcs( ad09["Mpeak_MERGERS"    ] ./ ad09["M2"      ], geq=0.1)
csw09[:mergerDM_acc_g5_l10      ] = findcs( ad09["Mpeak_MERGERS"    ] ./ ad09["M2"      ], geq=0.05, lt=0.1)
csw09[:mergerDM_acc_g0_l5       ] = findcs( ad09["Mpeak_MERGERS"    ] ./ ad09["M2"      ], gt=0.0, lt=0.05)
csw09[:start_thr                ] = findcs( as09["M"                ] .- as09["ΔM"      ], geq=2e10)
csw09[:start_disk               ] = findcs( as09["BVAL_0"           ] .- as09["ΔBVAL_0" ], geq=b_disk)
csw09[:start_int                ] = findcs( as09["BVAL_0"           ] .- as09["ΔBVAL_0" ], gt=b_ell, lt=b_disk)
csw09[:start_ell                ] = findcs( as09["BVAL_0"           ] .- as09["ΔBVAL_0" ], leq=b_ell)
csw09[:j_hi                     ] = findcs( log10.(pyplottable( as09["j_main"]  .- as09["Δj_main"] )), gt=jlim2 )
csw09[:j_md                     ] = findcs( log10.(pyplottable( as09["j_main"]  .- as09["Δj_main"] )), gt=jlim1 ,leq=jlim2 )
csw09[:j_lo                     ] = findcs( log10.(pyplottable( as09["j_main"]  .- as09["Δj_main"] )), leq=jlim1 )
csw09[:J_hi                     ] = findcs( log10.(pyplottable( as09["J_main"]  .- as09["ΔJ_main"] )), gt=Jlim2 )
csw09[:J_md                     ] = findcs( log10.(pyplottable( as09["J_main"]  .- as09["ΔJ_main"] )), gt=Jlim1 ,leq=Jlim2 )
csw09[:J_lo                     ] = findcs( log10.(pyplottable( as09["J_main"]  .- as09["ΔJ_main"] )), leq=Jlim1 )

csw09[:Nmerger_0   ] = findcs( as09["N_MERGERS"], eq=0)
csw09[:Nmerger_1_2 ] = findcs( as09["N_MERGERS"], gt=0,leq=2)
csw09[:Nmerger_3_10] = findcs( as09["N_MERGERS"], gt=2,leq=10)
csw09[:Nmerger_11_ ] = findcs( as09["N_MERGERS"], gt=10)


csw09[:transition               ] = idxclude(idxclude(csw09[:start_thr], csw09[:switches]), findcs(norm.(as09["Δlookbacktime"]), eq=0.0))

a = as09["M"][csw09[:transition]] .- as09["ΔM"][csw09[:transition]]
median_logMstarpre09  = median(pyplottable(log10.(pyplottable(a)))[findall(x->x .!== NaN, pyplottable(log10.(pyplottable(a))))])
pctl68up_logMstarpre09  = pctl68upper(pyplottable(log10.(pyplottable(a)))[findall(x->x .!== NaN, pyplottable(log10.(pyplottable(a))))])
pctl68lo_logMstarpre09  = pctl68lower(pyplottable(log10.(pyplottable(a)))[findall(x->x .!== NaN, pyplottable(log10.(pyplottable(a))))])
csw09[:Mpre_lo                  ] = idxmatch(csw09[:transition],findcs( log10.(as09["M"] .- as09["ΔM"]), leq=pctl68lo_logMstarpre09 ))
csw09[:Mpre_md                  ] = idxmatch(csw09[:transition],findcs( log10.(as09["M"] .- as09["ΔM"]), gt=pctl68lo_logMstarpre09, lt=pctl68up_logMstarpre09 ))
csw09[:Mpre_hi                  ] = idxmatch(csw09[:transition],findcs( log10.(as09["M"] .- as09["ΔM"]), geq=pctl68up_logMstarpre09 ))



### No switch

#asn03    = load("/home/moon/sfortune/spinevo/data/assembly_spinmap_noswitch_03Gyr.jld", "assembly_STARS")
#println("Loaded asn03")
#nsw03   = Dict{Symbol, Any}()
#nsw03[:halo3212                 ] = findcs( asn03["ID_ISUB"          ], eq=3212 )
#nsw03[:halo1414                 ] = findcs( asn03["ID_ISUB"          ], eq=1414 )
#nsw03[:z_0                     ]  = findcs( asn03["redshift"         ], eq=minimum(asn03["redshift"]) )
#nsw03[:z_2                     ]  = findcs( asn03["redshift"         ], geq=1.5, leq=2.5 )
#nsw03[:z_hi                     ] = findcs( asn03["snapNR"           ], leq=52 ) # borders: snap52=z1.18, snap100=z0.42
#nsw03[:z_md                     ] = findcs( asn03["snapNR"           ], gt=52, lt=100 ) # borders: snap52=z1.18, snap100=z0.42
#nsw03[:z_lo                     ] = findcs( asn03["snapNR"           ], geq=100 ) # borders: snap52=z1.18, snap100=z0.42
#nsw03[:snap136                  ] = findcs( asn03["snapNR"           ], eq=136 )
#nsw03[:snap36                   ] = findcs( asn03["snapNR"           ], eq=36 )
#nsw03[:switches                 ] = findcs( asn03["switch"           ], gt=0 )
#nsw03[:fakes                    ] = findcs( asn03["FAKEFLIP"         ], eq=1 )
#nsw03[:result_disk              ] = findcs( asn03["BVAL"             ], geq=b_disk )
#nsw03[:result_int               ] = findcs( asn03["BVAL"             ], gt=b_ell, lt=b_disk )
#nsw03[:result_ell               ] = findcs( asn03["BVAL"             ], leq=b_ell )
#nsw03[:z_welker                 ] = findcs( asn03["redshift"         ], geq=1.2, leq=3.8 )
#nsw03[:z_lower                  ] = findcs( asn03["redshift"         ], lt=1.2 )
#nsw03[:z_higher                 ] = findcs( asn03["redshift"         ], gt=3.8 )
#nsw03[:z_prepeak                ] = findcs( asn03["redshift"         ], geq=z_peak )
#nsw03[:z_postpeak               ] = findcs( asn03["redshift"         ], lt=z_peak )
#nsw03[:z_hirsch_0               ] = findcs( asn03["redshift"         ], leq=0.25 )
#nsw03[:z_hirsch_05              ] = findcs( asn03["redshift"         ], gt=0.25, lt=0.75 )
#nsw03[:z_hirsch_1               ] = findcs( asn03["redshift"         ], geq=0.75, lt=1.5 )
#nsw03[:z_hirsch_2               ] = findcs( asn03["redshift"         ], geq=1.5, lt=2.5 )
#nsw03[:z_hirsch_3               ] = findcs( asn03["redshift"         ], geq=2.5, lt=3.5 )
#nsw03[:z_hirsch_4               ] = findcs( asn03["redshift"         ], geq=3.5 )
#nsw03[:mergerSTARS_acc_eq0      ] = findcs( asn03["Mpeak_MERGERS"    ], eq=0.0 )
#nsw03[:mergerSTARS_acc_g10      ] = findcs( asn03["Mpeak_MERGERS"    ] ./ asn03["M2"      ], geq=0.1)
#nsw03[:mergerSTARS_acc_g5_l10   ] = findcs( asn03["Mpeak_MERGERS"    ] ./ asn03["M2"      ], geq=0.05, lt=0.1)
#nsw03[:mergerSTARS_acc_g0_l5    ] = findcs( asn03["Mpeak_MERGERS"    ] ./ asn03["M2"      ], gt=0.0, lt=0.05)
#nsw03[:start_thr                ] = findcs( asn03["M"                ] .- asn03["ΔM"      ], geq=2e10)
#nsw03[:start_disk               ] = findcs( asn03["BVAL_0"           ] .- asn03["ΔBVAL_0" ], geq=b_disk)
#nsw03[:start_int                ] = findcs( asn03["BVAL_0"           ] .- asn03["ΔBVAL_0" ], gt=b_ell, lt=b_disk)
#nsw03[:start_ell                ] = findcs( asn03["BVAL_0"           ] .- asn03["ΔBVAL_0" ], leq=b_ell)
#nsw03[:j_hi                     ] = findcs( log10.(pyplottable( asn03["j_main"]  .- asn03["Δj_main"] )), gt=jlim2 )
#nsw03[:j_md                     ] = findcs( log10.(pyplottable( asn03["j_main"]  .- asn03["Δj_main"] )), gt=jlim1 ,leq=jlim2 )
#nsw03[:j_lo                     ] = findcs( log10.(pyplottable( asn03["j_main"]  .- asn03["Δj_main"] )), leq=jlim1 )
#nsw03[:J_hi                     ] = findcs( log10.(pyplottable( asn03["J_main"]  .- asn03["ΔJ_main"] )), gt=Jlim2 )
#nsw03[:J_md                     ] = findcs( log10.(pyplottable( asn03["J_main"]  .- asn03["ΔJ_main"] )), gt=Jlim1 ,leq=Jlim2 )
#nsw03[:J_lo                     ] = findcs( log10.(pyplottable( asn03["J_main"]  .- asn03["ΔJ_main"] )), leq=Jlim1 )
#
#nsw03[:Nmerger_0   ] = findcs( asn03["N_MERGERS"], eq=0)
#nsw03[:Nmerger_1_2 ] = findcs( asn03["N_MERGERS"], gt=0,leq=2)
#nsw03[:Nmerger_3_10] = findcs( asn03["N_MERGERS"], gt=2,leq=10)
#nsw03[:Nmerger_11_ ] = findcs( asn03["N_MERGERS"], gt=10)
#
#
#nsw03[:transition               ] = idxclude(nsw03[:start_thr], nsw03[:switches])
#
#
## Dark matter 0.3 Gyr
##ad      = load("/home/moon/sfortune/spinevo/data/assembly_spinmap_noswitch_03Gyr.jld", "assembly_DM")
#
#
#
## 0.9Gyr data set
#asn09    = load("/home/moon/sfortune/spinevo/data/assembly_spinmap_noswitch_09Gyr.jld", "assembly_STARS")
#println("Loaded asn09")
#nsw09 = Dict{Symbol, Any}()
#nsw09[:halo3212                 ] = findcs( asn09["ID_ISUB"          ], eq=3212 )
#nsw09[:halo1414                 ] = findcs( asn09["ID_ISUB"          ], eq=1414 )
#nsw09[:z_0                     ]  = findcs( asn09["redshift"         ], eq=minimum(asn09["redshift"]) )
#nsw09[:z_2                     ]  = findcs( asn09["redshift"         ], geq=1.5, leq=2.5 )
#nsw09[:z_hi                     ] = findcs( asn09["snapNR"           ], leq=52 ) # borders: snap52=z1.18, snap100=z0.42
#nsw09[:z_md                     ] = findcs( asn09["snapNR"           ], gt=52, lt=100 ) # borders: snap52=z1.18, snap100=z0.42
#nsw09[:z_lo                     ] = findcs( asn09["snapNR"           ], geq=100 ) # borders: snap52=z1.18, snap100=z0.42
#nsw09[:snap136                  ] = findcs( asn09["snapNR"           ], eq=136 )
#nsw09[:snap36                   ] = findcs( asn09["snapNR"           ], eq=36 )
#nsw09[:switches                 ] = findcs( asn09["switch"           ], gt=0 )
#nsw09[:fakes                    ] = findcs( asn09["FAKEFLIP"         ], eq=1 )
#nsw09[:result_disk              ] = findcs( asn09["BVAL"             ], geq=b_disk )
#nsw09[:result_int               ] = findcs( asn09["BVAL"             ], gt=b_ell, lt=b_disk )
#nsw09[:result_ell               ] = findcs( asn09["BVAL"             ], leq=b_ell )
#nsw09[:z_welker                 ] = findcs( asn09["redshift"         ], geq=1.2, leq=3.8 )
#nsw09[:z_lower                  ] = findcs( asn09["redshift"         ], lt=1.2 )
#nsw09[:z_higher                 ] = findcs( asn09["redshift"         ], gt=3.8 )
#nsw09[:z_prepeak                ] = findcs( asn09["redshift"         ], geq=z_peak )
#nsw09[:z_postpeak               ] = findcs( asn09["redshift"         ], lt=z_peak )
#nsw09[:z_hirsch_0               ] = findcs( asn09["redshift"         ], leq=0.25 )
#nsw09[:z_hirsch_05              ] = findcs( asn09["redshift"         ], gt=0.25, lt=0.75 )
#nsw09[:z_hirsch_1               ] = findcs( asn09["redshift"         ], geq=0.75, lt=1.5 )
#nsw09[:z_hirsch_2               ] = findcs( asn09["redshift"         ], geq=1.5, lt=2.5 )
#nsw09[:z_hirsch_3               ] = findcs( asn09["redshift"         ], geq=2.5, lt=3.5 )
#nsw09[:z_hirsch_4               ] = findcs( asn09["redshift"         ], geq=3.5 )
#nsw09[:mergerSTARS_acc_eq0      ] = findcs( asn09["Mpeak_MERGERS"    ], eq=0.0 )
#nsw09[:mergerSTARS_acc_g10      ] = findcs( asn09["Mpeak_MERGERS"    ] ./ asn09["M2"      ], geq=0.1)
#nsw09[:mergerSTARS_acc_g5_l10   ] = findcs( asn09["Mpeak_MERGERS"    ] ./ asn09["M2"      ], geq=0.05, lt=0.1)
#nsw09[:mergerSTARS_acc_g0_l5    ] = findcs( asn09["Mpeak_MERGERS"    ] ./ asn09["M2"      ], gt=0.0, lt=0.05)
#nsw09[:start_thr                ] = findcs( asn09["M"                ] .- asn09["ΔM"      ], geq=2e10)
#nsw09[:start_disk               ] = findcs( asn09["BVAL_0"           ] .- asn09["ΔBVAL_0" ], geq=b_disk)
#nsw09[:start_int                ] = findcs( asn09["BVAL_0"           ] .- asn09["ΔBVAL_0" ], gt=b_ell, lt=b_disk)
#nsw09[:start_ell                ] = findcs( asn09["BVAL_0"           ] .- asn09["ΔBVAL_0" ], leq=b_ell)
#nsw09[:j_hi                     ] = findcs( log10.(pyplottable( asn09["j_main"]  .- asn09["Δj_main"] )), gt=jlim2 )
#nsw09[:j_md                     ] = findcs( log10.(pyplottable( asn09["j_main"]  .- asn09["Δj_main"] )), gt=jlim1 ,leq=jlim2 )
#nsw09[:j_lo                     ] = findcs( log10.(pyplottable( asn09["j_main"]  .- asn09["Δj_main"] )), leq=jlim1 )
#nsw09[:J_hi                     ] = findcs( log10.(pyplottable( asn09["J_main"]  .- asn09["ΔJ_main"] )), gt=Jlim2 )
#nsw09[:J_md                     ] = findcs( log10.(pyplottable( asn09["J_main"]  .- asn09["ΔJ_main"] )), gt=Jlim1 ,leq=Jlim2 )
#nsw09[:J_lo                     ] = findcs( log10.(pyplottable( asn09["J_main"]  .- asn09["ΔJ_main"] )), leq=Jlim1 )
#
#nsw09[:Nmerger_0   ] = findcs( asn09["N_MERGERS"], eq=0)
#nsw09[:Nmerger_1_2 ] = findcs( asn09["N_MERGERS"], gt=0,leq=2)
#nsw09[:Nmerger_3_10] = findcs( asn09["N_MERGERS"], gt=2,leq=10)
#nsw09[:Nmerger_11_ ] = findcs( asn09["N_MERGERS"], gt=10)
#
#
#nsw09[:transition               ] = idxclude(nsw09[:start_thr], nsw09[:switches])

# Dark matter 0.3 Gyr
#ad09    = load("/home/moon/sfortune/spinevo/data/assembly_spinmap_noswitch_09Gyr.jld", "assembly_DM")




# extra 

endstate = load("/home/moon/sfortune/spinevo/plots/endstate.jld")
# bval_end,
# nflips30,
# nflips45,
# nflips90,
# nflips135,
# G09_nflips30,
# G09_nflips45,
# G09_nflips90,
# G09_nflips135,
# nmergers2e8,
# nmergers2e9,
# nmergers2e10,
# G09_nmergers2e8,
# G09_nmergers2e9,
# G09_nmergers2e10,
# nsw_bval_end,
# nsw_nflips30,
# nsw_nflips45,
# nsw_nflips90,
# nsw_nflips135,
# nsw_G09_nflips30,
# nsw_G09_nflips45,
# nsw_G09_nflips90,
# nsw_G09_nflips135,
# nsw_nmergers2e8,
# nsw_nmergers2e9,
# nsw_nmergers2e10,
# nsw_G09_nmergers2e8,
# nsw_G09_nmergers2e9,
# nsw_G09_nmergers2e10


function sfrd_madau(z)
    # M_⊙ * yr-1 * Mpc-3 * h3
    return 0.015*((1+z)^2.7)/(1+((1+z)/2.9)^5.6)
end

# from Magneticum
#mgntcm  = load("/home/moon/sfortune/spinevo/plots/magneticum.jld")
global_sfr_coll     = load("/home/moon/sfortune/spinevo/plots/magneticum.jld", "global_sfr_coll")#mgntcm["global_sfr_coll"]
sub1e10_sfr_coll    = load("/home/moon/sfortune/spinevo/plots/magneticum.jld", "sub1e10_sfr_coll")#mgntcm["sub1e10_sfr_coll"]
z_uhr               = load("/home/moon/sfortune/spinevo/plots/magneticum.jld", "z_uhr")#mgntcm["uhr"]
sfr_uhr             = load("/home/moon/sfortune/spinevo/plots/magneticum.jld", "sfr_uhr")#mgntcm["uhr"]
z_hr                = load("/home/moon/sfortune/spinevo/plots/magneticum.jld", "z_hr")#mgntcm["hr"]
sfr_hr              = load("/home/moon/sfortune/spinevo/plots/magneticum.jld", "sfr_hr")#mgntcm["hr"]

z_madau     = LinRange(0, 12, 1000)
sfr_madau   = sfrd_madau.(z_madau)



sfr_coll    = Dict{String, Any}()
sfr_coll["subIDs"]          = Dict{Int64, Any}()
sfr_coll["redshift"]        = Array{Float64}(undef, 0)
sfr_coll["lookbacktime"]    = Array{Float64}(undef, 0)
sfr_coll["totSFR"]          = Array{Float64}(undef, 0)
sfr_coll["mean_sSFR"]       = Array{Float64}(undef, 0)
sfr_coll["N_halos"]         = Array{Int64}(undef, 0)
sfr_coll["totM"]            = Array{Float64}(undef, 0)
sfr_coll["boxsize"]         = Array{Float64}(undef, 0)
#for i in sort(unique(as03["snapNR"])) # snaps
for i in 1:length(global_sfr_coll["redshift"])
    idcs    = findcs(as03["redshift"], eq=global_sfr_coll["redshift"][i])#, comparewith=start_thr)
    #head    = read_header("$current_dir_simbox/groups_$(lpad(i, 3, "0"))/sub_$(lpad(i, 3, "0"))")
    #idcs    = Array{Int64}(undef, 0)
    #for ii in idcs_raw
    #    if !in(as03["subID"][ii], as03["subID"][idcs])
    #        idcs    = vcat(idcs, ii)
    #    end
    #end
    if length(idcs) > 0
        println("           Accepted z = $(global_sfr_coll["redshift"][i]) with $(length(idcs)) entries")
        sfr_coll["N_halos"]     = vcat(sfr_coll["N_halos"],     length(idcs) )
        sfr_coll["redshift"]    = vcat(sfr_coll["redshift"],    mean(as03["redshift"  ][idcs]))
        sfr_coll["lookbacktime"]= vcat(sfr_coll["lookbacktime"],mean(as03["lookbacktime"][idcs]))
        #sfr_coll["boxsize"]     = vcat(sfr_coll["boxsize"],     prod(convert_units_physical([48000.0, 48000.0, 48000.0], :pos, head)) ) #kpc^3
        sfr_coll["boxsize"]     = vcat(sfr_coll["boxsize"],     global_sfr_coll["boxsize"][i] ) #kpc^3
        sfr_coll["totSFR"]      = vcat(sfr_coll["totSFR"],      sum( as03["SFR"       ][idcs]))
        sfr_coll["mean_sSFR"]   = vcat(sfr_coll["mean_sSFR"],   mean(as03["SFR"       ][idcs] ./ as03["M"][idcs]))
        sfr_coll["totM"]        = vcat(sfr_coll["totM"],        sum( as03["M"         ][idcs]))
        sfr_coll["subIDs"][length(sfr_coll["N_halos"])]  = as03["subID"][idcs]
    else
        println("Rejected z = $(global_sfr_coll["redshift"][i]) with $(length(idcs)) entries")
    end
end


# check for problems
println("\nChecking for sfr_coll duplicates...")
for i in 1:length(sfr_coll["redshift"])
    if length(unique(sfr_coll["subIDs"][i])) != length(sfr_coll["subIDs"][i])
        println()
        println("\n\n$i   ---   $(length(unique(sfr_coll["subIDs"][i])))   ---   $(length(sfr_coll["subIDs"][i]))\n")
    else
        print("$i ")
    end
end


sthr_coll    = Dict{String, Any}()
sthr_coll["subIDs"]          = Dict{Int64, Any}()
sthr_coll["redshift"]        = Array{Float64}(undef, 0)
sthr_coll["lookbacktime"]    = Array{Float64}(undef, 0)
sthr_coll["totSFR"]          = Array{Float64}(undef, 0)
sthr_coll["mean_sSFR"]       = Array{Float64}(undef, 0)
sthr_coll["N_halos"]         = Array{Int64}(undef, 0)
sthr_coll["totM"]            = Array{Float64}(undef, 0)
sthr_coll["boxsize"]         = Array{Float64}(undef, 0)
#for i in sort(unique(as03["snapNR"])) # snaps
for i in 1:length(global_sfr_coll["redshift"])
    idcs    = idxclude(findcs(as03["redshift"], eq=global_sfr_coll["redshift"][i], comparewith=csw03[:transition]), csw03[:fakes])
    #head    = read_header("$current_dir_simbox/groups_$(lpad(i, 3, "0"))/sub_$(lpad(i, 3, "0"))")
    #idcs    = Array{Int64}(undef, 0)
    #for ii in idcs_raw
    #    if !in(as03["subID"][ii], as03["subID"][idcs])
    #        idcs    = vcat(idcs, ii)
    #    end
    #end
    if length(idcs) > 0
        println("           Accepted z = $(global_sfr_coll["redshift"][i]) with $(length(idcs)) entries")
        sthr_coll["N_halos"]     = vcat(sthr_coll["N_halos"],     length(idcs) )
        sthr_coll["redshift"]    = vcat(sthr_coll["redshift"],    mean(as03["redshift"  ][idcs]))
        sthr_coll["lookbacktime"]= vcat(sthr_coll["lookbacktime"],mean(as03["lookbacktime"][idcs]))
        #sthr_coll["boxsize"]     = vcat(sthr_coll["boxsize"],     prod(convert_units_physical([48000.0, 48000.0, 48000.0], :pos, head)) ) #kpc^3
        sthr_coll["boxsize"]     = vcat(sthr_coll["boxsize"],     global_sfr_coll["boxsize"][i] ) #kpc^3
        sthr_coll["totSFR"]      = vcat(sthr_coll["totSFR"],      sum( as03["SFR"       ][idcs]))
        sthr_coll["mean_sSFR"]   = vcat(sthr_coll["mean_sSFR"],   mean(as03["SFR"       ][idcs] ./ as03["M"][idcs]))
        sthr_coll["totM"]        = vcat(sthr_coll["totM"],        sum( as03["M"         ][idcs]))
        sthr_coll["subIDs"][length(sthr_coll["N_halos"])]  = as03["subID"][idcs]
    else
        println("Rejected z = $(global_sfr_coll["redshift"][i]) with $(length(idcs)) entries")
    end
end

idxmax_sfr_coll         = findcs(sfr_coll["totSFR"],  eq=maximum(sfr_coll["totSFR"]))[end]
idxmax_global_sfr_coll  = findcs(global_sfr_coll["totSFR"],   eq=maximum(global_sfr_coll["totSFR"]))[end]
idxmax_sub1e10_sfr_coll = findcs(sub1e10_sfr_coll["totSFR"],   eq=maximum(sub1e10_sfr_coll["totSFR"]))[end]
idxmax_madau            = findcs(sfr_madau,  eq=maximum(sfr_madau))[end]
idxmax_uhr              = findcs(sfr_uhr,  eq=maximum(sfr_uhr))[end]
idxmax_hr               = findcs(sfr_hr,  eq=maximum(sfr_hr))[end]


#csw03[:stay_disk       = start_disk[  start_disk  .∈  Ref(Set(result_disk     ))]
#csw03[:stay_int        = start_int[    start_int   .∈  Ref(Set(result_int      ))]
#csw03[:stay_ell        = start_ell[    start_ell   .∈  Ref(Set(result_ell      ))]
#csw03[:result_disk_fakes    = fakes[   fakes   .∈  Ref(Set(result_disk     ))]
#csw03[:result_int_fakes     = fakes[   fakes   .∈  Ref(Set(result_int      ))]
#csw03[:result_ell_fakes     = fakes[   fakes   .∈  Ref(Set(result_ell      ))]
#csw03[:result_disk_nfakes   ]   = n_fakes[ n_fakes .∈  Ref(Set(result_disk     ))]
#csw03[:result_int_nfakes    = n_fakes[ n_fakes .∈  Ref(Set(result_int      ))]
#csw03[:result_ell_nfakes    = n_fakes[ n_fakes .∈  Ref(Set(result_ell      ))]
#csw03[:start_disk_fakes    = fakes[   fakes   .∈  Ref(Set(start_disk     ))]
#csw03[:start_int_fakes     = fakes[   fakes   .∈  Ref(Set(start_int      ))]
#csw03[:start_ell_fakes     = fakes[   fakes   .∈  Ref(Set(start_ell      ))]
#csw03[:start_disk_nfakes   = n_fakes[ n_fakes .∈  Ref(Set(start_disk     ))]
#csw03[:start_int_nfakes    = n_fakes[ n_fakes .∈  Ref(Set(start_int      ))]
#csw03[:start_ell_nfakes    = n_fakes[ n_fakes .∈  Ref(Set(start_ell      ))]
#csw03[:stay_disk_fakes     = fakes[   fakes   .∈  Ref(Set(stay_disk     ))]
#csw03[:stay_int_fakes      = fakes[   fakes   .∈  Ref(Set(stay_int      ))]
#csw03[:stay_ell_fakes      = fakes[   fakes   .∈  Ref(Set(stay_ell      ))]
#csw03[:stay_disk_nfakes    = n_fakes[ n_fakes .∈  Ref(Set(stay_disk     ))]
#csw03[:stay_int_nfakes     = n_fakes[ n_fakes .∈  Ref(Set(stay_int      ))]
#csw03[:stay_ell_nfakes     = n_fakes[ n_fakes .∈  Ref(Set(stay_ell      ))]
#csw03[:ell_to_disk_nfakes  = start_ell_nfakes[start_ell_nfakes     .∈ Ref(Set(result_disk_nfakes))]
#csw03[:ell_to_int_nfakes   = start_ell_nfakes[start_ell_nfakes     .∈ Ref(Set(result_int_nfakes))]
#csw03[:int_to_disk_nfakes  = start_int_nfakes[start_int_nfakes     .∈ Ref(Set(result_disk_nfakes))]
#csw03[:int_to_ell_nfakes   = start_int_nfakes[start_int_nfakes     .∈ Ref(Set(result_ell_nfakes))]
#csw03[:disk_to_ell_nfakes  = start_disk_nfakes[start_disk_nfakes   .∈ Ref(Set(result_ell_nfakes))]
#csw03[:disk_to_int_nfakes  = start_disk_nfakes[start_disk_nfakes   .∈ Ref(Set(result_int_nfakes))]
#csw03[:mergersDM_1to50         = findcs(ad03["M2"]./ad03["Mpeak_MM"], geq=1, leq=50, comparewith=start_thr)











#G09_start_thr   = findcs(log10.(abs.(as09["M"] .- as09["ΔM"])), geq=log10.(2e10))
## borders: snap52=z1.18, snap100=z0.42
#G09_z_hi        = findcs(as09["snapNR"      ], leq=52, comparewith=G09_start_thr)
#G09_z_md        = findcs(as09["snapNR"      ], gt=52, lt=100, comparewith=G09_start_thr)
#G09_z_lo        = findcs(as09["snapNR"      ], geq=100, comparewith=G09_start_thr)
#G09_snap136     = findcs(as09["snapNR"      ], eq=136, comparewith=G09_start_thr)
#G09_snap36      = findcs(as09["snapNR"      ], eq=36, comparewith=G09_start_thr)
#
#
#G09_halo3212    = findcs(as09["ID_ISUB"      ], eq=3212, comparewith=G09_start_thr)
#G09_halo1414    = findcs(as09["ID_ISUB"      ], eq=1414, comparewith=G09_start_thr)
#
#G09_switches    = findcs(as09["switch"      ], eq=1, comparewith=G09_start_thr)
#G09_n_switches  = findcs(as09["switch"      ], eq=0, comparewith=G09_start_thr)
#G09_fakes       = findcs(as09["FAKEFLIP"    ], eq=1, comparewith=G09_start_thr)
#G09_n_fakes     = findcs(as09["FAKEFLIP"    ], eq=0, comparewith=G09_start_thr)
#
#G09_result_disk     = findcs(as09["BVAL"  ], geq=b_disk, comparewith=G09_start_thr)
#G09_result_int      = findcs(as09["BVAL"  ], gt=b_ell, lt=b_disk, comparewith=G09_start_thr)
#G09_result_ell      = findcs(as09["BVAL"  ], leq=b_ell, comparewith=G09_start_thr)
#
#G09_start_disk      = findcs(as09["BVAL_0"  ] .- as09["ΔBVAL_0"    ], geq=b_disk, comparewith=G09_start_thr)
#G09_start_int       = findcs(as09["BVAL_0"  ] .- as09["ΔBVAL_0"    ], gt=b_ell, lt=b_disk, comparewith=G09_start_thr)
#G09_start_ell       = findcs(as09["BVAL_0"  ] .- as09["ΔBVAL_0"    ], leq=b_ell, comparewith=G09_start_thr)
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
#G09_mergerSTARS_acc_g10     = findcs(         as09["Mpeak_MERGERS"] ./ as09["M2"],       geq=0.1, comparewith=G09_start_thr)
#G09_mergerSTARS_acc_g5_l10  = findcs(         as09["Mpeak_MERGERS"] ./ as09["M2"],       geq=0.05, lt=0.1, comparewith=G09_start_thr)
#G09_mergerSTARS_acc_g0_l5   = findcs(         as09["Mpeak_MERGERS"] ./ as09["M2"],       gt=0.0, lt=0.05, comparewith=G09_start_thr)
#G09_mergerSTARS_acc_eq0     = findcs(replace( as09["Mpeak_MERGERS"], missing => 0.0),  eq=0.0, comparewith=G09_start_thr)



#mergerDM_acc_g10     = findcs( ad03["Mpeak_MERGERS"] ./ ad03["M2"],       geq=0.1, comparewith=start_thr)
#mergerDM_acc_g5_l10  = findcs( ad03["Mpeak_MERGERS"] ./ ad03["M2"],       geq=0.05, lt=0.1, comparewith=start_thr)
#mergerDM_acc_g0_l5   = findcs( ad03["Mpeak_MERGERS"] ./ ad03["M2"],       gt=0.0, lt=0.05, comparewith=start_thr)
#mergerDM_acc_eq0     = findcs( ad03["Mpeak_MERGERS"], missing => 0.0),  eq=0.0, comparewith=start_thr)
#G09_mergerDM_acc_g10     = findcs(         ad09["Mpeak_MERGERS"] ./ ad09["M2"],       geq=0.1, comparewith=G09_start_thr)
#G09_mergerDM_acc_g5_l10  = findcs(         ad09["Mpeak_MERGERS"] ./ ad09["M2"],       geq=0.05, lt=0.1, comparewith=G09_start_thr)
#G09_mergerDM_acc_g0_l5   = findcs(         ad09["Mpeak_MERGERS"] ./ ad09["M2"],       gt=0.0, lt=0.05, comparewith=G09_start_thr)
#G09_mergerDM_acc_eq0     = findcs(replace( ad09["Mpeak_MERGERS"], missing => 0.0),  eq=0.0, comparewith=G09_start_thr)



#z_pre047   = findcs(as03["redshift"    ], geq=0.47)
#z_post047  = findcs(as03["redshift"    ], lt=0.47)
#z_pre064   = findcs(as03["redshift"    ], geq=0.64)
#z_post064  = findcs(as03["redshift"    ], lt=0.64)




as          = [as03, as09]
csw         = [csw03, csw09]
as_label1   = ["03", "09"]
as_label2   = ["0.3 Gyr", "1 Gyr"]
as_label3   = ["Instant", "Long-term"]

mt          = [Int.(LinRange(14,20,7)), Int.(LinRange(23,29,7))]
mt_label1   = ["imm", "earl"]
mt_label2   = ["Immediate", "Earlier"]
