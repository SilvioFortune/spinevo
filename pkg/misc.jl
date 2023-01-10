

@doc """
DESCRIPTION:\n
INPUT:\n
OUTPUT:\n
""" ->
function endstate_quickndirty()
    stories         = readdir("/home/moon/sfortune/spinevo/data/silvio_stories_spinmap")
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

    bval_max                = Array{Float64}(undef, 3, 0) # 1 b-val, 2 z, 3 mass
    maxd_nflips30           = Array{Int64}(undef, 2, 0) # 1 pre-max, 2 post-max
    maxd_nflips45           = Array{Int64}(undef, 2, 0) # 1 pre-max, 2 post-max
    maxd_nflips90           = Array{Int64}(undef, 2, 0) # 1 pre-max, 2 post-max
    maxd_nflips135          = Array{Int64}(undef, 2, 0) # 1 pre-max, 2 post-max
    maxd_G09_nflips30       = Array{Int64}(undef, 2, 0) # 1 pre-max, 2 post-max
    maxd_G09_nflips45       = Array{Int64}(undef, 2, 0) # 1 pre-max, 2 post-max
    maxd_G09_nflips90       = Array{Int64}(undef, 2, 0) # 1 pre-max, 2 post-max
    maxd_G09_nflips135      = Array{Int64}(undef, 2, 0) # 1 pre-max, 2 post-max
    maxd_nmergers2e8        = Array{Int64}(undef, 2, 0) # 1 pre-max, 2 post-max
    maxd_nmergers2e9        = Array{Int64}(undef, 2, 0) # 1 pre-max, 2 post-max
    maxd_nmergers2e10       = Array{Int64}(undef, 2, 0) # 1 pre-max, 2 post-max
    maxd_G09_nmergers2e8    = Array{Int64}(undef, 2, 0) # 1 pre-max, 2 post-max
    maxd_G09_nmergers2e9    = Array{Int64}(undef, 2, 0) # 1 pre-max, 2 post-max
    maxd_G09_nmergers2e10   = Array{Int64}(undef, 2, 0) # 1 pre-max, 2 post-max

    sstacked03    = Array{Int64}(undef, 0)
    istacked03    = Array{Int64}(undef, 0)
    istacked09    = Array{Int64}(undef, 0)
    sstacked09    = Array{Int64}(undef, 0)

    #nsw_bval_end        = Array{Float64}(undef, 0)
    #nsw_nflips30        = Array{Int64}(undef, 0)
    #nsw_nflips45        = Array{Int64}(undef, 0)
    #nsw_nflips90        = Array{Int64}(undef, 0)
    #nsw_nflips135       = Array{Int64}(undef, 0)
    #nsw_G09_nflips30    = Array{Int64}(undef, 0)
    #nsw_G09_nflips45    = Array{Int64}(undef, 0)
    #nsw_G09_nflips90    = Array{Int64}(undef, 0)
    #nsw_G09_nflips135   = Array{Int64}(undef, 0)
    #nsw_nmergers2e8     = Array{Int64}(undef, 0)
    #nsw_nmergers2e9     = Array{Int64}(undef, 0)
    #nsw_nmergers2e10    = Array{Int64}(undef, 0)
    #nsw_G09_nmergers2e8 = Array{Int64}(undef, 0)
    #nsw_G09_nmergers2e9 = Array{Int64}(undef, 0)
    #nsw_G09_nmergers2e10= Array{Int64}(undef, 0)
    println("Getting endstate data for $(length(stories)) entries:")
    flush(stdout)
    for i in 1:length(stories)
        print("$(i)", ifelse(i % 30 == 0, "\n", " "))
        flush(stdout)
        mcs03 = load(joinpath("/home/moon/sfortune/spinevo/data/silvio_stories_spinmap", stories[i]), "merger_collection_STARS")
        mcs09 = load(joinpath("/home/moon/sfortune/spinevo/data/silvio_stories_spinmap_09Gyr", stories[i]), "merger_collection_STARS")
        noswitch03  = findcs(mcs03["switch"         ], eq=0)
        noswitch09  = findcs(mcs09["switch"         ], eq=0)
        start_thr03 = findcs(mcs03["M"] .- mcs03["ΔM"], geq=2e10)
        start_thr09 = findcs(mcs09["M"] .- mcs09["ΔM"], geq=2e10)
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
        print(" -- ", length(mcs03["BVAL_0"])," ", length(mcs09["BVAL_0"]))
        slct03 = idxmatch(idxmatch(idxmatch(noswitch03,nofake03),end_thr03), start_thr03)
        slct09 = idxmatch(idxmatch(idxmatch(noswitch09,nofake09),end_thr09), start_thr09)
        println(" --- ", length(slct03)," ", length(slct09), " ", maximum(mcs03["BVAL_0"][slct03]), maximum(mcs09["BVAL_0"][slct09]))
        maxbvalpos03  = findcs(mcs03["BVAL_0"         ], eq=maximum(mcs03["BVAL_0"][slct03]))[end]
        maxbvalpos09  = findcs(mcs09["BVAL_0"         ], eq=maximum(mcs09["BVAL_0"][slct09]))[end]
        sstacked03   = vcat( sstacked03, mcs03["snapNR"][slct03] )
        sstacked09   = vcat( sstacked09, mcs09["snapNR"][slct09] )
        istacked09   = vcat( istacked09, mcs09["subID" ][slct09] )
        istacked03   = vcat( istacked03, mcs03["subID" ][slct03] )

        mmi03       = Array{Int64}(undef, 0)
        premaxdmmi03   = Array{Int64}(undef, 0) 
        postmaxdmmi03   = Array{Int64}(undef, 0) 
        mmi         = mergermap_indices(mcs03)
        for ii in 1:length(mmi["main"])
            if mmi["main"][ii] ∈ slct03
                mmi03   = vcat(mmi03, mmi["first"][ii]:mmi["last"][ii])
                if mmi["main"][ii] ≤ maxbvalpos03
                    premaxdmmi03   = vcat(premaxdmmi03, mmi["first"][ii]:mmi["last"][ii])
                else
                    postmaxdmmi03   = vcat(postmaxdmmi03, mmi["first"][ii]:mmi["last"][ii])
                end
            end
        end
        mmi09       = Array{Int64}(undef, 0)
        premaxdmmi09   = Array{Int64}(undef, 0) 
        postmaxdmmi09   = Array{Int64}(undef, 0) 
        mmi         = mergermap_indices(mcs09)
        for ii in 1:length(mmi["main"])
            if mmi["main"][ii] ∈ slct09
                mmi09   = vcat(mmi09, mmi["first"][ii]:mmi["last"][ii])
                if mmi["main"][ii] ≤ maxbvalpos03
                    premaxdmmi09   = vcat(premaxdmmi09, mmi["first"][ii]:mmi["last"][ii])
                else
                    postmaxdmmi09   = vcat(postmaxdmmi09, mmi["first"][ii]:mmi["last"][ii])
                end
            end
        end

        #nsw_mcs03 = load(joinpath("/home/moon/sfortune/spinevo/data/silvio_stories_spinmap_noswitch", stories[i]), "merger_collection_STARS")
        #nsw_mcs09 = load(joinpath("/home/moon/sfortune/spinevo/data/silvio_stories_spinmap_noswitch_09Gyr", stories[i]), "merger_collection_STARS")
        bval_end        = vcat( bval_end        , mcs03["BVAL_0"][end])
        nflips30        = vcat( nflips30        , length(findcs(mcs03["ϕ_flip"    ][slct03],geq=30   ))     )
        nflips45        = vcat( nflips45        , length(findcs(mcs03["ϕ_flip"    ][slct03],geq=45   ))     )
        nflips90        = vcat( nflips90        , length(findcs(mcs03["ϕ_flip"    ][slct03],geq=90   ))     )
        nflips135       = vcat( nflips135       , length(findcs(mcs03["ϕ_flip"    ][slct03],geq=135  ))     )
        G09_nflips30    = vcat( G09_nflips30    , length(findcs(mcs09["ϕ_flip"    ][slct09],geq=30   ))     )
        G09_nflips45    = vcat( G09_nflips45    , length(findcs(mcs09["ϕ_flip"    ][slct09],geq=45   ))     )
        G09_nflips90    = vcat( G09_nflips90    , length(findcs(mcs09["ϕ_flip"    ][slct09],geq=90   ))     )
        G09_nflips135   = vcat( G09_nflips135   , length(findcs(mcs09["ϕ_flip"    ][slct09],geq=135  ))     )
        nmergers2e8     = vcat( nmergers2e8     , length(findcs(mcs03["merger map"][9,mmi03],geq=2e8 ))     )
        nmergers2e9     = vcat( nmergers2e9     , length(findcs(mcs03["merger map"][9,mmi03],geq=2e9 ))     )
        nmergers2e10    = vcat( nmergers2e10    , length(findcs(mcs03["merger map"][9,mmi03],geq=2e10))     )
        G09_nmergers2e8 = vcat( G09_nmergers2e8 , length(findcs(mcs09["merger map"][9,mmi09],geq=2e8 ))     )
        G09_nmergers2e9 = vcat( G09_nmergers2e9 , length(findcs(mcs09["merger map"][9,mmi09],geq=2e9 ))     )
        G09_nmergers2e10= vcat( G09_nmergers2e10, length(findcs(mcs09["merger map"][9,mmi09],geq=2e10))     )

        bval_max                = hcat( bval_max, [mcs03["BVAL_0"][maxbvalpos03], mcs03["redshift"][maxbvalpos03], mcs03["M"][maxbvalpos03]] )
        maxd_nflips30           = hcat( maxd_nflips30           , [ length(idxmatch(findcs(mcs03["ϕ_flip"    ][slct03],geq=30   ), findcs(mcs03["redshift"][slct03],geq=mcs03["redshift"][maxbvalpos03]))), length(idxmatch(findcs(mcs03["ϕ_flip"    ][slct03],geq=30   ), findcs(mcs03["redshift"][slct03],lt=mcs03["redshift"][maxbvalpos03]))) ] )
        maxd_nflips45           = hcat( maxd_nflips45           , [ length(idxmatch(findcs(mcs03["ϕ_flip"    ][slct03],geq=45   ), findcs(mcs03["redshift"][slct03],geq=mcs03["redshift"][maxbvalpos03]))), length(idxmatch(findcs(mcs03["ϕ_flip"    ][slct03],geq=45   ), findcs(mcs03["redshift"][slct03],lt=mcs03["redshift"][maxbvalpos03]))) ] )
        maxd_nflips90           = hcat( maxd_nflips90           , [ length(idxmatch(findcs(mcs03["ϕ_flip"    ][slct03],geq=90   ), findcs(mcs03["redshift"][slct03],geq=mcs03["redshift"][maxbvalpos03]))), length(idxmatch(findcs(mcs03["ϕ_flip"    ][slct03],geq=90   ), findcs(mcs03["redshift"][slct03],lt=mcs03["redshift"][maxbvalpos03]))) ] )
        maxd_nflips135          = hcat( maxd_nflips135          , [ length(idxmatch(findcs(mcs03["ϕ_flip"    ][slct03],geq=135  ), findcs(mcs03["redshift"][slct03],geq=mcs03["redshift"][maxbvalpos03]))), length(idxmatch(findcs(mcs03["ϕ_flip"    ][slct03],geq=135  ), findcs(mcs03["redshift"][slct03],lt=mcs03["redshift"][maxbvalpos03]))) ] )
        maxd_G09_nflips30       = hcat( maxd_G09_nflips30       , [ length(idxmatch(findcs(mcs09["ϕ_flip"    ][slct09],geq=30   ), findcs(mcs09["redshift"][slct09],geq=mcs03["redshift"][maxbvalpos03]))), length(idxmatch(findcs(mcs09["ϕ_flip"    ][slct09],geq=30   ), findcs(mcs09["redshift"][slct09],lt=mcs03["redshift"][maxbvalpos03]))) ] )
        maxd_G09_nflips45       = hcat( maxd_G09_nflips45       , [ length(idxmatch(findcs(mcs09["ϕ_flip"    ][slct09],geq=45   ), findcs(mcs09["redshift"][slct09],geq=mcs03["redshift"][maxbvalpos03]))), length(idxmatch(findcs(mcs09["ϕ_flip"    ][slct09],geq=45   ), findcs(mcs09["redshift"][slct09],lt=mcs03["redshift"][maxbvalpos03]))) ] )
        maxd_G09_nflips90       = hcat( maxd_G09_nflips90       , [ length(idxmatch(findcs(mcs09["ϕ_flip"    ][slct09],geq=90   ), findcs(mcs09["redshift"][slct09],geq=mcs03["redshift"][maxbvalpos03]))), length(idxmatch(findcs(mcs09["ϕ_flip"    ][slct09],geq=90   ), findcs(mcs09["redshift"][slct09],lt=mcs03["redshift"][maxbvalpos03]))) ] )
        maxd_G09_nflips135      = hcat( maxd_G09_nflips135      , [ length(idxmatch(findcs(mcs09["ϕ_flip"    ][slct09],geq=135  ), findcs(mcs09["redshift"][slct09],geq=mcs03["redshift"][maxbvalpos03]))), length(idxmatch(findcs(mcs09["ϕ_flip"    ][slct09],geq=135  ), findcs(mcs09["redshift"][slct09],lt=mcs03["redshift"][maxbvalpos03]))) ] )
        maxd_nmergers2e8        = hcat( maxd_nmergers2e8        , [ length(findcs(mcs03["merger map"][9,premaxdmmi03],geq=2e8 )), length(findcs(mcs03["merger map"][9,postmaxdmmi03],geq=2e8 )) ] )
        maxd_nmergers2e9        = hcat( maxd_nmergers2e9        , [ length(findcs(mcs03["merger map"][9,premaxdmmi03],geq=2e9 )), length(findcs(mcs03["merger map"][9,postmaxdmmi03],geq=2e9 )) ] )
        maxd_nmergers2e10       = hcat( maxd_nmergers2e10       , [ length(findcs(mcs03["merger map"][9,premaxdmmi03],geq=2e10)), length(findcs(mcs03["merger map"][9,postmaxdmmi03],geq=2e10)) ] )
        maxd_G09_nmergers2e8    = hcat( maxd_G09_nmergers2e8    , [ length(findcs(mcs09["merger map"][9,premaxdmmi09],geq=2e8 )), length(findcs(mcs09["merger map"][9,postmaxdmmi09],geq=2e8 )) ] )
        maxd_G09_nmergers2e9    = hcat( maxd_G09_nmergers2e9    , [ length(findcs(mcs09["merger map"][9,premaxdmmi09],geq=2e9 )), length(findcs(mcs09["merger map"][9,postmaxdmmi09],geq=2e9 )) ] )
        maxd_G09_nmergers2e10   = hcat( maxd_G09_nmergers2e10   , [ length(findcs(mcs09["merger map"][9,premaxdmmi09],geq=2e10)), length(findcs(mcs09["merger map"][9,postmaxdmmi09],geq=2e10)) ] )

        #nsw_bval_end        = vcat( nsw_bval_end        , nsw_mcs03["BVAL_0"][end])
        #nsw_nflips30        = vcat( nsw_nflips30        , length(findcs(nsw_mcs03["ϕ_flip"         ],geq=30  )) )
        #nsw_nflips45        = vcat( nsw_nflips45        , length(findcs(nsw_mcs03["ϕ_flip"         ],geq=45  )) )
        #nsw_nflips90        = vcat( nsw_nflips90        , length(findcs(nsw_mcs03["ϕ_flip"         ],geq=90  )) )
        #nsw_nflips135       = vcat( nsw_nflips135       , length(findcs(nsw_mcs03["ϕ_flip"         ],geq=135 )) )
        #nsw_G09_nflips30    = vcat( nsw_G09_nflips30    , length(findcs(nsw_mcs09["ϕ_flip"         ],geq=30  )) )
        #nsw_G09_nflips45    = vcat( nsw_G09_nflips45    , length(findcs(nsw_mcs09["ϕ_flip"         ],geq=45  )) )
        #nsw_G09_nflips90    = vcat( nsw_G09_nflips90    , length(findcs(nsw_mcs09["ϕ_flip"         ],geq=90  )) )
        #nsw_G09_nflips135   = vcat( nsw_G09_nflips135   , length(findcs(nsw_mcs09["ϕ_flip"         ],geq=135 )) )
        #nsw_nmergers2e8     = vcat( nsw_nmergers2e8     , length(findcs(nsw_mcs03["merger map"][9,:],geq=2e8 )) )
        #nsw_nmergers2e9     = vcat( nsw_nmergers2e9     , length(findcs(nsw_mcs03["merger map"][9,:],geq=2e9 )) )
        #nsw_nmergers2e10    = vcat( nsw_nmergers2e10    , length(findcs(nsw_mcs03["merger map"][9,:],geq=2e10)) )
        #nsw_G09_nmergers2e8 = vcat( nsw_G09_nmergers2e8 , length(findcs(nsw_mcs09["merger map"][9,:],geq=2e8 )) )
        #nsw_G09_nmergers2e9 = vcat( nsw_G09_nmergers2e9 , length(findcs(nsw_mcs09["merger map"][9,:],geq=2e9 )) )
        #nsw_G09_nmergers2e10= vcat( nsw_G09_nmergers2e10, length(findcs(nsw_mcs09["merger map"][9,:],geq=2e10)) )
    end
    save("/home/moon/sfortune/spinevo/plots/endstate.jld", 
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
        "G09_nmergers2e10", G09_nmergers2e10,
        "bval_max"                , bval_max                ,
        "maxd_nflips30"           , maxd_nflips30           ,
        "maxd_nflips45"           , maxd_nflips45           ,
        "maxd_nflips90"           , maxd_nflips90           ,
        "maxd_nflips135"          , maxd_nflips135          ,
        "maxd_G09_nflips30"       , maxd_G09_nflips30       ,
        "maxd_G09_nflips45"       , maxd_G09_nflips45       ,
        "maxd_G09_nflips90"       , maxd_G09_nflips90       ,
        "maxd_G09_nflips135"      , maxd_G09_nflips135      ,
        "maxd_nmergers2e8"        , maxd_nmergers2e8        ,
        "maxd_nmergers2e9"        , maxd_nmergers2e9        ,
        "maxd_nmergers2e10"       , maxd_nmergers2e10       ,
        "maxd_G09_nmergers2e8"    , maxd_G09_nmergers2e8    ,
        "maxd_G09_nmergers2e9"    , maxd_G09_nmergers2e9    ,
        "maxd_G09_nmergers2e10"   , maxd_G09_nmergers2e10   






















        #"nsw_bval_end",         nsw_bval_end,
        #"nsw_nflips30",         nsw_nflips30,
        #"nsw_nflips45",         nsw_nflips45,
        #"nsw_nflips90",         nsw_nflips90,
        #"nsw_nflips135",        nsw_nflips135,
        #"nsw_G09_nflips30",     nsw_G09_nflips30,
        #"nsw_G09_nflips45",     nsw_G09_nflips45,
        #"nsw_G09_nflips90",     nsw_G09_nflips90,
        #"nsw_G09_nflips135",    nsw_G09_nflips135,
        #"nsw_nmergers2e8",      nsw_nmergers2e8,
        #"nsw_nmergers2e9",      nsw_nmergers2e9,
        #"nsw_nmergers2e10",     nsw_nmergers2e10,
        #"nsw_G09_nmergers2e8",  nsw_G09_nmergers2e8,
        #"nsw_G09_nmergers2e9",  nsw_G09_nmergers2e9,
        #"nsw_G09_nmergers2e10", nsw_G09_nmergers2e10
    )

    println("\n\n\n---------------------------\n\nNow witness the firepower of this fully armed and operational battle station!")
    return nothing
end

print("'endstate_quickndirty'   ")