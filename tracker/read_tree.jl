"""
Translation of read_tree.pro @ 
"/home/moon/dolag/IDL/read_tree.pro"
#"""

using Printf

# sfc: need to understand
function IS_DEF(x)
    #aux = size(x)
    return 
end

function read_tree(
                    snr, t, t_haloes=t_haloes, swap_endian=swap_endian, 
                    sub=".", notab=notab, orig="treedata", ifile=0, quiet=quiet
                    ) # orig="treedata.orig"
    fext = string(ifile, )
    if ifile > 9
        ifile = 0
    end
    if !IS_DEF(ifile)
        ifile = 0
    end
    if !IS_DEF(ifile)
        ifile = 0
    end

    

    return nothing
end