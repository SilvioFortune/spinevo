println("   Physical...")


@doc """
DESCRIPTION:\n
    - Calculate the orbital angular momentum with respect to the central subhalo
    - sfc: WARNING! This does not feel right. instead, it should be necessary to calculate everything with respect to the common center of mass
INPUT:\n
    - subID:    sub halo ID
    - centID:   central halo ID
    - snapNR:   snap number
    - boxNR:    box number
OUTPUT:\n
    - mass-independent orbital angular momentum vector
""" ->
function orbit_j(subID, centID, snapNR; boxNR=4) # Find central sub maybe using SUBFIND -> centID = FSUB[GRNR[subID+1]+1]
    if boxNR == 4
        pathtofile  = "/HydroSims/Magneticum/Box4/uhr_test/groups_$(@sprintf("%03i", snapNR))/sub_$(@sprintf("%03i", snapNR))"
    else
        error("boxNR not known.")
    end
    # Read subhalo features
    head        = read_header(pathtofile)
    spos        = convert_units_physical(read_subfind(pathtofile, "SPOS"), :pos, head)
    svel        = read_subfind(pathtofile, "SVEL") # sfc: no conversion since already physical units??
    # Return with masses added
    return ((spos[:,subID+1] .- spos[:,centID+1]) × (svel[:,subID+1] .- svel[:,centID+1]))
end

print("'orbit_j'   ")


@doc """
DESCRIPTION:\n
    - Rotation Matrix to align x with reference vector
INPUT:\n
    - x:    cartesian vector
    - ref:  reference vector
OUTPUT:\n
    - rotmat:   rotation matrix to align x with ref
""" ->
function caker(x; ref=[0,0,1])
    x       = x ./ norm(x)
    ref     = ref ./ norm(ref)
    vx      = x × ref
    cx      = transpose(x) * ref
    Vx      = [0 -vx[3] vx[2]; vx[3] 0 -vx[1]; -vx[2] vx[1] 0]
    return I + Vx + (Vx * Vx ./ (1+cx))
end

print("'caker'   ")


println()