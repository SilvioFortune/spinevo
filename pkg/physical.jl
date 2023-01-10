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
function orbit_j(subID, centID, snapNR; simbox=current_dir_simbox) # Find central sub maybe using SUBFIND -> centID = FSUB[GRNR[subID+1]+1]
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
INPUT:\n
OUTPUT:\n
""" ->
function relative_j(gsub, gcent)
    return ( (read_galaxy_prop(gsub, "SPOS", :physical) .-  read_galaxy_prop(gcent, "SPOS", :physical)) × ( read_galaxy_prop(gsub, "SVEL", :physical) .-  read_galaxy_prop(gcent, "SVEL", :physical)) )
end

print("'relative_j_fromtree'   ")


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


@doc """
DESCRIPTION:\n
INPUT:\n
OUTPUT:\n
""" ->
function get_radial_velocity(pos,vel)
    return dot(normalize(pos), vel)
end

print("'get_radial_velocity'   ")





@doc """
DESCRIPTION:\n
INPUT:\n
OUTPUT:\n
""" ->
function min_dj_v_flipangled(x)
    if 0 ≤ x ≤ 90
            return sind(x)
    elseif 90 < x ≤ 180
            return 1
    else
            error("x = $(x) is not a value between 0 and 180")
    end
end

print("'min_dj_v_flipangled'   ")



@doc """
DESCRIPTION:\n
INPUT:\n
OUTPUT:\n
""" ->
function min_dj_v_flipangle(x)
    if 0 ≤ x ≤ π/2
            return sin(x)
    elseif π/2 < x ≤ π
            return 1
    else
            error("x = $(x) is not a value between 0 and π")
    end
end

print("'min_dj_v_flipangle'   ")





@doc """
DESCRIPTION:\n
INPUT:\n
OUTPUT:\n
""" ->
function max_flipangled_v_dj(x)
    #return ifelse(-1 ≤ x ≤ 1, asind(x), NaN)
    if -1 ≤ x ≤ 1
            return asind(x)
    else
            return NaN
    end
end

print("'max_flipangled_v_dj'   ")



@doc """
DESCRIPTION:\n
INPUT:\n
OUTPUT:\n
""" ->
function max_flipangle_v_dj(x)
    #return ifelse(-1 ≤ x ≤ 1, asin(x), NaN)
    if -1 ≤ x ≤ 1
            return asin(x)
    else
            return NaN
    end
end

print("'max_flipangle_v_dj'   ")




##########################################################################
# average
@doc """
DESCRIPTION:\n
INPUT:\n
OUTPUT:\n
""" ->
function avg_dj_v_flipangled(x)
    return sqrt( (cosd(x))^(-2) - 1)
end

print("'avg_dj_v_flipangled'   ")


@doc """
DESCRIPTION:\n
INPUT:\n
OUTPUT:\n
""" ->
function avg_dj_v_flipangle(x)
    return sqrt( (cos(x))^(-2) - 1)
end

print("'avg_dj_v_flipangle'   ")





@doc """
DESCRIPTION:\n
INPUT:\n
OUTPUT:\n
""" ->
function avg_flipangled_v_dj(x)
    return acosd(1/sqrt(1+(x*x)))
end

print("'avg_flipangled_v_dj'   ")

@doc """
DESCRIPTION:\n
INPUT:\n
OUTPUT:\n
""" ->
function avg_flipangle_v_dj(x)
    return acos(1/sqrt(1+(x*x)))
end

print("'avg_flipangle_v_dj'   ")


println()