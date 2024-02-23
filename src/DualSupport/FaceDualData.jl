export FaceDualData, update_face_dual_data_value!, update_face_dual_data!,
        get_face_dual, get_m1_dual, get_00_dual


"""
    struct FaceDualData{TWONVARSP1, THREENVARSP1, TNUMBER}

Definition of FaceDualData, that holds the information needed to construct partial derivatives wrt its own properties
(the face) as well as its neighbors.
Parametric in types `TWONVARSP1`, two times the number of independent variables plus one, `THREENVARSP1`, three times the number of
independe variables plus one, and `TNUMBER`, the type of the number used for the calculations (usually floats, but
can be duals themselves).
"""
struct FaceDualData{TWONVARSP1, THREENVARSP1, TNUMBER}
    diff_cache_face::StarDiffCache{TWONVARSP1, TNUMBER}
    diff_cache_m1::StarDiffCache{THREENVARSP1, TNUMBER}
    diff_cache_00::StarDiffCache{THREENVARSP1, TNUMBER}
end

"""
    function CellDualData(nvars::Int, ::Type{TNUMBER}; is_ind_var=false, ind_var_i=0)

Instantiates an object of type FaceDualData, that holds the information needed to construct partial derivatives wrt its
own properties as well as its neighbors.
"""
function FaceDualData(nvars::Int, ::Type{TNUMBER}) where {TNUMBER}
    diff_cache_face = StarDiffCache(2*nvars, TNUMBER)
    diff_cache_m1 = StarDiffCache(3*nvars, TNUMBER)
    diff_cache_00 = StarDiffCache(3*nvars, TNUMBER)
    fd = FaceDualData{2*nvars+1, 3*nvars+1, TNUMBER}(diff_cache_face, 
                                diff_cache_00, diff_cache_m1)
    return fd
end

function Base.zero(::Type{FaceDualData{SIZE1,SIZE2,TNUMBER}}) where {SIZE1, SIZE2, TNUMBER}
    return FaceDualData((SIZE1-1)÷2, TNUMBER)
end

function Base.convert(::Type{FaceDualData{SIZE1, SIZE2, TN1}}, x::TN2) where {SIZE1, SIZE2, TN1<:Number, TN2<:Number} 
    cd = zero(FaceDualData{SIZE1,SIZE2,TN1})
    update_face_dual_data_value!(cd, x)
    return cd
end

function update_face_dual_data_value!(fd::FaceDualData, value)
    fd.diff_cache_face.dual_data[1] = value
    fd.diff_cache_m1.dual_data[1] = value
    fd.diff_cache_00.dual_data[1] = value
end

function update_face_dual_data!(fd::FaceDualData{SIZE1, SIZE2, TNUMBER}, dual::TDSC) where {SIZE1, SIZE2, TNUMBER, TDSC}
    update_face_dual_data_value!(fd, dual.value)
    twonvars = (SIZE1-1)
    nvars = twonvars÷2
    for i in 1:twonvars
        fd.diff_cache_face.dual_data[1+i] = dual.partials[i]
        fd.diff_cache_m1.dual_data[1+i] = dual.partials[i]
        fd.diff_cache_00.dual_data[1+nvars+i] = dual.partials[i]
    end
end

function get_face_dual(fd::FaceDualData)
    return get_dual(fd.diff_cache_face)
end

function get_m1_dual(fd::FaceDualData)
    return get_dual(fd.diff_cache_m1)
end

function get_00_dual(fd::FaceDualData)
    return get_dual(fd.diff_cache_00)
end