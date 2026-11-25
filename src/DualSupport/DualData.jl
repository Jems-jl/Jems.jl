export DualData, get_value, update_dual_data_value!, update_dual_data_local!,
        update_dual_data_mixed_p1!, update_dual_data_mixed_m1!, update_dual_data_mixed_full!

"""
    struct DualData{NVARSP1, THREENVARSP1, TNUMBER}

Definition of DualData, that holds the information needed to construct partial derivatives wrt its own properties
as well as its neighbors.
Parametric in types `NVARSP1`, the number of independent variables plus one, `THREENVARSP1`, three times the number of
independe variables plus one, and `TNUMBER`, the type of the number used for the calculations (usually floats, but
can be duals themselves).
"""
struct DualData{NVARSP1,THREENVARSP1,TNUMBER,DUAL_TAG}
    diff_cache_local::StarDiffCache{NVARSP1,TNUMBER,DUAL_TAG}
    diff_cache_m1::StarDiffCache{THREENVARSP1,TNUMBER,DUAL_TAG}
    diff_cache_00::StarDiffCache{THREENVARSP1,TNUMBER,DUAL_TAG}
    diff_cache_p1::StarDiffCache{THREENVARSP1,TNUMBER,DUAL_TAG}
end

"""
    function DualData(nvars::Int, ::Type{TNUMBER}; is_ind_var=false, ind_var_i=0)

Instantiates an object of type DualData, that holds the information needed to construct partial derivatives wrt its
own properties as well as its neighbors.
Use `is_ind_var=True` and `ind_var_i=i` to instantiate a DualData of a base independent variable,
with ones assigned in the appropriate spots
"""
function DualData(nvars::Int, ::Type{TNUMBER}, ::Type{DUAL_TAG}; 
                        is_ind_var=false, ind_var_i=0) where{TNUMBER,DUAL_TAG<:ForwardDiff.Tag}
    diff_cache_local = StarDiffCache(nvars, TNUMBER, DUAL_TAG)
    diff_cache_m1 = StarDiffCache(3*nvars, TNUMBER, DUAL_TAG)
    diff_cache_00 = StarDiffCache(3*nvars, TNUMBER, DUAL_TAG)
    diff_cache_p1 = StarDiffCache(3*nvars, TNUMBER, DUAL_TAG)
    dd = DualData{nvars+1,3*nvars+1,TNUMBER,DUAL_TAG}(diff_cache_local, 
                                diff_cache_m1, diff_cache_00, diff_cache_p1)
    if !is_ind_var
        return dd
    end

    if ind_var_i < 1 || ind_var_i > nvars
        throw(ArgumentError("ind_var_i=$ind_var_i must be larger or equal to one and not greater than nvars=$nvars"))
    end

    # The first entry of dual_du is the value of the variable, rest of entries are partials
    dd.diff_cache_local.dual_data[1+ind_var_i] = one(TNUMBER)
    # for these, the number of partials is nvars*3, we need to place it accordingly
    dd.diff_cache_m1.dual_data[1+ind_var_i] = one(TNUMBER)
    dd.diff_cache_00.dual_data[1+nvars+ind_var_i] = one(TNUMBER)
    dd.diff_cache_p1.dual_data[1+2*nvars+ind_var_i] = one(TNUMBER)

    return dd
end

"""
    function Base.zero(::Type{DualData{SIZE1,SIZE2,TNUMBER,DUAL_TAG}}) where {SIZE1, SIZE2, TNUMBER, DUAL_TAG}

Instantiates a DualData with zero entries (the neutral element for duals).
"""
function Base.zero(::Type{DualData{SIZE1,SIZE2,TNUMBER,DUAL_TAG}}) where {SIZE1,SIZE2,TNUMBER,DUAL_TAG}
    return DualData(SIZE1-1, TNUMBER, DUAL_TAG)
end

"""
    function Base.convert(::Type{DualData{SIZE1, SIZE2, TN1, DUAL_TAG}}, x::TN2) where {SIZE1, SIZE2, TN1<:Number, TN2<:Number, DUAL_TAG} 

Convert `x` of type `TN2` to a DualData object of types `SIZE1`, `SIZE2` and `TN1`.
"""
function Base.convert(::Type{DualData{SIZE1,SIZE2,TN1,DUAL_TAG}}, x::TN2) where {SIZE1,SIZE2,TN1<:Number,TN2<:Number,DUAL_TAG<:ForwardDiff.Tag} 
    dd = zero(DualData{SIZE1,SIZE2,TN1,DUAL_TAG})
    update_dual_data_value!(dd, x)
    return dd
end

function get_value(dd::DualData)
    return dd.diff_cache_00.dual_data[1]
end

"""
    function update_local_dual_data_value!(dd::DualData, value)

Updates all data of the DualData object to the given value.
"""
function update_dual_data_value!(dd::DualData, value)
    dd.diff_cache_local.dual_data[1] = value
    dd.diff_cache_m1.dual_data[1] = value
    dd.diff_cache_00.dual_data[1] = value
    dd.diff_cache_p1.dual_data[1] = value
end

"""
    function update_local_dual_data!(dd::DualData{SIZE1, SIZE2, TNUMBER, DUAL_TAG}, dual::TDSC) where {SIZE1, SIZE2, TNUMBER, DUAL_TAG<:ForwardDiff.Tag, TDSC}

Updates all data of the DualData object to the data of a given dual number.
"""
function update_dual_data_local!(dd::DualData{SIZE1,SIZE2,TNUMBER,DUAL_TAG}, dual::TDSC) where {SIZE1,SIZE2,TNUMBER,DUAL_TAG<:ForwardDiff.Tag,TDSC}
    update_dual_data_value!(dd, dual.value)
    nvars = SIZE1-1
    for i in 1:nvars
        dd.diff_cache_local.dual_data[1+i] = dual.partials[i]
        dd.diff_cache_m1.dual_data[1+i] = dual.partials[i]
        dd.diff_cache_00.dual_data[1+nvars+i] = dual.partials[i]
        dd.diff_cache_p1.dual_data[1+2*nvars+i] = dual.partials[i]
    end
end
function update_dual_data_mixed_p1!(dd::DualData{SIZE1,SIZE2,TNUMBER,DUAL_TAG}, dual::TDSC) where {SIZE1,SIZE2,TNUMBER,DUAL_TAG<:ForwardDiff.Tag,TDSC}
    # ensure we cannot cast accidentally this in the wrong place by setting NaNs
    dd.diff_cache_local.dual_data[1] = one(TNUMBER)+NaN
    dd.diff_cache_m1.dual_data[1] = dual.value
    dd.diff_cache_00.dual_data[1] = dual.value
    dd.diff_cache_p1.dual_data[1] = one(TNUMBER)+NaN
    # we only need 2/3 of the partials, 
    nvars = SIZE1-1
    for i in 1:2*nvars
        dd.diff_cache_m1.dual_data[1+i] = dual.partials[i+nvars]
        dd.diff_cache_00.dual_data[1+nvars+i] = dual.partials[i+nvars]
    end
    # remaining entries are zero
    for i in 1:nvars
        dd.diff_cache_m1.dual_data[1+i+2*nvars] = zero(TNUMBER)
        dd.diff_cache_00.dual_data[1+i] = zero(TNUMBER)
    end
end
function update_dual_data_mixed_m1!(dd::DualData{SIZE1,SIZE2,TNUMBER,DUAL_TAG}, dual::TDSC) where {SIZE1,SIZE2,TNUMBER,DUAL_TAG<:ForwardDiff.Tag,TDSC}
    # ensure we cannot cast accidentally this in the wrong place by setting NaNs
    dd.diff_cache_local.dual_data[1] = one(TNUMBER)+NaN
    dd.diff_cache_m1.dual_data[1] = one(TNUMBER)+NaN
    dd.diff_cache_00.dual_data[1] = dual.value
    dd.diff_cache_p1.dual_data[1] = dual.value
    # we only need 2/3 of the partials, 
    nvars = SIZE1-1
    for i in 1:2*nvars
        dd.diff_cache_00.dual_data[1+i] = dual.partials[i]
        dd.diff_cache_p1.dual_data[1+nvars+i] = dual.partials[i]
    end
    # remaining entries are zero
    for i in 1:nvars
        dd.diff_cache_00.dual_data[1+i+2*nvars] = dual.partials[i]
        dd.diff_cache_p1.dual_data[1+i] = dual.partials[i]
    end
end
function update_dual_data_mixed_full!(dd::DualData{SIZE1,SIZE2,TNUMBER,DUAL_TAG}, dual::TDSC) where {SIZE1,SIZE2,TNUMBER,DUAL_TAG<:ForwardDiff.Tag,TDSC}
    # same as above, but only 00 dual is available
    dd.diff_cache_local.dual_data[1] = one(TNUMBER)+NaN
    dd.diff_cache_m1.dual_data[1] = one(TNUMBER)+NaN
    dd.diff_cache_00.dual_data[1] = dual.value
    dd.diff_cache_p1.dual_data[1] = one(TNUMBER)+NaN
    # we need all partials
    nvars = SIZE1-1
    for i in 1:3*nvars
        dd.diff_cache_00.dual_data[1+i] = dual.partials[i]
    end
end

function get_local_dual(dd::DualData)
    return get_dual(dd.diff_cache_local)
end

function get_m1_dual(dd::DualData)
    return get_dual(dd.diff_cache_m1)
end

function get_00_dual(dd::DualData)
    return get_dual(dd.diff_cache_00)
end

function get_p1_dual(dd::DualData)
    return get_dual(dd.diff_cache_p1)
end