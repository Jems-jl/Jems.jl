export LocalDualData, update_local_dual_data_value!, update_local_dual_data!,
        get_local_dual, get_m1_dual, get_00_dual, get_p1_dual,
        get_mixed_00_dual, get_mixed_p1_dual


"""
    struct LocalDualData{NVARSP1, THREENVARSP1, TNUMBER}

Definition of LocalDualData, that holds the information needed to construct partial derivatives wrt its own properties
as well as its neighbors.
Parametric in types `NVARSP1`, the number of independent variables plus one, `THREENVARSP1`, three times the number of
independe variables plus one, and `TNUMBER`, the type of the number used for the calculations (usually floats, but
can be duals themselves).
"""
struct LocalDualData{NVARSP1,THREENVARSP1,TNUMBER,DUAL_TAG}
    diff_cache_local::StarDiffCache{NVARSP1,TNUMBER,DUAL_TAG}
    diff_cache_m1::StarDiffCache{THREENVARSP1,TNUMBER,DUAL_TAG}
    diff_cache_00::StarDiffCache{THREENVARSP1,TNUMBER,DUAL_TAG}
    diff_cache_p1::StarDiffCache{THREENVARSP1,TNUMBER,DUAL_TAG}
end

"""
    function LocalDualData(nvars::Int, ::Type{TNUMBER}; is_ind_var=false, ind_var_i=0)

Instantiates an object of type LocalDualData, that holds the information needed to construct partial derivatives wrt its
own properties as well as its neighbors.
Use `is_ind_var=True` and `ind_var_i=i` to instantiate a LocalDualData of a base independent variable,
with ones assigned in the appropriate spots
"""
function LocalDualData(nvars::Int, ::Type{TNUMBER}, ::Type{DUAL_TAG}; 
                        is_ind_var=false, ind_var_i=0) where{TNUMBER,DUAL_TAG<:ForwardDiff.Tag}
    diff_cache_local = StarDiffCache(nvars, TNUMBER, DUAL_TAG)
    diff_cache_m1 = StarDiffCache(3*nvars, TNUMBER, DUAL_TAG)
    diff_cache_00 = StarDiffCache(3*nvars, TNUMBER, DUAL_TAG)
    diff_cache_p1 = StarDiffCache(3*nvars, TNUMBER, DUAL_TAG)
    cd = LocalDualData{nvars+1,3*nvars+1,TNUMBER,DUAL_TAG}(diff_cache_local, 
                                diff_cache_m1, diff_cache_00, diff_cache_p1)
    if !is_ind_var
        return cd
    end

    if ind_var_i < 1 || ind_var_i > nvars
        throw(ArgumentError("ind_var_i=$ind_var_i must be larger or equal to one and not greater than nvars=$nvars"))
    end

    # The first entry of dual_du is the value of the variable, rest of entries are partials
    cd.diff_cache_local.dual_data[1+ind_var_i] = one(TNUMBER)
    # for these, the number of partials is nvars*3, we need to place it accordingly
    cd.diff_cache_m1.dual_data[1+ind_var_i] = one(TNUMBER)
    cd.diff_cache_00.dual_data[1+nvars+ind_var_i] = one(TNUMBER)
    cd.diff_cache_p1.dual_data[1+2*nvars+ind_var_i] = one(TNUMBER)

    return cd
end

"""
    function Base.zero(::Type{LocalDualData{SIZE1,SIZE2,TNUMBER,DUAL_TAG}}) where {SIZE1, SIZE2, TNUMBER, DUAL_TAG}

Instantiates a LocalDualData with zero entries (the neutral element for duals).
"""
function Base.zero(::Type{LocalDualData{SIZE1,SIZE2,TNUMBER,DUAL_TAG}}) where {SIZE1,SIZE2,TNUMBER,DUAL_TAG}
    return LocalDualData(SIZE1-1, TNUMBER, DUAL_TAG)
end

"""
    function Base.convert(::Type{LocalDualData{SIZE1, SIZE2, TN1, DUAL_TAG}}, x::TN2) where {SIZE1, SIZE2, TN1<:Number, TN2<:Number, DUAL_TAG} 

Convert `x` of type `TN2` to a LocalDualData object of types `SIZE1`, `SIZE2` and `TN1`.
"""
function Base.convert(::Type{LocalDualData{SIZE1,SIZE2,TN1,DUAL_TAG}}, x::TN2) where {SIZE1,SIZE2,TN1<:Number,TN2<:Number,DUAL_TAG<:ForwardDiff.Tag} 
    cd = zero(LocalDualData{SIZE1,SIZE2,TN1,DUAL_TAG})
    update_local_dual_data_value!(cd, x)
    return cd
end

"""
    function update_local_dual_data_value!(cd::LocalDualData, value)

Updates all data of the LocalDualData object to the given value.
"""
function update_local_dual_data_value!(cd::LocalDualData, value)
    cd.diff_cache_local.dual_data[1] = value
    cd.diff_cache_m1.dual_data[1] = value
    cd.diff_cache_00.dual_data[1] = value
    cd.diff_cache_p1.dual_data[1] = value
    # These are some attempts at speeding this up
    #@inbounds cd.diff_cache_local.dual_data[1] = value
    #@inbounds cd.diff_cache_m1.dual_data[1] = value
    #@inbounds cd.diff_cache_00.dual_data[1] = value
    #@inbounds cd.diff_cache_p1.dual_data[1] = value
end

"""
    function update_local_dual_data!(cd::LocalDualData{SIZE1, SIZE2, TNUMBER, DUAL_TAG}, dual::TDSC) where {SIZE1, SIZE2, TNUMBER, DUAL_TAG<:ForwardDiff.Tag, TDSC}

Updates all data of the LocalDualData object to the data of a given dual number.
"""
function update_local_dual_data!(cd::LocalDualData{SIZE1,SIZE2,TNUMBER,DUAL_TAG}, dual::TDSC) where {SIZE1,SIZE2,TNUMBER,DUAL_TAG<:ForwardDiff.Tag,TDSC}
    update_local_dual_data_value!(cd, dual.value)
    nvars = SIZE1-1
    for i in 1:nvars
        cd.diff_cache_local.dual_data[1+i] = dual.partials[i]
        cd.diff_cache_m1.dual_data[1+i] = dual.partials[i]
        cd.diff_cache_00.dual_data[1+nvars+i] = dual.partials[i]
        cd.diff_cache_p1.dual_data[1+2*nvars+i] = dual.partials[i]
    end
    #these are some attempts at speeding this up
    #@inbounds @views cd.diff_cache_local.dual_data[2:1+nvars] .= dual.partials
    #@inbounds @views cd.diff_cache_m1.dual_data[2:1+nvars] .= dual.partials
    #@inbounds @views cd.diff_cache_00.dual_data[2+nvars:1+2*nvars] .= dual.partials
    #@inbounds @views cd.diff_cache_p1.dual_data[2+2*nvars:1+3*nvars] .= dual.partials
end

function get_value(cd::LocalDualData)
    return cd.diff_cache_local.dual_data[1]
end

function get_local_dual(cd::LocalDualData)
    return get_dual(cd.diff_cache_local)
end

function get_m1_dual(cd::LocalDualData)
    return get_dual(cd.diff_cache_m1)
end

function get_00_dual(cd::LocalDualData)
    return get_dual(cd.diff_cache_00)
end

function get_p1_dual(cd::LocalDualData)
    return get_dual(cd.diff_cache_p1)
end

function get_mixed_00_dual(cd::LocalDualData)
    return get_mixed_dual(cd.diff_cache_m1)
end

function get_mixed_p1_dual(cd::LocalDualData)
    return get_mixed_dual(cd.diff_cache_00)
end
