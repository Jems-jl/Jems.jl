export MixedDualData, update_mixed_dual_data_value!, update_mixed_dual_data!,
        get_mixed_dual, get_m1_dual, get_00_dual


"""
    struct MixedDualData{TWONVARSP1, THREENVARSP1, TNUMBER}

Definition of MixedDualData, that holds the information needed to construct partial derivatives for a value that mixes local variables
from two different cells. One example could be a face pressure defined from local values of pressure defined at a cell center, although
the definition of local versus mixed does not need to correlate with face valued or cell center valued quantities.
Parametric in types `TWONVARSP1`, two times the number of independent variables plus one, `THREENVARSP1`, three times the number of
independe variables plus one, and `TNUMBER`, the type of the number used for the calculations (usually floats, but
can be duals themselves).
"""
struct MixedDualData{TWONVARSP1,THREENVARSP1,TNUMBER,DUAL_TAG}
    diff_cache_mixed::StarDiffCache{TWONVARSP1,TNUMBER,DUAL_TAG}
    diff_cache_m1::StarDiffCache{THREENVARSP1,TNUMBER,DUAL_TAG}
    diff_cache_00::StarDiffCache{THREENVARSP1,TNUMBER,DUAL_TAG}
end

"""
    function MixedDualData(nvars::Int, ::Type{TNUMBER}) where {TNUMBER}

Instantiates an object of type MixedDualData, that holds the information needed to construct partial derivatives wrt its
own properties as well as its neighbors.
"""
function MixedDualData(nvars::Int, ::Type{TNUMBER}, ::Type{DUAL_TAG}) where {TNUMBER,DUAL_TAG}
    diff_cache_mixed = StarDiffCache(2*nvars, TNUMBER, DUAL_TAG)
    diff_cache_m1 = StarDiffCache(3*nvars, TNUMBER, DUAL_TAG)
    diff_cache_00 = StarDiffCache(3*nvars, TNUMBER, DUAL_TAG)
    fd = MixedDualData{2*nvars+1,3*nvars+1,TNUMBER,DUAL_TAG}(diff_cache_mixed, 
                                diff_cache_00, diff_cache_m1)
    return fd
end

function Base.zero(::Type{MixedDualData{SIZE1,SIZE2,TNUMBER,DUAL_TAG}}) where {SIZE1, SIZE2, TNUMBER, DUAL_TAG}
    return MixedDualData((SIZE1-1)÷2, TNUMBER, DUAL_TAG)
end

function Base.convert(::Type{MixedDualData{SIZE1, SIZE2, TN1, DUAL_TAG}}, x::TN2) where {SIZE1, SIZE2, TN1<:Number, TN2<:Number, DUAL_TAG} 
    cd = zero(MixedDualData{SIZE1,SIZE2,TN1,DUAL_TAG})
    update_mixed_dual_data_value!(cd, x)
    return cd
end

function update_mixed_dual_data_value!(fd::MixedDualData, value)
    fd.diff_cache_mixed.dual_data[1] = value
    fd.diff_cache_m1.dual_data[1] = value
    fd.diff_cache_00.dual_data[1] = value
    # these are some attempts at speeding this up
    #@inbounds fd.diff_cache_mixed.dual_data[1] = value
    #@inbounds fd.diff_cache_m1.dual_data[1] = value
    #@inbounds fd.diff_cache_00.dual_data[1] = value
end

function update_mixed_dual_data!(fd::MixedDualData{SIZE1, SIZE2, TNUMBER, DUAL_TAG}, dual::TDSC) where {SIZE1, SIZE2, TNUMBER, DUAL_TAG, TDSC}
    update_mixed_dual_data_value!(fd, dual.value)
    twonvars = (SIZE1-1)
    nvars = twonvars÷2
    for i in 1:twonvars
        fd.diff_cache_mixed.dual_data[1+i] = dual.partials[i]
        fd.diff_cache_m1.dual_data[1+i] = dual.partials[i]
        fd.diff_cache_00.dual_data[1+nvars+i] = dual.partials[i]
    end
    # these are some attempts at speeding this up
    #@inbounds @views fd.diff_cache_mixed.dual_data[2:1+twonvars] .= dual.partials
    #@inbounds @views fd.diff_cache_m1.dual_data[2:1+twonvars] .= dual.partials
    #@inbounds @views fd.diff_cache_00.dual_data[2+nvars:1+nvars+twonvars] .= dual.partials
end

function get_value(fd::MixedDualData)
    return fd.diff_cache_mixed.dual_data[1]
end

function get_mixed_dual(fd::MixedDualData)
    return get_dual(fd.diff_cache_mixed)
end

function get_m1_dual(fd::MixedDualData)
    return get_dual(fd.diff_cache_m1)
end

function get_00_dual(fd::MixedDualData)
    return get_dual(fd.diff_cache_00)
end
