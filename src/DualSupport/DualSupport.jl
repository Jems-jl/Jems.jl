module DualSupport

using ForwardDiff
using StaticArrays

export CellDualData, update_cell_dual_data_value!, update_cell_dual_data!,
        get_cell_dual, get_m1_dual, get_00_dual, get_p1_dual

# Inspired by DiffCache from PreallocationTools (https://github.com/SciML/PreallocationTools.jl)
struct StarDiffCache{SIZE, TNUMBER}
    dual_data::MVector{SIZE,TNUMBER}
end

function StarDiffCache(nvars::Int, ::Type{TNUMBER}) where {TNUMBER}
    StarDiffCache{nvars+1, TNUMBER}(zeros(TNUMBER, nvars+1))
end

## This uses reinterpret
#function get_dual(sdc::StarDiffCache{SIZE, TNUMBER}) where{SIZE,TNUMBER}
#    reinterpret(ForwardDiff.Dual{Nothing, TNUMBER, SIZE}, sdc.dual_data)[1]
#end

# kudos to user Mason Protter from discourse.julia.com
# beware of caveats
# https://discourse.julialang.org/t/reinterpret-vector-into-single-struct/107709
function get_dual(sdc::StarDiffCache{SIZE, TNUMBER}) where {SIZE,TNUMBER}
    p::Ptr{ForwardDiff.Dual{Nothing, TNUMBER, SIZE-1}} = pointer(sdc.dual_data)
    unsafe_load(p)         # Load the first element from that pointer
end

struct CellDualData{NVARSP1, THREENVARSP1, TNUMBER}
    diff_cache_cell::StarDiffCache{NVARSP1, TNUMBER}
    diff_cache_m1::StarDiffCache{THREENVARSP1, TNUMBER}
    diff_cache_00::StarDiffCache{THREENVARSP1, TNUMBER}
    diff_cache_p1::StarDiffCache{THREENVARSP1, TNUMBER}
end

function CellDualData(nvars::Int, ::Type{TNUMBER}; is_ind_var=false, ind_var_i=0) where{TNUMBER}
    diff_cache_cell = StarDiffCache(nvars, TNUMBER)
    diff_cache_m1 = StarDiffCache(3*nvars, TNUMBER)
    diff_cache_00 = StarDiffCache(3*nvars, TNUMBER)
    diff_cache_p1 = StarDiffCache(3*nvars, TNUMBER)
    cd = CellDualData{nvars+1, 3*nvars+1, TNUMBER}(diff_cache_cell, 
                                diff_cache_00, diff_cache_m1, diff_cache_p1)
    if !is_ind_var
        return cd
    end

    if ind_var_i < 1 || ind_var_i > nvars
        throw(ArgumentError("ind_var_i=$ind_var_i must be larger or equal to one and not greater than nvars=$nvars"))
    end

    # The first entry of dual_du is the value of the variable, rest of entries are partials
    cd.diff_cache_cell.dual_data[1+ind_var_i] = one(TNUMBER)
    # for these, the number of partials is nvars*3, we need to place it accordingly
    cd.diff_cache_m1.dual_data[1+ind_var_i] = one(TNUMBER)
    cd.diff_cache_00.dual_data[1+nvars+ind_var_i] = one(TNUMBER)
    cd.diff_cache_p1.dual_data[1+2*nvars+ind_var_i] = one(TNUMBER)

    return cd
end

function Base.zero(::Type{CellDualData{SIZE1,SIZE2,TNUMBER}}) where {SIZE1, SIZE2, TNUMBER}
    return CellDualData(SIZE1-1, TNUMBER)
end

function Base.convert(::Type{CellDualData{SIZE1, SIZE2, TN1}}, x::TN2) where {SIZE1, SIZE2, TN1<:Number, TN2<:Number} 
    cd = zero(CellDualData{SIZE1,SIZE2,TN1})
    update_cell_dual_data_value!(cd, x)
    return cd
end

function update_cell_dual_data_value!(cd::CellDualData, value)
    cd.diff_cache_cell.dual_data[1] = value
    cd.diff_cache_m1.dual_data[1] = value
    cd.diff_cache_00.dual_data[1] = value
    cd.diff_cache_p1.dual_data[1] = value
end

function update_cell_dual_data!(cd::CellDualData{SIZE1, SIZE2, TNUMBER}, dual::TDSC) where {SIZE1, SIZE2, TNUMBER, TDSC}
    update_cell_dual_data_value!(cd, dual.value)
    nvars = SIZE1-1
    for i in 1:nvars
        cd.diff_cache_cell.dual_data[1+i] = dual.partials[i]
        cd.diff_cache_m1.dual_data[1+i] = dual.partials[i]
        cd.diff_cache_00.dual_data[1+nvars+i] = dual.partials[i]
        cd.diff_cache_p1.dual_data[1+2*nvars+i] = dual.partials[i]
    end
end

function get_cell_dual(cd::CellDualData)
    return get_dual(cd.diff_cache_cell)
end

function get_m1_dual(cd::CellDualData)
    return get_dual(cd.diff_cache_m1)
end

function get_00_dual(cd::CellDualData)
    return get_dual(cd.diff_cache_00)
end

function get_p1_dual(cd::CellDualData)
    return get_dual(cd.diff_cache_p1)
end

end
