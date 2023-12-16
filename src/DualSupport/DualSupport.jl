module DualSupport

using ForwardDiff

export CellDualData, update_cell_dual_data_value!, update_cell_dual_data!,
        get_cell_dual, get_m1_dual, get_00_dual, get_p1_dual

# Inspired by DiffCache from PreallocationTools (https://github.com/SciML/PreallocationTools.jl)
struct StarDiffCache{SIZE, TN}
    dual_data::Vector{TN}
end

function StarDiffCache(nvars::Int, ::Type{T}) where {T}
    StarDiffCache{nvars, T}(zeros(T, nvars+1))
end

function get_dual(sdc::StarDiffCache{SIZE, TN}) where{SIZE,TN}
    reinterpret(ForwardDiff.Dual{Nothing, TN, SIZE}, sdc.dual_data)[1]
end

struct CellDualData{NVARS, THREENVARS, TN}
    diff_cache_cell::StarDiffCache{NVARS, TN}
    diff_cache_m1::StarDiffCache{THREENVARS, TN}
    diff_cache_00::StarDiffCache{THREENVARS, TN}
    diff_cache_p1::StarDiffCache{THREENVARS, TN}
end

function CellDualData(nvars::Int, ::Type{TN}; is_ind_var=false, ind_var_i=0) where{TN}
    diff_cache_cell = StarDiffCache(nvars, TN)
    diff_cache_m1 = StarDiffCache(3*nvars, TN)
    diff_cache_00 = StarDiffCache(3*nvars, TN)
    diff_cache_p1 = StarDiffCache(3*nvars, TN)
    cd = CellDualData{nvars, 3*nvars, TN}(diff_cache_cell, 
                                diff_cache_00, diff_cache_m1, diff_cache_p1)
    if !is_ind_var
        return cd
    end

    if ind_var_i < 1 || ind_var_i > nvars
        throw(ArgumentError("ind_var_i=$ind_var_i must be larger than one and not greater than nvars=$nvars"))
    end

    # The first entry of dual_du is the value of the variable, rest of entries are partials
    cd.diff_cache_cell.dual_data[1+ind_var_i] = 1.0
    # for these, the number of partials is nvars*3, we need to place it accordingly
    cd.diff_cache_m1.dual_data[1+ind_var_i] = 1.0
    cd.diff_cache_00.dual_data[1+nvars+ind_var_i] = 1.0
    cd.diff_cache_p1.dual_data[1+2*nvars+ind_var_i] = 1.0

    return cd
end

function Base.zero(::Type{CellDualData{SIZE1,SIZE2,TN}}) where {SIZE1, SIZE2, TN}
    return CellDualData(SIZE1, TN)
end

function Base.convert(::Type{CellDualData{SIZE1, SIZE2, TN1}}, x::TN2) where {SIZE1, SIZE2, TN1<:Number, TN2<:Number} 
    cd = zero(CellDualData{SIZE1,SIZE2,TN1})
    update_cell_dual_data_value!(cd, x)
    return cd
end

function update_cell_dual_data_value!(cd::CellDualData, value)
    cd.diff_cache_cell.dual_data[1] = value
    cd.diff_cache_00.dual_data[1] = value
    cd.diff_cache_m1.dual_data[1] = value
    cd.diff_cache_p1.dual_data[1] = value
end

function update_cell_dual_data!(cd::CellDualData{SIZE1, SIZE2, TN}, dual::TDSC) where {SIZE1, SIZE2, TN, TDSC}
    update_cell_dual_data_value!(cd, dual.value)
    for i in 1:SIZE1
        cd.diff_cache_cell.dual_data[1+i] = dual.partials[i]
        cd.diff_cache_m1.dual_data[1+i] = dual.partials[i]
        cd.diff_cache_00.dual_data[1+SIZE1+i] = dual.partials[i]
        cd.diff_cache_p1.dual_data[1+2*SIZE1+i] = dual.partials[i]
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
