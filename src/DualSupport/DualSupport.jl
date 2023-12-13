module DualSupport

using ForwardDiff, PreallocationTools

struct CellDualData{TDC<:DiffCache, TDSC<:ForwardDiff.Dual, TDSF<:ForwardDiff.Dual}
    nvars::Int
    diff_cache_cell::TDC
    dual_sample_cell::TDSC
    diff_cache_m1::TDC
    diff_cache_00::TDC
    diff_cache_p1::TDC
    dual_sample_full::TDSF
end

function CellDualData(dual_sample_cell::Dual1, dual_sample_full::Dual2;
        is_ind_var=false, ind_var_i=0) where{Dual1<:ForwardDiff.Dual, Dual2<:ForwardDiff.Dual}
    # nvars is contained in parametric typing of dual
    nvars = typeof(dual_sample_cell).parameters[3]
    diff_cache_cell = DiffCache(zeros(1), nvars)
    diff_cache_m1 = DiffCache(zeros(1), 3*nvars)
    diff_cache_00 = DiffCache(zeros(1), 3*nvars)
    diff_cache_p1 = DiffCache(zeros(1), 3*nvars)
    cd = CellDualData(nvars, diff_cache_cell, dual_sample_cell,
            diff_cache_00, diff_cache_m1, diff_cache_p1, dual_sample_full)
    if !is_ind_var
        return cd
    end

    if ind_var_i < 1 || ind_var_i > nvars
        throw(ArgumentError("ind_var_i=$ind_var_i must be larger than one and not greater than nvars=$nvars"))
    end

    # The first entry of dual_du is the value of the variable, rest of entries are partials
    cd.diff_cache_cell.dual_du[1+ind_var_i] = 1.0
    # for these, the number of partials is nvars*3, we need to place it accordingly
    cd.diff_cache_m1.dual_du[1+ind_var_i] = 1.0
    cd.diff_cache_00.dual_du[1+nvars+ind_var_i] = 1.0
    cd.diff_cache_p1.dual_du[1+2*nvars+ind_var_i] = 1.0

    return cd
end

function Base.zero(::Type{CellDualData{TDC, TDSC, TDSF}}) where {TDC, TDSC, TDSF}
    nvars = TDSC.parameters[3]
    d1 = ForwardDiff.Dual(0.0, (zeros(nvars)...))
    d2 = ForwardDiff.Dual(0.0, (zeros(3*nvars)...))
    return CellDualData(d1,d2)
end

function Base.convert(::Type{CellDualData{TDC, TDSC, TDSF}}, x::TN) where {TDC, TDSC, TDSF, TN<:Number} 
    cd = zero(CellDualData{TDC, TDSC, TDSF})
    update_cell_dual_data_value!(cd, x)
    return cd
end

function update_cell_dual_data_value!(cd::CellDualData, value)
    cd.diff_cache_cell.du[1] = value
    cd.diff_cache_00.du[1] = value
    cd.diff_cache_m1.du[1] = value
    cd.diff_cache_p1.du[1] = value
    
    cd.diff_cache_cell.dual_du[1] = value
    cd.diff_cache_00.dual_du[1] = value
    cd.diff_cache_m1.dual_du[1] = value
    cd.diff_cache_p1.dual_du[1] = value
end

function update_cell_dual_data!(cd::CellDualData{TDC, TDSC, TDSF}, dual::TDSC) where {TDC, TDSC, TDSF}
    update_cell_dual_data_value!(cd, dual.value)
    for i in 1:cd.nvars
        cd.diff_cache_cell.dual_du[1+i] = dual.partials[i]
        cd.diff_cache_m1.dual_du[1+i] = dual.partials[i]
        cd.diff_cache_00.dual_du[1+cd.nvars+i] = dual.partials[i]
        cd.diff_cache_p1.dual_du[1+2*cd.nvars+i] = dual.partials[i]
    end
end

function get_cell_dual(cd::CellDualData)
    return get_tmp(cd.diff_cache_cell, cd.dual_sample_cell)[1]
end

function get_m1_dual(cd::CellDualData)
    return get_tmp(cd.diff_cache_m1, cd.dual_sample_full)[1]
end

function get_00_dual(cd::CellDualData)
    return get_tmp(cd.diff_cache_00, cd.dual_sample_full)[1]
end

function get_p1_dual(cd::CellDualData)
    return get_tmp(cd.diff_cache_p1, cd.dual_sample_full)[1]
end

end
