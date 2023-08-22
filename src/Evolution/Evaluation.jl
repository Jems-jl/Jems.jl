"""
    eval_cell_eqs(sm::StellarModel, k::Int, ind_vars_view::Vector{<:TT}) where{TT<:Real}

Evaluates the stellar structure equations of the stellar model, `sm`, at cell `k`, given the view of the independent
variables, `ind_vars_view`.
"""
function eval_cell_eqs!(sm::StellarModel, k::Int, result_view::AbstractVector{<:TT},
                        varm1::AbstractVector{<:TT}, var00::AbstractVector{<:TT}, varp1::AbstractVector{<:TT},
                        eosm1::AbstractVector{<:TT}, eos00::AbstractVector{<:TT}, eosp1::AbstractVector{<:TT},
                        κm1::TT, κ00::TT, κp1::TT) where {TT<:Real}
    for i = 1:(sm.nvars)
        result_view[i] = sm.structure_equations[i](sm, k, varm1, var00, varp1, eosm1, eos00, eosp1, κm1, κ00,
                                                   κp1)
    end
end

"""
    eval_eqs(sm::StellarModel)

Evaluates the stellar structure equations of the stellar model, `sm`, for all cells.
"""
function eval_eqs!(sm::StellarModel)
    for k = 1:(sm.nz)
        ki = 0
        kf = 0
        if k == 1
            ki = sm.nvars * (k - 1) + 1
            kf = sm.nvars * (k + 1)
        elseif k == sm.nz
            ki = sm.nvars * (k - 2) + 1
            kf = sm.nvars * (k)
        else
            ki = sm.nvars * (k - 2) + 1
            kf = sm.nvars * (k + 1)
        end
        resultview = view(sm.eqs_nums, ki:kf)
        cacheview = view(sm.diff_caches, k, :)
        varm1 = get_tmp(cacheview[1], resultview[1])
        var00 = get_tmp(cacheview[2], resultview[1])
        varp1 = get_tmp(cacheview[3], resultview[1])
        # collect eos and κ info
        eos00 = [dual.value for dual in view(sm.eos_results, k, :)]
        κ00 = sm.opacity_results[k].value
        if k != 1
            eosm1 = [dual.value for dual in view(sm.eos_results, k - 1, :)]
            κm1 = sm.opacity_results[k - 1].value
        else
            eosm1 = [NaN for dual in 1:length(sm.eos_results[k, :])]
            κm1 = NaN
        end
        if k != sm.nz
            eosp1 = [dual.value for dual in view(sm.eos_results, k + 1, :)]
            κp1 = sm.opacity_results[k + 1].value
        else
            eosp1 = [Nan for dual in 1:length(sm.eos_results[k, :])]
            κp1 = NaN
        end
        eval_cell_eqs!(sm, k, resultview, varm1, var00, varp1, eosm1, eos00, eosp1, κm1, κ00, κp1)

    end
end

"""
    eval_jacobian_row!(sm::StellarModel, k::int)

Evaluates row `k` of the Jacobian matrix of the given StellarModel `sm`.
"""
function eval_jacobian_row!(sm::StellarModel, k::Int)
    # ranges of ind_vars vector that needs to be considered, needs special cases for k=1 or nz
    # Jacobian has tridiagonal block structure:
    # x x - - - -
    # x x x - - -
    # - x x x - -
    # - - x x x -
    # - - - x x x
    # - - - - x x
    # In the first row, the first block corresponds to the derivatives of the structure equations with respect to the
    # variables at the first cell. The block to the right corresponds to the derivatives with respect to the variables
    # at the second cell. In the second row until row nz-1, the first block contains the derivatives of the structure
    # equations wrt. the variables of the previous cell, the second block wrt. to the varibles in the current cell, and
    # the third wrt. the variables of the next cell. The last row only contains the derivatives wrt. the penultimate
    # cell and the last cell.
    _set_diff_cache!(sm, k)
    ki = 0
    kf = 0
    if k == 1
        ki = sm.nvars * (k - 1) + 1
        kf = sm.nvars * (k + 1)
    elseif k == sm.nz
        ki = sm.nvars * (k - 2) + 1
        kf = sm.nvars * (k)
    else
        ki = sm.nvars * (k - 2) + 1
        kf = sm.nvars * (k + 1)
    end
    jac_view = view(sm.jacobian, (sm.nvars * (k - 1) + 1):(sm.nvars * k), ki:kf)
    resultview = view(sm.eqs, k, :)
    cacheview = view(sm.diff_caches, k, :)
    varm1 = get_tmp(cacheview[1], resultview[1])
    var00 = get_tmp(cacheview[2], resultview[1])
    varp1 = get_tmp(cacheview[3], resultview[1])

    # collect eos and κ info
    eos00 = view(sm.eos_results, k, :)
    κ00 = sm.opacity_results[k]
    if k != 1
        eosm1 = view(sm.eos_results, k - 1, :)
        κm1 = sm.opacity_results[k - 1]
    else
        eosm1 = Vector{typeof(eos00[1])}(undef, size(sm.eos_results[k, :]))
        κm1 = typeof(κ00)(NaN)
    end
    if k != sm.nz
        eosp1 = view(sm.eos_results, k + 1, :)
        κp1 = sm.opacity_results[k + 1]
    else
        eosp1 = Vector{typeof(eos00[1])}(undef, size(sm.eos_results[k, :]))
        κp1 = typeof(κ00)(NaN)
    end
    eval_cell_eqs!(sm, k, resultview, varm1, var00, varp1, eosm1, eos00, eosp1, κm1, κ00, κp1)

    # FOR SOME REASON, STORING RESULTS IN THE JACOBIAN IS CRAZY MEMORY EXPENSIVE
    # store results
    for i = 1:(sm.nvars)
        for j = 1:(k == 1 || k == sm.nz ? 2 * sm.nvars : 3 * sm.nvars)
            if (k == 1)  # for k==1 the correct derivatives are displaced!
                jac_view[i, j] = resultview[i].partials[j + sm.nvars]
            else
                jac_view[i, j] = resultview[i].partials[j]
            end
        end
    end
end

"""
    eval_jacobian!(sm::StellarModel)

Evaluates the whole Jacobian matrix of the given StellarModel `sm`.
"""
function eval_jacobian!(sm::StellarModel)
    Threads.@threads for k = 1:(sm.nz)
        eval_jacobian_row!(sm, k)
    end
end

function eval_info!(sm::StellarModel)
    tt = sm.eqs[1]  # get example of type of dual
    eval_eos!(sm, tt)
    eval_opacity!(sm, tt)
    eval_conv!(sm, tt)
    eval_∇!(sm, tt)
end

"""
    eval_eos!(sm::StellarModel)

Evaluates the equation of state of the current stellar model
"""
function eval_eos!(sm::StellarModel, eval_type::TT) where {TT<:Real}
    species_names = sm.varnames[(sm.nvars - sm.nspecies + 1):end]
    Threads.@threads for k = 1:(sm.nz)
        _set_diff_cache!(sm, k)
        var00 = get_tmp(sm.diff_caches[k, 2], eval_type)
        sm.eos_results[k, :] .= get_EOS_resultsTP(sm.eos, sm.isotope_data, var00[sm.vari[:lnT]],
                                                  var00[sm.vari[:lnP]],
                                                  var00[(sm.nvars - sm.nspecies + 1):(sm.nvars)], species_names)
    end
end

"""
    eval_opacity!(sm::StellarModel)

Evaluates the opacity of the current stellar model
"""
function eval_opacity!(sm::StellarModel, eval_type::TT) where {TT<:Real}
    species_names = sm.varnames[(sm.nvars - sm.nspecies + 1):end]
    Threads.@threads for k = 1:(sm.nz)
        _set_diff_cache!(sm, k)
        var00 = get_tmp(sm.diff_caches[k, 2], eval_type)
        sm.opacity_results[k] = get_opacity_resultsTP(sm.opacity, sm.isotope_data, var00[sm.vari[:lnT]],
                                                      var00[sm.vari[:lnP]],
                                                      var00[(sm.nvars - sm.nspecies + 1):(sm.nvars)], species_names)
    end
end

"""
    eval_mlt!(sm::StellarModel)

Evaluates the mlt info of the current stellar model
"""
function eval_conv!(sm::StellarModel, eval_type::TT) where {TT<:Real}
    Threads.@threads for k = 1:(sm.nz)
        _set_diff_cache!(sm, k)
        sm.conv_results[k, :] .= get_conv_results(sm.convection, sm.opt.convection.alpha_mlt, sm.eos_results[k, 7])
    end
end

function eval_∇!(sm::StellarModel, eval_type::TT) where {TT<:Real}
    Threads.@threads for k = 1:(sm.nz - 1)
        _set_diff_cache!(sm, k)
        var00 = get_tmp(sm.diff_caches[k, 2], eval_type)
        varp1 = get_tmp(sm.diff_caches[k, 3], eval_type)
        sm.∇[k] = get_∇ᵣ(sm, k, sm.opacity_results, var00, varp1)  # eval'd at face
        α, β = get_face_weights(sm, k)
        ∇ₐface = α * sm.eos_results[k, 7] + β * sm.eos_results[k + 1, 7]
        if ∇ₐface <= sm.∇[k]  # check if we're convective
            sm.∇[k] = sm.conv_results[k, 1]
        end
    end
end

function get_∇ᵣ(sm::StellarModel, k::Int, κ::AbstractVector{TT}, var00::AbstractVector{TT},
                varp1::AbstractVector{TT}) where {TT<:Real}
    α, β = get_face_weights(sm, k)
    κface = α * κ[k] + β * κ[k + 1]
    Tface = exp(α * var00[sm.vari[:lnT]] + β * varp1[sm.vari[:lnT]])
    Pface = exp(α * var00[sm.vari[:lnP]] + β * varp1[sm.vari[:lnP]])
    return 3κface * var00[sm.vari[:lum]] * LSUN * Pface / (16π * CRAD * CLIGHT * CGRAV * sm.m[k] * Tface^4)
end

function get_face_weights(sm::StellarModel, k::Int)
    if k == sm.nz
        print("cannot take face weight at outer boundary")
        return
    end
    α = sm.dm[k] / (sm.dm[k] + sm.dm[k + 1])
    β = 1 - α
    return α, β
end

function _set_diff_cache!(sm::StellarModel, k::Int)
    if k == 1
        for i = 1:(sm.nvars)
            sm.diff_caches[1, 2].du[i] = sm.ind_vars[sm.nvars * (k - 1) + i]
            sm.diff_caches[1, 3].du[i] = sm.ind_vars[sm.nvars * (k - 1) + sm.nvars + i]
            sm.diff_caches[1, 2].dual_du[(i - 1) * (3 * sm.nvars + 1) + 1] = sm.ind_vars[sm.nvars * (k - 1) + i]
            sm.diff_caches[1, 3].dual_du[(i - 1) * (3 * sm.nvars + 1) + 1] = sm.ind_vars[sm.nvars * (k - 1) + sm.nvars + i]
        end
    elseif k == sm.nz
        for i = 1:(sm.nvars)
            sm.diff_caches[sm.nz, 1].du[i] = sm.ind_vars[sm.nvars * (k - 2) + i]
            sm.diff_caches[sm.nz, 2].du[i] = sm.ind_vars[sm.nvars * (k - 2) + sm.nvars + i]
            sm.diff_caches[sm.nz, 1].dual_du[(i - 1) * (3 * sm.nvars + 1) + 1] = sm.ind_vars[sm.nvars * (k - 2) + i]
            sm.diff_caches[sm.nz, 2].dual_du[(i - 1) * (3 * sm.nvars + 1) + 1] = sm.ind_vars[sm.nvars * (k - 2) + sm.nvars + i]
        end
    else
        for i = 1:(sm.nvars)
            sm.diff_caches[k, 1].du[i] = sm.ind_vars[sm.nvars * (k - 2) + i]
            sm.diff_caches[k, 2].du[i] = sm.ind_vars[sm.nvars * (k - 2) + sm.nvars + i]
            sm.diff_caches[k, 3].du[i] = sm.ind_vars[sm.nvars * (k - 2) + 2 * sm.nvars + i]
            sm.diff_caches[k, 1].dual_du[(i - 1) * (3 * sm.nvars + 1) + 1] = sm.ind_vars[sm.nvars * (k - 2) + i]
            sm.diff_caches[k, 2].dual_du[(i - 1) * (3 * sm.nvars + 1) + 1] = sm.ind_vars[sm.nvars * (k - 2) + sm.nvars + i]
            sm.diff_caches[k, 3].dual_du[(i - 1) * (3 * sm.nvars + 1) + 1] = sm.ind_vars[sm.nvars * (k - 2) + 2 * sm.nvars + i]
        end
    end
end
