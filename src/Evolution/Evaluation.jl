"""
    eval_cell_eqs(sm::StellarModel, k::Int, ind_vars_view::Vector{<:TT}) where{TT<:Real}

Evaluates the stellar structure equations of the stellar model, `sm`, at cell `k`, given the view of the independent
variables, `ind_vars_view`.
"""
function eval_cell_eqs!(sm::StellarModel, k::Int, ind_vars_view::Vector{<:TT},
                        result_view::Vector{<:TT}) where {TT<:Real}
    # initialize as undefined, since m1 and p1 values are not defined at edges
    eosm1 = Vector{TT}(undef, sm.eos.num_results)
    eos00 = Vector{TT}(undef, sm.eos.num_results)
    eosp1 = Vector{TT}(undef, sm.eos.num_results)
    κm1::TT = NaN
    κ00::TT = NaN
    κp1::TT = NaN
    varm1::Vector{TT} = []
    var00::Vector{TT} = []
    varp1::Vector{TT} = []

    # collect required variables
    if k > 1 && k < sm.nz
        varm1 = ind_vars_view[1:(sm.nvars)]
        var00 = ind_vars_view[(sm.nvars + 1):(2 * sm.nvars)]
        varp1 = ind_vars_view[(2 * sm.nvars + 1):(3 * sm.nvars)]
    elseif k == 1
        var00 = ind_vars_view[1:(sm.nvars)]
        varp1 = ind_vars_view[(sm.nvars + 1):(2 * sm.nvars)]
    elseif k == sm.nz
        varm1 = ind_vars_view[1:(sm.nvars)]
        var00 = ind_vars_view[(sm.nvars + 1):(2 * sm.nvars)]
    end

    # collect eos and κ info
    eos00 = view(sm.eos_results, k, :)
    κ00 = sm.opacity_results[k]
    if k != 1
        eosm1 = view(sm.eos_results, k - 1, :)
        κm1 = sm.opacity_results[k - 1]
    end
    if k != sm.nz
        eosp1 = view(sm.eos_results, k + 1, :)
        κp1 = sm.opacity_results[k + 1]
    end

    # evaluate all equations!
    for i = 1:(sm.nvars)
        result_view[i] = sm.structure_equations[i](sm, k, varm1, var00, varp1, eosm1, eos00, eosp1, κm1, κ00, κp1)
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
        resultview = view(sm.eqs, (sm.nvars * (k - 1) + 1):(sm.nvars * k))
        eval_cell_eqs!(sm, k, view(sm.ind_vars, ki:kf), resultview)
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
    resultview = view(sm.eqs, (sm.nvars * (k - 1) + 1):(sm.nvars * k))
    eval_cell_eqs!(sm, k, view(sm.ind_vars, ki:kf), resultview)
    jac_view = view(sm.jacobian, (sm.nvars * (k - 1) + 1):(sm.nvars * k), ki:kf)
    for j=1:sm.nvars
        for i=1:sm.nvars
            jac_view[j, i] = resultview[j].partials[i]
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
    eval_eos!(sm)
    eval_opacity!(sm)
    eval_conv!(sm)
    eval_∇!(sm)
end

"""
    eval_eos!(sm::StellarModel)

Evaluates the equation of state of the current stellar model
"""
function eval_eos!(sm::StellarModel)
    species_names = sm.varnames[(sm.nvars - sm.nspecies + 1):end]
    Threads.@threads for k = 1:(sm.nz)
        var_here = sm.ind_vars[((k - 1) * sm.nvars + 1):(k * sm.nvars)]
        sm.eos_results[k, :] .= get_EOS_resultsTP(sm.eos, sm.isotope_data, var_here[sm.vari[:lnT]],
                                                  var_here[sm.vari[:lnP]],
                                                  var_here[(sm.nvars - sm.nspecies + 1):(sm.nvars)], species_names)
    end
end

"""
    eval_opacity!(sm::StellarModel)

Evaluates the opacity of the current stellar model
"""
function eval_opacity!(sm::StellarModel)
    species_names = sm.varnames[(sm.nvars - sm.nspecies + 1):end]
    Threads.@threads for k = 1:(sm.nz)
        var_here = sm.ind_vars[((k - 1) * sm.nvars + 1):(k * sm.nvars)]
        sm.opacity_results[k] = get_opacity_resultsTP(sm.opacity, sm.isotope_data, var_here[sm.vari[:lnT]],
                                                      var_here[sm.vari[:lnP]],
                                                      var_here[(sm.nvars - sm.nspecies + 1):(sm.nvars)], species_names)
    end
end

"""
    eval_mlt!(sm::StellarModel)

Evaluates the mlt info of the current stellar model
"""
function eval_conv!(sm::StellarModel)
    Threads.@threads for k = 1:(sm.nz)
        sm.conv_results[k, :] .= get_conv_results(sm.convection, sm.opt.convection.alpha_mlt, sm.eos_results[k, 7])
    end
end

function eval_∇!(sm::StellarModel)
    Threads.@threads for k = 1:(sm.nz - 1)
        sm.∇[k] = get_∇ᵣ(sm, k)  # eval'd at face
        α, β = get_face_weights(sm, k)
        ∇ₐface = α * sm.eos_results[k, 7] + β * sm.eos_results[k + 1, 7]
        if ∇ₐface <= sm.∇[k]  # check if we're convective
            sm.∇[k] = sm.conv_results[k, 1]
        end
    end
end

function get_∇ᵣ(sm::StellarModel, k::Int)
    α, β = get_face_weights(sm, k)
    κface = α * sm.opacity_results[k] + β * sm.opacity_results[k + 1]
    Tface = exp(α * sm.csi.lnT[k] + β * sm.csi.lnT[k + 1])
    Pface = exp(α * sm.csi.lnP[k] + β * sm.csi.lnT[k + 1])
    L = sm.csi.L[k] * LSUN
    return 3κface * L * Pface / (16π * CRAD * CLIGHT * CGRAV * sm.m[k] * Tface^4)
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
