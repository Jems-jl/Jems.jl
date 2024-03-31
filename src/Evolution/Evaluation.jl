"""
    eval_cell_eqs(m::StellarModel, k::Int, ind_vars_view::Vector{<:TT}) where{TT<:Real}

Evaluates the stellar structure equations of the model, `m`, at cell `k`, given the view of the independent
variables, `ind_vars_view`.
"""
function eval_cell_eqs!(m::AbstractModel, k::Int)
    # evaluate all equations! (except composition)
    for i = 1:(m.nvars - m.network.nspecies)
        m.solver_data.eqs_duals[k, i] = m.structure_equations[i].func(m, k)
    end
    # evaluate all composition equations
    for i = 1:(m.network.nspecies)
        m.solver_data.eqs_duals[k, m.nvars - m.network.nspecies + i] = equation_composition(m, k,
                                                                                            m.network.species_names[i])
    end
end

"""
    eval_jacobian_eqs_row!(m::AbstractModel, k::int)

Evaluates row `k` of the Jacobian matrix of the given Model `m`.
"""
function eval_jacobian_eqs_row!(m::AbstractModel, k::Int)
    #=
    Jacobian has tridiagonal block structure:
    x x - - - -
    x x x - - -
    - x x x - -
    - - x x x -
    - - - x x x
    - - - - x x
    In the first row, the first block corresponds to the derivatives of the structure equations with respect to the
    variables at the first cell. The block to the right corresponds to the derivatives with respect to the variables
    at the second cell. In the second row until row nz-1, the first block contains the derivatives of the structure
    equations wrt. the variables of the previous cell, the second block wrt. to the varibles in the current cell, and
    the third wrt. the variables of the next cell. The last row only contains the derivatives wrt. the penultimate
    cell and the last cell.

    Because this function acts on duals `eqs_duals`, we immediately evaluate the equations also.
    =#
    eval_cell_eqs!(m, k)  # evaluate equations on the duals, so we have jacobian also
    # populate the jacobian with the relevant entries
    jacobian_Lk = m.solver_data.jacobian_L[k]
    jacobian_Dk = m.solver_data.jacobian_D[k]
    jacobian_Uk = m.solver_data.jacobian_U[k]
    eqs_duals = m.solver_data.eqs_duals
    for i = 1:(m.nvars)
        if k == 1
            for j = 1:(m.nvars)
                jacobian_Lk[i, j] = 0
                jacobian_Dk[i, j] = eqs_duals[k, i].partials[j + m.nvars]
                jacobian_Uk[i, j] = eqs_duals[k, i].partials[j + 2 * m.nvars]
            end
        elseif k == m.props.nz
            for j = 1:(m.nvars)
                jacobian_Lk[i, j] = eqs_duals[k, i].partials[j]
                jacobian_Dk[i, j] = eqs_duals[k, i].partials[j + m.nvars]
                jacobian_Uk[i, j] = 0
            end
        else
            for j = 1:(m.nvars)
                jacobian_Lk[i, j] = eqs_duals[k, i].partials[j]
                jacobian_Dk[i, j] = eqs_duals[k, i].partials[j + m.nvars]
                jacobian_Uk[i, j] = eqs_duals[k, i].partials[j + 2 * m.nvars]
            end
        end
        # populate the eqs_numbers with relevant entries (will be RHS for linear solver)
        m.solver_data.eqs_numbers[(k - 1) * m.nvars + i] = eqs_duals[k, i].value
    end
end

"""
    eval_jacobian_eqs!(m::StellarModel)

Evaluates the whole Jacobian matrix and equations of the given StellarModel `m`.
"""
function eval_jacobian_eqs!(m::AbstractModel)
    Threads.@threads for k = 1:(m.props.nz)
        eval_jacobian_eqs_row!(m, k)
    end
end
