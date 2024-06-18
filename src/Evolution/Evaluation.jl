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
        m.solver_data.eqs_duals[k, m.nvars - m.network.nspecies + i] = m.composition_equation.func(m, k,
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
        eqs_dual = eqs_duals[k, i]
        partials = eqs_dual.partials
        if k == 1
            @inbounds @views jacobian_Lk[i, :] .= 0
            @inbounds @views jacobian_Dk[i, :] .= partials[(m.nvars+1):(2*m.nvars)]
            @inbounds @views jacobian_Uk[i, :] .= partials[(2 * m.nvars + 1):(3*m.nvars)]
        elseif k == m.props.nz
            @inbounds @views jacobian_Lk[i, :] .= partials[1:m.nvars]
            @inbounds @views jacobian_Dk[i, :] .= partials[(m.nvars+1):(2*m.nvars)]
            @inbounds @views jacobian_Uk[i, :] .= 0
        else
            @inbounds @views jacobian_Lk[i, :] .= partials[1:m.nvars]
            @inbounds @views jacobian_Dk[i, :] .= partials[(m.nvars+1):(2*m.nvars)]
            @inbounds @views jacobian_Uk[i, :] .= partials[(2 * m.nvars + 1):(3*m.nvars)]
        end
        # populate the eqs_numbers with relevant entries (will be RHS for linear solver)
        m.solver_data.eqs_numbers[(k - 1) * m.nvars + i] = eqs_dual.value
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
