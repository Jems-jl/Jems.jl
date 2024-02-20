using LinearAlgebra

##
# Thomas algorithm, taken from
# BCYCLIC: A parallel block tridiagonal matrix cyclic solver
# Hishman et al. 2010
# This is not BCYCLIC, but their description of the Thomas
# algorithm for block tridiagonal systems. Beware there are typos
# in their equations (3a) and (3b)
function thomas_algorithm!(sm)
    eqs_numbers = sm.solver_data.eqs_numbers
    solver_LU = sm.solver_data.solver_LU
    jacobian_D = sm.solver_data.jacobian_D
    jacobian_L = sm.solver_data.jacobian_L
    jacobian_U = sm.solver_data.jacobian_U
    solver_tmp1 = sm.solver_data.solver_tmp1
    solver_tmp2 = sm.solver_data.solver_tmp2
    solver_x = sm.solver_data.solver_x
    solver_β = sm.solver_data.solver_β
    solver_corr = sm.solver_data.solver_corr

    ##crappy pre-conditioning
    #for i in 1:sm.nz
    #    max_L = maximum(abs, jacobian_L[i], dims=2)
    #    max_D = maximum(abs, jacobian_D[i], dims=2)
    #    max_U = maximum(abs, jacobian_U[i], dims=2)
    #    if i ==1
    #        max_vals = max.(max_D, max_U)
    #    elseif i==sm.nz
    #        max_vals = max.(max_L, max_D)
    #    else
    #        max_vals = max.(max_L, max_D, max_U)
    #    end
    #    for j in 1:sm.nvars
    #        for l in 1:sm.nvars
    #            sm.solver_data.jacobian_L[i][j,l] /= max_vals[j]
    #            sm.solver_data.jacobian_D[i][j,l] /= max_vals[j]
    #            sm.solver_data.jacobian_U[i][j,l] /= max_vals[j]
    #        end
    #    end
    #    sm.solver_data.eqs_numbers[(sm.nvars*(i-1)+1):(sm.nvars*(i-1)+sm.nvars)] ./= max_vals
    #end

    # We store Δ_i in the diagonal
    b_1 = @view eqs_numbers[1:sm.nvars]
    solver_β[1] .= .- b_1  # this has an inverse sign for b_1 for our case
    # Δ_0 is just the inverse of the first diagonal block
    D_1 = jacobian_D[1]
    solver_LU[1] = lu!(D_1)
    for i=2:sm.props.nz
        # update β first
        L_i = jacobian_L[i]
        LU_Δ_im1 = solver_LU[i-1] # we have stored the LU factorization of Δ for previous zone here
        b_i = @view eqs_numbers[(i-1)*sm.nvars+1:i*sm.nvars]
        β_i = solver_β[i]
        β_im1 = solver_β[i-1]
        ldiv!(solver_x[i], LU_Δ_im1, β_im1) #sm.solver_x[i] is only used as placeholder here
        mul!(β_i, L_i, solver_x[i]) # here also β_i is a placeholder
        for j in 1:sm.nvars
            β_i[j] = -b_i[j] - β_i[j]
        end

        # update Δ next
        D_i = jacobian_D[i]
        U_im1 = jacobian_U[i-1]
        ldiv!(solver_tmp1[i], LU_Δ_im1, U_im1)
        mul!(solver_tmp2[i], L_i, solver_tmp1[i])
        D_i .= D_i .- solver_tmp2[i]
        solver_LU[i] = lu!(D_i)
    end

    x_N = solver_x[sm.props.nz]
    LU_Δ_N = solver_LU[sm.props.nz]
    β_N = solver_β[sm.props.nz] 
    ldiv!(x_N, LU_Δ_N, β_N)
    # backwards sweep
    for i=sm.props.nz-1:-1:1
        x_i = solver_x[i]
        x_ip1 = solver_x[i+1]
        β_i = solver_β[i]
        U_i = jacobian_U[i]
        LU_Δ_i = solver_LU[i]

        mul!(x_i,U_i,x_ip1)
        β_i .= β_i .- x_i
        ldiv!(x_i, LU_Δ_i, β_i)
    end
    # unload result into solver_corr
    for i=1:sm.props.nz
        for j=1:sm.nvars
            solver_corr[(i-1)*sm.nvars+j] = solver_x[i][j]
        end
    end
    return
end