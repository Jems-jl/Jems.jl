using LinearAlgebra

##
# Thomas algorithm, taken from
# BCYCLIC: A parallel block tridiagonal matrix cyclic solver
# Hishman et al. 2010
# This is not BCYCLIC, but their description of the Thomas
# algorithm for block tridiagonal systems. Beware there are typos
# in their equations (3a) and (3b)
function thomas_algorithm!(sm)
    # We store Δ_i in the diagonal
    b_1 = @view sm.eqs_numbers[1:sm.nvars]
    sm.solver_β[1] .= .- b_1 # this has an inverse sign for b_1 for our case
    # Δ_0 is just the inverse of the first diagonal block
    D_1 = sm.jacobian_D[1]
    D_1 .= inv(D_1)
    for i=2:sm.nz
        # update β first
        L_i = sm.jacobian_L[i]
        Δ_im1 = sm.jacobian_D[i-1] # we have stored Δ for previous zone here
        b_i = @view sm.eqs_numbers[(i-1)*sm.nvars+1:i*sm.nvars]
        β_i = sm.solver_β[i]
        β_im1 = sm.solver_β[i-1]
        mul!(sm.solver_x[i], Δ_im1, β_im1) #sm.solver_x[i] is only used as placeholder here
        mul!(β_i, L_i, sm.solver_x[i]) # here also β_i is a placeholder
        for j in 1:sm.nvars
            β_i[j] = .-b_i[j] - β_i[j]
        end

        # update Δ next
        tmp = sm.jacobian_tmp[1]
        D_i = sm.jacobian_D[i]
        U_im1 = sm.jacobian_U[i-1]
        mul!(tmp, L_i, Δ_im1)
        mul!(L_i, tmp, U_im1)
        D_i .= D_i .- L_i
        D_i .= inv(D_i)
    end

    x_N = sm.solver_x[sm.nz]
    Δ_N = sm.jacobian_D[sm.nz]
    β_N = sm.solver_β[sm.nz] 
    mul!(x_N,Δ_N,β_N)
    # backwards sweep
    for i=sm.nz-1:-1:1
        x_i = sm.solver_x[i]
        x_ip1 = sm.solver_x[i+1]
        β_i = sm.solver_β[i]
        U_i = sm.jacobian_U[i]
        Δ_i = sm.jacobian_D[i]

        mul!(x_i,U_i,x_ip1)
        β_i .= β_i .- x_i
        mul!(x_i,Δ_i,β_i)
    end
    #unload result into solver_corr
    for i=1:sm.nz
        for j=1:sm.nvars
            sm.solver_corr[(i-1)*sm.nvars+j] = sm.solver_x[i][j]
        end
    end
    return
end