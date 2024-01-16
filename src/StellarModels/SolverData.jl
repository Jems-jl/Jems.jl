using StaticArrays

@kwdef mutable struct SolverData{TNUMBER, TDUALFULL, TMATRIX, TLU, TVECTOR}
    eqs_numbers::Vector{TNUMBER}  # Stores the results of the equation evaluations (as numbers), size nz * nvars
    eqs_duals::Matrix{TDUALFULL}  # Stores the dual results of the equation evaluation, shape (nz, nvars)
    jacobian_D::Vector{TMATRIX}
    jacobian_U::Vector{TMATRIX}
    jacobian_L::Vector{TMATRIX}
    solver_LU::Vector{TLU}
    solver_tmp1::Vector{TMATRIX}
    solver_tmp2::Vector{TMATRIX}
    solver_β::Vector{TVECTOR}
    solver_x::Vector{TVECTOR}
    solver_corr::Vector{TNUMBER}
    newton_iters::Int
end

function SolverData(nvars, nz, nextra, use_static_arrays, number_type)
    # create the equation results matrix, holding dual numbers (for automatic differentiation, AD)
    dual_sample = ForwardDiff.Dual(zero(number_type), (zeros(number_type, 3 * nvars)...))
    eqs_duals = Matrix{typeof(dual_sample)}(undef, nz+nextra, nvars)
    for k = 1:(nz + nextra)
        for i = 1:nvars
            eqs_duals[k, i] = ForwardDiff.Dual(zero(number_type), (zeros(number_type, 3 * nvars)...))
        end
    end

    # create jacobian matrix (we have the diagonal and the upper and lower blocks)
    # we use static arrays, provided by StaticArrays. These are faster than regular
    # arrays for small nvars
    if use_static_arrays
        jacobian_D = [(@MMatrix zeros(number_type, nvars,nvars)) for i=1:(nz+nextra)]
        jacobian_U = [(@MMatrix zeros(number_type, nvars,nvars)) for i=1:(nz+nextra)]
        jacobian_L = [(@MMatrix zeros(number_type, nvars,nvars)) for i=1:(nz+nextra)]
        solver_LU = Vector{LU{number_type, MMatrix{nvars, nvars, number_type, nvars^2}, Vector{Int64}}}(undef,nz+nextra)
        solver_tmp1 = [(@MMatrix zeros(number_type,nvars,nvars)) for i=1:(nz+nextra)]
        solver_tmp2 = [(@MMatrix zeros(number_type,nvars,nvars)) for i=1:(nz+nextra)]
        solver_β = [(@MVector zeros(number_type, nvars)) for i=1:(nz+nextra)]
        solver_x = [(@MVector zeros(number_type, nvars)) for i=1:(nz+nextra)]
    else
        jacobian_D = [(zeros(number_type, nvars,nvars)) for i=1:(nz+nextra)]
        jacobian_U = [(zeros(number_type, nvars,nvars)) for i=1:(nz+nextra)]
        jacobian_L = [(zeros(number_type, nvars,nvars)) for i=1:(nz+nextra)]
        solver_LU = Vector{LU{number_type, Matrix{number_type}, Vector{Int64}}}(undef,nz+nextra)
        solver_tmp1 = [(zeros(number_type,nvars,nvars)) for i=1:(nz+nextra)]
        solver_tmp2 = [(zeros(number_type,nvars,nvars)) for i=1:(nz+nextra)]
        solver_β = [(zeros(number_type,nvars)) for i=1:(nz+nextra)]
        solver_x = [(zeros(number_type,nvars)) for i=1:(nz+nextra)]
        #When using regulars arrays one gets better performance by setting the BLAS
        #threads to one like this:
        #using LinearAlgebra
        #LinearAlgebra.BLAS.set_num_threads(1) # this allows for a faster linear solve
    end
    solver_corr = zeros(number_type, nvars*(nz+nextra))

    # create the equation results vector for the solver (holds plain numbers instead of duals)
    eqs_numbers = ones(number_type, nvars * (nz+nextra))

    SolverData(eqs_numbers = eqs_numbers,
               eqs_duals = eqs_duals,
               jacobian_D = jacobian_D,
               jacobian_U = jacobian_U,
               jacobian_L = jacobian_L,
               solver_LU = solver_LU,
               solver_tmp1 = solver_tmp1,
               solver_tmp2 = solver_tmp2,
               solver_β = solver_β,
               solver_x = solver_x,
               solver_corr = solver_corr,
               newton_iters = 0)
end
