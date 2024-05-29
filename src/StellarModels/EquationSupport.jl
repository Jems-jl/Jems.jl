using FunctionWrappers

"""
    struct TypeStableStructureEquation{TM,TD<:Real}

Structure that wraps a stellar structure equation into a type stable object, using FunctionWrappers.jl. This requires
that the stellar structure equations have the following signature:

    ```
    function structure_equation(::TM, ::Int)::TD
    ```

For typical usage, TM is the concrete type of a Model (say, a StellarModel), and TD the type of dual number being used
for automatic differentiation. The function must return an object of type TD, the result of the equation.
"""
struct TypeStableStructureEquation{TM,TD<:Real}
    func::FunctionWrappers.FunctionWrapper{TD,Tuple{TM,Int}}
end

"""
    struct TypeStableCompositionEquation{TM,TD<:Real}

Behaves similarly to TypeStableStructureEquation, but since this equation has a different calling signature, we need a 
new type.
"""
struct TypeStableCompositionEquation{TM,TD<:Real}
    func::FunctionWrappers.FunctionWrapper{TD,Tuple{TM,Int,Symbol}}
end
