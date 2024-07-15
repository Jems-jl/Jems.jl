module Rotation
using ..DualSupport

function get_i_rot(sm, k::Int)
    return 2 / 3 * exp(get_cell_dual(sm.props.lnr[k]))^2
end

function get_ω(sm, k::Int)
    return get_cell_dual(sm.props.j_rot[k]) / get_cell_dual(sm.props.i_rot[k])
end

function get_ν_ω(sm, k::Int)
    x::typeof(sm.props).parameters[2] = 1.0  # type of numbers used in the jac.
    return x  # cm^2 s^-1
end

end  # module Rotation