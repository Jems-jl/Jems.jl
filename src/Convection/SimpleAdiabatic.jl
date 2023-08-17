export get_conv_results


struct AdiabaticConvection <: AbstractConvection
    num_results::Int
    AdiabaticConvection() = new(1)
end


function get_conv_results(mlt::AdiabaticConvection, alpha_mlt, ∇ₐ)
    return ∇ₐ
end