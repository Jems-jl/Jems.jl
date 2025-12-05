struct TRhoProfile{TOBS,TOBS2,TAXIS}
    logRho::TOBS
    logT::TOBS
    mixing::TOBS2
    axis::TAXIS
end

function _read_trho_data(file)
    ρs::Vector{Float64} = []
    Ts::Vector{Float64} = []
    open(file) do io
        for line in eachline(io)
            ρ, T = split(line)
            push!(ρs, parse(Float64, ρ))
            push!(Ts, parse(Float64, T))
        end
    end
    return ρs, Ts
end

function TRhoProfile(grid_pos)
    axis = Axis(grid_pos, xlabel=L"\log_{10}(T/[\text{K}])", ylabel=L"\log_{10}(\rho/[\text{g\,cm^{-3}}])")

    h_burn_ρ, h_burn_T = _read_trho_data(pkgdir(Plotting, "data/PlotData", "hydrogen_burn.data"))
    lines!(axis, h_burn_ρ, h_burn_T, color=:gray, linestyle=:dash, linewidth=2)
    he_burn_ρ, he_burn_T = _read_trho_data(pkgdir(Plotting, "data/PlotData", "helium_burn.data"))
    lines!(axis, he_burn_ρ, he_burn_T, color=:gray, linestyle=:dash, linewidth=2)
    e_degen_ρ, e_degen_T = _read_trho_data(pkgdir(Plotting, "data/PlotData", "psi4.data"))
    lines!(axis, e_degen_ρ, e_degen_T, color=:gray, linestyle=:dash, linewidth=2)

    pgas_ρ = [-8, 5]
    pgas_T = log10(3.2e7) .+ (pgas_ρ .- log10(0.7e0))./3.0
    lines!(axis, pgas_ρ, pgas_T, color=:gray, linestyle=:dash, linewidth=2)

    mixing = Observable(mixing_colors[get.(Ref(mixing_map), [:convection], RGBAf(0.1, 0.1, 0.1))])
    logRho = Observable(zeros(Float64,1))
    logT = Observable(zeros(Float64,1))
    lines!(axis, logRho, logT, color=mixing)

    nom = LineElement(color=mixing_colors[mixing_map[:no_mixing]])
    conv = LineElement(color=mixing_colors[mixing_map[:convection]])
    axislegend(axis, [nom, conv], ["no mixing", "convection"],position=:lt)
    return TRhoProfile(logRho, logT, mixing, axis)
end

function update_plot!(p::TRhoProfile, m)
    p.logRho[] = get_value.(m.props.lnρ[1:(m.props.nz)]) .* log10(ℯ)
    p.logT[] = get_value.(m.props.lnT[1:(m.props.nz)]) .* log10(ℯ)
    p.mixing[] = mixing_colors[get.(Ref(mixing_map), (@view m.props.mixing_type[1:(m.props.nz)]), RGBAf(0.1, 0.1, 0.1))]
end
