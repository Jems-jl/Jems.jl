using DataFrames
using ForwardDiff, DataInterpolations
import ForwardDiff.Dual
using Jems.StellarModels
using Jems.DualExtrapolation; const de = DualExtrapolation
using Interpolations, Dierckx
using HDF5
using Jems.Constants, Statistics
using DataInterpolations: CubicSpline
using CairoMakie, LaTeXStrings, MathTeXEngine, Makie.Colors, PlotUtils
basic_theme = Theme(fonts=(regular=texfont(:text), bold=texfont(:bold),
                           italic=texfont(:italic), bold_italic=texfont(:bolditalic)),
                    fontsize=30, size=(1000, 750), linewidth=7,
                    Axis=(xlabelsize=40, ylabelsize=40, titlesize=40, xgridvisible=false, ygridvisible=false,
                          spinewidth=2.5, xminorticksvisible=true, yminorticksvisible=true, xtickalign=1, ytickalign=1,
                          xminortickalign=1, yminortickalign=1, xticksize=14, xtickwidth=2.5, yticksize=14,
                          ytickwidth=2.5, xminorticksize=7, xminortickwidth=2.5, yminorticksize=7, yminortickwidth=2.5,
                          xticklabelsize=35, yticklabelsize=35, xticksmirrored=true, yticksmirrored=true),
                    Legend=(patchsize=(70, 10), framevisible=false, patchlabelgap=20, rowgap=10))
set_theme!(basic_theme)
struct Difference_with_JEMS
    logM
    diff_logL
    diff_logT
    diff_pyth
    diff_tot
    zeta
end
function Difference_with_JEMS(interpolTrack::de.InterpolTrack, modeltrack::de.Track)
    if interpolTrack.logM != modeltrack.logM
       @warn "logM of modeltrack (=$modeltrack.logM) and interpoltrack (=$interpolTrack.logM) not the same"
    end
    if length(interpolTrack.zeta) != length(modeltrack.zeta)
        throw(ArgumentError("Not the same zeta sampling"))
    end
    diff_logL = modeltrack.logL_val - interpolTrack.logL_val; logL_norm = maximum(abs.(diff_logL)) 
    diff_logT = modeltrack.logT_val - interpolTrack.logT_val; logT_norm = maximum(abs.(diff_logT))
    diff_pyth = sqrt.((diff_logL/logL_norm).^2 .+ (diff_logT/logT_norm).^2); 
    diff_tot = sum(diff_pyth)
    weight_L = 1 / abs(modeltrack.logL_val[end] - modeltrack.logL_val[1])
    weight_T = 1 / abs(modeltrack.logT_val[end] - modeltrack.logT_val[1])
    factor_L = weight_L / (weight_L + weight_T); factor_T = weight_T / (weight_L + weight_T)
    diff_pyth = sqrt.((diff_logL*factor_L).^2 .+ (diff_logT*factor_T).^2); diff_tot = sum(diff_pyth)
    return Difference_with_JEMS(interpolTrack.logM, diff_logL, diff_logT, diff_pyth, diff_tot, interpolTrack.zeta)
end
##
gridpath = "DualRuns/DualGrid"
gridpath = "DualRuns/DualGrid2"
gridpath = "DualRuns/DualGrid3"

path = "DualRuns/DualGrid/logM_-0.1_X_0.7381_.history.hdf5"
get_logM(path) = parse(Float64, split(split(path, "logM_")[2], "_")[1])

#get ALL filepaths in gridpath
historypaths = filter(x -> occursin(".history.hdf5", x), readdir(gridpath))
profilepaths = filter(x -> occursin(".profiles.hdf5", x), readdir(gridpath))
N = length(historypaths)
println("$N history files found")
models = Dict()
modeltracks = Dict()
X_dual         = ForwardDiff.Dual{}(0.7381,  0.0,1.0,0.0,0.0,0.0)
Z_dual         = ForwardDiff.Dual{}(0.0134,  0.0,0.0,1.0,0.0,0.0)
Dfraction_dual = ForwardDiff.Dual{}(0.000312,0.0,0.0,0.0,1.0,0.0)
R_dual         = ForwardDiff.Dual{}(100*RSUN,0.0,0.0,0.0,0.0,1.0)

inititial_params_names = [:logM, :X, :Z, :Dfraction, :R]
all_logMs = []
i=8
only_use_these_logMs = [0.0,0.02]
only_use_these_logMs = []
for i in 1:N
    historypath = joinpath(gridpath, historypaths[i])
    profilepath = joinpath(gridpath, profilepaths[i])
    @assert get_logM(historypath) == get_logM(profilepath)
    logM = get_logM(historypath)
    if !isempty(only_use_these_logMs) if logM ∉ only_use_these_logMs; continue; end; end
    history_dual, profiles_dual = de.bookkeeping(historypath, profilepath,3)
    logM_dual = ForwardDiff.Dual{}(logM, 1.0,0.0,0.0,0.0,0.0)
    initial_params = [logM_dual, X_dual, Z_dual, Dfraction_dual, R_dual]
    model = de.Model_constructor(history_dual, profiles_dual, initial_params, inititial_params_names)
    models[logM] = model
    track = nothing
    try 
        X_init = 0.4; X_end = 0.0000001
        X_init = 0.99*model.history_value.X_center[1]; X_end = 0.01
        track = de.Track(model,X_init, X_end, 1000)
    catch 
        println(" Track FAILED for logM = $logM")
        continue
    end
    println(" Track OK for logM = $logM")
    modeltracks[logM] = track
    push!(all_logMs, logM)
end
##
interpol_asked = [0.02,0.03,0.04]
interpol_asked = [0.05,0.1,0.15]
interpolTrackxs = Dict()
modeltrack1, modeltrack2 = modeltracks[0.0], modeltracks[0.2]
extrapolGrid1 = de.ExtrapolGrid(modeltrack1, interpol_asked .- modeltrack1.logM);
extrapolGrid2 = de.ExtrapolGrid(modeltrack2, interpol_asked .- modeltrack2.logM);

for logM in interpol_asked
    interpolTrack = de.InterpolTrack(modeltrack1, modeltrack2,logM);
    interpolTrackxs[logM] = interpolTrack
end
## EERSTE PLOT INTERPOLATIE IN THESIS
fig = Figure(size=(1500,1200))
ax = Axis(fig[1:2,1], xlabel=L"$\log (T_{\text{eff}} / K)$", ylabel=L"$\log (L/L_\odot)$",xreversed=true)
down = de.plot!(modeltrack1, ax; color=:red, scatter=false,linestyle=:dash, linewidth=5)
up = de.plot!(modeltrack2, ax; color=:blue, scatter=false,linestyle=:dot, linewidth=5)
jems_compare = nothing; interpol = nothing; extrapol_down = nothing; extrapol_up = nothing
for logM in interpol_asked
    interpolTrack = interpolTrackxs[logM]
    interpol = de.plot!(interpolTrack, ax; color=:lightgreen, scatter=false,linewidth=5)
    jems_compare = de.plot!(modeltracks[logM], ax; color=:black, scatter=false,linewidth=2)
    text!(ax, interpolTrack.logT_val[1]+0.01,interpolTrack.logL_val[1]-0.01; text = string(logM),fontsize=40)
end
for extrapolTrack in extrapolGrid1.extrapoltracks
    extrapol_down = de.plot!(extrapolTrack, ax; color=:red, scatter=false,linewidth=2,linestyle = :dash)
end
for extrapolTrack in extrapolGrid2.extrapoltracks
    extrapol_up = de.plot!(extrapolTrack, ax; color=:blue, scatter=false,linewidth=2, linestyle = :dot)
end
axislegend(ax,[up, extrapol_up, jems_compare, interpol, extrapol_down, down],["JEMS log M = " * string(modeltrack2.logM),"Extrapolation from above","JEMS tracks for comparison", "Interpolated tracks","Extrapolation from below","JEMS log M = " * string(modeltrack1.logM)])
axup = Axis(fig[0,1],xreversed=true, ylabel=L"$\log (L/L_\odot)$") ; axup.xticklabelsvisible = true

interpolTrack = interpolTrackxs[0.1]
de.plot!(interpolTrack, axup; color=:lightgreen, scatter=false,linewidth=5)
text!(axup, interpolTrack.logT_val[1],interpolTrack.logL_val[1]+0.05; text = string(0.1),fontsize=50)
de.plot!(modeltracks[0.1], axup; color=:black, scatter=false,linewidth=2)
de.plot!(extrapolGrid1.extrapoltracks[2], axup; color=:red, scatter=false,linewidth=2,linestyle = :dash)
de.plot!(extrapolGrid2.extrapoltracks[2], axup; color=:blue, scatter=false,linewidth=2, linestyle = :dot)
fig
##
figpath = "Figures/DSE_interpolation_first_example.png"; save(figpath, fig, px_per_unit=3); @show figpath
##

interpol_asked = [0.25,0.3,0.35]
interpolTrackxs = Dict()
modeltrack1, modeltrack2 = modeltracks[0.2], modeltracks[0.4]
extrapolGrid1 = de.ExtrapolGrid(modeltrack1, interpol_asked .- modeltrack1.logM);
extrapolGrid2 = de.ExtrapolGrid(modeltrack2, interpol_asked .- modeltrack2.logM);

for logM in interpol_asked
    interpolTrack = de.InterpolTrack(modeltrack1, modeltrack2,logM);
    interpolTrackxs[logM] = interpolTrack
end

## PLOT TWEE FIGUREN, SOORT VAN KOPIE LINKS
fig = Figure(size=(2000,1000))
ax = Axis(fig[1,1], xlabel=L"$\log (T_{\text{eff}} / K)$", ylabel=L"$\log (L/L_\odot)$", title=L"Interpolation between $\log M = 0$ and $0.2",xreversed=true)
de.plot!(modeltrack1, ax; color=:red, scatter=false, label = "JEMS log M = " * string(modeltrack1.logM),linestyle=:dash, linewidth=5)
de.plot!(modeltrack2, ax; color=:blue, scatter=false,label = "JEMS log M = " * string(modeltrack2.logM),linestyle=:dot, linewidth=5)
for logM in interpol_asked
    interpolTrack = interpolTrackxs[logM]
    de.plot!(interpolTrack, ax; color=:lightgreen, scatter=true,linewidth=5)
    de.plot!(modeltracks[logM], ax; color=:black, scatter=false,linewidth=2)

end
for extrapolTrack in extrapolGrid1.extrapoltracks
    de.plot!(extrapolTrack, ax; color=:red, scatter=false,linewidth=2,linestyle = :dash)
end
for extrapolTrack in extrapolGrid2.extrapoltracks
    de.plot!(extrapolTrack, ax; color=:blue, scatter=false,linewidth=2, linestyle = :dot)
end
axislegend(ax,position = :lb)
ax2 = Axis(fig[1,2], xlabel=L"$\log (T_{\text{eff}} / K)$", title=L"Interpolation of $\log M = 0.05,0.10,0.15$",xreversed=true)
for logM in interpol_asked
    interpolTrack = interpolTrackxs[logM]
    de.plot!(interpolTrack, ax2; color=:lightgreen, scatter=false,label="Interpol",linewidth=10)
    de.plot!(modeltracks[logM], ax2; color=:black, scatter=false,label="JEMS",linewidth=2)
end
for extrapolTrack in extrapolGrid1.extrapoltracks
    de.plot!(extrapolTrack, ax2; color=:red, scatter=false,label="Extrapol from below",linewidth=2,linestyle = :dash)
end
for extrapolTrack in extrapolGrid2.extrapoltracks
    de.plot!(extrapolTrack, ax2; color=:blue, scatter=false,label="Extrapol from above",linewidth=2, linestyle = :dot)
end
leg = Legend(fig[1,3],ax2;tellwidth=true)
fig


## ################################"
interpol_asked = [-0.05,0.0,0.05,0.1,0.15,0.2,0.25]; interpolTrackxs = Dict(); diffs_with_JEMS = Dict()
modeltrack1, modeltrack2 = modeltracks[-0.1], modeltracks[0.3]
extrapolGrid1 = de.ExtrapolGrid(modeltrack1, interpol_asked .- modeltrack1.logM);
extrapolGrid2 = de.ExtrapolGrid(modeltrack2, interpol_asked .- modeltrack2.logM);
for logM in interpol_asked
    interpolTrack = de.InterpolTrack(modeltrack1, modeltrack2,logM);
    interpolTrackxs[logM] = interpolTrack; diffs_with_JEMS[logM] = Difference_with_JEMS(interpolTrack, modeltracks[logM])
end
##
fig = Figure(size=(2000,800)); logM1 = modeltrack1.logM; logM2 = modeltrack2.logM; title = "Interpolation between $logM1 and $logM2"
ax = Axis(fig[1,1], xlabel=L"$\log (T_{\text{eff}} / K)$", ylabel=L"$\log (L/L_\odot)$", title=title,xreversed=true)
de.plot!(modeltrack2, ax; color=:blue, scatter=false,label = "JEMS log M = " * string(modeltrack2.logM),linestyle=:dot, linewidth=5)
de.plot!(modeltrack1, ax; color=:red, scatter=false, label = "JEMS log M = " * string(modeltrack1.logM),linestyle=:dash, linewidth=5)
colormap = :viridis; nb_colors = length(interpol_asked); colors = palette(colormap); vals = collect(1:nb_colors).*255 ./ nb_colors
colors = [colors[round(Int,vals[i])] for i in 1:nb_colors]
for (logM,color) in zip(interpol_asked,colors)
    interpolTrack = interpolTrackxs[logM]
    de.plot!(interpolTrack, ax; color=color, scatter=true,linewidth=5)
    de.plot!(modeltracks[logM], ax; color=:black, scatter=false,linewidth=2)
end
for extrapolTrack in extrapolGrid1.extrapoltracks
    de.plot!(extrapolTrack, ax; color=:red, scatter=false,linewidth=2,linestyle = :dash)
end
for extrapolTrack in extrapolGrid2.extrapoltracks
    de.plot!(extrapolTrack, ax; color=:blue, scatter=false,linewidth=2, linestyle = :dot)
end
#de.plot_arrows!(ax, modeltrack1, 0:0.1:1,0.1);de.plot_arrows!(ax, modeltrack2, 0:0.1:1,-0.1)
axislegend(ax,position = :lb)
ax2 = Axis(fig[1,2], xlabel=L"$\log (T_{\text{eff}} / K)$", ylabel=L"$\log (L/L_\odot)$", title="$interpol_asked",xreversed=true)
for (logM,color) in zip(interpol_asked,colors)
    interpolTrack = interpolTrackxs[logM]
    de.plot!(interpolTrack, ax2; color=color, scatter=true,linewidth=5)
    de.plot!(modeltracks[logM], ax2; color=:black, scatter=false,linewidth=2)
end
for extrapolTrack in extrapolGrid1.extrapoltracks
    de.plot!(extrapolTrack, ax2; color=:red, scatter=false,label="Extrapol from below",linewidth=2,linestyle = :dash)
end
for extrapolTrack in extrapolGrid2.extrapoltracks
    de.plot!(extrapolTrack, ax2; color=:blue, scatter=false,label="Extrapol from above",linewidth=2, linestyle = :dot)
end
leg = Legend(fig[1,3],ax2;tellwidth=true)
fig
##
save("Figures/master_MSInterpo.png", fig, px_per_unit=3)
##
## ################################"
interpol_asked = 0.05:0.05:0.45; interpolTrackxs = Dict(); diffs_with_JEMS = Dict()
modeltrack1, modeltrack2 = modeltracks[0.0], modeltracks[0.5]
extrapolGrid1 = de.ExtrapolGrid(modeltrack1, interpol_asked .- modeltrack1.logM);
extrapolGrid2 = de.ExtrapolGrid(modeltrack2, interpol_asked .- modeltrack2.logM);
for logM in interpol_asked
    interpolTrack = de.InterpolTrack(modeltrack1, modeltrack2,logM);
    interpolTrackxs[logM] = interpolTrack; diffs_with_JEMS[logM] = Difference_with_JEMS(interpolTrack, modeltracks[logM])
end
## ANDERE MASSA'S, horende bij zeta diagnostics
interpol_asked = 0.05:0.05:0.65; interpolTrackxs = Dict(); diffs_with_JEMS = Dict()
modeltrack1, modeltrack2 = modeltracks[0.0], modeltracks[0.7]
extrapolGrid1 = de.ExtrapolGrid(modeltrack1, interpol_asked .- modeltrack1.logM);
extrapolGrid2 = de.ExtrapolGrid(modeltrack2, interpol_asked .- modeltrack2.logM);
for logM in interpol_asked
    interpolTrack = de.InterpolTrack(modeltrack1, modeltrack2,logM);
    interpolTrackxs[logM] = interpolTrack; diffs_with_JEMS[logM] = Difference_with_JEMS(interpolTrack, modeltracks[logM])
end
##
fig = Figure(size=(1000,600)); logM1 = modeltrack1.logM; logM2 = modeltrack2.logM; title = "Interpolation between $logM1 and $logM2"
ax = Axis(fig[1,1], xlabel=L"$\log (T_{\text{eff}} / K)$", ylabel=L"$\log (L/L_\odot)$",xreversed=true)
de.plot!(modeltrack2, ax; color=:blue, scatter=false,label = "JEMS log M = " * string(modeltrack2.logM),linestyle=:dot, linewidth=5)
scatter!(modeltrack2.logT_val[1], modeltrack2.logL_val[1]; color=:blue, markersize=20)
de.plot!(modeltrack1, ax; color=:red, scatter=false, label = "JEMS log M = " * string(modeltrack1.logM),linestyle=:dash, linewidth=5)
scatter!(modeltrack1.logT_val[1], modeltrack1.logL_val[1]; color=:red, markersize=20)

#lines!(ax, [modeltrack1.logT_val[1],modeltrack2.logT_val[1]], [modeltrack1.logL_val[1],modeltrack2.logL_val[1]]; color=:black, linestyle=:dash,linewidth=2)
zamslogT = Vector{Float64}(); zamslogL = Vector{Float64}(); 
for logM in range(logM1, logM2, 100)
    interpolTrack = de.InterpolTrack(modeltrack1, modeltrack2,logM);
    push!(zamslogT, interpolTrack.logT_val[1]); push!(zamslogL, interpolTrack.logL_val[1])
end
lines!(ax, zamslogT, zamslogL; color=:black, linewidth=2, linestyle=:dash,label="Interpolated ZAMS")


colormap = :viridis; nb_colors = length(interpol_asked); colors = palette(colormap); vals = collect(1:nb_colors).*255 ./ nb_colors
colors = [colors[round(Int,vals[i])] for i in 1:nb_colors]
for (logM,color) in zip(interpol_asked,colors)
    interpolTrack = interpolTrackxs[logM]
    de.plot!(interpolTrack, ax; color=color, scatter=false,linewidth=10,alpha=0.5)
    scatter!(interpolTrack.logT_val[1], interpolTrack.logL_val[1]; color=color, markersize=20)
    de.plot!(modeltracks[logM], ax; color=:black, scatter=false,linewidth=4)
    scatter!(modeltracks[logM].logT_val[1], modeltracks[logM].logL_val[1]; color=:black, markersize=10)
    text!(ax, interpolTrack.logT_val[1]+0.02,interpolTrack.logL_val[1]-0.08; text = string(logM),fontsize=30)
end
ax.xgridvisible = true; ax.ygridvisible = true
for extrapolTrack in extrapolGrid1.extrapoltracks
    #de.plot!(extrapolTrack, ax; color=:red, scatter=false,linewidth=2,linestyle = :dash)
end
for extrapolTrack in extrapolGrid2.extrapoltracks
    #de.plot!(extrapolTrack, ax; color=:blue, scatter=false,linewidth=2, linestyle = :dot)
end
#de.plot_arrows!(ax, modeltrack1, 0:0.1:1,0.1);de.plot_arrows!(ax, modeltrack2, 0:0.1:1,-0.1)
axislegend(ax,position = :lb)

#leg = Legend(fig[1,3],ax2;tellwidth=true)
fig
##
figpath = "Figures/DSE_Interpolation_0.0_0.6.png"; save(figpath, fig, px_per_unit=3); @show figpath
##
fig = Figure(size=(2000,1000))
ax1 = Axis(fig[1,1],aspect=2, xlabel=L"\zeta", ylabel = L"\Delta \log (L/L_\odot)"); 
ax2 = Axis(fig[2,1],aspect=2, xlabel=L"\zeta", ylabel = L"\Delta \log (T_{\text{eff}} / K)");
ax3 = Axis(fig[3,1],aspect=2, xlabel=L"\zeta", ylabel = L"\Delta");
ax4 = Axis(fig[3,2],aspect=1.8,xlabel="log M", ylabel="Max/Mean")
ax5 = Axis(fig[3,3],aspect=1.8,xlabel="log M", )
ax6 = Axis(fig[3,4],aspect=1.8,xlabel="log M", )
axA = Axis(fig[1,2],aspect=1.8,                ylabel="Max/Mean")
axB = Axis(fig[1,3],aspect=1.8,                )
axC = Axis(fig[1,4],aspect=1.8,                )
axD = Axis(fig[2,2],aspect=1.8,                ylabel="Max/Mean")
axE = Axis(fig[2,3],aspect=1.8,                )
axF = Axis(fig[2,4],aspect=1.8,                )
ax1.xlabelvisible = false; ax2.xlabelvisible = false; ax1.xticklabelsvisible=false; ax2.xticklabelsvisible=false
for ax in [ax1,ax2,ax3,ax4,ax5,ax6, axA,axB,axC,axD,axE,axF]
    ax.xticklabelsvisible=false; ax.yticklabelsvisible = false
end
ax3.xticklabelsvisible = true; ax4.xticklabelsvisible = true; ax5.xticklabelsvisible = true; ax6.xticklabelsvisible = true
ax1.yticklabelsvisible = true; ax2.yticklabelsvisible = true; ax3.yticklabelsvisible = true
colormap = :viridis; nb_colors = length(interpol_asked); colors = palette(colormap); vals = collect(1:nb_colors).*255 ./ nb_colors
colors = [colors[round(Int,vals[i])] for i in 1:nb_colors]
for (logM,color) in zip(interpol_asked, colors)
    interpolTrack = interpolTrackxs[logM]; diff = diffs_with_JEMS[logM]; logM = diff.logM
    scatter!(ax1, interpolTrack.zeta, diff.diff_logL, label = "$logM" ,color=color )
    scatter!(ax2, interpolTrack.zeta, diff.diff_logT, label = "$logM" ,color=color )
    scatter!(ax3, interpolTrack.zeta, diff.diff_pyth,  label = "$logM",color=color )
    
    axA.title = L"\zeta \in[0.0,0.5]"
    axB.title = L"\zeta \in[0.5,1.0]"
    axC.title = L"Total $\zeta \in[0.0,1.0]$"
    split = 500
    diff_Ltot = mean(abs.(diff.diff_logL[1:split])   )   ; scatter!(axA, logM, diff_Ltot,markersize=30  ,color=color )
    diff_Ltot = mean(abs.(diff.diff_logL[split:1000]))   ; scatter!(axB, logM, diff_Ltot,markersize=30  ,color=color )
    diff_Ltot = mean(abs.(diff.diff_logL[1:1000])  )   ; scatter!(axC, logM, diff_Ltot,markersize=30  ,color=color )
    diff_Ttot = mean(abs.(diff.diff_logT[1:split])   )   ; scatter!(axD, logM, diff_Ttot,markersize=30  ,color=color )
    diff_Ttot = mean(abs.(diff.diff_logT[split:1000]))   ; scatter!(axE, logM, diff_Ttot,markersize=30  ,color=color )
    diff_Ttot = mean(abs.(diff.diff_logT[1:1000])  )   ; scatter!(axF, logM, diff_Ttot,markersize=30  ,color=color )
    diff_tot =  mean(abs.(diff.diff_pyth[1:split])   )   ; scatter!(ax4, logM, diff_tot,markersize=30  ,color=color )
    diff_tot =  mean(abs.(diff.diff_pyth[split:1000]))   ; scatter!(ax5, logM, diff_tot,markersize=30  ,color=color )
    diff_tot =  mean(abs.(diff.diff_pyth[1:1000])  )   ; scatter!(ax6, logM, diff_tot,markersize=30  ,color=color )
    
    diff_Ltot = maximum(abs.(diff.diff_logL[1:split])   )   ; scatter!(axA, logM, diff_Ltot,markersize=30  ,color=color , marker=:cross )
    diff_Ltot = maximum(abs.(diff.diff_logL[split:1000]))   ; scatter!(axB, logM, diff_Ltot,markersize=30  ,color=color , marker=:cross )
    diff_Ltot = maximum(abs.(diff.diff_logL[1:1000])  )   ; scatter!(axC, logM, diff_Ltot,markersize=30  ,color=color   , marker=:cross )
    diff_Ttot = maximum(abs.(diff.diff_logT[1:split])   )   ; scatter!(axD, logM, diff_Ttot,markersize=30  ,color=color , marker=:cross )
    diff_Ttot = maximum(abs.(diff.diff_logT[split:1000]))   ; scatter!(axE, logM, diff_Ttot,markersize=30  ,color=color , marker=:cross )
    diff_Ttot = maximum(abs.(diff.diff_logT[1:1000])  )   ; scatter!(axF, logM, diff_Ttot,markersize=30  ,color=color   , marker=:cross )
    diff_tot =  maximum(abs.(diff.diff_pyth[1:split])   )   ; scatter!(ax4, logM, diff_tot,markersize=30  ,color=color  , marker=:cross )
    diff_tot =  maximum(abs.(diff.diff_pyth[split:1000]))   ; scatter!(ax5, logM, diff_tot,markersize=30  ,color=color  , marker=:cross )
    diff_tot =  maximum(abs.(diff.diff_pyth[1:1000])  )   ; scatter!(ax6, logM, diff_tot,markersize=30  ,color=color    , marker=:cross )    
end
for ax in [ax4,ax5,ax6, axA,axB,axC,axD,axE,axF]
    vlines!(ax, [modeltrack1.logM, modeltrack2.logM], color=:black, linestyle=:dot,alpha=0.4)
end
xlims!(ax1, [0,1.0]); 
xlims!(ax2, [0,1.0]); 
xlims!(ax3, [0,1.0])
linkyaxes!(ax4,ax5,ax6); linkyaxes!(axA,axB,axC); linkyaxes!(axD,axE,axF)
#ylims!(ax1, [-0.0025,0.010]); ylims!(ax2, [-0.01,0.01]); ylims!(ax3, [0,0.5])
fig 
##
figpath = "Figures/DSE_Interpolation_0.0_0.6_diffs_diagnostics.png"; save(figpath, fig, px_per_unit=3); @show figpath






interpol_asked = [0.01]
interpolTrackxs = Dict()
modeltrack1, modeltrack2 = modeltracks[0.0], modeltracks[0.02]
extrapolGrid1 = de.ExtrapolGrid(modeltrack1, interpol_asked .- modeltrack1.logM);
extrapolGrid2 = de.ExtrapolGrid(modeltrack2, interpol_asked .- modeltrack2.logM);
for logM in interpol_asked
    interpolTrack = de.InterpolTrack(modeltrack1, modeltrack2,logM);
    interpolTrackxs[logM] = interpolTrack
end
##
fig = Figure(size=(2000,800))
ax = Axis(fig[1,1], xlabel=L"$\log (T_{\text{eff}} / K)$", ylabel=L"$\log (L/L_\odot)$", title=L"Interpolation between $\log M = 0.0$ and $0.02",xreversed=true)
de.plot!(modeltrack1, ax; color=:red, scatter=false, label = "JEMS log M = " * string(modeltrack1.logM),linestyle=:dash, linewidth=5)
de.plot!(modeltrack2, ax; color=:blue, scatter=false,label = "JEMS log M = " * string(modeltrack2.logM),linestyle=:dot, linewidth=5)
for logM in interpol_asked
    interpolTrack = interpolTrackxs[logM]
    de.plot!(interpolTrack, ax; color=:lightblue, scatter=true,linewidth=5)
    de.plot!(modeltracks[logM], ax; color=:black, scatter=false,linewidth=2)
end
for extrapolTrack in extrapolGrid1.extrapoltracks
    de.plot!(extrapolTrack, ax; color=:red, scatter=false,linewidth=2,linestyle = :dash)
end
for extrapolTrack in extrapolGrid2.extrapoltracks
    de.plot!(extrapolTrack, ax; color=:blue, scatter=false,linewidth=2, linestyle = :dot)
end
de.plot_arrows!(ax, modeltrack1, 0:0.01:1,0.01)
de.plot_arrows!(ax, modeltrack2, 0:0.01:1,-0.01)
axislegend(ax,position = :lb)
ax2 = Axis(fig[1,2], xlabel=L"$\log (T_{\text{eff}} / K)$", ylabel=L"$\log (L/L_\odot)$",xreversed=true)
for logM in interpol_asked
    interpolTrack = interpolTrackxs[logM]
    de.plot!(interpolTrack, ax2; color=:lightblue, scatter=true,label="Interpol",linewidth=5)
    de.plot!(modeltracks[logM], ax2; color=:black, scatter=false,label="JEMS",linewidth=1)
end
for extrapolTrack in extrapolGrid1.extrapoltracks
    de.plot!(extrapolTrack, ax2; color=:red, scatter=false,label="Extrapol from below",linewidth=2,linestyle = :dash)
end
for extrapolTrack in extrapolGrid2.extrapoltracks
    de.plot!(extrapolTrack, ax2; color=:blue, scatter=false,label="Extrapol from above",linewidth=2, linestyle = :dot)
end
leg = Legend(fig[1,3],ax2;tellwidth=true)
fig
##
save("Figures/MSInterpo_mini.png", fig)
##


interpol_asked = [0.05]
interpolTrackxs = Dict()
modeltrack1, modeltrack2 = modeltracks[0.0], modeltracks[0.1]
extrapolGrid1 = de.ExtrapolGrid(modeltrack1, interpol_asked .- modeltrack1.logM);
extrapolGrid2 = de.ExtrapolGrid(modeltrack2, interpol_asked .- modeltrack2.logM);
for logM in interpol_asked
    interpolTrack = de.InterpolTrack(modeltrack1, modeltrack2,logM);
    interpolTrackxs[logM] = interpolTrack
end
##
fig = Figure(size=(2000,800))
ax = Axis(fig[1,1], xlabel=L"$\log (T_{\text{eff}} / K)$", ylabel=L"$\log (L/L_\odot)$", title=L"Interpolation between $\log M = 0.0$ and $0.1",xreversed=true)
de.plot!(modeltrack1, ax; color=:red, scatter=false, label = "JEMS log M = " * string(modeltrack1.logM),linestyle=:dash, linewidth=5)
de.plot!(modeltrack2, ax; color=:blue, scatter=false,label = "JEMS log M = " * string(modeltrack2.logM),linestyle=:dot, linewidth=5)
for logM in interpol_asked
    interpolTrack = interpolTrackxs[logM]
    de.plot!(interpolTrack, ax; color=:lightblue, scatter=true,linewidth=5)
    de.plot!(modeltracks[logM], ax; color=:black, scatter=false,linewidth=2)
end
for extrapolTrack in extrapolGrid1.extrapoltracks
    de.plot!(extrapolTrack, ax; color=:red, scatter=false,linewidth=2,linestyle = :dash)
end
for extrapolTrack in extrapolGrid2.extrapoltracks
    de.plot!(extrapolTrack, ax; color=:blue, scatter=false,linewidth=2, linestyle = :dot)
end
de.plot_arrows!(ax, modeltrack1, 0:0.01:1,0.05)
de.plot_arrows!(ax, modeltrack2, 0:0.01:1,-0.05)
axislegend(ax,position = :lb)
ax2 = Axis(fig[1,2], xlabel=L"$\log (T_{\text{eff}} / K)$", ylabel=L"$\log (L/L_\odot)$",xreversed=true)
for logM in interpol_asked
    interpolTrack = interpolTrackxs[logM]
    de.plot!(interpolTrack, ax2; color=:lightblue, scatter=true,label="Interpol",linewidth=5)
    de.plot!(modeltracks[logM], ax2; color=:black, scatter=false,label="JEMS log M = 0.05",linewidth=1)
end
for extrapolTrack in extrapolGrid1.extrapoltracks
    de.plot!(extrapolTrack, ax2; color=:red, scatter=false,label="Extrapol from below",linewidth=2,linestyle = :dash)
end
for extrapolTrack in extrapolGrid2.extrapoltracks
    de.plot!(extrapolTrack, ax2; color=:blue, scatter=false,label="Extrapol from above",linewidth=2, linestyle = :dot)
end
leg = Legend(fig[1,3],ax2;tellwidth=true)
fig
##
save("Figures/MSInterpo_mini.png", fig)


