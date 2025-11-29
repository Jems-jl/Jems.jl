
# TODO: Probably good idea to make these hardcoded options more flexible in the future
kipp_mixing_colors = [RGBAf(0.5, 0.5, 0.5), RGBAf(0, 0, 1)]
kipp_mixing_colors[1] = RGBAf(1, 1, 1)  # make no_mixing white in KippenLine diagram to avoid clutter
burning_colors = cgrad(:linear_wyor_100_45_c55_n256)
const mixing_map = Dict(:no_mixing => 1,
                        :convection => 2)

function basic_theme()
    return Theme(fonts=(regular=texfont(:text), bold=texfont(:bold),
                           italic=texfont(:italic), bold_italic=texfont(:bolditalic)),
                    fontsize=20, size=(1000, 750), linewidth=7,
                    Axis=(xlabelsize=30, ylabelsize=30, titlesize=30, xgridvisible=false, ygridvisible=false,
                          spinewidth=2.5, xminorticksvisible=true, yminorticksvisible=true, xtickalign=1, ytickalign=1,
                          xminortickalign=1, yminortickalign=1, xticksize=14, xtickwidth=2.5, yticksize=14,
                          ytickwidth=2.5, xminorticksize=7, xminortickwidth=2.5, yminorticksize=7, yminortickwidth=2.5,
                          xticklabelsize=25, yticklabelsize=25, xticksmirrored=true, yticksmirrored=true),
                    Legend=(patchsize=(70, 10), framevisible=false, patchlabelgap=20, rowgap=10))
end
