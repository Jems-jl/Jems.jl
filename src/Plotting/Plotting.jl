module Plotting

using Makie, LaTeXStrings, MathTeXEngine, Jems.DualSupport, Jems.Constants

include("Theming.jl")
include("Plotter.jl")
include("HRPlot.jl")
include("TRhoProfile.jl")
include("KippenLine.jl")
include("Abundances.jl")
include("HistoryPlot.jl")
include("ProfilePlot.jl")

end  # end module Plotting
