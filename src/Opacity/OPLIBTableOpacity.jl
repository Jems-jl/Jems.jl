using ForwardDiff
using Jems.Interpolations

export get_opacity_resultsTρ,RTTableOpacity, OpacityTableCollector, CompositeOpacity,
    get_opacity_table_collection
"""
    RTTableOpacity <: AbstractOpacity

Structure representing a single pre-calculated opacity table interpolation coefficients for a specific chemical composition.

# Fields
- `X::Float64`: Hydrogen mass fraction of this specific table.
- `Z::Float64`: Metallicity mass fraction of this specific table.
- `interpolator::BicubicInterpolation`: The bicubic interpolator built from the raw opacity data 
  (typically mapped over `log10(R)` and `log10(T)`) for this specific `(X, Z)` composition.
"""
struct RTTableOpacity <: AbstractOpacity 
    X :: Float64 
    Z :: Float64 
    interpolator :: BicubicInterpolation
end 

"""
    OpacityTableCollector <: AbstractOpacity

A grid of opacity tables spanning different Hydrogen (`X`) and Metallicity (`Z`) mass fractions.
This structure is used to interpolate opacities across different chemical compositions.

# Fields
- `Xs::Vector{Float64}`: Array of Hydrogen mass fractions defining the grid's X-axis.
- `Zs::Vector{Float64}`: Array of Metallicity mass fractions defining the grid's Z-axis.
- `tables::Matrix{RTTableOpacity}`: A 2D array of `RTTableOpacity` objects corresponding to 
  the `(X, Z)` coordinates.
"""
struct OpacityTableCollector <: AbstractOpacity
    Xs :: Vector{Float64}
    Zs :: Vector{Float64}
    tables :: Matrix{RTTableOpacity}
end 

"""
    CompositeOpacity <: AbstractOpacity

A composite structure that manages the transition between low-temperature and high-temperature 
opacity regimes (e.g., blending Ferguson/Alexander tables with OPAL tables).

# Fields
- `low_T_collector::OpacityTableCollector`: The collection of opacity tables used for the low-temperature regime.
- `high_T_collector::OpacityTableCollector`: The collection of opacity tables used for the high-temperature regime.
- `trans_logT_min::Float64`: The `log10(T)` boundary below which only the `low_T_collector` is used.
- `trans_logT_max::Float64`: The `log10(T)` boundary above which only the `high_T_collector` is used.
  (Temperatures falling between `min` and `max` are smoothly blended between both collectors).
"""
struct CompositeOpacity <: AbstractOpacity 
    low_T_collector :: OpacityTableCollector
    high_T_collector :: OpacityTableCollector
    trans_logT_min :: Float64
    trans_logT_max :: Float64 
end

"""
    RTTableOpacity(filepath)

Reads and parses a structured opacity data file to construct an `RTTableOpacity` object.
This constructor expects a specific tabular format (commonly used in OPAL or Ferguson/Alexander 
opacity tables). It dynamically searches the file for metadata (composition and grid sizes), 
extracts the `log10(R)` and `log10(T)` axes, reads the 2D opacity data, and builds the 
underlying bicubic interpolator.
"""
function RTTableOpacity(filepath :: String)
    local X_val::Float64
    local Z_val::Float64
    local num_Rs::Int
    local num_Ts::Int
    logRs = Float64[]

    local logTs::Vector{Float64}
    local data::Matrix{Float64}

    open(filepath, "r") do file 

        # Metadata Search  
        for line in eachline(file)
            if occursin(r"^\s*1\s+\d+", line)
                parts = split(line)
                X_val = parse(Float64, parts[3])
                Z_val = parse(Float64, parts[4])
                num_Rs = parse(Int, parts[5])
                num_Ts = parse(Int, parts[8])
                break 
            end
        end
        
        # grid start search 
        for line in eachline(file)
            row_vals = split(line)
            if length(row_vals) == num_Rs && all(x -> tryparse(Float64, x) !== nothing, row_vals)
                logRs = parse.(Float64, row_vals)
                break
            end
        end 

        local first_data_parts
        for line in eachline(file)
            row_vals = split(line)
            if length(row_vals) == (num_Rs + 1) && all(x -> tryparse(Float64, x) !== nothing, row_vals)
                first_data_parts = row_vals
                break
            end
        end

        # Pre-allocate 2D matrix directly
        logTs = Vector{Float64}(undef, num_Ts)
        data = Matrix{Float64}(undef, num_Rs, num_Ts)

        # Parse first temperature row
        logTs[1] = parse(Float64, first_data_parts[1])
        for R_idx in 1:num_Rs
            data[R_idx, 1] = parse(Float64, first_data_parts[R_idx + 1])
        end

        # Stream the rest of the lines
        T_idx = 2
        for line in eachline(file)
            if T_idx > num_Ts
                break
            end       
            row_parts = split(line)     
            # Skip empty lines
            if isempty(row_parts) 
                continue 
            end

            if length(row_parts) != (num_Rs + 1)
                error("Format Error at row $T_idx: Expected $(num_Rs+1) columns, found $(length(row_parts))")
            end

            logTs[T_idx] = parse(Float64, row_parts[1])
            
            for R_idx in 1:num_Rs
                val = parse(Float64, row_parts[R_idx + 1])
                data[R_idx, T_idx] = val
            end
            
            T_idx += 1
        end
    end   
    raw_data_3d = reshape(data, (1, size(data, 1), size(data, 2)))
    bicubic_interp = build_bicubic_interpolator(logRs, logTs, raw_data_3d)
    return RTTableOpacity(X_val, Z_val, bicubic_interp)
end

"""
    OpacityTableCollector(directory, mixture_key)

Scans a directory for opacity data files, filters them by a specific heavy-element mixture, 
and organizes them into a 2D rectangular grid spanning Hydrogen (`X`) and Metallicity (`Z`).
This constructor automates the bulk-loading of dozens or hundreds of individual `RTTableOpacity` 
objects. It extracts the unique `X` and `Z` coordinates from the loaded files and places each 
table into its corresponding `(X, Z)` slot in a pre-allocated matrix.
"""
function OpacityTableCollector(directory :: String, mixture_key :: String)
    files = readdir(directory; join= true)
    opacity_files = filter(f -> endswith(f, ".data"), files) # A safety check if the directory contains other things

    if isempty(opacity_files)
        error("No .data files found in $directory")
    end

    temp_tables = RTTableOpacity[]

    for file in opacity_files
        if occursin(mixture_key, basename(file))
            try
                push!(temp_tables, RTTableOpacity(file))
            catch e 
                println("Warning: Skipping file $file due to error: $e")
            end
        end
    end

    unique_Xs = sort(unique([t.X for t in temp_tables]))
    unique_Zs = sort(unique([t.Z for t in temp_tables]))

    grid = Matrix{RTTableOpacity}(undef, length(unique_Xs), length(unique_Zs)) #Creating a ractangular grid required for interpolation scheme 

    for t in temp_tables
        x_idx = searchsortedfirst(unique_Xs, t.X)
        z_idx = searchsortedfirst(unique_Zs, t.Z)
        grid[x_idx, z_idx] = t 
    end 
    return OpacityTableCollector(unique_Xs, unique_Zs, grid)

end 

"""
    get_opacity_table_collection(collection, lnT, lnρ, xa, species)

Performs a 4D interpolation to find the Rosseland mean opacity.
- Computes the composition (X, Z) and thermodynamic state (logT, logR).
- Locates the 4 bounding (X, Z) tables in the collection.
- Performs a bicubic interpolation in (logR, logT) space inside each of the 4 tables.
- Performs a bilinear interpolation in (X, Z) space to combine the 4 bicubic results.
"""
function get_opacity_table_collection(collection::OpacityTableCollector, lnT::TT, lnρ::TT, xa::AbstractVector{<:TT}, species::Vector{Symbol}) where {TT<:Real}
    
    inv_ln10 = 0.4342944819032518
    logT = lnT * inv_ln10
    logρ = lnρ * inv_ln10
    logR = logρ - 3*logT + 18
    
    # Evaluate composition 
    iH1  = findfirst(==(:H1), species)
    iH2  = findfirst(==(:H2), species)
    iHe3 = findfirst(==(:He3), species)
    iHe4 = findfirst(==(:He4), species)
    
    val_H1  = isnothing(iH1)  ? zero(TT) : xa[iH1]
    val_H2  = isnothing(iH2)  ? zero(TT) : xa[iH2]
    val_He3 = isnothing(iHe3) ? zero(TT) : xa[iHe3]
    val_He4 = isnothing(iHe4) ? zero(TT) : xa[iHe4]
   
    X_val = val_H1 + val_H2
    Y_val = val_He3 + val_He4 
    Z_val = 1.0 - X_val - Y_val
    
    # Safe value extraction for searchsortedlast
    val_X = ForwardDiff.value(X_val)
    val_Z = ForwardDiff.value(Z_val)
    
    Nx = length(collection.Xs)
    Nz = length(collection.Zs)

    # Find which (X,Z) box we're in
    ix = clamp(searchsortedlast(collection.Xs, val_X), 1, Nx-1)
    iz = clamp(searchsortedlast(collection.Zs, val_Z), 1, Nz - 1)

    # Bilinear weights for (X,Z) interpolation
    u_X = (X_val - collection.Xs[ix]) / (collection.Xs[ix+1] - collection.Xs[ix])
    u_Z = (Z_val - collection.Zs[iz]) / (collection.Zs[iz+1] - collection.Zs[iz])

    # Get the 4 bounding tables for (X,Z)
    ix0 = ix
    ix1 = clamp(ix + 1, 1, Nx)
    iz0 = iz
    iz1 = clamp(iz + 1, 1, Nz)

    t00 = isassigned(collection.tables, ix0, iz0) ? collection.tables[ix0, iz0] : collection.tables[ix, iz]
    t10 = isassigned(collection.tables, ix1, iz0) ? collection.tables[ix1, iz0] : collection.tables[ix, iz]
    t01 = isassigned(collection.tables, ix0, iz1) ? collection.tables[ix0, iz1] : collection.tables[ix, iz]
    t11 = isassigned(collection.tables, ix1, iz1) ? collection.tables[ix1, iz1] : collection.tables[ix, iz]

    # Get table grid bounds for clamping
    min_T = max(t00.interpolator.grid_y[1], t10.interpolator.grid_y[1], 
                t01.interpolator.grid_y[1], t11.interpolator.grid_y[1])
    max_T = min(t00.interpolator.grid_y[end], t10.interpolator.grid_y[end], 
                t01.interpolator.grid_y[end], t11.interpolator.grid_y[end])
    
    min_R = max(t00.interpolator.grid_x[1], t10.interpolator.grid_x[1], 
                t01.interpolator.grid_x[1], t11.interpolator.grid_x[1])
    max_R = min(t00.interpolator.grid_x[end], t10.interpolator.grid_x[end], 
                t01.interpolator.grid_x[end], t11.interpolator.grid_x[end])
    
    # Clamp logT and logR to valid range (preserves derivatives with Dual numbers)
    eff_logT = clamp(logT, min_T, max_T)
    eff_logR = clamp(logR, min_R, max_R)

    # Data postition in the (R,T) tables
    i, j, u, v = get_data_position(t00.interpolator.grid_x, t00.interpolator.grid_y, eff_logR, eff_logT)

    # Bicubic interpolation in (R,T)
    k00 = evaluate_interp(t00.interpolator, i, j, u, v)
    k10 = evaluate_interp(t10.interpolator, i, j, u, v)
    k01 = evaluate_interp(t01.interpolator, i, j, u, v)
    k11 = evaluate_interp(t11.interpolator, i, j, u, v)

    # Bilinear interpolation between the 4 opacity values
    k_z0 = k00 * (1 - u_X) + k10 * u_X
    k_z1 = k01 * (1 - u_X) + k11 * u_X
    log_kappa_final = k_z0 * (1 - u_Z) + k_z1 * u_Z

    return 10^log_kappa_final
end

"""
    smooth_step_func(x, floor, ceil)

Calculates a smooth Hermite interpolation weight between 0.0 and 1.0
to ensure continuous and differentiable transitions between physical regimes.
Used for blending low-T and high-T opacities
"""
@inline function smooth_step_func(x::T, floor::Float64, ceil::Float64) where T
    if x <= floor 
        return zero(T)
    elseif x >= ceil
        return one(T)
    else
        t = (x - floor)/(ceil - floor)

        return t * t * (3.0 -2.0 * t)
    end 

end

"""
    get_opacity_resultsTρ(composite::CompositeOpacity, lnT, lnρ, xa, species)

Evaluates the local temperature to determine if the physical state falls in the 
low-temperature regime, the high-temperature regime, or the transition zone. 
In the transition zone, it performs a mathematically smooth blend of both tables in log-space. 
Finally, calls get_opacity_table_collection to calculate the opacity.
"""
function get_opacity_resultsTρ(composite :: CompositeOpacity, lnT::TT, lnρ::TT, xa::AbstractVector{<:TT}, species::Vector{Symbol})::TT where {TT <: Real}
    inv_ln10 = 0.4342944819032518
    logT = lnT * inv_ln10
    val_logT = ForwardDiff.value(logT)

    # calculating the weight based on log

    if val_logT >= composite.trans_logT_max 
        return get_opacity_table_collection(composite.high_T_collector, lnT, lnρ, xa, species)

    elseif val_logT <= composite.trans_logT_min
        return get_opacity_table_collection(composite.low_T_collector, lnT, lnρ, xa, species)

    else 
        κ_low  = get_opacity_table_collection(composite.low_T_collector, lnT, lnρ, xa, species)
        κ_high = get_opacity_table_collection(composite.high_T_collector, lnT, lnρ, xa, species)

        #calculating the weight 
        w = smooth_step_func(logT, composite.trans_logT_min, composite.trans_logT_max)

        log_κ_low = log10(κ_low)
        log_κ_high = log10(κ_high)
        smooth_log_κ = (1-w) * log_κ_low + w * log_κ_high

        return 10^smooth_log_κ

    end
end
