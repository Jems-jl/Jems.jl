using Jems.Interpolations
using ForwardDiff

export set_EOS_resultsTρ!, EOSTableCollector


"""
    TρTableCollector <: AbstractEOS

Represents a single rectangular Equation of State (EOS) data table mapped over the 
density-temperature parameter space (`logQ` and `logT`) for a specific chemical composition.

# Fields
- `X::Float64`: The Hydrogen mass fraction of this specific table.
- `Z::Float64`: The Metallicity (heavy element) mass fraction of this specific table.
- `interpolator::Vector{BicubicInterpolation}`: A vector of pre-computed bicubic interpolators. 
  Each index corresponds to a specific thermodynamic quantity (e.g., pressure, internal energy, entropy) 
  defined over the (`logQ`, `logT`) grid.
- `col_names::Vector{String}`: The physical names of the variables corresponding to the `interpolator` indices.
"""
struct TρTableCollector <: AbstractEOS
    X :: Float64 
    Z :: Float64
    interpolator :: Vector{BilinearInterpolation}
    col_names :: Vector{String}
end 


"""
    EOS_Z_Slice

A structural grouping of `TρTableCollector` objects that all share the exact same Metallicity (Z). 
This acts as a 1D lookup table across different Hydrogen mass fractions (X).

# Fields
- `Z::Float64`: The shared Metallicity mass fraction for all tables in this slice.
- `Xs::Vector{Float64}`: A sorted array of the available Hydrogen mass fractions.
- `tables::Vector{TρTableCollector}`: The individual EOS tables corresponding to each X value, 
  sorted to perfectly match the order of `Xs`.
"""
struct EOS_Z_Slice
    Z :: Float64
    Xs :: Vector{Float64}
    tables :: Vector{TρTableCollector}
end

"""
    EOSTableCollector <: AbstractEOS

Struct handling the full 4D Equation of State interpolation grid. 
It organizes tabulated data into a nested structure to facilitate fast interpolation across 
Metallicity (Z), Hydrogen fraction (X), log10Q, and log10T.

# Fields
- `Zs::Vector{Float64}`: A sorted array of the available Metallicity grid points.
- `slices::Vector{EOS_Z_Slice}`: The Z-grouped slices of EOS tables, sorted to match `Zs`.
- `include_radiation::Bool`: Flag to determine whether analytical radiation pressure 
  (P_rad = a T^4 / 3) should be explicitly added to the gas pressure evaluated from the tables.
"""
struct EOSTableCollector <: AbstractEOS
    Zs :: Vector{Float64}
    slices :: Vector{EOS_Z_Slice}
    include_radiation :: Bool
end

"""
    TρTableCollector(filepath)

Constructor that parses an EOS data file from the provided `filepath`. 
Extracts metadata (version, X, Z, grid dimensions), streams the data block into a matrix, 
and builds bicubic interpolators for every column/variable present in the table.
"""
function TρTableCollector(filepath :: String)
    local X_val::Float64
    local Z_val::Float64
    local num_Ts::Int
    local num_Qs::Int
    local num_vars::Int
    local logTs::Vector{Float64}
    local logQs::Vector{Float64}
    local col_names::Vector{String}
    
    # READ ENTIRE FILE INTO MEMORY ONCE 
    file_str = read(filepath, String)
    len = sizeof(file_str)
    cursor = 1

    
    @inline function next_line_end(str, start_idx)
        idx = start_idx
        while idx <= len && str[idx] != '\n' && str[idx] != '\r'
            idx = nextind(str, idx)
        end
        return idx
    end

    
    while cursor <= len
        line_end = next_line_end(file_str, cursor)
        line = SubString(file_str, cursor, line_end - 1)
        
        if occursin("version", line)
            # Move cursor to the next line to read values
            cursor = line_end
            while cursor <= len && (file_str[cursor] == '\n' || file_str[cursor] == '\r')
                cursor = nextind(file_str, cursor)
            end
            
            line_end = next_line_end(file_str, cursor)
            val_line = SubString(file_str, cursor, line_end - 1)
            parts = split(val_line)
            
            X_val = parse(Float64, parts[2])
            Z_val = parse(Float64, parts[3])
            num_Ts   = parse(Int, parts[4])
            logT_min = parse(Float64, parts[5])
            logT_del = parse(Float64, parts[7])
            
            num_Qs    = parse(Int, parts[8])
            logQ_min  = parse(Float64, parts[9])
            logQ_del  = parse(Float64, parts[11])
            
            logTs = collect(range(logT_min, step = logT_del, length = num_Ts))
            logQs = collect(range(logQ_min, step = logQ_del, length = num_Qs))
            
        elseif occursin("logE", line)
            col_names = String.(split(line))
            num_vars = length(col_names)
            
            # Advance cursor past the header line and exit metadata search
            cursor = line_end
            while cursor <= len && (file_str[cursor] == '\n' || file_str[cursor] == '\r')
                cursor = nextind(file_str, cursor)
            end
            break
        end
        
        # Advance to next line
        cursor = line_end
        while cursor <= len && (file_str[cursor] == '\n' || file_str[cursor] == '\r')
            cursor = nextind(file_str, cursor)
        end
    end 
    
    
    eos_data_flat = Matrix{Float64}(undef, num_vars, num_Ts * num_Qs) 
    data_idx = 1

    while cursor <= len && data_idx <= (num_Ts * num_Qs)
        line_end = next_line_end(file_str, cursor)
        
        is_numeric_row = true
        col_idx = 1
        line_cursor = cursor
        
        while line_cursor < line_end
            # Skip spaces
            while line_cursor < line_end && isspace(file_str[line_cursor])
                line_cursor = nextind(file_str, line_cursor)
            end
            if line_cursor >= line_end
                break
            end

            # Find token end
            token_start = line_cursor
            while line_cursor < line_end && !isspace(file_str[line_cursor])
                line_cursor = nextind(file_str, line_cursor)
            end
            token_end = prevind(file_str, line_cursor)

            # Use tryparse to cleanly reject non-numeric text and skip the line
            val = tryparse(Float64, SubString(file_str, token_start, token_end))
            
            if isnothing(val)
                is_numeric_row = false
                break
            end
            
            if col_idx <= num_vars
                eos_data_flat[col_idx, data_idx] = val
            end
            col_idx += 1
        end

        # Verify it was a real data row before advancing data_idx
        if is_numeric_row && (col_idx - 1) == num_vars
            data_idx += 1
        end
        
        # Advance to next line
        cursor = line_end
        while cursor <= len && (file_str[cursor] == '\n' || file_str[cursor] == '\r')
            cursor = nextind(file_str, cursor)
        end
    end
    
    eos_data_3D = reshape(eos_data_flat, (num_vars, num_Ts, num_Qs))
    
    # 4. BUILD CONCRETE BICUBIC INTERPOLATORS (Using zero-allocation views)

    first_slice = reshape(view(eos_data_3D, 1, :, :), (1, num_Ts, num_Qs))
    first_interp = build_bilinear_interpolator(logTs, logQs, first_slice)
    
    interpolators = [first_interp]
    
    for k in 2:num_vars 
        data_slice = reshape(view(eos_data_3D, k, :, :), (1, num_Ts, num_Qs))
        interp = build_bilinear_interpolator(logTs, logQs, data_slice)
        push!(interpolators, interp)
    end 
    
    return TρTableCollector(X_val, Z_val, interpolators, col_names)
end
# The lines below are for bicubic interpolation
#     first_slice = reshape(view(eos_data_3D, 1, :, :), (1, num_Ts, num_Qs))
#     first_interp = build_bicubic_interpolator(logTs, logQs, first_slice)
    
#     interpolators = [first_interp]
    
#     for k in 2:num_vars 
#         data_slice = reshape(view(eos_data_3D, k, :, :), (1, num_Ts, num_Qs))
#         interp = build_bicubic_interpolator(logTs, logQs, data_slice)
#         push!(interpolators, interp)
#     end 
    
#     return TρTableCollector(X_val, Z_val, interpolators, col_names)
# end
"""
    EOSTableCollector(directory,include_radiation)

Constructor that scans a given `directory` for `.data` files, parses each into a `TρTableCollector`, 
and groups them sequentially into `EOS_Z_Slice` objects to construct the full 4D EOS grid.
"""
function EOSTableCollector(directory :: String; include_radiation::Bool = true)
    files = readdir(directory; join = true)
    candidate_files = filter(f -> endswith(f, ".data"), files)
    #Load all tables
    all_tables = [TρTableCollector(f) for f in candidate_files]

    unique_Zs = sort(unique([t.Z for t in all_tables]))

    slices = EOS_Z_Slice[]

    for z_val in unique_Zs
        tables_at_z = filter(t -> t.Z == z_val, all_tables)

        #Sort by X in each Z slice
        sort!(tables_at_z, by = t -> t.X)
        #X values in each Z slice 
        xs_at_z = [t.X for t in tables_at_z]

        push!(slices, EOS_Z_Slice(z_val, xs_at_z, tables_at_z))
    end

    return EOSTableCollector(unique_Zs, slices, include_radiation)
end 

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
    set_EOS_resultsTρ!(eos, r, lnT, lnρ, xa, species)

Main evaluation function to populate an `EOSResults` object `r` for a given Temperature (`lnT`) 
and Density (`lnρ`). Smoothly blends between table lookups (for standard stellar interiors) 
and analytical Ideal Gas results (for boundary/low-T regimes). 
Handles bilinear interpolation across composition (X, Z) and bicubic interpolation across Q and T.
"""
function set_EOS_resultsTρ!(eos::EOSTableCollector, r::EOSResults{TT}, lnT::TT, lnρ::TT,
                           xa::AbstractVector{<:TT}, species::Vector{Symbol}) where {TT<:Real}

   
    INV_LN10 = 0.4342944819032518
    logT = lnT * INV_LN10
    r.lnT, r.lnρ = lnT, lnρ
    r.T, r.ρ = exp(lnT), exp(lnρ)

    # Weights for transition regimes
    w_low = smooth_step_func(logT, 3.1, 3.2)
    w_high = 1.0 - smooth_step_func(logT, 7.0, 7.4)
    w = w_low * w_high
    iw = 1.0 - w

    
    if w <= 0.0001
        i_μ = get_μ_IdealEOS(xa, species)
        r.μ = i_μ
        r.Prad = eos.include_radiation ? (CRAD * r.T^4 / 3.0) : 0.0
        P_gas = CGAS * r.T * r.ρ / i_μ
        r.P = P_gas + r.Prad
        i_β = P_gas / r.P
        r.β, r.χ_ρ = i_β, i_β
        r.χ_T = 4.0 - 3.0 * i_β
        r.α, r.δ = 1.0 / r.χ_ρ, r.χ_T / r.χ_ρ
        r.u = CGAS * r.T / i_μ * (1.5 + 3.0 * (1.0 - i_β) / i_β)
        r.cₚ = CGAS / i_μ * (1.5 + 3.0 * (4.0 + i_β) * (1.0 - i_β) * r.α^2 + r.δ * r.α)
        r.∇ₐ = CGAS * r.δ / (i_β * i_μ * r.cₚ)
        r.Γ₁ = 1.0 / (r.α - r.δ * r.∇ₐ)
        r.lnP = log(r.P)
        return
    end

    
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

    
    logQ = (lnρ * INV_LN10) - 2.0*logT + 12.0

    
    sample_interp = eos.slices[1].tables[1].interpolator[1]
    i, j, u, v = get_data_position(sample_interp.grid_x, sample_interp.grid_y, logT, logQ)

    # Z-direction lookup
    val_Z = ForwardDiff.value(Z_val)
    iz = clamp(searchsortedlast(eos.Zs, val_Z), 1, length(eos.Zs) - 1)
    u_Z = (Z_val - eos.Zs[iz]) / (eos.Zs[iz+1] - eos.Zs[iz])

    slice0, slice1 = eos.slices[iz], eos.slices[iz+1]
    
    # X-direction lookup 
    val_X = ForwardDiff.value(X_val)
    ix0 = clamp(searchsortedlast(slice0.Xs, val_X), 1, length(slice0.Xs) - 1)
    ix1 = clamp(searchsortedlast(slice1.Xs, val_X), 1, length(slice1.Xs) - 1)

    u_X0 = (X_val - slice0.Xs[ix0]) / (slice0.Xs[ix0+1] - slice0.Xs[ix0])
    u_X1 = (X_val - slice1.Xs[ix1]) / (slice1.Xs[ix1+1] - slice1.Xs[ix1])

 
    function get_blended_var(idx)
        v00 = evaluate_interp(slice0.tables[ix0].interpolator[idx], i, j, u, v)
        v10 = evaluate_interp(slice0.tables[ix0+1].interpolator[idx], i, j, u, v)
        v01 = evaluate_interp(slice1.tables[ix1].interpolator[idx], i, j, u, v)
        v11 = evaluate_interp(slice1.tables[ix1+1].interpolator[idx], i, j, u, v)
        
        vz0 = v00 * (1 - u_X0) + v10 * u_X0
        vz1 = v01 * (1 - u_X1) + v11 * u_X1
        return vz0 * (1 - u_Z) + vz1 * u_Z
    end

    # Extract needed variables
    t_P_gas = 10^get_blended_var(2)
    t_u     = 10^get_blended_var(3)
    t_χ_ρ   = get_blended_var(5)
    t_χ_T   = get_blended_var(6)
    t_cₚ    = get_blended_var(7)
    t_μ     = get_blended_var(12)
    t_Γ₁    = get_blended_var(14)
    t_∇ₐ    = get_blended_var(16)

    # Blending and Assignment
    r.Prad = eos.include_radiation ? (CRAD * r.T^4 / 3.0) : 0.0

    if w >= 0.9999
        r.μ, r.u, r.χ_ρ, r.χ_T, r.cₚ, r.Γ₁, r.∇ₐ = t_μ, t_u, t_χ_ρ, t_χ_T, t_cₚ, t_Γ₁, t_∇ₐ
        r.P = t_P_gas + r.Prad
    else
        i_μ = get_μ_IdealEOS(xa, species)
        i_P_gas = CGAS * r.T * r.ρ / i_μ
        i_P_tot = i_P_gas + r.Prad
        i_β = i_P_gas / i_P_tot
        i_χ_ρ, i_χ_T = i_β, 4.0 - 3.0 * i_β
        i_α, i_δ = 1.0 / i_χ_ρ, i_χ_T / i_χ_ρ
        i_u  = CGAS * r.T / i_μ * (1.5 + 3.0 * (1.0 - i_β) / i_β)
        i_cₚ = CGAS / i_μ * (1.5 + 3.0 * (4.0 + i_β) * (1.0 - i_β) * i_α^2 + i_δ * i_α)
        i_∇ₐ = CGAS * i_δ / (i_β * i_μ * i_cₚ)
        i_Γ₁ = 1.0 / (i_α - i_δ * i_∇ₐ)
        
        r.μ   = w * t_μ     + iw * i_μ
        P_gas = w * t_P_gas + iw * i_P_gas
        r.P   = P_gas + r.Prad
        r.u   = w * t_u     + iw * i_u
        r.χ_ρ = w * t_χ_ρ   + iw * i_χ_ρ
        r.χ_T = w * t_χ_T   + iw * i_χ_T
        r.cₚ  = w * t_cₚ    + iw * i_cₚ
        r.∇ₐ  = w * t_∇ₐ    + iw * i_∇ₐ
        r.Γ₁  = w * t_Γ₁    + iw * i_Γ₁
    end

    r.lnP = log(r.P)
    r.β   = 1.0 - r.Prad / r.P
    r.α   = 1.0 / r.χ_ρ
    r.δ   = r.χ_T / r.χ_ρ
end

