using Jems 
using ForwardDiff
using DelimitedFiles

# export RT_table_opacity

"""
Structure to store data from a single OPLIB data file and store it as a 1D array of R,T alon with some metadata 
"""
struct RT_table_opacity <: Jems.Opacity.AbstractOpacity 
    X :: Float64 
    Z :: Float64 
    logTs :: Vector{Float64}
    logRs :: Vector{Float64}
    kap_data :: Vector{Float64} 
end 

"""
Structure to store a matrix with X,Z being the rows and columns, and each element being a RT opacity table 
"""
struct Opacity_table_collector <: Jems.Opacity.AbstractOpacity
    Xs :: Vector{Float64}
    Zs :: Vector{Float64}
    tables :: Matrix{RT_table_opacity} #Grid containing each opacity table at X,Z
end 

"""
Structure to store opacaity tables of X,Z,R,T from different data set, also the meta data required to smooth interopolate between them 
using a cubic spline 
"""
struct CompositeOpacity <: Jems.Opacity.AbstractOpacity 
    low_T_collector :: Opacity_table_collector # Holds the matrix of X,Z with each element being RT_table (for low T)
    high_T_collector :: Opacity_table_collector #Holds the matrix of X,Z with each element being RT_table (for high T)

    # Holds the transition points 
    trans_logT_min :: Float64
    trans_logT_max :: Float64 

end 
    
function RT_table_opacity(filepath :: String)
    
    #Reading the whole file into a string array 
    lines = readlines(filepath) 

    meta_idx = findfirst(l -> occursin(r"^\s*1\s+\d+", l), lines) # Looks for 1. for the start of metadata
    if meta_idx === nothing 
        error("Could not find Metadata line (starting with ' 1 '). File format unknown.")
    end

    parts = split(lines[meta_idx]) # Split the string into structured data and removing the spaces

    X_val = parse(Float64, parts[3])
    Z_val = parse(Float64, parts[4])
    num_Rs = parse(Int, parts[5])
    num_Ts = parse(Int, parts[8])

    grid_line_idx = 0 
    logRs = Float64[]

    for i in (meta_idx + 1):min(meta_idx + 20, length(lines))
        row_vals = split(lines[i])
        # Check if line is all numbers and has correct length
        if length(row_vals) == num_Rs && all(x -> tryparse(Float64, x) !== nothing, row_vals)
            grid_line_idx = i
            logRs = parse.(Float64, row_vals)
            break
        end
    end

    if grid_line_idx == 0
        error("Could not find Grid Line! Searched 20 lines after metadata for a row with $num_Rs numbers.")
    end
    # println("Grid (logRs) found at Line $grid_line_idx.")

    data_start_idx = 0
    
    for i in (grid_line_idx + 1):min(grid_line_idx + 10, length(lines))
        row_vals = split(lines[i])
        if length(row_vals) == (num_Rs + 1) && all(x -> tryparse(Float64, x) !== nothing, row_vals)
            data_start_idx = i
            break
        end
    end

    if data_start_idx == 0
        error("Could not find Data Start! Searched lines after grid for row with $(num_Rs + 1) numbers.")
    end
    # println("Data starts at Line $data_start_idx.")

    ## pre-allocating memory 

    logTs = Vector{Float64}(undef, num_Ts)
    data = Vector{Float64}(undef, num_Rs * num_Ts)


    for T_idx in 1:num_Ts
        current_line_idx = data_start_idx + T_idx - 1 
        if current_line_idx > length(lines)
            error("Unexpected End of File at line $line_idx")
        end
        # Format: "3.750  -0.6258  -0.7084 ..."
        row_parts = split(lines[current_line_idx])
        #Safety check
        if length(row_parts) != (num_Rs + 1)
            error("Format Error at line $current_line_idx: Expected $(num_Rs+1) columns, found $(length(row_parts))")
        end
        logTs[T_idx] = parse(Float64, row_parts[1]) #Storing the Temperature

        for R_idx in 1:num_Rs
            val = parse(Float64, row_parts[R_idx + 1])
            flat_idx = R_idx + (T_idx - 1) * num_Rs
            data[flat_idx] = val 

        end
    end 
    return RT_table_opacity(X_val, Z_val, logTs, logRs, data)
end 

function Opacity_table_collector(directory :: String, mixture_key :: String)
    files = readdir(directory; join= true)
    opacity_files = filter(f -> endswith(f, ".data"), files) # A safety check if the directory contains other things

    if isempty(opacity_files)
        error("No .data files found in $directory")
    end

    temp_tables = RT_table_opacity[]

    for file in opacity_files
        if occursin(mixture_key, basename(file))
            try
                push!(temp_tables, RT_table_opacity(file))
            catch e 
                println("Warning: Skipping file $file due to error: $e")
            end
        end
    end

    unique_Xs = sort(unique([t.X for t in temp_tables]))
    unique_Zs = sort(unique([t.Z for t in temp_tables]))

    grid = Matrix{RT_table_opacity}(undef, length(unique_Xs), length(unique_Zs))

    for t in temp_tables
        x_idx = searchsortedfirst(unique_Xs, t.X)
        z_idx = searchsortedfirst(unique_Zs, t.Z)
        grid[x_idx, z_idx] = t 
    end 
    println("Opacity tables loaded")
    return Opacity_table_collector(unique_Xs, unique_Zs, grid)

end 

@inline function damp_slope(x::T) where T
        if x > 1.0
            # Start at edge (1.0) + tiny fraction of the excess
            return 1.0 + 0.5 * (x - 1.0)
        elseif x < 0.0
            # Start at edge (0.0) + tiny fraction of the deficit
            return 0.0 + 0.5 * x
        else
            return x
        end
    end


# function get_log_kappa_per_table(table::RT_table_opacity, val_logT::Float64, val_logR::Float64, logT::T, logR::T) where T
 
#     # The following part was for clamping the value of T,R at the edges so that value of κ outside the boundary is constant and stable, but this constant value 
#     # causes issues with duals, so now the value is extrapolated smoothly outside the boundary using a very small slope
 
#     # # Clamping data 
#     # min_T, max_T = table.logTs[1], table.logTs[end]
#     # min_R, max_R = table.logRs[1], table.logRs[end]

#     # eff_logT = clamp(logT, min_T, max_T)
#     # eff_logR = clamp(logR, min_R, max_R)

#     T_idx = searchsortedlast(table.logTs, val_logT)
#     R_idx = searchsortedlast(table.logRs, val_logR)

#     ## Boundary cases (For now clamping to the edge)
#     T_idx = clamp(T_idx, 1, length(table.logTs)-1)
#     R_idx = clamp(R_idx, 1, length(table.logRs)-1)

#     #logR and logT values
#     R0 = table.logRs[R_idx]
#     R1 = table.logRs[R_idx+1]
#     T0 = table.logTs[T_idx]
#     T1 = table.logTs[T_idx+1]

#     # Calculating the slopes 
#     u = (logR - R0)/ (R1 - R0)
#     v = (logT - T0)/ (T1 - T0)

    

#     u_damped = damp_slope(u)
#     v_damped = damp_slope(v)
    
    
#     # Now we will look into the 1D array since we have the indexes

#     N_R = length(table.logRs)

#     # calculating the kappa values 
#     #Bottom left value for interpolation
#     idx_00 = R_idx + (T_idx - 1) * N_R
#     val_00 = table.kap_data[idx_00]
    
#     #Bottom right value for interpolation
#     idx_10 = (R_idx + 1)+ (T_idx - 1) * N_R
#     val_10 = table.kap_data[idx_10]

#     #Top left value for interpolation
#     idx_01 = R_idx + T_idx * N_R
#     val_01 = table.kap_data[idx_01]

#     #Top right value for interpolation
#     idx_11 = (R_idx + 1) + T_idx * N_R
#     val_11 = table.kap_data[idx_11]


#     #Bilinear Interpolation
#     # f(x,y) = (1-u)(1-v)f00 + u(1-v)f10 + (1-u)v f01 + uv f11
#     w00 = (1 - u_damped) * (1 - v_damped)
#     w10 = u_damped * (1 - v_damped)
#     w01 = (1 - u_damped) * v_damped
#     w11 = u_damped * v_damped

#     interpolated_log_kappa = (w00 * val_00) + (w10 * val_10) + 
#                              (w01 * val_01) + (w11 * val_11)

#     return interpolated_log_kappa

# end




# The following function implements Bicubic interpolation using Catmull-Rom Cubic Basis Functions (for unit square)

@inline function cubic_weights(t :: T) where T
    t2 = t * t
    t3 = t2 * t 

    # Standard Catmull-Rom weights: [1 t t^2 t^3] M_SH (M_SH is the  Catmull-Rom Cubic Basis matrix for 1D)
    # https://en.wikipedia.org/wiki/Catmull–Rom_spline
    w_m1 = -0.5 * t3 + 1.0 * t2 - 0.5 * t
    w_0  =  1.5 * t3 - 2.5 * t2 + 1.0
    w_1  = -1.5 * t3 + 2.0 * t2 + 0.5 * t
    w_2  =  0.5 * t3 - 0.5 * t2
    
    return w_m1, w_0, w_1, w_2
end

function get_log_kappa_per_table(table::RT_table_opacity, val_logT::Float64, val_logR::Float64, logT::T, logR::T) where T
    # Get table boundaries
    min_T, max_T = table.logTs[1], table.logTs[end]
    min_R, max_R = table.logRs[1], table.logRs[end]

    # Clamp the Dual Numbers (Preserves derivatives inside, zeros them outside)
    eff_logT = clamp(logT, min_T, max_T)
    eff_logR = clamp(logR, min_R, max_R)

    # Also clamp the Float values for index searching
    eff_val_logT = clamp(val_logT, min_T, max_T)
    eff_val_logR = clamp(val_logR, min_R, max_R)

    i_T = searchsortedlast(table.logTs, eff_val_logT)
    i_R = searchsortedlast(table.logRs, eff_val_logR)
    
    #length of the tables
    N_T = length(table.logTs)
    N_R = length(table.logRs)

    #clamping the indexes ()
    i_T = clamp(i_T, 1, N_T - 1)
    i_R = clamp(i_R, 1, N_R - 1)

    T_min, T_max = table.logTs[i_T], table.logTs[i_T+1]
    R_min, R_max = table.logRs[i_R], table.logRs[i_R+1]

    u = (eff_logT - T_min) / (T_max - T_min)
    v = (eff_logR - R_min) / (R_max - R_min)

    wt_m1, wt_0, wt_1, wt_2 = cubic_weights(u)
    wr_m1, wr_0, wr_1, wr_2 = cubic_weights(v)

    weights_T = (wt_m1, wt_0, wt_1, wt_2)
    weights_R = (wr_m1, wr_0, wr_1, wr_2)

    log_kappa = zero(T)

    for (idx_T, offset_T) in enumerate(-1:2)
        for (idx_R, offset_R) in enumerate(-1:2)
            curr_T_idx = i_T + offset_T
            curr_R_idx = i_R + offset_R

            real_idx_T = clamp(curr_T_idx, 1, N_T)
            real_idx_R = clamp(curr_R_idx, 1, N_R)

            flat_idx = real_idx_R + (real_idx_T - 1) * N_R
            
            val = table.kap_data[flat_idx]

            # --- ACCUMULATION ---
            # Add: Value * Weight_T * Weight_R
            log_kappa += val * weights_T[idx_T] * weights_R[idx_R]
        end
    end 
    return log_kappa
end 



function get_opacity_table_collection(collection :: Opacity_table_collector, lnT::TT, lnρ::TT, xa::AbstractVector{<:TT}, species::Vector{Symbol})::TT where {TT<:Real}
    
    inv_ln10 = 0.4342944819032518
    logT = lnT * inv_ln10
    logρ = lnρ * inv_ln10
    logR = logρ - 3*logT + 18
    
    val_logT = ForwardDiff.value(logT)
    val_logR = ForwardDiff.value(logR)

    # X_dual = zero(TT)
    # Y_dual = zero(TT)

    # for i in eachindex(species)
    #     name = species[i]
    #     # Sum Hydrogens (H1, H2) for X
    #     if name == :H1 
    #         X_dual += xa[i]
    #     # Sum Heliums (He3, He4) for Y
    #     elseif name == :He4
    #         Y_dual += xa[i]
    #     end
    # end
    iH1 = findfirst(==(:H1), species)
    X_dual = xa[iH1]
    iHe4 = findfirst(==(:He4), species)
    Y_dual = xa[iHe4]
    Z_dual = 1 - X_dual -Y_dual

    val_X = ForwardDiff.value(X_dual)
    val_Z = ForwardDiff.value(Z_dual)
    # println("Target Composition: X = $(round(val_X, digits=5)), Z = $(round(val_Z, digits=5))")
    Nx = length(collection.Xs)
    Nz = length(collection.Zs)

    ix = clamp(searchsortedlast(collection.Xs, val_X), 1, Nx-1)
    iz = clamp(searchsortedlast(collection.Zs, val_Z), 1, Nz - 1)

    # X0 = collection.Xs[ix]
    # X1 = collection.Xs[ix+1]
    # Z0 = collection.Zs[iz]
    # Z1 = collection.Zs[iz+1]
    
    # # wx = (X_dual - X0) / (X1 - X0)
    # # wz = (Z_dual - Z0) / (Z1 - Z0)

    # # wx = clamp(wx, 0.0, 1.0)
    # # wz = clamp(wz, 0.0, 1.0)

    # #Now using dampslope instead of clamping 
    # wx = damp_slope((X_dual - X0) / (X1 - X0))
    # wz = damp_slope((Z_dual - Z0) / (Z1 - Z0))


    # t00 = collection.tables[ix,   iz]   # Low X, Low Z
    # t10 = collection.tables[ix+1, iz]   # High X, Low Z
    # t01 = collection.tables[ix,   iz+1] # Low X, High Z
    # t11 = collection.tables[ix+1, iz+1] # High X, High Z

    # k00 = get_log_kappa_per_table(t00, val_logT, val_logR, logT, logR)
    # k10 = get_log_kappa_per_table(t10, val_logT, val_logR, logT, logR)
    # k01 = get_log_kappa_per_table(t01, val_logT, val_logR, logT, logR)
    # k11 = get_log_kappa_per_table(t11, val_logT, val_logR, logT, logR)

    # # Interpolate X first (at fixed Zs)
    # k_z0 = k00 * (1 - wx) + k10 * wx
    # k_z1 = k01 * (1 - wx) + k11 * wx

    # # Interpolate Z
    # log_kappa_final = k_z0 * (1 - wz) + k_z1 * wz

    """
    Bicubic interpolation
    """
    u_X = (X_dual - collection.Xs[ix]) / (collection.Xs[ix+1] - collection.Xs[ix])
    u_Z = (Z_dual - collection.Zs[iz]) / (collection.Zs[iz+1] - collection.Zs[iz])

    wX = cubic_weights(u_X)
    wZ = cubic_weights(u_Z)
    # @show wX, wZ, collection.Xs[ix], collection.Zs[iz]
    log_kappa_final = zero(TT)
    for (itx, offset_x) in enumerate(-1:2)
        for (itz, offset_z) in enumerate(-1:2)
            
            # Boundary handling for the table grid indices
            real_ix = clamp(ix + offset_x, 1, Nx)
            real_iz = clamp(iz + offset_z, 1, Nz)

            # Check if table exists (handles the X=1.0, Z>0 empty spots you mentioned)
            if !isassigned(collection.tables, real_ix, real_iz)
                #  if a corner of the stencil is missing, use the central table
                # In a real OPLIB grid, you'd want a more robust neighbor search here
                table = collection.tables[ix, iz]
            else
                table = collection.tables[real_ix, real_iz]
            end

            # Get log_kappa for this specific table at (T, R) using Bicubic
            k_val = get_log_kappa_per_table(table, val_logT, val_logR, logT, logR)

            # Accumulate using composition weights
            log_kappa_final += k_val * wX[itx] * wZ[itz]
        end
    end
    # # Final Safety Cap
    # limit = 50.0
    # val_final = ForwardDiff.value(log_kappa_final)
    # if val_final > limit
    #      # Limit + tiny slope
    #      log_kappa_final = limit + 1e-6 * (log_kappa_final - limit)
    # elseif val_final < -limit
    #      log_kappa_final = -limit + 1e-6 * (log_kappa_final + limit)
    # end

    return 10^log_kappa_final
end

# Smoothing function to interpolate smoothly between the two datasets 
# Using the function f(t)=3t^2 - 2t^3
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

function Jems.Opacity.get_opacity_resultsTρ(composite :: CompositeOpacity, lnT::TT, lnρ::TT, xa::AbstractVector{<:TT}, species::Vector{Symbol})::TT where {TT <: Real}
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
