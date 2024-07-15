"""
    struct RateCache{TT}

structure that saves frequently used powers/functions of temperature values, so as to avoid repeated calculations
"""
mutable struct RateCache{TT}
    T9::TT
    logT9::TT  # natural log
    cbrtT9::TT
    RateCache{TT}() where {TT} = new(zero(TT), zero(TT), zero(TT))
end


function update_rate_cache!(cache::RateCache{TT}, Temp::TT) where {TT}
    cache.T9 = Temp / 1e9
    cache.logT9 = log(cache.T9)
    cache.cbrtT9 = cbrt(cache.T9)
end