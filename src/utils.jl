
function rmax_rsum(r::V, w::V) where {V}
    rmax = mapreduce(abs, max, r)
    rsum = sum(w[i] * r[i]^2 for i=1:length(r))
    return rmax, rsum
end

