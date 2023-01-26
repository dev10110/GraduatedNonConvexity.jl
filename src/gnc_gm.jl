
## GM: μ: Inf -> 1

function ρ_GM(r, μ, c̄)
    if μ == Inf
        return r^2
    end
    
    return μ * c̄^2 * r^2 / (μ * c̄^2 + r^2)
end


function Φ_GM(w, μ, c̄)
    return μ*c̄^2 * (sqrt(w) - 1)^2
end


function weight_update_GM(r̂, μ, c̄)
    return ((μ * c̄^2) / (r̂^2 + μ * c̄^2) )^2
end


function rmax_rsum(r::V, w::V) where {F, V<: Vector{F}}
    rmax = zero(F)
    rsum = zero(F)
    for i=1:length(r)
        rmax = max(rmax, abs(r[i]))
        rsum += w[i] * r[i]^2
    end
    return rmax, rsum
end

function GNC_GM(N, data, LSQ_fn, RES_fn, c̄;
    max_iterations = 1000,
    μ_factor = 1.4,
    verbose=false,
    rtol=1e-6)
    
    # obtain the unweighted solution
    w = ones(N)
    x = LSQ_fn(w, data)
    rs = RES_fn(x, data)
    rmax, rsum = rmax_rsum(rs, w)
    
    if verbose 
        @printf " ## Graduated Nonconvexity Solver \n"
        @printf " ## Devansh Agrawal, 2023 \n\n"
        @printf "    iter |           μ |        rmax |        rsum\n"
        @printf "%8i |         Inf | %9.5e | %9.5e\n" 0 rmax rsum
    end

    # initialize μ
    μ = 2 * rmax^2 / c̄^2
    
    for iter = 1:max_iterations
        
        # update weights
        for i=1:length(rs)
            w[i] = weight_update_GM(rs[i], μ, c̄)
        end
        
        # solve LSQ
        x = LSQ_fn(w, data)
        rs = RES_fn(x, data)
        rmax, rsum_new = rmax_rsum(rs, w)
        
        # check for convergence
        if iter > 1 && abs(rsum_new - rsum) <= rtol
            if verbose
                println("rsum <= rtol. Done. \n")
            end
            break
        end
        rsum = rsum_new
        
        if verbose
            @printf "%8i | %9.5e | %9.5e | %9.5e\n" iter μ rmax rsum
        end
        
        # update μ
        if μ == 1
            if verbose
                print("μ == 1. Done. \n")
            end
            break
        else
            μ = max(1, μ/ μ_factor)
        end
        
    end
    
    return x
end


