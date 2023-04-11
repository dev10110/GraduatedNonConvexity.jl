
## TLS: μ: 0 -> Inf

function Φ_TLS(w, μ, c̄)
    return (μ * (1 - w) * c̄^2) / (μ + w)
end

function ρ_TLS(r, μ, c̄)
    
    if μ == 0
        return r^2
    end
    
    rsq = r^2
    
    if rsq <= μ/(μ+1) * c̄^2
        return rsq
    elseif rsq <= (μ+1)/μ c̄^2
        return 2 * c̄ * abs(r) * sqrt(μ*(μ+1)) - μ*(c̄^2 + rsq)
    else
        return c̄^2
    end
end

function weight_update_TLS(r̂, μ, c̄)
    
    r̂sq = r̂^2
    
    if r̂sq <= (μ / (μ + 1 )) * c̄^2
        return 1
    elseif r̂sq <= (μ + 1) * c̄^2 / μ
        return c̄ * sqrt(μ * (μ + 1)) / abs(r̂) - μ
    else
        return 0
    end

end

"""
    GNC_TLS(N, data, LSQ_fn, RES_fn, c)

Perform robust least squares given `data` using the Truncated Least Squares cost function. 

Inputs:
- `N`: the number of data points (equal to the number of residuals)
- `data`: the data which is to be fit
- `LSQ_fn`: Assumes `LSQ_fn(w, data)` returns the weighted least squares solution using the weights `w`
- `RES_fn`: Assumes `RES_fn(x, data)` returns a vector with the residuals given the candidate solution `x`. 
- `c`: the maximum residual of an inlier

Parameters:
- `max_iterations=1000`
- `μ_factor=1.4`: factor to increase μ by each iteration.
- `verbose=true`
- `rtol=1e-6`

Terminates when `rtol` is reached.

"""
function GNC_TLS(N, data, LSQ_fn, RES_fn, c̄;
    max_iterations = 1000,
    μ_factor = 1.4,
    verbose=true,
    rtol = 1e-6)

    # obtain the unweighted solution
    w = ones(N)
    x = LSQ_fn(w, data)
    rs = RES_fn(x, data)
    rmax, rsum = rmax_rsum(rs, w)
    
    if verbose 
        @printf " ## Graduated Nonconvexity Solver \n"
        @printf " ## Devansh Agrawal, 2023 \n\n"
        @printf "    iter |           μ |        rmax |        rsum \n"
        @printf "---------|-------------|-------------|-------------\n"
        @printf "%8i | %9.5e | %9.5e | %9.5e\n" 0 0 rmax rsum
    end

    # initialize μ
    μ = c̄^2 / ( max( 1e-12, 2 * rmax^2  - c̄^2) )
    
    for iter = 1:max_iterations
        
        # update weights
        for i=1:length(rs)
            w[i] = weight_update_TLS(rs[i], μ, c̄)
        end
        
        # solve LSQ
        x = LSQ_fn(w, data)
        rs = RES_fn(x, data)
        rmax, rsum_new = rmax_rsum(rs, w)
        
        # check for convergence
        if iter > 1 && abs(rsum_new - rsum) <= rtol
            if verbose
                println("Δrsum <= rtol. Done. \n")
            end
            break
        end
        rsum = rsum_new
        
        if verbose
            @printf "%8i | %9.5e | %9.5e | %9.5e\n" iter μ rmax rsum
        end
        
        # update μ
        μ  = μ * μ_factor
        
    end
    
    return x
end

"""
    GNC_TLS!(x, w, rs, data, LSQ_fn!, RES_fn!, c)

In-place version of `GNC_TLS`

Inputs:
- `x`: starting point for the robust least squares
- `w`: initial weights vector (data will be overriden) 
- `rs`: residuals vector (data will be overriden). Both `w` and `rs` need to be of the same length. 
- `data`: the data which is to be fit
- `LSQ_fn!`: Assumes `LSQ_fn!(x, w, data)` updates `x` in-place with the weighted least squares solution using weights `w`
- `RES_fn!`: Assumes `RES_fn!(rs, x, data)` updates `rs` in-place with the returns residuals given the candidate solution `x`. 
- `c`: the maximum residual of an inlier
Parameters:
- (same as `GNC_TLS`)


"""
function GNC_TLS!(x, w, rs, data, LSQ_fn!, RES_fn!, c̄;
    max_iterations = 1000,
    μ_factor = 1.4,
    verbose=true,
    rtol = 1e-6)

    # obtain the unweighted solution
    RES_fn!(rs, x, data)
    rmax, rsum = rmax_rsum(rs, w)
    
    if verbose 
        @printf " ## Graduated Nonconvexity Solver \n"
        @printf " ## Devansh Agrawal, 2023 \n\n"
        @printf "    iter |           μ |        rmax |        rsum \n"
        @printf "---------|-------------|-------------|-------------\n"
        @printf "%8i | %9.5e | %9.5e | %9.5e\n" 0 0 rmax rsum
    end

    # initialize μ
    μ = c̄^2 / (2 * rmax^2  - c̄^2)
    
    for iter = 1:max_iterations
        
        # update weights
        for i=1:length(rs)
            w[i] = weight_update_TLS(rs[i], μ, c̄)
        end
        
        # solve LSQ
        LSQ_fn!(x, w, data)
        RES_fn!(rs, x, data)
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
        μ  = μ * μ_factor
        
    end
    
    return
end
