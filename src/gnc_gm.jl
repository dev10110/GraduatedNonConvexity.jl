
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


"""
    GNC_GM(N, data, LSQ_fn, RES_fn, c)

Perform robust least squares given `data` using the German-McClure cost function. 

Inputs:
- `N`: the number of data points (equal to the number of residuals)
- `data`: the data which is to be fit
- `LSQ_fn`: Assumes `LSQ_fn(w, data)` returns the weighted least squares solution using the weights `w`
- `RES_fn`: Assumes `RES_fn(x, data)` returns a vector with the residuals given the candidate solution `x`. 
- `c`: the maximum residual of an inlier

Parameters:
- `max_iterations=1000`
- `μ_factor=1.4`: factor to decrease μ by each iteration.
- `verbose=true`
- `rtol=1e-6`

Terminates when `rtol` is reached, or `μ=1`.

"""
function GNC_GM(N, data, LSQ_fn, RES_fn, c̄;
    max_iterations = 1000,
    μ_factor = 1.4,
    verbose=true,
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


"""
    GNC_GM!(x, w, rs, data, LSQ_fn!, RES_fn!, c)

In-place version of `GNC_GM`

Inputs:
- `x`: starting point for the robust least squares
- `w`: initial weights vector (data will be overriden) 
- `rs`: residuals vector (data will be overriden). Both `w` and `rs` need to be of the same length. 
- `data`: the data which is to be fit
- `LSQ_fn!`: Assumes `LSQ_fn!(x, w, data)` updates `x` in-place with the weighted least squares solution using weights `w`
- `RES_fn!`: Assumes `RES_fn!(rs, x, data)` updates `rs` in-place with the returns residuals given the candidate solution `x`. 
- `c`: the maximum residual of an inlier
Parameters:
- (same as `GNC_GM`)


"""
function GNC_GM!(x, w, rs, data, LSQ_fn!, RES_fn!, c̄;
    max_iterations = 1000,
    μ_factor = 1.4,
    verbose=true,
    rtol=1e-6)
   
    # obtain the unweighted solution
    RES_fn!(rs, x, data)
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
        if μ == 1
            if verbose
                print("μ == 1. Done. \n")
            end
            break
        else
            μ = max(1, μ / μ_factor)
        end
        
    end
    
    return
end
