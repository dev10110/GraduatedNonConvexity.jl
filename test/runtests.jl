using GraduatedNonConvexity
using Test
using Random
using LinearAlgebra

# make sure it always runs the same test
Random.seed!(1)

@testset "GraduatedNonConvexity.jl - 2% outliers, no noise" begin


    ## construct dataset
    N=1000
    
    # this is the ground truth parameter we are trying to solve for
    β_gt = randn(2)
    
    x = randn(N, 2)
    x[:,2] .= 1
    
    # create some noise-free measurements
    y = x * β_gt
    
    # add some outliers to 2% of data
    for i=1:N
        rand() < 0.02 ? y[i] += 1.0 + rand()  : continue
    end
    
    # collect all the necessary data
    data = (x, y)
    
    ## define weighted least squares solver
    function least_sq_solver(w, data)
        
        X = data[1]
        y = data[2]
        
        W = diagm(w)
    
        # analytic expression for weighted least squares
        return (X' * W * X) \ (X' * W * y)
    end
    
    ## define residual function
    function residual_fn(β, data)
      X = data[1]
      y = data[2]
      return y - X * β
    end
    
    # define maximum residual of inliers
    c = 0.3
    
    ## Solve using GNC
    β_gm = GNC_GM(N, data, least_sq_solver, residual_fn, c) # or
    β_tls = GNC_TLS(N, data, least_sq_solver, residual_fn, c)
    
    ##  compare performance
    #@show norm(β_gm - β_gt)
    #@show norm(β_tls - β_gt)
    
    @test norm(β_gm - β_gt) <= 1e-4
    @test norm(β_tls - β_gt) <= 1e-4
end

@testset "GraduatedNonConvexity.jl - 20% outliers" begin

    ## construct dataset
    N=1000
    
    # this is the ground truth parameter we are trying to solve for
    β_gt = randn(2)
    
    x = randn(N, 2)
    x[:,2] .= 1
    
    # create some noisy measurements
    y = x * β_gt + 0.2*(2*rand(N) .-1) 
    
    # add some outliers to 20% of data
    for i=1:N
        rand() < 0.2 ? y[i] += 1.0 + rand()  : continue
    end
    
    # collect all the necessary data
    data = (x, y)
    
    ## define weighted least squares solver
    function least_sq_solver(w, data)
        
        X = data[1]
        y = data[2]
        
        W = diagm(w)
    
        # analytic expression for weighted least squares
        return (X' * W * X) \ (X' * W * y)
    end
    
    ## define residual function
    function residual_fn(β, data)
      X = data[1]
      y = data[2]
      return y - X * β
    end
    
    # define maximum residual of inliers
    c = 0.3
    
    ## Solve using GNC
    β_gm = GNC_GM(N, data, least_sq_solver, residual_fn, c) # or
    β_tls = GNC_TLS(N, data, least_sq_solver, residual_fn, c)
    
    ##  compare performance
    #@show norm(β_gm - β_gt)
    #@show norm(β_tls - β_gt)
    
    @test norm(β_gm - β_gt) <= 1e-2
    @test norm(β_tls - β_gt) <= 1e-2
end
