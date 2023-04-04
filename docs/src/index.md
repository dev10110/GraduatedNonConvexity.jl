```@meta
CurrentModule = GraduatedNonConvexity
```

# GraduatedNonConvexity

Documentation for [GraduatedNonConvexity](https://github.com/dev10110/GraduatedNonConvexity.jl).

## Introduction

This library is built to solve optimization problems of the form

```math
\min_{x} \sum_{i=1}^N \rho_\mu ( r ( y_i, x) )
```
where $$\rho_\mu$$ is the surrogate of a (non-convex) robust cost function $$\rho$$. 

To do this, it instead solves two problems repeatedly in a loop:
 
```math
x^{k+1} \gets \operatorname{argmin}_x \sum_{i=1}^N w_i^{k} r^2(y_i, x)
```  
```math
w^{k+1} \gets \operatorname{argmin}_w \sum_{i=1}^N w_i r^2(y_i, x^{k+1}) + \Phi_\rho(w_i)
```

i.e., it solves the non-convex robust optimization as a sequence of weighted least-squares problems. This method is particulary good when step 1 has a closed-form solution or a good solver, and step 2 has analytic solutions. 

This library assumes you have a function `LSQ_fn` that can perform the optimization in step 1. Then, it performs step 2 analytically, for some specific robust cost functions. 

The supported robust cost functions are:
- Geman McClure:
```math
\rho(r) = \frac{\bar{c}^2 r^2}{\bar{c}^2  + r^2}
```
- Truncated Least Squares:
```math
\rho(r) = \begin{cases}
r^2 & \text{if } r^2 \in [0, \bar{c}^2]\\
\bar{c}^2 & \text{if } r^2 \in [\bar{c}^2, \infty)
\end{cases}
```

## Quick Start
See the example below on how to use the library:

First, we construct the data
```@example main
using LinearAlgebra

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
```

Next, we define the weighted least squares solver, and the residual calculation functions

```@example main
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
```

And now solve:
```@example main
using GraduatedNonConvexity

# define maximum residual of inliers
c = 0.3

## Solve using GNC
β_gm = GNC_GM(N, data, least_sq_solver, residual_fn, c; verbose=false) # or
β_tls = GNC_TLS(N, data, least_sq_solver, residual_fn, c; verbose=false)
```

To save allocations, you can also provide in-place versions of the least squares solvers and the residual functions:
```
function least_sq_solver!(x, w, data)
    ...
end 

function residual_fn!(res, β, data)
    ...
end
```

and you can call the library in place:
```
β_gnc = zeros(2) # provide the initial guess for the solution
GNC_GM!(β_gnc, N, data, least_sq_solver!, residual_fn!, c; verbose=false) # or
GNC_TLS!(β_gnc, N, data, least_sq_solver!, residual_fn!, c; verbose=false)
```

## References
See the following paper for details on the method and implementation:
```
H. Yang, P. Antonante, V. Tzoumas, and L. Carlone,
“Graduated Non-Convexity for Robust Spatial Perception: 
From Non-Minimal Solvers to Global Outlier Rejection”,
IEEE Robotics and Automation Letters (RA-L), 2020.
```

## Documentation
```@index
```

```@autodocs
Modules = [GraduatedNonConvexity]
```
