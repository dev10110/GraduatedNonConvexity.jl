# GraduatedNonConvexity

<!-- [![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://dev10110.github.io/GraduatedNonConvexity.jl/stable/) -->
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://dev10110.github.io/GraduatedNonConvexity.jl/dev/)
[![Build Status](https://github.com/dev10110/GraduatedNonConvexity.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/dev10110/GraduatedNonConvexity.jl/actions/workflows/CI.yml?query=branch%3Amain)


## Quickstart

Suppose we are given data `x, y`

```
data = (x, y)
```

Suppose the model to fit is `y = x * m`, where we want to find `m`. 


Implement a weighted least squares solver 
```
function solver(w, data)
  # implement your weighted least squares optimizer
  # w is a N-vector of weights 
end
```

and a function to calculate the residuals
```
function residuals(c, data)
  # implement a function that returns the residuals for each data point given the candidate solution m
  # should return a N-vector of residuals
end
```

define maximum residual of inliers
```
c = 0.2
```

Solve using GNC
```
using GraduatedNonConvexity
m_gm = GNC_GM(N, data, solver, residuals, c) # or
m_tls = GNC_TLS(N, data, solver, residuals, c)
```

## Complete Example
see `examples/1d.jl`

## Documentation

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://dev10110.github.io/GraduatedNonConvexity.jl/dev/)
