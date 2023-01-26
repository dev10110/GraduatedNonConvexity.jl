```@meta
CurrentModule = GraduatedNonConvexity
```

# GraduatedNonConvexity

Documentation for [GraduatedNonConvexity](https://github.com/dev10110/GraduatedNonConvexity.jl).


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

See the example below on how to use the library:




## References
See the following paper for details on the method and implementation:
```
H. Yang, P. Antonante, V. Tzoumas, and L. Carlone,“Graduated Non-Convexity for Robust Spatial Perception: From Non-Minimal Solvers to Global Outlier Rejection”,IEEE Robotics and Automation Letters (RA-L), 2020.
```

```@index
```

```@autodocs
Modules = [GraduatedNonConvexity]
```
