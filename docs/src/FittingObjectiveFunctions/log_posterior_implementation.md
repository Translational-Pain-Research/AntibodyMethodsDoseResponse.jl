# Logarithmic posterior probability: How to implement

Consider the data and model from [Simple-example](@ref):
```@example 1
using FittingObjectiveFunctions, Plots #hide

X = collect(1:10)
Y = [1.0, 1.78, 3.64, 3.72, 5.33, 2.73, 7.52, 9.19, 6.73, 8.95]
ΔY = [0.38, 0.86, 0.29, 0.45, 0.66, 2.46, 0.39, 0.45, 1.62, 1.54]
model = ModelFunctions((x,λ)-> λ*x)

nothing #hide
```

## Log-likelihood and log-posterior

The general procedure to obtain the log-likelihood and log-posterior functions is the same as described in [Posterior-probability:-How-to-implement](@ref). However, the distributions and the prior need to be in logarithmic form, such that the default distributions of the `FittingData` constructor do not work. 
Thus, we need to define a `FittingData` object with logarithmic distributions:
``` @example 1
using Distributions
data_log_dists = FittingData(X,Y,ΔY,distributions = (y,m,Δy)-> logpdf(Normal(m,Δy),y))
nothing #hide
```
The log-likelihood function can be obtained by using the `log_posterior_objective` function:
``` @example 1
log_likelihood = log_posterior_objective(data_log_dists,model)
nothing #hide
```
As described in [Using-priors](@ref), the likelihood is obtained by using the prior `λ-> 1` (or in the logarithmic case `λ-> 0`). Again, this is what happens in the background. To use prior, it just needs to be passed as third argument:
``` @example 1
log_posterior = log_posterior_objective(data_log_dists,model, λ-> logpdf(Normal(1,0.1),λ))
nothing #hide
```
The resulting functions can be compared by adjusting the constant offset (see [Logarithmic scale](@ref))
``` @example 1
large_scope = plot(x-> log_likelihood(x) + 1.105, xlims = [0.9,1.2], label = "log_likelihood", legend = :topleft) #hide
plot!(log_posterior, label = "log_posterior")  #hide
small_scope = plot(x-> log_likelihood(x) + 1.105, xlims = [1.065,1.085], legend = :none) #hide
plot!(log_posterior)  #hide
plot(large_scope,small_scope, layout = (1,2), size = (800,300)) #hide
```

## Application: Regularized least squares

Recall the weighted least squares objective from [LSQ:-How-to-implement](@ref)
``` @example 1
lsq = lsq_objective(data_log_dists, model)
nothing #hide
```
Note that the `distributions` field of `data_log_dists` has no effect on least squares objectives. 

To replicate the least squares objective, unnormalized logarithmic distributions can be used:
``` @example 1
data_lsq_dists = FittingData(X,Y,ΔY, distributions = (y,m,Δy)-> -(y-m)^2/Δy^2)
lsq_likelihood = log_posterior_objective(data_lsq_dists,model)
nothing #hide
```

``` @example 1
plot(lsq, label = "lsq") #hide
plot!(lsq_likelihood, label = "lsq_likelihood") #hide
```
The `lsq_likelihood` is the same function as `lsq`, but with the opposite sign (because it is a logarithmic unnormalized posterior probability density). This could either be fixed by using `λ -> -lsq_likelihood(λ)`, or in this case by using `distributions = (y,m,Δy)-> -(y-m)^2/Δy^2`.

Now, a regularization can be implemented by using a corresponding logarithmic prior:
``` @example 1
lsq_posterior = log_posterior_objective(data_lsq_dists,model, λ -> - λ^2)
nothing #hide
```

## Derivatives

Analytical derivatives can be obtained almost in the same way as described in [LSQ:-partial-derivatives-and-gradients](@ref):

### Partial derivatives

``` @example 1
@doc log_posterior_partials #hide
```

### Gradient

``` @example 1
@doc log_posterior_gradient #hide
```