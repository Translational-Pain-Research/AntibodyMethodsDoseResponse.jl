# LSQ: How to implement

Consider the data and model from [Simple-example](@ref):
```@example 1
using FittingObjectiveFunctions, Plots #hide

X = collect(1:10)
Y = [1.0, 1.78, 3.64, 3.72, 5.33, 2.73, 7.52, 9.19, 6.73, 8.95]
ΔY = [0.38, 0.86, 0.29, 0.45, 0.66, 2.46, 0.39, 0.45, 1.62, 1.54]
data = FittingData(X,Y,ΔY)
model = ModelFunctions((x,λ)-> λ*x)

nothing #hide
```

## LSQ: objective functions

The uncertainty-weighted least squares objective function can be obtained as follows:
```@example 1
weighted_lsq = lsq_objective(data,model)
nothing # hide
```
`lsq_objective` returns a function that takes the model parameters `λ` as argument.

To obtain the standard least squares objective, the errors must be set to `1`. Recall the shortened constructor (see [FittingData](@ref)):

``` @example 1
data_no_errors = FittingData(X,Y)
standard_lsq = lsq_objective(data_no_errors,model)
w_lsq_plot = plot(weighted_lsq, label = "weighted lsq", xlims = [0,2], xlabel = "λ", ylabel = "lsq") #hide
plot!(standard_lsq, label = "standard lsq", xlims = [0,2], xlabel = "λ", ylabel = "lsq", color = 2) #hide
```


## LSQ: partial derivatives and gradients

### Examples

Redefine `model` to obtain analytical derivatives (see [ModelFunctions](@ref)):

``` @example 1
model = ModelFunctions((x,λ)->λ*x , partials = [(x,λ)-> x])
∂_weighted_lsq = lsq_partials(data,model)
```
``` @example 1
∂_weighted_lsq[1](1.1)
```
Note that `lsq_partials` returns a vector of abstract functions with `λ` as argument (one for each partial derivative), even in the 1-dimensional case.

``` @example 1
∇_weighted_lsq = lsq_gradient(data,model)
```
``` @example 1
∇_weighted_lsq([0.0],1.1) 
```
On the other hand `lsq_gradient` directly returns the gradient function, but with a different signature `(grad_vector,λ)`. The argument `grad_vector` must be a vector of appropriate type and length, that can be mutated.

!!! info "Why mutation for gradient function?"
	In some optimization algorithms, the gradient function is called multiple times during each iteration. Mutating an array allows to reduce the memory allocation overhead of creating a new gradient array every time.

### Partial derivatives

``` @example 1
@doc lsq_partials #hide
```

### Gradient

``` @example 1
@doc lsq_gradient #hide
```