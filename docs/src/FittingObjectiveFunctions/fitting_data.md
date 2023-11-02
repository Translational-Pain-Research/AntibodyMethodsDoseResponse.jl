# FittingData and ModelFunctions

```@example 1
using FittingObjectiveFunctions #hide
model = ModelFunctions((x,λ) -> λ*x) #hide
X = collect(1:10) #hide
Y = [1.0, 1.78, 3.64, 3.72, 5.33, 2.73, 7.52, 9.19, 6.73, 8.95] #hide
ΔY = [0.38, 0.86, 0.29, 0.45, 0.66, 2.46, 0.39, 0.45, 1.62, 1.54] #hide
nothing #hide
```

In the [Simple-example](@ref), it was already mentioned that the data needs to be summarized in a `FittingData` object and that information about the model needs to be summarized in a `ModelFunctions` object. The details about these data types are discussed here.

## FittingData
`FittingData` objects have the following general constructor:

```julia
data = FittingData(X,Y,ΔY, distributions = distribution_functions)
```
The resulting object has the following fields:
```julia
data.independent == X
data.dependent == Y
data.errors == ΔY
data.distributions == distributions_functions
```


!!! tip "Tip: shortened constructors"
	If no measurement errors `ΔY` are available, one can use the shortened constructor
	```julia
	FittingData(X,Y, distributions = distribution_functions)
	```
	In this case, the default errors `ones(length(X))` are used (leading e.g. to the standard least squares objective function: [Background:-LSQ](@ref))

!!! tip "Tip: distributions"
	The optional `distributions` keyword is used to specify the likelihood distributions (see [Background:-Posterior-probability](@ref)). In case distributions are not specified, the constructor defaults to normal distributions.

	The distributions can be specified as an array of functions, or a single function (if the same distribution shall be used for all data points) with the signature `(y,m,Δy)`, where

	* `y` is the measured dependent variable
	* `Δy` is the corresponding error
	* `m` are values the model function returns, when the parameters are varied, i.e. `m(x,λ)`.
	


## ModelFunctions

`ModelFunctions` objects have the following general constructor:
```julia
model = ModelFunctions(model_function, partial_derivatives = [derivative_functions...])
```
The resulting object has the following fields:
```julia
model.model == model_function
model.partials == [derivative_functions...]
```

The model function (and the partial derivatives) must have the signature `(x,λ)`, where `x` is the independent variable and `λ` is the parameter(array).

!!! tip "Tip: Partial derivatives"
	The keyword `partials` is optional, but is required for analytical partial derivatives and analytical gradient functions. The partial derivatives w.r.t. to the parameter are defined as array of functions.
	```math
	\text{m}(x,\lambda) = \lambda_1 x + \lambda_2 \ ,\quad \frac{\partial \text{m}(x,\lambda)}{\partial \lambda_1}  = x\ , \quad \frac{\partial \text{m}(x,\lambda)}{\partial \lambda_1} = 1
	```
	```julia
	ModelFunctions((x,λ)-> λ[1]*x + λ[2], partials = [(x,λ)-> x, (x,λ)-> 1])
	```

## Additional remarks

When a `FittingObject`/`ModelFunctions` object is created, some rudimentary consistency checks are made, e.g. that all arrays have the same lengths.
```@repl 1
FittingData([1,2,3],[1,2],[1,1,1])
```
Since the objects are mutable, the same consistency checks are repeated before objective functions are created. However, it can be useful to make a comprehensive consistency check. For this, an exemplary parameter (array) is needed:
```@example 1
data = FittingData(X,Y)
model = ModelFunctions((x,λ)-> λ*x)
consistency_check(data,model,1)
```
If everything works, `nothing` is returned, i.e. nothing happens. However, in case of a problem, an error is thrown:
```@repl 1
consistency_check(data,model,"1")
```

!!! note "Mutability of objects"
	Both `FittingData` objects and `ModelFunction` objects are mutable for convenience, as they are not performance relevant (only their fields are). However, when objective functions are created, the object fields are copied and enclosed in the objective function, to avoid accidental mutation. 
	
	I.e. once an objective function is created, it does not change, even if the `FittingData`/`ModelFunctions` objects are changed. To apply changes, a new objective function has to be created. 