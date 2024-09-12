# [Fitting](@id fitting)

The model generators return a [`ModelFunctions`](https://translational-pain-research.github.io/FittingObjectiveFunctions-documentation/API/#FittingObjectiveFunctions.ModelFunctions) object. The `model` field of a [`ModelFunctions`](https://translational-pain-research.github.io/FittingObjectiveFunctions-documentation/API/#FittingObjectiveFunctions.ModelFunctions) object is a pure Julia function, allowing to implement the model fitting from scratch.

```@example Fitting
using AntibodyMethodsDoseResponseConvenience
model, params = accumulation_model(create_grid([1,2,3]))
typeof(model.model) <: Function
```


Yet, since the model generators create a [`ModelFunctions`](https://translational-pain-research.github.io/FittingObjectiveFunctions-documentation/API/#FittingObjectiveFunctions.ModelFunctions) object, it is convenient to construct the fitting objective with [`FittingObjectiveFunctions.jl`](https://translational-pain-research.github.io/FittingObjectiveFunctions-documentation/). Then, only the minimization/maximization of the objective function remains to be implemented. Because these steps are always the same, the [`adaptive_dose_response_fit`](@ref) function summarizes the creation of the model function and the objective function, requiring only the implementation of a function minimizer.

!!! tip "Reminder: Convenience workflow"
	If there is no reason to avoid the dependencies of [`AntibodyMethodsDoseResponseConvenience.jl`](https://github.com/Translational-Pain-Research/AntibodyMethodsDoseResponseConvenience.jl), the workflow as described in the [quick start guide](@ref quick_start) should be used. [`FittingCondition`](@ref) and [`fit_condition`](@ref) expose the same options that are described here.

## The setting

We consider the following measurement data `(concentrations, responses, errors)`:
```@example Fitting
adaptive_result, data, replicates =  load_results("examples/default_concentrations") # hide
concentrations = data.independent # hide
responses = data.dependent # hide
errors = data.errors # hide
scatter(concentrations,responses, yerror = errors, xaxis = :log, legend = :none)
```

## Simple model fitting
 
 The data needs to be summarized in a [`FittingData`](https://translational-pain-research.github.io/FittingObjectiveFunctions-documentation/API/#FittingObjectiveFunctions.FittingData) object, as described in [Models](@ref binding_models):

```@example Fitting
data = FittingData(concentrations, responses, errors)
```
Next, a `OneDimGrid` needs to be created, ideally covering the concentration range:

```@example Fitting
grid = create_grid(LogRange(1e-10,1e-2,40))
```

Finally, a function minimizer needs to be implemented. For this, we use [`Optim.jl`](https://julianlsolvers.github.io/Optim.jl/stable/), here:

```@example Fitting
function minimizer(f,∇f,init)
	lower = zeros(length(init))
	upper = [Inf for i in 1:length(init)]
	return optimize(f,lower,upper, init, Fminbox(NelderMead()),
			Optim.Options(g_tol = 1e-12, iterations =2000)).minimizer
end
```
The implemented `minimizer` must take the objective function `f`, its gradient function `∇f`, if applicable, and an initial parameter array `init` as arguments and return the minimizing parameters. Furthermore, since ``K_\tau \geq 0``, the optimization domain should be limited.

Now, a model can be fitted to the `data` with [`adaptive_dose_response_fit`](@ref):
```julia
result = adaptive_dose_response_fit(grid,data,
		minimizer, 
		options = AdaptiveOptions(model = accumulation_model)
	)
```
```@example Fitting
result =  deserialize("examples/direct_fit/results.jld") # hide
```
`adaptive_dose_response_fit` returns an [`AdaptiveResult`](@ref) object, that has the following fields:

* `optimizer`: The estimated parameters (result of model fitting).
* `objective_value`: The objective function value for the estimated parameters.
* `grid`: A grid containing the estimated parameters as grid weights.
* `result`: The [`DoseResponseResult`](@ref) object corresponding to the grid (and the offset parameter).
* `time`: The elapsed time for the model fit (in seconds).

```@example Fitting
scatter(data, xaxis = :log, legend = :topleft, label = "data")
plot!(result.result, label = "fit result")
```

```@example Fitting
plot(DensityPlot(result.grid), xaxis = :log, color = 2, fill = 0, 
	fillalpha = 0.5, label = "fit result")
```

## [Adaptive model fitting](@id adaptive_model_fitting)

The model fit above can be improved in two areas. First, a regularization could be used. Second, the adaptive density approximation from [`AdaptiveDensityApproximation.jl`](https://translational-pain-research.github.io/AdaptiveDensityApproximation-documentation/) could be used to reduce the number of parameters. Here, we recreate the default optimization from the [`AntibodyMethodsDoseResponseConvenience.jl`](https://github.com/Translational-Pain-Research/AntibodyMethodsDoseResponseConvenience.jl) package as described in the [quick start guide](@ref quick_start) to illustrate some of the available options.

### Setting up the objective function properties

For the regularization, a log-posterior objective can be used, where a smoothing prior defines the regularization. 

Log-posterior objectives differ from posterior objectives only by taking the logarithm of all functions. Mathematically, there is no difference, but for the computation of tiny probabilities, taking the logarithm upfront is numerically beneficial. Accordingly, the prior also needs to be defined as log-prior, i.e. as the logarithm of the prior.

 While [`FittingObjectiveFunctions.jl`](https://translational-pain-research.github.io/FittingObjectiveFunctions-documentation/) expects standard functions as prior/log-prior, [`adaptive_dose_response_fit`](@ref) requires a prior-generating function (Here, the `500` is used to reproduce the `scale = 500` from the [quick start guide](@ref quick_start)):

```@example Fitting
function log_prior_generator(centers, volumes, offset)
	ℓV = log.(centers .+ volumes/2) .- log.(centers .- volumes/2)
	if isnothing(offset)
		return λ ->  -500*(sum((λ[i]/ℓV[i]-λ[i+1]/ℓV[i+1])^2 for i in 1:length(λ)-1))/length(λ)^2
	else
		return λ -> -500*(sum((λ[i]/ℓV[i]-λ[i+1]/ℓV[i+1])^2 for i in 1:length(λ)-2) + λ[end]^2)/length(λ)^2
	end
end
```

The prior generator above will create a new log-prior function for each step of the adaptive fitting. The returned log-prior function reads:

```math
\text{log-prior}(\lambda) = -\frac{500}{\text{length}(\lambda)^2} \left(\lambda_{\text{offset}} + \sum_{i=1}^{n-1} \left( \frac{\lambda_i}{\ell V_i} - \frac{\lambda_{i+1}}{\ell V_{i+1}} \right)^2 \right)
```

Since the default density-visualization of grids rescales the weights by using the visual interval lengths in a logarithmic scale `ℓV` (see [Background: log-volume normalization](@ref log_volume_normalization)), the smoothing should be applied to the rescaled parameters `λ`. Hence, the `λ[i]` are divided by `ℓV[i]`.

Next, observe that the log-prior is just the logarithm of a normal distribution (up to a missing normalization) for the difference of the rescaled parameters. Essentially, the prior assumes that there is no difference between neighboring parameters where the scale ``\frac{500}{\text{length}(\lambda)^2}`` expresses the strength/importance of this assumption. This is the aforementioned smoothing.

!!! info "Why prior-generator functions?"
	Using log-prior generating functions seems unnecessarily complicated, at first. However, during the adaptive fit, the underlying grid approximation changes, leading to different (visual) interval lengths. Defining a fixed function for the prior could not take the change of interval lengths into account, i.e. the parameters could not be rescaled properly. Hence, the prior needs to be recalculated after every change of the grid, which requires a function that generates the prior from the grid properties.


Without additional information about the measurement errors, a normal distribution is a sensible choice for the uncertainty distribution. Since the goal is a log-posterior objective, the logarithm of a normal distribution `(y,m,Δy)-> -(y-m)^2/Δy^2` must be used:

```@example Fitting
data = FittingData(concentrations,responses, errors, distributions = (y,m,Δy)-> -(y-m)^2/Δy^2)
``` 

Here, `y` denotes the data point value, `m` denotes the model value calculated from the parameters and `Δy` denotes the measurement error. Uncertainty distributions must take the arguments in this order `(y,m,Δy)` and must return the distribution / log-distribution value.

### Setting up the minimizers

As before, a minimizer needs to be defined:

```@example Fitting
function minimizer(f,∇f,init)
	lower = zeros(length(init))
	upper = [Inf for i in 1:length(init)]
	return optimize(f,lower,upper, init, Fminbox(NelderMead()),
		Optim.Options(g_tol = 1e-12, iterations =2000)).minimizer
end
```


!!! info "Minimization and log-posterior?"
	Likelihood and posterior objectives usually need to be maximized. However, `Optim.jl` only provides minimizers, as do some other optimization packages, expecting from the user to flip the sign of the function for a maximization. [`adaptive_dose_response_fit`](@ref) flips the sign of the posterior and log-posterior objectives automatically.

But, to recreate the default fitting from the [quick start guide](@ref quick_start), a second minimizer is needed (using the `LBFGS` algorithm):
```@example Fitting
function post_minimizer(f,∇f,init)
	lower = zeros(length(init))
	upper = [Inf for i in 1:length(init)]
	return optimize(f,lower,upper, init, Fminbox(LBFGS()),
		Optim.Options(g_tol = 1e-12, iterations =2000)).minimizer
end
```

The second minimizer will be applied after the adaptive optimization to fine-tune the results with a gradient based minimizer. 

### Setting up the options

Before fitting the data, the fitting options ([`AdaptiveOptions`](@ref)) and the initial grid need to be defined. Among others, the objective and the prior-generator defined above need to be selected (the additional options in the example are needed to obtain the defaults from the [quick start guide](@ref quick_start)):

```@example Fitting
adaptive_options = AdaptiveOptions(objective = :log_posterior, 
	iterations = 30, 
	offset = eps(), 
	prior_generator = log_prior_generator
)
```

The idea of the adaptive approximation is to start with a coarse grid, containing only 2 intervals (3 interval edges):
```@example Fitting
grid = create_grid(LogRange(1e-10,1e-2,3)) 
``` 

!!! info "Adaptive fitting"
	With only 2 parameters, common optimizers find a good minimum even for suboptimal initial parameters. Then, the grid can be refined in regions of interest. Next, the previous result can be used as good initial point for the optimization. This process can be repeated several times, increasing the number intervals for regions of interest while not wasting computation time with small intervals for uninteresting regions.


### Fitting 

Now, the data can be fitted with [`adaptive_dose_response_fit`](@ref):
```julia
adaptive_result = adaptive_dose_response_fit(grid,data,minimizer, options = fitting_options)
```

To fine-tune the result with the gradient-based minimizer from above, [`adaptive_dose_response_fit`](@ref) needs to be called again with different options (i.e. no iterations to obtain a single fit and using the offset from the previous result):

```julia
post_options = AdaptiveOptions(objective = :log_posterior, 
	offset = result.optimizer[end], 
	prior_generator = log_prior_generator
)

adaptive_result = adaptive_dose_response_fit(adaptive_result.grid,data,
		post_minimizer, 
		options = post_options
	)
```


```@example Fitting
adaptive_result # hide
```

As before, the fitted dose-response curve describes the data well:
```@example Fitting
scatter(data, xaxis = :log, legend = :topleft, label = "data")
plot!(adaptive_result.result, label = "fit result", color = 3)
```

But comparing the estimated densities from the simple model fit and the adaptive model reveals the difference:

```@example Fitting
plot(DensityPlot(result.grid), xaxis = :log, color = 2, 
	fill = 0, fillalpha = 0.5, legend = :topleft, label = "simple fit")

plot!(DensityPlot(adaptive_result.grid), color = 3, 
	fill = 0, fillalpha = 0.5, label = "adaptive fit")
```

## The same options in the convenience workflow

As mentioned above, the convenience workflow exposes the same options as [`adaptive_dose_response_fit`](@ref). Although the [adaptive model fitting](@ref adaptive_model_fitting) section uses options that are already the defaults of the convenience workflow, this section creates these options explicitly to illustrate how to modify the options and to illustrate the convenience gain.

First, the log-prior-generator (with the scale `500`) can be obtained with [`scaled_log_volume_prior`](@ref):

```@example FittingConvenience
using AntibodyMethodsDoseResponseConvenience #hide
adaptive_result, data, replicates =  load_results("examples/default_concentrations") # hide
concentrations = data.independent # hide
responses = data.dependent # hide
errors = data.errors # hide
log_prior_generator = scaled_log_volume_prior(500)
```

The [`FittingData`](https://translational-pain-research.github.io/FittingObjectiveFunctions-documentation/API/#FittingObjectiveFunctions.FittingData) object is created as before.

```@example FittingConvenience
data = FittingData(concentrations,responses, errors, distributions = (y,m,Δy)-> -(y-m)^2/Δy^2)
``` 

The `Optim.jl` minimizers can be obtained with less boilerplate code, using [`minimizer_generator`](@ref):

```@example FittingConvenience
minimizer = minimizer_generator(NelderMead())
post_minimizer = minimizer_generator(LBFGS())
``` 

In general, the grid is created automatically, based on the concentration range of the dose-response data. However, if one concentration is `0`, the automatic grid should not be used.

```@example FittingConvenience
grid = create_grid(LogRange(1e-10,1e-2,3)) 
``` 

!!! info "automatic grids and the zero concentration"
	The automatic grids are created with the [`LogRange`](@ref) function and are intended for logarithmic plots. A zero-concentration would lead to a `DomainError`, as [`LogRange`](@ref) demands positive numbers. To avoid an error, the automatic grid generator then substitutes `eps()` for `0` (only for the auto-generated grid) and raise a warning. However, too large grid domains lead to poor fit results. Hence, it is recommended to create the grid manually in those cases.

The same options as before are used:

```@example FittingConvenience
adaptive_options = AdaptiveOptions(objective = :log_posterior, iterations = 30, 
			offset = eps(), prior_generator = log_prior_generator)

post_options = AdaptiveOptions(objective = :log_posterior, offset = eps(), 
			prior_generator = log_prior_generator)
``` 

Finally, the objects, functions and options from above need to be summarized in a [`FittingCondition`](@ref) object:

```@example FittingConvenience
condition = FittingCondition(data,
	grid = grid,
	options_1 = adaptive_options,
	options_2 = post_options,
	minimizer_1 = minimizer,
	minimizer_2 = post_minimizer
	)
``` 

Fitting the `condition` is done by calling [`fit_condition`](@ref):

```julia
result = fit_condition(condition)
plot(DensityPlot(result.grid), xaxis = :log, fill = 0, fillalpha = 0.5)
```
```@example FittingConvenience
plot(DensityPlot(adaptive_result.grid), xaxis = :log, fill = 0, fillalpha = 0.5) #hide
```


## Multi-threaded fitting

The fitting of a single condition cannot be parallelized. However, fitting more than one condition can be done in parallel. Consider an array of [`FittingCondition`](@ref) objects:
```julia
conditions = [condition_1,condition_2,condition_3]
```
The conditions can be fitted in parallel with [`fit_conditions`](@ref) (observe the additional `s`):
```julia
fit_conditions(conditions)
```

!!! info "Multi-threading"
	The number of threads cannot be changed after Julia is launched. The number of threads needs to be set in before, e.g. with `julia -t 8` to run Julia with 8 threads. The number of threads can be checked with `Threads.nthreads()`.

## Changing the initial values


In each case discussed here, no initial parameter values were set for the model fitting. This is, because the wights of the gird are used as initial values. By default, [`create_grid`](https://translational-pain-research.github.io/AdaptiveDensityApproximation-documentation/api/#AdaptiveDensityApproximation.create_grid) sets all weights to `1`:

```@example FittingConvenience
grid = create_grid(LogRange(1e-10,1e-2,3)) 
export_weights(grid)
```

Changing the weights of the gird before running the fitting-process allows to specify the initial values for the model fitting. E.g. for the manual definition of adaptive fitting:

```julia
import_weights!(grid,[eps(),1])
adaptive_result = adaptive_dose_response_fit(grid,data,minimizer, options = fitting_options)
```

Or for the convenience workflow, e.g. after the `FittingCondition` has already been defined:

```julia
import_weights!(condition.grid,[eps(),1])
result = fit_condition(condition)
```