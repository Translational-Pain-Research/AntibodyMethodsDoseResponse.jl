# Results and Simulations

## Why `OneDimGrid` objects?

The model generators requires a `OneDimGrid` objects from the [AdaptiveDensityApproximation](@ref AdaptiveDensityApproximation) package (see [Models](@ref binding_models)). For this, the actual weights of the grid do not matter. The model generator returns a new parameter array with the goal to estimate these parameters. At this point it is also possible to choose parameter values for this array and to simulate dose-response curves with the corresponding model function.

But the [AdaptiveDensityApproximation](@ref AdaptiveDensityApproximation) offers additional methods, e.g. the name-giving adaptive approximation of densities or additional density conversion / analysis tools. For this reason, internal methods use `OneDimGrid` objects and construct the model functions from scratch for the calculation of dose-response curves. This is, among others, the reason why most functions (e.g. the convenience functions of [`AntibodyMethodsDoseResponseConvenience`](@ref api_convenience)) return or expect those grids. 


!!! tip
	Although not necessary for all analyses, it is highly recommended to use the `OneDimGrid` workflow, as advanced methods of this package will expect grids.


## Simulating the ``K_\tau``-density

To simulate a dose-response curve with a binding model, a density ``g(K_\tau)`` needs to be defined. This can be done with any Julia function, e.g. a probability density function:

```@example ResultsAndSimulations
using AntibodyMethodsDoseResponseConvenience
using Distributions

p(x) = pdf(Normal(1e-5,3e-6),x)
plot(p, xaxis = :log, xlims = [1e-10,1e-2], fill = 0, label = "density")
```
The `plot` function can be used without importing `Plots.jl`, as [`AntibodyMethodsDoseResponseConvenience`](@ref api_convenience) does this in the background automatically.

Next, a `OndDimGrid` needs to be defined to approximate the density with:

```@example ResultsAndSimulations
grid = create_grid(LogRange(1e-10,1e-2,50))
approximate_density!(grid, p, volume_normalization = true)
```

!!! info "volume_normalization = true"
	Without `volume_normalization = true` (the default), the density function itself is approximated (e.g. by using the function values at the centers of the respective intervals). The volume normalization approximates the area under the curve by using `function value × interval length`. This is a convenient choice for simulations, as the area under a density curve corresponds to the number of epitopes in the given interval, which is exactly what the weights of the grid should be.

To inspect the grid, the plotting recipe from the `AdaptiveDensityApproximationRecipes` package can be used (as already described in [Models](@ref binding_models)):


```@example ResultsAndSimulations
using AdaptiveDensityApproximationRecipes
plot(grid, xaxis = :log)
```

!!! info "y-axis scale"
	Since the grid approximates the area under the density curve (weights correspond to the height of the bars), the scale of the ``y``-axis is different, compared to the density plot.



## Simulating the dose-response curve

With the `OneDimeGrid` object `grid`, a dose-response can be simulated, using [`DoseResponseResult`](@ref):

```@example ResultsAndSimulations
concentrations = LogRange(1e-8,1e-2,16)
simulation_result = DoseResponseResult(grid,concentrations; model = accumulation_model, offset = 0)
scatter(simulation_result, xaxis = :log)
```

!!! info "DoseResponseResult"
	`DoseResponseResult` objects are used to construct and store dose-response curves, resulting from grids. This can be either because a grid was used to simulate a ``K_\tau``-density or because the result of a model-fit was stored in a grid. The package `AntibodyMethodsDoseResponseRecipes` provides a plotting recipe for `DoseResponseResult` objects.

!!! tip "filter_zeros"
	The same `filter_zeros` keyword that is used for [`FittingData`](@ref) plots can also be used for [`DoseResponseResult`](@ref) plots (see [Measurement Data plotting](@ref measurement_data_plotting)).

## Importing fitting results

Instead of simulating a dose-response curve, the `OneDimGrid` and [`DoseResponseResult`](@ref) types can also be used to create dose-response curves from a fitting result. Since no data has been fitted yet, we construct a result from the weights of the simulation grid above:

```@example ResultsAndSimulations
parameters = export_weights(grid)
push!(parameters, 0.2) # Add offset parameter
```

We assume now, that the `gird` was used to obtain the model function and that the `parameters` are the result of the model fit. In this case, the `offset` keyword was used. To import the result into the grid, the [`import_weights!`](@ref) functions can be used:

```@example ResultsAndSimulations
import_weights!(grid, parameters[1:end-1])
```

!!! info "Why end-1 ?"
	Since the `offset` keyword was used for the creation of model functions (in our scenario), the `parameters` array has one additional element at the end, the offset, which does not belong to the grid. If no `offset` is used, all parameters should be imported:
	```julia
	import_weights!(grid, parameters)
	```
	In fact, using the wrong selection of parameters usually leads to a `DimensionMismatch` error, as the number of parameters and the number of intervals in the grid do not match.


```@repl ResultsAndSimulations
import_weights!(grid, parameters) # Wrong parameter selection should fail.
```

After importing the `parameters` into the grid, the dose-response curve can be constructed as before (for any concentrations, not only the data concentrations):


```@example ResultsAndSimulations
concentrations = LogRange(1e-8,1e-2,16)
result = DoseResponseResult(grid,concentrations; model = accumulation_model, offset = parameters[end])
scatter(result, xaxis = :log, ylims = [0,1.4])
```
Observe the `offset = parameters[end]` part at the end of the [`DoseResponseResult`](@ref) call. In our scenario, the `parameters` array contains the offset as last element.


