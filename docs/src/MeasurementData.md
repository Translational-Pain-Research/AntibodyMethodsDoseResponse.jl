# Measurement data

To work with measurement data, it needs to be stored in a [`FittingData`](https://antibodypackages.github.io/FittingObjectiveFunctions-documentation/API/#FittingObjectiveFunctions.FittingData) object from [`FittingObjectiveFunctions.jl`](https://antibodypackages.github.io/FittingObjectiveFunctions-documentation/).

```julia
using FittingObjectiveFunctions
fitting_data = FittingData(concentrations, mean_responses, uncertainties)
```

## Check data

Good dose-response data must satisfy the following properties:

* The concentrations and the responses must be positive.
* All values must be real numbers.
* No `NaN` or `Inf`.

These properties can be checked, using [`dose_response_check`](@ref):

```@example MeasurementData
using AntibodyMethodsDoseResponseConvenience #hide
nothing #hide
```
```@repl MeasurementData
dose_response_check(FittingData([1,2],[-1,2]))
```

If [`dose_response_check`](@ref) does not throw an error, the [`FittingData`](https://antibodypackages.github.io/FittingObjectiveFunctions-documentation/API/#FittingObjectiveFunctions.FittingData) object contains proper dose-response data.

!!! info
	[`dose_response_check`](@ref) is used internally by some functions to ensure proper data transformation. E.g. data normalization or depletion corrections require proper dose-response data.

## Normalize data

Dose-response curves from different sources may require normalization. For example, replicates from different experiments could have used different exposure times. There are two methods to normalize the data. The first method requires a reference signal (scanned during the experiments) to which the other responses can be normalized. The second method normalizes the dose-response curve s.t. the strongest response is `1`.

To normalize the data to `1` use the [`normalize_data`](@ref) function without additional keywords.

```@example MeasurementData
data = FittingData([1,2,3],[1,2,3])
normalize_data(data)
```

To normalize the data to a reference signal, i.e. `2`, the `reference` keyword can be used:

```@example MeasurementData
normalize_data(data, reference = 2)
```

Another common normalization step is the removal of a constant offset, e.g. caused by autofluorescence. Since both offset removal and reference scaling are applied, it is important to understand the order of the corrections:

* The offset is subtracted from all response values
    + If the `reference` keyword is used, the offset is also subtracted from the reference value. This is, because reference measurements often suffer from the same offset.

* The responses are rescaled	
    + If the `reference` keyword was used, all responses are divided by the modified reference value.
    + Otherwise, the responses are rescaled, s.t. the largest response is `1`.

In the following example, `0.5` is subtracted from the response values `[1,2,3]->[0.5,1.5,2.5]`. Then the responses are divided by the largest response `2.5` to rescale the responses:
```@example MeasurementData
normalize_data(data, offset = 0.5)
```

!!! info "Measurement errors"
	The measurement errors are also normalized, following standard Gaussian error propagation (`reference` and `offset` are without errors). Thus, the measurement errors are just rescaled.

Without the `reference` keyword, it is possible to have the offset value be determined automatically, i.e. the dose-response curve is shifted s.t. the smallest response is zero:

```julia
normalize_data(data, offset = Inf)
```

Internally, any offset that is too large is replaced by the smallest response value. When the `reference` keyword is used, it is not recommended to use arbitrary offsets, as they are subtracted from the reference value, too. Furthermore, offset values that are larger than the reference value *(after being replaced by the smallest response)* throw an error:

```julia
normalize_data(data, offset = 2, reference = 0.1)
```
```@example MeasurementData
println("ERROR: DomainError with 1:") #hide
println("The offset value must be smaller than the reference value.") #hide
```

!!! tip "Mutation of data"
	[`normalize_data`](@ref) does not mutate the original [`FittingData`](https://antibodypackages.github.io/FittingObjectiveFunctions-documentation/API/#FittingObjectiveFunctions.FittingData) object, but returns a normalized copy. To mutate the original object, use [`normalize_data!`](@ref).


## Simple depletion correction

To keep the models computationally simple, e.g. by using analytical solutions, antibody depletion is not accounted for. Yet, depending on the experimental protocol, depletion can be unavoidable. Fortunately, there is simple approximation to obtain new concentrations that define lower bounds for the real binding process that is subject to depletion:

Let ``\{(a_i,r_i)\}`` be data points where, ``a_i`` denotes the initial antibody concentrations and ``r_i`` denotes the corresponding responses. Then lower-bound concentrations are 

```math
b_i = a_i - \widehat{\beta} r_i \qquad \text{where} \qquad \widehat{\beta} = \max \{\beta \in \mathbb{R}_0 \mid a_i -\beta r_i \geq 0 \ \forall i\}
```

This approximation holds true both for the accumulation model as well as the Langmuir isotherm model, albeit for different reasons (cf [Edwards et al](https://doi.org/10.1006/abio.1998.2814) and [https://arxiv.org/abs/2407.06052](https://arxiv.org/abs/2407.06052)). The maximization of the parameter ``\widehat{\beta}`` is necessary to obtain a worst-case scenario if the actual proportionality factor between the response signal and the actual number of bound complexes is unknown.


To obtain this "scale bound" ``\widehat{\beta}`` the [`scale_bound`](@ref) function can be used:

```@example MeasurementData
β = scale_bound(FittingData([1,2,3],[1,1.5,1.9]))
```

To obtain the lower-bound concentrations ``b_i``, the [`simple_depletion_correction`](@ref) function can be used:

```@example MeasurementData
simple_depletion_correction(FittingData([1,2,3],[1,1.5,1.9]), β)
```

Note that a new [`FittingData`](https://antibodypackages.github.io/FittingObjectiveFunctions-documentation/API/#FittingObjectiveFunctions.FittingData) object is returned, containing the lower-bound concentrations `b_i` instead of the initial concentrations `a_i`. It is possible to use different scales for [`simple_depletion_correction`](@ref). But for the [`scale_bound`](@ref) there is a shortcut (which internally calls [`scale_bound`](@ref)):

```@example MeasurementData
simple_depletion_correction(FittingData([1,2,3],[1,1.5,1.9]))
```

!!! tip "Mutation of data"
	Again, [`simple_depletion_correction`](@ref) does not mutate the original [`FittingData`](https://antibodypackages.github.io/FittingObjectiveFunctions-documentation/API/#FittingObjectiveFunctions.FittingData) object, but returns a corrected copy. To mutate the original object, use [`simple_depletion_correction!`](@ref).

## [Plotting](@id measurement_data_plotting)

Of course, [`FittingData`](https://antibodypackages.github.io/FittingObjectiveFunctions-documentation/API/#FittingObjectiveFunctions.FittingData) objects can be plotted with any plotting library by calling the fields:

```@example MeasurementData
using Plots
data = FittingData([0,1,2,3,4],[0,1,1.5,1.9,2.1])
plot(data.independent, data.dependent, yerrors = data.errors)
```

For convenience, [`AntibodyMethodsDoseResponseRecipes.jl`](https://github.com/AntibodyPackages/AntibodyMethodsDoseResponseRecipes.jl) contains plotting recipes for [`FittingData`](https://antibodypackages.github.io/FittingObjectiveFunctions-documentation/API/#FittingObjectiveFunctions.FittingData) objects:

```@example MeasurementData
using Plots, AntibodyMethodsDoseResponseRecipes
data = FittingData([0,1,2,3,4],[0,1,1.5,1.9,2.1])
plot(data)
```

Observe that the data point with `concentration = 0` is missing. Since dose-response curves are commonly plotted in a logarithmic scale and since [`Plots.jl`](https://docs.juliaplots.org/stable/) does not automatically remove zero values from logarithmic plots (which become ``-\infty``), the plotting recipe introduces the keyword `filter_zeros = [true,false]`. The `filter_zeros` keyword expects two `Bool` values. If the first value is `true`, all data points with `concentration = 0` are removed from the plot. Accordingly, if the second value is `true`, all data points with `response = 0` are removed from the plot.


The plotting recipe defines how [`Plots.jl`](https://docs.juliaplots.org/stable/) handles the data inside a [`FittingData`](https://antibodypackages.github.io/FittingObjectiveFunctions-documentation/API/#FittingObjectiveFunctions.FittingData) object. Thus, the usual keywords remain usable and can be combined with the `filter_zeros` keyword:
```@example MeasurementData
plot(data, color = :red, filter_zeros = [false,false])
```

## [`FittingData`](https://antibodypackages.github.io/FittingObjectiveFunctions-documentation/API/#FittingObjectiveFunctions.FittingData) and [`FittingCondition`](@ref)

Internally [`FittingCondition`](@ref) objects expect measurement data to be stored in [`FittingData`](https://antibodypackages.github.io/FittingObjectiveFunctions-documentation/API/#FittingObjectiveFunctions.FittingData) objects. Yet, the [quick start guide](@ref quick_start) mentions [`FittingData`](https://antibodypackages.github.io/FittingObjectiveFunctions-documentation/API/#FittingObjectiveFunctions.FittingData) objects only as optional method to implement specific measurement errors. This is, because [`FittingCondition`](@ref) offers constructors that automatically create the [`FittingData`](https://antibodypackages.github.io/FittingObjectiveFunctions-documentation/API/#FittingObjectiveFunctions.FittingData) objects.

The default constructor uses [`FittingData`](https://antibodypackages.github.io/FittingObjectiveFunctions-documentation/API/#FittingObjectiveFunctions.FittingData) objects:

```@example MeasurementData
data = FittingData([1,2,3],[1,1.5,1.9], [0.4,0.4,0.4], distributions = (y,m,Δy)-> -abs(y-m))
condition = FittingCondition(data)
condition.data 
```

Using the [`FittingData`](https://antibodypackages.github.io/FittingObjectiveFunctions-documentation/API/#FittingObjectiveFunctions.FittingData) constructor allows to specify different measurement errors, as mentioned in the tip, but also different uncertainty distributions for posterior based objectives.

The convenience constructor bypasses the need to construct the [`FittingData`](https://antibodypackages.github.io/FittingObjectiveFunctions-documentation/API/#FittingObjectiveFunctions.FittingData) explicitly, expecting only the concentrations and the responses:

```@example MeasurementData
condition = FittingCondition([1,2,3],[1,1.5,1.9])
condition.data 
```
Note that the measurement errors are set to `1`.

Finally, the convenience constructor allows to pass multiple response arrays:

```@example MeasurementData
condition = FittingCondition([1,2,3],[1,1.5,1.9],[0.8,1.7,2.8])
condition.data 
```

If multiple response arrays are provided to the constructor, it calculates the mean values for the responses and uses the standard deviation for the measurement errors. In addition, it adds the individual responses as `replicates`:

```@example MeasurementData
condition.replicates 
```