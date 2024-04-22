# [Quick start](@id quick_start)

The most convenient, and recommended, way to analyze dose-response data is to use the [`AntibodyMethodsDoseResponseConvenience`](@ref api_convenience) package. This tutorial will not cover all details (e.g. the internal procedure), but present a short introduction to get started. To see all details and available options, please have a look at the [AnitbodyMethodsDoseResponseConvenience API](@ref api_convenience).

## Starting point

Suppose, a dose-response experiment with 3 replicates was conducted, leading to the following dose-response curves:

```@example QuickStart
using AntibodyMethodsDoseResponseConvenience, Measures # hide
results, data, replicates =  load_results("examples/default_concentrations") # hide

scatter(data.independent, replicates[1].dependent, xaxis = :log, label = "replicate 1", legend = :topleft, xlabel = "concentration", ylabel = "response") # hide
scatter!(data.independent, replicates[2].dependent, xaxis = :log, label = "replicate 2", markershape = :rect) # hide
scatter!(data.independent, replicates[3].dependent, xaxis = :log, label = "replicate 3", markershape = :utriangle) # hide
```

## Loading Data

Importing general data into Julia is not the scope of the `AntibodyMethods` packages. For this, use e.g. `DelimitedFiles.jl` or `CSV.jl` in conjunction with `DataFrames.jl`. Nevertheless, as short introduction, assume that the data is stored in a csv file, where the columns are (concentrations, replicate 1, ..., replicate 3). Using `DelimitedFiles.jl` the data can be imported as follows:


```julia
using DelimitedFiles
csv_data = readdlm("path_to_file")
```
```@example QuickStart
csv_data = hcat(data.independent, replicates[1].dependent,replicates[2].dependent,replicates[3].dependent) # hide
```

To separate the different columns into independent arrays, as is needed in the following, use:

```@example QuickStart
conc  = csv_data[:,1] # concentrations
rep_1 = csv_data[:,2] # replicate 1
rep_2 = csv_data[:,3] # replicate 2
rep_3 = csv_data[:,4] # replicate 3
nothing #hide
```

## Defining a fitting condition

Before the data can be analyzed (model fitting), it needs to be summarized in a [`FittingCondition`](@ref) object. A [`FittingCondition`](@ref) object contains both the data and the fitting instructions, which can be specified with keywords. The (recommended) default analysis can be obtained by only passing the data, and the `scale` keyword.

```@example QuickStart
fitting_condition = FittingCondition(conc, rep_1,rep_2,rep_3, scale = 500, path = "path_to_store_results")
nothing # hide
```

!!! info "The scale keyword"
	To avoid overfitting of sparse, noisy dose-response data, some sort of regularization constraint is necessary. The `scale` keyword determines the scale of a penalty term to prevent jagged/spiky histograms. Not using the `scale` keyword results in a weighted least squares objective.
	
!!! warning "Inappropriate scale value"	
	Too small `scale` values only lead to an inefficient least squares objective. Too large scale values, however, can prevent the model to replicate the data at all. In cases, where the fitting result does not correspond to the data at all, try significantly smaller (orders of magnitude) `scale` values.

!!! info "The path keyword"
	Since curve-fitting can be a time-consuming process, it is recommended to save the results into files. This allows to re-plot the results at a later time, without having to re-run the fitting process. The `path` keyword defines the directory for the result files. If `path=""`, which is the default option, no files are saved.

!!! tip "Different measurement errors"
	Constructing a [`FittingCondition`](@ref) object by passing the different replicate responses will default to the standard deviation of the data points for the measurement error. If only a single response is used, the replicates field will be empty (`nothing`) and the measurement errors are set to `±1`. Different errors can be used by constructing the [`FittingData`](@ref) object manually:
	```julia
	errors = 0.1 .* responses
	data = FittingData(concentrations,responses, errors)
	fitting_condition = FittingCondition(data, scale = 500, path = "path_to_store_results")
	```

## Fitting a condition

After defining a [`FittingCondition`](@ref), the data can be analyzed, i.e. fitting the antibody-binding model to the data to obtain the estimated ``K_\tau``-density. Since the fitting instructions are already contained in the [`FittingData`](@ref) object, fitting the data is just a simple function call:

```julia
results, data, replicates = fit_condition(fitting_condition)
```

[`fit_condition`](@ref) returns an [`AdaptiveResult`](@ref) object `results`, the `data` as [`FittingData`](@ref) object and the `replicates` as array of [`FittingData`](@ref) objects for the respective replicates. If a `path` is specified for the [`FittingCondition`](@ref), the returned objects are also saved into files.

## Loading results*

In general, it is not necessary to load the results, as [`fit_condition`](@ref) returns all results. However, since the fitting process is time-consuming, it is a good idea to save the results into files (see `path` keyword above). After e.g. a restart of Julia, the results can be loaded as follows:


```julia
results, data, replicates = load_results("path_to_stored_results")
```

This allows to re-plot the results at a later time, without the need to re-run the analysis.

!!! warning
	The `results`, `data` and `replicates` are saved with `Serialization.jl`. Loading objects into a new instance of Julia does not re-instantiate referenced functions. In other words, the distribution functions of loaded `FittingData` objects do not work. Hence, loaded `FittingData` objects cannot be used for fitting. Instead, one should define new `FittingData` functions from scratch, e.g. by

		new_data = FittingData(data.independent, data.dependent, data.errors, distributions = ....)

	The default distributions created by the `FittingCondition` constructor are

		(y,m,Δy)-> -(y-m)^2/Δy^2

## Plotting the results

To visualize the results, different plotting functions are provided. The most straightforward one is [`bin_analysis_plot`](@ref), which returns 2 plots, the dose-response plot and the density plot:

```@example QuickStart
using Measures
dr_plot, density_plot = bin_analysis_plot(results,data,replicates)
plot(dr_plot, density_plot, layout = (1,2), size = (800,300), margins = 4mm)
```

The `Measures` package is used here, to define the margins between the individual plots.

!!! info "Plot objects"
	[`AntibodyMethodsDoseResponseConvenience`](@ref api_convenience) uses `Plots.jl` in the background to generate the plots. All methods and options of `Plots.jl` are immediately available; it is not necessary to import the `Plots` with  `using Plots`. The returned plots `dr_plot` and `density_plot` are full plot objects, that can be saved directly with `savefig(dr_plot,"file_name")`. `Plots.jl` also allows to re-draw plot objects, allowing to compose images (see [Partial plotting and the combination of plots](@ref partial_plotting)).


The second plotting function is [`peak_analysis_plot`](@ref), which visualizes the effect of different peaks on the dose-response curve:

```@example QuickStart
individual_dr_plot, cumulative_dr_plot, peak_plot = peak_analysis_plot(results, data)
plot(individual_dr_plot, cumulative_dr_plot, peak_plot, layout = (2,2), size = (800,600))
```

The top-left plot shows the individual dose-response curves of the different peaks (lower-left plot) as they would appear if only the respective peak was present alone. The cumulative plot (top-right) shows the contribution of each peak to the total dose-response curve.


!!! info "Selection of peaks"
	Without additional keywords, the peaks are determined automatically using [`peak_detection`](@ref). Peak regions can also be manually selected. E.g. two isolated ranges ``[10^{-7},10^{-5}]`` and ``[10^{-4}, 10^{-3}]`` could be selected using `bins = [[1e-7,1e-5],[1e-4,1e-3]]`.


## [Partial plotting and the combination of plots](@id partial_plotting)

Note that the peak analysis plots did not contain the replicates. This is by design, as the different plots can be combined. For this, we create a partial [`bin_analysis_plot`](@ref), where only the data points and the replicates are plotted:

```@example QuickStart
dr_data_plot, empty_density_plot = bin_analysis_plot(nothing,data, replicates)
plot(dr_data_plot)
```
Setting any of the arguments to `nothing` will omit the argument in the plot. E.g. to just plot the resulting curve, without the data `bin_analysis_plot(results,nothing,nothing)` would be used.

The cumulative dose-response, i.e. the contribution of the peaks to the total dose-response curve, can now be plotted on top of the data plot `dr_data_plot`, using the additional keyword `cumulative_dr_plot = dr_data_plot`. Since the data points are already contained in the `dr_data_plot`, they should not be plotted again with the [`peak_analysis_plot`](@ref), hence the `nothing`:

```@example QuickStart
individual_dr_plot, cumulative_dr_plot, peak_plot = peak_analysis_plot(results, nothing, cumulative_dr_plot = dr_data_plot)
plot(cumulative_dr_plot, peak_plot, layout = (1,2), size = (800,300), margins = 4mm)
```

## [Plotting options](@id plotting_options)

The plotting functions [`bin_analysis_plot`](@ref), [`peak_analysis_plot`](@ref) and [`uncertainty_plot`](@ref) use keyword-generating functions to modify the plotting options. The keyword-generating functions generate keywords and the default keyword arguments for the plot, but allow to overwrite the arguments for the individual keywords. For example the [`data_options`](@ref) function generate the following keywords:

```@example QuickStart
data_options()
```
The individual keywords can by changed when calling the function. Furthermore, it is possible to add new keywords:

```@example QuickStart
data_options(label = "new label text", markershape = :utriangle)
```

There is no need to save the resulting keyword arguments. The keyword generating function can be called (as keyword argument) in the plotting function:

```@example QuickStart
dr_plot, empty_plot = bin_analysis_plot(results,data, replicates,
	data_arguments = data_options(label = "new label text", markershape = :utriangle)
	)
plot(dr_plot)
```

!!! info "Why keyword generating functions?"
	The first idea to allow plot modifications involves global keywords, that could then be changed individually. However, different objects of the plots (e.g. the result and the data) have coinciding keywords, e.g. the color. But writing new, unique keywords for all the existing `Plots.jl` keywords is impossible.

	Hence tuples of keywords that correspond to the different parts of the plot (e.g. the result, the data and the replicates) are used. This allows to use the same keyword but with different options. E.g. `data_arguments` expects a tuple of keywords, which then only affect the plotting of the data points.

	To define the default behavior, a default tuple of keywords could be used. However, then if only one keyword needed to be changed, the user would have to create a new tuple, copying over the unchanged keywords, which is inconvenient. Thus, keyword-(tuple) generating functions are provided that produce the default keyword tuple but allow to change individual keywords. In fact `data_arguments = data_options()` is set as default internally.


Another example where generating functions are used to modify the plot, are the base-plot functions. In the previous section ([Partial plotting and the combination of plots](@ref partial_plotting)) it was used that [`bin_analysis_plot`](@ref), [`peak_analysis_plot`](@ref) and [`uncertainty_plot`](@ref) plot on top of preexisting plots. In fact, even when e.g. `dr_plot` is not explicitly changed, `dr_plot = dr_base_plot()` is used to create an empty plot (with some options pre-configured). Calling [`dr_base_plot`](@ref) allows to change those pre-configured options:

```@example QuickStart
dr_plot, empty_plot = bin_analysis_plot(results,data, replicates,
		dr_plot = dr_base_plot(xlabel = "concentrations", legend = :none)
		)
plot(dr_plot)
```

A complete list of these generating functions and how they are used in the different plotting functions can be found in the [AnitbodyMethodsDoseResponseConvenience API](@ref api_convenience).