# Density plots for `OneDimGrid`

Consider the following density:

```@example ResultsAndSimulations
using AntibodyMethodsDoseResponseConvenience
using Distributions

p(x) = pdf(Normal(1e-5,3e-6),x)
plot(p, xaxis = :log, xlims = [1e-10,1e-2], fill = 0, label = "density")
```

In the previous tutorials it was explained how this density can be approximated with a `OneDimGrid` and how the `OneDimGrid` can be plotted:

```@example ResultsAndSimulations
using AdaptiveDensityApproximationRecipes
grid = create_grid(LogRange(1e-10,1e-2,50))
approximate_density!(grid, p, volume_normalization = true)
plot(grid, xaxis = :log)
```

While the grid-plot can be helpful, it can also be distracting for an assessment of the density approximation. 

## The `DensityPlot` recipe

The [`AntibodyMethodsDoseResponseRecipes.jl`](https://github.com/AntibodyPackages/AntibodyMethodsDoseResponseRecipes.jl) package offers a plotting recipe to plot the density that a `OneDimGrid` approximates (this package is automatically imported by [`AntibodyMethodsDoseResponseConvenience.jl`](@ref api_convenience)):

```@example ResultsAndSimulations
plot(DensityPlot(grid), xaxis = :log, fill = 0, volume_normalization = :linear, label = ":linear")
```

Observe that the ``y``-axis now has the same scale as in the density plot.

!!! info "volume_normalization"
	The `OneDimGrid` weights correspond to the number of epitopes ``\lambda_j`` with ``K_\tau`` in a given interval ``I_j``, not to the density value ``g_j = \frac{\lambda_j}{\text{length}(I_j)}``. Thus it is necessary to divide the weights by the interval lengths (`volume_normalization = :linear`) to obtain the true density shape. However, using the visual volumes of the logarithmic scale (`volume_normalization = :log`) represents the contribution of peaks to the dose-response curve more accurately (see [Background: log-volume normalization](@ref log_volume_normalization)). Hence, it is the default normalization despite not being the "true" density plot.


```@example ResultsAndSimulations
plot(DensityPlot(grid), xaxis = :log, volume_normalization = :none, 
	fill = 0, fillapha = 0.5, label = ":none = grid plot", legend = :topleft)

plot!(DensityPlot(grid), volume_normalization = :log, fill = 0, fillalpha = 0.5, label = ":log")
```

With `DensityAnnotations` it is possible to display the total number of epitopes within selected annotation bins, on top of an already existing density plot:

```@example ResultsAndSimulations
plot(DensityPlot(grid), label = "grid")
plot!(DensityAnnotations(grid), xaxis = :log, 
	volume_normalization = :log,
	annotation_bins = [[1e-10,1e-6], [1e-6,1e-4],[1e-4,1e-2]],
	annotation_size = 10,
	annotation_offset = [0.2,0.1,0.2],
	hover_points = false,
	annotation_color = :black
	)
```

### Keywords and default arguments

* `volume_normalization = :log`: This keyword is needed to estimate the height of the annotation markers. Hence it should match the choice for the density plot.
* `annotation_bins = []`: The selected bins as vector of vectors. The bins are selected based on the ``K_\tau`` value, not based on the grid interval index. 
* `annotation_size = 10`: Font size for the annotation labels (that show the number of epitopes within the bin).
* `annotation_offset = 0.05`: Offset of the annotation labels from the top of the plot (as fraction of the plot height). When a single `value` is provided, the offsets alternate between `0` and the provided `value`. Alternatively, the offsets can be provided independently for each annotation bin, by providing a length-matched array of offsets.
* `hover_points = false`: Only recommended if the [Plotly/PlotlyJs backend](https://docs.juliaplots.org/latest/backends/) is used. If `true`, scatter points with the annotation label as tooltip are added, as the [Plotly/PlotlyJs backend](https://docs.juliaplots.org/latest/backends/) handles annotations differently.
* `annotation_color = :black`: The color for the annotation labels and bin lines.


!!! info "Number of epitopes"
	The "number of epitopes" is calculated in the dimension/unit that the response values are provided in. If the response value is the total number of bound antibody-epitope complexes, the calculated "number of epitopes" is indeed the number of epitopes. If the responses are measured as quantity that is only proportional to the number of bound complexes, the "number of epitopes" is only proportional to the actual number of epitopes, with the same proportionality constant.

## [Background: log-volume normalization](@id log_volume_normalization)

!!! warning "log-volume and standard volume"
	**log-volume** means the visual volume (visual interval length) in a plot with logarithmic scale, **not the logarithm of the volume!**. It corresponds to `volume_normalization = :log`. **Standard volume** is the true volume (interval length), corresponding to `volume_normalization = :linear`.

To illustrate the benefit of the log-volume normalization for the analysis of dose-response data, we may consider 2 peaks:

```@example ResultsAndSimulations
p_A(x) =  pdf(Normal(1e-6, 3e-7),x)
p_B(x) = pdf(Normal(1e-5,3e-6),x)
p(x) = p_A(x) + p_B(x)

plot(p, xaxis = :log, xlims = [1e-10,1e-2], fill = 0, label = "density")
```

As before, the density is approximated with a `OneDimGrid`:


```@example ResultsAndSimulations
grid = create_grid(LogRange(1e-10,1e-2,100))
approximate_density!(grid, p, volume_normalization = true)
```

Next, we compare the preferred (default) log-volume normalization with the standard volume normalization (which is true to the actual density). To identify the individual effect of the peaks, we approximate `p_A(x)` and `p_B(x)` individually:

```@example ResultsAndSimulations

grid_A = create_grid(LogRange(1e-10,1e-2,100))
approximate_density!(grid_A, p_A, volume_normalization = true)

grid_B = create_grid(LogRange(1e-10,1e-2,100))
approximate_density!(grid_B, p_B, volume_normalization = true)

lv_plot = plot(DensityPlot(grid_A), xaxis = :log, fill = 0, title = "log-volume", label = "peak A")
plot!(DensityPlot(grid_B), fill = 0, label = "peak B")


v_plot = plot(DensityPlot(grid_A), xaxis = :log, volume_normalization = :linear, fill = 0, title = "volume", label = "peak A")
plot!(DensityPlot(grid_B), volume_normalization = :linear, fill = 0, label = "peak B")

plot(lv_plot, v_plot, layout = (1,2), size = (800,300))
```

With the log-volume normalization both peaks have the same size, suggesting that each peak has the same effect strength w.r.t. to the dose-response curve, albeit at different locations. With the standard volume normalization, however, peak A is significantly larger than peak B, indicating that peak B affects the dose-response curve only marginally. 

Since we have approximated grids for each peak individually, as well as for the whole density, it can easily be checked how strong the effect of each peak is, by simulating the corresponding dose-response curves:


```@example ResultsAndSimulations
concentrations = LogRange(1e-10,1e-2, 100)

dr_A = DoseResponseResult(grid_A, concentrations)
dr_B = DoseResponseResult(grid_B, concentrations)
dr_total = DoseResponseResult(grid, concentrations)

plot(dr_A, label = "peak A", xaxis = :log, legend = :topleft, linewidth = 3)
plot!(dr_B, label = "peak B", linewidth = 3)
plot!(dr_total, label = "both peaks", linewidth = 3)
```

As the log-volume normalized plot suggested, both peaks have the same effect strength w.r.t. to the dose-response curve. The only difference is the location of the effect.

!!! remark "Visualization"
	While the standard volume normalization shows the true density, the log-volume normalization is more useful for the analysis of dose-response curves.