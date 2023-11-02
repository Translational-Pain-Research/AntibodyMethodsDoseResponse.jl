# [API - Convenience package](@id api_convenience)


## Fitting data

```@docs
scaled_log_volume_prior
```

```@docs
minimizer_generator
```

```@docs
FittingCondition
```

```@docs
fit_condition
```

```@docs
fit_conditions
```




## Loading data

```@docs
load_results
```

## Modify default plotting options

!!! info
	In the plotting functions, the default options are not set up as tuples of keyword arguments, but instead as function calls. The following functions generate the default keyword arguments, but allow to selectively change individual keywords or pass new `Plots.jl` keywords. In this way, changing the keywords does not require to manually repeat the unchanged default keyword arguments.


```@docs
dr_base_plot
```

```@docs
density_base_plot
```

```@docs
data_options
```

```@docs
replicate_options
```

```@docs
fit_options
```

```@docs
density_options
```

```@docs
eu_options
```

```@docs
du_options
```


## Plotting functions

```@docs
bin_analysis_plot
```

```@docs
peak_analysis_plot
```


```@docs
uncertainty_plot
```