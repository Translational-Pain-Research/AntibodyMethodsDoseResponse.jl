# API

## Create new grids


```@docs
create_grid
```

## Approximate density functions

```@docs
approximate_density!
```


## Import and Export


```@docs
export_weights
```

```@docs
export_all
```

```@docs
import_weights!
```




## Refine a grid


```@docs 
refine!
```


```@docs
subdivide!
```


## Restrict the grid domain

```@docs
restrict_domain!
```

```@docs
select_indices
```

## Simple calculations

```@docs
sum(::Union{AdaptiveDensityApproximation.Grid, AdaptiveDensityApproximation.OneDimGrid})
```

```@docs
sum(::Function,::Union{AdaptiveDensityApproximation.Grid, AdaptiveDensityApproximation.OneDimGrid})
```

```@docs
prod(::Union{AdaptiveDensityApproximation.Grid, AdaptiveDensityApproximation.OneDimGrid})
```

```@docs
prod(::Function,::Union{AdaptiveDensityApproximation.Grid, AdaptiveDensityApproximation.OneDimGrid})
```

```@docs
integrate
```

## Integral models



```@docs
integral_model
```

## Numeric PDF and CDF

```@docs
get_pdf
```

```@docs
get_cdf
```

## Grid slices

```@docs
get_slice
```

```@docs
dimension
```
