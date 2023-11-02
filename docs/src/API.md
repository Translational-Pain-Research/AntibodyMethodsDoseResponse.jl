# API


## Measurement data processing

```@docs
dose_response_check
```

```@docs
normalize_data
```

```@docs
normalize_data!
```


## Depletion correction of measurement data


```@docs
scale_bound
```

```@docs
simple_depletion_correction
```

```@docs
simple_depletion_correction!
```

## Models

```@docs
accumulation_model
```

```@docs
langmuir_model
```

```@docs
accumulation_inv_const_model 
```

```@docs
langmuir_inv_const_model
```


## Result type

```@docs
DoseResponseResult
```

## Adaptive fitting

```@docs
AdaptiveOptions
```

```@docs
area_scaled_variation
```

```@docs
log_area_scaled_variation
```

```@docs
adaptive_dose_response_fit
```

```@docs
AdaptiveResult
```


## Convenience methods

```@docs
LogRange
```

```@docs
peak_detection
```



## Uncertainty estimation


```@docs
EpitopeUncertainty
```

```@docs
DoseResponseUncertainty
```