# Background: LSQ

For data points ``\{(x_i,y_i,\Delta y_i)\}_{i=1^N}``, the (weighted) least squares objective function is

```math
\text{lsq}(\lambda) = \sum_{i=1}^N \frac{(y_i - m(x_i,\lambda))^2}{\Delta y_i^2}
```
For the standard least squares objective function, one sets ``Δy_i = 1`` for all ``i = 1,\ldots,n``.

The optimal parameters ``\lambda``, given the data ``\{(x_i,y_i,\Delta y_i)\}`` and the model function ``m(x,\lambda)``, are those that minimize the least squares objective function ``\text{lsq}(\lambda)``. An explanation for this statement can be found in [Background:-Posterior-probability](@ref)
