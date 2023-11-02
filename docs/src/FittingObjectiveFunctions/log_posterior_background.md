# Background: Logarithmic posterior probability

The general setting is the same as in [Background:-Posterior-probability](@ref). The starting point is

```math
p(\lambda \mid \{x_i\}_{i=1}^N, \{y_i\}_{i=1}^N, m) \propto  p_0(\lambda\mid \{x_i\}_{i=1}^N, m) \prod_{i_1}^n \ell_i(y_i\mid \lambda, x_i,m)\ ,
```
assuming the likelihoods ``\ell_i`` are known.

## Product of small numbers

In the formula for the posterior likelihood, it can happen that many small numbers (close to zero) need to be multiplied together. Because [floating point numbers](https://en.wikipedia.org/wiki/Floating-point_arithmetic#Representable_numbers,_conversion_and_rounding) can only represent numbers up to a certain precision, such products, though theoretically non-zero, tend to be rounded to zero.

For example, consider the following array as the likelihood values:
``` @example 1
using Distributions, BenchmarkTools
small_values = [pdf(Normal(0,1),10+i) for i in 1:10]
```
Although the values are non-zero, the product is rounded to zero:
``` @example 1
prod(small_values)
```
One could use floating point types with higher precision:
``` @example 1 
small_values_high_precision = BigFloat.(small_values)
prod(small_values_high_precision)
```
However, this entails a huge performance loss together with increased memory usage: 

``` @example 1
@benchmark prod(small_values)
```

``` @example 1
@benchmark prod(small_values_high_precision)
```

## Logarithmic scale

Since posterior probabilities are often unnormalized anyways, one is not interested in the particular values, but only in relative differences. But then, any strictly monotonic function can be applied to compare relative differences. A convenient choice for such a strictly monotonic (increasing) function is the logarithm.

Note that because of proportionality, there is a constant ``\alpha > 0`` such that
```math
p(\lambda \mid \{x_i\}_{i=1}^N, \{y_i\}_{i=1}^N, m) =  \alpha \cdot  p_0(\lambda\mid \{x_i\}_{i=1}^N, m) \prod_{i_1}^n \ell_i(y_i\mid \lambda, x_i,m)\ .
```
Applying the natural logarithm leads to:

```math
\begin{aligned}
\ln (p(\lambda \mid \{x_i\}_{i=1}^N, \{y_i\}_{i=1}^N, m)) &=  \ln \left(\alpha \cdot  p_0(\lambda\mid \{x_i\}_{i=1}^N, m) \prod_{i_1}^n \ell_i(y_i\mid \lambda, x_i,m) \right) \\ 
&= \ln(p_0(\lambda\mid \{x_i\}_{i=1}^N, m)) + \sum_{i=1}^N \ln(\ell_i(y_i\mid \lambda, x_i,m)) + \ln(\alpha)
\end{aligned}
```
Using the logarithm allowed to exchange the multiplication of small numbers for an addition in the logarithmic scale, at the cost of having to calculate the logarithm of every value. However, the cost of calculating the logarithm is the worst case scenario. In many cases, it is possible if not easier to implement logarithms of the involved densities (e.g. for the normal distribution, laplace distribution, etc.).


To shorten the notation, denote the logarithms of the distributions by ``L_p = \ln \circ\ p`` for the posterior, ``L_i= \ln \circ\ \ell_i`` for the likelihoods and ``L_0 =\ln \circ\ p_0`` for the prior: 
```math
L_p(\lambda \mid \{x_i\}_{i=1}^N, \{y_i\}_{i=1}^N, m) =   L_0(\lambda\mid \{x_i\}_{i=1}^N, m) +  \sum_{i_1}^n L_i(y_i\mid \lambda, x_i,m) + \text{const.}
```

## Effect of Logarithmic scale

Using the logarithm has two effects. First of all, the product becomes a sum. This alone would suffice to prevent the rounding to zero problem:

``` @example 1
sum(small_values)
```
In addition, the logarithm has the effect of compressing the number scale for numbers larger than 1 and to stretch out the number scale for numbers between 0 and 1:
``` @example 1
log_values = log.(small_values)
```
Of course, the sum is still non-zero:
``` @example 1
sum(log_values)
```
While the use of higher precision floating point numbers (`BigFloat`) meant a huge performance loss, the log scale method dose not impair performance:
``` @example 1
@benchmark sum(log_values)
```

## Logarithmic density example

In [Logarithmic-scale](@ref) it was mentioned that some distributions are even easier to be implemented in a logarithmic scale. This is not only the case for the definition of densities from scratch, but also applies for `Distributions.jl`. Observe that "far" away from the mean, the pdf of a normal distribution is rounded to zero:
``` @example 1
pdf(Normal(0,1),100)
```
Obviously, the logarithm cannot be applied to this. However, `Distributions.jl` offers a `logpdf` function:
``` @example 1
logpdf(Normal(0,1),100)
```
This allows values even further away from the mean:
``` @example 1
logpdf(Normal(0,1),10^20)
```