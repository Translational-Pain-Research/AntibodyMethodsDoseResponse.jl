# Background: Posterior probability


The posterior objective (and the [Log posterior objective](#Log-posterior-objective) which is numerically favorable) allows to define more general objective functions. Form a Bayesian perspective, one is interested in the probability density for a particular parameter ``\lambda`` given the data ``\{x_i\}_{i=1}^N, \{y_i\}_{i=1}^N`` and the model ``m(x,\lambda)``:
``` math
p(\lambda \mid \{x_i\}_{i=1}^N, \{y_i\}_{i=1}^N, m)
```

## Applying Bayes' theorem

Using Bayes' theorem, the probability density can be rewritten as:
``` math
p(\lambda \mid \{x_i\}_{i=1}^N, \{y_i\}_{i=1}^N, m) = \frac{\ell(\{y_i\}_{i=1}^N \mid \{x_i\}_{i=1}^N , \lambda, m )\cdot p_0(\lambda\mid \{x_i\}_{i=1}^N, m)}{p(\{y_i\}_{i=1}^N \mid \{x_i\}_{i=1}^N , m)}
```
The denominator is but a normalization constant, that does not depend on ``\lambda``, i.e. can be ignored for optimization problems (and MCMC sampling):
``` math
p(\lambda \mid \{x_i\}_{i=1}^N, \{y_i\}_{i=1}^N, m) \propto \ell(\{y_i\}_{i=1}^N \mid \{x_i\}_{i=1}^N , \lambda, m )\cdot p_0(\lambda\mid \{x_i\}_{i=1}^N, m)
```
Because of the proportionality, one may refer to the right hand side as **unnormalized posterior**.

* ``\ell`` is a proper probability distribution for ``\{y_i\}_{i=1}^N`` given ``\{x_i\}_{i=1}^N, \lambda, m``. However, it can also be regarded as function of ``\lambda``, for fixed ``\{x_i\}_{i=1}^N``, ``\{y_i\}_{i=1}^N`` and ``m`` (which is needed, since the data is fixed, but different parameters need to be tested for model fitting). In this case, one calls it the **likelihood** function of ``\lambda``. It is no longer a proper probability density (still positive but no longer normalized).

* ``p_0`` is the so called **prior** distribution. It determines the probability of the parameters, before the data was obtained. This is sometimes called *belief in parameters* or *initial knowledge*.

!!! default "The prior and objectivity"
	A common critique is, that the prior is not objective. While the choice of prior can be subjective, it must be explicitly stated making all assumptions transparent. This allows for an objective comparison of the different approaches.

	In fact, there are two common types of priors in least squares fitting.

	1. ``p_0(\lambda\mid \{x_i\}_{i=1}^N, m) = 1``, i.e. a uniform prior. Since one usually uses a computer, there is a largest number ``b <\infty`` and a smallest number ``a > -\infty`` that the computer can use. Then one may choose the uniform distribution ``p_0(\lambda \mid \{x_i\}_{i=1}^N, m) = \frac{1}{b-a}``. Sine the posterior probability is only considered up to proportionality, one can simply use ``p_0(\lambda\mid \{x_i\}_{i=1}^N, m) = 1``. This leads to a **maximum likelihood** objective.

	2. In ill-defined problems, it is common practice to use some kind of regularization. In some cases, these regularizations correspond to certain priors. For example, the Tikhonov regularization essentially uses the prior ``p_0(\lambda\mid \{x_i\}_{i=1}^N, m) \propto \exp(-||\Gamma \lambda ||^2)``.

## Independent data points
A common assumption is that the data points are independent. While this is not a necessity, writing general likelihood functions is usually not trivial. If the data points are independent, the likelihood function becomes a product of likelihood functions for the individual data point likelihoods:
```math
\ell(\{y_i\}_{i=1}^N \mid \lambda, \{x_i\}_{i=1}^N m ) = \prod_{i=1}^N \ell_i(y_i\mid \lambda, x_i, m)
```
Note that the likelihoods can differ for the different data points, denoted by ``\ell_i`` here.  Thus the posterior probability / the objective function becomes
```math
p(\lambda \mid \{x_i\}_{i=1}^N, \{y_i\}_{i=1}^N, m) \propto  p_0(\lambda\mid \{x_i\}_{i=1}^N, m) \prod_{i_1}^n \ell_i(y_i\mid \lambda, x_i,m)
```

In general, ``x_i`` is only the measured value, while the true value ``\mathcal{X}_i`` is unknown. If the distribution ``p_i(\mathcal{X}_i\mid \lambda, x_i, m)`` is known, marginalization can be used to express the likelihood
```math
\ell(y_i\mid \lambda, x_i, m) = \int \ell(y_i \mid \mathcal{X}_i, \lambda, x_i, m)\cdot p(\mathcal{X}_i\mid \lambda, x_i,m) \ d\mathcal{X}_i
```
It can happen, that this integral can only be solved numerically. Since this is computationally expensive, and needs to be redone for every new value of ``\lambda``, the resulting posterior distribution is often not suited for optimization/sampling purposes. Another approach could be data-augmentation, e.g. to sample ``\ell(y_i, \mathcal{X}_i\mid \lambda, x_i, m)``, which is not the scope of this package.

The likelihood ``\ell(y_i\mid \mathcal{X}_i, \lambda, x_i, m)`` is essentially given by the probability distribution ``q_i(y_i\mid \mathcal{Y}_i)`` to measure ``y_i`` when the true value is ``\mathcal{Y}_i``, since ``\mathcal{Y}_i = m(\mathcal{X}_i,\lambda)`` by assumption of the model:
```math
\ell(y_i \mid \mathcal{X}_i, \lambda, x_i, m) = q_i(y_i\mid m(\mathcal{X}_i, \lambda))
```

## No ``x``-uncertainty
A convenient situation is, when the distinction between ``x_i`` and ``\mathcal{X}_i`` can be neglected, e.g. because the independent variable can be measured with high precision. Then ``p_i(\mathcal{X}_i\mid \lambda, x_i, m)`` becomes a Dirac distribution, and
```math
\ell(y_i\mid \lambda, x_i, m) = q(y_i\mid m(x_i,\lambda))
```
Hence, the posterior probability reads 
```math
p(\lambda \mid \{x_i\}_{i=1}^N, \{y_i\}_{i=1}^N, m) \propto  p_0(\lambda\mid \{x_i\}_{i=1}^N, m) \prod_{i_1}^n q_i(y_i\mid m(x_i,\lambda))
```

## Retrieving the LSQ objective
Using the aforementioned uniform prior ``p_0(\lambda\mid \{x_i\}_{i=1}^N, m) = 1`` and assuming normal distributions for ``q_i`` with standard deviations ``\Delta y_i`` leads to

```math
\begin{aligned}
p(\lambda \mid \{x_i\}_{i=1}^N, \{y_i\}_{i=1}^N, m) &\propto   \prod_{i_1}^n \frac{1}{\sqrt{2\pi}\Delta y_i}\exp\left(- \frac{(y_i - m(x_i,\lambda))^2}{2\Delta y_i}\right) \\ 
&\propto \prod_{i_1}^n \exp\left(- \frac{(y_i - m(x_i,\lambda))^2}{2\Delta y_i}\right)\\ 
& \quad = \exp\left(- \sum_{i=1}^N \frac{(y_i - m(x_i,\lambda))^2}{2\Delta y_i}\right)
\end{aligned}
```
Maximizing this function is equivalent to minimizing
```math
\sum_{i=1}^N \frac{(y_i - m(x_i,\lambda))^2}{2\Delta y_i}
```

which is the weighted least squares objective (up to a factor ``\frac{1}{2}``) function (see [Background: LSQ](@ref)).

