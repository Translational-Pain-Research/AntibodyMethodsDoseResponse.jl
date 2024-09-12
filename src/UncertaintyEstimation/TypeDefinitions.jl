####################################################################################################
# Define uncertainty structs.
####################################################################################################

export EpitopeUncertainty, DoseResponseUncertainty

"""
	struct EpitopeUncertainty

Data type to store uncertainty estimates for the `weights` of a K_τ `grid`.

**Fields**

* `levels`: List of uncertainty levels. They can be sample quantiles or fractions of the best objective value, depending on the constructor.
* `lower`: Matrix of estimated lower bounds for the weights at the corresponding uncertainty levels (dimension order : [level, grid parameter index]).
* `upper`: Matrix of estimated upper bounds for the wights at the corresponding uncertainty levels (dimension order : [level, gird parameter index]).
* `lower_offset`: Vector of lower bounds for the offset parameter corresponding to the uncertainty levels.
* `upper_offset`: Vector of upper bounds for the offset parameter corresponding to the uncertainty levels.

----

----

**Default constructor**

	EpitopeUncertainty(levels, lower, upper, lower_offset= nothing, upper_offset = nothing)

----

----

**Construction from bin-wise shifting**

	EpitopeUncertainty(data::FittingData,grid::OneDimGrid, bins = collect(1:length(grid)); 
		keywords...)

Estimate uncertainty by shifting all grid weights uniformly, one bin at a time, while keeping the other bins fixed. Admissible weights (for a given level) are determined by calculating the objective-function value (objective function automatically generated) for the shifted weights.

The bins can be defined as vector of indices, e.g. `[[1,2,3],[4,5,6]]` or `[1,3,5]` which is converted to `[[1],[3],[5]]`. To obtain bin indices from the grid-domain, use [`select_indices`](https://translational-pain-research.github.io/AdaptiveDensityApproximation-documentation/api/#AdaptiveDensityApproximation.select_indices).

The following keywords are available:

* `levels = collect(0.1:0.1:1)`: The uncertainty levels as fractions of the best objective value.
* `steps::Integer = 10^4`: Number of intermediate shifts to be tested. The maximal shift range is determined automatically.
* `bisections::Integer = 10^2`: Number of interval-bisection steps to determine the maximal shift range.
* `volume_normalization = :none`: Changes the scaling of the weight shifts for each individual weight in a bin. Use `:none` to apply the same shift for each weight. Use `:linear` to scale the shifts with the interval volumes corresponding to the weights.  And use `:log` to scale the sifts with the visual interval volumes (as they appear in a logarithmic plot) corresponding to the weights.
* `options::AdaptiveOptions = AdaptiveOptions()`: The objective function is automatically generated, using the same construction as [`adaptive_dose_response_fit`](@ref). See [`AdaptiveOptions`](@ref) and [`adaptive_dose_response_fit`](@ref) for the details. To use the offset estimation from a fit result, add the estimated value to the [`AdaptiveOptions`](@ref) with the `offset` keyword: `options = AdaptiveOptions(other_options..., offset = estimated_offset)`.




----

----

**Construction from samples**

	EpitopeUncertainty(samples; keywords... ) 

Estimate the uncertainty as credibility intervals (symmetric interval around the median) from `samples` (drawn from a posterior distribution). The samples can either be passed as array of samples (parameters) or as matrix (order: [parameter index, sample index]).

The following keywords are available:

* `levels = collect(0.1:0.1:1)`: The uncertainty levels as quantiles.
* `offset::Bool = false`: If `true` the last parameter element is treated as offset parameter, otherwise, all parameters are treated as grid weights.
"""
struct EpitopeUncertainty
	levels::Vector{T} where T <: Real
	lower::Matrix{T} where T <: Real
	upper::Matrix{T} where T <: Real
	lower_offset::Union{Nothing,Vector{T}} where T <: Real
	upper_offset::Union{Nothing,Vector{T}} where T <: Real

	# No need for keyword arguments. If upper_offset is to be specified, lower_offset must be specified as well.
	function EpitopeUncertainty(levels,lower,upper, lower_offset=nothing, upper_offset=nothing)

		if size(lower,2) != size(upper,2)
			throw(DimensionMismatch("Parameter length for upper and lower differs!"))
		end

		if isnothing(lower_offset) || isnothing(upper_offset)
			if !(length(levels)==size(lower,1)==size(upper,1))
				throw(DimensionMismatch("Numbers of levels, of lower bound estimates and of upper bound estimates do not match!"))
			end
			return new(levels,lower,upper, nothing,nothing)
		else
			if !(length(levels)==size(lower,1)==size(upper,1)==length(lower_offset)==length(upper_offset))
				throw(DimensionMismatch("Numbers of levels, of lower bound estimates, of upper bound estimates, of lower offsets and of upper offsets do not match!"))
			end
			return new(levels,lower,upper, lower_offset,upper_offset)
		end
	end

end









"""
	struct DoseResponseUncertainty

Data type to store uncertainty estimates for the response values of a [`DoseResponseResult`](@ref) object.

**Fields**

* `levels`: List of uncertainty levels, corresponding to the uncertainty levels of a [`EpitopeUncertainty`](@ref) object (if constructed from an `EpitopeUncertainty` object).
* `concentration`: Vector of concentrations (no uncertainty).
* `lower`: Matrix of estimated lower bounds for the responses at the corresponding uncertainty level (dimension order : [level, concentration/response index]).
* `upper`: Matrix of estimated upper bounds for the responses at the corresponding uncertainty level (dimension order : [level, concentration/response index]).

----

----

**Default constructor**

	DoseResponseUncertainty(levels,concentrations,lower,upper)

----

----

**Construction from an EpitopeUncertainty object**

	DoseResponseUncertainty(grid::OneDimGrid,
		eu::EpitopeUncertainty,
		concentrations::AbstractVector; 
		keywords...
	)

Estimate the dose-response uncertainty from an [`EpitopeUncertainty`](@ref) object `eu` for the provided `concentrations`. The `grid` should be the grid that was used to create the `EpitopeUncertainty` object `eu`.

The following keywords are available:

* `bins = [collect(1:length(grid))]`: The response bounds are calculated as point-wise minima/maxima of responses created from the grid weights, where one bin at a time is replaced with the [`EpitopeUncertainty`](@ref) lower and upper bound, while keeping the other weights fixed. For the minima/maxima all response values, iterating over all bins, are considered. Ideally, the bins should correspond to the bins that were used to construct the [`EpitopeUncertainty`](@ref) object `eu`.
* `model::Function = accumulation_model`: The model that is used to calculate the response values. The available model functions are [`accumulation_model`](@ref), [`accumulation_inv_const_model`](@ref), [`langmuir_model`](@ref) and [`langmuir_inv_const_model`](@ref).

There is no `offset` keyword, as the offsets are determined by the [`EpitopeUncertainty`](@ref) object.
"""
struct DoseResponseUncertainty
	levels::Vector{T} where T <: Real
	concentrations::Vector{T} where T <: Real
	lower::Matrix{T} where T <: Real
	upper::Matrix{T} where T <: Real

	function DoseResponseUncertainty(levels,concentrations,lower,upper)

		if !(length(concentrations) == size(lower,2) == size(upper,2))
			throw(DimensionMismatch("Lengths of concentrations, of lower bound responses and of upper bound responses do not match!"))
		end

		if !(length(levels)==size(lower,1)==size(upper,1))
			throw(DimensionMismatch("Number of levels, of lower bound estimates and of upper bound estimates do not match!"))
		end

		return new(levels,concentrations,lower,upper)
	end
end