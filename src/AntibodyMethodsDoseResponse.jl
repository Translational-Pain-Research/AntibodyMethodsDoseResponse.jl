module AntibodyMethodsDoseResponse

	# Dependencies from standard registry.
	using Statistics

	# Dependencies not in the standard registry.
	using FittingObjectiveFunctions
	using AdaptiveDensityApproximation


	include("Data.jl")
	include("Models.jl")


	include("Convenience/GeneralMethods.jl")
	include("Convenience/AdaptiveFitting.jl")

	include("UncertaintyEstimation/InternalMethods.jl")
	include("UncertaintyEstimation/TypeDefinitions.jl")
	include("UncertaintyEstimation/EpitopeUncertaintyConstructors.jl")
	include("UncertaintyEstimation/DoseResponseUncertaintyConstructors.jl")

end
