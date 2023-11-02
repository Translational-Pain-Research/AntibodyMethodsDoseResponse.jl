using AntibodyMethodsDoseResponse
import AntibodyMethodsDoseResponse as DR
using FittingObjectiveFunctions
using AdaptiveDensityApproximation
using Test
using Statistics

@testset "AntibodyMethodsDoseResponse.jl" begin

    include("Data.jl")
    include("Models.jl")

    @testset "Convenience methods" begin
        include("Convenience/GeneralMethods.jl")
        include("Convenience/AdaptiveFitting.jl")
    end

    @testset "Uncertainty estimation" begin
        include("UncertaintyEstimation/EpitopeUncertaintyAuxiliaryFunctions.jl")
        include("UncertaintyEstimation/EpitopeUncertainty.jl")
        include("UncertaintyEstimation/DoseResponseUncertainty.jl")
    end

end
