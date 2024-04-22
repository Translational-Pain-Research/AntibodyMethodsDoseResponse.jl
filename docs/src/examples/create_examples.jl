using AntibodyMethodsDoseResponseConvenience
using Distributions, Random

cd(@__DIR__)


density(k) = pdf(Normal(1e-6,5e-7),k) + pdf(Normal(1e-4,5e-5),k)
simulation_grid = create_grid(LogRange(1e-10,1e-2,100))
approximate_density!(simulation_grid,density, volume_normalization = true)

concentrations = LogRange(1e-10,1e-2,16)
responses = DoseResponseResult(simulation_grid,concentrations).responses

Random.seed!(1234)
deviations = rand(Truncated(Normal(0,0.1),0,Inf),16)
replicate_1 = responses .+ deviations .* responses

Random.seed!(4321)
deviations = rand(Truncated(Normal(0,0.1),0,Inf),16)
replicate_2 = responses .+ deviations .* responses

Random.seed!(1122)
deviations = rand(Truncated(Normal(0,0.1),0,Inf),16)
replicate_3 = responses .+ deviations .* responses



condition = FittingCondition(concentrations, replicate_1,replicate_2,replicate_3, path = pwd()*"/default_concentrations", scale = 500)


results = fit_condition(condition)



u_concentrations = LogRange(extrema(condition.data.independent)...,100)

bins, ranges = peak_detection(results.grid,0.01)
eu = EpitopeUncertainty(condition.data,results.grid,bins, levels = [1e-10,1e-5,0.1,0.25,0.5,0.75,0.9,1])
du = DoseResponseUncertainty(results.grid,eu, u_concentrations, bins = bins)

mkpath("uncertainty")
serialize("uncertainty/eu.jld", eu)
serialize("uncertainty/du.jld", du)



bins, ranges = peak_detection(results.grid,0.01)
eu = EpitopeUncertainty(condition.data,results.grid,bins, levels = [1e-10,1e-5,0.1,0.25,0.5,0.75,0.9,1], volume_normalization = :linear)
du = DoseResponseUncertainty(results.grid,eu, u_concentrations, bins = bins)

mkpath("uncertainty_volume")
serialize("uncertainty_volume/eu.jld", eu)
serialize("uncertainty_volume/du.jld", du)



bins, ranges = peak_detection(results.grid,0.01)
eu = EpitopeUncertainty(condition.data,results.grid,bins, levels = [1e-10,1e-5,0.1,0.25,0.5,0.75,0.9,1], volume_normalization = :log)
du = DoseResponseUncertainty(results.grid,eu, u_concentrations, bins = bins)

mkpath("uncertainty_log_volume")
serialize("uncertainty_log_volume/eu.jld", eu)
serialize("uncertainty_log_volume/du.jld", du)


temp_bins, temp_ranges = peak_detection(results.grid,0.01)
bins = [temp_bins[2]]
eu = EpitopeUncertainty(condition.data,results.grid,bins, levels = [1e-10,1e-5,0.1,0.25,0.5,0.75,0.9,1], options = condition.options_2)
du = DoseResponseUncertainty(results.grid,eu, u_concentrations, bins = bins)

mkpath("uncertainty_single")
serialize("uncertainty_single/eu.jld", eu)
serialize("uncertainty_single/du.jld", du)


temp_bins, temp_ranges = peak_detection(results.grid,0.01)
bins = [elm for elm in temp_bins[2]]
eu = EpitopeUncertainty(condition.data,results.grid,bins, levels = [1e-10,1e-5,0.1,0.25,0.5,0.75,0.9,1], options = condition.options_2)
du = DoseResponseUncertainty(results.grid,eu, u_concentrations, bins = bins)

mkpath("uncertainty_split")
serialize("uncertainty_split/eu.jld", eu)
serialize("uncertainty_split/du.jld", du)


bins = collect(1:length(results.grid))
eu = EpitopeUncertainty(condition.data,results.grid,bins, levels = [1e-10,1e-5,0.1,0.25,0.5,0.75,0.9,1], options = condition.options_2)
du = DoseResponseUncertainty(results.grid,eu, u_concentrations, bins = bins)

mkpath("uncertainty_total")
serialize("uncertainty_total/eu.jld", eu)
serialize("uncertainty_total/du.jld", du)


bins = [collect(1:length(results.grid))]
eu = EpitopeUncertainty(condition.data,results.grid,bins, levels = [1e-10,1e-5,0.1,0.25,0.5,0.75,0.9,1], options = condition.options_2)
du = DoseResponseUncertainty(results.grid,eu, u_concentrations, bins = bins)

mkpath("uncertainty_global")
serialize("uncertainty_global/eu.jld", eu)
serialize("uncertainty_global/du.jld", du)


grid = create_grid(LogRange(1e-10,1e-2,40))

default_minimizer_1 = function(f,∇f,init)
	lower = zeros(length(init))
	upper = [Inf for i in 1:length(init)]
	return optimize(f,lower,upper, init, Fminbox(NelderMead()),Optim.Options(g_tol = 1e-12, iterations =2000)).minimizer
end

result = adaptive_dose_response_fit(grid, condition.data, default_minimizer_1)

mkpath("direct_fit")
serialize("direct_fit/results.jld", result)


