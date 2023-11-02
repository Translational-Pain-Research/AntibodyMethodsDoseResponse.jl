using Documenter, FittingObjectiveFunctions, AdaptiveDensityApproximation, AntibodyMethodsDoseResponse, AntibodyMethodsDoseResponseConvenience

FOF = "FittingObjectiveFunctions"
ADA = "AdaptiveDensityApproximation"

makedocs(sitename="AntibodyMethodsDoseResponse", pages = [
"Introduction"=>"index.md" ,
"Quick start" => "QuickStart.md",
"Detailed explanations" => ["Measurement data" => "MeasurementData.md",
							"Models" => "Models.md",
							"Results and simulations" => "ResultsAndSimulations.md",
							"Density plots for `OneDimGrid`" => "DensityPlots.md",
							"Fitting" => "Fitting.md",
							"Uncertainty estimation" => "Uncertainty.md"
							],

"FittingObjectiveFunctions" => [
	"Introduction"=>"$FOF/index.md" ,
	"FittingData and ModelFunctions"=>"$FOF/fitting_data.md",
	"Least squares objective"=>["Background"=>"$FOF/lsq_background.md","How to implement" => "$FOF/lsq_implementation.md"],
	"Posterior probability"=> ["Background"=>"$FOF/posterior_background.md", "How to implement"=>"$FOF/posterior_implementation.md"],
	"Logarithmic posterior probability"=>["Background"=>"$FOF/log_posterior_background.md", "How to implement"=>"$FOF/log_posterior_implementation.md"],
],

"AdaptiveDensityApproximation"=>"$ADA/index.md" ,

"API"=> ["AntibodyMethodsDoseResponse" => "API.md" ,
		"AntibodyMethodsDoseResponseConvenience" => "API_Convenience.md",
		"FittingObjectiveFunctions"=>"FittingObjectiveFunctions/API.md" ,
		"AdaptiveDensityApproximation"=>"AdaptiveDensityApproximation/api.md"] ,
])