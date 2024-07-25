using Documenter, AntibodyMethodsDoseResponse, AntibodyMethodsDoseResponseConvenience


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



"API"=> ["AntibodyMethodsDoseResponse" => "API.md" ,
		"AntibodyMethodsDoseResponseConvenience" => "API_Convenience.md",
		] ,
])