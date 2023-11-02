@testset "General methods" begin

	@testset "LogRange" begin
		# All arguments must be larger than 0.
		@test_throws DomainError LogRange(-1,10,5)
		@test_throws DomainError LogRange(1,-10,5)
		@test_throws DomainError LogRange(1,10,0)
		@test_throws DomainError LogRange(1,10,1,0)

		# LogRange(x,y,n) should create n values.
		@test length(LogRange(1,10,10)) == 10
		# LogRange(x,y,n) should start with x and end with y.
		@test LogRange(1,10,10)[1] == 1 && LogRange(1,10,10)[end] == 10
		# Default base is 10.
		@test LogRange(1,100,3) == [1.0,10.0,100.0]
		# Test base 2.
		@test LogRange(1,4,3,2) == [1.0,2.0,4.0]
	end










	@testset "peak_detection" begin
		grid = create_grid(LinRange(1,101,101))

		# Relative threshold must be in (0,1].
		@test_throws DomainError peak_detection(grid,0)
		@test_throws DomainError peak_detection(grid,-1)
		@test_throws DomainError peak_detection(grid,1.1)

		# fill = true despite all blocks belong to the peak (nothing left to be filled should not raise an error).
		@test peak_detection(grid, fill = true, volume_normalization = :none)[1] == [collect(1:100)]


		@testset "No volume normalization" begin
			# [0,...,0 , 0.25 , 0.5 , 1.0 , 0.5 , 0.25 , 0,...,0 , 0.25 , 0.5 , 1.0 , 0.5 , 0.25 , 0, ..., 0]
			weights = zeros(length(grid))
			weights[[25,75]] .= 1
			weights[[24,26,74,76]] .= 0.5
			weights[[23,27,73,77]] .= 0.25
			import_weights!(grid,weights)
			
			# Threshold 0.1 (find all non-zero):
			peak_indices, peak_domains = peak_detection(grid,0.1, fill = false, volume_normalization = :none)
			@test peak_indices == [collect(23:27), collect(73:77)]
			@test peak_domains == [[23.0,28.0], [73.0,78.0]]

			# Threshold 0.3 (find only [0.5,1,0.5]):
			peak_indices, peak_domains = peak_detection(grid,0.3, fill = false, volume_normalization = :none)
			@test peak_indices == [collect(24:26), collect(74:76)]
			@test peak_domains == [[24.0,27.0], [74.0,77.0]]

			# fill = true:
			peak_indices, peak_domains = peak_detection(grid,0.3, fill = true, volume_normalization = :none)
			# Using sets, as gaps are filled after peak groups have been added.
			@test Set(peak_indices) == Set([collect(1:23),collect(24:26), collect(27:73), collect(74:76), collect(77:100)])
			@test Set(peak_domains) == Set([[1.0,24.0],[24.0,27.0], [27.0,74.0], [74.0,77.0],[77.0,101.0]])
		end



		grid = create_grid(LogRange(1,101,101))
		centers, volumes, default_weights = export_all(grid)
		edges = [[centers[i]-0.5*volumes[i], centers[i] + 0.5*volumes[i]] for i in eachindex(centers)]


		@testset "Linear volume normalization" begin
			# [0,...,0 , 0.25 , 0.5 , 1.0 , 0.5 , 0.25 , 0,...,0 , 0.25 , 0.5 , 1.0 , 0.5 , 0.25 , 0, ..., 0]
			weights = zeros(length(grid))
			# Multiply weights by volumes to counteract the volume_normalization for the test.
			weights[[25,75]] .= 1 .* volumes[[25,75]]
			weights[[24,26,74,76]] .= 0.5 .* volumes[[24,26,74,76]]
			weights[[23,27,73,77]] .= 0.25 .* volumes[[23,27,73,77]]
			import_weights!(grid,weights)
			
			# Threshold 0.1 (find all non-zero):
			peak_indices, peak_domains = peak_detection(grid,0.1, fill = false, volume_normalization = :linear)
			@test peak_indices == [collect(23:27), collect(73:77)]
			@test peak_domains == [[edges[23][1],edges[27][end]], [edges[73][1],edges[77][end]]]

			# Threshold 0.3 (find only [0.5,1,0.5]):
			peak_indices, peak_domains = peak_detection(grid,0.3, fill = false, volume_normalization = :linear)
			@test peak_indices == [collect(24:26), collect(74:76)]
			@test peak_domains == [[edges[24][1],edges[26][end]], [edges[74][1],edges[76][end]]]

			# fill = true:
			peak_indices, peak_domains = peak_detection(grid,0.3, fill = true, volume_normalization = :linear)
			# Using sets, as gaps are filled after peak groups have been added.
			@test Set(peak_indices) == Set([collect(1:23),collect(24:26), collect(27:73), collect(74:76), collect(77:100)])
			@test Set(peak_domains) == Set([[edges[1][1],edges[23][end]],[edges[24][1],edges[26][end]], [edges[27][1],edges[73][end]], [edges[74][1],edges[76][end]],[edges[77][1],edges[100][end]]])
		end


		@testset "Logarithmic volume normalization" begin
			# [0,...,0 , 0.25 , 0.5 , 1.0 , 0.5 , 0.25 , 0,...,0 , 0.25 , 0.5 , 1.0 , 0.5 , 0.25 , 0, ..., 0]
			weights = zeros(length(grid))

			log_volumes = [log(edges[i][end])-log(edges[i][1]) for i in eachindex(centers)]
			# Multiply weights by log_volumes to counteract the volume_normalization for the test.
			weights[[25,75]] .= 1 .* log_volumes[[25,75]]
			weights[[24,26,74,76]] .= 0.5 .* log_volumes[[24,26,74,76]]
			weights[[23,27,73,77]] .= 0.25 .* log_volumes[[23,27,73,77]]
			import_weights!(grid,weights)
			
			# Threshold 0.1 (find all non-zero):
			peak_indices, peak_domains = peak_detection(grid,0.1, fill = false, volume_normalization = :log)
			@test peak_indices == [collect(23:27), collect(73:77)]
			@test peak_domains == [[edges[23][1],edges[27][end]], [edges[73][1],edges[77][end]]]

			# Threshold 0.3 (find only [0.5,1,0.5]):
			peak_indices, peak_domains = peak_detection(grid,0.3, fill = false, volume_normalization = :log)
			@test peak_indices == [collect(24:26), collect(74:76)]
			@test peak_domains == [[edges[24][1],edges[26][end]], [edges[74][1],edges[76][end]]]

			# fill = true:
			peak_indices, peak_domains = peak_detection(grid,0.3, fill = true, volume_normalization = :log)
			# Using sets, as gaps are filled after peak groups have been added.
			@test Set(peak_indices) == Set([collect(1:23),collect(24:26), collect(27:73), collect(74:76), collect(77:100)])
			@test Set(peak_domains) == Set([[edges[1][1],edges[23][end]],[edges[24][1],edges[26][end]], [edges[27][1],edges[73][end]], [edges[74][1],edges[76][end]],[edges[77][1],edges[100][end]]])
		end
	end
end