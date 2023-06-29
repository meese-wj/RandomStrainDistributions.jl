using Test
using RandomStrainDistributions

@testset "RandomStrainDistributions.jl" begin

    let 
        @testset "Testing Random Dislocations" begin
            ςon = 1/50
            Lx = 100
            Ly = Lx
            rsd = RandomDislocationDistribution(; concentration = ςon, Lx = Lx, Ly = Ly )
            all_dis = collect_dislocations(rsd)
            
            @testset "Topological charge sums to zero" begin
                # this is tested by making sure the burgers vectors sum to zero
                @test sum( burgersvector, all_dis[:] ) == zerovector( burgersvector(all_dis[1]) )
            end

            @testset "Burgers Vector selection" begin
                # this tests that the only Burgers vectors selected are those provided
                # in the tetragonal_burgers_vectors
                those_provided = true
                for rand_bv ∈ burgersvector.( all_dis[:] )
                    count = zero(Int)
                    for provided_bv ∈ tetragonal_burgers_vectors
                        count += rand_bv == provided_bv ? one(count) : zero(count)
                    end
                    
                    if count != one(count)
                        those_provided = false
                    end
                end

                @test those_provided
            end

            @testset "All origins are unique" begin
                # this makes sure no two dislocations sit on top of one another
                for idx ∈ eachindex(all_dis), jdx ∈ eachindex(all_dis)
                    if dislocationorigin(all_dis[idx]) == dislocationorigin(all_dis[jdx]) && idx != jdx
                        @show (idx, jdx)
                        @show all_dis[idx]
                        @show all_dis[jdx]
                    end
                end

                @test length( unique( dislocationorigin.(all_dis) ) ) == length(all_dis)
            end

            @testset "All origins are within the grid" begin
                # this makes sure no dislocations were generated outside the system
                # and that no dislocations are centered on the grid itself.
                all_inside = true
                no_grid_touching = true
                for origin ∈ dislocationorigin.(all_dis[:])
                    if !( 0 < xcomponent(origin) < Lx + one(Lx) && 0 < ycomponent(origin) < Ly + one(Ly) )
                        all_inside = false
                    end
                    if ( xcomponent(origin) - floor(xcomponent(origin)) ≈ zero(xcomponent(origin)) ) && ( ycomponent(origin) - floor(ycomponent(origin)) ≈ zero(ycomponent(origin)) )
                        no_grid_touching == false
                    end
                end
                @test all_inside
                @test no_grid_touching
            end


        end
    end

    let 
        @time @testset "Testing Non-Random Number of Dislocations" begin
            con = 0.001
            Lx = 100
            Ly = 100

            rdd = RandomDislocationDistribution(; concentration = con, Lx = Lx, Ly = Ly, random_defect_number = false)
            num_dislocs = Int[]
            for _ ∈ 1:10_000
                dislocs = collect_dislocations(rdd)
                push!(num_dislocs, length(dislocs))
            end
                
            @test length(unique(num_dislocs)) == 1
            @test unique(num_dislocs)[1] == con * Lx * Ly
        end        
    end

end