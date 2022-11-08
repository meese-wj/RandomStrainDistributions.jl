
function unique_origins(origins)
    distinct_origins = eltype(origins)[]
    for (idx, origin) ∈ enumerate(origins)
        counts = zero(Int)
        for distinct ∈ distinct_origins
            boo = isequal(origin, distinct)
            counts += boo ? one(counts) : zero(counts)
        end
        if counts == zero(counts)
            push!(distinct_origins, Vector2D(origin))
        end
    end
    return distinct_origins
end

@testset "RandomStrainDistributions.jl" begin

    let 
        @testset "Testing Random Dislocations" begin
            ςon = 1/50
            Lx = 100
            Ly = Lx
            rsd = RandomDislocationDistribution(; expected_num_dislocations = round(Int, ςon * Lx^2), Lx = Lx, Ly = Ly )
            all_dis = collect_dislocations(rsd)
            
            @testset "Topological charge sums to zero" begin
                # this is tested by making sure the burgers vectors sum to zero
                @test isequal( sum( all_dis[Int(BurgersVector), :] ), Vector2D(0., 0.) )
            end

            @testset "Burgers Vector selection" begin
                # this tests that the only Burgers vectors selected are those provided
                # in the tetragonal_burgers_vectors
                those_provided = true
                for rand_bv ∈ all_dis[Int(BurgersVector), :]
                    count = zero(Int)
                    for provided_bv ∈ tetragonal_burgers_vectors
                        count += isequal(rand_bv, provided_bv) ? one(count) : zero(count)
                    end
                    
                    if count != one(count)
                        those_provided = false
                    end
                end

                @test those_provided
            end

            @testset "All origins are unique" begin
                # this makes sure no two dislocations sit on top of one another
                @test length( unique_origins( all_dis[Int(DislocationOrigin), :] ) ) == size(all_dis)[2]
            end

            @testset "All origins are within the grid" begin
                # this makes sure no dislocations were generated outside the system
                # and that no dislocations are centered on the grid itself.
                all_inside = true
                no_grid_touching = true
                for origin ∈ all_dis[Int(DislocationOrigin), :]
                    if !( 0 < origin.vec[1] < Lx + one(Lx) && 0 < origin.vec[2] < Ly + one(Ly) )
                        all_inside = false
                    end
                    if ( origin.vec[1] - floor(origin.vec[1]) ≈ zero(origin.vec[1]) ) && ( origin.vec[2] - floor(origin.vec[2]) ≈ zero(origin.vec[2]) )
                        no_grid_touching == false
                    end
                end
                @test all_inside
                @test no_grid_touching
            end


        end
    end

end