
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
                @test isequal( sum( burgersvector.( all_dis[:] ) ), Vector2D(0., 0.) )
            end

            @testset "Burgers Vector selection" begin
                # this tests that the only Burgers vectors selected are those provided
                # in the tetragonal_burgers_vectors
                those_provided = true
                for rand_bv ∈ burgersvector.( all_dis[:] )
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
                @test length( unique_origins( dislocationorigin.(all_dis[:]) ) ) == length(all_dis)
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

end