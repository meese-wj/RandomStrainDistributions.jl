
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
            burgers_vectors = [ Vector2D(1., 0), Vector2D(-1., 0), Vector2D(0, 1.), Vector2D(0, -1.) ]
            ςon = 1/50
            Lx = 100
            rsd = RandomDislocationDistribution(; expected_num_dislocations = round(Int, ςon * Lx^2), Lx = Lx, burgers_vectors = burgers_vectors )
            all_dis = collect_dislocations(rsd)
            
            @testset "Topological charge sums to zero" begin
                # this is tested by making sure the burgers vectors sum to zero
                @test isequal( sum( all_dis[Int(BurgersVector), :] ), Vector2D(0., 0.) )
            end

            @testset "All origins are unique" begin
                # this makes sure no two dislocations sit on top of one another
                @test length( unique_origins( all_dis[Int(DislocationOrigin), :] ) ) == size(all_dis)[2]
            end

        end
    end

end