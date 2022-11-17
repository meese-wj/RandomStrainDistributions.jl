
function iter_time(num = Int(1e6))
    yval::Float64 = one(Float64)
    for idx ∈ UnitRange(-num, num)
        yval += rand()
    end
    yval
end

@show @allocated iter_time()
@time iter_time()
