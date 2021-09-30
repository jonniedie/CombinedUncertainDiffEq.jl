using CombinedUncertainDiffEq
using CombinedUncertainDiffEq: Uniform, to_distribution, to_interval
using Test

@testset "Conversions" begin
    i = 1..2
    d = Uniform(1, 2)

    # Individual
    @test to_distribution(i) == d
    @test to_distribution(d) == d
    @test to_interval(i) == i
    @test to_interval(d) == i

    # Round trip
    @test to_interval(to_distribution(i)) == i
    @test to_distribution(to_interval(d)) == d
end
