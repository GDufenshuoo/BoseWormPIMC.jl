using BoseWormPIMC
using Test

@testset "BoseWormPIMC.jl" begin
    @test BoseWormPIMC.pimc("MC.in")
end
