using Aqua, FastAlmostBandedMatrices, Test

@testset "Aqua Q/A" begin
    Aqua.test_all(FastAlmostBandedMatrices; ambiguities = false)
end
