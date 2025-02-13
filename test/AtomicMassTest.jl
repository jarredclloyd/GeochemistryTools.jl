@testset "Atomic Mass" begin
    @test get_atomicmass('H').val ≈ 1.008
    @test get_atomicmass("He").val ≈ 4.002602
    @test get_atomicmass("Lithium").val ≈ 6.94
    @test get_atomicmass("beryllium").val ≈ 9.0121831
    @test get_atomicmass(:B).val ≈ 10.81
    @test get_atomicmass(6).val ≈ 12.011
end
