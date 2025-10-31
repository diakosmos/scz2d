using Test 
using SCZ2D
using SCZ2D.MyCheb 

@testset "MyCheb.jl" begin

    a = SCZ2D.MyCheb.mycheb(5)
    @test a[1] == 1 
    @test a[2] ≈ 1/√2
    @test a[3] == 0.0
    @test a[4] ≈ -1/√2 
    @test a[5] == -1
    #=
       @testset "Basic functionality" begin
           @test my_function(2, 3) == 5
           @test_throws DomainError my_function(-1, 0)
       end
       
       @testset "Edge cases" begin
           @test isapprox(my_other_function(1e-10), 0.0, atol=1e-9)
       end
    =#

end

