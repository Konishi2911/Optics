using Test
import Optics

@testset "Rays.Segment" begin
    ray = Optics.SingleRay([0.0, 0.0], [1.0, 1.0], nothing)
    seg = Optics.Segment([1.0, 0.0], [0.0, 1.0])
    @test Optics.intersect_point(ray, seg) == [0.5, 0.5]

    ray2 = Optics.SingleRay([1.0, 1.0], [1.0, 1.0], nothing)
    @test Optics.intersect_point(ray2, seg) === nothing
end

@testset "Rays.Arc" begin
    ray = Optics.SingleRay([0.0, 0.0], [1.0, 0.0], nothing)
    arc = Optics.Arc([0.0, 0.0], 1.0, -π, π/2)

    @test Optics.intersect_point(ray, arc) == [1.0, 0.0]
    @test Optics.find_local(arc, [1.0, 0.0]) == 2/3
end