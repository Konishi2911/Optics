using Optim

abstract type AbstractGeometricElement end

struct Segment <: AbstractGeometricElement
    p1::Vector{Float64}
    p2::Vector{Float64}
end

struct Arc <: AbstractGeometricElement
    center::Vector{Float64}
    radius::Float64
    start_angle::Float64
    end_angle::Float64
end


function geom2d(seg::Segment)
    fn = (t) -> begin
        x = seg.p1[1] + (seg.p2[1] - seg.p1[1]) * t
        y = seg.p1[2] + (seg.p2[2] - seg.p1[2]) * t
        return x, y
    end
    return fn    
end

function geom2d(arc::Arc)
    fn = (t) -> begin
        theta = arc.start_angle + (arc.end_angle - arc.start_angle) * t
        x = arc.center[1] + arc.radius * cos(theta)
        y = arc.center[2] + arc.radius * sin(theta)
        return x, y
    end
    return fn    
end

function normal(seg::Segment, t::Float64)
    v = seg.p2 - seg.p1
    v = v / norm(v)
    return [-v[2], v[1]]
end

function normal(arc::Arc, t::Float64)
    theta = arc.start_angle + (arc.end_angle - arc.start_angle) * t
    nx = cos(theta)
    ny = sin(theta)
    return [nx, ny]
end

function find_local(geom::AbstractGeometricElement, x::Vector{Float64})
    f(t) = (x[1] - geom2d(geom)(t)[1])^2 + (x[2] - geom2d(geom)(t)[2])^2
    res = Optim.optimize(f, 0.0, 1.0)    
    
    return Optim.minimizer(res)
end