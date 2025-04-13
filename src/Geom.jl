using Optim
using ForwardDiff

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

struct AsphericCurve <: AbstractGeometricElement
    diameter::Float64
    offset::Float64
    coeffs::Vector
    conic_constant::Float64
    radius::Float64

    is_mirrored::Bool
end

function AsphericCurve(diameter::Float64, coeffs::Vector{Float64}, conic_constant::Float64, radius::Float64; is_mirrored::Bool = false)
    return AsphericCurve(diameter, 0.0, coeffs, conic_constant, radius, is_mirrored)
end

function AsphericCurve(curve::AsphericCurve, offset::Float64; is_mirrored::Bool = curve.is_mirrored)
    return AsphericCurve(curve.diameter, offset, curve.coeffs, curve.conic_constant, curve.radius, is_mirrored)
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

function geom2d(curve::AsphericCurve)
    fn = (t) -> begin
        r = curve.radius
        k = curve.conic_constant
        s = (t - 0.5) * curve.diameter

        y = s
        x = s^2 / (r * (1 + sqrt(1 - (1 + k) * s^2 / r^2)))
        for i in 1:length(curve.coeffs)
            x += curve.coeffs[i] * s^(2*(i+1))
        end

        if curve.is_mirrored
            return -x + curve.offset, y
        else
            return x + curve.offset, y
        end
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

function normal(curve::AsphericCurve, t::Float64)
    f = geom2d(curve)
    fx(t) = f(t)[1]
    fy(t) = f(t)[2]

    gradx(t) = ForwardDiff.derivative(fx, t)
    grady(t) = ForwardDiff.derivative(fy, t)
    tan = [gradx(t), grady(t)] / norm([gradx(t), grady(t)])
    return [-tan[2], tan[1]]
end

function find_local(geom::AbstractGeometricElement, x::Vector{Float64})
    f(t) = (x[1] - geom2d(geom)(t)[1])^2 + (x[2] - geom2d(geom)(t)[2])^2
    res = Optim.optimize(f, 0.0, 1.0)    
    
    return Optim.minimizer(res)
end