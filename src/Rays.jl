using LinearAlgebra
using Optim
using Plots
using Makie

struct SingleRay
    origin::Vector{Float64}
    direction::Vector{Float64}

    medium::Union{AbstractLens, Nothing}
end

struct PointSource
    position::Vector{Float64}
end

function intersect_point(ray::SingleRay, seg::Segment)
    q1 = seg.p1 - ray.origin
    q2 = seg.p2 - ray.origin
    n = normal(seg, 0.5)

    dotr = dot(ray.direction, n)
    dot1 = dot(q1, n)
    if dot1 * dotr < 1e-10
        return nothing
    end

    cross1 = ray.direction[1] * q1[2] - ray.direction[2] * q1[1]
    cross2 = ray.direction[1] * q2[2] - ray.direction[2] * q2[1]

    if cross1 * cross2 > 0
        return nothing
    else
        t = abs(cross1) / (abs(cross1) + abs(cross2))
        return seg.p1 + t * (seg.p2 - seg.p1)
    end
end

function intersect_point(ray::SingleRay, arc::Arc)
    p = ray.origin - arc.center
    pp = dot(p, p)
    qq = dot(ray.direction, ray.direction)
    pq = dot(p, ray.direction)

    det = pq^2 - qq * (pp - arc.radius^2)
    if det < 0
        return nothing
    end
    sqrt_det = sqrt(det)
    t1 = (-pq + sqrt_det) / qq
    t2 = (-pq - sqrt_det) / qq

    pi1 = ray.origin + t1 * ray.direction
    pi2 = ray.origin + t2 * ray.direction
    theta1 = atan(pi1[2] - arc.center[2], pi1[1] - arc.center[1])
    theta2 = atan(pi2[2] - arc.center[2], pi2[1] - arc.center[1])

    if arc.start_angle > pi || arc.end_angle > pi 
        if theta1 < 0
            theta1 += 2 * pi
        end
        if theta2 < 0
            theta2 += 2 * pi
        end
    end

    if t1 > 0 && t2 > 0
        if arc.start_angle < theta1 < arc.end_angle && 
            arc.start_angle < theta2 < arc.end_angle

            if norm(pi1 - ray.origin) < norm(pi2 - ray.origin)
                return pi1
            else
                return pi2
            end
        elseif arc.start_angle < theta1 < arc.end_angle
            return pi1
        elseif arc.start_angle < theta2 < arc.end_angle
            return pi2
        else
            return nothing
        end
    elseif t1 > 0
        # Only t1 is valid
        if arc.start_angle < theta1 < arc.end_angle
            return pi1
        else
            return nothing
        end
    elseif t2 > 0
        # Only t2 is valid
        if arc.start_angle < theta2 < arc.end_angle
            return pi2
        else
            return nothing
        end
    else
        return nothing
    end
end

function intersect_point(ray::SingleRay, curve::AsphericCurve)
    f = geom2d(curve)
    g(s) = norm(f(s[1]) - (ray.origin + s[2] * ray.direction))
    res = Optim.optimize(g, [0.0, 0.0])
    if res.f_val < 1e-8
        return res.minimizer[1]
    else
        return nothing
    end
end

function plot!(rays::Vector{SingleRay}, extend_length::Float64 = 100.0)
    for i in 1:length(rays)-1
        Plots.plot!(
            [rays[i].origin[1], rays[i+1].origin[1]], 
            [rays[i].origin[2], rays[i+1].origin[2]], 
            color=:blue,
            legend=false
        )
    end
    Plots.plot!(
        [rays[end].origin[1], rays[end].origin[1] + extend_length * rays[end].direction[1]], 
        [rays[end].origin[2], rays[end].origin[2] + extend_length * rays[end].direction[2]], 
        color=:blue,
        legend=false
    )
end

function plot!(axis, rays::Vector{SingleRay}, extend_length::Float64 = 100.0; linewidth=2)
    for i in 1:length(rays)-1
        Makie.lines!(
            axis,
            [rays[i].origin[1], rays[i+1].origin[1]], 
            [rays[i].origin[2], rays[i+1].origin[2]], 
            color=:blue,
            linewidth=linewidth,
        )
    end
    Makie.lines!(
        axis,
        [rays[end].origin[1], rays[end].origin[1] + extend_length * rays[end].direction[1]], 
        [rays[end].origin[2], rays[end].origin[2] + extend_length * rays[end].direction[2]], 
        color=:blue,
        linewidth=linewidth,
    )
end