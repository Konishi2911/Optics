using LinearAlgebra
using Plots

struct SingleRay
    origin::Vector{Float64}
    direction::Vector{Float64}

    medium::Union{AbstractLens, Nothing}
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

    if t1 > 0 && t2 > 0
        # Sort the angles and points
        if arc.start_angle < theta1 < arc.end_angle && 
            arc.start_angle < theta2 < arc.end_angle

            if norm(pi1 - ray.origin) < norm(pi2 - ray.origin)
                return pi1
            else
                return pi2
            end
        elseif theta1 > arc.start_angle
            return pi1
        elseif theta2 < arc.end_angle
            return pi2
        else
            return nothing
        end
    elseif t1 > 0
        # Only t1 is valid
        if theta1 > arc.start_angle && theta1 < arc.end_angle
            return pi1
        else
            return nothing
        end
    elseif t2 > 0
        # Only t2 is valid
        if theta2 > arc.start_angle && theta2 < arc.end_angle
            return pi2
        else
            return nothing
        end
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