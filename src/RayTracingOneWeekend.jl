import FromFile: @from
using Grassmann, ProgressMeter, StaticArrays, LinearAlgebra
using IterTools
Point, Color = (SVector{3, Float64}, SVector{3, Float64})
red, white, blue = Color(1, 0, 0), Color(1, 1, 1), Color(0.5, 0.7, 1.0)
PPM_FILE = "/Users/alokbeniwal/.julia/dev/RayTracingOneWeekend/target/raytraced.ppm"

normsq(x) = norm(x)^2

# hittable
abstract type Solid end

struct Sphere{R <: Real} <: Solid
    center::Point
    radius::R
end

struct HitRecord
    pt::Point
    normal::Point
    t::Float64
    front_face::Bool
end

# Ray = Tuple{Point, Point} # start, end
struct Ray
    origin::Point
    direction::Point # really a vector, affine space
end

at(ray::Ray, t::Real) = Ray(ray.origin, t * ray.direction)

# TODO: each ray can be done totally parallel to the others
function color(
    world::AbstractArray{<:Solid},
    ray::Ray;
    t_min = eps(Float64),
    t_max = Inf,
)::Color
    hit_record = hit(world, ray; t_min = t_min, t_max = t_max)
    if hit_record !== nothing
        return (hit_record.normal + white) / 2
    end
    t = (normalize(ray.direction).y + 1.0) / 2
    (1 - t) * white + t * blue
end

function test_img(h, w)::Matrix{Color}
    mat = zeros(Color, (h, w))

    for row in h:-1:1, col in 1:w # Column major, so outer loop is *first*.
        r, g, b = (col - 1) / (w - 1), (row - 1) / (h - 1), 0.25
        mat[row, col] = Color(r, g, b)
    end
    mat
end

# return record or nothing
function hit(sphere::Sphere, ray::Ray; t_min = -Inf, t_max = Inf)::Union{HitRecord, Nothing}
    oc = ray.origin - sphere.center
    a = normsq(ray.direction)
    half_b = oc ⋅ ray.direction
    c = normsq(oc) - sphere.radius^2
    discriminant = half_b^2 - a * c
    if discriminant < 0
        return nothing
    end

    root = -(half_b + sqrt(discriminant)) / a
    if root < t_min || root > t_max
        # try other root
        root = -(half_b - sqrt(discriminant)) / a
        if root < t_min || root > t_max
            return nothing
        end
    end

    println(ray)
    println(at(ray,root))
    println(sphere)
    pt = at(ray, root).direction
    normal = (pt - sphere.center) / sphere.radius
    println(normal)
    front_face = normal ⋅ ray.direction < 0

    HitRecord(pt, front_face ? normal : -normal, root, front_face)
end

function hit(
    objects::AbstractArray{<:Solid},
    ray::Ray;
    t_min = -Inf,
    t_max = Inf,
)::Union{HitRecord, Nothing}
    closest, closest_t = nothing, t_max
    for object in objects
        hit_record = hit(object, ray; t_min = t_min, t_max = closest_t)
        if hit_record !== nothing && hit_record.t < closest_t
            closest_t, closest = hit_record.t, hit_record
        end
    end

    closest
end

function main()
    # Image
    aspect_ratio = 16 / 9
    img_width = 400
    img_height = Int(img_width ÷ aspect_ratio)
    viewport_height = 2.0
    viewport_width = viewport_height * aspect_ratio
    focal_len = 1.0

    # World
    world = [Sphere(Point(0, 0, -1), 0.5), Sphere(Point(0, -100.5, -1), 100)]

    # Camera
    origin = zero(Point)

    horizontal = Point([viewport_width, 0, 0])
    vertical = Point([0, viewport_height, 0])
    inwards = Point(0, 0, focal_len)

    lower_left_corner = origin - horizontal / 2 - vertical / 2 - inwards

    function write_ppm(matrix::Matrix{Color}, filename = PPM_FILE)
        rows, cols = size(matrix)
        function write_color(io, c)
            r, g, b = Int.(floor.(255.999 * c))
            write(io, "$r $g $b \n")
        end
        open(filename, "w") do io
            write(
                io,
                """
          P3
          $(cols) $(rows)
          255
          """,
            )
            for row in img_height:-1:1, col in 1:img_width
                # TODO: check (u,v), since not 0 indexed
                u, v = (col - 1) / (img_width - 1), (row - 1) / (img_height - 1)
                ray = Ray(origin, lower_left_corner + u * horizontal + v * vertical)
                c = color(world, ray)
                write_color(io, c)
            end
        end
    end

    write_ppm(test_img(img_height, img_width))
end

main()

# TODO: intersect ray and viewport (line meet bivector)
# TODO: can broadcast to color in parallel
