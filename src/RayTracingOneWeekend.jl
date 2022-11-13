import FromFile: @from
using Grassmann, ProgressMeter, StaticArrays, LinearAlgebra
using IterTools
Point, Color = (SVector{3, Float64}, SVector{3, Float64})
PPM_FILE = "/Users/alokbeniwal/.julia/dev/RayTracingOneWeekend/target/raytraced.ppm"

normsq(x) = norm(x)^2

# start, end
Ray = Tuple{Point, Point}

at(ray, t) = Ray((ray[1], t * ray[2]))

function color(ray::Ray)
    t = (normalize(ray[2])[2] + 1.0) / 2 # unit y direction
    # TODO: ones is white, other is blue
    white, blue = ones(Color), Color([0.5, 0.7, 1.0])
    (1 - t) * white + t * blue
end

# TODO: intersect ray and viewport (line meet bivector)
function main()
    aspect_ratio = 16 / 9
    img_width = 400
    img_height = Int(img_width ÷ aspect_ratio)
    viewport_height = 2.0
    viewport_width = viewport_height * aspect_ratio
    focal_len = 1.0
    # rand_img = rand(Color, (img_height, img_width))
    img_width,img_height=256,256
    initial_scene = let
        mat = zeros(Color, (img_height, img_width))
        # XXX: Julia is column major, so outer loop is *first*.
        for col in 1:img_width, row in 1:img_height
            r, g, b = (col-1) / (img_width-1), (row-1) / (img_height-1), 0.25
            mat[row, col] = Color([r, g, b])
        end
        mat
    end

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
            for row in (rows:-1:1), col in 1:cols
                # TODO: check (u,v), since not 0 indexed
                # u, v = (row - 1) / (rows - 1), (col - 1) / (cols - 1)
                # ray = Ray((zero(Point), lower_left_corner + u * horizontal + v * vertical))
                # c = color(ray)
                write_color(io, matrix[row, col])
            end
        end
    end

    write_ppm(initial_scene)
end

main()
