import FromFile: @from
using Grassmann

const PPM_FILE = "/Users/alokbeniwal/.julia/dev/RayTracingOneWeekend/target/raytraced.ppm"


function rand_img()
    rand_color() = rand(0:255)
    img = map(_ -> (rand_color(), rand_color(), rand_color()), zeros(256, 256))
    img
end

function write_ppm(matrix, filename = PPM_FILE)
    rows, cols = size(matrix)
    open(filename, "w") do io
        write(
            io,
            """
      P3
      $(cols) $(rows)
      """,
        )
        for (r, g, b) in matrix # col major reading, vertical
            write(io, "$r $g $b \n")
        end
    end
end

write_ppm(rand_img())
