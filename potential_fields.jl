module Potential_fields
export potential_oscillator, potential_barrier
    
    Base.@kwdef mutable struct Potential_barrier_params
        init_params = init_params
        magnitude::Float64
        x_0::Float64
        width::Float64
        x_axis::Array{Float64} = [init_params.x_min, x_0 - 0.5width, x_0 - 0.5width, x_0 + 0.5width, x_0 + 0.5width, init_params.x_max]
        y_axis::Array{Float64} = [0, 0, 0.5*(magnitude/20), 0.5*(magnitude/20), 0, 0]
    end

    function potential_oscillator(init_params, mag)
        v = zeros(init_params.Nx)
        v .= 0.5*mag*(init_params.x).^2 
        return v
    end

    function potential_barrier(init_params, barriers)
        v = zeros(init_params.Nx)
        x_axis, y_axis = Float64[], Float64[] 
        for i in 1:length(barriers)
            barrier_init_params = Potential_barrier_params(init_params = init_params, magnitude = barriers[i][1], x_0 = barriers[i][2], width=barriers[i][3])

            total_length = init_params.x_max - init_params.x_min
            position_start = barrier_init_params.x_0 - .5*barrier_init_params.width - init_params.x_min
            position_stop = barrier_init_params.x_0 + .5*barrier_init_params.width - init_params.x_min
            coord_start = floor(position_start/total_length * init_params.Nx)
            coord_stop  = floor(position_stop/total_length * init_params.Nx)
            v[Int64(coord_start):Int64(coord_stop)] .+= barrier_init_params.magnitude
            append!(x_axis, barrier_init_params.x_axis)
            append!(y_axis, barrier_init_params.y_axis)
            println(x_axis)
        end
        return v, x_axis, y_axis
    end

end
