include("potential_fields.jl")
include("FD_methods.jl")
using GLMakie
#=

User defined domain and wavefunction parameters.

The following parameters must be explicitly defined: Nx, Nt, x_min, x_max, t_max, k, σ, x_0

The spatial domain is initialized by defining the number of spatial step sizes, Nx, and the boundary/end points of 
the domain, x_min and x_max. The temporal domain is established by setting the total number of time steps, Nt, and the total 
simulated time, t. 
The time and spatial step sizes, dt and dx, and the 1D time and space arrays, t and x, which denote all the discrete temporal 
and spatial grid points are generated from the domain parameters.

The gaussian wavelet, ψ, is a 1D complex vector/array generated from the k, σ and x_0 wavefunction 
parameters and the x domain parameter.
k   = wavenumber
σ   = spatial width/spread of the packet
x_0 = coordinate the packet is centered at
=#
println("Hello, world")
Base.@kwdef mutable struct Params
    Nx::Int64
    Nt::Int64
    x_min::Float64
    x_max::Float64
    t_max::Float64
    
    dt::Float64       = t_max/(Nt-1)
    dx::Float64       = (abs(x_min)+abs(x_max))/(Nx-1)
    t::Array{Float64} = collect(0:dt:t_max)
    x::Array{Float64} = collect(x_min:dx:x_max)

    k::Float64
    σ::Float64
    x_0::Float64

    # Initial wavefunction shape at t=0
    ψ::Array{ComplexF64}    = (((2π*(σ^2))^(-0.25)).*exp.((-(x .- x_0).^2)/(4*σ^2)+(im*k.*x)))
    ψ_rk::Array{ComplexF64} = (((2π*(σ^2))^(-0.25)).*exp.((-(x .- x_0).^2)/(4*σ^2)+(im*k.*x)))
    
    potential_field::Vector{Float64} = zeros(Int64, Nx)
end

#################################################################################################
#################################################################################################
#################################################################################################

fig = Figure(backgroundcolor = RGBf0(0.98, 0.98, 0.98),resolution = (1200, 800))
ax1 = Axis(fig[1, 1])
ax2 = Axis(fig[2, 1])
ax3 = Axis(fig[2, 2])
ax_node = Axis(fig[2,3], yscale=log10, title="Reflected", xlabel = "Time(s)") 
ax_max = Axis(fig[1, 3])

ylims!(ax1, -1, 1), ylims!(ax2, -1, 1), ylims!(ax3, -.05, 1.5)


scene, layout = layoutscene()

layout[1, 1] = Label(fig, "x,x")
layout[1, 2] = x_range = IntervalSlider(fig, range = -30:1:30, startvalues = (-5, 5), width=250)
layout[1, 3] = Label(fig, @lift string($(x_range.interval)))

x_0_slider   = labelslider!(fig, "x_0:", -30:1:30; format = x -> "$(x)", startvalue=4)
k_slider     = labelslider!(fig, "k:", -15:0.1:15; format = x -> "$(x)")
σ_slider     = labelslider!(fig, "σ:", 0.1:0.1:10; format = x -> "$(x)")
Nx_slider    = labelslider!(fig, "Nx:", 10:1:2048; format = x -> "$(x)")
Nt_slider    = labelslider!(fig, "Nt:", 10:1:512; format = x -> "$(x)")
t_max_slider = labelslider!(fig, "t:", 0.1:0.1:10; format = x -> "$(x)")

fig[1,2] = vgrid!(layout,
                x_0_slider.layout,
                k_slider.layout,
                σ_slider.layout,
                Nx_slider.layout,
                Nt_slider.layout,
                t_max_slider.layout)

menu = Menu(fig[3, 1], options = ["Crank Nicolson - ABC on", "Crank Nicolson - ABC off", "RK4"])
# menu_potentials = Menu(fig[3, 2], options = ["Oscillator", "Barrier"])
# menu_barrier = Menu(fig[3, 3], options = ["Oscillator", "Barrier"])

toggle = Toggle(fig[3,2])
fig[3,2] = grid!(hcat(Label(fig, "Initialise") , toggle ,  Label(fig, "Animate") ))

init_params =  @lift Params(Nx=Int64($(Nx_slider.slider.value)), 
                            Nt=Int64($(Nt_slider.slider.value)), 
                            σ=$(σ_slider.slider.value),
                            x_min=$(x_range.interval)[1], 
                            x_max=$(x_range.interval)[2], 
                            x_0=$(x_0_slider.slider.value),
                            k=$(k_slider.slider.value), 
                            t_max=$(t_max_slider.slider.value))      

@lift begin
    empty!(ax1), empty!(ax2), empty!(ax3), empty!(ax_max)
    lines!(fig[1,1], init_params[].x, real.($(init_params).ψ), color=:red,  axis = (xlabel="x", ylabel = "y",))
    lines!(fig[2,1], init_params[].x, imag.($(init_params).ψ), color=:blue,  axis = (xlabel="x", ylabel = "y",))
    lines!(fig[2,2], init_params[].x, abs2.($(init_params).ψ), color=:black,  axis = (xlabel="x", ylabel = "y",))
    lines!(fig[1,3], init_params[].x, real.($(init_params).ψ), color=:red,  axis = (xlabel="x", ylabel = "y",))
end

function animate_plot()

    if menu.selection[] == "Crank Nicolson - ABC on"
        BC_par =  FDM.CN_boundary_conditions(init_params = init_params[], ABC = "true", v1=init_params[].k, v2=init_params[].k)
        CN = FDM.CN_matrix(init_params[], BC_par)
    elseif menu.selection[] == "Crank Nicolson - ABC off"
        BC_par =  FDM.CN_boundary_conditions(init_params = init_params[], ABC = "false", v1=init_params[].k, v2=init_params[].k)
        CN = FDM.CN_matrix(init_params[], BC_par)
    end

    real_node = Node(real.(init_params[].ψ))
    imag_node = Node(imag.(init_params[].ψ))
    abs2_node = Node(abs2.(init_params[].ψ))

    lines!(fig[1,1], init_params[].x, real_node, color=:red,  axis = (xlabel="x", ylabel = "y",))
    lines!(fig[1,3], init_params[].x, real_node, color=:red,  axis = (xlabel="x", ylabel = "y",))
    lines!(fig[2,1], init_params[].x, imag_node, color=:blue,  axis = (xlabel="x", ylabel = "y",))
    lines!(fig[2,2], init_params[].x, abs2_node, color=:black,  axis = (xlabel="x", ylabel = "y",))


    norm_init = sum(abs2.(init_params[].ψ[1:end]))
    norm_node = Node(Point2f0[(0,1)])
    lines!(fig[2, 3], norm_node)


    record(fig, "test.mp4") do io
        for i in init_params[].t
            if menu.selection[] == "Crank Nicolson - ABC on" || menu.selection[] == "Crank Nicolson - ABC off"
                init_params[].ψ = CN*init_params[].ψ
            elseif menu.selection[] == "RK4"
                init_params[].ψ = FDM.rk4(init_params[])
            end

            current_norm = sum(abs2.(init_params[].ψ[2:end-1]))/norm_init

            update_norm = Point2f0(i, current_norm)
            norm_node[] = push!(norm_node[], update_norm)
            ylims!(ax_node, floor(current_norm, digits=20)/1.1, 1.1)
            xlims!(ax_node, 0, i+0.05)
            tmp_max = maximum([abs(minimum(real.(init_params[].ψ[2:end-1]))), maximum(real.(init_params[].ψ[2:end-1]))])
            ylims!(ax_max,-tmp_max,tmp_max)
            real_node[] = real.(init_params[].ψ)
            imag_node[] = imag.(init_params[].ψ)
            abs2_node[] = abs2.(init_params[].ψ)
            ax3.title = "t = $(round(i, digits=4))"
            recordframe!(io)
        end
    end
end

@lift begin
    if $(toggle.active) == true
        empty!(ax1), empty!(ax2), empty!(ax3), empty!(ax_max)
        @async animate_plot()
    end
end

