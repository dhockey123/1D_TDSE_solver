include("potential_fields.jl")
include("Solvers.jl")
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
println("If you ran this from terminal it probably won't work. Run it in VSCode with Ctrl+a then Ctrl+Enter")

Base.@kwdef mutable struct Params
    Nx::Int64
    Nt::Int64
    x_min::Float64
    x_max::Float64
    t_max::Float64
    
    dt::Float64             = t_max/(Nt-1)
    dx::Float64             = (abs(x_min)+abs(x_max))/(Nx-1)
    t::Array{Float64}       = collect(0:dt:t_max)
    x::Array{Float64}       = collect(x_min:dx:x_max)
    k_space::Array{Float64} = [ collect(0:((Nx-1)/2)) ; collect(-(Nx)/2:-1) ] *2π/(Nx*dx)

    k::Float64
    σ::Float64
    x_0::Float64

    # Initial wavefunction shape at t=0
    ψ::Array{ComplexF64}    = (((2π*(σ^2))^(-0.25)).*exp.((-(x .- x_0).^2)/(4*σ^2)+(im*k.*x)))
    ψ_rk::Array{ComplexF64} = (((2π*(σ^2))^(-0.25)).*exp.((-(x .- x_0).^2)/(4*σ^2)+(im*k.*x)))

    norm_init::Float64 = sum(abs2.(ψ[1:end]))
    
    potential_field::Vector{Float64} = zeros(Int64, Nx)
end

#################################################################################################
#################################### DRAW GRIDS #################################################
#################################################################################################

fig = Figure(backgroundcolor = RGBf0(0.98, 0.98, 0.98),resolution = (1200, 900))
ax1 = Axis(fig[1, 1], title = "ℝ(Ψ) ")
ax2 = Axis(fig[2, 1], title = "𝕀(Ψ)")
ax3 = Axis(fig[2, 2], title = "|Ψ|²")
ax_node = Axis(fig[2,3], yscale=log10, title="‖Ψ‖ ", xlabel = "Time(s)") 
ax_max = Axis(fig[1, 3], title = "ℝ(Ψ) - Dynamic limits")

################################################################################################
################################  INTERACTIVE ##################################################
################################################################################################

ylims!(ax1, -1, 1), ylims!(ax2, -1, 1), ylims!(ax3, -.2, 1.5)

s1, interval_layout = layoutscene(padding=0)
s2, U_barrier_layout = layoutscene(padding=0)

x_0_slider   = labelslider!(fig, "x_0:", -30:1:30; format = x -> "$(x)", startvalue=4)
k_slider     = labelslider!(fig, "k:", -15:0.1:15; format = x -> "$(x)")
σ_slider     = labelslider!(fig, "σ:", 0.1:0.1:10; format = x -> "$(x)")
Nx_slider    = labelslider!(fig, "Nx:", 10:2:2048; format = x -> "$(x)")
Nt_slider    = labelslider!(fig, "Nt:", 10:1:512; format = x -> "$(x)")
t_max_slider = labelslider!(fig, "t:", 0.1:0.1:10; format = x -> "$(x)")

Δt = @lift round($(t_max_slider.slider.value)/$(Nt_slider.slider.value), digits=5)
Δx = @lift round((($(x_range.interval)[2]-$(x_range.interval)[1]))/$(Nx_slider.slider.value), digits=5)
Δxx = @lift round($(Δx)^2, digits = 5)
Δt_label = Label(fig,@lift string("Δt=", $(Δt)))
Δx_label = Label(fig,@lift string("Δx=", $(Δx)))
Δx2_label = Label(fig,@lift string("Δx² ", $(Δxx)))

interval_layout[2,1] = Label(fig, "x,x")
interval_layout[2,2] = x_range = IntervalSlider(fig, range = -30:1:30, startvalues = (-5, 5), tellwidth=true)
interval_layout[2,3] = Label(fig, @lift string($(x_range.interval)))
interval_layout[1,1] = Δt_label
interval_layout[1,2] = Δx_label
interval_layout[1,3] = Δx2_label



fig[1,2] = vgrid!(interval_layout,
                x_0_slider.layout,
                k_slider.layout,
                σ_slider.layout,
                Nx_slider.layout,
                Nt_slider.layout,
                t_max_slider.layout, tellwidth=true)


solver_menu = Menu(fig[3, 1], options = ["Crank Nicolson - ABC on", "Crank Nicolson - ABC off", "RK4", "Pseudospectral_FFT"])
potentials_menu = Menu(fig[3, 3], options = ["Oscillator", "Barrier"])

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




#################################################################################################
######################################## ANIMATE PLOT ###########################################
#################################################################################################

function animate_plot()
    # [[mag_1, pos_1, width_1], [mag_2, pos_2, width_2], [mag_3, pos_3, width_3]]
    barriers = [[100, -4, 1], [100, 4, 1]]

    if potentials_menu.selection[] == "Oscillator"
        init_params[].potential_field = Potential_fields.potential_oscillator(init_params[], 3)
        lines!(fig[1,1], init_params[].x, (init_params[].x .- init_params[].x_0).^2/50, color=:black, linestyle=:dash)
    end
    if potentials_menu.selection[] == "Barrier"
        init_params[].potential_field, drawx, drawy = Potential_fields.potential_barrier(init_params[], barriers)
        lines!(fig[1,1], drawx, drawy, color=:black, linestyle=:dash)
        lines!(fig[2,1], drawx, drawy, color=:black, linestyle=:dash)
        lines!(fig[2,2], drawx, drawy, color=:black, linestyle=:dash)
        lines!(fig[1,3], drawx, drawy, color=:black, linestyle=:dash)
    end


    if solver_menu.selection[] == "Crank Nicolson - ABC on"
        BC_par =  Solver.CN_boundary_conditions(init_params = init_params[], ABC = "true", v1=init_params[].k, v2=init_params[].k)
        CN = Solver.CN_matrix(init_params[], BC_par)
    elseif solver_menu.selection[] == "Crank Nicolson - ABC off"
        BC_par =  Solver.CN_boundary_conditions(init_params = init_params[], ABC = "false", v1=init_params[].k, v2=init_params[].k)
        CN = Solver.CN_matrix(init_params[], BC_par)
    end

    real_node = Node(real.(init_params[].ψ))
    imag_node = Node(imag.(init_params[].ψ))
    abs2_node = Node(abs2.(init_params[].ψ))
    norm_node = Node(Point2f0[(0,1)])

    lines!(fig[1,1], init_params[].x, real_node, color=:red,  axis = (xlabel="x", ylabel = "y",))
    lines!(fig[1,3], init_params[].x, real_node, color=:red,  axis = (xlabel="x", ylabel = "y",))
    lines!(fig[2,1], init_params[].x, imag_node, color=:blue,  axis = (xlabel="x", ylabel = "y",))
    lines!(fig[2,2], init_params[].x, abs2_node, color=:black,  axis = (xlabel="x", ylabel = "y",))
    lines!(fig[2,3], norm_node)

    record(fig, "test.mp4") do io
        for i in init_params[].t

            if solver_menu.selection[] == "Crank Nicolson - ABC on" || solver_menu.selection[] == "Crank Nicolson - ABC off"
                init_params[].ψ = CN*init_params[].ψ
            elseif solver_menu.selection[] == "RK4"
                init_params[].ψ = Solver.rk4(init_params[])
            elseif solver_menu.selection[] == "Pseudospectral_FFT"
                init_params[].ψ = Solver.Pseudospectral_FFT(init_params[])
            end

            axis_limits = maximum([abs(minimum(real.(init_params[].ψ[2:end-1]))), maximum(real.(init_params[].ψ[2:end-1]))])
            current_norm = sum(abs2.(init_params[].ψ[2:end-1]))/init_params[].norm_init
            update_norm = Point2f0(i, current_norm)
            
            ylims!(ax_node, floor(current_norm, digits=20)/1.1, 1.1)
            xlims!(ax_node, 0, i+0.05)
            ylims!(ax_max,-axis_limits,axis_limits)

            norm_node[] = push!(norm_node[], update_norm)
            real_node[] = real.(init_params[].ψ)
            imag_node[] = imag.(init_params[].ψ)
            abs2_node[] = abs2.(init_params[].ψ)
            ax3.title = "|Ψ|²     t = $(round(i, digits=4))"
            recordframe!(io)
        end
        println("\n\nYou CAN set toggle to 'Initialise' and run as new. The less buggy option is hitting Ctrl + a and Ctrl + Enter and running it again.. \n\n")
    end
    println("\n\nYou CAN set toggle to 'Initialise' and run as new. The less buggy option is hitting Ctrl + a and Ctrl + Enter and running it again.. \n\n")
end

@lift begin
    if $(toggle.active) == true
        empty!(ax1), empty!(ax2), empty!(ax3), empty!(ax_max)
        @async animate_plot()
    end
    empty!(ax1), empty!(ax2), empty!(ax3), empty!(ax_max), empty!(ax_node)
end
