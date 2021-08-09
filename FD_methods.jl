module FDM
export CN_matrix, CN_boundary_conditions, rk4

    const m, ħ = 1, 1
    using LinearAlgebra
    #=
    This defines the parameters for the absorbing boundary conditions (ABC) and, based on how the code is currently structured, is
    required to generate the Crank Nicolson propagation matrix.

    For implementing the one-way boundary conditions the user defined parameters are; ABC, a1, a2.
    e.g BC = CN_boundary_conditions(ABC="true", a1=3, a2=2)

    ABC = "true"/"false" and determines whether or not the one-way boundaries will be implemented into the CN matrix
    v1, v2 = parameters used to define the absorbed group velocities of the wave packet(s). Since ħ, m = 1, 1; 
    For best results v1, v2 = k, k. Where k is the user defined wavenumber in the "Params" struct. See sections 2.3 and 3.3.

    The construction of the absorbing boundary points, the Φ's , are shown in Section 3.3
    =#

    Base.@kwdef mutable struct CN_boundary_conditions
        init_params             = init_params
        ABC::String     = "false"
        v1::Float64     = 0
        v2::Float64     = 0
        c1::Float64     = 2/(abs(v1)+abs(v2))
        c2::Float64     = (abs(v1)*abs(v2))/(abs(v1)+abs(v2))
        Φ_1::ComplexF64 =  im/(2*init_params.dx) + (im*c1)/(2*init_params.dt) + c2/4
        Φ_2::ComplexF64 = -im/(2*init_params.dx) + (im*c1)/(2*init_params.dt) + c2/4
        Φ_3::ComplexF64 = -im/(2*init_params.dx) + (im*c1)/(2*init_params.dt) - c2/4
        Φ_4::ComplexF64 =  im/(2*init_params.dx) + (im*c1)/(2*init_params.dt) - c2/4
    end


    #=
    Generates the 2nd order accurate in space difference matrix (of size Nx * Nx) for the kinetic part of the hamiltonian operator. See T from Eq 3.3.
    Even though this is for the Crank Nicolson matrix this is done slightly differently than as shown in section 3.2, but ultimately 
    the results are the same, see Eq 2.31 for proof. 
    =#
    function diff_matrix(init_params)
        A = zeros(ComplexF64, init_params.Nx, init_params.Nx)
        A[diagind(A, 0)] .= -2
        A[diagind(A, -1)] .= A[diagind(A, 1)] .= 1
        return A
    end

    # Generates the difference matrix for the potential part of the hamiltonian operator. See V from Eq 3.3.
    function voltage_diff_matrix(init_params)
        A = zeros(ComplexF64, init_params.Nx, init_params.Nx)
        A[diagind(A, 0)] .= init_params.potential_field
        return A
    end

    #= Generates the Crank Nicolson matrix from the forward time central space and backward time central space matricepar. 
    The BTCS and FTCS matrices are analogous to the A and B matrices, respectively, as shown in Section 3.2

    The one-way boundaries are implemented by editing the the top left and bottom right corners of the FTCS and BTCS matrices, see Section 3.3
    =#

    function CN_matrix(init_params, BC_par)
        #Combines kinetic and potential operators to produce hamiltonian operator
        T = (-ħ^2)/(2*m*init_params.dx^2).*diff_matrix(init_params)
        H = T .+ voltage_diff_matrix(init_params)
        #Produces the A(BTCS) and B(FTCS) matrices shown in Eq 3.10
        #For some clarity in the logic used in constructing the FTCS and BTCS see Eq 2.31 !!!!
        FTCS = Matrix(I, init_params.Nx, init_params.Nx) - (im*init_params.dt/2ħ).*H
        BTCS = Matrix(I, init_params.Nx, init_params.Nx) + (im*init_params.dt/2ħ).*H

        # Implements the one-way / Absorbing boundary conditions.
        if BC_par.ABC == "true"
            BTCS[1,1] = BTCS[end,end]   = BC_par.Φ_1
            BTCS[1,2] = BTCS[end,end-1] = BC_par.Φ_2
            FTCS[1,1] = FTCS[end,end]   = BC_par.Φ_3
            FTCS[1,2] = FTCS[end,end-1] = BC_par.Φ_4
        end
        # Returns CN propagation matrix. CN*ψ_rk^(n+1) = (BTCS)^(-1) * (FTCS) * ψ_rk^(n)
        return FTCS / BTCS
    end

    function rk4(init_params)
        A = (im*ħ)/(2*m*init_params.dx^2).*diff_matrix(init_params) .- (im/ħ)*voltage_diff_matrix(init_params)
        k1 = A*init_params.ψ_rk
        k2 = A*(init_params.ψ_rk + k1*init_params.dt/2)
        k3 = A*(init_params.ψ_rk + k2*init_params.dt/2)
        k4 = A*(init_params.ψ_rk + k3*init_params.dt)
        init_params.ψ_rk = init_params.ψ_rk + (init_params.dt/6)*(k1 + 2*k2 + 2*k3 + k4)
    end

end

#Clean RK4 and CN constants A.