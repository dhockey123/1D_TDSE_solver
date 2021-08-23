
1D TDSE Solver
===============================
# About

A one-dimensional time-dependent Schrödinger equation visualiser that illustrates the dynamics of a gaussian wavepacket in the presence of different potential fields and also examines the effectiveness of the one-way boundaries devised by *Fevens et. al*

 Numerical solutions to the TDSE can be obtained by using the Crank Nicolson, Pseudospectral FFT or Explicit Runge Kutte method. The choice of one-way boundary conditions is also given when using the CN method.

# Usage

Using Ctrl+a and Ctrl+Enter in VSCode from "main.jl" is the recommended way to run the code. Does NOT work directly from terminal ( This might just be my system though ). 

User interaction is enabled for specifying the following,

### Domain parameters
```
Nx    : Number of grid points
Nt    : Number of time steps
t_max : Simulated runtime
x, x  : Boundary points of domain
```

### Wave parameters
```
σ   : Initial spatial spread 
k   : Initial wavenumber
x_0 : Centered position of wave packet 
```

### Numerical solver
```
Crank Nicolson - ABC on  (One-way/Absorbing boundaries)
Crank Nicolson - ABC off (Reflective boundaries)
Pseudospectral FFT
4th order Runge Kutte
```
### Potential Field 
```
Potential well
Barriers
No choice = None
```

## Constructing potential barriers
In-app construction of potential barriers is currently not available.
Barriers are constructed by defining the potential magnitude, position barrier is centered at and width of barrier in the array called "barriers"
Can construct as many barriers as needed, just ensure the barriers don't exist outside solution domain.
### example
```
# barriers = [[mag_1, pos_1, width_1], [mag_2, pos_2, width_2]]
barriers = [[100, 0, 2], [-33, 4, 4]]
```
## Changing initial conditions after animating plot

To plot again with different initial conditions wait till t reachs t_max, toggle back to "Initialise" and set up as before. The less buggy option is to rerun again from VSCode ( its much faster the 2nd time) .

# Examples

![gitgif_CN](https://user-images.githubusercontent.com/81137805/130281259-af0bbe3f-78b3-4c0f-add0-437719d5d141.gif)
![gitgif_oscillator](https://user-images.githubusercontent.com/81137805/130282258-8f24555d-58f2-4c70-8c7a-67183b6564b0.gif)


# Version Info
```
Julia   v1.6.2
GLMakie v0.4.4
FFTW    v1.4.3
```


# Task list
- [] Implement higher-order finite difference methods
- [] Produce error tests and verify numerical methods
- [] Enable in-app construction of potential barriers
- [] Implement Chebyshev solutin method 

# References for one-way boundaries

```
Fevens T, Jiang H. Absorbing Boundary Conditions for the Schrödinger Equation. SIAM Journal on Scientific Computing. 1999;21(1):255-282.

Kuska J. Absorbing boundary conditions for the Schrödinger equation on finite intervals. Physical Review B. 1992;46(8):5000-5003.

Shibata T. Absorbing boundary conditions for the finite-difference time-domain calculation of the one-dimensional Schrödinger equation. Physical Review B. 1991;43(8):6760-6763.
```







