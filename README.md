# DE Solver
This code uses Runge Kutta algorithm and shooting method for solving the Differential Equation

### ode/ode_coupled_rk.f95
This code is used for solving coupled and **second order non-eigen** ODE's using 4th order Runge Kutta method. For the latter, change `ode2` function to look like `z1=y`.

### ode/ode_evp_rk_ord1.f95
This code is used for solving **first order eigen value** ODE's using 4th order Runge Kutta method. However, you are required to put two estimates to get the eigen values.

### ode/ode_evp_rk_ord2.f95
This code is used for solving **second order eigen value** ODE's using 4th order Runge Kutta method. However, you are required to put two estimates to get the eigen values. Please keep an eye on normalization while solving.

### ode/ode_ord2_ev_det.f95
This code is used to determine the possible **eigen values** of a **second order ODE**. It basically plots y<sub>exp</sub>-y<sub>calc</sub> vs different values of the eigen value parameter (**m** in this code) and the zeros of this function are the eigen values. This code is a half-finished one and can be considered as an extension of *ode/ode_evp_rk_ode2.f95*.
