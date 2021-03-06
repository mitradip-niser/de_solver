# DE Solver
This code uses Runge Kutta algorithm and shooting method for solving the Differential Equation

### [ode/ode_coupled_rk.f95](https://github.com/mitradip-niser/de_solver/blob/master/ode/ode_coupled_rk.f95)
This code is used for solving coupled and **second order non-eigen** ODE's using 4th order Runge Kutta method. For the latter, change `ode2` function to look like `z1=y`.

### [ode/ode_evp_rk_ord1.f95](https://github.com/mitradip-niser/de_solver/blob/master/ode/ode_evp_rk_ord1.f95)
This code is used for solving **first order eigen value** ODE's using 4th order Runge Kutta method. However, you are required to put two estimates to get the eigen values.

### [ode/ode_evp_rk_ord2.f95](https://github.com/mitradip-niser/de_solver/blob/master/ode/ode_evp_rk_ord2.f95)
This code is used for solving **second order eigen value** ODE's using 4th order Runge Kutta method. However, you are required to put two estimates to get the eigen values. Please keep an eye on normalization while solving.

### [ode/ode_ord2_ev_det.f95](https://github.com/mitradip-niser/de_solver/blob/master/ode/ode_ord2_ev_det.f95)
This code is used to determine the possible **eigen values** of a **second order ODE**. It basically plots y<sub>exp</sub>-y<sub>calc</sub> vs different values of the eigen value parameter (**m** in this code) and the zeros of this function are the eigen values. This code is a half-finished one and can be considered as an extension of *ode/ode_evp_rk_ode2.f95*.

### [ode/ode_ord2_all_eval_range.f95](https://github.com/mitradip-niser/de_solver/blob/master/ode/ode_ord2_all_eval_range.f95)
This code is used to determine **all the eigen values and eigen vectors** in a given range for a second order ODE. There are options for normalization of the eigen vector, which is generally done in case of solving quantum mechanical systems. However, one has to be a bit careful while setting the `stepsize`,i.e. it should not exceed the least difference between two eigen values.


#### COMPILATION AND USAGE:
These codes have been complied using **gfortran 5.4.0** and tested in **Linux Mint 18.1**. It is expected that the codes will compile in the same or higher version of gfortran. Since there are quite a few interactions with files, the user is requested to download and complie the code rather than using any online compiler.
