      subroutine ode1(x,y,z,y1)
!     This is an ode of form y'=f(x,y,z)
      real::x,y,z,y1
      y1=z
      return
      end

      subroutine ode2(x,y,z,z1)
!     This is an ode of form z'=g(x,y,z)
      real::x,y,z,z1
      z1=y
      return
      end

      subroutine ode_solve(x0,xend,n,y0,z0,xarr,yarr,zarr)
!     Solve coupled ODE using RK4 method
      real::x0,xend,y0,z0,xarr(100),yarr(100),zarr(100)
      integer::n,i
!     Variable Declaration
!     --------------------
!     x0                Start point
!     xend              End Point
!     n                 No of steps to take
!     y0                y0=y(x0)
!     z0                z0=z(x0)
!     xarr,yarr,zarr    Array to return the values of x,y and z as array
!     Rest of the variables are for internal work
      real::h,k(4),l(4),x,y,z,f
      h=(xend-x0)/n
      x=x0
      y=y0
      z=z0

      do i=1,n
      call ode1(x,y,z,f)
      k(1)=h*f
      call ode2(x,y,z,f)
      l(1)=h*f
      call ode1(x+h/2,y+k(1)/2,z+l(1)/2,f)
      k(2)=h*f
      call ode2(x+h/2,y+k(1)/2,z+l(1)/2,f)
      l(2)=h*f
      call ode1(x+h/2,y+k(2)/2,z+l(2)/2,f)
      k(3)=h*f
      call ode2(x+h/2,y+k(2)/2,z+l(2)/2,f)
      l(3)=h*f
      call ode1(x+h,y+k(3),z+l(3),f)
      k(4)=h*f
      call ode2(x+h,y+k(3),z+l(3),f)
      l(4)=h*f
      x=x+h
      y=y+(k(1)+2*k(2)+2*k(3)+k(4))/6
      z=z+(l(1)+2*l(2)+2*l(3)+l(4))/6
      xarr(i)=x
      yarr(i)=y
      zarr(i)=z
      enddo
      return
      end


!     Main Code begins here
      program ode_coupled_rk
      implicit none
      real::xarr(100),yarr(100),zarr(100)
      integer::i
      call ode_solve(0.,5.,100,1.,1.,xarr,yarr,zarr)

      open(10,file='c_od.txt')
      do i=1,100
      write(10,*)xarr(i),yarr(i),zarr(i)
      enddo
      stop
      end
