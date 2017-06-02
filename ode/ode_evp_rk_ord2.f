      subroutine ode1(x,y,z,w,y1)
!     This is an ode of form y'=f(x,y,z)
      real::x,y,z,w,y1
      y1=z
      return
      end

      subroutine ode2(x,y,z,w,z1)
      real::x,y,z,w,z1
!     Please change this ODE according to your needs.
!     'w' is the eigen value
      z1=-1.*w*y
      return
      end

      subroutine ode3(x,y,z,w,w1)
      real::x,y,z,w,w1
      w1=0
      return
      end

      subroutine ode_solve(x0,xend,n,y0,z0,w0,xarr,yarr,zarr,warr)
!     Solve coupled ODE using RK4 method
      real::x0,xend,y0,z0,w0,xarr(100),yarr(100),zarr(100),warr(100)
      integer::n,i
!     Variable Declaration
!     --------------------
!     x0                Start point
!     xend              End Point
!     n                 No of steps to take
!     y0                y0=y(x0)
!     z0                z0=z(x0)
!     w0                w0=w(x0)
!     xarr,yarr,zarr    Array to return the values of x,y and z as array
!     Rest of the variables are for internal work
      real::h,k(4),l(4),m(4),x,y,z,f,w
      h=(xend-x0)/n
      x=x0
      y=y0
      z=z0
      w=w0

      do i=1,n
      call ode1(x,y,z,w,f)
      k(1)=h*f
      call ode2(x,y,z,w,f)
      l(1)=h*f
      call ode3(x,y,z,w,f)
      m(1)=h*f

      call ode1(x+h/2,y+k(1)/2,z+l(1)/2,w+m(1)/2,f)
      k(2)=h*f
      call ode2(x+h/2,y+k(1)/2,z+l(1)/2,w+m(1)/2,f)
      l(2)=h*f
      call ode3(x+h/2,y+k(1)/2,z+l(1)/2,w+m(1)/2,f)
      m(2)=h*f

      call ode1(x+h/2,y+k(2)/2,z+l(2)/2,w+m(2)/2,f)
      k(3)=h*f
      call ode2(x+h/2,y+k(2)/2,z+l(2)/2,w+m(2)/2,f)
      l(3)=h*f
      call ode3(x+h/2,y+k(2)/2,z+l(2)/2,w+m(2)/2,f)
      m(3)=h*f

      call ode1(x+h,y+k(3),z+l(3),w+m(3),f)
      k(4)=h*f
      call ode2(x+h,y+k(3),z+l(3),w+m(3),f)
      l(4)=h*f
      call ode3(x+h,y+k(3),z+l(3),w+m(3),f)
      m(4)=h*f

      x=x+h
      y=y+(k(1)+2*k(2)+2*k(3)+k(4))/6
      z=z+(l(1)+2*l(2)+2*l(3)+l(4))/6
      w=w+(m(1)+2*m(2)+2*m(3)+m(4))/6
      xarr(i)=x
      yarr(i)=y
      zarr(i)=z
      warr(i)=w
      enddo
      return
      end


      subroutine eigen_val(x0,y0,z0,xn,yn,m1,m2,convg,ev,xarr,yarr)
!     Subroutine for determining eigen value using shooting method
!     Variable Declaration
!     --------------------
!     x0,xn     Limits
!     y0        y0=y(x0)
!     yn        yn=y(xn)
!     m1,m2     Guess values for eigen value
!     convg     Convergence Criteria
!     ev        For sending back eigen value
!     xarr,yarr The eigen vector
!     All other variables are for internal use
      real::x0,y0,z0,xn,yn,z1,z2,ev,convg
      real::xarr(100),yarr(100),zarr(100),warr(100)
      real::x,y,z,b1,b2,b3,z3,m1,m2
      integer::counter
      counter=0
      z1=m1
      z2=m2

      call ode_solve(x0,xn,100,y0,z0,z1,xarr,yarr,zarr,warr)
      if (abs(yarr(100)-yn)<convg) then
      ev=z1
      goto 20
      end if
      b1=yarr(100)

      call ode_solve(x0,xn,100,y0,z0,z2,xarr,yarr,zarr,warr)
      if (abs(yarr(100)-yn)<convg) then
      ev=z2
      goto 20
      end if
      b2=yarr(100)

10    if (b1/=b2) then
      z3=z2+(((z2-z1)*(yn-b2))/(b2-b1))
      else
      z3=z2+((z2-z1)*(yn-b2))
      endif
      counter=counter+1

      call ode_solve(x0,xn,100,y0,z0,z3,xarr,yarr,zarr,warr)
      z1=z2
      z2=z3
      b1=b2
      b2=yarr(100)
      write(*,*)counter,b2,z2
      if (abs(b2-yn)>=convg) then
      goto 10
      end if

      ev=z2
20    return
      end

!     Main Code begins here
      program ode_evp_rk_ord1
      implicit none
      real::xarr(100),yarr(100),eval
      integer::i
      call eigen_val(0.,0.,5.,3.142,0.,5.,30.,0.0001,eval,xarr,yarr)

      open(10,file='eval2_trial.txt')
      write(10,*)'Eigen Value=',eval
      do i=1,100
      write(10,*)xarr(i),yarr(i)
      enddo
      stop
      end
