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

      subroutine normalize(xarr,yarr,num,mode)
!     Normalization Mode
!     1     Area Normalization (As done in Quantum Mechanical System)
!     2     Maxima Normalization (Absolute Value)
      real::xarr(100),yarr(100),norm
      integer::num,mode,i
      norm=0.
      if (mode==1) then
        do i=1,num-1
          norm=norm+(((abs(yarr(i))+abs(yarr(i+1)))/2.)**2)*(xarr(i+1)-xarr(i))
        end do
        write(*,*)norm
        do i=1,num
        yarr(i)=yarr(i)/sqrt(norm)
      end do
      elseif (mode==2) then
      do i=1,num
      if (norm<=abs(yarr(i))) then
        norm=abs(yarr(i))
      end if
      end do
      do i=1,num
        yarr(i)=yarr(i)/norm
      end do
      end if
      end

      subroutine eigen_diff_val(x0,y0,z0,xn,yn,m1,m2,stepsize,convg,ev,xarr,yarr)
!     Subroutine for determining eigen value using shooting method
!     Variable Declaration
!     --------------------
!     x0,xn     Limits
!     y0        y0=y(x0)
!     yn        yn=y(xn)
!     m1,m2     Guess values for eigen value
!     stepsize  The size of step taken for next iteration
!     convg     Convergence Criteria
!     ev        For sending back eigen value
!     xarr,yarr The eigen vector
!     All other variables are for internal use
      real::x0,y0,z0,xn,yn,z1,z2,ev,convg,stepsize
      real::xarr(100),yarr(100),zarr(100),warr(100)
      real::x,y,z,b1,b2,b3,z3,m1,m2,df,dk,dg
      integer::i
      df=m1
      open(15,file='all_eval.txt')
      call ode_solve(x0,xn,100,y0,z0,df,xarr,yarr,zarr,warr)
      b1=yarr(100)-yn
20    df=df+stepsize
      call ode_solve(x0,xn,100,y0,z0,df,xarr,yarr,zarr,warr)
      b2=yarr(100)-yn
      if ((b1*b2)<=0) then
      dg=df-stepsize
30    dk=(df+dg)/2.
      call ode_solve(x0,xn,100,yn,z0,dk,xarr,yarr,zarr,warr)
      b3=yarr(100)-yn
      if (abs(b3)>=convg) then
        if (b3*b1<=0) then
          b2=b3
          df=dk
        else
          b1=b3
          dg=dk
        end if
        goto 30
      end if
      call normalize(xarr,yarr,100,1)
      write(15,*)' '
      write(15,*)'Eigenvalue=',dk
      do i=1,100
      write(15,*)xarr(i),yarr(i)
      enddo
      endif
      b1=b2
      if (df<=m2) then
      goto 20
      endif

      return
      end

!     Main Code begins here
      program ode_evp_rk_ord1
      implicit none
      real::xarr(100),yarr(100),eval
      integer::i
      call eigen_diff_val(0.,0.,1.,3.142,0.,0.,100.,1.25,0.00005,eval,xarr,yarr)
      stop
      end
