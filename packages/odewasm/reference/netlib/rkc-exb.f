
c                           Example B
c
c  This is a simplification of Example 1 of B.P. Sommeijer, L.F. Shampine, 
c  and J.G. Verwer, RKC: an Explicit Solver for Parabolic PDEs that shows
c  the use of RKC on a substantial problem.  Semi-discretization of the 
c  heat equation in three space variables results in 19**3 = 6859 equations.
c  The inhomogeneous term and the boundary conditions have been specified so
c  that there is an analytical solution evaluated in sol(x,y,z,t).  The 
c  maximum error of the numerical solution at TEND is measured by comparison
c  to a reference solution computed to high accuracy, so the error reported 
c  is the error of the time integration, not the difference between the 
c  solutions of the ODEs and the PDE.  Some statistics about the integration 
c  are also displayed.
c
c  WARNING: This program expects the file exb.ref containing the reference
c           solution to be present in the same directory.
c
      integer          ndim
      parameter       (ndim=19*19*19)
      integer          info(4),idid
      double precision t,tend,rtol,atol
      double precision y(ndim), work(8+4*ndim)
      integer          neqn,i
      double precision yref(ndim),error
      integer          nfe,nsteps,naccpt,nrejct,nfesig,maxm
      common /rkcdid/  nfe,nsteps,naccpt,nrejct,nfesig,maxm
      integer          nx,ny,nz
      common    /grid/ nx,ny,nz
      external         f
c
      t = 0d0
      tend = 0.7d0
c-----------------------------------------
c  Define the mesh and the number of ODEs.
c  Define the initial values.
c-----------------------------------------      
      nx = 19
      ny = 19
      nz = 19
      neqn = nx*ny*nz
      call exact(neqn,t,y)
c--------------------------------------
c  Load the reference solution at TEND.
c-------------------------------------- 
      open(10,file='exb.ref')
      read(10,*) yref
c---------------------------------------------------------------
c  info(1) = 1 -- compute a solution at TEND only.
c  info(2) = 1 -- SPCRAD returns a bound on the spectral radius.
c  info(3) = 1 -- the Jacobian is constant.
c  info(4) = 0 -- ATOL is a scalar.
c---------------------------------------------------------------
      info(1) = 1
      info(2) = 1
      info(3) = 1
      info(4) = 0
c
      rtol = 1d-2
      atol = rtol
c
      idid = 0
      call rkc(neqn,f,y,t,tend,rtol,atol,info,work,idid)
c---------------------------------
c  Was the integration successful?
c---------------------------------
      if(idid .ne. 1) then
        write(*,*) ' Failed at t = ',t,' with idid = ',idid
        stop
      endif
c
      error = 0d0
      do 10 i = 1,neqn
          error = max(error,abs(y(i) - yref(i)))
10      continue
      write(*,'(/a,d8.1,a,f6.1,a,d8.2)') ' With rtol = atol =',rtol, 
     &  ', the maximum error at tend =',tend,' was',error
      write(*,'(a,i5,a)') ' The integration cost',nfe,
     &  ' function evaluations.'
      write(*,'(a,i4,a,i3,a)') ' There were',nsteps,' steps (',
     &  nrejct,' rejected).'
      write(*,'(a,i4/)') ' The maximum number of stages used was',
     &  maxm
      end

      subroutine exact(neqn,t,y)
      integer          neqn
      double precision t,y(neqn)
      integer          i,j,k,l
      double precision dx,dy,dz,sol
      integer          nx,ny,nz
      common    /grid/ nx,ny,nz
c
      dx = 1d0/(nx+1)
      dy = 1d0/(ny+1)
      dz = 1d0/(nz+1)
      do 30 i = 1,nx
        do 20 j = 1,ny
          do 10 k = 1,nz
            l = i + (j-1)*nx + (k-1)*nx*ny
            y(l) = sol(i*dx,j*dy,k*dz,t)
10        continue
20      continue            
30    continue
      return
      end

      double precision function sol(x,y,z,t)
      double precision x,y,z,t
      double precision arg
      arg = 5d0*(x + 2d0*y + 1.5d0*z - 0.5d0 - t)
      sol = tanh(arg)
      return
      end

      double precision function spcrad(neqn,t,y)
      integer          neqn
      double precision t,y(neqn)
      integer          nx,ny,nz
      common    /grid/ nx,ny,nz
      spcrad = 4d0*((nx+1)**2 + (ny+1)**2 + (nz+1)**2)
      return
      end

      subroutine f(neqn,t,y,dydt)
      integer          neqn
      double precision t,y(neqn),dydt(neqn)
      integer          i,j,k,l
      double precision u(0:20,0:20,0:20),dx,dy,dz,dxsq,dysq,dzsq,
     &                 arg,sh,ch,sol
      integer          nx,ny,nz
      common    /grid/ nx,ny,nz
c
      dx = 1d0/(nx+1)
      dy = 1d0/(ny+1)
      dz = 1d0/(nz+1)
      dxsq = dx*dx
      dysq = dy*dy
      dzsq = dz*dz
      do 30 i = 1,nx
        do 20 j = 1,ny
          do 10 k = 1,nz
            u(i,j,k) = y(i + (j-1)*nx + (k-1)*nx*ny)
10        continue
20      continue   
30    continue   
c
      do 50 i = 1,nx
        do 40 j = 1,ny
          u(i,j,0) = sol(i*dx,j*dy,0d0,t)
          u(i,j,nz+1) = sol(i*dx,j*dy,1d0,t)
40      continue
50    continue
c
      do 70 i = 1,nx
        do 60 k = 1,nz
          u(i,0,k) = sol(i*dx,0d0,k*dz,t)
          u(i,ny+1,k) = sol(i*dx,1d0,k*dz,t)
60      continue
70    continue
c
      do 90 j = 1,ny
        do 80 k = 1,nz
          u(0,j,k) = sol(0d0,j*dy,k*dz,t)
          u(nx+1,j,k) = sol(1d0,j*dy,k*dz,t)
80      continue          
90    continue
c
      do 120 i = 1,nx
        do 110 j = 1,ny
          do 100 k = 1,nz
            arg = 5d0*(i*dx + 2d0*j*dy + 1.5d0*k*dz - 0.5d0 - t)
            sh = sinh(arg)
            ch = cosh(arg)
            l = i + (j-1)*nx + (k-1)*nx*ny
            dydt(l) = (u(i-1,j,k) - 2d0*u(i,j,k) + u(i+1,j,k))/dxsq +
     &                (u(i,j-1,k) - 2d0*u(i,j,k) + u(i,j+1,k))/dysq +
     &                (u(i,j,k-1) - 2d0*u(i,j,k) + u(i,j,k+1))/dzsq +
     &                (-5d0*ch + 362.5d0*sh)/(ch**3)
100       continue
110     continue
120   continue
c
      return
      end
