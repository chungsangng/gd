c23456-----------------------------------------------------------------2
      program gd_3_B
c
ccccccccccccccccccccccccccccccccccccccccc
c                                       c
c     Author: C. S. Ng                  c
c                                       c
c     Last change: 2022-01-11           c
c                                       c
ccccccccccccccccccccccccccccccccccccccccc
c
c  This program calculate a 3D BGK mode by Hankel transform using
c  explicit J_0 with a guassian function
c  Next order of the galactic model by reading the previous order
c     and output fields for the next order
c  Use arrays to store the Bessel functions
c
c  calculate B and B field lines from solution from gd_3
c
      integer amax
      PARAMETER (amax=4096)
      real bf,bd,zeta,tau,zd,hd,hdi,hdf,tol,kd,zetad,taud,zmax,s
      real ka,ka2,cp,ca,e1,bd2,zt,zt1,p0,a0,ztd,kzd1,kzd2,fb
      real x,xmax,dx,dx2,d2x,k,kmax,dk,dxk,z,dz,dz2,d2z,kz,kz2
      real fp(0:amax),fa(0:amax),ip(0:amax),ia(0:amax)
      real psi(0:amax,0:amax),aphi(0:amax,0:amax)
      real gp(0:amax,0:amax)
      real psik(0:amax,0:amax)
      real ga(0:amax,0:amax),ppk,apk
      real a(0:amax),b(0:amax),c(0:amax),f(0:amax),t(0:amax)
      real dpdr,d2pdr,dpdz,d2pdz,prhs,dadr,d2adr,dadz,d2adz,gm
      real psum2,p0sum2,pdif2,asum2,a0sum2,adif2,p0m,a0m
c
      integer nmax,n,i,j,l,ij,n1,ip1,im1,jp1,jm1,nhd,ihd
c
      namelist/in/bf,bd,zeta,tau,zd,hdi,hdf,nhd,tol,kd,zetad,
     .   taud,zmax,s,nmax
      read(*,in)
c
      write(*,in)
      write(*,*)bf
      write(*,*)bd
      write(*,*)zeta
      write(*,*)tau
      write(*,*)zd
      write(*,*)hdi
      write(*,*)hdf
      write(*,*)nhd
      write(*,*)tol
      write(*,*)kd
      write(*,*)zetad
      write(*,*)taud
      write(*,*)zmax
      write(*,*)s
      write(*,*)nmax
c
      n = 2**nmax
      n1 = n - 1
      write(*,*)n
      ka2 = 1./kd/(zetad*bf - 2.*kd*taud)
c
      if (ka2.le.0.) then
        write(*,*)"ka2 <= 0", ka2
        STOP
      endif
c
      ka = sqrt(ka2)
      xmax = s*ka
      dx = xmax/n
      write(*,*)xmax
      write(*,*)dx
      dx2 = dx*dx
      d2x = 2.*dx
      kmax = 2.*s/ka
      dk = kmax/n
      dxk = dx*dk
      write(*,*)zmax
      dz = zmax/n
      write(*,*)dz
      dz2 = dz*dz
      d2z = 2.*dz
c
c ka is the width of the guassian function
c nmax is the maximum power of 2 such that # of grid is 2^nmax
c
      cp = -zd*hd/2.
      zt = zeta/tau
      zt1 = 1. + zt
      bd2 = bd*bd
      ca = bd2*zd*hd*kd
      ztd = zetad/taud
      kzd1 = 2.*kd*zetad
      kzd2 = kd*(zetad*bf - 2.*kd*taud)
c
      open(unit=10,file='gd_3-psi.txt',status='old')
      do j = 0,n
        do l = 0,n
          read(10,*) psi(j,l)
        enddo
      enddo
      close(unit=10)
      open(unit=20,file='gd_3-aphi.txt',status='old')
      do j = 0,n
        x = j*dx
        do l = 0,n
          read(20,*) aphi(j,l)
          gp(j,l) = (x*(aphi(j,l) + x*bf/2.))**0.5
c         gp(j,l) = x*(aphi(j,l) + x*bf/2.)
        enddo
      enddo
      close(unit=20)
c
      do j = 1,n1
        jp1 = j + 1
        jm1 = j - 1
        do i = 1,n1
          x = i*dx
          ip1 = i + 1
          im1 = i - 1
          dadr = (aphi(ip1,j)-aphi(im1,j))/d2x
          dadz = (aphi(i,jp1)-aphi(i,jm1))/d2z
          ga(i,j) = sqrt((bf + dadr + aphi(i,j)/x)**2 + dadz**2)
        enddo
      enddo
c
      j = n
      jm1 = n1
      do i = 1,n1
        x = i*dx
        ip1 = i + 1
        im1 = i - 1
        dadr = (aphi(ip1,j)-aphi(im1,j))/d2x
        dadz = (3.*aphi(i,j)-4.*aphi(i,jm1)+aphi(i,jm1-1))/d2z
        ga(i,j) = sqrt((bf + dadr + aphi(i,j)/x)**2 + dadz**2)
      enddo
c
      i = n
      im1 = n1
      x = n*dx
      do j = 1,n1
        jp1 = j + 1
        jm1 = j - 1
        dadz = (aphi(i,jp1)-aphi(i,jm1))/d2z
        dadr = (3.*aphi(i,j)-4.*aphi(im1,j)+aphi(im1-1,j))/d2x
        ga(i,j) = sqrt((bf + dadr + aphi(i,j)/x)**2 + dadz**2)
      enddo
c
      i = 0
      ip1 = 1
      x = 0.
      do j = 1,n1
        jp1 = j + 1
        jm1 = j - 1
        dadz = (aphi(i,jp1)-aphi(i,jm1))/d2z
        dadr = (-aphi(ip1+1,j)+4.*aphi(ip1,j))/d2x
        ga(i,j) = sqrt((bf + dadr)**2 + dadz**2)
      enddo
c
      j = 0
      jp1 = 1
      do i = 1,n1
        x = i*dx
        ip1 = i + 1
        im1 = i - 1
        p0 = psi(i,0)
        a0 = aphi(i,0)
        e1 = exp(-kzd1*x*a0-kzd2*x*x-ztd*p0)
        ppk = cp*e1
        apk = ca*x*e1
        dadz = -(3.*aphi(i,j)-4.*aphi(i,jp1)+aphi(i,jp1+1))/d2z
        dadr = (aphi(ip1,j)-aphi(im1,j))/d2x
        ga(i,j) = sqrt((bf + dadr + aphi(i,j)/x)**2 + apk**2)
      enddo
c
      i = n
      im1 = n1
      x = n*dx
      j = n
      jm1 = n1
      dadz = (3.*aphi(i,j)-4.*aphi(i,jm1)+aphi(i,jm1-1))/d2z
      dadr = (3.*aphi(i,j)-4.*aphi(im1,j)+aphi(im1-1,j))/d2x
      ga(i,j) = sqrt((bf + dadr + aphi(i,j)/x)**2 + dadz**2)
c
      i = 0
      ip1 = 1
      j = n
      jm1 = n1
      dadz = (3.*aphi(i,j)-4.*aphi(i,jm1)+aphi(i,jm1-1))/d2z
      dadr = (-aphi(ip1+1,j)+4.*aphi(ip1,j))/d2x
      ga(i,j) = sqrt((bf + dadr)**2 + dadz**2)
c
      i = n
      im1 = n1
      x = n*dx
      j = 0
      jp1 = 1
      p0 = psi(i,0)
      a0 = aphi(i,0)
      e1 = exp(-kzd1*x*a0-kzd2*x*x-ztd*p0)
      ppk = cp*e1
      apk = ca*x*e1
      dadz = -(3.*aphi(i,j)-4.*aphi(i,jp1)+aphi(i,jp1+1))/d2z
      dadr = (3.*aphi(i,j)-4.*aphi(im1,j)+aphi(im1-1,j))/d2x
      ga(i,j) = sqrt((bf + dadr)**2 + apk**2)
c
      i = 0
      ip1 = 1
      x = 0.
      j = 0
      jp1 = 1
      p0 = psi(i,0)
      a0 = aphi(i,0)
      e1 = exp(-kzd1*x*a0-kzd2*x*x-ztd*p0)
      ppk = cp*e1
      apk = ca*x*e1
      dadr = (-aphi(ip1+1,j)+4.*aphi(ip1,j))/d2x
      dadz = -(3.*aphi(i,j)-4.*aphi(i,jp1)+aphi(i,jp1+1))/d2z
      ga(i,j) = sqrt((bf + dadr)**2 + apk**2)
c
      open(unit=10,file='gd_3-BL.txt',status='new')
      do j = 0,n
        do l = 0,n
          write(10,*) gp(j,l)
        enddo
      enddo
      close(unit=10)
      open(unit=20,file='gd_3-BS.txt',status='new')
      do j = 0,n
        do l = 0,n
          write(20,*) ga(j,l)
        enddo
      enddo
      close(unit=20)
c
1000  continue
      STOP
      END     
c23456-----------------------------------------------------------------2
