c23456-----------------------------------------------------------------2
      program gd_3
c
ccccccccccccccccccccccccccccccccccccccccc
c                                       c
c     Author: C. S. Ng                  c
c                                       c
c     Last change: 2022-01-07           c
c                                       c
ccccccccccccccccccccccccccccccccccccccccc
c
c  This program calculate a 3D BGK mode by Hankel transform using
c  explicit J_0 and J_1
c  Next order of the galactic model by reading the previous order
c     and output fields for the next order
c  Use arrays to store the Bessel functions
c
c average the output to improve convergence
c iteration from hdi to hdf in nhd steps
c
      integer amax
      PARAMETER (amax=4096)
      real bf,bd,zeta,tau,zd,hd,hdi,hdf,tol,kd,zetad,taud,zmax,s
      real ka,ka2,cp,ca,e1,bd2,zt,zt1,p0,a0,ztd,kzd1,kzd2,fb
      real x,xmax,dx,k,kmax,dk,dxk,z,dz,dz2,kz,kz2,fhd
      real fp(0:amax),fa(0:amax),ip(0:amax),ia(0:amax)
      real j0(0:amax*amax),j1(0:amax*amax)
      real psi(0:amax,0:amax),aphi(0:amax,0:amax)
      real psi0(0:amax,0:amax),aphi0(0:amax,0:amax)
      real psik(0:amax,0:amax)
      real ga(0:amax,0:amax),ppk(0:amax),apk(0:amax)
      real a(0:amax),b(0:amax),c(0:amax),f(0:amax),t(0:amax)
      real bessj0, bessj1
      real psum2,p0sum2,pdif2,asum2,a0sum2,adif2,dif
      integer nmax,n,i,j,l,ij,nhd,ihd
c
      namelist/in/bf,bd,zeta,tau,zd,hdi,hdf,nhd,tol,kd,zetad,
     .   taud,zmax,s,nmax
      read(*,in)
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
      kmax = 2.*s/ka
      dk = kmax/n
      dxk = dx*dk
      write(*,*)zmax
      dz = zmax/n
      write(*,*)dz
      dz2 = dz*dz
c
c ka is the width of the guassian function
c nmax is the maximum power of 2 such that # of grid is 2^nmax
c
      hd = hdi
      zt = zeta/tau
      zt1 = 1. + zt
      bd2 = bd*bd
      ztd = zetad/taud
      kzd1 = 2.*kd*zetad
      kzd2 = kd*(zetad*bf - 2.*kd*taud)
c
      do j = 0,n
        do i = 0,n
          j0(i*j) = 100.
        enddo
      enddo
c
      do j = 0,n
        do i = 0,n
          ij = i*j
          if (j0(ij).eq.100.) then
            j0(ij) = bessj0(ij*dxk)
            j1(ij) = bessj1(ij*dxk)
          endif
        enddo
      enddo
c
      if (nhd.lt.0) then
c
        nhd = -nhd
        cp = zd*hd*ka2/4.
        ca = -bd2*zd*hd*kd*ka2*ka2/4.
        do i = 0,n
          k = i*dk
          e1 = exp(-(k*ka/2.)**2)
          fp(i) = cp*k*e1/sqrt(k*k + zt1)
          fa(i) = ca*k*e1
        enddo
c
        do l = 0,n
          z = l*dz
          do j = 0,n
            do i = 0,n
              ij = i*j
              k = i*dk
              ip(i) = fp(i)*j0(ij)*exp(-sqrt(k*k+zt1)*z)
              ia(i) = fa(i)*j1(ij)*exp(-k*z)
            enddo
            do i=1,nmax+1
              call trapzd(ip,0.,kmax,psi0(j,l),i,n)
              call trapzd(ia,0.,kmax,aphi0(j,l),i,n)
            enddo
          enddo
        enddo
c
      else
        open(unit=10,file='gd_3-psi0.txt',status='old')
        do j = 0,n
          do l = 0,n
            read(10,*) psi0(j,l)
          enddo
        enddo
        close(unit=10)
        open(unit=20,file='gd_3-aphi0.txt',status='old')
        do j = 0,n
          do l = 0,n
            read(20,*) aphi0(j,l)
          enddo
        enddo
        close(unit=20)
      endif
c
      fhd = (hdf/hdi)**(1./nhd)
c
      do ihd = 1,nhd
        hd = hd*fhd
        cp = -zd*hd/2.
        ca = -bd2*zd*hd*kd
        write(*,*)"hd = ", hd
        do j = 0,n
          do l = 0,n
            psi0(j,l) = psi0(j,l)*fhd
            aphi0(j,l) = aphi0(j,l)*fhd
          enddo
        enddo
        dif = 1.
        do while (dif.gt.tol)
c
c - forward transform
c
          do l = 0,n
            z = l*dz
            do j = 0,n
              do i = 0,n
                ij = i*j
                x = i*dx
                p0 = psi0(i,l)
                ip(i) = x*j0(ij)*(exp(-zt*p0)-exp(p0)+zt1*p0)
              enddo
              do i=1,nmax+1
                call trapzd(ip,0.,xmax,ga(j,l),i,n)
              enddo
            enddo
          enddo
c
          do j = 0,n
            do i = 0,n
              ij = i*j
              x = i*dx
              p0 = psi0(i,0)
              a0 = aphi0(i,0)
              e1 = exp(-kzd1*x*a0-kzd2*x*x-ztd*p0)
              ip(i) = cp*x*j0(ij)*e1
              ia(i) = ca*x*x*j1(ij)*e1
            enddo
            do i=1,nmax+1
              call trapzd(ip,0.,xmax,ppk(j),i,n)
              call trapzd(ia,0.,xmax,apk(j),i,n)
            enddo
          enddo
c
          do j = 0,n
            k = j*dk
            kz2 = k*k + zt1
            kz = sqrt(kz2)
            fb = -2. - kz2*dz2
            do i = 0,n
              a(i) = 1.
              c(i) = 1.
              b(i) = fb
              f(i) = -ga(j,i)*dz2
            enddo
            c(0) = 2.
            f(0) = f(0) + 2.*dz*ppk(j)
            b(n) = b(n) + exp(-kz*dz)
            call tridia(a,b,c,f,n+1,t)
            do i = 0,n
              psik(j,i) = t(i)
            enddo
          enddo
c
c - backward transform
c
          do l = 0,n
            z = l*dz
            do j = 0,n
              do i = 0,n
                ij = i*j
                k = i*dk
                ip(i) = k*j0(ij)*psik(i,l)
                ia(i) = apk(i)*j1(ij)*exp(-k*z)
              enddo
              do i=1,nmax+1
                call trapzd(ip,0.,kmax,psi(j,l),i,n)
                call trapzd(ia,0.,kmax,aphi(j,l),i,n)
              enddo
            enddo
          enddo
c
          psum2 = 0.
          p0sum2 = 0.
          pdif2 = 0.
          asum2 = 0.
          a0sum2 = 0.
          adif2 = 0.
          do j = 0,n
            do l = 0,n
              psum2 = psum2 + psi(j,l)**2
              p0sum2 = p0sum2 + psi0(j,l)**2
              pdif2 = pdif2 + (psi(j,l)-psi0(j,l))**2
              psi0(j,l) = (psi(j,l)+psi0(j,l))/2.
              asum2 = asum2 + aphi(j,l)**2
              a0sum2 = a0sum2 + aphi0(j,l)**2
              adif2 = adif2 + (aphi(j,l)-aphi0(j,l))**2
              aphi0(j,l) = (aphi(j,l)+aphi0(j,l))/2.
            enddo
          enddo
          dif = sqrt(2.*pdif2/(psum2 + p0sum2))
     .            + sqrt(2.*adif2/(asum2 + a0sum2))
          write(*,*) dif
          call flush()
        enddo
      enddo
c
      open(unit=10,file='gd_3-psi.txt',status='new')
      do j = 0,n
        do l = 0,n
          write(10,*) psi0(j,l)
        enddo
      enddo
      close(unit=10)
      open(unit=20,file='gd_3-aphi.txt',status='new')
      do j = 0,n
        do l = 0,n
          write(20,*) aphi0(j,l)
        enddo
      enddo
      close(unit=20)
c
      write(*,*)"p0sum2 = "
      write(*,*) p0sum2
      write(*,*)"psum2 = "
      write(*,*) psum2
      write(*,*) "pdif ="
      write(*,*) sqrt(2.*pdif2/(psum2 + p0sum2))
      write(*,*)"a0sum2 = "
      write(*,*) a0sum2
      write(*,*)"asum2 = "
      write(*,*) asum2
      write(*,*) "adif ="
      write(*,*) sqrt(2.*adif2/(asum2 + a0sum2))
c
1000  continue
      STOP
      END     
c23456-----------------------------------------------------------------2
c The following are from the Numerical Recipes -- tested and modified.
c
      FUNCTION bessj0(x)
      REAL bessj0,x
      REAL ax,xx,z
      DOUBLE PRECISION p1,p2,p3,p4,p5,q1,q2,q3,q4,q5,r1,r2,r3,r4,r5,r6,
     *s1,s2,s3,s4,s5,s6,y
      SAVE p1,p2,p3,p4,p5,q1,q2,q3,q4,q5,r1,r2,r3,r4,r5,r6,s1,s2,s3,s4,
     *s5,s6
      DATA p1,p2,p3,p4,p5/1.d0,-.1098628627d-2,.2734510407d-4,
     *-.2073370639d-5,.2093887211d-6/, q1,q2,q3,q4,q5/-.1562499995d-1,
     *.1430488765d-3,-.6911147651d-5,.7621095161d-6,-.934945152d-7/
      DATA r1,r2,r3,r4,r5,r6/57568490574.d0,-13362590354.d0,
     *651619640.7d0,-11214424.18d0,77392.33017d0,-184.9052456d0/,s1,s2,
     *s3,s4,s5,s6/57568490411.d0,1029532985.d0,9494680.718d0,
     *59272.64853d0,267.8532712d0,1.d0/
      if(abs(x).lt.8.)then
        y=x**2
        bessj0=(r1+y*(r2+y*(r3+y*(r4+y*(r5+y*r6)))))/(s1+y*(s2+y*(s3+y*
     *(s4+y*(s5+y*s6)))))
      else
        ax=abs(x)
        z=8./ax
        y=z**2
        xx=ax-.785398164
        bessj0=sqrt(.636619772/ax)*(cos(xx)*(p1+y*(p2+y*(p3+y*(p4+y*
     *p5))))-z*sin(xx)*(q1+y*(q2+y*(q3+y*(q4+y*q5)))))
      endif
      return
      END
c23456-----------------------------------------------------------------2
      FUNCTION bessj1(x)
      REAL bessj1,x
      REAL ax,xx,z
      DOUBLE PRECISION p1,p2,p3,p4,p5,q1,q2,q3,q4,q5,r1,r2,r3,r4,r5,r6,
     *s1,s2,s3,s4,s5,s6,y
      SAVE p1,p2,p3,p4,p5,q1,q2,q3,q4,q5,r1,r2,r3,r4,r5,r6,s1,s2,s3,s4,
     *s5,s6
      DATA r1,r2,r3,r4,r5,r6/72362614232.d0,-7895059235.d0,
     *242396853.1d0,-2972611.439d0,15704.48260d0,-30.16036606d0/,s1,s2,
     *s3,s4,s5,s6/144725228442.d0,2300535178.d0,18583304.74d0,
     *99447.43394d0,376.9991397d0,1.d0/
      DATA p1,p2,p3,p4,p5/1.d0,.183105d-2,-.3516396496d-4,
     *.2457520174d-5,-.240337019d-6/, q1,q2,q3,q4,q5/.04687499995d0,
     *-.2002690873d-3,.8449199096d-5,-.88228987d-6,.105787412d-6/
      if(abs(x).lt.8.)then
        y=x**2
        bessj1=x*(r1+y*(r2+y*(r3+y*(r4+y*(r5+y*r6)))))/(s1+y*(s2+y*(s3+
     *y*(s4+y*(s5+y*s6)))))
      else
        ax=abs(x)
        z=8./ax
        y=z**2
        xx=ax-2.356194491
        bessj1=sqrt(.636619772/ax)*(cos(xx)*(p1+y*(p2+y*(p3+y*(p4+y*
     *p5))))-z*sin(xx)*(q1+y*(q2+y*(q3+y*(q4+y*q5)))))*sign(1.,x)
      endif
      return
      END
c23456-----------------------------------------------------------------2
      SUBROUTINE trapzd(func,a,b,s,n,nmax)
      INTEGER n,nmax
      REAL a,b,s
      real func(0:nmax)
      INTEGER it,j,dn,nx
      REAL del,sum,tnm,x
      if (n.eq.1) then
        s=0.5*(b-a)*(func(0)+func(nmax))
      else
        it=2**(n-2)
        tnm=it
        dn=nmax/tnm
        nx=dn/2
        sum=0.
        do j=1,it
          sum=sum+func(nx)
          nx=nx+dn
        enddo
        s=0.5*(s+(b-a)*sum/tnm)
      endif
      return
      END
c23456-----------------------------------------------------------------2
c tridiagonal matrix algorithm
c
c copied from Y. Jaluria K. Torrance ,computational heat transfer pp336
c
c a,b and c are the three elements in each row. With b at the diagonal
c f is the constant on the right-hand side of each equation. n is the
c number of equations and t is the variable to be computed.
c
      subroutine tridia(a,b,c,f,n,t)
c
      integer nmax
      parameter (nmax = 4097)
      integer n,i,j
c
      real a(nmax),b(nmax),c(nmax),f(nmax),t(nmax),d
c
c Elimiante the a's by Gaussian elimination and determine
c the new coefficients.
c
      do i = 2,n
        d = a(i)/b(i-1)
        b(i) = b(i) - c(i-1)*d
        f(i) = f(i) - f(i-1)*d
      enddo
c
c back subsititution
c
      t(n) = f(n)/b(n)
      do i = 1,n-1
        j = n-i
        t(j) = (f(j) - c(j)*t(j+1))/b(j)
      enddo
      return
      end
c----------------------------------------------------------------------2
