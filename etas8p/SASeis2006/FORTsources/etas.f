      program etasap
c-----------------------------------------------------------------------
c Maximum likelihood procedure for the ETAS point process model.
c Subroutine FUNC4 corresponds to the exact likelihood and FUNC9 to the
c approximate version: see one of the references below.
c 
c     References:
c     Ogata, Y. (1988). J. Amer. Statist. Soc., 83, pp. 9-27.
c       ---     (1989). Tectonophysics, 169, pp. 159-174.
c       ---     (1992). J. Geophys. Res. 97, pp. 19845-19871.
c     Ogata, Y., Matsu'ura S.R., Katsura, K. (1993). submitted to 
c                Geophys. Res. Letters.
c-----------------------------------------------------------------------
      implicit real * 8 (a-h,o-z)
      parameter(ldata=17777, npara=5)
      common/xyod/xx(ldata),xmg(ldata),xmag0
      common/param/xini(npara),n
      common /range/tstart,ntstar
      common /kkxy/kkx,kky,kkt
      common t,nn,mm,iappr,nfunct
c
      call input
c
      do 10 i=2,nn
      if(xx(i).ge.xx(i-1)) go to 10
      write(6,*) 'reverse occurrence time'
      write(6,*) i,xx(i),xx(i-1),xmg(i),xmg(i-1)
   10 continue
      write(6,8) nfunct
    8 format(1h ,' funct = ',i5)
      write(6,9) t,nn,mm
      write(6,5) xmag0
    5 format(1h ,'reference magnitudes; xmag0',5x,f10.4)
    3 format(1h ,10f12.4/(1h ,10f12.4))
    9 format(1h ,'t,nn,mm',5x,f10.4,2i6)
    2 format(f10.2,2i10)
    1 format(8f10.2)
c
      call finout
c
   20 continue
      stop
      end
c***********************************************************************
      subroutine input
c
c       Reading parameters
c
      implicit real * 8 (a-h,o-z)
      parameter(ldata=17777, npara=5)
      common/param/xini(npara),n
      dimension z(ldata),amg(ldata)
      character*80 hypodata
      common /xyod/xx(ldata),xmg(ldata),xmag0
      common /kkxy/kkx,kky,kkt
      common /fukasa/dep(ldata)
      common t,nn,mm,iappr,nfunct
      common /range/tstart,ntstar
      open(unit=1,file='etas.open')
c     open(unit=1,file='./etaspc.open')
c     open(unit=1,file='siminit.dat')
*     read(1,112) hypodata
  112 format(a)
      amx2=10.0
      read(1,*) nfunct,iappr
      read(1,*) zts,zte,tstart
      read(1,*) amx1,xmag0
      write(6,*) nfunct,iappr
      write(6,*) zts,zte,tstart
      write(6,*) amx1,xmag0
      n=5
      read(1,*) (xini(i),i=1,n)
      close(unit=1)
      write(6,*) '      minmag, dep1, ts, tend,tstart'
      write(6,5) amx1,zts,zte,tstart
    1 format(1h ,20a4)
    4 format(8f10.3)
    5 format(1h ,8f10.3)
    9 format(3i10)
    2 format(f10.0,i10)
    3 format(20a4)
c
      call zisin(t,nd,z,amg,dep,hypodata)
c
      t=zte-zts
      tstart=tstart-zts
      nnn=0
      nn=0
      ntstar=0
      smg=0.0
      do 10 i=1,nd
      if(amg(i).lt.amx1.or.amg(i).gt.amx2) go to 10
c     if(dep(i).lt.depx1.or.dep(i).gt.depx2) go to 10
      if(z(i).lt.zts.or.z(i).gt.zte) go to 10
      nn=nn+1
      write(6,1002) nn,i,amg(i),z(i),dep(i)
 1002 format(2i7,4f12.5,f5.0)
      if(z(i).lt.tstart) ntstar=nn
      xx(nn)=z(i)-zts
      xmg(nn)=amg(i)
      dep(nn)=dep(i)
      amct=amx1
      if(xmg(nn).ge.amct) smg=smg+(xmg(nn)-amct+.05)
      if(xmg(nn).ge.amct) nnn=nnn+1
   10 continue
  111 format(3(i5,f12.5,f4.1))
      mm=nd
      bvl=nnn/log(10.)/smg
      write(6,*) 'read #data; selected #data; #tstart; b-value; M_0'
      write(6,*)  nd, nn, ntstar, bvl, amct
      if(zte.ne.0.0) t=zte
      t=zte-zts
      return
      end
c***********************************************************************
      subroutine zisin(t,nd,z,amg,dep,hypodata)
c
c            Reading hypocenter data
c
      implicit real * 8 (a-h,o-z)
      parameter(ldata=17777, npara=5)
      dimension z(ldata),amg(ldata),dep(ldata),fmt(20)
      character*80 hypodata
*     open(unit=10,file=hypodata)
      open(unit=10,file='work.etas')
      read(10,1001) fmt
 1001 format(20a4)
      write(6,1001) fmt
      i=1
   10 read(10,*,end=20) idum,e,an,amg(i),z(i),dep(i)
 1002 format(i5,4f12.5,f5.0)
      i=i+1
      go to 10
   20 nd=i-1
      close(unit=10)
      return
      end
c***********************************************************************
      subroutine finout
      implicit real * 8 (a-h,o-z)
      parameter(ldata=17777,npara=5)
      external func4,func9
      common/xyod/xx(ldata),xmg(ldata),xmag0
      common t,nn,mm,iappr,nfunct
      common/param/xini(npara),n
      common /ddd/ f,aic2
      common /kkxy/kkx,kky,kkt
      dimension x(npara)
c
      do 10 i=1,nn
   10 xmg(i)=xmg(i)-xmag0
c
      do 60 i=1,n
   60 x(i)=xini(i)
      write(6,1020)   n
      write(6,1030)  (x(i),i=1,n)
c
      do 70 i=1,n
   70 x(i)=sqrt(x(i))
c
      do 30 ii=1,1 
c
      if(nfunct.eq.4) call davidn(x,n,4,func4)
c
      if(nfunct.eq.9) call davidn(x,n,9,func9)
c
   30 continue
c
      do 80 i=1,n
   80 x(i)=x(i)**2
c
      write(6,1040) f,(x(i),i=1,n)
      aic2=f+n
 1005 format(e20.10,3i10)
      write(6,1001) aic2
 1001 format(1h ,'aic/2 =',e20.10)
   20 continue
      return
 1000 format(3i10,2f15.6)
 1010 format(8f10.4)
 1020 format(1h ,'n=',i3)
 1030 format(1h ,'x=',6e12.5)
 1040 format(1h , 'neg max lklhd=',1 e16.7
     3      /' max lklhd est.=',9e12.5/('                 ',9e12.5))
 1050 format(4d20.13)
 1060 format(e25.15)
      end
c***********************************************************************
      subroutine  davidn( x,n,ihes,funct )
c
c          minimization by davidon-fletcher-powell procedure
c
c       the following subroutines are directly called by this subroutine
c             funct
c             hesian
c             linear
c          inputs:
c             x:       vector of initial values
c             k:       dimension of the vector x
c
c          output:
c             x:       vector of minimizing solution
c
      implicit  real * 8  ( a-h , o-z )
      parameter(npara=5)
      external funct
      dimension  x(npara) , dx(npara) , g(npara) , g0(npara) , y(npara)
      dimension  h(npara,npara) , wrk(npara) , s(npara)
      common /ccc/ isw,ipr
      common /ddd/ f,aic2
      data  tau1 , tau2  /  1.0d-5 , 1.0d-5  /
      data  eps1 , eps2  / 1.0d-5 , 1.0d-5  /
      ramda = 0.5d0
      const1 = 1.0d-70
c
c          initial estimate of inverse of hessian
c
cgk
      iccc=0
22221 continue
      iccc=iccc+1
cgk
      do  20   i=1,n
      do  10   j=1,n
   10 h(i,j) = 0.0d00
      s(i) = 0.0d00
      dx(i) = 0.0d00
   20 h(i,i) = 1.0d00
      isw = 0
c
      call  funct( n,x,xm,g,ig )
c
      write( 6,340 )     xm
c
      icc = 0
c      iteration
11110 continue
      icc = icc + 1
      do  11111   ic=1,n
      if( ic .eq. 1 .and. icc .eq. 1 )     go to 120
c
      do  40   i=1,n
   40 y(i) = g(i) - g0(i)
      do  60   i=1,n
      sum = 0.0d00
      do  50   j=1,n
   50 sum = sum + y(j) * h(i,j)
   60 wrk(i) = sum
      s1 = 0.0d00
      s2 = 0.0d00
      do  70   i=1,n
      s1 = s1 + wrk(i) * y(i)
   70 s2 = s2 + dx(i) * y(i)
      if( s1.le.const1 .or. s2.le.const1 )  go to 900
      if( s1 .le. s2 )     go to 100
c
c          update the inverse of hessian matrix
c
c               ---  davidon-fletcher-powell type correction  ---
c
      do  90   i=1,n
      do  90   j=i,n
      h(i,j) = h(i,j) + dx(i)*dx(j)/s2 - wrk(i)*wrk(j)/s1
   90 h(j,i) = h(i,j)
      go to  120
c
c               ---  fletcher type correction  ---
c
  100 continue
      stem = s1 / s2 + 1.0d00
      do  110   i=1,n
      do  110   j=i,n
      h(i,j) = h(i,j)- (dx(i)*wrk(j)+wrk(i)*dx(j)-dx(i)*dx(j)*stem)/s2
  110 h(j,i) = h(i,j)
c
c
c
  120 continue
      ss = 0.0d00
      do  150   i=1,n
      sum = 0.0d00
      do  140   j=1,n
  140 sum = sum + h(i,j)*g(j)
      ss = ss + sum * sum
  150 s(i) = -sum
c
c
      s1 = 0.0d00
      s2 = 0.0d00
      do  170   i=1,n
      s1 = s1 + s(i)*g(i)
  170 s2 = s2 + g(i)*g(i)
      ds2 = dsqrt(s2)
      gtem = dabs(s1) / ds2
c     write(6,610)gtem,ds2
      if( gtem .le. tau1  .and.  ds2 .le. tau2 )     go to  900
      if( s1 .lt. 0.0d00 )     go to  200
      do  190   i=1,n
      do  180   j=1,n
  180 h(i,j) = 0.0d00
      h(i,i) = 1.0d00
  190 s(i) = -s(i)
  200 continue
c
      ed = xm
c
c          linear  search
c
      call  linear( x,s,ramda,ed,n,ig,funct )
c
      write( 6,330 )     ramda , ed , s1 , s2
c
      s1 = 0.0d00
      do  210   i=1,n
      dx(i) = s(i) * ramda
      s1 = s1 + dx(i) * dx(i)
      g0(i) = g(i)
  210 x(i) = x(i) + dx(i)
      xmb = xm
      isw = 0
c
      call  funct( n,x,xm,g,ig )
c
      s2 = 0.d0
      do  220     i=1,n
  220 s2 = s2 + g(i)*g(i)
      if( dsqrt(s2) .gt. tau2 )   go to  11111
      if( xmb/xm-1.d0 .lt. eps1  .and.  dsqrt(s1) .lt. eps2 )  go to 900
11111 continue
      if( icc .ge. 5 )     go to 900
      go to 11110
  900 continue
      write( 6,600 )
      write( 6,610 )     (x(i)**2,i=1,n)
      write( 6,601 )
      write( 6,610 )     (g(i)**2,i=1,n)
      return
  330 format( 1h ,'lambda =',d15.7,5x,'-LL =',d23.15,2x,d9.2,2x,d9.2)
  340 format( 1h ,4x,'initial (-1)*Log-Likelihood =',d23.15)
  600 format( 1h ,'-----  x  -----' )
  601 format( 1h ,'***  gradient  ***' )
  610 format( 1h ,10d13.5 )
      end
c***********************************************************************
      subroutine  linear( x,h,ram,ee,k,ig,funct )
c
c     this subroutine performs the linear search along the direction spe
c     by the vector h
c----------------------------------------------------------------------
c       the following subroutine is directly called by this subroutine:
c             funct
c----------------------------------------------------------------------
c
c     inputs:
c        x:       vector of position
c        h:       search direction
c        k:       dimension of vector x
c
c     outputs:
c        ram:     optimal step width
c        e2:      minimum function value
c        ig:      error code
c
      implicit  real  *8 ( a-h,o-z )
      parameter(npara=5)
      external funct
      integer  return,sub
      dimension  x(npara) , h(npara) , x1(npara)
      dimension  g(npara)
      common /ccc/ isw,ipr
c
      isw = 1
      ipr=7
      if( ram .le. 1.0d-30 )  ram = 0.01d0
      const2 = 1.0d-60
      hnorm = 0.d0
      do 10  i=1,k
   10 hnorm = hnorm + h(i)**2
      hnorm = dsqrt( hnorm )
c
      ram2 = ram
      e1 =ee
      ram1 = 0.d0
c
      do 20  i=1,k
   20 x1(i) = x(i) + ram2*h(i)
      call  funct( k,x1,e2,g,ig )
c     if(ipr.ge.7)  write(6,2)  ram2,e2
      if(ipr.ge.7)  write(6,8)  e2,(x1(i)**2,i=1,k)
    8 format(' -ll=',d13.5,1x,5d12.5)
c
      if( ig .eq. 1 )  go to  50
      if( e2 .gt. e1 )  go to 50
   30 ram3 = ram2*2.d0
      do 40  i=1,k
   40 x1(i) = x(i) + ram3*h(i)
      call  funct( k,x1,e3,g,ig )
      if( ig.eq.1 )  go to  500
c     if( ipr.ge.7 )  write(6,3)  ram3,e3
      if(ipr.ge.7)  write(6,8)  e3,(x1(i)**2,i=1,k)
      if( e3 .gt. e2 )  go to 70
      ram1 = ram2
      ram2 = ram3
      e1 = e2
      e2 = e3
      go to 30
c
   50 ram3 = ram2
      e3 = e2
      ram2 = ram3*0.1d0
      if( ram2*hnorm .lt. const2 )  go to  400
      do 60  i=1,k
   60 x1(i) = x(i) + ram2*h(i)
      call  funct( k,x1,e2,g,ig )
c     if(ipr.ge.7)  write(6,4)  ram2,e2
      if(ipr.ge.7)  write(6,8)  e2,(x1(i)**2,i=1,k)
      if( e2.gt.e1 )  go to 50
c
   70 assign 80 to return
      go to 200
c
   80 do 90  i=1,k
   90 x1(i) = x(i) + ram*h(i)
      call  funct( k,x1,ee,g,ig )
c     if(ipr.ge.7)  write(6,5)  ram,ee
      if(ipr.ge.7)  write(6,8)  ee,(x1(i)**2,i=1,k)
c
      ifg = 0
      assign  300 to  sub
      assign 200 to sub
   95 assign 130 to return
      if( ram .gt. ram2 )  go to 110
      if( ee .ge. e2 )  go to 100
      ram3 = ram2
      ram2 = ram
      e3 =e2
      e2 =ee
      go to  sub,( 200,300 )
c
  100 ram1 = ram
      e1 = ee
      go to  sub,( 200,300 )
c
  110 if( ee .le. e2 )  go to 120
      ram3 = ram
      e3 = ee
      go to  sub,( 200,300 )
c
  120 ram1 = ram2
      ram2 = ram
      e1 = e2
      e2 = ee
      go to  sub,( 200,300 )
c
  130 do 140  i=1,k
  140 x1(i) = x(i) + ram*h(i)
      call  funct( k,x1,ee,g,ig )
c     if( ipr.ge.7 )  write(6,6)  ram,ee
      if(ipr.ge.7)  write(6,8)  ee,(x1(i)**2,i=1,k)
      assign 200 to sub
      ifg = ifg+1
      ifg = 0
      if( ifg .eq. 1 )  go to 95
c
      if( e2 .lt. ee )  ram = ram2
      return
c
c      -------  internal subroutine sub1  -------
  200 a1 = (ram3-ram2)*e1
      a2 = (ram1-ram3)*e2
      a3 = (ram2-ram1)*e3
      b2 = (a1+a2+a3)*2.d0
      b1 = a1*(ram3+ram2) + a2*(ram1+ram3) + a3*(ram2+ram1)
      if( b2 .eq. 0.d0 )  go to 210
      ram = b1 /b2
      go to return ,( 80,130 )
c
  210 ig = 1
      ram = ram2
      return
c
c      -------  internal subroutine sub2  -------
c
  300 if( ram3-ram2 .gt. ram2-ram1 )  go to 310
      ram = (ram1+ram2)*0.5d0
      go to return ,( 80,130 )
c
  310 ram = (ram2+ram3)*0.5d0
      go to return ,( 80,130 )
c ------------------------------------------------------------
c
  400 ram = 0.d0
      return
c ------------------------------------------------------------
c
  500 ram = (ram2+ram3)*0.5d0
  510 do 520  i=1,k
  520 x1(i) = x(i) + ram*h(i)
      call  funct( k,x1,e3,g,ig )
c     if( ipr.ge.7 )  write(6,7)  ram,e3
      if(ipr.ge.7)  write(6,8)  e3,(x1(i)**2,i=1,k)
      if( ig.eq.1 )  go to 540
      if( e3.gt.e2 )  go to 530
      ram1 = ram2
      ram2 = ram
      e1 = e2
      e2 = e3
      go to 500
c
  530 ram3 = ram
      go to 70
c
  540 ram = (ram2+ram)*0.5d0
      go to 510
c
c ------------------------------------------------------------
    1 format( 1h ,'lambda =',d18.10, 10x,'e1 =',d25.17 )
    2 format( 1h ,'lambda =',d18.10, 10x,'e2 =',d25.17 )
    3 format( 1h ,'lambda =',d18.10, 10x,'e3 =',d25.17 )
    4 format( 1h ,'lambda =',d18.10, 10x,'e4 =',d25.17 )
    5 format( 1h ,'lambda =',d18.10, 10x,'e5 =',d25.17 )
    6 format( 1h ,'lambda =',d18.10, 10x,'e6 =',d25.17 )
    7 format( 1h ,'lambda =',d18.10, 10x,'e7 =',d25.17 )
      end
c***********************************************************************
      subroutine func4(n,b,f,h,ifg)
c----------------------------------------------------------------------
c     likelihood function of modified oomori type
c     hawkes' point process
c     the optimization w.r.t. a5 is not yet succeeded
c     at the date of 26th oct. 1981.
c     -----this is succeeded at 23 dec.1983-----
c----------------------------------------------------------------------
      implicit real * 8 (a-h,o-z)
      parameter(ldata=17777, npara=5)
      common/xyod/xx(ldata),xmg(ldata),xmag0
      common t,nn,mm,iappr
      common /ddd/fff,aic2
      common /range/tstart,ntstar
      common /kkxy/kkx,kky,kkt
      dimension b(npara),h(npara),ggt(npara),gt(npara),at(npara)
      ifg=0
      a1=b(1)**2
      a2=b(2)**2
      a3=b(3)**2
      a4=b(4)**2
      a5=b(5)**2
c
      if(a4.gt.10.d0) go to 50
      if(a5.gt.3.d0) go to 50
c
        ff=0.0
      d1ff=0.0
      d2ff=0.0
      d3ff=0.0
      d4ff=0.0
      d5ff=0.0
      do 110 k=1,kkt
  110 at(k)=b(k+5)
      do 120 k=1,kkt
  120 ggt(k)=0.0
      rmdti1=0.0
c
      if(xx(1).lt.tstart) go to 130
c
      do 140 k=1,kkt
  140 rmdti1=rmdti1+at(k)*(xx(1)/t)**k
      ramdai=a1+rmdti1
      if(ramdai.le.0.0) go to 50
      ff=log(ramdai)
      d1ff=1.d0/ramdai
      do 150 k=1,kkt
      ggt(k)=ggt(k)+(xx(1)/t)**k/ramdai
  150 continue
c
  130 continue
c
      do 10 i=2,nn
c
      rmdi=0.0
      rmd1i=0.0
      rmdmi=0.0
      rmdli=0.0
c
      do 20 j=1,i-1
      rmdi=rmdi+exp(a4*xmg(j))/(xx(i)-xx(j)+a3)**a5
      rmd1i=rmd1i+exp(a4*xmg(j))/(xx(i)-xx(j)+a3)**(a5+1)
      rmdmi=rmdmi+xmg(j)*exp(a4*xmg(j))/(xx(i)-xx(j)+a3)**a5
      rmdli=rmdli-exp(a4*xmg(j))/(xx(i)-xx(j)+a3)**a5
     &                          *log(xx(i)-xx(j)+a3)
   20 continue
c
      if(xx(i).lt.tstart) go to 10
c
      rmdti1=0.0
      if(kkt.eq.0) go to 160
      do 170 k=1,kkt
  170 rmdti1=rmdti1+at(k)*(xx(i)/t)**k
  160 continue
c
      ramdai=a1+a2*rmdi+rmdti1
      if(ramdai.le.0.0) go to 50
      ff=ff+log(ramdai)
      d1ff=d1ff+1.d0/ramdai
      d2ff=d2ff+rmdi/ramdai
      d3ff=d3ff-a2*a5*rmd1i/ramdai
      d4ff=d4ff+a2*rmdmi/ramdai
      if(a5.ne.1.d0) d5ff=d5ff+a2*rmdli/ramdai
      do 180 k=1,kkt
      ggt(k)=ggt(k)+(xx(i)/t)**k/ramdai
  180 continue
c
   10 continue
c
      d3ft=0.0
      d4ft=0.0
      d5ft=0.0
      ft=0.0
c
      if(a5.eq.1.d0) then
c
      do 30 i=ntstar+1,nn
      ft=ft+1.d0*(log(t-xx(i)+a3)-log(a3))*exp(a4*xmg(i))
      d3ft=d3ft+((t-xx(i)+a3)**(-a5)-a3**(-a5))*exp(a4*xmg(i))
      d4ft=d4ft+(log(t-xx(i)+a3)-log(a3))*xmg(i)*exp(a4*xmg(i))
   30 continue
c
      do 40 i=1,ntstar
      ft=ft+1.d0*(log(t-xx(i)+a3)-log(a3))*exp(a4*xmg(i))
      ft=ft-1.d0*(log(tstart-xx(i)+a3)-log(a3))*exp(a4*xmg(i))
      d3ft=d3ft+((t-xx(i)+a3)**(-a5)-a3**(-a5))*exp(a4*xmg(i))
      d3ft=d3ft-((tstart-xx(i)+a3)**(-a5)-a3**(-a5))*exp(a4*xmg(i))
      d4ft=d4ft+(log(t-xx(i)+a3)-log(a3))*xmg(i)*exp(a4*xmg(i))
      d4ft=d4ft-(log(tstart-xx(i)+a3)-log(a3))*xmg(i)*exp(a4*xmg(i))
   40 continue
c
      else
c
      do 60 i=ntstar+1,nn
      ft=ft+1.d0/(1.d0-a5)*((t-xx(i)+a3)**(1.d0-a5)-a3**(1.d0-a5))
     &     *exp(a4*xmg(i))
      d3ft=d3ft+((t-xx(i)+a3)**(-a5)-a3**(-a5))*exp(a4*xmg(i))
      d4ft=d4ft+1.d0/(1.d0-a5)*((t-xx(i)+a3)**(1.d0-a5)
     &                 -a3**(1.d0-a5))*xmg(i)*exp(a4*xmg(i))
      d5ft=d5ft+1.d0/(1.d0-a5)*(-(t-xx(i)+a3)**(1.d0-a5)
     &      *log(t-xx(i)+a3)+a3**(1.d0-a5)*log(a3))*exp(a4*xmg(i))
   60 continue
c
      do 70 i=1,ntstar
      ft=ft+1.d0/(1.d0-a5)*((     t-xx(i)+a3)**(1.d0-a5)-a3**(1.d0-a5))
     &     *exp(a4*xmg(i))
      ft=ft-1.d0/(1.d0-a5)*((tstart-xx(i)+a3)**(1.d0-a5)-a3**(1.d0-a5))
     &     *exp(a4*xmg(i))
      d3ft=d3ft+((     t-xx(i)+a3)**(-a5)-a3**(-a5))*exp(a4*xmg(i))
      d3ft=d3ft-((tstart-xx(i)+a3)**(-a5)-a3**(-a5))*exp(a4*xmg(i))
      d4ft=d4ft+1.d0/(1.d0-a5)*((     t-xx(i)+a3)**(1.d0-a5)
     &                 -a3**(1.d0-a5))*xmg(i)*exp(a4*xmg(i))
      d4ft=d4ft-1.d0/(1.d0-a5)*((tstart-xx(i)+a3)**(1.d0-a5)
     &                 -a3**(1.d0-a5))*xmg(i)*exp(a4*xmg(i))
      d5ft=d5ft+1.d0/(1.d0-a5)*(-(     t-xx(i)+a3)**(1.d0-a5)
     &      *log(     t-xx(i)+a3)+a3**(1.d0-a5)*log(a3))*exp(a4*xmg(i))
      d5ft=d5ft-1.d0/(1.d0-a5)*(-(tstart-xx(i)+a3)**(1.d0-a5)
     &      *log(tstart-xx(i)+a3)+a3**(1.d0-a5)*log(a3))*exp(a4*xmg(i))
c
   70 continue
c
      endif
c
      fs=a1*(t-tstart)+a2*ft
      if(a1.gt.0.0) d1fs=(t-tstart)
      d2fs=ft
      d3fs=a2*d3ft
      d4fs=a2*d4ft
      d5fs=0.0
      if(a5.ne.1.d0) d5fs=a2*(ft/(1.d0-a5)+d5ft)
c
      stsum=0.0
      do 80 k=1,kkt
      stsum=stsum+at(k)*t*(     t/t)**(k+1)/(k+1)
      stsum=stsum-at(k)*t*(tstart/t)**(k+1)/(k+1)
   80 continue
c
      f=ff-fs-stsum
      g1=d1ff-d1fs
      g2=d2ff-d2fs
      g3=d3ff-d3fs
      g4=d4ff-d4fs
      g5=d5ff-d5fs
c
      do 90 k=1,kkt
      gt(k)=ggt(k)-t*(     t/t)**(k+1)/(k+1)
      gt(k)=ggt(k)+t*(tstart/t)**(k+1)/(k+1)
   90 continue
c
      f=-f
      h(1)=-g1*2.0d0*b(1)
      h(2)=-g2*2.0d0*b(2)
      h(3)=-g3*2.0d0*b(3)
      h(4)=-g4*2.0d0*b(4)
      h(5)=-g5*2.0d0*b(5)
      do 100 k=1,kkt
  100 h(k+5)=-gt(k)
      if(a5.eq.1.d0) h(5)=0.0
      fff=f
c     write(6,1030) fff,a1,a2,a3,a4,a5
 1030 format(1h ,'f=',e12.5,'; x=',6e12.5)
    3 format(1h ,110x,d18.10)
    1 format(1h ,7d18.10)
      return
   50 continue
      ifg=1
      f=1.0d30
      return
      end
c***********************************************************************
      subroutine func9(n,b,f,h,ifg)
c-----------------------------------------------------------------------
c     likelihood function of "etas" model with the approximation
c     method: this is copied fron func2 and modified (1992 october)
c     improved for the effective use of 'iap' processor.
c     when is=16, the gradient were not converge because approximation
c     for a3 did not fit well. this subrourine has overcomed this
c     by making use of the direct diffrerential of lamdai: see pi#(.).
c-----------------------------------------------------------------------
      implicit real * 8 (a-h,o-z)
      parameter(ldata=17777, npara=5)
      common/xyod/xx(ldata),xmg(ldata),xmag0
      common t,nn,mm,iappr
      common /ddd/fff,aic2
      common /range/tstart,ntstar
      dimension b(npara),h(npara)
      dimension xi1(144),xi2(144),wx1(144),wx2(144)
      dimension fi1(144),fi2(144),alf1(144),alf2(144),ci1(144),ci2(144)
      dimension rmd(ldata),rmdc(ldata),rmdm(ldata),rmdl(ldata)
      data ixhiab /0/
      save ixhiab,delta0,xi0,xi1,xi2,wx0,wx1,wx2
c
      if(ixhiab.eq.0) call hiab(delta0,xi0,xi1,xi2,wx0,wx1,wx2)
c
      pi2=1.570796326794397d0
c take 'is' below as one of the 1,2,4,8 or 16:
c the larger is the faster but the less accurate.
      is=iappr
      im=is
c
      delta=delta0*im
      if(ixhiab.eq.0) write(6,*) 'is=',is
      ixhiab=1
      ifg=0
      a1=b(1)**2
      a2=b(2)**2
      a3=b(3)**2
      a4=b(4)**2
      a5=b(5)**2
c
      if(a5.gt.10.d0) go to 50
c
      ff=0.0
      d1ff=0.0
      d2ff=0.0
      d3ff=0.0
      d4ff=0.0
      d5ff=0.0
c
      if(ixhiab.eq.0) call hiab(delta,xi0,xi1,xi2,wx0,wx1,wx2)
      ixhiab=1
c
c the folloings are not redundant: consider the case where i=1.
      rmdi=0.0
      rmdci=0.0
      rmdmi=0.0
      rmdli=0.0
      fi0=0.0
      ci0=0.0
      alf0=0.0
      do 90 j=is,144,im
      fi1(j)=0.0
      fi2(j)=0.0
      ci1(j)=0.0
      ci2(j)=0.0
      alf1(j)=0.0
      alf2(j)=0.0
   90 continue
      cpg=a3**(-a5)/gam(0,a5)
      cpg3=-a5*a3**(-a5-1)/gam(0,a5)
      cpg5=-cpg*(log(a3)+gam(1,a5)/gam(0,a5))
      do 10 i=1,nn
      if(xx(i).eq.0.0.and.i.eq.1) go to 10
      if(i.eq.1) go to 10
      if (a4*xmg(i-1).gt.170) go to 50
      if (xx(i)-xx(i-1).lt.0.0) write(6,*) 'reverse',i,xx(i),xmg(i)
c
      ci0=exp(-xi0/a3*(xx(i)-xx(i-1)))*
     &    (ci0+xi0/a3**2*(xx(i)-xx(i-1))*(fi0+exp(a4*xmg(i-1))))
      qi0=ci0*xi0**(a5-1)*exp(-xi0)
      fi0=(fi0+exp(a4*xmg(i-1)))*exp(-xi0/a3*(xx(i)-xx(i-1)))
      gi0=fi0*xi0**(a5-1)*exp(-xi0)
      hi0=gi0*log(xi0)
      alf0=(alf0+xmg(i-1)*exp(a4*xmg(i-1)))*exp(-xi0/a3*(xx(i)-xx(i-1)))
      blf0=alf0*xi0**(a5-1)*exp(-xi0)
      rmd(i)=gi0*wx0
      rmdc(i)=qi0*wx0
      rmdm(i)=blf0*wx0
      rmdl(i)=hi0*wx0
   10 continue
c
      do 20 j=is,144,im
      do 25 i=2,nn
c
      ci1(j)=exp(-xi1(j)/a3*(xx(i)-xx(i-1)))*
     &   (ci1(j)+xi1(j)/a3**2*(xx(i)-xx(i-1))*(fi1(j)+exp(a4*xmg(i-1))))
      qi1=ci1(j)*xi1(j)**(a5-1)*exp(-xi1(j))
      ci2(j)=exp(-xi2(j)/a3*(xx(i)-xx(i-1)))*
     &   (ci2(j)+xi2(j)/a3**2*(xx(i)-xx(i-1))*(fi2(j)+exp(a4*xmg(i-1))))
      qi2=ci2(j)*xi2(j)**(a5-1)*exp(-xi2(j))
      rmdc(i)=rmdc(i)+qi1*wx1(j)+qi2*wx2(j)
      fi1(j)=(fi1(j)+exp(a4*xmg(i-1)))*exp(-xi1(j)/a3*(xx(i)-xx(i-1)))
      gi1=fi1(j)*xi1(j)**(a5-1)*exp(-xi1(j))
      hi1=gi1*log(xi1(j))
      fi2(j)=(fi2(j)+exp(a4*xmg(i-1)))*exp(-xi2(j)/a3*(xx(i)-xx(i-1)))
      gi2=fi2(j)*xi2(j)**(a5-1)*exp(-xi2(j))
      hi2=gi2*log(xi2(j))
      rmd(i)=rmd(i)+gi1*wx1(j)+gi2*wx2(j)
      rmdl(i)=rmdl(i)+hi1*wx1(j)+hi2*wx2(j)
c
      alf1(j)=(alf1(j)+xmg(i-1)*exp(a4*xmg(i-1))) *
     &      exp(-xi1(j)/a3*(xx(i)-xx(i-1)))
      blf1=alf1(j)*xi1(j)**(a5-1)*exp(-xi1(j))
      alf2(j)=(alf2(j)+xmg(i-1)*exp(a4*xmg(i-1))) *
     &      exp(-xi2(j)/a3*(xx(i)-xx(i-1)))
      blf2=alf2(j)*xi2(j)**(a5-1)*exp(-xi2(j))
      rmdm(i)=rmdm(i)+blf1*wx1(j)+blf2*wx2(j)
   25 continue
   20 continue
c
c     do 15 i=1,nn
      do 15 i=ntstar+1,nn
      if(i.eq.1) go to 80
      rmdi=rmd(i)*pi2*delta
      rmdci=rmdc(i)*pi2*delta
      rmdmi=rmdm(i)*pi2*delta
      rmdli=rmdl(i)*pi2*delta
   80 continue
      ramdai=a1+a2*cpg*rmdi
      ff=ff+log(ramdai)
      d1ff=d1ff+1.d0 /ramdai
      d2ff=d2ff+1.d0 /ramdai* cpg *rmdi
      d3ff=d3ff+a2   /ramdai*(cpg3*rmdi+cpg*rmdci)
      d4ff=d4ff+a2   /ramdai* cpg *rmdmi
      if(a5.ne.1.d0)
     &d5ff=d5ff+a2   /ramdai*(cpg5*rmdi+cpg*rmdli)
   15 continue
c
      d3ft=0.0
      d4ft=0.0
      d5ft=0.0
      ft=0.0
c
      if(a5.eq.1.d0) then
c
      do 30 i=ntstar+1,nn
      effmag=exp(a4*xmg(i))
      ft=ft+1.d0*(log(t-xx(i)+a3)-log(a3))*effmag
      d3ft=d3ft+((t-xx(i)+a3)**(-a5)-a3**(-a5))*effmag
      d4ft=d4ft+(log(t-xx(i)+a3)-log(a3))*xmg(i)*effmag
   30 continue
c
      do 40 i=1,ntstar
      effmag=exp(a4*xmg(i))
      ft=ft+1.d0*(log(t-xx(i)+a3)-log(a3))*effmag
      ft=ft-1.d0*(log(tstart-xx(i)+a3)-log(a3))*effmag
      d3ft=d3ft+((t-xx(i)+a3)**(-a5)-a3**(-a5))*effmag
      d3ft=d3ft-((tstart-xx(i)+a3)**(-a5)-a3**(-a5))*effmag
      d4ft=d4ft+(log(t-xx(i)+a3)-log(a3))*xmg(i)*effmag
      d4ft=d4ft-(log(tstart-xx(i)+a3)-log(a3))*xmg(i)*effmag
   40 continue
c
      else
c
      do 60 i=ntstar+1,nn
      effmag=exp(a4*xmg(i))
      ft=ft+1.d0/(1.d0-a5)*((t-xx(i)+a3)**(1.d0-a5)-a3**(1.d0-a5))
     &     *effmag
      d3ft=d3ft+((t-xx(i)+a3)**(-a5)-a3**(-a5))*effmag
      d4ft=d4ft+1.d0/(1.d0-a5)*((t-xx(i)+a3)**(1.d0-a5)
     &                 -a3**(1.d0-a5))*xmg(i)*effmag
      d5ft=d5ft+1.d0/(1.d0-a5)**2 *
     &           ((t-xx(i)+a3)**(1.d0-a5)-a3**(1.d0-a5)) * effmag
     &         +1.d0/(1.d0-a5)    * ( -(t-xx(i)+a3)**(1.d0-a5)
     &      *log(t-xx(i)+a3)+a3**(1.d0-a5)*log(a3) ) * effmag
   60 continue
c
      do 70 i=1,ntstar
      effmag=exp(a4*xmg(i))
      ft=ft+1.d0/(1.d0-a5)*((t-xx(i)+a3)**(1.d0-a5)-a3**(1.d0-a5))
     &     *effmag
      ft=ft-1.d0/(1.d0-a5)*((tstart-xx(i)+a3)**(1.d0-a5)-a3**(1.d0-a5))
     &     *effmag
      d3ft=d3ft+((t-xx(i)+a3)**(-a5)-a3**(-a5))*effmag
      d3ft=d3ft-((tstart-xx(i)+a3)**(-a5)-a3**(-a5))*effmag
      d4ft=d4ft+1.d0/(1.d0-a5)*((t-xx(i)+a3)**(1.d0-a5)
     &                 -a3**(1.d0-a5))*xmg(i)*effmag
      d4ft=d4ft-1.d0/(1.d0-a5)*((tstart-xx(i)+a3)**(1.d0-a5)
     &                 -a3**(1.d0-a5))*xmg(i)*effmag
      d5ft=d5ft+1.d0/(1.d0-a5)**2 *
     &           ((t-xx(i)+a3)**(1.d0-a5)-a3**(1.d0-a5)) * effmag
     &         +1.d0/(1.d0-a5)    * ( -(t-xx(i)+a3)**(1.d0-a5)
     &      *log(t-xx(i)+a3)+a3**(1.d0-a5)*log(a3) ) * effmag
      d5ft=d5ft-1.d0/(1.d0-a5)**2 *
     &           ((tstart-xx(i)+a3)**(1.d0-a5)-a3**(1.d0-a5)) * effmag
     &         -1.d0/(1.d0-a5)    * ( -(tstart-xx(i)+a3)**(1.d0-a5)
     &      *log(tstart-xx(i)+a3)+a3**(1.d0-a5)*log(a3) ) * effmag
   70 continue
c
      endif
c
      fs=a1*(t-tstart)+a2*ft
      d1fs=0.0
      if(a1.gt.0.0) d1fs=t-tstart
      d2fs=ft
      d3fs=a2*d3ft
      d4fs=a2*d4ft
      d5fs=0.0
      if(a5.ne.1.d0) d5fs=a2*d5ft
c
      f=ff-fs
      g1=d1ff-d1fs
      g2=d2ff-d2fs
      g3=d3ff-d3fs
      g4=d4ff-d4fs
      g5=d5ff-d5fs
c
      f=-f
      h(1)=-g1*2.0d0*b(1)
      h(2)=-g2*2.0d0*b(2)
      h(3)=-g3*2.0d0*b(3)
      h(4)=-g4*2.0d0*b(4)
      h(5)=-g5*2.0d0*b(5)
      if(a5.eq.1.d0) h(5)=0.0
      fff=f
c     write(6,1030) fff,a1,a2,a3,a4,a5
 1030 format(1h ,'f=',e12.5,'; x=',6e12.5)
    3 format(1h ,110x,d18.10)
    1 format(1h ,7d13.6)
      return
   50 continue
      ifg=1
      f=1.0d30
      return
      end
c***********************************************************************
      subroutine hiab(h,a0,a1,a2,b0,b1,b2)
      implicit real*8(a-h,o-z)
      dimension a1(144),a2(144),b1(144),b2(144)
      pi2=1.570796326794397d0
      h=1.d0/32
      eh=exp(h)
      ehi=1/eh
      eni=1.d0
      in=144
      inc=in
      a0=1.d0/exp(pi2)
      b0=2.d0*a0
      do 20 i=1,in
      eni=ehi*eni
      en=1.d0/eni
      s1=h*i
      a1(i)=1.d0/exp(pi2*(s1+en))
      b1(i)=a1(i)*(1.d0+en)
      a2(i)=exp(pi2*(s1-eni))
      b2(i)=a2(i)*(1.d0+eni)
   20 continue
      return
      end
c***********************************************************************
c     real*8 function gam(id,q)
      real*8 function dbgam(id,q)
c hitac    real function dbgam*8(id,q)
      implicit real*8(a-h, o-z)
      dimension a(10)
      data a/ 0.99999 99998 71452d0, 0.42278 43615 29813d0,
     1        0.41183 94326 49605d0, 0.08158 87915 49927d0,
     2        0.07416 87114 09713d0, 0.00004 89152 06125d0,
     3        0.01038 02945 70428d0, -.00162 85524 78086d0,
     4        0.00082 46883 39196d0, -.00000 66427 76723d0 /
      fact=1
      dfac=0.0
      d2fac=0.0
      p=q
   30 if(p.ge.1.d0 .and. p.le.2.d0) then
        d2fac=d2fac*p+2*dfac
        dfac=dfac*p+fact
        fact=fact*p
        x=p-1
      else if (p.lt.1) then
        d2fac=d2fac*p+2*dfac
        dfac=dfac*p+fact
        fact=fact*p
        p=p+1
        go to 30
      else
        d2fac=(d2fac*(p-1))/(p-1)**2
     &        -2*(dfac*(p-1)-fact)/(p-1)**3
        dfac=(dfac*(p-1)-fact)/(p-1)**2
        dfac=(dfac*(p-1)-fact)/(p-1)**2
        fact=fact/(p-1)
        p=p-1
        go to 30
      end if
c
      gamm=0.0
      gam1=0.0
      dgam=0.0
      d2gam=0.0
      eps=0.1d-6
      do 10 i=1,10
      gamm=gamm+a(i)*x**(i-1)
      dgam=dgam+(i-1)*a(i)*x**(i-2)
      d2gam=d2gam+(i-1)*(i-2)*a(i)*x**(i-3)
   10 continue
      if(id.eq.0) gam=gamm/fact
      if(id.eq.1) gam=(dgam*fact-gamm*dfac)/fact**2
      if(id.eq.2) gam=(d2gam*fact-gamm*d2fac)/fact**2
     &             -2*(dgam*fact-gamm*dfac)/fact**3 *dfac
      dbgam=gam
      return
      end
c***********************************************************************
c     real*8 function qbgam(id,qq)
      real*8 function gam(id,qq)
c hitac     real function gam*8(id,q)
c     implicit real*16(a-h, o-z)
      implicit real* 8(a-h, o-z)
      real*8 qq
      dimension a(11),b(11)
      data a/ -2 98354.32785 74342 13883 04376 59         d0,
     1        -2 38495.39700 18198 87246 87344 23         d0,
     2        -1 17049.47601 21780 68840 38544 45         d0,
     3        -  39494.45048 30157 19364 21824 091        d0,
     4        -  10466.99423 82752 14053 30650 531        d0,
     5        -   2188.21811 00718 16369 39479 5998       d0,
     6        -    380.51122 08641 73465 75849 22631      d0,
     7        -     52.83123 75563 58453 83718 97838 2    d0,
     8        -      6.12857 17637 04498 30688 94282 12   d0,
     9        -       .50280 18054 41681 24673 64198 75   d0,
     t        -       .03343 06032 23305 95274 51566 0112 d0 /
      data b/ -2 98354.32785 74342 13883 04385 24         d0,
     1        -1 12355.86087 48644 91134 23064 08         d0,
     2           53327.16689 11814 21574 85686 311        d0,
     3            8571.16049 89070 43851 96114 7763       d0,
     4        -   4734.86597 70282 11706 55681 977        d0,
     5             196.04976 12885 58583 89970 39621      d0,
     6             125.77333 67869 88864 59666 47426      d0,
     7        -     20.53126 15310 06727 64513 92906 7    d0,
     8               1.                                   d0,
     9               0.0                                  d0,
     t               0.0                                  d0 /
      fact=1
      dfac=0.0
      d2fac=0.0
      p=qq
   30 if(p.ge.1.d0 .and. p.le.2.d0) then
        d2fac=d2fac*p+2*dfac
        dfac=dfac*p+fact
        fact=fact*p
        x=p-1
      else if (p.lt.1.d0) then
        d2fac=d2fac*p+2*dfac
        dfac=dfac*p+fact
        fact=fact*p
        p=p+1
        go to 30
      else
        d2fac=d2fac/(p-1)+2*dfac/(p-1)**2*(-1)
     &        +fact/(p-1)**3*(-1)*(-2)
        dfac=dfac/(p-1)+fact/(p-1)**2*(-1)
        fact=fact/(p-1)
        p=p-1
        go to 30
      end if
c
      gam1=a(1)
      gam2=b(1)
      do 10 i=2,11
      gam1=gam1+a(i)*x**(i-1)
   10 gam2=gam2+b(i)*x**(i-1)
      dga1=a(2)
      dga2=b(2)
      do 40 i=3,11
      dga1=dga1+(i-1)*a(i)*x**(i-2)
      dga2=dga2+(i-1)*b(i)*x**(i-2)
   40 continue
      d2ga1=2*a(3)
      d2ga2=2*b(3)
      do 50 i=4,11
      d2ga1=d2ga1+(i-1)*(i-2)*a(i)*x**(i-3)
      d2ga2=d2ga2+(i-1)*(i-2)*b(i)*x**(i-3)
   50 continue
c
      if(id.eq.0) gam=gam1/gam2/fact
      if(id.eq.1) gam=(dga1*gam2*fact-gam1*dga2*fact-gam1*gam2*dfac)
     & / (gam2*fact)**2
      if(id.eq.2) gam=(d2ga1*gam2*fact+dga1*dga2*fact+dga1*gam2*dfac
     &                -dga1*dga2*fact-gam1*d2ga2*fact-gam1*dga2*dfac
     &                -dga1*gam2*dfac-gam1*dga2*dfac-gam1*gam2*d2fac)
     & / (gam2*fact)**2
     &             -2*(dga1*gam2*fact-gam1*dga2*fact-gam1*gam2*dfac)
     & / (gam2*fact)**3 * (dga2*fact+gam2*dfac)
      eps=0.1d-4
      return
      end
