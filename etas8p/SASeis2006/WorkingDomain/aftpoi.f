c83-11-11-15:06:19/87-03-25-10:21 nhpossn fort
      program momori
c-----------------------------------------------------------------------
c     this program performs the maximum likelihood calculation and
c     estimate of the parameters of modified omori formula assuming
c     piecewise non-stationary poisson process.
c
c     this is made based on the 'nhpossn4'.
c
c     func5---------exponential decay poisson process
c     func6---------the modified omori's model for a interval (t0,t1)
c     func9---------model for detecting the 2nd aftershock
c     func10--------model for a pair of 2nd aftershock sequences
c
c     this program is desined by y.ogata and programed by y.ogata and
c     k.katsura.  see ogata (1983; j.p.e., vol.31,pp.115-124).
c-----------------------------------------------------------------------
      implicit real * 8 (a-h,o-z)
      real*4 time
      common/xyod/xdumy,xx(19999)
      common/y/y(50)
      common t,nn,nfunct
      common/range/t0,t1,t2,t3
      dimension xxxx(50),x(50),ahaic(30)
      dimension axxx(19999),amag(19999)
      xdumy=0.0
      call input(nnnn,axxx,amag)
      kaisu=1
      iend=0
c     if(nfunct.eq.6) kaisu=100
      do 9753 ijkl=1,kaisu
      if(ijkl.ne.1) write(6,9751)
 9751 format(1h )
      if(nn.eq.0) go to 9753
      call repara(ijkl,nnn,xxxx,x)
      call dav(nnn,x,ahaic)
      call output(nnn,x,ahaic)
      if(nfunct.ne.5) call sizes(nnn,x)
c     call clock(time)
c      write(6,9752) time
 9752 format(1h ,'time= ',f10.3)
 9753 continue
 9754 continue
      stop
      end
      subroutine zisin(t,nd,z,amg,dep,xp,yp,hypodata)
c
c            reading hypocenter data
c
      implicit real * 8 (a-h,o-z)
c     parameter(ldata=7999, npara=5)
      dimension z(19999),amg(19999),dep(19999)
      dimension xp(19999),yp(19999)
      character*60 hypodata,fmt
*     write(6,112) hypodata
  112 format(1h ,a)
c     open(unit=10,file=hypodata)
      open(unit=10,file='work.etas')
      read(10,1001) fmt
 1001 format(a)
      i=1
*  10 read(10,1002,end=20) e,an,amg(i),z(i),dep(i)
   10 read(10,*,end=20) inum,xp(i),yp(i),amg(i),z(i),dep(i)
 1002 format(5x,4f12.5,f5.0,5x)
      i=i+1
      go to 10
   20 nd=i-1
      close(unit=10)
      return
      end
      subroutine input(nnnn,xxxx,amag)
c-----------------------------------------------------------------------
c     read card format
card 1: nfunct (i10); function selection
c             =5 exponential decay model
c             =6 simplest modified omori formula
c             =9 case where a secondary aftershock sequence exists.
c            =10 case where two secondary aftershock sequences exist.
card 2: fmt (20a4);   format of cards no.4 and 5
card 3: t,nn (f10.2,i10); whole time span (t>or=t1), # of events
card 4: (xx(k),k=1,nn); occurred times of events, in fmt of card 2
card 5: amag(k),k=1,nn; corresponding magnitudes (not necessary)
c *** case where nfunct is not = 5 ***
card 6: t0,t1,t2,bmag,t3 (5f10.3)
c       (t0,t1) = observed time span,
c       t2 = starting time of a secondary aftershocks
c       t3 = starting time of another secondary aftershocks
c       bmag = cut off magnitude
c
c you have to add another two input cards (see subroutine repara).
c
c-----------------------------------------------------------------------
      implicit real * 8 (a-h,o-z)
      common/xyod/xdumy,xx(19999)
      common/y/y(50)
      common t,nn,nfunct
      common/range/t0,t1,t2,t3
      dimension fmt(20)
      dimension xxxx(19999)
      dimension amag(19999),ag(19999)
      dimension dep(19999),xp(19999),yp(19999)
      character*60 hypodata
      open(unit=1,file='aftpoi.open')
*     read(1,112) hypodata
  112 format(a)
      xdumy=0.0
      read(1,*) nfunct
    5 format(20a4)
      write(6,6) nfunct
    6 format(1h ,' funct = ',i5)
      read(1,*) zts,zte,tstart
      read(1,*) amx1,xmag0
      call zisin(t,nd,xxxx,amag,dep,xp,yp,hypodata)
      zts=tstart
      t=zte-zts
ckk      tstart=tstart-zts
      tstart=tstart
      nnn=0
      nn=0
      ntstar=0
      do 10 i=1,nd
      if(amag(i).lt.amx1) go to 10
      if(xxxx(i).lt.zts.or.xxxx(i).gt.zte) go to 10
      nn=nn+1
      if(xxxx(i).lt.tstart) ntstar=nn
      xxxx(nn)=xxxx(i)
      amag(nn)=amag(i)
      dep(nn)=dep(i)
   10 continue
      do 1122 k=1,nn
 1122 xx(k)=xxxx(k)
      if(nfunct.lt.6) go to 155
      t0=zts
      t1=zte
      bmag=amx1
  155 continue
      write(6,9) t,nn
    7 format(1h ,' t0 = ',f10.4,5x,' t1 = ',f10.4,' t2 = ',f10.4,
     &               ' bmag = ',f10.4,' t3 =',f10.4)
      write(6,3) (xx(k),k=1,nn)
    3 format(/1h ,'input data'/1h ,5f12.4/(1h ,5f12.4))
    9 format(1h ,' t =',f10.3,' nn =',i10)
    2 format(f10.0,2i10)
    1 format(8f10.3)
    4 format(2i10)
 9753 continue
 9754 continue
      return
      end
      subroutine repara(ijkl,nnn,xxxx,x)
c-----------------------------------------------------------------------
card 7 nnn (i10) # of parameters
card 8 (x(i),i=1,n) (8f10.2) initial estimates
c-----------------------------------------------------------------------
      implicit real * 8 (a-h,o-z)
      common /range/t0,t1,t2,t3
      common t,nn,nfunct
      dimension xxxx(50),x(50),xini(50)
      n=5
      read(1,*) (xini(i),i=1,n)
 1010 format(7f10.4)
      nnn=n-1
      do 41 i=1,nnn
      xxxx(i)=xini(i)
      if(i.eq.nnn) xxxx(i)=xini(n)
   41 x(i)=xxxx(i)
      do 45 i=1,nnn
      if(nfunct.eq.5.or.nfunct.eq.6) x(i)=sqrt(x(i))
      if(nfunct.eq.9.and.x(i).ne.0.0) x(i)=log(x(i))
      if(nfunct.eq.10.and.x(i).ne.0.0) x(i)=log(x(i))
   45 continue
      return
      end
c
c
c
c
c
      subroutine dav(n,x,ahaic)
      implicit real * 8 (a-h,o-z)
      external func5,func6,func9,func10
      common t,nn,nfunct
      common/range/t0,t1,t2,t3
      common/y/y(50)
      common/ddd/f,aic
      dimension x(50),ahaic(30)
      write(6,1020) n
      write(6,1030)  (x(i),i=1,n)
      do 30 ii=1,2
      if(nfunct.eq.5) call davidn(x,n,0,func5)
      if(nfunct.eq.6) call davidn(x,n,0,func6)
      if(nfunct.eq.9) call davidn(x,n,0,func9)
      if(nfunct.eq.10) call davidn(x,n,0,func10)
   30 continue
   80 continue
      ahaic(1)=aic
      return
 1020 format(1h ,3x,'input data'/1h ,5x,'n=',i3)
 1030 format(1h ,                   5x,'x=',6e16.7)
      end
      subroutine output(n,x,ahaic)
      implicit real * 8 (a-h,o-z)
      common t,nn,nfunct
      common/range/t0,t1,t2,t3
      common/ddd/f,aic
      common/y/y(50)
      dimension x(50),ahaic(30)
      do 70 i=1,n
      if(nfunct.ne.9.and.nfunct.ne.10) x(i)=x(i)**2
      if(x(i).eq.0.0) go to 70
      if(nfunct.eq.9.or.nfunct.eq.10) x(i)=exp(x(i))
   70 continue
ckk      write(6,1040) f,(x(i),i=1,n)
      x0=0.0
      write(6,1040) f,(x(i),i=1,3),x0,x(4)
      ncount=1
      do 110 iii=1,ncount
      write(6,1080) iii,ahaic(iii)
  110 continue
 1080 format(1h ,i10,d20.10)
      return
 1000 format(3i10,2f15.6)
 1010 format(8f10.4)
 1020 format(1h ,3x,'input data'/1h ,5x,'n=',i3)
 1030 format(1h ,                   5x,'x=',6e16.7)
 1040 format(
     2      /1h ,'neg max lklhd=',1 e16.7
     3    /1h ,'max lklhd est.=',10e12.5/('                 ',10e12.5))
 1050 format(4d20.13)
 1100 format(i10)
 1060 format(e25.15)
 1070 format(1h ,'  c = ',e25.15)
      end
      subroutine  linear( x,h,ram,ee,k,ig,funct )
c
c     this subroutine performs the linear search along the direction spe
c     by the vector h
c       ----------------------------------------------------------------
c       the following subroutine is directly called by this subroutine:
c             funct
c       ----------------------------------------------------------------
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
      integer  return,sub
      dimension  x(1) , h(1) , x1(50)
      dimension  g(50)
      external funct
      common     / ccc /  isw , ipr
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
      if(ipr.ge.7)  write(6,8)  ram2,(x1(i)**2,i=1,k)
    8 format(1h ,'-ll=',d13.5,1x,4d12.5)
c
      if( ig .eq. 1 )  go to  50
      if( e2 .gt. e1 )  go to 50
   30 ram3 = ram2*2.d0
      do 40  i=1,k
   40 x1(i) = x(i) + ram3*h(i)
      call  funct( k,x1,e3,g,ig )
      if( ig.eq.1 )  go to  500
c     if( ipr.ge.7 )  write(6,3)  ram3,e3
      if(ipr.ge.7)  write(6,8)  ram3,(x1(i)**2,i=1,k)
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
      if(ipr.ge.7)  write(6,8)  ram2,(x1(i)**2,i=1,k)
      if( e2.gt.e1 )  go to 50
c
   70 assign 80 to return
      go to 200
c
   80 do 90  i=1,k
   90 x1(i) = x(i) + ram*h(i)
      call  funct( k,x1,ee,g,ig )
c     if(ipr.ge.7)  write(6,5)  ram,ee
      if(ipr.ge.7)  write(6,8)  ram,(x1(i)**2,i=1,k)
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
      if(ipr.ge.7)  write(6,8)  ram,(x1(i)**2,i=1,k)
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
      if(ipr.ge.7)  write(6,8)  ram,(x1(i)**2,i=1,k)
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
      e n d
      subroutine  davidn( x,n,ihes,funct )
c
c          minimization by davidon-fletcher-powell procedure
c
c       ----------------------------------------------------------------
c       the following subroutines are directly called by this subroutine
c             funct
c             hesian
c             linear
c       ----------------------------------------------------------------
c          inputs:
c             x:       vector of initial values
c             k:       dimension of the vector x
c             ihes:    =0   inverse of hessian matrix is not available
c                      =1   inverse of hessian matrix is available
c
c          output:
c             x:       vector of minimizing solution
c
      implicit  real * 8  ( a-h , o-z )
      dimension  x(50) , dx(50) , g(50) , g0(50) , y(50)
      dimension  h(50,50) , wrk(50) , s(50)
      dimension  ht(50,50),hf(50,50)
      external funct
      common     / ccc /  isw , ipr
      common     / ddd /   f , aic
      common t,nn,nfunct
      data  tau1 , tau2  /  1.0d-5 , 1.0d-5  /
      data  eps1 , eps2  / 1.0d-5 , 1.0d-5  /
      ramda = 0.5d0
      const1 = 1.0d-70
c
c          initial estimate of inverse of hessian
c
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
c          inverse of hessian computation (if available)
c
c     if( ihes .eq. 1 )   call  hesian( x,n,h )
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
      write( 6,330 )     ramda , f , s1 , s2
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
      write( 6,610 )     (x(i),i=1,n)
      write( 6,601 )
      write( 6,610 )     (g(i),i=1,n)
c     the followings are inserted by y.o.(13/10/82)
      write(6,602)
      do 651 i=1,n
  651 write(6,610) (h(i,j),j=1,n)
      if(nfunct.ne.6) return
      call fisher(x,n,hf)
      write(6,605)
  605 format(/1h ,'***  fisher matrix  ***')
      do 654 i=1,n
  654 write(6,604) (hf(i,j),j=1,n)
      call invdet(hf,hdet,n,50)
      write(6,606)
  606 format(/1h ,'***  inverse fisher  ***')
      do 655 i=1,n
  655 write(6,604) (hf(i,j),j=1,n)
  602 format(/1h ,'***  estimated inverse hessian  ***')
  604 format(1h ,8d15.5/(1h ,8d15.5))
  603 format(/1h ,'***  inverse hessian  ***')
c-----------------------------------------------------------------------
      return
  330 format( 1h ,'lambda =',d15.7,5x,'-LL =',d23.15,2x,d9.2,2x,d9.2)
  340 format( 1h ,4x,'initial (-1)*Log-Likelihood =',d23.15)
  600 format(/1h ,'-----  x  -----' )
  601 format(/1h ,'***  gradient  ***' )
  610 format( 1h ,10d13.5 )
      end
c
c
c
      subroutine invdet(x,xdet,mm,mj)
c-----------------------------------------------------------------------
c       the inverse and determinant of x computation
c----------------------------------------------------------------l 11890
c       inputs:
c          x:     mm*mm square matrix
c          mm:    dimension of x
c          mj:    absolute dimension of x in the main program
c
c       outputs:
c          x:     inverse of x
c          xdet:  determinant of x
c
      implicit  real * 8  ( a-h , o-z )
      dimension x(mj,mj)
      dimension  ids(100)
      xdet = 1.0d00
      do 10 l=1,mm
c     pivoting at l-th stage
      xmaxp=0.10000d-10
      maxi=0
      do 110 i=l,mm
    1 if( abs(xmaxp) .ge. abs(x(i,l)) )     go to 110
      xmaxp=x(i,l)
      maxi=i
  110 continue
      ids(l)=maxi
      if(maxi.eq.l) go to 120
      if(maxi.gt.0) go to 121
      xdet = 0.0d00
      go to 140
c     row interchange
  121 do 14 j=1,mm
      xc=x(maxi,j)
      x(maxi,j)=x(l,j)
   14 x(l,j)=xc
      xdet=-xdet
c 120 xdet=xdet*xmaxp
  120 continue
      xc = 1.0d00 / xmaxp
      x(l,l)=1.0d00
      do 11 j=1,mm
   11 x(l,j)=x(l,j)*xc
      do 12 i=1,mm
      if(i.eq.l) go to 12
      xc=x(i,l)
      x(i,l) = 0.0d00
      do 13 j=1,mm
   13 x(i,j)=x(i,j)-xc*x(l,j)
   12 continue
   10 continue
      if(mm.gt.1) go to 123
      go to 140
c     column interchange
  123 mm1=mm-1
      do 130 j=1,mm1
      mmj=mm-j
      jj=ids(mmj)
      if(jj.eq.mmj) go to 130
      do 131 i=1,mm
      xc=x(i,jj)
      x(i,jj)=x(i,mmj)
  131 x(i,mmj)=xc
  130 continue
  140 return
      end
c
c
      subroutine func5(n,b,f,h,ifg)
c-----------------------------------------------------------------------
c     likelihood function of exp-decay poisson process
c     lammbda = a1 + a2*exp(-a3*t)
c-----------------------------------------------------------------------
      implicit real * 8 (a-h,o-z)
      common/xyod/xdumy,xx(19999)
      common/y/y(50)
      common t,nn,nfunct
      common/ddd/ff,aic
      dimension b(50),h(50),g(50)
      ifg=0
      a1=b(1)**2
      a2=b(2)**2
      a3=b(3)**2
      f1=0.0
      gg1=0.0
      gg2=0.0
      gg3=0.0
      do 20 i=1,nn
      ramdai=a1+a2*exp(-a3*xx(i))
      if(ramdai.le.0.0) go to 50
      f1=f1+log(ramdai)
      uni=1.0d00
      gg1=gg1+uni/ramdai
      gg2=gg2+exp(-a3*xx(i))/ramdai
      gg3=gg3-a2*xx(i)*exp(-a3*xx(i))/ramdai
   20 continue
c
c
      sasump=(1.0d0-exp(-a3*t))/a3
      gs1=t
      gs2=sasump
      gs3=-sasump*a2/a3+a2/a3*t*exp(-a3*t)
      go to 240
c
c
   50 continue
      ifg=1
      f=1.0d30
      return
c
c
  240 continue
      f=f1-a1*t-a2*sasump
      g(1)=gg1-gs1
      g(2)=gg2-gs2
      g(3)=gg3-gs3
      f=-f
      h(1)=-g(1)
      h(1)=h(1)*2.0d00*b(1)
      h(2)=-g(2)
      h(2)=h(2)*2.0d00*b(2)
      h(3)=-g(3)
      h(3)=h(3)*2.0d00*b(3)
      ff=f
      na=0
      do 800 i=1,n
  800 if(b(i).ne.0.0) na=na+1
      aic=ff+na
    3 format(1h ,110x,d18.10)
    1 format(1h ,7d18.10)
      return
      end
c
c
      subroutine func6(n,b,f,h,ifg)
c-----------------------------------------------------------------------
c     likelihood function of the omori's poisson process
c     lammbda = a1 + a2/(a3+t)**a4
c-----------------------------------------------------------------------
      implicit real * 8 (a-h,o-z)
      common/xyod/xdumy,xx(19999)
      common/ddd/ff,aic
      common/y/y(50)
      common/range/t0,t1,t2,t3
      common/grad/g(50)
      common t,nn,nfunct
      dimension b(50),h(50)
      ifg=0
      a1=b(1)**2
      a2=b(2)**2
      a3=b(3)**2
      a4=b(4)**2
      if(a4.gt.5.0d00) go to 119
      if(a3.gt.10000.) go to 119
      f1=0.0
      gg1=0.0
      gg2=0.0
      gg3=0.0
      gg4=0.0
      do 20 i=1,nn
      ramdai=a1+a2/(a3+xx(i))**a4
      if(ramdai.le.0.0) go to 50
      f1=f1+log(ramdai)
      uni=1.0d00
      gg1=gg1+uni/ramdai
      gg2=gg2+uni/ramdai/(a3+xx(i))**a4
      gg3=gg3-a2*a4/ramdai/(a3+xx(i))**(a4+uni)
      gg4=gg4-a2*log(a3+xx(i))/ramdai/(a3+xx(i))**a4
   20 continue
c
c
      if(b(4).eq.uni) sasump=log(t1+a3)-log(t0+a3)
      if(b(4).gt.uni) sasump=
     & (uni/(t1+a3)**(a4-uni)-uni/(t0+a3)**(a4-uni))/(uni-a4)
      if(b(4).lt.uni) sasump=
     & ((t1+a3)**(uni-a4)-(t0+a3)**(uni-a4))/(uni-a4)
      gs1=t1-t0
      gs2=sasump
      if(b(4).ne.uni) gs3=a2*(uni/(t1+a3)**a4-uni/(t0+a3)**a4)
      if(b(4).eq.uni) gs3=a2*(uni/(t1+a3)-uni/(t0+a3))
      if(b(4).gt.uni) gs4=
     & a2/(uni-a4)**2
     &  *(uni/(t1+a3)**(a4-uni)-uni/(t0+a3)**(a4-uni))
     & +a2/(uni-a4)
     &  *(-log(t1+a3)/(t1+a3)**(a4-uni)+log(t0+a3)/(t0+a3)**(a4-uni))
      if(b(4).lt.uni) gs4=
     & a2/(uni-a4)**2
     &  *((t1+a3)**(uni-a4)-(t0+a3)**(uni-a4))
     & +a2/(uni-a4)
     &  *(-(t1+a3)**(uni-a4)*log(t1+a3)+(t0+a3)**(uni-a4)*log(t0+a3))
      go to 240
c
c
   50 continue
      ifg=1
      f=1.0d30
      return
c
c
  240 continue
      f=f1-a1*(t1-t0)-a2*sasump
      g(1)=gg1-gs1
      g(2)=gg2-gs2
      g(3)=gg3-gs3
      g(4)=gg4-gs4
      if(b(4).eq.1.0d00) g(4)=0.0
      f=-f
      h(1)=-g(1)
      h(1)=h(1)*2.0d00*b(1)
      h(2)=-g(2)
      h(2)=h(2)*2.0d00*b(2)
      h(3)=-g(3)
      h(3)=h(3)*2.0d00*b(3)
      h(4)=-g(4)
      h(4)=h(4)*2.0d00*b(4)
      ff=f
      na=0
      do 800 i=1,n
  800 if(b(i).ne.0.0) na=na+1
      aic=ff+na
    3 format(1h ,110x,d18.10)
    1 format(1h ,7d18.10)
      return
  119 continue
      f=1.0d50
      ifg=1
      return
      end
c
      subroutine func9(n,b,f,h,ifg)
c-----------------------------------------------------------------------
c     likelihood function of the omori's poisson process
c     lammbda = a1 + a2/(a3+t)**a4 + a5/(a6+t-t2)**a7 * i(t.gt.t2)
c     a model for detecting the second aftershock
c     this is not succeeded at the moment (82/11/22)
c-----------------------------------------------------------------------
      implicit real * 8 (a-h,o-z)
      common/xyod/xdumy,xx(19999)
      common/ddd/ff,aic
      common/y/y(50)
      common/range/t0,t1,t2,t3
      common/grad/g(50)
      common t,nn,nfunct
      dimension b(50),h(50)
      ifg=0
      a1=0.0
      a1=b(1)**2
      if(b(1).ne.0.0) a1=exp(b(1))
      a2=0.0
      a2=b(2)**2
      if(b(2).ne.0.0) a2=exp(b(2))
      a3=0.0
      a3=b(3)**2
      if(b(3).ne.0.0) a3=exp(b(3))
      a4=0.0
      a4=b(4)**2
      if(b(4).ne.0.0) a4=exp(b(4))
      a5=0.0
      a5=b(5)**2
      if(b(5).ne.0.0) a5=exp(b(5))
      a6=0.0
      a6=b(6)**2
      if(b(6).ne.0.0) a6=exp(b(6))
      a7=0.0
      a7=b(7)**2
      if(b(7).ne.0.0) a7=exp(b(7))
      if(b(7).eq.0.0) a7=exp(b(4))
      if(b(6).eq.0.0) a6=exp(b(3))
c     if(a4.gt.5.0d00) go to 119
      if(a3.gt.10000.) go to 119
      uni=1.0d00
      f1=0.0
      gg1=0.0
      gg2=0.0
      gg3=0.0
      gg4=0.0
      gg5=0.0
      gg6=0.0
      gg7=0.0
      if(a7*log(a6+t1-t2).gt.150.) go to 50
      if(a4*log(a3).lt.-150.) go to 50
      if(a4*log(a3+t1).gt.150.) go to 50
      do 20 i=1,nn
      ramdai=a1+a2/(a3+xx(i))**a4
      if(xx(i).gt.t2)
     &ramdai=a1+a2/(a3+xx(i))**a4+a5/(a6+xx(i)-t2)**a7
      if(ramdai.le.0.0) go to 50
      gg1=gg1+uni/ramdai
      gg2=gg2+uni/ramdai/(a3+xx(i))**a4
      gg3=gg3-a2*a4/ramdai/(a3+xx(i))**(a4+uni)
      gg4=gg4-a2*log(a3+xx(i))/ramdai/(a3+xx(i))**a4
c
      if(xx(i).le.t2) go to 10
      if(ramdai.le.0.0) go to 50
      gg5=gg5+uni/ramdai/(a6+xx(i)-t2)**a7
      gg6=gg6-a5*a7/ramdai/(a6+xx(i)-t2)**(a7+uni)
      gg7=gg7-a5*log(a6+xx(i)-t2)/ramdai/(a6+xx(i)-t2)**a7
c
   10 continue
      f1=f1+log(ramdai)
   20 continue
c
c
      if(b(4).eq.uni) sasump=log(t1+a3)-log(t0+a3)
      if(b(4).gt.uni) sasump=
     & (uni/(t1+a3)**(a4-uni)-uni/(t0+a3)**(a4-uni))/(uni-a4)
      if(b(4).lt.uni) sasump=
     & ((t1+a3)**(uni-a4)-(t0+a3)**(uni-a4))/(uni-a4)
      gs1=t1-t0
      gs2=sasump
      if(b(4).ne.uni) gs3=a2*(uni/(t1+a3)**a4-uni/(t0+a3)**a4)
      if(b(4).eq.uni) gs3=a2*(uni/(t1+a3)-uni/(t0+a3))
      if(b(4).gt.uni) gs4=
     & a2/(uni-a4)**2
     &  *(uni/(t1+a3)**(a4-uni)-uni/(t0+a3)**(a4-uni))
     & +a2/(uni-a4)
     &  *(-log(t1+a3)/(t1+a3)**(a4-uni)+log(t0+a3)/(t0+a3)**(a4-uni))
      if(b(4).lt.uni) gs4=
     & a2/(uni-a4)**2
     &  *((t1+a3)**(uni-a4)-(t0+a3)**(uni-a4))
     & +a2/(uni-a4)
     &  *(-(t1+a3)**(uni-a4)*log(t1+a3)+(t0+a3)**(uni-a4)*log(t0+a3))
c
      t3=t1-t2
      t4=0.0
      if(b(7).eq.uni) sasumq=log(t3+a6)-log(t4+a6)
      if(b(7).gt.uni) sasumq=
     & (uni/(t3+a6)**(a7-uni)-uni/(t4+a6)**(a7-uni))/(uni-a7)
      if(b(7).lt.uni) sasumq=
     & ((t3+a6)**(uni-a7)-(t4+a6)**(uni-a7))/(uni-a7)
      gs5=sasumq
      if(b(7).ne.uni) gs6=a5*(uni/(t3+a6)**a7-uni/(t4+a6)**a7)
      if(b(7).eq.uni) gs6=a5*(uni/(t3+a6)-uni/(t4+a6))
      if(b(7).gt.uni) gs7=
     & a5/(uni-a7)**2
     &  *(uni/(t3+a6)**(a7-uni)-uni/(t4+a6)**(a7-uni))
     & +a5/(uni-a7)
     &  *(-log(t3+a6)/(t3+a6)**(a7-uni)+log(t4+a6)/(t4+a6)**(a7-uni))
      if(b(7).lt.uni) gs7=
     & a5/(uni-a7)**2
     &  *((t3+a6)**(uni-a7)-(t4+a6)**(uni-a7))
     & +a5/(uni-a7)
     &  *(-(t3+a6)**(uni-a7)*log(t3+a6)+(t4+a6)**(uni-a7)*log(t4+a6))
c
      go to 240
c
c
   50 continue
      ifg=1
      f=1.0d30
      return
c
c
  240 continue
      f=f1-a1*(t1-t0)-a2*sasump-a5*sasumq
      g(1)=gg1-gs1
      g(2)=gg2-gs2
      g(3)=gg3-gs3
      g(4)=gg4-gs4
      g(5)=gg5-gs5
      g(6)=gg6-gs6
      g(7)=gg7-gs7
      f=-f
      h(1)=-g(1)
      h(1)=h(1)*a1
      if(b(1).eq.0.0) h(1)=0.0
      h(2)=-g(2)
      h(2)=h(2)*a2
      if(b(2).eq.0.0) h(2)=0.0
      h(3)=-g(3)
      if(b(6).eq.0.0) h(3)=-(g(3)+g(6))
      h(3)=h(3)*a3
      if(b(3).eq.0.0) h(3)=0.0
      h(4)=-g(4)
      if(b(7).eq.0.0) h(4)=-(g(4)+g(7))
      h(4)=h(4)*a4
      if(b(4).eq.0.0) h(4)=0.0
      h(5)=-g(5)
      h(5)=h(5)*a5
      if(b(5).eq.0.0) h(5)=0.0
      h(6)=-g(6)
      h(6)=h(6)*a6
      if(b(6).eq.0.0) h(6)=0.0
      h(7)=-g(7)
      h(7)=h(7)*a7
      if(b(7).eq.0.0) h(7)=0.0
      ff=f
      na=0
      do 800 i=1,n
  800 if(b(i).ne.0.0) na=na+1
      aic=ff+na
    3 format(1h ,110x,d18.10)
    1 format(1h ,7d18.10)
      return
  119 continue
      f=1.0d50
      ifg=1
      return
      end
c
c
c
      subroutine func10(n,b,f,h,ifg)
c-----------------------------------------------------------------------
c     likelihood function of the omori's poisson process
c     lammbda = a1 + a2/(a3+t)**a4 + a5/(a6+t-t2)**a7 * i(t.gt.t2)
c                  + a8/(a9+t-t3)**a10 * i(t.gt.t3)
c     a model for detecting the third aftershock
c     this is not succeeded at the moment (82/11/22)
c-----------------------------------------------------------------------
      implicit real * 8 (a-h,o-z)
      common/xyod/xdumy,xx(19999)
      common/ddd/ff,aic
      common/y/y(50)
      common/range/t0,t1,t2,t3
      common/grad/g(50)
      common t,nn,nfunct
      dimension b(50),h(50)
      do 111 i=1,10
      if(b(i).gt.270.d0) go to 50
  111 continue
      ifg=0
      a1=0.0
      a1=b(1)**2
      if(b(1).ne.0.0) a1=exp(b(1))
      a2=0.0
      a2=b(2)**2
      if(b(2).ne.0.0) a2=exp(b(2))
      a3=0.0
      a3=b(3)**2
      if(b(3).ne.0.0) a3=exp(b(3))
      a4=0.0
      a4=b(4)**2
      if(b(4).eq.0.0) a4=exp(b(10))
      if(b(4).ne.0.0) a4=exp(b(4))
      a5=0.0
      a5=b(5)**2
      if(b(5).ne.0.0) a5=exp(b(5))
      a6=0.0
      a6=b(6)**2
      if(b(6).ne.0.0) a6=exp(b(6))
      a7=0.0
      a7=b(7)**2
      if(b(7).ne.0.0) a7=exp(b(7))
      a8=0.0
      a8=b(8)**2
      if(b(8).ne.0.0) a8=exp(b(8))
      a9=0.0
      a9=b(9)**2
      if(b(9).ne.0.0) a9=exp(b(9))
      a10=0.0
      a10=b(10)**2
      if(b(10).ne.0.0) a10=exp(b(10))
      if(b(7).eq.0.0) a7=exp(b(4))
      if(b(10).eq.0.0) a10=exp(b(4))
      if(b(6).eq.0.0) a6=exp(b(3))
      if(b(9).eq.0.0) a9=exp(b(3))
      if(a3.gt.10000.) go to 119
      uni=1.0d00
      f1=0.0
      gg1=0.0
      gg2=0.0
      gg3=0.0
      gg4=0.0
      gg5=0.0
      gg6=0.0
      gg7=0.0
      gg8=0.0
      gg9=0.0
      gg10=0.0
      if(a7*log(a6+t1-t2).gt.150.) go to 50
      if(a4*log(a3).lt.-150.) go to 50
      if(a4*log(a3+t1).gt.150.) go to 50
      do 20 i=1,nn
      ramdai=a1+a2/(a3+xx(i))**a4
      if(xx(i).gt.t2)
     &ramdai=a1+a2/(a3+xx(i))**a4+a5/(a6+xx(i)-t2)**a7
      if(xx(i).gt.t3)
     &ramdai=a1+a2/(a3+xx(i))**a4+a5/(a6+xx(i)-t2)**a7
     &         +a8/(a9+xx(i)-t3)**a10
      if(ramdai.le.0.0) go to 50
      gg1=gg1+uni/ramdai
      gg2=gg2+uni/ramdai/(a3+xx(i))**a4
      gg3=gg3-a2*a4/ramdai/(a3+xx(i))**(a4+uni)
      gg4=gg4-a2*log(a3+xx(i))/ramdai/(a3+xx(i))**a4
c
      if(xx(i).le.t2) go to 10
      if(ramdai.le.0.0) go to 50
      gg5=gg5+uni/ramdai/(a6+xx(i)-t2)**a7
      gg6=gg6-a5*a7/ramdai/(a6+xx(i)-t2)**(a7+uni)
      gg7=gg7-a5*log(a6+xx(i)-t2)/ramdai/(a6+xx(i)-t2)**a7
c
      if(xx(i).le.t3) go to 10
      if(ramdai.le.0.0) go to 50
      gg8=gg8+uni/ramdai/(a9+xx(i)-t3)**a10
      gg9=gg9-a8*a10/ramdai/(a9+xx(i)-t3)**(a10+uni)
      gg10=gg10-a8*log(a9+xx(i)-t3)/ramdai/(a9+xx(i)-t3)**a10
c
   10 continue
      f1=f1+log(ramdai)
   20 continue
c
c
      if(b(4).eq.uni) sasump=log(t1+a3)-log(t0+a3)
      if(b(4).gt.uni) sasump=
     & (uni/(t1+a3)**(a4-uni)-uni/(t0+a3)**(a4-uni))/(uni-a4)
      if(b(4).lt.uni) sasump=
     & ((t1+a3)**(uni-a4)-(t0+a3)**(uni-a4))/(uni-a4)
      gs1=t1-t0
      gs2=sasump
      if(b(4).ne.uni) gs3=a2*(uni/(t1+a3)**a4-uni/(t0+a3)**a4)
      if(b(4).eq.uni) gs3=a2*(uni/(t1+a3)-uni/(t0+a3))
      if(b(4).gt.uni) gs4=
     & a2/(uni-a4)**2
     &  *(uni/(t1+a3)**(a4-uni)-uni/(t0+a3)**(a4-uni))
     & +a2/(uni-a4)
     &  *(-log(t1+a3)/(t1+a3)**(a4-uni)+log(t0+a3)/(t0+a3)**(a4-uni))
      if(b(4).lt.uni) gs4=
     & a2/(uni-a4)**2
     &  *((t1+a3)**(uni-a4)-(t0+a3)**(uni-a4))
     & +a2/(uni-a4)
     &  *(-(t1+a3)**(uni-a4)*log(t1+a3)+(t0+a3)**(uni-a4)*log(t0+a3))
c
      t12=t1-t2
      t02=0.0
      if(b(7).eq.uni) sasumq=log(t12+a6)-log(t02+a6)
      if(b(7).gt.uni) sasumq=
     & (uni/(t12+a6)**(a7-uni)-uni/(t02+a6)**(a7-uni))/(uni-a7)
      if(b(7).lt.uni) sasumq=
     & ((t12+a6)**(uni-a7)-(t02+a6)**(uni-a7))/(uni-a7)
      gs5=sasumq
      if(b(7).ne.uni) gs6=a5*(uni/(t12+a6)**a7-uni/(t02+a6)**a7)
      if(b(7).eq.uni) gs6=a5*(uni/(t12+a6)-uni/(t02+a6))
      if(b(7).gt.uni) gs7=
     & a5/(uni-a7)**2
     &  *(uni/(t12+a6)**(a7-uni)-uni/(t02+a6)**(a7-uni))
     & +a5/(uni-a7)
     & *(-log(t12+a6)/(t12+a6)**(a7-uni)+log(t02+a6)/(t02+a6)**(a7-uni))
      if(b(7).lt.uni) gs7=
     & a5/(uni-a7)**2
     &  *((t12+a6)**(uni-a7)-(t02+a6)**(uni-a7))
     & +a5/(uni-a7)
     & *(-(t12+a6)**(uni-a7)*log(t12+a6)+(t02+a6)**(uni-a7)*log(t02+a6))
c
      t13=t1-t3
      t03=0.0
      if(b(10).eq.uni) sasumr=log(t13+a9)-log(t03+a9)
      if(b(10).gt.uni) sasumr=
     & (uni/(t13+a9)**(a10-uni)-uni/(t03+a9)**(a10-uni))/(uni-a10)
      if(b(10).lt.uni) sasumr=
     & ((t13+a9)**(uni-a10)-(t03+a9)**(uni-a10))/(uni-a10)
      gs8=sasumr
      if(b(10).ne.uni) gs9=a8*(uni/(t13+a9)**a10-uni/(t03+a9)**a10)
      if(b(10).eq.uni) gs9=a8*(uni/(t13+a9)-uni/(t03+a9))
      if(b(10).gt.uni) gs10=
     & a8/(uni-a10)**2
     &  *(uni/(t13+a9)**(a10-uni)-uni/(t03+a9)**(a10-uni))
     & +a8/(uni-a10)
     & *(-log(t13+a9)/(t13+a9)**(a10-uni)+log(t03+a9)/(t03+a9)**(a10-uni
     &))
      if(b(10).lt.uni) gs10=
     & a8/(uni-a10)**2
     &  *((t13+a9)**(uni-a10)-(t03+a9)**(uni-a10))
     & +a8/(uni-a10)
     & *(-(t13+a9)**(uni-a10)*log(t13+a9)+(t03+a9)**(uni-a10)*log(t03+a9
     &))
c
      go to 240
c
c
   50 continue
      ifg=1
      f=1.0d30
      return
c
c
  240 continue
      f=f1-a1*(t1-t0)-a2*sasump-a5*sasumq-a8*sasumr
      g(1)=gg1-gs1
      g(2)=gg2-gs2
      g(3)=gg3-gs3
      g(4)=gg4-gs4
      g(5)=gg5-gs5
      g(6)=gg6-gs6
      g(7)=gg7-gs7
      g(8)=gg8-gs8
      g(9)=gg9-gs9
      g(10)=gg10-gs10
      if(b(4).eq.1.0d00) g(4)=0.0
      if(b(7).eq.1.0d00) g(7)=0.0
      f=-f
      h(1)=-g(1)
      h(1)=h(1)*a1
      if(b(1).eq.0.0) h(1)=0.0
      h(2)=-g(2)
      h(2)=h(2)*a2
      if(b(2).eq.0.0) h(2)=0.0
      h(3)=-g(3)
      if(b(6).eq.0.0) h(3)=-(g(3)+g(6))
      if(b(6).eq.0.0.and.b(9).eq.0.0) h(3)=-(g(3)+g(6)+g(9))
      h(3)=h(3)*a3
      if(b(3).eq.0.0) h(3)=0.0
      h(4)=-g(4)
      if(b(7).eq.0.0.and.b(10).ne.0.0) h(4)=-(g(4)+g(7))
      if(b(7).ne.0.0.and.b(10).eq.0.0) h(4)=-(g(4)+g(10))
      if(b(7).eq.0.0.and.b(10).eq.0.0) h(4)=-(g(4)+g(7)+g(10))
      h(4)=h(4)*a4
      if(b(4).eq.0.0) h(4)=0.0
      h(5)=-g(5)
      h(5)=h(5)*a5
      if(b(5).eq.0.0) h(5)=0.0
      h(6)=-g(6)
      h(6)=h(6)*a6
      if(b(6).eq.0.0) h(6)=0.0
      h(7)=-g(7)
      if(b(10).eq.0.0) h(7)=-(g(7)+g(10))
      h(7)=h(7)*a7
      if(b(7).eq.0.0) h(7)=0.0
      h(8)=-g(8)
      h(8)=h(8)*a8
      if(b(8).eq.0.0) h(8)=0.0
      h(9)=-g(9)
      h(9)=h(9)*a9
      if(b(9).eq.0.0) h(9)=0.0
      h(10)=-g(10)
      if(b(4).eq.0) h(10)=-(g(4)+g(10))
      h(10)=h(10)*a10
      if(b(10).eq.0.0) h(10)=0.0
      ff=f
      na=0
      do 800 i=1,n
  800 if(b(i).ne.0.0) na=na+1
      aic=ff+na
    3 format(1h ,110x,d18.10)
    1 format(1h ,7d18.10)
      go to 300
  119 continue
      f=1.0d50
      ifg=1
  300 continue
      return
      end
c
c
c
      function sf1(x,q)
      implicit real * 8 (a-h,o-z)
      sf1=x**(1.d0-q)/(1.d0-q)
      return
      end
      function sf11(x,q)
      implicit real * 8 (a-h,o-z)
      sf11=log(x)
      return
      end
      function sf2(x,q)
      implicit real * 8 (a-h,o-z)
      sf2=sf1(x,q)*(log(x)-1.d0/(1.d0-q))
      return
      end
      function sf21(x,q)
      implicit real * 8 (a-h,o-z)
      sf21=(log(x))**2/2.d0
      return
      end
      function sf3(x,q)
      implicit real * 8 (a-h,o-z)
      sf3=sf1(x,q)*(log(x))**2-2.d0/(1.d0-q)*sf2(x,q)
      return
      end
      function sf31(x,q)
      implicit real * 8 (a-h,o-z)
      sf31=(log(x))**3/3.d0
      return
      end
c
c     real function gm*8 (x,q,c)
      real*8 function gm (x,q,c)
      implicit real * 8 (a-h,o-z)
      gm=0.0
      if(x.eq.c) go to 20
      gmi=x**(-q)
      do 10 i=1,50
      i1=i-1
      if(i1.eq.0) i1=1
c     gmi=gmi*x/i1
      gmi=gmi*(x-c)/i1
      gm=gm+(-1)**(i-1)*gmi/(i-q)
      if(gmi/gm.lt.1.d-13) go to 20
   10 continue
   20 return
      end
c
c     real function dgm*8 (x,q,c)
      real*8 function dgm (x,q,c)
      implicit real*8(a-h,o-z)
      dgm=0.0
      if(x.eq.c) go to 30
      dgmi=x**(-q)
      do 10 i=1,50
      i1=i-1
      if(i1.eq.0) i1=1
c     dgmi=dgmi*x/i1
      dgmi=dgmi*(x-c)/i1
      dgm=dgm+(-1)**i*dgmi/(i-q)**2
      if(dgmi/dgm.lt.1.d-13) go to 20
   10 continue
c  20 dgm=-dgm-gm(x,q)*log(x)
   20 dgm=-dgm-gm(x,q,c)*log(x)
   30 continue
      return
      end
c
c
c
      subroutine fisher(b,n,h)
c-----------------------------------------------------------------------
c     fisher information matrix of the modified omori model
c     lammbda = a1 + a2/(a3+t)**a4
c-----------------------------------------------------------------------
      implicit real * 8 (a-h,o-z)
      common/range/t0,t1,t2,t3
      common t,nn,nfunct
      dimension b(50),h(50,50)
      a1=b(1)**2
      a2=b(2)**2
      a3=b(3)**2
      a4=b(4)**2
c---- careful about q = 1 in the functions ---
      if(a1.ne.0.0) h(1,1)=(t1-t0)/a1-log(t1-t0)/a1
      if(a1.eq.0.0) h(1,1)=1.0d0
      h(1,2)=0.0
      h(1,3)=0.0
      h(1,4)=0.0
      h(2,2)=(sf1(t1+a3,a4)-sf1(t0+a3,a4))/a2
      h(2,3)=-a4*(sf1(t1+a3,a4+1)-sf1(t0+a3,a4+1))
      h(2,4)=-(sf2(t1+a3,a4)-sf2(t0+a3,a4))
      h(3,3)=a2*a4**2*(sf1(t1+a3,a4+2)-sf1(t0+a3,a4+2))
      h(3,4)=a2*a4*(sf2(t1+a3,a4+1)-sf2(t0+a3,a4+1))
      h(4,4)=a2*(sf3(t1+a3,a4)-sf3(t0+a3,a4))
      do 10 i=1,4
      do 10 j=i,4
   10 h(j,i)=h(i,j)
      return
      end
c
c
c
      subroutine sizes(n,x)
      implicit real * 8 (a-h,o-z)
      dimension x(50),cls(20),ti(20)
      dimension ak(20),p(20),c(20)
      common/range/t0,t1,t2,t3
      ti(1)=t2
      ti(2)=t3
      kn=(n-1)/3
      do 10 k=1,kn
      ak(k)=x(3*k-1)
      c(k)= x(3*k)
      if(c(k).eq.0.0) c(k)=c(k-1)
      p(k)= x(3*k+1)
      if(p(k).eq.0.0) p(k)=p(k-1)
   10 continue
      cls(1)=ak(1)*((t1+c(1))**(1-p(1))-c(1)**(1-p(1)))/(1-p(1))
      if(p(1).eq.1.d0) cls(1)=ak(1)*(log(t1+c(1))-log(c(1)))
      if(kn.eq.1) go to 25
      do 20 k=2,kn
      if(p(k).eq.1.d0) then
      cls(k)=ak(k)*(log(t1-ti(k-1)+c(k))-log(c(k)))
      else
      cls(k)=ak(k)*((t1-ti(k-1)+c(k))**(1-p(k))-c(k)**(1-p(k)))/(1-p(k))
      end if
   20 continue
   25 continue
      write(6,3)
      write(6,4) t0,ak(1),c(1),p(1),cls(1)
      if(kn.eq.1) go to 35
      do 30 k=2,kn
      write(6,4) ti(k-1),ak(k),c(k),p(k),cls(k)
   30 continue
   35 continue
    3 format(1h ,'    ti         k          c         p         cls')
    4 format(1h ,5e11.4)
      return
      end
      
