      program eptren
c
c     this program carries out maximum likelihood estimates of
c   intensity rates of either exponential plynomial or exponential
c   fourier series of non-stationary poisson process models.
c   the optimal oreders are selected automatically by the use of
c   akaike information criterion(aic). a subroutine for plotting
c   a graph of the estimated lines are attached.
c
c     structure of the program
c
c          eptren
c             i---inputs
c             i---reduc1
c             i---reduc2
c             i---davidn----------funct
c             i      i---hesian
c             i      +---linear---funct
c             i---fincal
c             +---output---printr---trenfn
c                             +-----cyclfn
c
c     this program is designed by y. ogata, and programed by y. ogata
c   and k. katsura, inst. statist. math., tokyo. (31/01/85)
c
c     reference of the model, see
c   maclean,c.j.(1974) "estimation and testing og an exponential
c   polynomial rate function within the non-stationary poisson
c   process." biometrika 61, pp.81-86.
c
      implicit real * 8 (a-h,o-z)
      dimension xx(9000),x(100),aic(20),xa(100,20)
      dimension amg(9000)
      external funct1,funct2
      call inputs(xx,amg,t,nn,nfunct,ipl,x,nb,ni,cycle)
      if(nfunct.eq.1) call reduc1(t,xx,nn,nb,ni)
      if(nfunct.eq.2) call reduc2(t,xx,nn,nb,ni,cycle)
      do 10 i=1,nb
      ii=i
      if(nfunct.eq.2) ii=2*i-1
      do 20 j=1,ii
      x(j)=0.0
   20 continue
      if(nfunct.eq.1) call davidn(x,ii,ihes,funct1)
      if(nfunct.eq.2) call davidn(x,ii,ihes,funct2)
      call fincal(ii,x,aic(i),xa(1,i),t,nfunct)
   10 continue
      call output(xx,amg,nn,t,xa,aic,nb,ipl,nfunct,cycle)
      stop
      end
      subroutine inputs(xx,amg,t,nn,nfunct,ipl,x,nb,ni,cycle)
c
c  inputs
c
c     nfunct=1:  exponential polynomial trend fitting.
c           =0:  exponential fourier series fitting.
c     ipl=1:  xy-plotter is given, otherwise not.
c     in =5:  card reading of data and input parameters.
c     fmt:  reading format of point process data.
c     t:  length of observed time interval of points.
c     nn:  number of observed points in (0,t).
c     xx(i):  point process data.
c     nb:  number of parameters.
c     ni:  number of subdivisions in either (0,t) or (0,cycle)
c          for the numerical integration of intensity.
c     cycle:  periodicity to be investigated.
c
c     it should be remarked that no initial estimates for the parameters
c   are necessary for the input.  for the reason for this see
c   y. ogata (1978), "the asymptotic behaviour of maximum likelihood
c   estimators for stationary point processes." ann. inst. statist.
c   math. vol. 30, page 255.
c
      implicit real * 8 (a-h,o-z)
      dimension xx(1),x(1),amg(1)
c     dimension fmt(10)
      character*8 fmt(10)
      open(4,file='eptren.open')
      read(4,*) nfunct,ires,am1
      read(4,*) ipl
      read(4,*) in
    1 format(8i10)
c     read(4,2) fmt
    2 format(10a8)
      if(ires.ne.1) open(in,file='work.etas')
      if(ires.eq.1) open(in,file='work.res')
c     read(in,3) t,nn
c     read(in,*) t,nn
      if(ires.ne.1) read(in,2) fmt
    3 format(f10.0,i10)
c     write(6,*) fmt
c     read(in,fmt) (xx(i),i=1,nn)
      i=1
   10 continue
      if(ires.ne.1) read(10,*,end=20) idum,e,an,amg(i),xx(i),dep    
      if(ires.eq.1) read(10,*,end=20) idum,e,an,amg(i),zi,dep,xx(i)    
 1002 format(i5,4f12.5,f5.0)
      i=i+1
      go to 10
   20 n1=i-1
      t=xx(n1)
      nn=0
      do 30 i=1,n1
      if(amg(i).lt.am1) go to 30
      if(xx(i).lt.0.0) go to 30
      nn=nn+1
      xx(nn)=xx(i)
      amg(nn)=amg(i)
   30 continue
c     read(in,*) (xx(i),i=1,nn)
      ipl=0
c     write(6,4) nfunct,ipl,in,fmt
c   4 format('funct =',i5,5x,'ipl =',i5,5x,' in =',i5//
c    &           'fmt = ',5x,10a8)
      write(6,4) nfunct,ipl,in
    4 format('funct =',i5,5x,'ipl =',i5,5x,' in =',i5/)
      write(6,5) t,nn
    5 format(/'t =',f15.5,'     nn =',i10)
      write(6,6) (xx(i),i=1,nn)
    6 format(//'data set (xx(i),i=1,nn)'/(5f12.5))
      read(4,1) nb,ni
      if(nfunct.eq.2) read(4,*) cycle
      do 50 i=1,nb
   50 x(i)=0.0
    7 format(8f10.0)
      write(6,8) nb,ni
    8 format(//' no. of parameters ',i5,5x,'no. of subdiv.',i10///)
      if(nfunct.eq.2) write(6,9) cycle
    9 format(//' cycle =',f15.5)
      return
      end
      subroutine reduc1(t,xx,nn,nb,ni)
c
c     reduction of the data for the (exp)polynomial trend analysis
c
c     input
c           t: length of the observed time interval
c           xx: data of events
c           nn: number of data
c           nb: the largest number of parameters
c           ni: number of subdivisions in (0,t)
c     output
c           rxz:  values of polynomial elements at each events
c           sxz:  values of polynomial at each selected points
c           delta:  length of subdivision by the selected points
      implicit real * 8(a-h,o-z)
      common /rd1fn1/delta,rxz(20),sxz(4001,20),ns
      dimension xx(1)
      ns=ni
      do 10 j=1,nb
      rxz(j)=0.0
      do 10 i=1,nn
      rxz(j)=rxz(j)+(xx(i)/t)**(j-1)
   10 continue
      delta=1.0d0/ni
      sxz(1,1)=1.d0
      do 30 j=2,nb
   30 sxz(1,j)=0.0
      do 20 i=2,ni+1
      do 20 j=1,nb
      sxz(i,j)=((i-1)*delta)**(j-1)
   20 continue
      return
      end
      subroutine funct1(n,a,f,g,ifg)
c
c     negative log likelihood function and gradients for log linearly
c     parametrized intensity (exponential polynomial trend)
c
c   input
c        n: number of parameters
c        ns: number of subdivision in (0,1)
c        a: parameters
c        delta: length of subdivisions in (0,1)
c        rxz: values of elements at each events
c        sxz: values of elements at each selected points
c   output
c        f: value of the function
c        g: gradient vecter
c        ifg: index for restrictions
c
      implicit real * 8 (a-h,o-z)
      common /rd1fn1/delta,rxz(20),sxz(4001,20),ns
      common     / ddd /  r , ff , aic , sd
      dimension g(1),gs(30),a(1)
      ifg=0
      fxx=0.0
      do 20 j=1,n
   20 fxx=fxx+a(j)*rxz(j)
c
      ssxx=0.0
      do 10 j=1,n
   10 gs(j)=0.0
      do 30 i=1,ns+1
      rmd=0.0
      do 40 j=1,n
   40 rmd=rmd+a(j)*sxz(i,j)
      if(rmd.gt.100.d0) go to 80
      exprm=exp(rmd)
      exprmd=exprm
      if(i.eq.1.or.i.eq.ns+1) exprmd=exprm/2
      ssxx=ssxx+exprmd
      do 50 j=1,n
   50 gs(j)=gs(j)+sxz(i,j)*exprmd
   30 continue
      f=-fxx+ssxx*delta
      ff=f
      do 60 j=1,n
   60 g(j)=-rxz(j)+gs(j)*delta
      return
   80 continue
      f=1.d30
      ifg=1
      return
      end
      subroutine reduc2(t,xx,nn,nb,ni,cycle)
c
c   reduction of the data for the exponential fourier trend analysis
c
c     input
c           t: length of the observed time interval
c           cycle: periodicity
c           xx: data of events
c           nn: number of data
c           nb: the largest number of parameters
c           ni: number of subdivisions in (0,cycle)
c     output
c           rxz:  values of fourier elements at each events
c           sxz:  values of fourier at each selected points
c           delta:  length of subdivision by the selected points
c
      implicit real * 8(a-h,o-z)
      common /rd2fn2/delta,rxc(20),sxc(4001,20),rxs(20),sxs(4001,20),tr,
     &               it,ns,nnd
      dimension xx(1)
      data pi/3.14159265358979d0/
      ns=ni
      nnd=nn
c     cycle=365.d0
      it=int(t/cycle)
      tr=t-it*cycle
      do 10 k=1,nb
      rxc(k)=0.0
      rxs(k)=0.0
      do 10 i=1,nn
      rxc(k)=rxc(k)+cos(2*pi*k*xx(i)/cycle)
      rxs(k)=rxs(k)+sin(2*pi*k*xx(i)/cycle)
   10 continue
      delta=cycle/ni
      do 30 k=2,nb
      sxc(1,k)=0.0
   30 sxs(1,k)=0.0
      do 20 i=1,ni+1
      do 20 k=1,nb
c     sxc(i,k)=cos(2*pi*k*(i-0.5d0)*delta/cycle)
      sxc(i,k)=cos(2*pi*k*(i-1)*delta/cycle)
c     sxs(i,k)=sin(2*pi*k*(i-0.5d0)*delta/cycle)
      sxs(i,k)=sin(2*pi*k*(i-1)*delta/cycle)
   20 continue
      return
      end
      subroutine funct2(n,a,f,g,ifg)
c
c     negative log likelihood function and gradients for log linearly
c     parametrized intensity (periodicity)
c
c     input
c        n: number of parameters
c        nnd: number of events
c        ns: number of subdivision in (0,cycle)
c        a: parameters
c        delta: length of subdivisions in (0,cycle)
c        rxz: values of elements at each events
c        sxz: values of elements at each selected points
c     output
c        f: value of the function
c        g: gradient vecter
c        ifg: index for restrictions
c
      implicit real * 8 (a-h,o-z)
      common /rd2fn2/delta,rxc(20),sxc(4001,20),rxs(20),sxs(4001,20),tr,
     &               it,ns,nnd
      common     / ddd /  r , ff , aic , sd
      dimension g(20),a(1)
      dimension gs(20),gc(20),gcp(20),gsp(20)
      ifg=0
      fxx=a(1)*nnd
      n2=(n-1)/2
      if(n2.eq.0) go to 25
      do 20 j=1,n2
   20 fxx=fxx+a(2*j)*rxc(j)+a(2*j+1)*rxs(j)
   25 continue
c
      ssxx=0.0
      ssxxp=0.0
      g(1)=1.0d0
      if(n2.eq.0) go to 15
      do 10 j=1,n2
      gc(j)=0.0
      gcp(j)=0.0
      gs(j)=0.0
   10 gsp(j)=0.0
   15 continue
      do 30 i=1,ns+1
      rmd=a(1)
      if(n2.eq.0) go to 45
      do 40 j=1,n2
   40 rmd=rmd+a(2*j)*sxc(i,j)+a(2*j+1)*sxs(i,j)
   45 continue
      if(rmd.gt.100.d0) go to 80
      exprm=exp(rmd)
      exprmd=exprm
      if(i.eq.1.or.i.eq.ns+1) exprmd=exprm/2
      ssxx=ssxx+exprmd
      if(i*delta.le.tr) ssxxp=ssxx
      if(n2.eq.0) go to 55
      do 50 j=1,n2
      gc(j)=gc(j)+sxc(i,j)*exprmd
      if(i*delta.le.tr) gcp(j)=gc(j)
      gs(j)=gs(j)+sxs(i,j)*exprmd
      if(i*delta.le.tr) gsp(j)=gs(j)
   50 continue
   55 continue
   30 continue
c
      f=-fxx+(it*ssxx+ssxxp)*delta
      g(1)=-nnd+(it*ssxx+ssxxp)*delta
      if(n2.eq.0) go to 65
      do 60 j=1,n2
      g(2*j)=-rxc(j)+(it*gc(j)+gcp(j))*delta
   60 g(2*j+1)=-rxs(j)+(it*gs(j)+gsp(j))*delta
   65 continue
      ff=f
      return
   80 continue
      f=1.d30
      ifg=1
      return
      end
      subroutine  davidn( x,n,ihes,funct )
c
c     minimization by davidon-fletcher-powell procedure
c
c     this subroutine was copied from timsac 78.
c
c     ----------------------------------------------------------------
c     the following subroutines are directly called by this subroutine
c          funct
c          hesian
c          linear
c     ----------------------------------------------------------------
c     inputs:
c            x: vector of initial values
c            k: dimension of the vector x
c            ihes: =0   inverse of hessian matrix is not available
c                  =1   inverse of hessian matrix is available
c
c     output:
c            x: vector of minimizing solution
c
      implicit  real * 8  ( a-h , o-z )
      dimension  x(82) , dx(82) , g(82) , g0(82) , y(82)
      dimension  h(82,82) , wrk(82) , s(82)
      common     / ccc /  isw,ipr
      common     / ddd /  r , f , aic , sd
      external funct
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
c     write( 6,340 )     xm , sd , aic
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
      ds2 = sqrt(s2)
      gtem = abs(s1) / ds2
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
      write( 6,330 )     ramda , f
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
      if( sqrt(s2) .gt. tau2 )   go to  11111
      if( xmb/xm-1.d0 .lt. eps1  .and.  sqrt(s1) .lt. eps2 )  go to 900
11111 continue
      if( icc .ge. 5 )     go to 900
      go to 11110
  900 continue
      write( 6,600 )
      write( 6,610 )     (x(i),i=1,n)
      write( 6,601 )
      write( 6,610 )     (g(i),i=1,n)
      return
  330 format( 'lambda =',d15.7,5x,'log likelihood =',d23.15)
  340 format( 28x,'log-likelihood =',d23.15,5x,'sd =',d22.15,5x,
     *  'aic =',d23.15 )
  600 format( /'-----  x  -----' )
  601 format( /'***  gradient  ***' )
  610 format( 5d13.5 )
      end
c
c
c
      subroutine  linear( x,h,ram,ee,k,ig,funct )
c
c     this subroutine was copied from timsac 78.
c
c     this subroutine performs the linear search along the direction
c     specified by the vector h
c   -----------------------------------------------------------------
c     the following subroutine is directly called by this subroutine:
c         funct
c   -----------------------------------------------------------------
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
      dimension  x(100) , h(100) , x1(82)
      dimension  g(82)
      common     / ccc /  isw , ipr
c
      isw = 1
      ipr = 7
      if( ram .le. 1.0d-30 )  ram = 0.01d0
      const2 = 1.0d-60
      hnorm = 0.d0
      do 10  i=1,k
   10 hnorm = hnorm + h(i)**2
      hnorm = sqrt( hnorm )
c
      ram2 = ram
      e1 =ee
      ram1 = 0.d0
c
      do 20  i=1,k
   20 x1(i) = x(i) + ram2*h(i)
      call  funct( k,x1,e2,g,ig )
      if(ipr.ge.7)  write(6,2)  ram2,e2
c
      if( ig .eq. 1 )  go to  50
      if( e2 .gt. e1 )  go to 50
   30 ram3 = ram2*2.d0
      do 40  i=1,k
   40 x1(i) = x(i) + ram3*h(i)
      call  funct( k,x1,e3,g,ig )
      if( ig.eq.1 )  go to  500
      if( ipr.ge.7 )  write(6,3)  ram3,e3
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
      if(ipr.ge.7)  write(6,4)  ram2,e2
c     write(6,8) e1,e2,e3,ee
    8 format(1h ,5d20.10)
      if( e2.gt.e1 )  go to 50
c
   70 assign 80 to return
      go to 200
c
   80 do 90  i=1,k
   90 x1(i) = x(i) + ram*h(i)
      call  funct( k,x1,ee,g,ig )
      if(ipr.ge.7)  write(6,5)  ram,ee
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
      if( ipr.ge.7 )  write(6,6)  ram,ee
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
      if( ipr.ge.7 )  write(6,7)  ram,e3
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
      subroutine fincal(n,x,aic,xa,t,nfunct)
      implicit real * 8 (a-h,o-z)
      common     / ddd /  r , xm , aicc , sd
      dimension x(1),xa(1)
      if(nfunct.eq.2) go to 20
      xa(1)=x(1)-log(t)
      if(n.eq.1) go to 40
      do 10 i=2,n
      xa(i)=x(i)/t**(i-1)
   10 continue
   40 continue
      aic=2*xm+2*n
      return
   20 do 30 i=1,n
      xa(i)=x(i)
   30 continue
      aic=2*xm+2*n
      return
      end
      subroutine output( xx,amg,nn,t,xa,aic,n,ipl,nfunct,cycle )
      implicit real * 8 (a-h,o-z)
      dimension xx(1),amg(1),aic(1),xa(100,20)
      acmin=1.d10
      do 10 i=1,n
      write(6,1) aic(i)
    1 format(////' a i c ',f15.2)
      ii=i
      if(nfunct.eq.2) ii=2*i-1
      write(6,2) (xa(j,i),j=1,ii)
    2 format(//' parameters '/(2x,5d13.5))
      if(acmin.lt.aic(i)) go to 10
      acmin=aic(i)
      imin=i
   10 continue
      iimin=imin
      tt=t
      if(nfunct.eq.2) iimin=2*imin-1
      if(nfunct.eq.2) tt=cycle
      if(nfunct.eq.1) open(7,file='out.eptrend1')
      if(nfunct.eq.1) write(7,*) 'eptren1'
      if(nfunct.eq.2) open(7,file='out.epcycle1')
      if(nfunct.eq.2) write(7,*) 'epcycle1'
c     x1=1.0
      do 20 i=1,nn
      xd=xx(i)
      if(nfunct.eq.2) xd=xx(i)-int(xx(i)/tt)*tt
      x1=amg(i)
      write(7,*) xd,x1
   20 continue
      close(7)
      call printr(tt,xa(1,imin),iimin,nfunct)
      write(6,3) iimin,(xa(i,imin),i=1,iimin)
    3 format(i10/(3d21.13))
      return
      end
      subroutine trenfn(xa,x,y,n)
      real * 8 xa(1),yy,x,y
      yy=xa(1)
      if(n.eq.1) go to 20
      do 10 i=2,n
       yy=yy+xa(i)*x**(i-1)
   10 continue
   20 continue
      y=exp(yy)
      return
      end
      subroutine cyclfn(xa,x,y,n)
      real * 8 xa(1),yy,x,y,pi
      data pi/3.14159265358979d0/
      yy=xa(1)
      if(n.eq.1) go to 20
      n2=(n-1)/2
      do 10 i=1,n2
       yy=yy+xa(2*i)*cos(2*pi*i*x)+xa(2*i+1)*sin(2*pi*i*x)
   10 continue
   20 continue
      y=exp(yy)
      return
      end
      subroutine printr(t,xa,n,nfunct)
      real * 8 t,xa,xx,yy
      dimension x(2000),y(2000),xa(1)
      character*1 xo
      data xo/'o'/
      nn=101
      ymin=0.0
      ymax=0.0
      do 10 i=1,nn
      x(i)=t*(i-1)/(nn-1)
      xx=x(i)
      if(nfunct.eq.2) xx=1.d0*(i-1)/nn
      if(nfunct.eq.1) call trenfn(xa,xx,yy,n)
      if(nfunct.eq.2) call cyclfn(xa,xx,yy,n)
      y(i)=yy
      if(ymin.gt.y(i)) ymin=y(i)
      if(ymax.lt.y(i)) ymax=y(i)
   10 continue
      x1=ymin
      x2=ymin+(ymax-ymin)*2/4
      x3=ymax
      write(6,1) x1,x2,x3
    1 format(//'x-values  f-values  ',f8.5,2(10x,e11.4)/
     &       21x,'+',4('---------+'))
      if(nfunct.eq.1) open(7,file='out.eptrend2')
      if(nfunct.eq.1) write(7,*) 'eptrend2'
      if(nfunct.eq.2) open(7,file='out.epcycle2')
      if(nfunct.eq.2) write(7,*) 'epcycle2'
      do 20 i=1,nn
      l=(y(i)-ymin)*40/(ymax-ymin)+0.5
      if(l.gt.0) write(6,2) x(i),y(i),(xo,j=1,l)
      if(l.le.0) write(6,2) x(i),y(i)
    2 format(2(1x,d9.4),' i',40a1)
      write(7,*) x(i),y(i)
   20 continue
      close(7)
      return
      end
