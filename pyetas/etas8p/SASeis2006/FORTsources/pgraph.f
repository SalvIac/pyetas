      program pgraph
c
c     this provides the following several graphical outputs for the
c     point process data set:
c
c     1. cumulative numbers of events versus time, and positions of
c        spikes ( subroutine cumlat ).
c     2. series of counted events in the moving interval
c        with a fixed length.  this is a kind of moving average.
c        dotted lines indicate i*(standard error), i=1,2,3, assuming
c        the stationary poisson process ( subroutine count1 or count2 ).
c     3. log survivor curve with i*(standard error), i=1,2,3, assuming
c        the stationary poisson process, and the similar graph in
c        which (x,y) plots are rotated and shifted in such a way
c        that the standard error lines and expectation lines are
c        parallel ( subroutine surviv ).
c     4. empirical distribution of u(i)=exp(-c*x(i)) where x(i)
c        is the i-th interval length between consecutive points, and
c        lines of .95 & .99 significance bands of the two sided
c        kolmogorov-smirnov test assuming the uniform distribution are
c        given.  the related graph of ( u(i),u(i+1) ) plots are also
c        carried out ( subroutine unifrm ).
c     5. estimated intensity  mf(t)  under the palm probability.
c        this is related to the covariance density  c(t)  by the
c        relation of   c(t) = m*(mf(t)-m), where m is the mean
c        intensity of the point process.   0.95 & 0.99 error bands
c        are provided assuming the stationary poisson process
c        ( subroutine palmxy, palmpr ).
c     6. estimation of variance time curve with 0.95 & 0.99 error
c        lines assuming the stationary poisson process
c        ( subroutine vtcxyp, vtcprt ).
c
c   structure of the program
c
c        pgraph
c           i---input
c           i---count1---shimiz
c           i---surviv---errbr2---plsinv
c           i      +-------unifrm---unitsq
c           i---palmpr
c           i---vtc
c           i---vtcprt
c
c     this program was designed by y. ogata, and programmed by
c        k. katsura, inst. statist. math., tokyo. (31/01/85)
c
c   references
c
c   y.ogata & k.shimazaki (1984). "transition from aftershock to normal
c      activity: the 1965 rat islands earthquake aftershock sequence."
c      bulletin of the seismological society of america, vol. 74, no.
c      5, pp. 1757-1765.
c   y. ogata (1985). in preparation.
c
      dimension z(3000),xmg(3000),x(3000),erres(3000)
      dimension errest(3000)
      dimension w(3000),x1(3000),xt(300),sigma(3000)
      data xmg/3000*0.0/
      real*8 a(20),b,bm,zd(3000),xmgd(3000),xmgmax,wop,a55
      real*8 ax(20),ay(20),ac(20),at(20),xx(3000),yy(3000),ymg(3000)
c
      call input(nfunct,isi,icnt,z,xmg,xmgd,yb,t,t0,t1,t2,a,ak,c1,p,
     &           ak2,c2,p2,kkx,xmgmax,nn,ipoint,days,h,zd,delta,dmax,ns,
     &           xx,yy,ymg,c,cm,d,dm,ax,ay,ac,at,kky,kkc,kkt,mm,cycle)
c
      call tmchg0(x,xt,zd,xmgd,tm,ttt,nn,inn,xmin,xmax,t)
c ----------------------------------------------
c
      call count1(x,nn,h,nfunct,xmg,xmgmax,xt,inn,tm)
c
      call surviv(x,nn,nfunct,ttt,isi)
c
      tt=t
      if(nfunct.eq.0) tt=ttt
c
      call palmpr(z,nn,tt,days,ipoint)
c
      if(nfunct.ne.0) call palmpr(x,nn,float(nn),days*nn/tt,ipoint)
c
      call vtc(z,nn,delta,dmax,tt,sigma,k,erres,errest)
c
c
      call vtcprt(sigma,k,delta,nn,tt,erres,errest)
c
      if(nfunct.eq.0) go to 20
      call vtc(x,nn,delta*nn/tt,dmax*nn/tt,float(nn),sigma,k,erres,erres
     &t)
      call vtcprt(sigma,k,delta*nn/tt,nn,float(nn),erres,errest)
   20 continue
c  ---------------------------------------------
      stop
      end
      subroutine input(nfunct,isi,icnt,z,xmg,xmgd,yb,t,t0,t1,t2,a,
     &                 ak,c1,p,ak2,c2,p2,kkx,xmgmax,nn,ipoint,days,h,zd,
     &                 delta,dmax,ns,xx,yy,ymg,c,cm,d,dm,ax,ay,ac,at,
     &                 kky,kkc,kkt,mm,cycle)
c
c  inputs
c
c    nfunct:  model option for data transformation (noise analysis)
c           =0  data set is directly displayed.
c           =1  the modified omori aftershock (incudes 2ndary one) trend
c           =2  the original and modified omori type response self-
c                 exciting series with effection of the magnitudes
c           =3  the laguerre polynomial type response self exciting with
c                 effection of the magnitudes, inputs point process
c                 data, trend and cyclic component (see program lintes)
c           =4  the exponential polynomial trend (see program eptren).
c     iyy:  option whether the input data is read (=1) or not (=0).
c     mags:  option for data reading
c          =0  simple point process data only are read; otherwise
c                magnitude data of earthquake (eqk) attached are read.
c     idep:  option for data reading.
c          =0  depths of eqk are not read; otherwise read.
c     in:  option for reading data set
c        =5  input is made by the card reading.
c     fmt:  format for reading the point process data set.
c     t:  length of the total time interval in which  points occur
c     nn:  number of points observed in (0,t)
c     zd(i):  occurring data of points in (0,t)
c     xmgd(i):  magnitud data of eqk if available.
c     dep(i):  depths of eqk data if available.
c
c                             * * *
c     ak,c,p,ak2,c2,p2:  coefficients of the modified omori formula;
c                          suffix 2 stands for secondary aftershoks.
c     t0,t1,t2;  t2 is the time period of initiating secondary shocks in
c                  the observed time interval (t0,t1) when the time of
c                  main shock is set to be 0.0.
c     yb:  scaling proportion of y to x in displaying (x,y)
c            coordinates.
c                             * * *
c     a(i):  coefficients of the polynomial trend.
c     xmgmax:  xmgmax = max(xmgd(i); i=1,nn)
c                             * * *
c     kkx:  kkx-1 is the order of the laguerre polynomial for the
c             self-exciting response  function.
c     kky:  kky-1 is the order of the laguerre polynomial for the
c             input response  function.
c     kkc:  kkc is the order of the fourier series
c     kkt:  kky-1 is the order of the polynomial for thetrend.
c     ax(i),ay(i),ac(i),at(i):
c               coefficients of the polynomials corresponding the above;
c               stationary poisson rate is a(1).
c     c:  the exponential coefficient of the self-exciting laguerre
c           response function
c     cm:  the exponential coefficient to measure the scale of the
c            magnitude affection in self-exciting part
c     d:  the exponential coefficient of the input laguerre
c           response function
c     dm:  the exponential coefficient to measure the scale of the
c            magnitude affection in input part
c     cycle:  periodicity
c     yy(i),ymg(i):  input data and their magnitudes
c                             * * *
c     ns:  number of subdivisions of the obsrved interval (0,t).
c     kkx:  number of parameters in the exponential polynomials.
c     a(i):  coefficients of the exponential polynomials.
c                             * * *
c
c     icnt:  option from two graphs count1 or count2.
c          =0  count1, otherwise count2; count1 gives more datail than
c              count2.
c     h:  time length of the moving interval in which points are
c           counted to show the graph count1.
c     isi:  option of displaying the error bands in the log survivor
c             function.
c         =0  displayed; otherwise not.
c     ipoint:  number of subintervals in (0,days) to estimate a non
c                parametric intensity under the palm probability measure
c     days:  length  of interval to display the intensity estimate
c              under the palm probability.
c     delta:  length of a subinterval unit in (0,dmax) to compute the
c               variance time curve.
c     dmax:  time length of a interval to display the variance time
c              curve; this is usually less than t/4.
c
c     dimension z(3000),xmg(3000),fmt(20),dep(3000)
      dimension z(3000),xmg(3000),dep(3000)
      character*8 fmt(10)
      real*8 a(20),zd(3000),xmgd(3000),xx(3000),yy(3000),ymg(3000)
      real*8 ax(20),ay(20),ac(20),at(20)
      real*8 xmgmax
      do 40 i=1,3000
      xmg(i)=0.0
      xmgd(i)=0.0
      dep(i)=0.0
   40 continue
      yb=0.0
      xmgmax=0.0
      xmg1=0.0
      xmg2=10.0
      d1=0.0
      d2=0.0
      xz1=0.0
      xz2=0.0
      open(4,file='pgraph.open')
      nfunct=0
      mgs=1
      iyy=0
      idep=0
c     read(4,4) nfunct,iyy,mgs,idep
      read(4,*) ires,am1
      read(4,*) in
c     read(4,19) fmt
      if(ires.ne.1) open(in,file='work.etas')
      if(ires.eq.1) open(in,file='work.res')
      if(ires.ne.1) read(in,19) fmt
      if(ires.ne.1) write(6,*) fmt
c     read(in,fmt) (xx(i),i=1,nn)
      i=1
   10 continue
      if(ires.ne.1) read(10,*,end=20) idum,e,an,xmgd(i),zd(i),dept    
      if(ires.eq.1) read(10,*,end=20) idum,e,an,xmgd(i),zi,dept,zd(i)    
      i=i+1
      go to 10
   20 n1=i-1
      t=zd(n1)
      nn=0
      do 30 i=1,n1
      if(xmgd(i).lt.am1) go to 30
      if(zd(i).lt.0.0) go to 30
      nn=nn+1
      zd(nn)=zd(i)
      xmgd(nn)=xmgd(i)
   30 continue
c     read(in,1) t,nn,mm
c     read(in,fmt) (zd(i),i=1,nn)
c     if(mgs.ne.0) read(4,fmt) (xmgd(i),i=1,nn)
c     if(idep.ne.0) read(4,fmt) (dep(i),i=1,nn)
c     read(in,*) (zd(i),i=1,nn)
c     if(mgs.ne.0) read(4,*) (xmgd(i),i=1,nn)
c     if(idep.ne.0) read(4,*) (dep(i),i=1,nn)
      read(4,*) icnt,h
      read(4,*) isi
      read(4,*) ipoint,days
      read(4,*) delta,dmax
      if(delta.ne.0.0.and.dmax.ne.0.0) go to 103
      dmax=t/4
      delta=dmax/100
  103 continue
      if(yb.eq.0.0) yb=1.0
      if(nfunct.ne.1) t1=t
      if(nfunct.ne.1) t0=0.0
      ii=0
      do 101 i=1,nn
      if(xmg1.gt.xmgd(i).or.xmg2.lt.xmgd(i)) go to 101
      if(d1.gt.dep(i).or.d2.lt.dep(i)) go to 101
      if(xz2.eq.0) go to 102
      if(xz1.gt.zd(i).or.xz2.lt.zd(i)) go to 101
  102 continue
      ii=ii+1
      z(ii)=zd(i)
      xx(ii)=zd(i)
      zd(ii)=zd(i)
      xmg(ii)=xmgd(i)
      xmgd(ii)=xmgd(i)
      if(xmg(ii).eq.0.0) xmg(ii)=6.0
  101 continue
      nn=ii
c     write(6,6) nfunct,mgs,idep,in,fmt
      write(6,6) nfunct,mgs,idep,in
      if(nfunct.ne.3) write(6,7) t,nn
      if(nfunct.eq.3) write(6,23) t,nn,mm
      write(6,8) (zd(i),i=1,nn)
      if(mgs.ne.0) write(6,8) (xmgd(i),i=1,nn)
      write(6,15) icnt,h
      write(6,16) isi
      write(6,17) ipoint,days
      write(6,18) delta,dmax
    1 format(f10.0,2i10)
    2 format(8f10.0)
    3 format(i10,f10.0)
    4 format(8i10)
    5 format(i10/(8f10.0))
c   6 format(1h ,' funct =',i5,5x,'mgs =',i5,5x,'idep =',i5
c    &       //1h ,' in =',i5//1h ,' fmt = ',20a4)
    6 format(1h ,' funct =',i5,5x,'mgs =',i5,5x,'idep =',i5
     &       //1h ,' in =',i5/)
    7 format(/1h ,' t =',f12.5,5x,' nn =',i10//1h ,' input data ')
    8 format(/1h ,5f12.5/(1h ,5f12.5))
    9 format(/1h ,' ak =',f12.5,3x,' c =',f12.5,3x,' p =',f12.5,3x,
     &             ' ak2 =',f12.5,3x,' c2 =',f12.5,3x,' p2 =',f12.5)
c  10 format(/1h ,' t0 =',f12.5,3x,' t1 =',f12.5,3x,'t2 =',f12.5)
   11 format(/1h ,' yb =',f12.5)
   12 format(/1h ,' a(i),i=1,10'/1h ,10f12.5)
   13 format(/1h ,' kkx =',i5,3x,' kky =',i5,3x,' kkc =',i5,3x,
     &             ' kkt =',i5)
   14 format(/1h ,' xmgmax =',f12.5)
   15 format(/1h ,' icnt =',i5,3x,' h =',f12.5)
   16 format(/1h ,' isi =',i5)
   17 format(/1h ,' ipoint =',i5,3x,' days =',f12.5)
   18 format(/1h ,' delta =',f12.5,3x,'dmax =',f12.5)
   19 format(10a8)
   21 format(//1h ,' ns =',i5,5x,'kkx =',i5//1h ,' a(i),i=1,kkx'/
     &        (1h ,10d13.5))
   22 format(i10/(4d20.14))
   23 format(//1h ,' t =',f12.5,5x,' nn =',i10,5x,' mm =',i10//
     &         1h ,'input data')
   24 format(//1h ,' t0 =',f12.5,5x,' c =',f12.5,5x,' cm =',f12.5,5x,
     &             ' d =',f12.5,5x,' dm =',f12.5)
   25 format(//1h ,10d12.5/(1h ,10d12.5))
      return
      end
      subroutine count1(z,nn,h,nfunct,xmg,xmgmax,xt,inn,tm)
      dimension z(1),xmg(1),x(6000),ix(6000),y(6000),xt(1)
      dimension xnm(2000)
      real*8 xmgmax
      data xnm/2000*0.0/
      i=1
      j=1
      k=0
      n=0
      nmax=0
      xx=0.0
   10 if(z(i)-xx.gt.z(j)+h-xx) go to 20
      k=k+1
      x(k)=z(i)
      n=n+1
      ix(k)=n
      xnm(n+1)=xnm(n+1)+(x(k)-x(k-1))
      if(nmax.lt.n+1) nmax=n+1
      i=i+1
      if(i.gt.nn) go to 30
      go to 10
   20 k=k+1
      xx=z(j)+h
      x(k)=z(j)+h
      n=n-1
      ix(k)=n
      xnm(n+1)=xnm(n+1)+(x(k)-x(k-1))
      j=j+1
      go to 10
   30 continue
      kn=k
    1 format((1h ,10(f7.3,i4,1x)))
      open(7,file='out.pgPTnum')
      write(7,*) h,tm
      do 40 i=1,kn
      y(i)=shimiz(ix(i),h)
      write(7,*) x(i),y(i)
   40 continue
      close(7)
      return
      end
      function shimiz(ix,h)
c   the following transformation provides a normal approximation n(0,1)
c   of poisson variables with mean v; see the equation (50) in
c     r. shimizu(1984). "normal approximation for asymmetric distribu-
c     tions (in japanese)".  proc. inst. statist. math., vol. 32, no.2.
c
      shimiz=(33.d0*ix+29.d0-h-(32.d0*ix+31.d0)*(h/(ix+1.d0))**
     &       (1.d0/4.d0))/(9.d0*sqrt(ix+1.d0))
c
      return
      end
      subroutine surviv(z,n,nfunct,ttt,isi)
c
c     logarithm of cumulative plots with respect to ordered interval
c     lengths. lines of error bounds are given. the differences with
c     the theoretical stationary poisson process are shown in a
c     separated graph with lines of error bounds.
c
      dimension z(1),x(3000),y(3000),w(3000)
      do 100 i=2,n
      x(i-1)=z(i)-z(i-1)
      if(nfunct.eq.0) x(i-1)=x(i-1)*ttt/n
      if(x(i-1).lt.0.0) write(6,2) i,z(i-1),z(i)
    2 format(1h ,' error ',i5,2f15.5)
      if(x(i-1).lt.0.0) x(i-1)=0.0
      w(i-1)=x(i-1)
  100 continue
      n1=n-1
      do 10 i=1,n1-1
      xmin=x(i)
      jj=i
      do 20 j=i+1,n1
      if(x(j).gt.xmin) go to 20
      jj=j
      xmin=x(j)
   20 continue
      x(jj)=x(i)
      x(i)=xmin
   10 continue
      do 30 i=1,n1
      y(i)=n1+1-i
   30 continue
    1 format(1h ,10f13.5)
      do 50 i=1,n1+2
      if(nfunct.eq.0) x(i)=x(i)*n/ttt
   50 continue
      if(isi.eq.0) call errbr2(x,y,n1,xw,n,ttt)
      call errplt(x,n1,xw)
      do 60 i=1,n1+2
      if(nfunct.eq.0) x(i)=x(i)*ttt/n
   60 continue
      call unifrm(x,n1,xw,yw,ttt,w)
      return
      end
      function plsinv(n,k,z,isw)
      real*8 a,b,c,d,z,uu,ud
      a=1.d0-1.d0/(9.d0*(n-k+1))
      b=1.d0-1.d0/(9.d0*k)
      c=1.d0/(9.d0*(n-k+1))
      d=1.d0/(9.d0*k)
      ud=k/(n-k+1.d0)*((a*b-sqrt((a*b)**2-(a**2-c*z**2)*(b**2-d*z**2)))
     &    /(a**2-c*z**2))**3
      uu=k/(n-k+1.d0)*((a*b+sqrt((a*b)**2-(a**2-c*z**2)*(b**2-d*z**2)))
     &    /(a**2-c*z**2))**3
      if(isw.eq.1) plsinv=log(1.d0+ud)
      if(isw.eq.2) plsinv=log(1.d0+uu)
      return
      end
      subroutine errbr2(x,y,n,xw,nx,ttt)
c   this relates to the subroutine errplt.  error bars are drawed by
c   the use of the inversion of the paulson's approximation.
      dimension x(1),y(1)
      dimension xx(3000,6),yy(3000),nc(6)
      real* 8 stderr(6)/0.15866d0,0.84134d0,0.022750d0,0.977250d0,
     &                  0.0013499d0,0.9986501d0/
      do 50 i=1,6
      xx(1,i)=-log(stderr(i))/n
   50 continue
      do 10 k=2,n-1
      xx(k,1)=plsinv(n,k,1.d0,1)
      xx(k,2)=plsinv(n,k,1.d0,2)
      xx(k,3)=plsinv(n,k,2.d0,1)
      xx(k,4)=plsinv(n,k,2.d0,2)
      xx(k,5)=plsinv(n,k,3.d0,1)
      xx(k,6)=plsinv(n,k,3.d0,2)
   10 continue
      do 60 i=1,6
      xx(n,i)=-log(1.d0-stderr(i)**(1.d0/n))
   60 continue
      open(7,file='out.pgSurviv')
      write(7,*) 'surviv'
      do 40 k=1,k
      write(7,*) x(k)*ttt/nx,(xx(k,i)*ttt/nx,i=1,6),n-k+1
   40 continue
      close(7)
      return
      end
      subroutine errplt(x,n,xw)
c   a transformation from order statistic of the exponential random
c   variables, x(k), to the normal n(0,1).  let u(k)=exp(a*x(k)),then
c   we have order statistics of uniform random variables u(k), k=1,n.
c   a transformation is made through the normal approximation of the
c   random variables,  beta(k,n-k+1) = u(k)/(1-u(k)), which is called
c   by name of the paulson's approximation.
      dimension x(1),y(3000)
      plson(n,k,v)=
     &    -((1.-1./(9*k))-((n-k+1.)/k*v)**(1./3.)*(1.-1./(9*(n-k+1))))
     &     /sqrt((1./(9*k))+((n-k+1.)/k*v)**(2./3.)*(1./(9*(n-k+1))))
      open(7,file='out.pgSurDev')
      write(7,*)'errplt'
      do 10 k=1,n
      v=exp(x(k))-1.
      y(k)=plson(n,k,v)
      write(7,*) k,y(k)
   10 continue
      close(7)
    1 format(1h ,10f13.5)
      return 
      end
      subroutine unifrm(x,n,xw,yw,ttt,w)
      dimension x(1),xx(5000),y(5000),w(1),ww(5000)
      rmd=(n+1)/ttt
      open(7,file='out.pgInter1')
      x1=1.35810*sqrt(float(n))/n
      x2=1.62762*sqrt(float(n))/n
      write(7,*) x1,x2
      do 10 i=1,n
      xx(i)=exp(-x(i)*rmd)
      ww(i)=exp(-w(i)*rmd)
   10 continue
      do 30 i=1,n
      y(i)=(n+1.-i)/n
      write(7,*) xx(i),y(i)
   30 continue
      close(7)
      call unitsq(ww,n,xw,yw)
      return
      end
      subroutine unitsq(x,nn,xw,yw)
      dimension x(1),xx(3000),yy(3000)
      n=nn-1
      open(7,file='out.pgInter2')
      write(7,*)'unitsq'
      do 10 i=1,n
      xx(i)=x(i)
      yy(i)=x(i+1)
      write(7,*) xx(i),yy(i)
   10 continue
      close(7)
      return
      end
      subroutine palmpr(x,n,t,t1,n1)
      dimension x(1),xx(1000),p(4)
      character*1 xl(101),xmi,xii,xtt,xst,bl
      dimension xp(4)
      data p/-2.57583,-1.95996,1.95996,2.57583/
      data xmi,xii,xtt,xst,bl/'-','i','|','*',' '/
c     write(6,4)
    4 format(1h )
      do 10 i=1,n1
      xx(i)=0.0
   10 continue
      nc=0
      do 20 i=1,n-1
      if(x(i).gt.t-t1) go to 20
      nc=nc+1
      do 30 j=i+1,n
      if(x(j)-x(i).gt.t1)go to 30
      ii=(x(j)-x(i))*n1/t1+1
      if(ii.lt.1) ii=1
      xx(ii)=xx(ii)+1
   30 continue
   20 continue
      dlt=t1/n1
      rmd=n/t
      rmd1=rmd*dlt*nc
      xp1=rmd1-1./2+2*rmd1*((1-1/(36*rmd1)+p(1)/(6*sqrt(rmd1)))**3-1)
      xp1=xp1/dlt/nc
      xp2=rmd1-1./2+2*rmd1*((1-1/(36*rmd1)+p(2)/(6*sqrt(rmd1)))**3-1)
      xp2=xp2/dlt/nc
      xp3=rmd1-1./2+2*rmd1*((1-1/(36*rmd1)+p(3)/(6*sqrt(rmd1)))**3-1)
      xp3=xp3/dlt/nc
      xp4=rmd1-1./2+2*rmd1*((1-1/(36*rmd1)+p(4)/(6*sqrt(rmd1)))**3-1)
      xp4=xp4/dlt/nc
      open(7,file='out.pgPalm')
      write(7,*) xp1,xp2,xp3,xp4
      do 40 i=1,n1
      xx(i)=xx(i)/dlt/nc
      write(7,*) t1/n1*i,xx(i)
   40 continue
      close(7)
      xmax=2*rmd
      xmin=0.0
      do 55 i=1,n1
      if(xmin.gt.xx(i)) xmin=xx(i)
      if(xmax.lt.xx(i)) xmax=xx(i)
   55 continue
      rmd=rmd*dlt*nc
      do 50 i=1,4
      xp(i)=rmd-1./2.+2*rmd*((1.-1./(36*rmd)+p(i)/(6*sqrt(rmd)))**3-1)
      xp(i)=xp(i)/dlt/nc
      if(xmin.gt.xp(i)) xmin=xp(i)
      if(xmax.lt.xp(i)) xmax=xp(i)
   50 continue
      xm=(xmax+xmin)/2
      x0=0.0
      write(6,1) xmin,xm,xmax
    1 format('time span ',1x,e10.4,13x,e10.4,15x,e10.4)
      do 60 i=1,101
      xl(i)=bl
   60 continue
      xl(1)=xtt
      xl(26)=xtt
      xl(51)=xtt
      write(6,2) (xl(i),i=1,51)
    2 format(13x,51a1)
      do 70 i=1,101
      xl(i)=xmi
   70 continue
      xl(1)=xii
c     write(6,2) (xl(i),i=1,101)
c     write(6,3) x0
c   3 format(1h+,5x,e13.4,'-')
      write(6,3) x0,(xl(i),i=1,51)
    3 format(e12.4,' ',51a1)
      rmd=n/t
      do 80 i=1,n1
      do 90 j=1,101
      xl(j)=bl
   90 continue
      xl(1)=xii
      do 100 j=1,4
      ii=(xp(j)-xmin)/(xmax-xmin)*50+1
      xl(ii)=xtt
  100 continue
      ii=(rmd-xmin)/(xmax-xmin)*50+1
      xl(ii)=xtt
      ii=xx(i)/xmax*50+1
      xl(ii)=xst
      xt1=t1*i/n1
      write(6,5) xt1,(xl(j),j=1,51)
    5 format(d12.4,1x,101a1)
   80 continue
      return
      end
      subroutine vtc(x,n,delta,dmax,t,sigma,k,erres,errest)
c
c   calculation of the variance time curve; see pp. 115-118 in
c   d.r.cox & p.a.w.lewis (1966). the statistical analysis of series
c   of events. methuen, london.
c
c     inputs
c            x:  data
c            n:  number of data
c            t:  length of the obeserved interval.
c            delta:  small interval length no interval length in the
c                      series contain no more then two or three events.
c            dmax:  maximum length og displayed interval
c            sigma:  sum of the products of successive non-dverlappuig
c                      entries in the column.
c            k:  number of culculated variance time curve values.
c
      real*8 a,sigm1,sigm2,ak,avar,amean,r,rmd,xx,rt,sig
      dimension x(3000),a(3000),sigm1(3000),sigm2(3000),sigc0(3000)
     &        ,ak(3000),avar(3000),sigma(3000),xx(3000),amean(3000)
      dimension erres(3000),errest(3000)
      nn=t/delta
      do 10 i=1,nn
   10 xx(i)=0.0
      do 20 i=1,n
      ii=x(i)/delta+1
   20 xx(ii)=xx(ii)+1
      a(1)=nn
      a(2)=nn-1
      sigm1(1)=0.0
      sigm2(1)=0.0
      sigm1(2)=0.0
      sigm2(2)=0.0
      do 30 i=1,nn
      sigm1(1)=sigm1(1)+xx(i)
      sigm2(1)=sigm2(1)+xx(i)**2
      if(i.ne.nn) sigm1(2)=sigm1(2)+xx(i)+xx(i+1)
      if(i.ne.nn) sigm2(2)=sigm2(2)+(xx(i)+xx(i+1))**2
   30 continue
      rmd=n/t
      amean(1)=sigm1(1)/a(1)
      amean(2)=sigm1(2)/a(2)
      sigc0(1)=sigm2(1)-amean(1)*sigm1(1)
      sigc0(2)=sigm2(2)-amean(2)*sigm1(2)
      ak(1)=3*a(1)/(3*a(1)*(a(1)-1)+1-1)
      ak(2)=3*a(2)/(3*a(2)*(a(2)-2)+2**2-1)
      avar(1)=ak(1)*sigc0(1)
      avar(2)=ak(2)*sigc0(2)
      sigma(1)=avar(1)
      erres(1)=1*delta*rmd/(3*a(1))*(4*1**2*delta*rmd+3*1+2*delta*rmd)
      erres(1)=sqrt(erres(1))
      r=1
      rt=rmd*r*delta
      t0=r*delta/t
      errest(1)=(2./3+4./3/r)*rt**2*t0+rt*t0
      errest(1)=sqrt(errest(1))
      k=1
      sigma(2)=avar(2)
      erres(2)=2*delta*rmd/(3*a(2))*(4*2**2*delta*rmd+3*2+2*delta*rmd)
      erres(2)=sqrt(erres(2))
      r=2
      rt=rmd*r*delta
      t0=r*delta/t
      errest(2)=(2./3+4./3/r)*rt**2*t0+rt*t0
      errest(2)=sqrt(errest(2))
      k=2
      k=2
   40 k=k+1
      k4=4*(k-2)
      n4=nn-k4+1
      a(k)=n4
      sigm1(k)=0.0
      sigm2(k)=0.0
      do 50 i=1,n4
      sig=0.0
      do 60 j=1,k4
      ij=i-1+j
      sigm1(k)=sigm1(k)+xx(ij)
      sig=sig+xx(ij)
   60 continue
      sigm2(k)=sigm2(k)+sig**2
   50 continue
      amean(k)=sigm1(k)/a(k)
      sigc0(k)=sigm2(k)-amean(k)*sigm1(k)
      r=k4
      ak(k)=3*a (k)/(3*a(k)*(a(k)-r)+r**2-1)
      avar(k)=ak(k)*sigc0(k)
      sigma(k)=avar(k)
      erres(k)=r*delta*rmd/(3*a(k))*(4*r**2*delta*rmd+3*r+2*delta*rmd)
      erres(k)=sqrt(erres(k))
      rt=rmd*r*delta
      t0=r*delta/t
      errest(k)=(2./3+4./3/r)*rt**2*t0+rt*t0
      errest(k)=sqrt(errest(k))
      if(dmax.ne.0.0.and.dmax.lt.(k4+4)*delta) return
      if((k4+4)*delta.lt.t/4) go to 40
      return
      end
      subroutine vtcprt(sigma,n,delta,nn,t,erres,errest)
      character*1 xl(101),xmi,xii,xtt,xst,bl,xo
      dimension sigma(1),erres(1),errest(1)
      dimension x(2000),y(5)
      data xmi,xii,xtt,xst,bl,xo/'-','|','|','*',' ','o'/
      write(6,4)
    4 format(1h )
      xmax=4*delta*(n-2)
      n1=4*(n-2)+2
      smin=0.0
      smax=0.0
      x(1)=delta
      x(2)=delta*2
      do 20 i=3,n
      x(i)=4*delta*(i-2)
   20 continue
      open(7,file='out.pgVTC')
      write(7,*) 'vtc'
      write(7,*) n,n1
      do 10 i=1,n
      if(smin.gt.sigma(i)) smin=sigma(i)
      if(smax.lt.sigma(i)) smax=sigma(i)
   10 continue
      sm=(smin+smax)/2
      x0=0
      write(6,1) smin,sm,smax
    1 format(19x,e10.4,13x,e10.4,15x,e10.4)
      do 60 i=1,101
      xl(i)=bl
   60 continue
      xl(1)=xtt
      xl(26)=xtt
      xl(51)=xtt
      write(6,2) (xl(i),i=1,51)
    2 format(20x,101a1)
      do 70 i=1,101
      xl(i)=xmi
   70 continue
      xl(1)=xii
c     write(6,2) (xl(i),i=1,101)
c     write(6,3) x0
c   3 format(1h+,5x,e13.4,'-')
      write(6,3) x0,(xl(i),i=1,51)
    3 format(1h ,5x,e13.4,' ',101a1)
      im=(n/t-smin)/(smax-smin)*50+1
      do 80 i=1,n1
      do 90 j=1,101
      xl(j)=bl
   90 continue
      xl(1)=xii
      xl(im)=xtt
      i1=0
      if(i.eq.1) i1=1
      if(i.eq.2) i1=2
      if(mod(i,4).eq.0) i1=i/4+2
c     if(i1.eq.0) go to 85
      if(i1.eq.0) go to 80
      ii=(i*delta*nn/t-smin)/(smax-smin)*50+1
      if(ii.gt.0.and.ii.le.50) xl(ii)=xo
      y(1)=i*delta*nn/t
      do 100 j=1,4
      if(j.eq.1) k=-3
      if(j.eq.2) k=-2
      if(j.eq.3) k=2
      if(j.eq.4) k=3
      ii=(i*delta*nn/t+erres(i1)*k-smin)/(smax-smin)*50+1
      if(ii.gt.0.and.ii.le.100) xl(ii)=xii
      y(j+1)=i*delta*nn/t+erres(i1)*k
  100 continue
      ii=(sigma(i1)-smin)/(smax-smin)*50+1
      xl(ii)=xst
      write(7,*) x(i1),sigma(i1),(y(j),j=1,5)
   85 continue
      if(i.ne.n1) write(6,2) (xl(j),j=1,51)
      if(i.eq.n1) write(6,3) xmax,(xl(j),j=1,51)
   80 continue
      close(7)
c     write(6,3) xmax
      return
      end
      subroutine tmchg0(x,xt,zd,xmgd,tm,ttt,nn,inn,xmin,xmax,t)
      dimension x(3000),xt(300)
      real*8 zd(3000),x1,x2,xmgd(3000)
      do 11 i=1,200
      yi=365.25*i
      if(yi.gt.t) go to 12
      xt(i)=yi*nn/t
   11 continue
   12 inn=i
      tm=nn
      ttt=t
      open(7,file='out.pgCumMT')
      x1=1.35810*sqrt(float(nn))
      x2=1.62762*sqrt(float(nn))
      write(7,*) t,nn
      write(7,*) x1,x2
      do 160 i=1,nn
      x(i)=zd(i)/t*nn
      write(7,*) zd(i),xmgd(i)
  160 continue
      close(7)
      t=nn
      xmin=0.
      xmax=x(nn)
      return
      end
