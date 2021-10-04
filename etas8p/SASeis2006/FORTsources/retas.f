      program etarpp
c-----------------------------------------------------------------------
c Subroutine FUNC4 corresponds to the exact version and FUNC9 to the
c approximate version: see the references below.
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
      common t,nn,mm,iappr,nfunct
c
      call input
c
      do 10 i=2,nn
      if(xx(i).ge.xx(i-1)) go to 10
      write(6,*) 'reverse occurrence time'
      write(6,*) i,xx(i),xx(i-1),xmg(i),xmg(i-1)
   10 continue
      write(6,6) nfunct
      write(6,4) t,nn,mm
      write(6,5) xmag0
      write(6,*)
    6 format(1h ,' funct = ',i5)
    5 format(1h ,'reference magnitudes; xmag0',5x,f10.4)
    4 format(1h ,'(T,nn,mm) =',5x,f10.4,2i6)
    3 format(1h ,10f12.4/(1h ,10f12.4))
    2 format(f10.2,2i10)
    1 format(8f10.2)
c
      call residual
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
c     common/epi/ xp(ldata),yp(ldata)
      common/hyp/ xp(ldata),yp(ldata),dep(ldata)
      common /xyod/xx(ldata),xmg(ldata),xmag0
      common t,nn,mm,iappr,nfunct
      common /range/tstart,ntstar
      character*60 hypodata
      open(unit=1,file='etas.open')
*     open(unit=1,file='./etaspc.open')
*     read(1,112) hypodata
  112 format(a)
      read(1,*) nfunct,iappr
      read(1,*) zts,zte,tstart
      read(1,*) amx1,xmag0
      n=5
      read(1,*) (xini(i),i=1,n)
      close(unit=1)
      write(6,*) ' cutoff-mag/  ref-mag/    t_0  /     T   / dead-time.'
      write(6,5) amx1,xmag0,zts,zte,tstart
    1 format(1h ,20a4)
    4 format(8f10.3)
    5 format(1h ,8f10.3)
    6 format(3i10)
    2 format(f10.0,i10)
    3 format(20a4)
c
      call zisin(t,nd,z,amg,dep,xp,yp,hypodata)
c
      t=zte-zts
      tstart=tstart-zts
      nnn=0
      nn=0
      ntstar=0
      do 10 i=1,nd
      if(amg(i).lt.amx1) go to 10
      if(z(i).lt.zts.or.z(i).gt.zte) go to 10
      nn=nn+1
      if(z(i).lt.tstart) ntstar=nn
      xx(nn)=z(i)-zts
      xmg(nn)=amg(i)
      dep(nn)=dep(i)
      xp(nn)=xp(i)
      yp(nn)=yp(i)
      dep(nn)=dep(i)
   10 continue
      write(6,*) 'input hypocenter data'
      write(6,7) (i,xx(i),xmg(i),i=1,nn)
    7 format(3(i5,f12.5,f4.1))
      mm=nd
      write(6,*) 'read #data; selected #data; #tstart; M_0'
      write(6,8)  nd, nn, ntstar, amx1
    8 format(i10,i16,i9,f6.2)
      if(zte.ne.0.0) t=zte
      t=zte-zts
      return
      end
c***********************************************************************
      subroutine zisin(t,nd,z,amg,dep,xp,yp,hypodata)
c
c            Reading hypocenter data
c
      implicit real * 8 (a-h,o-z)
      parameter(ldata=17777, npara=5)
      dimension z(ldata),amg(ldata),dep(ldata),fmt(20)
      dimension xp(ldata),yp(ldata)
      character*40 hypodata
*     write(6,112) hypodata
  112 format(1h ,a)
*     open(unit=10,file=hypodata)
      open(unit=10,file='./work.etas')
      read(10,1001) fmt
 1001 format(20a4)
      i=1
   10 read(10,*,end=20) idum,xp(i),yp(i),amg(i),z(i),dep(i)
 1002 format(5x,4f12.5,f5.0,5x)
      i=i+1
      go to 10
   20 nd=i-1
      close(unit=10)
      return
      end
c***********************************************************************
      subroutine residual
      implicit real * 8 (a-h,o-z)
      parameter(ldata=17777, npara=5)
      common/param/a(npara),n
      common /range/tstart,ntstar
      common t,nn,mm,iappr,nfunct
      common/xyod/xx(ldata),xmg(ldata),xmag0
c     common/epi/ xp(ldata),yp(ldata)
      common/hyp/ xp(ldata),yp(ldata),dep(ldata)
      dimension x(ldata),xmg0(ldata)
      func41(t,tx,xm,a3,a4)=(log(t-tx+a3)-log(a3))*exp(a4*xm)
      func4p(t,tx,xm,a3,a4,a5)=1.d0/(1.d0-a5)*((t-tx+a3)**(1.d0-a5)
     &      -a3**(1.d0-a5))*exp(a4*xm)
c
      do 40 i=1,nn
   40 xmg0(i)=xmg(i)-xmag0
c
      chtsta=a(1)*tstart
      ft=0.0
      do 30 j=1,ntstar
      if(a(5).eq.1.d0) ft=ft+func41(tstart,xx(j),xmg0(j),a(3),a(4))
      if(a(5).ne.1.d0) ft=ft+func4p(tstart,xx(j),xmg0(j),a(3),a(4),a(5))
   30 continue
      chtsta=chtsta+a(2)*ft
c
      x(1)=a(1)*xx(1)-chtsta
      do 10 i=2,nn
      ft=0.0
      do 20 j=1,i-1
      if(a(5).eq.1.d0) ft=ft+func41(xx(i),xx(j),xmg0(j),a(3),a(4))
      if(a(5).ne.1.d0) ft=ft+func4p(xx(i),xx(j),xmg0(j),a(3),a(4),a(5))
   20 continue
      x(i)=a(1)*xx(i)+a(2)*ft-chtsta
   10 continue
c
c     write(6,1002) (i-ntstar,xmg(i),x(i),x(i)-x(i-1),i=1,nn)
      write(6,1002) (i-ntstar,xmg(i),x(i),x(i)-x(i-1),xx(i),i=1,nn)
      open(unit=1,file='work.res')
*     write(1,1001) (i-ntstar,xmg(i),x(i),i=1,nn)
c     write(1,1003) (i-ntstar,xp(i),yp(i),xmg(i),x(i),i=1,nn)
*     write(1,1003) (i-ntstar,xp(i),yp(i),xmg(i),x(i),xx(i),i=1,nn)
c     write(1,1004) (i-ntstar,xp(i),yp(i),xmg(i),x(i),xx(i),i=1,nn)
      write(1,1005) (i-ntstar,xp(i),yp(i),xmg(i),xx(i),
     &                                         dep(i),x(i),i=1,nn)
c1003 format(i5,4f12.5,5x)
 1003 format(i5,5f12.5)
 1004 format(i6,2f12.5,f6.2,2f15.5)
 1005 format(i5,2f12.5,f12.1,f12.5,f8.2,2x,f12.5)
      close(unit=1)
 1001 format(i5,24x,2f12.5,5x)
c1002 format(i5,24x,3f12.5)
 1002 format(i5,4f12.5)
      return
      end
