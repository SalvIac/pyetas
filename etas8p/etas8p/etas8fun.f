
c********************************************************************
c
c    function for computing the background rate at the position of
C    all the events, stored in a common shared vector zbkgd       
c   
c********************************************************************     
 
      subroutine bkgdcalc()
      
      implicit real*8 (a-h,o-z)
      include 'mpif.h'
      include 'common.inc'
      
    
      external pgauss,pcauchy
      real*8 w(6)

c      write(*,*)'background started'
c      write(*,*)myrank,ista2,iend2,'s='
      do i=ista2,iend2
         s=0d0
         do j=1,nn
           r0=sqrt((xx(j)-xx(i))**2+(yy(j)-yy(i))**2)
           w(1)=zbandw(j)
           s=s+zprob(j)*dgauss(r0,w)
c           write(*,*) i, j, r0, w(1), s
         enddo
         zbkgd(i)=s/(tz-tstart)
c         if(zbkgd(i).lt.1d-20)zbkgd(i)=1d-20
c         write(*,*)myrank,i,ista2,iend2,zbandw(i),zprob(i), 's=',s, zbkgd(i)
      enddo
      
c	write(*,*)'passing here'


      do i=0, nprocs-1
        is1=nn*i/nprocs+1
        ie1=nn*(i+1)/nprocs
        call mpi_bcast(zbkgd(is1),ie1-is1+1,mpi_double_precision,i,
     &      mpi_comm_world, ierr)
      enddo 


        s=0d0
        do i=ista2,iend2
          w(1)=zbandw(i) 
c		write(*,*) 'passing here 2' ,i
          s=s+zprob(i)*polyint(pgauss,w,tx,ty,npoly,xx(i),yy(i))
C          s=s+zprob(i)*rectint(pgauss,w,0d0,tx,0d0,ty,xx(i),yy(i))
c          write(*,*) 'passing here 3',i
        enddo 
c 	write(*,*)'passing here 4'
         
          
    
         s1=0
c        write(*,*)myrank,s,s1
         call mpi_reduce(s,s1,1,mpi_double_precision,mpi_sum,0,
     &     mpi_comm_world, ierr)
         call mpi_bcast(s1,1,mpi_double_precision,0,
     &     mpi_comm_world, ierr)
c        write(*,*)myrank,s1
          
        xint0=s1

c        write(*,*)"tz=", tz, "tz-tstart=", tz-tstart
c        write(*,*)"xint0=", xint0
c      write(*,*)'background finished'
      return 
      end
C********************************************************************



C********************************************************************
C   pgauss calculates the integral  of the gaussian kernel function
C     with bandwidth w(1) from 0 to r
C
C********************************************************************   
       real*8 function pcauchy(r,w)

       implicit real*8 (a-h,o-z)
       real*8 w(6)
       parameter(PI=3.14159265358979323846264338327d0)
           
       
       pcauchy=1/2.0d0/PI*(1d0-(1+r**2/w(1)**2)**0.5)
      
       return
       end

 
C********************************************************************
C   dgauss calculates the value of the gaussian kernel function
C     at r with bandwidth w(1)
C
C********************************************************************   
       real*8 function dcauchy(r,w)

       implicit real*8 (a-h,o-z)
       real*8 w(6)
          
       parameter(PI=3.14159265358979323846264338327d0)

       dcauchy=PI/2.5/w(1)**2/(1+r**2/w(1)**2)**1.5

       return
       end

                   

C********************************************************************
C   pgauss calculates the integral  of the gaussian kernel function
C     with bandwidth w(1) from 0 to r
C
C********************************************************************   
       real*8 function pgauss(r,w)

       implicit real*8 (a-h,o-z)
       real*8 w(6)
       parameter(PI=3.14159265358979323846264338327d0)

       pgauss=1/2d0/PI*(1d0-exp(-r**2/2d0/w(1)**2))
      
       return
       end

  
C********************************************************************
C   dgauss calculates the value of the gaussian kernel function
C     at r with bandwidth w(1)
C
C********************************************************************   
       real*8 function dgauss(r,w)

       implicit real*8 (a-h,o-z)
       real*8 w(6)
          
       parameter(PI=3.14159265358979323846264338327d0)

       dgauss=1/2d0/PI/w(1)**2*exp(-r**2/2d0/w(1)**2)
       return
       end

                   

C****************************************************************
C     function to calculate the conditional function value at
C      each the time of each event
c      i: sequence no of events
c      b(7): sqroot of the paratemers
c****************************************************************
      subroutine xlamb(i,b, fv1, g1)
        
      implicit real*8 (a-h,o-z)

      real*8 b(8), g1(50)
      include 'common1.inc'
     
      xmu=b(1)**2
      a2=b(2)**2
      c=b(3)**2
      alfa=b(4)**2
      p=b(5)**2
      d=b(6)**2
      q=b(7)**2
      gamma=b(8)**2
      
      s=xmu*zbkgd(i)
      sg1=zbkgd(i)
      sg2=0d0
      sg3=0d0
      sg4=0d0
      sg5=0d0
      sg6=0d0
      sg7=0d0
      sg8=0d0

      if(i.gt.1)then
       do j=1, i-1
         delt=zz(i)-zz(j)

         pr1=exp(alfa*zmg(j))
         pr2=(p-1)/c*(1d0+delt/c)**(-p)

         ssig=d*exp(gamma*zmg(j))
         bbb=(q-1)/ssig/pi      
         dist2=(xx(i)-xx(j))**2+(yy(i)-yy(j))**2
         pr3=bbb*(dist2/ssig+1d0)**(-q)  
       
         s=s+a2*pr1*pr2*pr3
         sg2=sg2+pr1*pr2*pr3

         pr2_c=pr2*(-1d0/c-p/(c+delt)+p/c)
         sg3=sg3+a2*pr1*pr2_c*pr3

         pr1_alfa=pr1*zmg(j)
         sg4=sg4+a2*pr1_alfa*pr2*pr3

         pr2_p=pr2*(1d0/(p-1)-log(1d0+delt/c))
         sg5=sg5+a2*pr1*pr2_p*pr3

         pr3_d=pr3/d*(-1d0+q*(1d0-1d0/(1d0+dist2/ssig)))
         sg6=sg6+a2*pr1*pr2*pr3_d

         pr3_q= pr3*(1d0/(q-1)-log(1+dist2/ssig))
         sg7=sg7+a2*pr1*pr2*pr3_q
        
         pr3_gamma=pr3*(-zmg(j)+q*zmg(j)*(1d0-1d0/(1d0+dist2/ssig)))
         sg8=sg8+a2*pr1*pr2*pr3_gamma

       enddo
      endif
      fv1=s
      g1(1)=sg1*2*b(1)
      g1(2)=sg2*2*b(2)
      g1(3)=sg3*2*b(3)
      g1(4)=sg4*2*b(4)
      g1(5)=sg5*2*b(5)
      g1(6)=sg6*2*b(6)
      g1(7)=sg7*2*b(7)
      g1(8)=sg8*2*b(8)

      return
      end


   
C  &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
C     Function to calculate integral of the contribution of a
C       previous events 
c          xint(i,b)  
C     
C  *************************************************************
        subroutine  xint(i,b,fv,h)
        implicit real*8 (a-h,o-z)
        include 'common1.inc'
        real*8 b(8),h(8)

        real*8 w(7)      
      
        external fr,dgamma_fr,dd_fr, dq_fr


        xmu=b(1)**2
        a2=b(2)**2
        c=b(3)**2
        alfa=b(4)**2
        p=b(5)**2
        d=b(6)**2
        q=b(7)**2
        gamma=b(8)**2
               
        gi=0d0
        gic=0d0
        gip=0d0
          
        if(zz(i).gt.tstart)then
           ttemp=tz-zz(i)

           gi=1d0-(1d0+ttemp/c)**(1-p)           
           gic=-(1d0-gi)*(1d0-p)*(1d0/(c+ttemp)-1d0/c) 
           gip=-(1d0-gi)*(log(c)-log(c+ttemp))
        else
           ttemp2=tz-zz(i)
           ttemp1=tstart-zz(i)
           gi2=1-(1d0+ttemp2/c)**(1-p)           
           gic2=-(1d0-gi2)*(1d0-p)*(1d0/(c+ttemp2)-1d0/c) 
           gip2=-(1d0-gi2)*(log(c)-log(c+ttemp2))
           gi1=1-(1d0+ttemp1/c)**(1-p)           
           gic1=-(1d0-gi1)*(1d0-p)*(1d0/(c+ttemp1)-1d0/c) 
           gip1=-(1d0-gi1)*(log(c)-log(c+ttemp1))

           gi=gi2-gi1
           gic=gic2-gic1
           gip=gip2-gip1     
        endif

        
        w(1)=gamma
        w(2)=d
        w(3)=q
        w(4)=zmg(i)

         
        si=polyint(fr,w,tx,ty,npoly,xx(i),yy(i))
        sid=polyint(dd_fr,w,tx,ty,npoly,xx(i),yy(i)) 
        siq=polyint(dq_fr,w,tx,ty,npoly,xx(i),yy(i))
        sigamma=polyint(dgamma_fr,w,tx,ty,npoly,xx(i),yy(i))

      
C         if(myrank.eq.0)print*, 'gi',si,'sia',sia,'alfa',sialfa,'d',sid
C        print *, tz, zz(i),ttemp, c, p



        sk=a2*exp(alfa*zmg(i))

        fv=sk*gi*si
      
        h(1)=0d0
        h(2)=sk*gi*si/a2*2d0*b(2)
        h(3)=sk*gic*si*2*b(3)
        h(4)=sk*gi*si*zmg(i)*2*b(4)
        h(5)=sk*gip*si*2*b(5)
        h(6)=sk*gi*sid*2*b(6)              
        h(7)=sk*gi*siq*2*b(7)              
        h(8)=sk*gi*sigamma*2*b(8)              

c        write(*,*)(h(i),i=1,6)
        end


           
C*******************************************************************
C       fr(r)=\int f(r)r dr
c*******************************************************************
              
        real*8 function fr(r,w)

        implicit real*8 (a-h, o-z)

        real*8 w(6)
        parameter(PI=3.14159265358979323846264338327d0)

        gamma=w(1)
        d=w(2)
        q=w(3)
        xmag=w(4)
        
        ssig=d*exp(gamma*xmag)        
        

        fr=(1d0-(1+r*r/ssig)**(1d0-q))/pi/2d0
        
c        print *,'fr=',fr
        return
        end


C*******************************************************************
C       dgamma_fr(r)=d\int f(r)r dr \over  d\alfa
c*******************************************************************
              
        real*8 function dgamma_fr(r,param)

        implicit real*8 (a-h, o-z)
        parameter(PI=3.14159265358979323846264338327d0)
 
        
        real*8 param(6)

         
        gamma=param(1)
        d=param(2)
        q=param(3)
        xmag=param(4)
        
        ssig=d*exp(gamma*xmag)

        dgamma_fr=(-q+1)*(1+r*r/ssig)**(-q)*xmag*r*r/ssig/2d0/pi
        
        return
        end

C*******************************************************************
C       dd_fr(r)=d\int f(r)r dr \over  d\d
c*******************************************************************
              
        real*8 function dd_fr(r,param)

        implicit real*8 (a-h, o-z)
        parameter(PI=3.14159265358979323846264338327d0)

        real*8 param(6)

      
        gamma=param(1)
        d=param(2)
        q=param(3)
        xmag=param(4)
        
        ssig=d*exp(gamma*xmag)

        dd_fr=(-q+1)*(1+r*r/ssig)**(-q)/d*r*r/ssig/2d0/pi
        
      
        return
        end


C*******************************************************************
C       dq_fr(r)=d\int f(r)r dr \over  dq
c*******************************************************************
              
        real*8 function dq_fr(r,param)

        implicit real*8 (a-h, o-z)



        parameter(PI=3.14159265358979323846264338327d0)

        real*8 param(6)

        
        gamma=param(1)
        d=param(2)
        q=param(3)
        xmag=param(4)
       
C        write(*,*) a2, alfa, d, xmag

        ssig=d*exp(gamma*xmag)        

        
        dq_fr=(1d0+r*r/ssig)**(-q+1)*log(1d0+r*r/ssig)/pi/2d0
          
        return
        end











