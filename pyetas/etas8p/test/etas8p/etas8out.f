CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   The subroutine outback outputs the background rate on a lattice of
C      mx*my grids
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine outrates(file1)
      implicit real*8 (a-h,o-z)
      character*80 file1
c      include 'mpif.h'
      include 'common.inc'
      real*8 rtx1, rtx2, rty1, rty2
      common /range/ rtx1, rtx2, rty1, rty2

      real*8 w(10)

          xmu=x(1)**2
          a2=x(2)**2
          c=x(3)**2
          alfa=x(4)**2
          p=x(5)**2
          d=x(6)**2
          q=x(7)**2
          gamma=x(8)**2

      if(myrank.eq.0) then
         open(34,file=file1)
         write(34,*)mx
         write(34,*)my
  
         tx1=rtx1
         tx2=rtx2
         ty1=rty1
         ty2=rty2

  
         do i=1,mx
            x0=tx1+(tx2-tx1)/(mx-1)*(i-1)
            do j=1,my
               y0=ty1+(ty2-ty1)/(my-1)*(j-1)
               s=0d0     
               s1=0d0
               do ii=1,nn
                  r0=sqrt((xx(ii)-x0)**2+(yy(ii)-y0)**2)
                  w(1)=zbandw(ii)
                  s=s+zprob(ii)*dgauss(r0,w)
                  s1=s1+dgauss(r0,w)
               enddo
               xlamb = 0.0 + xmu * s/(tz-tstart)
               do ii=1, nn
                 delt=tz - zz(ii)
                 pr1=exp(alfa*zmg(ii))
                 ssig=d*exp(gamma*zmg(ii))
                 bbb=(q-1)/pi/ssig      
                 dist2=(xx(ii)-x0)**2+(yy(ii)-y0)**2
                 pr2=(p-1)/c*(1d0+delt/c)**(-p)
                 pr3=bbb*(dist2/ssig+1d0)**(-q)         
                 xlamb = xlamb + a2*pr1*pr2*pr3
               enddo           
               write(34,991)s/(tz-tstart),s1/(tz-tstart), xlamb
             enddo
         enddo
      endif
      return 
 991  format(1x,3f20.10)
      end


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   The subroutine outprob outputs the probability of each event
C      being a background event or not
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

        subroutine outprob(file1)
        implicit real*8 (a-h,o-z)
        character*80 file1
c        include 'mpif.h'
        include 'common.inc'

        
        if(myrank.eq.2)then
             open(35,file=file1)
             do i=1,nn
                write(35,991)i,zprob(i),zbandw(i),zbkgd(i)
             enddo
        endif

   
       return 

 990  format(6(1x,i5,f7.4))

 991  format(1x,i8,3f20.10)
        end





