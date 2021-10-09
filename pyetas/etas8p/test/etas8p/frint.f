                         
      
        real*8 function frint(func,para,x1,y1,x2,y2,xc,yc)
        
        
        implicit real*8 (a-h,o-z)
        external func

        real*8 para(6)    
        integer id

        id=1

        det=(x1*y2+y1*xc+x2*yc)-(x1*yc+x2*y1+xc*y2)

        if(det.lt.0.0)id=-1
        if(abs(det).lt.1e-10)then
           frint=0
           return
        endif

        r1=dist(x1,y1,xc,yc)
        r2=dist(x2,y2,xc,yc)
        r12=dist(x1,y1,x2,y2)

c        write(*,*)'r12=',r12,x1,y1,x2,y2      

c        print *, (r1*r1+r2*r2-r12*r12)/(2.0*r1*r2)

        theta=(r1*r1+r2*r2-r12*r12)/(2.0*r1*r2)
        if(abs(theta).gt.1d0)theta=1d0-1e-10
        theta=acos(theta)

        if(r1+r2.gt.1e-20)then        
          x0=x1+r1/(r1+r2)*(x2-x1)
          y0=y1+r1/(r1+r2)*(y2-y1)
         else 
          frint=0
          return
         endif

        r0=dist(x0,y0,xc,yc)

c        write(*,*)'frint r=', r1,r0,r2,x0,y0,xc,yc
       
        f1=func(r1,para)
        f2=func(r0,para)
        f3=func(r2,para)
c        write(*,*)'in frint',f1,f2,f3,r1,r2,r0,theta,
c     &      id*(f1/6+f2*2.0/3.0+f3/6.0)*theta
        frint=id*(f1/6+f2*2.0/3.0+f3/6.0)*theta
c        print *,'ff', det,frint,f1,f2,f3,theta
        return
        end


        real*8 function dist(x1,y1,x2,y2)
        
        implicit real*8 (a-h,o-z)
        t1=x1-x2
        t2=y1-y2
        dist=sqrt(t1*t1+t2*t2)
        return
        end

        real*8 function rectint(func,para,xmin,xmax,ymin,ymax,xc,yc)
       
        implicit real*8 (a-h,o-z)
        real*8 para(6)
        external func
         
        ndiv=1000
        
c        write(*,*)'hello', ndiv          

        sum=0
        do i=1, ndiv
          x1=xmin+(xmax-xmin)/ndiv*(i-1)
          y1=ymin
          x2=xmin+(xmax-xmin)/ndiv*i
          y2=ymin      
          sum=sum+frint(func,para,x1,y1,x2,y2,xc,yc)
c         write(*,*)i, 1, sum
c          write(*,*)frint(func,para,x1,y1,x2,y2,xc,yc)
        enddo  
 
c        write(*,*)sum
        do i=1, ndiv
          x1=xmax
          y1=ymin+(ymax-ymin)/ndiv*(i-1)
          x2=xmax
          y2=ymin+(ymax-ymin)/ndiv*i      
c          write(*,*)x1,y1,x2,y2
          sum=sum+frint(func,para,x1,y1,x2,y2,xc,yc)
c         write(*,*)i, 2, sum
        enddo 
c          write(*,*)sum 

        do i=1, ndiv
          x1=xmax+(xmin-xmax)/ndiv*(i-1)
          y1=ymax
          x2=xmax+(xmin-xmax)/ndiv*i
          y2=ymax      
          sum=sum+frint(func,para,x1,y1,x2,y2,xc,yc)
c         write(*,*)i, 3, sum
        enddo  

c        write(*,*)sum

         do i=1, ndiv
          x1=xmin
          y1=ymax+(ymin-ymax)/ndiv*(i-1)
          x2=xmin
          y2=ymax+(ymin-ymax)/ndiv*i      
          sum=sum+frint(func,para,x1,y1,x2,y2,xc,yc) 
c         write(*,*)i, 4, sum    
        enddo  

c          write(*,*)sum
         rectint=sum
         end






        real*8 function polyint(func,para,xp,yp,np, xc,yc)
       
        implicit real*8 (a-h,o-z)
        real*8 para(6), xp(np+1), yp(np+1)
        external func



c        write(*,*)' polyint started'         
        ndiv=1000
        

       
c        write(*,*)'hello', ndiv          

        sum=0


        do j=1, np
          do i=1, ndiv
              x1=xp(j)+(xp(j+1)-xp(j))/ndiv*(i-1)
              y1=yp(j)+(yp(j+1)-yp(j))/ndiv*(i-1)
              x2=xp(j)+(xp(j+1)-xp(j))/ndiv*i
              y2=yp(j)+(yp(j+1)-yp(j))/ndiv*i
              sum=sum+frint(func,para,x1,y1,x2,y2,xc,yc)
c             write(*,*)i, 1, sum
c             write(*,*)frint(func,para,x1,y1,x2,y2,xc,yc)
           enddo  
        enddo
 

c          write(*,*)'polyint finished'
         polyint=sum
	
         end

         
C
C    This programme is edit by Zhaung Jiancang in 2000.6.18 SUN.

c        real*8 function func1(r,para)
c        implicit real*8 (a-h,o-z)
c        real*8 para(6)
c  
c        parameter(pi=3.1415926535897)c
c
c         sigma2=para(1)*exp(para(2)*para(3))
c          sigma2=1
  
c         write(*,*)r,para(1),para(2),para(3)     
c         write(*,*)sigma2,r,exp(-r*r/2/sigma2)
c         write(*,*)'r=',r
c          func1=1/pi/2*(1-exp(-r*r/2/sigma2))
          
c        end  
