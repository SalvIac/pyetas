        subroutine polyse(px,py,n,xx,yy,m,flag)
        implicit real*8 (a-h, o-z)
        double precision px(1:n+1),py(1:n+1),x,y,temp,esp
        double precision xx(1:m),yy(1:m)
        integer te,n,ip,sgn,is,m,flag(1:m)
        esp=1e-8
        px(n+1)=px(1)
        py(n+1)=py(1)
c        write(*,*) n,m
c        do k=1,n
c           write(*,*) px(k), py(k)
c        enddo
        do k=1,m
        x=xx(k)
        y=yy(k)
        ip=0
        te=-1
        do i=1,n
c          write(*,*)py(i),y,py(i+1)
c          write(*,*)'sgn',sgn(y-py(i))*sgn(py(i+1)-y)
          is=sgn(y-py(i))*sgn(py(i+1)-y) 
         if(is.ge.0) then
           if(is.gt.0) then
              temp=px(i+1)*abs(y-py(i))+px(i)*abs(py(i+1)-y)
     &              -abs(py(i+1)-py(i))*x
c                write(*,*) 'temp=',temp
                if(temp.gt.0) then 
                    ip=ip+1
c                    write(*,*)i,'ip=',ip
                  elseif(temp.eq.0) then
                    te=0
                    goto 20
                  else
                endif
            elseif(y.eq.py(i).and.py(i+1).gt.y)  then
               if (x-px(i).lt.0 ) then 
                     ip=ip+1
c                    write(*,*)i,'ip=',ip
               elseif (x-px(i).eq.0) then
                     te=0
                     goto 20
               endif
            elseif((y.eq.py(i+1)).and.py(i).gt.y)  then
               if (x-px(i+1).lt.0 ) then 
                     ip=ip+1
cc                     write(*,*)i,'ip=',ip
               elseif (x-px(i+1).eq.0) then
                     te=0
                     goto 20
               endif
            elseif(y.eq.py(i).and.py(i+1).lt.y)  then
               if (x-px(i).eq.0) then
                     te=0
                     goto 20
               endif
            elseif((y.eq.py(i+1)).and.py(i).lt.y)  then
               if (x-px(i+1).eq.0) then
                     te=0
                     goto 20
               endif
            elseif(y.eq.py(i+1).and.py(i).eq.y)  then
                if(sgn(x-px(i))*sgn(px(i+1)-x).gt.0) then
                     te=0
                     goto 20
                endif
             else
          endif
        endif
        enddo
cc        write(*,*)i,'ip=',ip
        if (ip-ip/2*2.eq.1) then 
           te=1
         else 
           te=-1
        endif
 20     continue
        flag(k)=te
        enddo
        return
        end

        integer function sgn(x)
        double precision x
          if (x .gt. 0.0d0) then 
              sgn=1
            elseif (x .lt. 0.0d0) then
              sgn=-1
            else 
              sgn=0
           endif
          return
         end

        function ycentroid(tx,ty, np)
        implicit real*8 (a-h, o-z)

        real *8 tx(np+1), ty(np+1)
        
        tx(np+1)=tx(1)
        ty(np+1)=ty(1)
        area=0
        ymoment=0

        do i=1, np
          temp=tx(i)*ty(i+1)-tx(i+1)*ty(i)
          area=area+temp
          ymoment=ymoment+(ty(i)+ty(i+1))*temp
        enddo
     
        area=area/2
C        write(*,*) area, ymoment
        ycentroid=ymoment/6d0/area
        return
        end
        
