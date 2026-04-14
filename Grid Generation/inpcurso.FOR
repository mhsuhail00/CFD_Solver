C     Last change:  FD   11 Jul 2005    5:34 pm
      PARAMETER (NP1=1000,NP2=1000)
      implicit REAL *8 (a-h,o-z)
      dimension x(2,NP1,NP2),dxi(2),alph(np1,np2),beta(np1,np2),gamma(np
     :1,np2),n(2),ajac(np1,np2)
      DIMENSION dxdxi(np1,np2),dydxi(np1,np2),dxde(np1,np2),dyde(np1,np2
     :),dex(np1,np2),dey(np1,np2),dxix(np1,np2),dxiy(np1,np2)

      DIMENSION xnox(np1),xnix(np1),xnoy(np1),xniy(np1)

      OPEN(2,FILE='e121.dat',STATUS='old')
      read(2,*)n(1),n(2)
      READ(2,*)dxi(1),dxi(2)
         do i=1,n(1)
          do j=1,n(2)
           read(2,*)x(1,i,j),x(2,i,j)
           enddo
          enddo
      CLOSE(2)
c-----------------------------------------------------
c        reading complete
c-----------------------------------------------------
       do i =2,n(1)-1
       do j =2,n(2)-1

       dxdxi(i,j)=(x(1,i+1,j)-x(1,i-1,j))/(2*dxi(1))
       dydxi(i,j)=(x(2,i+1,j)-x(2,i-1,j))/(2*dxi(1))
       dxde(i,j) =(x(1,i,j+1)-x(1,i,j-1))/(2*dxi(2))
       dyde(i,j) =(x(2,i,j+1)-x(2,i,j-1))/(2*dxi(2))

	enddo
	enddo
	
	i=1
	j=1

       dxdxi(i,j)=(x(1,i+1,j)-x(1,n(1)-1,j))/(2*dxi(1))
       dydxi(i,j)=(x(2,i+1,j)-x(2,n(1)-1,j))/(2*dxi(1))
       dxde(i,j) =(x(1,i,j+1)-x(1,i,j))/dxi(2)
       dyde(i,j) =(x(2,i,j+1)-x(2,i,j))/dxi(2)

	i=1
	do j=2,n(2)-1

       dxdxi(i,j)=(x(1,i+1,j)-x(1,n(1)-1,j))/(2*dxi(1))
       dydxi(i,j)=(x(2,i+1,j)-x(2,n(1)-1,j))/(2*dxi(1))
       dxde(i,j) =(x(1,i,j+1)-x(1,i,j-1))/(2*dxi(2))
       dyde(i,j) =(x(2,i,j+1)-x(2,i,j-1))/(2*dxi(2))

	enddo

	i=1
	j=n(2)

       dxdxi(i,j)=(x(1,i+1,j)-x(1,n(1)-1,j))/(2*dxi(1))
       dydxi(i,j)=(x(2,i+1,j)-x(2,n(1)-1,j))/(2*dxi(1))
       dxde(i,j) =(x(1,i,j)-x(1,i,j-1))/dxi(2)
       dyde(i,j) =(x(2,i,j)-x(2,i,j-1))/dxi(2)
	
	j=n(2)
	do i=2,n(1)-1

       dxdxi(i,j)=(x(1,i+1,j)-x(1,i-1,j))/(2*dxi(1))
       dydxi(i,j)=(x(2,i+1,j)-x(2,i-1,j))/(2*dxi(1))
       dxde(i,j) =(x(1,i,j)-x(1,i,j-1))/dxi(2)
       dyde(i,j) =(x(2,i,j)-x(2,i,j-1))/dxi(2)

	enddo	

	j=1

	do i=2,n(1)-1
	dxdxi(i,j)=(x(1,i+1,j)-x(1,i-1,j))/(2*dxi(1))
        dydxi(i,j)=(x(2,i+1,j)-x(2,i-1,j))/(2*dxi(1))
        dxde(i,j) =(x(1,i,j+1)-x(1,i,j))/dxi(2)
        dyde(i,j) =(x(2,i,j+1)-x(2,i,j))/dxi(2)

	enddo

	i=n(1)
	do j=1,n(2)
	
       dxdxi(i,j)=dxdxi(1,j)
       dydxi(i,j)=dydxi(1,j)
       dxde(i,j) =dxde(1,j) 
       dyde(i,j) =dyde(1,j) 

	enddo
       do i =1,n(1)
       do j =1,n(2)

       ajac(i,j)=(dxdxi(i,j)*dyde(i,j)-dxde(i,j)*dydxi(i,j))

       enddo
       enddo
C---------------------------------------------------------------------------
C     CALCULATING INNER AND OUTER NORMALS (XNIX,XNIY) & (XNOX,XNOY) RESPECTIVELY
C---------------------------------------------------------------------------

      j=1
        do i=1,n(1)

        xnix(i)=x(1,i,j)/sqrt(x(1,i,j)**2+x(2,i,j)**2)
        xniy(i)=x(2,i,j)/sqrt(x(1,i,j)**2+x(2,i,j)**2)

c      WRITE(*,*)XNIX(I),I,XNIY(I)
c      PAUSE

        end do
      j=N(2)
        do i=1,n(1)

        xnox(i)=-x(1,i,j)/sqrt(x(1,i,j)**2+x(2,i,j)**2)
        xnoy(i)=-x(2,i,j)/sqrt(x(1,i,j)**2+x(2,i,j)**2)

c      WRITE(*,*)XNoX(I),I,XNoY(I)
c      PAUSE

        end do


c     Calculating DXI/DX,DXI/DY,DE/DX,DE/DY
c      do i=1,n(1)
c        do j=1,n(2)
c          if(x(1,i,j).eq.20.0)then
c           if(x(2,i,j).eq.5.0)then
c           iur=i
c           jur=j
c          endif
c          if(x(2,i,j).eq.-5.0)then
c          ilr=i
c          jlr=j
c          endif
c        endif
c        if(x(1,i,j).eq.-5.0)then
c          if(x(2,i,j).eq.5.0)then
c          iul=i
c         jul=j
c          endif
c          if(x(2,i,j).eq.-5.0)then
c          ill=i
c          jll=j
c          endif
c        endif
c        enddo
c        enddo
c       WRITE(*,*)iur,jur,ilr,jlr,iul,jul,ill,jll
c       pause
c----------------------------------------------------------------
c      FILE WRITING
C----------------------------------------------------------------
      OPEN(1,FILE='inp_600.dat',STATUS='unknown')

      WRITE(1,*)n(1),n(2),dxi(1),dxi(2)
c      WRITE(1,*)iur,jur,ilr,jlr,iul,jul,ill,jll

      do j =1,n(2)
       do i =1,n(1)

        WRITE(1,*)x(1,i,j),x(2,i,j)

       enddo
       WRITE(1,*)
      enddo


      do j =1,n(2)
       do i =1,n(1)

         dxix(i,j)=dyde(i,j)/ajac(i,j)
         dxiy(i,j)=-dxde(i,j)/ajac(i,j)
         dex(i,j) =-dydxi(i,j)/ajac(i,j)
         dey(i,j) =dxdxi(i,j)/ajac(i,j)
c         (x(2,i,j+1)-x(2,i,j-1))/(2*dxi(2))
c      IF(i.eq.1.and.j.eq.1)WRITE(*,*)i,j,dyde(i,j),x(2,i,j+1),x(2,i,j-1)

      WRITE(1,*)dxix(i,j),dxiy(i,j),dex(i,j),dey(i,j)


       enddo
       WRITE(1,*)
      enddo

      do j =1,n(2)
       do i =1,n(1)

       alph(i,j)=dxix(i,j)**2+dxiy(i,j)**2

       beta(i,j)=dxix(i,j)*dex(i,j)+dxiy(i,j)*dey(i,j)

       gamma(i,j)=dex(i,j)**2+dey(i,j)**2


        WRITE(1,*)alph(i,j),beta(i,j),gamma(i,j)

       enddo
       WRITE(1,*)
      enddo

      do j =1,n(2)
       do i =1,n(1)
        WRITE(1,*)ajac(i,j)
       enddo
      enddo

       do i =1,n(1)
        WRITE(1,*)xnix(i),xniy(i),xnox(i),xnoy(i)

       enddo

       CLOSE(1)

       WRITE(*,*)'done'





       stop
       end


