C     Last change:  FD   11 Jul 2005    5:19 pm
c       program for o grid
c       declare the variables
      PARAMETER (np1=1000,np2=1000)
      implicit REAL *8 (a-h,o-z)
      dimension x(2,np1,np2),xi(2,np2),a(np1,np2),xo(2,np1,np2)
      COMMON/as/ae(np1,np2),aw(np1,np2),as(np1,np2),an(np1,np2)
     :,ase(np1,np2),ane(np1,np2),asw(np1,np2),anw(np1,np2)
     :,ap(np1,np2)
      dimension alph(np1,np2),beta(np1,np2),gamma(np1,np2),dyde(np1,np2)
      dimension dxi(2),dxdxi(np1,np2),dxde(np1,np2),dydxi(np1,np2),
     :sumx(2),phis(np1,np2),q1(np1,np2),q(2,np1,np2)
      DIMENSION dxix(np1,np2),dxiy(np1,np2),dex(np1,np2),dey(np1,np2),
     :ajac(np1,np2)

      COMMON /bee/n(2)

c-------input.dat ---INPUT
c        open(unit=1,file='c2164.txt',status='unknown')
c-------ogrid3.dat---result file after one iteration
c        open(unit=5,file='ogrid3.dat',status='unknown')

        itsum=0
        w=0.9
	OPEN(1,FILE='cylinder_points.txt',STATUS='unknown')
        read (1,*)n(1),n(2)
        do k=1,2

        do i= 1,n(1)
        do j= 1,n(2)

           q(k,i,j)=0

        end do
        end do

        end do
        dxi(1)=1.0/(n(1)-1)
        dxi(2)=1.0/(n(2)-1)
        itermax=2000
        tol=5e-3
c        read(*,*)itermax

c       scanning xi space
        do k=1,2
         xi(k,1)=0
          do i=2,n(k)
           xi(k,i)=xi(k,i-1)+dxi(k)
          enddo
        enddo

c---------------------------------
c       setting bc
c---------------------------------
      
       do i =1,n(1)
      
        read(1,*)x(1,i,1),x(2,i,1),x(1,i,n(2)),x(2,i,n(2))

   
       enddo
       CLOSE(1)

c-------------------------------
c          setting ic
c-------------------------------

        
        do j=2,n(2)-1

         x(1,1,j)=x(1,1,j-1)+(x(1,1,n(2))-x(1,1,1))/(n(2)-1)
         x(1,n(1),j)=x(1,1,j)
         x(2,1,j)=0.0
         x(2,n(1),j)=x(2,1,j)         



        enddo

        do i=2,n(1)-1
         do j=2,n(2)-1
          do k=1,2
       x(k,i,j)=((j-1)*x(k,i,n(2))+(n(2)-j)*x(k,i,1))/(n(2)-1)

        end do
        enddo
        enddo
       
	

c       x(1,i,j)=x(1,i,1)*xi(2,j)+x(1,i,n(2))*(1-xi(2,j))
c   : +x(1,1,j)*xi(2,i)+x(1,n(1),j)*(1-xi(2,i))
c       x(2,i,j)=x(2,i,1)*xi(2,j)+x(2,i,n(2))*(1-xi(2,j))
c   : +x(2,1,j)*xi(2,i)+x(2,n(1),j)*(1-xi(2,i))

        OPEN(2,FILE='e121i.dat',STATUS='unknown')
        do j =1,n(2)
        do i =1,n(1)

        WRITE(2,*)x(1,i,j),x(2,i,j),'1'
       enddo
       WRITE(2,*)
c       WRITE(2,*)
       enddo
       CLOSE(2)

        OPEN(2,FILE='ini.dat',STATUS='unknown')
       do i =1,n(1)
       do j =1,n(2),n(2)-1
        WRITE(2,*)x(1,i,j),x(2,i,j),'1'
       enddo
       WRITE(2,*)
       enddo
       CLOSE(2)
c	pause

C------------------------------------------
c        SETTING OF BC AND IC COMPLETED
c------------------------------------------
c       setting up the alph,beta,gmma in the interior
c------------------------------------------
c------------------------------------------
c        iteration starts
c------------------------------------------
        call cpu_time(start_time)
        DO III=1,itermax

       do k=1,2
        do i=1,n(1)
         do j=1,n(2)
          xo(k,i,j)=x(k,i,j)
         enddo
       enddo
      enddo

c       do k=1,2

       do i=1,n(1)-1
       do j=1,n(2)-1
c-----------------------------------------
            if (i.eq.1.or.j.eq.1) then

           dxde(i,j)= (x(1,i,j+1)-x(1,i,j))/dxi(2)
           dyde(i,j)= (x(2,i,j+1)-x(2,i,j))/dxi(2)
           dxdxi(i,j)=(x(1,i+1,j)-x(1,i,j))/dxi(1)
           dydxi(i,j)=(x(2,i+1,j)-x(2,i,j))/dxi(1)

             dxde(n(1),j)=dxde(1,j)
             dyde(n(1),j)=dyde(1,j)
             dxde(i,n(2))=dxde(i,1)
             dyde(i,n(2))=dyde(i,1)

             dxdxi(n(1),j)=dxdxi(1,j)
             dydxi(n(1),j)=dxdxi(1,j)
             dxdxi(i,n(2))=dxdxi(i,1)
             dydxi(i,n(2))=dxdxi(i,1)
           else


           dxde(i,j)= (x(1,i,j+1)-x(1,i,j-1))/(2*dxi(2))
           dyde(i,j)= (x(2,i,j+1)-x(2,i,j-1))/(2*dxi(2))
           dxdxi(i,j)=(x(1,i+1,j)-x(1,i-1,j))/(2*dxi(1))
           dydxi(i,j)=(x(2,i+1,j)-x(2,i-1,j))/(2*dxi(1))

          endif
c           write(*,*)dxde(i,j),i.j
        enddo
        enddo



         do i=1,n(1)
         do j=1,n(2)


           alph(i,j)=dxde(i,j)**2+dyde(i,j)**2
           beta(i,j)=dxdxi(i,j)*dxde(i,j)+dydxi(i,j)*dyde(i,j)
           gamma(i,j)=dxdxi(i,j)**2+dydxi(i,j)**2



         enddo
         enddo
C--------SETTING A AND Q MATRIX IN THE INTERIER
c------K LOOP START
        do k=1,2

        do j=2,n(2)-1
        DO i=2,n(1)-1


           q(k,i,j) =0.0
           ae(i,j)  =alph(i,j)/dxi(1)**2
           aw(i,j)  =alph(i,j)/dxi(1)**2
           as(i,j)  =gamma(i,j)/dxi(2)**2
           an(i,j)  =gamma(i,j)/dxi(2)**2
           asw(i,j) =-0.5*beta(i,j)/(dxi(1)*dxi(2))
           anw(i,j) =0.5*beta(i,j)/(dxi(1)*dxi(2))
           ase(i,j) =0.5*beta(i,j)/(dxi(1)*dxi(2))
           ane(i,j) =-0.5*beta(i,j)/(dxi(1)*dxi(2))
           ap(i,j)  =-2*((alph(i,j)/dxi(1)**2)+(gamma(i,j)/dxi(2)**2))

        end do
        END DO
c---------setting row colounm deletion

        DO i=2,n(1)-1
          j=2
c          IF (j.eq.2) THEN
         q(k,i,j)=q(k,i,j)-as(i,j)*x(k,i,j-1)-asw(i,j)*x(k,i-1,j-1)-
     :ase(i,j)*x(k,i+1,j-1)
          as(i,j)=0
         ase(i,j)=0
          asw(i,j)=0
c          END IF

          j=n(2)-1
c          IF(j.eq.n(2)-1)then
         q(k,i,j)=q(k,i,j)-an(i,j)*x(k,i,j+1)-anw(i,j)*x(k,i-1,j+1)-
     :ane(i,j)*x(k,i+1,j+1)
          an(i,j)=0
          ane(i,j)=0
          anw(i,j)=0
c          END IF

        enddo

        do i=1,n(1)
        do j=1,n(2)

         phis(i,j)=x(k,i,j)
         q1(i,j)=q(k,i,j)

        end do
        end do

C---------CALLING SIP MODIFIED-------------
        CALL SIP9P (PHIS,Q1,its)

        do i=1,n(1)
        do j=1,n(2)

         x(k,i,j)=phis(i,j)

        end do
        end do


        itsum=itsum+its

        sumx(k)=0

        do i=1,n(1)
         do j=1,n(2)
c          sumx(k)=(sumx(k)+(x(k,i,j)-xo(k,i,j))**2)/(n(1)*n(2))
	sumx(k)=(sumx(k)+(x(k,i,j)-xo(k,i,j))**2)
          enddo
        enddo
        sumx(k)=sumx(k)/(n(1)*n(2))

        sumx(k)=sumx(k)**0.5

        do i=2,n(1)-1
        do j=2,n(2)-1
          x(k,i,j)=w*x(k,i,j)+(1-w)*xo(k,i,j)
        ENDDO
        enddo

        write(*,*)III,sumx(k),k
        IF(sumx(k).lt.tol)goto 30
          do J=2,n(2)-1
C          X(K,1,J)=(2*X(K,1,J)+(X(K,2,J)+X(K,n(1)-1,J))/2.0)/4.0
          X(1,1,J)=(X(1,2,J)+X(1,n(1)-1,J))/2.0
          X(1,n(1),J)=X(1,1,J)
        end do     
      end do

C-------K LOOP END-----------------------------
      ENDDO
  30  continue
c-----iteration loop endded

c       DO K=1,2

c        do J=2,n(2)-1
c          X(K,1,J)=(X(K,2,J)+X(K,n(1)-1,J))/2.0
c          X(K,n(1),J)=X(K,1,J)
c        end do

c       END DO



 	do i=1,1

       do j=2,n(2)-1
       if (i.eq.1) then
       inn=n(1)-1
       ipp=i+1
       else
       inn=i-1
       ipp=i+1
       end if

       jnn=j-1
       jpp=j+1

       x(1,i,j)=0.25*(x(1,ipp,j)+x(1,inn,j)+x(1,i,jpp)+x(1,i,jnn))

       IF(i.eq.1)x(1,n(1),j)=x(1,i,j)

       end do
       end do

       OPEN(2,FILE='cf_grid.dat',STATUS='unknown')
       	WRITE(2,*) 'zone'
	      WRITE(2,*) 'I=',n(1)
	      WRITE(2,*) 'J=',n(2)
       do j =1,n(2)
       do i =1,n(1)

        WRITE(2,*)x(1,i,j),x(2,i,j),'1'
       enddo
       enddo
       CLOSE(2)

c-----CALCULATION COMPLETED-----------------------


c-----parameter for file input into solver------------------
       OPEN(2,FILE='e121.dat',STATUS='unknown')
C       WRITING N1,N2 AND X,Y

       write(2,*)n(1),n(2)
       WRITE(2,*)dxi(1),dxi(2)
       do i =1,n(1)
         do j =1,n(2)
          WRITE(2,*)x(1,i,j),x(2,i,j)
         enddo
        WRITE(2,*)
       enddo

       do i =1,n(1)
         do j =1,n(2)
           WRITE(2,*)alph(i,j),beta(i,j),gamma(i,j)
        enddo
        WRITE(2,*)
       enddo




       CLOSE(2)



      stop
      end

       subroutine sip9p(phi,q,iter)
      PARAMETER (np1=1000,np2=1000)
      implicit REAL *8 (a-h,o-z)
      COMMON/as/ae(np1,np2),aw(np1,np2),as(np1,np2),an(np1,np2),
     :ase(np1,np2),ane(np1,np2),asw(np1,np2),anw(np1,np2)
     :,ap(np1,np2)
      COMMON /bee/n(2)
      dimension BP(Np1,Np2),BS(Np1,Np2), BN(Np1,Np2), BSE(Np1,Np2),
     :BE(np1,np2), BNE(np1,np2), BW(np1,np2), BSW(np1,np2), BNW(np1,np2)
      DIMENSION RES(np1,np2),Qp(np1,np2),del(np1,np2),phi(np1,np2),
     :q(np1,np2),phio(np1,np2)


C       tol=1e-1
       maxiter=1
        alp=0.92

      do j=1,n(2)
       do i=1,n(1)
      
       bsw(i,j)=0
       bn(i,j)=0
       bs(i,j)=0
       bse(i,j)=0
       bnw(i,j)=0
       bne(i,j)=0
       be(i,j)=0
       bw(i,j)=0
       bp(i,j)=0
       qp(i,j)=0
       del(i,j)=0
      enddo
      enddo


      do j=2,n(2)-1
       do i=2,n(1)-1
      bsw(i,j)=asw(i,j)

      bw(i,j)=(aw(i,j)+alp*anw(i,j)-bsw(i,j)*bn(i-1,j-1))/
     :(1+alp*bn(i-1,j))

      bs(i,j)=(as(i,j)+alp*ase(i,j)-bsw(i,j)*be(i-1,j-1))/
     :(1+alp*be(i,j-1))

      ad=anw(i,j)+ase(i,j)-bs(i,j)*be(i,j-1)-bw(i,j)*bn(i-1,j)
      bp(i,j)=ap(i,j)-alp*ad-bs(i,j)*bn(i,j-1)-bw(i,j)*be(i-1,j)-
     :bsw(i,j)*bne(i-1,j-1)


      bn(i,j)=(an(i,j)+alp*anw(i,j)-alp*bw(i,j)*bn(i-1,j)-
     : bw(i,j)*bne(i-1,j))/bp(i,j)

      be(i,j)=(ae(i,j)+alp*ase(i,j)-alp*bs(i,j)*be(i,j-1)-
     :bs(i,j)*bne(i,j-1))/bp(i,j)

      bne(i,j)=ane(i,j)/bp(i,j)
       end do
      end do

      DO iter=1,MAXITER

      do i=1,n(1)
       do j=1,n(2)

       phio(i,j)=phi(i,j)

       end do
      end do

       SSUM=0.0
       do j=2,n(2)-1
         do i=2,n(1)-1

      res(i,j)=q(i,j)-ap(i,j)*phi(i,j)-ae(i,j)*phi(i+1,j)-an(i,j)*
     :phi(i,j+1)-as(i,j)*phi(i,j-1)-aw(i,j)*phi(i-1,j)-anw(i,j)*
     :phi(i-1,j+1)-ane(i,j)*phi(i+1,j+1)-asw(i,j)*phi(i-1,j-1)
     :-ase(i,j)*phi(i+1,j-1)

       SSUM=ssum+ABS(res(i,j))

      qp(i,j)=(res(i,j)-bs(i,j)*qp(i,j-1)-bw(i,j)*qp(i-1,j)
     :-bsw(i,j)*qp(i-1,j-1))/bp(i,j)

       enddo
      enddo

       IF(iter.eq.1)sumnor=ssum
       sumav=ssum/(sumnor*n(1)*n(2))
c	sumav=ssum/sumnor
c       WRITE(*,*)sumnor

       do j=n(2)-1,2,-1
       do i=n(1)-1,2,-1


       del(i,j)=qp(i,j)-bn(i,j)*del(i,j+1)-be(i,j)*del(i+1,j)-
     :bne(i,j)*del(i+1,j+1)
       phi(i,j)=phi(i,j)+del(i,j)

       enddo
       enddo
c      WRITE(*,*)iter

C        IF(sumav.lt.tol)GOTO 20



       END DO
C   20  continue

       iter=iter-1
      return
      end
