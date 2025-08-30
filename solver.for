C     Last change:  ALI  21 Feb 2011    3:09 pm
      parameter (np1=350,np2=570)
      implicit REAL *8 (a-h,o-z)
      DIMENSION ae(np1,np2),aw(np1,np2),as(np1,np2),an(np1,np2),
     :ase(np1,np2),ane(np1,np2),asw(np1,np2),anw(np1,np2)
     :,ap(np1,np2)
      COMMON /bee/n(2)

      DIMENSION alph(NP1,NP2),beta(NP1,NP2),gamma(NP1,NP2)
      character filnam(100)*12,resfile*40

      DIMENSION aue(np1,np2),auw(np1,np2),aun(np1,np2),aus(np1,np2),
     :aune(np1,np2),ause(np1,np2),ausw(np1,np2),aunw(np1,np2),
     :aup(np1,np2),ate(np1,np2),atw(np1,np2),atn(np1,np2),ats(np1,np2),
     :atne(np1,np2),atse(np1,np2),atsw(np1,np2),atnw(np1,np2),
     :atp(np1,np2),bus(np1),buse(np1),busw(np1),bts(np1),btse(np1),
     :btsw(np1),bun(np1),bune(np1),bunw(np1),btn(np1),btne(np1)
     :,btnw(np1),
     :aunn(np1,np2),auss(np1,np2),auee(np1,np2),auww(np1,np2),
     :aunnee(np1,np2),aunnww(np1,np2),aussee(np1,np2),aussww(np1,np2),
     :aunne(np1,np2),aunnw(np1,np2),ausse(np1,np2),aussw(np1,np2),
     :aunee(np1,np2),aunww(np1,np2),ausee(np1,np2),ausww(np1,np2),
     :auup(np1,np2),
     :atnn(np1,np2),atss(np1,np2),atee(np1,np2),atww(np1,np2),
     :atnnee(np1,np2),atnnww(np1,np2),atssee(np1,np2),atssww(np1,np2),
     :atnne(np1,np2),atnnw(np1,np2),atsse(np1,np2),atssw(np1,np2),
     :atnee(np1,np2),atnww(np1,np2),atsee(np1,np2),atsww(np1,np2),
     :atup(np1,np2)

      DIMENSION ajac(np1,np2),dxix(np1,np2),dxiy(np1,np2),dex(np1,np2),
     :dey(np1,np2),x(2,np1,np2),dxi(2),q(np1,np2),si(np1,np2),dil(np1,np
     :2),xnox(np1),xnix(np1),xnoy(np1),xniy(np1),xnixi(np1),xnoxi(np1),
     :xniet(np1),xnoet(np1),vr(2,np1),vth(2,np1),qup(np1,np2),
     :qvp(np1,np2),qu(np1,np2),qv(np1,np2),qt(np1,np2),p1(np1,np2),
     :q1(np1,np2),sol(np1,np2)

      DIMENSION u(3,np1,np2),us(2,np1,np2),h(3,np1,np2),pcor(np1,np2),p(
     :np1,np2),d2u(3),conv(3),uxi(np1,np2),uet(np1,np2),up(3,np1,np2),
     :vdotn(np1),uold(3,np1,np2),vort(np1,np2),thi(np1),alc(3)
      REAL*8 Nuss,p_grid,a_grid,ar,aaa,sgn,f_ar

        Ri=0.0

        F=0.0
        
        Pr=0.71
        Pi=ACOS(-1.0)

        thetamax=pi/12

        speed_amp=thetamax*2.0*pi*F  
        accn_amp=2.0*pi*F*speed_amp          
 
        alpha=82.0	!angle from gravity vector
        uinf=SIN(alpha*pi/180.0)
	vinf=COS(alpha*pi/180.0)
	Re=1000
        ubar=0.05
        dt=0.01e-02
        norm=0
	MAXSTEP=5000000
        restart=1
        eps=1e-2
	nsnap=0
        maxsnap=100
	iflag=1

        OPEN(2,FILE='INP.DAT',STATUS='old')
         read(2,*)n(1),n(2),dxi(1),dxi(2)
         read(2,*)p_grid,a_grid,ar 
c         READ(2,*)dxi(1),dxi(2)
          READ(2,*)ic1,ic2,ic3,ic4

         do j=1,n(2)
          do i=1,n(1)
           read(2,*)aaa,bbb,x(1,i,j),x(2,i,j)
          enddo
         enddo

         do j=1,n(2)
          do i=1,n(1)
           read(2,*)dxix(i,j),dxiy(i,j),dex(i,j),dey(i,j)
          enddo
         enddo

         do j=1,n(2)
          do i=1,n(1)
           read(2,*)alph(i,j),beta(i,j),gamma(i,j)
          enddo
         enddo

         do j =1,n(2)
          do i =1,n(1)
           read(2,*)ajac(i,j)
          enddo
         enddo

         do i =1,n(1)
          read(2,*)xnix(i),xniy(i),xnox(i),xnoy(i)
         enddo
          
         do j=1,n(2)
          do i=1,n(1)
c            read(2,*)P1(i,j),Q1(i,j)
            P1(i,j)=0.0
            Q1(i,j)=0.0
          enddo
         enddo
         
       irem=0
       n(2)=n(2)-irem
       if(irem.ne.0)then
        do i=1,n(1)
          xnox(i)=-dex(i,n(2))/sqrt(gamma(i,n(2)))
          xnoy(i)=-dey(i,n(2))/sqrt(gamma(i,n(2)))
        end do
       endif 

c     --------------------------------------------------------
c       generating filenames for saving the snapshots
c     --------------------------------------------------------
        DO I=1,maxsnap

        FILNAM(I)='SNAP000.DAT'

        END DO

        DO K=1,MAXSNAP

        I3=K/100

        I2=(K-100*I3)/10

        I1=K-I2*10-I3*100

        FILNAM(K)(5:5)=char(i3+48)

        FILNAM(K)(6:6)=char(i2+48)

        FILNAM(K)(7:7)=char(i1+48)

        END DO

c--------------------------------------------------------
c     CALCULATING NXi AND Net AT OUTER AND INNER POINTS
C--------------------------------------------------------
c       at inner first
        j=1
        do i=1,n(1)

        xnixi(i)=dxix(i,j)*xnix(i)+dxiy(i,j)*xniy(i)
        xniet(i)=dex(i,j)*xnix(i)+dey(i,j)*xniy(i)

        end do

        j=n(2)
        do i=1,n(1)

        xnoxi(i)=dxix(i,j)*xnox(i)+dxiy(i,j)*xnoy(i)
        xnoet(i)=dex(i,j)*xnox(i)+dey(i,j)*xnoy(i)

        end do

c=============================
c	n(2)=245
c	dxi(2)=1.0/(n(2)-1)
c=============================

c-------------------------------------------------------------------
       open (1,FILE='bound.dat',status='unknown')

       do j=1,n(2),n(2)-1
       do i=1,n(1)

       WRITE(1,*)i,j,x(1,i,j),x(2,i,j),'1'

       end do

       WRITE(1,*)

       end do
       CLOSE(1)

c---------------------------------------
c         Aplying Initial conditions
c---------------------------------------
       if (restart.eq.0) then
       loop=1
       time=0
          do j=1,n(2)
          do i=1,n(1)

          U(1,i,j)=uinf
          U(2,i,j)=vinf
          U(3,i,j)=0.0
          
          Uxi(i,j)=0
          Uet(i,j)=0
          P(i,j)=0
          Up(1,i,j)=uinf
          up(2,i,j)=vinf
          pcor(i,j)=0
          si(i,j)=0
          enddo
          enddo

       else

       open (1,FILE='spa100.dat',status='unknown',form='unformatted')
       read (1)loop,time,dmax

       read(1)x,si,u,p

       CLOSE(1)
       endif

	iiflag=0
	iflag=0
	t_period=100.0
	if (iflag.eq.1) then 

	icycles= time/t_period
	tstart=(icycles+1)*t_period
	t_incr=t_period/maxsnap
	i_loop=t_incr/dt
	loop_snap=loop+(tstart-time)/dt
	iflag=0
	iiflag=1	
	nsnap=1
	write (*,*)tstart,time,loop_snap,i_loop,loop  

        endif
	
c--------------------------------------------------------
c     CALCULATING V dot N at outer to determine inward or outward flow
c     IF V DOT N IS POSITIVE INWARD AND OUTWARD FLOW and storing indices in iino,iinl
C--------------------------------------------------------
c         kk=0
c         j=n(2)
c         do i=1,n(1)
c         vdotn(i)=u(1,i,j)*xnox(i)+u(2,i,j)*xnoy(i)
c
c         if (vdotn(i).gt.0) then
c          kk=kk+1
c           if (kk.eq.1) then
c           iino=i
c          end if
c         end if
c         end do

c        iinl=iino+kk-1

c----------------------------------------------------
c       APPLYING BOUNDARY CONDITION
c---------setting boundary conditions----------------
c---------solid-fluid boundary
          j=1
          do k=1,2
          do i=1,n(1)

           IF(k.eq.1)then
            U(k,i,j)=-Speed_amp*x(2,i,j)
            up(k,i,j)=U(k,i,j)
           else
            U(k,i,j)=Speed_amp*x(1,i,j)
            up(k,i,j)=U(k,i,j)
           endif

          end do
          end do

          j=1
          do i=1,n(1)
            U(3,i,j)=1.0
          end do

c----------------------------------------------------
c       setting bc at infinity
c----------------------------------------------------
        j=n(2)

        do i=1,n(1)-1

        vnn=u(1,i,j)*xnox(i)+u(2,i,j)*xnoy(i)

        IF(vnn.ge.0)then

c        inflow dirichlet conditions

           u(1,i,j)=uinf

           u(2,i,j)=vinf
           u(3,i,j)=0.0
           up(1,i,j)=u(1,i,j)
           up(2,i,j)=u(2,i,j)

        else

c---------------------Neuman condition---------------
          inn=i-1
          ipp=i+1

         if (i.eq.1) then
          inn=n(1)-1
         endif
         jnn=j-1

        U(1,i,j)=u(1,i,jnn)
        U(2,i,j)=u(2,i,jnn)
        U(3,i,j)=u(3,i,jnn)

        if (i.eq.1) then
         u(1,n(1),j)=u(1,i,j)
         u(2,n(1),j)=u(2,i,j)
         u(3,n(1),j)=u(3,i,j)
        end if

       endif

      end do


c     ------------------------------------------------------
c       forming coeff matrix for velocity
c     ------------------------------------------------------
      do j=2,n(2)-1
       do i=1,n(1)-1

       if (i.eq.1) then
        inn=n(1)-1
        ipp=i+1
        else
        inn=i-1
        ipp=i+1
        end if
        jpp=j+1
        jnn=j-1

        if (j.eq.2.or.j.eq.n(2)-1) then

        aue(i,j)=-dt*(alph(i,j)/(dxi(1)**2)+P1(i,j)/(2.0*dxi(1)))/Re
        auw(i,j)=-dt*(alph(i,j)/(dxi(1)**2)-P1(i,j)/(2.0*dxi(1)))/Re
        aun(i,j)=-dt*(gamma(i,j)/(dxi(2)**2)+Q1(i,j)/(2.0*dxi(2)))/Re
        aus(i,j)=-dt*(gamma(i,j)/(dxi(2)**2)-Q1(i,j)/(2.0*dxi(2)))/Re

        aune(i,j)=dt*beta(I,J)/(2.0*dxi(1)*dxi(2)*Re)
	ausw(i,j)=dt*beta(I,J)/(2.0*dxi(1)*dxi(2)*Re)
        aunw(i,j)=-dt*beta(I,J)/(2.0*dxi(1)*dxi(2)*Re)
        ause(i,j)=-dt*beta(I,J)/(2.0*dxi(1)*dxi(2)*Re)
        aup(i,j)=1+dt*2.0*(alph(i,j)/(dxi(1)**2)+gamma(i,j)/
     :  (dxi(2)**2))/Re

c       -----------------------------------------------------
c       forming coeff matrix for temperature
c       -----------------------------------------------------
        ate(i,j)=-dt*(alph(i,j)/(dxi(1)**2)+P1(i,j)/(2.0*dxi(1)))/
     :           (Re*Pr)
        atw(i,j)=-dt*(alph(i,j)/(dxi(1)**2)-P1(i,j)/(2.0*dxi(1)))/
     :           (Re*Pr)
        atn(i,j)=-dt*(gamma(i,j)/(dxi(2)**2)+Q1(i,j)/(2.0*dxi(2)))/
     :           (Re*Pr)
        ats(i,j)=-dt*(gamma(i,j)/(dxi(2)**2)-Q1(i,j)/(2.0*dxi(2)))/
     :           (Re*Pr)

        atne(i,j)=dt*(beta(I,J)/(2.0*dxi(1)*dxi(2)))/(Re*Pr)
        atsw(i,j)=dt*(beta(I,J)/(2.0*dxi(1)*dxi(2)))/(Re*Pr)
        atnw(i,j)=-dt*(beta(I,J)/(2.0*dxi(1)*dxi(2)))/(Re*Pr)
        atse(i,j)=-dt*(beta(I,J)/(2.0*dxi(1)*dxi(2)))/(Re*Pr)
        atp(i,j)=1+dt*2.0*(alph(i,j)/(dxi(1)**2)+gamma(i,j)/
     :  (dxi(2)**2))/(Re*Pr)

C       ----------------------------------------------------
C       FORTH ORDER COEFFICIENT MATRIX FOR VELOCITY
c       ----------------------------------------------------
        else
        aue(i,j)=(-dt)*((4.0*alph(i,j))/(3.0*(dxi(1)**2))+
     :  (2.0*P1(i,j))/(3.0*dxi(1)))/Re
        auw(i,j)=(-dt)*((4.0*alph(i,j))/(3.0*(dxi(1)**2))-
     :  (2.0*P1(i,j))/(3.0*dxi(1)))/Re
        aun(i,j)=(-dt)*((4.0*gamma(i,j))/(3.0*(dxi(2)**2))+
     :  (2.0*Q1(i,j))/(3*dxi(2)))/Re
        aus(i,j)=(-dt)*((4.0*gamma(i,j))/(3.0*(dxi(2)**2))-
     :  (2.0*Q1(i,j))/(3*dxi(2)))/Re
        aune(i,j)=(-dt)*(-8.0*beta(I,J)/(9.0*dxi(1)*dxi(2)))/Re
	aunw(i,j)=(-dt)*(8.0*beta(I,J)/(9.0*dxi(1)*dxi(2)))/Re
        ause(i,j)=(-dt)*(8.0*beta(I,J)/(9.0*dxi(1)*dxi(2)))/Re
        ausw(i,j)=(-dt)*(-8.0*beta(I,J)/(9.0*dxi(1)*dxi(2)))/Re

        aunn(i,j)=(-dt)*(-gamma(i,j)/(12.0*(dxi(2)**2))-
     :  Q1(i,j)/(12.0*dxi(2)))/Re
        auss(i,j)=(-dt)*(-gamma(i,j)/(12.0*(dxi(2)**2))+
     :  Q1(i,j)/(12.0*dxi(2)))/Re
        auee(i,j)=(-dt)*(-alph(i,j)/(12.0*(dxi(1)**2))-
     :  P1(i,j)/(12.0*dxi(1)))/Re
        auww(i,j)=(-dt)*(-alph(i,j)/(12.0*(dxi(1)**2))+
     :  P1(i,j)/(12.0*dxi(1)))/Re
        aunnee(i,j)=(-dt)*(-beta(I,J)/(72.0*dxi(1)*dxi(2)))/Re
        aunnww(i,j)=(-dt)*(beta(I,J)/(72.0*dxi(1)*dxi(2)))/Re
        aussee(i,j)=(-dt)*(beta(I,J)/(72.0*dxi(1)*dxi(2)))/Re
        aussww(i,j)=(-dt)*(-beta(I,J)/(72.0*dxi(1)*dxi(2)))/Re
        aunne(i,j)=(-dt)*(beta(I,J)/(9.0*dxi(1)*dxi(2)))/Re
        aunnw(i,j)=(-dt)*(-beta(I,J)/(9.0*dxi(1)*dxi(2)))/Re
        ausse(i,j)=(-dt)*(-beta(I,J)/(9.0*dxi(1)*dxi(2)))/Re
        aussw(i,j)=(-dt)*(beta(I,J)/(9.0*dxi(1)*dxi(2)))/Re
        aunee(i,j)=(-dt)*(beta(I,J)/(9.0*dxi(1)*dxi(2)))/Re
        aunww(i,j)=(-dt)*(-beta(I,J)/(9.0*dxi(1)*dxi(2)))/Re
        ausee(i,j)=(-dt)*(-beta(I,J)/(9.0*dxi(1)*dxi(2)))/Re
        ausww(i,j)=(-dt)*(beta(I,J)/(9.0*dxi(1)*dxi(2)))/Re

        aup(i,j)=1+dt*(5.0*alph(i,j)/(2.0*(dxi(1)**2))+
     :  5.0*gamma(i,j)/(2.0*(dxi(2)**2)))/Re

C       ----------------------------------------------------
C       FORTH ORDER COEFFICIENT MATRIX FOR TEMPERATURE
c       ----------------------------------------------------
        ate(i,j)=aue(i,j)/Pr
        atw(i,j)=auw(i,j)/Pr
        atn(i,j)=aun(i,j)/Pr
        ats(i,j)=aus(i,j)/Pr
        atne(i,j)=aune(i,j)/Pr
	atnw(i,j)=aunw(i,j)/Pr
        atse(i,j)=ause(i,j)/Pr
        atsw(i,j)=ausw(i,j)/Pr

        atnn(i,j)=aunn(i,j)/Pr
        atss(i,j)=auss(i,j)/Pr
        atee(i,j)=auee(i,j)/Pr
        atww(i,j)=auww(i,j)/Pr
        atnnee(i,j)=aunnee(i,j)/Pr
        atnnww(i,j)=aunnww(i,j)/Pr
        atssee(i,j)=aussee(i,j)/Pr
        atssww(i,j)=aussww(i,j)/Pr
        atnne(i,j)=aunne(i,j)/Pr
        atnnw(i,j)=aunnw(i,j)/Pr
        atsse(i,j)=ausse(i,j)/Pr
        atssw(i,j)=aussw(i,j)/Pr
        atnee(i,j)=aunee(i,j)/Pr
        atnww(i,j)=aunww(i,j)/Pr
        atsee(i,j)=ausee(i,j)/Pr
        atsww(i,j)=ausww(i,j)/Pr

        atp(i,j)=1+dt*(5.0*alph(i,j)/(2.0*(dxi(1)**2))+
     :  5.0*gamma(i,j)/(2.0*(dxi(2)**2)))/(Re*Pr)

        end if

       if (j.eq.2) then
       bus(i)=aus(i,j)
       buse(i)=ause(i,j)
       busw(i)=ausw(i,j)
       bts(i)=ats(i,j)
       btse(i)=atse(i,j)
       btsw(i)=atsw(i,j)

       aus(i,j)=0
       ause(i,j)=0
       ausw(i,j)=0
       ats(i,j)=0
       atse(i,j)=0
       atsw(i,j)=0

       end if

       if (j.eq.n(2)-1) then
       bun(i)=aun(i,j)
       bune(i)=aune(i,j)
       bunw(i)=aunw(i,j)
       btn(i)=atn(i,j)
       btne(i)=atne(i,j)
       btnw(i)=atnw(i,j)

       aun(i,j)=0
       aune(i,j)=0
       aunw(i,j)=0
       atn(i,j)=0
       atne(i,j)=0
       atnw(i,j)=0

       end if

       if (i.eq.1) then
       aue(n(1),j) =aue(i,j)
       auw(n(1),j) =auw(i,j)
       aun(n(1),j) =aun(i,j)
       aus(n(1),j) =aus(i,j)
       aune(n(1),j)=aune(i,j)
       ause(n(1),j)=ause(i,j)
       ausw(n(1),j)=ausw(i,j)
       aunw(n(1),j)=aunw(i,j)
       aup(n(1),j) =aup(i,j)

       aunn(n(1),j)=aunn(i,j)
       aunnee(n(1),j)=aunnee(i,j)
       aunnww(n(1),j)=aunnww(i,j)
       aunne(n(1),j)=aunne(i,j)
       aunnw(n(1),j)=aunnw(i,j)
       aunee(n(1),j)=aunee(i,j)
       aunww(n(1),j)=aunww(i,j)
       auss(n(1),j)=auss(i,j)
       aussee(n(1),j)=aussee(i,j)
       aussww(n(1),j)=aussww(i,j)
       ausse(n(1),j)=ausse(i,j)
       aussw(n(1),j)=aussw(i,j)
       ausee(n(1),j)=ausee(i,j)
       ausww(n(1),j)=ausww(i,j)
       auee(n(1),j)=auee(i,j)
       auww(n(1),j)=auww(i,j)

       ate(n(1),j) =ate(i,j)
       atw(n(1),j) =atw(i,j)
       atn(n(1),j) =atn(i,j)
       ats(n(1),j) =ats(i,j)
       atne(n(1),j)=atne(i,j)
       atse(n(1),j)=atse(i,j)
       atsw(n(1),j)=atsw(i,j)
       atnw(n(1),j)=atnw(i,j)
       atp(n(1),j) =atp(i,j)

       atnn(n(1),j)=atnn(i,j)
       atnnee(n(1),j)=atnnee(i,j)
       atnnww(n(1),j)=atnnww(i,j)
       atnne(n(1),j)=atnne(i,j)
       atnnw(n(1),j)=atnnw(i,j)
       atnee(n(1),j)=atnee(i,j)
       atnww(n(1),j)=atnww(i,j)
       atss(n(1),j)=atss(i,j)
       atssee(n(1),j)=atssee(i,j)
       atssww(n(1),j)=atssww(i,j)
       atsse(n(1),j)=atsse(i,j)
       atssw(n(1),j)=atssw(i,j)
       atsee(n(1),j)=atsee(i,j)
       atsww(n(1),j)=atsww(i,j)
       atee(n(1),j)=atee(i,j)
       atww(n(1),j)=atww(i,j)
       end if

      end do
      end do

c-------------------------------------------------------------
c-----Forming A matrix for pressure poisson Equation
c-------------------------------------------------------------
      do j=2,n(2)-1
       do i=1,n(1)-1

        if (i.eq.1) then
        inn=n(1)-1
        ipp=i+1
        else
        inn=i-1
        ipp=i+1
        end if
        jpp=j+1
        jnn=j-1

c         EAST COMPONENT(I+1,J)
       aae=(dxix(i,j)/(2.*dxi(1)**2))*(dxix(i,j)+dxix(ipp,j))
       bbe=(dex(i,j)/(8.*dxi(1)*dxi(2)))*(dxix(i,jpp)-dxix(i,jnn))
       cce=(dxiy(i,j)/(2.*dxi(1)**2))*(dxiy(i,j)+dxiy(ipp,j))
       dde=(dey(i,j)/(8.*dxi(1)*dxi(2)))*(dxiy(i,jpp)-dxiy(i,jnn))

       ae(i,j)=aae+bbe+cce+dde

c         WEST COMPONENT(I-1,J)
       aaw=(dxix(i,j)/(2.*dxi(1)**2))*(dxix(i,j)+dxix(inn,j))
       bbw=(dex(i,j)/(8.*dxi(1)*dxi(2)))*(dxix(i,jnn)-dxix(i,jpp))
       ccw=(dxiy(i,j)/(2.*dxi(1)**2))*(dxiy(i,j)+dxiy(inn,j))
       ddw=(dey(i,j)/(8.*dxi(1)*dxi(2)))*(dxiy(i,jnn)-dxiy(i,jpp))





       aw(i,j)=aaw+bbw+ccw+ddw

c        NORTH COMPONENT(I,J+1)
       aan=(dxix(i,j)/(8.*dxi(1)*dxi(2)))*(dex(ipp,j)-dex(inn,j))
       bbn=(dex(i,j)/(2.*dxi(2)**2))*(dex(i,j)+dex(i,jpp))
       ccn=(dxiy(i,j)/(8.*dxi(1)*dxi(2)))*(dey(ipp,j)-dey(inn,j))
       ddn=(dey(i,j)/(2.*dxi(2)**2))*(dey(i,j)+dey(i,jpp))

       an(i,j)=aan+bbn+ccn+ddn

c        SOUTH COMPONENNT(I,J-1)
       aas=(dxix(i,j)/(8.*dxi(1)*dxi(2)))*(dex(inn,j)-dex(ipp,j))
       bbs=(dex(i,j)/(2.*dxi(2)**2))*(dex(i,j)+dex(i,jnn))
       ccs=(dxiy(i,j)/(8.*dxi(1)*dxi(2)))*(dey(inn,j)-dey(ipp,j))
       dds=(dey(i,j)/(2.*dxi(2)**2))*(dey(i,j)+dey(i,jnn))

       as(i,j)=aas+bbs+ccs+dds

c       NORTH EAST COMPONENT(I+1,J+1)
        aane=(dxix(i,j)/(8.*dxi(1)*dxi(2)))*(dex(i,j)+dex(ipp,j))
        bbne=(dex(i,j)/(8.*dxi(1)*dxi(2)))*(dxix(i,j)+dxix(i,jpp))
        ccne=(dxiy(i,j)/(8.*dxi(1)*dxi(2)))*(dey(i,j)+dey(ipp,j))
        ddne=(dey(i,j)/(8.*dxi(1)*dxi(2)))*(dxiy(i,j)+dxiy(i,jpp))

        ane(i,j)=aane+bbne+ccne+ddne

c       SOUTH WEST COMPONENT(I-1,J-1)
        aasw=(dxix(i,j)/(8.*dxi(1)*dxi(2)))*(dex(i,j)+dex(inn,j))
        bbsw=(dex(i,j)/(8.*dxi(1)*dxi(2)))*(dxix(i,j)+dxix(i,jnn))
        ccsw=(dxiy(i,j)/(8.*dxi(1)*dxi(2)))*(dey(i,j)+dey(inn,j))
        ddsw=(dey(i,j)/(8.*dxi(1)*dxi(2)))*(dxiy(i,j)+dxiy(i,jnn))

        asw(i,j)=aasw+bbsw+ccsw+ddsw

c       NORTH WEST(I-1,J+1)
        aanw=-(dxix(i,j)/(8.*dxi(1)*dxi(2)))*(dex(i,j)+dex(inn,j))
        bbnw=-(dex(i,j)/(8.*dxi(1)*dxi(2)))*(dxix(i,j)+dxix(i,jpp))
        ccnw=-(dxiy(i,j)/(8.*dxi(1)*dxi(2)))*(dey(i,j)+dey(inn,j))
        ddnw=-(dey(i,j)/(8.*dxi(1)*dxi(2)))*(dxiy(i,j)+dxiy(i,jpp))

        anw(i,j)=aanw+bbnw+ccnw+ddnw

c       SOUTH EAST COMPONENTS(I+1,J-1)
        aase=-(dxix(i,j)/(8.*dxi(1)*dxi(2)))*(dex(i,j)+dex(ipp,j))
        bbse=-(dex(i,j)/(8.*dxi(1)*dxi(2)))*(dxix(i,j)+dxix(i,jnn))
        ccse=-(dxiy(i,j)/(8.*dxi(1)*dxi(2)))*(dey(i,j)+dey(ipp,j))
        ddse=-(dey(i,j)/(8.*dxi(1)*dxi(2)))*(dxiy(i,j)+dxiy(i,jnn))

        ase(i,j)=aase+bbse+ccse+ddse

c       node itself P
      pxi=(1.0/(2.*dxi(1)**2))
      pet=(1.0/(2.*dxi(2)**2))
      aap=-dxix(i,j)*(2.*dxix(i,j)+dxix(inn,j)+dxix(ipp,j))
      bbp=-dex(i,j)*(2.*dex(i,j)+dex(i,jnn)+dex(i,jpp))
      ccp=-dxiy(i,j)*(2.*dxiy(i,j)+dxiy(inn,j)+dxiy(ipp,j))
      ddp=-dey(i,j)*(2.*dey(i,j)+dey(i,jnn)+dey(i,jpp))

      ap(i,j)=(aap*pxi+bbp*pet+ccp*pxi+ddp*pet)

      if (i.eq.1) then

      ae(n(1),j) =ae(i,j)
      aw(n(1),j) =aw(i,j)
      an(n(1),j) =an(i,j)
      as(n(1),j) =as(i,j)
      ane(n(1),j)=ane(i,j)
      ase(n(1),j)=ase(i,j)
      asw(n(1),j)=asw(i,j)
      anw(n(1),j)=anw(i,j)
      ap(n(1),j) =ap(i,j)

      end if

      end do
      end do

C----------------------------------------------------------
c     START OF TIME LOOP
C----------------------------------------------------------
      DO loop=loop,maxstep

      time=time+dt
c      pause
C--------FLOW FIELD INSIDE DOMAIN
c      U in xi and eta
        do i=1,n(1)
        do j=1,n(2)

        uxi(i,j)=dxix(i,j)*U(1,i,j)+dxiy(i,j)*u(2,i,j)
	uet(i,j)=dex(i,j)*U(1,i,j)+dey(i,j)*u(2,i,j)
        uold(3,i,j)=u(3,i,j)

        enddo
        enddo


C----------------------------------------------------------
C     CONVECTION TERM
C----------------------------------------------------------
C         K LOOP STARTS
        do j=2,n(2)-1
        do i=1,n(1)-1
      if (i.eq.1.or.i.eq.2.or.i.eq.(n(1)-1)) then

	if (i.eq.1)then
        inn=n(1)-1
        ipp=i+1
        inn2=n(1)-2
        ipp2=i+2
	endif

	if (i.eq.2)then
        inn=i-1
        ipp=i+1
        inn2=n(1)-1
        ipp2=i+2
	endif

	if (i.eq.n(1)-1)then
        inn=i-1
        ipp=i+1
        inn2=i-2
        ipp2=2
	endif

      else

        inn=i-1
        ipp=i+1
        inn2=i-2
        ipp2=i+2

      end if

	jpp=j+1
        jnn=j-1
        jpp2=j+2
        jnn2=j-2

       DO K=1,3

        IF(k.le.2)then
	pec1=uxi(i,j)*re*dxi(1)/alph(i,j)

	pec2=uet(i,j)*re*dxi(2)/gamma(i,j)
        else
	pec1=uxi(i,j)*re*Pr*dxi(1)/alph(i,j)

	pec2=uet(i,j)*re*Pr*dxi(2)/gamma(i,j)
        endif

c======================================================================
c     CONVECTIVE TERM -THIRD ORDER ASYMMETRIC UPWIND DIFFERENCING IN
C     CENTER AND CENTRAL AT BOUNDARY + HYBRID DIFFERENCING
C======================================================================	
	if (j.ge.3.and.j.le.n(2)-2) then

	IF (pec1.le.2.and.pec1.gt.-2) THEN	   
C	CENTRAL 4TH ORDER

       xpp=8.0*(u(k,ipp,j)-u(k,inn,j))
       xnn=u(k,ipp2,j)-u(k,inn2,j)

       du_xi=(1.0/12.0)*(xpp-xnn)/dxi(1)

	ELSE

C	UPWIND 3RD ORDER	
        ak1=uxi(i,j)*(-u(k,ipp2,j)+8*u(k,ipp,j)-8*u(k,inn,j)
     :                 +u(k,inn2,j))/(12.0*dxi(1))
        ak2=dabs(uxi(i,j))*(u(k,ipp2,j)-4*u(k,ipp,j)+6*u(k,i,j)
     :                 -4*u(k,inn,j)+u(k,inn2,j))/(4.0*dxi(1))   
        ak1=ak1+ak2
        du_xi=ak1/uxi(i,j) 

	endif

     	else

C	NEAR BOUNDARY ALWAYS CENTRAL	
        xpp=8.0*(u(k,ipp,j)-u(k,inn,j))
        xnn=u(k,ipp2,j)-u(k,inn2,j)

        du_xi=(1.0/12.0)*(xpp-xnn)/dxi(1)

	endif 

	if (j.ge.3.and.j.le.n(2)-2) then

	if (pec2.le.2.and.pec2.gt.-2) then

C	CENTRAL 4TH ORDER
       ypp=8.0*(u(k,i,jpp)-u(k,i,jnn))
       ynn=u(k,i,jpp2)-u(k,i,jnn2)

       du_et=(1.0/12.0)*(ypp-ynn)/dxi(2)

      else

c	upwind 3RD ORDER
        ak3=uet(i,j)*(-u(k,i,jpp2)+8*u(k,i,jpp)-8*u(k,i,jnn)
     :                 +u(k,i,jnn2))/(12.0*dxi(2))
        ak4=dabs(uet(i,j))*(u(k,i,jpp2)-4*u(k,i,jpp)+6*u(k,i,j)
     :                 -4*u(k,i,jnn)+u(k,i,jnn2))/(4.0*dxi(2))   
        ak3=ak3+ak4 

        du_et=ak3/uet(i,j)

       endif

	ELSE 
C	NEAR BOUNDARY ALWAYS CENTRAL
	du_et=0.5*(u(k,i,jpp)-u(k,i,jnn))/dxi(2)  

      end if
	
	conv(k)=uxi(i,j)*du_xi+uet(i,j)*du_et

        END do! end of k-loop

C-------------------------------------------------------
C		DIFFUSION
C-------------------------------------------------------

c-------------------------------------------------
c	Guessed velocity field (star )
c-------------------------------------------------
       dp_dxi=(p(ipp,j)-p(inn,j))/(2.*dxi(1))
       dp_de =(p(i,jpp)-p(i,jnn))/(2.*dxi(2))
       dp_dx=(dxix(i,j)*dp_dxi+dex(i,j)*dp_de)
       dp_dy=(dxiy(i,j)*dp_dxi+dey(i,j)*dp_de)

       qu(i,j)=dt*(-conv(1)-dp_dx)+u(1,i,j)
       qv(i,j)=dt*(-conv(2)-dp_dy+Ri*u(3,i,j))+u(2,i,j)
       qt(i,j)=-dt*conv(3)+u(3,i,j)

       qup(i,j)=qu(i,j)+dt*dp_dx
       qvp(i,j)=qv(i,j)+dt*dp_dy

       IF(j.eq.2) then
       sumu=bus(i)*u(1,i,jnn)+buse(i)*u(1,ipp,jnn)+busw(i)*u(1,inn,jnn)
       qu(i,j)=qu(i,j)-sumu
       sumv=bus(i)*u(2,i,jnn)+buse(i)*u(2,ipp,jnn)+busw(i)*u(2,inn,jnn)
       qv(i,j)=qv(i,j)-sumv
       sumt=bts(i)*u(3,i,jnn)+btse(i)*u(3,ipp,jnn)+btsw(i)*u(3,inn,jnn)
       qt(i,j)=qt(i,j)-sumt

       sumu=bus(i)*up(1,i,jnn)+buse(i)*up(1,ipp,jnn)+busw(i)*
     + up(1,inn,jnn)
       qup(i,j)=qup(i,j)-sumu
       sumv=bus(i)*up(2,i,jnn)+buse(i)*up(2,ipp,jnn)+busw(i)*
     + up(2,inn,jnn)
       qvp(i,j)=qvp(i,j)-sumv
       END if

       if (j.eq.n(2)-1) then
       sumu=bun(i)*u(1,i,jpp)+bune(i)*u(1,ipp,jpp)+bunw(i)*u(1,inn,jpp)
       qu(i,j)=qu(i,j)-sumu
       sumv=bun(i)*u(2,i,jpp)+bune(i)*u(2,ipp,jpp)+bunw(i)*u(2,inn,jpp)
       qv(i,j)=qv(i,j)-sumv
       sumt=btn(i)*u(3,i,jpp)+btne(i)*u(3,ipp,jpp)+btnw(i)*u(3,inn,jpp)
       qt(i,j)=qt(i,j)-sumt

       sumu=bun(i)*up(1,i,jpp)+bune(i)*up(1,ipp,jpp)+bunw(i)*
     + up(1,inn,jpp)
       qup(i,j)=qup(i,j)-sumu
       sumv=bun(i)*up(2,i,jpp)+bune(i)*up(2,ipp,jpp)+bunw(i)*
     + up(2,inn,jpp)
       qvp(i,j)=qvp(i,j)-sumv

       end if

        IF(i.eq.1)then
        qu(n(1),j)=qu(1,j)
        qv(n(1),j)=qv(1,j)
        qt(n(1),j)=qt(1,j)
        qup(n(1),j)=qup(1,j)
        qvp(n(1),j)=qvp(1,j)
        endif

        enddo
        enddo  ! end of space scan

c       'solving u-vel'
        do i=1,n(1)
        do j=1,n(2)
        sol(i,j)=u(1,i,j)
        end do
        end do

        CALL gauss(aup,aue,aus,aun,auw,ause,ausw,aune,aunw,auss,aussee,
     :  aussww,ausse,aussw,ausee,ausww,aunn,aunnee,aunnww,aunne,aunnw,
     :  aunee,aunww,auee,auww,sol,qu)

        do j=2,n(2)-1

        do i=1,n(1)-1
        us(1,i,j)=sol(i,j)
        if (i.eq.1) then
        us(1,n(1),j)=sol(i,j)
        end if
        end do
        end do

c       'solving v-vel'
        do i=1,n(1)
        do j=1,n(2)
        sol(i,j)=u(2,i,j)
        end do
        end do

        CALL gauss(aup,aue,aus,aun,auw,ause,ausw,aune,aunw,auss,aussee,
     :  aussww,ausse,aussw,ausee,ausww,aunn,aunnee,aunnww,aunne,aunnw,
     :  aunee,aunww,auee,auww,sol,qv)

        do j=2,n(2)-1
        do i=1,n(1)-1
        us(2,i,j)=sol(i,j)
        if (i.eq.1) then
        us(2,n(1),j)=sol(i,j)
        end if
        end do
        end do

c       'solving T'
        do i=1,n(1)
        do j=1,n(2)
        sol(i,j)=u(3,i,j)
        end do
        end do

        CALL gauss(atp,ate,ats,atn,atw,atse,atsw,atne,atnw,atss,atssee,
     :  atssww,atsse,atssw,atsee,atsww,atnn,atnnee,atnnww,atnne,atnnw,
     :  atnee,atnww,atee,atww,sol,qt)

        do j=2,n(2)-1
        do i=1,n(1)-1
        u(3,i,j)=sol(i,j)
        if (i.eq.1) then
        u(3,n(1),j)=sol(i,j)
        end if
        end do
        end do

c       'solving up-vel'

        do i=1,n(1)
        sol(i,1)=up(1,i,1)
        end do
        
        do i=1,n(1)
        do j=2,n(2)
        sol(i,j)=0.0
        end do
        end do

        CALL gauss(aup,aue,aus,aun,auw,ause,ausw,aune,aunw,auss,aussee,
     :  aussww,ausse,aussw,ausee,ausww,aunn,aunnee,aunnww,aunne,aunnw,
     :  aunee,aunww,auee,auww,sol,qup)

        do j=2,n(2)-1
        do i=1,n(1)-1
        up(1,i,j)=sol(i,j)
        if (i.eq.1) then
        up(1,n(1),j)=sol(i,j)
        end if
        end do
        end do

c       'solving vp-vel'

        do i=1,n(1)
        sol(i,1)=up(2,i,1)
        end do

        do i=1,n(1)
        do j=2,n(2)
        sol(i,j)=0.0
        end do
        end do

        CALL gauss(aup,aue,aus,aun,auw,ause,ausw,aune,aunw,auss,aussee,
     :  aussww,ausse,aussw,ausee,ausww,aunn,aunnee,aunnww,aunne,aunnw,
     :  aunee,aunww,auee,auww,sol,qvp)

        do j=2,n(2)-1
        do i=1,n(1)-1
        up(2,i,j)=sol(i,j)
        if (i.eq.1) then
        up(2,n(1),j)=sol(i,j)
        end if
        end do
        end do

c------------------------------------------------------
c      updating the bc for up
c------------------------------------------------------
        j=n(2)
        do i=1,n(1)-1

c        vnn=u(1,i,j)*xnox(i)+u(2,i,j)*xnoy(i)
        vnn=uinf*xnox(i)+vinf*xnoy(i)

        IF(vnn.ge.0)then
        UP(1,I,J)=U(1,I,J)
        UP(2,I,J)=U(2,I,J)
        else

        inn=i-1
        ipp=i+1
        if(i.eq.1)inn=n(1)-1

        jnn=j-1

        up(1,i,j)=(5.*up(1,i,jnn)-4.*up(1,i,jnn-1)+up(1,i,jnn-2))/2.0
        up(2,i,j)=(5.*up(2,i,jnn)-4.*up(2,i,jnn-1)+up(2,i,jnn-2))/2.0
        endif

        if (i.eq.1) then
        up(1,n(1),j)=up(1,i,j)
        up(2,n(1),j)=up(2,i,j)
        end if

        end do

c----------------------------------------------------------
c     calculation of star velocities at i+-1/2 and j+-1/2
c----------------------------------------------------------
        do j=2,n(2)-1
        do i=1,n(1)-1
         if (i.eq.1) then
        inn=n(1)-1
        ipp=i+1
        else
        inn=i-1
        ipp=i+1
        end if
        jpp=j+1
        jnn=j-1

          dpdxi_ip=(p(ipp,j)-p(i,j))/dxi(1)
          dpde_ip=(p(ipp,jpp)+p(i,jpp)-p(i,jnn)-p(ipp,jnn))/(4*dxi(2))

          dpdxi_in=(p(i,j)-p(inn,j))/dxi(1)
          dpde_in=(p(i,jpp)+p(inn,jpp)-p(i,jnn)-p(inn,jnn))/(4*dxi(2))

          dpdxi_jp=(p(ipp,jpp)-p(inn,jpp)+p(ipp,j)-p(inn,j))/(4*dxi(1))
          dpde_jp=(p(i,jpp)-p(i,j))/dxi(2)

          dpdxi_jn=(p(ipp,j)-p(inn,j)+p(ipp,jnn)-p(inn,jnn))/(4*dxi(1))
          dpde_jn=(p(i,j)-p(i,jnn))/dxi(2)

      us_ip=0.5*(up(1,i,j)+up(1,ipp,j))-0.5*dt*((dxix(i,j)+dxix(ipp,j))
     :*dpdxi_ip+(dex(i,j)+dex(ipp,j))*dpde_ip)

      us_in=0.5*(up(1,i,j)+up(1,inn,j))-0.5*dt*((dxix(i,j)+dxix(inn,j))
     :*dpdxi_in+(dex(i,j)+dex(inn,j))*dpde_in)

      us_jp=0.5*(up(1,i,j)+up(1,i,jpp))-0.5*dt*((dxix(i,j)+dxix(i,jpp))
     :*dpdxi_jp+(dex(i,j)+dex(i,jpp))*dpde_jp)

      us_jn=0.5*(up(1,i,j)+up(1,i,jnn))-0.5*dt*((dxix(i,j)+dxix(i,jnn))
     :*dpdxi_jn+(dex(i,j)+dex(i,jnn))*dpde_jn)

      vs_ip=0.5*(up(2,i,j)+up(2,ipp,j))-0.5*dt*((dxiy(i,j)+dxiy(ipp,j))
     :*dpdxi_ip+(dey(i,j)+dey(ipp,j))*dpde_ip)

      vs_in=0.5*(up(2,i,j)+up(2,inn,j))-0.5*dt*((dxiy(i,j)+dxiy(inn,j))
     :*dpdxi_in+(dey(i,j)+dey(inn,j))*dpde_in)

      vs_jp=0.5*(up(2,i,j)+up(2,i,jpp))-0.5*dt*((dxiy(i,j)+dxiy(i,jpp))
     :*dpdxi_jp+(dey(i,j)+dey(i,jpp))*dpde_jp)

      vs_jn=0.5*(up(2,i,j)+up(2,i,jnn))-0.5*dt*((dxiy(i,j)+dxiy(i,jnn))
     :*dpdxi_jn+(dey(i,j)+dey(i,jnn))*dpde_jn)

       dusdxi=(us_ip-us_in)/dxi(1)
       dusde=(us_jp-us_jn)/dxi(2)
       dvsdxi=(vs_ip-vs_in)/dxi(1)
       dvsde=(vs_jp-vs_jn)/dxi(2)

      q(i,j)=(dxix(i,j)*dusdxi)+(dex(i,j)*dusde)+(dxiy(i,j)*dvsdxi)+
     :(dey(i,j)*dvsde)

      q(i,j)=q(i,j)/dt

         enddo

        enddo

C      INITIALIZING THE PCORR

        do i=1,n(1)
        do j=1,n(2)

        pcor(i,j)=0
        uold(1,i,j)=u(1,i,j)
        uold(2,i,j)=u(2,i,j)

        end do
        end do

c----------------------------------------------------
c         performing  Guass Siedel iterations
c----------------------------------------------------
c      call goss9p(pcor,q)
       call sip9p(ap,ae,as,an,aw,ase,asw,ane,anw,pcor,q)


c------apply boundary condition on Pcor

        if (norm.eq.1) then

        WRITE(*,*)'hello'

        else

c       --------------solid-boundary-------------------

        j=1

	do i=1,n(1)-1

        pcor(i,j)=pcor(i,j+1)

        if (i.eq.1)pcor(n(1),j)=pcor(i,j)

        enddo


c       ----------------artificial boundary--------------
        j=n(2)

        do i=1,n(1)-1

        vnn=uinf*xnox(i)+vinf*xnoy(i)

        pcor(i,j)=0
        IF(vnn.ge.0)pcor(i,j)=pcor(i,j-1)


        if (i.eq.1) then
         pcor(n(1),j)=pcor(i,j)
        end if

        enddo

        endif
c--------------------------------------------------------
c-----updating U and V from Pcor in the interior
c---------------------------------------------------------
      do i=1,n(1)-1
      do j=2,n(2)-1

        if (i.eq.1) then
        inn=n(1)-1
        ipp=i+1
        else
        inn=i-1
        ipp=i+1
        end if
        jpp=j+1
        jnn=j-1

      dpcor_dxi=0.5*(pcor(ipp,j)-pcor(inn,j))/dxi(1)

      dpcor_de=0.5*(pcor(i,jpp)-pcor(i,jnn))/dxi(2)

      u(1,i,j)=us(1,i,j)-dt*(dxix(i,j)*dpcor_dxi+dex(i,j)*dpcor_de)

      u(2,i,j)=us(2,i,j)-dt*(dxiy(i,j)*dpcor_dxi+dey(i,j)*dpcor_de)
c
      if (i.eq.1) then

      u(1,n(1),j)=u(1,i,j)
      u(2,n(1),j)=u(2,i,j)

      end if

      enddo
      enddo

         do i=1,n(1)-1
         do j=2,n(2)-1

         p(i,j)=P(i,j)+pcor(i,j)

         if (i.eq.1) then

         p(n(1),j)=p(i,j)

         end if

         end do
         end do

c==========================================================
c       Evaluating Vr and Vth from U and V velocity just
c       before the outer plane in vr,vth index 1 is n(2)-1
c==========================================================
        j=n(2)-1
        do i=1,n(1)-1

	costh=x(1,i,j)/(x(1,i,j)**2+x(2,i,j)**2)**0.5
        sinth=x(2,i,j)/(x(1,i,j)**2+x(2,i,j)**2)**0.5

        vr(1,i)=u(1,i,j)*costh+u(2,i,j)*sinth
        vth(1,i)=-u(1,i,j)*sinth+u(2,i,j)*costh

        if (i.eq.1) then
        vr(1,n(1))=vr(1,i)
	vth(1,n(1))=vth(1,i)
        end if

         end do

c===========================================================
c       Calculating circulation at the 2nd last level in jth
c===========================================================
        circ=0.0
        j=n(2)-1
        do i=1,n(1)-1
        de=1.0/(n(1)-1)
      f1=(u(1,i,j)*dey(i,j)-u(2,i,j)*dex(i,j))*ABS(ajac(i,j))
      f2=(u(1,i+1,j)*dey(i+1,j)-u(2,i+1,j)*dex(i+1,j))*ABS(ajac(i+1,j))

        circ=circ+de*0.5*(f1+f2)

        end do

c=========================================================
c       Predicting values for vr and vth at outer
c=========================================================
        j=n(2)
        do i=1,n(1)-1
        eps=1e-2
      cr=(x(1,i,j-1)**2+x(2,i,j-1)**2)**0.5/(x(1,i,j)**2+x(2,i,j)**2)**
     :0.5
        costh=x(1,i,j)/((x(1,i,j)**2+x(2,i,j)**2)**0.5)
	sinth=x(2,i,j)/((x(1,i,j)**2+x(2,i,j)**2)**0.5)

	vrinf=uinf*costh+vinf*sinth
	vtinf=-uinf*sinth+vinf*costh
	
	if (ABS(circ).gt.eps) then
        kk=1
        else
        kk=2
        end if

        vr(2,i)=vr(1,i)*cr**2+vrinf*(1-cr**2)
        vth(2,i)=vth(1,i)*cr**kk+vtinf*(1-cr**kk)

        if (i.eq.1) then
        vr(2,n(1))=vr(2,i)
	vth(2,n(1))=vth(2,i)
       
        end if

        end do

c--------------------------------------------------
c       updating the bc of U And V
c---------------------------------------------------


c	-----------------cylinder_oscillation--------------

          j=1
          do k=1,2
          do i=1,n(1)

           IF(k.eq.1)then
            U(k,i,j)=-speed_amp*cos(2.0*pi*F*time)*x(2,i,j)      ! line edited
            up(k,i,j)=U(k,i,j)
           else
            U(k,i,j)=speed_amp*cos(2.0*pi*F*time)*x(1,i,j)       ! line edited
            up(k,i,j)=U(k,i,j)
           endif

          end do
          end do                        

        j=n(2)

        do i=1,n(1)-1

        vnn=uinf*xnox(i)+vinf*xnoy(i)
        IF(vnn.ge.0)then
         u(1,i,j)=uinf
         U(2,i,j)=vinf
         u(3,i,j)=0.0

        else

         costh=x(1,i,j)/(x(1,i,j)**2+x(2,i,j)**2)**0.5
         sinth=x(2,i,j)/(x(1,i,j)**2+x(2,i,j)**2)**0.5

         U(1,i,j)=costh*vr(2,i)-sinth*vth(2,i)
         U(2,i,j)=sinth*vr(2,i)+costh*vth(2,i)
         U(3,i,j)=uold(3,i,j)-(uet(i,j)*dt/dxi(2))
     :                       *(uold(3,i,j)-uold(3,i,j-1))
        endif

        if (i.eq.1) then
         U(1,n(1),j)=U(1,1,j)
         U(2,n(1),j)=U(2,1,j)
         U(3,n(1),j)=U(3,1,j)
        end if

        enddo


c=============================
c	apply BE for updating pressure
c=============================

c========================================================================
c	APPLYING MOMENTUM EQUATION ON inlet AND SOLID BOUNDARY
c       and Gresho's condition at outflow
C========================================================================
c     obtaining the new uxi and uet
       do j=1,n(2)
        do i=1,n(1)
         uxi(i,j)=dxix(i,j)*U(1,i,j)+dxiy(i,j)*u(2,i,j)
	 uet(i,j)=dex(i,j)*U(1,i,j)+dey(i,j)*u(2,i,j)
        end do
       end do  	

c     at solid boundary
	J=1
	do i=1,n(1)-1

	do k=1,2

	conv(k)=0
	d2u(k)=0
	alc(k)=0

	if (i.eq.1) then

	ipp=i+1
	inn=n(1)-1
      else
      ipp=i+1
      inn=i-1
	endif

	jpp=j+1
        jpp2=j+2

c	diffusive
	aa=alph(i,j)*(U(k,ipp,j)+U(k,inn,j)-2*U(k,i,j))/(dxi(1)**2)

      gg=gamma(i,j)*(U(k,i,jpp+1)+U(k,i,j)-2*U(k,i,jpp))/(dxi(2)**2)

      bb=beta(i,j)*(U(k,ipp,jpp)+U(k,inn,j)-U(k,inn,jpp)-U(k,ipp,j))
     :/(2*dxi(1)*dxi(2))

       qqq=Q1(i,j)*(-3*U(k,i,j)+4*U(k,i,jpp)-U(k,i,jpp2))/(2*dxi(2))


       d2u(k)=aa+gg-2*bb+qqq

c	convective
      conv(k)=uxi(i,j)*0.5*(u(k,ipp,j)-u(k,inn,j))/dxi(1)
      
      conv(k)=conv(k)+uet(i,j)*(u(k,i,jpp)-u(k,i,j))/dxi(2)

c     local

c	alc(k)=(u(k,i,j)-uold(k,i,j))/dt
	

        IF(k.eq.1)then
            alc(k)=accn_amp*sin(2.0*pi*F*time)*x(2,i,j)      ! line edited
        else
            alc(k)=-accn_amp*sin(2.0*pi*F*time)*x(1,i,j)       ! line edited
        endif


	if (k.eq.1) dp_dx=1.0*d2u(k)/re-conv(k)-alc(k)
	if (k.eq.2) dp_dy=1.0*d2u(k)/re-conv(k)-alc(k)+Ri*u(3,i,j)

	enddo

	p(i,j)=p(i,j+1)-(dp_dx*(-dxiy(i,j)*ajac(i,j))
     :+dp_dy*(dxix(i,j)*ajac(i,j)))*dxi(2)

	if(i.eq.1) p(n(1),j)=p(i,j)

	enddo

c     at exit boundary
	J=N(2)

	do i=1,n(1)-1
          vnn=uinf*xnox(i)+vinf*xnoy(i)
          IF(vnn.ge.0)then
c       -------------momentum equation----------------------------------
           do k=1,2
	   conv(k)=0
	   d2u(k)=0
	   alc(k)=0

           ipp=i+1
           inn=i-1
	   if(i.eq.1)inn=n(1)-1

           jnn=j-1
           jnn2=j-2

c	diffusive
	aa=alph(i,j)*(U(k,ipp,j)+U(k,inn,j)-2*U(k,i,j))/(dxi(1)**2)

      gg=gamma(i,j)*(U(k,i,j)+U(k,i,jnn-1)-2*U(k,i,jnn))/(dxi(2)**2)

      bb=beta(i,j)*(U(k,ipp,j)+U(k,inn,jnn)-U(k,ipp,jnn)-U(k,inn,j))
     :/(2*dxi(1)*dxi(2))

      qqq=Q1(i,j)*(3*U(k,i,j)-4*U(k,i,jnn)+U(k,i,jnn2))/(2*dxi(2))

           d2u(k)=aa+gg-2*bb+qqq

c	convective
      conv(k)=uxi(i,j)*0.5*(u(k,ipp,j)-u(k,inn,j))/dxi(1)

      conv(k)=conv(k)+uet(i,j)*(3.*u(k,i,j)-4*u(k,i,jnn)
     :                  +u(k,i,jnn2))/(2*dxi(2))

c     local
	alc(k)=(u(k,i,j)-uold(k,i,j))/dt

	if (k.eq.1) dp_dx=1.*d2u(k)/re-conv(k)-alc(k)
	if (k.eq.2) dp_dy=1.*d2u(k)/re-conv(k)-alc(k)+Ri*U(3,i,j)

	enddo   ! k-loop

	p(i,j)=p(i,j-1)+(dp_dx*(-dxiy(i,j)*ajac(i,j))
     :+dp_dy*(dxix(i,j)*ajac(i,j)))*dxi(2)

         ELSE

c       -------------gresho's condition---------------------------------
       p(i,j)=0.5*(1.0/re)*((3*uet(i,j)-4*uet(i,j-1)+uet(i,j-2))/dxi(2))
         ENDIF

         if(i.eq.1)p(n(1),j)=p(i,j)

	enddo

c----------------------------------
c-----calculation of si
c----------------------------------
      j=1
      do i=1,n(1)
       si(i,j)=0
      end do

      do i=1,n(1)
      do j=2,n(2)

      ca=(dxix(i,j)*u(1,i,j)*ABS(ajac(i,j))+dxix(i,j-1)*u(1,i
     :,j-1)*ABS(ajac(i,j-1)))
      cb=(dxiy(i,j)*u(2,i,j)*ABS(ajac(i,j))+dxiy(i,j-1)*u(2,i
     :,j-1)*ABS(ajac(i,j-1)))

      si(i,j)=si(i,j-1)+(ca+cb)*0.5*dxi(2)

      enddo
      enddo

C----------------------------
c      DILATION AND VORTICITY
C----------------------------
        dmax=0
        do i=1,n(1)-1
        do j=2,n(2)-1

        if (i.eq.1) then
        inn=n(1)-1
        ipp=i+1
        else
        inn=i-1
        ipp=i+1
        end if
        jpp=j+1
        jnn=j-1

          dil(i,j)=dxix(i,j)*(U(1,ipp,j)-U(1,inn,j))/(2*dxi(1))+dex(i,j)
     :*(U(1,i,jpp)-U(1,i,jnn))/(2*dxi(2))+dey(i,j)*(U(2,i,jpp)-U(2,i,jnn
     :))/(2*dxi(2))+dxiy(i,j)*(U(2,ipp,j)-U(2,inn,j))/(2*dxi(1))

        dv_dxi=0.5/dxi(1)*(u(2,ipp,j)-u(2,inn,j))
        dv_det=0.5/dxi(2)*(u(2,i,jpp)-u(2,i,jnn))

        dv_dx=dxix(i,j)*dv_dxi + dex(i,j)*dv_det

        du_dxi=0.5/dxi(1)*(u(1,ipp,j)-u(1,inn,j))
        du_det=0.5/dxi(2)*(u(1,i,jpp)-u(1,i,jnn))

        du_dy=dxiy(i,j)*du_dxi + dey(i,j)*du_det

        vort(i,j)=dv_dx-du_dy

       if (i.eq.1) then
         dil(n(1),j)=dil(i,j)
         vort(n(1),j)=vort(i,j)
       end if

        if (dil(i,j).gt.dmax) then
           dmax=dil(i,j)
        end if

        end do
        end do

        do j=1,n(2),n(2)-1
        do i=1,n(1)-1
        if (i.eq.1) then
        inn=n(1)-1
        ipp=i+1
        else
        inn=i-1
        ipp=i+1
        end if
        jpp=j+1
        jnn=j-1

        dv_dxi=0.5/dxi(1)*(u(2,ipp,j)-u(2,inn,j))
        IF(j.eq.1)dv_det=1.0/dxi(2)*(u(2,i,jpp)-u(2,i,j))
        IF(j.eq.n(2))dv_det=1.0/dxi(2)*(u(2,i,j)-u(2,i,jnn))

        dv_dx=dxix(i,j)*dv_dxi + dex(i,j)*dv_det

        du_dxi=0.5/dxi(1)*(u(1,ipp,j)-u(1,inn,j))
        IF(j.eq.1)du_det=1.0/dxi(2)*(u(1,i,jpp)-u(1,i,j))
        IF(j.eq.n(2))du_det=1.0/dxi(2)*(u(1,i,j)-u(1,i,jnn))

        du_dy=dxiy(i,j)*du_dxi+dey(i,j)*du_det

        vort(i,j)=dv_dx-du_dy

        if (i.eq.1) then
	vort(n(1),j)=vort(i,j)
        end if

        end do
        end do

        WRITE(*,*)loop,Dmax

C=========================================================
C       Calculation of lift,drag,moment and Nusselt number
C=========================================================

c       ----------------------------------------------------
c       calculating pressure and vorticity surface integrals
c       for forces
c       ----------------------------------------------------
        j=1

        pr_x=0.0
        pr_y=0.0
        vor_x=0.0
        vor_y=0.0

	do i=1,n(1)-1

        ip=i+1

        PJ1=p(i,j)*ajac(i,j)
        pj2=p(ip,j)*ajac(ip,j)

        VJ1=vort(i,j)*ajac(i,j)
        VJ2=vort(ip,j)*ajac(ip,j)

        fp1_x=PJ1*dex(i,j)
        fp2_x=PJ2*dex(ip,j)

        fp1_y=PJ1*dey(i,j)
        fp2_y=PJ2*dey(ip,j)

        fv1_x=VJ1*dey(i,j)
        fv2_x=VJ2*dey(ip,j)

        fv1_y=VJ1*dex(i,j)
        fv2_y=VJ2*dex(ip,j)

        pr_x=pr_x+0.5*dxi(1)*(fp1_x+fp2_x)
        pr_y=pr_y+0.5*dxi(1)*(fp1_y+fp2_y)

	vor_x=vor_x+0.5*dxi(1)*(fv1_x+fv2_x)
        vor_y=vor_y+0.5*dxi(1)*(fv1_y+fv2_y)

	end do
        cx=2*pr_x+(2.0/Re)*vor_x
        cy=2*pr_y-(2.0/Re)*vor_y




         CL_pr=2*pr_y*SIN(alpha*pi/180.0)-2*pr_x*COS(alpha*pi/180.0)
         CD_pr=2*pr_y*COS(alpha*pi/180.0)+2*pr_x*SIN(alpha*pi/180.0)
         CL_vor=-(2.0/Re)*vor_y*SIN(alpha*pi/180.0)
     :   -(2.0/Re)*vor_x*COS(alpha*pi/180.0)         
         CD_vor=-(2.0/Re)*vor_y*COS(alpha*pi/180.0)
     :   +(2.0/Re)*vor_x*SIN(alpha*pi/180.0)


        cl=cy*SIN(alpha*pi/180.0)-cx*COS(alpha*pi/180.0)
        cd=cy*COS(alpha*pi/180.0)+cx*SIN(alpha*pi/180.0)

c       -------------------------------------------------------
c       calculating surface pressure,vorticity and temp. integrals
c       for moment coefficient and Nusselt number
c       -------------------------------------------------------
        press_i=0.0
        vor_i=0.0
        temp_i=0.0

	do i=1,n(1)-1

        ip=i+1

        PJ1=p(i,j)*ajac(i,j)
        pj2=p(ip,j)*ajac(ip,j)

        VJ1=vort(i,j)*ajac(i,j)
        VJ2=vort(ip,j)*ajac(ip,j)

        TJ1=ajac(i,j)*(dex(i,j)**2+dey(i,j)**2)
        TJ2=ajac(ip,j)*(dex(ip,j)**2+dey(ip,j)**2)

        fp1=PJ1*(x(1,i,j)*dey(i,j)-x(2,i,j)*dex(i,j))
        fp2=PJ2*(x(1,ip,j)*dey(ip,j)-x(2,ip,j)*dex(ip,j))

        fv1=VJ1*(x(1,i,j)*dex(i,j)+x(2,i,j)*dey(i,j))
        fv2=VJ2*(x(1,ip,j)*dex(ip,j)+x(2,ip,j)*dey(ip,j))

        fh1=TJ1*(4*U(3,i,j+1)-3*U(3,i,j)-U(3,i,j+2))/(2*dxi(2))
        fh2=TJ2*(4*U(3,ip,j+1)-3*U(3,ip,j)-U(3,ip,j+2))/(2*dxi(2))

        press_i=press_i+0.5*dxi(1)*(fp1+fp2)

	vor_i=vor_i+0.5*dxi(1)*(fv1+fv2)

        temp_i=temp_i+0.5*(fh1+fh2)*dxi(1)

	end do

        cm=2*press_i-(2.0/Re)*vor_i
         Nuss=(2*temp_i)/(Pi*(3*(1+(1./ar))-((3+(1./ar))*((3./ar)+1))**
     :  0.5))

c----------------------------------------------------------
c       FILE WRITING
C----------------------------------------------------------
      if(MOD(loop,100).eq.0) then

      open (1,FILE='spt100.dat',status='unknown')
	WRITE(1,*) 'zone'
	WRITE(1,*) 'I=',n(1)
	WRITE(1,*) 'J=',n(2)
      do j=1,n(2)
       do i=1,n(1)

      WRITE(1,"(2(f15.9,1x),6(e20.13,2x))")x(1,i,j),x(2,i,j)
     :,u(1,i,j),u(2,i,j),u(3,i,j),p(i,j),si(i,j),vort(i,j)

       enddo
      WRITE(1,*)

      enddo
      CLOSE(1)

      open (1,FILE='spa100.dat',status='unknown',form='unformatted')
      WRITE(1)loop,time,dmax

        write(1)x,si,u,p 

       CLOSE(1)

      open (1,FILE='COEFF_HIS.dat',status='unknown',ACCESS='append')

      WRITE(1,"(5(f15.8,2x))")time,CL,CD,CM,NUSS

      CLOSE(1)

    
      open (1,FILE='COEFF_HIS_pr_vor.dat'
     : ,status='unknown',ACCESS='append')

      WRITE(1,"(5(f15.8,2x))")time,CL_pr,CD_pr,CL_vor,CD_vor

      CLOSE(1)


c================================================================
c       local nusselt number profile on cylinder
c================================================================
      open (5,FILE='SURF_DIST.dat',status='unknown')

       do i=1,n(1)
       dthdn=-(4*U(3,i,2)-3*U(3,i,1)-U(3,i,3))/(2*dxi(2))
       dthdn=dthdn*sqrt(dex(i,1)**2+dey(i,1)**2)

       WRITE(5,*)(i-1)*dxi(1),p(i,1),vort(i,1),dthdn
       end do

      CLOSE(5)

      end if

      if (iiflag.eq.1) then

	if (loop.eq.loop_snap) then
	nsnap=nsnap+1

	if (nsnap.eq.(maxsnap+1)) goto 101

	open (1,FILE=filnam(nsnap),status='unknown')

      do j=1,n(2)
       do i=1,n(1)

      WRITE(1,"(2(f13.9,1x),3(e12.5,2x))")x(1,i,j),x(2,i,j),si(i,j)
     :,u(3,i,j),vort(i,j)

       enddo
      WRITE(1,*)

      enddo
      CLOSE(1)

	loop_snap=loop_snap+i_loop

	endif

      endif

101    continue
       END DO
C      END OF TIME LOOP

  100   stop
        end

      subroutine sip9p(ap,ae,as,an,aw,ase,asw,ane,anw,phi,q)

      PARAMETER (NP1=350,NP2=570)
      implicit REAL *8 (a-h,o-z)

      COMMON /bee/n(2)
      dimension be(np1,np2),bw(np1,np2),bs(np1,np2),bn(np1,np2),
     :bse(np1,np2),bne(np1,np2),bnw(np1,np2),bsw(np1,np2),bp(np1,np2)

      DIMENSION RES(np1,np2),Qp(np1,np2),del(np1,np2),phi(np1,np2),
     :q(np1,np2),phio(np1,np2),ae(np1,np2),aw(np1,np2),as(np1,np2),
     :an(np1,np2),ase(np1,np2),ane(np1,np2),asw(np1,np2),anw(np1,np2)
     :,ap(np1,np2)

       tol=0.75e-2
       maxiter=100000
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
      enddo
      enddo

      do j=2,n(2)-1
       do i=1,n(1)-1

       IF (i.EQ.1) THEN
        inn=n(1)-1
        ipp=i+1

        else

        inn=i-1
        ipp=i+1

       END IF
        jpp=j+1
        jnn=j-1

      bsw(i,j)=asw(i,j)

      bw(i,j)=(aw(i,j)+alp*anw(i,j)-bsw(i,j)*bn(inn,jnn))/
     :(1+alp*bn(inn,j))

      bs(i,j)=(as(i,j)+alp*ase(i,j)-bsw(i,j)*be(inn,jnn))/
     :(1+alp*be(i,jnn))

      ad=anw(i,j)+ase(i,j)-bs(i,j)*be(i,jnn)-bw(i,j)*bn(inn,j)

      bp(i,j)=ap(i,j)-alp*ad-bs(i,j)*bn(i,jnn)-bw(i,j)*be(inn,j)-
     :bsw(i,j)*bne(inn,jnn)

      bn(i,j)=(an(i,j)+alp*anw(i,j)-alp*bw(i,j)*bn(inn,j)-
     : bw(i,j)*bne(inn,j))/bp(i,j)

      be(i,j)=(ae(i,j)+alp*ase(i,j)-alp*bs(i,j)*be(i,jnn)-
     :bs(i,j)*bne(i,jnn))/bp(i,j)

      bne(i,j)=ane(i,j)/bp(i,j)

       if (i.eq.1) then

       bsw(n(1),j)=bsw(i,j)
       bn(n(1),j)=bn(i,j)
       bs(n(1),j)=bs(i,j)
       bse(n(1),j)=bse(i,j)
       bnw(n(1),j)=bnw(i,j)
       bne(n(1),j)=bne(i,j)
       be(n(1),j)=be(i,j)
       bw(n(1),j)=bw(i,j)
       bp(n(1),j)=bp(i,j)

      end if

       end do
      end do

      do j=1,n(2)
       do i=1,n(1)

       qp(i,j)=0
       del(i,j)=0
      enddo
      enddo

      DO iter=1,MAXITER

      do i=1,n(1)
       do j=1,n(2)

       phio(i,j)=phi(i,j)

       end do
      end do

       SSUM=0.0
       do j=2,n(2)-1
         do i=1,n(1)-1

        IF (i.eq.1) THEN
        inn=n(1)-1
        ipp=i+1

        else

        inn=i-1
        ipp=i+1

       END IF
        jpp=j+1
        jnn=j-1

      res(i,j)=q(i,j)-ap(i,j)*phi(i,j)-ae(i,j)*phi(ipp,j)-an(i,j)*
     :phi(i,jpp)-as(i,j)*phi(i,jnn)-aw(i,j)*phi(inn,j)-anw(i,j)*
     :phi(inn,jpp)-ane(i,j)*phi(ipp,jpp)-asw(i,j)*phi(inn,jnn)
     :-ase(i,j)*phi(ipp,jnn)

       SSUM=ssum+ABS(res(i,j))

      qp(i,j)=(res(i,j)-bs(i,j)*qp(i,jnn)-bw(i,j)*qp(inn,j)
     :-bsw(i,j)*qp(inn,jnn))/bp(i,j)

       if (i.eq.1) then

       res(n(1),j)=res(i,j)
       qp(n(1),j)=qp(i,j)

       end if

       enddo
      enddo

       IF(iter.eq.1)then
       if(ssum.ne.0.0)then 
       sumnor=ssum
       else
       sumnor=1.0
       endif
       endif


       sumav=ssum/sumnor

       do j=n(2)-1,2,-1
       do i=n(1)-1,1,-1

       IF (i.eq.1) THEN
        inn=n(1)-1
        ipp=i+1

        else

        inn=i-1
        ipp=i+1

       END IF
        jpp=j+1
        jnn=j-1

       del(i,j)=qp(i,j)-bn(i,j)*del(i,jpp)-be(i,j)*del(ipp,j)-
     :bne(i,j)*del(ipp,jpp)
       phi(i,j)=phi(i,j)+del(i,j)

       if (i.eq.1) then
        phi(n(1),j)=phi(i,j)
       end if


       enddo
       enddo

        IF(sumav.lt.tol)GOTO 20

      END DO
   20  continue

      return
      end

      subroutine gauss(ap,ae,as,an,aw,ase,asw,ane,anw,ass,assee,
     :  assww,asse,assw,asee,asww,ann,annee,annww,anne,annw,
     :  anee,anww,aee,aww,phi,q)

      PARAMETER (NP1=350,NP2=570)
      implicit REAL *8 (a-h,o-z)

      COMMON /bee/n(2)

      DIMENSION RES(np1,np2),phi(np1,np2),q(np1,np2),phio(np1,np2),
     :ae(np1,np2),aw(np1,np2),as(np1,np2),an(np1,np2),ase(np1,np2),
     :ane(np1,np2),asw(np1,np2),anw(np1,np2),ap(np1,np2),
     :ann(np1,np2),ass(np1,np2),aee(np1,np2),aww(np1,np2),
     :annee(np1,np2),annww(np1,np2),assee(np1,np2),assww(np1,np2),
     :anne(np1,np2),annw(np1,np2),asse(np1,np2),assw(np1,np2),
     :anee(np1,np2),anww(np1,np2),asee(np1,np2),asww(np1,np2)

       tol=0.75e-2
       maxiter=100000
      DO iter=1,MAXITER

      do i=1,n(1)
       do j=1,n(2)

       phio(i,j)=phi(i,j)

       end do
      end do

       SSUM=0.0
       do j=2,n(2)-1
         do i=1,n(1)-1

        inn=i-1
        inn2=i-2
        ipp=i+1
        ipp2=i+2

        jpp=j+1
        jpp2=j+2
        jnn=j-1
        jnn2=j-2

        IF (i.eq.1) THEN
        inn=n(1)-1
        inn2=n(1)-2
        END IF

        IF(i.eq.2)then
        inn2=n(1)-1
        endif

        IF(i.eq.n(1)-1)then
        ipp2=2
        endif

        IF(j.eq.2.or.j.eq.n(2)-1)then

        res(i,j)=q(i,j)-ap(i,j)*phi(i,j)-ae(i,j)*phi(ipp,j)-an(i,j)*
     :  phi(i,jpp)-as(i,j)*phi(i,jnn)-aw(i,j)*phi(inn,j)-anw(i,j)*
     :  phi(inn,jpp)-ane(i,j)*phi(ipp,jpp)-asw(i,j)*phi(inn,jnn)
     :  -ase(i,j)*phi(ipp,jnn)

        else

c       forth order
        res(i,j)=q(i,j)-ap(i,j)*phi(i,j)-ae(i,j)*phi(ipp,j)-an(i,j)*
     :  phi(i,jpp)-as(i,j)*phi(i,jnn)-aw(i,j)*phi(inn,j)-anw(i,j)*
     :  phi(inn,jpp)-ane(i,j)*phi(ipp,jpp)-asw(i,j)*phi(inn,jnn)
     : -ase(i,j)*phi(ipp,jnn)-aee(i,j)*phi(ipp2,j)-aww(i,j)*phi(inn2,j)
     : -annee(i,j)*phi(ipp2,jpp2)-anee(i,j)*phi(ipp2,jpp)
     : -asee(i,j)*phi(ipp2,jnn)-assee(i,j)*phi(ipp2,jnn2)
     : -anne(i,j)*phi(ipp,jpp2)-asse(i,j)*phi(ipp,jnn2)
     : -annw(i,j)*phi(inn,jpp2)-assw(i,j)*phi(inn,jnn2)
     : -annww(i,j)*phi(inn2,jpp2)-anww(i,j)*phi(inn2,jpp)
     : -asww(i,j)*phi(inn2,jnn)-assww(i,j)*phi(inn2,jnn2)
     : -ann(i,j)*phi(i,jpp2)-ass(i,j)*phi(i,jnn2)

        endif
       SSUM=ssum+ABS(res(i,j))
       if (i.eq.1) then

       res(n(1),j)=res(i,j)

       end if

       enddo
      enddo

       IF(iter.eq.1)then
       if(ssum.ne.0.0)then 
       sumnor=ssum
       else
       sumnor=1.0
       endif
       endif

       sumav=ssum/sumnor

       do j=2,n(2)-1
       do i=1,n(1)-1

        inn=i-1
        inn2=i-2
        ipp=i+1
        ipp2=i+2

        jpp=j+1
        jpp2=j+2
        jnn=j-1
        jnn2=j-2

        IF (i.eq.1) THEN
        inn=n(1)-1
        inn2=n(1)-2
        END IF

        IF(i.eq.2)then
        inn2=n(1)-1
        endif

        IF(i.eq.n(1)-1)then
        ipp2=2
        endif

        IF(j.eq.2.or.j.eq.n(2)-1)then

        phi(i,j)=(q(i,j)-ae(i,j)*phi(ipp,j)-an(i,j)*
     :  phi(i,jpp)-as(i,j)*phi(i,jnn)-aw(i,j)*phi(inn,j)-anw(i,j)*
     :  phi(inn,jpp)-ane(i,j)*phi(ipp,jpp)-asw(i,j)*phi(inn,jnn)
     :  -ase(i,j)*phi(ipp,jnn))/ap(i,j)
        else

c       forth order

        phi(i,j)=(q(i,j)-ae(i,j)*phi(ipp,j)-an(i,j)*
     :  phi(i,jpp)-as(i,j)*phi(i,jnn)-aw(i,j)*phi(inn,j)-anw(i,j)*
     :  phi(inn,jpp)-ane(i,j)*phi(ipp,jpp)-asw(i,j)*phi(inn,jnn)
     : -ase(i,j)*phi(ipp,jnn)-aee(i,j)*phi(ipp2,j)-aww(i,j)*phi(inn2,j)
     : -annee(i,j)*phi(ipp2,jpp2)-anee(i,j)*phi(ipp2,jpp)
     : -asee(i,j)*phi(ipp2,jnn)-assee(i,j)*phi(ipp2,jnn2)
     : -anne(i,j)*phi(ipp,jpp2)-asse(i,j)*phi(ipp,jnn2)
     : -annw(i,j)*phi(inn,jpp2)-assw(i,j)*phi(inn,jnn2)
     : -annww(i,j)*phi(inn2,jpp2)-anww(i,j)*phi(inn2,jpp)
     : -asww(i,j)*phi(inn2,jnn)-assww(i,j)*phi(inn2,jnn2)
     : -ann(i,j)*phi(i,jpp2)-ass(i,j)*phi(i,jnn2))/ap(i,j)

        endif
       if (i.eq.1) then
        phi(n(1),j)=phi(i,j)
 
       end if

       enddo
       enddo

       IF(sumav.lt.tol)GOTO 20


      END DO
   20 continue

      return
      end

