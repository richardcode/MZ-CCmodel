      program gnanadesikan
c     
c     function: multilayer Gnanadesikan model
c     April 2013
c
c     this version includes hacked-in heat flux diagnostics at 3 depths
c       and fixes the temperature levels independent of surface warming
c
      implicit none
c
      integer nz,nt,nout,nmax
c
      parameter(nz=100,nt=3,nout=3001,nmax=50000001)
c
      real*8 h(nz),d(0:nz),theta(0:nz)
      real*8 qek(nz),qeddy(nz),qn(nz),qu(nz)
      real*8 qtot(nz,nt),ab(nt),fanth(nmax),famoc(nmax)
      real*8 dout(0:nz,nout)
      real*8 qekout(0:nz,nout),qeddyout(0:nz,nout)
      real*8 qnout(0:nz,nout),quout(0:nz,nout),qtotout(0:nz,nout)
      real*8 heatout(nout),heatupout(0:nout+1),thetasout(nout)
      real*8 heatekout(0:nout+1),heatuout(0:nout+1)
      real*8 heateddyout(0:nout+1),heatnout(0:nout+1)
      real*8 pi,year,tmax,dt,t1,t2,t3,t4,depth,dek,qek0,ddrake,socean
      real*8 theta0,thetas,delta,surface,tnadw10,tnadw20,a
      real*8 dwind,qek00,kappagm0,deddy,dtnadw,tnadw1,tnadw2
      real*8 thetas0,dthetas,dtheta,rho0,cp,heat,heatuptake
      real*8 hmin,kappav,kappagm,lx,ly,qn0,dqn,qn00
      real*8 heatn,heatu,heateddy,heatek
      real*8 heatn1,heatu1,heateddy1,heatek1,heatuptake1
      real*8 heatn2,heatu2,heateddy2,heatek2,heatuptake2
      real*8 heatn3,heatu3,heateddy3,heatek3,heatuptake3
      real*8 heatup1out(0:nout+1)
      real*8 heatek1out(0:nout+1),heatu1out(0:nout+1)
      real*8 heateddy1out(0:nout+1),heatn1out(0:nout+1)
      real*8 heatup2out(0:nout+1)
      real*8 heatek2out(0:nout+1),heatu2out(0:nout+1)
      real*8 heateddy2out(0:nout+1),heatn2out(0:nout+1)
      real*8 heatup3out(0:nout+1)
      real*8 heatek3out(0:nout+1),heatu3out(0:nout+1)
      real*8 heateddy3out(0:nout+1),heatn3out(0:nout+1)
      real*8 hhmax,hhmin,qekmax,qekmin,qumax,qumin
      real*8 qedmax,qedmin,qnmax,qnmin

      real*8 Tatm,TT(nmax),TTmax,c,lambda
      real*8 To,co,gamma
      real*8 forcing(nmax),amp,init_for
      real*8 rhs(nt),rhso(nt)
      real*8 time(nout),Tout(nout),forcingout(nout),Toout(nout)
      
      real*8 heattout(nout)

c		
      integer nstop,n,k,kmin,n1,n2,n3,n4,ndump,ndumpp


      pi=3.14159265358979
      Tatm=0.                   ! global mean atm temp perturbation at t0
      To=0.                     ! global mean surface temp perturbation at t0
      amp=4.                    ! forcing amplitude
      init_for=0.               ! initial forcing
      c=1000000.                ! heat capacity of the atm
      co=1000000000.            ! heat capacity of the ocean
      lambda=2.                 ! feedback rate
      gamma=1.                  ! coupling strength btw atm and ocean
      

c     Southern ocean switch (0 for no Drake Passage, 1 for Drake Passage)
        socean=1.

c     eddy diffusivities (placed here to scale dt with kappav)
C            fractional change in eddy transport (0.1 = +10%)
        kappav=1.e-5           
        kappagm0=1.e3*socean   
        deddy=0.           

c     Southern Ocean Ekman transport, Ekman depth, Drake Passage depth,
c            fractional change in Ekman transport (0.1 = +10%)
        qek0=30.e6*socean      
        dek=100.           
        depth=5.e3           
        ddrake=0.8*depth     
        dwind=0.              

c     NADW formation, fractional reduction in NADW formation (0.5 = -50%),
c            temperature range of NADW formation, NADW warming 
        qn0=20.e6        
        dqn=0.              
        tnadw10=6.           
        tnadw20=2.             
        dtnadw=0.              

c     surface temperature, surface interface number (fractional)
        thetas0=15. 
        thetas=thetas0      
c
c ------
c
c     length of year, integration time, time step, no steps, output times
        year=31557600.    			
	tmax=1.2e4*year    			
        dt=0.25e-2*year       			
        if (kappav.gt.1.e-5) dt=dt*1.e-5/kappav
        if (qek0.lt.15.e6) dt=0.1e-2*year
        nstop=int(tmax/dt)
c
c     output counter
        ndump=0
        ndumpp=0

c
c     define anthropogenic forcing time series: 
c       (a) antropogenic forcing; (b) AMOC collapse
        t1=tmax-3000.*year 
        t2=t1+200.*year  
        t3=t1+100.*year  
        t4=t1+200.*year
        n1=int(t1/dt)
        n2=int(t2/dt)
        n3=int(t3/dt)
        n4=int(t4/dt)
       
c ------


        do n=1,nstop

C     Forcing F(CO2)
           if (n.le.n1) then
              forcing(n)=init_for
           elseif (n.gt.n1.and.n.lt.n4) then
              forcing(n)=(n*dt-t1)/(t4-t1)
              forcing(n)=(sin(0.5*pi*forcing(n)))**2
              forcing(n)=init_for+amp*forcing(n)
           else
              forcing(n)=init_for+amp
           endif


c     MOC collapse, over 100 years (Meridional Overturning Circulation)
          if (n.le.n3) then
            famoc(n)=0.
          elseif (n.gt.n3.and.n.lt.n2) then  ! MOC transition (in sine^2)
            famoc(n)=(n*dt-t3)/(t2-t3)       ! from 100 to 200 years
            famoc(n)=(sin(0.5*pi*famoc(n)))**2 
          else 
            famoc(n)=1.
          endif

          !if (n.ge.n1) print*,forcing(n)
          
        enddo


c     constants required for heat content
        rho0=1027.
        cp=3992.

c     surface area north of ACC, length and width of ACC
        a=2.e14 	  
        lx=2.e7   
        ly=2.e6  


c ------

c     adams bashforth parameters (integration parameters, (3) multistep method)
        ab(1)=(23./12.)*dt/a   
        ab(2)=-(16./12.)*dt/a  
        ab(3)=(5./12.)*dt/a    
        do k=1,nz
           qtot(k,2)=0.         !Defined lower (=qtot(k,1))
           qtot(k,3)=0.         !Defined lower (=qtot(k,2))
        enddo
        
        rhs(2)=0.
        rhs(3)=0.
        rhso(2)=0.
        rhso(3)=0.

        print*,'step1: OK'
c ------

        do n=1,nstop

           rhs(1)=-lambda*Tatm+forcing(n)-gamma*(Tatm-To)
           Tatm=Tatm+a/c*(ab(1)*rhs(1)+ab(2)*rhs(2)+ab(3)*rhs(3))
           rhs(3)=rhs(2)
           rhs(2)=rhs(1)
           
           rhso(1)=gamma*(Tatm-To)
           To=To+a/co*(ab(1)*rhso(1)+ab(2)*rhso(2)+ab(3)*rhso(3))
           rhso(3)=rhso(2)
           rhso(2)=rhso(1)

           TT(n)=To

           if (n.ge.n1.and.mod(n-n1,(nstop-n1)/(nout-1)).eq.0) then
              ndump=ndump+1 
              time(ndump)=ndump
              Tout(ndump)=Tatm
              Toout(ndump)=TT(n)
              forcingout(ndump)=forcing(n)
              !print*,Tout(ndump),Toout(ndump)
           end if

         enddo

         TTmax=maxval(TT)+thetas0
         surface=real(nz)*(TTmax-thetas)/TTmax

c     temperature and initial depths of interfaces
        do k=0,nz
          theta(k)=TTmax*(real(nz-k)/real(nz))
          dtheta=TTmax/real(nz)               
          if (theta(k).ge.thetas) then         
            d(k)=0.0                        
          else
            d(k)=depth*(real(k)-surface)/(real(nz)-surface) 
          endif
          !print *,int(d(k)),k,theta(k)
        enddo
 
c     minimum upper layer thickness (numerical parameter)
        hmin=1.e-3

        ndump=0.
        ndumpp=0.

        print*,'step2: OK'

c -------------------------------------------------------------

c     main loop
     
      do n=1,nstop
                  
c       set thetas
          thetas=thetas0+TT(n)
          surface=real(nz)*(TTmax-thetas)/TTmax 

c       find surface layer and fraction not outcropped
          kmin=int(surface+1)        
    	  delta=real(kmin)-surface
          
c       calculate layer thicknesses
          do k=1,nz
            h(k)=max(d(k)-d(k-1),hmin) 
          end do

          if (n.eq.1) then
             print*,'step2.0: OK',TTmax
          endif


        do k=kmin,nz-1  	

          qn00=qn0*(1.0-dqn*famoc(n))   		!NADW  (control solution *(1-variation))
          tnadw1=tnadw10+dtnadw*famoc(n)       		!Temp of NADW formation (upper)
          tnadw2=tnadw20+dtnadw*famoc(n)   		!Temp of NADW formation (lower)
          if (theta(k).le.thetas.and.theta(k).gt.tnadw1) then  ! cf. (10) in the paper
            qn(k)=(thetas-theta(k))/(thetas-tnadw1)
            qn(k)=qn00*sin(0.5*pi*qn(k))
          elseif (theta(k).le.tnadw1.and.theta(k).gt.tnadw2) then 
            qn(k)=(tnadw1-theta(k))/(tnadw1-tnadw2)
            qn(k)=qn00*(cos(0.5*pi*qn(k)))**2
          else
            qn(k)=0.0
          endif 
c
          kappagm=kappagm0*(1.+deddy*fanth(n))   	! eddy diff
          qeddy(k)=kappagm*d(k)*(lx/ly)/		! cf. (6) and (7) in the paper
     -             max((thetas-theta(k))/thetas,1.e-4)	!? /max ?
          if (d(k).gt.ddrake) then 
            qeddy(k)=qeddy(k)*(depth-d(k))/(depth-ddrake)	
          endif
c     
c          print *,k,kappagm*d(k)*(lx/ly), 
c     -             max((thetas-theta(k))/thetas,1.e-3)
c
          qek00=qek0*(1.+dwind*fanth(n))   		! Ekman transport
          if (d(k).lt.dek) then 			! cf. (5) in the paper
            qek(k)=qek00*d(k)/dek
            qeddy(k)=qeddy(k)*d(k)/dek
          elseif (d(k).gt.ddrake) then 
            qek(k)=qek00*(depth-d(k))/(depth-ddrake)
          else
            qek(k)=qek00
          endif
c	
c
          if (k.eq.kmin) then                      ! Diapycanal upwelling
             qu(k)=a*kappav*(delta/h(k)-1./h(k+1)) ! cf. (9) in the paper
c     print *,delta/h(k),1./h(k+1)
          else
             qu(k)=a*kappav*(1/h(k)-1./h(k+1))
          endif
c
c
c	Equation (2) leading the model
          qtot(k,1)=qek(k)+qu(k)-qn(k)-qeddy(k)
c

       enddo
      
       
       if (n.eq.1) then
          print*,'step2.1: OK'
       endif
c
c
c       heat content and uptake (including decomposition)
c         subscripts 1-3 refer to >10, 5-10, <5 degrees
c         (NB: temperature classes hard-wired here for now) 
c
          heatek=0. 	!Ekman 
          heatek1=0.
          heatek2=0.
          heatek3=0.
          heatn=0.  	!NADW
          heatn1=0.
          heatn2=0.
          heatn3=0.
          heatu=0.  	!Diapyc upwelling
          heatu1=0.
          heatu2=0.
          heatu3=0.
          heateddy=0. 	!Eddy diff
          heateddy1=0.
          heateddy2=0.
          heateddy3=0.
          heatuptake=0. !uptake
          heatuptake1=0.
          heatuptake2=0.
          heatuptake3=0.
          heat=rho0*cp*dtheta*d(kmin) ! cf. (3) (heat over the surface layer (k=kmin))
c          print*,d(kmin)

          do k=kmin+1,nz-1
            heatek=heatek+rho0*cp*dtheta*qek(k)
            heatn=heatn+rho0*cp*dtheta*qn(k)
            heatu=heatu+rho0*cp*dtheta*qu(k)
            heateddy=heateddy+rho0*cp*dtheta*qeddy(k)
            heatuptake=heatuptake+rho0*cp*dtheta*qtot(k,1)
            heat=heat+rho0*cp*dtheta*a*d(k)
            if (theta(k).ge.10.0) then
              heatek1=heatek1+rho0*cp*dtheta*qek(k)
              heatn1=heatn1+rho0*cp*dtheta*qn(k)
              heatu1=heatu1+rho0*cp*dtheta*qu(k)
              heateddy1=heateddy1+rho0*cp*dtheta*qeddy(k)
              heatuptake1=heatuptake1+rho0*cp*dtheta*qtot(k,1)
            elseif (theta(k).ge.5.0.and.theta(k).lt.10.0) then
              heatek2=heatek2+rho0*cp*dtheta*qek(k)
              heatn2=heatn2+rho0*cp*dtheta*qn(k)
              heatu2=heatu2+rho0*cp*dtheta*qu(k)
              heateddy2=heateddy2+rho0*cp*dtheta*qeddy(k)
              heatuptake2=heatuptake2+rho0*cp*dtheta*qtot(k,1)
            else
              heatek3=heatek3+rho0*cp*dtheta*qek(k)
              heatn3=heatn3+rho0*cp*dtheta*qn(k)
              heatu3=heatu3+rho0*cp*dtheta*qu(k)
              heateddy3=heateddy3+rho0*cp*dtheta*qeddy(k)
              heatuptake3=heatuptake3+rho0*cp*dtheta*qtot(k,1)
            endif
          enddo
c

c     step forward (with 3-steps Adams-Bashforth linear explicit method) 
        do k=0,nz-1
           if (k.lt.kmin) then
              d(k)=0.
           elseif (k.ge.kmin) then
              d(k)=d(k)+ab(1)*qtot(k,1)+ab(2)*qtot(k,2)+ab(3)*qtot(k,3)
           endif
           if (k.eq.kmin) d(k)=max(d(k),hmin)
        enddo
c
        
        if (n.eq.1) then
           print*,'step2.2: OK'
        endif



        if (mod(n,(nstop/(nout-1))).eq.0) then
           ndumpp=ndumpp+1
           heattout(ndumpp)=heat
        endif

        if (n.ge.n1.and.mod(n-n1,(nstop-n1)/(nout-1)).eq.0) then

c         save data 

            ndump=ndump+1 	
            do k=0,nz
              if (k.lt.kmin.or.k.eq.nz) then
                qekout(k,ndump)=0. 
                qeddyout(k,ndump)=0. 
                qnout(k,ndump)=0. 
                quout(k,ndump)=0.
                qtotout(k,ndump)=0. 
                dout(k,ndump)=0.
              else 
                qekout(k,ndump)=qek(k)*1.e-6
                qeddyout(k,ndump)=qeddy(k)*1.e-6
                qnout(k,ndump)=qn(k)*1.e-6
                quout(k,ndump)=qu(k)*1.e-6
                qtotout(k,ndump)=qtot(k,1)*1.e-6
                dout(k,ndump)=d(k)       
              endif
              dout(nz,ndump)=d(nz)
              heatout(ndump)=heat
              heatupout(ndump)=heatuptake
              heatekout(ndump)=heatek
              heatuout(ndump)=heatu
              heatnout(ndump)=heatn
              heateddyout(ndump)=heateddy
              heatup1out(ndump)=heatuptake1
              heatek1out(ndump)=heatek1
              heatu1out(ndump)=heatu1
              heatn1out(ndump)=heatn1
              heateddy1out(ndump)=heateddy1
              heatup2out(ndump)=heatuptake2
              heatek2out(ndump)=heatek2
              heatu2out(ndump)=heatu2
              heatn2out(ndump)=heatn2
              heateddy2out(ndump)=heateddy2
              heatup3out(ndump)=heatuptake3
              heatek3out(ndump)=heatek3
              heatu3out(ndump)=heatu3
              heatn3out(ndump)=heatn3
              heateddy3out(ndump)=heateddy3
              thetasout(ndump)=thetas
            enddo
        endif
c
c       shuffle forcing
        do k=kmin,nz-1
          qtot(k,3)=qtot(k,2)
          qtot(k,2)=qtot(k,1)
        end do

      ! if (mod(n,40).eq.0) then
      !   print*,Tatm,To
      !endif


      end do 	  
c     

      print*,'step3: OK'

c     end of main (temporal) loop
c
c
c --------------------------------------------------------------------------
c
c
c     write data to files
c        
        open(unit=9,file='qek.dat',status='unknown') 
        do k=0,nz
         do n=1,ndump
          write(9,*) qekout(k,n)
         end do
        end do
        close(9)
c
        open(unit=9,file='qeddy.dat',status='unknown') 
        do k=0,nz
         do n=1,ndump
          write(9,*) qeddyout(k,n)
         end do
        end do
        close(9)
c
        open(unit=9,file='qn.dat',status='unknown') 
        do k=0,nz
         do n=1,ndump
          write(9,*) qnout(k,n)
         end do
        end do
        close(9)
c
        open(unit=9,file='qu.dat',status='unknown') 
        do k=0,nz
         do n=1,ndump
          write(9,*) quout(k,n)
         end do
        end do
        close(9)
c
        open(unit=9,file='qtot.dat',status='unknown') 
        do k=0,nz
         do n=1,ndump
          write(9,*) qtotout(k,n)
         end do
        end do
        close(9)
c
        open(unit=9,file='d.dat',status='unknown') 
        do k=0,nz
         do n=1,ndump
          write(9,*) dout(k,n)!
         end do
        end do
        close(9)
c
        open(unit=9,file='heat.dat',status='unknown') 
        do n=1,ndump
          write(9,*) heatout(n)
        end do
        close(9)
c
        open(unit=9,file='heatt.dat',status='unknown') 
        do n=1,ndump
          write(9,*) heattout(n)
        end do
        close(9)
c
        heatupout(0)=heatupout(1)
        heatupout(nout+1)=heatupout(nout)
        open(unit=9,file='heatup.dat',status='unknown') 
        do n=1,ndump
          write(9,*) 0.25*(2.*heatupout(n)
     -                   +heatupout(n+1)+heatupout(n-1))/a
        end do
        close(9)
c
        heatekout(0)=heatekout(1)
        heatekout(nout+1)=heatekout(nout)
        open(unit=9,file='heatek.dat',status='unknown') 
        do n=1,ndump
          write(9,*) 0.25*(2.*heatekout(n)
     -                   +heatekout(n+1)+heatekout(n-1))/a
        end do
        close(9)
c
        heatuout(0)=heatuout(1)
        heatuout(nout+1)=heatuout(nout)
        open(unit=9,file='heatu.dat',status='unknown') 
        do n=1,ndump
          write(9,*) 0.25*(2.*heatuout(n)
     -                   +heatuout(n+1)+heatuout(n-1))/a
        end do
        close(9)
c
        heateddyout(0)=heateddyout(1)
        heateddyout(nout+1)=heateddyout(nout)
        open(unit=9,file='heateddy.dat',status='unknown') 
        do n=1,ndump
          write(9,*) 0.25*(2.*heateddyout(n)
     -                   +heateddyout(n+1)+heateddyout(n-1))/a
        end do
        close(9)
c
        heatnout(0)=heatnout(1)
        heatnout(nout+1)=heatnout(nout)
        open(unit=9,file='heatn.dat',status='unknown') 
        do n=1,ndump
          write(9,*) 0.25*(2.*heatnout(n)
     -                   +heatnout(n+1)+heatnout(n-1))/a
        end do
        close(9)
c
        heatup2out(0)=heatup2out(1)
        heatup2out(nout+1)=heatup2out(nout)
        open(unit=9,file='heatup2.dat',status='unknown') 
        do n=1,ndump
          write(9,*) 0.25*(2.*heatup2out(n)
     -                   +heatup2out(n+1)+heatup2out(n-1))/a
        end do
        close(9)
c
        heatek2out(0)=heatek2out(1)
        heatek2out(nout+1)=heatek2out(nout)
        open(unit=9,file='heatek2.dat',status='unknown') 
        do n=1,ndump
          write(9,*) 0.25*(2.*heatek2out(n)
     -                   +heatek2out(n+1)+heatek2out(n-1))/a
        end do
        close(9)
c
        heatu2out(0)=heatu2out(1)
        heatu2out(nout+1)=heatu2out(nout)
        open(unit=9,file='heatu2.dat',status='unknown') 
        do n=1,ndump
          write(9,*) 0.25*(2.*heatu2out(n)
     -                   +heatu2out(n+1)+heatu2out(n-1))/a
        end do
        close(9)
c
        heateddy2out(0)=heateddy2out(1)
        heateddy2out(nout+1)=heateddy2out(nout)
        open(unit=9,file='heateddy2.dat',status='unknown') 
        do n=1,ndump
          write(9,*) 0.25*(2.*heateddy2out(n)
     -                   +heateddy2out(n+1)+heateddy2out(n-1))/a
        end do
        close(9)
c
        heatn2out(0)=heatn2out(1)
        heatn2out(nout+1)=heatn2out(nout)
        open(unit=9,file='heatn2.dat',status='unknown') 
        do n=1,ndump
          write(9,*) 0.25*(2.*heatn2out(n)
     -                   +heatn2out(n+1)+heatn2out(n-1))/a
        end do
        close(9)
c
        heatup3out(0)=heatup3out(1)
        heatup3out(nout+1)=heatup3out(nout)
        open(unit=9,file='heatup3.dat',status='unknown') 
        do n=1,ndump
          write(9,*) 0.25*(2.*heatup3out(n)
     -                   +heatup3out(n+1)+heatup3out(n-1))/a
        end do
        close(9)
c
        heatek3out(0)=heatek3out(1)
        heatek3out(nout+1)=heatek3out(nout)
        open(unit=9,file='heatek3.dat',status='unknown') 
        do n=1,ndump
          write(9,*) 0.25*(2.*heatek3out(n)
     -                   +heatek3out(n+1)+heatek3out(n-1))/a
        end do
        close(9)
c
        heatu3out(0)=heatu3out(1)
        heatu3out(nout+1)=heatu3out(nout)
        open(unit=9,file='heatu3.dat',status='unknown') 
        do n=1,ndump
          write(9,*) 0.25*(2.*heatu3out(n)
     -                   +heatu3out(n+1)+heatu3out(n-1))/a
        end do
        close(9)
c
        heateddy3out(0)=heateddy3out(1)
        heateddy3out(nout+1)=heateddy3out(nout)
        open(unit=9,file='heateddy3.dat',status='unknown') 
        do n=1,ndump
          write(9,*) 0.25*(2.*heateddy3out(n)
     -                   +heateddy3out(n+1)+heateddy3out(n-1))/a
        end do
        close(9)
c
        heatn3out(0)=heatn3out(1)
        heatn3out(nout+1)=heatn3out(nout)
        open(unit=9,file='heatn3.dat',status='unknown') 
        do n=1,ndump
          write(9,*) 0.25*(2.*heatn3out(n)
     -                   +heatn3out(n+1)+heatn3out(n-1))/a
        end do
        close(9)
c
        heatup1out(0)=heatup1out(1)
        heatup1out(nout+1)=heatup1out(nout)
        open(unit=9,file='heatup1.dat',status='unknown') 
        do n=1,ndump
          write(9,*) 0.25*(2.*heatup1out(n)
     -                   +heatup1out(n+1)+heatup1out(n-1))/a
        end do
        close(9)
c
        heatek1out(0)=heatek1out(1)
        heatek1out(nout+1)=heatek1out(nout)
        open(unit=9,file='heatek1.dat',status='unknown') 
        do n=1,ndump
          write(9,*) 0.25*(2.*heatek1out(n)
     -                   +heatek1out(n+1)+heatek1out(n-1))/a
        end do
        close(9)
c
        heatu1out(0)=heatu1out(1)
        heatu1out(nout+1)=heatu1out(nout)
        open(unit=9,file='heatu1.dat',status='unknown') 
        do n=1,ndump
          write(9,*) 0.25*(2.*heatu1out(n)
     -                   +heatu1out(n+1)+heatu1out(n-1))/a
        end do
        close(9)
c
        heateddy1out(0)=heateddy1out(1)
        heateddy1out(nout+1)=heateddy1out(nout)
        open(unit=9,file='heateddy1.dat',status='unknown') 
        do n=1,ndump
          write(9,*) 0.25*(2.*heateddy1out(n)
     -                   +heateddy1out(n+1)+heateddy1out(n-1))/a
        end do
        close(9)
c
        heatn1out(0)=heatn1out(1)
        heatn1out(nout+1)=heatn1out(nout)
        open(unit=9,file='heatn1.dat',status='unknown') 
        do n=1,ndump
          write(9,*) 0.25*(2.*heatn1out(n)
     -                   +heatn1out(n+1)+heatn1out(n-1))/a
        end do
        close(9)
c
        open(unit=9,file='thetas.dat',status='unknown') 
        do n=1,ndump
          write(9,*) thetasout(n)
        end do
        close(9)     
c     
        open(unit=9,file='forcing.dat',status='unknown') 
        do n=1,ndump
          write(9,*) forcingout(n)
        end do
        close(9)  
c
        open(unit=9,file='Tatm.dat',status='unknown') 
        do n=1,ndump
          write(9,*) Tout(n)
        end do
        close(9)  
c
        open(unit=9,file='To.dat',status='unknown') 
        do n=1,ndump
          write(9,*) Toout(n)
        end do
        close(9)  
c
        open(unit=9,file='Param.dat',status='unknown') 
        write(9,*) c,co,lambda,gamma,amp
        close(9)  
c
        stop 
        end

