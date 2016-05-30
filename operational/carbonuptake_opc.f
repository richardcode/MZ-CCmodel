      program Ocean Carbon Uptake
    
c     September 2013
c     - A program designed to model the Atlantic Ocean's uptake and
c     circulation of carbon.

      implicit none
      integer nz,nt,nout,nmax
      parameter(nz=100,nt=3,nout=1001,nmax=100000001)

      real*8 dc(0:nz),h(nz),z(nz),d(0:nz),theta(0:nz),w(0:nz),phi0(0:nz)
      real*8 qek(nz),qeddy(nz),qn(nz),qu(nz),pco2(nmax),dDICdt(0:nz,nt)
      real*8 qtot(nz,nt),phitot(nz,nt),ab(nt),fanth(nmax),famoc(nmax)
      real*8 dout(0:nz,nout),dcout(0:nz,nout),g(0:nz,nt),f(0:nz)
      real*8 qekout(0:nz,nout),qeddyout(0:nz,nout),DIC(0:nz)
      real*8 qnout(0:nz,nout),quout(0:nz,nout),qtotout(0:nz,nout)
      real*8 heatout(nout),heatupout(0:nout+1),thetasout(nout)
      real*8 heatekout(0:nout+1),heatuout(0:nout+1),zz(nz),phiMOC
      real*8 heateddyout(0:nout+1),heatnout(0:nout+1),phiek,phied
      real*8 pi,year,tmax,dt,t1,t2,t3,depth,dek,qek0,ddrake,socean
      real*8 theta0,thetas,delta,surface,tnadw10,tnadw20,a,pco20
      real*8 dwind,qek00,kappagm0,deddy,dtnadw,tnadw1,tnadw2,dc0
      real*8 thetas0,dthetas,dtheta,rho0,cp,heat,heatuptake
      real*8 hmin,kappav,kappagm,lx,ly,qn0,dqn,qn00,phi00
      real*8 heatn,heatu,heateddy,heatek,DIC0,phidif,phi4(nz)
      real*8 heatn1,heatu1,heateddy1,heatek1,heatuptake1
      real*8 heatn2,heatu2,heateddy2,heatek2,heatuptake2
      real*8 heatn3,heatu3,heateddy3,heatek3,heatuptake3
      real*8 heatup1out(0:nout+1),phi1(nz),phi2(nz),phi3(nz)
      real*8 heatek1out(0:nout+1),heatu1out(0:nout+1)
      real*8 heateddy1out(0:nout+1),heatn1out(0:nout+1)
      real*8 heatup2out(0:nout+1),Beta,Alpha
      real*8 heatek2out(0:nout+1),heatu2out(0:nout+1)
      real*8 heateddy2out(0:nout+1),heatn2out(0:nout+1)
      real*8 heatup3out(0:nout+1),pco2out(0:nout+1)
      real*8 heatek3out(0:nout+1),heatu3out(0:nout+1)
      real*8 heateddy3out(0:nout+1),heatn3out(0:nout+1)
		
	  integer nstop,n,k,j,kmin,n1,n2,n3,ndump

c     Constants required for heat content
        rho0=1027.
        cp=3992.
	pco20=28.35	! Carbon pressure
        dc0=9.74e-5	! Carbon anomaly
        phi00=5.e-6	! air-sea flux
	DIC0=2.e-3	! initial DIC

      pi=3.14159265358979

c     Southern ocean switch (0 for no Drake Passage, 1 for Drake Passage)
        socean=1.
c
c     Eddy diffusivities (placed here to scale dt with kappav)
c            fractional change in eddy transport (0.1 = +10%)
        kappav=1.e-5
        kappagm0=1.e3*socean
        deddy=0.
c
c     Southern Ocean Ekman transport, Ekman depth, Drake Passage depth,
c            fractional change in Ekman transport (0.1 = +10%)
        qek0=30.e6*socean
        dek=100.
        depth=5.e3
        ddrake=0.8*depth
        dwind=0.
c
c     NADW formation, fractional reduction in NADW formation (0.5 = -50%),
c            temperature range of NADW formation, NADW warming 
        qn0=20.e6
        dqn=0.
        tnadw10=6.
        tnadw20=2.
        dtnadw=0.
c                
c     Surface temperature, surface interface number (fractional)
c       Note: have to modify theta0 if dthetas > 5
c         (theta0 hard-wired to keep levels the same across each scenario;
c          small number added to ensure data output times do not coincide 
c          with a new layer appearing at the surface => noisy diagnostics!)
        thetas0=15.
        dthetas=0.
        theta0=15.0+1.45e-3
        thetas=thetas0
        surface=real(nz)*(theta0-thetas)/theta0	  
c
c     Length of year, integration time, time step, no steps, output times
	    year=31557600.
	    tmax=1.1e4*year
        dt=0.25e-2*year
        if (kappav.gt.1.e-5) dt=dt*1.e-5/kappav
        if (qek0.lt.15.e6) dt=0.1e-2*year
        nstop=int(tmax/dt)
c
c     Output counter
        ndump=0
c
c     Define anthropogenic forcing time series: 
c       (a) anthropogenic forcing; (b) AMOC collapse;
        t1=tmax-1000.*year
c       Set to 140 for ramp to 4xCO2
        t2=t1+140.*year
        t3=t1+100.*year
        n1=int(t1/dt)
        n2=int(t2/dt)
        n3=int(t3/dt)
        do n=1,nstop
c          if (n.le.n1) then
            fanth(n)=0.
c          elseif (n.gt.n1.and.n.lt.n2) then 
c            fanth(n)=(n*dt-t1)/(t2-t1)
c           fanth(n)=(sin(0.5*pi*fanth(n)))**2
c          else 
c            fanth(n)=1.
c          endif
c          if (n.le.n3) then
            famoc(n)=0.
c          elseif (n.gt.n3.and.n.lt.n2) then 
c            famoc(n)=(n*dt-t3)/(t2-t3)
c            famoc(n)=(sin(0.5*pi*famoc(n)))**2
c          else 
c           famoc(n)=1.
c         endif
        enddo

c     Surface area north of ACC, length and width of ACC
        a=2.e14	  
        lx=2.e7
        ly=2.e6
c
c     Adams Bashforth parameters
  	    ab(1)=(23./12.)*dt/a
	    ab(2)=-(16./12.)*dt/a
        ab(3)=(5./12.)*dt/a
        do k=1,nz
          qtot(k,2)=0.
          qtot(k,3)=0.
        enddo
c
c     Temperature and initial depths of interfaces
        do k=0,nz
          theta(k)=theta0*(real(nz-k)/real(nz))
          dtheta=theta0/real(nz)
          if (theta(k).ge.thetas) then 
            d(k)=0.0
          else
            d(k)=depth*(real(k)-surface)/(real(nz)-surface)
          endif
c          print *,d(k),k,theta(k)
        enddo

	kmin=int(surface+1)

!     Initial air-Sea fluxes	
	do k=0,nz
        if (k.lt.kmin) then
          phi0(k)=phi00*(real(kmin-k))	!linear distribution over the 'surface' first layers. 
					!(to be modified? -> considering the outcropping surface)
	else
	  phi0(k)=0.0
	endif
	enddo

!     Initial amounts of carbon anomaly
	do k=0,nz
          if (k.lt.kmin) then
            dc(k)=0.0
	  else
	    dc(k)=dc0*(real(k-kmin))/real(nz-kmin)
	  endif
	enddo

!     Initial amounts of DIC
	do k=0,nz
	DIC(k)=DIC0+400*(real(k)/100.)
	enddo

!     pCO2 forcing
      print *, 'n1 = ',n1
      print *, 'n2 = ',n2
        do n=0,nstop
          if (n.le.n1) then
            pco2(n)=28.35
          elseif (n.gt.n1.and.n.lt.n2) then 
c            pco2(n)=((n*dt-t1)/(t2-t1))
c            pco2(n)=28.35+(pco20-28.35)*(sin(0.5*pi*pco2(n)))**2
c            Increase partial pressure of CO2 by 1%/yr to 4xCO2 stabilisation
             pco2(n)=28.35*exp(log(1.01)*(n*dt-t1)/year)
          elseif (n.ge.n2) then
            pco2(n)=28.35*exp(log(1.01)*(n2*dt-t1)/year)
          endif
        enddo
c      print *, pco2

c     Minimum upper layer thickness (numerical parameter)
        hmin=1.e-3

c     Main loop
      print *, nstop
c
      do n=1,nstop
c
c       Set thetas
          thetas=thetas0+dthetas*fanth(n)
          surface=real(nz)*(theta0-thetas)/theta0
c
c       Find surface layer and fraction not outcropped
          kmin=int(surface+1)
    	  delta=real(kmin)-surface
c          
c       Calculate layer thicknesses
          do k=1,nz
            h(k)=max(d(k)-d(k-1),hmin)
          end do
c	
          do k=kmin,nz-1
c
          qn00=qn0*(1.0-dqn*famoc(n))
          tnadw1=tnadw10+dtnadw*famoc(n)
          tnadw2=tnadw20+dtnadw*famoc(n)
          if (theta(k).le.thetas.and.theta(k).gt.tnadw1) then 
            qn(k)=(thetas-theta(k))/(thetas-tnadw1)
            qn(k)=qn00*sin(0.5*pi*qn(k))
          elseif (theta(k).le.tnadw1.and.theta(k).gt.tnadw2) then 
            qn(k)=(tnadw1-theta(k))/(tnadw1-tnadw2)
            qn(k)=qn00*(cos(0.5*pi*qn(k)))**2
          else
            qn(k)=0.0
          endif 
c
          kappagm=kappagm0*(1.+deddy*fanth(n))
          qeddy(k)=kappagm*d(k)*(lx/ly)/
     -             max((thetas-theta(k))/thetas,1.e-4)
          if (d(k).gt.ddrake) then 
            qeddy(k)=qeddy(k)*(depth-d(k))/(depth-ddrake)
          endif
c     
          qek00=qek0*(1.+dwind*fanth(n))
          if (d(k).lt.dek) then 
            qek(k)=qek00*d(k)/dek
            qeddy(k)=qeddy(k)*d(k)/dek
          elseif (d(k).gt.ddrake) then 
            qek(k)=qek00*(depth-d(k))/(depth-ddrake)
          else
            qek(k)=qek00
          endif
c
          if (k.eq.kmin) then
            qu(k)=a*kappav*(delta/h(k)-1./h(k+1))
          else
            qu(k)=a*kappav*(1/h(k)-1./h(k+1))
          endif
c
          qtot(k,1)=qek(k)+qu(k)-qn(k)-qeddy(k)
c
c          if ((n.le.n1.and.mod(n,nstop/10).eq.0).or.
c     -        (n.ge.n1.and.mod(n,nstop/500).eq.0)) then 
c            if (k.eq.kmin) print *,'time:',real(n)/real(nstop)
c           print *,theta(k),d(k),1.e-6*qek(k),-1.e-6*qeddy(k),
c   -            -1.e-6*qn(k),1.e-6*qu(k),1.e-6*qtot(k,1)

c            if (k.eq.nz-1) print *,'heat flux:',heatuptake/a
c            if (k.eq.nz-1) print *,' '
c
        enddo
c
c       Heat content and uptake (including decomposition)
c         subscripts 1-3 refer to >10, 5-10, <5 degrees
c         (NB: temperature classes hard-wired here for now) 
c
          heatek=0.
          heatek1=0.
          heatek2=0.
          heatek3=0.
          heatn=0.
          heatn1=0.
          heatn2=0.
          heatn3=0.
          heatu=0.
          heatu1=0.
          heatu2=0.
          heatu3=0.
          heateddy=0.
          heateddy1=0.
          heateddy2=0.
          heateddy3=0.
          heatuptake=0.
          heatuptake1=0.
          heatuptake2=0.
          heatuptake3=0.
          heat=rho0*cp*dtheta*d(kmin)
c
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


!       The new additional algebra in the diffusion equation!
       
	do j=1,3
	do k=kmin,nz-1
	  dDICdt(k,j)=(1./101325.)*(1.54818e-4)*pco2(n)*thetas*qtot(k,j)/a	!Pressure = 101325 Pa, 
          dDICdt(k,j)=dDICdt(k,j)+18.3*(real(pco2(n)-pco2(n-1))/dt)
	  g(k,j)=((dc(k-2)-dc(k))/(h(k)+h(k-2)))
	  g(k,j)=g(k,j)-((dc(k)-dc(k+2))/(h(k)+h(k+2)))
	  g(k,j)=kappav*g(k,j)
	  g(k,j)=g(k,j)+0.5*(h(k+1)+h(k-1))*dDICdt(k,j)
          g(k,j)=g(k,j)-(phi0(k+1)-phi0(k-1))
	  g(k,j)=g(k,j)*dc(k)
	  f(k)=d(k)*(dc(k-1)-dc(k+1))-0.5*(dc(k)*(h(k+1)+h(k-1)))
        enddo
	enddo

!     Define phitot
	do j=1,3
	do k=kmin,nz-1
 	phitot(k,j)=-real((qu(k)*dc(k))*(dc(k-1)-dc(k+1)))+(g(k,j)/f(k))
       	enddo
	enddo

!     Time-Stepping:
!     a) DIC; b) Carbon anomaly; c) Interface Depth;
	do k=kmin,nz-1
	  DIC(k)=DIC(k)+ab(1)*a*dDICdt(k,1)+ab(2)*a*dDICdt(k,2)
	  DIC(k)=DIC(k)+ab(3)*a*dDICdt(k,3)
	enddo
    	
	do k=kmin,nz-1
          dc(k)=dc(k)+ab(1)*phitot(k,1)+ab(2)*phitot(k,2)
          dc(k)=dc(k)+ab(3)*phitot(k,3)
	enddo

	do k=kmin,nz-1
          d(k)=d(k)+ab(1)*qtot(k,1)+ab(2)*qtot(k,2)+ab(3)*qtot(k,3)
        if (k.eq.kmin) d(k)=max(d(k),hmin)
        enddo
   
       if (n.ge.n1.and.mod(n-n1,(nstop-n1)/(nout-1)).eq.0) then

c         Save data 
            ndump=ndump+1 
c           print *,n,n1,n-n1,nstop-n1,nout-1,(nstop-n1)/(nout-1)
c           print *,ndump,nout,n,nstop
            do k=0,nz
              if (k.lt.kmin.or.k.eq.nz) then
                qekout(k,ndump)=0. 
                qeddyout(k,ndump)=0. 
                qnout(k,ndump)=0. 
                quout(k,ndump)=0.
                qtotout(k,ndump)=0. 
                dout(k,ndump)=0.
                dcout(k,ndump)=0.
              else 
                qekout(k,ndump)=qek(k)*1.e-6
                qeddyout(k,ndump)=qeddy(k)*1.e-6
                qnout(k,ndump)=qn(k)*1.e-6
                quout(k,ndump)=qu(k)*1.e-6
                qtotout(k,ndump)=qtot(k,1)*1.e-6
                dout(k,ndump)=d(k)
		dcout(k,ndump)=dc(k)      
              endif
              dout(nz,ndump)=d(nz)
	      dcout(nz,ndump)=dc(nz)
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
c             Add the forcing pCO2 to the output data
              pco2out(ndump) = pco2(n)
	enddo
	endif

c
c     Shuffle forcing	
        do k=kmin,nz-1
          qtot(k,3)=qtot(k,2)
          qtot(k,2)=qtot(k,1)
        end do
c	print*,f(50)

c     Preliminary Data Writing!
	do k=1,nz
	z(k) = h(k)*dc(k)
	zz(k) = h(k)*DIC(k)
c	phi1(k)=(rho0*qek(k)*h(k)*dc(k))/(5000.)
c	phi2(k)=(rho0*qeddy(k)*h(k)*dc(k))/(5000.)
c	phi3(k)=(rho0*qn(k)*h(k)*dc(k))/(5000.)
c	phi4(k)=rho0*d(k)*h(k)*(g(k,j)/f(k))/(5000.)
	enddo
	Beta = 0.
	Alpha = 0.
c	phiek = 0.
c	phied = 0.
c	phiMOC = 0.
c	phidif = 0.
	do k=1,99
	Beta = Beta + z(k)
	Alpha = Alpha + zz(k)
c	phiek = phiek + phi1(k)
c	phied = phied + phi2(k)
c	phiMOC = phiMOC + phi3(k)
c	phidif = phidif + phi4(k)
	enddo
	open(unit=9,file='dc50.dat',status='unknown')        
	write(9,*) Alpha/(5000.)

        
      enddo
      close(9)	  
c
c     End of main loop
c
c     Write data to files
c
        print *, 'ndump =',ndump
        open(unit=9,file='qek.dat',status='unknown') 
        do k=0,nz
         do n=1,ndump
          write(9,*) qekout(k,n)
         enddo
        enddo
        close(9)
c
        open(unit=9,file='qeddy.dat',status='unknown') 
        do k=0,nz
         do n=1,ndump
          write(9,*) qeddyout(k,n)
         enddo
        enddo
        close(9)
c
        open(unit=9,file='qn.dat',status='unknown') 
        do k=0,nz
         do n=1,ndump
          write(9,*) qnout(k,n)
         enddo
        enddo
        close(9)
c
        open(unit=9,file='qu.dat',status='unknown') 
        do k=0,nz
         do n=1,ndump
          write(9,*) quout(k,n)
         enddo
        enddo
        close(9)
c
        open(unit=9,file='qtot.dat',status='unknown') 
        do k=0,nz
         do n=1,ndump
          write(9,*) qtotout(k,n)
         enddo
        enddo
        close(9)
c
        open(unit=9,file='d.dat',status='unknown') 
         do n=1,ndump
          write(9,*) dout(k,n)
         enddo
        close(9)
c
        open(unit=9,file='heat.dat',status='unknown') 
        do n=1,ndump
          write(9,*) heatout(n)
        enddo
        close(9)
c
        open(unit=9,file='pco2.dat',status='unknown')
        do n=1,ndump
          write(9,*) pco2out(n)
        enddo
        close(9)
c
        open(unit=9,file='carbonanomaly.dat',status='unknown')
        do k=0,nz
	  do n=1,ndump-1
          write(9,"(d20.14,2x)",advance="no") dcout(k,n)
        enddo
	write(9,"(d20.14,2x)") dcout(k,n)
        enddo
        close(9)
c
      stop 
      end
