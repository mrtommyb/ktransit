CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine transitmodel(nfit,nplanet,nplanetmax,sol,nmax,npt,time,
     .  itime,ntt,tobs,omc,tmodel,dtype)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
C Here are the f2py-specific comments.
Cf2py intent(in) :: nfit,nplanet,nplanetmax,sol,nmax,npt,time,itime,ntt
Cf2py intent(in) :: tobs,omc,dtype
Cf2py intent(out) :: tmodel

      integer nfit,npt,i,j,nintg,dtype(npt),ii,nplanet,nplanetmax,nmax,
     .   caltran
      parameter(nintg=11)
      double precision sol(nfit),per,epoch,b,RpRs,tmodel(npt),
     .  time(npt),phi,adrs,bt(nintg),bs2,Pi,tpi,c1,c2,c3,c4,
     .  itime(npt),t,tflux(nintg),dnintg,dnintgm1,pid2,zpt,eccn,w,Eanom,
     .  Tanom,trueanomaly,phi0,vt(nintg),K,voff,drs,distance,dil,G,
     .  esw,ecw,Manom,Ted,ell,Ag,Cs,fDB,tide(nintg),alb(nintg),
     .  albedomod,phase,ratio,ab,tm,mu(nintg),y2,x2,incl,mulimbf(nintg),
     .  occ(nintg),bp(nintg)
      integer ntt(nplanetmax)
      double precision tobs(nplanetmax,nmax),omc(nplanetmax,nmax),ttcor

      Pi=acos(-1.d0)!define Pi and 2*Pi
      tPi=2.0d0*Pi
      pid2=Pi/2.0d0
      G=6.674d-11 !N m^2 kg^-2  Gravitation constant
      Cs=2.99792458e8 !Speed of light
      fDB=5.46 !Doppler Boosting factor

      c1=sol(2)      !non-linear limb-darkening
      c2=sol(3)
      c3=sol(4)
      c4=sol(5)
      dil=sol(6)     !dilution parameter (model scaling)
      voff=sol(7)    !velocity zero point
      zpt=sol(8)     !flux zero point.

      do 17 i=1,npt
        tmodel(i)=0.0d0
 17   continue

      do 16 ii=1,nplanet

        per=sol(10*(ii-1)+8+2)     !Period (days)
c        sol(10*(ii-1)+8+3)=abs(sol(10*(ii-1)+8+3)) !make b positive
c        bs2=abs(sol(10*(ii-1)+8+3))
        bs2=sol(10*(ii-1)+8+3)*sol(10*(ii-1)+8+3)
        b=sqrt(bs2)       !impact parameter
c        write(0,*) bs2,b
        RpRs=abs(sol(10*(ii-1)+8+4))    !Rp/R*

        ecw=sol(10*(ii-1)+8+6)
        esw=sol(10*(ii-1)+8+5)
        eccn=sqrt(ecw*ecw+esw*esw) !eccentricity
        if(eccn.ge.1.0) eccn=0.99
        if(eccn.eq.0.0d0)then
            w=0.0d0
        else
            if(ecw.eq.0.0d0)then
                w=Pi/2.0d0
            else
                w=atan(esw/ecw)
            endif
            if((ecw.gt.0.0d0).and.(esw.lt.0.0d0))then
                w=tPi+w
            elseif((ecw.lt.0.0d0).and.(esw.ge.0.0d0))then
                w=Pi+w
            elseif((ecw.le.0.0d0).and.(esw.lt.0.0d0))then
                w=Pi+w
            endif
        endif

c        write(0,*) sol(7),sol(8),w
c        write(0,*) "w:",acos(sol(7)/eccn),asin(sol(8)/eccn)
c        read(5,*)

C       a/R*
c        adrs=sol(5)*per/tpi*sqrt(1-sol(3))*(1+sol(8))/sqrt(1-eccn*eccn)
        adrs=1000.0*sol(1)*G*(Per*86400.0d0)**2/(3.0d0*Pi)
        adrs=adrs**(1.0d0/3.0d0)
c        write(0,*) "a/R*:",adrs

C       Find inclination
        incl=acos(b/adrs)


c        K=abs(sol(10*(ii-1)+8+7))
        K=sol(10*(ii-1)+8+7)

        ted=sol(10*(ii-1)+8+8)/1.0d6 !Occultation Depth
        ell=sol(10*(ii-1)+8+9)/1.0d6 !Ellipsoidal variations
        ag=sol(10*(ii-1)+8+10)/1.0d6 !Phase changes

        dnintg=dble(nintg) !convert integer to double
        dnintgm1=2.0*dnintg-2.0

C     Find phase at centre of transit
        epoch=sol(10*(ii-1)+8+1)   !center of transit time (days)
c        phi1=(epoch/per-int(epoch/per))*twopi
        Eanom=tan(w/2.0d0)/sqrt((1.0d0+eccn)/(1.0d0-eccn)) !mean anomaly
        Eanom=2.0d0*atan(Eanom)
        phi0=Eanom-eccn*sin(Eanom)

        do 10 i=1,npt
            call lininterp(tobs,omc,nplanetmax,nmax,ii,ntt,time(i),
     .          ttcor)
c            write(0,*) ii,time(i),ttcor
c            read(5,*)
            do 11 j=1,nintg
                tflux(j)=0.0 !initialize model
C               sample over integration time

c                t=t-ttcor

                t=time(i)+itime(i)*(2.0*dble(j)-dnintg-1.0)/dnintgm1-
     .              epoch-ttcor


c                write(0,*) itime(i)
C               get orbital position (mean anomaly)
                phi=t/per-floor(t/per)
                phi=phi*tPi+phi0
                Manom=phi
                if(Manom.gt.tPi) Manom=Manom-tPi
                if(Manom.lt.0.0d0) Manom=Manom+tPi
                call kepler(Manom,Eanom,eccn)
                Tanom=trueanomaly(eccn,Eanom)
                if(phi.gt.Pi) phi=phi-tPi
                drs=distance(adrs,eccn,Tanom)
C               Added this (2014/04/23)
                incl=acos(b/drs)
                x2=drs*Sin(Tanom-w)
                y2=drs*Cos(Tanom-w)*cos(incl)

c                x2=drs*Cos(Tanom+w)
c                y2=drs*Sin(Tanom+w)*cos(incl)


c                bt(j)=sqrt(bs2+(drs*sin(Tanom-w))**2)
                bt(j)=sqrt(x2*x2+y2*y2)
C               Correct for light-travel time!
c            if((abs(bt(j))-RpRs.le.1.0d0).and.(abs(phi).gt.Pid2))then
c              t=time(i)-ltt+itime(i)*(2.0*dble(j)-dnintg-1.0)/dnintgm1
c     .              -epoch
c              phi=(t-ltt)/per-floor((t-ltt)/per)
c              phi=phi*tPi
c              Manom=phi+w
c              if(Manom.gt.tPi) Manom=Manom-tPi
c              if(Manom.lt.0.0d0) Manom=Manom+tPi
c              call kepler(Manom,Eanom,eccn)
c              Tanom=trueanomaly(eccn,Eanom)
c                if(phi.gt.Pi) phi=phi-tPi
c                drs=distance(adrs,eccn,Tanom)
c                bt(j)=sqrt(bs2+(drs*sin(Tanom-phi0))**2)
c            endif
c                vt(j)=K*(cos(Pid2+Tanom-phi0)+eccn*cos(w))

                vt(j)=K*(cos(Tanom-w+pid2)+eccn*cos(-w+pid2))
c                vt(j)=-K*(cos(Tanom+w)+ecw)

                tide(j)=ell*(drs/adrs)**(1.0d0/3.0d0)*
     .              cos(2.0d0*(Pid2+Tanom-w))
c                tide(j)=ell*(drs/adrs)**(1.0d0/3.0d0)*
c     .              cos(2.0d0*(Pid2+phi))

                alb(j)=albedomod(Pi,ag,Tanom-w)*adrs/drs
c                alb(j)=albedomod(Pi,ag,phi)*adrs/drs

                if(j.eq.nintg/2+1)then
c                    phase=Tanom-w!phi(nintg/2+1)
c                    if(phase.gt.Pi) phase=phase-tPi
c                    if(phase.lt.-Pi) phase=phase+tPi
c                    write(6,*) time(i),x2/adrs,y2/adrs/cos(incl)
c                    write(6,*) time(i),vt(j),w
c                     write(6,*) time(i),vt(j),w
                endif

 11         continue
            if(dtype(i).eq.0)then
c                if(abs(phase).lt.Pid2)then
                if(y2.ge.0.0d0)then
C       If we have a transit
                    caltran=0 !if zero, there is no transit
                    do 18 j=1,nintg
                        if(bt(j).lt.1.0d0+RpRs)then
                           caltran=1
                        endif
 18                 continue
                    if(caltran.eq.1) then
                       if((c3.eq.0.0).and.(c4.eq.0.0))then
                         call occultquad(bt,c1,c2,RpRs,tflux,mu,nintg)
                       else
                       call occultsmall(RpRs,c1,c2,c3,c4,nintg,bt,tflux)
c                         call occultnl(RpRs,c1,c2,c3,c4,bt,tflux,
c     .                     mulimbf,nintg)
c                         write(0,*) RpRs,abs(sol(10*(ii-1)+8+4))
c                         write(6,550) RpRs,(bt(j),tflux(j),j=1,nintg)
 550                     format(30(F8.5,1X))
                       endif
                    else
                        do 19 j=1,nintg
                           tflux(j)=1.0d0
 19                     continue
                    endif
                    tm=0.0d0
                    do 12 j=1,nintg
                        if(RpRs.le.0.0)tflux(j)=1.0d0
C                   model=transit+doppler+ellipsodial
                        tm=tm+tflux(j)-fDB*vt(j)/Cs+tide(j)+alb(j)
 12                 continue
                    tm=tm/dnintg
                else
C       We have an eclipse
                    tm=0.0d0
                    do 20 j=1,nintg
                      bp(j)=bt(j)/RpRs
 20                 continue
                    call occultuniform(bp,1.0/RpRs,occ,nintg)
                    do 14 j=1,nintg
                        ratio=1.0d0-occ(j)

C                      Old estimate, replaced by analytic function
c                        ratio=1.0d0
c                        ab=dabs(bt(j))
c                        if((ab.ge.1.0d0).and.(ab-RpRs.le.1.0d0))then
c                            ratio=(1.0d0+RpRs-ab)/(2.0d0*RpRs)
c                        elseif((ab.lt.1.0d0).and.(ab+RpRs.ge.1.0d0))then
c                            ratio=(RpRs+1.0d0-ab)/(2.0d0*RpRs)
c                        elseif(ab-RpRs.gt.1.0d0)then
c                            ratio=0.0d0
c                        endif
c                        write(0,*) bt(j),ratio,rationew
c                        read(5,*)
                        if(RpRs.le.0.0d0) ratio=0.0d0
                        tm=tm+(1.0d0-ted*ratio)
     .                      -fDB*vt(j)/Cs+tide(j)+alb(j)
 14                 continue
                    tm=tm/dnintg
                endif
                tm=tm+(1.0d0-tm)*dil-1.0d0!add dilution
            else
                tm=0.0d0
                do 13 j=1,nintg
                    tm=tm+vt(j)
 13             continue
                tm=tm/dnintg
c            write(0,*) "rv:",tmodel(i)
c            read(5,*)
            endif
            tmodel(i)=tmodel(i)+tm
 10     continue

C     Need to add zero points (voff, zpt)

c        do 9 i=1,npt
c            write(6,*) time(i),tmodel(i)
c 9      continue

 16   continue

      do 15 i=1,npt
        if(dtype(i).eq.0)then
            tmodel(i)=tmodel(i)+zpt+1.0d0
        else
            tmodel(i)=tmodel(i)+voff
        endif
 15   continue

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      double precision function albedomod(Pi,ag,phi)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      double precision Pi,phi,alpha,phase,ag

      phi=phi+Pi
      if(phi.gt.2.0*Pi) phi=phi-2.0*Pi


      alpha=abs(phi)
c      alpha=2.0*Pi*t/Per+phi
      alpha=alpha-2.0*Pi*int(alpha/(2.0*Pi))
      if(alpha.gt.Pi) alpha=abs(alpha-2.0*pi)
c      write(6,*) t,alpha
c      phase=(1.0d0+cos(alpha))/2.0d0
      phase=(sin(alpha)+(Pi-alpha)*cos(alpha))/Pi  !Lambertian Sphere

      albedomod=ag*phase

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      double precision function distance(asep,eccn,Tanom)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      double precision asep,eccn,Tanom

      distance=asep*(1.0d0-eccn*eccn)/(1+eccn*cos(Tanom))

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      double precision function trueanomaly(eccn,Eanom)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      double precision eccn,Eanom,temp(2)

      temp(1)=sqrt((1.0d0+eccn)/(1.0d0-eccn))
      temp(2)=tan(Eanom/2.)
      trueanomaly=2.0d0*atan(temp(1)*temp(2))

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine kepler(Manom,Eanom,eccn)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer i,itmax
      parameter(itmax=100)
      double precision Manom,Eanom,Eold,eccn,diff,thres
      thres=1.0d-8

      Eold=Eanom
      Eanom=Manom+eccn*sin(Eanom)
      diff=abs(1.0d0-Eanom/Eold)
      Eold=Eanom
      i=0
      do while ((diff.gt.thres).and.(i.lt.itmax))
        Eanom=Manom+eccn*sin(Eanom)
        diff=abs(1.0d0-Eanom/Eold)
        Eold=Eanom
        i=i+1
      enddo

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine invkepler(Eanom,Manom,eccn)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer i,itmax
      parameter(itmax=100)
      double precision Manom,Eanom,Mold,eccn,diff,thres
      thres=1.0d-6

      Mold=Manom
      Manom=Eanom-eccn*sin(Manom)
      diff=abs(1.0d0-Manom/Mold)
      Mold=Manom
      i=0
      do while ((diff.gt.thres).and.(i.lt.itmax))
        Manom=Eanom-eccn*sin(Manom)
        diff=abs(1.0d0-Manom/Mold)
        Mold=Manom
        i=i+1
      enddo
c      if(i.ge.itmax) write(0,*) "invkepler itmax"

      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine readttfile(nunit,nplanetmax,nmax,nplanet,ntt,tobs,omc)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nunit,ntt(nplanet),i,nmax,nplanetmax,nplanet
      double precision tobs(nplanetmax,nmax),omc(nplanetmax,nmax),err

      i=1
 10   read(nunit,*,end=11) tobs(nplanet,i),omc(nplanet,i),err
        if(err.eq.0.0d0) goto 10
        i=i+1
      goto 10
 11   continue
      ntt(nplanet)=i-1

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine lininterp(x,y,npmax,nmax,np,npt,xin,yout)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer i,npmax,nmax,npt(npmax),np
      double precision x(npmax,nmax),y(npmax,nmax),xin,yout

C     Default is zero
      yout=0.0d0
      if(npt(np).eq.0) return

c      write(0,*) (rhoin(i),i=1,npt)
c      write(0,*) (rhoierr(i),i=1,npt)
      if(xin.lt.x(np,1))then
        yout=y(np,1)+(xin-x(np,1))/(x(np,2)-x(np,1))*(y(np,2)-y(np,1))
      elseif(xin.gt.x(np,npt(np)))then
        yout=y(np,npt(np))+(xin-x(np,npt(np)))/(x(np,npt(np))-
     .      x(np,npt(np)-1))*(y(np,npt(np))-y(np,npt(np)-1))
      else
        do 10 i=1,npt(np)-1
            if((xin.gt.x(np,i)).and.(xin.le.x(np,i+1)))then
                yout=y(np,i)+(xin-x(np,i))/(x(np,i+1)-x(np,i))*
     .              (y(np,i+1)-y(np,i))
            endif
 10     continue
      endif
c      write(0,*) drho,dsig
c      read(5,*)

      return
      end

      subroutine occultquad(z0,u1,u2,p,muo1,mu0,nz)
C  This routine computes the lightcurve for occultation
C  of a quadratically limb-darkened source without microlensing.
C  Please cite Mandel & Agol (2002) if you make use of this routine
C  in your research.  Please report errors or bugs to agol@tapir.caltech.edu
      implicit none
      integer i,nz
      double precision z0(nz),u1,u2,p,muo1(nz),mu0(nz),
     &       mu(nz),lambdad(nz),etad(nz),lambdae(nz),lam,
     &       pi,x1,x2,x3,z,omega,kap0,kap1,q,Kk,Ek,Pk,n,ellec,ellk,rj
      if(abs(p-0.5d0).lt.1.d-3) p=0.5d0
C
C Input:
C
C rs   radius of the source (set to unity)
C z0   impact parameter in units of rs
C p    occulting star size in units of rs
C u1   linear    limb-darkening coefficient (gamma_1 in paper)
C u2   quadratic limb-darkening coefficient (gamma_2 in paper)
C
C Output:
C
C muo1 fraction of flux at each z0 for a limb-darkened source
C mu0  fraction of flux at each z0 for a uniform source
C
C Limb darkening has the form:
C  I(r)=[1-u1*(1-sqrt(1-(r/rs)^2))-u2*(1-sqrt(1-(r/rs)^2))^2]/(1-u1/3-u2/6)/pi
C
C To use this routine
C
C Now, compute pure occultation curve:
      omega=1.d0-u1/3.d0-u2/6.d0
      pi=acos(-1.d0)
C Loop over each impact parameter:
      do i=1,nz
C substitute z=z0(i) to shorten expressions
        z=z0(i)
        x1=(p-z)**2
        x2=(p+z)**2
        x3=p**2-z**2
C the source is unocculted:
C Table 3, I.
        if(z.ge.1.d0+p) then
          lambdad(i)=0.d0
          etad(i)=0.d0
          lambdae(i)=0.d0
          goto 10
        endif
C the  source is completely occulted:
C Table 3, II.
        if(p.ge.1.d0.and.z.le.p-1.d0) then
          lambdad(i)=1.d0
          etad(i)=1.d0
          lambdae(i)=1.d0
          goto 10
        endif
C the source is partly occulted and the occulting object crosses the limb:
C Equation (26):
        if(z.ge.abs(1.d0-p).and.z.le.1.d0+p) then
          kap1=acos(min((1.d0-p*p+z*z)/2.d0/z,1.d0))
          kap0=acos(min((p*p+z*z-1.d0)/2.d0/p/z,1.d0))
          lambdae(i)=p*p*kap0+kap1
          lambdae(i)=(lambdae(i)-0.5d0*sqrt(max(4.d0*z*z-
     &               (1.d0+z*z-p*p)**2,0.d0)))/pi
        endif
C the occulting object transits the source star (but doesn't
C completely cover it):
        if(z.le.1.d0-p) lambdae(i)=p*p
C the edge of the occulting star lies at the origin- special
C expressions in this case:
        if(abs(z-p).lt.1.d-4*(z+p)) then
C Table 3, Case V.:
          if(z.ge.0.5d0) then
            lam=0.5d0*pi
            q=0.5d0/p
            Kk=ellk(q)
            Ek=ellec(q)
C Equation 34: lambda_3
            lambdad(i)=1.d0/3.d0+16.d0*p/9.d0/pi*(2.d0*p*p-1.d0)*Ek-
     &                 (32.d0*p**4-20.d0*p*p+3.d0)/9.d0/pi/p*Kk
C Equation 34: eta_1
            etad(i)=1.d0/2.d0/pi*(kap1+p*p*(p*p+2.d0*z*z)*kap0-
     &              (1.d0+5.d0*p*p+z*z)/4.d0*sqrt((1.d0-x1)*(x2-1.d0)))
            if(p.eq.0.5d0) then
C Case VIII: p=1/2, z=1/2
              lambdad(i)=1.d0/3.d0-4.d0/pi/9.d0
              etad(i)=3.d0/32.d0
            endif
            goto 10
          else
C Table 3, Case VI.:
            lam=0.5d0*pi
            q=2.d0*p
            Kk=ellk(q)
            Ek=ellec(q)
C Equation 34: lambda_4
            lambdad(i)=1.d0/3.d0+2.d0/9.d0/pi*(4.d0*(2.d0*p*p-1.d0)*Ek+
     &                 (1.d0-4.d0*p*p)*Kk)
C Equation 34: eta_2
            etad(i)=p*p/2.d0*(p*p+2.d0*z*z)
            goto 10
          endif
        endif
C the occulting star partly occults the source and crosses the limb:
C Table 3, Case III:
        if((z.gt.0.5d0+abs(p-0.5d0).and.z.lt.1.d0+p).or.(p.gt.0.5d0.
     &      and.z.gt.abs(1.d0-p)*1.0001d0.and.z.lt.p)) then
          lam=0.5d0*pi
          q=sqrt((1.d0-(p-z)**2)/4.d0/z/p)
          if(q.gt.1.0d0) q=0.99999d0
          Kk=ellk(q)
          Ek=ellec(q)
          n=1.d0/x1-1.d0
          Pk=Kk-n/3.d0*rj(0.d0,1.d0-q*q,1.d0,1.d0+n)
C Equation 34, lambda_1:
          lambdad(i)=1.d0/9.d0/pi/sqrt(p*z)*(((1.d0-x2)*(2.d0*x2+
     &        x1-3.d0)-3.d0*x3*(x2-2.d0))*Kk+4.d0*p*z*(z*z+
     &        7.d0*p*p-4.d0)*Ek-3.d0*x3/x1*Pk)
          if(z.lt.p) lambdad(i)=lambdad(i)+2.d0/3.d0
C Equation 34, eta_1:
          etad(i)=1.d0/2.d0/pi*(kap1+p*p*(p*p+2.d0*z*z)*kap0-
     &          (1.d0+5.d0*p*p+z*z)/4.d0*sqrt((1.d0-x1)*(x2-1.d0)))
          goto 10
        endif
C the occulting star transits the source:
C Table 3, Case IV.:
        if(p.le.1.d0.and.z.le.(1.d0-p)*1.0001d0) then
          lam=0.5d0*pi
          q=sqrt((x2-x1)/(1.d0-x1))
          if(q.gt.1.0d0) q=0.99999d0
          Kk=ellk(q)
          Ek=ellec(q)
          n=x2/x1-1.d0
          Pk=Kk-n/3.d0*rj(0.d0,1.d0-q*q,1.d0,1.d0+n)
C Equation 34, lambda_2:
          lambdad(i)=2.d0/9.d0/pi/sqrt(1.d0-x1)*((1.d0-5.d0*z*z+p*p+
     &         x3*x3)*Kk+(1.d0-x1)*(z*z+7.d0*p*p-4.d0)*Ek-3.d0*x3/x1*Pk)
          if(z.lt.p) lambdad(i)=lambdad(i)+2.d0/3.d0
          if(abs(p+z-1.d0).le.1.d-4) then
            lambdad(i)=2/3.d0/pi*acos(1.d0-2.d0*p)-4.d0/9.d0/pi*
     &            sqrt(p*(1.d0-p))*(3.d0+2.d0*p-8.d0*p*p)
          endif
C Equation 34, eta_2:
          etad(i)=p*p/2.d0*(p*p+2.d0*z*z)
        endif
 10     continue
C Now, using equation (33):
        muo1(i)=1.d0-((1.d0-u1-2.d0*u2)*lambdae(i)+(u1+2.d0*u2)*
     &      lambdad(i)+u2*etad(i))/omega
C Equation 25:
        mu0(i)=1.d0-lambdae(i)
      enddo
      return
      end

      FUNCTION rc(x,y)
      REAL*8 rc,x,y,ERRTOL,TINY,SQRTNY,BIG,TNBG,COMP1,COMP2,THIRD,C1,C2,
     *C3,C4
      PARAMETER (ERRTOL=.04d0,TINY=1.69d-38,SQRTNY=1.3d-19,BIG=3.d37,
     *TNBG=TINY*BIG,COMP1=2.236d0/SQRTNY,COMP2=TNBG*TNBG/25.d0,
     *THIRD=1.d0/3.d0,C1=.3d0,C2=1.d0/7.d0,C3=.375d0,C4=9.d0/22.d0)
      REAL*8 alamb,ave,s,w,xt,yt
      if(x.lt.0..or.y.eq.0..or.(x+abs(y)).lt.TINY.or.(x+
     *abs(y)).gt.BIG.or.(y.lt.-COMP1.and.x.gt.0..and.x.lt.COMP2))pause
     *'invalid arguments in rc'
      if(y.gt.0.d0)then
        xt=x
        yt=y
        w=1.
      else
        xt=x-y
        yt=-y
        w=sqrt(x)/sqrt(xt)
      endif
1     continue
        alamb=2.d0*sqrt(xt)*sqrt(yt)+yt
        xt=.25d0*(xt+alamb)
        yt=.25d0*(yt+alamb)
        ave=THIRD*(xt+yt+yt)
        s=(yt-ave)/ave
      if(abs(s).gt.ERRTOL)goto 1
      rc=w*(1.d0+s*s*(C1+s*(C2+s*(C3+s*C4))))/sqrt(ave)
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software

      FUNCTION rj(x,y,z,p)
      REAL*8 rj,p,x,y,z,ERRTOL,TINY,BIG,C1,C2,C3,C4,C5,C6,C7,C8
      PARAMETER (ERRTOL=.05d0,TINY=2.5d-13,BIG=9.d11,C1=3.d0/14.d0,
     *C2=1.d0/3.d0,C3=3.d0/22.d0,C4=3.d0/26.d0,C5=.75d0*C3,
     *C6=1.5d0*C4,C7=.5d0*C2,C8=C3+C3)
CU    USES rc,rf
      REAL*8 a,alamb,alpha,ave,b,beta,delp,delx,dely,delz,ea,eb,ec,ed,
     *ee,fac,pt,rcx,rho,sqrtx,sqrty,sqrtz,sum,tau,xt,yt,zt,rc,rf
      if(min(x,y,z).lt.0..or.min(x+y,x+z,y+z,abs(p)).lt.TINY.or.max(x,y,
     *z,abs(p)).gt.BIG)then
         write(0,*) x,y,z,abs(p)
         pause 'invalid arguments in rj'
      endif
      sum=0.d0
      fac=1.d0
      if(p.gt.0.d0)then
        xt=x
        yt=y
        zt=z
        pt=p
      else
        xt=min(x,y,z)
        zt=max(x,y,z)
        yt=x+y+z-xt-zt
        a=1.d0/(yt-p)
        b=a*(zt-yt)*(yt-xt)
        pt=yt+b
        rho=xt*zt/yt
        tau=p*pt/yt
        rcx=rc(rho,tau)
      endif
1     continue
        sqrtx=sqrt(xt)
        sqrty=sqrt(yt)
        sqrtz=sqrt(zt)
        alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
        alpha=(pt*(sqrtx+sqrty+sqrtz)+sqrtx*sqrty*sqrtz)**2
        beta=pt*(pt+alamb)**2
        sum=sum+fac*rc(alpha,beta)
        fac=.25d0*fac
        xt=.25d0*(xt+alamb)
        yt=.25d0*(yt+alamb)
        zt=.25d0*(zt+alamb)
        pt=.25d0*(pt+alamb)
        ave=.2d0*(xt+yt+zt+pt+pt)
        delx=(ave-xt)/ave
        dely=(ave-yt)/ave
        delz=(ave-zt)/ave
        delp=(ave-pt)/ave
      if(max(abs(delx),abs(dely),abs(delz),abs(delp)).gt.ERRTOL)goto 1
      ea=delx*(dely+delz)+dely*delz
      eb=delx*dely*delz
      ec=delp**2
      ed=ea-3.d0*ec
      ee=eb+2.d0*delp*(ea-ec)
      rj=3.d0*sum+fac*(1.d0+ed*(-C1+C5*ed-C6*ee)+eb*(C7+delp*
     *(-C8+delp*C4))+delp*ea*(C2-delp*C3)-C2*delp*ec)/(ave*sqrt(ave))
      if (p.le.0.d0) rj=a*(b*rj+3.d0*(rcx-rf(xt,yt,zt)))
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software

      function ellec(k)
      implicit none
      double precision k,m1,a1,a2,a3,a4,b1,b2,b3,b4,ee1,ee2,ellec
C Computes polynomial approximation for the complete elliptic
C integral of the second kind (Hasting's approximation):
      m1=1.d0-k*k
      a1=0.44325141463d0
      a2=0.06260601220d0
      a3=0.04757383546d0
      a4=0.01736506451d0
      b1=0.24998368310d0
      b2=0.09200180037d0
      b3=0.04069697526d0
      b4=0.00526449639d0
      ee1=1.d0+m1*(a1+m1*(a2+m1*(a3+m1*a4)))
      ee2=m1*(b1+m1*(b2+m1*(b3+m1*b4)))*log(1.d0/m1)
      ellec=ee1+ee2
      return
      end

      function ellk(k)
      implicit none
      double precision a0,a1,a2,a3,a4,b0,b1,b2,b3,b4,ellk,
     &       ek1,ek2,k,m1
C Computes polynomial approximation for the complete elliptic
C integral of the first kind (Hasting's approximation):
      m1=1.d0-k*k
      a0=1.38629436112d0
      a1=0.09666344259d0
      a2=0.03590092383d0
      a3=0.03742563713d0
      a4=0.01451196212d0
      b0=0.5d0
      b1=0.12498593597d0
      b2=0.06880248576d0
      b3=0.03328355346d0
      b4=0.00441787012d0
      ek1=a0+m1*(a1+m1*(a2+m1*(a3+m1*a4)))
      ek2=(b0+m1*(b1+m1*(b2+m1*(b3+m1*b4))))*log(m1)
      ellk=ek1-ek2
      return
      end

      FUNCTION rf(x,y,z)
      REAL*8 rf,x,y,z,ERRTOL,TINY,BIG,THIRD,C1,C2,C3,C4
      PARAMETER (ERRTOL=.08d0,TINY=1.5d-38,BIG=3.d37,THIRD=1.d0/3.d0,
     *C1=1.d0/24.d0,C2=.1d0,C3=3.d0/44.d0,C4=1.d0/14.d0)
      REAL*8 alamb,ave,delx,dely,delz,e2,e3,sqrtx,sqrty,sqrtz,xt,yt,zt
      if(min(x,y,z).lt.0.d0.or.min(x+y,x+z,y+z).lt.TINY.or.max(x,y,
     *z).gt.BIG)pause 'invalid arguments in rf'
      xt=x
      yt=y
      zt=z
1     continue
        sqrtx=sqrt(xt)
        sqrty=sqrt(yt)
        sqrtz=sqrt(zt)
        alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
        xt=.25d0*(xt+alamb)
        yt=.25d0*(yt+alamb)
        zt=.25d0*(zt+alamb)
        ave=THIRD*(xt+yt+zt)
        delx=(ave-xt)/ave
        dely=(ave-yt)/ave
        delz=(ave-zt)/ave
      if(max(abs(delx),abs(dely),abs(delz)).gt.ERRTOL)goto 1
      e2=delx*dely-delz**2
      e3=delx*dely*delz
      rf=(1.d0+(C1*e2-C2-C3*e3)*e2+C4*e3)/sqrt(ave)
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software 0NL&WR2.

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      double precision function mandelagol(nintg,R1,R2,x1,x2,y1,y2,c,
     .  b0,mu,mulimb0,mulimbf,dist)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Adapted from Mandel and Agol, 2002, ApJ 580, L171
      implicit none
      integer i,nintg,sflag
      double precision R1,R2,x1,x2(nintg),y1,y2(nintg),c(4),
     .  c1,c2,c3,c4,mulimb0(nintg),mulimbf(nintg,5),rl,b0(nintg),
     .  mu(nintg),dist(nintg)

      mu=0

      c1=c(1)
      c2=c(2)
      c3=c(3)
      c4=c(4)
      rl=R2/R1
      sflag=0
      do 10 i=1,nintg
       dist(i)=Sqrt((x2(i)-x1)*(x2(i)-x1)+(y2(i)-y1)*(y2(i)-y1))/(R1+R2)
c        if(dist(i).ge.1.0d0)then
c            sflag=sflag+1
c            b0(i)=2.0
c        else
            b0(i)=(R1+R2)*dist(i)/R1
c        endif
 10   continue

c      write(6,500) "hello",(b0(i),i=1,nintg)
 500  format(A5,11(1X,F5.3))
      if(sflag.lt.nintg) then
c          call occultnl(rl,c1,c2,c3,c4,b0,mulimb0,mulimbf,nintg)
c          call occultsmall(rl,c1,c2,c3,c4,nintg,b0,mulimb0)
      endif
c      if(b0(1).le.1.0) write(6,*)b0(1),mu

C      mandelagol=mulimb0(1)
      mandelagol=0.0
      do 11 i=1,nintg
c        mandelagol=mandelagol+mu(i)
c        if(dist(i).ge.1.0d0) mulimb0(i)=1.0d0
        mandelagol=mandelagol+mulimb0(i)
 11   continue
      mandelagol=mandelagol/dble(nintg)
c      write(6,*) "hello2",mandelagol

      return
      end

      subroutine occultsmall(p,c1,c2,c3,c4,nz,z,mu)
      implicit none
      integer i,nz
c      parameter (nz=201)
      real*8 p,c1,c2,c3,c4,z(nz),mu(nz),i1,norm,
     &       x,tmp,iofr,pi
C This routine approximates the lightcurve for a small
C planet. (See section 5 of Mandel & Agol (2002) for
C details):
C Input:
C  p      ratio of planet radius to stellar radius
C  c1-c4  non-linear limb-darkening coefficients
C  z      impact parameters (positive number normalized to stellar
C        radius)- this is an array which MUST be input to the routine
C  NOTE:  nz must match the size of z & mu in calling routine
C Output:
C  mu     flux relative to unobscured source for each z
C
      pi=acos(-1.d0)
      norm=pi*(1.d0-c1/5.d0-c2/3.d0-3.d0*c3/7.d0-c4/2.d0)
      i1=1.d0-c1-c2-c3-c4
      do i=1,nz
        mu(i)=1.d0
        if(z(i).gt.1.d0-p.and.z(i).lt.1.d0+p) then
          x=1.d0-(z(i)-p)**2
          tmp=(1.d0-c1*(1.d0-0.8d0*x**0.25)
     &             -c2*(1.d0-2.d0/3.d0*x**0.5)
     &             -c3*(1.d0-4.d0/7.d0*x**0.75)
     &             -c4*(1.d0-0.5d0*x))
          mu(i)=1.d0-tmp*(p**2*acos((z(i)-1.d0)/p)
     &        -(z(i)-1.d0)*sqrt(p**2-(z(i)-1.d0)**2))/norm
        endif
        if(z(i).le.1.d0-p.and.z(i).ne.0.d0) then
          mu(i)=1.d0-pi*p**2*iofr(c1,c2,c3,c4,z(i),p)/norm
        endif
        if(z(i).eq.0.d0) then
          mu(i)=1.d0-pi*p**2/norm
        endif
      enddo
      return
      end

      function iofr(c1,c2,c3,c4,r,p)
      implicit none
      real*8 r,p,c1,c2,c3,c4,sig1,sig2,iofr
      sig1=sqrt(sqrt(1.d0-(r-p)**2))
      sig2=sqrt(sqrt(1.d0-(r+p)**2))
      iofr=1.d0-c1*(1.d0+(sig2**5-sig1**5)/5.d0/p/r)
     &         -c2*(1.d0+(sig2**6-sig1**6)/6.d0/p/r)
     &         -c3*(1.d0+(sig2**7-sig1**7)/7.d0/p/r)
     &         -c4*(p**2+r**2)
      return
      end

      subroutine occultnl(rl,c1,c2,c3,c4,b0,mulimb0,mulimbf,nb)
c; Please cite Mandel & Agol (2002) if making use of this routine.
      implicit none
      integer i,j,nb,nr,i1,i2,nmax
      parameter (nmax=2**16)
      real*8 mulimbf(nb,5),pi,c1,c2,c3,c4,rl,bt0(nb),b0(nb),mulimb0(nb),
     &       mulimb(nb),mulimbp(nb),dt,t(nmax+1),th(nmax+1),r(nmax+1),
     &       sig,mulimb1(nb),mulimbhalf(nb),mulimb3half(nb),mulimb2(nb),
     &       sig1,sig2,omega,dmumax,fac,mu(nb),f1,f2
      pi=acos(-1.d0)
C  This routine uses the results for a uniform source to
C  compute the lightcurve for a limb-darkened source
C  (5-1-02 notes)
C Input:
C   rl        radius of the lens   in units of the source radius
C   c1-c4     limb-darkening coefficients
C   b0        impact parameter normalized to source radius
C Output:
C  mulimb0 limb-darkened magnification
C  mulimbf lightcurves for each component
C
C  First, make grid in radius:
C  Call magnification of uniform source:
      call occultuniform(b0,rl,mulimb0,nb)
      i1=nb
      i2=1
      fac=0.d0
      do i=1,nb
        bt0(i)=b0(i)
        mulimbf(i,1)=1.d0
        mulimbf(i,2)=0.8d0
        mulimbf(i,3)=2.d0/3.d0
        mulimbf(i,4)=4.d0/7.d0
        mulimbf(i,5)=0.5d0
        mulimb(i)=mulimb0(i)
        if(mulimb0(i).ne.1.d0) then
          i1=min(i1,i)
          i2=max(i2,i)
        endif
        fac=max(fac,abs(mulimb0(i)-1.d0))
      enddo
C print,rl
      omega=4.*((1.d0-c1-c2-c3-c4)/4.+c1/5.+c2/6.+c3/7.+c4/8.)
      nr=2
      dmumax=1.d0
c      write(6,*) 'i1,i2 ',i1,i2
      do while (dmumax.gt.fac*1.d-3)
        do i=i1,i2
          mulimbp(i)=mulimb(i)
        enddo
        nr=nr*2
c        write(6,*) 'nr ',nr
        dt=0.5d0*pi/dble(nr)
        if(nr+1.gt.nmax) then
         write(0,*) "M&A Seg: ",nr+1,b0(1)
         return
        endif
        do j=1,nr+1
          t(j) =dt*dble(j-1)
          th(j)=t(j)+0.5d0*dt
          r(j)=sin(t(j))
        enddo
        sig=sqrt(cos(th(nr)))
        do i=i1,i2
          mulimbhalf(i) =sig**3*mulimb0(i)/(1.d0-r(nr))
          mulimb1(i)    =sig**4*mulimb0(i)/(1.d0-r(nr))
          mulimb3half(i)=sig**5*mulimb0(i)/(1.d0-r(nr))
          mulimb2(i)    =sig**6*mulimb0(i)/(1.d0-r(nr))
        enddo
        do j=2,nr
          do i=1,nb
            b0(i)=bt0(i)/r(j)
          enddo
C  Calculate uniform magnification at intermediate radii:
          call occultuniform(b0,rl/r(j),mu,nb)
C  Equation (29):
          sig1=sqrt(cos(th(j-1)))
          sig2=sqrt(cos(th(j)))
          dmumax=0.d0
          do i=i1,i2
            f1=r(j)*r(j)*mu(i)/(r(j)-r(j-1))
            f2=r(j)*r(j)*mu(i)/(r(j+1)-r(j))
            mulimbhalf(i) =mulimbhalf(i) +f1*sig1**3-f2*sig2**3
            mulimb1(i)    =mulimb1(i)    +f1*sig1**4-f2*sig2**4
            mulimb3half(i)=mulimb3half(i)+f1*sig1**5-f2*sig2**5
            mulimb2(i)    =mulimb2(i)    +f1*sig1**6-f2*sig2**6
            mulimb(i)=((1.d0-c1-c2-c3-c4)*mulimb0(i)+c1*mulimbhalf(i)*dt
     &        +c2*mulimb1(i)*dt+c3*mulimb3half(i)*dt+c4*mulimb2(i)*dt)
     &        /omega
            if(mulimb(i)+mulimbp(i).ne.0.d0) then
              dmumax=max(dmumax,abs(mulimb(i)-mulimbp(i))/(mulimb(i)+
     &               mulimbp(i)))
            endif
          enddo
        enddo
      enddo
      do i=i1,i2
        mulimbf(i,1)=mulimb0(i)
        mulimbf(i,2)=mulimbhalf(i)*dt
        mulimbf(i,3)=mulimb1(i)*dt
        mulimbf(i,4)=mulimb3half(i)*dt
        mulimbf(i,5)=mulimb2(i)*dt
        mulimb0(i)=mulimb(i)
      enddo
      do i=1,nb
        b0(i)=bt0(i)
      enddo
      return
      end

      subroutine occultuniform(b0,w,muo1,nb)
      implicit none
      integer i,nb
      real*8 muo1(nb),w,b0(nb),z,pi,lambdae,kap0,kap1
      if(abs(w-0.5d0).lt.1.d-3) w=0.5d0
      pi=acos(-1.d0)
C  This routine computes the lightcurve for occultation
C  of a uniform source without microlensing  (Mandel & Agol 2002).
C Input:
C
C  rs   radius of the source (set to unity)
C  b0   impact parameter in units of rs
C  w    occulting star size in units of rs
C
C Output:
C  muo1 fraction of flux at each b0 for a uniform source
C
C  Now, compute pure occultation curve:
      do i=1,nb
C  substitute z=b0(i) to shorten expressions
        z=b0(i)
C  the source is unocculted:
C  Table 3, I.
        if(z.ge.1.d0+w) then
          muo1(i)=1.d0
          goto 1
        endif
C  the  source is completely occulted:
C  Table 3, II.
        if(w.ge.1.d0.and.z.le.w-1.d0) then
          muo1(i)=0.d0
          goto 1
        endif
C  the source is partly occulted and the occulting object crosses the limb:
C  Equation (26):
        if(z.ge.abs(1.d0-w).and.z.le.1.d0+w) then
          kap1=acos(min((1.d0-w*w+z*z)/2.d0/z,1.d0))
          kap0=acos(min((w*w+z*z-1.d0)/2.d0/w/z,1.d0))
          lambdae=w*w*kap0+kap1
          lambdae=(lambdae-0.5d0*sqrt(max(4.d0*z*z-(1.d0+z*z-w*w)**2,
     &            0.d0)))/pi
          muo1(i)=1.d0-lambdae
        endif
C  the occulting object transits the source star (but doesn't
C  completely cover it):
        if(z.le.1.d0-w) muo1(i)=1.d0-w*w
 1      continue
      enddo
C muo1=1.d0-lambdae
      return
      end

