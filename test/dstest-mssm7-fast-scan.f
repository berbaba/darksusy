      program dstest
c
c     This program tests some DarkSUSY routines
c
c     In addition to the DarkSUSY programs and data files, this test
c     program chooses a few mSUGRA models randomly and uses Isasugra to
c     calculate the low energy models. These are then fed into DarkSUSY
c     and the relic density is calculated.
c     If you want to calculate rates and stuff, use the calls as
c     outlined in dstest.f

c     
c-----This line is 72 columns long--------------------------------------
c
      implicit none
      real*8 oh2,xf,dsrdomega                            ! relic density
      character*80 message,scr
      integer ii,jj,idum,nfc,iwar,ierr,iend,unphys,hwarning,excl,nevts,
     &  ipass 
      real*8 sigsip,sigsin,sigsdp,sigsdn           ! cross-sections
      real*8 eth,thmax,rateea,energy,theta  	   ! neutrino telescopes
      real*8 gaugfr,ratesunu,ratesumu,convrate 	   ! neutrino telescopes
      integer m,ch,rtype,ptype,istat               ! neutrino telescopes
      real*8 Emin,Emax                             ! neutrino telescopes
      integer nbins, masswimp			   ! neutrino telescopes
      real*8 tabE,tabYieldneutrino                 ! neutrino telescopes
      character channel*20			   ! neutrino telescopes
      character outputfile*40	 	 	   ! neutrino telescopes     
      real*8 dsrndlin,dsrndsgn,dsabsq
c
c     Here we include the file dssusy.h which contains common block
c     definitions for particle masses, susy parameters, and results.
c     We also include dsidtag.h which gives the possibility to tag
c     every model with a 12 character id tag (variable idtag). This
c     id tag is printed in case of errors for that model.
c
      include 'dssusy.h'
      include 'dsidtag.h'

c      idum = -363469   ! seed for random number generator
      
      write(*,*) 'How many model points to scan?'
      read(5,*) nevts
      write(*,*) 'Which SEED?'
      read(5,*) idum
      if (nevts.eq.0.or.idum.eq.0) goto 100  ! exit program
c
c     This call initializes the DarkSUSY package. This call should
c     be the first call in any program using DarkSUSY.
c
      call dsinit
c      call aldata   ! to make sure isajet common block datas are set
                    ! on most compilers, this is not needed

c
c     The amount of output onto standard output can be controlled by the
c     variable prtlevel. Setting prtlevel=0 suppresses all messages,
c     except fatal errors that terminate the program. Setting prtlevel=1
c     displays informational messages. Setting prtlevel to higher values
c     displays more and more information. In this example, we set
c     prtlevel=1 to have a message from dsreadpar. We reset it to 0
c     later. Had we set prtlevel=2 here, we would also have messages for
c     each variable which is (re)set in dsreadpar. Caution: if prtlevel
c     is reset in dsreadpar, the new value applies to subsequent
c     messages.
c
      prtlevel=0
c
c     This is an example of redefining DarkSUSY parameters reading them
c     from an input file, in this case the file dsdefault.dat which
c     contains the default values (this is done for test purposes only,
c     because the default values are already included in the code). The
c     user can copy the file dsdefault.dat to a new file, change the
c     values of the variables, and read in the new file.  Not all
c     variables need to be assigned in the new file, but only those
c     whose value differs from the default value. Note: dsdefault.dat
c     resets prtlevel to 0, ie no messages.
c     Note: the dsdefault.dat is not upto date, but a call like the
c     following can be used to reset parameters to your own choices.
c
c      open (unit=10,file='dsdefault.dat',status='old')
c      call dsreadpar(10)
c      close (10)


c     Now we are ready to scan the supersymmetric parameter space. 
c     We here read some models from dstest-isasugra.mod, which
c     are a few of the Benchmark points in hep-ph/0306219
c
      

c      open(unit=10,file='dstest-isasugra.mod',status='old',
c     &  form='formatted')
c      read(10,'(a)') scr ! read header lines
c      read(10,'(a)') scr ! read header lines
c      read(10,'(a)') scr ! read header lines


c     these are minimal SUGRA
c 20   read(10,1001,end=100) idtag,m0,mhf,a0,sgnmu,tgb
c 1001 format(1x,A12,3(1x,E14.8),1x,e7.1,5(1x,e14.8))
      
      write(outputfile,'(a25)')'output_Sun_mssm7_scan.txt'
      open (unit=41,file=outputfile)
      
      iend=0
      ipass=0
      
 20   call random_model(nevts,idum,ipass,iend)

      if (iend.eq.1) then
        goto 100  ! exit program
      endif	
      
c      write(*,*) ' '
c      write(*,*) 'Model: ',idtag
c      write(*,*) 'Parameters: ',m0,mhf,a0,sgnmu,tgb

c      Now define the model, see dsgive_model_isasugra how we transfer
c      these varibles to DarkSUSY common blocks
c      call dsgive_model_isasugra(m0,mhf,a0,sgnmu,tgb)
         

c      Now set up the model by a call to dssusy_isasugra. dssusy_isasugra
c      calls isasugra and transfers the results to DarkSUSY common blocks.
c      DarkSUSY is then set-up for this model. Check the flags unphys
c      and hwarning which are non-zero if the model is not good.
      call dssusy(unphys,hwarning)

c      write(*,*) '   hwarning = ',hwarning
c      write(*,*) '     unphys = ',unphys
      if (hwarning.eq.0.and.unphys.eq.0) then  ! Model OK
      
c      call dsacset('pdg2002f')      
      call dsacbnd(excl)
      
      if (excl.eq.0) then ! Accelerators exclusions OK

c...Now calculate the relic density with all coannihilations
c       Calculating omega h^2, fast method...'
        oh2=dsrdomega(1,1,xf,ierr,iwar,nfc)

c     WMAP 7yrs limit: D. Larson, J. Dunkley, et. al.
c     ApJS 192 (2011) 14
c     Omega(WMAP-7yrs)h² = 0.1120 +- 0.0056 (1sigma) 
c			 = 0.1120 +- 0.0112 (2sigma)

      if(oh2.le.0.1232d0)then

      ipass=1
      
c
c     We now have all the masses and couplings calculated. Let's print
c     some out to see what they are. dsabsq is just the square of the
c     absolute value of it's complex*16 argument.

c        write(*,*) '  Neutralino mass = ',mass(kn(1))
        gaugfr=dsabsq(neunmx(1,1))+dsabsq(neunmx(1,2))
c        write(*,*) '  Gaugino fraction = ',
c     &    gaugfr
c        write(*,*) '  H1 mass =  ',mass(kh1),width(kh1)
c        write(*,*) '  H2 mass =  ',mass(kh2),width(kh2)
c        write(*,*) '  H3 mass =  ',mass(kh3),width(kh3)
c        write(*,*) '  H+- mass = ',mass(khc),width(khc)

c...Now calculate the relic density with all coannihilations
c        write(*,*) '   Calculating omega h^2, please be patient...'
c        oh2=dsrdomega(1,1,xf,ierr,iwar,nfc)
c        write(*,*) '   Oh2 = ',oh2

c...Go on calculating rates etc, as shown in dstest.f if we want to
      
c     Neutrino telescope rates
c
c     To calculate the rates in neutrino telescopes we call dsntrates. We
c     can either calculate the flux of neutrinos, the conversion rate
c     (per volume element) of neutrinos to muons or the muon flux. This
c     is determined by the argument rtype. We can also choose the
c     energy threshold and the maximal half-aperture angle from the
c     center of the Earth/Sun we are interested in. The rates for both
c     the earth and the sun are returned in units of km^-2 yr^-1 for
c     rtype=1,3 and in units of km^-3 yr^-1 for rtype=2. If some
c     warnings were issued, the flag istat is non-zero.
c
c     The defualt calculation method is to use the full expressions by
c     Gould and numerically integrate them over the velocity distribution
c     as specified by the halo profile, e.g. a gaussian for an
c     isothermal sphere. This numerical integration is rather slow and
c     there is thus an option to use tables (read from disk, recreated
c     if absent) instead. This is the default. Changing to numerical 
c     integration directly is handled by a call to dsntset (see dsntset.f
c     for details). Also other (approximate) formulae or the Damour
c     Krauss population are available and can be chosen by a call to dsntset.
c     Note: For the earth, which captures from a population of WIMPs bound
c     in the solar system, a new estimate (Lundberg and Edsjo, 
c     astro-ph/0401113) is used as a default.
c
c     

c      write(*,*) 'Calculating rates in neutrino telescopes'
      eth=1.0d0      ! energy threshold (of neutrino/muon), GeV
      thmax=10.0d0   ! the maximum half-aperture angle, degrees
      rtype=1        ! 1=neutrino flux
                     ! 2=neutrino to muon conversion rate
                     ! 3=muon flux
      ptype=3        ! 1=particles only
                     ! 2=anti-particles only
                     ! 3=sum of particle and anti-particle rates

      call dsntrates(eth,thmax,rtype,ptype,rateea,ratesunu,istat)
      rtype=3
      call dsntrates(eth,thmax,rtype,ptype,rateea,ratesumu,istat)
      rtype=2
      call dsntrates(eth,thmax,rtype,ptype,rateea,convrate,istat)

c      write(*,*) '  Flux from the Earth = ',rateea,
c     &  ' km^-2 yr^-1'
c      write(*,*) '  Flux from the Sun =   ',ratesu,
c     &  ' km^-2 yr^-1'            

c
c     Let's start with scattering
c     cross sections for direct detection experiments by calling
c     dsddneunuc which returns the spin-independent and the spin-dependent
c     scattering cross sections off protons and neutrons. The cross
c     sections are returned in units of cm^2.
c
c     The next line sets the spin independent form factors to be the best 
c     available, i.e. in order Fourier-Bessel, Sum-of-gaussians, Fermi, 
c     Lewin-Smith. This is also the default setting.
c     Use call dsrdset('help',' ') and dsrdset('si','help') 
c     to print out the other possibilities.
c
      call dsddset('si','best')
c
c     The next line sets the spin dependent form factors to be the best 
c     available, i.e. in order interacting shell model, odd group model, 
c     single particle shell model. This is also the default setting.
c     Use call dsrdset('help',' ') and dsrdset('sd','help') 
c     to print out the other possibilities.
c
      call dsddset('sd','best')

c      write (*,*) 'Calculating scattering cross sections...'
      call dsddneunuc(sigsip,sigsin,sigsdp,sigsdn)
c      write(*,*) '  sigsip (pb) = ',sigsip*1.0d36
c      write(*,*) '  sigsin (pb) = ',sigsin*1.0d36
c      write(*,*) '  sigsdp (pb) = ',sigsdp*1.0d36
c      write(*,*) '  sigsdn (pb) = ',sigsdn*1.0d36
c      call dsddsigma(1,1,1,sigsip,sigsdp)
c      write(*,*) ' proton: sigsi (pb) = ',sigsip*1.0d36,
c     &     ' sigsd (pb) = ',sigsdp*1.0d36

c     mass2q(3) == amsq**2
      write(41,51) mu,' ',m2,' ',ma,' ',tanbe,' ',mass2q(3),
     &   ' ',asoftu(3),' ',asoftd(3), 
     &   ' ',oh2,' ',mass(kn(1)),' ',gaugfr,' ',mass(kgluin),
     &   ' ',mass(kh1),' ',width(kh1),' ',mass(kh2),' ',width(kh2),
     &   ' ',mass(kh3),' ',width(kh3),' ',mass(khc),' ',width(khc),
     &   ' ',ratesunu,' ',ratesumu,' ',convrate,
     &   ' ',sigsdp*1.0d36,' ',sigsip*1.0d36,
     &   ' ',sigsdn*1.0d36,' ',sigsin*1.0d36
 51   format(f,a,f,a,f,a,f,a,f,a,f,a,f,a,e,a,f,a,e,a,f,a,f,a,e,a,f,a,
     &   e,a,f,a,e,a,f,a,e,a,e,a,e,a,e,a,e,a,e,a,e,a,e)
      
      else
	ipass=0
c     WMAP 7yrs - 2sigma exclusion
      endif
      
      else
        ipass=0
c     Accelerators exclusions
      endif
      
      else
	ipass=0
c     Models exclusions
      endif
c
c     end of loop
c
      goto 20

 100  continue

c      close(10)
      close(41)
      stop
      end
      
      
      
      subroutine random_model(nevts,idum,ipass,iend)
c
c     To generate model parameters in a random way
c
      implicit none
      include 'dssusy.h'
      include 'dsidtag.h'
      real*8 dsrndlog,dsrndlin,dsrndsgn
      integer first,n,i,iend,idum,nevts,ipass,nint,iint
      real*8 amumin,amumax,am2min,am2max,amamin,amamax,
     &  atbmin,atbmax,amsqmin,amsqmax,atmmin,atmmax,
     &  abmmin,abmmax,
     &  amu,am2,ama,atb,amsq,atm,abm
      real*8 amutmp,am2tmp,amatmp,atbtmp,amsqtmp,atmtmp,
     &  abmtmp
      data first/0/, n/0/,amutmp/0/,am2tmp/0/,amatmp/0/,
     &  atbtmp/0/,amsqtmp/0/,atmtmp/0/,abmtmp/0/,
     &  nint/0/,iint/0/
      save first,n,amutmp,am2tmp,amatmp,atbtmp,amsqtmp,
     &  atmtmp,abmtmp,nint,iint
c... First time set random number seed
      if (first.eq.0) then
         n=0
         first=1
	 amutmp=0
	 am2tmp=0
	 amatmp=0
	 atbtmp=0
	 amsqtmp=0
	 atmtmp=0
	 abmtmp=0
      endif
c...  Limits of the Mu-M2-MA-tan(beta)-Msq-At-Ab region explored
      amumin = 10.d0
      amumax = 10000.d0
      am2min = 10.d0
      am2max = 10000.d0
      amamin = 60.d0
      amamax = 1000.d0
      atbmin = 1.001d0
      atbmax = 60.d0
      amsqmin = 50.d0 
      amsqmax = 3000.d0
      atmmin = -3.d0
      atmmax = 3.d0
      abmmin = -3.d0
      abmmax = 3.d0
c... random point
      if (ipass.eq.1) then
         n=n+1
	 write (*,'(a,i8.8)') 'Model: rndm',n  
	 if(nint.eq.0)then
	   iint=1
	 endif
      endif

c      write (idtag,'(a,i8.8)') 'rndm',n  
      if (n.lt.nevts) then	
        if (ipass.eq.0.or.iint.eq.0) then
          amu = dsrndsgn(idum)*dsrndlog(idum,amumin,amumax)
	  amutmp=amu
          am2 = dsrndsgn(idum)*dsrndlog(idum,am2min,am2max)
	  am2tmp=am2
	  ama = dsrndlog(idum,amamin,amamax)	
	  amatmp=ama	   		 
          atb = dsrndlog(idum,atbmin,atbmax)
	  atbtmp=atb
          amsq = dsrndlog(idum,amsqmin,amsqmax) 
	  amsqtmp=amsq	   
          atm = dsrndlin(idum,atmmin,atmmax) 
	  atmtmp=atm
	  abm = dsrndlin(idum,abmmin,abmmax) 
	  abmtmp=abm
	else
	  amu=dsrndlin(idum,amutmp,amutmp+100.d0)
	  am2=dsrndlin(idum,am2tmp,am2tmp+100.d0)
	  ama=amatmp
	  atb=atbtmp
	  amsq=amsqtmp
	  atm=atmtmp
	  abm=abmtmp
	  nint=nint+1
	  if(nint.gt.100) then
	    iint=0
	    nint=0
	    amutmp=amutmp+100.d0
	    am2tmp=am2tmp+100.d0
	  endif  
	endif
						 	 
c       Now define the model, see dsgive_model_isasugra how we transfer 
c       these varibles to DarkSUSY common blocks		 	 
        call dsgive_model(amu,am2,ama,atb,amsq,atm,abm) 
c	write(*,*) ' '
c       write(*,*) 'Model: ',idtag
c       write(*,*) 'Parameters: ',amu,am2,ama,atb,amsq,atm,abm 	 
      else 
        iend = 1
      endif      
      return
      end
