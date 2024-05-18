module heatm
   ! provides heatr for NJOY2016
   use locale
   implicit none
   private
   public heatr

   ! global variables

   ! units for heatr
   integer::nendf,nin,nout,nplot

   ! heatr global parameters
   integer,parameter::nbuf=1000
   integer,parameter::nqamax=30
   integer,parameter::ilmax=200
   integer,parameter::maxmf6=320

   ! global variables
   integer::matd,mgam,npk,mtp(28),nqa,mta(nqamax),mt303,mt19
   integer::ne,kchk,iprint,lqs(nqamax)
   real(kr),dimension(:),allocatable::qbar
   real(kr)::qa(nqamax),efirst,elast,za,awr,elist(ilmax)
   real(kr)::break, break_i
   real(kr)::qdel,etabmax
   integer::mt103,mt104,mt105,mt106,mt107,mt16
   integer::miss4(250),nmiss4
   integer::i6,mt6(maxmf6),i6g,mt6no(maxmf6)
   integer::i6p,mt6yp(maxmf6)
   integer::jp,jpn,jpp
   real(kr)::q,zat,awrt,zap,awp
   integer::lct
   integer::idame
   real(kr)::ebot,etop
   integer::mt458,nply,lfc,ifc1,ifc2,ifc3,ifc4,ifc5,ifc6
   real(kr),dimension(:),allocatable::c458,cpoly,hpoly,afr,anp,agp
   real(kr)::emc2,tm,rtm
   integer::icntrl(12)

   ! target
   integer::izat

   ! projectile
   integer::izap

contains

   subroutine heatr
   !-------------------------------------------------------------------
   !
   ! Compute heating kerma (kinetic energy release in material)
   ! and radiation damage energy production.
   !
   ! The prompt kerma is computed pointwise on the grid of the
   ! total cross section from the input pendf tape and written
   ! onto the output PENDF tape at infinite dilution using the
   ! 300 series of MT numbers.  All temperatures on the input PENDF
   ! tape for the desired material are processed.  The dictionary
   ! is revised.  Reaction Q values are obtained from the ENDF
   ! tape unless the user enters his own value.  Partial kermas
   ! can be requested for self-shielding calculations or other
   ! purposes.  The code uses the energy balance method where
   ! photon files are available and deposits all photon energy
   ! locally when files are not available.  This assures
   ! consistency between neutron heating and energy deposition by
   ! subsequent photon interactions.  An exception is made for
   ! capture where recoil is computed by momentum conservation.
   ! Photon files are used to estimate the average photon momentum
   ! when available.  A diagnostic message is printed if the
   ! momentum calculation leads to a significant error in
   ! energy conservation.
   !
   ! If desired, the energy-balance kerma factors can be compared
   ! with conservative kinematic limits (set iprint=2).
   ! A plot file for viewr can be automatically prepared.
   !
   ! Damage energy is computed using the Lindhard electronic
   ! screening damage function with a displacement threshold
   ! from a table of default values for important elements
   ! or a value provided by the user.
   !
   !---input specifications (free format)---------------------------
   !
   ! card 1
   !    nendf    unit for endf tape
   !    nin      unit for input pendf tape
   !    nout     unit for output pendf tape
   !    nplot    unit for graphical check output
   ! card 2
   !    matd     material to be processed
   !    npk      number of partial kermas desired (default=0)
   !    nqa      number of user q values (default=0)
   !    ntemp    number of temperatures to process
   !             (default=0, meaning all on pendf)
   !    local    0/1=gamma rays transported/deposited locally
   !             (default=0)
   !    iprint   print (0 min, 1 max, 2 check) (default=0)
   !    ed       displacement energy for damage
   !             (default from built-in table)
   ! card 3      for npk gt 0 only
   !    mtk      mt numbers for partial kermas desired
   !             total (mt301) will be provided automatically.
   !             partial kerma for reaction mt is mt+300
   !             and may not be properly defined unless
   !             a gamma file for mt is on endf tape.
   !             special values allowed--
   !               303   non-elastic (all but mt2)
   !               304   inelastic (mt51 thru 91)
   !               318   fission (mt18 or mt19, 20, 21, 38)
   !               401   disappearance (mt102 thru 120)
   !               442   total photon ev-barns
   !               443   total kinematic kerma (high limit)
   !             damage energy production values--
   !               444   total
   !               445   elastic (mt2)
   !               446   inelastic (mt51 thru 91)
   !               447   disappearance (mt102 thru 120)
   !          cards 4 and 5 for nqa gt 0 only
   ! card 4
   !    mta      mt numbers for users q values
   ! card 5
   !    qa       user specified q values (ev)
   !               (if qa.ge.99.e6, read in variable qbar
   !                  for this reaction)
   ! card 5a     variable qbar (for reactions with qa flag only)
   !    qbar      tab1 record giving qbar versus e (1000 words max)
   !
   !-------------------------------------------------------------------
   use mainio ! provides nsysi,nsyso
   ! internals
   integer::ntemp,local,npkk,i,nz,j,isave,iz,loc,nzq,nz0
   integer::nscr,nend4,nend6,iold,inew,itemp,idone,nb,nw
   integer::Damage_settings, k, l
   real(kr)::time,flag
   integer,parameter::npkmax=28
   integer::mtk(2+npkmax)
   real(kr)::z(17)
   real(kr),dimension(:),allocatable::tmp
   real(kr),parameter::qtest=99.e6_kr
   real(kr),parameter::qflag=-1.e-9_kr
   real(kr),parameter::zero=0
   integer,parameter::maxqbar=10000

   !--start

   !--read user input.
   nplot=0
   read(nsysi,*) nendf,nin,nout,nplot
   read(nsysi,*) matd,npk,nqa,ntemp,local,iprint,break, Damage_settings

   !-- Damage settings card
   if (Damage_settings.gt.0) then
         read(nsysi,*) icntrl(1), icntrl(2), icntrl(3)
   end if
   
   !-- Threshold displacement energy treatments
   select case (icntrl(1))
      
      case(1)  
         break_i = 0.0 ! No threshold
      case(2)  
         break_i = 2.0*break ! original 3-level Kinchin-Pease damage energy 
      case(3)
         break_i = 2.0*break/0.8 ! NRT damage energy
      case(4)
         break_i = break ! sharp transition Kinchin-Pease damage energy 

   end select

   !-- Efficiency/Partition mode
   select case (icntrl(2))
      
   case(0)  

   case(1)  

   case(3)

   case(4)

   end select

   do i = 1,10
      do j=1,10
         do k=1,10

         print*, df(REAL(i*10,8),REAL(j,8),REAL(j,8),REAL(k,8),REAL(k,8))

         end do
      end do
   end do

   end subroutine heatr

   real(kr) function df(e,zr,ar,zl,al)
   !-------------------------------------------------------------------
   ! Damage function using the Lindhard partition of
   ! energy between atomic and electronic motion.
   ! Call with e=0 for each reaction to precompute the constants.
   !-------------------------------------------------------------------
   ! externals
   real(kr)::e,zr,ar,zl,al
   ! internals
   real(kr)::el,rel,denom,fl,ep,dam,damage_energy
   real(kr)::threshold_factor, critical_mass
   real(kr),parameter::twothd=.666666667e0_kr
   real(kr),parameter::threeq=.75e0_kr
   real(kr),parameter::sixth=.166666667e0_kr
   real(kr),parameter::onep5=1.5e0_kr
   real(kr),parameter::c1=30.724e0_kr
   real(kr),parameter::c2=.0793e0_kr
   real(kr),parameter::c3=3.4008e0_kr
   real(kr),parameter::c4=.40244e0_kr
   real(kr),parameter::zero=0
   save rel,fl

   if (zr.eq.zero) then
      df=0
   else if (e.le.zero) then
      el=c1*zr*zl*sqrt(zr**twothd+zl**twothd)*(ar+al)/al
      rel=1/el
      denom=(zr**twothd+zl**twothd)**threeq*ar**onep5*sqrt(al)
      fl=c2*zr**twothd*sqrt(zl)*(ar+al)**onep5/denom
      df=0
   else if (e.lt.break) then
      df=0
   else
      ep=e*rel
      dam=e/(1+fl*(c3*ep**sixth+c4*ep**threeq+ep))
      df=dam
   endif





   !option to inhibit or isolate low recoil mass contributions to dam
   if ( icntrl(3) .gt. 0) then
      critical_mass = al - icntrl(3)*1.0
      if (ar .le. critical_mass) then
         dam = 0.0
         df = 0.0
      endif
   else if ( icntrl(3) .lt. 0) then
      critical_mass = al + icntrl(3)*1.0
      if (ar .ge. critical_mass) then
         dam = 0.0
         df = 0.0
      endif
   endif

   damage_energy = df

   !-- Applying various damage threshold treatments

   select case(icntrl(1))
      
      case(0) ! Default

      case(1) ! No threshold

      case(2) ! Original 3-level Kinchin-Pease damage energy 
         if ( damage_energy .lt. break) then
            dam = 0.0
            df  = 0.0
            threshold_factor = 1.0
         elseif ( damage_energy .lt. 2.0*break) then
            if ( dam .le. 0.0) then
              threshold_factor = 0.0
            else
              threshold_factor = 2.0*break/dam
            endif
            df = threshold_factor*dam
         else
            threshold_factor = 1.0
            df = threshold_factor*dam
         endif
      
      case(3) ! NRT damage energy

         if ( damage_energy .lt. break) then
            dam = 0.0
            df  = 0.0
            threshold_factor = 1.0
         elseif ( damage_energy .lt. break_i) then
            if ( dam .le. 0.0) then
              threshold_factor = 1.0
              df = 0.0
              dam = 0.0
            else
              threshold_factor = 2.0*break/0.8/dam
            endif
            df = threshold_factor*dam
            dam = df
         else 
            threshold_factor = 1.0
            df = threshold_factor*dam
            dam = df
         endif

      case(4) ! sharp transition Kinchin-Pease damage energy 
         if ( damage_energy .lt. break) then
            dam = 0.0
            df  = 0.0
            threshold_factor = 1.0
         else
            threshold_factor = 1.0
            df = threshold_factor*dam
         endif

   end select

   return
   end function df

end module heatm