! Equations module for dark energy with constant equation of state parameter w
! allowing for perturbations based on a quintessence model
! by Antony Lewis (http://cosmologist.info/)

! Dec 2003, fixed (fatal) bug in tensor neutrino setup
! Changes to tight coupling approximation
! June 2004, fixed problem with large scale polarized tensors; support for vector modes
! Generate vector modes on their own. The power spectrum is taken from the scalar parameters.
! August 2004, fixed reionization term in lensing potential
! Nov 2004, change massive neutrino l_max to be consistent with massless if light
! Apr 2005, added DoLateRadTruncation option
! June 2006, added support for arbitrary neutrino mass splittings
! Nov 2006, tweak to high_precision transfer function accuracy at lowish k
! June 2011, improved radiation approximations from arXiv: 1104.2933; Some 2nd order tight coupling terms
!            merged fderivs and derivs so flat and non-flat use same equations; more precomputed arrays
!            optimized neutrino sampling, and reorganised neutrino integration functions
! Feb 2012, updated PPF version but now only simple case for w, w_a (no anisotropic stresses etc)
! Feb 2013: fixed various issues with accuracy at larger neutrino masses
! Oct 2013: fix PPF, consistent with updated equations_cross

module LambdaGeneral
  use precision
  use ModelParams
  implicit none

  real(dl)  :: w_lam = -1_dl !p/rho for the dark energy (assumed constant)
  ! w_lam is now w0
  !comoving sound speed. Always exactly 1 for quintessence
  !(otherwise assumed constant, though this is almost certainly unrealistic)
  real(dl) :: cs2_lam = 1_dl
  !cs2_lam now is ce^2

  logical :: use_tabulated_w = .false.
  real(dl) :: wa_ppf = 0._dl
  real(dl) :: c_Gamma_ppf = 0.4_dl
  integer, parameter :: nwmax = 5000, nde = 2000
  integer :: nw_ppf
  real(dl) w_ppf(nwmax), a_ppf(nwmax), ddw_ppf(nwmax)
  real(dl) rde(nde),ade(nde),ddrde(nde)
  real(dl), parameter :: amin = 1.d-9
  logical :: is_cosmological_constant

  private nde,ddw_ppf,rde,ade,ddrde,amin
contains

  subroutine DarkEnergy_ReadParams(Ini)
    use IniFile
    Type(TIniFile) :: Ini
    character(LEN=Ini_max_string_len) wafile
    integer i

    if (Ini_HasKey_File(Ini,'usew0wa')) then
       stop 'input variables changed from usew0wa: now use_tabulated_w or w, wa'
    end if

    use_tabulated_w = Ini_Read_Logical_File(Ini,'use_tabulated_w',.false.)
    if(.not. use_tabulated_w)then
       w_lam = Ini_Read_Double_File(Ini,'w', -1.d0)
       wa_ppf = Ini_Read_Double_File(Ini,'wa', 0.d0)
       if (Feedback >0) write(*,'("(w0, wa) = (", f8.5,", ", f8.5, ")")') w_lam,wa_ppf
    else
       wafile = Ini_Read_String_File(Ini,'wafile')
       open(unit=10,file=wafile,status='old')
       nw_ppf=0
       do i=1,nwmax+1
          read(10,*,end=100)a_ppf(i),w_ppf(i)
          a_ppf(i)=dlog(a_ppf(i))
          nw_ppf=nw_ppf+1
       enddo
       write(*,'("Note: ", a, " has more than ", I8, " data points")') trim(wafile), nwmax
       write(*,*)'Increase nwmax in LambdaGeneral'

       stop
100    close(10)
       write(*,'("read in ", I8, " (a, w) data points from ", a)') nw_ppf, trim(wafile)
       call setddwa
       call interpolrde
    endif
    cs2_lam = Ini_Read_Double_File(Ini,'cs2_lam',1.d0)

    call setcgammappf

  end subroutine DarkEnergy_ReadParams


  subroutine setddwa
    real(dl), parameter :: wlo=1.d30, whi=1.d30

    call spline(a_ppf,w_ppf,nw_ppf,wlo,whi,ddw_ppf) !a_ppf is lna here

  end subroutine setddwa


  function w_de(a)
    real(dl) :: w_de, al
    real(dl), intent(IN) :: a

    if(.not. use_tabulated_w) then
       w_de=w_lam+wa_ppf*(1._dl-a)
    else
       al=dlog(a)
       if(al.lt.a_ppf(1)) then
          w_de=w_ppf(1)                   !if a < minimum a from wa.dat
       elseif(al.gt.a_ppf(nw_ppf)) then
          w_de=w_ppf(nw_ppf)              !if a > maximus a from wa.dat
       else
          call cubicsplint(a_ppf,w_ppf,ddw_ppf,nw_ppf,al,w_de)
       endif
    endif
  end function w_de  ! equation of state of the PPF DE


  function drdlna_de(al)
    real(dl) :: drdlna_de, a
    real(dl), intent(IN) :: al

    a=dexp(al)
    drdlna_de=3._dl*(1._dl+w_de(a))

  end function drdlna_de


  subroutine interpolrde
    real(dl), parameter :: rlo=1.d30, rhi=1.d30
    real(dl) :: atol, almin, al, rombint, fint
    integer :: i
    external rombint
    atol=1.d-5
    almin=dlog(amin)
    do i=1,nde
       al=almin-almin/(nde-1)*(i-1)    !interpolate between amin and today
       fint=rombint(drdlna_de, al, 0._dl, atol)+4._dl*al
       ade(i)=al
       rde(i)=dexp(fint) !rho_de*a^4 normalize to its value at today
    enddo
    call spline(ade,rde,nde,rlo,rhi,ddrde)
  end subroutine interpolrde

  function grho_de(a)  !8 pi G a^4 rho_de
    real(dl) :: grho_de, al, fint
    real(dl), intent(IN) :: a
    external rombint

    if(.not. use_tabulated_w) then
       grho_de=grhov*a**(1._dl-3.*w_lam-3.*wa_ppf)*exp(-3.*wa_ppf*(1._dl-a))
    else
       if(a.eq.0.d0)then
          grho_de=0.d0      !assume rho_de*a^4-->0, when a-->0, OK if w_de always <0.
       else
          al=dlog(a)
          if(al.lt.ade(1))then
             fint=rde(1)*(a/amin)**(1.-3.*w_de(amin))    !if a<amin, assume here w=w_de(amin)
          else              !if amin is small enough, this extrapolation will be unnecessary.
             call cubicsplint(ade,rde,ddrde,nde,al,fint)
          endif
          grho_de=grhov*fint
       endif
    endif
  end function grho_de


  !-------------------------------------------------------------------
  SUBROUTINE cubicsplint(xa,ya,y2a,n,x,y)
    INTEGER n
    real(dl) x,y,xa(n),y2a(n),ya(n)
    INTEGER k,khi,klo
    real(dl)a,b,h
    klo=1
    khi=n
1   if (khi-klo.gt.1) then
       k=(khi+klo)/2
       if(xa(k).gt.x)then
          khi=k
       else
          klo=k
       endif
       goto 1
    endif
    h=xa(khi)-xa(klo)
    if (h.eq.0.) stop 'bad xa input in splint'
    a=(xa(khi)-x)/h
    b=(x-xa(klo))/h
    y=a*ya(klo)+b*ya(khi)+&
         ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.d0
  END SUBROUTINE cubicsplint
  !--------------------------------------------------------------------


  subroutine setcgammappf

    c_Gamma_ppf=0.4d0*sqrt(cs2_lam)

  end subroutine setcgammappf


end module LambdaGeneral

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!Return OmegaK - modify this if you add extra fluid components
!RL commented out 032724- this is physically incorrect and we should modify H0 instead. The curvature should always be determined at the input level
!function GetOmegak()
!  use precision
!  use ModelParams
  !DG May 25 2015 Added so that massless neutrinos can be self consistently included
  !negligible for most cosmologies but included for completenes
!  use constants 
!  real(dl) nnu
!  real(dl)  GetOmegak,rhocrit,omegah2_rad

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !DG May 25 2015 Added so that massless neutrinos can be self consistently included
  !negligible for most cosmologies but included for completeness
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !Self-consistent inclusion of massive neutrinos in computation of curvature
  !different from regular CAMB
  !RL: Although this routine is in equations_ppf, from the DG comments here it seems that this particular function is an add-on (i.e. it's not CAMB-intrinsic, so I still need to fix)
!  rhocrit=(8.0d0*const_pi*G*1.d3/(3.0d0*((1.d7/(MPC_in_sec*c*1.d2))**(2.0d0))))**(-1.0d0)
!  omegah2_rad=((CP%TCMB**4.0d0)/(rhocrit))/(c**2.0d0) !RL replaced the COBE value with the input value
!  omegah2_rad=omegah2_rad*a_rad*1.d1/(1.d4)
!  nnu=CP%Num_Nu_massless
!  omegah2_rad=omegah2_rad+(nnu*grhor*(c**2.0d0)/((1.d5**2.0d0)))/3.0d0 
!  GetOmegak = 1.0d0 - (CP%omegab+CP%omegac+CP%omegav+CP%omegan+CP%omegaax)-omegah2_rad/((CP%H0/1.d2)**2.0d0)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!end function GetOmegak

subroutine init_background
  use LambdaGeneral
  use ModelParams !RL
  integer :: i_check
  !This is only called once per model, and is a good point to do any extra initialization.
  !It is called before first call to dtauda, but after massive neutrinos are initialized and after GetOmegak
  is_cosmological_constant = .not. use_tabulated_w .and. w_lam==-1_dl .and. wa_ppf==0._dl
  
  CP%tau0=TimeOfz(0._dl) !RL 061924 moved tau0 here so to initialte tau_osc

  if (CP%a_osc .gt. 1._dl) then
     CP%tau_osc=CP%tau0 + 1._dl !To make sure no switch happens before the present day. 
  else
     CP%tau_osc=DeltaTime(0._dl, CP%a_osc, in_tol = 1.0d-8)
  end if
  !RL 062624 - add scenarios  potential problems that the background switches but the perturbation does not
  if (CP%a_osc .le. 1._dl .and. CP%tau_osc .gt. CP%tau0) then
     write(*, *) 'a_osc <= 1 and tau_osc > tau0. Rarely happens, but likely a_osc is too &
&close to 1 that there is an issue related to the tolerance of the tau integration. &
&Taking tau_osc as the smaller value (tau0) instead.'
     CP%tau_osc = min(CP%tau_osc, CP%tau0)
  end if

end  subroutine init_background


!Background evolution
  !RL 07312023 - constructing a background evolution function
function grhoax_frac(a_in)
  !rho_ax(a)/rho_crit.
  use precision
  use constants
  use ModelParams
  implicit none
  real(dl) grhoax_frac
  real(dl), intent(IN) :: a_in
  real(dl) a, a_min
  real(dl) a2, v1_bg, v2_bg, grhoaxh2_ov_grhom, grhoaxh2_ov_grhom_test, wcorr_coeff !RL added grhoax_kg, wcorr_coeff
  !RL 010625
  a = a_in
  a_min = 10._dl**(loga_table(1))

  a2 = a**2.0d0
  
  if (a .lt. CP%a_osc) then
      if (a .lt. a_min) then
         grhoaxh2_ov_grhom = rhoaxh2ovrhom_logtable(1)
      else         
         !Note that in the background the spline table is log10(rho)
         call spline_out(loga_table,rhoaxh2ovrhom_logtable,rhoaxh2ovrhom_logtable_buff,ntable,dlog10(a),grhoaxh2_ov_grhom)
      end if
      
     grhoax_frac = (10._dl**grhoaxh2_ov_grhom)/(CP%H0**2.0d0/1.0d4)        

  else
     !RL: Initialize the unused variables to be safe
     v1_bg = 0.0d0
     v2_bg = 0.0d0
     grhoaxh2_ov_grhom = 0.0d0
     !w correction to the background
     !wcorr_coeff = CP%ah_osc*CP%a_osc/((CP%ma/CP%H0_eV)*(CP%H0/100.0d0))
     wcorr_coeff = CP%ahosc_ETA*CP%a_osc/((CP%ma/CP%H0_eV)*(CP%H0/100.0d0)) !RL082924

     grhoax_frac=(CP%rhorefp_ovh2)*((CP%a_osc/a)**3.0d0)*dexp((wcorr_coeff**2.0d0)*3.0d0*CP%wEFA_c*(1.0d0/(a2**2.0d0) &
          &- 1.0d0/(CP%a_osc**4.0d0))/4.0d0)
  endif
  
end function grhoax_frac


function dtauda(a)
  !get d tau / d a
  use precision
  use constants
  use ModelParams
  use MassiveNu
  use LambdaGeneral
  implicit none
  real(dl) dtauda, grhoax_frac !RL
  real(dl), intent(IN) :: a
  real(dl) rhonu,grhoa2, a2, v1_test, v2_test, grhotest 
  integer i
  integer nu_i
  external grhoax_frac

  a2=a**2._dl

  !  8*pi*G*rho*a**4.
  grhoa2=grhok*a2+(grhoc+grhob)*a+grhog+grhornomass
  if (is_cosmological_constant) then
     grhoa2=grhoa2+grhov*a2**2._dl
  else
     grhoa2=grhoa2+ grho_de(a)
  end if

  if (CP%Num_Nu_massive /= 0) then
     !Get massive neutrino density relative to massless
     do nu_i = 1, CP%nu_mass_eigenstates
        call Nu_rho(a*nu_masses(nu_i),rhonu)
        grhoa2=grhoa2+rhonu*grhormass(nu_i)
     end do
  end if
  grhoa2 = grhoa2 + grhoax_frac(a)*grhom*(a2**2._dl)
  dtauda=sqrt(3._dl/grhoa2)


end function dtauda

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!Gauge-dependent perturbation equations

module GaugeInterface
  use precision
  use ModelParams
  use MassiveNu
  use LambdaGeneral
  use Errors
  implicit none
  public

  !Description of this file. Change if you make modifications.
  character(LEN=*), parameter :: Eqns_name = 'equations_ppf-Oct13'

  integer, parameter :: basic_num_eqns = 5

  logical :: DoTensorNeutrinos = .false.

  logical :: DoLateRadTruncation = .true.!, RL testing with .false. 07/11/2023 - it seems that the setting here is overwritten by what's defined by inidriver read from the inifile
  !if true, use smooth approx to radiation perturbations after decoupling on
  !small scales, saving evolution of irrelevant osciallatory multipole equations

  logical, parameter :: second_order_tightcoupling = .true.

  real(dl) :: Magnetic = 0._dl
  !Vector mode anisotropic stress in units of rho_gamma
  real(dl) :: vec_sig0 = 1._dl
  !Vector mode shear
  integer, parameter :: max_l_evolve = 512 !Maximum l we are ever likely propagate

  !Supported scalar initial condition flags
  integer, parameter :: initial_adiabatic=1, initial_iso_CDM=2, &
       initial_iso_baryon=3,  initial_iso_neutrino=4, initial_iso_neutrino_vel=5, initial_vector = 0, &
                                !!!!!!!!!!!!!!!!!!!!!!!!!
                                !!Isocurvature initial condition flags
       initial_iso_axion=6 
  integer, parameter :: initial_nummodes =  initial_iso_axion !DM: added axion isocurvature
!!!!!!!!!!!!!!!!!!!!
  type EvolutionVars
     real(dl) q, q2
     real(dl) k_buf,k2_buf ! set in initial

     integer w_ix !Index of two quintessence equations
     integer r_ix !Index of the massless neutrino hierarchy
     integer g_ix !Index of the photon neutrino hierarchy

!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!Axion variables
!!!!!Axion fluid indices

     integer a_ix !Index of the two axion fluid equations
     integer a_kg_ix !RL: Index of the two axion field equations - and will change to fluid equations after the switch   

     integer q_ix !index into q_evolve array that gives the value q
     logical TransferOnly

     !       nvar  - number of scalar (tensor) equations for this k
     integer nvar,nvart, nvarv

     !Max_l for the various hierarchies
     integer lmaxg,lmaxnr,lmaxnu,lmaxgpol,MaxlNeeded
     integer lmaxnrt, lmaxnut, lmaxt, lmaxpolt, MaxlNeededt
     logical EvolveTensorMassiveNu(max_nu)
     integer lmaxnrv, lmaxv, lmaxpolv

     integer polind  !index into scalar array of polarization hierarchy

     !array indices for massive neutrino equations
     integer nu_ix(max_nu), nu_pert_ix
     integer nq(max_nu), lmaxnu_pert
     logical has_nu_relativistic

     !Initial values for massive neutrino v*3 variables calculated when switching
     !to non-relativistic approx
     real(dl) G11(max_nu),G30(max_nu)
     !True when using non-relativistic approximation
     logical MassiveNuApprox(max_nu)
     real(dl) MassiveNuApproxTime(max_nu)

     !True when truncating at l=2,3 when k*tau>>1 (see arXiv:1104.2933)
     logical high_ktau_neutrino_approx

     !Massive neutrino scheme being used at the moment
     integer NuMethod

     !True when using tight-coupling approximation (required for stability at early times)
     logical TightCoupling, TensTightCoupling
     real(dl) TightSwitchoffTime

     !Numer of scalar equations we are propagating
     integer ScalEqsToPropagate
     integer TensEqsToPropagate
     !beta > l for closed models
     integer FirstZerolForBeta
     !Tensor vars
     real(dl) aux_buf

     real(dl) pig, pigdot !For tight coupling
     real(dl) poltruncfac

     !PPF parameters
     real(dl) dgrho_e_ppf, dgq_e_ppf
     real(dl) dgrhoec_ppf, dgqec_ppf, vTc_ppf 

     logical no_nu_multpoles, no_phot_multpoles
     integer lmaxnu_tau(max_nu)  !lmax for massive neutinos at time being integrated
     logical nu_nonrelativistic(max_nu)
     logical oscillation_started !RL: whether the axion oscillation has started yet
     logical output_done !RL testing (to enable output to the screen only once, right after oscillation started)
     real(dl) metric_delta(2) !RL 090323 adding boundary condition delta function
     real(dl) renorm_c !RL 050724

     real(dl) denlk(max_l_evolve),denlk2(max_l_evolve), polfack(max_l_evolve)
     real(dl) Kf(max_l_evolve)

     integer E_ix, B_ix !tensor polarization indices
     real(dl) denlkt(4,max_l_evolve),Kft(max_l_evolve)
  end type EvolutionVars

  !precalculated arrays
  real(dl) polfac(max_l_evolve),denl(max_l_evolve),vecfac(max_l_evolve),vecfacpol(max_l_evolve)

  real(dl), parameter :: ep0=1.0d-2
  integer, parameter :: lmaxnu_high_ktau=3

  real(dl) epsw
  real(dl) nu_tau_notmassless(nqmax0+1,max_nu), nu_tau_nonrelativistic(max_nu),nu_tau_massive(max_nu)
contains

  subroutine GaugeInterface_ScalEv(EV,y,tau,tauend,tol1,ind,c,w)
    type(EvolutionVars) EV
    real(dl) c(24),w(EV%nvar,9), y(EV%nvar), tol1, tau, tauend, tau_osc
    integer ind
    !write(*, *), 'GaugeInterface_ScalEv, EV%q, tau, tauend, c(1)', EV%q, tau, tauend, c(1)
    call dverk(EV,EV%ScalEqsToPropagate,derivs,tau,y,tauend,tol1,ind,c,EV%nvar,w) !
    !!write(*, *) 'GaugeInterface_ScalEv, EV%nvar, EV%a_kg_ix, EV%a_kg_ix + 1, y(EV%a_kg_ix), y(EV%a_kg_ix+1), y(1:48)', EV%nvar, EV%a_kg_ix, EV%a_kg_ix + 1, y(EV%a_kg_ix), y(EV%a_kg_ix+1), y(1:48)
    
    !RL added the 1e-3 trying to see if decreasing the tolerance will suppress the spike due to integrating to tau_switch
    !write(*, *) 'Is c(1:9) always zero?', c(1:9) == 0._dl
    !RL 050224
    !!if (tau .lt. TimeSteps%points(2) .and. c(1) == 0) then
    !!   write(*, *) 'c(1) reset to 0 but still before CMB output, c(1), tau', c(1), tau
    !!   ind=2
    !!   c(1) = 2
    !!end if
    
    if (ind==-3) then
       call GlobalError('Dverk error -3: the subroutine was unable  to  satisfy  the  error ' &
            //'requirement  with a particular step-size that is less than or * ' &
            //'equal to hmin, which may mean that tol is too small' &
            //'--- but most likely you''ve messed up the y array indexing; ' &
            //'compiling with bounds checking may (or may not) help find the problem.',error_evolution)
    end if
  end subroutine GaugeInterface_ScalEv

  function next_nu_nq(nq) result (next_nq)
    integer, intent(in) :: nq
    integer q, next_nq

    if (nq==0) then
       next_nq=1
    else
       q = nu_q(nq)
       if (q>=10) then
          next_nq = nqmax
       else
          next_nq = nq+1
       end if
    end if

  end function next_nu_nq

  recursive subroutine GaugeInterface_EvolveScal(EV,tau,y,tauend,tol1,ind,c,w)
    use ThermoData
    !use lvalues
    !use constants
    use ModelData
    type(EvolutionVars) EV, EVout
    real(dl) c(24),w(EV%nvar,9), y(EV%nvar), yout(EV%nvar), tol1, tau, tauend
    integer ind, nu_i
    real(dl) cs2, opacity, dopacity
    real(dl) tau_switch_ktau, tau_switch_nu_massless, tau_switch_nu_massive, next_switch
    real(dl) tau_switch_no_nu_multpoles, tau_switch_no_phot_multpoles,tau_switch_nu_nonrel
    real(dl) tau_switch_oscillation !RL adding axion oscillation switch time
    real(dl) noSwitch, smallTime
    real(dl) sources(CTransScal%NumSources)
    real(dl) yprimetest(EV%nvar) !RL testing 09012023
    !!write(*, *) 'In GaugeInterface_EvolveScal, EV%q, tau, tauend, c(1)', EV%q, tau, tauend, c(1)

    noSwitch= CP%tau0+1
    !smallTime =  min(tau, 1/EV%k_buf)/100 !The original smallTime - turned on again on 11/07/2022 for testing
    smallTime = 0._dl !RL killing smallTime 

    tau_switch_ktau = noSwitch
    tau_switch_no_nu_multpoles= noSwitch
    tau_switch_no_phot_multpoles= noSwitch

    !Massive neutrino switches
    tau_switch_nu_massless = noSwitch
    tau_switch_nu_nonrel = noSwitch
    tau_switch_nu_massive= noSwitch

    !RL: axion oscillation switches - first it's noSwitch
    tau_switch_oscillation = noSwitch

    !Evolve equations from tau to tauend, performing switches in equations if necessary.

    !RL: deciding whether oscillation has already started. If it has not, then the tau_switch will be the tau_osc. If it has (i.e. we have passed the oscillation switch), then it remains noSwitch (note that this ScalarEv subroutine is recursive and is called over and over again until some condition is met)
    !RL: tip: checked that in EvolutionVars, if you first assign a variable it's 0 (False) by default until you assign it to True
    if (.not. EV%oscillation_started) then
       tau_switch_oscillation = CP%tau_osc  
    end if
    
    !RL hacking to temprarily comment this out 07/11/2023-----------
    if (.not. EV%high_ktau_neutrino_approx .and. .not. EV%no_nu_multpoles ) then
       tau_switch_ktau=  max(20, EV%lmaxnr-4)/EV%k_buf
    end if
    !RL-----------

    if (CP%Num_Nu_massive /= 0) then
       do nu_i = 1, CP%Nu_mass_eigenstates
          if (EV%nq(nu_i) /= nqmax) then
             tau_switch_nu_massless = min(tau_switch_nu_massless,nu_tau_notmassless(next_nu_nq(EV%nq(nu_i)),nu_i))
          else if (.not. EV%nu_nonrelativistic(nu_i)) then
             tau_switch_nu_nonrel = min(nu_tau_nonrelativistic(nu_i),tau_switch_nu_nonrel)
          else if (EV%NuMethod==Nu_trunc .and..not. EV%MassiveNuApprox(nu_i)) then
             tau_switch_nu_massive = min(tau_switch_nu_massive,EV%MassiveNuApproxTime(nu_i))
          end if
       end do
    end if

    if (DoLateRadTruncation) then

       if (.not. EV%no_nu_multpoles) & !!.and. .not. EV%has_nu_relativistic .and. tau_switch_nu_massless ==noSwitch)  &
            tau_switch_no_nu_multpoles=max(15/EV%k_buf*AccuracyBoost,min(taurend,matter_verydom_tau))

       if (.not. EV%no_phot_multpoles .and. (.not.CP%WantCls .or. EV%k_buf>0.03*AccuracyBoost)) &
            tau_switch_no_phot_multpoles =max(15/EV%k_buf,taurend)*AccuracyBoost
    end if

    next_switch = min(tau_switch_oscillation, tau_switch_ktau, tau_switch_nu_massless,EV%TightSwitchoffTime, &
         tau_switch_nu_massive, tau_switch_no_nu_multpoles, tau_switch_no_phot_multpoles, &
         tau_switch_nu_nonrel, noSwitch)

    if (next_switch < tauend) then
       if (next_switch > tau+smallTime) then
          call GaugeInterface_ScalEv(EV, y, tau,next_switch,tol1,ind,c,w)
          if (global_error_flag/=0) return
       end if

       EVout=EV

       if (next_switch == tau_switch_oscillation) then
          EVout%oscillation_started = .true.
          call SetupScalarArrayIndices(EVout)
          !write(*, *) 'In GaugeInterface, right before CopyScalarVariableArray, tau, CP%tau_osc, their fractional difference'
          !write(*, '(36e52.42,\)') tau, CP%tau_osc, tau/CP%tau_osc - 1.0d0
          call CopyScalarVariableArray(y,yout, EV, EVout)
          EV = EVout
          y = yout
          
          ind = 1 
          
       else if (next_switch == EV%TightSwitchoffTime) then
          !TightCoupling
          EVout%TightCoupling=.false.
          EVout%TightSwitchoffTime = noSwitch
          call SetupScalarArrayIndices(EVout)
          call CopyScalarVariableArray(y,yout, EV, EVout)
          EV=EVout
          y=yout
          ind=1
          !Set up variables with their tight coupling values
          y(EV%g_ix+2) = EV%pig
          call thermo(tau,cs2,opacity,dopacity)

          if (second_order_tightcoupling) then
             ! Francis-Yan Cyr-Racine November 2010

             y(EV%g_ix+3) = (3._dl/7._dl)*y(EV%g_ix+2)*(EV%k_buf/opacity)*(1._dl+dopacity/opacity**2) + &
                  (3._dl/7._dl)*EV%pigdot*(EV%k_buf/opacity**2)*(-1._dl)

             y(EV%polind+2) = EV%pig/4 + EV%pigdot*(1._dl/opacity)*(-5._dl/8._dl- &
                  (25._dl/16._dl)*dopacity/opacity**2) + &
                  EV%pig*(EV%k_buf/opacity)**2*(-5._dl/56._dl)
             y(EV%polind+3) = (3._dl/7._dl)*(EV%k_buf/opacity)*y(EV%polind+2)*(1._dl + &
                  dopacity/opacity**2) + (3._dl/7._dl)*(EV%k_buf/opacity**2)*((EV%pigdot/4._dl)* &
                  (1._dl+(5._dl/2._dl)*dopacity/opacity**2))*(-1._dl)
          else
             y(EV%g_ix+3) = 3./7*y(EV%g_ix+2)*EV%k_buf/opacity
             y(EV%polind+2) = EV%pig/4
             y(EV%polind+3) =y(EV%g_ix+3)/4
          end if
          !write(*,*) 'c(1) reset 3', c(1)
       else if (next_switch==tau_switch_ktau) then
          !k tau >> 1, evolve massless neutrino effective fluid up to l=2
          !write(*, *) 'RL: Is EV%high_ktau_neutrino_approx ever .true.?', EVout%high_ktau_neutrino_approx
          EVout%high_ktau_neutrino_approx=.true.
          EV%nq(1:CP%Nu_mass_eigenstates) = nqmax
          call SetupScalarArrayIndices(EVout)
          call CopyScalarVariableArray(y,yout, EV, EVout)
          y=yout
          EV=EVout
          !write(*, *) 'RL: Is EV%high_ktau_neutrino_approx .true. after this?', EV%high_ktau_neutrino_approx
          !write(*,*) 'c(1) reset 4', c(1)
       else if (next_switch == tau_switch_nu_massless) then
          !Mass starts to become important, start evolving next momentum mode
          do nu_i = 1, CP%Nu_mass_eigenstates
             if (EV%nq(nu_i) /= nqmax .and. &
                  next_switch == nu_tau_notmassless(next_nu_nq(EV%nq(nu_i)),nu_i)) then
                EVOut%nq(nu_i) = next_nu_nq(EV%nq(nu_i))
                call SetupScalarArrayIndices(EVout)
                call CopyScalarVariableArray(y,yout, EV, EVout)
                EV=EVout
                y=yout
                exit
             end if
          end do
          !write(*,*) 'c(1) reset 5', c(1)
       else if (next_switch == tau_switch_nu_nonrel) then
          !Neutrino becomes non-relativistic, don't need high L
          do nu_i = 1, CP%Nu_mass_eigenstates
             if (.not. EV%nu_nonrelativistic(nu_i) .and.  next_switch==nu_tau_nonrelativistic(nu_i) ) then
                EVout%nu_nonrelativistic(nu_i)=.true.
                call SetupScalarArrayIndices(EVout)
                call CopyScalarVariableArray(y,yout, EV, EVout)
                EV=EVout
                y=yout
                exit
             end if
          end do
          !write(*,*) 'c(1) reset 6', c(1)
       else if (next_switch == tau_switch_nu_massive) then
          !Very non-relativistic neutrinos, switch to truncated velocity-weight hierarchy
          do nu_i = 1, CP%Nu_mass_eigenstates
             if (.not. EV%MassiveNuApprox(nu_i) .and. next_switch== EV%MassiveNuApproxTime(nu_i) ) then
                call SwitchToMassiveNuApprox(EV,y, nu_i)
                exit
             end if
          end do
          !write(*,*) 'c(1) reset 7', c(1)
       else if (next_switch==tau_switch_no_nu_multpoles) then
          !Turn off neutrino hierarchies at late time where slow and not needed.
          ind=1
          EVout%no_nu_multpoles=.true.
          EVout%nq(1:CP%Nu_mass_eigenstates ) = nqmax
          call SetupScalarArrayIndices(EVout)
          call CopyScalarVariableArray(y,yout, EV, EVout)
          y=yout
          EV=EVout
          !write(*, *) 'Is EV%no_nu_multpoles still activated?', EV%no_nu_multpoles
          !write(*,*) 'c(1) reset 8', c(1)
       else if (next_switch==tau_switch_no_phot_multpoles) then
          !Turn off photon hierarchies at late time where slow and not needed.
          ind=1
          EVout%no_phot_multpoles=.true.
          call SetupScalarArrayIndices(EVout)
          call CopyScalarVariableArray(y,yout, EV, EVout)
          y=yout
          EV=EVout
         ! write(*,*) 'c(1) reset 9', c(1)
       end if
       !write(*, *) 'Before GaugeInterface_EvolveScal called again recursively, c(1)', c(1)
       call GaugeInterface_EvolveScal(EV,tau,y,tauend,tol1,ind,c,w)
       return
    end if

    !write(*, *) 'Before GI_ScalEv call, c(1)', c(1)
    call GaugeInterface_ScalEv(EV,y,tau,tauend,tol1,ind,c,w)
    !write(*, *) 'After GI_ScalEv call, c(1)', c(1)
    !if (EV%oscillation_started == .true. .and. EV%oscillation_output_done == .false.) then
    !write(*, *) 'What is ind after one call of ScalEv?', ind
    !write(*, *) 'What is c here when ind is not reset to 1?'
    !write(*, *) c
    !EV%oscillation_output_done = .true.
    !end if

    !write(*, *) 'Are the sources the same from that in cmbmain?', sources
  end subroutine GaugeInterface_EvolveScal

  subroutine GaugeInterface_EvolveTens(EV,tau,y,tauend,tol1,ind,c,w)
    use ThermoData
    type(EvolutionVars) EV, EVout
    real(dl) c(24),w(EV%nvart,9), y(EV%nvart),yout(EV%nvart), tol1, tau, tauend
    integer ind
    real(dl) opacity, cs2

    !<<<<<<<<<< TENSOR TAU_OSC SWITCH FIX -- BLOCK 1 of 2 -- BEGIN <<<<<<<<<<
    !Retained for the tensor support planned in AxiECAMB 2.0; unreachable in this
    !version, which refuses get_tensor_cls = T. The original subroutine body is
    !preserved verbatim between the two stamped blocks, so deleting both blocks
    !restores the pre-fix code exactly. The fix only affects get_tensor_cls = T runs
    !with an oscillating axion, which previously crashed (Error in dverk at tau_osc):
    !only the scalar evolution set EV%oscillation_started, so derivst kept reading the
    !KG background tables (which end at a_osc) past the switch. Here we stop dverk
    !exactly at CP%tau_osc (tau-space, consistent with the scalar switch), set the flag,
    !and restart integration (ind = 1) across the discontinuous RHS. No variable
    !remapping is needed: the axion enters the tensor equations only through the
    !background density.
    !BLOCK 1 handles tau_osc arriving before the tight-coupling switch-off
    !(or with no tight-coupling switch pending in this interval).
    if (.not. EV%oscillation_started .and. CP%tau_osc < tauend .and. &
         .not. (EV%TensTightCoupling .and. tauend > EV%TightSwitchoffTime &
         .and. EV%TightSwitchoffTime <= CP%tau_osc)) then
       if (CP%tau_osc > tau) then
          call dverk(EV,EV%TensEqsToPropagate, derivst,tau,y,CP%tau_osc,tol1,ind,c,EV%nvart,w)
       end if
       EV%oscillation_started = .true.
       ind = 1
    end if
    !>>>>>>>>>> TENSOR TAU_OSC SWITCH FIX -- BLOCK 1 of 2 -- END >>>>>>>>>>

    if (EV%TensTightCoupling .and. tauend > EV%TightSwitchoffTime) then
       if (EV%TightSwitchoffTime > tau) then
          !write(*, *) 'dverk in GaugeInterface_EvolveTens is called' !RL testing
          call dverk(EV,EV%TensEqsToPropagate, derivst,tau,y,EV%TightSwitchoffTime,tol1,ind,c,EV%nvart,w)
       end if
       EVOut=EV
       EVOut%TensTightCoupling = .false.
       call SetupTensorArrayIndices(EVout)
       call CopyTensorVariableArray(y,yout,Ev, Evout)
       Ev = EvOut
       y=yout
       call thermo(tau,cs2,opacity)
       y(EV%g_ix+2)= 32._dl/45._dl*EV%k_buf/opacity*y(3)
       y(EV%E_ix+2) = y(EV%g_ix+2)/4
    end if

    !<<<<<<<<<< TENSOR TAU_OSC SWITCH FIX -- BLOCK 2 of 2 -- BEGIN <<<<<<<<<<
    !BLOCK 2 handles tau_osc arriving after the tight-coupling switch-off
    !(the tight-coupling block above has then already integrated up to its
    !switch time, so we continue from there to tau_osc).
    if (.not. EV%oscillation_started .and. CP%tau_osc < tauend) then
       if (CP%tau_osc > tau) then
          call dverk(EV,EV%TensEqsToPropagate, derivst,tau,y,CP%tau_osc,tol1,ind,c,EV%nvart,w)
       end if
       EV%oscillation_started = .true.
       ind = 1
    end if
    !>>>>>>>>>> TENSOR TAU_OSC SWITCH FIX -- BLOCK 2 of 2 -- END >>>>>>>>>>

    call dverk(EV,EV%TensEqsToPropagate, derivst,tau,y,tauend,tol1,ind,c,EV%nvart,w)


  end subroutine GaugeInterface_EvolveTens

  function DeltaTimeMaxed(a1,a2, tol) result(t)
    real(dl) a1,a2,t
    real(dl), optional :: tol
    if (a1>1._dl) then
       t=0
    elseif (a2 > 1._dl) then
       t = DeltaTime(a1,1.01_dl, tol)
    else
       t= DeltaTime(a1,a2) 
    end if
  end function DeltaTimeMaxed

  subroutine GaugeInterface_Init
    !Precompute various arrays and other things independent of wavenumber
    integer j, nu_i
    real(dl) a_nonrel, a_mass,a_massive, time, nu_mass

    epsw = 100/CP%tau0

    if (CP%WantScalars) then
       do j=2,max_l_evolve
          polfac(j)=real((j+3)*(j-1),dl)/(j+1)
       end do
    end if

    if (CP%WantVectors) then
       do j=2,max_l_evolve
          vecfac(j)=real((j+2),dl)/(j+1)
          vecfacpol(j)=real((j+3)*j,dl)*(j-1)*vecfac(j)/(j+1)**2
       end do
    end if

    do j=1,max_l_evolve
       denl(j)=1._dl/(2*j+1)
    end do

    do nu_i=1, CP%Nu_Mass_eigenstates
       nu_mass = max(0.1_dl,nu_masses(nu_i))
       a_mass =  1.e-1_dl/nu_mass/lAccuracyBoost
       !if (HighAccuracyDefault) a_mass=a_mass/4
       time=DeltaTime(0._dl,nu_q(1)*a_mass)
       nu_tau_notmassless(1, nu_i) = time
       do j=2,nqmax
          !times when each momentum mode becomes signficantly nonrelativistic
          time= time + DeltaTimeMaxed(nu_q(j-1)*a_mass,nu_q(j)*a_mass, 0.01_dl)
          nu_tau_notmassless(j, nu_i) = time
       end do

       a_nonrel =  2.5d0/nu_mass*AccuracyBoost !!!Feb13tweak
       nu_tau_nonrelativistic(nu_i) =DeltaTimeMaxed(0._dl,a_nonrel)
       a_massive =  17.d0/nu_mass*AccuracyBoost
       nu_tau_massive(nu_i) =nu_tau_nonrelativistic(nu_i) + DeltaTimeMaxed(a_nonrel,a_massive)
    end do

  end subroutine GaugeInterface_Init


  subroutine SetupScalarArrayIndices(EV, max_num_eqns)
    !Set up array indices after the lmax have been decided
    use MassiveNu
    !Set the number of equations in each hierarchy, and get total number of equations for this k
    type(EvolutionVars) EV
    integer, intent(out), optional :: max_num_eqns
    integer neq, maxeq, nu_i
    integer neq_dummynu !RL testing

    neq=basic_num_eqns
    maxeq=neq
    if (.not. EV%no_phot_multpoles) then
       !Photon multipoles
       EV%g_ix=basic_num_eqns+1
       if (EV%TightCoupling) then
          neq=neq+2
       else
          neq = neq+ (EV%lmaxg+1)
          !Polarization multipoles
          EV%polind = neq -1 !polind+2 is L=2, for polarizationthe first calculated
          neq=neq + EV%lmaxgpol-1
       end if
    end if
    
    neq_dummynu = neq
    if (.not. EV%no_nu_multpoles) then
       !Massless neutrino multipoles
       EV%r_ix= neq+1
       if (EV%high_ktau_neutrino_approx) then
          neq=neq + 3
       else
          neq=neq + (EV%lmaxnr+1)
       end if
    end if
    !write(*, *) 'EV%no_nu_multpoles, EV%lmaxnr, neq_nu'
    !write(*, *) EV%no_nu_multpoles, EV%lmaxnr, neq -neq_dummynu 
    !write(*, *) 'k'
    !write(*, '(36e52.42)') EV%q
    maxeq = maxeq +  (EV%lmaxg+1)+(EV%lmaxnr+1)+EV%lmaxgpol-1

    !Dark energy
    if (.not. is_cosmological_constant) then
       EV%w_ix = neq+1
       neq=neq+1 !ppf
       maxeq=maxeq+1
    else
       EV%w_ix=0
    end if

    !Massive neutrinos
    if (CP%Num_Nu_massive /= 0) then
       EV%has_nu_relativistic = any(EV%nq(1:CP%Nu_Mass_eigenstates)/=nqmax)
       if (EV%has_nu_relativistic) then
          EV%lmaxnu_pert=EV%lmaxnu
          EV%nu_pert_ix=neq+1
          neq = neq+ EV%lmaxnu_pert+1
          maxeq=maxeq+ EV%lmaxnu_pert+1
       else
          EV%lmaxnu_pert=0
       end if

       do nu_i=1, CP%Nu_Mass_eigenstates
          if (EV%high_ktau_neutrino_approx) then
             if (HighAccuracyDefault .and. CP%WantTransfer .and. EV%q < 1.d0) then
                EV%lmaxnu_tau(nu_i)=max(4,lmaxnu_high_ktau)
             else
                EV%lmaxnu_tau(nu_i)=lmaxnu_high_ktau
             end if
          else
             EV%lmaxnu_tau(nu_i) =max(min(nint(0.8_dl*EV%q*nu_tau_nonrelativistic(nu_i)*lAccuracyBoost),EV%lmaxnu),3)
!!!Feb13tweak
             if (EV%nu_nonrelativistic(nu_i)) EV%lmaxnu_tau(nu_i)=min(EV%lmaxnu_tau(nu_i),nint(4*lAccuracyBoost))
          end if
          EV%lmaxnu_tau(nu_i)=min(EV%lmaxnu,EV%lmaxnu_tau(nu_i))

          EV%nu_ix(nu_i)=neq+1
          if (EV%MassiveNuApprox(nu_i)) then
             neq = neq+4
          else
             neq = neq+ EV%nq(nu_i)*(EV%lmaxnu_tau(nu_i)+1)
          endif
          maxeq = maxeq + nqmax*(EV%lmaxnu+1)
       end do
    else
       EV%has_nu_relativistic = .false.
    end if

    !Axions - adding on to the DE parameters - RL reading: this seems to be putting the axion EoMs to the last of the array. 
    EV%a_ix=neq+1 !RL reading: adding EV%a_ix to the end of the entire array, and since there are two equations for the axion EoM, neq is added by 2. This is repeatedly done after all the other equation components are reshuffled (neutrinos, etc.)
    EV%a_kg_ix = neq + 3 !RL adding the new KG equation solvers
    !neq=neq+2
    neq = neq + 4 !RL
    !maxeq=maxeq+2
    maxeq = maxeq + 4 !RL

    !         print*, 'a_ix =', EV%

    EV%ScalEqsToPropagate = neq
    if (present(max_num_eqns)) then
       max_num_eqns=maxeq
    end if

    !!write(*, *) 'end of SetupScalarArrayIndices, EV%ScalEqsToPropagate', EV%ScalEqsToPropagate 
  end subroutine SetupScalarArrayIndices

  subroutine CopyScalarVariableArray(y,yout, EV, EVout)
    use ModelData !RL for the integration by parts boundary condition
    type(EvolutionVars) EV, EVout
    real(dl), intent(in) :: y(EV%nvar)
    real(dl), intent(out) :: yout(EVout%nvar)
    real(dl) :: yprime(EV%nvar), yprimeout(EV%nvar) !RL for EFA
    integer lmax,i, nq
    integer nnueq,nu_i, ix_off, ix_off2, ind, ind2
    real(dl) q, pert_scale
    real(dl) v1_bg, v2_bg, drhoax_kg, grhoax_kg, a, a2, k, k2!, tau !RL inserted for the need of switching equations
    real(dl) tU, tV, tW, tdvarphi_c, tdvarphi_cp, tdvarphi_s, tdvarphi_sp, tdrho_ef, u_ax_ef, weight, u_ax_efa, wcorr_coeff, w_ax!RL for EFA , clxax_efa; weight is the weight for a smooth transition in ktau for addressing the quasistatic equilibrium issue (and the sigma-eta boundary discontinuity issue)
    real(dl) dtauda !RL added to check with BG
    real(dl) tdP_ef_test, csquared_ax_test, kamnorm_test !RL 042524 - test delta P
    real(dl) pnu, rhonu, dHsqdmt_term_pert, A_coeff_pert!RL testing
    real(dl) tdvarphi_c_altest, tdvarphi_cp_altest, tdvarphi_s_altest, tdvarphi_sp_altest !RL testing
    real(dl) dgpi_out, dgpiout_out
    real(dl) sources_temp(CTransScal%NumSources) !RL for the integration by parts boundary condition (though not very useful)
    real(dl) metricdelta_test(2) !RL 102723
    !RL: locals for runtime "no-ETA" criterion (see comment block below)
    real(dl) :: aux_num, aux_den, w_ax_p1_inst, u_ax_kg_inst
    logical  :: use_inst

    yprime = 0
    yout=0
    yout(1:basic_num_eqns) = y(1:basic_num_eqns)
    !RL adding a, a2 and k for computing fluid variables at the switch
    a=y(1)
    a2=a*a
    k=EV%k_buf
    k2=EV%k2_buf

    ! New code: NOTE that there is only one more index here
    if (.not. is_cosmological_constant) then
       yout(EVout%w_ix)=y(EV%w_ix)
    end if

    if (.not. EV%no_phot_multpoles .and. .not. EVout%no_phot_multpoles) then
       if (.not. EV%oscillation_started .and. EVout%oscillation_started) then
       end if
       
       if (EV%TightCoupling .or. EVout%TightCoupling) then
          lmax=1
       else
          lmax = min(EV%lmaxg,EVout%lmaxg)
       end if
       yout(EVout%g_ix:EVout%g_ix+lmax)=y(EV%g_ix:EV%g_ix+lmax)
       if (.not. EV%TightCoupling .and. .not. EVout%TightCoupling) then
          lmax = min(EV%lmaxgpol,EVout%lmaxgpol)
          yout(EVout%polind+2:EVout%polind+lmax)=y(EV%polind+2:EV%polind+lmax)
       end if
    end if

    if (.not. EV%no_nu_multpoles .and. .not. EVout%no_nu_multpoles) then
       if (.not. EV%oscillation_started .and. EVout%oscillation_started) then
       end if
       if (EV%high_ktau_neutrino_approx .or. EVout%high_ktau_neutrino_approx) then
          lmax=2
       else
          lmax = min(EV%lmaxnr,EVout%lmaxnr)
       end if
       yout(EVout%r_ix:EVout%r_ix+lmax)=y(EV%r_ix:EV%r_ix+lmax)
    end if

    if (CP%Num_Nu_massive /= 0) then
       if (.not. EV%oscillation_started .and. EVout%oscillation_started) then
       end if
       do nu_i=1,CP%Nu_mass_eigenstates
          ix_off=EV%nu_ix(nu_i)
          ix_off2=EVOut%nu_ix(nu_i)
          if (EV%MassiveNuApprox(nu_i) .and. EVout%MassiveNuApprox(nu_i)) then
             nnueq=4
             yout(ix_off2:ix_off2+nnueq-1)=y(ix_off:ix_off+nnueq-1)
          else if (.not. EV%MassiveNuApprox(nu_i) .and. .not. EVout%MassiveNuApprox(nu_i)) then
             lmax=min(EV%lmaxnu_tau(nu_i),EVOut%lmaxnu_tau(nu_i))
             nq = min(EV%nq(nu_i), EVOut%nq(nu_i))
             do i=1,nq
                ind= ix_off + (i-1)*(EV%lmaxnu_tau(nu_i)+1)
                ind2=ix_off2+ (i-1)*(EVOut%lmaxnu_tau(nu_i)+1)
                yout(ind2:ind2+lmax) = y(ind:ind+lmax)
             end do
             do i=nq+1, EVOut%nq(nu_i)
                lmax = min(EVOut%lmaxnu_tau(nu_i), EV%lmaxnr)
                ind2=ix_off2+ (i-1)*(EVOut%lmaxnu_tau(nu_i)+1)
                !yout(ind2:ind2+lmax) = y(EV%r_ix:EV%r_ix+lmax) !RL commented this line out to implement MX's fix when EV%no_nu_multpoles is F and EV%high_ktau_neutrino_approx is T
                !------------------ RL copying MX's bugfix
                if ((.not. EV%high_ktau_neutrino_approx) .or. lmax<=2) then
                   yout(ind2:ind2+lmax) = y(EV%r_ix:EV%r_ix+lmax)
                else
                   yout(ind2:ind2+2) = y(EV%r_ix:EV%r_ix+2)
                   yout(ind2+3:ind2+lmax) = 0._dl
                end if
                !------------------

                !Add leading correction for the mass
                q=nu_q(i)
                pert_scale=(nu_masses(nu_i)/q)**2/2
                lmax = min(lmax,EV%lmaxnu_pert)
                yout(ind2:ind2+lmax) = yout(ind2:ind2+lmax) &
                     + y(EV%nu_pert_ix:EV%nu_pert_ix+lmax)*pert_scale
             end do
          end if
       end do

       if (EVOut%has_nu_relativistic .and. EV%has_nu_relativistic) then
          if (.not. EV%oscillation_started .and. EVout%oscillation_started) then
       end if
          lmax = min(EVOut%lmaxnu_pert, EV%lmaxnu_pert)
          yout(EVout%nu_pert_ix:EVout%nu_pert_ix+lmax)=  y(EV%nu_pert_ix:EV%nu_pert_ix+lmax)
       end if
    end if

        ! Axions 
    yout(EVout%a_ix)=y(EV%a_ix)
    yout(EVout%a_ix+1)=y(EV%a_ix+1)
    ! RL adding KG
    if (.not. EV%oscillation_started .and. EVout%oscillation_started) then !switch, where EV%oscillation_started is false but EVout%oscillation_started is true
       
       call spline_out(loga_table,phinorm_table,phinorm_table_ddlga,ntable,dlog10(a),v1_bg)
       call spline_out(loga_table,phidotnorm_table,phidotnorm_table_ddlga,ntable,dlog10(a),v2_bg)
       !Add 1+w for the momentum term
       !w_ax_p1 = (2.0d0*(v2_bg**2.0d0)/a2)/((v2_bg**2.0d0)/a2 + (CP%m_ovH0*v1_bg)**2.0d0)
       !Originally, dv1 (i.e. delta_v1) is y(EV%a_kg_ix), and dv2 (i.e. delta_v2) is y(EV%a_kg_ix+1)
       drhoax_kg = (v2_bg*y(EV%a_kg_ix+1)/a2 + (CP%m_ovH0**2.0d0)*v1_bg*y(EV%a_kg_ix)*EV%renorm_c)*2.0d0 !RL 050324
       grhoax_kg = (v2_bg)**2.0d0/a2+(CP%m_ovH0*v1_bg)**2.0d0
       !For EFA, we will need to know hLdot, i.e. 2*k*z here, so some extra variables will be computed again
       call derivs(EV,EV%ScalEqsToPropagate,CP%tau_osc,y,yprime)
             
       !RL adding w correction to the background 
       !wcorr_coeff = CP%ah_osc*CP%a_osc/(CP%m_ovH0*(CP%H0/100.0d0))
       wcorr_coeff = CP%ahosc_ETA*CP%a_osc/(CP%m_ovH0*(CP%H0/100.0d0)) !RL082924
       !!write(*, *) 'In CopyScalarVariableArrays, CP%wEFA_c', CP%wEFA_c
       w_ax = ((wcorr_coeff/a2)**2.0d0)*CP%wEFA_c !(Not worry about h_in or h_out for low mass axions yet since this construction only works in RD and while the effect is small in late times, the w_ax correction itself is incorrect already)
       
       !Compute the U coefficient for constructing the pert EFA variables
       tU = ((-CP%tvarphi_c + CP%tvarphi_sp)*2.0d0*yprime(3) &
            &- 2.0d0*k2*y(EV%a_kg_ix+1)/(a2*(CP%m_ovH0**2.0d0)*CP%H0_in_Mpc_inv))/(a*CP%m_ovH0*CP%H0_in_Mpc_inv)&
            &+ 6.0d0*CP%ah_osc*y(EV%a_kg_ix)*EV%renorm_c/(a*CP%m_ovH0*CP%H0/100.0d0)
       
       !V coefficient
       tV = (-(CP%tvarphi_s + CP%tvarphi_cp)*2.0d0*yprime(3) &
            &+ 2.0d0*k2*y(EV%a_kg_ix)*EV%renorm_c/(a*CP%m_ovH0*CP%H0_in_Mpc_inv))/(a*CP%m_ovH0*CP%H0_in_Mpc_inv) &
            &+ 6.0d0*CP%ah_osc*y(EV%a_kg_ix+1)/(a2*(CP%m_ovH0**2.0d0)*CP%H0/100.0d0)

       !W coefficient
       tW = CP%A_coeff**2.0d0 + 3.0d0*CP%A_coeff*CP%ah_osc/(a*CP%m_ovH0*CP%H0/100.0d0) &
            &+ 2.0d0*k2/((a*CP%m_ovH0*CP%H0_in_Mpc_inv)**2.0d0) + 4.0d0
       !RL testing alternative tW with A_alt-------------
!RL0430       tW = CP%A_coeff_alt**2.0d0 + 3.0d0*CP%A_coeff_alt*CP%ah_osc/(a*CP%m_ovH0*CP%H0/100.0d0) + 2.0d0*k2/((a*CP%m_ovH0*CP%H0_in_Mpc_inv)**2.0d0) + 4.0d0
       !------------------------
       !The pert boundary conditions
        tdvarphi_c = y(EV%a_kg_ix)*EV%renorm_c 
        tdvarphi_cp = (-2.0d0*tU - (CP%A_coeff + 3.0d0*CP%ah_osc/(a*CP%m_ovH0*CP%H0/100.0d0))*tV)/(2.0d0*tW)
        tdvarphi_s = y(EV%a_kg_ix + 1)/(a*CP%m_ovH0) - tdvarphi_cp
        tdvarphi_sp = (CP%A_coeff*tU - (2.0d0 + k2/((a*CP%m_ovH0*CP%H0_in_Mpc_inv)**2.0d0))*tV)/(2.0d0*tW)
       !RL using the alt test boundary conditions
!RL0430       tdvarphi_c = y(EV%a_kg_ix)
!RL0430       tdvarphi_cp = (-2.0d0*tU - (CP%A_coeff_alt + 3.0d0*CP%ah_osc/(a*CP%m_ovH0*CP%H0/100.0d0))*tV)/(2.0d0*tW)
!RL0430       tdvarphi_s = y(EV%a_kg_ix + 1)/(a*CP%m_ovH0) - tdvarphi_cp
!RL0430       tdvarphi_sp = (CP%A_coeff_alt*tU - (2.0d0 + k2/((a*CP%m_ovH0*CP%H0_in_Mpc_inv)**2.0d0))*tV)/(2.0d0*tW)

       !RL: runtime "no-ETA" criterion (paper-appendix-C breakdown in deep
       !sub-Jeans iso modes).  When
       !  (delta-phi_c'^2 + delta-phi_s'^2) > X^2 * (delta-phi_c^2 + delta-phi_s^2),
       !the (C8) bilinear projects to a runaway delta_a; use instantaneous
       !KG values for the post-switch initial conditions instead.
       !Anchors:  iso-broken ratio ~5e3,  m>>Heq iso-OK ~0.1,  adi-OK ~3e-3.
       !X=10 chosen as a global validity check on the EFA averaging: bypasses
       !the runaway when the perturbation timescale is much shorter than 1/m.
       !Adi never triggers; iso modes in the deep sub-Jeans regime do, where
       !the iso TF is small anyway.  Edit the literal 100.0_dl below to retune.
       use_inst = .false.
       aux_num = tdvarphi_cp**2 + tdvarphi_sp**2
       aux_den = tdvarphi_c**2  + tdvarphi_s**2
       if (aux_den > 0._dl) then
          use_inst = (aux_num > aux_den * 100.0_dl)  ! X^2 = 100.0  (X = 10.0)
       end if

       !Normalized deltarho_ef
        tdrho_ef = (CP%m_ovH0**2.0d0)*(CP%tvarphi_s*tdvarphi_cp - CP%tvarphi_c*tdvarphi_sp &
             &+ CP%tvarphi_cp*tdvarphi_cp + CP%tvarphi_sp*tdvarphi_sp &
             &+ tdvarphi_s*(2.0d0*CP%tvarphi_s + CP%tvarphi_cp) &
             &+ tdvarphi_c*(2.0d0*CP%tvarphi_c - CP%tvarphi_sp))
       !Normalized deltaP_ef test, RL 042524 - note the factor of 2 due to the normalization of the variables. See Eqn. (46) of candidacy notes
       tdP_ef_test = tdrho_ef - 2._dl*(CP%m_ovH0**2.0d0)*(tdvarphi_s*CP%tvarphi_s + tdvarphi_c*CP%tvarphi_c)

       kamnorm_test = k2/(a2*((CP%m_ovH0*CP%H0_in_Mpc_inv)**2.0_dl))
      if (kamnorm_test .lt. 1.e-14_dl) then
          !RL dealing with machine precision issue - Taylor expand to leading order
          csquared_ax_test = kamnorm_test/4.0_dl + 5.0_dl*((1/(a*dtauda(a)))**2.0_dl)/(4.0_dl*(k2/kamnorm_test))
          !!write(*, *) 'machine precision, kamnorm, csquared_ax_test', kamnorm_test, csquared_ax_test
    else
       csquared_ax_test = (sqrt(1.0_dl + kamnorm_test) - 1.0_dl)**2.0_dl/(kamnorm_test) &
            &+ 5.0_dl*((1/(a*dtauda(a)))**2.0_dl)/(4.0_dl*(k2/kamnorm_test))       
       !!csquared_ax_test = (sqrt(1.0_dl + kamnorm_test) - 1.0_dl)**2.0_dl/(kamnorm_test) + 1.1_dl*((1/(a*dtauda(a)))**2.0_dl)/((k2/kamnorm_test))
      end if
       
       !Normalized u_ax_ef = (1+w)thetaax_ef/k. This is a bit long so declared a variable name for it
       !First just compute the RHS of the equation (without background rho+P)
      u_ax_ef = k*CP%m_ovH0*(tdvarphi_c*(CP%tvarphi_s + CP%tvarphi_cp) &
           &+ tdvarphi_s*(-CP%tvarphi_c + CP%tvarphi_sp))/(a*CP%H0_in_Mpc_inv)
       
       !Then take into account the background omaxh2_ef
       u_ax_ef = u_ax_ef/(CP%rhorefp_ovh2*(CP%H0**2.0d0/1.0d4))
       !RL 07/07/2023: add u_ax_efa taking into account the metric
  
       weight = (k/(CP%ahosc_ETA*CP%H0_in_Mpc_inv/(CP%H0/100.0d0)))**2._dl/&
            &(3._dl + (k/(CP%ahosc_ETA*CP%H0_in_Mpc_inv/(CP%H0/100.0d0)))**2._dl) !RL 082024

       
       u_ax_efa = u_ax_ef*(1._dl + w_ax)/(1._dl + (CP%Prefp/(CP%rhorefp_ovh2*(CP%H0**2.0d0/1.0d4)))) !- k*(CP%tau_osc**2._dl)*(CP%Prefp/(CP%rhorefp_hsq*(CP%H0**2.0d0/1.0d4)) - w_ax)*yprime(3)/((k*CP%tau_osc)**2._dl + (sup_C)**2._dl)
       
!       
       !Now the LHS of the two EoMs are assigned clxax_kg and u_ax_kg
       !RL: change for efa
       if (use_inst) then
          ! Instantaneous KG initial conditions.  No (u_ax_ef - u_ax_efa)/k
          ! correction here -- that belongs to the ETA reconstruction.  The
          ! metric_delta block below still runs and captures the rho_a
          ! discontinuity through (yprime - yprimeout) automatically.
          !
          ! RL 2026-05-21 Option B (advisor PS reinterpretation): enforce
          ! continuity in (delta_rho_a, rho_a*u_a) at the switch by writing
          !   y(a_kg_ix)   = delta_rho_KG  / rho_EFA
          !   y(a_kg_ix+1) = rho_KG*u_KG   / rho_EFA
          ! so the EFA evolution "sees"
          !   rho_EFA * y(a_kg_ix)   = delta_rho_KG
          !   rho_EFA * y(a_kg_ix+1) = rho_KG * u_KG.
          ! Previous form  yout = drhoax_kg/grhoax_kg  and  u_ax_kg_inst
          ! left an O(H/m) rho_KG/rho_EFA mismatch in the imprinted IC,
          ! which the stiff post-switch evolution amplified beyond O(H/m).
          w_ax_p1_inst = (2.0_dl * v2_bg**2 / a2) / &
                         (v2_bg**2 / a2 + (CP%m_ovH0 * v1_bg)**2)
          u_ax_kg_inst = w_ax_p1_inst * k * y(EV%a_kg_ix) * EV%renorm_c &
                         / (CP%H0_in_Mpc_inv * v2_bg)
          yout(EVout%a_kg_ix)   = drhoax_kg &
                                  / (CP%rhorefp_ovh2 * (CP%H0**2.0d0/1.0d4))
          yout(EVout%a_kg_ix+1) = u_ax_kg_inst * grhoax_kg &
                                  / (CP%rhorefp_ovh2 * (CP%H0**2.0d0/1.0d4))
       else
          yout(EVout%a_kg_ix) = tdrho_ef/(CP%rhorefp_ovh2*(CP%H0**2.0d0/1.0d4)) !deltaax_ef
          !RL 103023: deltaax_efa needs to change in order to keep v_ef = v_efa and sigma_ef = sigma_efa
          yout(EVout%a_kg_ix) = yout(EVout%a_kg_ix) + (3._dl*CP%ah_osc*CP%H0_in_Mpc_inv/(CP%H0/100.0d0))*(u_ax_ef-u_ax_efa)/k

          yout(EVout%a_kg_ix+1) = u_ax_efa!u_ax_ef !uax_ef, taken to be continuous, i.e. u_ax_efa
       end if
       
       !Now we have finished switching to EFA for EVout, call derivs again to get variables on the EFA side
       
       EVout%output_done = .false.
       call derivs(EVout,EVout%ScalEqsToPropagate,CP%tau_osc,yout,yprimeout)
       
       !We have y and yout variables. Construct the corresponding source boundary values
       if (CP%flat) then
          !Constructing adotoa*sigma out of the eta*k equation and k*z in derivs
          !Note that k and Kf(1) are the same between EV and EVout
          
          EVout%metric_delta(1)= ((-yprime(3)/k+3._dl*yprime(2)/k2)*(2._dl*(yprime(1)/a))/k &
               &- y(2)/k) - ((-yprimeout(3)/k+3._dl*yprimeout(2)/k2)*(2._dl*(yprimeout(1)/a))/k - yout(2)/k)
          !!EVout%metric_delta(1)= ((-yprime(3)/k+3._dl*yprime(2)/k2)*(2._dl*(yprime(1)/a)-CP%opac_tauosc)/k - y(2)/k + dgpi_out/k2) - ((-yprimeout(3)/k+3._dl*yprimeout(2)/k2)*(2._dl*(yprimeout(1)/a)-CP%opac_tauosc)/k - yout(2)/k + dgpiout_out/k2)
          EVout%metric_delta(2)= (-yprime(3)/k+3._dl*yprime(2)/k2)/k - &
               &(-yprimeout(3)/k+3._dl*yprimeout(2)/k2)/k
          
       else
          
          EVout%metric_delta(1)= ((-yprime(3)/k+3._dl*(yprime(2)-CP%curv*(-yprime(3)/k))/k2)*&
               &(2._dl*(yprime(1)/a))/(k*EV%Kf(1)) - y(2)/(k*EV%Kf(1))) -&
               & ((-yprimeout(3)/k+3._dl*(yprimeout(2) - CP%curv*(-yprimeout(3)/k))/k2)*&
               &(2._dl*(yprimeout(1)/a))/(k*EVout%Kf(1)) - yout(2)/(k*EV%Kf(1)))
          !!EVout%metric_delta(1)= ((-yprime(3)/k+3._dl*(yprime(2)-CP%curv*(-yprime(3)/k))/k2)*(2._dl*(yprime(1)/a)-CP%opac_tauosc)/(k*EV%Kf(1))+ dgpi_out/k2) - ((-yprimeout(3)/k+3._dl*(yprimeout(2) - CP%curv*(-yprimeout(3)/k))/k2)*(2._dl*(yprimeout(1)/a)-CP%opac_tauosc)/(k*EVout%Kf(1)) - yout(2)/(k*EV%Kf(1)) + dgpiout_out/k2)
          EVout%metric_delta(2)= (-yprime(3)/k+3._dl*(yprime(2)-CP%curv*(-yprime(3)/k))/k2)/(k*EV%Kf(1)) - &
               &(-yprimeout(3)/k+3._dl*(yprimeout(2) - CP%curv*(-yprimeout(3)/k))/k2)/(k*EV%Kf(1))
          !--------
       end if
       
       yout(2) = (EVout%metric_delta(2)*yprimeout(1)/a)*(k*EV%Kf(1))*weight + y(2) !yout(2) is etaTEFAnew*k
       
        EVout%metric_delta(1) = EVout%metric_delta(1) - weight*(yprimeout(1)/a)*EVout%metric_delta(2)
          
        EVout%metric_delta(2) = EVout%metric_delta(2)*(1._dl -weight)


    else 
       yout(EVout%a_kg_ix) = y(EV%a_kg_ix)
       yout(EVout%a_kg_ix+1) = y(EV%a_kg_ix+1)
    end if
    
    
  end subroutine CopyScalarVariableArray


  subroutine SetupTensorArrayIndices(EV, maxeq)
    type(EvolutionVars) EV
    integer nu_i, neq
    integer, optional, intent(out) :: maxeq
    neq=3
    EV%g_ix = neq-1 !EV%g_ix+2 is quadrupole
    if (.not. EV%TensTightCoupling) then
       EV%E_ix = EV%g_ix + (EV%lmaxt-1)
       EV%B_ix = EV%E_ix + (EV%lmaxpolt-1)
       neq = neq+ (EV%lmaxt-1)+(EV%lmaxpolt-1)*2
    end if
    if (present(maxeq)) then
       maxeq =3 + (EV%lmaxt-1)+(EV%lmaxpolt-1)*2
    end if
    EV%r_ix = neq -1
    if (DoTensorNeutrinos) then
       neq = neq + EV%lmaxnrt-1
       if (present(maxeq)) maxeq = maxeq+EV%lmaxnrt-1
       if (CP%Num_Nu_massive /= 0 ) then
          do nu_i=1, CP%nu_mass_eigenstates
             EV%EvolveTensorMassiveNu(nu_i) = nu_tau_nonrelativistic(nu_i) < 0.8*tau_maxvis*AccuracyBoost
             if (EV%EvolveTensorMassiveNu(nu_i)) then
                EV%nu_ix(nu_i)=neq-1
                neq = neq+ nqmax*(EV%lmaxnut-1)
                if (present(maxeq)) maxeq = maxeq + nqmax*(EV%lmaxnut-1)
             end if
          end do
       end if
    end if

    EV%TensEqsToPropagate = neq

  end  subroutine SetupTensorArrayIndices

  subroutine CopyTensorVariableArray(y,yout, EV, EVout)
    type(EvolutionVars) EV, EVOut
    real(dl), intent(in) :: y(EV%nvart)
    real(dl), intent(out) :: yout(EVout%nvart)
    integer lmaxpolt, lmaxt, nu_i, ind, ind2, i

    yout=0
    yout(1:3) = y(1:3)
    if (.not. EVOut%TensTightCoupling .and. .not.EV%TensTightCoupling) then
       lmaxt = min(EVOut%lmaxt,EV%lmaxt)
       yout(EVout%g_ix+2:EVout%g_ix+lmaxt)=y(EV%g_ix+2:EV%g_ix+lmaxt)
       lmaxpolt = min(EV%lmaxpolt, EVOut%lmaxpolt)
       yout(EVout%E_ix+2:EVout%E_ix+lmaxpolt)=y(EV%E_ix+2:EV%E_ix+lmaxpolt)
       yout(EVout%B_ix+2:EVout%B_ix+lmaxpolt)=y(EV%B_ix+2:EV%B_ix+lmaxpolt)
    end if
    if (DoTensorNeutrinos) then
       lmaxt=min(EV%lmaxnrt,EVOut%lmaxnrt)
       yout(EVout%r_ix+2:EVout%r_ix+lmaxt)=y(EV%r_ix+2:EV%r_ix+lmaxt)
       do nu_i =1, CP%nu_mass_eigenstates
          if (EV%EvolveTensorMassiveNu(nu_i)) then
             lmaxt=min(EV%lmaxnut,EVOut%lmaxnut)
             do i=1,nqmax
                ind= EV%nu_ix(nu_i) + (i-1)*(EV%lmaxnut-1)
                ind2=EVOut%nu_ix(nu_i)+ (i-1)*(EVOut%lmaxnut-1)
                yout(ind2+2:ind2+lmaxt) = y(ind+2:ind+lmaxt)
             end do
          end if
       end do
    end if

  end subroutine CopyTensorVariableArray

  subroutine GetNumEqns(EV)
    use MassiveNu
    !Set the number of equations in each hierarchy, and get total number of equations for this k
    type(EvolutionVars) EV
    real(dl) scal, max_nu_mass
    integer nu_i,q_rel,j
    
    if (CP%Num_Nu_massive == 0) then
       EV%lmaxnu=0
       max_nu_mass=0
    else
       max_nu_mass = maxval(nu_masses(1:CP%Nu_mass_eigenstates))
       do nu_i = 1, CP%Nu_mass_eigenstates
          !Start with momentum modes for which t_k ~ time at which mode becomes non-relativistic
          q_rel=0
          do j=1, nqmax
             !two different q's here EV%q ~k
             if (nu_q(j) > nu_masses(nu_i)*adotrad/EV%q) exit
             q_rel = q_rel + 1
          end do

          if (q_rel>= nqmax-2) then
             EV%nq(nu_i)=nqmax
          else
             EV%nq(nu_i)=q_rel
          end if
          !q_rel = nint(nu_masses(nu_i)*adotrad/EV%q) !two dffierent q's here EV%q ~k
          !EV%nq(nu_i)=max(0,min(nqmax0,q_rel)) !number of momentum modes to evolve intitially
          EV%nu_nonrelativistic(nu_i) = .false.
       end do

       EV%NuMethod = CP%MassiveNuMethod
       if (EV%NuMethod == Nu_Best) EV%NuMethod = Nu_Trunc
       !l_max for massive neutrinos
       if (CP%Transfer%high_precision) then
          EV%lmaxnu=nint(25*lAccuracyBoost)
       else
          !!EV%lmaxnu=max(3,nint(10*lAccuracyBoost)) !Default
          EV%lmaxnu=max(3,nint(10*lAccuracyBoost)) !RL 07/20/2023
          if (max_nu_mass>700) EV%lmaxnu=max(3,nint(15*lAccuracyBoost)) !Feb13 tweak
          !RL 07/13/2023-wanted to understand what transfer_high_precision does-------------------
          !!EV%lmaxnu=nint(25*lAccuracyBoost)
          !RL 07/13/2023-------------------
       endif
    end if

    if (CP%closed) then
       EV%FirstZerolForBeta = nint(EV%q*CP%r)
    else
       EV%FirstZerolForBeta=l0max !a large number
    end if

    EV%oscillation_started = .false. !RL
    EV%metric_delta(:) = 0._dl !RL 090323
    EV%high_ktau_neutrino_approx = .false.
    if (CP%WantScalars) then
       EV%TightCoupling=.true.
       EV%no_phot_multpoles =.false.
       EV%no_nu_multpoles =.false.
       EV%MassiveNuApprox=.false.

       if (HighAccuracyDefault .and. CP%AccuratePolarization) then
          EV%lmaxg  = max(nint(11*lAccuracyBoost),3)
       else
          EV%lmaxg  = max(nint(8*lAccuracyBoost),3)
       end if
       EV%lmaxnr = max(nint(14*lAccuracyBoost),3)
       if (max_nu_mass>700 .and. HighAccuracyDefault) EV%lmaxnr = max(nint(32*lAccuracyBoost),3) !Feb13 tweak

       EV%lmaxgpol = EV%lmaxg
       if (.not.CP%AccuratePolarization) EV%lmaxgpol=max(nint(4*lAccuracyBoost),3)

       !!write(*, *) 'Before large scale adjustment, lmaxnr?', EV%lmaxnr
       !if (EV%q < 0.05) then !Default
       if (EV%q < 0.05) then !RL 07/19/2023
          !Large scales need fewer equations
          scal  = 1          
          if (CP%AccuratePolarization) scal = 4  !But need more to get polarization right
          !write(*, *) 'In GetNumEqns, CP%AccuratePolarization', CP%AccuratePolarization
          EV%lmaxgpol=max(3,nint(min(8,nint(scal* 150* EV%q))*lAccuracyBoost))
          !!EV%lmaxnr=max(3,nint(min(7,nint(sqrt(scal)* 150 * EV%q))*lAccuracyBoost)) !Default
          EV%lmaxg=max(3,nint(min(8,nint(sqrt(scal) *300 * EV%q))*lAccuracyBoost))
          
          !write(*, *) 'In GetNumEqns, scal, EV%lmaxnr', scal, EV%lmaxnr
          !!EV%lmaxnr=EV%lmaxg !RL 07/13/2023
          !!EV%lmaxnr=max(3,nint(min(14,nint(sqrt(scal)* 280 * EV%q))*lAccuracyBoost)) !RL from WH smoother lmaxnr
          EV%lmaxnr=max(3,nint(min(8,nint(sqrt(scal)* 450 * EV%q))*lAccuracyBoost)) !RL modifying WH smoother lmaxnr
          !!EV%lmaxnr = 3 !RL hacking to test 07/18/2023
          if (EV%lmaxnr < EV%lmaxnu) then
                ! Nov 2020 change following Pavel Motloch report (RL added from newest CAMB)
             EV%lmaxnr = EV%lmaxnu
                !EV%lmaxnu = min(EV%lmaxnu, EV%lmaxnr) ! may be better but have not tested and makes small result changes
          endif
          
          if (CP%AccurateReionization) then
             EV%lmaxg=EV%lmaxg*4
             EV%lmaxgpol=EV%lmaxgpol*2
          end if
       end if

       if (EV%TransferOnly) then
          EV%lmaxgpol = min(EV%lmaxgpol,nint(5*lAccuracyBoost))
          EV%lmaxg = min(EV%lmaxg,nint(6*lAccuracyBoost))
       end if

       if (CP%Transfer%high_precision) then
          if (HighAccuracyDefault) then
             EV%lmaxnr=max(nint(45*lAccuracyBoost),3)
          else
             EV%lmaxnr=max(nint(30*lAccuracyBoost),3)
          endif
          if (EV%q > 0.04 .and. EV%q < 0.5) then !baryon oscillation scales
             EV%lmaxg=max(EV%lmaxg,10)
          end if
       end if

       if (CP%closed) then
          EV%lmaxnu=min(EV%lmaxnu, EV%FirstZerolForBeta-1)
          EV%lmaxnr=min(EV%lmaxnr, EV%FirstZerolForBeta-1)
          EV%lmaxg=min(EV%lmaxg, EV%FirstZerolForBeta-1)
          EV%lmaxgpol=min(EV%lmaxgpol, EV%FirstZerolForBeta-1)
       end if

       EV%poltruncfac=real(EV%lmaxgpol,dl)/max(1,(EV%lmaxgpol-2))
       EV%MaxlNeeded=max(EV%lmaxg,EV%lmaxnr,EV%lmaxgpol,EV%lmaxnu)
       if (EV%MaxlNeeded > max_l_evolve) stop 'Need to increase max_l_evolve; if using openmp, may also need to increase stack size'
       !write(*, *) 'In GetNumEqns, calling SetupScalarArrayIndices once if WantScalars, EV%q_ix, EV%q', EV%q_ix, EV%q
       call SetupScalarArrayIndices(EV,EV%nvar)
       !write(*, *) 'In GetNumEqns, after SetupScalarArrayIndices, EV%q_ix, EV%q', EV%q_ix, EV%q
       if (CP%closed) EV%nvar=EV%nvar+1 !so can reference lmax+1 with zero coefficient
       EV%lmaxt=0
    else
       EV%nvar=0
    end if

    if (CP%WantTensors) then
       EV%TensTightCoupling = .true.
       EV%lmaxt=max(3,nint(8*lAccuracyBoost))
       EV%lmaxpolt = max(3,nint(4*lAccuracyBoost))
       ! if (EV%q < 1e-3) EV%lmaxpolt=EV%lmaxpolt+1
       if (DoTensorNeutrinos) then
          EV%lmaxnrt=nint(6*lAccuracyBoost)
          EV%lmaxnut=EV%lmaxnrt
       else
          EV%lmaxnut=0
          EV%lmaxnrt=0
       end if
       if (CP%closed) then
          EV%lmaxt=min(EV%FirstZerolForBeta-1,EV%lmaxt)
          EV%lmaxpolt=min(EV%FirstZerolForBeta-1,EV%lmaxpolt)
          EV%lmaxnrt=min(EV%FirstZerolForBeta-1,EV%lmaxnrt)
          EV%lmaxnut=min(EV%FirstZerolForBeta-1,EV%lmaxnut)
       end if
       EV%MaxlNeededt=max(EV%lmaxpolt,EV%lmaxt, EV%lmaxnrt, EV%lmaxnut)
       if (EV%MaxlNeededt > max_l_evolve) stop 'Need to increase max_l_evolve'
       call SetupTensorArrayIndices(EV, EV%nvart)
    else
       EV%nvart=0
    end if


    if (CP%WantVectors) then
       EV%lmaxv=max(10,nint(8*lAccuracyBoost))
       EV%lmaxpolv = max(5,nint(5*lAccuracyBoost))

       EV%nvarv=(EV%lmaxv)+(EV%lmaxpolv-1)*2+3

       EV%lmaxnrv=nint(30*lAccuracyBoost)

       EV%nvarv=EV%nvarv+EV%lmaxnrv
       if (CP%Num_Nu_massive /= 0 ) then
          stop 'massive neutrinos not supported for vector modes'
       end if
    else
       EV%nvarv=0
    end if

    !write(*, *) 'End of GetNumEqns, k, EV%lmaxnr', EV%q, EV%lmaxnr !RL
  end subroutine GetNumEqns

  !cccccccccccccccccccccccccccccccccc
  subroutine SwitchToMassiveNuApprox(EV,y, nu_i)
    !When the neutrinos are no longer highly relativistic we use a truncated
    !energy-integrated hierarchy going up to third order in velocity dispersion
    type(EvolutionVars) EV, EVout
    integer, intent(in) :: nu_i

    real(dl) a,a2,pnu,clxnu,dpnu,pinu,rhonu
    real(dl) qnu
    real(dl) y(EV%nvar), yout(EV%nvar)

    a=y(1)
    a2=a*a
    EVout=EV
    EVout%MassiveNuApprox(nu_i)=.true.
    call SetupScalarArrayIndices(EVout)
    call CopyScalarVariableArray(y,yout, EV, EVout)

    !Get density and pressure as ratio to massles by interpolation from table
    call Nu_background(a*nu_masses(nu_i),rhonu,pnu)

    !Integrate over q
    call Nu_Integrate_L012(EV, y, a, nu_i, clxnu,qnu,dpnu,pinu)
    !clxnu_here  = rhonu*clxnu, qnu_here = qnu*rhonu
    dpnu=dpnu/rhonu
    qnu=qnu/rhonu
    clxnu = clxnu/rhonu
    pinu=pinu/rhonu

    yout(EVout%nu_ix(nu_i))=clxnu
    yout(EVout%nu_ix(nu_i)+1)=dpnu
    yout(EVout%nu_ix(nu_i)+2)=qnu
    yout(EVout%nu_ix(nu_i)+3)=pinu

    call Nu_Intvsq(EV,y, a, nu_i, EVout%G11(nu_i),EVout%G30(nu_i))
    !Analytic solution for higher moments, proportional to a^{-3}
    EVout%G11(nu_i)=EVout%G11(nu_i)*a2*a/rhonu
    EVout%G30(nu_i)=EVout%G30(nu_i)*a2*a/rhonu

    EV=EVout
    y=yout

  end subroutine SwitchToMassiveNuApprox

  subroutine MassiveNuVarsOut(EV,y,yprime,a,grho,gpres,dgrho,dgq,dgpi, gdpi_diff,pidot_sum)
    implicit none
    type(EvolutionVars) EV
    real(dl) :: y(EV%nvar), yprime(EV%nvar),a, grho,gpres,dgrho,dgq,dgpi, gdpi_diff,pidot_sum 
    !grho = a^2 kappa rho
    !gpres = a^2 kappa p
    !dgrho = a^2 kappa \delta\rho
    !dgp =  a^2 kappa \delta p
    !dgq = a^2 kappa q (heat flux)
    !dgpi = a^2 kappa pi (anisotropic stress)
    !dgpi_diff = a^2 kappa (3*p -rho)*pi

    integer nu_i
    real(dl) pinudot,grhormass_t, rhonu, pnu,  rhonudot
    real(dl) adotoa, grhonu_t,gpnu_t
    real(dl) clxnu, qnu, pinu, dpnu
    real(dl) dtauda

    do nu_i = 1, CP%Nu_mass_eigenstates
       grhormass_t=grhormass(nu_i)/a**2

       !Get density and pressure as ratio to massless by interpolation from table
       call Nu_background(a*nu_masses(nu_i),rhonu,pnu)

       if (EV%MassiveNuApprox(nu_i)) then
          clxnu=y(EV%nu_ix(nu_i))
          dpnu=y(EV%nu_ix(nu_i)+1) !RL added back for testing only, by consulting the corresponding line in SwitchToMassiveNuApprox (Not used in the original CAMB since it was commented out)
          qnu=y(EV%nu_ix(nu_i)+2)
          pinu=y(EV%nu_ix(nu_i)+3)
          pinudot=yprime(EV%nu_ix(nu_i)+3)
       else
          !Integrate over q
          call Nu_Integrate_L012(EV, y, a, nu_i, clxnu,qnu,dpnu,pinu)
          !clxnu_here  = rhonu*clxnu, qnu_here = qnu*rhonu
          dpnu=dpnu/rhonu !RL commented back
          !write(*, *) 'Rayne what is rhonu???', rhonu
          !write(*, *) 'Rayne, normalizing dpnu, what is it?', dpnu
          qnu=qnu/rhonu
          clxnu = clxnu/rhonu
          pinu=pinu/rhonu
          adotoa = 1/(a*dtauda(a))
          rhonudot = Nu_drho(a*nu_masses(nu_i),adotoa,rhonu)

          call Nu_pinudot(EV,y, yprime, a,adotoa, nu_i,pinudot)
          pinudot=pinudot/rhonu - rhonudot/rhonu*pinu
       endif
       !write(*, *) 'Rayne, EV%MassiveNuApprox(nu_i)', EV%MassiveNuApprox(nu_i)
       !write(031523, '(36e52.42,\)'), a, dpnu, rhonu, clxnu !RL

       grhonu_t=grhormass_t*rhonu
       gpnu_t=grhormass_t*pnu

       grho = grho  + grhonu_t
       gpres= gpres + gpnu_t

       dgrho= dgrho + grhonu_t*clxnu
       !write(*, *) 'Rayne, dgp, grhormass_t*dpnu', dgp, grhormass_t*dpnu
       !dgp = dgp + grhormass_t*dpnu !RL removed for cleanup on June 14 25 !RL added this line: (dpnu/(rhonu*clxnu))(sound speed) * grhonu_t*clxnu = (dpnu/rhonu)*grhonu_t = dpnu*grhormass_t. dpnu and rhonu are both in units of the mean density of one eigenstate of massless neutrinos, hence it works to convert the conventions this way. (Also can do sanity check against the two lines above that compute grhonu_t and gpnu_t)
       dgq  = dgq   + grhonu_t*qnu
       dgpi = dgpi  + grhonu_t*pinu
       gdpi_diff = gdpi_diff + pinu*(3*gpnu_t-grhonu_t)
       pidot_sum = pidot_sum + grhonu_t*pinudot
    end do

  end subroutine MassiveNuVarsOut

  subroutine Nu_Integrate_L012(EV,y,a,nu_i,drhonu,fnu,dpnu,pinu)
    type(EvolutionVars) EV
    !  Compute the perturbations of density and energy flux
    !  of one eigenstate of massive neutrinos, in units of the mean
    !  density of one eigenstate of massless neutrinos, by integrating over
    !  momentum.
    integer, intent(in) :: nu_i
    real(dl), intent(in) :: a, y(EV%nvar)
    real(dl), intent(OUT) ::  drhonu,fnu
    real(dl), optional, intent(OUT) :: dpnu,pinu
    real(dl) tmp, am, aq,v, pert_scale
    integer iq, ind

    !  q is the comoving momentum in units of k_B*T_nu0/c.

    drhonu=0
    fnu=0
    if (present(dpnu)) then
       dpnu=0
       pinu=0
    end if
    am=a*nu_masses(nu_i)
    ind=EV%nu_ix(nu_i)
    do iq=1,EV%nq(nu_i)
       aq=am/nu_q(iq)
       v=1._dl/sqrt(1._dl+aq*aq)
       drhonu=drhonu+ nu_int_kernel(iq)* y(ind)/v
       fnu=fnu+nu_int_kernel(iq)* y(ind+1)
       if (present(dpnu)) then
          dpnu=dpnu+  nu_int_kernel(iq)* y(ind)*v
          pinu=pinu+ nu_int_kernel(iq)*y(ind+2)*v
       end if
       ind=ind+EV%lmaxnu_tau(nu_i)+1
    end do
    ind = EV%nu_pert_ix
    do iq=EV%nq(nu_i)+1,nqmax
       !Get the rest from perturbatively relativistic expansion
       aq=am/nu_q(iq)
       v=1._dl/sqrt(1._dl+aq*aq)
       pert_scale=(nu_masses(nu_i)/nu_q(iq))**2/2
       tmp = nu_int_kernel(iq)*(y(EV%r_ix)  + pert_scale*y(ind)  )
       drhonu=drhonu+ tmp/v
       fnu=fnu+nu_int_kernel(iq)*(y(EV%r_ix+1)+ pert_scale*y(ind+1))
       if (present(dpnu)) then
          dpnu=dpnu+ tmp*v
          pinu = pinu+ nu_int_kernel(iq)*(y(EV%r_ix+2)+ pert_scale*y(ind+2))*v
       end if
    end do

    if (present(dpnu)) then
       dpnu = dpnu/3
    end if

  end subroutine Nu_Integrate_L012

  subroutine Nu_pinudot(EV,y, ydot, a,adotoa, nu_i,pinudot)
    type(EvolutionVars) EV
    integer, intent(in) :: nu_i
    real(dl), intent(in) :: a,adotoa, y(EV%nvar), ydot(EV%nvar)

    !  Compute the time derivative of the mean density in massive neutrinos
    !  and the shear perturbation.
    real(dl) pinudot
    real(dl) aq,q,v,aqdot,vdot
    real(dl) psi2,psi2dot
    real(dl) am, pert_scale
    integer iq,ind

    !  q is the comoving momentum in units of k_B*T_nu0/c.
    pinudot=0._dl
    ind=EV%nu_ix(nu_i)+2
    am=a*nu_masses(nu_i)
    do iq=1,EV%nq(nu_i)
       q=nu_q(iq)
       aq=am/q
       aqdot=aq*adotoa
       v=1._dl/sqrt(1._dl+aq*aq)
       vdot=-aq*aqdot/(1._dl+aq*aq)**1.5d0
       pinudot=pinudot+nu_int_kernel(iq)*(ydot(ind)*v+y(ind)*vdot)
       ind=ind+EV%lmaxnu_tau(nu_i)+1
    end do
    ind = EV%nu_pert_ix+2
    do iq=EV%nq(nu_i)+1,nqmax
       q=nu_q(iq)
       aq=am/q
       aqdot=aq*adotoa
       pert_scale=(nu_masses(nu_i)/q)**2/2
       v=1._dl/sqrt(1._dl+aq*aq)
       vdot=-aq*aqdot/(1._dl+aq*aq)**1.5d0
       psi2dot=ydot(EV%r_ix+2)  + pert_scale*ydot(ind)
       psi2=y(EV%r_ix+2)  + pert_scale*y(ind)
       pinudot=pinudot+nu_int_kernel(iq)*(psi2dot*v+psi2*vdot)
    end do

  end subroutine Nu_pinudot

  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  function Nu_pi(EV,y, a, nu_i) result(pinu)
    type(EvolutionVars) EV
    integer, intent(in) :: nu_i
    real(dl), intent(in) :: a, y(EV%nvar)
    real(dl) :: am
    real(dl) pinu,q,aq,v
    integer iq, ind

    if (EV%nq(nu_i)/=nqmax) stop 'Nu_pi nq/=nqmax0'
    pinu=0
    ind=EV%nu_ix(nu_i)+2
    am=a*nu_masses(nu_i)
    do iq=1, EV%nq(nu_i)
       q=nu_q(iq)
       aq=am/q
       v=1._dl/sqrt(1._dl+aq*aq)
       pinu=pinu+nu_int_kernel(iq)*y(ind)*v
       ind =ind+EV%lmaxnut+1
    end do

  end function Nu_pi

  !cccccccccccccccccccccccccccccccccccccccccccccc
  subroutine Nu_Intvsq(EV,y, a, nu_i, G11,G30)
    type(EvolutionVars) EV
    integer, intent(in) :: nu_i
    real(dl), intent(in) :: a, y(EV%nvar)
    real(dl), intent(OUT) ::  G11,G30

    !  Compute the third order variables (in velocity dispersion)
    !by integrating over momentum.
    real(dl) aq,q,v, am
    integer iq, ind

    !  q is the comoving momentum in units of k_B*T_nu0/c.
    am=a*nu_masses(nu_i)
    ind=EV%nu_ix(nu_i)
    G11=0._dl
    G30=0._dl
    if (EV%nq(nu_i)/=nqmax) stop 'Nu_Intvsq nq/=nqmax0'
    do iq=1, EV%nq(nu_i)
       q=nu_q(iq)
       aq=am/q
       v=1._dl/sqrt(1._dl+aq*aq)
       G11=G11+nu_int_kernel(iq)*y(ind+1)*v**2
       if (EV%lmaxnu_tau(nu_i)>2) then
          G30=G30+nu_int_kernel(iq)*y(ind+3)*v**2
       end if
       ind = ind+EV%lmaxnu_tau(nu_i)+1
    end do

  end subroutine Nu_Intvsq


  subroutine MassiveNuVars(EV,y,a,grho,gpres,dgrho,dgq, wnu_arr)
    implicit none
    type(EvolutionVars) EV
    real(dl) :: y(EV%nvar), a, grho,gpres,dgrho,dgq
    real(dl), intent(out), optional :: wnu_arr(max_nu)
    !grho = a^2 kappa rho
    !gpres = a^2 kappa p
    !dgrho = a^2 kappa \delta\rho
    !dgp =  a^2 kappa \delta p
    !dgq = a^2 kappa q (heat flux)
    integer nu_i
    real(dl) grhormass_t, rhonu, qnu, clxnu, grhonu_t, gpnu_t, pnu

    do nu_i = 1, CP%Nu_mass_eigenstates
       grhormass_t=grhormass(nu_i)/a**2

       !Get density and pressure as ratio to massless by interpolation from table
       call Nu_background(a*nu_masses(nu_i),rhonu,pnu)

       if (EV%MassiveNuApprox(nu_i)) then
          clxnu=y(EV%nu_ix(nu_i))
          qnu=y(EV%nu_ix(nu_i)+2)
       else
          !Integrate over q
          call Nu_Integrate_L012(EV, y, a, nu_i, clxnu,qnu)
          !clxnu_here  = rhonu*clxnu, qnu_here = qnu*rhonu
          qnu=qnu/rhonu
          clxnu = clxnu/rhonu
       endif

       grhonu_t=grhormass_t*rhonu
       gpnu_t=grhormass_t*pnu

       grho = grho  + grhonu_t
       gpres= gpres + gpnu_t
       dgrho= dgrho + grhonu_t*clxnu
       dgq  = dgq   + grhonu_t*qnu

       if (present(wnu_arr)) then
          wnu_arr(nu_i) =pnu/rhonu
       end if
    end do

  end subroutine MassiveNuVars


  subroutine GrowthRate(EV,y,tau,k,a,growth,clxtot)

    ! DM: This subroutine computes the growth rate as a function of (k,a) for a model
    ! with an axion dark matter component.
    ! Modified from axion-cosmomc-2011 using the full scalar field case
    ! Note this gives the full growth rate with axions but DOES NOT INCLUDE MASSIVE NEUTRINOS
    ! Including these correctly requires outputting delta_total 
    ! and numerically taking derivatives.

    ! Could put this in derivs, so the integrator doesn't need to be called
    ! again every time you need the growth rate?

    implicit none
    type(EvolutionVars) EV
    real(dl) :: y(EV%nvar),yprime(EV%nvar),a,a2,k,k2,growth
    real(dl) :: deltadot,adotoa,clxtot,grho,tau
    real(dl) :: grhoc_dot,grhoc_t,clxc_dot,clxc,dgrhoc_dot
    real(dl) :: v1_bg, v2_bg, LHS_bg, dv1, dv2, drhoax_kg, grhoax_kg, clxax_kg, clxax_kg_dot, u_ax_kg, u_ax_kg_dot !RL adding the KG variables
    real(dl) :: grhoax_dot,clxax_dot,grhoax_t,clxax,dgrhoax_dot
    real(dl) :: w_ax, w_ax_p1, wcorr_coeff, csquared_ax,cad2 !RL added w_ax_p1, wcorr_coeff 
    real(dl) :: grhob_t,grhob_dot,clxb,clxb_dot,dgrhob_dot
    real(dl) :: grhor_t,grhog_t,grhov_t
    real(dl) :: grhomat,dgrhomat,grhomat_dot,dgrhomat_dot
    real(dl) dorp
    real(dl) gr
    integer i
    

    ! DM: note g(density_var_dot) = a^2(density_var_dot) NOT (a^2 density_var)_dot
    !!write(*, *) 'GrowthRate is called' !RL testing
    !!write(*, *) 'Rayne, in GrowthRate, tau, tauosc', tau, CP%tau_osc
    call derivs(EV,EV%ScalEqsToPropagate,tau,y,yprime)
    !!if (tau .eq. CP%tau_osc .and. .not. EV%oscillation_started) then
    !!   write(*, *) 'Rayne, in GrowthRate on the KG side, tau = tauosc, derivs is called, yprime(2), yprime(3)', yprime(2), yprime(3)
    !!end if
    a = y(1)
    a2 = a**2.
    k=EV%k_buf
    k2=EV%k2_buf

    ! Get LCDM variables
    grhob_t=grhob/a
    grhoc_t=grhoc/a
    grhor_t=grhornomass/a2
    grhog_t=grhog/a2
    grhov_t=grhov*a**(-1-3*w_lam)
    clxc = y(3)
    clxb = y(4)
    clxc_dot = yprime(3)
    clxb_dot = yprime(4)
    !Axion variables
    !update axion variables 9/19 to include interpolated densities and eq of state from Dan's w_evolve module




    clxax = y(EV%a_ix)
    clxax_dot = yprime(EV%a_ix)
    !RL commenting out 07312023
    !RL adding KG variables 
    if (.not. EV%oscillation_started) then
       dv1 = y(EV%a_kg_ix)
       dv2 = y(EV%a_kg_ix + 1)
       !Spline to get the phi and phidot at this point
       
       call spline_out(loga_table,phinorm_table,phinorm_table_ddlga,ntable,dlog10(a),v1_bg)
       call spline_out(loga_table,phidotnorm_table,phidotnorm_table_ddlga,ntable,dlog10(a),v2_bg)
       !RL getting 1+w for momentum
       w_ax_p1 = (2.0d0*(v2_bg**2.0d0)/a2)/((v2_bg**2.0d0)/a2 + (CP%m_ovH0*v1_bg)**2.0d0)
       !!w_ax = w_ax_p1 - 1.0_dl
       drhoax_kg = (v2_bg*dv2/a2 + (CP%m_ovH0**2.0d0)*v1_bg*dv1*EV%renorm_c)*2.0d0 !RL 050324
       !Get the grhoax from field variables too!
       grhoax_kg = (v2_bg)**2.0d0/a2+(CP%m_ovH0*v1_bg)**2.0d0
       !!if (v1_bg .eq. v1_bg .and. v2_bg .eq. v2_bg) then
          !CAMB normalization
       !!   dorp = grhom*grhoax_kg/(CP%H0**2.0d0/1.0d4)
       !!else
       !!   dorp=0.0d0
       !!endif
       !write(*, *) 'Rayne, compare whether the grhoax_kg matches what the GDM calculates', grhoax_kg/(1.0d1**(dble(gr)) * (CP%H0**2.0d0/1.0d4)) - 1.0d0 !RL: this one is tried and true
       clxax_kg = drhoax_kg/grhoax_kg !(original grhoax table)(1.0d1**(dble(gr)) * (CP%H0**2.0d0/1.0d4))
       !clxax_kg = drhoax_kg/(1.0d1**(dble(gr)) * (CP%H0**2.0d0/1.0d4))
       u_ax_kg = w_ax_p1*k*dv1*EV%renorm_c/(CP%H0_in_Mpc_inv*v2_bg)
    else
       !RL adding w correction to the background 
       !wcorr_coeff = CP%ah_osc*CP%a_osc/(CP%m_ovH0*(CP%H0/100.0d0))
       wcorr_coeff = CP%ahosc_ETA*CP%a_osc/(CP%m_ovH0*(CP%H0/100.0d0)) !RL082924
       !!write(*, *) 'In GrowthRate, CP%wEFA_c', CP%wEFA_c
       w_ax = ((wcorr_coeff/a2)**2.0d0)*CP%wEFA_c
       w_ax_p1 = 1.0_dl + w_ax
       dorp=grhom*CP%rhorefp_ovh2*((CP%a_osc/a)**3.0d0)*dexp((wcorr_coeff**2.0d0)*3.0d0&
            &*CP%wEFA_c*(1.0d0/(a2**2.0d0) - 1.0d0/(CP%a_osc**4.0d0))/4.0d0) !RL 110923
       !dorp=grhom*CP%rhorefp_hsq*((CP%a_osc/a)**3.0d0)
       !past tauosc, directly obtain clxax and u_ax
       clxax_kg = y(EV%a_kg_ix)
       u_ax_kg = y(EV%a_kg_ix+1)
    end if

    ! Compute the Hubble rate
    grhoax_t = dorp*a2
    grho = grhoax_t+grhoc_t+grhob_t+grhov_t+grhor_t+grhog_t
    adotoa = sqrt(grho/3.)

    ! Compute the time derivatives of densities and perturbations

    ! LCDM variables
    grhoc_dot = -3.*adotoa*grhoc_t
    grhob_dot = -3.*adotoa*grhob_t
    dgrhoc_dot = grhoc_t*clxc_dot+clxc*grhoc_dot
    dgrhob_dot = grhob_t*clxb_dot+clxb*grhob_dot	

    ! Axions
    !make change to 1+w_ax, is this correct? 9/22
    grhoax_dot = -3.*adotoa*(w_ax_p1)*grhoax_t !RL used the w_ax_p1 from kg
    dgrhoax_dot = grhoax_t*clxax_dot+clxax*grhoax_dot !RL 012323: No I still have to derive clxax_dot. GrowthRate doesn't seem to be called anywhere apart from in the transfer function and even so it doesn't seem to interfere with the transfer function so LET'S SAVE IT TO LATER
    !dgrhoax_dot = grhoax_t*clxax_kg_dot+clxax_kg*grhoax_dot !RL used pert from kg instead


    ! Compute and output delta

    grhomat = grhoc_t+grhob_t+grhoax_t
    dgrhomat = clxc*grhoc_t+clxb*grhob_t+clxax_kg*grhoax_t !RL used pert from kg instead
    clxtot = dgrhomat/grhomat 

    ! Compute deltadot

    grhomat_dot = grhoc_dot+grhob_dot+grhoax_dot
    dgrhomat_dot = dgrhoc_dot+dgrhob_dot+dgrhoax_dot

    deltadot = -clxtot*(grhomat_dot/grhomat)+(dgrhomat_dot/grhomat)

    ! Output the growth rate f = d log delta/ d log a		
    growth = deltadot/(adotoa*clxtot)
    !updated for correct axion homogeneous equations of motion
  end subroutine GrowthRate

  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine output(EV,y,j,tau,sources,dgpi_out)
    use ThermoData
    use lvalues
    use constants
    use ModelData
    implicit none
    integer j
    type(EvolutionVars) EV
    real(dl), target :: y(EV%nvar),yprime(EV%nvar)
    real(dl), dimension(:),pointer :: ypol,ypolprime

    real(dl) dgq,grhob_t,grhor_t,grhoc_t,grhog_t,grhov_t,sigma,polter
    real(dl) qgdot,pigdot,pirdot,vbdot,dgrho, dgp !RL added dgp for testing
    real(dl), optional, intent(out) :: dgpi_out !RL 091023 for switch boundary condition of source integration
    real(dl) a,a2,dz,z,clxc,clxb,vb,clxg,qg,pig,clxr,qr,pir

    real(dl) tau,x,divfac
    real(dl) dgpi_diff, pidot_sum
    real(dl), target :: pol(3),polprime(3)
    !dgpi_diff = sum (3*p_nu -rho_nu)*pi_nu

    real(dl) k,k2  ,adotoa, grho, gpres,etak,phi,dgpi
    real(dl)  diff_rhopi, octg, octgprime,vq
    real(dl) sources(CTransScal%NumSources)
    !        real(dl) t4,t92
    real(dl) ISW
    real(dl) w_eff
    real(dl) hdotoh,ppiedot

    real(dl) :: v1_bg, v2_bg, LHS_bg, dv1, dv2, drhoax_kg, grhoax_kg, clxax_kg, clxax_kg_dot, u_ax_kg, u_ax_kg_dot !RL adding the KG variables
    !real(dl) :: dv2_analytical, dv1_analytical, clxax_analytical !RL adding the analytical solution for testing purposes
    !real(dl) :: v1_bg_test, v2_bg_test, v1_bg_mepsilon, v1_bg_pepsilon, v2_bg_mepsilon, v2_bg_pepsilon !RL background check cont.d
    real(dl) :: epsilon_bgcheck !RL used to check the BG
    real(dl) grhoax_t,clxax,v_ax, w_ax, w_ax_p1, wcorr_coeff, gpres_ax, cad2! Axion variables !RL added w_ax_p1, wcorr_coeff
    real(dl) growth,clxtot ! DM: Growth Rate
    real(dl) dorp
    real(dl) gr, gr_delogged
    integer i
    real(dl) dv1_quasitest, hLddot_forquasitest, v2dot_forquasitest, dv2_quasitest, clxax_quasitest
    real(dl) dpax_kg, dgpaxa2_kg, dgpnumass_temp, dgrhonumass_temp !RL testing gr_minus2x, gr_minus2x_delogged, gr_minus, gr_minus_delogged, gr_plus, gr_plus_delogged, gr_plus2x, gr_plus2x_delogged, dlnaepsilon, dgrdlna, dgrdlna_minus, dgrdlna_plus, gP, gP_minus, gP_plus, dgPdlna, dgPlnadgrlna
    real(dl) delta_rho_sync_quasitest, delta_p_sync_quasitest, delta_rho_rest_quasitest, delta_p_rest_quasitest
    real(dl) kamnorm, csquared_ax!RL testing
    real(dl) opacity_output, dopacity_output, opacity_use, dopacity_use, cs2_output !RL added to call thermo to compute opacity test
    real(dl) H0_in_Mpc_inv !RL testing
    real(dl) dtauda
    external dtauda
    real(dl) wef_test, gpresaxef_test !RL 08152023
    real(dl) sigma_test, v_ax_test !RL
    !!real(dl) :: leadingterms(5) !RL 08142023
    H0_in_Mpc_inv = dble(CP%H0)/dble(c/1.0d3) !RL

    yprime = 0

    !EV%output_done = .false.
    call derivs(EV,EV%ScalEqsToPropagate,tau,y,yprime)

    if (EV%TightCoupling .or. EV%no_phot_multpoles) then
       pol=0
       polprime=0
       ypolprime => polprime
       ypol => pol
    else
       ypolprime => yprime(EV%polind+1:)
       ypol => y(EV%polind+1:)
    end if

    k=EV%k_buf
    k2=EV%k2_buf

    a   =y(1)
    a2  =a*a
    etak=y(2)
    clxc=y(3)
    clxb=y(4)
    vb  =y(5)
    vbdot =yprime(5)

    !Compute expansion rate from: grho 8*pi*rho*a**2

    grhob_t=grhob/a
    grhoc_t=grhoc/a
    grhor_t=grhornomass/a2
    grhog_t=grhog/a2

    !Axion variables
    !update axion variables 9/19 to include interpolated densities and eq of state from Dan's w_evolve module


    !RL adding KG variables
    if (.not. EV%oscillation_started) then
       dv1 = y(EV%a_kg_ix)
       !write(*, *) 'Check that the dv2s are consistent:', yprime(EV%a_kg_ix)/y(EV%a_kg_ix + 1) - 1.0d0
       dv2 = y(EV%a_kg_ix+1)
       !spline to get the phi and phidot at this point
       !!if (EV%q_ix .eq. 68) then
         !write(*, *) '68thk in output, y(EV%a_kg_ix)*EV%renorm_c, y(EV%a_kg_ix+1)*EV%renorm_c, y(EV%a_kg_ix), y(EV%a_kg_ix+1):', k, y(EV%a_kg_ix)*EV%renorm_c, y(EV%a_kg_ix+1)*EV%renorm_c, y(EV%a_kg_ix), y(EV%a_kg_ix+1)
         !! write(042824, '(36E52.42E3)') k, y(EV%a_kg_ix)*EV%renorm_c, y(EV%a_kg_ix+1)*EV%renorm_c, CP%tvarphi_c, CP%tvarphi_s, CP%tvarphi_cp, CP%tvarphi_sp, tdvarphi_c, tdvarphi_s, tdvarphi_cp, tdvarphi_sp, tdrho_ef, tdP_ef_test, csquared_ax_test
       !!end if
       
       call spline_out(loga_table,phinorm_table,phinorm_table_ddlga,ntable,dlog10(a),v1_bg)
       call spline_out(loga_table,phidotnorm_table,phidotnorm_table_ddlga,ntable,dlog10(a),v2_bg)
       !obtain 1+w for momentum 
       w_ax_p1 = (2.0d0*(v2_bg**2.0d0)/a2)/((v2_bg**2.0d0)/a2 + (CP%m_ovH0*v1_bg)**2.0d0)
       w_ax = w_ax_p1 - 1.0_dl
       !cad2 is not assigned here since adotoa is not assigned yet (and not needed when we have KG)
       !drhoax and clxax calculated from KG
       drhoax_kg = (v2_bg*dv2/a2 + (CP%m_ovH0**2.0d0)*v1_bg*dv1*EV%renorm_c)*2.0d0!RL 050324
       !grhoax from background KG
       grhoax_kg = ((v2_bg)**2.0d0/a2+(CP%m_ovH0*v1_bg)**2.0d0)
       !renormalize to CAMB convension
       if (v1_bg .eq. v1_bg .and. v2_bg .eq. v2_bg) then !RL inherited DG flag
          dorp = grhom*grhoax_kg/(CP%H0**2.0d0/1.0d4)
       else
          dorp=0.0d0
       endif
       clxax_kg = drhoax_kg/grhoax_kg 
       u_ax_kg = w_ax_p1*k*dv1*EV%renorm_c/(CP%H0_in_Mpc_inv*v2_bg)
       v_ax_test = k*dv1*EV%renorm_c/(CP%H0_in_Mpc_inv*v2_bg)
       !dpax_kg = (v2_bg*dv2/a2 - (CP%m_ovH0**2.0d0)*v1_bg*dv1)*2.0d0 !RL for quasistatic test
       !write(*, *) 'Rayne, dpax_kg, drhoax_kg, their ratio', dpax_kg, drhoax_kg, dpax_kg/drhoax_kg
       !dgpaxa2_kg = grhom*dpax_kg*a2/(CP%H0**2.0d0/1.0d4) !RL switching to CAMB normalization
       !RL testing
       !!call spline_out(loga_table,CP%wef_table_test,CP%wef_table_test_buff,ntable,dlog10(a),wef_test)
    else
       !past tauosc, directly obtain clxax and u_ax
       !write(*, *) 'Rayne, after the switch, a, aosc, what is this mysterious value of v1?', a, CP%a_osc, v1_bg
       v1_bg = 0.0_dl
       v2_bg = 0.0_dl
       drhoax_kg = 0.0_dl
       grhoax_kg = 0.0_dl
       !RL adding w correction to the background
       !wcorr_coeff = CP%ah_osc*CP%a_osc/(CP%m_ovH0*(CP%H0/100.0d0))
       wcorr_coeff = CP%ahosc_ETA*CP%a_osc/(CP%m_ovH0*(CP%H0/100.0d0)) !RL082924
       w_ax = ((wcorr_coeff/a2)**2.0d0)*CP%wEFA_c
       w_ax_p1 = 1.0_dl + w_ax
       !I can have cad2 here because I need it and it only depends on the H at aosc per my selection of approximation
       !!cad2=w_ax*((w_ax + 7.0_dl/3.0_dl)/w_ax_p1)
       !write(*, *) 'Rayne, cad2 from output', cad2
       !cad2 = 0._dl
       dorp=grhom*CP%rhorefp_ovh2*((CP%a_osc/a)**3.0d0)*dexp((wcorr_coeff**2.0d0)*3.0d0*&
            &CP%wEFA_c*(1.0d0/(a2**2.0d0) - 1.0d0/(CP%a_osc**4.0d0))/4.0d0) !RL 110923
       !!write(*, *) 'In output, CP%wEFA_c', CP%wEFA_c
       !dorp=grhom*CP%rhorefp_hsq*((CP%a_osc/a)**3.0d0)
       clxax_kg = y(EV%a_kg_ix)
       u_ax_kg = y(EV%a_kg_ix+1)
       v_ax_test = u_ax_kg/w_ax_p1
    end if
    !!write(*, *) 'a, tau, w_ax, CP%H0_eV*10, ((wcorr_coeff/a2)**2.0d0), CP%wEFA_c', a, tau, w_ax, CP%H0_eV*10, ((wcorr_coeff/a2)**2.0d0), CP%wEFA_c


    grhoax_t=dorp*a2

    !Axion pert variables from DG
    clxax=y(EV%a_ix)
    v_ax=y(EV%a_ix+1)

    !if (.not. EV%oscillation_started) then !RL
    !   call spline_out(loga_table,CP%wax_table,CP%wax_table_buff,ntable,dlog10(a),w_ax)
    !else
    !   w_ax=0.0d0
    !endif
    !write(*,*) 'Rayne, output w_ax', w_ax
    gpres_ax=w_ax*grhoax_t
    !RL 08152023
    !!gpresaxef_test=wef_test*grhoax_t
    !write(*, *) 'Rayne, w_ax, wef_test', w_ax, wef_test


    grho=grhob_t+grhoc_t+grhor_t+grhog_t+grhoax_t ! RH removed the lambda term for later
    gpres=(grhog_t+grhor_t)/3+gpres_ax ! RH removed the lambda pressure term for later

    !  8*pi*a*a*SUM[rho_i*clx_i] add radiation later

    dgrho=grhob_t*clxb+grhoc_t*clxc+grhoax_t*clxax_kg !RL replaced with axion pert from kg

    !RL added dgp: 8*pi*a*a*SUM[p_i] only axions for now, baryons have negligible pressure, cdm has no pressure------------------
    !!dgpaxa2_kg = grhoax_t*clxax_kg*dpax_kg/drhoax_kg
    !!dgp = dgpaxa2_kg
    !--------------------------------

    !  8*pi*a*a*SUM[(rho_i+p_i)*v_i]
    dgq=grhob_t*vb+u_ax_kg*grhoax_t ! only baryons and axions for now

    if (is_cosmological_constant) then
       w_eff = -1_dl
       grhov_t=grhov*a2
       grho = grho+grhov_t ! add in the DE

    else
       !ppf
       w_eff=w_de(a)   !effective de
       grhov_t=grho_de(a)/a2
       grho = grho+grhov_t
       dgrho=dgrho+EV%dgrho_e_ppf !RL: DE has no pressure pert (Wayne taught me...but I'm not fully grasping it until I derive it) so I'm not adding dgp here.
       dgq=dgq+EV%dgq_e_ppf
    end if


    gpres=gpres+grhov_t*w_eff ! RH add in the DE now

    dgpi= 0
    dgpi_diff = 0
    pidot_sum = 0


    if (CP%Num_Nu_Massive /= 0) then
       call MassiveNuVarsOut(EV,y,yprime,a,grho,gpres,dgrho,dgq,dgpi, dgpi_diff,pidot_sum) 
    end if
    

    adotoa=sqrt((grho+grhok)/3)

    if (EV%no_nu_multpoles) then
       z=(0.5_dl*dgrho/k + etak)/adotoa
       dz= -adotoa*z - 0.5_dl*dgrho/k
       clxr=-4*dz/k
       qr=-4._dl/3*z
       pir=0
       pirdot=0
    else
       clxr=y(EV%r_ix)
       qr  =y(EV%r_ix+1)
       pir =y(EV%r_ix+2)
       pirdot=yprime(EV%r_ix+2)
    end if

    
    !_______________RL The following are for working mode (and for the debug mode that uses all timesteps the same as defined in TimeSteps%npoints)
    opacity_use = opac(j)
    dopacity_use = dopac(j)
    !_______________
    !RL also put comment flags at each place opac(j) and dopac(j) are replaced with the dummy variable

    if (EV%no_phot_multpoles) then
       z=(0.5_dl*dgrho/k + etak)/adotoa
       dz= -adotoa*z - 0.5_dl*dgrho/k
       clxg=-4*dz/k -4/k*opacity_use*(vb+z) !RL replaced opac(j) with opacity_use
       qg=-4._dl/3*z
       pig=0
       pigdot=0
       octg=0
       octgprime=0
       qgdot = -4*dz/3
    else
       if (EV%TightCoupling) then
          pig = EV%pig
          !pigdot=EV%pigdot
          if (second_order_tightcoupling) then
             octg = (3._dl/7._dl)*pig*(EV%k_buf/opacity_use) !RL replaced opac(j) with opacity_use
             ypol(2) = EV%pig/4 + EV%pigdot*(1._dl/opacity_use)*(-5._dl/8._dl) !RL replaced opac(j) with opacity_use
             ypol(3) = (3._dl/7._dl)*(EV%k_buf/opacity_use)*ypol(2) !RL replaced opac(j) with opacity_use
          else
             ypol(2) = EV%pig/4
             octg=0
          end if
          octgprime=0
       else
          pig =y(EV%g_ix+2)
          pigdot=yprime(EV%g_ix+2)
          octg=y(EV%g_ix+3)
          octgprime=yprime(EV%g_ix+3)
       end if
       clxg=y(EV%g_ix)
       qg  =y(EV%g_ix+1)
       qgdot =yprime(EV%g_ix+1)
    end if

    dgrho = dgrho + grhog_t*clxg+grhor_t*clxr
    !RL added dgp for the photons and massless neutrinos---------
    !!dgp = dgp + (grhog_t*clxg+grhor_t*clxr)/3._dl
    !------------------
    dgq   = dgq   + grhog_t*qg+grhor_t*qr
    dgpi  = dgpi  + grhor_t*pir + grhog_t*pig

    !  Get sigma (shear) and z from the constraints
    !  have to get z from eta for numerical stability
    z=(0.5_dl*dgrho/k + etak)/adotoa
    sigma=(z+1.5_dl*dgq/k2)/EV%Kf(1)

    if (is_cosmological_constant) then
       ppiedot=0
    else
       hdotoh=(-3._dl*grho-3._dl*gpres -2._dl*grhok)/6._dl/adotoa
       ppiedot=3._dl*EV%dgrho_e_ppf+EV%dgq_e_ppf*(12._dl/k*adotoa+k/adotoa-3._dl/k*(adotoa+hdotoh))+ &
            grhov_t*(1+w_eff)*k*z/adotoa -2._dl*k2*EV%Kf(1)*(yprime(EV%w_ix)/adotoa-2._dl*y(EV%w_ix))
       ppiedot=ppiedot*adotoa/EV%Kf(1)
    end if

    polter = 0.1_dl*pig+9._dl/15._dl*ypol(2)

    if (CP%flat) then
       x=k*(CP%tau0-tau)
       divfac=x*x
    else
       x=(CP%tau0-tau)/CP%r
       divfac=(CP%r*rofChi(x))**2*k2
    end if


    if (EV%TightCoupling) then
       if (second_order_tightcoupling) then
          pigdot = EV%pigdot
          ypolprime(2)= (pigdot/4._dl)*(1+(5._dl/2._dl)*(dopacity_use/opacity_use**2))!RL replaced opac(j) with opacity_use, dopac(j) with dopacity_use
       else
          pigdot = -dopacity_use/opacity_use*pig + 32._dl/45*k/opacity_use*(-2*adotoa*sigma  &
               +etak/EV%Kf(1)-  dgpi/k +vbdot ) !RL replaced opac(j) with opacity_use, dopac(j) with dopacity_use
          ypolprime(2)= pigdot/4
       end if
    end if

    pidot_sum =  pidot_sum + grhog_t*pigdot + grhor_t*pirdot
    diff_rhopi = pidot_sum - (4*dgpi+ dgpi_diff )*adotoa + ppiedot


    !Maple's fortran output - see scal_eqs.map
    !2phi' term (\phi' + \psi' in Newtonian gauge)
    ISW = (4.D0/3.D0*k*EV%Kf(1)*sigma+(-2.D0/3.D0*sigma-2.D0/3.D0*etak/adotoa)*k &
         -diff_rhopi/k**2-1.D0/adotoa*dgrho/3.D0+(3.D0*gpres+5.D0*grho)*sigma/k/3.D0 &
         -2.D0/k*adotoa/EV%Kf(1)*etak)*expmmu(j)
    !RL 08152023 test change
    !!ISW = (4.D0/3.D0*k*EV%Kf(1)*sigma+(-2.D0/3.D0*sigma-2.D0/3.D0*etak/adotoa)*k &
    !!     -diff_rhopi/k**2-1.D0/adotoa*dgrho/3.D0+(3.D0*(gpres-gpres_ax + gpresaxef_test)+5.D0*grho)*sigma/k/3.D0 &
    !!     -2.D0/k*adotoa/EV%Kf(1)*etak)*expmmu(j) !
    !e.g. to get only late-time ISW
    !  if (1/a-1 < 30) ISW=0

    !The rest, note y(9)->octg, yprime(9)->octgprime (octopoles)
    !RL replaced opac(j) with opacity_use, dopac(j) with dopacity_use
     !write(*,*) 'Rayne, before sources(1), vis(j)', vis(j)
    sources(1)= ISW +  ((-9.D0/160.D0*pig-27.D0/80.D0*ypol(2))/k**2*opacity_use+ &
         (11.D0/10.D0*sigma- 3.D0/8.D0*EV%Kf(2)*ypol(3)+vb-9.D0/80.D0*EV%Kf(2)*octg+3.D0/40.D0*qg)/k- &
         (-180.D0*ypolprime(2)-30.D0*pigdot)/k**2/160.D0)*dvis(j) + &
         (-(9.D0*pigdot+ 54.D0*ypolprime(2))/k**2*opacity_use/160.D0+pig/16.D0+clxg/4.D0+3.D0/8.D0*ypol(2) + &
         (-21.D0/5.D0*adotoa*sigma-3.D0/8.D0*EV%Kf(2)*ypolprime(3) + &
         vbdot+3.D0/40.D0*qgdot- 9.D0/80.D0*EV%Kf(2)*octgprime)/k + &
         (-9.D0/160.D0*dopacity_use*pig-21.D0/10.D0*dgpi-27.D0/80.D0*dopacity_use*ypol(2))/k**2)*vis(j) + &
         (3.D0/16.D0*ddvis(j)*pig+9.D0/8.D0*ddvis(j)*ypol(2))/k**2+21.D0/10.D0/k/EV%Kf(1)*vis(j)*etak

    ! Doppler term
    !   sources(1)=  (sigma+vb)/k*dvis(j)+((-2.D0*adotoa*sigma+vbdot)/k-1.D0/k**2*dgpi)*vis(j) &
    !         +1.D0/k/EV%Kf(1)*vis(j)*etak

    !Equivalent full result
    !    t4 = 1.D0/adotoa
    !    t92 = k**2
    !    sources(1) = (4.D0/3.D0*EV%Kf(1)*expmmu(j)*sigma+2.D0/3.D0*(-sigma-t4*etak)*expmmu(j))*k+ &
    !        (3.D0/8.D0*ypol(2)+pig/16.D0+clxg/4.D0)*vis(j)
    !    sources(1) = sources(1)-t4*expmmu(j)*dgrho/3.D0+((11.D0/10.D0*sigma- &
    !         3.D0/8.D0*EV%Kf(2)*ypol(3)+vb+ 3.D0/40.D0*qg-9.D0/80.D0*EV%Kf(2)*y(9))*dvis(j)+(5.D0/3.D0*grho+ &
    !        gpres)*sigma*expmmu(j)+(-2.D0*adotoa*etak*expmmu(j)+21.D0/10.D0*etak*vis(j))/ &
    !        EV%Kf(1)+(vbdot-3.D0/8.D0*EV%Kf(2)*ypolprime(3)+3.D0/40.D0*qgdot-21.D0/ &
    !        5.D0*sigma*adotoa-9.D0/80.D0*EV%Kf(2)*yprime(9))*vis(j))/k+(((-9.D0/160.D0*pigdot- &
    !        27.D0/80.D0*ypolprime(2))*opac(j)-21.D0/10.D0*dgpi -27.D0/80.D0*dopac(j)*ypol(2) &
    !        -9.D0/160.D0*dopac(j)*pig)*vis(j) - diff_rhopi*expmmu(j)+((-27.D0/80.D0*ypol(2)-9.D0/ &
    !        160.D0*pig)*opac(j)+3.D0/16.D0*pigdot+9.D0/8.D0*ypolprime(2))*dvis(j)+9.D0/ &
    !        8.D0*ddvis(j)*ypol(2)+3.D0/16.D0*ddvis(j)*pig)/t92


    if (x > 0._dl) then
       !E polarization source
       sources(2)=vis(j)*polter*(15._dl/8._dl)/divfac
       !factor of four because no 1/16 later
    else
       sources(2)=0
    end if

    if (CTransScal%NumSources > 2) then
       !Get lensing sources
       !Can modify this here if you want to get power spectra for other tracer
       if (tau>taurend .and. CP%tau0-tau > 0.1_dl) then
          !phi_lens = Phi - 1/2 kappa (a/k)^2 sum_i rho_i pi_i
          !Neglect pi contributions because not accurate at late time anyway
          phi = -(dgrho +3*dgq*adotoa/k)/(k2*EV%Kf(1)*2)
          ! - dgpi/k2/2

          sources(3) = -2*phi*f_K(tau-tau_maxvis)/(f_K(CP%tau0-tau_maxvis)*f_K(CP%tau0-tau))

          !         sources(3) = -2*phi*(tau-tau_maxvis)/((CP%tau0-tau_maxvis)*(CP%tau0-tau))
          !We include the lensing factor of two here
       else
          sources(3) = 0
       end if
    end if

    !    if (abs(k-0.5)<1e-1) then
    !      write(889, *) a, grho/const_rhocrit
    !  end if

    !if (abs(a/2.934062d-4 - 1._dl) < 1.d-4) then
    !   write(*, *) 'Rayne, estimating a, adotoa at aeq:', a, adotoa
    !end if

    !RL: test write out the variables
    !write(001, *) a, k, etak, clxc, clxb, vb, clxg, qg, pig, clxr, qr, pir, pigdot, pirdot, clxax

    !if (EV%oscillation_started == .true. .and. EV%oscillation_output_done == .false.) then
    !   write(*, *) 'Rayne what is a, a_osc, their fractional difference, CP%a_osc, m/H, m/H over dfac - 1, right at the switch?', a, a_osc, a/a_osc - 1._dl, CP%a_osc, CP%m_ovH0*CP%H0_in_Mpc_inv/(adotoa/a), (CP%m_ovH0*CP%H0_in_Mpc_inv/(adotoa/a))/10._dl - 1._dl
    !write(001, '(15e25.15,\)') k, a, a_osc, tau, CP%tau_osc, clxax, clxc 
    !   write(*, *) 'Oscillation started, here is the last y(EV%a_ix) and y(EV%a_ix+1):'
    !   write(*, '(36e52.42,\)') y(EV%a_ix), y(EV%a_ix+1)
    !   write(*, '(36e52.42,\)') k, a, a_osc, tau, CP%tau_osc, clxax, clxc
    !end if


    !RL testing the quasistatic evolution (have to be after z is defined)
    !dv1_quasitest = -z*v2_bg*(CP%H0_in_Mpc_inv)/k
    !hLddot_forquasitest = -2._dl*k*z*adotoa - dgrho - 3._dl*dgp
    !dv2/dtau (with normalization), from the background KG
    !v2dot_forquasitest = -2.0d0*adotoa*v2_bg - a2*(CP%m_ovH0**2.0d0)*CP%H0_in_Mpc_inv*v1_bg
    !dv2_quasitest = -hLddot_forquasitest*v2_bg/(2.0d0*k2) - z*v2dot_forquasitest/k
    !delta_rho_sync_quasitest = 2.0d0*(v2_bg*dv2_quasitest/a2 + (CP%m_ovH0**2.0d0)*v1_bg*dv1_quasitest)
    !delta_p_sync_quasitest = 2.0d0*(v2_bg*dv2_quasitest/a2 - (CP%m_ovH0**2.0d0)*v1_bg*dv1_quasitest)
    !clxax_quasitest = delta_rho_sync_quasitest/((v2_bg)**2.0d0/a2+(CP%m_ovH0*v1_bg)**2.0d0)
    !cad2 = 1.0d0 + 2.0d0 * (CP%m_ovH0**2.0d0) * a2 * v1_bg/(3 * (adotoa/CP%H0_in_Mpc_inv) * v2_bg)
    !delta_rho_rest_quasitest = delta_rho_sync_quasitest + 2.0d0*(3.0d0*adotoa*v2_bg*dv1_quasitest/(a2*CP%H0_in_Mpc_inv))
    !delta_p_rest_quasitest = delta_p_sync_quasitest + 2.0d0*(3.0d0*adotoa*cad2*v2_bg*dv1_quasitest/(a2*CP%H0_in_Mpc_inv))
    !mnorm=2.0d0*CP%m_ovH0*a*CP%H0*1.d3/(c) !for output cs2
    !write(030123, '(36e52.42,\)') k, a, CP%a_osc, tau, CP%tau_osc, adotoa, CP%m_ovH0*CP%H0_in_Mpc_inv/(adotoa/a), v1_bg, v2_bg, dorp, w_ax_p1, cad2, dv1_quasitest, dv1, dv2, clxax_kg, thetaax_kg, z, dgrho, clxc !LHS_bg ,  clxax_quasitest, z, dgrho, clxg
    !!write(032823, '(36e52.42,\)') CP%m_ovH0, k, a, CP%a_osc, tau, CP%tau_osc, adotoa, hLddot_forquasitest, z, v1_bg, v2_bg, v2dot_forquasitest, dv1, dv1_quasitest, dv2, dv2_quasitest, clxax_kg, clxax_quasitest, 2.0d0*v2_bg*dv2/a2, 2.0d0*(CP%m_ovH0**2.0d0)*v1_bg*dv1 !dgrho, dgp, clxax_kg, dgpaxa2_kg, dgpnumass_temp, grhoax_t*clxax_kg, grhoc_t*clxc, grhob_t*clxb, grhog_t*clxg, grhor_t*clxr, dgrhonumass_temp
    !write(033123, '(36e52.42,\)') CP%m_ovH0, k, k/(CP%m_ovH0*CP%H0_in_Mpc_inv), a, CP%a_osc, tau, CP%tau_osc, adotoa, CP%m_ovH0*CP%H0_in_Mpc_inv/(adotoa/a), v1_bg, v2_bg, dv1, dv1_quasitest, dv2, dv2_quasitest, clxax_kg, clxax_quasitest, delta_rho_sync_quasitest, delta_p_sync_quasitest, delta_rho_rest_quasitest, delta_p_rest_quasitest, delta_p_rest_quasitest/delta_rho_rest_quasitest, ((k/mnorm)**2.0d0)/(1.0d0+(k/mnorm)**2.0d0)
    kamnorm = k2/(a2*((CP%m_ovH0*CP%H0_in_Mpc_inv)**2.0_dl)) !RL testing
  if (kamnorm .lt. 1.e-13_dl) then
    !!write(*, *) 'Rayne, machine precision', kamnorm
    !RL dealing with machine precision issue - Taylor expand to leading order
   csquared_ax = kamnorm/4.0_dl + 5.0_dl*(adotoa**2.0_dl)/(4.0_dl*(k2/kamnorm))   
  else  
     csquared_ax = (sqrt(1.0_dl + kamnorm) - 1.0_dl)**2.0_dl/(kamnorm) + 5.0_dl*(adotoa**2.0_dl)/(4.0_dl*(k2/kamnorm))     
     !csquared_ax = (sqrt(1.0_dl + kamnorm) - 1.0_dl)**2.0_dl/(kamnorm) + 1.1_dl*(adotoa**2.0_dl)/((k2/kamnorm))
   end if
   !!write(081423, '(36e52.42,\)') CP%m_ovH0, k, k/(CP%m_ovH0*CP%H0_in_Mpc_inv), a, CP%a_osc, tau, CP%tau_osc, 2.0_dl*k*z, etak, clxax_kg, sources(1), sources(2), sources(3), ISW, pig, pigdot, clxg, ypol(2), ypolprime(2), ypol(3), ypolprime(3), opacity_use, dopacity_use, sigma, EV%Kf(1),  EV%Kf(2), vb, vbdot, octg, octgprime, qg, qgdot, dgpi, vis(j), dvis(j), ddvis(j), diff_rhopi, dgrho, gpres, grho, expmmu(j) !(4.D0/3.D0*k*EV%Kf(1)*sigma)*expmmu(j), (-2.D0/3.D0*sigma-2.D0/3.D0*etak/adotoa)*k*expmmu(j), -diff_rhopi*expmmu(j)/k**2-1.D0/adotoa*dgrho/3.D0, (3.D0*(gpres-gpres_ax + gpresaxef_test)+5.D0*grho)*sigma/k/3.D0, (-2.D0/k*adotoa/EV%Kf(1)*etak)*expmmu(j) 
    !write(*,'(36e52.42,\)') CP%H0_in_Mpc_inv
    !write(043023, '(36e52.42,\)') CP%m_ovH0, k, k/(CP%m_ovH0*CP%H0_in_Mpc_inv), a, CP%a_osc, tau, CP%tau_osc, adotoa, CP%m_ovH0*CP%H0_in_Mpc_inv/(adotoa/a), v1_bg, v2_bg, dorp*(CP%H0**2.0d0/1.0d4)/grhom, dv1, dv2, drhoax_kg, clxc, clxax_kg
    !write(*, *) 'k, clxc, clxax_kg in output', k, clxc, clxax_kg !RL
    !Testing the analytical solution
    !dv2_analytical = (CP%m_ovH0**2.0d0) * CP%H0_in_Mpc_inv * k * z * a2 * v1_bg/(35.0d0 * (adotoa**2.0d0))
    !write(*, *) 'Rayne, what is the first entry of CP%m_ovH0, CP%H0_in_Mpc_inv, k * z, a2, v1_bg, adotoa, dv2_analytical in the output analytical solution?', CP%m_ovH0, CP%H0_in_Mpc_inv, k * z, a2, v1_bg, adotoa, dv2_analytical
    !dv1_analytical = ((CP%m_ovH0* CP%H0_in_Mpc_inv)**2.0d0) * k * z * a2 * v1_bg/(210.0d0 * (adotoa**3.0d0))
    !clxax_analytical = (v2_bg*dv2_analytical/a2 + (CP%m_ovH0**2.0d0)*v1_bg*dv1_analytical)*2.0d0/grhoax_kg 
    !write(111001, '(36e52.42,\)') k, a, a_osc, tau, CP%tau_osc, adotoa, dv1, dv1_analytical, dv2, dv2_analytical, clxax, clxax_kg, clxax_analytical!, v1_bg, v2_bg,

   !!    write(073123, '(36e52.42,\)') a, CP%a_osc, tau, CP%tau_osc, adotoa/H0_in_Mpc_inv, dorp*(CP%H0**2.0d0/1.0d4)/grhom, clxax_kg
   
    !!write(*, *) 'Rayne, in output, j, expmmu(j)', j, expmmu(j)
    if (present(dgpi_out)) then !RL 091023
       dgpi_out = dgpi
    !!else  !RL 092123
       !write(092023, '(36e52.42)') a, CP%a_osc, tau, CP%tau_osc, adotoa/H0_in_Mpc_inv, dorp*(CP%H0**2.0d0/1.0d4)/grhom, clxax_kg
    end if
    
 
  end subroutine output


  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine outputt(EV,yt,n,j,tau,dt,dte,dtb)
    !calculate the tensor sources for open and closed case
    use ThermoData

    implicit none
    integer j,n
    type(EvolutionVars) :: EV
    real(dl), target :: yt(n), ytprime(n)
    real(dl) tau,dt,dte,dtb,x,polterdot,polterddot,prefac
    real(dl) pig, pigdot, octg, aux, polter, shear, adotoa,a
    real(dl) sinhxr,cothxor
    real(dl) k,k2
    real(dl), dimension(:),pointer :: E,Bprime,Eprime
    real(dl), target :: pol(3),polEprime(3), polBprime(3)
    real(dl) dtauda

    call derivst(EV,EV%nvart,tau,yt,ytprime)

    k2=EV%k2_buf
    k=EV%k_buf
    aux=EV%aux_buf
    shear = yt(3)

    x=(CP%tau0-tau)/CP%r

    !  And the electric part of the Weyl.
    if (.not. EV%TensTightCoupling) then
       !  Use the full expression for pigdt
       pig=yt(EV%g_ix+2)
       pigdot=ytprime(EV%g_ix+2)
       E => yt(EV%E_ix+1:)
       Eprime=> ytprime(EV%E_ix+1:)
       Bprime => ytprime(EV%B_ix+1:)
       octg=ytprime(EV%g_ix+3)
    else
       !  Use the tight-coupling approximation
       a =yt(1)
       adotoa = 1/(a*dtauda(a))
       pigdot=32._dl/45._dl*k/opac(j)*(2._dl*adotoa*shear+ytprime(3))
       pig = 32._dl/45._dl*k/opac(j)*shear
       pol=0
       polEprime=0
       polBprime=0
       E=>pol
       EPrime=>polEPrime
       BPrime=>polBPrime
       E(2)=pig/4._dl
       EPrime(2)=pigdot/4
       octg=0
    endif

    sinhxr=rofChi(x)*CP%r

    if (EV%q*sinhxr > 1.e-8_dl) then
       prefac=sqrt(EV%q2*CP%r*CP%r-CP%Ksign)
       cothxor=cosfunc(x)/sinhxr

       polter = 0.1_dl*pig + 9._dl/15._dl*E(2)
       polterdot=9._dl/15._dl*Eprime(2) + 0.1_dl*pigdot
       polterddot = 9._dl/15._dl*(-dopac(j)*(E(2)-polter)-opac(j)*(  &
            Eprime(2)-polterdot) + k*(2._dl/3._dl*Bprime(2)*aux - 5._dl/27._dl*Eprime(3)*EV%Kft(2))) &
            +0.1_dl*(k*(-octg*EV%Kft(2)/3._dl + 8._dl/15._dl*ytprime(3)) - &
            dopac(j)*(pig - polter) - opac(j)*(pigdot-polterdot))

       dt=(shear*expmmu(j) + (15._dl/8._dl)*polter*vis(j)/k)*CP%r/sinhxr**2/prefac

       dte=CP%r*15._dl/8._dl/k/prefac* &
            ((ddvis(j)*polter + 2._dl*dvis(j)*polterdot + vis(j)*polterddot)  &
            + 4._dl*cothxor*(dvis(j)*polter + vis(j)*polterdot) - &
            vis(j)*polter*(k2 -6*cothxor**2))

       dtb=15._dl/4._dl*EV%q*CP%r/k/prefac*(vis(j)*(2._dl*cothxor*polter + polterdot) + dvis(j)*polter)
    else
       dt=0._dl
       dte=0._dl
       dtb=0._dl
    end if

  end subroutine outputt

  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine outputv(EV,yv,n,j,tau,dt,dte,dtb)
    !calculate the vector sources
    use ThermoData

    implicit none
    integer j,n
    type(EvolutionVars) :: EV
    real(dl), target :: yv(n), yvprime(n)
    real(dl) tau,dt,dte,dtb,x,polterdot
    real(dl) vb,qg, pig, polter, sigma
    real(dl) k,k2
    real(dl), dimension(:),pointer :: E,Eprime

    call derivsv(EV,EV%nvarv,tau,yv,yvprime)

    k2=EV%k2_buf
    k=EV%k_buf
    sigma = yv(2)
    vb  = yv(3)
    qg  = yv(4)
    pig = yv(5)


    x=(CP%tau0-tau)*k

    if (x > 1.e-8_dl) then
       E => yv(EV%lmaxv+3:)
       Eprime=> yvprime(EV%lmaxv+3:)

       polter = 0.1_dl*pig + 9._dl/15._dl*E(2)
       polterdot=9._dl/15._dl*Eprime(2) + 0.1_dl*yvprime(5)

       if (yv(1) < 1e-3) then
          dt = 1
       else
          dt =0
       end if
       dt= (4*(vb+sigma)*vis(j) + 15._dl/2/k*( vis(j)*polterdot + dvis(j)*polter) &
            + 4*(expmmu(j)*yvprime(2)) )/x

       dte= 15._dl/2*2*polter/x**2*vis(j) + 15._dl/2/k*(dvis(j)*polter + vis(j)*polterdot)/x

       dtb= -15._dl/2*polter/x*vis(j)
    else
       dt=0
       dte=0
       dtb=0
    end if

  end subroutine outputv


  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine initial(EV, y, tau)
    !  Initial conditions.
    use constants
    use ThermoData
    implicit none

    type(EvolutionVars) EV
    real(dl) y(EV%nvar)
    real(dl) Rp15,tau,x,x2,x3,om,omtau, &
         Rc,Rb,Rv,Rg,grhonu,chi
    real(dl) k,k2
    real(dl) a,a2, dtauda, adotoa, dgrho, v1_bg, z, iqg, rhomass,a_massive, ep !RL adding adotoa, dgrho, v1_bg, v2_bg for the axions pert IC
    integer l,i, nu_i, j, ind
    integer, parameter :: i_clxg=1,i_clxr=2,i_clxc=3, i_clxb=4, &
         i_qg=5,i_qr=6,i_vb=7,i_pir=8, i_eta=9, i_aj3r=10,i_clxq=11,i_vq=12, &
         i_clxax=13,i_v_ax=14, i_dphi_ax=15, i_dphidot_ax=16 !RL adding KG variables

    integer, parameter :: i_max = i_dphidot_ax !i_max = i_v_ax - RL changed
    real(dl) initv(6,1:i_max), InitVec(1:i_max)
    real(dl) Ra,omr,AA,FF,frac,aeq !DM: quantities for axion isocurvature
    !variables for finding a_osc/aeq and which is first
    double precision lnamin,lnamax,dlna
    double precision, allocatable::a_arr(:),lna_arr(:),f_arr(:),fax_arr(:)
    double precision, allocatable::fmatt_arr(:),frad_arr(:)
    double precision,allocatable::lna_arr_small(:),a_arr_small(:)
    double precision, allocatable::fmatt_arr_small(:)
    double precision, allocatable::fax_arr_small(:),frad_arr_small(:)
    double precision, allocatable:: f_arr_small(:),a_eq_buff(:)
    double precision, allocatable::aeq_buff(:)
    integer llo_out,lhi_out
    double precision ho_out,a0_out,b0_out,Konstant
    !double precision H_eV

    !write(*, *) 'We are at initial settings, what is my logical switch?', EV%oscillation_started
    !write(*, *) 'That is funny. What is this other ktau switch?', EV%high_ktau_neutrino_approx
    !get hubble in units of eV in terms of standarf formula
    !H_eV=1.d14*6.5821d-27*dble(CP%H0)/(MPC_in_sec*c)
    !!write(*, *) 'Rayne, H_eV', H_eV
    


    !print*, 'mass Renee', CP%ma, 'CP%m_ovH0', CP%m_ovH0
    !renormalize mass in terms of hubble and use throughout boltzmann code as needed and store in big internal data structure eV
    !print*,CP%m_ovH0
    if (CP%flat) then
       EV%k_buf=EV%q
       EV%k2_buf=EV%q2
       EV%Kf(1:EV%MaxlNeeded)=1._dl
    else
       EV%k2_buf=EV%q2-CP%curv
       EV%k_buf=sqrt(EV%k2_buf)

       do l=1,EV%MaxlNeeded
          EV%Kf(l)=1._dl-CP%curv*(l*(l+2))/EV%k2_buf
       end do
    end if

    k=EV%k_buf
    k2=EV%k2_buf

    do j=1,EV%MaxlNeeded
       EV%denlk(j)=denl(j)*k*j
       EV%denlk2(j)=denl(j)*k*EV%Kf(j)*(j+1)
       EV%polfack(j)=polfac(j)*k*EV%Kf(j)*denl(j)
    end do

    !Get time to switch off tight coupling
    !The numbers here are a bit of guesswork
    !The high k increase saves time for very small loss of accuracy
    !The lower k ones are more delicate. Nead to avoid instabilities at same time
    !as ensuring tight coupling is accurate enough
    if (EV%k_buf > epsw) then
       if (EV%k_buf > epsw*5) then
          ep=ep0*5/AccuracyBoost
          if (HighAccuracyDefault) ep = ep*0.65
       else
          ep=ep0
       end if
    else
       ep=ep0
    end if
    if (second_order_tightcoupling) ep=ep*2
    EV%TightSwitchoffTime = min(tight_tau,Thermo_OpacityToTime(EV%k_buf/ep))


    y=0

    !  k*tau, (k*tau)**2, (k*tau)**3
    x=k*tau
    x2=x*x
    x3=x2*x
    rhomass =  sum(grhormass(1:CP%Nu_mass_eigenstates))
    grhonu=rhomass+grhornomass

    om = (grhob+grhoc)/sqrt(3*(grhog+grhonu))  ! DM: for a_i<<a_osc, axions do not contribute in background expansion 
    ! DM: therefore adiabtic initial conditions completely unchanged in this case.    
    omtau=om*tau
    Rv=grhonu/(grhonu+grhog)

    Rg = 1-Rv
    Rc=(CP%omegac)/(CP%omegac+CP%omegab) ! DM: no axions in here for same reason
    Rb=1-Rc
    Rp15=4*Rv+15

    if (CP%Scalar_initial_condition > initial_nummodes) &
         stop 'Invalid initial condition for scalar modes'

    a=tau*adotrad*(1+omtau/4)
    a2=a*a
    adotoa = 1.0_dl/(a*dtauda(a))!1.0_dl/(tau*(1+omtau/4)) !RL added, which is just adotrad/a - updated to using dtauda
    call spline_out(loga_table,phinorm_table,phinorm_table_ddlga,ntable,dlog10(a),v1_bg)
    !write(*, *) 'Rayne, is the spline in IC successful?', v1_bg
    initv=0

    !  Set adiabatic initial conditions

    chi=1  !Get transfer function for chi
    initv(1,i_clxg)=-chi*EV%Kf(1)/3*x2*(1-omtau/5)
    initv(1,i_clxr)= initv(1,i_clxg)
    initv(1,i_clxb)=0.75_dl*initv(1,i_clxg)
    initv(1,i_clxc)=initv(1,i_clxb)
    initv(1,i_qg)=initv(1,i_clxg)*x/9._dl
    initv(1,i_qr)=-chi*EV%Kf(1)*(4*Rv+23)/Rp15*x3/27
    initv(1,i_vb)=0.75_dl*initv(1,i_qg)
    initv(1,i_pir)=chi*4._dl/3*x2/Rp15*(1+omtau/4*(4*Rv-5)/(2*Rv+15))
    initv(1,i_aj3r)=chi*4/21._dl/Rp15*x3
    initv(1,i_eta)=-chi*2*EV%Kf(1)*(1 - x2/12*(-10._dl/Rp15 + EV%Kf(1)))
    initv(1,i_clxax)=0
    initv(1,i_v_ax)=0
    !RL tested that adotoa is not pre-defined
    !write(*, *) 'Rayne, adotrad, 1/dtauda(a), their fractional difference', adotrad, 1/dtauda(a), adotrad*dtauda(a) - 1.0_dl  - RL tested that the fractional difference is ~3e-8 and hence deleted dtauda since adotrad is already available !UPDATE: used  1/dtauda(a) since that's supposed to be more accurate
    !write(*, *) 'Rayne, testing the initial adotoa' !RL: when using 1/dtauda(a), it give exactly the same adotoa as the first entry in derivs
    !write(*, '(36e52.42,\)') adotoa

    !total perturbations: radiation, nutrinos (all massless in RD), and matter terms UNFINISHED
    !  8*pi*a*a*SUM[rho_i*clx_i]
    !  grho = 8*pi*rho*a**2
    dgrho = (grhonu + grhog)*initv(1,i_clxg)/a2 + (grhob + grhoc)*initv(1,i_clxb)/a !RL added, adiabatic initial pert conditions
    !write(*, *) 'Rayne, what is the fraction of the matter contribution in dgrho?', ((grhob + grhoc)*initv(1,i_clxb)/a)/dgrho
    z = (0.5_dl*dgrho/k -initv(1,i_eta)*k/2)/adotoa !RL (we're only specifying adiabatic IC here, the sign change can be left later to be done collectively)
    !write(*, *) 'Rayne, what is your dgrho?', dgrho !Why is there an extra negative sign...
    !write(*, *) 'What is grhob, grhoc, initv(1,i_clxb)?', grhob, grhoc, initv(1,i_clxb)
    !dv2_analytical = 
    !dv1_analytical = 
    !clxax_analytical = (v2_bg*dv2_analytical/a2 + (CP%m_ovH0**2.0d0)*v1_bg*dv1_analytical)*2.0d0/grhoax_kg
    !write(*, *) 'Rayne, is your initial dgrho, etak, adotoa, z correct?', dgrho, initv(1,i_eta)*k/2, adotoa, z !RL tested z to be ~7e-9 fractionally from the first entry of z in derivs
    !write(*, *) 'Rayne, what are your calculated k*z, clxcdot, and their fractional differences?', k*z, -chi*EV%Kf(1)*(1-omtau/5)*k*x/2, -chi*EV%Kf(1)*(1-omtau/5)*k*x/2/(k*z) - 1.0_dl !RL tested to be ~1e-8 fractional different
    initv(1,i_dphi_ax) = ((CP%m_ovH0*CP%H0_in_Mpc_inv)**2.0d0) * &
         &(chi*EV%Kf(1)*(1-omtau/5)*k*x/2) * a2 * v1_bg/(210.0d0 * (adotoa**3.0d0)) !RL - clxcdotmethod
    initv(1,i_dphidot_ax) = (CP%m_ovH0**2.0d0) * CP%H0_in_Mpc_inv &
         &* (chi*EV%Kf(1)*(1-omtau/5)*k*x/2) * a2 * v1_bg/(35.0d0 * (adotoa**2.0d0)) !RL- clxcdotmethod
    
    !initv(1,i_dphi_ax) = ((CP%m_ovH0*CP%H0_in_Mpc_inv)**2.0d0) * (k * z) * a2 * v1_bg/(210.0d0 * (adotoa**3.0d0)) !RL - kzmethod
    !initv(1,i_dphidot_ax) = (CP%m_ovH0**2.0d0) * CP%H0_in_Mpc_inv * (k * z) * a2 * v1_bg/(35.0d0 * (adotoa**2.0d0)) !RL- kzmethod
    !write(*, *) 'Rayne what are your CP%m_ovH0, CP%H0_in_Mpc_inv, chi*EV%Kf(1)*(1-omtau/5)*k*x/2, a2, v1_bg, adotoa, initv(1,i_dphidot_ax) in the init?', CP%m_ovH0, CP%H0_in_Mpc_inv, chi*EV%Kf(1)*(1-omtau/5)*k*x/2, a2, v1_bg, adotoa, initv(1,i_dphidot_ax) !Should be right if there's a sign flip

    if (CP%Scalar_initial_condition/= initial_adiabatic) then
       !CDM isocurvature

       initv(2,i_clxg)= Rc*omtau*(-2._dl/3 + omtau/4)
       initv(2,i_clxr)=initv(2,i_clxg)
       initv(2,i_clxb)=initv(2,i_clxg)*0.75_dl
       initv(2,i_clxc)=1+initv(2,i_clxb)
       initv(2,i_qg)=-Rc/9*omtau*x
       initv(2,i_qr)=initv(2,i_qg)
       initv(2,i_vb)=0.75_dl*initv(2,i_qg)
       initv(2,i_pir)=-Rc*omtau*x2/3/(2*Rv+15._dl)
       initv(2,i_eta)= Rc*omtau*(1._dl/3 - omtau/8)*EV%Kf(1)
       initv(2,i_aj3r)=0
       !RL 072926: baryon iso is the baryon-CDM compensation (row 3 = CDM row * omega_b/omega_c) and
       !does not exist at omega_c = 0, which axfrac = 1 hits exactly. Build it only where defined;
       !elsewhere it stays zero (initv=0 above) so initial_vector avoids 0*NaN. The refusal for modes
       !that actually read row 3 is in CAMB_ValidateParams -- this routine runs inside an OMP k-loop.
       if (Rc > 0._dl) then
          initv(3,:) = initv(2,:)*(Rb/Rc)
          initv(3,i_clxc) = initv(3,i_clxb)
          initv(3,i_clxb)= initv(3,i_clxb)+1
       end if

       !neutrino isocurvature density mode

       initv(4,i_clxg)=Rv/Rg*(-1 + x2/6)
       initv(4,i_clxr)=1-x2/6
       initv(4,i_clxc)=-omtau*x2/80*Rv*Rb/Rg
       initv(4,i_clxb)= Rv/Rg/8*x2
       iqg = - Rv/Rg*(x/3 - Rb/4/Rg*omtau*x)
       initv(4,i_qg) =iqg
       initv(4,i_qr) = x/3
       initv(4,i_vb)=0.75_dl*iqg
       initv(4,i_pir)=x2/Rp15
       initv(4,i_eta)=EV%Kf(1)*Rv/Rp15/3*x2

       !neutrino isocurvature velocity mode

       initv(5,i_clxg)=Rv/Rg*x - 2*x*omtau/16*Rb*(2+Rg)/Rg**2
       initv(5,i_clxr)=-x -3*x*omtau*Rb/16/Rg
       initv(5,i_clxc)=-9*omtau*x/64*Rv*Rb/Rg
       initv(5,i_clxb)= 3*Rv/4/Rg*x - 9*omtau*x/64*Rb*(2+Rg)/Rg**2
       iqg = Rv/Rg*(-1 + 3*Rb/4/Rg*omtau+x2/6 +3*omtau**2/16*Rb/Rg**2*(Rg-3*Rb))
       initv(5,i_qg) =iqg
       initv(5,i_qr) = 1 - x2/6*(1+4*EV%Kf(1)/(4*Rv+5))
       initv(5,i_vb)=0.75_dl*iqg
       initv(5,i_pir)=2*x/(4*Rv+5)+omtau*x*6/Rp15/(4*Rv+5)
       initv(5,i_eta)=2*EV%Kf(1)*x*Rv/(4*Rv+5) + omtau*x*3*EV%Kf(1)*Rv/32*(Rb/Rg - 80/Rp15/(4*Rv+5))
       initv(5,i_aj3r) = 3._dl/7*x2/(4*Rv+5)

       !quintessence isocurvature mode

       !axion isocurvature mode
       ! DM: See pdf notes on axion theory
       ! Derived by Dan Grin using a Matrix ODE formalism
       ! All assume tau_i<<tau_osc
       initv(6,:)=0.0d0
       Ra=grhoax/(grhog+grhonu)/a_osc**3.0d0 !why is there a grhoax that's not declared in the subroutine? also grhog and a_osc, let's ask Dan when sending a follow up !RL 050224: There's a grhoax declared in modules. It's different than the dorp used at the reference point. This notation is really confusing
       omr=(grhog+grhonu)/(grhob+grhoc)
       frac=grhoax/(grhoax+grhob+grhoc) ! Fraction of matter in axions	





       ! print*,a_osc,CP%a_osc
       !calculate normalization factor for isocurvature perturbation	 
       if(a_osc.le.CP%aeq) then
          FF=1.0d0 ! DM: F(f) sets a_0 from a_osc FF=1 if aosc<aeq
          Konstant=1.0d0-frac
       else		 											
          FF=((1.0d0-frac)+frac*(CP%aeq/a_osc)**3.0d0)**(-1.0d0)
          Konstant=(1.0d0-frac)/((1.0d0-frac)+frac*((CP%aeq/a_osc)**3.0d0))
       end if


       !RL reading: noticed that there's this axion isocurvature initial condition specifications that we might need to incorporate into the field variables too in the future. But stay tuned for now
       initv(6,i_clxax)= 0._dl
       !1.0d0-(x*x)/(10.0d0)&
       !     -Konstant*(x*x)*omtau/(180.0d0)+(1.0d0/600.0d0)*(x**4.0d0)&
       !     +(53.0d0*(Konstant**2.0d0)*(x**2.0d0)*(omtau**2.0d0))/(140.0d0*5.0d0*16.0d0)
       ! Normalise to \delta_a=1. 
       initv(6,i_v_ax)= 0._dl
       !-x/5.0d0+(x*x*x)/30.0d0+Konstant*x*omtau/(30.0d0)&
       !     -x*(Konstant**2.0d0)*(omtau**2.0d0)/(84.d0)&
       !     +(11.0d0*(Konstant**3.0d0)*x*(omtau**(3.0d0)))/(42.0d0*64.0d0)&
       !     -Konstant*(x**3.0d0)*omtau/(120.0d0)-(x**5.0d0)/(600.0d0)&
       !     -(11.0d0*(Konstant**4.0d0)*x*(omtau**4.0d0))/(63.0d0*128.0d0)&
       !     +(929.0d0*(Konstant**2.0d0))*(x**3.0d0)*(omtau**2.0d0)/(3780.0d0*16.0d0*5.0d0)&
       !     +((CP%m_ovH0**2.0d0)*x*(omtau**4.0d0)/(225.0d0*5.0d0*64.0d0*CP%omegar&
       !     *(((CP%omegab+CP%omegac+CP%omegaax)/(4.0d0*CP%omegar))**4.0d0)))*(FF**4.0d0)
       initv(6,i_dphi_ax) = 1._dl !RL 011626 simpler field normalization where delta v_ini = 1
       !RL 111225: k correction (1st term) as well as m/H correction (2nd term)
       initv(6,i_dphidot_ax) = - initv(6,i_dphi_ax)*(k2/(3._dl*CP%H0_in_Mpc_inv) + &
            &a2 * (CP%m_ovH0**2.0d0) * CP%H0_in_Mpc_inv/5._dl)/adotoa


       !		initv(6,i_clxax)=1.0d0-(x*x)/(10.0d0)
       !		initv(6,i_v_ax)=-x/5.0d0

       ! Normalisation of *gamma* perturbation

       AA=Ra*(omr**4.0d0)*(FF**4.0d0)
       !		print*,AA,omtau

       initv(6,i_clxg)= -AA/3.*omtau**4.
       initv(6,i_clxr)=initv(6,i_clxg)
       initv(6,i_clxb)=0.75_dl*initv(6,i_clxg)
       initv(6,i_clxc)=initv(6,i_clxb)
       initv(6,i_eta)=-(0.50*initv(6,i_clxg))
       !DM: this eta is -2eta_s, see notes and below setting y(2)

       initv(6,i_qg)=initv(6,i_clxg)*x/15.0d0
       initv(6,i_qr)=initv(6,i_clxg)*x/15.0d0
       initv(6,i_vb)=0.75_dl*initv(6,i_qg)
       initv(6,i_pir)= (-9.0d0*Rb)*omtau*initv(6,i_clxg)*(1.0d0-Konstant)/(5.0d0*(75.0d0+4.0d0*Rv))
       !DM: seems pir=2 sigma_nu
       initv(6,i_aj3r)= initv(6,i_pir)*x/(7.0d0)
       !DM: assuming this is exactly F_{\nu 3}.
       !	 print*,initv(6,:)
       !		print*,AA
    end if



    !write(*, *) 'Rayne, what is y(2) before this initialization?', -InitVec(i_eta)*k/2
    !write(*, *) 'Rayne, what is the scalar initial conditions then?', -initv(1,i_eta)*k/2


    if (CP%Scalar_initial_condition==initial_vector) then
       InitVec = 0
       do i=1,initial_nummodes
          InitVec = InitVec+ initv(i,:)*CP%InitialConditionVector(i)
       end do
    else
       InitVec = initv(CP%Scalar_initial_condition,:)
       if (CP%Scalar_initial_condition==initial_adiabatic) InitVec = -InitVec !RL note: sign change due to change in gauge conventions (eta = +1 or -1)
       !So we start with chi=-1 as before
    end if

    y(1)=a
    y(2)= -InitVec(i_eta)*k/2
    !get eta_s*k, where eta_s is synchronous gauge variable
    !  CDM
    y(3)=InitVec(i_clxc)

    !  Baryons
    y(4)=InitVec(i_clxb)
    y(5)=InitVec(i_vb)


    ! Axions

    y(EV%a_ix)=InitVec(i_clxax)
    y(EV%a_ix+1)=InitVec(i_v_ax)
    !RL adding KG variables
    !y(EV%a_kg_ix) = InitVec(i_dphi_ax)
    EV%renorm_c = sqrt(CP%rhorefp_ovh2)*CP%H0*((k*CP%tau_osc)**2._dl)/CP%m_ovH0
    y(EV%a_kg_ix) = InitVec(i_dphi_ax)/EV%renorm_c !RL 050324
    y(EV%a_kg_ix+1) = InitVec(i_dphidot_ax) !RL 050324

    !  Photons
    y(EV%g_ix)=InitVec(i_clxg)
    y(EV%g_ix+1)=InitVec(i_qg)

    if (.not. is_cosmological_constant) then
       y(EV%w_ix) = InitVec(i_clxq) !ppf: Gamma=0, i_clxq stands for i_Gamma
    end if

    !  Neutrinos
    y(EV%r_ix)=InitVec(i_clxr)
    y(EV%r_ix+1)=InitVec(i_qr)
    y(EV%r_ix+2)=InitVec(i_pir)

    if (EV%lmaxnr>2) then
       y(EV%r_ix+3)=InitVec(i_aj3r)
    endif

    if (CP%Num_Nu_massive == 0) return

    do nu_i = 1, CP%Nu_mass_eigenstates
       EV%MassiveNuApproxTime(nu_i) = Nu_tau_massive(nu_i)
       a_massive =  20000*k/nu_masses(nu_i)*AccuracyBoost*lAccuracyBoost
       if (a_massive >=0.99) then
          EV%MassiveNuApproxTime(nu_i)=CP%tau0+1
       else if (a_massive > 17.d0/nu_masses(nu_i)*AccuracyBoost) then
          EV%MassiveNuApproxTime(nu_i)=max(EV%MassiveNuApproxTime(nu_i),DeltaTime(0._dl,a_massive, 0.01_dl))
       end if
       ind = EV%nu_ix(nu_i)
       do  i=1,EV%nq(nu_i)
          y(ind:ind+2)=y(EV%r_ix:EV%r_ix+2)
          if (EV%lmaxnu_tau(nu_i)>2) y(ind+3)=InitVec(i_aj3r)
          ind = ind + EV%lmaxnu_tau(nu_i)+1
       end do
    end do

  end subroutine initial


  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine initialt(EV,yt,tau)
    !  Initial conditions for tensors
    use ThermoData
    implicit none
    real(dl) bigR,tau,x,aj3r,elec, pir, rhomass
    integer l
    type(EvolutionVars) EV
    real(dl) k,k2 ,a, omtau
    real(dl) yt(EV%nvart)
    real(dl) tens0, ep, tensfac

    if (CP%flat) then
       EV%aux_buf=1._dl
       EV%k2_buf=EV%q2
       EV%k_buf=EV%q
       EV%Kft(1:EV%MaxlNeededt)=1._dl !initialize for flat case
    else
       EV%k2_buf=EV%q2-3*CP%curv
       EV%k_buf=sqrt(EV%k2_buf)
       EV%aux_buf=sqrt(1._dl+3*CP%curv/EV%k2_buf)
    endif

    k=EV%k_buf
    k2=EV%k2_buf

    do l=1,EV%MaxlNeededt
       if (.not. CP%flat) EV%Kft(l)=1._dl-CP%curv*((l+1)**2-3)/k2
       EV%denlkt(1,l)=k*denl(l)*l !term for L-1
       tensfac=real((l+3)*(l-1),dl)/(l+1)
       EV%denlkt(2,l)=k*denl(l)*tensfac*EV%Kft(l) !term for L+1
       EV%denlkt(3,l)=k*denl(l)*tensfac**2/(l+1)*EV%Kft(l) !term for polarization
       EV%denlkt(4,l)=k*4._dl/(l*(l+1))*EV%aux_buf !other for polarization
    end do

    if (k > 0.06_dl*epsw) then
       ep=ep0
    else
       ep=0.2_dl*ep0
    end if

    !    finished_tightcoupling = ((k/opacity > ep).or.(1._dl/(opacity*tau) > ep))
    EV%TightSwitchoffTime = min(tight_tau,Thermo_OpacityToTime(EV%k_buf/ep))

    a=tau*adotrad
    rhomass =  sum(grhormass(1:CP%Nu_mass_eigenstates))
    omtau = tau*(grhob+grhoc)/sqrt(3*(grhog+rhomass+grhornomass))

    if (DoTensorNeutrinos) then
       bigR = (rhomass+grhornomass)/(rhomass+grhornomass+grhog)
    else
       bigR = 0._dl
    end if

    x=k*tau

    yt(1)=a
    tens0 = 1

    yt(2)= tens0
    !commented things are for the compensated mode with magnetic fields; can be neglected
    !-15/28._dl*x**2*(bigR-1)/(15+4*bigR)*Magnetic*(1-5./2*omtau/(2*bigR+15))

    elec=-tens0*(1+2*CP%curv/k2)*(2*bigR+10)/(4*bigR+15) !elec, with H=1

    !shear
    yt(3)=-5._dl/2/(bigR+5)*x*elec
    !          + 15._dl/14*x*(bigR-1)/(4*bigR+15)*Magnetic*(1 - 15./2*omtau/(2*bigR+15))

    yt(4:EV%nvart)=0._dl

    !  Neutrinos
    if (DoTensorNeutrinos) then
       pir=-2._dl/3._dl/(bigR+5)*x**2*elec
       !           + (bigR-1)/bigR*Magnetic*(1-15./14*x**2/(15+4*bigR))
       aj3r=  -2._dl/21._dl/(bigR+5)*x**3*elec !&
       !           + 3._dl/7*x*(bigR-1)/bigR*Magnetic
       yt(EV%r_ix+2)=pir
       yt(EV%r_ix+3)=aj3r
       !Should set up massive too, but small anyway..
    end if

  end subroutine initialt

  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine initialv(EV,yv,tau)
    !  Initial conditions for vectors

    implicit none
    real(dl) bigR,Rc,tau,x,pir
    type(EvolutionVars) EV
    real(dl) k,k2 ,a, omtau
    real(dl) yv(EV%nvarv)

    if (CP%flat) then
       EV%k2_buf=EV%q2
       EV%k_buf=EV%q
    else
       stop 'Vectors not supported in non-flat models'
    endif

    k=EV%k_buf
    k2=EV%k2_buf

    omtau = tau*(grhob+grhoc)/sqrt(3*(grhog+grhornomass))

    a=tau*adotrad*(1+omtau/4)

    x=k*tau

    bigR = (grhornomass)/(grhornomass+grhog)
    Rc=CP%omegac/(CP%omegac+CP%omegab)

    yv(1)=a


    yv(2)= vec_sig0*(1- 15._dl/2*omtau/(4*bigR+15)) + 45._dl/14*x*Magnetic*(BigR-1)/(4*BigR+15)
    !qg
    yv(4)= vec_sig0/3* (4*bigR + 5)/(1-BigR)*(1  -0.75_dl*omtau*(Rc-1)/(bigR-1)* &
         (1 - 0.25_dl*omtau*(3*Rc-2-bigR)/(BigR-1))) &
         -x/2*Magnetic
    yv(3)= 3._dl/4*yv(4)

    yv(5:EV%nvarv) = 0

    !        if (.false.) then
    !         yv((EV%lmaxv-1+1)+(EV%lmaxpolv-1)*2+3+1) = vec_sig0/6/bigR*x**2*(1+2*bigR*omtau/(4*bigR+15))
    !         yv((EV%lmaxv-1+1)+(EV%lmaxpolv-1)*2+3+2) = -2/3._dl*vec_sig0/bigR*x*(1 +3*omtau*bigR/(4*bigR+15))
    !         yv((EV%lmaxv-1+1)+(EV%lmaxpolv-1)*2+3+3) = 1/4._dl*vec_sig0/bigR*(5+4*BigR)
    !         yv((EV%lmaxv-1+1)+(EV%lmaxpolv-1)*2+3+4) =1/9.*x*vec_sig0*(5+4*bigR)/bigR
    !         yv(4) = 0
    !         yv(3)= 3._dl/4*yv(4)
    !          return
    !        end if

    !  Neutrinos
    !q_r
    yv((EV%lmaxv-1+1)+(EV%lmaxpolv-1)*2+3+1) = -1._dl/3*vec_sig0*(4*BigR+5)/bigR &
         + x**2*vec_sig0/6/BigR +0.5_dl*x*(1/bigR-1)*Magnetic
    !pi_r
    pir=-2._dl/3._dl*x*vec_sig0/BigR - (1/bigR-1)*Magnetic
    yv((EV%lmaxv-1+1)+(EV%lmaxpolv-1)*2+3+1 +1)=pir
    yv((EV%lmaxv-1+1)+(EV%lmaxpolv-1)*2+3+1 +2)=3._dl/7*x*Magnetic*(1-1/BigR)

  end subroutine initialv


  subroutine outtransf(EV, y, tau,Arr)
    !write out clxc, clxb, clxg, clxn, clxax
    use Transfer
    implicit none
    type(EvolutionVars) EV

    real(dl) clxc, clxb, clxg, clxr, k,k2
    real(dl) grho,gpres,dgrho,dgq,a,a2
    real Arr(:)
    real(dl) y(EV%nvar)
    real(dl) growth,clxtot !DM: growth vars, clxtot is unnecessary, just needs a name
    real(dl) tau
    real(dl) dorp
    real(dl) gr
    real(dl) clxax, dv1, dv2, v1_bg, v2_bg, wcorr_coeff, drhoax_kg, grhoax_kg, clxax_kg, grhoax_t !RL added clxax_kg w_ax_p1 not needed for transfer function

    integer i	

    a    = y(1)
    !write(*, *) 'In transf, a', a
    a2 = a*a
    clxc = y(3)
    clxb = y(4)
    !write(*, *) 'Rayne, in transf, k, EV%no_nu_multpoles, EV%no_phot_multpoles', EV%k_buf, EV%no_nu_multpoles, EV%no_phot_multpoles
    if (EV%no_nu_multpoles) then
       clxr=0
    else
       clxr = y(EV%r_ix)
    end if

    if (EV%no_phot_multpoles) then
       clxg=0
    else
       clxg = y(EV%g_ix)
    end if

    k    = EV%k_buf
    k2   = EV%k2_buf
    !write(*, *) 'Rayne, we''re in outtransf for k = ', k
    ! Axions
    clxax=y(EV%a_ix)

    !if (.not. EV%oscillation_started) then !RL
    !   call spline_out(loga_table,rhoaxh2ovrhom_logtable,rhoaxh2ovrhom_logtable_buff,ntable,dlog10(a),gr)
    !write(*,*) a, gr
    !compute log10 of density use of interpolation
    !print*,a,grhoc,grhom*(1.0d1**(dble(gr)))/((grhoc+grhob)/(a**3.0d0))
    !   if (gr .eq. gr) then
    !delog it and multiply by physical units
    !      dorp=grhom*(1.0d1**(dble(gr)))
    !   else
    !      dorp=0.0d0
    !   endif
    !else
    !dorp=grhom*CP%rhorefp_hsq*((CP%a_osc/a)**3.0d0)
    !endif

    !RL added clxax_kg
    if (.not. EV%oscillation_started) then
       dv1 = y(EV%a_kg_ix)
       !write(*, *) 'Check that the dv2s are consistent:', yprime(EV%a_kg_ix)/y(EV%a_kg_ix + 1) - 1.0d0
       dv2 = y(EV%a_kg_ix+1)
       !spline to get the phi and phidot at this point
       call spline_out(loga_table,phinorm_table,phinorm_table_ddlga,ntable,dlog10(a),v1_bg)
       call spline_out(loga_table,phidotnorm_table,phidotnorm_table_ddlga,ntable,dlog10(a),v2_bg)
       !obtain 1+w for momentum
       !w_ax_p1 = (2.0_dl*(v2_bg**2.0_dl)/a2)/((v2_bg**2.0_dl)/a2 + (CP%m_ovH0*v1_bg)**2.0_dl)
       !drhoax and clxax calculated from KG
       drhoax_kg = (v2_bg*dv2/a2 + (CP%m_ovH0**2.0d0)*v1_bg*dv1*EV%renorm_c)*2.0d0 !RL 050324
       !grhoax from background KG
       grhoax_kg = (v2_bg)**2.0d0/a2+(CP%m_ovH0*v1_bg)**2.0d0
       !translate to CAMB normalization
       if (v1_bg .eq. v1_bg .and. v2_bg .eq. v2_bg) then
          dorp = grhom*grhoax_kg/(CP%H0**2.0d0/1.0d4)
       else
          dorp=0.0d0
       end if

       clxax_kg = drhoax_kg/grhoax_kg 
       !Don't need momentum for transfer function
    else
       !past tauosc, directly obtain clxax
       !RL adding w correction to the background 
       !wcorr_coeff = CP%ah_osc*CP%a_osc/(CP%m_ovH0*(CP%H0/100.0d0))
       wcorr_coeff = CP%ahosc_ETA*CP%a_osc/(CP%m_ovH0*(CP%H0/100.0d0)) !RL082924
       !!write(*, *) 'In outtransf, CP%wEFA_c', CP%wEFA_c
       !w_ax_p1 = 1.0_dl + 3.0d0*((wcorr_coeff/a2)**2.0d0)/2.0d0
       dorp=grhom*CP%rhorefp_ovh2*((CP%a_osc/a)**3.0d0)*dexp((wcorr_coeff**2.0d0)*3.0d0*&
            &CP%wEFA_c*(1.0d0/(a2**2.0d0) - 1.0d0/(CP%a_osc**4.0d0))/4.0d0)
       !dorp=grhom*CP%rhorefp_hsq*((CP%a_osc/a)**3.0d0)
       clxax_kg = y(EV%a_kg_ix)
    end if

    grhoax_t=dorp*(a**2.0d0)
    Arr(Transfer_kh) = k/(CP%h0/100._dl)
    if (CP%axion_isocurvature) then
       ! RL 051525 + RL 06**26: in axion-iso mode all species' perturbations
       ! were sourced by the same IC `delta-phi-tilde_i = 1`, so every output
       ! transfer is in "per-delta-phi-tilde_i" units.  Multiply each output
       ! by phi_i/2 = phinorm_table(1)/2 to convert to "per-delta_a(tau_i)=1"
       ! units, matching the IC role that delta_c(tau_i)=1 plays in public-
       ! CAMB CDM-iso.  Single overall multiplier per output -- preserves
       ! the dgrho compensation (clxax_kg, clxc, clxb stay raw in the dgrho
       ! sum below; the Transfer_tot output picks up the same phi_i/2 factor).
       Arr(Transfer_cdm)   = (clxc/k2)     * (phinorm_table(1)/2._dl)
       Arr(Transfer_b)     = (clxb/k2)     * (phinorm_table(1)/2._dl)
       Arr(Transfer_g)     = (clxg/k2)     * (phinorm_table(1)/2._dl)
       Arr(Transfer_r)     = (clxr/k2)     * (phinorm_table(1)/2._dl)
       Arr(Transfer_axion) = (clxax_kg/k2) * (phinorm_table(1)/2._dl)
    else
       Arr(Transfer_cdm)   = clxc/k2
       Arr(Transfer_b)     = clxb/k2
       Arr(Transfer_g)     = clxg/k2
       Arr(Transfer_r)     = clxr/k2
       Arr(Transfer_axion) = clxax_kg/k2
    end if

    !write(*, *) 'Rayne, transfer function for axions, transfer function for cdm', Arr(Transfer_axion), Arr(Transfer_cdm)
    !write(111608, '(36e52.42,\)') k, a, CP%a_osc, tau, CP%tau_osc, Arr(Transfer_axion), Arr(Transfer_cdm), Arr(Transfer_b), Arr(Transfer_g)

    dgrho = 0
    grho =  0

    if (CP%Num_Nu_Massive > 0) then
       call MassiveNuVars(EV,y,a,grho,gpres,dgrho,dgq)
       if (CP%axion_isocurvature) then
          Arr(Transfer_nu) = (dgrho/grho/k2) * (phinorm_table(1)/2._dl)
       else
          Arr(Transfer_nu) = dgrho/grho/k2
       end if
    else
       Arr(Transfer_nu) = 0
    end if

    !If we want DE perturbations to get \delta\rho/\rho_m    
    !dgrho=dgrho+y(EV%w_ix)*grhov*a**(-1-3*w_lam)
    !Arr(Transfer_r) = y(EV%w_ix)/k2

    !dgrho = dgrho+(clxc*grhoc + clxb*grhob)/a 
    !grho =  grho+(grhoc+grhob)/a + grhov*a**(-1-3*w_lam)

    ! DM: only include axions if they are clustered: a_osc kluge for bias as in Hlozek et al 2015. Be careful when a_osc.ge.1
    !! DM16: this modification with ma\geq 1e-25 eV to treat axions in non-linear lensing S4 fiducial model.
    ! See additional changes in halofit, and check for consistency.

    !if (a_osc.lt.CP%aeq) then
    !if (CP%ma.ge.1.e-25) then
    if (CP%ma .ge. 1.44e-32_dl) then !RL 070324
       dgrho = dgrho+(clxc*grhoc + clxb*grhob)/a+grhoax_t*clxax_kg !RL replaced axion pert with kg
       grho =  grho+(grhoc+grhob)/a+grhoax_t
    else
       !	write(*,*)'axions are not clustering in P(k)'
       dgrho = dgrho+(clxc*grhoc + clxb*grhob)/a
       grho =  grho+(grhoc+grhob)/a
    end if

    ! DM: output the growth rate

    !!write(*, *) 'RL, in outtransf, before calling GrowthRate, EV%oscillation_started, CP%a_osc', EV%oscillation_started, CP%a_osc
    !call GrowthRate(EV,y,tau,k,a,growth,clxtot)
    !!write(*, *) 'RL, in outtransf, after calling GrowthRate'
    
    growth = 0._dl
    Arr(Transfer_f) = growth
    if (CP%axion_isocurvature) then
       Arr(Transfer_tot) = (dgrho/grho/k2) * (phinorm_table(1)/2._dl)
    else
       Arr(Transfer_tot) = dgrho/grho/k2
    end if

  end subroutine outtransf

  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine derivs(EV,n,tau,ay,ayprime)
    !  Evaluate the time derivatives of the perturbations
    !  ayprime is not necessarily GaugeInterface.yprime, so keep them distinct
    use ThermoData
    use MassiveNu
    use constants
    implicit none
    type(EvolutionVars) EV

    integer n,nu_i
    real(dl) ay(n),ayprime(n)
    real(dl) tau,w
    real(dl) k,k2

    !  Internal variables.

    real(dl) opacity
    real(dl) photbar,cs2,pb43,grho,slip,clxgdot, &
         clxcdot,clxbdot,adotdota,gpres,clxrdot,etak
    real(dl) q,aq,v
    real(dl) G11_t,G30_t, wnu_arr(max_nu)

    real(dl) dgq,grhob_t,grhor_t,grhoc_t,grhog_t,grhov_t,sigma,polter
    real(dl) qgdot,qrdot,pigdot,pirdot,vbdot,dgrho,adotoa
    real(dl) a,a2,dz,z,clxc,clxb,vb,clxg,qg,pig,clxr,qr,pir
    real(dl) E2, dopacity
    integer l,i,ind, ind2, off_ix, ix
    real(dl) dgs,sigmadot !, ddz
    !non-flat vars
    real(dl) cothxor !1/tau in flat case
    !ppf
    real(dl) Gamma,S_Gamma,ckH,Gammadot,Fa,dgqe,dgrhoe, vT
    real(dl) w_eff, grhoT
    ! Axions
    real(dl) :: v1_bg, v2_bg, dv1, dv2, drhoax_kg, grhoax_kg, clxax_kg, clxax_kg_dot, u_ax_kg, u_ax_kg_dot !RL adding the KG variables
    real(dl) grhoax_t, clxax, gpres_ax
    real(dl) clxaxdot,v_ax,v_axdot
    real(dl) w_ax, w_ax_p1, wcorr_coeff, csquared_ax,cad2 !RL added w_ax_p1 for DE era still keeping w_ax for DM era
    real(dl) dorp
    real(dl) gr,kamnorm !RL modified mnorm
    real(dl) dv1_quasitest, clxax_quasitest !RL testing

    k=EV%k_buf
    k2=EV%k2_buf

    a=ay(1)
    a2=a*a

    etak=ay(2)

    !  CDM variables
    clxc=ay(3)

    !  Baryon variables
    clxb=ay(4)
    vb=ay(5)

    ! Axion variables
    clxax=ay(EV%a_ix)
    v_ax=ay(EV%a_ix+1)
    !RL adding KG variables
    !(!!)cs2 is not computed here since adotoa is not assigned yet. It is computed later!
    if (.not. EV%oscillation_started) then
       dv1 = ay(EV%a_kg_ix) !RL 050324: not renormalizing here, saving repeated multiplication and division in the RK process - instead renormalize the delta rho etc.
       dv2 = ay(EV%a_kg_ix+1)
       !spline to get the phi and phidot at this point
       call spline_out(loga_table,phinorm_table,phinorm_table_ddlga,ntable,dlog10(a),v1_bg)
       call spline_out(loga_table,phidotnorm_table,phidotnorm_table_ddlga,ntable,dlog10(a),v2_bg)
       !obtain 1+w for momentum
       w_ax_p1 = (2.0_dl*(v2_bg**2.0_dl)/a2)/((v2_bg**2.0_dl)/a2 + (CP%m_ovH0*v1_bg)**2.0_dl)
       w_ax = w_ax_p1 - 1.0_dl
       !drhoax and clxax calculated from KG
       drhoax_kg = (v2_bg*dv2/a2 + (CP%m_ovH0**2.0d0)*v1_bg*dv1*EV%renorm_c)*2.0d0 !RL 050324
       !grhoax from background KG
       grhoax_kg = (v2_bg)**2.0d0/a2+(CP%m_ovH0*v1_bg)**2.0d0
       if (v1_bg .eq. v1_bg .and. v2_bg .eq. v2_bg) then !RL inherited DG flag
          dorp = grhom*grhoax_kg/(CP%H0**2.0d0/1.0d4)
       else
          dorp=0.0d0
       end if
       clxax_kg = drhoax_kg/grhoax_kg 
       u_ax_kg = w_ax_p1*k*dv1*EV%renorm_c/(CP%H0_in_Mpc_inv*v2_bg) !RL 050324
    else !past tauosc
       !RL: EFA 
       v1_bg = 0._dl
       v2_bg = 0._dl
       drhoax_kg = 0._dl
       grhoax_kg = 0._dl
       !RL adding w correction to the background 
       !wcorr_coeff = CP%ah_osc*CP%a_osc/(CP%m_ovH0*(CP%H0/100.0d0))
       wcorr_coeff = CP%ahosc_ETA*CP%a_osc/(CP%m_ovH0*(CP%H0/100.0d0)) !RL082924
       w_ax = ((wcorr_coeff/a2)**2.0d0)*CP%wEFA_c
       w_ax_p1 = 1.0_dl + w_ax
       dorp=grhom*CP%rhorefp_ovh2*((CP%a_osc/a)**3.0d0)*dexp((wcorr_coeff**2.0d0)*3.0d0*&
            &CP%wEFA_c*(1.0d0/(a2**2.0d0) - 1.0d0/(CP%a_osc**4.0d0))/4.0d0) !RL 110923
       !write(*, *) 'In derivs, CP%wEFA_c', CP%wEFA_c
       !RL: now I need cad2 and from my approximation scheme it only depend on H at aosc
       !cad2 = w_ax*((w_ax + 7.0_dl/3.0_dl)/w_ax_p1) !RL 081324 - define cad2 after the w_total is obtained
       !dorp=grhom*CP%rhorefp_hsq*((CP%a_osc/a)**3.0d0)
       clxax_kg = ay(EV%a_kg_ix)
       u_ax_kg = ay(EV%a_kg_ix+1)
    end if

    !  Compute expansion rate from: grho 8*pi*rho*a**2

    grhob_t=grhob/a
    grhoc_t=grhoc/a
    grhor_t=grhornomass/a2
    grhog_t=grhog/a2
    if (is_cosmological_constant) then
       grhov_t=grhov*a2
       w_eff = -1_dl
    else
       !ppf
       w_eff=w_de(a)   !effective de
       grhov_t=grho_de(a)/a2
    end if



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Axion variables, updated 9/19/2013 to include DG spline implmentation
!!!!!!!updated 9/1/2015	to include simple cubic scaling when a>> aosc
    !if (tau .le. CP%tau_osc) then
    !if (tau .lt. CP%tau_osc) then
    !if (.not. EV%oscillation_started) then !RL
    !   call spline_out(loga_table,rhoaxh2ovrhom_logtable,rhoaxh2ovrhom_logtable_buff,ntable,dlog10(a),gr)
    !write(*,*) a, gr
    !compute log10 of density use of interpolation
    !print*,a,grhoc,grhom*(1.0d1**(dble(gr)))/((grhoc+grhob)/(a**3.0d0))
    !   if (gr .eq. gr) then
    !delog it and multiply by physical units
    !      dorp=grhom*(1.0d1**(dble(gr)))
    !   else
    !      dorp=0.0d0
    !   endif
    !else
    !   dorp=grhom*CP%rhorefp_hsq*((CP%a_osc/a)**3.0d0)     
    !endif
    grhoax_t=dorp*a2  !RL changed, originally (a**2.0d0)

    gpres_ax=w_ax*grhoax_t
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 



    !  Get sound speed and ionisation fraction.
    if (EV%TightCoupling) then
       call thermo(tau,cs2,opacity,dopacity)
    else
       call thermo(tau,cs2,opacity)
    end if


    gpres=gpres_ax
    grho=grhob_t+grhoc_t+grhor_t+grhog_t+grhov_t+grhoax_t

    !RL comments out this mnorm and, and computes the corrected csquared_ax after adotoa is assigned
    !mnorm=2.0d0*CP%m_ovH0*a*CP%H0*1.d3/(c)
    !csquared_ax=((k/mnorm)**2.0d0)/(1.0d0+(k/mnorm)**2.0d0)

    !total perturbations: matter terms first, then add massive nu, de and radiation 
    !  8*pi*a*a*SUM[rho_i*clx_i]
    dgrho=grhob_t*clxb+grhoc_t*clxc+grhoax_t*clxax_kg !RL replaced axion pert with kg
    !  8*pi*a*a*SUM[(rho_i+p_i)*v_i]
    !switch 9/22 v_ax*(1+w_ax) as the axion momentum variable to avoid infinities
    dgq=grhob_t*vb+grhoax_t*u_ax_kg !RL replaced axion pert with kg


    if (CP%Num_Nu_Massive > 0) then
       call MassiveNuVars(EV,ay,a,grho,gpres,dgrho,dgq, wnu_arr)
    end if

    if (CP%flat) then
       adotoa=sqrt(grho/3._dl)
       cothxor=1._dl/tau
    else
       adotoa=sqrt((grho+grhok)/3._dl)
       cothxor=1._dl/tanfunc(tau/CP%r)/CP%r
    end if

    ! if (w_lam /= -1 .and. w_Perturb) then
    !    clxq=ay(EV%w_ix)
    !    vq=ay(EV%w_ix+1)
    !    dgrho=dgrho + clxq*grhov_t
    !    dgq = dgq + vq*grhov_t*(1+w_lam)
    !end if

    if (EV%no_nu_multpoles) then
       !RSA approximation of arXiv:1104.2933, dropping opactity terms in the velocity
       !Approximate total density variables with just matter terms
       z=(0.5_dl*dgrho/k + etak)/adotoa
       dz= -adotoa*z - 0.5_dl*dgrho/k
       clxr=-4*dz/k
       qr=-4._dl/3*z
       pir=0
    else
       !  Massless neutrinos
       clxr=ay(EV%r_ix)
       qr  =ay(EV%r_ix+1)
       pir =ay(EV%r_ix+2)
    endif

    if (EV%no_phot_multpoles) then
       if (.not. EV%no_nu_multpoles) then
          z=(0.5_dl*dgrho/k + etak)/adotoa
          dz= -adotoa*z - 0.5_dl*dgrho/k
          clxg=-4*dz/k-4/k*opacity*(vb+z)
          qg=-4._dl/3*z
       else
          clxg=clxr-4/k*opacity*(vb+z)
          qg=qr
       end if
       pig=0
    else
       !  Photons
       clxg=ay(EV%g_ix)
       qg=ay(EV%g_ix+1)
       if (.not. EV%TightCoupling) pig=ay(EV%g_ix+2)
    end if

    !  8*pi*a*a*SUM[rho_i*clx_i] - radiation terms
    dgrho=dgrho + grhog_t*clxg+grhor_t*clxr

    !  8*pi*a*a*SUM[(rho_i+p_i)*v_i]
    dgq=dgq + grhog_t*qg+grhor_t*qr

    !  Photon mass density over baryon mass density
    photbar=grhog_t/grhob_t
    pb43=4._dl/3*photbar

    ayprime(1)=adotoa*a

    if (.not. is_cosmological_constant) then
       !ppf
       grhoT = grho - grhov_t
       vT= dgq/(grhoT+gpres)
       Gamma=ay(EV%w_ix)

       !sigma for ppf
       sigma = (etak + (dgrho + 3*adotoa/k*dgq)/2._dl/k)/EV%kf(1) - k*Gamma
       sigma = sigma/adotoa

       S_Gamma=grhov_t*(1+w_eff)*(vT+sigma)*k/adotoa/2._dl/k2
       ckH=c_Gamma_ppf*k/adotoa
       Gammadot=S_Gamma/(1+ckH*ckH)- Gamma -ckH*ckH*Gamma
       Gammadot=Gammadot*adotoa
       ayprime(EV%w_ix)=Gammadot

       if(ckH*ckH.gt.3.d1)then
          Gamma=0
          Gammadot=0.d0
          ayprime(EV%w_ix)=Gammadot
       endif

       Fa=1+3*(grhoT+gpres)/2._dl/k2/EV%kf(1)
       dgqe=S_Gamma - Gammadot/adotoa - Gamma
       dgqe=-dgqe/Fa*2._dl*k*adotoa + vT*grhov_t*(1+w_eff)
       dgrhoe=-2*k2*EV%kf(1)*Gamma-3/k*adotoa*dgqe
       dgrho=dgrho+dgrhoe
       dgq=dgq+dgqe

       EV%dgrho_e_ppf=dgrhoe
       EV%dgq_e_ppf=dgqe
       EV%dgrhoec_ppf=dgrhoe+3._dl*(1+w_eff)*grhov_t*adotoa*vT/k
       EV%dgqec_ppf=dgqe+(1+w_eff)*grhov_t*sigma
       EV%vTc_ppf=vT+sigma
    end if

    !  Get sigma (shear) and z from the constraints
    ! have to get z from eta for numerical stability
    z=(0.5_dl*dgrho/k + etak)/adotoa
    if (CP%flat) then
       !eta*k equation
       sigma=(z+1.5_dl*dgq/k2)
       ayprime(2)=0.5_dl*dgq
    else
       sigma=(z+1.5_dl*dgq/k2)/EV%Kf(1)
       ayprime(2)=0.5_dl*dgq + CP%curv*z
    end if

    !if (w_lam /= -1 .and. w_Perturb) then
    !
    !   ayprime(EV%w_ix)= -3*adotoa*(cs2_lam-w_lam)*(clxq+3*adotoa*(1+w_lam)*vq/k) &
    !       -(1+w_lam)*k*vq -(1+w_lam)*k*z
    !
    !   ayprime(EV%w_ix+1) = -adotoa*(1-3*cs2_lam)*vq + k*cs2_lam*clxq/(1+w_lam)
    !
    !end if
    !
    !  CDM equation of motion
    !write(*, *) 'Rayne, what is the first entry of dgrho, etak, adotoa, z in derivs?', dgrho, etak, adotoa, z
    clxcdot=-k*z
    ayprime(3)=clxcdot

    !  Baryon equation of motion.
    clxbdot=-k*(z+vb)
    ayprime(4)=clxbdot
    !write(120708, '(36e52.42,\)') tau, a, adotoa !RL

    !RL 081324 adding gpres
    gpres=(grhog_t+grhor_t)/3._dl+grhov_t*w_eff+gpres_ax
    ! --------------------------------------
    ! AXIONS
    ! Axion equation of motion updated 9/22/2012 to include inite.note
    !9/22 switch to (1+w_ax)v_ax as momentum variable for axions
    !in this epoch (this part of if statement w_ax)=0
    !9/22 DG update to evolve(1+w)*v_ax to avoid infinities when correct gauged EOM are used, see DG notes

    !RL adding KG variables
    !Also, the dv1 here is actually sqrt(4 pi G/3) dphi, dv2 here is sqrt(4 pi G/3H0^2) dphidot, following the v_vec internal version of axion_background. Derivation see equations_ppf tampering notes. Will change variable names after verifying that they work since this is confusing
    if (.not. EV%oscillation_started) then
       !write(*, *) 'Rayne, in derivs, before oscillation starts, k*z'
       !write(*, '(36e52.42,\)') k*z
       !The spline of v1_bg and v2_bg is already done at the beginning to obtain clxax_kg
       !RL testing by turning both the metric terms off for both KG and EFA
       !ayprime(EV%a_kg_ix) = dv2 * CP%H0_in_Mpc_inv
       !ayprime(EV%a_kg_ix+1) = -2 * adotoa * dv2 - k2*dv1/(CP%H0_in_Mpc_inv) - a2*(CP%m_ovH0**2.0_dl)*dv1*CP%H0_in_Mpc_inv !- k*z*v2_bg
       !The untouched KG---------------
       ayprime(EV%a_kg_ix) = dv2 * CP%H0_in_Mpc_inv/EV%renorm_c 
       ayprime(EV%a_kg_ix+1) = -2 * adotoa * dv2 - k2*dv1*EV%renorm_c/(CP%H0_in_Mpc_inv) &
            &- a2*(CP%m_ovH0**2.0_dl)*dv1*EV%renorm_c*CP%H0_in_Mpc_inv - k*z*v2_bg !RL 050324
    else !Now past the oscillation phase, use EFA
       !w_ax_p1 is defined at the start of the subroutine
       cad2 = w_ax*((1._dl+ gpres/grho)/w_ax_p1 +1._dl) !RL 081324 gpres is "changed" again in a later chunk but that's ppf and doens't matter
       kamnorm = k2/(a2*((CP%m_ovH0*CP%H0_in_Mpc_inv)**2.0_dl))
     if (kamnorm .lt. 1.e-14_dl) then 
          !write(*, *) 'Rayne, machine precision', kamnorm
          !RL dealing with machine precision issue - Taylor expand to leading order
        csquared_ax = kamnorm/4.0_dl + 5.0_dl*(adotoa**2.0_dl)/(4.0_dl*(k2/kamnorm))
        !csquared_ax = kamnorm/4.0_dl + 3.0_dl*(adotoa**2.0_dl)/(2.0_dl*(k2/kamnorm))
      else
         csquared_ax = (sqrt(1.0_dl + kamnorm) - 1.0_dl)**2.0_dl/(kamnorm) + 5.0_dl*(adotoa**2.0_dl)/(4.0_dl*(k2/kamnorm))
         !csquared_ax = (sqrt(1.0_dl + kamnorm) - 1.0_dl)**2.0_dl/(kamnorm) + 3.0_dl*(adotoa**2.0_dl)/(2.0_dl*(k2/kamnorm))
!           csquared_ax = (sqrt(1.0_dl + kamnorm) - 1.0_dl)**2.0_dl/(kamnorm) + 1.1_dl*(adotoa**2.0_dl)/((k2/kamnorm))
     end if
       !RL replaced w_ax with w_ax_p1 - 1.0_dl, 1.0_dl-w_ax with 2.0_dl-w_axp1, 1.0_dl+w_ax with w_ax_p1; all the clxax, etc. below are from the KG     
       !clxax_kg_dot = -w_ax_p1*(thetaax_kg+k*z) - 3.0_dl*(csquared_ax-w_ax)*adotoa*clxax_kg - 9.0_dl*w_ax_p1*(csquared_ax-cad2)*(adotoa**2.0_dl)*thetaax_kg/k2
       !thetaax_kg_dot= -(1.0_dl-3.0_dl*csquared_ax)*adotoa*thetaax_kg + csquared_ax*k2*clxax_kg/w_ax_p1
       !---------------------------
       !Rl testing to put the correction to w_ax and cad2 out of the pert EoM
       !w_ax = 3.0_dl*(adotoa**2.0_dl)/(2.0_dl*(k2/kamnorm))
       !w_ax_p1 = 1.0_dl + w_ax
       !cad2 = w_ax*((w_ax + 7.0_dl/3.0_dl)/w_ax_p1)
       !write(*, *) 'Rayne, cad2', cad2
       !clxax_kg_dot = -k*(u_ax_kg + z)-3.0_dl*(csquared_ax)*adotoa*clxax_kg-9.0_dl*(csquared_ax)*(adotoa**2.0_dl)*u_ax_kg/k
       !u_ax_kg_dot=-adotoa*u_ax_kg+3.0_dl*csquared_ax*adotoa*u_ax_kg+k*csquared_ax*clxax_kg+3.0_dl*(w_ax)*adotoa*u_ax_kg
       !RL testing while turning the metric off (kz) for both KG and EFA
       !clxax_kg_dot = -k*(u_ax_kg )-3.0_dl*(csquared_ax-w_ax)*adotoa*clxax_kg-9.0_dl*(csquared_ax-cad2)*(adotoa**2.0_dl)*u_ax_kg/k
       !u_ax_kg_dot=-adotoa*u_ax_kg+3.0_dl*csquared_ax*adotoa*u_ax_kg+k*csquared_ax*clxax_kg+3.0_dl*(w_ax-cad2)*adotoa*u_ax_kg
       !--------------------------- untouched EFA EoM
     clxax_kg_dot = -k*(u_ax_kg + z*w_ax_p1)-3.0_dl*(csquared_ax-w_ax)*adotoa*clxax_kg-&
          &9.0_dl*(csquared_ax-cad2)*(adotoa**2.0_dl)*u_ax_kg/k
     u_ax_kg_dot=-adotoa*u_ax_kg+3.0_dl*csquared_ax*adotoa*u_ax_kg+k*csquared_ax*clxax_kg+&
          &3.0_dl*(w_ax-cad2)*adotoa*u_ax_kg
       ayprime(EV%a_kg_ix) = clxax_kg_dot
       ayprime(EV%a_kg_ix+1) = u_ax_kg_dot
    end if
    ! --------------------------------------

    !Equations combined piecewise parameters     (DG GDM version - should only put clxax_kg in the metric, note the RHS equation)
    !if (.not. EV%oscillation_started) then

    !   csquared_ax_use=1.0d0 
    !cad2_use=dble(cad2) !axionCAMB original
    

    !   cad2_use = 1.0d0 + 2.0d0 * (CP%m_ovH0**2.0d0) * a2 * v1_bg/(3 * (adotoa/CP%H0_in_Mpc_inv) * v2_bg) !RL testing
    !write(*, *) 'Rayne, cad2, cad2_use, cad2_vs_-7/3, cad2_use_vs-7/3, their fractional difference', cad2, cad2_use, cad2/(-7.0_dl/3.0_dl) - 1.0_dl, cad2_use/(-7.0_dl/3.0_dl) - 1.0_dl, cad2_use/cad2 - 1.0_dl
    !write(*, *) 'Rayne, adotoa_bg, adotoa, their fractional difference', adotoa_bg, adotoa, adotoa/adotoa_bg - 1.0_dl
    !write(110308, '(36e52.42,\)') tau, a, cad2, cad2_use, adotoa, adotoa_bg

    !RL for debugging purposes used adotoa_bg instead of adotoa
    !   w_axp1_use = (2.0d0*(v2_bg**2.0d0)/a2)/((v2_bg**2.0d0)/a2 + (CP%m_ovH0*v1_bg)**2.0d0)
    !write(*, *) 'Rayne, tau, w_axp1, w_axp1_use, their fractional difference', tau, w_axp1, w_axp1_use, w_axp1/w_axp1_use - 1.0_dl
    !RL replaced w_ax with w_axp1 - 1.0d0, 1.0d0-w_ax with 2.0d0-w_axp1, 1.0d0+w_ax with w_axp1
    !RL using w_axp1use calculated from the field variables
    !   clxaxdot=-3.0d0*adotoa*(csquared_ax_use)*(2.0d0-w_axp1_use)*clxax-(k*v_ax+k*z*(w_axp1_use))&
    !        &-9.0d0*((adotoa)**2.0d0)*csquared_ax_use*v_ax*(1.0d0-cad2_use)/k
    !   v_axdot=(-adotoa*v_ax+k*csquared_ax_use*clxax)+v_ax*3.0d0*adotoa*(w_axp1_use - 1.0d0-cad2_use)&
    !        &+csquared_ax_use*3.0d0*adotoa*v_ax
    !else
    !write(*, *) 'a = ', a, 'after oscillation' !RL trying to see when the functions are called
    !   csquared_ax_use=csquared_ax
    !write(*, *) 'csquared_ax_use is zero exactly?', csquared_ax_use == 0._dl, csquared_ax_use
    !   cad2_use=0.0d0
    !   w_ax=0.0d0
    !   w_axp1 = 1.0d0
    !RL replaced w_ax with w_axp1 - 1.0d0, 1.0d0-w_ax with 2.0d0-w_axp1, 1.0d0+w_ax with w_axp1
    !   clxaxdot=-3.0d0*adotoa*(csquared_ax_use)*(2.0d0-w_axp1)*clxax-(k*v_ax+k*z*(w_axp1))&
    !        &-9.0d0*(adotoa**2.0d0)*csquared_ax_use*v_ax*(1.0d0-cad2_use)/k
    !   v_axdot=(-adotoa*v_ax+k*csquared_ax_use*clxax)+v_ax*3.0d0*adotoa*(w_axp1 - 1.0d0-cad2_use) &  	
    !        &+csquared_ax_use*3.0d0*adotoa*v_ax
    !endif
    !write(082401, '(36e52.42,\)') k, a, a_osc, tau, CP%tau_osc, w_ax, 1.0_dl + w_ax, csquared_ax_use, cad2_use 

    !axionCAMB's GDM variables
    clxaxdot = 0. !RL: already useless, but didn't get rid of due to isocurvature
    v_axdot = 0.
    ayprime(EV%a_ix)=clxaxdot
    ayprime(EV%a_ix+1)=v_axdot



    !  Photon equation of motion
    clxgdot=-k*(4._dl/3._dl*z+qg)

    ! old comment:Small k: potential problem with stability, using full equations earlier is NOT more accurate in general
    ! Easy to see instability in k \sim 1e-3 by tracking evolution of vb

    !  Use explicit equation for vb if appropriate

    if (EV%TightCoupling) then
       !  ddota/a
       !gpres=gpres + grhov_t*w_eff !RL 081324 moving this to the outside since I need gpres
       adotdota=(adotoa*adotoa-gpres)/2

       pig = 32._dl/45/opacity*k*(sigma+vb)

       !  First-order approximation to baryon-photon splip
       slip = - (2*adotoa/(1+pb43) + dopacity/opacity)* (vb-3._dl/4*qg) &
            +(-adotdota*vb-k/2*adotoa*clxg +k*(cs2*clxbdot-clxgdot/4))/(opacity*(1+pb43))

       if (second_order_tightcoupling) then
          ! by Francis-Yan Cyr-Racine simplified (inconsistently) by AL assuming flat
          !AL: First order slip seems to be fine here to 2e-4

          !  8*pi*G*a*a*SUM[rho_i*sigma_i]
          dgs = grhog_t*pig+grhor_t*pir

          ! Define shear derivative to first order
          sigmadot = -2*adotoa*sigma-dgs/k+etak

          !Once know slip, recompute qgdot, pig, pigdot
          qgdot = k*(clxg/4._dl-pig/2._dl) +opacity*slip

          pig = 32._dl/45/opacity*k*(sigma+3._dl*qg/4._dl)*(1+(dopacity*11._dl/6._dl/opacity**2)) &
               + (32._dl/45._dl/opacity**2)*k*(sigmadot+3._dl*qgdot/4._dl)*(-11._dl/6._dl)

          pigdot = -(32._dl/45._dl)*(dopacity/opacity**2)*k*(sigma+3._dl*qg/4._dl)*(1 + &
               dopacity*11._dl/6._dl/opacity**2 ) &
               + (32._dl/45._dl/opacity)*k*(sigmadot+3._dl*qgdot/4._dl)*(1+(11._dl/6._dl) &
               *(dopacity/opacity**2))

          EV%pigdot = pigdot
       end if

       !  Use tight-coupling approximation for vb
       !  zeroth order approximation to vbdot + the pig term
       vbdot=(-adotoa*vb+cs2*k*clxb  &
            +k/4*pb43*(clxg-2*EV%Kf(1)*pig))/(1+pb43)

       vbdot=vbdot+pb43/(1+pb43)*slip

       EV%pig = pig
    else
       vbdot=-adotoa*vb+cs2*k*clxb-photbar*opacity*(4._dl/3*vb-qg)
    end if

    ayprime(5)=vbdot

    if (.not. EV%no_phot_multpoles) then
       !  Photon equations of motion
       ayprime(EV%g_ix)=clxgdot
       qgdot=4._dl/3*(-vbdot-adotoa*vb+cs2*k*clxb)/pb43 &
            +EV%denlk(1)*clxg-EV%denlk2(1)*pig
       ayprime(EV%g_ix+1)=qgdot

       !  Use explicit equations for photon moments if appropriate
       if (.not. EV%tightcoupling) then
          E2=ay(EV%polind+2)
          polter = pig/10+9._dl/15*E2 !2/15*(3/4 pig + 9/2 E2)
          ix= EV%g_ix+2
          if (EV%lmaxg>2) then
             pigdot=EV%denlk(2)*qg-EV%denlk2(2)*ay(ix+1)-opacity*(pig - polter) &
                  +8._dl/15._dl*k*sigma
             ayprime(ix)=pigdot
             do  l=3,EV%lmaxg-1
                ix=ix+1
                ayprime(ix)=(EV%denlk(l)*ay(ix-1)-EV%denlk2(l)*ay(ix+1))-opacity*ay(ix)
             end do
             ix=ix+1
             !  Truncate the photon moment expansion
             ayprime(ix)=k*ay(ix-1)-(EV%lmaxg+1)*cothxor*ay(ix) -opacity*ay(ix)
          else !closed case
             pigdot=EV%denlk(2)*qg-opacity*(pig - polter) +8._dl/15._dl*k*sigma
             ayprime(ix)=pigdot
          endif
          !  Polarization
          !l=2
          ix=EV%polind+2
          if (EV%lmaxgpol>2) then
             ayprime(ix) = -opacity*(ay(ix) - polter) - k/3._dl*ay(ix+1)
             do l=3,EV%lmaxgpol-1
                ix=ix+1
                ayprime(ix)=-opacity*ay(ix) + (EV%denlk(l)*ay(ix-1)-EV%polfack(l)*ay(ix+1))
             end do
             ix=ix+1
             !truncate
             ayprime(ix)=-opacity*ay(ix) + &
                  k*EV%poltruncfac*ay(ix-1)-(EV%lmaxgpol+3)*cothxor*ay(ix)
          else !closed case
             ayprime(ix) = -opacity*(ay(ix) - polter)
          endif
       end if
    end if

    if (.not. EV%no_nu_multpoles) then
       !  Massless neutrino equations of motion.
       clxrdot=-k*(4._dl/3._dl*z+qr)
       ayprime(EV%r_ix)=clxrdot
       qrdot=EV%denlk(1)*clxr-EV%denlk2(1)*pir
       ayprime(EV%r_ix+1)=qrdot
       if (EV%high_ktau_neutrino_approx) then
          !ufa approximation for k*tau>>1, more accurate when there are reflections from lmax
          !Method from arXiv:1104.2933
          !                if (.not. EV%TightCoupling) then
          !                 gpres=gpres+ (grhog_t+grhor_t)/3 +grhov_t*w_lam
          !                 adotdota=(adotoa*adotoa-gpres)/2
          !                end if
          !                ddz=(2*adotoa**2 - adotdota)*z  &
          !                  + adotoa/(2*k)*( 6*(grhog_t*clxg+grhor_t*clxr) + 2*(grhoc_t*clxc+grhob_t*clxb) ) &
          !                   - 1._dl/(2*k)*( 2*(grhog_t*clxgdot+grhor_t*clxrdot) + grhoc_t*clxcdot + grhob_t*clxbdot )
          !                dz= -adotoa*z - 0.5_dl*dgrho/k
          !                pirdot= -3*pir*cothxor + k*(qr+4._dl/3*z)
          pirdot= -3*pir*cothxor - clxrdot
          ayprime(EV%r_ix+2)=pirdot

          !                pirdot=k*(0.4_dl*qr-0.6_dl*ay(EV%lmaxg+10)+8._dl/15._dl*sigma)
          !                ayprime(EV%lmaxg+9)=pirdot
          !                ayprime(3+EV%lmaxg+7)=k*ay(3+EV%lmaxg+6)- &
          !                                      (3+1)*cothxor*ay(3+EV%lmaxg+7)
          !               ayprime(3+EV%lmaxg+7+1:EV%lmaxnr+EV%lmaxg+7)=0
       else
          ix=EV%r_ix+2
          if (EV%lmaxnr>2) then
             pirdot=EV%denlk(2)*qr- EV%denlk2(2)*ay(ix+1)+8._dl/15._dl*k*sigma
             ayprime(ix)=pirdot
             do l=3,EV%lmaxnr-1
                ix=ix+1
                ayprime(ix)=(EV%denlk(l)*ay(ix-1) - EV%denlk2(l)*ay(ix+1))
             end do
             !  Truncate the neutrino expansion
             ix=ix+1
             ayprime(ix)=k*ay(ix-1)- (EV%lmaxnr+1)*cothxor*ay(ix)
          else
             pirdot=EV%denlk(2)*qr +8._dl/15._dl*k*sigma
             ayprime(ix)=pirdot
          end if
       end if
    end if ! no_nu_multpoles

    !  Massive neutrino equations of motion.
    if (CP%Num_Nu_massive == 0) return

    do nu_i = 1, CP%Nu_mass_eigenstates
       if (EV%MassiveNuApprox(nu_i)) then
          !Now EV%iq0 = clx, EV%iq0+1 = clxp, EV%iq0+2 = G_1, EV%iq0+3=G_2=pinu
          !see astro-ph/0203507
          G11_t=EV%G11(nu_i)/a/a2
          G30_t=EV%G30(nu_i)/a/a2
          off_ix = EV%nu_ix(nu_i)
          w=wnu_arr(nu_i)
          ayprime(off_ix)=-k*z*(w+1) + 3*adotoa*(w*ay(off_ix) - ay(off_ix+1))-k*ay(off_ix+2)
          ayprime(off_ix+1)=(3*w-2)*adotoa*ay(off_ix+1) - 5._dl/3*k*z*w - k/3*G11_t
          ayprime(off_ix+2)=(3*w-1)*adotoa*ay(off_ix+2) - k*(2._dl/3*EV%Kf(1)*ay(off_ix+3)-ay(off_ix+1))
          ayprime(off_ix+3)=(3*w-2)*adotoa*ay(off_ix+3) + 2*w*k*sigma - k/5*(3*EV%Kf(2)*G30_t-2*G11_t)
       else
          ind=EV%nu_ix(nu_i)
          do i=1,EV%nq(nu_i)
             q=nu_q(i)
             aq=a*nu_masses(nu_i)/q
             v=1._dl/sqrt(1._dl+aq*aq)

             ayprime(ind)=-k*(4._dl/3._dl*z + v*ay(ind+1))
             ind=ind+1
             ayprime(ind)=v*(EV%denlk(1)*ay(ind-1)-EV%denlk2(1)*ay(ind+1))
             ind=ind+1
             if (EV%lmaxnu_tau(nu_i)==2) then
                ayprime(ind)=-ayprime(ind-2) -3*cothxor*ay(ind)
             else
                ayprime(ind)=v*(EV%denlk(2)*ay(ind-1)-EV%denlk2(2)*ay(ind+1)) &
                     +k*8._dl/15._dl*sigma
                do l=3,EV%lmaxnu_tau(nu_i)-1
                   ind=ind+1
                   ayprime(ind)=v*(EV%denlk(l)*ay(ind-1)-EV%denlk2(l)*ay(ind+1))
                end do
                !  Truncate moment expansion.
                ind = ind+1
                ayprime(ind)=k*v*ay(ind-1)-(EV%lmaxnu_tau(nu_i)+1)*cothxor*ay(ind)
             end if
             ind = ind+1
          end do
       end if
    end do

    if (EV%has_nu_relativistic) then
       ind=EV%nu_pert_ix
       ayprime(ind)=+k*a2*qr -k*ay(ind+1)
       ind2= EV%r_ix
       do l=1,EV%lmaxnu_pert-1
          ind=ind+1
          ind2=ind2+1
          ayprime(ind)= -a2*(EV%denlk(l)*ay(ind2-1)-EV%denlk2(l)*ay(ind2+1)) &
               +   (EV%denlk(l)*ay(ind-1)-EV%denlk2(l)*ay(ind+1))
       end do
       ind=ind+1
       ind2=ind2+1
       ayprime(ind)= k*(ay(ind-1) -a2*ay(ind2-1)) -(EV%lmaxnu_pert+1)*cothxor*ay(ind)
    end if

    !RL testing the quasistatic evolution (have to be after z is defined)
    !dv1_quasitest = -z*v2_bg*(CP%H0_in_Mpc_inv)/k
    !clxax_quasitest = ((CP%m_ovH0**2.0d0)*v1_bg*dv1_quasitest)/((v2_bg)**2.0d0/a2+(CP%m_ovH0*v1_bg)**2.0d0)
    !write(03021923, '(36e52.42,\)') k, a, CP%a_osc, tau, CP%tau_osc, adotoa, CP%m_ovH0*CP%H0_in_Mpc_inv/(adotoa/a), dorp, w_ax_p1, cad2, clxax_kg, thetaax_kg, clxax_quasitest, dv1, dv1_quasitest, z, v2_bg, dgrho, clxg

    !!if (tau .eq. CP%tau_osc .and. .not. EV%oscillation_started) then
    !!   write(*, *) 'Rayne, in derivs on the KG side, tau = tauosc, ayprime(2), ayprime(3),adotoa/H0, adotoa*sigma/k', ayprime(2), ayprime(3), adotoa/CP%H0_in_Mpc_inv, adotoa*sigma/k
    !!end if
    
    !!if (tau .eq. CP%tau_osc .and. EV%oscillation_started) then
    !!   write(*, *) 'Rayne, in derivs on the EFA side, tau = tauosc, ayprime(2), ayprime(3), adotoa/H0, adotoa*sigma/k', ayprime(2), ayprime(3), adotoa/CP%H0_in_Mpc_inv, adotoa*sigma/k
    !!end if
    !if (.not. EV%output_done) then     
    !   write(*, '(A, 36E52.42)') 'Rayne, in derivs, etak', etak
    !   EV%output_done = .true.
    !end if
    
    
  end subroutine derivs



  subroutine derivsv(EV,n,tau,yv,yvprime)
    !  Evaluate the time derivatives of the vector perturbations, flat case
    use ThermoData
    use MassiveNu
    implicit none
    type(EvolutionVars) EV
    integer n,l
    real(dl), target ::  yv(n),yvprime(n)
    real(dl) ep,tau,grho,rhopi,cs2,opacity,gpres
    logical finished_tightcoupling
    real(dl), dimension(:),pointer :: neut,neutprime,E,B,Eprime,Bprime
    real(dl)  grhob_t,grhor_t,grhoc_t,grhog_t,grhov_t,polter
    real(dl) sigma, qg,pig, qr, vb, rhoq, vbdot, photbar, pb43
    real(dl) k,k2,a,a2, adotdota
    real(dl) pir,adotoa
    ! ppf
    real(dl) w_eff
    ! Axions
    real(dl) v1_bg, v2_bg, grhoax_kg !RL added background variables
    real(dl) grhoax_t,gpres_ax
    real(dl) dorp
    real(dl) gr,w_ax, w_ax_p1, wcorr_coeff !RL added w_ax_p1

    stop 'ppf not implemented for vectors'

    k2=EV%k2_buf
    k=EV%k_buf

    !E and B start at l=2. Set up pointers accordingly to fill in y arrays
    E => yv(EV%lmaxv+3:)
    Eprime=> yvprime(EV%lmaxv+3:)
    B => E(EV%lmaxpolv:)
    Bprime => Eprime(EV%lmaxpolv:)
    neutprime => Bprime(EV%lmaxpolv+1:)
    neut => B(EV%lmaxpolv+1:)

    a=yv(1)

    sigma=yv(2)

    a2=a*a

    !  Get sound speed and opacity, and see if should use tight-coupling

    call thermo(tau,cs2,opacity)
    if (k > 0.06_dl*epsw) then
       ep=ep0
    else
       ep=0.2_dl*ep0
    end if

    finished_tightcoupling = &
         ((k/opacity > ep).or.(1._dl/(opacity*tau) > ep .and. k/opacity > 1d-4))


    ! Compute expansion rate from: grho=8*pi*rho*a**2
    ! Also calculate gpres: 8*pi*p*a**2
    grhob_t=grhob/a
    grhoc_t=grhoc/a
    grhor_t=grhornomass/a2
    grhog_t=grhog/a2
    grhov_t=grhov*a**(-1-3*w_lam)

    ! Axion variables, updated 9/19 to include dan's spline stuff

    !if (tau .le. CP%tau_osc) then
    !if (tau .lt. CP%tau_osc) then
    !if (a .lt. CP%a_osc) then
    !if (a .le. CP%a_osc) then !RL
    if (.not. EV%oscillation_started) then !RL
       !call spline_out(loga_table,rhoaxh2ovrhom_logtable,rhoaxh2ovrhom_logtable_buff,ntable,dlog10(a),gr)
       !write(*,*) a, gr
       !compute log10 of density use of interpolation
       !print*,a,grhoc,grhom*(1.0d1**(dble(gr)))/((grhoc+grhob)/(a**3.0d0))
       !if (gr .eq. gr) then
       !delog it and multiply by physical units
       !   dorp=grhom*(1.0d1**(dble(gr)))
       !RL
       
       call spline_out(loga_table,phinorm_table,phinorm_table_ddlga,ntable,dlog10(a),v1_bg)
       call spline_out(loga_table,phidotnorm_table,phidotnorm_table_ddlga,ntable,dlog10(a),v2_bg)
       if (v1_bg .eq. v1_bg .and. v2_bg .eq. v2_bg) then
          !Get the grhoax from field variables 
          grhoax_kg = (v2_bg)**2.0_dl/a2+(CP%m_ovH0*v1_bg)**2.0_dl
          dorp = grhom*grhoax_kg/(CP%H0**2.0d0/1.0d4)
       else
          dorp=0.0d0
       endif
       !obtain 1+w 
       w_ax_p1 = (2.0_dl*(v2_bg**2.0_dl)/a2)/((v2_bg**2.0_dl)/a2 + (CP%m_ovH0*v1_bg)**2.0_dl)
       w_ax = w_ax_p1 - 1.0_dl
    else
       !RL adding w correction to the background 
       !wcorr_coeff = CP%ah_osc*CP%a_osc/(CP%m_ovH0*(CP%H0/100.0d0))
       wcorr_coeff = CP%ahosc_ETA*CP%a_osc/(CP%m_ovH0*(CP%H0/100.0d0)) !RL082924
       w_ax = ((wcorr_coeff/a2)**2.0d0)*CP%wEFA_c
       w_ax_p1 = 1.0_dl + w_ax
       dorp=grhom*CP%rhorefp_ovh2*((CP%a_osc/a)**3.0d0)*dexp((wcorr_coeff**2.0d0)*3.0d0*&
            &CP%wEFA_c*(1.0d0/(a2**2.0d0) - 1.0d0/(CP%a_osc**4.0d0))/4.0d0)
       !dorp=grhom*CP%rhorefp_hsq*((CP%a_osc/a)**3.0d0)

    endif
    grhoax_t=dorp*(a**2.0d0)

    gpres_ax=w_ax_p1*grhoax_t 

    grho=grhob_t+grhoc_t+grhor_t+grhog_t+grhov_t+grhoax_t
    gpres=(grhog_t+grhor_t)/3._dl+grhov_t*w_eff+gpres_ax
    adotoa=sqrt(grho/3._dl)
    adotdota=(adotoa*adotoa-gpres)/2

    photbar=grhog_t/grhob_t
    pb43=4._dl/3*photbar

    yvprime(1)=adotoa*a

    vb = yv(3)
    qg = yv(4)
    qr = neut(1)

    !  8*pi*a*a*SUM[(rho_i+p_i)*v_i]
    rhoq=grhob_t*vb+grhog_t*qg+grhor_t*qr
    !  sigma = 2*rhoq/k**2
    !for non-large k this expression for sigma is unstable at early times
    !so propagate sigma equation separately (near total cancellation in rhoq)
    ! print *,yv(2),2*rhoq/k**2

    if (finished_tightcoupling) then
       !  Use explicit equations:

       pig = yv(5)

       polter = 0.1_dl*pig + 9._dl/15._dl*E(2)

       vbdot = -adotoa*vb-photbar*opacity*(4._dl/3*vb-qg) - 0.5_dl*k*photbar*Magnetic

       !  Equation for the photon heat flux stress

       yvprime(4)=-0.5_dl*k*pig + opacity*(4._dl/3*vb-qg)

       !  Equation for the photon anisotropic stress
       yvprime(5)=k*(2._dl/5*qg -8/15._dl*yv(6))+8._dl/15._dl*k*sigma  &
            -opacity*(pig - polter)
       ! And for the moments
       do  l=3,EV%lmaxv-1
          yvprime(l+3)=k*denl(l)*l*(yv(l+2)-   &
               vecfac(l)*yv(l+4))-opacity*yv(l+3)
       end do
       !  Truncate the hierarchy
       yvprime(EV%lmaxv+3)=k*EV%lmaxv/(EV%lmaxv-1._dl)*yv(EV%lmaxv+2)- &
            (EV%lmaxv+2._dl)*yv(EV%lmaxv+3)/tau-opacity*yv(EV%lmaxv+3)

       !E equations

       Eprime(2) = - opacity*(E(2) - polter) + k*(1/3._dl*B(2) - 8._dl/27._dl*E(3))
       do l=3,EV%lmaxpolv-1
          Eprime(l) =-opacity*E(l) + k*(denl(l)*(l*E(l-1) - &
               vecfacpol(l)*E(l+1)) + 2._dl/(l*(l+1))*B(l))
       end do
       !truncate
       Eprime(EV%lmaxpolv)=0._dl

       !B-bar equations

       do l=2,EV%lmaxpolv-1
          Bprime(l) =-opacity*B(l) + k*(denl(l)*(l*B(l-1) - &
               vecfacpol(l)*B(l+1)) - 2._dl/(l*(l+1))*E(l))
       end do
       !truncate
       Bprime(EV%lmaxpolv)=0._dl
    else
       !Tight coupling expansion results

       pig = 32._dl/45._dl*k/opacity*(vb + sigma)

       EV%pig = pig

       vbdot=(-adotoa*vb  -3._dl/8*pb43*k*Magnetic  -3._dl/8*k*pb43*pig &
            - pb43/(1+pb43)/opacity*(0.75_dl*k*adotoa*pb43**2/(pb43+1)*Magnetic + vb*&
            ( 2*pb43*adotoa**2/(1+pb43) + adotdota)) &
            )/(1+pb43)

       !  Equation for the photon heat flux
       ! Get drag from vbdot expression
       yvprime(4)=-0.5_dl*k*pig - &
            (vbdot+adotoa*vb)/photbar - 0.5_dl*k*Magnetic

       !  Set the derivatives to zero
       yvprime(5:n)=0._dl
       yv(5)=pig
       E(2)=  pig/4
    endif

    yvprime(3) = vbdot

    !  Neutrino equations:

    !  Massless neutrino anisotropic stress
    pir=neut(2)
    neutprime(1)= -0.5_dl*k*pir
    neutprime(2)=2._dl/5*k*qr -8._dl/15._dl*k*neut(3)+ 8._dl/15._dl*k*sigma
    !  And for the moments
    do  l=3,EV%lmaxnrv-1
       neutprime(l)=k*denl(l)*l*(neut(l-1)- vecfac(l)*neut(l+1))
    end do

    !  Truncate the hierarchy
    neutprime(EV%lmaxnrv)=k*EV%lmaxnrv/(EV%lmaxnrv-1._dl)*neut(EV%lmaxnrv-1)-  &
         (EV%lmaxnrv+2._dl)*neut(EV%lmaxnrv)/tau


    !  Get the propagation equation for the shear

    rhopi=grhog_t*pig+grhor_t*pir+ grhog_t*Magnetic

    yvprime(2)=-2*adotoa*sigma -rhopi/k

  end subroutine derivsv



  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine derivst(EV,n,tau,ayt,aytprime)
    !  Evaluate the time derivatives of the tensor perturbations.
    use ThermoData
    use MassiveNu
    implicit none
    type(EvolutionVars) EV
    integer n,l,i,ind, nu_i
    real(dl), target ::  ayt(n),aytprime(n)
    real(dl) tau,grho,rhopi,cs2,opacity,pirdt
    real(dl), dimension(:),pointer :: neut,neutprime,E,B,Eprime,Bprime
    real(dl) q,aq,v
    real(dl)  grhob_t,grhor_t,grhoc_t,grhog_t,grhov_t,polter
    real(dl) Hchi,pinu, pig
    real(dl) k,k2,a,a2
    real(dl) pir, adotoa, rhonu, shear

    real(dl) cothxor
    ! Axions       
    real(dl) v1_bg, v2_bg, grhoax_kg, grhoax_t !RL added background field & grhoax
    real(dl) dorp, wcorr_coeff !RL added w correction
    real(dl) gr


    k2=EV%k2_buf
    k= EV%k_buf

    a=ayt(1)

    Hchi=ayt(2)

    shear=ayt(3)

    a2=a*a

    ! Compute expansion rate from: grho=8*pi*rho*a**2
    ! Also calculate gpres: 8*pi*p*a**2
    grhob_t=grhob/a
    grhoc_t=grhoc/a
    grhor_t=grhornomass/a2
    grhog_t=grhog/a2
    if (is_cosmological_constant) then
       grhov_t=grhov*a2
    else
       grhov_t=grho_de(a)/a2
    end if
    ! Axion variables, updated 9/19 to include dan's spline stuff
    !if (tau .le. CP%tau_osc) then
    !if (tau .lt. CP%tau_osc) then
    !if (a .lt. CP%a_osc) then
    !if (a .le. CP%a_osc) then !RL
    !RL 071026: EV%oscillation_started is now set for tensors by the tau_osc
    !switch event in GaugeInterface_EvolveTens (it used to be scalar-only,
    !which left this branch on the KG tables past a_osc and crashed dverk).
    if (.not. EV%oscillation_started) then !RL
       !call spline_out(loga_table,rhoaxh2ovrhom_logtable,rhoaxh2ovrhom_logtable_buff,ntable,dlog10(a),gr)
       !write(*,*) a, gr
       !compute log10 of density use of interpolation
       !print*,a,grhoc,grhom*(1.0d1**(dble(gr)))/((grhoc+grhob)/(a**3.0d0))
       !if (gr .eq. gr) then
       !delog it and multiply by physical units
       !   dorp=grhom*(1.0d1**(dble(gr)))
       
       call spline_out(loga_table,phinorm_table,phinorm_table_ddlga,ntable,dlog10(a),v1_bg)
       call spline_out(loga_table,phidotnorm_table,phidotnorm_table_ddlga,ntable,dlog10(a),v2_bg)
       if (v1_bg .eq. v1_bg .and. v2_bg .eq. v2_bg) then
          !Get the grhoax from field variables 
          grhoax_kg = (v2_bg)**2.0_dl/a2+(CP%m_ovH0*v1_bg)**2.0_dl
          dorp = grhom*grhoax_kg/(CP%H0**2.0d0/1.0d4)
       else
          dorp=0.0d0
       endif
    else
       !RL adding w correction to the background 
       !wcorr_coeff = CP%ah_osc*CP%a_osc/(CP%m_ovH0*(CP%H0/100.0d0))
       wcorr_coeff = CP%ahosc_ETA*CP%a_osc/(CP%m_ovH0*(CP%H0/100.0d0)) !RL082924
       !w_ax_p1 = 1.0_dl + 3.0d0*((wcorr_coeff/a2)**2.0d0)/2.0d0 !Not declared, not needed
       dorp=grhom*CP%rhorefp_ovh2*((CP%a_osc/a)**3.0d0)*dexp((wcorr_coeff**2.0d0)*3.0d0*&
            &CP%wEFA_c*(1.0d0/(a2**2.0d0) - 1.0d0/(CP%a_osc**4.0d0))/4.0d0)
       !dorp=grhom*CP%rhorefp_hsq*((CP%a_osc/a)**3.0d0)
    endif

    grhoax_t=dorp*(a**2.0d0)


    grho=grhob_t+grhoc_t+grhor_t+grhog_t+grhov_t+grhoax_t

    !Do massive neutrinos
    if (CP%Num_Nu_Massive >0) then
       do nu_i=1,CP%Nu_mass_eigenstates
          call Nu_rho(a*nu_masses(nu_i),rhonu)
          grho=grho+grhormass(nu_i)*rhonu/a2
       end do
    end if

    if (CP%flat) then
       cothxor=1._dl/tau
       adotoa=sqrt(grho/3._dl)
    else
       cothxor=1._dl/tanfunc(tau/CP%r)/CP%r
       adotoa=sqrt((grho+grhok)/3._dl)
    end if

    aytprime(1)=adotoa*a

    call thermo(tau,cs2,opacity)

    if (.not. EV%TensTightCoupling) then
       !  Don't use tight coupling approx - use explicit equations:
       !  Equation for the photon anisotropic stress


       !E and B start at l=2. Set up pointers accordingly to fill in ayt arrays
       E => ayt(EV%E_ix+1:)
       B => ayt(EV%B_ix+1:)
       Eprime=> aytprime(EV%E_ix+1:)
       Bprime => aytprime(EV%B_ix+1:)

       ind = EV%g_ix+2

       !  Photon anisotropic stress
       pig=ayt(ind)
       polter = 0.1_dl*pig + 9._dl/15._dl*E(2)

       if (EV%lmaxt > 2) then
          aytprime(ind)=-EV%denlkt(2,2)*ayt(ind+1)+k*8._dl/15._dl*shear  &
               -opacity*(pig - polter)

          do l=3, EV%lmaxt -1
             ind = ind+1
             aytprime(ind)=EV%denlkt(1,L)*ayt(ind-1)-EV%denlkt(2,L)*ayt(ind+1)-opacity*ayt(ind)
          end do

          !Truncate the hierarchy
          ind=ind+1
          aytprime(ind)=k*EV%lmaxt/(EV%lmaxt-2._dl)*ayt(ind-1)- &
               (EV%lmaxt+3._dl)*cothxor*ayt(ind)-opacity*ayt(ind)

          !E and B-bar equations

          Eprime(2) = - opacity*(E(2) - polter) + EV%denlkt(4,2)*B(2) - &
               EV%denlkt(3,2)*E(3)

          do l=3, EV%lmaxpolt-1
             Eprime(l) =(EV%denlkt(1,L)*E(l-1)-EV%denlkt(3,L)*E(l+1) + EV%denlkt(4,L)*B(l)) &
                  -opacity*E(l)
          end do
          l= EV%lmaxpolt
          !truncate: difficult, but setting l+1 to zero seems to work OK
          Eprime(l) = (EV%denlkt(1,L)*E(l-1) + EV%denlkt(4,L)*B(l)) -opacity*E(l)

          Bprime(2) =-EV%denlkt(3,2)*B(3) - EV%denlkt(4,2)*E(2)  -opacity*B(2)
          do l=3, EV%lmaxpolt-1
             Bprime(l) =(EV%denlkt(1,L)*B(l-1) -EV%denlkt(3,L)*B(l+1) - EV%denlkt(4,L)*E(l)) &
                  -opacity*B(l)
          end do
          l=EV%lmaxpolt
          !truncate
          Bprime(l) =(EV%denlkt(1,L)*B(l-1) - EV%denlkt(4,L)*E(l))  -opacity*B(l)

       else !lmax=2

          aytprime(ind)=k*8._dl/15._dl*shear-opacity*(pig - polter)
          Eprime(2) = - opacity*(E(2) - polter) + EV%denlkt(4,2)*B(2)
          Bprime(2) = - EV%denlkt(4,2)*E(2)  -opacity*B(2)
       end if

    else  !Tight coupling
       pig = 32._dl/45._dl*k/opacity*shear
    endif

    rhopi=grhog_t*pig


    !  Neutrino equations:
    !  Anisotropic stress
    if (DoTensorNeutrinos) then
       neutprime => aytprime(EV%r_ix+1:)
       neut => ayt(EV%r_ix+1:)

       !  Massless neutrino anisotropic stress
       pir=neut(2)

       rhopi=rhopi+grhor_t*pir

       if (EV%lmaxnrt>2) then
          pirdt=-EV%denlkt(2,2)*neut(3) + 8._dl/15._dl*k*shear
          neutprime(2)=pirdt
          !  And for the moments
          do  l=3, EV%lmaxnrt-1
             neutprime(l)= EV%denlkt(1,L)*neut(l-1) -EV%denlkt(2,L)*neut(l+1)
          end do

          !  Truncate the hierarchy
          neutprime(EV%lmaxnrt)=k*EV%lmaxnrt/(EV%lmaxnrt-2._dl)*neut(EV%lmaxnrt-1)-  &
               (EV%lmaxnrt+3._dl)*cothxor*neut(EV%lmaxnrt)
       else
          pirdt= 8._dl/15._dl*k*shear
          neutprime(2)=pirdt
       end if

       !  Massive neutrino equations of motion and contributions to anisotropic stress.
       if (CP%Num_Nu_massive > 0) then
          do nu_i=1,CP%Nu_mass_eigenstates
             if (.not. EV%EvolveTensorMassiveNu(nu_i)) then
                rhopi=rhopi+ grhormass(nu_i)/a2*pir !- good approx, note no rhonu weighting
             else
                ind=EV%nu_ix(nu_i)+2

                pinu= Nu_pi(EV, ayt(ind),a, nu_i)
                rhopi=rhopi+ grhormass(nu_i)/a2*pinu

                do i=1,nqmax
                   q=nu_q(i)
                   aq=a*nu_masses(nu_i)/q
                   v=1._dl/sqrt(1._dl+aq*aq)
                   if (EV%lmaxnut>2) then
                      aytprime(ind)=-v*EV%denlkt(2,2)*ayt(ind+1)+8._dl/15._dl*k*shear
                      do l=3,EV%lmaxnut-1
                         ind=ind+1
                         aytprime(ind)=v*(EV%denlkt(1,L)*ayt(ind-1)-EV%denlkt(2,L)*ayt(ind+1))
                      end do
                      ind = ind+1
                      !Truncate moment expansion.
                      aytprime(ind)=k*v*EV%lmaxnut/(EV%lmaxnut-2._dl)*ayt(ind-1)-(EV%lmaxnut+3)*cothxor*ayt(ind)
                   else
                      aytprime(ind)=8._dl/15._dl*k*shear
                   end if
                   ind=ind+1
                end do
             end if
          end do
       end if
    end if

    !  Get the propagation equation for the shear

    if (CP%flat) then
       aytprime(3)=-2*adotoa*shear+k*Hchi-rhopi/k
    else
       aytprime(3)=-2*adotoa*shear+k*Hchi*(1+2*CP%curv/k2)-rhopi/k
    endif

    aytprime(2)=-k*shear

  end subroutine derivst



  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

end module GaugeInterface
