!!!!!!!!!!!!!!!! axion_background !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This module computes the evolution of the background in the presence of an axion field
! We use exact equations up until m=3H, and then treat as a w=0 fluid
! For more details, see Hlozek et al, 2014, arXiv:1410.2896
! 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Input cosmological parameters, and receive evolution of equation of state, density, 
!adiabatic sound speed
! Also receive initial field value, oscillation scale factor, equality scale factor
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This code solves the Klein-Gordon equation for a scalar field in an expanding universe (including
!photons, massive and massless neutrinos, cold dark matter, and baryons)
! \ddot{\phi}+2H \phi+m^2 a^2\phi=0 is transformed using the variable definitions
!where dots denote derivatives with respect to conformal time and H
!is the conformal Hubble parameter and a is the cosmological scale factor
!phi=\sqrt(3/(4\pi G)) v_{1}^{twiddle}
!phi_dot=u_2, u_2=H_{0} u_{2}^twiddle, dimensionless conformal time t^twiddle=H_{0} t
!H_{0} is the usual Hubble parameter, and finally u_{2}^{twiddle}=\sqrt{3/(4\pi G)} v_{2}^{twiddle}
!A dimensionless Hubble parameter littleh is defined as H/(100 km /s/ Mpc)
!as is a dimensionless axion mass maxion_twiddle=m/(H_0), where h
!is the dimensionless Hubble parameter
!This transforms the second order ODE into a pair of first order ODES
!In terms of these variables, there is a nice expression for the adiabatic
! sound speed c_ad^{2}=Pdot/\rho_dot=1+(2/3) maxion_twiddle^{2} v_{1}^{twiddle}a^{2}/(v_{2}^{twiddle}*littleh)





module axion_background

contains

  subroutine w_evolve(Params, badflag)
    !Params is the overall parameter structure to be passed back to camb, including the estimated a_osc (RL added tau_osc)
    !value, table for interpolation of the adiabatic sound speed, axion EOS and axion energy density
    !flag for failed histories where H becomes negative (axion drives collapse). Flag is raised in derivs_bg routine or
    !in w_evolve if a history has a turnaround within 10 log(a) bin steps of the last time step
    !Give subroutine access to model parameters from camb, constants from constants.f90, massive neutrino routines

    use ModelParams
    use constants
    use Precision
    use MassiveNu
    implicit none
    type(CAMBparams) :: Params 
    integer badflag ! flag for failed histories where H becomes negative (axion drives collapse). Flag is raised in derivs_bg routine.
    integer i 
    integer k,j !general indices



    !Internal: Neutrino-related constants from CAMB input
    !contribution to H^2/(100 km /s/ Mpc)^2 of massive and massless neutrinos, respectively
    real(dl) lhsqcont_massive(Params%Nu_mass_eigenstates),lhsqcont_massless
    !zeta(3), conversion factor to get neutrino masses in eV for check comparison with rest of CAMB
    real(dl) zeta3, conv
    !Correction factor from massless to massive neutrino energy densities, constant to go from massive neutrino
    !mass fractions to massive neutrino massive values
    real(dl)  rhonu,nu_constant
    real(dl) :: omaxh2_convergencethresh !RL 040926
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! non-trivial aspects to sharing neutrino data structures in new subroutines of CAMB
    ! so we recompute somethings twice, but these are single numbers
    ! should not be a significant slow down



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Internal:: Cosmological Parameters and useful algebraic combinations
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !omega h^2 of various species
    !Matter, Cosmological constant, massive neutrino density
    real(dl) omegah2_m, omegah2_lambda,omnuh2
    !baryons + cold dark matter, Hubble parameter in units of eV, omega_curvature
    real(dl) omegah2_regm,H_ev,rhocrit,omk 
    !baryons, dark matter separately, axions, axion mass in units of Hubble, dimensionless Hubble Parameter
    real(dl) omegah2_b,omegah2_dm,omegah2_ax,maxion_twiddle,hnot
    !Dimensionless Hubble^2
    real (dl) hsq
    !Desired axion fraction omega_a/(omega_a+omega_b+omega_c)
    !real(dl) fax
    ! scale factor at equality estimated under assumption axions redshift purely as CDM
    real(dl) regzeq

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Internal: Evolution control parameters (other than those defined in modules.f90)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    !dfac sets m=dfac*H condition required for code to switch from evolving KG+Einstein equations recast in terms of fluid variables
    !to evolving WKB equations
    !Details in Hlozek et al 2014. arXiv:1410.2896
    real(dl) dfac
    !how much past a_osc to go in the final a integration at the end (a_osc is defined later as the time of scalar field oscillation)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Internal: Evolution variable arrays
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(dl), allocatable :: a_arr(:), v_vec(:,:), littlehfunc(:), littlehfunc_buff(:), diagnostic(:)
    !scalar field and its derivative (in appropriate units), conformal Hubble parameter /(100 km /s/ Mpc) as a function of scale factor
    !m/(dfac*littlehfunc*a) diagnostic parameter that is used to identify approximate aosc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Internal: Fluid quantities of interest, computed after evolution generated
    !Scalar field energy density
    real(dl), allocatable :: rhoaxh2_ov_rhom(:)
    !log_e(m/3H) used for array of this quantity to find a_osc, the time of scalar field oscillation
    !!arrays used to find aosc once final scalar field initial condition is determined
    real(dl), allocatable :: f_arr(:)
    !Array used to find estimate of scale factor at matter-radiation equality (axions counter as matter)
    real(dl), allocatable ::  eq_arr(:)
!!!!!!!!!!!!!!!!!!!!!!!!

    ! NB some of these quantities a_arr, and rhoaxh2_ov_rhom are here so that
    ! we don't have to keep taking logarithms and exponentials of lists
    ! ultimately we just want the logs of these for cubic spline interpolation in the rest of CAMB
    ! but to compute Hubble, field derivatives, we need the quantities themselves and not their logs
    ! By keeping these quantities we can compute them once inside the scalar field routine (per initial condition)
    ! and keep the information to use for when it is needed
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Internal storage arrays for spline interpolation of scalar field solution
    !for fixed scalar initial condition
    real(dl), allocatable :: v_buff(:) ! spline buffer scalar field and its derivative with respect to a
    real (dl), allocatable :: abuff(:) !spline buffer for scale factor interpolation (ie finding 'when
    ! various quantities like rho_rad/rho_m, m/3H, etc, reach values of interest)
    !
    real (dl), allocatable ::  eq_arr_buff(:)
    !nstop is a control index that excludes really bad scalar field evolution histories and
    !prevents them from screwing up the spline
    integer  nstop
!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Quantities of interest that are generated for each scalar field evolution 
    ! and sometimes read into arrays
    ! energy fraction output (omega_ax/(omega_ax+omega_c+omega_b))
    !real(dl) fout
    !Estimated scalar factor when m=dfac*H, scalar field v(1,i) when this value is reached,
    !Estimated scalar field derivative at same point i time
    real(dl) phiosc,phidosc
    !log(aosc), aosc
    real(dl) laosc
    !output value of scalar field value needed to get desired axion density
    !obtained via cubic splines, final fout/fax for the corresponding scalar field history
    !RL real(dl) vtw,final_fout_check
    !critical densities, joining value of rho when axion starts to act like CDM, better log initial value so derivatives can be taken t offset, v2 value at aosc joining point
    !used to set full density evolution of final chosen scalar field history
    real(dl) rhorefp, Prefp !RL added Prefp
    real(dl) v1_ref, v2_ref !RL testing, rhorefp_field
    real(dl) littlehauxi, wcorr_coeff, omaxh2_wcorr, A_coeff, tvarphi_c,tvarphi_cp,tvarphi_s,tvarphi_sp !RL for EFA
    real(dl) lh_skip !RL 020224 for skipping recombination
    real(dl) y_phase !RL 120624 for phase-finding skips
!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Internal storage arrays to sweep through quantities of interest for different 
    !scalar field evolutions corresponding to different initial conditions
    !
    ! omega_ax/(omega_ax+omega_c+omega_b) for different scalar field initial conditions
    !RL real (dl) fout_arr(nphi)
    !array of aosc values for different initial conditions
    ! for the scalar field
    !RL real(dl) aosc_arr(nphi)
    !log of aosc values for different scalar field histories
    !RL real(dl) laosc_arr(nphi)
    !array of fout_arr(nphi)/fax (actual axion mass fraction vs desired axion mass fraction)
    !logarithm later taken for root finding
    !RL real(dl) fout_check_arr(nphi)
    ! stepsize in log of initial scalar field value
    real(dl) dlogvtwiddle_init
    !array of initial v1 values
    !RL real(dl) vtwiddle_init_arr(nphi)
    !initial guess of v1 and their corresponding omaxh2s (RL)
    real(dl) v1_initguess(3), omaxh2_guess(3), aosc_guess(3) !RL 082924
    integer iter_c
    logical bisec_bracketed !RL 033124
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    !WORK NEEDED ON COMMENTS
!!!!!! Simple analytic estimates for required scalar field initial condition
!!! used to set a broad range of scalar initial conditions 
    !Scalar field evolution is generated for each initial condition and
    !cubic spline interpolation is then used to find the one which generates
    !desired axion energy fraction today
    !initial value max and min for phi_init, and final output 
    real(dl) vtwiddle_initmax,vtwiddle_initmin,vtwiddle_init
    !list of initial values for phi_init
    real(dl) vtwiddle_initlist(3)
!!!!!!!!
    !WORK NEEDED ON COMMENTS


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Internal storage arrays for spline interpolation of quantities of interest
    ! for different scalar initial conditions
    !array of axion energy fractions, buffer array for spline fitting in space of initial phi values
    !RL real (dl) vtwiddle_init_arr_buff(nphi)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Spline derivatives at endpoints (used for all splines)
    real(dl) d1,d2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !RL testing the spline for phi
    real(dl) d1v1, d2v1, d1v2, d2v2
!!!!!!!!!!!!!
    !Timing variables
    !real :: clock_start, clock_stop   
!!!!!!!!!!

    !WORK NEEDED ON OMMENTS
    ! numbers for initial and final value of a, and stepsize
    ! initial needed to ensure axion sub-dominance to everything else at a_init
    real(dl) a_init,a_m, a_lambda, a_rel,as_scalar,as_rad,as_matt,a_final,dloga,log_a_final,log_a_init
    !!real(dl) kfinal(1:2),svec(1:16),avec(1:16) 
    real(dl), allocatable :: cmat(:,:), kvec(:,:), kfinal(:), svec(:), avec(:)
    !!    real(dl) nu_massless_degeneracy, fractional_number,actual_massless,neff_i !RL added DG fix 06/27/2023
    !END WORK NEEDED ON COMMENTS


    !RL 041124: deallocate global arrays before the code can be run again with new parameters
    if (allocated(loga_table)) deallocate(loga_table)
    if (allocated(phinorm_table)) deallocate(phinorm_table)
    if (allocated(phidotnorm_table)) deallocate(phidotnorm_table)
    if (allocated(phinorm_table_ddlga)) deallocate(phinorm_table_ddlga)
    if (allocated(phidotnorm_table_ddlga)) deallocate(phidotnorm_table_ddlga)
    if (allocated(rhoaxh2ovrhom_logtable)) deallocate(rhoaxh2ovrhom_logtable)
    if (allocated(rhoaxh2ovrhom_logtable_buff)) deallocate(rhoaxh2ovrhom_logtable_buff)
    !Time Code
    !clock_start = 0.0
    !call cpu_time(clock_start)
    !Set value of neutrino constants as elsewhre in camb
    zeta3=1.2020569031595942853997d0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! DM cosmological parameters etc. in units for integrator
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    hnot = Params%H0/100.d0
    hsq=hnot**2.0d0
    omegah2_dm=Params%omegac*(hnot**2.0d0)
    omegah2_b=Params%omegab*(hnot**2.0d0) 
    omnuh2=Params%omegan*(hnot**2.0d0)   

    omegah2_lambda=Params%omegav*(hnot**2.0d0)  
    omegah2_ax=Params%omegaax*(hnot**2.0d0)   
    maxion_twiddle = Params%m_ovH0    
    !Total densities
    omegah2_regm=omegah2_dm+omegah2_b
    omegah2_m=omegah2_regm+omegah2_ax
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Initialize massive and massless neutrino variables
    
    !!    call Nu_init !RL moved to inidriver_axion
    !06/27/2023, RL incorporating 2nd batch of DG's neutrino fix-------------------
    !!    nu_massless_degeneracy = Params%Num_Nu_massless !N_eff for massless neutrinos
    !!    if (Params%Num_nu_massive > 0) then
    !!       if (Params%Nu_mass_eigenstates==0) stop 'Have Num_nu_massive>0 but no nu_mass_eigenstates'
    !!       if (Params%Nu_mass_eigenstates==1 .and. Params%Nu_mass_numbers(1)==0) Params%Nu_mass_numbers(1) = Params%Num_Nu_Massive
    !!       if (all(Params%Nu_mass_numbers(1:Params%Nu_mass_eigenstates)==0)) Params%Nu_mass_numbers=1 !just assume one for all
    !!       if (Params%share_delta_neff) then
    !!          !default case of equal heating of all neutrinos
    !!          fractional_number = Params%Num_Nu_massless + Params%Num_Nu_massive
    !!          actual_massless = int(Params%Num_Nu_massless + 1e-6_dl)
    !!          neff_i = fractional_number/(actual_massless + Params%Num_Nu_massive)
    !!          nu_massless_degeneracy = neff_i*actual_massless
    !!         Params%Nu_mass_degeneracies(1:Params%Nu_mass_eigenstates) = Params%Nu_mass_numbers(1:Params%Nu_mass_eigenstates)*neff_i
    !!       end if
    !!    else !DG 6/6/2023 Accidentally wasn't counting massless neutrinos unless
    !massive neutrinos were on in first edition of Spring 2023 bugfix now
    !fixed, much simpler, very little degeneracy remultiplication needed
    !Don't need to reshare fractionality of remaining neutrinos so 
    !input parameters transparent
    !!       nu_massless_degeneracy = Params%Num_Nu_massless
    !!    endif
    !   print*,Params%Nu_massless_degeneracy,nu_massless_degeneracy
    !print*,'nu_massless_degneracy',nu_massless_degeneracy
    !RL 06/27/2023--------------------
    

    grhom = 3.0d0*(hsq*1.d10)/(c**2.0d0) 
    grhog = ((kappa/(c**2.0d0)*4.0d0*sigma_boltz)/(c**3.0d0))*(Params%TCMB**4.0d0)*(Mpc**2.0d0) !RL replaced the COBE_CMBTemp
    grhor = (7.0d0/8.0d0)*((4.0d0/11.0d0)**(4.0d0/3.0d0))*grhog 
    !calculate critical density
    rhocrit=(8.0d0*const_pi*G*1.d3/(3.0d0*((1.d7/(MPC_in_sec*c*1.d2))**(2.0d0))))**(-1.0d0)
    Params%omegah2_rad=((Params%TCMB**4.0d0)/(rhocrit))/(c**2.0d0) !RL replaced the COBE_CMBTemp
    Params%omegah2_rad=Params%omegah2_rad*a_rad*1.d1/(1.d4)
    !calculate omega rad using standard formula
    !Contribution of photons and massless neutrinos to H/(100 km /s/Mpc)
    !lhsqcont_massless=(Params%Num_Nu_massless*grhor*(c**2.0d0)/((1.d5**2.0d0)))/3.0d0
    !RL pasted DG's fix or lhsqcont_massless
    lhsqcont_massless=(Params%Nu_massless_degeneracy*grhor*(c**2.0d0)/((1.d5**2.0d0)))/3.0d0

    Params%omegah2_rad=Params%omegah2_rad+lhsqcont_massless
    !const from modules.f90
    nu_constant=(7.0d0/120.0d0)*(const_pi**4.0d0)/(zeta3*1.5d0)
    nu_constant=nu_constant*omnuh2*(grhom/grhor)/hsq

    do k=1,Params%Nu_mass_eigenstates,1
       !Compute neutrino masses corresponding to input mass fractions and degeneracies
       Nu_masses(k)=nu_constant*Params%Nu_mass_fractions(k)/(Params%Nu_mass_degeneracies(k))
       !Compute contribution of massive neutrinos to Hubble parameter/100 km/s/Mpc
       lhsqcont_massive(k)=Params%Nu_mass_degeneracies(k)*(grhor*(c**2.0d0)/((1.d5**2.0d0)))/3.0d0
       !Print some useful checks for neutrinos
       !!       call Nu_rho(Nu_Masses(k),rhonu) !RL FLAG 06/27/2023 DG uncommented this line
       !	print*,Nu_masses(k),omegah2_rad+lhsqcont_massive(1)*rhonu
       !	conv = k_B*(8.0d0*grhor/grhog/7.0d0)**0.25d0*COBE_CMBTemp/elecV * &
       !   &(Params%Nu_mass_degeneracies(k)/dble(Params%Nu_mass_numbers(k)))**0.25d0!
       !	call Nu_rho(1.0d0,rhonu)
       !	print*,Nu_masses(k)*conv
    enddo



    !print*,lhsqcont_massive(1)/Params%Nu_mass_degeneracies(1)

    !! Params%omegah2_rad=Params%omegah2_rad+

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Switch between  evolution regimes
    !definition of m=nH pmatch point, use 3 for convention
    dfac=Params%dfac
    !!if (Params%ma .lt. 1.e-29) then !RL 110823
    Params%wEFA_c = 9._dl/8._dl 
    !!Params%wEFA_c = 1.22_dl !RL testing 040324
    !!else
!!!Params%wEFA_c = 3._dl/2._dl
    !!end if

!!!!!!!!!!!!!
    !!    write(*, *) 'RL 06/27/2023 try printing lhsqcont_massive', lhsqcont_massive

!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Constants for integrator from
    !This is an implementation of the 8th order formulas on page 75 of 
    !Classical Fifth-, Sixth-, Seventh-, and Eighth-Order Runge-Kutta Formulas with Stepsize Control
    !by E Fehlberg
    !http://hdl.handle.net/2060/19680027281
    !(NASA Huntsville, 1968)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Constant matrix
    allocate(cmat(16,16), kvec(2,16)) !RL011025 
    allocate(kfinal(2), svec(16), avec(16))
    avec=0.0d0
    kvec=0.0d0
    avec(1)=0.4436894037649818d0
    avec(2)=0.6655341056474727d0
    avec(3)=0.9983011584712091d0
    avec(4)=0.31550d0
    avec(5)=0.5054410094816906d0
    avec(6)=0.1714285714285714d0
    avec(7)=0.8285714285714285d0
    avec(8)=0.6654396612101156d0
    avec(9)=0.2487831796806265d0
    avec(10)=0.1090d0
    avec(11)=0.8910d0
    avec(12)=0.3995d0
    avec(13)=0.6005d0
    avec(14)=1.0d0
    avec(15)=0.0d0
    avec(16)=1.0d0
    cmat=0.0d0
    cmat(1,1)=avec(1)
    cmat(2,1)=0.1663835264118681d0
    cmat(2,2)=0.49915057d0
    cmat(3,1)=0.24957528d0
    cmat(3,2)=0.0d0
    cmat(3,3)=0.74872586d0
    cmat(4,1)=0.20661891d0
    cmat(4,2)=0.0d0
    cmat(4,3)=0.17707880d0
    cmat(4,4)=-0.68197715d-1
    cmat(5,1)=0.10927823d0
    cmat(5,2)=0.0d0
    cmat(5,3)=0.0d0
    cmat(5,4)=0.40215962d-2
    cmat(5,5)=0.39214118d0
    cmat(6,1)=0.98899281d-1
    cmat(6,2)=0.0d0
    cmat(6,3)=0.0d0
    cmat(6,4)=0.35138370d-1
    cmat(6,5)=0.12476099d0
    cmat(6,6)=-0.55745546d-1
    cmat(7,1)=-0.36806865d0
    cmat(7,2)=0.0d0
    cmat(7,3)=0.0d0
    cmat(7,4)=0.0d0
    cmat(7,5)=-0.22273897d1
    cmat(7,6)=0.13742908d1
    cmat(7,7)=0.20497390d1
    cmat(8,1)=0.45467962d-1
    cmat(8,2)=0.0d0
    cmat(8,3)=0.0d0
    cmat(8,4)=0.0d0
    cmat(8,5)=0.0d0
    cmat(8,6)=0.32542131d0
    cmat(8,7)=0.28476660d0
    cmat(8,8)=0.97837801d-2
    cmat(9,1)=0.60842071d-1
    cmat(9,2)=0.0d0
    cmat(9,3)=0.0d0
    cmat(9,4)=0.0d0
    cmat(9,5)=0.0d0
    cmat(9,6)=-0.21184565d-1
    cmat(9,7)=0.19596557d0
    cmat(9,8)=-0.42742640d-2
    cmat(9,9)=0.17434365d-1
    cmat(10,1)=0.54059783d-1
    cmat(10,2)=0.0d0
    cmat(10,3)=0.0d0
    cmat(10,4)=0.0d0
    cmat(10,5)=0.0d0
    cmat(10,6)=0.0d0
    cmat(10,7)=.11029325d0
    cmat(10,8)=-.12565008d-2
    cmat(10,9)=0.36790043d-2
    cmat(10,10)=-.57780542d-1
    cmat(11,1)=.12732477d0
    cmat(11,2)=0.0d0
    cmat(11,3)=0.0d0
    cmat(11,4)=0.0d0
    cmat(11,5)=0.0d0
    cmat(11,6)=0.0d0
    cmat(11,7)=0.0d0
    cmat(11,8)=0.11448805
    cmat(11,9)=0.28773020
    cmat(11,10)=0.50945379d0
    cmat(11,11)=-0.14799682d0
    cmat(12,1)=-0.36526793d-2
    cmat(12,2)=0.0d0
    cmat(12,3)=0.0d0
    cmat(12,4)=0.0d0
    cmat(12,5)=0.0d0
    cmat(12,6)=0.81629896d-1
    cmat(12,7)=-0.38607735d0
    cmat(12,8)=0.30862242d-1
    cmat(12,9)=-0.58077254d-1
    cmat(12,10)=0.33598659d0
    cmat(12,11)=0.41066880d0
    cmat(12,12)=-0.11840245d-1
    cmat(13,1)=-0.12375357d1
    cmat(13,2)=0.0d0
    cmat(13,3)=0.0d0
    cmat(13,4)=0.0d0
    cmat(13,5)=0.0d0
    cmat(13,6)=-0.24430768d2
    cmat(13,7)=0.54779568d0
    cmat(13,8)=-0.44413863d1
    cmat(13,9)=0.10013104d2
    cmat(13,10)=-0.14995773d2
    cmat(13,11)=0.58946948d1
    cmat(13,12)=0.17380377d1
    cmat(13,13)=0.27512330d2
    cmat(14,1)=-0.35260859d0
    cmat(14,2)=0.0d0
    cmat(14,3)=0.0d0
    cmat(14,4)=0.0d0
    cmat(14,5)=0.0d0
    cmat(14,6)=-0.18396103d0
    cmat(14,7)=-0.65570189d0
    cmat(14,8)=-.39086144d0
    cmat(14,9)=0.26794646d0
    cmat(14,10)=-0.10383022d1
    cmat(14,11)=0.16672327d1
    cmat(14,12)=0.49551925d0
    cmat(14,13)=.11394001d1
    cmat(14,14)=0.51336696d-1
    cmat(15,1)=0.10464847d-2
    cmat(15,2)=0.0d0
    cmat(15,3)=0.0d0
    cmat(15,4)=0.0d0
    cmat(15,5)=0.0d0
    cmat(15,6)=0.0d0
    cmat(15,7)=0.0d0
    cmat(15,8)=0.0d0
    cmat(15,9)=-0.67163886d-2
    cmat(15,10)=0.81828762d-2
    cmat(15,11)=-0.42640342d-2
    cmat(15,12)=0.280090294741d-3
    cmat(15,13)=-0.87835333d-2
    cmat(15,14)=0.10254505d-1
    cmat(15,15)=0.0d0
    cmat(16,1)=-0.13536550d1
    cmat(16,2)=0.0d0
    cmat(16,3)=0.0d0
    cmat(16,4)=0.0d0
    cmat(16,5)=0.0d0
    cmat(16,6)=-0.18396103d0
    cmat(16,7)=-0.65570189d0
    cmat(16,8)=-0.39086144d0
    cmat(16,9)=0.27466285d0
    cmat(16,10)=-0.10464851d1
    cmat(16,11)=0.16714967d1
    cmat(16,12)=0.49523916d0
    cmat(16,13)=0.11481836d1
    cmat(16,14)=0.41082191d-1
    cmat(16,15)=0.0d0
    cmat(16,16)=1.0d0

    svec(1)=0.32256083d-1
    svec(2)=0.0d0
    svec(3)=0.0d0
    svec(4)=0.0d0
    svec(5)=0.0d0
    svec(6)=0.0d0
    svec(7)=0.0d0
    svec(8)=0.0d0
    svec(9)=0.25983725d0
    svec(10)=0.92847805d-1
    svec(11)=.16452330d0
    svec(12)=0.176659510d0
    svec(13)=0.23920102d0
    svec(14)=0.39484274d-2
    svec(15)=0.3072649547580d-1 !RL noting: DG fixed svec(15) third digit

    !END NEEDS WORK
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! begin initialization procedure, shoot for best phi_init
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!NEEDS WORK IN COMMNENTS
    !Initial aosc guess to tell computer to try to find aosc, also a value that flags when axion doesnt start oscillating by today
    aosc_guess = (/15.0d0, 15.0d0, 15.0d0/)
    !compute curvature parameter for this set of omegas
    omk=1.0d0-(omegah2_m+Params%omegah2_rad+omegah2_lambda+omnuh2)/hsq

!!!!!!!!!
    ! Determine a range of epochs to be well before in starting axion field evolution
    !scale factor at time of axions becoming irrelevant in competition with different species
    !and other important transitionary epochs
    as_matt=(omegah2_regm/(maxion_twiddle**2.0d0))**(1.0d0/3.0d0)
    as_rad=(Params%omegah2_rad/(maxion_twiddle**2.0d0))**(1.0d0/4.0d0)
    !scale factor at equality of other matter with radiation
    a_m=(Params%omegah2_rad/(omegah2_regm))
    !This times 10^-7 will be another maximum initial scale factor (see below)
    a_rel=10.0d0
    !Subdominance of dark energy to ordinary radiation
    a_lambda=(Params%omegah2_rad/omegah2_lambda)**(0.25d0)
    !aosc (start well before coherent oscillation of the axion field)
    !RL: is the precise number of a_osc defined yet?
    as_scalar=(omegah2_ax/(maxion_twiddle**2.0d0))**(1.0d0/3.0d0)
    !write(*, *) 'is the precise number of a_osc defined yet?', Params%a_osc, as_scalar
    !find safe initial a (take analytic estimates for when axions are negligible to everything else
    ! and look at scalar factors 7 orders of magnitude smaller as initial scale factors!!!
    a_init=min(a_rel,a_lambda,a_m,as_matt,as_rad,as_scalar)*1.d-8
    !maximum a
    a_final=1.0d0
    !log of various a values, dlog a, step back one
    log_a_init=dlog(a_init)
    log_a_final=dlog(a_final)   
    !set log a bin size
    dloga=(log_a_final-log_a_init)/(dble(ntable-1)) ! RH: DG's original
    allocate(a_arr(ntable)) !RL 111623
    a_arr(1)=dexp(log_a_init)
    !Feed to output parameter table
    !initial loga and a value
    allocate(loga_table(ntable)) !RL 112823
    loga_table(1)=log_a_init
    !Generate arrays of log scale factor for output and internal scale factor
    forall(i=2:ntable)
       loga_table(i)=log_a_init+dloga*dble(i-1)
       a_arr(i)=dexp(loga_table(i))
    end forall
    !RL 033124 - it is too much work to take the entire solve for background thing out for a bracketing part, hence I'll monitor it with a logical
    bisec_bracketed = .false.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

!!!!!!!
    !assuming axion acts like cosmological constant then matter, 
    !analytically obtain
    !initial v_1^twiddle value to get right density today,
    ! under different assumptions
    ! Expressions will agree with arXiv:1410.2896 in appropriate asymptotic limits
    !**************
    !During Axion domination
    !vtwiddle_initlist(1)=fax
    !Oscillation starts during matter domination
    !vtwiddle_initlist(2)=fax*omegah2_m/(dsqrt(maxion_twiddle)*((Params%omegah2_rad**0.750d0)))
    !Oscillation starts during radiation domination
    !vtwiddle_initlist(3)=fax*omegah2_m/(dsqrt(maxion_twiddle)*(Params%omegah2_rad))
    !Concatenate these initial v_1^twiddle
    !vtwiddle_initlist=dsqrt(vtwiddle_initlist)
    !**************

    !-------------------RL replacing this spline procedure to get the initial conditions with bisection
    !First guesses of the initial v_1. This is done with m/H_* = 3 regardless of what the actual m/H_* we're using since this is good enough for the approximation and also we used the assumption that v_1(m/H_* = 3) ~ v_1,ini. The true v1 will be taken care of in the bisection process.
    if (maxion_twiddle .lt. 3.0d0) then !The case where axions never oscillate
       v1_initguess(2) = dsqrt(omegah2_ax)/maxion_twiddle
       !In this case, !maximum a isn't changing and remains 1.0d0
    else if ((maxion_twiddle**2.0d0)/9.0d0 .lt. &
         &(omegah2_m**4.0d0)/(Params%omegah2_rad**3.0d0) + (Params%omegah2_rad**3.0d0)/(omegah2_m**4.0d0)) then !Oscillation starts approximately in matter-domination
       v1_initguess(2) = dsqrt(omegah2_ax)/(3.0d0*dsqrt(omegah2_m)/hnot)
    else !Oscillation starts approximately in radiation-domination
       v1_initguess(2) = dsqrt(omegah2_ax)/(((9.0d0*Params%omegah2_rad/hsq)**0.375d0)*(maxion_twiddle**0.25d0))
    end if
    !RL assigning the initial guesses for bisection
    !Make the range that spans 2x or /2 of the initial guess
    v1_initguess(1) = v1_initguess(2)/2.0d0
    v1_initguess(3) = v1_initguess(2)*2.0d0
    !Reassign the middle entry to be the average of the two ends of the guess range
    v1_initguess(2) = (v1_initguess(1) + v1_initguess(3))/2.0d0

    !Initiate the omaxh2_guesses to be unreasonable values to signal that we haven't computed it
    omaxh2_guess = (/42.0d0, 42.0d0, 42.0d0/)
    !RL 113023: initialize ah_osc as well, to be used as the iterating H^EF. Initially, set to an impossible value, under which case the instantaneous Hubble is computed. Once this is done, littlehauxi will be updated to H^EF and will be used in auxiIC
!!!!littlehauxi = -1._dl

!!!!


    !------------------------------------
    !**************
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Set range of initial scalar field initial conditions to try
    !bracket range of initial vtwiddle 1 values (basically phi 1)
    !Spanning many orders of magnitude beyond the simples guesses above to bracket a wide range of initial conditions
    !and make sure to choose the correct one
    !vtwiddle_initmin=min(vtwiddle_initlist(1),vtwiddle_initlist(2),vtwiddle_initlist(3))/1.0d2
    !vtwiddle_initmax=max(vtwiddle_initlist(1),vtwiddle_initlist(2),vtwiddle_initlist(3))*100.0d1
    !set size of log step in v1 initial
    !dlogvtwiddle_init=(dlog(vtwiddle_initmax)-dlog(vtwiddle_initmin))/(dble(nphi)-1.0d0)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Generate array of initial phi values
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !forall(j=1:nphi)
    !   vtwiddle_init_arr(j)=dexp(dble(j-1)*dlogvtwiddle_init+dlog(vtwiddle_initmin))
    !end forall
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !**************



!!!!!!!!!!!!!!!!!LOOP OVER SCALAR FIELD INITIAL CONDITIONS
    !**************
    !do j=1,nphi,1************RL commenting this line out to replace with a while loop
    !**************
    !Initialize iteration
    iter_c = 0
    allocate(v_vec(2,ntable))
    allocate(littlehfunc(ntable), littlehfunc_buff(ntable), diagnostic(ntable), f_arr(ntable), v_buff(ntable), abuff(ntable))
    do while (iter_c < nphi) !RL reusing nphi to be the max of iterations - iter_c should really be just a handful and not supposed to come even close to 150
       do j = 1,3,1 !Go over the 3 variables, evaluate the corresponding omaxh2
          !Control bad value of aosc -- if a history never gets to m=dfac*H the value of this parameter alerts the code
          ! that this is a history for which coherent oscillation never begins
          if (omaxh2_guess(j) .gt. 1.0d0) then !unreasonable final omaxh2 (including crashed histories or reassigned empty values (42)
             aosc_guess(j)=15.0d0
!!!!initial phi value
             !********************
             !vtwiddle_init=vtwiddle_init_arr(j)
             !********************
             vtwiddle_init=v1_initguess(j) !************ RL

             !start axion at rest at top of hill, with phi guessed at current guess
             v_vec(1,1)=vtwiddle_init
             v_vec(2,1)=0.0d0 !axionCAMB original (1st place, initial)
             !calculate dimensionless hubble including standard matter and axion potential and kinetic energy
             call lh(omegah2_regm,Params%omegah2_rad,omegah2_lambda,omk,hsq,&
                  &maxion_twiddle,a_arr(1),v_vec(1:2,1),littlehfunc(1),badflag,&
                  lhsqcont_massless,lhsqcont_massive,Params%Nu_mass_eigenstates,Nu_masses)
             !RL modifying v_vec(2,1) after obtaining the initial Hubble
             v_vec(2,1)= - vtwiddle_init * (a_arr(1)**2.0d0) * (maxion_twiddle**2.0d0) * hnot/(5.0d0 * littlehfunc(1)) 

             !calculate initial value of dphi/da and dphidot/da
             !Compute first step parameters for the pair of first-order ODEs being solved
             kvec=0.0d0
             kfinal=0.0d0
             call next_step(a_arr(1),v_vec(1:2,1),kvec(1:2,1:16),kfinal(1:2),avec(1:16),&
                  &omegah2_regm,Params%omegah2_rad,&
                  &omegah2_lambda,omk,hsq,&
                  &maxion_twiddle,badflag,dloga,16,cmat(1:16,1:16),lhsqcont_massless,lhsqcont_massive,&
                  &Params%Nu_mass_eigenstates,Nu_masses)

             diagnostic(1)=dfac*littlehfunc(1)/(a_arr(1)*hnot) !RL fixed by adding hnot

!!!!

!!!!!!!!!!!!!!!!!!!!!!
             !Integration performed using eight order Runge-Kutta method
             !Derivatives at 8 points in advance of the one in question (with trial function values)
             !Added to current value using optimal quadrature coefficients under certain assumptions
             !This is an implementation of the 8th order formulas on page 75 of 
             !Classical Fifth-, Sixth-, Seventh-, and Eighth-Order Runge-Kutta Formulas with Stepsize Control
             !by E Fehlberg
             !http://hdl.handle.net/2060/19680027281
             !(NASA Huntsville, 1968)
             !It turns out a reasonably accurate integrator is required to accurately obtain the adiabatic sound speed at earlier times
             !Using super-horizon analytic solutions, we found that this integrator + 5000 grid points
             !was necessary to avoid exciting a non-physically large low-l ISW effect
             ! DM please add comments on which condition led us to this realiztion???



             do i=2,ntable,1
!!!!integrate ODE using 16 pt (8th order Runge-Kutta) rule
                !increment fluid (homogeneous values) using precomputed steps
                v_vec(:,i)=v_vec(:,i-1)+(svec(1)*kvec(:,1)+svec(9)*kvec(:,9)+svec(10)*kvec(:,10)&
                     &+svec(11)*kvec(:,11)+svec(12)*kvec(:,12)+svec(13)*kvec(:,13)&
                     &+svec(14)*kvec(:,14)+svec(15)*kvec(:,15))
!!!!
                !calculate hubble for next step
                call lh(omegah2_regm,Params%omegah2_rad,omegah2_lambda,omk,hsq,&
                     &maxion_twiddle,a_arr(i),v_vec(1:2,i),littlehfunc(i),badflag,&
                     &lhsqcont_massless,lhsqcont_massive,Params%Nu_mass_eigenstates,Nu_masses)

                kvec=0.0d0
                kfinal=0.0d0

                !Compute next steps in scalar field and its derivative
                call next_step(a_arr(i),v_vec(1:2,i),kvec(1:2,1:16),kfinal(1:2),&
                     &avec(1:16),omegah2_regm,Params%omegah2_rad,omegah2_lambda,omk,hsq,&
                     &maxion_twiddle,badflag,dloga,16,cmat(1:16,1:16),&
                     &lhsqcont_massless,lhsqcont_massive,Params%Nu_mass_eigenstates,Nu_masses)

                !compare m to nH and identify a reasonable guess for when a given history crosses the m=3H 
                ! condition (and thus when the code will switch from pure scalar evolution to coherent oscillation)
                diagnostic(i)=dfac*littlehfunc(i)/(a_arr(i)*hnot) !RL fixed by adding hnot

                ! find first guess for aosc 
                if (i .gt. 1) then
                   if (aosc_guess(j) .eq. 15.0d0) then
                      if (maxion_twiddle .lt.diagnostic(i-1)) then
                         if (maxion_twiddle .ge.diagnostic(i)) then
                            aosc_guess(j)=a_arr(i)
                         endif
                      endif
                   endif
                   !RL 020124 - also find rough guess for dfac at z = 800. There's no need to be very precise since it's just a threshold and the conservativeness of the skip will be taken care of in the inidriver_axion.F90 iteration
                   !Note that I can't use the condition aosc .eq. 15.0d0 but have to use the condition iter_c = 0, since a_arr(ntable) will be changed according to the aosc found in this iteration (but aosc might still be 15.0 due to the reassignment of omaxh2 in the end, since we had wanted to recalculate aosc), and you'll run into trouble when a_arr(ntable) < z=800
                endif
                !RL printing to see how the ODE solver is called
                !write(*, *) 'from the Runge-Kutta block'
             enddo
!!!!now there's an array of dlog(diagnostic), use a simple spline to find aosc

             f_arr(1:ntable)=maxion_twiddle/(diagnostic(1:ntable))
             f_arr=dlog(f_arr)
                


             if (maxion_twiddle .lt. 3._dl .or. Params%ma .lt. 1.44e-32_dl) then
                !(case of no oscillation and m/H_* < 10, just use result from full scalar field evolution) !RL 061124
                !write(*, *) 'maxion_twiddle', maxion_twiddle
                !RL f_arr=dexp(f_arr)
                !RL fout=(v_vec(2,ntable)/a_arr(ntable))**2.0d0+(maxion_twiddle*v_vec(1,ntable))**2.0d0
                omaxh2_guess(j) = (v_vec(2,ntable)/a_arr(ntable))**2.0d0+(maxion_twiddle*v_vec(1,ntable))**2.0d0
                
                !aosc_arr(j)=aosc
             else !RL 061124

                if (aosc_guess(j) .ne. 15._dl) then

                   !d1=(loga_table(2)-loga_table(1))/(f_arr(2)-f_arr(1))
                   !d2=(loga_table(ntable)-loga_table(ntable-1))/(f_arr(ntable)-f_arr(ntable-1))
                   !RL switched to natural spline which will be triggered when d1 and d2 are both super large
                   d1 = 1.0d50
                   d2 = 1.0d50
                   call spline(f_arr(1:ntable),(loga_table(1:ntable)),ntable,d1,d2,abuff(1:ntable))
                   !print*, 'Rayne is natural spline correctly implemented?', abuff(1), abuff(ntable)
                   call spline_out(f_arr(1:ntable),loga_table(1:ntable),&
                        &abuff(1:ntable)&
                        &,ntable,0.0d0,laosc)
                   !write(*, *) 'Splining iter_c, j, laosc', iter_c, j, laosc

                   !use a simple spline to calculate phi at time oscillation begins (important to get relic density later, which will be off by a factor of order unity depending on how aosc is defined)
                   !Finds root of equation log(m/3H)=0 using a cubic spline
                   !d1=(v_vec(1,2)-v_vec(1,1))/(loga_table(2)-loga_table(1))
                   !d2=(v_vec(1,ntable)-v_vec(1,ntable-1))/(loga_table(ntable)-loga_table(ntable-1))

                   !RL 010925 comment out the NaN check here. The isnan function here is not part of the standard library, and this block isn't particularly useful
!!!!!!                   do i=1,ntable-10,1
                      !!       print*,i,isnan(f_arr(i))
                      !
!!!!!!                      do k=0,10,1
!!!!!!                         call lh(omegah2_regm,Params%omegah2_rad,omegah2_lambda,omk,hsq,&
!!!!!!                              &maxion_twiddle,a_arr(i+k),v_vec(1:2,i+k),littlehfunc(i+k),badflag,&
!!!!!!                              &lhsqcont_massless,lhsqcont_massive,Params%Nu_mass_eigenstates,Nu_masses)
!!!!!!                         if (((isnan(f_arr(i+k))).or.(isnan(v_vec(1,i+k))) ).or.(isnan(v_vec(2,i+k))))then
!!!!!!                            print*,i,k,littlehfunc(i+k),v_vec(1,i+k),v_vec(2,i+k),a_arr(i+k),dexp(laosc)
!!!!!!                         endif
!!!!!!                      enddo
!!!!!!                   enddo

                else 
                   !!write(*, *) 'Switch is after the present day, switching at a = 1.0 - 1e-3 instead'!RL 061124
                   aosc_guess(j) = 1.0_dl - 1.e-3_dl
                   laosc = dlog(aosc_guess(j)) !RL 061124
                end if
                
                !RL replacing it with the analytical derivatives - note here it's natural log, not log 10!!----------------
                d1 = v_vec(2,1)*hnot/littlehfunc(1) 
                d2 = v_vec(2,ntable)*hnot/littlehfunc(ntable)
                !---------------------

                !Find scalar field value at moment of m=nH (fiducial onset of coherent oscillation)
                call spline(loga_table(1:ntable),(v_vec(1,1:ntable)),ntable,d1,d2,&
                     v_buff(1:ntable))
                call spline_out(loga_table(1:ntable),(v_vec(1,1:ntable)),&
                     &v_buff(1:ntable)&
                     &,ntable,laosc,phiosc)  

                !Find scalar field derivative with scale factor at moment of m=nH (fiducal onset of coherent
                !oscillation)  
                !use a cubic spline to calculate phidot at time rolling begins, we can see perfect equipartition between kinetic and potential scalar field energy is not yet achieved at m=3H, but this is good enough for government work
                !d1=(v_vec(2,2)-v_vec(2,1))/(loga_table(2)-loga_table(1))
                !d2=(v_vec(2,ntable)-v_vec(2,ntable-1))/(loga_table(ntable)-loga_table(ntable-1))
                !RL replacing it with the analytical derivatives - note here it's natural log, not log 10!!----------------
                d1 = -(2.0d0*v_vec(2,1) + &
                     & ((maxion_twiddle*dexp(loga_table(1)))**2.0d0)*v_vec(1,1)*hnot/littlehfunc(1))
                d2 = -(2.0d0*v_vec(2,ntable) + &
                     &((maxion_twiddle*dexp(loga_table(ntable)))**2.0d0)*v_vec(1,ntable)*hnot/littlehfunc(ntable))
                !---------------------
                call spline(loga_table(1:ntable),(v_vec(2,1:ntable)),ntable,d1,d2,v_buff(1:ntable))   
                call spline_out(loga_table(1:ntable),(v_vec(2,1:ntable)),&
                     &v_buff(1:ntable)&
                     &,ntable,laosc,phidosc)
                aosc_guess(j)=dexp(laosc)
                !RL: flag lest my upper bound for aosc is too small - unlikely to happen, but still
                ![Commented out now, might not be necessary due to the aosc = 15 initialization
                !if (aosc .gt. a_final) then
                !   print* 'warning: end value of a too small that aosc is not covered'
                !   stop
                !end if

                !RL: Compute the auxiliary initial conditions
                call auxiIC(Params, omegah2_regm,Params%omegah2_rad,omegah2_lambda,omk,hnot,&
                     &maxion_twiddle,aosc_guess(j), (/phiosc, phidosc/), badflag,lhsqcont_massless,&
                     &lhsqcont_massive,Params%Nu_mass_eigenstates,Nu_masses, littlehauxi, &
                     &Params%ahosc_ETA, A_coeff, tvarphi_c,tvarphi_cp,tvarphi_s,tvarphi_sp, rhorefp, Prefp)
                !RL: rhorefp now is the effective fluid version - the omax(a)h2 at aosc
                !rhorefp = (maxion_twiddle**2.0d0)*(tvarphi_c**2.0d0 + tvarphi_s**2.0d0 + (tvarphi_cp**2.0d0 + tvarphi_sp**2.0d0)/2.0d0 - tvarphi_c*tvarphi_sp + tvarphi_s*tvarphi_cp)

                !also the effective fluid version of the pressure
                !Prefp = (maxion_twiddle**2.0d0)*(tvarphi_cp**2.0d0/2.0d0 + tvarphi_sp**2.0d0/2.0d0 - tvarphi_c*tvarphi_sp + tvarphi_s*tvarphi_cp)
                !calculate axion energy density fraction (compared to total matter+axion)
                !parcel these numbers of interest into an array which sweeps through different v1 values  
                !Beyond a>aosc, (RL:modifying the a^-3 scaling with H/m corrections)
                !wcorr_coeff = littlehauxi*aosc/(maxion_twiddle*hnot)
                wcorr_coeff = Params%ahosc_ETA*aosc_guess(j)/(maxion_twiddle*hnot) !RL082924
                !Here a is taken to be 1.0 to compute the present-day omaxh2
                !RL 110923
                !!write(*, *) 'In bisection, Params%wEFA_c', Params%wEFA_c
                omaxh2_wcorr = rhorefp*(aosc_guess(j)**3.0d0)*dexp((wcorr_coeff**2.0d0)*3.0d0*&
                     &Params%wEFA_c*(1.0d0 - 1.0d0/(aosc_guess(j)**4.0d0))/4.0d0)                  

                omaxh2_guess(j) = omaxh2_wcorr
                !RL fout=omaxh2_wcorr/(omegah2_regm+omaxh2_wcorr)
                
                !RL 020124 obtain Hubble at z = 800 if that's later than aosc
                !RL 022624 modified to obtain Hubble at z = 800 at all times since there's a window to skip before the recombination window which depends on aeq, and aeq can only be evaluated later
                !RL 010425 modified again to restrain the evaluation to later than aosc since now we use aeq_LCDM
!!!!!!!!if (aosc .lt. Params%a_skip .and. aosc .gt. Params%a_skipst) then
                !Obtain the ah from EFA at the recombination jump redshift
                !RL 010425: notice here I used v_vec(1:2,1), which is dummy and not used in lh() since rho_EFA is present
                if (aosc_guess(j) .lt. Params%a_skip .and. aosc_guess(j) .ge. Params%a_skipst) then
                   call lh(omegah2_regm,Params%omegah2_rad,omegah2_lambda,omk,hsq,&
                        &maxion_twiddle,Params%a_skip,v_vec(1:2,1),lh_skip,badflag,&
                        &lhsqcont_massless,lhsqcont_massive,Params%Nu_mass_eigenstates,Nu_masses,&
                        &rhorefp*((aosc_guess(j)/Params%a_skip)**3.0d0)*dexp((wcorr_coeff**2.0d0)*3.0d0*&
                        &Params%wEFA_c*(1.0d0/(Params%a_skip**4.0d0) - 1.0d0/(aosc_guess(j)**4.0d0))/4.0d0))
                   !!write(*, *) 'Rayne, dfac, littlehauxi/aosc, Params%ahosc_ETA/aosc', Params%dfac, littlehauxi/aosc, Params%ahosc_ETA/aosc
                   !write(*, *) 'Rayne, dfac, Params%ahosc_ETA/aosc, m/<H>', Params%dfac, Params%ahosc_ETA/aosc, Params%dfac*(littlehauxi/Params%ahosc_ETA)
                   Params%dfac_skip = min(littlehauxi/Params%ahosc_ETA, 1._dl)*lh_skip
                   !Immediately reuse to get the skip dfac
                   Params%dfac_skip = (maxion_twiddle*Params%a_skip*hnot/Params%dfac_skip)

                   !1.01 is a conservative scaling to prevent infinite loops in inidriver_axion.F90
                   !!write(*, *) 'Rayne, in background, aosc, Params%a_skip, Params%dfac, Params%dfac_skip, rhof', aosc, Params%a_skip, Params%dfac, Params%dfac_skip, rhorefp*((aosc/Params%a_skip)**3.0d0)*dexp((wcorr_coeff**2.0d0)*3.0d0*Params%wEFA_c*(1.0d0/(Params%a_skip**4.0d0) - 1.0d0/(aosc**4.0d0))/4.0d0)
!!!!!!!else
!!!!!!!   Params%dfac_skip = 0._dl
!!!!!!!end if
                end if                

             end if
          end if
          

          !Convergence criterion for the phi_ini shooting loop: the trial abundance must
          !match the requested one to this fractional precision. Tightening it further has
          !negligible cost but no measurable effect on observables.
          omaxh2_convergencethresh = 1.0d-6

          if (abs(omaxh2_guess(j)/omegah2_ax - 1.0d0) .lt. omaxh2_convergencethresh) then

             !Diagnostic for background abundance-convergence problems. Left in place,
             !commented, so it can be switched back on quickly when the shooting loop
             !misbehaves for a new mass/abundance combination.
             !print*, 'AxiECAMB abundance converged: j =', j, ' omaxh2(trial), omaxh2(target), frac. diff. =', &
             !     omaxh2_guess(j), omegah2_ax, abs(omaxh2_guess(j)/omegah2_ax - 1.0d0)
             vtwiddle_init = v1_initguess(j)
             !Don't forget to assign aosc too
             Params%a_osc= aosc_guess(j)
             !if (Params%a_osc .ge. 1.0d0) then
             !   Params%a_osc=1.0d0
             !endif
             iter_c = -1
             exit
          end if
       enddo
          !write(*, *) 'Rayne, iter_c, aosc_guess(1), aosc_guess(2), aosc_guess(3), omaxh2_guess(1), omaxh2_guess(2), omaxh2_guess(3), v1_initguess(1), v1_initguess(2), v1_initguess(3)', iter_c, aosc_guess(1), aosc_guess(2), aosc_guess(3), omaxh2_guess(1), omaxh2_guess(2), omaxh2_guess(3), v1_initguess(1), v1_initguess(2), v1_initguess(3)
          !!write(040324, '(36e52.42)') iter_c, aosc, Params%wEFA_c, omaxh2_guess(1), omaxh2_guess(2), omaxh2_guess(3), v1_initguess(1), v1_initguess(2), v1_initguess(3)
          
          if (iter_c .eq. -1) exit !The previous exit is to exit that if statement if I found the solution, but I still need to exit the while loop
          if (iter_c .eq. nphi - 1) then !Warning sign if iter_c exceeds maximum iteration
             print*, 'Warning: exceeding the maximum number of iteration for bisection: &
                  &number of iterations:', iter_c, 'omaxh2 result:', omaxh2_guess(2), 'omaxh2 input:', &
                  &omegah2_ax, 'fractional error:', omaxh2_guess(2)/omegah2_ax - 1.0d0,&
                  &'aosc:', aosc_guess(2),'. The code will proceed - possibly wEF was not set to converge &
&accurately enough which is ok, but please do check for unreasonable inputs.'
             vtwiddle_init = v1_initguess(2)
             !Don't forget to assign aosc too
             Params%a_osc= aosc_guess(2)
             if (Params%a_osc .ge. 1.0d0) then
                Params%a_osc=1.0d0
             endif
             exit
             !RL: note we didn't go ahead and assign the vtwiddle_init and corresponding aosc, hence the code will crash if this error appears (tested with a very small nphi). But under normal circumstances this should not happen
          end if
          !Applies to the initial cases, if the two initially guessed v_ini don't straddle the target
          if ((omaxh2_guess(1) - omegah2_ax)*(omaxh2_guess(3) - omegah2_ax) > 0 .and. .not. bisec_bracketed) then
             !!write(*, *) 'Initial phase, still bracketing'
             if (omaxh2_guess(3) < omegah2_ax) then
                v1_initguess(3) = v1_initguess(3)*2.0d0
                omaxh2_guess(3) = 42.0d0
             else if (omaxh2_guess(1) > omegah2_ax) then
                v1_initguess(1) = v1_initguess(1)/2.0d0
                omaxh2_guess(1) = 42.0d0
             end if
             !If this doesn't apply, we can do bisection
          else if ((omaxh2_guess(1) - omegah2_ax)*(omaxh2_guess(3) - omegah2_ax) < 0 .and. .not. bisec_bracketed) then
             !!write(*, *) 'Bracket found, enter bisection'
             if  ((omaxh2_guess(1) - omegah2_ax)*(omaxh2_guess(2) - omegah2_ax) < 0) then
                v1_initguess(3) = v1_initguess(2)
                omaxh2_guess(3) = omaxh2_guess(2)
             else if ((omaxh2_guess(3) - omegah2_ax)*(omaxh2_guess(2) - omegah2_ax) < 0) then
                v1_initguess(1) = v1_initguess(2)
                omaxh2_guess(1) = omaxh2_guess(2)
             end if
             bisec_bracketed = .true.
          else if ((omaxh2_guess(1) - omegah2_ax)*(omaxh2_guess(3) - omegah2_ax) < 0 .and. bisec_bracketed) then
             !!write(*, *) 'Bisection mid-process, iter_c, omaxh2_guess(1), omaxh2_guess(2), omaxh2_guess(3)'
             !!write(*, *) iter_c
             !!write(*, '(36e52.42)') omaxh2_guess(1), omaxh2_guess(2), omaxh2_guess(3)
             !!write(*, *) 'Bisection mid-process, v1_initguess(1), v1_initguess(2), v1_initguess(3)'
             !!write(*, '(36e52.42)') v1_initguess(1), v1_initguess(2), v1_initguess(3)
             !!write(*, *) 'Bisection mid-process, aosc_guess(1), aosc_guess(2), aosc_guess(3)'
             !!write(*, '(36e52.42)') aosc_guess(1), aosc_guess(2), aosc_guess(3)
             if  ((omaxh2_guess(1) - omegah2_ax)*(omaxh2_guess(2) - omegah2_ax) < 0) then
                v1_initguess(3) = v1_initguess(2)
                omaxh2_guess(3) = omaxh2_guess(2)
                aosc_guess(3) = aosc_guess(2)
             else if ((omaxh2_guess(3) - omegah2_ax)*(omaxh2_guess(2) - omegah2_ax) < 0) then
                v1_initguess(1) = v1_initguess(2)
                omaxh2_guess(1) = omaxh2_guess(2)
                aosc_guess(1) = aosc_guess(2)
             end if
             !The most unlikely case (and forbidden) that the bracket is lost after bisection begins. If that's the case, Huston we have a problem and the code should be stopped
          else if ((omaxh2_guess(1) - omegah2_ax)*(omaxh2_guess(3) - omegah2_ax) > 0 .and. bisec_bracketed) then
             write(*, *) 'Bisection already started, bracket lost, stop.'
             stop
          end if

          !Regardless of finding the right initial guesses or doing bisection, the reassignment of the middle value is the same
          !!write(*, *) 'Before v1guess(2) reassignment, v1_initguess(1), v1_initguess(2), v1_initguess(3)'
          !!write(*, '(36e52.42)') v1_initguess(1), v1_initguess(2), v1_initguess(3)
          v1_initguess(2) = (v1_initguess(1) + v1_initguess(3))/2.0_dl
          !!write(*, *) 'After v1guess(2) reassignment, v1_initguess(1), v1_initguess(2), v1_initguess(3)'
          !!write(*, '(36e52.42)') v1_initguess(1), v1_initguess(2), v1_initguess(3)

          omaxh2_guess(2) = 42.0d0
          aosc_guess(2) = 15.0d0
          !Adjust the end value of a if the higher omaxh2 guess is close to the true value by approx. 10%
          if (omaxh2_guess(3)/omegah2_ax .lt. 1.1d0 .and. a_final .eq. 1.0d0 .and. aosc_guess(3)*1.1d0 .lt. 1.0d0) then

             !!write(*, *) 'Rayne threshold aosc*1.1d0', aosc*1.1d0
             a_final=aosc_guess(3)*1.1d0 !multiply that aosc by 1.1 just to be safe
             !reassign the log_a_final
             log_a_final=dlog(a_final)   
             !reset log a bin size
             dloga=(log_a_final-log_a_init)/(dble(ntable-1)) ! RH: DG's original
             !Feed to output parameter table
             loga_table(1)=log_a_init
             !Generate arrays of log scale factor for output and internal scale factor
             forall(i=2:ntable)
                loga_table(i)=log_a_init+dloga*dble(i-1)
                a_arr(i)=dexp(loga_table(i))
             end forall
          end if
          !if (aosc .lt. 1.0_dl) then
          !Update H used in auxiIC with H^EF of this iteration
!!!!write(*, *) 'Rayne, littlehauxi in bisection before reassignment', littlehauxi
!!!!call lh(omegah2_regm,Params%omegah2_rad,omegah2_lambda,omk,hsq,maxion_twiddle,aosc,(/phiosc, phidosc/),littlehauxi,badflag,lhsqcont_massless,lhsqcont_massive,Params%Nu_mass_eigenstates,Nu_masses, rhorefp)
!!!!write(*, *) 'Rayne, littlehauxi in bisection after reassignment', littlehauxi
          !Update wEFA with wEF of this iteration, AFTER the reassignment of H
          !!write(*, *) 'Rayne, wEFA in bisection before reassignment', Params%wEFA_c
          !   Params%wEFA_c = (Prefp/rhorefp)/((littlehauxi/(maxion_twiddle*hnot*aosc))**2._dl)
          !!write(*, *) 'Rayne, wEFA in bisection after reassignment', Params%wEFA_c
          !end if
          !Proceed one step in the iteration
          iter_c = iter_c + 1
       end do
       deallocate(f_arr, v_buff, abuff) !RL 111623
       deallocate(cmat, kvec, kfinal, svec, avec) !RL 112223



       !Compute m/3H as a function of scale factor, as well as the scalar field equation of state
       !and appropriately normalized energy density
       allocate(rhoaxh2_ov_rhom(ntable))
       allocate(phinorm_table(ntable), phidotnorm_table(ntable))
       allocate(phinorm_table_ddlga(ntable), phidotnorm_table_ddlga(ntable))
       forall(i=1:ntable)
          diagnostic(i)=dfac*littlehfunc(i)/(a_arr(i)*hnot) !RL fixed by adding hnot
          !RL: trying to save phi and phidot too - but with a factor of sqrt(3/4 pi G) and sqrt(3/4 pi G) H0, as specified for v_vec
          phinorm_table(i) = v_vec(1,i)
          phidotnorm_table(i) = v_vec(2,i)

          !!tabulate axion energy density
          rhoaxh2_ov_rhom(i)=(v_vec(2,i)/a_arr(i))**2.0d0+(maxion_twiddle*v_vec(1,i))**2.0d0
       end forall
       deallocate(diagnostic)
       !write(*, *) 'Rayne, the initial v1 and rhorefp/(maxion_twiddle**2.0d0) in the output table'
       !write(*, '(36e52.42,\)') v_vec(1,1), rhorefp/(maxion_twiddle**2.0d0)
       !write(*, *) 'Rayne, the final grhoaxtableinternal in the output table', rhoaxh2_ov_rhom(ntable)
       !write(*, *) 'Rayne, the last entry of rhoaxh2_ov_rhom', rhoaxh2_ov_rhom(ntable)
       !write(*, *) 'Rayne, constructing grhoax from the last entry of a_arr, v_vec(1), v_vec(2)', (v_vec(2,ntable)/a_arr(ntable))**2.0_dl+ (maxion_twiddle * v_vec(1,ntable))**2.0_dl
       !This is the axion energy density *h^2/(3H_0^2/8\pi G), where h is the dimensionless Hubble Parameter
       !so in camb's definitions grhoa2_axion=grhom*grhox_table_internal(i)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! Spline final arrays

       !RL rewriting the whole thing----------------------------
       !First change the basis of loga_table, since we're doing it anyway
       !write(*, *) 'Rayne, dexp(loga_table(ntable)) before changing basis', dexp(loga_table(ntable))
       loga_table=dlog10(dexp(loga_table))
       
       !The above change of basis test passed, not this problem
       !Then spline the background field variables
       !RL using analytical derivative at the boundary - v1; switched to log10 base

       d1 = dlog(10._dl)*phidotnorm_table(1)*hnot/littlehfunc(1) 
       d2 = dlog(10._dl)*phidotnorm_table(ntable)*hnot/littlehfunc(ntable)
       !d1v1 = d1
       !d2v1 = d2
       call spline(loga_table(1:ntable),phinorm_table,ntable,d1,d2,phinorm_table_ddlga)
       !Then v2
       d1 = -dlog(10._dl) * (2._dl*phidotnorm_table(1) + &
            &((maxion_twiddle*(10._dl**(loga_table(1))))**2._dl)*phinorm_table(1)*hnot/littlehfunc(1))
       d2 = -dlog(10._dl) * (2._dl*phidotnorm_table(ntable) + &
            &((maxion_twiddle*(10._dl**(loga_table(ntable))))**2._dl)*phinorm_table(ntable)*hnot/littlehfunc(ntable))
       !d1v2 = d1
       !d2v2 = d2

       call spline(loga_table(1:ntable),phidotnorm_table,ntable,d1,d2,phidotnorm_table_ddlga)

       !Then, with complete history in hand, spline to find axion energy density at a=aosc
       call spline_out(loga_table,phinorm_table,phinorm_table_ddlga,ntable,dlog10(Params%a_osc),v1_ref)
       call spline_out(loga_table,phidotnorm_table,phidotnorm_table_ddlga,ntable,dlog10(Params%a_osc),v2_ref)
       !Compute the coefficient the auxiliary initial conditions
       !Since this touches all cosmological components a subroutine is written
       !write(*, *) 'Rayne, before the final auxiIC is called, Params%a_osc, v1_ref, v2_ref', Params%a_osc, v1_ref, v2_ref
       call auxiIC(Params, omegah2_regm,Params%omegah2_rad,omegah2_lambda,omk,hnot,&
            &maxion_twiddle,Params%a_osc, (/v1_ref, v2_ref/), badflag,lhsqcont_massless,&
            &lhsqcont_massive,Params%Nu_mass_eigenstates,Nu_masses, littlehauxi, &
            &Params%ahosc_ETA, A_coeff, tvarphi_c,tvarphi_cp,tvarphi_s,tvarphi_sp, rhorefp, Prefp)
       

       !Store the background EFA field variables at the switch
       Params%ah_osc = littlehauxi
       !write(*, *) 'Rayne, check the last littlehauxi: Params%ah_osc/Params%a_osc, Params%ahosc_ETA/Params%a_osc', Params%ah_osc/Params%a_osc, Params%ahosc_ETA/Params%a_osc
       Params%A_coeff = A_coeff
       !RL testing another A_coeff - computed here but only used for the perturbation 
       Params%A_coeff_alt =  A_coeff + 2.0d0*littlehauxi/(Params%a_osc*hnot*maxion_twiddle)
       !write(*, *) 'Rayne, Params%A_coeff, Params%A_coeff_alt', Params%A_coeff, Params%A_coeff_alt
       Params%tvarphi_c = tvarphi_c
       Params%tvarphi_s = tvarphi_s
       Params%tvarphi_cp = tvarphi_cp
       Params%tvarphi_sp = tvarphi_sp
       !rhorefp now is the effective fluid version
       !!rhorefp = (maxion_twiddle**2.0d0)*(tvarphi_c**2.0d0 + tvarphi_s**2.0d0 + (tvarphi_cp**2.0d0 + tvarphi_sp**2.0d0)/2.0d0 - tvarphi_c*tvarphi_sp + tvarphi_s*tvarphi_cp)
       !also the effective fluid version of the pressure
       !!Prefp = (maxion_twiddle**2.0d0)*(tvarphi_cp**2.0d0/2.0d0 + tvarphi_sp**2.0d0/2.0d0 - tvarphi_c*tvarphi_sp + tvarphi_s*tvarphi_cp)
       !rhorefp = (v2_ref)**2.0_dl/(Params%a_osc**2.0_dl)+(maxion_twiddle*v1_ref)**2.0_dl
       !Rescale for CAMB using hsq and store to derived type (CP)
       Params%rhorefp_ovh2=rhorefp/hsq
       Params%Prefp = Prefp !RL: I don't want to rescale with hsq anymore...
       !Temporarily using the background energy density here with finite difference spline boundary conditions
       !Switch to w=0 at aosc
       !Final approach guideline: analytical boundary derivatives for the splined variable (e.g. dlog(rhoax)/dlna), table extends to sufficient "a" that no problems occur in splining at a=aosc, never use the lookup table to define density after aosc.

       !RL: since we don't use the table beyond a_osc, and we need spline tables to be an extension of KG, no rescaling used here
       !do i=1,ntable,1
       !   if (a_arr(i) .gt. Params%a_osc) then
       !      !RL changing the scaling relation from w=0 to an H-varying w with H\propto a^-2 (RD condition)
       !      wcorr_coeff = Params%ah_osc*Params%a_osc/(maxion_twiddle*hnot)
       !      rhoaxh2_ov_rhom(i)=rhorefp*((Params%a_osc/a_arr(i))**3.0d0)*dexp((wcorr_coeff**2.0d0)*9.0d0*(1.0d0/(a_arr(i)**4.0d0) - 1.0d0/(Params%a_osc**4.0d0))/8.0d0)
       !   endif
       !enddo
       allocate(rhoaxh2ovrhom_logtable(ntable), rhoaxh2ovrhom_logtable_buff(ntable))
       rhoaxh2ovrhom_logtable=dlog10(rhoaxh2_ov_rhom)
       !RL replaced with analytical boundary condition (note that here it's dlog10(rho)/dlog10(a))

       d1 = -6.0_dl*(phidotnorm_table(1)/(10._dl**(loga_table(1))))**2._dl/rhoaxh2_ov_rhom(1)
       d2 = -6.0_dl*(phidotnorm_table(ntable)/(10._dl**(loga_table(ntable))))**2._dl/rhoaxh2_ov_rhom(ntable)
    
       call spline(loga_table,rhoaxh2ovrhom_logtable,ntable,d1,d2,rhoaxh2ovrhom_logtable_buff)
       

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!calculate axion adiabatic sound speed Pdot/rhodot, asymptoting to constant (physical) value at low a, when small machine numbers start to behave badly
       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!
       !use same spline method used in rest of subroutine to find new aeq value in this cosmological history (will be needed in initial condition code to normalize, time, scale factor, etc, in particular for isocurvature mode
       allocate(eq_arr(ntable), eq_arr_buff(ntable))
       forall(i=1:ntable)
          eq_arr(i)=((omegah2_regm/(a_arr(i)**3.0d0)+rhoaxh2_ov_rhom(i))/(Params%omegah2_rad/(a_arr(i)**4.0d0)))
       end forall
       eq_arr=dlog(eq_arr)
       d1=(loga_table(2)-loga_table(1))/(eq_arr(2)-eq_arr(1))
       d2=(loga_table(ntable)-loga_table(ntable-1))/(eq_arr(ntable)-eq_arr(ntable-1))      
       call spline(eq_arr(1:ntable),(loga_table(1:ntable)),ntable,d1,d2,eq_arr_buff(1:ntable))
       call spline_out(eq_arr(1:ntable),loga_table(1:ntable),&
            &eq_arr_buff(1:ntable)&
            &,ntable,0.0d0,Params%aeq)
       Params%aeq=10._dl**(Params%aeq) !RL fixed 012524
       !Sometimes this spline breaks if a_osc<a_eq, in that case simpler expressions can be used
       regzeq=(Params%omegah2_rad+sum(lhsqcont_massive))/(omegah2_b+omegah2_dm+omegah2_ax)
       aeq_LCDM = (((Params%TCMB**4.0d0)/(rhocrit))/(c**2.0d0)*a_rad*1.d1/(1.d4))*&
            &(1._dl + (Params%nu_massless_degeneracy + &
            &sum(Params%Nu_mass_degeneracies(1:Params%Nu_mass_eigenstates)))*(7._dl/8._dl)*&
            &((4._dl/11._dl)**(4._dl/3._dl)))/(omegah2_b+omegah2_dm+omegah2_ax) !RL 031724, the approximate aeq assuming a matter scaling of axions, used for phase-finding
       !print*, 'RL, omegah2_b+omegah2_dm+omegah2_ax + (((Params%TCMB**4.0d0)/rhocrit)/(c**2.0d0)*a_rad*1.d1/(1.d4)) *sum(Params%Nu_mass_degeneracies(1:Params%Nu_mass_eigenstates)) * (7._dl/8._dl)*((4._dl/11._dl)**(4._dl/3._dl))', omegah2_b+omegah2_dm+omegah2_ax + (((Params%TCMB**4.0d0)/rhocrit)/(c**2.0d0)*a_rad*1.d1/(1.d4)) *sum(Params%Nu_mass_degeneracies(1:Params%Nu_mass_eigenstates)) * (7._dl/8._dl)*((4._dl/11._dl)**(4._dl/3._dl))
       
       !!write(*, *) 'regzeq, aeq_LCDM, their fractional difference', regzeq, aeq_LCDM, regzeq/aeq_LCDM - 1._dl
       if (Params%a_osc.lt.regzeq) then
          Params%aeq=regzeq
       endif
!!!!!!!
       deallocate(eq_arr, eq_arr_buff) !RL 111623
       


       !!feed out real (dl) valued version of all this stuff for output
       !put scalar field initial condition in planck units
       !Params%phiinit=vtwiddle_init*sqrt(6.0d0)
       !print*, 'vtwiddle_init', vtwiddle_init
       !if (Params%use_axfrac) then
       !   Params%axfrac = Params%axfrac
       !else
       !   Params%axfrac = fax
       !end if

       !The isocurvature amplitude amp_i is built in inidriver_axion.F90, after this routine
       !returns, where phinorm_table(1) supplies v1 and all three parameterizations are
       !handled together. The background solver no longer knows about isocurvature amplitudes.

       !The converged initial field value phi_ini/Mpl is reported once from
       !inidriver_axion.F90, after the abundance iteration. Reporting it here would print
       !one line per call, and this routine is invoked repeatedly while movH_switch is
       !re-solved.

       !RL 070926: frozen-field initial condition requires m/H < 0.1 at the start of
       !the background integration (paper Sec. "Isocurvature parameters and interface").
       !The same threshold is enforced at the perturbation start in GetTauStart (cmbmain).
       if (maxion_twiddle * a_arr(1) * hnot / littlehfunc(1) .ge. 0.1_dl) then
          write(*,*) 'Error: m_ax/H_ini =', maxion_twiddle*a_arr(1)*hnot/littlehfunc(1), &
               ' >= 0.1; the frozen-field initial condition is violated.'
          stop
       end if

       !output omega_r
       Params%omegar=Params%omegah2_rad/hsq
       
       ! Write full background: KG table up to a_osc, then analytic a^-3


       deallocate(a_arr, v_vec, littlehfunc, littlehfunc_buff, rhoaxh2_ov_rhom)
     end subroutine w_evolve

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!
     ! Begin derivative routine
!!!!!!!!!!!!!!!!!!!!!!!!

     subroutine derivs_bg(a,v,dvt_dloga,omegah2_regm,omegah2_rad,omegah2_lambda,omk,hsq,maxion_twiddle,badflag,&
          &lhsqcont_massless,lhsqcont_massive,Nu_mass_eigenstates,Nu_masses)
       use constants
       use Precision
       implicit none
       integer badflag,Nu_mass_eigenstates
       real(dl) a
       real(dl) dvt_dloga(1:2),dvt_da(1:2),lhr
       real(dl) v(1:2)
       real(dl) omegah2_regm,omegah2_rad,omegah2_lambda
       real(dl) maxion_twiddle,Nu_masses(Nu_mass_eigenstates)
       real(dl) omk,hsq,lhsqcont_massless,lhsqcont_massive(Nu_mass_eigenstates)

       !write(*, *) 'derivs_bg in the background is called' !RL testing
       !Compute hubble dimensionless 
       call lh(omegah2_regm,omegah2_rad,omegah2_lambda,omk,hsq,maxion_twiddle,a,v,lhr,badflag,&
            &lhsqcont_massless,lhsqcont_massive,Nu_mass_eigenstates,Nu_masses)

       !calculate fluid derivatives (d/dloga )for next step
       !Solving Equation described and defined in top part of this fortran file
       dvt_da(1)=v(2)*dsqrt(hsq)/(a*(lhr)) !RL added the solution proposed by astralsight5 *dsqrt(hsq)
       dvt_da(2)=-2.0d0*v(2)/(a)-(maxion_twiddle**2.0d0)*a*v(1)*dsqrt(hsq)/(lhr)! RL added the solution proposed by astralsight5 *dsqrt(hsq)
       dvt_dloga(1:2)=a*dvt_da(1:2)
     end subroutine derivs_bg

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !RL 120624: optional subroutine: obtain the analytical approximation of the mass oscillation phase
     subroutine get_phase_info(Params, y_beta, beta_coeff, movHETA, beta2x)
       use ModelParams
       use constants
       implicit none
       type(CAMBparams) :: Params
       real(dl) y_beta, movHETA, beta_coeff, beta2x
       y_beta = Params%a_osc/aeq_LCDM
       movHETA = Params%dfac*Params%ah_osc/Params%ahosc_ETA
       beta_coeff = (4._dl*(y_beta**2 - y_beta - 2.0_dl + 2.0_dl*sqrt(1.0_dl + y_beta)))/(3._dl*(y_beta**2))
       beta2x = movHETA*beta_coeff - const_pi*3._dl*(1.0_dl + y_beta)/(4.0_dl + 3.0_dl*y_beta)
     
     end subroutine get_phase_info
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !RL added: supply the initial conditions at aosc 
     subroutine auxiIC(Params, omegah2_regm,omegah2_rad,omegah2_lambda,omk,hnot,maxion_twiddle,&
          &a, v, badflag,lhsqcont_massless,lhsqcont_massive,Nu_mass_eigenstates,Numasses, &
          &littlehauxi, lhETA, A_coeff, tvarphi_c, tvarphi_cp, tvarphi_s, tvarphi_sp, rhorefp, Prefp)
       use ModelParams
       use constants
       use Precision
       use MassiveNu
       implicit none
       type(CAMBparams) :: Params 
       integer badflag,Nu_mass_eigenstates,i,iter_EFA,maxiter !RL 041024 added iter_EFA
       real(dl) omegah2_regm,omegah2_rad,omegah2_lambda,omk,hnot,maxion_twiddle,a, a2
       real(dl) v(1:2),littlehauxi,lhETA, lhETA_upd, lhsqcont_massless,lhsqcont_massive(Nu_mass_eigenstates)
       real(dl) mass_correctors(Nu_mass_eigenstates), w_nu(Nu_mass_eigenstates), rhonu, pnu 
       real(dl) numasses(Nu_mass_eigenstates)
       real(dl) w_ax, wEFA_c_upd !RL
       real(dl) Hovm_ins, Hovm_ETA, dHsqdmt_term, A_coeff, A_denom, tvarphi_c, tvarphi_cp, tvarphi_s, tvarphi_sp, rhorefp, Prefp !RL added for outputting the initial conditions
       real(dl) tol_EFA !RL 041024: tolerance to determine that <H> and wEFA have converged
       !!real(dl) H_eV !RL 110823
       !Hubble parameter in eV
       !!H_eV=1.d14*6.5821d-25*hnot/(MPC_in_sec*c)
       !At this point, first compute dimensionless *conformal* hubble including standard matter and axion potential and kinetic energy
       call lh(omegah2_regm,omegah2_rad,omegah2_lambda,omk,hnot**2.0d0,maxion_twiddle,a,v,&
            &littlehauxi,badflag, lhsqcont_massless,lhsqcont_massive,Nu_mass_eigenstates,Numasses)
       !RL 041024: initial setups - we don't have rhorefp and Prefp yet, so the first step is using H(instantaneous) to construct EFA quantities
       lhETA = littlehauxi
       !write(*, *) 'In auxiIC, initial lhETA = littlehauxi', lhETA
       !write(*, *) 'Rayne, littleh at refp, H_refp in eV, m/H_refp', littlehauxi, 1.d14*6.5821d-25*littleh_ref/(Params%a_osc*MPC_in_sec*c), Params%ma/(1.d14*6.5821d-25*littlehauxi/(Params%a_osc*MPC_in_sec*c)) !(RL tested, the m/H is right at this point)
       !I also need dH^2/d(mt), which is evaluated using dH^2/dlna, see note1_stepplanning in EFA notes for details. Since I still need the energy density and pressure of all cosmological components, the following code referred heavily to lh
       !Compute corrections to neutrino energy density and pressure from massless to massive case
       do i=1,Nu_mass_eigenstates,1
          call Nu_background(a*numasses(i),rhonu,pnu)
          mass_correctors(i)=rhonu
          w_nu(i) = pnu/rhonu
       enddo
       !write(*, *) 'Rayne, i, w_nu(i), w_nu, w_nu(1:Nu_mass_eigenstates), (1.0d0 + w_nu(i)), (1.0d0 + w_nu)', i, w_nu(i), w_nu, w_nu(1:Nu_mass_eigenstates), (1.0d0 + w_nu(i)), (1.0d0 + w_nu)
       a2 = a**2.0d0
       !RL 041024 - initialize rhorefp for computing dH^2/dtau^2. The initial value is the instantaneous energy density ratio
       rhorefp = (v(2)/a)**2.0d0+(maxion_twiddle*v(1))**2.0d0
       tol_EFA = 1.e-7_dl
       maxiter = 30 !Maximum iteration for the <H> convergence
       Hovm_ins = littlehauxi/(a*hnot*maxion_twiddle)
       !Start the iteration process.
       do iter_EFA = 1, maxiter, 1
          !H/m = (H/H0)(H0/m) = (lh/(a*hnot))/m_twiddle  
          Hovm_ETA = lhETA/(a*hnot*maxion_twiddle)
          !Compute the eos of axions since it's a bit long
          !w_ax_p1 = (2.0d0*(v(2)**2.0d0)/a2)/((v(2)**2.0d0)/a2 + (maxion_twiddle*v(1))**2.0d0)
          !With the eos correction from 0
          !wcorr_coeff = littlehauxi*a/(maxion_twiddle*hnot)

          w_ax = (Hovm_ETA**2.0d0)*Params%wEFA_c !RL 110923

          !Sum up all the components together first - multiply by HH0^2/m later. The reason why I don't take out the -3.0 is I don't want to compute fractions for EoS for example 1/3. Note this expression carries an h^2, so we don't need to multiply by another hnot^2 afterwards
          dHsqdmt_term = omegah2_regm*(-3.0d0)/a+& !matter (a^-3*a^2, explained below)
               &omegah2_rad*(-4.0d0)/(a2)+& !radiation (a^-4*a^2, again explained below)
               &sum(lhsqcont_massive*mass_correctors*(-3.0d0*(1.0d0 + w_nu)),Nu_mass_eigenstates)/(a2)+&!massive neutrino correction (a^-4*a^2) 
                                !&omegah2_lambda+& !dark energy has w=-1 (?? Do we do anything with time-changing w_a? Maybe not according to how the axion_background lh code was done...)
               &rhorefp*(-3.0d0*(w_ax + 1.0d0))*(a2)+& !axions - note what's computed here is the actual omaxh2(a), not at today (a^2 again explained below). RL 041024: this should also apply when we substitute <H> for the dH^2/dmt term since the continuity equation doesn't rely on w being a constant
               &omk*(hnot**2.0d0)*(-2.0d0) !curvature (a^-2*a^2)

          !Now multiply by HH0^2/m. This multiplication carries all the dimension. But better yet we also have a m/H^3 term, so this overall gives H0^2/H^2 = (H0/H)^2. Since lh = ah, H0/H = hnot/h = a(hnot/lh), (H0/H)^2 = a^2 (hnot/lh)^2. This cancels out the a^2 term we evaluated before so we simply remove them from above. Also as explained above the h^2 term is already included so we just divide by lh^2
          dHsqdmt_term = dHsqdmt_term/(lhETA**2.0d0) !Now this term is dimensionless
          !Now we can compute the coefficient that links the 2nd order derivatives wrt (mt) to the 1st derivatives.
          A_coeff = (-Hovm_ETA/2.0d0)*(3.0d0 - dHsqdmt_term)
          !Also a denominator term in the IC terms appear more than once
          A_denom = A_coeff**2.0d0 + 3.0d0*A_coeff*Hovm_ins + 4.0d0
          !Now construct the 4 IC sin & cos terms (with normalization sqrt(4\pi G/3)h)
          tvarphi_c = v(1)
          tvarphi_cp = -3.0d0*Hovm_ins*(2.0d0*v(1) + (A_coeff + 3.0d0*Hovm_ins)*v(2)/(a*maxion_twiddle))/A_denom
          tvarphi_s = v(2)/(a*maxion_twiddle) - tvarphi_cp
          tvarphi_sp = 3.0d0*Hovm_ins*(A_coeff*v(1) - 2.0d0*v(2)/(a*maxion_twiddle))/A_denom
          !rhorefp now is the effective fluid version
          rhorefp = (maxion_twiddle**2.0d0)*(tvarphi_c**2.0d0 + tvarphi_s**2.0d0 + &
               &(tvarphi_cp**2.0d0 + tvarphi_sp**2.0d0)/2.0d0 - tvarphi_c*tvarphi_sp + tvarphi_s*tvarphi_cp)
          !also the effective fluid version of the pressure
          Prefp = (maxion_twiddle**2.0d0)*(tvarphi_cp**2.0d0/2.0d0 + tvarphi_sp**2.0d0/2.0d0 -&
               & tvarphi_c*tvarphi_sp + tvarphi_s*tvarphi_cp)

          !Compare the new H and wEFA coefficient
          wEFA_c_upd = (Prefp/rhorefp)/((lhETA/(maxion_twiddle*hnot*a))**2._dl)
          call lh(omegah2_regm,omegah2_rad,omegah2_lambda,omk,hnot**2.0d0,maxion_twiddle,&
               &a,v,lhETA_upd,badflag, lhsqcont_massless,lhsqcont_massive,Nu_mass_eigenstates,Numasses,rhorefp)
          !!write(*, *) 'auxiIC iteration, wEFA_c_upd, Params%wEFA_c, abs(wEFA_c_upd/Params%wEFA_c -1.0_dl), littlehauxi_upd, littlehauxi, littlehauxi_upd/littlehauxi -1.0_dl'
          !!   write(*, *) wEFA_c_upd, Params%wEFA_c, abs(wEFA_c_upd/Params%wEFA_c -1.0_dl), littlehauxi_upd, littlehauxi, littlehauxi_upd/littlehauxi -1.0_dl
          !if the new quantities are sufficiently close to the previous quantity, end the loop; else reassign, and loop again
          if (abs(wEFA_c_upd/Params%wEFA_c -1.0_dl) .lt. tol_EFA) then !!!-  .and. abs(lhETA_upd/lhETA -1.0_dl) .lt. tol_EFA) then
             !RL 041128 - note that we still need to get the instantaneous ah, or else this algorithm is outputting the <H> as the instantaneous ah!!

             !! write(*, *) 'auxiIC iteration converged, lhETA_upd, lhETA, lhETA_upd/lhETA -1.0_dl, instantaneous ah'
             !!write(*, *) lhETA_upd, lhETA, lhETA_upd/lhETA -1.0_dl, littlehauxi
             exit
          else
             Params%wEFA_c = wEFA_c_upd
             lhETA = lhETA_upd
          end if
       end do
       !write(*, *) 'Converged Params%wEFA_c, lhETA', Params%wEFA_c, lhETA

       if (iter_EFA .eq. maxiter) then !Warning sign if iter_EFA exceeds maximum iteration
          print '(A, I0, A, E13.5E3, A, E13.5E3)', &
               &'Warning: exceeding the maximum number of iteration (', iter_EFA, ') for bisection: fractional error in w = ', &
               &wEFA_c_upd/Params%wEFA_c -1.0_dl, ', fractional error in H_ETA = ', lhETA_upd/lhETA -1.0_dl
          !write(*, *) wEFA_c_upd, Params%wEFA_c, lhETA_upd, lhETA, lhETA_upd/lhETA -1.0_dl
          !RL: note we didn't go ahead and assign the vtwiddle_init and corresponding aosc, hence the code will crash if this error appears (tested with a very small nphi). But under normal circumstances this should not happen
       end if

     end subroutine auxiIC
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !Compute aH/(100 km/s/Mpc) at any moment in our cosmological history !RL changed to aH to correspond to what axionCAMB has actually done
     subroutine lh(omegah2_regm,omegah2_rad,omegah2_lambda,omk,hsq,maxion_twiddle,a,v,littlehfunc,&
          &badflag,&
          &lhsqcont_massless,lhsqcont_massive,Nu_mass_eigenstates,nu_masses,rho_f)
       use constants
       use Precision
       use MassiveNu
       implicit none
       integer badflag,Nu_mass_eigenstates,i
       real(dl) omegah2_regm,omegah2_rad,omegah2_lambda,omk,hsq,maxion_twiddle,a
       real(dl) v(1:2),littlehfunc,lhsqcont_massless,lhsqcont_massive(Nu_mass_eigenstates)
       real(dl) mass_correctors(Nu_mass_eigenstates),rhonu
       real(dl) nu_masses(Nu_mass_eigenstates)
       real(dl), optional :: rho_f !RL 113023 added an optional input for the time average energy density

       !Compute corrections to neutrino energy density from massless to massive case
       do i=1,Nu_mass_eigenstates,1
          call Nu_rho(a*nu_masses(i),rhonu)
          mass_correctors(i)=rhonu
       enddo
!!!!!!!


       !Compute H/(100 km /s/Mpc), contributions from normal stuff
       littlehfunc=(omegah2_regm/(a**3.0d0)+omegah2_rad/(a**4.0d0))+&
            &sum(lhsqcont_massive*mass_correctors,Nu_mass_eigenstates)/(a**4.0d0)

       !Contributions from cosmological constant and scalar field
       littlehfunc=littlehfunc+omegah2_lambda
       if (present(rho_f)) then !RL 113023
!!!write(*, *) 'Rayne, in lh, making sure v(1) and v(2) are not used'
          littlehfunc=littlehfunc+rho_f
       else
          littlehfunc=littlehfunc+(maxion_twiddle*v(1))**2.0d0+((v(2)/a)**2.0d0)
       end if

       littlehfunc=littlehfunc*(a**2.0d0)+omk*hsq !RL: note this line is the line that changes Hubble to conformal Hubble
       ! DM flag histories where h goes negative, i.e. do not want collapsing universe
       if (littlehfunc .le. 0.0d0) then
          badflag=1
       endif
       !if (isnan(littlehfunc)) then !RL 010925 isnan is not available here, change to old NaN check methods
       if (.not. littlehfunc == littlehfunc) then
          badflag = 1
       endif
       !

       littlehfunc=dsqrt(littlehfunc)

     end subroutine lh


!!!!! Take next step in integration using method described above

     subroutine next_step(a,v,kvec,kfinal,avec,omegah2_regm,omegah2_rad,omegah2_lambda,omk,&
          &hsq,maxion_twiddle,badflag,dloga,nstep,cmat,&
          &lhsqcont_massless,lhsqcont_massive,Nu_mass_eigenstates,Nu_masses)
       use constants
       use Precision
       implicit none
       integer nstep,cp,m,badflag
       real(dl) hsq,a,v(1:2),kvec(1:2,1:nstep),omegah2_regm,omegah2_rad,omegah2_lambda,maxion_twiddle,dloga,omk
       real(dl) vfeed(1:2),cmat(1:nstep,1:nstep),kfinal(1:2),avec(1:nstep)
       integer Nu_mass_eigenstates
       real(dl) lhsqcont_massless,lhsqcont_massive(Nu_mass_eigenstates)
       real(dl) Nu_masses(Nu_mass_eigenstates)

       kvec=0.0d0
       kfinal=0.0d0
       call derivs_bg(a,v(1:2),kvec(1:2,1),omegah2_regm,omegah2_rad,omegah2_lambda,omk,hsq,maxion_twiddle,badflag,&
            &lhsqcont_massless,lhsqcont_massive,Nu_mass_eigenstates,Nu_masses)
       kvec(1:2,1)=kvec(1:2,1)*dloga
       do m=1,nstep,1
          do cp=1,2,1
             vfeed(cp) = dot_product(cmat(m, 1:m), kvec(cp, 1:m)) !RL 011025 reverse kvec order

          enddo
          if (m .le. (nstep-1)) then
             call derivs_bg(a*dexp(dloga*avec(m)),v(1:2)+vfeed(1:2),&
                  &kvec(1:2,m+1),omegah2_regm,omegah2_rad,omegah2_lambda,omk,hsq&
                  &,maxion_twiddle,badflag,lhsqcont_massless,lhsqcont_massive,Nu_mass_eigenstates,Nu_masses)
             kvec(1:2,m+1)=kvec(1:2,m+1)*dloga
          else
             call derivs_bg(a*dexp(dloga*avec(m)),v(1:2)+vfeed(1:2),&
                  &kfinal(1:2),omegah2_regm,omegah2_rad,omegah2_lambda,omk,hsq&
                  &,maxion_twiddle,badflag,lhsqcont_massless,lhsqcont_massive,Nu_mass_eigenstates,Nu_masses)
             kfinal(1:2)=kfinal(1:2)*dloga
          endif
       enddo
       !write(*, *) 'Rayne, testing whether v1, 2 solves KG better than 1e-4 in next_step'
       !write(*, *) '1st term:', kvec(1,1:2)
     end subroutine next_step

!!!!!!!!
     

   end module axion_background
