!     Code for Anisotropies in the Microwave Background
!     by Antony Lewis (http://cosmologist.info/) and Anthony Challinor
!     See readme.html for documentation. This is a sample driver routine that reads
!     in one set of parameters and produdes the corresponding output.

program driver
  use IniFile
  use CAMB
  use LambdaGeneral
  use Lensing
  use AMLUtils
  use Transfer
  use constants
  use Precision
  use Bispectrum
  use CAMBmain
  use NonLinear ! RH made this change for axions

  !Tell CAMB to solve the Klein-Gordon equation for a background axion field
  !to determine what initial condition is needed to get desired relic abundance of axions
  !today

  use axion_background
#ifdef NAGF95
  use F90_UNIX
#endif
  implicit none

  Type(CAMBparams) P
  character(LEN=Ini_max_string_len) numstr, VectorFileName, &
       InputFile, ScalarFileName, TensorFileName, TotalFileName, LensedFileName,&
       LensedTotFileName, LensPotentialFileName,ScalarCovFileName
  integer i
  character(LEN=Ini_max_string_len) TransferFileNames(max_transfer_redshifts), &
       MatterPowerFileNames(max_transfer_redshifts), outroot, version_check
  real(dl) output_factor, nmassive,omnuh2,nu_massless_degeneracy,fractional_number
  real(dl) actual_massless,neff_i
  real clock_start, clock_stop ! RH timing RL reusing 07/18/2023	

  type (CAMBdata)  :: AxionIsoData ! Adding this for the iso stuff
  type (CAMBdata)  :: AxionAdiData ! Adding this for the iso stuff
#ifdef WRITE_FITS
  character(LEN=Ini_max_string_len) FITSfilename
#endif

  logical bad

  !Cosmo parameters for dg integrator
  !dimensionless Hubble
  real(dl) hnot
  !a_equality, omega_radiation, initial scalar field value,omega_rad H^2, rho_crit, number of massless neutrinos
  integer iter_dfac, iter_dfacETA
!!!real(dl) dfac_prev1, dfac_prev2, aosc_prev1, aosc_prev2 !RL adding to perform bisection for the recombination jump case
  real(dl) aeq,omegar,phiinit,omegah2_rad,rhocrit, nnu, rh_num_nu_massless
  !Timing variables
  real clock_totstart, clock_totstop ! RH timing
  integer reni ! RH
  !Control Flag
  integer badflag
  real(dl), allocatable :: RHCl_temp(:,:,:), RHCl_temp_lensed(:,:,:), RHCl_temp_tensor(:,:,:) !RL 111523 moving the RHCl_temp arrays here (might change later) to eliminate the stack overflow problem
  !!call cpu_time(clock_totstart) ! RH timing

  real(dl) twobeta_tgt, beta_coeff, y_phase, movHETA_beta, movHETA_new
  real(dl) twobeta_new, twobeta_old1, twobeta_old2, hETA_beta, hETA_new
  real(dl) hosc_new, hosc_old1, hosc_old2, hETA_old1, hETA_old2, beta_tol !RL 030624
  character(LEN=Ini_max_string_len) filenametest !RL 042824
  ! End axion stuff

  character(len=500) :: testthetafilename !RL022224 temporary
  character(len=32) :: dfac_str !RL022224 temporary
  character(len=6) :: abun_str, abun_str1 !RL022224 temporary
  character(len=5) :: H0_str, H0_str1 !RL 091924 temporary
  character(len=32) :: max_str !RL022224 temporary
  character(len=4) :: omk_str, omk_str1 !RL090324 temporary
  character(len=1) :: useaxfrac_str, dolateradtrunc_str !RL022224 temporary


  InputFile = ''
  if (GetParamCount() /= 0)  InputFile = GetParam(1)
  if (InputFile == '') stop 'No parameter input file'

  call Ini_Open(InputFile, 1, bad, .false.)
  if (bad) stop 'Error opening parameter file'

  Ini_fail_on_not_found = .false.

  outroot = Ini_Read_String('output_root')
  if (outroot /= '') outroot = trim(outroot) // '_'

  highL_unlensed_cl_template = Ini_Read_String_Default('highL_unlensed_cl_template',highL_unlensed_cl_template)

  call CAMB_SetDefParams(P)

  P%WantScalars = Ini_Read_Logical('get_scalar_cls')
  P%WantVectors = Ini_Read_Logical('get_vector_cls',.false.)
  P%WantTensors = Ini_Read_Logical('get_tensor_cls',.false.)

  P%OutputNormalization=outNone
  output_factor = Ini_Read_Double('CMB_outputscale',1.d0)

  P%WantCls= P%WantScalars .or. P%WantTensors .or. P%WantVectors

  P%PK_WantTransfer=Ini_Read_Logical('get_transfer')

  AccuracyBoost  = Ini_Read_Double('accuracy_boost',AccuracyBoost)
  lAccuracyBoost = Ini_Read_Real('l_accuracy_boost',lAccuracyBoost)
  HighAccuracyDefault = Ini_Read_Logical('high_accuracy_default',HighAccuracyDefault)

  P%NonLinear = Ini_Read_Int('do_nonlinear',NonLinear_none)
  if (P%NonLinear == 1 .or. P%NonLinear == 2) then
     !write(*, *) 'do_nonlinear is temporarily disabled in this version, setting back to 0.'
     !P%NonLinear = 0
     write(*, *) 'Warning: do_nonlinear options are not well-tested in this version'
  end if

  P%DoLensing = .false.
  if (P%WantCls) then
     if (P%WantScalars  .or. P%WantVectors) then
        P%Max_l = Ini_Read_Int('l_max_scalar')
        P%Max_eta_k = Ini_Read_Double('k_eta_max_scalar',P%Max_l*2._dl)
        if (P%WantScalars) then
           P%DoLensing = Ini_Read_Logical('do_lensing',.false.)
           if (P%DoLensing) lensing_method = Ini_Read_Int('lensing_method',1)
        end if
        if (P%WantVectors) then
           if (P%WantScalars .or. P%WantTensors) stop 'Must generate vector modes on their own'
           i = Ini_Read_Int('vector_mode')
           if (i==0) then
              vec_sig0 = 1
              Magnetic = 0
           else if (i==1) then
              Magnetic = -1
              vec_sig0 = 0
           else
              stop 'vector_mode must be 0 (regular) or 1 (magnetic)'
           end if
        end if
     end if

     if (P%WantTensors) then
        P%Max_l_tensor = Ini_Read_Int('l_max_tensor')
        P%Max_eta_k_tensor =  Ini_Read_Double('k_eta_max_tensor',Max(500._dl,P%Max_l_tensor*2._dl))
     end if
  endif



  !  Read initial parameters.

  call DarkEnergy_ReadParams(DefIni)

  P%H0     = Ini_Read_Double('hubble')
!!!write(H0_str1,'(F5.2)') P%H0
  !write(*, *) 'abun_str T', abun_str
!!!write(H0_str, '(A,A,A)') H0_str1(1:2), 'd', H0_str1(4:5)
  P%H0_in_Mpc_inv = dble(P%H0)/dble(c/1.0d3) !RL
  P%H0_eV = h_P*P%H0_in_Mpc_inv/(Mpc_in_sec*2._dl*const_pi*elecV) !RL

  P%omegak= Ini_Read_Double('omk')
  !RL 120124: I need radiation fraction for assigning the DE fraction of axions, which in turn needs the neutrino fraction. Hence I have to take out this assignment separately
  if (Ini_Read_Logical('use_physical',.false.)) then 
     P%omegan = Ini_Read_Double('omnuh2')/(P%H0/100)**2
  else
     P%omegan = Ini_Read_Double('omega_neutrino')
  end if

  !!P%dfac = Ini_Read_Double('movH_switch') !RL 092623 switch time
  P%dfac = 10._dl !RL 121924 making movH internal
  ntable = nint(P%dfac*100) + 1 !RL 111123

  P%tcmb   = Ini_Read_Double('temp_cmb',COBE_CMBTemp)
  P%yhe    = Ini_Read_Double('helium_fraction',0.24_dl)

  !Compute  some basic constants
  rhocrit=(8.0d0*const_pi*G*1.d3/(3.0d0*((1.d7/(MPC_in_sec*c*1.d2))**(2.0d0))))**(-1.0d0)
  omegah2_rad=((P%tcmb**4.0d0)/(rhocrit))/(c**2.0d0) !RL replaced the COBE temperature
  omegah2_rad=omegah2_rad*a_rad*1.d1/(1.d4)

  !DG May 25 2015
  !Neutrino stuff out of usual order so that OmegaK Can be self consistently computed
  !If you restructure this make sure that P%omegav is self-consistently computed including massive and massless
  P%Num_Nu_massless  = Ini_Read_Double('massless_neutrinos')

  P%Nu_mass_eigenstates = Ini_Read_Int('nu_mass_eigenstates',1)
  if (P%Nu_mass_eigenstates > max_nu) stop 'too many mass eigenstates'

  numstr = Ini_Read_String('massive_neutrinos')
  read(numstr, *) nmassive

  if (abs(nmassive-nint(nmassive))>1e-6) stop 'massive_neutrinos should now be integer (or integer array)'
  read(numstr,*, end=100, err=100) P%Nu_Mass_numbers(1:P%Nu_mass_eigenstates)
  P%Num_Nu_massive = sum(P%Nu_Mass_numbers(1:P%Nu_mass_eigenstates))

  if (P%Num_Nu_massive>0) then
     P%share_delta_neff = Ini_Read_Logical('share_delta_neff', .true.)
     numstr = Ini_Read_String('nu_mass_degeneracies')

     if (P%share_delta_neff) then
        if (numstr/='') write (*,*) 'WARNING: nu_mass_degeneracies ignored when share_delta_neff'
     else
        if (numstr=='') stop 'must give degeneracies for each eigenstate if share_delta_neff=F'
        read(numstr,*) P%Nu_mass_degeneracies(1:P%Nu_mass_eigenstates)
     end if

     numstr = Ini_Read_String('nu_mass_fractions')
     if (numstr=='') then
        if (P%Nu_mass_eigenstates >1) stop 'must give nu_mass_fractions for the eigenstates'
        P%Nu_mass_fractions(1)=1
     else
        read(numstr,*) P%Nu_mass_fractions(1:P%Nu_mass_eigenstates)
     end if
  end if

  !RL moved here 020625, tested Cl identical before and after the suite of changes
  if (P%Num_Nu_Massive /= sum(P%Nu_mass_numbers(1:P%Nu_mass_eigenstates))) then
     if (sum(P%Nu_mass_numbers(1:P%Nu_mass_eigenstates))/=0) stop 'Num_Nu_Massive is not sum of Nu_mass_numbers'
  end if

  if (P%Omegan == 0 .and. P%Num_Nu_Massive /=0) then
     print*, 'omeganuh2 is set to 0 but we still have massive neutrinos, resetting'
     if (P%share_delta_neff) then
        P%Num_Nu_Massless = P%Num_Nu_Massless + P%Num_Nu_Massive
     else
        P%Num_Nu_Massless = P%Num_Nu_Massless + sum(P%Nu_mass_degeneracies(1:P%Nu_mass_eigenstates))
     end if
     P%Num_Nu_Massive  = 0
     P%Nu_mass_numbers = 0
  end if
  P%Nu_massless_degeneracy = P%Num_Nu_massless !N_eff for massless neutrinos !RL 061425

!!!4/8 DG Error in original AxionCAMB
!!! massless neutrino contribution wrong
!!!When neutrinos are massive
  !!correcting H contributions already
  if (P%Num_nu_massive > 0) then
     if (P%Nu_mass_eigenstates==0) stop 'Have Num_nu_massive>0 but no nu_mass_eigenstates'
     if (P%Nu_mass_eigenstates==1 .and. P%Nu_mass_numbers(1)==0) P%Nu_mass_numbers(1) = P%Num_Nu_Massive
     if (all(P%Nu_mass_numbers(1:P%Nu_mass_eigenstates)==0)) P%Nu_mass_numbers=1 !just assume one for all
     if (P%share_delta_neff) then
        !default case of equal heating of all neutrinos
        fractional_number = P%Num_Nu_massless + P%Num_Nu_massive
        actual_massless = int(P%Num_Nu_massless + 1e-6_dl)
        neff_i = fractional_number/(actual_massless + P%Num_Nu_massive)
        nu_massless_degeneracy = neff_i*actual_massless
        P%Nu_massless_degeneracy=nu_massless_degeneracy
        P%Nu_mass_degeneracies(1:P%Nu_mass_eigenstates) = P%Nu_mass_numbers(1:P%Nu_mass_eigenstates)*neff_i
     end if
     if (abs(sum(P%Nu_mass_fractions(1:P%Nu_mass_eigenstates))-1) > 1e-4) &
          stop 'Nu_mass_fractions do not add up to 1'
  else
     P%Nu_mass_eigenstates = 0
  end if

  grhog= ((kappa/(c**2.0d0)*4.0d0*sigma_boltz)/(c**3.0d0))*(P%tcmb**4.0d0)*(Mpc**2.0d0) !RL replaced the COBE temperature
  P%grhor = (7.0d0/8.0d0)*((4.0d0/11.0d0)**(4.0d0/3.0d0))*grhog

  ! RH added this here because we are calculating omegav - but not trying to change anything globally

  if (P%omegan == 0 .and. P%Num_Nu_Massive /=0) then 
     !        print*, 'we are where omeganuh2=0 but we still have massive neutrinos'
     !	print*, P%Num_Nu_Massless, P%Num_Nu_Massive, 'here 1'
     if (P%share_delta_neff) then
        rh_num_Nu_Massless = P%Num_Nu_Massless + P%Num_Nu_Massive
     else 
        rh_Num_Nu_Massless = P%Num_Nu_Massless + sum(P%Nu_mass_degeneracies(1:P%Nu_mass_eigenstates)) 
     end if

     ! Note that this will be computed correctly later - this is not carrying through, but only for omegak
     !        CP%Num_Nu_Massive  = 0 
     !        CP%Nu_mass_numbers = 0 
  end if

  if (P%omegan ==0 .and. P%Num_Nu_Massive ==0)  rh_num_nu_massless = P%Num_Nu_Massless 

  ! RH again calculating this here so we make sure to add the neutrino density correctly to get omegak

  if (P%omegan > 0 .and. P%Num_Nu_massive > 0) then
     rh_num_nu_massless = P%Num_Nu_Massless ! only using the massless neutrinos, and adding in omnuh2 later
  end if


  !     print*, P%Num_Nu_Massless, P%Num_Nu_Massive, rh_num_nu_massless, 'here 2'
  omegah2_rad= omegah2_rad+(rh_Num_Nu_massless*P%grhor*(c**2.0d0)/((1.d5**2.0d0)))/3.0d0     
  P%omegah2_rad = omegah2_rad

  if (Ini_Read_Logical('use_physical',.false.)) then 
     P%omegab = Ini_Read_Double('ombh2')/(P%H0/100)**2

     P%use_axfrac = Ini_Read_Logical('use_axfrac',.false.)
     ! read in axion mass
     P%ma     = Ini_Read_Double('m_ax') !! RH axion mass
     if (P%ma < 0) P%ma = 10**P%ma ! RH making this exponential from the inidriver
     P%m_ovH0 = P%ma/P%H0_eV !RL 050724

     if (P%use_axfrac) then
        !! Compute axion fractions rather than densities
        P%omegada = Ini_Read_Double('omdah2')/(P%H0/100)**2
        ! Read in Axion faction and compute density
        P%axfrac = Ini_Read_Double('axfrac')
        if (P%m_ovH0 .ge. 10._dl) then !RL 120124 - DM case

           P%omegaax = P%axfrac*P%omegada 
           P%omegac = (1-P%axfrac)*P%omegada
        else !RL 120124 - DE case
           write(*, *) 'Note: m/H0 < 10, axfrac is ULA fraction in dark energy'
!!!write(*,*) 'P%omegac*h2 before assignment', P%omegac*((P%H0/1.d2)**2.0d0)
           P%omegac = P%omegada
           P%omegaax = P%axfrac*(1._dl-P%omegab-P%omegac - P%omegan -P%omegak - P%omegah2_rad/((P%H0/1.d2)**2.0d0))
!!!write(*,*) 'P%axfrac, P%omegaax', P%axfrac, P%omegaax
        end if

     else 
        ! read in axion densities and matter densities
        P%omegac = Ini_Read_Double('omch2')/(P%H0/100)**2
        P%omegaax = Ini_Read_Double('omaxh2')/(P%H0/100)**2
        if (P%m_ovH0 .ge. 10._dl) then !RL 120124 - DM case 
           P%axfrac = P%omegaax/(P%omegac+P%omegaax)
        else !RL 120124 - DE case
           P%axfrac = P%omegaax/(1.0d0-P%omegab-P%omegac - P%omegan -P%omegak - P%omegah2_rad/((P%H0/1.d2)**2.0d0))
        end if

     endif

!!!!!!
     !Compute value of cosmological constant including curvature and radiation (photons + massless neutrinos) self consistently	
     P%omegav = 1._dl-P%omegab-P%omegac - P%omegan -P%omegak-P%omegaax - P%omegah2_rad/((P%H0/1.d2)**2.0d0)
!!!write(*, *) 'Budget for DE:', P%axfrac*(1._dl-P%omegab-P%omegac - P%omegan -P%omegak - P%omegah2_rad/((P%H0/1.d2)**2.0d0))
!!!write(*, *) 'P%omegaax', P%omegaax

     P%Hinf = Ini_Read_Double('Hinf') ! H inflation in GeV 
     P%Hinf = (10**P%Hinf)/mplanck ! computing the ratio of Hinflation to Mplanck
     !       print*, 'This is Hinflation renee', P%Hinf
     P%axion_isocurvature = Ini_Read_Logical('axion_isocurvature', .true.)
     !RL 121924 disable isocurvature
     if (P%axion_isocurvature .eqv. .true.) then
        write(*, *) 'WARNING: axion isocurvature disabled in this release, proceeding without'
        P%axion_isocurvature = .false.
     end if


  else

     P%omegab = Ini_Read_Double('omega_baryon')
     P%omegac = Ini_Read_Double('omega_cdm')
     P%omegav = Ini_Read_Double('omega_lambda')
     P%omegaax = Ini_Read_Double('omega_axion')/(P%H0/100)**2
     P%ma     = Ini_Read_Double('m_ax')  

  end if

  !	print*, 'hi renee omk', P%omegak, 'omegav', P%omegav, 'grhog', grhog, 'P%grhor', (P%grhor*(c**2.0d0)/((1.d5**2.0d0)))/3.0d0, 'omegah2_rad', P%omegah2_rad 






  !JD 08/13 begin changes for nonlinear lensing of CMB + LSS compatibility
  !P%Transfer%redshifts -> P%Transfer%PK_redshifts and P%Transfer%num_redshifts -> P%Transfer%PK_num_redshifts
  !in the P%WantTransfer loop.
  if (((P%NonLinear==NonLinear_lens .or. P%NonLinear==NonLinear_both) .and. P%DoLensing) &
       .or. P%PK_WantTransfer) then
     P%Transfer%high_precision=  Ini_Read_Logical('transfer_high_precision',.false.)
  else
     P%transfer%high_precision = .false.
  endif

  if (P%NonLinear/= NonLinear_none) call NonLinear_ReadParams(DefIni) ! RH axions

  ! RH making changes here
  if (P%PK_WantTransfer)  then
     P%WantTransfer  = .true.
     P%transfer%kmax          =  Ini_Read_Double('transfer_kmax')
     P%transfer%k_per_logint  =  Ini_Read_Int('transfer_k_per_logint')
     P%transfer%PK_num_redshifts =  Ini_Read_Int('transfer_num_redshifts')

     transfer_interp_matterpower = Ini_Read_Logical('transfer_interp_matterpower ', transfer_interp_matterpower)
     transfer_power_var = Ini_read_int('transfer_power_var',transfer_power_var)
     if (P%transfer%PK_num_redshifts > max_transfer_redshifts) stop 'Too many redshifts'
     do i=1, P%transfer%PK_num_redshifts
        P%transfer%PK_redshifts(i)  = Ini_Read_Double_Array('transfer_redshift',i,0._dl)
        transferFileNames(i)     = Ini_Read_String_Array('transfer_filename',i)
        MatterPowerFilenames(i)  = Ini_Read_String_Array('transfer_matterpower',i)
        if (TransferFileNames(i) == '') then
           TransferFileNames(i) =  trim(numcat('transfer_',i))//'.dat'
        end if
        if (MatterPowerFilenames(i) == '') then
           MatterPowerFilenames(i) =  trim(numcat('matterpower_',i))//'.dat'
        end if
        if (TransferFileNames(i)/= '') &
             TransferFileNames(i) = trim(outroot)//TransferFileNames(i)
        if (MatterPowerFilenames(i) /= '') &
             MatterPowerFilenames(i)=trim(outroot)//MatterPowerFilenames(i)
     end do
  else
     P%Transfer%PK_num_redshifts = 1
     P%Transfer%PK_redshifts = 0
  end if



  if ((P%NonLinear==NonLinear_lens .or. P%NonLinear==NonLinear_both) .and. P%DoLensing) then
     P%WantTransfer  = .true.
     call Transfer_SetForNonlinearLensing(P%Transfer)
  end if

  call Transfer_SortAndIndexRedshifts(P%Transfer)
  !JD 08/13 end changes

  P%transfer%kmax=P%transfer%kmax*(P%h0/100._dl)

  Ini_fail_on_not_found = .false.

  DebugParam = Ini_Read_Double('DebugParam',DebugParam)
  ALens = Ini_Read_Double('Alens',Alens)

  call Reionization_ReadParams(P%Reion, DefIni)
  call InitialPower_ReadParams(P%InitPower, DefIni, P%WantTensors)
  call Recombination_ReadParams(P%Recomb, DefIni)
  if (Ini_HasKey('recombination')) then
     i = Ini_Read_Int('recombination',1)
     if (i/=1) stop 'recombination option deprecated'
  end if

  call Bispectrum_ReadParams(BispectrumParams, DefIni, outroot)

  if (P%WantScalars .or. P%WantTransfer) then
     P%Scalar_initial_condition = Ini_Read_Int('initial_condition',initial_adiabatic)

     if (P%Scalar_initial_condition == initial_vector) then
        P%InitialConditionVector=0
        numstr = Ini_Read_String('initial_vector',.true.)
        read (numstr,*) P%InitialConditionVector(1:initial_iso_axion)
     end if
     if (P%Scalar_initial_condition/= initial_adiabatic) use_spline_template = .false.
  end if

  if (P%WantScalars) then
     ScalarFileName = trim(outroot)//Ini_Read_String('scalar_output_file')
     LensedFileName =  trim(outroot) //Ini_Read_String('lensed_output_file')
     LensPotentialFileName =  Ini_Read_String('lens_potential_output_file')
     if (LensPotentialFileName/='') LensPotentialFileName = concat(outroot,LensPotentialFileName)
     ScalarCovFileName =  Ini_Read_String_Default('scalar_covariance_output_file','scalCovCls.dat',.false.)
     if (ScalarCovFileName/='') then
        has_cl_2D_array = .true.
        ScalarCovFileName = concat(outroot,ScalarCovFileName)
     end if
  end if
  if (P%WantTensors) then
     TensorFileName =  trim(outroot) //Ini_Read_String('tensor_output_file')
     if (P%WantScalars)  then
        TotalFileName =  trim(outroot) //Ini_Read_String('total_output_file')
        LensedTotFileName = Ini_Read_String('lensed_total_output_file')
        if (LensedTotFileName/='') LensedTotFileName= trim(outroot) //trim(LensedTotFileName)
     end if
  end if
  if (P%WantVectors) then
     VectorFileName =  trim(outroot) //Ini_Read_String('vector_output_file')
  end if

#ifdef WRITE_FITS
  if (P%WantCls) then
     FITSfilename =  trim(outroot) //Ini_Read_String('FITS_filename',.true.)
     if (FITSfilename /='') then
        inquire(file=FITSfilename, exist=bad)
        if (bad) then
           open(unit=18,file=FITSfilename,status='old')
           close(18,status='delete')
        end if
     end if
  end if
#endif


  Ini_fail_on_not_found = .false.

  !optional parameters controlling the computation

  P%AccuratePolarization = Ini_Read_Logical('accurate_polarization',.true.)
  P%AccurateReionization = Ini_Read_Logical('accurate_reionization',.false.)
  P%AccurateBB = Ini_Read_Logical('accurate_BB',.false.)
  P%DerivedParameters = Ini_Read_Logical('derived_parameters',.true.)

  version_check = Ini_Read_String('version_check')
  if (version_check == '') then
     !tag the output used parameters .ini file with the version of CAMB being used now
     call TNameValueList_Add(DefIni%ReadValues, 'version_check', version)
  else if (version_check /= version) then
     write(*,*) 'WARNING: version_check does not match this CAMB version'
  end if
  !Mess here to fix typo with backwards compatibility !RL laughing hard about this 07/11/2023
  if (Ini_HasKey('do_late_rad_trunction')) then
     DoLateRadTruncation = Ini_Read_Logical('do_late_rad_trunction',.true.)
     if (Ini_HasKey('do_late_rad_truncation')) stop 'check do_late_rad_xxxx'
  else
     DoLateRadTruncation = Ini_Read_Logical('do_late_rad_truncation',.true.)
  end if
  DoTensorNeutrinos = Ini_Read_Logical('do_tensor_neutrinos',DoTensorNeutrinos )
  FeedbackLevel = Ini_Read_Int('feedback_level',FeedbackLevel)

  P%MassiveNuMethod  = Ini_Read_Int('massive_nu_approx',Nu_best)

  ThreadNum      = Ini_Read_Int('number_of_threads',ThreadNum)
  use_spline_template = Ini_Read_Logical('use_spline_template',use_spline_template)

  DoTensorNeutrinos = DoTensorNeutrinos .or. HighAccuracyDefault
  if (do_bispectrum) then
     lSampleBoost   = 50
  else
     lSampleBoost   = Ini_Read_Double('l_sample_boost',lSampleBoost)
  end if
  if (outroot /= '') then
     if (InputFile /= trim(outroot) //'params.ini') then
        call Ini_SaveReadValues(trim(outroot) //'params.ini',1)
     else
        write(*,*) 'Output _params.ini not created as would overwrite input'
     end if
  end if

!!!!! RL testing 061924 - this has to be before ini close
!!!  write(dfac_str,'(F5.1)') P%dfac
!!!  write(omk_str1,'(F4.2)') P%omegak
!!!  write(omk_str, '(A,A,A)') omk_str1(1:1), 'd', omk_str1(3:4)
!!!  write(max_str,'(ES32.2E3)') P%ma
!!!  max_str = adjustl(max_str)
!!!  if (P%use_axfrac) then
!!!     useaxfrac_str = 'T'
!!!     write(*, *) 'axfracT, axfrac:', P%axfrac
!!!     if (P%axfrac .eq. 1.e-6_dl) then
!!!        abun_str = '1d0e-6'
!!!     else
!!!        write(abun_str1,'(F6.4)') P%axfrac
  !write(*, *) 'abun_str T', abun_str
!!!        write(abun_str, '(A,A,A)') abun_str1(1:1), 'd', abun_str1(3:6)
!!!        write(*, *) 'abun_str:', abun_str
!!!     end if
!!!  else
!!!     useaxfrac_str = 'F'
!!!     !write(*, *) 'axfrac?', P%axfrac
!!!     if (Ini_Read_Double('omaxh2') .eq. 1.e-6_dl) then
!!!        abun_str = '1d0e-6'
!!!    else
!!!        write(abun_str1,'(F6.4)') Ini_Read_Double('omaxh2')
!!!        write(*, *) 'abun_str1:', abun_str1
  ! write(*, *) 'abun_str F, abun_str(3:3):', abun_str, abun_str(3:3)
!!!        write(abun_str, '(A,A,A)') abun_str1(1:1), 'd', abun_str1(3:6)
!!!     end if
!!!  end if

!!!  abun_str = trim(adjustl(abun_str))
  !write(*, *) 'RL checkpoint, useaxfrac_str, teststr: ', useaxfrac_str, 'teststr'
!!!  if (abun_str == '0d00') abun_str = '1e-6'
!!!  if (DoLateRadTruncation) then
!!!     dolateradtrunc_str = 'T'
!!!  else
!!!     dolateradtrunc_str = 'F'
!!!  end if
!!!  write(*, *) 'RL checkpoint, dolateradtrunc_str:', dolateradtrunc_str

  call Ini_Close

  ! DM: The place axion evolution is called
  ! This computes axion parameters and also creates the lookup table for axions during slow-roll.
  ! Tables for density, equation of state and sound-speed. grhoax_table, wax_table, cs2_table.
  ! Sampled at loga_table values. Later splined to any a necessary.
  ! Sound speed is the sound speed only before oscillations (more precisely the adiabatic sound speed)
  !, afterwards fluid representation changes.
  ! For more details see Hlozek et al 2014. arXiv:1410.2896

  !call cpu_time(clock_start) ! RH timing
  ! Run axion background evolution and then with arrays in hand for interpolation, run the regular CAMB
  !!allocate(P%grhoax_table(ntable))
  !!allocate(P%grhoax_table_buff(ntable))
  call init_massive_nu(P%omegan /=0) !RL added 07/10/23
  P%a_skip = 1._dl/(800._dl + 1._dl) !RL for skipping
  P%a_skipst = 1._dl/(1300._dl + 1._dl) !lower threshold of recombination skip
  P%dfac_skip = 0._dl !RL initializing P%dfac_skip just to make sure it doesn't get assigned to some random numbers
  call w_evolve(P, badflag)

  !First check if we're in the window before recombination where the photon ETA is an issue
  if (P%dfac .lt. 23._dl .and. P%m_ovH0 .ge. 10._dl .and. P%ma .lt. 1.e-25_dl &
       &.and. P%a_osc*(P%omegaax/(P%omegac+P%omegaax))/aeq_LCDM .gt. 0.03_dl &
       &.and. P%a_osc .lt. P%a_skipst) then !RL 022624 - is an empirically tuned number 0.03_dl 
     !if (P%a_osc .lt. 1._dl) then
     twobeta_tgt = 7.08_dl*const_pi!22.25_dl!10.5_dl*const_pi!
     !P%dfac = twobeta_tgt + 0.75_dl*const_pi !First guess using radiation domination
     P%dfac = twobeta_tgt + 0.75_dl*const_pi - twobeta_tgt**2/&
          &(4._dl*(twobeta_tgt + 2._dl*P%m_ovH0*(aeq_LCDM**1.5_dl)/&
          &sqrt(2._dl*(P%omegac + P%omegab + P%omegan + P%omegaax)))) !First guess considering matter-radiation equality in LCDM
     ntable = nint(P%dfac*100) + 1     
     call w_evolve(P, badflag)
     call get_phase_info(P, y_phase, beta_coeff, movHETA_new, twobeta_new)

     !Rough guess of target ETA values using beta. This is useful in bracket finding
     movHETA_beta = (twobeta_tgt + const_pi*3._dl*(1.0_dl + y_phase)/(4.0_dl + 3.0_dl*y_phase))/beta_coeff
     hETA_beta = (P%dfac/movHETA_beta) * (P%ah_osc/P%a_osc)
     beta_tol = 2.e-2_dl*const_pi
     if (abs(twobeta_new - twobeta_tgt) .gt. beta_tol) then !RL: if so, we initiate the shooting and bisecting process

        !First we find the bracket
        hosc_old1 = P%ah_osc/P%a_osc
        hETA_old1 = P%ahosc_ETA/P%a_osc
        twobeta_old1 = twobeta_new
        hosc_new = 2._dl*(hETA_beta-hETA_old1)+hosc_old1
        P%dfac = P%dfac * ((P%ah_osc/P%a_osc)/hosc_new)
        ntable = nint(P%dfac*100) + 1
        call w_evolve(P, badflag)
        call get_phase_info(P, y_phase, beta_coeff, movHETA_new, twobeta_new)
        hosc_old2 = P%ah_osc/P%a_osc
        hETA_old2 = P%ahosc_ETA/P%a_osc
        twobeta_old2 = twobeta_new

        iter_dfacETA = 1
        do while (iter_dfacETA .lt. 500)
           if (abs(twobeta_new - twobeta_tgt) .lt. beta_tol) then 

              iter_dfacETA = -1
              exit
           else if ((twobeta_old2-twobeta_tgt)*(twobeta_old1-twobeta_tgt) .lt. 0._dl) then

              exit
           else
              hosc_new = 1._dl*(hETA_beta-hETA_old1)+hosc_old2 !Still use hETA_beta-hETA_old1 to make the shooting step larger, or else the shooting step diminishes each time
              P%dfac = P%dfac * ((P%ah_osc/P%a_osc)/hosc_new)

              ntable = nint(P%dfac*100) + 1
              call w_evolve(P, badflag)
              call get_phase_info(P, y_phase, beta_coeff, movHETA_new, twobeta_new)
              hosc_old2 = P%ah_osc/P%a_osc
              hETA_old2 = P%ahosc_ETA/P%a_osc
              twobeta_old2 = twobeta_new
           end if
        end do

        !If one of the brackets hit the target (not very likely), no need to proceed with the bisection. Else do bisection
        if (iter_dfacETA .ne. -1) then
           !After bracket is found, do bisection - reuse iter_dfacETA
           ! Calculate new guess as the midpoint of the bracket
           hosc_new = (hosc_old1 + hosc_old2) / 2._dl
           P%dfac = P%dfac * ((P%ah_osc/P%a_osc)/hosc_new)
           ntable = nint(P%dfac*100) + 1
           call w_evolve(P, badflag)
           call get_phase_info(P, y_phase, beta_coeff, movHETA_new, twobeta_new)
           ! Check which side of the target the new guess falls on and update the brackets

           iter_dfacETA = 1
           do while (iter_dfacETA .lt.500)
              if (abs(twobeta_new - twobeta_tgt) .lt. beta_tol) then 

                 exit
              else
                 if ((twobeta_new - twobeta_tgt) * (twobeta_old1 - twobeta_tgt) .lt. 0._dl) then
                    hosc_old2 = P%ah_osc/P%a_osc
                    hETA_old2 = P%ahosc_ETA/P%a_osc
                    twobeta_old2 = twobeta_new
                 else
                    hosc_old1 = P%ah_osc/P%a_osc
                    hETA_old1 = P%ahosc_ETA/P%a_osc
                    twobeta_old1 = twobeta_new
                 end if
                 ! Calculate new guess as the midpoint of the bracket
                 hosc_new = (hosc_old1 + hosc_old2) / 2._dl
                 P%dfac = P%dfac * ((P%ah_osc/P%a_osc)/hosc_new)

                 ntable = nint(P%dfac*100) + 1
                 call w_evolve(P, badflag)
                 call get_phase_info(P, y_phase, beta_coeff, movHETA_new, twobeta_new)

              end if
           end do
        end if
     else
     end if

  end if

  do iter_dfac = 1, 500
!!!write(*, *) 'P%a_osc, P%a_skip*(1._dl - 1.e-2_dl), P%a_skipst', P%a_osc, P%a_skip*(1._dl - 1.e-2_dl), P%a_skipst
     if (P%a_osc .lt. P%a_skip*(1._dl - 1.e-2_dl) .and. P%a_osc .ge. P%a_skipst) then !
        !RL 032024: 1e-2 is the tolerence to eliminate additional loops if we don't skip to exactly after a_skip due to numerical factors
        P%dfac = P%dfac_skip
        ntable = nint(P%dfac*100) + 1
        call w_evolve(P, badflag)
!!!write(*, *) 'Rayne, each call, P%a_osc, P%a_skip, P%a_skip*(1._dl - 1.e-2_dl), P%dfac, P%dfac_skip', P%a_osc, P%a_skip, P%a_skip*(1._dl - 1.e-2_dl), P%dfac, P%dfac_skip
     else
        exit
     end if
  end do

  if (iter_dfac .gt. 500 .and. P%a_osc .lt. P%a_skip) then
     print*, 'Warning: maximum iteration reached, but aosc still not skipped sufficiently: &
          &P%a_osc, P%a_skip', P%a_osc, P%a_skip
  end if


  !!write(*, *) CP%dfac, dfac_str
  !!write(*, *) 'omk_str', omk_str ', '_omk', omk_str'
  !!write(testthetafilename, '(A,A,A,A,A,A,A,A,A,A,A,A,I0,A,I0,A,A,A,I0,A,A,A)') '../Testdata/housecleaning3/AxiEF2_axf', useaxfrac_str, abun_str, '_omdah20d12_H0', H0_str, '_hyb2ktau2_Ct6tau*t12hw6_wEFAtowEF_m', max_str(1:1), 'd', max_str(3:4), 'e-', max_str(8:9), '_m-nuon_aB', int(AccuracyBoost), '_laB', int(lAccuracyBoost), '_Lratc', dolateradtrunc_str, '_dfac', int(P%dfac), 'd', dfac_str(5:5), '_zstar_rs-star_thetastar_sigma8_om0.dat' !_scalclsF_dfac
  !!write(testthetafilename, '(A,A,A,A,A,A,A,A,A,A,A,A,I0,A,I0,A,A,A,I0,A,A,A)') '../testdata_midway/Asm_axf', useaxfrac_str, abun_str, '_omdah20d12_H0', H0_str, '_hyb2ktau2_Ct6tau*t12hw6_wEFAtowEF_m', max_str(1:1), 'd', max_str(3:4), 'e-', max_str(8:9), '_m-nuon_aB', int(AccuracyBoost), '_laB', int(lAccuracyBoost), '_Lratc', dolateradtrunc_str, '_dfac', int(P%dfac), 'd', dfac_str(5:5), '_zstar_rs-star_thetastar_sigma8_om0.dat' !

  !!  call CreateTxtFile(filenametest,050924)
  !call cpu_time(clock_stop) ! RH timing 
  !print*, 'timing after dans routine', clock_stop - clock_start

  if (.not. CAMB_ValidateParams(P)) stop 'Stopped due to parameter error'

#ifdef RUNIDLE
  call SetIdle
#endif


  !call cpu_time(clock_start) ! RH timing
  !write(*, *) 'RL, clock_start'


!!!!! This is where we need to be renee, but where are the cls
!!!! regenerate the spectra here
  if (global_error_flag==0) then 

!!!write(*, *) 'RL, calling CAMB_GetResults'
     call CAMB_GetResults(P)


     if (P%axion_isocurvature) then 
        !          print*, 'computing isocurvature' 
        if (P%WantScalars) then
           allocate(RHCl_temp(size(Cl_scalar, 1), size(Cl_scalar, 2), size(Cl_scalar, 3))) !RL 111523
           !!print *, 'shape(Cl_scalar)', shape(Cl_scalar) !RL 111523: shape doesn't seem to work. But let's leave isocurvature for later
           RHCl_temp(lmin:P%Max_l,1,C_Temp:C_last) = Cl_scalar(lmin:P%Max_l,1,C_Temp:C_last)
        end if

        if (P%DoLensing) then
           allocate(RHCl_temp_lensed(size(Cl_lensed, 1), size(Cl_lensed, 2), size(Cl_lensed, 3))) !RL 111523
           RHCl_temp_lensed(lmin:P%Max_l,1,C_Temp:C_Cross) = Cl_lensed(lmin:P%Max_l,1,C_Temp:C_Cross)
        end if

        if (P%WantTensors) then
           allocate(RHCl_temp_tensor(size(Cl_tensor, 1), size(Cl_tensor, 2), size(Cl_tensor, 3)))  !RL 111523
           RHCl_temp_tensor(lmin:P%Max_l,1,C_Temp:C_Cross) = Cl_tensor(lmin:P%Max_l,1,C_Temp:C_Cross)
        end if

        P%Scalar_initial_condition = 6
        P%InitPower%rat(1) =  0
        P%InitPower%ant(1) = 0
        P%InitPower%ScalarPowerAmp(1) = P%amp_i
        P%InitPower%an(1)= 1-P%r_val/8.d0
        call CAMB_GetResults(P)
        if (P%WantScalars)  then 
           Cl_scalar(lmin:P%Max_l,1,C_Temp:C_last) = Cl_scalar(lmin:P%Max_l,1,C_Temp:C_last)  &
                +  RHCl_temp(lmin:P%Max_l,1,C_Temp:C_last)
        end if

        if (P%DoLensing) then 
           Cl_lensed(lmin:lmax_lensed,1,C_Temp:C_Cross) = Cl_lensed(lmin:lmax_lensed,1,C_Temp:C_Cross) &
                +  RHCl_temp_lensed(lmin:lmax_lensed,1,C_Temp:C_Cross) 
        end if

        if (P%WantTensors) then 
           Cl_tensor(lmin:P%Max_l_tensor,1,C_Temp:C_Cross) = Cl_tensor(lmin:P%Max_l_tensor,1,C_Temp:C_Cross)  &
                +  RHCl_temp_tensor(lmin:P%Max_l_tensor,1,C_Temp:C_Cross) 
        end if
     end if
  end if


  if (global_error_flag/=0) then
     write (*,*) 'Error result '//trim(global_error_message)
     stop
  endif
  !!

  !call cpu_time(clock_stop) ! RH timing 
  !print*, 'after getresults', clock_stop - clock_start
  if (P%PK_WantTransfer) then
     call Transfer_SaveToFiles(MT,TransferFileNames)
     call Transfer_SaveMatterPower(MT,MatterPowerFileNames)
     call Transfer_output_sig8(MT)
  end if
  !!
  if (P%WantCls) then
     call output_cl_files(ScalarFileName, ScalarCovFileName, TensorFileName, TotalFileName, &
          LensedFileName, LensedTotFilename, output_factor)
     !
     call output_lens_pot_files(LensPotentialFileName, output_factor)
     !
     if (P%WantVectors) then
        call output_veccl_files(VectorFileName, output_factor)
     end if
     !
#ifdef WRITE_FITS
     if (FITSfilename /= '') call WriteFitsCls(FITSfilename, CP%Max_l)
#endif
  end if

  !RL 111123-----
  if (allocated(loga_table)) deallocate(loga_table)
  if (allocated(phinorm_table)) deallocate(phinorm_table)
  if (allocated(phidotnorm_table)) deallocate(phidotnorm_table)
  if (allocated(phinorm_table_ddlga)) deallocate(phinorm_table_ddlga)
  if (allocated(phidotnorm_table_ddlga)) deallocate(phidotnorm_table_ddlga)
  if (allocated(rhoaxh2ovrhom_logtable)) deallocate(rhoaxh2ovrhom_logtable)
  if (allocated(rhoaxh2ovrhom_logtable_buff)) deallocate(rhoaxh2ovrhom_logtable_buff)
  if (allocated(RHCl_temp)) deallocate(RHCl_temp)
  if (allocated(RHCl_temp_lensed)) deallocate(RHCl_temp_lensed)
  if (allocated(RHCl_temp_tensor)) deallocate(RHCl_temp_tensor)

  call CAMB_cleanup
  !!call cpu_time(clock_totstop) ! RH timing	
  !!print*, 'Total time', clock_totstop - clock_totstart
  stop
  !
100 stop 'Must give num_massive number of integer physical neutrinos for each eigenstate'
end program driver


#ifdef RUNIDLE
!If in Windows and want to run with low priorty so can multitask
subroutine SetIdle
  USE DFWIN
  Integer dwPriority
  Integer CheckPriority

  dwPriority = 64 ! idle priority
  CheckPriority = SetPriorityClass(GetCurrentProcess(), dwPriority)

end subroutine SetIdle
#endif

