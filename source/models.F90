!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Module to relate fitted paramertes with different models
!
! Copyright - 2026 Ada Lovelace Centre (ALC)
!             Scientific Computing Department (SCD)
!             The Science and Technology Facilities Council (STFC)
!
! Author:        i.scivetti  March 2026
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Module models 

  Use constants,     Only : exp_ocf_lower_bound,&
                            calio2021,&
                            lipari1980,&
                            morrone2006,&
                            pi, &
                            fmt_mean_std, &
                            fmt_mean_only

  Use fitting_setup, Only : fit_type
  
  Use numprec,       Only : wi,&
                            wp 
  
  Use unit_output,   Only : error_stop, &
                            info

  Implicit None
  Private

  Type quantity
    Real(Kind=wp), Allocatable :: value(:)
    Real(Kind=wp) :: suma
    Real(Kind=wp) :: mean
    Real(Kind=wp) :: std
    Character(Len=256) :: label
  End Type
  
  Type, Public ::  model_type
    Type(quantity) :: A_inf
    Type(quantity) :: kf
    Type(quantity) :: kb
    Type(quantity) :: k2
    Type(quantity) :: inv_kf
    Type(quantity) :: inv_kb
    Type(quantity) :: inv_k2
    Type(quantity) :: D
    Type(quantity) :: Dw
    Type(quantity) :: theta_cone
    Type(quantity) :: theta_lib
    Type(quantity) :: Tau_cone
    Type(quantity) :: Tau_sf
    Type(quantity) :: Tau_fast
    Type(quantity) :: Tau_slow
    Type(quantity) :: Tau_lib
    Type(quantity) :: Tau_host
    Logical        :: caution
  Contains
    Private
      Procedure  :: alloc_kinetics    => allocate_transfer_kinetics
      Procedure  :: alloc_diffusion   => allocate_diffusion_array
      Procedure  :: alloc_wic         => allocate_wic_array
      Procedure  :: alloc_react_ocf   => allocate_react_ocf_array
      Final :: cleanup    
  End Type
  
  Public :: extract_physics_from_fitting

Contains

  Subroutine allocate_react_ocf_array(T, std_samples)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Allocate arrays to compute compute_time related
    ! to the orientation of reactive species
    !
    ! author    - i.scivetti March 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Class(model_type),   Intent(InOut)  :: T
    Integer(Kind=wi),    Intent(In   )  :: std_samples
    
    Integer(Kind=wi)    :: fail(3)
    Character(Len=256)  :: message
    
    fail=0
    Write (message,'(1x,1a)') '***ERROR: Allocation problems for kinetcs&
                                & (subroutine allocate_wic_array). Check the value of "std_samples"'
                                
    Allocate(T%Tau_sf%value(std_samples),   Stat=fail(1))
    Allocate(T%Tau_fast%value(std_samples), Stat=fail(2))
    Allocate(T%Tau_slow%value(std_samples), Stat=fail(3))
    
    If (Any(fail > 0)) Then
      Call error_stop(message)
    End If

  End Subroutine allocate_react_ocf_array  
  
  Subroutine allocate_wic_array(T, std_samples)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Allocate arrays to compute compute_wobbling_in_cone
    ! related parameters
    !
    ! author    - i.scivetti March 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Class(model_type),   Intent(InOut)  :: T
    Integer(Kind=wi),    Intent(In   )  :: std_samples    
    
    Integer(Kind=wi)    :: fail(3)
    Character(Len=256)  :: message
    
    fail=0
    Write (message,'(1x,1a)') '***ERROR: Allocation problems for kinetcs&
                                & (subroutine allocate_wic_array). Check the value of "std_samples"'
                                
    Allocate(T%Dw%value(std_samples), Stat=fail(1))
    Allocate(T%theta_cone%value(std_samples), Stat=fail(2))
    Allocate(T%theta_lib%value(std_samples), Stat=fail(3))
    
    If (Any(fail > 0)) Then
      Call error_stop(message)
    End If

  End Subroutine allocate_wic_array
  
  Subroutine allocate_transfer_kinetics(T, std_samples)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Allocate arrays to compute kinetic quantities
    !
    ! author    - i.scivetti March 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Class(model_type),   Intent(InOut)  :: T
    Integer(Kind=wi),    Intent(In   )  :: std_samples    
    
    Integer(Kind=wi)    :: fail(3)
    Character(Len=256)  :: message
    
    fail=0
    Write (message,'(1x,1a)') '***ERROR: Allocation problems for kinetcs&
                                & (subroutine allocate_transfer_kinetics). Check the value of "std_samples"'
                                
    Allocate(T%kf%value(std_samples), Stat=fail(1))
    Allocate(T%kb%value(std_samples), Stat=fail(2))
    Allocate(T%k2%value(std_samples), Stat=fail(3))
    
    If (Any(fail > 0)) Then
      Call error_stop(message)
    End If

  End Subroutine allocate_transfer_kinetics  

  Subroutine allocate_diffusion_array(T, std_samples)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Allocate arrays to compute diffusion
    !
    ! author    - i.scivetti March 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Class(model_type),   Intent(InOut)  :: T
    Integer(Kind=wi),    Intent(In   )  :: std_samples    
    
    Integer(Kind=wi)    :: fail
    Character(Len=256)  :: message
    
    fail=0
    Write (message,'(1x,1a)') '***ERROR: Allocation problems for kinetcs&
                                & (subroutine allocate_diffusion_array). Check the value of "std_samples"'
                                
    Allocate(T%D%value(std_samples), Stat=fail)
    
    If (fail > 0) Then
      Call error_stop(message)
    End If

  End Subroutine allocate_diffusion_array
  
 
  Subroutine cleanup(T)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Deallocate variables
    !
    ! author    - i.scivetti March 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(model_type) :: T

    If (Allocated(T%kf%value)) Then
      Deallocate(T%kf%value)
    End If 

    If (Allocated(T%kb%value)) Then
      Deallocate(T%kb%value)
    End If 
    
    If (Allocated(T%k2%value)) Then
      Deallocate(T%k2%value)
    End If 

    If (Allocated(T%D%value)) Then
      Deallocate(T%D%value)
    End If 

    If (Allocated(T%theta_cone%value)) Then
      Deallocate(T%theta_cone%value)
    End If 
    
    If (Allocated(T%theta_lib%value)) Then
      Deallocate(T%theta_lib%value)
    End If 

    If (Allocated(T%Dw%value)) Then
      Deallocate(T%Dw%value)
    End If    
    
    If (Allocated(T%Tau_sf%value)) Then
      Deallocate(T%Tau_sf%value)
    End If 
    
    If (Allocated(T%Tau_fast%value)) Then
      Deallocate(T%Tau_fast%value)
    End If 

    If (Allocated(T%Tau_slow%value)) Then
      Deallocate(T%Tau_slow%value)
    End If     
    
  End Subroutine cleanup

  Subroutine extract_physics_from_fitting(fit_data, model_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Extract physical quantities from the fitted parameters using models
    !
    ! author    - i.scivetti March 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(fit_type),   Intent(InOut) :: fit_data
    Type(model_type), Intent(InOut) :: model_data

    ! Define the time limit for evaluation of the computed, time-related quantities
    If (fit_data%end_time%fread) Then
      fit_data%t_limit=fit_data%end_time%value
    Else
      fit_data%t_limit=fit_data%input_data(fit_data%num_input_data)%time
    End If    

    ! Initialise flag to control the computed quantities
    model_data%caution=.False.
    
    If (Trim(fit_data%what_to_fit%type) == 'ocf') Then
      If (.Not. fit_data%reactive_species%stat) Then
        Call compute_wobbling_in_cone(fit_data, model_data)
      Else
        Call orientation_reactive_species(fit_data, model_data)
      End If
    End If
    
    If (Trim(fit_data%functional_form) == 'tcf_1exp_to_zero'     .Or. &
        Trim(fit_data%functional_form) == 'tcf_1exp_to_constant' .Or. &
        Trim(fit_data%functional_form) == 'tcf_2exp_to_zero'     .Or. &
        Trim(fit_data%functional_form) == 'spcf_1exp_to_zero'     .Or. &
        Trim(fit_data%functional_form) == 'spcf_1exp_to_constant' .Or. &
        Trim(fit_data%functional_form) == 'spcf_2exp_to_zero') Then          
      Call compute_rates(fit_data, model_data)
    End If

    If (Trim(fit_data%functional_form) == 'msd_linear') Then
      Call compute_diffusion_coefficient(fit_data, model_data)
    End If

  End Subroutine extract_physics_from_fitting

  Subroutine compute_ocf_react_times(fit_data, i, model_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to compute times from reactive ocf
    !
    ! author    - i.scivetti March 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(fit_type),   Intent(InOut) :: fit_data
    Integer(Kind=wi), Intent(In   ) :: i
    Type(model_type), Intent(InOut) :: model_data
    
    Real(Kind=wp)  :: Tau_sf, Tau_fast, Tau_slow
    Real(Kind=wp)  :: T1, T2, T3

    ! Initialise from fitting
    T3 =fit_data%param%T3%value(i)
    T2 =fit_data%param%T2%value(i)
    T1 =fit_data%param%T1%value(i)

    If (Trim(fit_data%functional_form) == 'ocf_2exp_to_constant' .Or. &
      & Trim(fit_data%functional_form) == 'ocf_2exp_to_zero') Then
      Tau_fast=1.0_wp/(1.0_wp/T2+1.0_wp/T3)
      Tau_slow=T3
      If (fit_data%fit_superfast_time%stat) Then
        Tau_sf=1.0_wp/(1.0_wp/T1+1.0_wp/T2+1.0_wp/T3)
      End If              
    Else If (Trim(fit_data%functional_form) == 'ocf_1exp_to_constant' .Or. &
             Trim(fit_data%functional_form) == 'ocf_1exp_to_zero') Then
      Tau_fast=T2
      Tau_slow=0.0_wp         
      If (fit_data%fit_superfast_time%stat) Then
        Tau_sf=1.0_wp/(1.0_wp/T1+1.0_wp/T2)
      End If      
    End If

    If (.Not. fit_data%fit_superfast_time%stat) Then
      Tau_sf=0.0_wp
    End If    
    
    model_data%Tau_sf%value(i)=Tau_sf
    model_data%Tau_fast%value(i)=Tau_fast
    model_data%Tau_slow%value(i)=Tau_slow
  
  End Subroutine compute_ocf_react_times    
  
  Subroutine orientation_reactive_species(fit_data, model_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to analyse the fitted parameters from the orientation of
    ! reactive species
    !
    ! author    - i.scivetti March 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(fit_type),   Intent(InOut)   :: fit_data
    Type(model_type), Intent(InOut)   :: model_data
    
    Character(Len=256) :: messages(5)
    Integer            :: i
 
    Call model_data%alloc_react_ocf(fit_data%std_samples%value)

    model_data%Tau_slow%label='Tau_slow'
    model_data%Tau_fast%label='Tau_fast'    
    model_data%Tau_sf%label='Tau*' 
    model_data%A_inf%label='A_inf'
    
    !reset sums
    model_data%Tau_sf%suma= 0.0_wp
    model_data%Tau_fast%suma= 0.0_wp
    model_data%Tau_slow%suma= 0.0_wp
    
    Do i = 1, fit_data%std_samples%value
      If (fit_data%param%valid(i)) Then
        fit_data%xstd=(i-1)*fit_data%delta_std-fit_data%max_times_std%value
        Call compute_ocf_react_times(fit_data, i, model_data)
        model_data%Tau_sf%suma=model_data%Tau_sf%suma+model_data%Tau_sf%value(i)*exp(-0.5*fit_data%xstd**2)
        model_data%Tau_fast%suma =model_data%Tau_fast%suma +model_data%Tau_fast%value(i)*exp(-0.5*fit_data%xstd**2)
        model_data%Tau_slow%suma=model_data%Tau_slow%suma+model_data%Tau_slow%value(i)*exp(-0.5*fit_data%xstd**2)
      End If
    End Do
 
    ! Define means
    model_data%Tau_sf%mean=model_data%Tau_sf%suma/fit_data%sum_weights    
    model_data%Tau_fast%mean =model_data%Tau_fast%suma/fit_data%sum_weights     
    model_data%Tau_slow%mean=model_data%Tau_slow%suma/fit_data%sum_weights

    !reset sums
    model_data%Tau_sf%suma=0.0_wp
    model_data%Tau_fast%suma=0.0_wp
    model_data%Tau_slow%suma=0.0_wp
    
    Do i = 1, fit_data%std_samples%value
      If (fit_data%param%valid(i)) Then
        fit_data%xstd=(i-1)*fit_data%delta_std-fit_data%max_times_std%value
        model_data%Tau_sf%suma=model_data%Tau_sf%suma+&
                              &(model_data%Tau_sf%value(i)-model_data%Tau_sf%mean)**2 * exp(-0.5*fit_data%xstd**2)
        model_data%Tau_fast%suma=model_data%Tau_fast%suma+&
                              &(model_data%Tau_fast%value(i)-model_data%Tau_fast%mean)**2 * exp(-0.5*fit_data%xstd**2)
        model_data%Tau_slow%suma=model_data%Tau_slow%suma+&
                              &(model_data%Tau_slow%value(i)-model_data%Tau_slow%mean)**2 * exp(-0.5*fit_data%xstd**2)
      End If
    End Do
    
    ! Define std
    model_data%Tau_sf%std=sqrt(model_data%Tau_sf%suma/fit_data%sum_weights)    
    model_data%Tau_fast%std=sqrt(model_data%Tau_fast%suma/fit_data%sum_weights)     
    model_data%Tau_slow%std=sqrt(model_data%Tau_slow%suma/fit_data%sum_weights)
    
    Write(messages(1), '(1x,a)') 'Description of relevant&
                                 & features related to the "orientation of the reactive species" ['&
                                 &//Trim(calio2021)//']:'
 
    If (Trim(fit_data%functional_form) == 'ocf_2exp_to_constant') Then 
       Write(messages(2), '(3x,a)') '* Remanent orientational fraction (A_inf); dimensionless.'
       Write(messages(3), '(3x,a)') '* Slow relaxation time (Tau_slow); in ps.'
       Write(messages(4), '(3x,a)') '* Fast relaxation time (Tau_fast); in ps.'
       Call info(messages, 4)
    Else If (Trim(fit_data%functional_form) == 'ocf_2exp_to_zero') Then  
       Write(messages(2), '(3x,a)') '* Slow relaxation time (Tau_slow); in ps.'
       Write(messages(3), '(3x,a)') '* Fast relaxation time (Tau_fast); in ps.'
       Call info(messages, 3)
    Else If (Trim(fit_data%functional_form) == 'ocf_1exp_to_constant') Then  
       Write(messages(2), '(3x,a)') '* Remanent orientational fraction (A_inf); dimensionless.'
       Write(messages(3), '(3x,a)') '* Net relaxation time (Tau_net); in ps.'
       Call info(messages, 3)
    Else If (Trim(fit_data%functional_form) == 'ocf_1exp_to_zero') Then  
       Write(messages(2), '(3x,a)') '* Net relaxation time (Tau_net); in ps.'
       Call info(messages, 2)
    End If

    If (fit_data%fit_superfast_time%stat) Then
      Write(messages(1), '(3x,a)') '* Relaxation time from "superfast" changes (Tau*); in ps.'
    Call info(messages, 1)
    End If
    
    Write(messages(1), '(1x,a)') 'The user is referred to the attached notes for details.'
    Write(messages(2), '(a)') ' ----------------------------------------'    
    If (fit_data%std_samples%value /=1) Then
      Write(messages(3), '(1x,a,12x,a,10x,a)') 'Quantity', 'MEAN', 'STD'
    Else
      Write(messages(3), '(1x,a,9x,a)') 'Quantity', 'MEAN (STD not computed)'
    End If    
      Write(messages(4), '(a)') ' ----------------------------------------'    
    Call info(messages, 4)       

    If (Trim(fit_data%functional_form) == 'ocf_2exp_to_constant') Then
      model_data%A_inf%mean=fit_data%param%A3%mean
      model_data%A_inf%std =fit_data%param%A3%std
      Call print_quantity(messages(1), fit_data%std_samples%value, 'A_inf   ', model_data%A_inf)
      Call print_quantity(messages(2), fit_data%std_samples%value, 'Tau_slow', model_data%Tau_slow)
      Call print_quantity(messages(3), fit_data%std_samples%value, 'Tau_fast', model_data%Tau_fast)
      Call info(messages, 3)
    Else If (Trim(fit_data%functional_form) == 'ocf_2exp_to_zero') Then
      Call print_quantity(messages(1), fit_data%std_samples%value, 'Tau_slow', model_data%Tau_slow)
      Call print_quantity(messages(2), fit_data%std_samples%value, 'Tau_fast', model_data%Tau_fast)    
      Call info(messages, 2)
    Else If (Trim(fit_data%functional_form) == 'ocf_1exp_to_constant') Then
      model_data%A_inf%mean=fit_data%param%A2%mean
      model_data%A_inf%std =fit_data%param%A2%std
      model_data%Tau_fast%label='Tau_net'
      Call print_quantity(messages(1), fit_data%std_samples%value, 'A_inf   ', model_data%A_inf)
      Call print_quantity(messages(2), fit_data%std_samples%value, 'Tau_net ', model_data%Tau_fast)    
      Call info(messages, 2)
    Else If (Trim(fit_data%functional_form) == 'ocf_1exp_to_zero') Then
      model_data%Tau_fast%label='Tau_net'
      Call print_quantity(messages(1), fit_data%std_samples%value, 'Tau_net ', model_data%Tau_fast)
      Call info(messages, 1)
    End If
    
    If (fit_data%fit_superfast_time%stat) Then
     Call print_quantity(messages(1), fit_data%std_samples%value, 'Tau*    ', model_data%Tau_sf)
     Call info(messages, 1)
    End If
   
    Call info(' ----------------------------------------', 1)
    
    Call compare_wrt_time_domain(model_data%Tau_fast, fit_data, model_data%caution)
   
    If (Trim(fit_data%functional_form) == 'ocf_2exp_to_constant' .Or. &
        Trim(fit_data%functional_form) == 'ocf_2exp_to_zero') Then    
        
      Call compare_mean_values(model_data%Tau_fast, model_data%Tau_slow, model_data%caution)
      Call compare_wrt_time_domain(model_data%Tau_slow, fit_data, model_data%caution)
      If (fit_data%std_samples%value > 1) Then
        Call compare_bounds(model_data%Tau_fast, model_data%Tau_slow, model_data%caution)
      End If
     
    End If

    If (fit_data%fit_superfast_time%stat) Then
      Call compare_mean_values(model_data%Tau_sf, model_data%Tau_fast, model_data%caution)
      If (fit_data%std_samples%value > 1) Then
        Call compare_bounds(model_data%Tau_sf, model_data%Tau_fast, model_data%caution)
      End If
    End If

    ! Evaluate remanent contribution
    If (Trim(fit_data%functional_form) == 'ocf_2exp_to_constant' .Or. &
        Trim(fit_data%functional_form) == 'ocf_1exp_to_constant') Then

       If (model_data%A_inf%mean<0) Then
         model_data%caution=.True.
         Write(messages(1), '(1x,a)') 'CAUTION: the mean value of "A_inf" is negative. THIS SHOULD NOT HAPPEN!'
         Call info(messages, 1)                
       End If
        
      If (fit_data%std_samples%value > 1) Then
        Call mean_vs_std(model_data%A_inf, model_data%caution)

        If ((model_data%A_inf%mean - model_data%A_inf%std) < 0.0_wp) Then
          model_data%caution=.True.
          Write(messages(1), '(1x,a)') 'CAUTION: the lower bound of "A_inf" is negative.'
          Call info(messages, 1)                
        End If        
      End If

    End If        
    
    If (model_data%caution) Then
      Call advise_user(fit_data)
    End If
    
  End Subroutine orientation_reactive_species
  
  Subroutine compute_wic_quantities(fit_data, i, model_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to compute wic quantities
    !
    ! author    - i.scivetti March 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(fit_type),   Intent(InOut) :: fit_data
    Integer(Kind=wi), Intent(In   ) :: i
    Type(model_type), Intent(InOut) :: model_data
    
    Character(Len=256)  :: messages(5)
    Real(Kind=wp)  :: A1, A2, T1, T2
    Real(Kind=wp)  :: cos_cone, cos_lib, cos0 
    Real(Kind=wp)  :: t_net, term1, term2 

    ! Initialise from fitting
    T2 =fit_data%param%T2%value(i)
    A1 =fit_data%param%A1%value(i)

    If (Trim(fit_data%functional_form) /= 'ocf_1exp_zero') Then  
      A2 =fit_data%param%A2%value(i)
    Else
      A2=0.0_wp
    End If

    If (fit_data%fit_superfast_time%stat) Then
      T1=fit_data%param%T1%value(i)
    Else
      T1=0.0_wp
    End If
    
    ! Values
    If (A1>=0) Then
      cos_lib=(-1.0_wp+sqrt(1+8.0_wp*sqrt(A1)))/2.0_wp
      model_data%theta_lib%value(i)=acos(cos_lib)/pi*180
    End If

    If (Trim(fit_data%functional_form) /= 'ocf_1exp_to_zero') Then  
      If (A2 >= 0.0_wp) Then
        cos_cone=(-1.0_wp+sqrt(1+8.0_wp*sqrt(A2)))/2.0_wp
      Else
        Write(messages(1), '(a)') '*** PROBLEMS with the wobbling-in-a-cone model'
        Write(messages(3), '(3x,a)') 'Review the input data for fitting.'
        Write(messages(5), '(3x,a)') 'If that does not work, try changing the option of "fitting_function".'
        If (fit_data%std_samples%value /= 1) Then
          If (i==(fit_data%std_samples%value-1)/2+1) Then
            Write(messages(2), '(a)') '*** Fitted parameter A2 for the mean set of values (STD=0) is negative!'
            Write(messages(4), '(3x,a)') 'Consider changing the time range for analysis.'
          Else
            Write(messages(2), '(a)') '*** At least one of the fitted values for A2 is negative!'
            Write(messages(4), '(3x,a)') 'Consider changing the time range for analysis and/or&
                                    & reducing the value of "max_times_std".'
          End If
        Else
          Write(messages(2), '(a)') '*** Fitted parameter A2 for the mean set of values (STD=0) is negative!'
          Write(messages(4), '(3x,a)') 'Consider changing the time range for analysis.'
        End If
        Call info(messages, 5)
        
        Call error_stop(' ')
        
      End If
    Else
      cos_cone=0.0_wp
    End If  
  
    cos0=cos_cone
    model_data%theta_cone%value(i)=acos(cos_cone)/pi*180    
    
    t_net=(1-A1)*(T1*T2)/(T1+T2)+(A1-A2)*T2
    term1=(cos0**2 *(1.0_wp+cos0)**2 * (log((1.0_wp+cos0)/2)+(1-cos0)/2))/(cos0-1)/2.0_wp
    term2=(1-cos0)*(6.0_wp+8.0_wp*cos0-cos0**2-12.0_wp*cos0**3-7.0_wp*cos0**4)/24.0_wp
    
    model_data%Dw%value(i)=(term1+term2)/t_net
  
  End Subroutine compute_wic_quantities  
  
  Subroutine compute_wobbling_in_cone(fit_data, model_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to compute quantities from the wobbling in a cone model
    !
    ! author    - i.scivetti March 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(fit_type),   Intent(InOut) :: fit_data
    Type(model_type), Intent(InOut) :: model_data

    Character(Len=256)  :: messages(5)
    Integer        :: i

    Call model_data%alloc_wic(fit_data%std_samples%value)

    !Define labels for quantities
    model_data%theta_lib%label ='theta_lib'
    model_data%theta_cone%label='theta_cone'  
    model_data%Tau_cone%label='Tau_cone'
    model_data%Tau_host%label='Tau_host'
    model_data%Tau_lib%label ='Tau_lib'
    model_data%Dw%label='Dw'

    ! Initialise relevant arrays
    model_data%theta_cone%suma=0.0_wp
    model_data%theta_lib%suma=0.0_wp
    model_data%Dw%suma=0.0_wp
    
    Do i = 1, fit_data%std_samples%value
      If (fit_data%param%valid(i)) Then
        fit_data%xstd=(i-1)*fit_data%delta_std-fit_data%max_times_std%value
        Call compute_wic_quantities(fit_data, i, model_data)
        model_data%theta_cone%suma=model_data%theta_cone%suma+model_data%theta_cone%value(i)*exp(-0.5*fit_data%xstd**2)
        model_data%theta_lib%suma =model_data%theta_lib%suma +model_data%theta_lib%value(i)*exp(-0.5*fit_data%xstd**2)
        model_data%Dw%suma=model_data%Dw%suma+model_data%Dw%value(i)*exp(-0.5*fit_data%xstd**2)
      End If
    End Do
 
    model_data%theta_cone%mean=model_data%theta_cone%suma/fit_data%sum_weights    
    model_data%theta_lib%mean =model_data%theta_lib%suma/fit_data%sum_weights     
    model_data%Dw%mean=model_data%Dw%suma/fit_data%sum_weights
    model_data%Tau_cone%mean=fit_data%param%T2%mean

    model_data%theta_cone%suma=0.0_wp
    model_data%theta_lib%suma=0.0_wp
    model_data%Dw%suma=0.0_wp
    
    Do i = 1, fit_data%std_samples%value
      If (fit_data%param%valid(i)) Then
        fit_data%xstd=(i-1)*fit_data%delta_std-fit_data%max_times_std%value
        model_data%theta_cone%suma=model_data%theta_cone%suma+&
                              &(model_data%theta_cone%value(i)-model_data%theta_cone%mean)**2 * exp(-0.5*fit_data%xstd**2)
        model_data%theta_lib%suma=model_data%theta_lib%suma+&
                              &(model_data%theta_lib%value(i)-model_data%theta_lib%mean)**2 * exp(-0.5*fit_data%xstd**2)
        model_data%Dw%suma=model_data%Dw%suma+(model_data%Dw%value(i)-model_data%Dw%mean)**2 * exp(-0.5*fit_data%xstd**2)
      End If
    End Do
    
    model_data%theta_cone%std=sqrt(model_data%theta_cone%suma/fit_data%sum_weights)    
    model_data%theta_lib%std=sqrt(model_data%theta_lib%suma/fit_data%sum_weights)     
    model_data%Dw%std=sqrt(model_data%Dw%suma/fit_data%sum_weights)
    model_data%Tau_cone%std=fit_data%param%T2%std

    Write(messages(1), '(1x,a)') 'The following quantities are derived from the wobbling-in-a-cone&
                                 & model ['//Trim(lipari1980)//']:'
    If (Trim(fit_data%functional_form) == 'ocf_2exp_to_constant') Then  
       Write(messages(2), '(3x,a)') '* Final cone angle (theta_cone) and librational angle (theta_lib); in degrees.'
       Write(messages(3), '(3x,a)') '* Diffussion associated with the final cone angle wobbling (Dw); in 1/ps.'
       Write(messages(4), '(3x,a)') '* Effective relaxation time for the wobbling cone (Tau_cone)&
                                    & and from the host (Tau_host); in ps.'
    Else If (Trim(fit_data%functional_form) == 'ocf_2exp_to_zero') Then  
       Write(messages(2), '(3x,a)') '* Final cone angle (theta_cone) and librational angle (theta_lib); in degrees.'
       Write(messages(3), '(3x,a)') '* Diffussion associated with the final cone angle wobbling (Dw); in 1/ps.'
       Write(messages(4), '(3x,a)') '* Effective relaxation time for the wobbling cone (Tau_cone)&
                                    & and from the host (Tau_host); in ps.'
    Else If (Trim(fit_data%functional_form) == 'ocf_1exp_to_constant') Then  
       Write(messages(2), '(3x,a)') '* Final cone angle (theta_cone) and librational angle (theta_lib); in degrees.'
       Write(messages(3), '(3x,a)') '* Diffussion associated with the final cone angle wobbling (Dw); in 1/ps.'
       Write(messages(4), '(3x,a)') '* Effective relaxation time for the wobbling cone (Tau_cone); in ps.'
    Else If (Trim(fit_data%functional_form) == 'ocf_1exp_to_zero') Then  
       Write(messages(2), '(3x,a)') '* Librational angle (theta_lib); in degrees. The final cone angle is equal to&
                                     & 90 degrees for this option of the "fitting_function" directive.'
       Write(messages(3), '(3x,a)') '* Diffussion associated with the final cone angle wobbling (Dw); in 1/ps.'
       Write(messages(4), '(3x,a)') '* Effective relaxation time for the wobbling cone (Tau_cone); in ps.'
    End If
    Call info(messages, 4)

    If (fit_data%fit_superfast_time%stat) Then
      Write(messages(1), '(3x,a)') '* Relaxation time for librational motion (Tau_lib); in ps.'
    Call info(messages, 1)
    End If

    Write(messages(1), '(1x,a)') 'The user is referred to the attached notes for details.'
    Write(messages(2), '(a)') ' ------------------------------------------'    
    If (fit_data%std_samples%value /=1) Then
      Write(messages(3), '(1x,a,13x,a,10x,a)') 'Quantity', 'MEAN', 'STD'
    Else
      Write(messages(3), '(1x,a,11x,a)') 'Quantity', 'MEAN (STD not computed)'
    End If    
      Write(messages(4), '(a)') ' ------------------------------------------'    
    Call info(messages, 4)    
    
       
    If (Trim(fit_data%functional_form) == 'ocf_2exp_to_constant' .Or. &
        Trim(fit_data%functional_form) == 'ocf_2exp_to_zero') Then
      model_data%Tau_host%mean=fit_data%param%T3%mean  
      model_data%Tau_host%std=fit_data%param%T3%std
      Call print_quantity(messages(1), fit_data%std_samples%value, 'theta_lib ', model_data%theta_lib)
      Call print_quantity(messages(2), fit_data%std_samples%value, 'theta_cone', model_data%theta_cone)
      Call print_quantity(messages(3), fit_data%std_samples%value, 'Dw        ', model_data%Dw)
      Call print_quantity(messages(4), fit_data%std_samples%value, 'Tau_cone  ', model_data%Tau_cone)
      Call print_quantity(messages(5), fit_data%std_samples%value, 'Tau_host  ', model_data%Tau_host)
      Call info(messages, 5)
    Else If (Trim(fit_data%functional_form) == 'ocf_1exp_to_constant') Then
      Call print_quantity(messages(1), fit_data%std_samples%value, 'theta_lib ', model_data%theta_lib)
      Call print_quantity(messages(2), fit_data%std_samples%value, 'theta_cone', model_data%theta_cone)
      Call print_quantity(messages(3), fit_data%std_samples%value, 'Dw        ', model_data%Dw)
      Call print_quantity(messages(4), fit_data%std_samples%value, 'Tau_cone  ', model_data%Tau_cone)
      Call info(messages, 4)
    Else If (Trim(fit_data%functional_form) == 'ocf_1exp_to_zero') Then
      Call print_quantity(messages(1), fit_data%std_samples%value, 'theta_lib ', model_data%theta_lib)
      Call print_quantity(messages(2), fit_data%std_samples%value, 'Dw        ', model_data%Dw)
      Call print_quantity(messages(3), fit_data%std_samples%value, 'Tau_cone  ', model_data%Tau_cone)
      Call info(messages, 3)
    End If
    
    If (fit_data%fit_superfast_time%stat) Then
      model_data%Tau_lib%mean=fit_data%param%T1%mean  
      model_data%Tau_lib%std=fit_data%param%T1%std
          
      Call print_quantity(messages(1), fit_data%std_samples%value, 'Tau_lib   ', model_data%Tau_lib)
     Call info(messages, 1)
    End If
   
    Call info(' ------------------------------------------', 1)

    ! Check the validity of computed quantities
    Call compare_wrt_time_domain(model_data%Tau_cone, fit_data, model_data%caution)
    
    If (fit_data%std_samples%value > 1) Then
      Call mean_vs_std(model_data%theta_lib, model_data%caution)
      Call mean_vs_std(model_data%tau_cone, model_data%caution)
    End If
    
    If (Trim(fit_data%functional_form) /= 'ocf_1exp_to_zero') Then
      Call compare_mean_values(model_data%theta_lib, model_data%theta_cone, model_data%caution)
      If (fit_data%std_samples%value > 1) Then
        Call mean_vs_std(model_data%theta_cone, model_data%caution)
        Call compare_bounds(model_data%theta_lib, model_data%theta_cone, model_data%caution)  
      End If
    End If
    
    If (Trim(fit_data%functional_form) == 'ocf_2exp_to_constant' .Or. &
        Trim(fit_data%functional_form) == 'ocf_2exp_to_zero') Then
      Call compare_wrt_time_domain(model_data%Tau_host, fit_data, model_data%caution)
      Call compare_mean_values(model_data%Tau_cone, model_data%Tau_host, model_data%caution)
      If (fit_data%std_samples%value > 1) Then
        Call mean_vs_std(model_data%Tau_host, model_data%caution)
        Call compare_bounds(model_data%Tau_cone, model_data%Tau_host, model_data%caution)
      End If
    End If

    If (fit_data%fit_superfast_time%stat) Then
      Call compare_wrt_time_domain(model_data%Tau_lib, fit_data, model_data%caution)
      Call compare_mean_values(model_data%Tau_lib, model_data%Tau_cone, model_data%caution)
      If (fit_data%std_samples%value > 1) Then
         Call mean_vs_std(model_data%tau_lib, model_data%caution)
         Call compare_bounds(model_data%Tau_lib, model_data%Tau_cone, model_data%caution)
      End If
    End If
    
    If (model_data%caution) Then
      Call advise_user(fit_data)
    Else
          Write(messages(1), '(1x,a)') '***NOTE: Derived quantities are expected to depend on whether the user chooses'
          Write(messages(2), '(1x,a)') '   a fitting function that tends to a non-zero constant or decays to zero.'
          Call info(messages, 2)
    
      If (Trim(fit_data%functional_form) == 'ocf_2exp_to_constant' .Or. &
          Trim(fit_data%functional_form) == 'ocf_2exp_to_zero') Then
          Write(messages(1), '(1x,a)') '***NOTE: care must be paid when using two exponentials to interpret OCF profiles.'
          Write(messages(2), '(1x,a)') '   Particularly, "Tau_host" must be much larger than "Tau_cone", otherwise results for "Dw"'
          Write(messages(3), '(1x,a)') '   are expected to differ significantly from those obtained using one exponential.'
          Call info(messages, 3)
      End If    
    End If    
    
  End Subroutine compute_wobbling_in_cone  
  
  Subroutine compute_diffusion_coefficient(fit_data, model_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to compute the average compute_diffusion coefficient and STD
    ! from fits.
    !
    ! author    - i.scivetti March 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(fit_type),   Intent(InOut) :: fit_data
    Type(model_type), Intent(InOut) :: model_data

    Character(Len=256)  :: messages(5), dim
    Integer        :: i
 
    Call model_data%alloc_diffusion(fit_data%std_samples%value)

    model_data%D%suma=0.0_wp
    
    Do i = 1, fit_data%std_samples%value
      If (fit_data%param%valid(i)) Then
        fit_data%xstd=(i-1)*fit_data%delta_std-fit_data%max_times_std%value
        model_data%D%value(i)=fit_data%param%A1%value(i)/2.0_wp/fit_data%msd%dim
        model_data%D%suma=model_data%D%suma+model_data%D%value(i)*exp(-0.5*fit_data%xstd**2)
      End If
    End Do
 
    model_data%D%mean=model_data%D%suma/fit_data%sum_weights    

    model_data%D%suma=0.0_wp
    
    Do i = 1, fit_data%std_samples%value
      If (fit_data%param%valid(i)) Then
        fit_data%xstd=(i-1)*fit_data%delta_std-fit_data%max_times_std%value
        model_data%D%suma=model_data%D%suma+(model_data%D%value(i)-model_data%D%mean)**2 * exp(-0.5*fit_data%xstd**2)
      End If
    End Do
    
    model_data%D%std=sqrt(model_data%D%suma/fit_data%sum_weights)
    
    ! Define labels
    model_data%D%label='D'

    Write (dim,'(i1)') fit_data%msd%dim
    Write(messages(1), '(1x,a)') "The diffusion coefficient (D) within the approximation of the Einstein's relation"
    Write(messages(2), '(1x,a)') 'is computed as D=A1/(2*dim), where dim in this case is equal to '//Trim(dim)//&
                                &' (see above).'
    Write(messages(3), '(1x,a)') '(D is un units of Angstrom^2/ps)'
    Call info(messages, 3)
    
    Write(messages(1), '(a)') ' ---------------------------------'    
    If (fit_data%std_samples%value /=1) Then
      Write(messages(2), '(1x,a,5x,a,8x,a)') 'Quantity', 'MEAN', 'STD'
    Else
      Write(messages(2), '(1x,a,2x,a)') 'Quantity', 'MEAN (STD not computed)'
    End If    
    Write(messages(3), '(a)') ' ---------------------------------'    
    Call print_quantity(messages(4), fit_data%std_samples%value, 'D', model_data%D)
    Write(messages(5), '(a)') ' ---------------------------------'
    Call info(messages, 5)

    If (fit_data%std_samples%value > 1) Then
      Call mean_vs_std(model_data%D, model_data%caution)
    End If
    
    If (model_data%caution) Then
      Call advise_user(fit_data)
    End If
    
  End Subroutine compute_diffusion_coefficient
  
  
  Subroutine compute_first_order_kinetics(fit_data, i, model_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to compute rates
    !
    ! author    - i.scivetti March 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(fit_type),   Intent(InOut) :: fit_data
    Integer(Kind=wi), Intent(In   ) :: i
    Type(model_type), Intent(InOut) :: model_data
    
    Real(Kind=wp)  :: A0, A1, A2

    A1 =fit_data%param%A1%value(i)

    If (Trim(fit_data%functional_form) == 'tcf_2exp_to_zero') Then  
      A0 =fit_data%param%A0%value(i)
      A2 =fit_data%param%A2%value(i)
    Else If (Trim(fit_data%functional_form) == 'tcf_1exp_to_constant') Then  
      A0 =fit_data%param%A0%value(i)
      A2 = 0.0_wp
    Else If (Trim(fit_data%functional_form) == 'tcf_1exp_to_zero') Then  
      A0 = 1.0_wp
      A2 = 0.0_wp
    Else If (Trim(fit_data%functional_form) == 'spcf_2exp_to_zero') Then  
      A0 =fit_data%param%A0%value(i)
      A2 =fit_data%param%A2%value(i)
    Else If (Trim(fit_data%functional_form) == 'spcf_1exp_to_constant') Then  
      A0 =fit_data%param%A0%value(i)
      A2 = 0.0_wp
    Else If (Trim(fit_data%functional_form) == 'spcf_1exp_to_zero') Then  
      A0 = 1.0_wp
      A2 = 0.0_wp      
    End If
  
    ! Values
    model_data%kf%value(i)=A0*A1+A2*(1.0_wp-A0)
    model_data%kb%value(i)=((A1-A2)**2 * A0 * (1.0_wp-A0))/model_data%kf%value(i)
    model_data%k2%value(i)=A1*A2/model_data%kf%value(i)
  
  End Subroutine compute_first_order_kinetics
  
  Subroutine compute_rates(fit_data, model_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to compute average rates and associated errors from the fitting
    !
    ! author    - i.scivetti March 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(fit_type),   Intent(InOut) :: fit_data
    Type(model_type), Intent(InOut) :: model_data

    Character(Len=256)  :: messages(6)
    Integer        :: i
 
    Call model_data%alloc_kinetics(fit_data%std_samples%value)

    ! Define labels
    model_data%kf%label=    'kf'
    model_data%kb%label=    'kb'
    model_data%k2%label=    'k2'   
    model_data%inv_kf%label='1/(kf)'
    model_data%inv_kb%label='1/(kb)'
    model_data%inv_k2%label='1/(k2)'
    
    model_data%kf%suma=0.0_wp
    model_data%kb%suma=0.0_wp
    model_data%k2%suma=0.0_wp
    
    Do i = 1, fit_data%std_samples%value
      If (fit_data%param%valid(i)) Then
        fit_data%xstd=(i-1)*fit_data%delta_std-fit_data%max_times_std%value
        Call compute_first_order_kinetics(fit_data, i, model_data)
        model_data%kf%suma=model_data%kf%suma+model_data%kf%value(i)*exp(-0.5*fit_data%xstd**2)
        model_data%kb%suma=model_data%kb%suma+model_data%kb%value(i)*exp(-0.5*fit_data%xstd**2)
        model_data%k2%suma=model_data%k2%suma+model_data%k2%value(i)*exp(-0.5*fit_data%xstd**2)
      End If
    End Do
 
    model_data%kf%mean=model_data%kf%suma/fit_data%sum_weights    
    model_data%kb%mean=model_data%kb%suma/fit_data%sum_weights     
    model_data%k2%mean=model_data%k2%suma/fit_data%sum_weights

    model_data%kf%suma=0.0_wp
    model_data%kb%suma=0.0_wp
    model_data%k2%suma=0.0_wp
    
    Do i = 1, fit_data%std_samples%value
      If (fit_data%param%valid(i)) Then
        fit_data%xstd=(i-1)*fit_data%delta_std-fit_data%max_times_std%value
        model_data%kf%suma=model_data%kf%suma+(model_data%kf%value(i)-model_data%kf%mean)**2 * exp(-0.5*fit_data%xstd**2)
        model_data%kb%suma=model_data%kb%suma+(model_data%kb%value(i)-model_data%kb%mean)**2 * exp(-0.5*fit_data%xstd**2)
        model_data%k2%suma=model_data%k2%suma+(model_data%k2%value(i)-model_data%k2%mean)**2 * exp(-0.5*fit_data%xstd**2)
      End If
    End Do
    
    model_data%kf%std=sqrt(model_data%kf%suma/fit_data%sum_weights)    
    model_data%kb%std=sqrt(model_data%kb%suma/fit_data%sum_weights)     
    model_data%k2%std=sqrt(model_data%k2%suma/fit_data%sum_weights)

    If (Trim(fit_data%what_to_fit%type) == 'tcf') Then
      Write(messages(1), '(1x,a)') 'Inverse transfer rates for the "'//Trim(fit_data%species_name%type)//'" sites&
                                  & (computed with the first-order kinetics approximation)'
    Else If (Trim(fit_data%what_to_fit%type) == 'spcf') Then
      Write(messages(1), '(1x,a)') 'Inverse transfer rates for the special pair formed between the&
                                  & "'//Trim(fit_data%species_name%type)//'" sites and neighbouring species&
                                  & (computed with the first-order kinetics approximation)'
    End If
    
    Write(messages(2), '(1x,a)') 'The following quantities are derived ['//Trim(morrone2006)//']:'
    
    If (Trim(fit_data%functional_form) == 'tcf_2exp_to_zero' .Or. Trim(fit_data%functional_form) == 'spcf_2exp_to_zero') Then
       Write(messages(3), '(3x,a)') '* Forward  rate constant (kf) in 1/ps; 1/(kf) in ps.'
       Write(messages(4), '(3x,a)') '* Backward rate constant (kb) in 1/ps; 1/(kb) in ps.'
       Write(messages(5), '(3x,a)') '* Second forward rate constant (k2) in 1/ps; 1/(k2) in ps.'       
       Call info(messages, 5)
    Else If (Trim(fit_data%functional_form) == 'tcf_1exp_to_constant' .Or. &
           & Trim(fit_data%functional_form) == 'spcf_1exp_to_constant') Then
       Write(messages(3), '(3x,a)') '* Forward  rate constant (kf) in 1/ps; 1/(kf) in ps.'
       Write(messages(4), '(3x,a)') '* Backward rate constant (kb) in 1/ps; 1/(kb) in ps.'
       Call info(messages, 4)    
    Else If (Trim(fit_data%functional_form) == 'tcf_1exp_to_zero' .Or. Trim(fit_data%functional_form) == 'spcf_1exp_to_zero') Then  
       Write(messages(3), '(3x,a)') '* Forward rate constant (kf) in 1/ps; 1/(kf) in ps.'
       Call info(messages, 3)    
    End If

    Write(messages(1), '(1x,a)') 'The user is referred to the attached notes for details.'
    Write(messages(2), '(a)') ' --------------------------------------'    
    If (fit_data%std_samples%value /=1) Then
      Write(messages(3), '(1x,a,10x,a,9x,a)') 'Quantity', 'MEAN', 'STD'
    Else
      Write(messages(3), '(1x,a,7x,a)') 'Quantity', 'MEAN (STD not computed)'
    End If    
    Write(messages(4), '(a)') ' --------------------------------------'    
    Call info(messages, 4)        

    ! Define inverse of kf
    model_data%inv_kf%mean=1.0_wp/model_data%kf%mean    
    model_data%inv_kf%std=model_data%kf%std/(model_data%kf%mean)**2
    
    If (Trim(fit_data%functional_form) == 'tcf_2exp_to_zero' .Or. Trim(fit_data%functional_form) == 'spcf_2exp_to_zero') Then  
      model_data%inv_kb%mean=1.0_wp/model_data%kb%mean    
      model_data%inv_k2%mean=1.0_wp/model_data%k2%mean     
      model_data%inv_kb%std=model_data%kb%std/(model_data%kb%mean)**2
      model_data%inv_k2%std=model_data%k2%std/(model_data%k2%mean)**2      
      Call print_quantity(messages(1), fit_data%std_samples%value, 'kf    ', model_data%kf)
      Call print_quantity(messages(2), fit_data%std_samples%value, 'kb    ', model_data%kb)
      Call print_quantity(messages(3), fit_data%std_samples%value, 'k2    ', model_data%k2)
      Call print_quantity(messages(4), fit_data%std_samples%value, '1/(kf)', model_data%inv_kf)
      Call print_quantity(messages(5), fit_data%std_samples%value, '1/(kb)', model_data%inv_kb)
      Call print_quantity(messages(6), fit_data%std_samples%value, '1/(k2)', model_data%inv_k2)
      Call info(messages, 6)
    Else If (Trim(fit_data%functional_form) == 'tcf_1exp_to_constant' .Or. &
           & Trim(fit_data%functional_form) == 'spcf_1exp_to_constant') Then  
      model_data%inv_kb%mean=1.0_wp/model_data%kb%mean    
      model_data%inv_kb%std=model_data%kb%std/(model_data%kb%mean)**2
      Call print_quantity(messages(1), fit_data%std_samples%value, 'kf    ', model_data%kf)
      Call print_quantity(messages(2), fit_data%std_samples%value, 'kb    ', model_data%kb)
      Call print_quantity(messages(3), fit_data%std_samples%value, '1/(kf)', model_data%inv_kf)
      Call print_quantity(messages(4), fit_data%std_samples%value, '1/(kb)', model_data%inv_kb)      
      Call info(messages, 4)
    Else If (Trim(fit_data%functional_form) == 'tcf_1exp_to_zero' .Or. Trim(fit_data%functional_form) == 'spcf_1exp_to_zero') Then  
      Call print_quantity(messages(1), fit_data%std_samples%value, 'kf    ', model_data%kf)
      Call print_quantity(messages(2), fit_data%std_samples%value, '1/(kf)', model_data%inv_kf)
      Call info(messages, 2)
    End If
    Call info(' --------------------------------------', 1)

    If (fit_data%std_samples%value > 1) Then
      Call mean_vs_std(model_data%inv_kf, model_data%caution)
    End If
 
    Call compare_wrt_time_domain(model_data%inv_kf, fit_data, model_data%caution)
    
    If (Trim(fit_data%functional_form) == 'tcf_2exp_to_zero'      .Or. &
        Trim(fit_data%functional_form) == 'tcf_1exp_to_constant'  .Or. &
        Trim(fit_data%functional_form) == 'spcf_2exp_to_zero'     .Or. &
        Trim(fit_data%functional_form) == 'spcf_1exp_to_constant') Then  

      Call compare_wrt_time_domain(model_data%inv_kb, fit_data, model_data%caution) 
      Call compare_mean_values(model_data%inv_kf, model_data%inv_kb, model_data%caution)
      
      
      If (fit_data%std_samples%value > 1) Then
        Call mean_vs_std(model_data%inv_kb, model_data%caution)
        Call compare_bounds(model_data%inv_kf, model_data%inv_kb, model_data%caution)
      End If
      
      If (Trim(fit_data%functional_form) == 'tcf_2exp_to_zero' .Or. Trim(fit_data%functional_form) == 'spcf_2exp_to_zero') Then 
        Call compare_wrt_time_domain(model_data%inv_k2, fit_data, model_data%caution) 
        Call compare_mean_values(model_data%inv_kf, model_data%inv_k2, model_data%caution)
        Call compare_mean_values(model_data%inv_kb, model_data%inv_k2, model_data%caution)
        
        If (fit_data%std_samples%value > 1) Then
          Call mean_vs_std(model_data%inv_k2, model_data%caution)
          Call compare_bounds(model_data%inv_kf, model_data%inv_k2, model_data%caution)
          Call compare_bounds(model_data%inv_kb, model_data%inv_k2, model_data%caution)
        End If
      End If
    End If    
    
    If (model_data%caution) Then
      Call advise_user(fit_data)
    End If
    
  End Subroutine compute_rates  

  Subroutine advise_user(fit_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to advise the the user if computed quantities exhibit issues
    !
    ! author    - i.scivetti March 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(fit_type),     Intent(InOut) :: fit_data  

    Character(Len=256)  :: messages(2)
    
    Write(messages(1), '(1x,a)') ' '    
    Write(messages(2), '(1x,a)') '***PROBLEMS*** Computed quantities appear to have issues.'
    Call info(messages, 2)
    Write(messages(1), '(1x,a)') 'ADVISE 1: review the choice of the formula to fit the data&
                                & (option of the "fitting_function" directive)'
    Write(messages(2), '(1x,a)') 'ADVISE 2: check the input data and the time domain for fitting'
    Call info(messages, 2)
  
    If (fit_data%std_samples%value /= 1) Then
      Write(messages(1), '(1x,a)') 'ADVISE 3: consider reducing the value of "max_times_std"'
      Call info(messages, 1)

      If (fit_data%input_parameters%fread) Then
         If (fit_data%valid_solutions /= fit_data%std_samples%value) Then
          Write(messages(1), '(a)') 'IMPORTANT: revise the values defined in the "&input_parameters" block'
          Call info(messages, 1)
         End If
      End If
    End If    
    
  End Subroutine advise_user

  Subroutine compare_wrt_time_domain(T, fit_data, caution)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to compare time-related quantity against the time limit of the 
    ! data fitted
    !
    ! author    - i.scivetti March 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(quantity),    Intent(In   ) :: T
    Type(fit_type),    Intent(In   ) :: fit_data
    Logical,           Intent(InOut) :: caution  

    Character(Len=256)  :: message

    If (T%mean > fit_data%t_limit) Then
      caution=.True.
      If (fit_data%end_time%fread) Then
        Write(message, '(1x,a)') 'CAUTION: the mean value of "'//Trim(T%label)//'" is larger than the "end_time" directive!&
                                     & THIS COMPROMISES THE VALITIY CONDITIONS!'
      Else
        Write(message, '(1x,a)') 'CAUTION: the mean value of "'//Trim(T%label)//'" is larger than the largest recorded time!&
                                     & THIS COMPROMISES THE VALITIY CONDITIONS!'
      End If
      Call info(message, 1)                
    End If    
    
  End Subroutine compare_wrt_time_domain 
 
  Subroutine compare_mean_values(T1, T2, caution)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to compare two mean values
    !
    ! author    - i.scivetti March 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(quantity),    Intent(In   ) :: T1
    Type(quantity),    Intent(In   ) :: T2
    Logical,           Intent(InOut) :: caution  

    Character(Len=256)  :: message

    If (T1%mean > T2%mean) Then
      caution=.True.
      Write(message, '(1x,a)') 'CAUTION: "'//Trim(T1%label)//'" is larger than "'//Trim(T2%label)//'"! THIS SHOULD NOT HAPPEN!'
      Call info(message, 1)                
    End If
    
  End Subroutine compare_mean_values
 
  Subroutine compare_bounds(T1, T2, caution)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to compare bounds between two quantities
    !
    ! author    - i.scivetti March 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(quantity),    Intent(In   ) :: T1
    Type(quantity),    Intent(In   ) :: T2
    Logical,           Intent(InOut) :: caution  

    Character(Len=256)  :: message

    If ((T1%mean+T1%std) > (T2%mean-T2%std)) Then
      caution=.True.
      Write(message, '(1x,a)') 'CAUTION: the upper bound for "'//Trim(T1%label)//'" is larger than lower bound for "'&
                               &//Trim(T2%label)//'"'
      Call info(message, 1)                
    End If
    
  End Subroutine compare_bounds
 
  Subroutine mean_vs_std(T, caution)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to compare the mean and the std
    !
    ! author    - i.scivetti March 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(quantity),    Intent(In   ) :: T
    Logical,           Intent(InOut) :: caution  
    
    Character(Len=256)  :: message

    If (T%mean < T%std) Then
      caution=.True.
      Write(message, '(1x,a)') 'CAUTION: the STD for "'//Trim(T%label)//'" is larger than the MEAN value!'
      Call info(message, 1)                
    End If
    
  End Subroutine mean_vs_std
 
  Subroutine print_quantity(line, nfits, descriptor, T) 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to print computed quantities (averages and stds)
    !
    ! author    - i.scivetti March 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
    Character(Len=256), Intent(  Out)  :: line
    Integer(Kind=wi),   Intent(In   )  :: nfits
    Character(*),       Intent(In   )  :: descriptor
    Type(quantity),     Intent(In   )  :: T
    
    If (nfits == 1) Then
      Write(line, fmt_mean_only) descriptor, T%mean
    Else
      Write(line, fmt_mean_std)  descriptor, T%mean, T%std
    End If
  
  End Subroutine print_quantity
  
End module models
