!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Module for the setting up of variables and procedures for fitting
!
! Copyright - 2026 Ada Lovelace Centre (ALC)
!             Scientific Computing Department (SCD)
!             The Science and Technology Facilities Council (STFC)
!
! Author:        i.scivetti  March 2026
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Module fitting_setup 

   Use constants,    Only : max_std_samples, &
                            fmt_mean_std, &
                            fmt_mean_only
   
                            
   Use input_types,  Only : in_integer, &
                            in_logic,   &
                            in_scalar,  &
                            in_param,   & 
                            in_string
 
   Use numprec,      Only : wi,&
                            wp 

  Use unit_output,   Only : error_stop, &
                            info

  Implicit None
  Private

  Type raw_type
    Real(Kind=wp) :: time
    Real(Kind=wp) :: avg
    Real(Kind=wp) :: std
  End Type
  
  Type param_type
    Real(Kind=wp), Allocatable :: value(:)
    Character(Len=256), Allocatable :: fit_error(:)
    Real(Kind=wp)  :: suma
    Real(Kind=wp)  :: mean
    Real(Kind=wp)  :: std
    Type(in_param) :: ini
  End Type

  Type set_param_type
    Type(param_type)  :: A3
    Type(param_type)  :: A2
    Type(param_type)  :: A1
    Type(param_type)  :: A0
    Type(param_type)  :: T1
    Type(param_type)  :: T2
    Type(param_type)  :: T3
    Real(Kind=wp), Allocatable :: ssr(:)
    Logical, Allocatable :: fit_ok(:)
    Logical, Allocatable :: valid(:)
  End Type

  !Type to describe the msd
  Type :: msd_type
    Type(in_string)  :: invoke
    Type(in_string)  :: units
    Type(in_string)  :: select
    Real(Kind=wp)    :: factor
    Integer(Kind=wi) :: dim
  End Type
  
  Type, Public :: fit_type
    Logical              :: caution_model
    Character(Len=256)   :: formula 
    Character(Len=256)   :: formula_show
    Character(Len=256)   :: quantity
    Character(Len=256)   :: set_fit
    Character(Len=256)   :: init_param
    Character(Len=256)   :: instruct_fitting
    Character(Len=256)   :: exec_gnuplot
    Character(Len=256)   :: check_fit
    Character(Len=256)   :: check_ndf
    Character(Len=256)   :: check_window
    Character(Len=256)   :: extract_results
    Character(Len=256)   :: via_params
    Character(Len=256)   :: dump_params
    Character(Len=256)   :: param_list
    Character(Len=256)   :: header
    Character(Len=256)   :: separation_line
    Character(Len=256)   :: ssr_string
    Character(Len=256)   :: results_string
    Character(Len=256)   :: xstd_char
    Character(Len=256)   :: functional_form
    Character(Len=256)   :: rm_gnuplot_output
    Character(Len=256)   :: rm_gnuplot_input
    Real(Kind=wp)        :: xstd
    Real(Kind=wp)        :: delta_std
    Real(Kind=wp)        :: sum_weights
    Real(Kind=wp)        :: t_limit
    Real(Kind=wp)        :: average_ssr
    Integer(Kind=wi)     :: ndf
    Integer(Kind=wi)     :: nparam
    Integer(Kind=wi)     :: valid_solutions
    Type(in_string)      :: filename
    Type(in_string)      :: input_parameters
    Type(in_string)      :: what_to_fit
    Type(in_string)      :: fitting_function
    Type(in_string)      :: species_name
    Type(in_param)       :: max_times_std
    Type(in_integer)     :: std_samples
    Type(in_param)       :: start_time
    Type(in_param)       :: end_time
    Type(in_integer)     :: max_iterations
    Type(in_logic)       :: reactive_species
    Type(in_logic)       :: plot_fittings
    Type(in_logic)       :: plot_raw_data
    Type(in_logic)       :: print_average_ssr
    Type(in_logic)       :: fit_superfast_time
    Type(set_param_type) :: param
    Type(msd_type)       :: msd
    Integer(Kind=wi)     :: num_input_data
    Type(in_integer)     :: time_column
    Type(in_integer)     :: avg_column
    Type(in_integer)     :: std_column
    Type(raw_type), Allocatable :: input_data(:)
    Real(Kind=wp), Allocatable  :: raw_data(:,:)
    Real(Kind=wp), Allocatable  :: profile(:,:)
  Contains
    Private
      Procedure, Public    :: alloc_parameters  => allocate_fitting_parameters
      Procedure, Public    :: alloc_data        => allocate_data_arrays_for_printing
      Final                :: cleanup
  End Type

  Public ::  print_parameter
  
Contains

  Subroutine allocate_fitting_parameters(T)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Allocate arrays for fitting
    !
    ! author    - i.scivetti March 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Class(fit_type),   Intent(InOut)  :: T
    
    Integer(Kind=wi)    :: fail(17)
    Character(Len=256)  :: message

    fail=0
    Write (message,'(1x,1a)') '***ERROR: Allocation problems for fitting parameters&
                                & (subroutine allocate_fitting_parameters). Check the value of "std_samples"'
                                
    Allocate(T%param%A3%value(T%std_samples%value), Stat=fail(1))
    Allocate(T%param%A2%value(T%std_samples%value), Stat=fail(2))
    Allocate(T%param%A1%value(T%std_samples%value), Stat=fail(3))
    Allocate(T%param%A0%value(T%std_samples%value), Stat=fail(4))
    Allocate(T%param%T1%value(T%std_samples%value), Stat=fail(5))
    Allocate(T%param%T2%value(T%std_samples%value), Stat=fail(6))
    Allocate(T%param%T3%value(T%std_samples%value), Stat=fail(7))
    
    Allocate(T%param%A3%fit_error(T%std_samples%value), Stat=fail(8))
    Allocate(T%param%A2%fit_error(T%std_samples%value), Stat=fail(9))
    Allocate(T%param%A1%fit_error(T%std_samples%value), Stat=fail(10))
    Allocate(T%param%A0%fit_error(T%std_samples%value), Stat=fail(11))
    Allocate(T%param%T1%fit_error(T%std_samples%value), Stat=fail(12))
    Allocate(T%param%T2%fit_error(T%std_samples%value), Stat=fail(13))
    Allocate(T%param%T3%fit_error(T%std_samples%value), Stat=fail(14))

    Allocate(T%param%fit_ok(T%std_samples%value), Stat=fail(15))
    Allocate(T%param%valid(T%std_samples%value),  Stat=fail(16))
    Allocate(T%param%ssr(T%std_samples%value),    Stat=fail(17))
    
    If (Any(fail > 0)) Then
      Call error_stop(message)
    Else 
      T%param%fit_ok=.True.
      T%param%A3%value=0.0_wp
      T%param%A2%value=0.0_wp
      T%param%A1%value=0.0_wp
      T%param%A0%value=0.0_wp
      T%param%T1%value=0.0_wp
      T%param%T2%value=0.0_wp
      T%param%T3%value=0.0_wp
    End If

  End Subroutine allocate_fitting_parameters

  Subroutine allocate_data_arrays_for_printing(T)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Allocate arrays for the raw data to fit
    !
    ! author    - i.scivetti March 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Class(fit_type),   Intent(InOut)  :: T
    
    Integer(Kind=wi)    :: fail
    Character(Len=256)  :: message

    fail=0
    
    Write (message,'(1x,1a)') '***ERROR: Allocation problems for kinetcs&
                                & (subroutine allocate_data_arrays_for_printing).&
                                & Check the amount of data of the input file'
                                
    Allocate(T%input_data(T%num_input_data), Stat=fail)
    
    If (fail > 0) Then
      Call error_stop(message)
    End If
    
    If (T%plot_raw_data%stat) Then
      Allocate(T%raw_data(T%num_input_data, max_std_samples), Stat=fail)
      If (fail > 0) Then
        Call error_stop(message)
      End If
    End If

    If (T%plot_fittings%stat) Then
      Allocate(T%profile(T%num_input_data, max_std_samples), Stat=fail)
      If (fail > 0) Then
        Call error_stop(message)
      End If
    End If
    
  End Subroutine allocate_data_arrays_for_printing
  
 
  Subroutine cleanup(T)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Deallocate variables
    !
    ! author    - i.scivetti March 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(fit_type) :: T

    If (Allocated(T%param%A3%value)) Then
      Deallocate(T%param%A3%value)
    End If 
    
    If (Allocated(T%param%A2%value)) Then
      Deallocate(T%param%A2%value)
    End If 

    If (Allocated(T%param%A1%value)) Then
      Deallocate(T%param%A1%value)
    End If 

    If (Allocated(T%param%A0%value)) Then
      Deallocate(T%param%A0%value)
    End If 

    If (Allocated(T%param%T1%value)) Then
      Deallocate(T%param%T1%value)
    End If 
    
    If (Allocated(T%param%T2%value)) Then
      Deallocate(T%param%T2%value)
    End If 
    
    If (Allocated(T%param%T3%value)) Then
      Deallocate(T%param%T3%value)
    End If 
    
    If (Allocated(T%param%A3%fit_error)) Then
      Deallocate(T%param%A3%fit_error)
    End If 
    
    If (Allocated(T%param%A2%fit_error)) Then
      Deallocate(T%param%A2%fit_error)
    End If 

    If (Allocated(T%param%A1%fit_error)) Then
      Deallocate(T%param%A1%fit_error)
    End If 

    If (Allocated(T%param%A0%fit_error)) Then
      Deallocate(T%param%A0%fit_error)
    End If 
    
    If (Allocated(T%param%T1%fit_error)) Then
      Deallocate(T%param%T1%fit_error)
    End If 
    
    If (Allocated(T%param%T2%fit_error)) Then
      Deallocate(T%param%T2%fit_error)
    End If 

    If (Allocated(T%param%T3%fit_error)) Then
      Deallocate(T%param%T3%fit_error)
    End If 

    If (Allocated(T%param%fit_ok)) Then
      Deallocate(T%param%fit_ok)
    End If 
    
    If (Allocated(T%param%valid)) Then
      Deallocate(T%param%valid)
    End If 
    
    If (Allocated(T%param%ssr)) Then
      Deallocate(T%param%ssr)
    End If 

    If (Allocated(T%input_data)) Then
      Deallocate(T%input_data)
    End If 

    If (Allocated(T%raw_data)) Then
      Deallocate(T%raw_data)
    End If 

    If (Allocated(T%profile)) Then
      Deallocate(T%profile)
    End If 

  End Subroutine cleanup

  Subroutine print_parameter(message, nfits, label, mean, std) 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to print averages and std from the fitted parameters
    !
    ! author    - i.scivetti March 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
    Character(Len=256), Intent(  Out)  :: message
    Integer(Kind=wi),   Intent(In   )  :: nfits
    Character(*),       Intent(In   )  :: label
    Real(Kind=wp),      Intent(In   )  :: mean
    Real(Kind=wp),      Intent(In   )  :: std
    
    If (nfits == 1) Then
      Write(message, fmt_mean_only) label, mean
    Else
      Write(message, fmt_mean_std) label, mean, std
    End If
    
  
  End Subroutine print_parameter  
  
End Module fitting_setup
  
