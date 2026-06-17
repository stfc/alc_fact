!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
!
! Copyright - 2026 Ada Lovelace Centre (ALC)
!             Scientific Computing Department (SCD)
!             The Science and Technology Facilities Council (STFC)
!
! Author:     i.scivetti  March 2026
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Module fitting_checks 

  Use constants,      Only : Bohr_to_A, &
                             exp_ocf_lower_bound,&
                             max_std_samples, &
                             max_iterations_fitting
                       
  Use fileset,        Only : file_type, &
                             FILE_SET
                             
   Use input_types,   Only : in_param
                      
   Use numprec,       Only : wp 
                      
   Use process_data,  Only : capital_to_lower_case 
   
   Use fitting_setup, Only : fit_type
    
   Use unit_output,   Only : error_stop, &
                             info

  Implicit None
  Private

  Public :: check_fitting_settings, check_validity_parameters
  Public :: check_times_involved, check_validity_mean_values

Contains

  Subroutine check_fitting_settings(files, fit_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to check the correctness of fitting-related directives
    !
    ! author    - i.scivetti March 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),    Intent(In   ) :: files(:)
    Type(fit_type),     Intent(InOut) :: fit_data

    Character(Len=256)  :: messages(6)
    Character(Len=64 )  :: error_set, max_samples, word
    Logical             :: filename_exist, freact

    error_set = '***ERROR in file '//Trim(files(FILE_SET)%filename)//' -'

    If (fit_data%what_to_fit%fread) Then
      If (fit_data%what_to_fit%fail) Then
        Write (messages(1),'(2(1x,a))') Trim(error_set), 'Wrong (or missing) settings for the "what_to_fit" directive.'
        Call info(messages, 1)
        Call error_stop(' ')
      Else
        If (Trim(fit_data%what_to_fit%type) /= 'tcf'  .And. &
            Trim(fit_data%what_to_fit%type) /= 'spcf' .And. &
            Trim(fit_data%what_to_fit%type) /= 'msd'  .And. &
            Trim(fit_data%what_to_fit%type) /= 'ocf') Then
            Write (messages(1),'(2(1x,a))')  Trim(error_set), 'Invalid input for the "what_to_fit" directive.'
            Call info(messages, 1)
            Write (messages(1),'(a)')        'Implemented options for "what_to_fit":'
            Write (messages(2),'(a)')        ' - ocf'
            Write (messages(3),'(a)')        ' - tcf'
            Write (messages(4),'(a)')        ' - spcf'
            Write (messages(5),'(a)')        ' - msd'
            Call info(messages, 5)
            Call info(' ', 1)
            Write (messages(1),'(1x,a)') 'What quantity do you want to fit?'
            Call info(messages, 1)
            Call error_stop(' ')
        End If
      End If
    Else
       Write (messages(1),'(2(1x,a))')  Trim(error_set), 'The user must define the "what_to_fit" directive!'
       Call info(messages, 1)
       Write (messages(1),'(a)')  'Implemented options for "what_to_fit":'
       Write (messages(2),'(a)')  ' - ocf'
       Write (messages(3),'(a)')  ' - tcf'
       Write (messages(4),'(a)')  ' - spcf'
       Write (messages(5),'(a)')  ' - msd'
       Call info(messages, 5)
       Call info(' ', 1)
       Write (messages(1),'(1x,a)') 'What quantity do you want to fit?'
       Call info(messages, 1)       
       Call error_stop(' ')
    End If    

    If (fit_data%reactive_species%fread) Then
      If (fit_data%reactive_species%fail) Then
        Write (messages(1),'(2(1x,a))') Trim(error_set), 'Wrong (or missing) setting for the "reactive_species" directive,&
                                      & which must be set either .True. or .False.'
        Call info(messages, 1)
        Call error_stop(' ')
      End If
    Else
      If(Trim(fit_data%what_to_fit%type) == 'ocf') Then
          Write (messages(1),'(2(1x,a))') Trim(error_set),&
                                & 'The option "ocf" for "what_to_fit" requires the "reactive_species" directive,&
                                & which must be set either .True. or .False.'
        Call info(messages, 1)
        Call error_stop(' ')      
      End If    
    End If    

    ! Define species_name
    If (fit_data%species_name%fread) Then
      If (fit_data%species_name%fail) Then
        Write (messages(1),'(2(1x,a))') Trim(error_set), 'Wrong (or missing) settings for the "species_name" directive.'
        Call info(messages, 1)
        Call error_stop(' ')
      End If
    Else
       Write (messages(1),'(2(1x,a))')  Trim(error_set), 'The user must define the "species_name" directive'
       Call info(messages, 1)
       Call error_stop(' ')
    End If        

    If (fit_data%filename%fread) Then
      If (fit_data%filename%fail) Then
        Write (messages(1),'(2(1x,a))') Trim(error_set), 'Wrong (or missing) settings for the "filename" directive.'
        Call info(messages, 1)
        Call error_stop(' ')
      Else
        Inquire(file=Trim(fit_data%filename%type), Exist=filename_exist) 
        If (.Not. filename_exist) Then
           Write (messages(1),'(2(1x,a))') Trim(error_set), 'File "'//Trim(fit_data%filename%type)//'" (defined for the&
                                         & "filename" directive) DOES NOT exist!' 
          Call info(messages, 1)
          Call error_stop(' ')
        End If
      End If
    Else
       Write (messages(1),'(2(1x,a))')  Trim(error_set), 'The user must define the "filename" directive&
                                      & with the name of the file that contains the data to be fitted (in columns).'
       Call info(messages, 1)
       Call error_stop(' ')
    End If
    
    
    If (fit_data%fitting_function%fread) Then
      If (fit_data%fitting_function%fail) Then
        Write (messages(1),'(2(1x,a))') Trim(error_set), 'Wrong (or missing) settings for the "fitting_function" directive.'
        Call info(messages, 1)
        Write (messages(1),'(a)')        'Implemented options for "fitting_function":'
        Write (messages(2),'(a)')        ' - 1exp_to_zero'
        Write (messages(3),'(a)')        ' - 1exp_to_constant'
        Write (messages(4),'(a)')        ' - 2exp_to_zero'
        Write (messages(5),'(a)')        ' - 2exp_to_constant'
        Write (messages(6),'(a)')        ' - linear'
        Call info(messages, 6)        
        Call error_stop(' ')
      End If
    Else
       Write (messages(1),'(2(1x,a))')  Trim(error_set), 'The user must define the "fitting_function" directive'
       Call info(messages, 1)
        Write (messages(1),'(a)')        'Implemented options for "fitting_function":'
        Write (messages(2),'(a)')        ' - 1exp_to_zero'
        Write (messages(3),'(a)')        ' - 1exp_to_constant'
        Write (messages(4),'(a)')        ' - 2exp_to_zero'
        Write (messages(5),'(a)')        ' - 2exp_to_constant'
        Write (messages(6),'(a)')        ' - linear'
        Call info(messages, 6)               
       Call error_stop(' ')
    End If

    If (Trim(fit_data%what_to_fit%type) == 'ocf') Then
      If (Trim(fit_data%fitting_function%type) /= '1exp_to_zero'     .And. &
          Trim(fit_data%fitting_function%type) /= '1exp_to_constant' .And. &
          Trim(fit_data%fitting_function%type) /= '2exp_to_zero'     .And. &
          Trim(fit_data%fitting_function%type) /= '2exp_to_constant') Then
        Write (messages(1),'(2(1x,a))')  Trim(error_set), 'Invalid option for the "fitting_function" directive.'
        Call info(messages, 1)
        Write (messages(1),'(a)')        'Implemented options for OCF:'
        Write (messages(2),'(a)')        ' - 1exp_to_zero'
        Write (messages(3),'(a)')        ' - 1exp_to_constant'
        Write (messages(4),'(a)')        ' - 2exp_to_zero'
        Write (messages(5),'(a)')        ' - 2exp_to_constant'
        Call info(messages, 5)
        Call error_stop(' ')      
      End If
    Else If (Trim(fit_data%what_to_fit%type) == 'tcf') Then
      If (Trim(fit_data%fitting_function%type) /= '1exp_to_zero'     .And. &
          Trim(fit_data%fitting_function%type) /= '1exp_to_constant' .And. &
          Trim(fit_data%fitting_function%type) /= '2exp_to_zero') Then
        Write (messages(1),'(2(1x,a))')  Trim(error_set), 'Invalid option for the "fitting_function" directive.'
        Call info(messages, 1)
        Write (messages(1),'(a)')        'Implemented options for TCF:'
        Write (messages(2),'(a)')        ' - 1exp_to_zero'
        Write (messages(3),'(a)')        ' - 1exp_to_constant'
        Write (messages(4),'(a)')        ' - 2exp_to_zero'
        Call info(messages, 4)
        Call error_stop(' ')      
      End If
    Else If (Trim(fit_data%what_to_fit%type) == 'spcf') Then
      If (Trim(fit_data%fitting_function%type) /= '1exp_to_zero'     .And. &
          Trim(fit_data%fitting_function%type) /= '1exp_to_constant' .And. &
          Trim(fit_data%fitting_function%type) /= '2exp_to_zero') Then
        Write (messages(1),'(2(1x,a))')  Trim(error_set), 'Invalid option for the "fitting_function" directive.'
        Call info(messages, 1)
        Write (messages(1),'(a)')        'Implemented options for SPCF:'
        Write (messages(2),'(a)')        ' - 1exp_to_zero'
        Write (messages(3),'(a)')        ' - 1exp_to_constant'
        Write (messages(4),'(a)')        ' - 2exp_to_zero'
        Call info(messages, 4)
        Call error_stop(' ')      
      End If      
    Else If (Trim(fit_data%what_to_fit%type) == 'msd') Then
      If (Trim(fit_data%fitting_function%type) /= 'linear') Then
        Write (messages(1),'(2(1x,a))')  Trim(error_set), 'Invalid option for the "fitting_function" directive.'
        Call info(messages, 1)
        Write (messages(1),'(a)')        'Implemented options for MSD:'
        Write (messages(2),'(a)')        ' - linear'
        Call info(messages, 2)
        Call error_stop(' ')      
      End If    
    End If        

    ! Define functional_form for fitting
    fit_data%functional_form=Trim(fit_data%what_to_fit%type)//'_'//Trim(fit_data%fitting_function%type)

    If (Trim(fit_data%functional_form) == 'ocf_1exp_to_zero'     .Or. &
        Trim(fit_data%functional_form) == 'ocf_1exp_to_constant' .Or. &
        Trim(fit_data%functional_form) == 'ocf_2exp_to_zero'     .Or. &
        Trim(fit_data%functional_form) == 'ocf_2exp_to_constant') Then
        If(fit_data%reactive_species%stat) Then
          fit_data%quantity='Orientational Correlation Function (OCF) of the reactive species "'&
                           &//Trim(fit_data%species_name%type)//'"'
        Else
          fit_data%quantity='Orientational Correlation Function (OCF) of the nonreactive species "'&
                           &//Trim(fit_data%species_name%type)//'"'
        End If
    End If

    If (Trim(fit_data%functional_form) == 'tcf_1exp_to_zero'     .Or. &
        Trim(fit_data%functional_form) == 'tcf_1exp_to_constant' .Or. &
        Trim(fit_data%functional_form) == 'tcf_2exp_to_zero') Then
        If (fit_data%reactive_species%fread) Then
           If (.Not.fit_data%reactive_species%stat) Then
             Write (messages(1),'(2(1x,a))') Trim(error_set), 'The user has requested to fit&
                                            & Transfer Correlation Function (TCF) data but the&
                                            & "reactive_species" directive was set to .FALSE.'
             Write (messages(2),'(1x,a)')   ' TCF implicitly refers to compute the transfer of reactive systems.&
                                            & Please resolve this inconsistency in the settings.'
             Call info(messages,2)
             Call error_stop(' ')
           Else
             freact=.True.
           End If
        Else
          freact=.True. 
        End If
      If (freact) Then
         fit_data%quantity='Transfer Correlation Function (TCF) of the reactive species "'&
                             &//Trim(fit_data%species_name%type)//'"'
         fit_data%reactive_species%stat=.True.                    
      End If
    End If

    If (Trim(fit_data%functional_form) == 'spcf_1exp_to_zero'     .Or. &
        Trim(fit_data%functional_form) == 'spcf_1exp_to_constant' .Or. &
        Trim(fit_data%functional_form) == 'spcf_2exp_to_zero') Then
        If (fit_data%reactive_species%fread) Then
           If (.Not.fit_data%reactive_species%stat) Then
             Write (messages(1),'(2(1x,a))') Trim(error_set), 'The user has requested to fit&
                                            & Special Pair Correlation Function (SPCF) data but the&
                                            & "reactive_species" directive was set to .FALSE.'
             Write (messages(2),'(1x,a)')   ' SPCF implicitly refers to compute the transfer of reactive systems.&
                                            & Please resolve this inconsistency in the settings.'
             Call info(messages,2)
             Call error_stop(' ')
           Else
             freact=.True.
           End If
        Else
          freact=.True. 
        End If
      If (freact) Then
         fit_data%quantity='Special Pair Correlation Function (SPCF) from the reactive species "'&
                             &//Trim(fit_data%species_name%type)//'"'
         fit_data%reactive_species%stat=.True.                    
      End If
    End If    
    
    If (Trim(fit_data%functional_form) == 'msd_linear') Then
        If (fit_data%reactive_species%fread) Then
           If (fit_data%reactive_species%stat) Then
             Write (messages(1),'(2(1x,a))') Trim(error_set), 'The user has requested to fit&
                                           & Mean Square Displacement (MSD) data but the&
                                           & "reactive_species" directive was set to .TRUE.'
             Write (messages(2),'(1x,a)')   ' MSD data implicitly refers to nonreactive species.&
                                           & Please resolve this inconsistency in the settings.'
             Call info(messages,2)
             Call error_stop(' ')
           Else
             freact=.False.
           End If
        Else
           freact=.False.
        End If
      If (.Not. freact) Then
        fit_data%quantity='Mean Square Displacement (MSD) of the nonreactive species "'//Trim(fit_data%species_name%type)//'"' 
        fit_data%reactive_species%stat=.False.
      End If
    End If

    ! Include superfast response in the fittings
    If (fit_data%fit_superfast_time%fread) Then
      If (fit_data%fit_superfast_time%fail) Then
        Write (messages(1),'(2(1x,a))') Trim(error_set), 'Missing (or wrong) specification for directive&
                                  & "fit_superfast_time" (choose either .True. or .False.)'
        Call info(messages,1)
        Call error_stop(' ')
      End If
      If (fit_data%fit_superfast_time%stat) Then
         If (Trim(fit_data%what_to_fit%type)/='ocf')Then
           Write (messages(1),'(2(1x,a))') Trim(error_set), 'Setting "fit_superfast_time" to .True. is only&
                                         & meaningful for the fitting of "ocf" data.'
           Call info(messages,1)
           Call error_stop(' ')         
         End If
      End If
    Else
       fit_data%fit_superfast_time%stat=.False.
    End If    
    
    If (fit_data%std_samples%fread) Then
      If (fit_data%std_samples%fail) Then
        Write (messages(1),'(2(1x,a))') Trim(error_set),&
                                & 'Wrong (or missing) settings for "std_samples" directive. It must be an odd, positive integer&
                                & with 1 the lowest possible value (i.e. 1, 3, 5, 7, etc.). If set to 1, only&
                                & the average values will be considered for fitting.'
        Call info(messages, 1)
        Call error_stop(' ')
      Else
        If (fit_data%std_samples%value < 1 .Or. (mod(fit_data%std_samples%value,2)==0)) Then
           Write (messages(1),'(1x,2a, i2,a)') Trim(error_set), &
                                &' The value for directive "std_samples" must be an odd, positive integer&
                                & with 1 the lowest possible value (i.e. 1, 3, 5, 7, etc.). If set to 1, only&
                                & the average values will be considered for fitting.'
          Call info(messages, 1)
          Call error_stop(' ')
        End If
          If (fit_data%std_samples%value > max_std_samples) Then
             Write(max_samples, '(i2)') max_std_samples
             Write (messages(1),'(1x,2a, i2,a)') Trim(error_set), ' For convenience, the value of "std_samples" must be lower&
                                               & than a pre-defined maximum of '//Trim(max_samples)
             Call info(messages, 1)
             Call error_stop(' ')
          End If
      End If    
    Else  
      fit_data%std_samples%value=1
    End If

    If (fit_data%max_iterations%fread) Then
      If (fit_data%max_iterations%fail) Then
        Write (messages(1),'(2(1x,a))') Trim(error_set),&
                                & 'Wrong (or missing) settings for "max_iterations" directive. It must be a positive integer.'
        Call info(messages, 1)
        Call error_stop(' ')
      Else
          If (fit_data%max_iterations%value < 0) Then
            Write (messages(1),'(1x,2a, i2,a)') Trim(error_set), &
                                    &' Input value for directive "max_iterations" must be a positive integer.'
            Call info(messages, 1)
            Call error_stop(' ')
          End If
          If (fit_data%max_iterations%value > max_iterations_fitting) Then
            Write(word,'(i4)') max_iterations_fitting   
            Write (messages(1),'(1x,2a, i2,a)') Trim(error_set), &
                                    &' Input value for directive "max_iterations" cannot exceed a predefined maximum limit of '&
                                    &//Trim(word)//' iterations'
            Call info(messages, 1)
            Call error_stop(' ')
          End If
       End If   
    Else  
      Write (messages(1),'(2(1x,a))')  Trim(error_set), 'The user must define the "max_iterations" directive, which&
                                    & indicates the maximum number of iterations for fitting.'
      Call info(messages, 1)
      Call error_stop(' ')
    End If

    ! max_times_std
    If (fit_data%max_times_std%fread) Then
      If (fit_data%max_times_std%fail) Then
         Write (messages(1),'(2(1x,a))') Trim(error_set),&
                                & 'Wrong (or missing) settings for "max_times_std" directive'
        Call info(messages, 1)
        Call error_stop(' ')
      Else
        If (fit_data%max_times_std%value<0.0_wp) Then
          Write (messages(1),'(2(1x,a))') Trim(error_set),&
                                & 'Value for directive "max_times_std" must be a positive'
          Call info(messages, 1)
          Call error_stop(' ')
        End If
      End If
    Else
      If (fit_data%std_samples%value /= 1) Then 
        Write (messages(1),'(2(1x,a))')  Trim(error_set), 'The user must define the "max_times_std" directive.&
                                    & Although it depends on the data (including STD) to be fitted, we recommend&
                                    & not to exceed the value of 1.'
        Call info(messages, 1)
        Call error_stop(' ')
      End If
    End If
    
    
    If (fit_data%plot_fittings%fread) Then
      If (fit_data%plot_fittings%fail) Then
        Write (messages(1),'(2(1x,a))') Trim(error_set), 'Wrong specification for directive&
                                  & "plot_fittings" (choose either .True. or .False.)'
        Call info(messages,1)
        Call error_stop(' ')
      End If
    Else
       fit_data%plot_fittings%stat=.True.
    End If
    
    If (fit_data%plot_raw_data%fread) Then
      If (fit_data%plot_raw_data%fail) Then
        Write (messages(1),'(2(1x,a))') Trim(error_set), 'Wrong specification for directive&
                                  & "plot_raw_data" (choose either .True. or .False.)'
        Call info(messages,1)
        Call error_stop(' ')
      End If
    Else
       fit_data%plot_raw_data%stat=.True.
    End If

 
    If (fit_data%time_column%fread) Then
      If (fit_data%time_column%fail) Then
        Write (messages(1),'(2(1x,a))') Trim(error_set),&
                                & 'Wrong (or missing) settings for "time_column_number" directive.'
        Call info(messages, 1)
        Call error_stop(' ')
      Else
        If (fit_data%time_column%value < 1) Then
           Write (messages(1),'(1x,2a, i2,a)') Trim(error_set), &
                                  &' The value for directive "time_column_number" must be a positive integer'
          Call info(messages, 1)
          Call error_stop(' ')
        End If
      End If    
    Else  
      Write (messages(1),'(2(1x,a))')  Trim(error_set), 'The user must define the "time_column_number" directive, which refers&
                                     & to the column number of the data file corresponding to the recorded time.'
      Call info(messages, 1)
      Call error_stop(' ')
    End If
    
    If (fit_data%avg_column%fread) Then
      If (fit_data%avg_column%fail) Then
        Write (messages(1),'(2(1x,a))') Trim(error_set),&
                                & 'Wrong (or missing) settings for "avg_column_number" directive.'
        Call info(messages, 1)
        Call error_stop(' ')
      Else
        If (fit_data%avg_column%value < 1) Then
           Write (messages(1),'(1x,2a, i2,a)') Trim(error_set), &
                                  &' The value for directive "avg_column_number" must be a positive integer'
          Call info(messages, 1)
          Call error_stop(' ')
        End If
      End If    
    Else  
      Write (messages(1),'(2(1x,a))')  Trim(error_set), 'The user must define the "avg_column_number" directive, which refers&
                                     & to the column number of the data file corresponding to the average values.'
      Call info(messages, 1)
      Call error_stop(' ')
    End If

    If (fit_data%std_column%fread) Then
      If (fit_data%std_column%fail) Then
        Write (messages(1),'(2(1x,a))') Trim(error_set),&
                                & 'Wrong settings for "std_column_number" directive.'
        Call info(messages, 1)
        Call error_stop(' ')
      Else
        If (fit_data%std_column%value < 1) Then
           Write (messages(1),'(1x,2a, i2,a)') Trim(error_set), &
                                  &' The value for directive "std_column_number" must be a positive integer'
          Call info(messages, 1)
          Call error_stop(' ')
        End If
      End If    
    Else  
      If (fit_data%std_samples%value /= 1) Then
        Write (messages(1),'(2(1x,a))')  Trim(error_set), 'The user must define the "std_column_number" directive, which refers&
                                       & to the column number of the data file corresponding to the STD values.'
        Call info(messages, 1)
        Call error_stop(' ')
      End If
    End If

    If (fit_data%time_column%value == fit_data%avg_column%value) Then
      Write (messages(1),'(2(1x,a))')  Trim(error_set), 'The values assigned to "time_column_number" and "avg_column_number"&
                                     & must be different!'
      Call info(messages, 1)
      Call error_stop(' ')
    End If
    
    If (fit_data%std_column%fread .And. fit_data%std_samples%value /= 1) Then
      If (fit_data%time_column%value == fit_data%std_column%value) Then
        Write (messages(1),'(2(1x,a))')  Trim(error_set), 'The values assigned to "time_column_number" and "std_column_number"&
                                       & must be different!'
        Call info(messages, 1)
        Call error_stop(' ')
      End If
      If (fit_data%avg_column%value == fit_data%std_column%value) Then
        Write (messages(1),'(2(1x,a))')  Trim(error_set), 'The values assigned to "avg_column_number" and "std_column_number"&
                                       & must be different!'
        Call info(messages, 1)
        Call error_stop(' ')
      End If      
    End If
    
    ! Check start_time
    Call check_time_directive(fit_data%start_time, 'start_time', error_set, .False.)
    ! Check end_time
    Call check_time_directive(fit_data%end_time, 'end_time', error_set, .False.)
    
    If (fit_data%end_time%fread) Then
       If((fit_data%end_time%value < fit_data%start_time%value) .Or. &
           Abs(fit_data%end_time%value-fit_data%start_time%value) < epsilon(1.0_wp) ) Then
           Write (messages(1),'(2(1x,a))') Trim(error_set),&
                                & 'the value for "start_time" must be lower than "end_time".&
                                & Make sure the units ensure this condition.'
          Call info(messages, 1)
          Call error_stop(' ')
       End If    
    End If

    If (Trim(fit_data%functional_form) == 'msd_linear') Then
      If (.Not. fit_data%msd%invoke%fread) Then
           Write (messages(1),'(2(1x,a))') Trim(error_set), 'The user has requested to fit the&
                                          & Mean Square Displacement (MSD) but the &msd block&
                                          & is not defined. Please add the block'
           Call info(messages,1)
           Call error_stop(' ')
      Else
        Call check_msd(files, fit_data)
        If (fit_data%end_time%fread) Then
          If(Abs(fit_data%start_time%value) < epsilon(1.0_wp)) Then
            Write (messages(1),'(2(1x,a))') Trim(error_set),&
                                & 'the value for "start_time" for MSD fitting must be larger than zero'
            Call info(messages, 1)
            Call error_stop(' ')
          End If
        Else
          Write (messages(1),'(2(1x,a))') Trim(error_set),&
                                & 'the user must define the "start_time" for MSD fitting, larger than zero'
          Call info(messages, 1)
          Call error_stop(' ')
        End If    
      End If
    End If
    
    If (fit_data%std_column%fread .Or.  fit_data%max_times_std%fread) Then
      If (.Not. fit_data%std_samples%fread) Then
        Write (messages(1),'(2(1x,a))') Trim(error_set),&
                                & 'To include the standard deviations in the fitting the user must&
                                & define the "std_samples" directive.'
        Call info(messages, 1)
        Call error_stop(' ')
      End If  
    End If    
    
    If (fit_data%print_average_ssr%fread) Then
      If (fit_data%print_average_ssr%fail) Then
        Write (messages(1),'(2(1x,a))') Trim(error_set), 'Missing (or wrong) specification for directive&
                                  & "print_average_ssr" (choose either .True. or .False.)'
        Call info(messages,1)
        Call error_stop(' ')
      End If
    Else
       fit_data%print_average_ssr%stat=.False.
    End If    

    ! Input parameters
    If (fit_data%input_parameters%fread) Then
      Call check_input_parameters(files, fit_data)
    Else      
      Call set_initial_parameters(fit_data)
    End If
    
  End Subroutine check_fitting_settings

  Subroutine check_time_directive(T, tag, error_set, kill)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to check time related directives
    !
    ! author    - i.scivetti March 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(in_param),          Intent(InOut)  :: T
    Character(Len=*),        Intent(In   )  :: tag
    Character(Len=*),        Intent(In   )  :: error_set
    Logical,                 Intent(In   )  :: kill

    Character(Len=256)  :: messages(2)
    
    If (T%fread) Then
      If (T%fail) Then
        Write (messages(1),'(2(1x,a))') Trim(error_set), 'Wrong (or missing) settings for the "'&
                                      &//Trim(T%tag)//'" directive.'
        Call info(messages, 1)
        Call error_stop(' ')
      Else
        If (T%value < 0.0_wp) Then
          Write (messages(1),'(2(1x,a))') Trim(error_set), &
                                    &'Input value for "'//Trim(T%tag)//&
                                    &'" MUST NOT be negative'
          Call info(messages, 1)
          Call error_stop(' ')
        End If
        Call capital_to_lower_case(T%units)
        If (Trim(T%units) /= 'fs' .And. &
           Trim(T%units) /= 'ps') Then
           Write (messages(1),'(2(1x,a))')  Trim(error_set),&
                                    & 'Units for directive "'//Trim(T%tag)//&
                                    &'" must be "fs" or "ps". Have you included the units?'
          Call info(messages, 1)
          Call error_stop(' ')
        End If
        ! Transform to fs
        If (Trim(T%units) == 'fs') Then
           T%value=T%value/1000.0_wp
        End If
      End If
    Else 
      If (kill)then
        Write (messages(1),'(2(1x,a))')  Trim(error_set), 'The user must define the "'//Trim(tag)//'" directive'
        Call info(messages, 1)
        Call error_stop(' ')
      End If
    End If
    
  End Subroutine check_time_directive
  
  Subroutine check_msd(files, fit_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to check the settings of the &MSD block
    !
    ! author    - i.scivetti March 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),    Intent(In   ) :: files(:)
    Type(fit_type),    Intent(InOut) :: fit_data

    Character(Len=256)  :: messages(2)
    Character(Len=64 )  :: error_set

    error_set = '***ERROR in the &MSD block of file '//Trim(files(FILE_SET)%filename)//' -'

    If (fit_data%msd%select%fread) Then
      If (fit_data%msd%select%fail) Then
        Write (messages(1),'(2(1x,a))') Trim(error_set), 'Wrong (or missing) settings for the "select" directive.'
        Call info(messages, 1)
        Call error_stop(' ')
      Else
        If (Trim(fit_data%msd%select%type)/='x'  .And. &
            Trim(fit_data%msd%select%type)/='y'  .And. &
            Trim(fit_data%msd%select%type)/='z'  .And. &
            Trim(fit_data%msd%select%type)/='xy' .And. &
            Trim(fit_data%msd%select%type)/='xz' .And. &
            Trim(fit_data%msd%select%type)/='yz' .And. &
            Trim(fit_data%msd%select%type)/='xyz') Then
             Write (messages(1),'(2(1x,a))') Trim(error_set), &
                                    &'Wrong input for "select". Valid options: "x", "y", "z", "xy",&
                                    & "xz", "yz" or "xyz"'
          Call info(messages, 1)
          Call error_stop(' ')
        End If
      End If
    Else
       Write (messages(1),'(2(1x,a))')  Trim(error_set), 'The user must define the "select" directive'
       Write (messages(2),'( (1x,a))') 'Valid options: "x", "y", "z", "xy", "xz", "yz" or "xyz"'
       Call info(messages, 2)
       Call error_stop(' ')
    End If
    
    If (Trim(fit_data%msd%select%type)=='x'  .Or. &
        Trim(fit_data%msd%select%type)=='y'  .Or. &
        Trim(fit_data%msd%select%type)=='z') Then
      fit_data%msd%dim=1
    Else If (Trim(fit_data%msd%select%type)=='xy' .Or. &
             Trim(fit_data%msd%select%type)=='xz' .Or. &
             Trim(fit_data%msd%select%type)=='yz') Then
      fit_data%msd%dim=2    
    Else If (Trim(fit_data%msd%select%type)=='xyz') Then
      fit_data%msd%dim=3
    End If

    If (fit_data%msd%units%fread) Then
      If (fit_data%msd%units%fail) Then
        Write (messages(1),'(2(1x,a))') Trim(error_set), 'Wrong (or missing) settings for the "units" directive.'
        Call info(messages, 1)
        Call error_stop(' ')
      Else
        If (Trim(fit_data%msd%units%type)/='angstrom^2' .And. Trim(fit_data%msd%units%type)/='bohr^2') Then
             Write (messages(1),'(2(1x,a))') Trim(error_set), &
                                    &'Wrong input for "units". Valid options: "Angstrom^2" and "Bohr^2"'
          Call info(messages, 1)
          Call error_stop(' ')
        Else
          If (Trim(fit_data%msd%units%type)=='angstrom^2') Then
            fit_data%msd%factor= 1.0_wp
          Else If (Trim(fit_data%msd%units%type)=='bohr^2') Then
            fit_data%msd%factor= Bohr_to_A**2
          End If
        End If
      End If
    Else
       Write (messages(1),'(2(1x,a))')  Trim(error_set), 'The user must define the "units" directive'
       Write (messages(2),'( (1x,a))') 'Valid options: "Angstrom^2" and "Bohr^2"'
       Call info(messages, 2)
       Call error_stop(' ')
    End If
    
  End Subroutine check_msd
 
  Subroutine check_input_parameters(files, fit_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to check settings of the &input_parameters block
    !
    ! author    - i.scivetti March 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),    Intent(In   ) :: files(:)
    Type(fit_type),     Intent(InOut) :: fit_data

    If (Trim(fit_data%functional_form) == 'ocf_2exp_to_constant') Then
      Call check_param_if_defined(files,fit_data%param%A1%ini,'A1')
      Call check_param_if_defined(files,fit_data%param%A2%ini,'A2')
      Call check_param_if_defined(files,fit_data%param%A3%ini,'A3')
      Call check_param_if_defined(files,fit_data%param%T2%ini,'T2')
      Call check_param_if_defined(files,fit_data%param%T3%ini,'T3')
    Else If (Trim(fit_data%functional_form) == 'ocf_2exp_to_zero') Then
      Call check_param_if_defined(files,fit_data%param%A1%ini,'A1')
      Call check_param_if_defined(files,fit_data%param%A2%ini,'A2')
      Call check_param_if_defined(files,fit_data%param%T2%ini,'T2')
      Call check_param_if_defined(files,fit_data%param%T3%ini,'T3')
    Else If (Trim(fit_data%functional_form) == 'ocf_1exp_to_constant') Then
      Call check_param_if_defined(files,fit_data%param%A1%ini,'A1')
      Call check_param_if_defined(files,fit_data%param%A2%ini,'A2')
      Call check_param_if_defined(files,fit_data%param%T2%ini,'T2')
    Else If (Trim(fit_data%functional_form) == 'ocf_1exp_to_zero') Then    
      Call check_param_if_defined(files,fit_data%param%A1%ini,'A1')
      Call check_param_if_defined(files,fit_data%param%T2%ini,'T2')
    Else If (Trim(fit_data%functional_form) == 'tcf_2exp_to_zero') Then
      Call check_param_if_defined(files,fit_data%param%A0%ini,'A0')
      Call check_param_if_defined(files,fit_data%param%A1%ini,'A1')
      Call check_param_if_defined(files,fit_data%param%A2%ini,'A2')
    Else If (Trim(fit_data%functional_form) == 'tcf_1exp_to_constant') Then
      Call check_param_if_defined(files,fit_data%param%A0%ini,'A0')
      Call check_param_if_defined(files,fit_data%param%A1%ini,'A1')
    Else If (Trim(fit_data%functional_form) == 'tcf_1exp_to_zero') Then
      Call check_param_if_defined(files,fit_data%param%A1%ini,'A1')
    Else If (Trim(fit_data%functional_form) == 'spcf_2exp_to_zero') Then
      Call check_param_if_defined(files,fit_data%param%A0%ini,'A0')
      Call check_param_if_defined(files,fit_data%param%A1%ini,'A1')
      Call check_param_if_defined(files,fit_data%param%A2%ini,'A2')
    Else If (Trim(fit_data%functional_form) == 'spcf_1exp_to_constant') Then
      Call check_param_if_defined(files,fit_data%param%A0%ini,'A0')
      Call check_param_if_defined(files,fit_data%param%A1%ini,'A1')
    Else If (Trim(fit_data%functional_form) == 'spcf_1exp_to_zero') Then
      Call check_param_if_defined(files,fit_data%param%A1%ini,'A1')      
    Else If (Trim(fit_data%functional_form) == 'msd_linear') Then
      Call check_param_if_defined(files,fit_data%param%A0%ini,'A0')
      Call check_param_if_defined(files,fit_data%param%A1%ini,'A1')
    End If

    If (Trim(fit_data%what_to_fit%type)=='ocf' .And. fit_data%fit_superfast_time%stat) Then
      Call check_param_if_defined(files,fit_data%param%T1%ini,'T1')
    End If    
  
  End Subroutine check_input_parameters

  Subroutine check_param_if_defined(files,T,tag)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to check time related directives
    !
    ! author    - i.scivetti March 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),  Intent(In   ) :: files(:)
    Type(in_param),   Intent(InOut) :: T
    Character(Len=*), Intent(In   ) :: tag

    Character(Len=256)  :: message
    Character(Len=64 )  :: error_set

    error_set = '***ERROR in the "&input_parameters" block (file '//Trim(files(FILE_SET)%filename)//'):'
    
    If (T%fread) Then
      If (T%fail) Then
         Write (message,'(2(1x,a))') Trim(error_set),&
                                & 'Wrong (or missing) settings for parameter "'//Trim(T%tag)//'"'
        Call info(message, 1)
        Call error_stop(' ')
      Else
        If (Abs(T%value)<epsilon(1.0_wp)) Then
          Write (message,'(2(1x,a))') Trim(error_set),&
                                & 'Input value for parameter "'//Trim(T%tag)//'" must be non-zero'
          Call info(message, 1)
          Call error_stop(' ')
        End If
      End If
    Else
      Write (message,'(1x,a)') Trim(error_set)
      Call info(message, 1)
      Write (message,'(1x,a)')  'Parameter "'//Trim(tag)//'" must be defined&
                               & for the selected fitting options".'
      Call info(message, 1)
      Call error_stop(' ')
    End If
    
  End Subroutine check_param_if_defined

  Subroutine set_initial_parameters(fit_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to define the initial value for the parameters
    !
    ! author    - i.scivetti March 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(fit_type),     Intent(InOut) :: fit_data

    If (Trim(fit_data%functional_form) == 'ocf_1exp_to_constant') Then
      If (fit_data%reactive_species%stat) Then
        fit_data%param%A2%ini%value=0.01_wp
        fit_data%param%A1%ini%value=0.10_wp
        fit_data%param%T2%ini%value=3.00_wp
      Else
        fit_data%param%A2%ini%value=0.10
        fit_data%param%A1%ini%value=0.70
        fit_data%param%T2%ini%value=11.0
      End If
    Else If (Trim(fit_data%functional_form) == 'ocf_1exp_to_zero') Then
      If (fit_data%reactive_species%stat) Then
        fit_data%param%A1%ini%value=0.1_wp
        fit_data%param%T2%ini%value=3.0_wp
      Else
        fit_data%param%A1%ini%value=0.7_wp
        fit_data%param%T2%ini%value=11.0_wp
      End If
    Else If (Trim(fit_data%functional_form) == 'ocf_2exp_to_constant') Then
      If (fit_data%reactive_species%stat) Then
        fit_data%param%A3%ini%value=0.01_wp
        fit_data%param%A2%ini%value=0.05_wp
        fit_data%param%A1%ini%value=0.10_wp
        fit_data%param%T2%ini%value=0.70_wp
        fit_data%param%T3%ini%value=4.00_wp        
      Else
        fit_data%param%A3%ini%value=0.05_wp
        fit_data%param%A2%ini%value=0.10_wp
        fit_data%param%A1%ini%value=0.70_wp
        fit_data%param%T2%ini%value=1.00_wp
        fit_data%param%T3%ini%value=10.0_wp        
      End If
    Else If (Trim(fit_data%functional_form) == 'ocf_2exp_to_zero') Then
        If (fit_data%reactive_species%stat) Then
          fit_data%param%A2%ini%value=0.01_wp
          fit_data%param%A1%ini%value=0.10_wp
          fit_data%param%T2%ini%value=1.00_wp
          fit_data%param%T3%ini%value=4.00_wp        
        Else
          fit_data%param%A2%ini%value=0.10_wp
          fit_data%param%A1%ini%value=0.70_wp
          fit_data%param%T2%ini%value=1.00_wp
          fit_data%param%T3%ini%value=10.0_wp        
        End If
    Else If (Trim(fit_data%functional_form) == 'tcf_2exp_to_zero') Then
          fit_data%param%A0%ini%value=0.50_wp
          fit_data%param%A1%ini%value=0.10_wp
          fit_data%param%A2%ini%value=0.01_wp
    Else If (Trim(fit_data%functional_form) == 'tcf_1exp_to_constant') Then
          fit_data%param%A0%ini%value=0.50_wp
          fit_data%param%A1%ini%value=0.10_wp
    Else If (Trim(fit_data%functional_form) == 'tcf_1exp_to_zero') Then
          fit_data%param%A1%ini%value=0.10_wp
    Else If (Trim(fit_data%functional_form) == 'spcf_2exp_to_zero') Then
          fit_data%param%A0%ini%value=0.50_wp
          fit_data%param%A1%ini%value=0.10_wp
          fit_data%param%A2%ini%value=0.01_wp
    Else If (Trim(fit_data%functional_form) == 'spcf_1exp_to_constant') Then
          fit_data%param%A0%ini%value=0.50_wp
          fit_data%param%A1%ini%value=0.10_wp
    Else If (Trim(fit_data%functional_form) == 'spcf_1exp_to_zero') Then
          fit_data%param%A1%ini%value=0.10_wp          
    Else If (Trim(fit_data%functional_form) == 'msd_linear') Then
          fit_data%param%A0%ini%value=0.2_wp
          fit_data%param%A1%ini%value=1.0_wp
    End If
    
    If (Trim(fit_data%what_to_fit%type)=='ocf' .And. fit_data%fit_superfast_time%stat) Then
     fit_data%param%T1%ini%value=fit_data%param%T2%ini%value/50.0_wp
    End If
    
  End Subroutine set_initial_parameters

  Subroutine check_validity_parameters(A0, A1, A2, A3, T1, T2, T3, functional_form, fok, fit_libration, error_msg)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to check the validity of the involved parameters 
    !
    ! author    - i.scivetti March 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Real(kind=wp),      Intent(In   ) :: A0
    Real(kind=wp),      Intent(In   ) :: A1
    Real(kind=wp),      Intent(In   ) :: A2
    Real(kind=wp),      Intent(In   ) :: A3
    Real(kind=wp),      Intent(In   ) :: T1
    Real(kind=wp),      Intent(In   ) :: T2
    Real(kind=wp),      Intent(In   ) :: T3
    Character(Len=256), Intent(In   ) :: functional_form
    Logical,            Intent(InOut) :: fok
    Logical,            Intent(In   ) :: fit_libration
    Character(Len=256), Optional,  Intent(  Out) :: error_msg

    Character(Len=256) :: head
    
    head='*** PROBLEMS with the "&input_parameters" block:'
    
    If (Trim(functional_form) == 'ocf_2exp_to_constant' .Or. &
       Trim(functional_form) == 'ocf_2exp_to_zero') Then
       If (T2 >= T3) Then
         fok=.False.
         If (Present(error_msg))then
           Write(error_msg, '(a)') Trim(head)//' T2 must be lower than T3!'
         End If
       End If
       If (T3 < 0.0_wp .Or. T2 < 0.0_wp) Then
         fok=.False.
         If (Present(error_msg))then
           Write(error_msg, '(a)') Trim(head)//' T2 and T3 must be larger than zero!'
         End If  
       End If
       If (A1 > 1.0_wp .Or. A1 < -1.0_wp) Then
         fok=.False.
         If (Present(error_msg))then
           Write(error_msg, '(a)') Trim(head)//' A1 must lie between -1.0 and 1.0!'
         End If  
       End If
       If (A2 > A1) Then
         fok=.False.
         If (Present(error_msg))then
           Write(error_msg, '(a)') Trim(head)//' A2 must be lower than A1!'
         End If 
       End If
       If (A2 > 1.0_wp .Or. A2 < -1.0_wp) Then
         fok=.False.
         If (Present(error_msg))then
           Write(error_msg, '(a)') Trim(head)//' A2 must lie between -1.0 and 1.0!'
         End If  
       End If
       
       If (Trim(functional_form) == 'ocf_2exp_to_constant') Then
          If (A3 > 1.0_wp .Or. A3 < -1.0_wp) Then
            fok=.False.
            If (Present(error_msg))then
              Write(error_msg, '(a)') Trim(head)//' A3 must lie between -1.0 and 1.0!'
            End If  
          End If
          If (A3 > A2) Then
           fok=.False.
           If (Present(error_msg))then
             Write(error_msg, '(a)') Trim(head)//' A3 must be lower than A2!'
           End If  
          End If
       End If
    End If   

    If (Trim(functional_form) == 'ocf_1exp_to_constant' .Or. &  
       Trim(functional_form) == 'ocf_1exp_to_zero') Then
       If (T2 < 0.0_wp) Then
         fok=.False.
         If (Present(error_msg))then
           Write(error_msg, '(a)') Trim(head)//' T2 must be larger than zero!'
         End If
       End If
       If (A1 > 1.0_wp .Or. A1 < -1.0_wp) Then
         fok=.False.
         If (Present(error_msg))then
           Write(error_msg, '(a)') Trim(head)//' A1 must lie between -1.0 and 1.0!'
         End If
       End If
       If (Trim(functional_form) == 'ocf_1exp_to_constant') Then
          If (A2 > 1.0_wp .Or. A2 < -1.0_wp) Then
            fok=.False.
            If (Present(error_msg))then
              Write(error_msg, '(a)') Trim(head)//' A2 must lie between -1.0 and 1.0!'
            End If
          End If
          If (A2 > A1) Then
            fok=.False.
            If (Present(error_msg))then
              Write(error_msg, '(a)') Trim(head)//' A2 must be lower than A1!'
            End If
          End If
       End If
    End If
    
    If (fit_libration) Then
      If (T1 >= T2) Then
        fok=.False.
        If (Present(error_msg))then
          Write(error_msg, '(a)') Trim(head)//' T1 must be lower than T2!'
        End If
      End If
      If (T1 < 0.0_wp) Then
        fok=.False.
        If (Present(error_msg))then
          Write(error_msg, '(a)') Trim(head)//' T1 must be larger than zero!'
        End If  
      End If    
    End If
    

    If (Trim(functional_form) == 'tcf_1exp_to_constant' .Or. &  
        Trim(functional_form) == 'tcf_1exp_to_zero') Then
       If (A1 < 0.0_wp) Then
         fok=.False.
         If (Present(error_msg))then
           Write(error_msg, '(a)') Trim(head)//' A1 must be larger than zero!'
         End If  
       End If
       If (Trim(functional_form) == 'tcf_1exp_to_constant') Then
          If (A0 > 1.0_wp .Or. A0 < 0.0_wp) Then
            fok=.False.
            If (Present(error_msg))then
              Write(error_msg, '(a)') Trim(head)//' A0 must positive and lower or equal to 1.0!'
            End If  
          End If
       End If
    End If

    If (Trim(functional_form) == 'tcf_2exp_to_zero') Then
       If (A1 < 0.0_wp .Or. A2 < 0.0_wp) Then
         fok=.False.
         If (Present(error_msg))then
           Write(error_msg, '(a)') Trim(head)//' A1 and A2 must be larger than zero!'
         End If
       End If
       If (A1 <= A2) Then
         fok=.False.
         If (Present(error_msg))then
           Write(error_msg, '(a)') Trim(head)//' A1 must be larger than A2!'
         End If  
       End If
       If (A0 > 1.0_wp .Or. A0 < 0.0_wp) Then
         fok=.False.
         If (Present(error_msg))then
           Write(error_msg, '(a)') Trim(head)//' A1 must positive and lower or equal to 1.0!'
         End If  
       End If
    End If

    If (Trim(functional_form) == 'spcf_1exp_to_constant' .Or. &  
        Trim(functional_form) == 'spcf_1exp_to_zero') Then
       If (A1 < 0.0_wp) Then
         fok=.False.
         If (Present(error_msg))then
           Write(error_msg, '(a)') Trim(head)//' A1 must be larger than zero!'
         End If  
       End If
       If (Trim(functional_form) == 'spcf_1exp_to_constant') Then
          If (A0 > 1.0_wp .Or. A0 < 0.0_wp) Then
            fok=.False.
            If (Present(error_msg))then
              Write(error_msg, '(a)') Trim(head)//' A0 must positive and lower or equal to 1.0!'
            End If  
          End If
       End If
    End If

    If (Trim(functional_form) == 'spcf_2exp_to_zero') Then
       If (A1 < 0.0_wp .Or. A2 < 0.0_wp) Then
         fok=.False.
         If (Present(error_msg))then
           Write(error_msg, '(a)') Trim(head)//' A1 and A2 must be larger than zero!'
         End If
       End If
       If (A1 <= A2) Then
         fok=.False.
         If (Present(error_msg))then
           Write(error_msg, '(a)') Trim(head)//' A1 must be larger than A2!'
         End If  
       End If
       If (A0 > 1.0_wp .Or. A0 < 0.0_wp) Then
         fok=.False.
         If (Present(error_msg))then
           Write(error_msg, '(a)') Trim(head)//' A0 must positive and lower or equal to 1.0!'
         End If  
       End If
    End If    
    
    If (Trim(functional_form) == 'msd_linear') Then
       If (A1 < 0.0_wp .Or. Abs(A1) < epsilon(1.0_wp)) Then
         fok=.False.
         If (Present(error_msg))then
           Write(error_msg, '(a)') Trim(head)//' A1 must be larger than zero!'
         End If  
       End If
    End If
    
  End Subroutine check_validity_parameters

  Subroutine check_times_involved(fit_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Evaluate time domain and compare with time settings (if any) 
    ! Warn the user of problems 
    !
    ! author    - i.scivetti March 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(fit_type),     Intent(InOut) :: fit_data
    
    Character(Len=256) :: message

  
    ! Check for Start_time
    If(fit_data%start_time%fread) Then
      If(fit_data%start_time%value < fit_data%input_data(1)%time) Then
        Write(message, '(1x,a)') '**** WARNING: the value set for "start_time" is lower than&
                             & the minimum time recorded. The fitting will consider such minimum value&
                             & as the lowest bound.'
        Call info(message, 1)                     
      End If
      If (Trim(fit_data%what_to_fit%type)=='ocf') Then      
        If(fit_data%fit_superfast_time%stat) Then
          If(fit_data%start_time%value > fit_data%input_data(1)%time) Then
            Write(message, '(1x,a)') '**** WARNING: the user is advised to remove/comment the&
                                  & "start_time" directive to fit superfast times.'
            Call info(message, 1)
          End If
          If(fit_data%start_time%value > exp_ocf_lower_bound) Then
            Write(message, '(1x,a)') '**** WARNING: the value of "start_time" is not adequate to&
                                  & fit the superfast time reliably. Please check.'
            Call info(message, 1)
          End If
        Else
          If(fit_data%start_time%value < exp_ocf_lower_bound   .And. &
           & fit_data%input_data(1)%time < exp_ocf_lower_bound) Then
             Write(message, '(1x,a)') '**** WARNING: the user is advised to analyse&
                                & the input data at short times. This might require setting&
                                & "fit_superfast_time" to .True. for the current value of "start_time".'
            Call info(message, 1)
          End If 
        End If
      End If
    Else
      If(fit_data%input_data(1)%time < 0.0_wp) Then
        Write(message, '(1x,a)') '**** WARNING: the minimum time recorded is negative! The user&
                            & is advised to set the "start_time" directive or remove&
                            & negative times from the input data.'
        Call info(message, 1)      
      End If
    End If

    ! Check for End_time
    If(fit_data%end_time%fread) Then
      If(fit_data%end_time%value > fit_data%input_data(fit_data%num_input_data)%time) Then
        Write(message, '(1x,a)') '**** WARNING: the value set for "end_time" is larger than&
                            & the maximum time recorded.'
        Call info(message, 1)                     
      End If
    End If    

    ! General check
    If (Trim(fit_data%what_to_fit%type)=='ocf') Then      
      If(fit_data%fit_superfast_time%stat) Then
        If(fit_data%input_data(1)%time > exp_ocf_lower_bound) Then
          Write(message, '(1x,a)') '**** WARNING: the lowest recorded time is not adequate to fit&
                               & superfast times.'
          Call info(message, 1)
        End If
      Else
        If(fit_data%input_data(1)%time < exp_ocf_lower_bound) Then
          Write(message, '(1x,a)') '**** WARNING: the value for the lowest recorded time suggests&
                               & setting "fit_superfast_time" to .True.'
          Call info(message, 1)
        End If          
      End If
    End If    
       
  End Subroutine check_times_involved   

  Subroutine check_validity_mean_values(fit_data, fvalid)  
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
    ! Subroutine to check the physical validity of the mean values End module fitting_checks
    !
    ! author    - i.scivetti March 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(fit_type),     Intent(In   ) :: fit_data
    Logical,            Intent(  Out) :: fvalid

    fvalid=.True.

    If (Trim(fit_data%what_to_fit%type) == 'ocf') Then
      If (fit_data%param%A1%mean > 1.0_wp) Then
        fvalid=.False.
      End If
    
      If (Trim(fit_data%functional_form) == 'ocf_2exp_to_constant') Then
         If (fit_data%param%A3%mean < -1.0_wp) Then
           fvalid=.False.
         End If
      End If
      
      If (Trim(fit_data%functional_form) == 'ocf_1exp_to_constant') Then
         If (fit_data%param%A2%mean < -1.0_wp) Then
           fvalid=.False.
         End If
      End If
    End If    
    
  End Subroutine check_validity_mean_values
  
End module fitting_checks  
