!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Module to extract errors from the fitting of correlations
!
! Copyright - 2026 Ada Lovelace Centre (ALC)
!             Scientific Computing Department (SCD)
!             The Science and Technology Facilities Council (STFC)
!
! Author:     i.scivetti  March 2026
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Module fitting 

  Use fileset,        Only : file_type, &
                             FILE_SET, &
                             FILE_FIT_PARAMS,&
                             FILE_DUMP, &
                             FILE_INPUT_DATA, &
                             FILE_RAW_DATA, &
                             FILE_FITTED_DATA, &
                             refresh_out
                      
  Use gnuplot_setup,  Only : delete_gnuplot_working_files,&
                             generate_gnuplot_input,&
                             extract_data_from_file,&
                             define_gnuplot_actions
                             
  Use numprec,        Only : wi,&
                             wp 
                      
  Use process_data,   Only : remove_symbols
                            
  Use fitting_setup,  Only : fit_type, &
                             print_parameter
  
  Use fitting_checks, Only : check_validity_parameters,&
                             check_times_involved, &
                             check_validity_mean_values
   
  Use unit_output,    Only : error_stop, &
                             info

  Implicit None
  Private

  Public :: obtain_parameters_from_fitting

Contains

  Subroutine obtain_parameters_from_fitting(files, fit_data)  
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to compute the error associated with fittings
    !
    ! author    - i.scivetti March 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),    Intent(InOut) :: files(:)
    Type(fit_type),     Intent(InOut) :: fit_data

    Call fit_data%alloc_parameters() 
    
    Call set_fitting_directives(files, fit_data)
    Call extract_data_from_file(files, fit_data)    
    Call print_fitting_settings(files, fit_data)
    Call check_validity_input_data(fit_data)
    Call execute_fitting(files, fit_data)
    Call print_fitting_results(fit_data)
    Call refresh_out(files)
    Call display_rawdata_and_fittings(files, fit_data)
    
  End Subroutine obtain_parameters_from_fitting

  Subroutine execute_fitting(files, fit_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Execute fitting in three steps
    ! 1) fit average data and variations w.r.t these averages considering the STD, if reported
    ! 2) compute average of fitted parameters
    ! 3) compute the errors (STD) of the fitted parameters
    !
    ! author    - i.scivetti March 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),    Intent(InOut) :: files(:)
    Type(fit_type),     Intent(InOut) :: fit_data 
 
    Call fit_vs_std(files, fit_data)
    Call compute_average_over_stds(fit_data)
    Call compute_std_over_stds(fit_data)
    
  End Subroutine execute_fitting
  
  Subroutine display_rawdata_and_fittings(files, fit_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Obtain raw data to be fitted. Print the raw data and fittigs in colums. 
    !
    ! author    - i.scivetti March 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),    Intent(InOut) :: files(:)
    Type(fit_type),     Intent(InOut) :: fit_data

    If (fit_data%plot_raw_data%stat) Then
      Call plot_raw_data(files, fit_data)
    End If

    If (fit_data%plot_fittings%stat) Then
      Call plot_fittings(files, fit_data)
    End If
  
  End Subroutine display_rawdata_and_fittings

  Subroutine plot_fittings(files, fit_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Plots only those fittings that are physically valid (and have converged)
    !
    ! author    - i.scivetti March 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),    Intent(InOut) :: files(:)
    Type(fit_type),     Intent(InOut) :: fit_data
  
    Integer(Kind=wi)   :: iunit, i, j, neff
    Character(Len=256) :: xstd(fit_data%std_samples%value)
    Integer(Kind=wi)   :: indx(fit_data%std_samples%value)
  
    neff=0
    Do i=1, fit_data%std_samples%value
      If (fit_data%param%valid(i) .And. fit_data%param%fit_ok(i)) Then
        neff=neff+1
        indx(neff)=i
        fit_data%xstd=(i-1)*fit_data%delta_std-fit_data%max_times_std%value
        Write(xstd(indx(neff)),'(f12.3)') fit_data%xstd
        If (fit_data%xstd < 0.0_wp) Then
          xstd(indx(neff))='AVG'//Trim(Adjustl(xstd(neff)))//'*STD'
        Else If (Abs(fit_data%xstd) < epsilon(1.0_wp)) Then
          xstd(indx(neff))='AVG+'//Trim(Adjustl(xstd(neff)))//'*STD'
        Else
          xstd(indx(neff))='AVG+'//Trim(Adjustl(xstd(neff)))//'*STD'
        End If
        Call obtain_fitted_values(fit_data, indx(neff)) 
      End If
    End Do
  
    If (neff /= 0) Then
      Open(Newunit=files(FILE_FITTED_DATA)%unit_no, File=files(FILE_FITTED_DATA)%filename, Status='Replace')
      iunit=files(FILE_FITTED_DATA)%unit_no
      Write(iunit, '(a,6x,(*(a,3x)))') "#     Time", (Trim(xstd(indx(i))), i=1, neff)
      Do j=1, fit_data%num_input_data
        Write(iunit, '(e13.4,(*(3x,e13.4)))') fit_data%input_data(j)%time, (fit_data%profile(j,indx(i)), i=1, neff)
      End Do
      
      Close(iunit) 
    End If
  
  End Subroutine plot_fittings
  
  Subroutine obtain_fitted_values(fit_data, neff)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Compute fitted value
    !
    ! author    - i.scivetti March 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(fit_type),     Intent(InOut) :: fit_data
    Integer(Kind=wi),   Intent(In   ) :: neff

    Real(Kind=wp)    :: A0, A1, A2, A3, T2, T3, T1, t, fit
    Integer(Kind=wi) :: j
    
    A3=fit_data%param%A3%value(neff)
    A2=fit_data%param%A2%value(neff)
    A1=fit_data%param%A1%value(neff) 
    A0=fit_data%param%A0%value(neff)
    T1=fit_data%param%T1%value(neff)
    T2=fit_data%param%T2%value(neff)
    T3=fit_data%param%T3%value(neff)
    
    Do j=1, fit_data%num_input_data
      t=fit_data%input_data(j)%time
      If (Trim(fit_data%functional_form) == 'ocf_2exp_to_constant') Then
        If (fit_data%fit_superfast_time%stat) Then
          fit= A3+(A2-A3)*exp(-t/T3)+(A1-A2)*exp(-t/T2-t/T3)+(1.0_wp-A1)*exp(-t/T1-t/T2-t/T3)
        Else
          fit= A3+(A2-A3)*exp(-t/T3)+(A1-A2)*exp(-t/T2-t/T3)
        End If
      Else If (Trim(fit_data%functional_form) == 'ocf_2exp_to_zero') Then
        If (fit_data%fit_superfast_time%stat) Then
          fit= (A2-A3)*exp(-t/T3)+(A1-A2)*exp(-t/T2-t/T3)+(1.0_wp-A1)*exp(-t/T1-t/T2-t/T3)
        Else
          fit= (A2-A3)*exp(-t/T3)+(A1-A2)*exp(-t/T2-t/T3)
        End If
      Else If (Trim(fit_data%functional_form) == 'ocf_1exp_to_constant') Then
        If (fit_data%fit_superfast_time%stat) Then
          fit=A2+(A1-A2)*exp(-t/T2)+(1.0_wp-A1)*exp(-t/T1-t/T2)
        Else
          fit=A2+(A1-A2)*exp(-t/T2)
        End If
      Else If (Trim(fit_data%functional_form) == 'ocf_1exp_to_zero') Then
        If (fit_data%fit_superfast_time%stat) Then
          fit=A1*exp(-t/T2)+(1.0_wp-A1)*exp(-t/T1-t/T2)
        Else
          fit=A1*exp(-t/T2)
        End If
      Else If (Trim(fit_data%functional_form) == 'tcf_2exp_to_zero') Then
        fit=A0*exp(-A1*t)+(1-A0)*exp(-A2*t)
      Else If (Trim(fit_data%functional_form) == 'tcf_1exp_to_constant') Then
        fit=A0*exp(-A1*t)+(1-A0)
      Else If (Trim(fit_data%functional_form) == 'tcf_1exp_to_zero') Then
        fit=exp(-A1*t)
      Else If (Trim(fit_data%functional_form) == 'spcf_2exp_to_zero') Then
        fit=A0*exp(-A1*t)+(1-A0)*exp(-A2*t)
      Else If (Trim(fit_data%functional_form) == 'spcf_1exp_to_constant') Then
        fit=A0*exp(-A1*t)+(1-A0)
      Else If (Trim(fit_data%functional_form) == 'spcf_1exp_to_zero') Then
        fit=exp(-A1*t)        
      Else If (Trim(fit_data%functional_form) == 'msd_linear') Then
        fit=A0 + A1*t
      End If
      fit_data%profile(j,neff)=fit
    End Do
    
  End Subroutine obtain_fitted_values
  
  Subroutine plot_raw_data(files, fit_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Plot all the raw data: deviations from the average value by multiples of
    ! the STD
    !
    ! author    - i.scivetti March 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),    Intent(InOut) :: files(:)
    Type(fit_type),     Intent(InOut) :: fit_data
  
    Integer :: iunit, i, j
    Character(Len=256) :: xstd(fit_data%std_samples%value)
  
    Do i=1, fit_data%std_samples%value 
      fit_data%xstd=(i-1)*fit_data%delta_std-fit_data%max_times_std%value
      Write(xstd(i),'(f12.3)') fit_data%xstd
      If (fit_data%xstd < 0.0_wp) Then
        xstd(i)='AVG'//Trim(Adjustl(xstd(i)))//'*STD'
      Else If (Abs(fit_data%xstd) < epsilon(1.0_wp)) Then
        xstd(i)='AVG+'//Trim(Adjustl(xstd(i)))//'*STD'
      Else
        xstd(i)='AVG+'//Trim(Adjustl(xstd(i)))//'*STD'
      End If
      Do j=1, fit_data%num_input_data
        fit_data%raw_data(j,i)=fit_data%input_data(j)%avg+fit_data%xstd*fit_data%input_data(j)%std
      End Do
    End Do
  
    Open(Newunit=files(FILE_RAW_DATA)%unit_no, File=files(FILE_RAW_DATA)%filename, Status='Replace')
    iunit=files(FILE_RAW_DATA)%unit_no
    Write(iunit, '(a,6x,(*(a,3x)))') "#     Time", (Trim(xstd(i)), i=1, fit_data%std_samples%value)
    Do j=1, fit_data%num_input_data
      Write(iunit, '(e13.4,(*(3x,e13.4)))') fit_data%input_data(j)%time, (fit_data%raw_data(j,i), i=1, fit_data%std_samples%value)
    End Do
  
    Close(iunit) 
  
  End Subroutine plot_raw_data

  Subroutine compute_average_over_stds(fit_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to compute the mean associated with fittings
    !
    ! author    - i.scivetti March 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(fit_type),     Intent(InOut) :: fit_data

    Integer(kind=wi)   :: i
    
    fit_data%param%A3%suma=0.0_wp
    fit_data%param%A2%suma=0.0_wp
    fit_data%param%A1%suma=0.0_wp
    fit_data%param%A0%suma=0.0_wp
    fit_data%param%T1%suma=0.0_wp
    fit_data%param%T2%suma=0.0_wp
    fit_data%param%T3%suma=0.0_wp
    fit_data%sum_weights=0.0_wp
    fit_data%average_ssr=0.0_wp
    
    fit_data%param%valid=fit_data%param%fit_ok
    fit_data%valid_solutions=0
    
    Do i = 1, fit_data%std_samples%value
      If (fit_data%param%fit_ok(i)) Then
        Call check_validity_parameters(fit_data%param%A0%value(i), fit_data%param%A1%value(i),&
                                      &fit_data%param%A2%value(i), fit_data%param%A3%value(i),&
                                      &fit_data%param%T1%value(i), fit_data%param%T2%value(i),&
                                      &fit_data%param%T3%value(i),&
                                      &fit_data%functional_form, fit_data%param%valid(i),&
                                      &fit_data%fit_superfast_time%stat)  
      End If
      If (fit_data%param%valid(i)) Then
        fit_data%valid_solutions=fit_data%valid_solutions+1
        fit_data%xstd=(i-1)*fit_data%delta_std-fit_data%max_times_std%value
        fit_data%sum_weights=fit_data%sum_weights+exp(-0.5*fit_data%xstd**2)
        fit_data%param%A3%suma=fit_data%param%A3%suma+fit_data%param%A3%value(i)*exp(-0.5*fit_data%xstd**2)
        fit_data%param%A2%suma=fit_data%param%A2%suma+fit_data%param%A2%value(i)*exp(-0.5*fit_data%xstd**2)
        fit_data%param%A1%suma=fit_data%param%A1%suma+fit_data%param%A1%value(i)*exp(-0.5*fit_data%xstd**2)
        fit_data%param%A0%suma=fit_data%param%A0%suma+fit_data%param%A0%value(i)*exp(-0.5*fit_data%xstd**2)
        fit_data%param%T2%suma=fit_data%param%T2%suma+fit_data%param%T2%value(i)*exp(-0.5*fit_data%xstd**2)
        fit_data%param%T3%suma=fit_data%param%T3%suma+fit_data%param%T3%value(i)*exp(-0.5*fit_data%xstd**2)
        If (fit_data%fit_superfast_time%stat) Then
          fit_data%param%T1%suma=fit_data%param%T1%suma+fit_data%param%T1%value(i)*exp(-0.5*fit_data%xstd**2)
        End If
        If (fit_data%print_average_ssr%stat) Then
          fit_data%average_ssr=fit_data%average_ssr+fit_data%param%ssr(i)
        End If
      End If  
    End Do
    
    If (fit_data%valid_solutions /= 0) Then 
      fit_data%param%A3%mean=fit_data%param%A3%suma/fit_data%sum_weights
      fit_data%param%A2%mean=fit_data%param%A2%suma/fit_data%sum_weights
      fit_data%param%A1%mean=fit_data%param%A1%suma/fit_data%sum_weights
      fit_data%param%A0%mean=fit_data%param%A0%suma/fit_data%sum_weights
      fit_data%param%T2%mean=fit_data%param%T2%suma/fit_data%sum_weights
      fit_data%param%T3%mean=fit_data%param%T3%suma/fit_data%sum_weights
      If (fit_data%fit_superfast_time%stat) Then
        fit_data%param%T1%mean=fit_data%param%T1%suma/fit_data%sum_weights  
      End If
      If (fit_data%print_average_ssr%stat) Then
        fit_data%average_ssr=fit_data%average_ssr/fit_data%valid_solutions
      End If      
    End If 
    
  End Subroutine compute_average_over_stds
  
  Subroutine compute_std_over_stds(fit_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to compute the STD associated with fittings
    !
    ! author    - i.scivetti March 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(fit_type),     Intent(InOut) :: fit_data

    Integer(kind=wi)   :: i
    
    fit_data%param%A3%suma=0.0_wp
    fit_data%param%A2%suma=0.0_wp
    fit_data%param%A1%suma=0.0_wp
    fit_data%param%A0%suma=0.0_wp
    fit_data%param%T1%suma=0.0_wp
    fit_data%param%T2%suma=0.0_wp
    fit_data%param%T3%suma=0.0_wp

    Do i = 1, fit_data%std_samples%value
      If (fit_data%param%valid(i)) Then
        fit_data%xstd=(i-1)*fit_data%delta_std-fit_data%max_times_std%value
        fit_data%param%A3%suma=fit_data%param%A3%suma+ &
        & (fit_data%param%A3%value(i)-fit_data%param%A3%mean)**2 * exp(-0.5*fit_data%xstd**2)
        fit_data%param%A2%suma=fit_data%param%A2%suma+ &
        & (fit_data%param%A2%value(i)-fit_data%param%A2%mean)**2 * exp(-0.5*fit_data%xstd**2)
        fit_data%param%A1%suma=fit_data%param%A1%suma+ &
        & (fit_data%param%A1%value(i)-fit_data%param%A1%mean)**2 * exp(-0.5*fit_data%xstd**2)
        fit_data%param%A0%suma=fit_data%param%A0%suma+ &
        & (fit_data%param%A0%value(i)-fit_data%param%A0%mean)**2 * exp(-0.5*fit_data%xstd**2)
        fit_data%param%T2%suma=fit_data%param%T2%suma+ &
        & (fit_data%param%T2%value(i)-fit_data%param%T2%mean)**2 * exp(-0.5*fit_data%xstd**2)
        fit_data%param%T3%suma=fit_data%param%T3%suma+ &
        & (fit_data%param%T3%value(i)-fit_data%param%T3%mean)**2 * exp(-0.5*fit_data%xstd**2)
        If (fit_data%fit_superfast_time%stat) Then
          fit_data%param%T1%suma=fit_data%param%T1%suma+ &
          & (fit_data%param%T1%value(i)-fit_data%param%T1%mean)**2 * exp(-0.5*fit_data%xstd**2)        
        End If
      End If
    End Do
    
    If (fit_data%valid_solutions /= 0) Then
      fit_data%param%A3%std=sqrt(fit_data%param%A3%suma)/sqrt(fit_data%sum_weights)
      fit_data%param%A2%std=sqrt(fit_data%param%A2%suma)/sqrt(fit_data%sum_weights)
      fit_data%param%A1%std=sqrt(fit_data%param%A1%suma)/sqrt(fit_data%sum_weights)
      fit_data%param%A0%std=sqrt(fit_data%param%A0%suma)/sqrt(fit_data%sum_weights)
      fit_data%param%T2%std=sqrt(fit_data%param%T2%suma)/sqrt(fit_data%sum_weights)
      fit_data%param%T3%std=sqrt(fit_data%param%T3%suma)/sqrt(fit_data%sum_weights)
      If (fit_data%fit_superfast_time%stat) Then
        fit_data%param%T1%std=sqrt(fit_data%param%T1%suma)/sqrt(fit_data%sum_weights)
      End If      
    End If
    
  End Subroutine compute_std_over_stds  
  
  Subroutine fit_vs_std(files, fit_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Fit the correlation for various multiples of its STD
    !
    ! author    - i.scivetti March 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),    Intent(InOut) :: files(:)
    Type(fit_type),     Intent(InOut) :: fit_data

    Integer(kind=wi)   :: i, j, io, Nfail, k(3,2)
    Character(Len=256) :: string
    Character(Len=256) :: message, word, messages(3)
    Character(Len=256) :: ndata, ndf
    
    Call delete_gnuplot_working_files(files)
 
    Nfail=0 
    If (fit_data%std_samples%value /= 1) Then
      fit_data%delta_std=2.0_wp*fit_data%max_times_std%value/(fit_data%std_samples%value - 1) 
    Else
      fit_data%delta_std=0.0_wp
    End If
    Call info(' ', 1) 
    Call info('Start of fitting', 1) 
    Call info('================', 1) 

    ! Check times
    Call check_times_involved(fit_data)
    
    ! start by fitting the average correlation to check if everything is ok(ish)
    fit_data%xstd=0.0_wp
    Write(fit_data%xstd_char,'(f6.3)') fit_data%xstd
    Call generate_gnuplot_input(files, fit_data)
    
    Call execute_command_line(Trim(fit_data%exec_gnuplot))

    ! Check if the AVG converges and does not fail for the functional form 
    Open(Newunit=files(FILE_DUMP)%unit_no, File=files(FILE_DUMP)%filename, Status='Replace')
    Call execute_command_line(Trim(fit_data%check_fit))
    Read (files(FILE_DUMP)%unit_no, Fmt=*, iostat=io) string
    Close(files(FILE_DUMP)%unit_no, Status='Delete')

    If (.Not. is_iostat_end(io)) Then
      Write(message, '(a)') 'PROBLEMS! Fitting for the average correlation has NOT CONVERGED/FAILED!&
                           & Please check the input data and definition of the settings.' 
      Open(Newunit=files(FILE_DUMP)%unit_no, File=files(FILE_DUMP)%filename, Status='Replace')
      Call execute_command_line(Trim(fit_data%check_window))
      Read (files(FILE_DUMP)%unit_no, Fmt=*, iostat=io) string
      Close(files(FILE_DUMP)%unit_no, Status='Delete')
      If (.Not. is_iostat_end(io)) Then
         Write(message, '(a)') 'PROBLEMS! There is no data within the specified time window&
                              & defined within the time domain. Please check&
                              & the "'//Trim(fit_data%filename%type)//'" as well as the numbers assigned to "column" directives.'
      End If 
      Call info (message, 1)
      If (is_iostat_end(io)) Then
        If (fit_data%input_parameters%fread) Then
          Write(message, '(a)') 'Moreover, the user should either i) revise the values of the "&input_parameters" block or&
                              & ii) try removing the "&input_parameters" block (if this has not been tried already).'
        Else
          Write(message, '(a)') 'If problems persist, the user should also consider using the "&input_parameters" block.'
        End If
        Call info (message, 1) 
      End If
      Call delete_gnuplot_working_files(files)
      Call error_stop(' ')
    Else
      ! Extract the number of degrees of freedom
      Call execute_command_line(Trim(fit_data%check_ndf))
      Open(Newunit=files(FILE_DUMP)%unit_no, File=files(FILE_DUMP)%filename, Status='Old')
      Read (files(FILE_DUMP)%unit_no, Fmt=*) (word, i=1,5), fit_data%ndf
      Close(files(FILE_DUMP)%unit_no, Status='Delete')
      Call execute_command_line(fit_data%extract_results)
      Open(Newunit=files(FILE_FIT_PARAMS)%unit_no, File=files(FILE_FIT_PARAMS)%filename, Status='Old')
      Call read_results_from_list(files(FILE_FIT_PARAMS)%unit_no, fit_data, (fit_data%std_samples%value+1)/2)
      Call redefine_initial_parameters(fit_data, (fit_data%std_samples%value+1)/2)
      Call delete_gnuplot_working_files(files)
    End If

    ! If the fitting for the AVG correlation worked, let's explore variations of it in times*std 
    k(1,1)= (fit_data%std_samples%value+1)/2-1
    k(2,1)= 1
    k(3,1)=-1
    k(1,2)= (fit_data%std_samples%value+1)/2+1
    k(2,2)= fit_data%std_samples%value
    k(3,2)= 1

    Write(ndata, Fmt=*) fit_data%ndf+fit_data%nparam
    Write(ndf, Fmt=*)   fit_data%ndf
    
    ! Things are ready to print headers
    Write (messages(1),'(1x,a)') 'Number of data points: '//Trim(Adjustl(ndata))  
    Write (messages(2),'(1x,a)') 'Number of degrees of freedom (NDF): '//Trim(Adjustl(ndf))  
    Write (messages(3),'(1x,a)') 'Fitting formula: '//Trim(fit_data%formula_show)
    Call info(messages, 3)
    
    Do j=1,2
      Do i = k(1,j), k(2,j), k(3,j)
        fit_data%xstd=(i-1)*fit_data%delta_std-fit_data%max_times_std%value
        Write(fit_data%xstd_char,'(f6.3)') fit_data%xstd
        Call generate_gnuplot_input(files, fit_data)
        Call execute_command_line(Trim(fit_data%exec_gnuplot))
        Open(Newunit=files(FILE_DUMP)%unit_no, File=files(FILE_DUMP)%filename, Status='Replace')
        Call execute_command_line(Trim(fit_data%check_fit))
        Read (files(FILE_DUMP)%unit_no, Fmt=*, iostat=io) string
        Close(files(FILE_DUMP)%unit_no, Status='Delete')
        If (.Not. is_iostat_end(io)) Then
          If (Abs(fit_data%xstd) < epsilon(1.0_wp)) Then
            Write(message, '(a)') 'WARNING! Fitting for the average correlation has NOT CONVERGED/FAILED! ACTION:&
                                & discarded from the error analysis'
          
          Else 
            Write(message, '(a)') 'WARNING! Fitting NOT CONVERGED/FAILED for the correlation corresponding to a deviation of "'&
                                &//Trim(Adjustl(fit_data%xstd_char))//'*STD" (with respect to the average correlation). ACTION:&
                                & discarded from the error analysis'
          End If  
          fit_data%param%fit_ok(i)=.False.
          Nfail=Nfail+1
          Call info (message, 1)
          Call execute_command_line(fit_data%rm_gnuplot_output)
        Else  
           Call execute_command_line(fit_data%extract_results)
           Open(Newunit=files(FILE_FIT_PARAMS)%unit_no, File=files(FILE_FIT_PARAMS)%filename, Status='Old')
           Call read_results_from_list(files(FILE_FIT_PARAMS)%unit_no, fit_data, i)
           Call redefine_initial_parameters(fit_data, i)
           Call delete_gnuplot_working_files(files)
        End If
        Call refresh_out(files)
      End Do
      Call redefine_initial_parameters(fit_data, (fit_data%std_samples%value+1)/2)
    End Do

    If (Nfail/=0) Then
       Call info(' ', 1)
       Write(message, '(a)') 'The user should check the input value for "max_times_std" as well as the values in the "'&
                             &//Trim(fit_data%filename%type)//'" file'
       Call info (message, 1)                      
    End If

    If(Nfail==fit_data%std_samples%value) Then
       Call info(' ', 1)
       Write(message, '(a)') '***ERROR: NONE of the fittings has converged! Please review the settings'
       Call info (message, 1)
       Call error_stop(' ')
    End If

    If (Nfail/=0) Then
       Call info(' ', 1)
       Write(message, '(a)') 'CAUTION: it is likely that the estimation of the mean values and std for&
                            & the parameters is NOT CORRECT'
       Call info (message, 1)
    End If
    
  End Subroutine fit_vs_std
   
  Subroutine set_fitting_directives(files, fit_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Set the directives for fitting
    !
    ! author    - i.scivetti March 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),    Intent(InOut) :: files(:)
    Type(fit_type),     Intent(InOut) :: fit_data
  
 
    If (Trim(fit_data%functional_form) == 'ocf_1exp_to_constant') Then
      If (fit_data%fit_superfast_time%stat) Then
        Write(fit_data%header,'(2(a,8x),10x,4(a,20x))') 'Deviation (xSTD)', 'SSR/NDF', 'A2','A1','T2', 'T1' 
        Write(fit_data%separation_line,Fmt=*) Repeat("=", 135)   
        fit_data%formula='f(x) = A2+(A1-A2)*exp(-x/T2)+(1-A1)*exp(-x/T1-x/T2)'
        fit_data%formula_show='f(t)=A2+(A1-A2)*exp(-t/T2)+(1-A1)*exp(-t/T1-t/T2)' 
        Write(fit_data%init_param, '(4(a,e12.4),a)')'A2=',fit_data%param%A2%ini%value,&
                                                & '; A1=',fit_data%param%A1%ini%value,&
                                                & '; T2=',fit_data%param%T2%ini%value,&
                                                & '; T1=',fit_data%param%T1%ini%value,';'
        fit_data%via_params='via A2, A1, T2, T1;'
        Write(fit_data%param_list,'(4(a,14x),a)') '\|A2', '=\|A1', '=\|T2', '=\|T1', '='
        fit_data%nparam=4      
      Else
        Write(fit_data%header,'(2(a,8x),10x,3(a,20x))') 'Deviation (xSTD)', 'SSR/NDF', 'A2','A1','T2' 
        Write(fit_data%separation_line,Fmt=*) Repeat("=", 110)   
        fit_data%formula='f(x) = A2+(A1-A2)*exp(-x/T2)'
        fit_data%formula_show='f(t)=A2+(A1-A2)*exp(-t/T2)' 
        Write(fit_data%init_param, '(3(a,e12.4),a)')'A2=',fit_data%param%A2%ini%value,&
                                                & '; A1=',fit_data%param%A1%ini%value,&
                                                & '; T2=',fit_data%param%T2%ini%value,';'
        fit_data%via_params='via A2, A1, T2;'
        Write(fit_data%param_list,'(3(a,14x),a)') '\|A2', '=\|A1', '=\|T2', '='
        fit_data%nparam=3
      End If
    Else If (Trim(fit_data%functional_form) == 'ocf_1exp_to_zero') Then
      If (fit_data%fit_superfast_time%stat) Then
        Write(fit_data%separation_line,Fmt=*) Repeat("=", 115)   
        Write(fit_data%header,'(2(a,8x),10x,3(a,20x))') 'Deviation (xSTD)', 'SSR/NDF', 'A1','T2', 'T1'
        fit_data%formula='f(x) = A1*exp(-x/T2)+(1-A1)*exp(-x/T1-x/T2)'
        fit_data%formula_show='f(t) = A1*exp(-t/T2)+(1-A1)*exp(-t/T1-t/T2)'
        Write(fit_data%init_param, '(3(a,e12.4),a)')'A1=',fit_data%param%A1%ini%value,&
                                                & '; T2=',fit_data%param%T2%ini%value,&
                                                & '; T1=',fit_data%param%T1%ini%value,';'
        fit_data%via_params='via A1, T2, T1;'
        Write(fit_data%param_list,'(3(a,14x),a)') '\|A1', '=\|T2', '=\|T1', '='
        fit_data%nparam=3
      Else
        Write(fit_data%separation_line,Fmt=*) Repeat("=", 86)   
        Write(fit_data%header,'(2(a,8x),10x,2(a,20x))') 'Deviation (xSTD)', 'SSR/NDF', 'A1','T2'
        fit_data%formula='f(x) = A1*exp(-x/T2)'
        fit_data%formula_show='f(t) = A1*exp(-t/T2)'
        Write(fit_data%init_param, '(2(a,e12.4),a)')'A1=',fit_data%param%A1%ini%value,&
                                                & '; T2=',fit_data%param%T2%ini%value,';'
        fit_data%via_params='via A1, T2;'
        Write(fit_data%param_list,'(2(a,14x),a)') '\|A1', '=\|T2', '='
        fit_data%nparam=2
      End If
    Else If (Trim(fit_data%functional_form) == 'ocf_2exp_to_constant') Then
      If (fit_data%fit_superfast_time%stat) Then
        Write(fit_data%separation_line,Fmt=*) Repeat("=", 180)
        Write(fit_data%header,'(2(a,8x),10x,6(a,20x))') 'Deviation (xSTD)', 'SSR/NDF', 'A3','A2','A1','T3','T2', 'T1'
        fit_data%formula='f(x)=A3+(A2-A3)*exp(-x/T3)+(A1-A2)*exp(-x/T2-x/T3)+(1-A1)*exp(-x/T1-x/T2-x/T3)'
        fit_data%formula_show='f(t) = A3+(A2-A3)*exp(-t/T3)+(A1-A2)*exp(-t/T2-t/T3)+(1-A1)*exp(-t/T1-t/T2-t/T3)'
        Write(fit_data%init_param, '(6(a,e12.4),a)')'A3=',fit_data%param%A3%ini%value,&
                                                & '; A2=',fit_data%param%A2%ini%value,&
                                                & '; A1=',fit_data%param%A1%ini%value,&
                                                & '; T3=',fit_data%param%T3%ini%value,&
                                                & '; T2=',fit_data%param%T2%ini%value,&
                                                & '; T1=',fit_data%param%T1%ini%value,';'
        fit_data%via_params='via A3, A2, A1, T3, T2, T1;'
        Write(fit_data%param_list,'(6(a,14x),a)') '\|A3','=\|A2', '=\|A1', '=\|T3', '=\|T2', '=\|T1', '='
        fit_data%nparam=6      
      Else
        Write(fit_data%separation_line,Fmt=*) Repeat("=", 152)
        Write(fit_data%header,'(2(a,8x),10x,5(a,20x))') 'Deviation (xSTD)', 'SSR/NDF', 'A3','A2','A1','T3','T2'
        fit_data%formula='f(x)=A3+(A2-A3)*exp(-x/T3)+(A1-A2)*exp(-x/T2-x/T3)'
        fit_data%formula_show='f(t) = A3+(A2-A3)*exp(-t/T3)+(A1-A2)*exp(-t/T2-t/T3)'
        Write(fit_data%init_param, '(5(a,e12.4),a)')'A3=',fit_data%param%A3%ini%value,&
                                                & '; A2=',fit_data%param%A2%ini%value,&
                                                & '; A1=',fit_data%param%A1%ini%value,&
                                                & '; T3=',fit_data%param%T3%ini%value,&
                                                & '; T2=',fit_data%param%T2%ini%value,';'
        fit_data%via_params='via A3, A2, A1, T3, T2;'
        Write(fit_data%param_list,'(5(a,14x),a)') '\|A3','=\|A2', '=\|A1', '=\|T3', '=\|T2', '='
        fit_data%nparam=5
      End If
    Else If (Trim(fit_data%functional_form) == 'ocf_2exp_to_zero') Then
      If (fit_data%fit_superfast_time%stat) Then
        Write(fit_data%separation_line,Fmt=*) Repeat("=", 160)   
        Write(fit_data%header,'(2(a,8x),10x,5(a,20x))') 'Deviation (xSTD)','SSR/NDF', 'A2','A1','T3','T2', 'T1'
        fit_data%formula='f(x)=A2*exp(-x/T3)+(A1-A2)*exp(-x/T2-x/T3)+(1-A1)*exp(-x/T1-x/T2-x/T3)'
        fit_data%formula_show='f(t) = A2*exp(-t/T3)+(A1-A2)*exp(-t/T2-t/T3)+(1-A1)*exp(-t/T1-t/T2-t/T3)'
        Write(fit_data%init_param, '(5(a,e12.4),a)')'A2=',fit_data%param%A2%ini%value,&
                                                & '; A1=',fit_data%param%A1%ini%value,&
                                                & '; T3=',fit_data%param%T3%ini%value,&
                                                & '; T2=',fit_data%param%T2%ini%value,&
                                                & '; T1=',fit_data%param%T1%ini%value,';'
        fit_data%via_params='via A2, A1, T3, T2, T1;'
        Write(fit_data%param_list,'(5(a,14x),a)') '\|A2', '=\|A1', '=\|T3', '=\|T2', '=\|T1', '='
        fit_data%nparam=5
      Else
        Write(fit_data%separation_line,Fmt=*) Repeat("=", 128)   
        Write(fit_data%header,'(2(a,8x),10x,4(a,20x))') 'Deviation (xSTD)','SSR/NDF', 'A2','A1','T3','T2'
        fit_data%formula='f(x)=A2*exp(-x/T3)+(A1-A2)*exp(-x/T2-x/T3)'
        fit_data%formula_show='f(t) = A2*exp(-t/T3)+(A1-A2)*exp(-t/T2-t/T3)'
        Write(fit_data%init_param, '(4(a,e12.4),a)')'A2=',fit_data%param%A2%ini%value,&
                                                & '; A1=',fit_data%param%A1%ini%value,&
                                                & '; T3=',fit_data%param%T3%ini%value,&
                                                & '; T2=',fit_data%param%T2%ini%value,';'
        fit_data%via_params='via A2, A1, T3, T2;'
        Write(fit_data%param_list,'(4(a,14x),a)') '\|A2', '=\|A1', '=\|T3', '=\|T2', '='
        fit_data%nparam=4
      End If
    Else If (Trim(fit_data%functional_form) == 'tcf_2exp_to_zero') Then
      Write(fit_data%header,'(2(a,8x),10x,3(a,20x))') 'Deviation (xSTD)', 'SSR/NDF', 'A0','A1','A2' 
      Write(fit_data%separation_line,Fmt=*) Repeat("=", 110)   
      fit_data%formula='f(x) = A0*exp(-A1*x)+(1-A0)*exp(-A2*x)'
      fit_data%formula_show='f(t)=A0*exp(-A1*t)+(1-A0)*exp(-A2*t)' 
      Write(fit_data%init_param, '(3(a,e12.4),a)')'A0=',fit_data%param%A0%ini%value,&
                                              & '; A1=',fit_data%param%A1%ini%value,&
                                              & '; A2=',fit_data%param%A2%ini%value,';'
      fit_data%via_params='via A0, A1, A2;'
      Write(fit_data%param_list,'(3(a,14x),a)') '\|A0', '=\|A1', '=\|A2', '='
      fit_data%nparam=3
    Else If (Trim(fit_data%functional_form) == 'tcf_1exp_to_constant') Then
      Write(fit_data%header,'(2(a,8x),10x,2(a,20x))') 'Deviation (xSTD)', 'SSR/NDF', 'A0','A1' 
      Write(fit_data%separation_line,Fmt=*) Repeat("=", 86)   
      fit_data%formula='f(x) = A0*exp(-A1*x)+(1-A0)'
      fit_data%formula_show='f(t)=A0*exp(-A1*t)+(1-A0)' 
      Write(fit_data%init_param, '(2(a,e12.4),a)')'A0=',fit_data%param%A0%ini%value,&
                                              & '; A1=',fit_data%param%A1%ini%value,';'
      fit_data%via_params='via A0, A1;'
      Write(fit_data%param_list,'(2(a,14x),a)') '\|A0', '=\|A1', '='
      fit_data%nparam=2
    Else If (Trim(fit_data%functional_form) == 'tcf_1exp_to_zero') Then
      Write(fit_data%header,'(2(a,8x),10x,1(a,20x))') 'Deviation (xSTD)', 'SSR/NDF', 'A1' 
      Write(fit_data%separation_line,Fmt=*) Repeat("=", 62)   
      fit_data%formula='f(x) = exp(-A1*x)'
      fit_data%formula_show='f(t)=exp(-A1*t)' 
      Write(fit_data%init_param, '(1(a,e12.4),a)')'A1=',fit_data%param%A1%ini%value,';'
      fit_data%via_params='via A1;'
      Write(fit_data%param_list,'(1(a,14x),a)') '\|A1', '='
      fit_data%nparam=1
    Else If (Trim(fit_data%functional_form) == 'spcf_2exp_to_zero') Then
      Write(fit_data%header,'(2(a,8x),10x,3(a,20x))') 'Deviation (xSTD)', 'SSR/NDF', 'A0','A1','A2' 
      Write(fit_data%separation_line,Fmt=*) Repeat("=", 110)   
      fit_data%formula='f(x) = A0*exp(-A1*x)+(1-A0)*exp(-A2*x)'
      fit_data%formula_show='f(t)=A0*exp(-A1*t)+(1-A0)*exp(-A2*t)' 
      Write(fit_data%init_param, '(3(a,e12.4),a)')'A0=',fit_data%param%A0%ini%value,&
                                              & '; A1=',fit_data%param%A1%ini%value,&
                                              & '; A2=',fit_data%param%A2%ini%value,';'
      fit_data%via_params='via A0, A1, A2;'
      Write(fit_data%param_list,'(3(a,14x),a)') '\|A0', '=\|A1', '=\|A2', '='
      fit_data%nparam=3
    Else If (Trim(fit_data%functional_form) == 'spcf_1exp_to_constant') Then
      Write(fit_data%header,'(2(a,8x),10x,2(a,20x))') 'Deviation (xSTD)', 'SSR/NDF', 'A0','A1' 
      Write(fit_data%separation_line,Fmt=*) Repeat("=", 86)   
      fit_data%formula='f(x) = A0*exp(-A1*x)+(1-A0)'
      fit_data%formula_show='f(t)=A0*exp(-A1*t)+(1-A0)' 
      Write(fit_data%init_param, '(2(a,e12.4),a)')'A0=',fit_data%param%A0%ini%value,&
                                              & '; A1=',fit_data%param%A1%ini%value,';'
      fit_data%via_params='via A0, A1;'
      Write(fit_data%param_list,'(2(a,14x),a)') '\|A0', '=\|A1', '='
      fit_data%nparam=2
    Else If (Trim(fit_data%functional_form) == 'spcf_1exp_to_zero') Then
      Write(fit_data%header,'(2(a,8x),10x,1(a,20x))') 'Deviation (xSTD)', 'SSR/NDF', 'A1' 
      Write(fit_data%separation_line,Fmt=*) Repeat("=", 62)   
      fit_data%formula='f(x) = exp(-A1*x)'
      fit_data%formula_show='f(t)=exp(-A1*t)' 
      Write(fit_data%init_param, '(1(a,e12.4),a)')'A1=',fit_data%param%A1%ini%value,';'
      fit_data%via_params='via A1;'
      Write(fit_data%param_list,'(1(a,14x),a)') '\|A1', '='
      fit_data%nparam=1      
    Else If (Trim(fit_data%functional_form) == 'msd_linear') Then
      Write(fit_data%header,'(2(a,8x),10x,2(a,20x))') 'Deviation (xSTD)', 'SSR/NDF', 'A0','A1' 
      Write(fit_data%separation_line,Fmt=*) Repeat("=", 86)   
      fit_data%formula='f(x) = A0+A1*x'
      fit_data%formula_show='f(t)=A0+A1*t' 
      Write(fit_data%init_param, '(2(a,e12.4),a)')'A0=',fit_data%param%A0%ini%value,&
                                              & '; A1=',fit_data%param%A1%ini%value,';'
      fit_data%via_params='via A0, A1;'
      Write(fit_data%param_list,'(2(a,14x),a)') '\|A0', '=\|A1', '='
      fit_data%nparam=2
    End If

    Call define_gnuplot_actions(files, fit_data)

  End Subroutine set_fitting_directives                            
  
  Subroutine read_results_from_list(iunit, fit_data, i)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to read parameters from the generated list
    !
    ! author    - i.scivetti March 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Integer(Kind=wi),    Intent(In   ) :: iunit
    Type(fit_type),      Intent(InOut) :: fit_data
    Integer(Kind=wi),    Intent(In   ) :: i
  
    Character(Len=256) :: word
    Integer(kind=wi)   :: j

    Read(iunit,Fmt=*) (word, j=1,7), fit_data%param%ssr(i)
    fit_data%param%ssr(i)=fit_data%param%ssr(i)/fit_data%ndf
    If (Trim(fit_data%functional_form) == 'ocf_2exp_to_constant') Then
      Call read_fitted_parameters(iunit, fit_data%param%A3%value(i), fit_data%param%A3%fit_error(i))
      Call read_fitted_parameters(iunit, fit_data%param%A2%value(i), fit_data%param%A2%fit_error(i))
      Call read_fitted_parameters(iunit, fit_data%param%A1%value(i), fit_data%param%A1%fit_error(i))
      Call read_fitted_parameters(iunit, fit_data%param%T3%value(i), fit_data%param%T3%fit_error(i))
      Call read_fitted_parameters(iunit, fit_data%param%T2%value(i), fit_data%param%T2%fit_error(i))
      If (fit_data%fit_superfast_time%stat)Then
        Call read_fitted_parameters(iunit, fit_data%param%T1%value(i), fit_data%param%T1%fit_error(i))
      End If
    Else If (Trim(fit_data%functional_form) == 'ocf_2exp_to_zero') Then
      Call read_fitted_parameters(iunit, fit_data%param%A2%value(i), fit_data%param%A2%fit_error(i))
      Call read_fitted_parameters(iunit, fit_data%param%A1%value(i), fit_data%param%A1%fit_error(i))
      Call read_fitted_parameters(iunit, fit_data%param%T3%value(i), fit_data%param%T3%fit_error(i))
      Call read_fitted_parameters(iunit, fit_data%param%T2%value(i), fit_data%param%T2%fit_error(i))
      If (fit_data%fit_superfast_time%stat)Then
        Call read_fitted_parameters(iunit, fit_data%param%T1%value(i), fit_data%param%T1%fit_error(i))
      End If      
    Else If (Trim(fit_data%functional_form) == 'ocf_1exp_to_constant') Then
      Call read_fitted_parameters(iunit, fit_data%param%A2%value(i), fit_data%param%A2%fit_error(i))
      Call read_fitted_parameters(iunit, fit_data%param%A1%value(i), fit_data%param%A1%fit_error(i))
      Call read_fitted_parameters(iunit, fit_data%param%T2%value(i), fit_data%param%T2%fit_error(i))
      If (fit_data%fit_superfast_time%stat)Then
        Call read_fitted_parameters(iunit, fit_data%param%T1%value(i), fit_data%param%T1%fit_error(i))
      End If      
    Else If (Trim(fit_data%functional_form) == 'ocf_1exp_to_zero') Then
      Call read_fitted_parameters(iunit, fit_data%param%A1%value(i), fit_data%param%A1%fit_error(i))
      Call read_fitted_parameters(iunit, fit_data%param%T2%value(i), fit_data%param%T2%fit_error(i))
      If (fit_data%fit_superfast_time%stat)Then
        Call read_fitted_parameters(iunit, fit_data%param%T1%value(i), fit_data%param%T1%fit_error(i))
      End If      
    Else If (Trim(fit_data%functional_form) == 'tcf_2exp_to_zero') Then
      Call read_fitted_parameters(iunit, fit_data%param%A0%value(i), fit_data%param%A0%fit_error(i))
      Call read_fitted_parameters(iunit, fit_data%param%A1%value(i), fit_data%param%A1%fit_error(i))
      Call read_fitted_parameters(iunit, fit_data%param%A2%value(i), fit_data%param%A2%fit_error(i))
    Else If (Trim(fit_data%functional_form) == 'tcf_1exp_to_constant') Then
      Call read_fitted_parameters(iunit, fit_data%param%A0%value(i), fit_data%param%A0%fit_error(i))
      Call read_fitted_parameters(iunit, fit_data%param%A1%value(i), fit_data%param%A1%fit_error(i))
    Else If (Trim(fit_data%functional_form) == 'tcf_1exp_to_zero') Then
      Call read_fitted_parameters(iunit, fit_data%param%A1%value(i), fit_data%param%A1%fit_error(i))
    Else If (Trim(fit_data%functional_form) == 'spcf_2exp_to_zero') Then
      Call read_fitted_parameters(iunit, fit_data%param%A0%value(i), fit_data%param%A0%fit_error(i))
      Call read_fitted_parameters(iunit, fit_data%param%A1%value(i), fit_data%param%A1%fit_error(i))
      Call read_fitted_parameters(iunit, fit_data%param%A2%value(i), fit_data%param%A2%fit_error(i))
    Else If (Trim(fit_data%functional_form) == 'spcf_1exp_to_constant') Then
      Call read_fitted_parameters(iunit, fit_data%param%A0%value(i), fit_data%param%A0%fit_error(i))
      Call read_fitted_parameters(iunit, fit_data%param%A1%value(i), fit_data%param%A1%fit_error(i))
    Else If (Trim(fit_data%functional_form) == 'spcf_1exp_to_zero') Then
      Call read_fitted_parameters(iunit, fit_data%param%A1%value(i), fit_data%param%A1%fit_error(i))      
    Else If (Trim(fit_data%functional_form) == 'msd_linear') Then
      Call read_fitted_parameters(iunit, fit_data%param%A0%value(i), fit_data%param%A0%fit_error(i))
      Call read_fitted_parameters(iunit, fit_data%param%A1%value(i), fit_data%param%A1%fit_error(i))
    End If
    
  End Subroutine read_results_from_list
  
  Subroutine redefine_initial_parameters(fit_data, indx)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to redefine_initial parameters using the fittd values of the 
    ! average correlation
    !
    ! author    - i.scivetti March 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(fit_type),     Intent(InOut) :: fit_data
    Integer(Kind=wi),   Intent(In   ) :: indx
  
    If (Trim(fit_data%functional_form) == 'ocf_2exp_to_constant') Then
      If (fit_data%fit_superfast_time%stat)Then
        Write(fit_data%init_param,'(6(a,e13.4),a)') 'A3=', fit_data%param%A3%value(indx),&
                                                & '; A2=', fit_data%param%A2%value(indx),&
                                                & '; A1=', fit_data%param%A1%value(indx),&
                                                & '; T3=', fit_data%param%T3%value(indx),&
                                                & '; T2=', fit_data%param%T2%value(indx),&
                                                & '; T1=', fit_data%param%T1%value(indx),';'
      Else
        Write(fit_data%init_param,'(5(a,e13.4),a)') 'A3=', fit_data%param%A3%value(indx),&
                                                & '; A2=', fit_data%param%A2%value(indx),&
                                                & '; A1=', fit_data%param%A1%value(indx),&
                                                & '; T3=', fit_data%param%T3%value(indx),&
                                                & '; T2=', fit_data%param%T2%value(indx),';'
      End If                                        
    Else If (Trim(fit_data%functional_form) == 'ocf_2exp_to_zero') Then
      If (fit_data%fit_superfast_time%stat) Then
        Write(fit_data%init_param,'(5(a,e13.4),a)') 'A2=', fit_data%param%A2%value(indx),&
                                                & '; A1=', fit_data%param%A1%value(indx),&
                                                & '; T3=', fit_data%param%T3%value(indx),&
                                                & '; T2=', fit_data%param%T3%value(indx),&
                                                & '; T1=', fit_data%param%T1%value(indx),';'      
      Else
        Write(fit_data%init_param,'(4(a,e13.4),a)') 'A2=', fit_data%param%A2%value(indx),&
                                                & '; A1=', fit_data%param%A1%value(indx),&
                                                & '; T3=', fit_data%param%T3%value(indx),&
                                                & '; T2=', fit_data%param%T2%value(indx),';'
      End If                                        
    Else If (Trim(fit_data%functional_form) == 'ocf_1exp_to_constant') Then
      If (fit_data%fit_superfast_time%stat) Then
        Write(fit_data%init_param,'(4(a,e13.4),a)') 'A2=', fit_data%param%A2%value(indx),&
                                                & '; A1=', fit_data%param%A1%value(indx),&
                                                & '; T2=', fit_data%param%T2%value(indx),&
                                                & '; T1=', fit_data%param%T1%value(indx),';'      
      Else
        Write(fit_data%init_param,'(3(a,e13.4),a)') 'A2=', fit_data%param%A2%value(indx),&
                                                & '; A1=', fit_data%param%A1%value(indx),&
                                                & '; T2=', fit_data%param%T2%value(indx),';'
      End If                                        
    Else If (Trim(fit_data%functional_form) == 'ocf_1exp_to_zero') Then
      If (fit_data%fit_superfast_time%stat) Then
        Write(fit_data%init_param,'(3(a,e13.4),a)') 'A1=', fit_data%param%A1%value(indx),&
                                                & '; T2=', fit_data%param%T2%value(indx),&
                                                & '; T1=', fit_data%param%T1%value(indx),';'      
      Else
        Write(fit_data%init_param,'(2(a,e13.4),a)') 'A1=', fit_data%param%A1%value(indx),&
                                                & '; T2=', fit_data%param%T2%value(indx),';'
      End If                                        
    Else If (Trim(fit_data%functional_form) == 'tcf_2exp_to_zero') Then
      Write(fit_data%init_param,'(3(a,e13.4),a)') 'A0=', fit_data%param%A0%value(indx),&
                                              & '; A1=', fit_data%param%A1%value(indx),&
                                              & '; A2=', fit_data%param%A2%value(indx),';'
    Else If (Trim(fit_data%functional_form) == 'tcf_1exp_to_constant') Then
      Write(fit_data%init_param,'(2(a,e13.4),a)') 'A0=', fit_data%param%A0%value(indx),&
                                              & '; A1=', fit_data%param%A1%value(indx),';'
    Else If (Trim(fit_data%functional_form) == 'tcf_1exp_to_zero') Then
      Write(fit_data%init_param,'(1(a,e13.4),a)') 'A1=', fit_data%param%A1%value(indx),';'
    Else If (Trim(fit_data%functional_form) == 'spcf_2exp_to_zero') Then
      Write(fit_data%init_param,'(3(a,e13.4),a)') 'A0=', fit_data%param%A0%value(indx),&
                                              & '; A1=', fit_data%param%A1%value(indx),&
                                              & '; A2=', fit_data%param%A2%value(indx),';'
    Else If (Trim(fit_data%functional_form) == 'spcf_1exp_to_constant') Then
      Write(fit_data%init_param,'(2(a,e13.4),a)') 'A0=', fit_data%param%A0%value(indx),&
                                              & '; A1=', fit_data%param%A1%value(indx),';'
    Else If (Trim(fit_data%functional_form) == 'spcf_1exp_to_zero') Then
      Write(fit_data%init_param,'(1(a,e13.4),a)') 'A1=', fit_data%param%A1%value(indx),';'      
    Else If (Trim(fit_data%functional_form) == 'msd_linear') Then
      Write(fit_data%init_param,'(2(a,e13.4),a)') 'A0=', fit_data%param%A0%value(indx),&
                                              & '; A1=', fit_data%param%A1%value(indx),';'
    End If
  
  End Subroutine redefine_initial_parameters
  
  Subroutine read_fitted_parameters(iunit, value, error)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Read parameters from each line of the dump list
    !
    ! author    - i.scivetti March 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
    Integer(kind=wi),   Intent(In   ) :: iunit
    Real(kind=wp),      Intent(  Out) :: value
    Character(Len=256), Intent(  Out) :: error
  
    Character(Len=256) :: symbol1, symbol2, symbol3
    Character(Len=256) :: line
    
    Read (iunit, Fmt='(a)') line
    Call remove_symbols(line,'/')
    Read (line, Fmt=*) symbol1, symbol2, value,&
                      & symbol1, symbol2, symbol3, error
  
  End Subroutine read_fitted_parameters
  
  Subroutine print_fitting_settings(files, fit_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to print fitting-related directives
    !
    ! author    - i.scivetti March 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),    Intent(In   ) :: files(:)
    Type(fit_type),     Intent(InOut) :: fit_data
  
    Logical             :: valid_params
    Character(Len=256)  :: messages(7)
    Character(Len=256)  :: std_max, np, nit, tstart, tend, dim, error_msg
    
    Call info(' ', 1) 
    Call info('Settings for fitting with GNUPLOT', 1) 
    Call info('=================================', 1) 
    
    Write(np, '(i6)')   fit_data%std_samples%value
    Write(nit, '(i6)')   fit_data%max_iterations%value
    If (fit_data%std_samples%value /= 1) Then
      Write(std_max, '(f6.3)') fit_data%max_times_std%value
    Else
      fit_data%max_times_std%value=0.0_wp
    End If

    Write (messages(1),'(1x,a)') '* Fitting the '//Trim(fit_data%quantity)//' using the formula "'&
                               &//Trim(fit_data%formula_show)//'"'
    Call info(messages, 1)
    If (Trim(fit_data%what_to_fit%type)=='ocf') Then
      If (fit_data%reactive_species%stat) Then
        If (fit_data%fit_superfast_time%stat) Then
          Write (messages(1),'(1x,a)') '* The user has requested to include the superfast time constant T1 in the fitting'
        Else
          Write (messages(1),'(1x,a)') '* The fitting formula implicitly assumes an instantaneous drop in the orientation'
        End If
        Call info(messages, 1)
      Else
        If (fit_data%fit_superfast_time%stat) Then
          Write (messages(1),'(1x,a)') '* The user has requested to include the fast, librational time constant T1 in the fitting'
        Else
          Write (messages(1),'(1x,a)') '* The fitting formula implicitly assumes an instantaneous librational motion'
        End If
        Call info(messages, 1)
      End If
    End If
    
    If (fit_data%std_samples%value /= 1) Then
      Write (messages(1),'(1x,a)') '* Each fitting will be optimised '//Trim(fit_data%via_params)//' a maximum number of '&
                               &//Trim(Adjustl(nit))//' iterations is set'
    Else
      Write (messages(1),'(1x,a)') '* The fitting will be optimised '//Trim(fit_data%via_params)//' a maximum number of '&
                               &//Trim(Adjustl(nit))//' iterations is set'
    End If    
                               
    Write (messages(2),'(1x,a)') '* The data to fit is read from file: '//Trim(fit_data%filename%type)
    
    If (fit_data%std_samples%value /= 1) Then
      Write (messages(3),'(1x,a)') '* A total of '//Trim(Adjustl(np))//&
                               &' deviations from the average values will be considered between +/- '&
                               &//Trim(Adjustl(std_max))//' times the reported STD'  
    Else
      If (fit_data%std_samples%fread) Then
        If (fit_data%max_times_std%fread) Then
          Write (messages(3),'(1x,a)') '* Since the "std_samples" directive was set to 1, only the average correlation&
                                       & will be fitted ("max_times_std" is redundant)'
        Else
          Write (messages(3),'(1x,a)') '* Since the "std_samples" directive was set to 1, only the average correlation&
                                       & will be fitted'
        End If
      Else
        If (fit_data%max_times_std%fread) Then
          Write (messages(3),'(1x,a)') '* Only the average correlation&
                                       & will be fitted ("max_times_std" is redundant)'
        Else
          Write (messages(3),'(1x,a)') '* Only the average correlation will be fitted'
        End If      
      End If
    End If
   
                             
    If (fit_data%input_parameters%fread) Then
      Write (messages(4),'(1x,a)') '* The initial parameters are defined through the "&input_parameters" block' 
      valid_params=.True.
      Call check_validity_parameters(fit_data%param%A0%ini%value, fit_data%param%A1%ini%value,&
                                    &fit_data%param%A2%ini%value, fit_data%param%A3%ini%value,&
                                    &fit_data%param%T1%ini%value, fit_data%param%T2%ini%value,&
                                    &fit_data%param%T3%ini%value, &
                                    &fit_data%functional_form, valid_params,&
                                    &fit_data%fit_superfast_time%stat, error_msg)
      If (.Not. valid_params) Then
         Call info(messages, 4)
         Call info(error_msg, 1)
         Call error_stop(' ')
      End If
    Else
      Write (messages(4),'(1x,a)') '* The initial parameters are defined internally' 
    End If
    Write (messages(5),'(1x,a)') '* The initial values for the parameters are: '//Trim(fit_data%init_param) 
    
    If (fit_data%start_time%fread) Then                           
      Write(tstart, '(f12.3)') fit_data%start_time%value
      If (fit_data%end_time%fread) Then
         Write(tend, '(f12.3)') fit_data%end_time%value
         Write (messages(6),'(1x,a)') '* The user has requested to fit the data between '//Trim(Adjustl(tstart))//&
                                     &' and '//Trim(Adjustl(tend))//' ps'
      Else
         If (Abs(fit_data%start_time%value) > epsilon(1.0_wp) ) Then
           Write (messages(6),'(1x,a)') '* The user has requested to fit data above '//Trim(Adjustl(tstart))//' ps'
         Else
           Write (messages(6),'(1x,a)') '* The fitting will consider all input data'
         End If
      End If
    Else
      If (fit_data%end_time%fread) Then
         Write(tend, '(f12.3)') fit_data%end_time%value
         Write (messages(6),'(1x,a)') '* The user has requested to fit data below '//Trim(Adjustl(tend))//' ps'    
      Else
         Write (messages(6),'(1x,a)') '* The fitting will consider all input data'
      End If
    End If
    Call info(messages, 6)

    If (Trim(fit_data%functional_form) == 'msd_linear') Then
        Write (dim,'(i1)') fit_data%msd%dim
        Write (messages(1),'(1x,a)') '* MSD data corresponds to values computed for the "'//Trim(fit_data%msd%select%type)//&
                                    &'" coordinates. Thus, the dimesion (dim) is equal to '//Trim(adjustl(dim)) 
        Call info(messages, 1)
    End If

    Call refresh_out(files) 
    
  End Subroutine print_fitting_settings
  
  Subroutine check_validity_input_data(fit_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Several checks to make sure the user has provided sensible input data for 
    ! fitting
    !
    ! author    - i.scivetti March 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(fit_type),     Intent(In   ) :: fit_data

    Integer(Kind=wi)   :: i
    Logical            :: time_ok, avg_ok, std_ok
    Character(Len=256)  :: messages(2)
    
    time_ok=.True.
    avg_ok=.True.
    std_ok=.True.
    
    Do i=1, fit_data%num_input_data
      If (fit_data%std_column%fread) Then
       If (fit_data%input_data(i)%std < 0.0_wp .And. std_ok) Then
         Call info('***PROBLEMS*** the column number corresponding to the STD values (value assigned to&
                & the "std_column_number" directive) contains at least one negative value!', 1)
         Call info ('Have you specified the right column number? The run is terminated here. Sorry!', 1)
         Call error_stop(' ')
       End If 
      End If
      If (fit_data%input_data(i)%time < 0.0_wp .And. time_ok) Then
        Call info('***WARNING*** the column number corresponding to the recorded times (value assigned to&
                & the "time_column_number" directive) contains at least one negative value!', 1)      
        time_ok=.False.
      End If
      If (avg_ok) Then
        If (Trim(fit_data%what_to_fit%type)=='ocf') Then
          If ((fit_data%input_data(i)%avg < -1.0_wp) .Or. (fit_data%input_data(i)%avg > 1.0_wp)) Then
            Write(messages(2),'(14x,a)') '!!! values must be between -1.0 and 1.0 !!!'
            avg_ok=.False.
          End If
        Else If (Trim(fit_data%what_to_fit%type)=='spcf' .Or. Trim(fit_data%what_to_fit%type)=='tcf') Then
          If ((fit_data%input_data(i)%avg < 0.0_wp) .Or. (fit_data%input_data(i)%avg > 1.0_wp)) Then
            Write(messages(2),'(14x,a)') '!!! values must be between  0.0 and 1.0 !!!'
            avg_ok=.False.
          End If
        Else If (Trim(fit_data%what_to_fit%type)=='msd') Then
          If (fit_data%input_data(i)%avg < 0.0_wp) Then
            Write(messages(2),'(14x,a)') '!!! values must be positive !!!'
            avg_ok=.False.
          End If        
        End If
      End If
    End Do    

    
    If (.Not. avg_ok) Then
      Write(messages(1),'(a)') '***WARNING*** the column number corresponding to the average values (value assigned to&
                & the "avg_column_number" directive) seems to have incorrect data for the selected option of "what_to_fit"'
      Call info(messages, 2)
    End If
    
    If ((.Not. time_ok) .Or. (.Not. avg_ok)) Then 
      Write(messages(1),'(14x,a)') 'Have you specified the right column number? The run will continue,&
                                  & but fitting issues are expected'
      Call info(messages, 1)
    End If
    
  End Subroutine check_validity_input_data
      
  Subroutine print_fitting_results(fit_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to print the results of fitting
    !
    ! author    - i.scivetti March 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(fit_type),     Intent(InOut) :: fit_data
    
    Character(Len=256)  :: message, header(5)
    Character(Len=256)  :: messages(5), char_units
    Integer(kind=wi)   :: i
    Logical            :: fvalid
    
    fit_data%valid_solutions=0
    
    If (Trim(fit_data%functional_form) == 'ocf_2exp_to_constant') Then
      If (fit_data%fit_superfast_time%stat) Then
        Write(message,'(1x,a)') '(A3, A2 and A1 are dimensionless; T1, T2 and T3 are in units of ps)' 
      Else
        Write(message,'(1x,a)') '(A3, A2 and A1 are dimensionless; T2 and T3 are in units of ps)' 
      End If
    Else If (Trim(fit_data%functional_form) == 'ocf_2exp_to_zero') Then
      If (fit_data%fit_superfast_time%stat) Then
        Write(message,'(1x,a)') '(A2 and A1 are dimensionless; T1, T2 and T3 are in units of ps)' 
      Else
        Write(message,'(1x,a)') '(A2 and A1 are dimensionless; T2 and T3 are in units of ps)' 
      End If    
    Else If (Trim(fit_data%functional_form) == 'ocf_1exp_to_constant') Then        
      If (fit_data%fit_superfast_time%stat) Then
        Write(message,'(1x,a)') '(A2 and A1 are dimensionless; T2 and T1 are in units of ps)'
      Else
        Write(message,'(1x,a)') '(A2 and A1 are dimensionless; T2 is in units of ps)' 
      End If    
    Else If (Trim(fit_data%functional_form) == 'ocf_1exp_to_zero') Then
      If (fit_data%fit_superfast_time%stat) Then
        Write(message,'(1x,a)') '(A1 is dimensionless; T2 and T1 are in units of ps)' 
      Else
        Write(message,'(1x,a)') '(A1 is dimensionless; T2 is in units of ps)' 
      End If    
    Else If (Trim(fit_data%functional_form) == 'tcf_2exp_to_zero') Then
      Write(message,'(1x,a)') '(A0 is dimensionless; A1 and A2 are in units of 1/ps)' 
    Else If (Trim(fit_data%functional_form) == 'tcf_1exp_to_constant') Then
      Write(message,'(1x,a)') '(A0 is dimensionless; A1 is in units of 1/ps)' 
    Else If (Trim(fit_data%functional_form) == 'tcf_1exp_to_zero') Then
      Write(message,'(1x,a)') '(A1 is in units of 1/ps)'
    Else If (Trim(fit_data%functional_form) == 'spcf_2exp_to_zero') Then
      Write(message,'(1x,a)') '(A0 is dimensionless; A1 and A2 are in units of 1/ps)' 
    Else If (Trim(fit_data%functional_form) == 'spcf_1exp_to_constant') Then
      Write(message,'(1x,a)') '(A0 is dimensionless; A1 is in units of 1/ps)' 
    Else If (Trim(fit_data%functional_form) == 'spcf_1exp_to_zero') Then
      Write(message,'(1x,a)') '(A1 is in units of 1/ps)'       
    Else If (Trim(fit_data%functional_form) == 'msd_linear') Then
      Write(message,'(1x,a)') '(A0 is in units of Angstrom^2; A1 is in units of Angstrom^2/ps)' 
    End If
    Call info(message, 1)
    char_units=message

    Call info(' ', 1) 
    Call info(' Table of fitted parameters (only converged fittings)', 1)
    Call info(fit_data%separation_line, 1)
    Write(message,'(1x,a)') Trim(fit_data%header)
    Call info(message, 1)
    Call info(fit_data%separation_line, 1)
    Do i = 1, fit_data%std_samples%value
      If (fit_data%param%fit_ok(i)) Then
        fit_data%xstd=(i-1)*fit_data%delta_std-fit_data%max_times_std%value
        If (Trim(fit_data%functional_form) == 'ocf_2exp_to_constant') Then
          If (fit_data%fit_superfast_time%stat) Then
            Write(message,'(f12.3,8x,e13.4,6x,12(e13.4,1x,a,1x))') fit_data%xstd, fit_data%param%ssr(i), &
                                                & fit_data%param%A3%value(i), Trim(fit_data%param%A3%fit_error(i)),&
                                                & fit_data%param%A2%value(i), Trim(fit_data%param%A2%fit_error(i)),&
                                                & fit_data%param%A1%value(i), Trim(fit_data%param%A1%fit_error(i)),&
                                                & fit_data%param%T3%value(i), Trim(fit_data%param%T3%fit_error(i)),&
                                                & fit_data%param%T2%value(i), Trim(fit_data%param%T2%fit_error(i)),& 
                                                & fit_data%param%T1%value(i), Trim(fit_data%param%T1%fit_error(i))
          Else
            Write(message,'(f12.3,8x,e13.4,6x,10(e13.4,1x,a,1x))') fit_data%xstd, fit_data%param%ssr(i), &
                                                & fit_data%param%A3%value(i), Trim(fit_data%param%A3%fit_error(i)),&
                                                & fit_data%param%A2%value(i), Trim(fit_data%param%A2%fit_error(i)),&
                                                & fit_data%param%A1%value(i), Trim(fit_data%param%A1%fit_error(i)),&
                                                & fit_data%param%T3%value(i), Trim(fit_data%param%T3%fit_error(i)),&
                                                & fit_data%param%T2%value(i), Trim(fit_data%param%T2%fit_error(i))
          End If                                    
        Else If (Trim(fit_data%functional_form) == 'ocf_2exp_to_zero') Then
          If (fit_data%fit_superfast_time%stat) Then
            Write(message,'(f12.3,8x,e13.4,6x,10(e13.4,1x,a,1x))') fit_data%xstd, fit_data%param%ssr(i), &
                                                & fit_data%param%A2%value(i), Trim(fit_data%param%A2%fit_error(i)),&
                                                & fit_data%param%A1%value(i), Trim(fit_data%param%A1%fit_error(i)),&
                                                & fit_data%param%T3%value(i), Trim(fit_data%param%T3%fit_error(i)),&
                                                & fit_data%param%T2%value(i), Trim(fit_data%param%T2%fit_error(i)),&
                                                & fit_data%param%T1%value(i), Trim(fit_data%param%T1%fit_error(i))
          Else
            Write(message,'(f12.3,8x,e13.4,6x,8(e13.4,1x,a,1x))') fit_data%xstd, fit_data%param%ssr(i), &
                                                & fit_data%param%A2%value(i), Trim(fit_data%param%A2%fit_error(i)),&
                                                & fit_data%param%A1%value(i), Trim(fit_data%param%A1%fit_error(i)),&
                                                & fit_data%param%T3%value(i), Trim(fit_data%param%T3%fit_error(i)),&
                                                & fit_data%param%T2%value(i), Trim(fit_data%param%T2%fit_error(i))
          End If                                    
        Else If (Trim(fit_data%functional_form) == 'ocf_1exp_to_constant') Then
          If (fit_data%fit_superfast_time%stat) Then
          Write(message,'(f12.3,8x,e13.4,6x,8(e13.4,1x,a,1x))') fit_data%xstd, fit_data%param%ssr(i), &
                                              & fit_data%param%A2%value(i), Trim(fit_data%param%A2%fit_error(i)),&
                                              & fit_data%param%A1%value(i), Trim(fit_data%param%A1%fit_error(i)),&
                                              & fit_data%param%T2%value(i), Trim(fit_data%param%T2%fit_error(i)),&
                                              & fit_data%param%T1%value(i), Trim(fit_data%param%T1%fit_error(i))
          Else
          Write(message,'(f12.3,8x,e13.4,6x,6(e13.4,1x,a,1x))') fit_data%xstd, fit_data%param%ssr(i), &
                                              & fit_data%param%A2%value(i), Trim(fit_data%param%A2%fit_error(i)),&
                                              & fit_data%param%A1%value(i), Trim(fit_data%param%A1%fit_error(i)),&
                                              & fit_data%param%T2%value(i), Trim(fit_data%param%T2%fit_error(i))
          End If                                    
        Else If (Trim(fit_data%functional_form) == 'ocf_1exp_to_zero') Then
          If (fit_data%fit_superfast_time%stat) Then
            Write(message,'(f12.3,8x,e13.4,6x,6(e13.4,1x,a,1x))') fit_data%xstd, fit_data%param%ssr(i), &
                                                & fit_data%param%A1%value(i), Trim(fit_data%param%A1%fit_error(i)),&
                                                & fit_data%param%T2%value(i), Trim(fit_data%param%T2%fit_error(i)),&
                                                & fit_data%param%T1%value(i), Trim(fit_data%param%T1%fit_error(i))
          Else
            Write(message,'(f12.3,8x,e13.4,6x,4(e13.4,1x,a,1x))') fit_data%xstd, fit_data%param%ssr(i), &
                                                & fit_data%param%A1%value(i), Trim(fit_data%param%A1%fit_error(i)),&
                                                & fit_data%param%T2%value(i), Trim(fit_data%param%T2%fit_error(i))
          End If                                    
        Else If (Trim(fit_data%functional_form) == 'tcf_2exp_to_zero') Then
          Write(message,'(f12.3,8x,e13.4,6x,6(e13.4,1x,a,1x))') fit_data%xstd, fit_data%param%ssr(i), &
                                              & fit_data%param%A0%value(i), Trim(fit_data%param%A0%fit_error(i)),&
                                              & fit_data%param%A1%value(i), Trim(fit_data%param%A1%fit_error(i)),&
                                              & fit_data%param%A2%value(i), Trim(fit_data%param%A2%fit_error(i))
        Else If (Trim(fit_data%functional_form) == 'tcf_1exp_to_constant') Then
          Write(message,'(f12.3,8x,e13.4,6x,4(e13.4,1x,a,1x))') fit_data%xstd, fit_data%param%ssr(i), &
                                              & fit_data%param%A0%value(i), Trim(fit_data%param%A0%fit_error(i)),&
                                              & fit_data%param%A1%value(i), Trim(fit_data%param%A1%fit_error(i))
        Else If (Trim(fit_data%functional_form) == 'tcf_1exp_to_zero') Then
          Write(message,'(f12.3,8x,e13.4,6x,2(e13.4,1x,a,1x))') fit_data%xstd, fit_data%param%ssr(i), &
                                              & fit_data%param%A1%value(i), Trim(fit_data%param%A1%fit_error(i))
        Else If (Trim(fit_data%functional_form) == 'spcf_2exp_to_zero') Then
          Write(message,'(f12.3,8x,e13.4,6x,6(e13.4,1x,a,1x))') fit_data%xstd, fit_data%param%ssr(i), &
                                              & fit_data%param%A0%value(i), Trim(fit_data%param%A0%fit_error(i)),&
                                              & fit_data%param%A1%value(i), Trim(fit_data%param%A1%fit_error(i)),&
                                              & fit_data%param%A2%value(i), Trim(fit_data%param%A2%fit_error(i))
        Else If (Trim(fit_data%functional_form) == 'spcf_1exp_to_constant') Then
          Write(message,'(f12.3,8x,e13.4,6x,4(e13.4,1x,a,1x))') fit_data%xstd, fit_data%param%ssr(i), &
                                              & fit_data%param%A0%value(i), Trim(fit_data%param%A0%fit_error(i)),&
                                              & fit_data%param%A1%value(i), Trim(fit_data%param%A1%fit_error(i))
        Else If (Trim(fit_data%functional_form) == 'spcf_1exp_to_zero') Then
          Write(message,'(f12.3,8x,e13.4,6x,2(e13.4,1x,a,1x))') fit_data%xstd, fit_data%param%ssr(i), &
                                              & fit_data%param%A1%value(i), Trim(fit_data%param%A1%fit_error(i))
        Else If (Trim(fit_data%functional_form) == 'msd_linear') Then
          Write(message,'(f12.3,8x,e13.4,6x,4(e13.4,1x,a,1x))') fit_data%xstd, fit_data%param%ssr(i), &
                                              & fit_data%param%A0%value(i), Trim(fit_data%param%A0%fit_error(i)),&
                                              & fit_data%param%A1%value(i), Trim(fit_data%param%A1%fit_error(i))
        End If
        Call info(message, 1)
      End If
    End Do
    Call info(fit_data%separation_line, 1)
    Call info(' (*)SSR/NDF: Sum of Square of Residuals/NDF', 1)
    
    ! Inform is the fitting is not valid 
    Do i = 1, fit_data%std_samples%value
      If (fit_data%param%fit_ok(i)) Then
        If (.Not. fit_data%param%valid(i)) Then
          fit_data%xstd=(i-1)*fit_data%delta_std-fit_data%max_times_std%value
          Write(fit_data%xstd_char,'(f6.3)') fit_data%xstd
          If (Abs(fit_data%xstd) < epsilon(1.0_wp)) Then
           Write(message, '(a)') 'EXTRA WARNING! Fitting for the average correlation is physically invalid! ACTION:&
                                & discarded from the analysis'
          
          Else 
            Write(message, '(a)') 'EXTRA WARNING! Fitting of correlation corresponding&
                                & to a deviation of "'//Trim(Adjustl(fit_data%xstd_char))//'*STD"&
                                & is physically invalid. ACTION: discarded from the analysis'
          End If
          Call info(message, 1)
        Else
          fit_data%valid_solutions=fit_data%valid_solutions+1
        End If  
      End If
    End Do
    
    If (fit_data%valid_solutions == 0) Then
       Call info (' ', 1) 
       If (fit_data%std_samples%value /= 1) Then
         Write(messages(1), '(a)') 'SERIOUS PROBLEMS: None of the fittings above corresponds to a physically&
                        & invalid solution. No meaningful results can be derived. The analysis stops here.'
         Write(messages(2), '(a)') 'The user should review the input data (plus the "column_number" directives),&
                            & the choice for "fitting_function", the time domain for fitting and the value of the&
                            & "max_times_std" directive.'
       Else
         Write(messages(1), '(a)') 'SERIOUS PROBLEMS: The fitting above corresponds to a physically&
                        & invalid solution. No meaningful result can be derived. The analysis stops here.'
         Write(messages(2), '(a)') 'The user should review the input data (plus the "column_number" directives),&
                            & the choice for "fitting_function" and the time domain for fitting.'
       End If
       Call info(messages, 2)
        If (fit_data%input_parameters%fread) Then
         Write(message, '(a)') '*** ADVISE: either i) revise the values of the "&input_parameters" block or&
                              & ii) try removing the "&input_parameters" block (if this has not been tried already).'
         Call info (message, 1) 
        End If 
        If (fit_data%fit_superfast_time%stat) Then
         Write(message, '(a)') '*** ADVISE: inspect if the use of "fit_superfast_time" is correct in this case'
         Call info (message, 1) 
        End If         
       Call error_stop(' ')
    Else
      If (fit_data%valid_solutions /= fit_data%std_samples%value) Then
         Write(message, '(a)') 'The user should review the input data, choice of the functional form&
                            & for fitting and the value of the "max_times_std" directive.'
         Call info(message, 1)
         If (fit_data%input_parameters%fread) Then
           Write(message, '(a)') '*** ADVISE either i) revise the values of the "&input_parameters" block or&
                               & ii) try removing the "&input_parameters" block.'
           Call info (message, 1) 
         End If 
      End If
      If (fit_data%print_average_ssr%fread) Then
         Call info(' ', 1)
         Write(message, '(1x,a,e13.4)') 'Average SSR/NDF (only converged fittings)=', fit_data%average_ssr 
         Call info (message, 1)         
      End If
      Call info(' ', 1)
      Call info(' ', 1)
      If (fit_data%std_samples%value /=1) Then
        Write(header(1), '(1x,a)') 'MEAN value and STD for each parameter'
        Write(header(3), '(a)') ' ----------------------------------'
        Write(header(4), '(1x,a,5x,a,9x,a)') 'Parameter', 'MEAN', 'STD'
        Write(header(5), '(a)') ' ----------------------------------'
      Else
        Write(header(1), '(1x,a)') 'MEAN value for each parameter (no STD for this analysis)'
        Write(header(3), '(a)') ' ---------------------'
        Write(header(4), '(1x,a,5x,a)') 'Parameter', 'MEAN'
        Write(header(5), '(a)') ' ---------------------'
      End If
      Write(header(2), '(a)') Trim(char_units)
      Call info(header, 5)
      
      
      If (Trim(fit_data%functional_form) == 'ocf_2exp_to_constant') Then
        Call print_parameter(messages(1), fit_data%std_samples%value, 'A3', fit_data%param%A3%mean, fit_data%param%A3%std)
        Call print_parameter(messages(2), fit_data%std_samples%value, 'A2', fit_data%param%A2%mean, fit_data%param%A2%std)
        Call print_parameter(messages(3), fit_data%std_samples%value, 'A1', fit_data%param%A1%mean, fit_data%param%A1%std)
        Call print_parameter(messages(4), fit_data%std_samples%value, 'T3', fit_data%param%T3%mean, fit_data%param%T3%std)
        Call print_parameter(messages(5), fit_data%std_samples%value, 'T2', fit_data%param%T2%mean, fit_data%param%T2%std)
        Call info(messages, 5)
      Else If (Trim(fit_data%functional_form) == 'ocf_2exp_to_zero') Then
        Call print_parameter(messages(1), fit_data%std_samples%value, 'A2', fit_data%param%A2%mean, fit_data%param%A2%std)
        Call print_parameter(messages(2), fit_data%std_samples%value, 'A1', fit_data%param%A1%mean, fit_data%param%A1%std)
        Call print_parameter(messages(3), fit_data%std_samples%value, 'T3', fit_data%param%T3%mean, fit_data%param%T3%std)
        Call print_parameter(messages(4), fit_data%std_samples%value, 'T2', fit_data%param%T2%mean, fit_data%param%T2%std)      
        Call info(messages, 4)
      Else If (Trim(fit_data%functional_form) == 'ocf_1exp_to_constant') Then
        Call print_parameter(messages(1), fit_data%std_samples%value, 'A2', fit_data%param%A2%mean, fit_data%param%A2%std)
        Call print_parameter(messages(2), fit_data%std_samples%value, 'A1', fit_data%param%A1%mean, fit_data%param%A1%std)
        Call print_parameter(messages(3), fit_data%std_samples%value, 'T2', fit_data%param%T2%mean, fit_data%param%T2%std)
        Call info(messages, 3)
      Else If (Trim(fit_data%functional_form) == 'ocf_1exp_to_zero') Then
        Call print_parameter(messages(1), fit_data%std_samples%value, 'A1', fit_data%param%A1%mean, fit_data%param%A1%std)
        Call print_parameter(messages(2), fit_data%std_samples%value, 'T2', fit_data%param%T2%mean, fit_data%param%T2%std)
        Call info(messages, 2)
      Else If (Trim(fit_data%functional_form) == 'tcf_2exp_to_zero' .Or.&
               Trim(fit_data%functional_form) == 'spcf_2exp_to_zero') Then
        Call print_parameter(messages(1), fit_data%std_samples%value, 'A0', fit_data%param%A0%mean, fit_data%param%A0%std)
        Call print_parameter(messages(2), fit_data%std_samples%value, 'A1', fit_data%param%A1%mean, fit_data%param%A1%std)
        Call print_parameter(messages(3), fit_data%std_samples%value, 'A2', fit_data%param%A2%mean, fit_data%param%A2%std)
        Call info(messages, 3)
      Else If (Trim(fit_data%functional_form) == 'tcf_1exp_to_constant' .Or. &
               Trim(fit_data%functional_form) == 'spcf_1exp_to_constant') Then  
        Call print_parameter(messages(1), fit_data%std_samples%value, 'A0', fit_data%param%A0%mean, fit_data%param%A0%std)
        Call print_parameter(messages(2), fit_data%std_samples%value, 'A1', fit_data%param%A1%mean, fit_data%param%A1%std)
        Call info(messages, 2)
      Else If (Trim(fit_data%functional_form) == 'tcf_1exp_to_zero' .Or. &
               Trim(fit_data%functional_form) == 'spcf_1exp_to_zero') Then
        Call print_parameter(messages(1), fit_data%std_samples%value, 'A1', fit_data%param%A1%mean, fit_data%param%A1%std)
        Call info(messages, 1)
      Else If (Trim(fit_data%functional_form) == 'msd_linear') Then  
        Call print_parameter(messages(1), fit_data%std_samples%value, 'A0', fit_data%param%A0%mean, fit_data%param%A0%std)
        Call print_parameter(messages(2), fit_data%std_samples%value, 'A1', fit_data%param%A1%mean, fit_data%param%A1%std)
        Call info(messages, 2)
      End If

      If (Trim(fit_data%what_to_fit%type)=='ocf') Then
        If (fit_data%fit_superfast_time%stat) Then
          Call print_parameter(messages(1), fit_data%std_samples%value, 'T1', fit_data%param%T1%mean, fit_data%param%T1%std)
          Call info(messages, 1)          
        End If
      End If
     
      Call check_validity_mean_values(fit_data, fvalid)
      
      If(.Not. fvalid) Then
       If (fit_data%input_parameters%fread) Then
          If (Trim(fit_data%functional_form) == 'ocf_1exp_to_zero'     .Or. &
              Trim(fit_data%functional_form) == 'ocf_1exp_to_constant' .Or. &
              Trim(fit_data%functional_form) == 'ocf_2exp_to_zero'     .Or. &
              Trim(fit_data%functional_form) == 'ocf_2exp_to_constant') Then
            Write(message, '(1x,a)') 'CAUTION: the user should corroborate the physical validity of the computed&
                         & MEAN values (and STD), and check initial parameters&
                         & set in the "&input_parameters" block'
            Call info(message, 1)
          End If
        Else   
          Write(message, '(a)') 'CAUTION: the user should corroborate the physical validity of the computed&
                                 & MEAN values (and STD) for the parameters.'
          Call info(message, 1)          
        End If
      End If
      
    End If
    
  End Subroutine print_fitting_results

 
End Module fitting
