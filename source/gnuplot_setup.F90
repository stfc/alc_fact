!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Module to manage files for GNUPLOT operation
!
! Copyright - 2026 Ada Lovelace Centre (ALC)
!             Scientific Computing Department (SCD)
!             The Science and Technology Facilities Council (STFC)
!
! Author:        i.scivetti  March 2026
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Module gnuplot_setup 

  Use constants,     Only : min_gnuplot_version,&
                            min_gnuplot_subversion

  Use fileset,       Only : file_type, &
                            FILE_GNUPLOT_INPUT, &
                            FILE_GNUPLOT_OUTPUT, &
                            FILE_FIT_PARAMS,&  
                            FILE_DUMP, &
                            FILE_INPUT_DATA

  Use fitting_setup, Only : fit_type
  
  Use numprec,       Only : wi, wp
  
  Use process_data,   Only : remove_symbols,&
                             check_for_symbols
  
  Use unit_output,   Only : error_stop, &
                            info

  Implicit None

  Private
 
  Public :: delete_gnuplot_working_files, generate_gnuplot_input
  Public :: check_gnuplot_availability, extract_data_from_file
  Public :: define_gnuplot_actions
  
Contains

  Subroutine check_gnuplot_availability(files)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Check if gnuplot is available
    !
    ! author    - i.scivetti March 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),    Intent(InOut) :: files(:)
    
    Character(len=256) :: string, messages(2)
    Character(Len=256) :: name, version_char
    Integer(Kind=wi)   :: io, version, sub_version
    Logical            :: flag
     
    Open(Newunit=files(FILE_DUMP)%unit_no, File=files(FILE_DUMP)%filename, Status='Replace')
    Call execute_command_line('which gnuplot > '//Trim(files(FILE_DUMP)%filename))
    Read (files(FILE_DUMP)%unit_no, Fmt='(a)', iostat=io) string
    Close(files(FILE_DUMP)%unit_no, Status='Delete')

    Write(messages(1), '(1x, a)') 'SERIOUS PROBLEM: this code uses the GNUPLOT software&
                                   & to execute the fittings.'    
    
    If (is_iostat_end(io)) Then
      Write(messages(2), '(1x, a)') 'GNUPLOT does not seem to be installed/loaded.&
                                   & Please make GNUPLOT available in this machine'
      Call info(' ', 1)
      Call info(messages, 2)
      Call error_stop(' ')
    Else
      flag=.True.
      Open(Newunit=files(FILE_DUMP)%unit_no, File=files(FILE_DUMP)%filename, Status='Replace')
      Call execute_command_line('gnuplot --version > '//Trim(files(FILE_DUMP)%filename))
      Read (files(FILE_DUMP)%unit_no, Fmt=*, iostat=io) name, version_char
      Call remove_symbols(version_char,'.')
      Read (version_char, Fmt=*, iostat=io) version, sub_version
      Close(files(FILE_DUMP)%unit_no, Status='Delete')
      If (version<min_gnuplot_version) Then
        flag=.False.
      Else
        If (version<min_gnuplot_subversion) Then
          flag=.False.
        End If
      End If
      If (.Not. flag) Then
        Write(messages(2), '(1x, a)') 'However, the minimum version of GNUPLOT required is "5.2".&
                                    & Please update your installed version of GNUPLOT.'
        Call info(' ', 1)
        Call info(messages, 2)
        Call error_stop(' ')      
      End If
    End if
    
  
  End Subroutine check_gnuplot_availability

  Subroutine delete_gnuplot_working_files(files)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Clean warking files
    !
    ! author    - i.scivetti March 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),    Intent(InOut) :: files(:)

    Character(Len=256) :: temp_file(3)
    Integer(kind=wi)   :: i
    Logical            :: file_exists

    temp_file(1)=Trim(files(FILE_FIT_PARAMS)%filename)
    temp_file(2)=Trim(files(FILE_GNUPLOT_INPUT)%filename) 
    temp_file(3)=Trim(files(FILE_GNUPLOT_OUTPUT)%filename)
    
    Do i=1, 3
      Inquire(File=Trim(temp_file(i)), EXIST=file_exists)
      If (file_exists) Then
        Call execute_command_line('rm '//Trim(temp_file(i)))
      End If
    End Do
    
  End Subroutine delete_gnuplot_working_files

  Subroutine extract_data_from_file(files, fit_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Extract data from the generated files
    !
    ! author    - i.scivetti March 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),    Intent(InOut) :: files(:)
    Type(fit_type),     Intent(InOut) :: fit_data

    Integer(Kind=wi)   :: iunit, i
    Character(Len=256) :: timecol, avgcol, stdcol, is_rubbish
    Logical            :: found

    ! print column numbers to strings
    Write(timecol,*) fit_data%time_column%value
    Write(avgcol,*)  fit_data%avg_column%value
    If (fit_data%std_column%fread) Then
      Write(stdcol,*)  fit_data%std_column%value
    End If
    
    Open(Newunit=files(FILE_GNUPLOT_INPUT)%unit_no, File=files(FILE_GNUPLOT_INPUT)%filename, Status='Replace')
    iunit=files(FILE_GNUPLOT_INPUT)%unit_no
    
    Write(iunit, '(a)') 'FILE = "'//Trim(fit_data%filename%type)//'"' 
    Write(iunit, '(a)') 'stats FILE u 0 nooutput'
    Write(iunit, '(a)') 'rowCount = STATS_records'
    Write(iunit, '(a)') 'array time[rowCount]; array AVG[rowCount]; array STD[rowCount]'
    Write(iunit, '(a)') 'stats FILE u (time[int(rowCount-$0)] = $'//Trim(Adjustl(timecol))//') nooutput'
    Write(iunit, '(a)') 'stats FILE u (AVG[int(rowCount-$0)]  = $'//Trim(Adjustl(avgcol))//') nooutput'
    If (fit_data%std_column%fread) Then
      Write(iunit, '(a)') 'stats FILE u (STD[int(rowCount-$0)]  = $'//Trim(Adjustl(stdcol))//') nooutput'
    End If
    Write(iunit, '(a)') 'set print "'//Trim(files(FILE_INPUT_DATA)%filename)//'"'
    Write(iunit, '(a)') 'print rowCount'
    If (fit_data%std_column%fread) Then
      Write(iunit, '(a)') 'do for [i=rowCount:1:-1] {print time[i], AVG[i], STD[i]}'
    Else
      Write(iunit, '(a)') 'do for [i=rowCount:1:-1] {print time[i], AVG[i]}'
    End If
    Close(iunit)
    Call execute_command_line(Trim(fit_data%exec_gnuplot))
    
    Call execute_command_line('rm '//Trim(files(FILE_GNUPLOT_INPUT)%filename))

    Open(Newunit=files(FILE_INPUT_DATA)%unit_no, File=files(FILE_INPUT_DATA)%filename, Status='Old')
    iunit=files(FILE_INPUT_DATA)%unit_no
    Read (iunit,*) fit_data%num_input_data
    found=.False.
    
    Read (iunit,'(a)') is_rubbish
    Call check_for_symbols(is_rubbish, 'undefined', found)
    If (found) Then
     Call info('***PROBLEMS*** Wrong definition of the "column_number" directives!&
             & Please check the number of columns of the data files provided.', 1)
     Close(iunit, Status='Delete') 
     Call error_stop(' ')
    Else
     Backspace iunit
    End If

    Call fit_data%alloc_data()

    Do i=1, fit_data%num_input_data
      Read (iunit,'(a)') is_rubbish
      Call check_for_symbols(is_rubbish, 'NaN', found)
      If (found) Then
       Call info('***PROBLEMS*** Non-numerical value(s) found as part of the relevant input columns!&
               & Please check the content of the input file "'//Trim(fit_data%filename%type)//'"', 1)
       Call error_stop(' ')
      Else
       Backspace iunit
      End If
      If (fit_data%std_column%fread) Then
       Read(iunit,*) fit_data%input_data(i)%time, fit_data%input_data(i)%avg, fit_data%input_data(i)%std
       If (Trim(fit_data%what_to_fit%type)=='msd') Then
         fit_data%input_data(i)%avg =  fit_data%msd%factor * fit_data%input_data(i)%avg 
         fit_data%input_data(i)%std =  fit_data%msd%factor * fit_data%input_data(i)%std
        End If
      Else
        Read(iunit,*) fit_data%input_data(i)%time, fit_data%input_data(i)%avg
       If (Trim(fit_data%what_to_fit%type)=='msd') Then
         fit_data%input_data(i)%avg =  fit_data%msd%factor * fit_data%input_data(i)%avg 
        End If        
      End If
    End Do
    Close(iunit, Status='Delete')
    
  End Subroutine extract_data_from_file
  
  Subroutine define_gnuplot_actions(files, fit_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Setting gnuplot instructions for fitting
    !
    ! author    - i.scivetti March 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),    Intent(InOut) :: files(:)
    Type(fit_type),     Intent(InOut) :: fit_data
    
    Character(Len=256) :: num_it
    Character(Len=256) :: tstart, tend, whattofit
    Character(Len=256) :: timecol, avgcol, stdcol
    
    Write(num_it, '(i6)') fit_data%max_iterations%value
    fit_data%dump_params=Trim(files(FILE_GNUPLOT_OUTPUT)%filename)//' > '//Trim(files(FILE_FIT_PARAMS)%filename)
    fit_data%set_fit='set fit  logfile "'//Trim(files(FILE_GNUPLOT_OUTPUT)%filename)//'" quiet maxiter '//Trim(Adjustl(num_it))
    fit_data%exec_gnuplot='gnuplot '//Trim(files(FILE_GNUPLOT_INPUT)%filename)
    fit_data%check_fit='grep "Fit stopped\|BREAK" '//Trim(files(FILE_GNUPLOT_OUTPUT)%filename)//&
                   &' > '//Trim(files(FILE_DUMP)%filename)
    fit_data%check_ndf='grep "degrees of freedom" '//Trim(files(FILE_GNUPLOT_OUTPUT)%filename)//&
                   &' > '//Trim(files(FILE_DUMP)%filename)
    fit_data%check_window='grep "No data to fit" '//Trim(files(FILE_GNUPLOT_OUTPUT)%filename)//&
                   &' > '//Trim(files(FILE_DUMP)%filename)

    fit_data%ssr_string='final sum of squares of residuals' 
    fit_data%results_string=Trim(fit_data%ssr_string)//Trim(fit_data%param_list)
    fit_data%extract_results='grep -e "'//Trim(fit_data%results_string)//'" '//Trim(fit_data%dump_params)
    
    ! print column numbers to strings
    Write(timecol,'(i3)') fit_data%time_column%value
    Write(avgcol,'(i3)')  fit_data%avg_column%value
    If (fit_data%std_column%fread) Then
      Write(stdcol,'(i3)')  fit_data%std_column%value
    End If
    
    If (fit_data%std_samples%value == 1) Then
      whattofit='($'//Trim(Adjustl(avgcol))//')'
    Else
      whattofit='($'//Trim(Adjustl(avgcol))//'+$'//Trim(Adjustl(stdcol))//'*xstd)'
    End If
    
    If (fit_data%start_time%fread) Then
      Write(tstart,'(f12.3)') fit_data%start_time%value
      If (fit_data%end_time%fread) Then
        Write(tend,'(f12.3)') fit_data%end_time%value
        fit_data%instruct_fitting="fit ["//Trim(Adjustl(tstart))//":"//Trim(Adjustl(tend))//"] f(x) '"&
                                &//Trim(fit_data%filename%type)//&
                                &"' using "//Trim(Adjustl(timecol))//":"//Trim(whattofit)//" "//Trim(fit_data%via_params)
      Else
      fit_data%instruct_fitting="fit ["//Trim(Adjustl(tstart))//":] f(x) '"//Trim(fit_data%filename%type)//&
                                &"' using "//Trim(Adjustl(timecol))//":"//Trim(whattofit)//" "//Trim(fit_data%via_params)
      End If
    Else
      If (fit_data%end_time%fread) Then
        Write(tend,'(f12.3)') fit_data%end_time%value
        fit_data%instruct_fitting="fit [:"//Trim(Adjustl(tend))//"] f(x) '"&
                                &//Trim(fit_data%filename%type)//&
                                &"' using "//Trim(Adjustl(timecol))//":"//Trim(whattofit)//" "//Trim(fit_data%via_params)
      Else
        fit_data%instruct_fitting="fit [:] f(x) '"//Trim(fit_data%filename%type)//&
                                &"' using "//Trim(Adjustl(timecol))//":"//Trim(whattofit)//" "//Trim(fit_data%via_params)
      End If
    End If
    
    fit_data%rm_gnuplot_output='rm '//Trim(files(FILE_GNUPLOT_OUTPUT)%filename)
    fit_data%rm_gnuplot_input='rm '//Trim(files(FILE_GNUPLOT_INPUT)%filename)
    
  End Subroutine define_gnuplot_actions

  Subroutine generate_gnuplot_input(files, fit_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Print input file for GNUPLOT fitting
    !
    ! author    - i.scivetti March 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),    Intent(InOut) :: files(:)
    Type(fit_type),     Intent(InOut) :: fit_data

    Integer(Kind=wi)    :: iunit
    
    Open(Newunit=files(FILE_GNUPLOT_INPUT)%unit_no, File=files(FILE_GNUPLOT_INPUT)%filename, Status='Replace')
    iunit=files(FILE_GNUPLOT_INPUT)%unit_no
    Write(iunit,'(a)') Trim(fit_data%formula)
    Write(iunit,'(a)') Trim(fit_data%set_fit)
    Write(iunit,'(a)') Trim(fit_data%init_param)
    Write(iunit,'(a)') 'xstd='//Trim(Adjustl(fit_data%xstd_char))
    Write(iunit,'(a)') Trim(fit_data%instruct_fitting)
    Close(iunit)
  
  End Subroutine generate_gnuplot_input  
  
End module gnuplot_setup
