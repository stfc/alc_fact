!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Module thar reads input data from the SETTINGS file
!
! Author        - i.scivetti   March 2026
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Module input_settings 

  Use fileset,           Only : file_type, &
                                FILE_SET, &  
                                FILE_OUT, & 
                                refresh_out
  Use numprec,           Only : wi, &
                                wp
  Use process_data,      Only : capital_to_lower_case, &
                                check_for_rubbish, &
                                get_word_length
  Use unit_output,       Only : error_stop,&
                                info
  Use fitting_setup,     Only : fit_type
  
  Use fitting_checks,    Only : check_fitting_settings  

  Implicit None
  
  Private
  Public :: read_settings

Contains

  Subroutine duplication_error(directive)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Aborts execution when duplication for
   ! a directive is found
   !
   ! author - i. scivetti  March 2026
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Character(Len=*), Intent(In   ) :: directive

    Character(Len=256)  :: message

    Write (message,'(4a)') '***ERROR - Directive "', Trim(directive), '" is duplicated!'
    Call error_stop(message)

  End Subroutine duplication_error  

  Subroutine set_read_status(word, io, fread, fail, string)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to:
    !  - prevent duplication
    !  - define input directive is read by setting fread=.True. 
    !  - test if there was a problem with reading a directive, indicated by io/=0. This sets fail=.True.
    !
    ! author    - i.scivetti March 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Character(Len=*), Intent(In   ) :: word
    Integer(Kind=wi), Intent(In   ) :: io
    Logical,          Intent(  Out) :: fread 
    Logical,          Intent(InOut) :: fail
    Character(Len=*), Optional, Intent(InOut) :: string

    If (fread)then
      Call duplication_error(word)
    Else
      fread=.True.
      If (io /= 0) Then
        fail=.True.
      End If
    End If

    If (present(string)) then
      Call capital_to_lower_case(string)
    End If

  End Subroutine set_read_status 

  Subroutine read_settings(files, fit_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to read settings from SETTINGS file.
    ! Lines starting with # are ignored and assumed as comments. 
    ! If a directive is identified during the reading of the file, subroutine "set_read_fail" 
    ! assigns fread=.True. On the contrary, the subroutine assigns fail=.True. (fail=.False.) 
    ! if the format/syntax for the directive is correct (incorrect)
    ! If the directive is repeated the execution is aborted via subroutine duplication 
    ! 
    ! author        - i.scivetti March 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),    Intent(InOut) :: files(:)
    Type(fit_type),     Intent(InOut) :: fit_data 
 
    Logical            :: safe
    Character(Len=256) :: word
    Integer(Kind=wi)   :: length, io, iunit
  
    Character(Len=256)  :: message

    Character(Len=32 )  :: set_file
    Character(Len=32 )  :: set_error

    set_file = Trim(files(FILE_SET)%filename)
    set_error = '***ERROR in the '//Trim(set_file)//' file.'

    ! Open the SETTINGS file with settings
    Inquire(File=files(FILE_SET)%filename, Exist=safe)
    
    If (.not.safe) Then
      Call info(' ', 1)
      Write (message,'(4(1x,a))') Trim(set_error), 'File', Trim(set_file), '(settings for analysis) not found'
      Call error_stop(message)
    Else
      Open(Newunit=files(FILE_SET)%unit_no, File=Trim(set_file), Status='old')
      iunit=files(FILE_SET)%unit_no 
    End If

     Read (iunit, Fmt=*, iostat=io) word
     ! If nothing is found, complain and abort
     If (is_iostat_end(io)) Then
       Write (message,'(3(1x,a))') Trim(set_error), Trim(set_file), 'file seems to be empty?. Please check'
       Call error_stop(message)
     End If
     ! Check header has "#" as the first character 
     If (word(1:1)/='#') Then
       Write (message,'(4(1x,a))') Trim(set_error), 'Heading comment in file', Trim(set_file), & 
                                  'is required and MUST be preceded with the symbol "#"'
       Call error_stop(message)
     End If

    Do
      Read (iunit, Fmt=*, iostat=io) word
      If (io /= 0) Then
        Exit
      end If
      Call check_for_rubbish(iunit, Trim(set_file)) 
      Call get_word_length(word,length)
      Call capital_to_lower_case(word)

      If (word(1:1) == '#' .Or. word(1:3) == '   ') Then
        ! Do nothing if line is a comment of we have an empty line
        Read (iunit, Fmt=*, iostat=io) word

      Else If (word(1:length) == 'filename') Then
         Read (iunit, Fmt=*, iostat=io) word, fit_data%filename%type
         Call set_read_status(word, io, fit_data%filename%fread, fit_data%filename%fail)

      Else If (word(1:length) == 'fitting_function') Then
         Read (iunit, Fmt=*, iostat=io) word, fit_data%fitting_function%type
         Call set_read_status(word, io, fit_data%fitting_function%fread, fit_data%fitting_function%fail,&
                            & fit_data%fitting_function%type)         
         
      Else If (word(1:length) == 'what_to_fit') Then
         Read (iunit, Fmt=*, iostat=io) word, fit_data%what_to_fit%type
         Call set_read_status(word, io, fit_data%what_to_fit%fread, fit_data%what_to_fit%fail,&
                            & fit_data%what_to_fit%type)

      Else If (word(1:length) == 'species_name') Then
         Read (iunit, Fmt=*, iostat=io) word, fit_data%species_name%type
         Call set_read_status(word, io, fit_data%species_name%fread, fit_data%species_name%fail)
                                     
      Else If (Trim(word)=='start_time') Then
        Read (iunit, Fmt=*, iostat=io) fit_data%start_time%tag, fit_data%start_time%value, fit_data%start_time%units
        Call set_read_status(word, io, fit_data%start_time%fread, fit_data%start_time%fail)
         
      Else If (Trim(word)=='end_time') Then
        Read (iunit, Fmt=*, iostat=io) fit_data%end_time%tag, fit_data%end_time%value, fit_data%end_time%units
        Call set_read_status(word, io, fit_data%end_time%fread, fit_data%end_time%fail)
        
      Else If (word(1:length) == 'max_times_std') Then 
        Read (iunit, Fmt=*, iostat=io) word, fit_data%max_times_std%value 
        Call set_read_status(word, io, fit_data%max_times_std%fread, fit_data%max_times_std%fail)

      Else If (word(1:length) == 'std_samples') Then 
        Read (iunit, Fmt=*, iostat=io) word, fit_data%std_samples%value 
        Call set_read_status(word, io, fit_data%std_samples%fread, fit_data%std_samples%fail)

      Else If (word(1:length) == 'max_iterations') Then 
        Read (iunit, Fmt=*, iostat=io) word, fit_data%max_iterations%value 
        Call set_read_status(word, io, fit_data%max_iterations%fread, fit_data%max_iterations%fail)

      Else If (word(1:length) == 'reactive_species') Then
        Read (iunit, Fmt=*, iostat=io) word, fit_data%reactive_species%stat
        Call set_read_status(word, io, fit_data%reactive_species%fread, fit_data%reactive_species%fail)

      Else If (word(1:length) == 'plot_fittings') Then
        Read (iunit, Fmt=*, iostat=io) word, fit_data%plot_fittings%stat
        Call set_read_status(word, io, fit_data%plot_fittings%fread, fit_data%plot_fittings%fail)

      Else If (word(1:length) == 'plot_raw_data') Then
        Read (iunit, Fmt=*, iostat=io) word, fit_data%plot_raw_data%stat
        Call set_read_status(word, io, fit_data%plot_raw_data%fread, fit_data%plot_raw_data%fail)

      Else If (word(1:length) == 'time_column_number') Then 
        Read (iunit, Fmt=*, iostat=io) word, fit_data%time_column%value 
        Call set_read_status(word, io, fit_data%time_column%fread, fit_data%time_column%fail)
        
      Else If (word(1:length) == 'avg_column_number') Then 
        Read (iunit, Fmt=*, iostat=io) word, fit_data%avg_column%value 
        Call set_read_status(word, io, fit_data%avg_column%fread, fit_data%avg_column%fail)

      Else If (word(1:length) == 'std_column_number') Then 
        Read (iunit, Fmt=*, iostat=io) word, fit_data%std_column%value 
        Call set_read_status(word, io, fit_data%std_column%fread, fit_data%std_column%fail)

      Else If (word(1:length) == '&msd') Then
        Read (iunit, Fmt=*, iostat=io) fit_data%msd%invoke%type
        Call set_read_status(word, io, fit_data%msd%invoke%fread, fit_data%msd%invoke%fail)
        !Read information inside the block
        Call read_msd(iunit, fit_data)

      Else If (word(1:length) == '&input_parameters') Then
        Read (iunit, Fmt=*, iostat=io) fit_data%input_parameters%type
        Call set_read_status(word, io, fit_data%input_parameters%fread, fit_data%input_parameters%fail)
        !Read information inside the block
        Call read_input_parameters(iunit, fit_data)
        
      Else If (word(1:length) == 'fit_superfast_time') Then
        Read (iunit, Fmt=*, iostat=io) word, fit_data%fit_superfast_time%stat
        Call set_read_status(word, io, fit_data%fit_superfast_time%fread, fit_data%fit_superfast_time%fail) 
        
      Else If (word(1:length) == 'print_average_ssr') Then
        Read (iunit, Fmt=*, iostat=io) word, fit_data%print_average_ssr%stat
        Call set_read_status(word, io, fit_data%print_average_ssr%fread, fit_data%print_average_ssr%fail)
        
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
      ! Directive not recognised. Inform and kill 
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
      Else
        If (word(1:1)=='&') Then
          Write (message,'(1x,4a)') Trim(set_error), ' Unknown directive found: "', Trim(word),&
                                  &'. Do you use "&" to define a block? If so,&
                                  & make sure the block is valid and has right syntax.'
        Else
          Write (message,'(1x,a)') Trim(set_error)//' Unknown directive found: "'//Trim(word)//'".&
                                  & Have you correctly defined the previous directives? Have you forgotten something maybe?'
        End If 
        Call error_stop(message)
      End If

    End Do
    ! Close file
    Close(files(FILE_SET)%unit_no)

  End Subroutine read_settings

  Subroutine read_input_parameters(iunit, fit_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to read input parameters from block
    !
    ! author    - i.scivetti March 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Integer(Kind=wi),  Intent(In   ) :: iunit
    Type(fit_type),    Intent(InOut)  :: fit_data 

    Integer(Kind=wi)   :: io, length
    Character(Len=256) :: message, word
    Character(Len=256) :: set_error
    
    set_error = '***ERROR in the &input_parameters block (SETTINGS file).'

    Do
      Read (iunit, Fmt=*, iostat=io) word
      If (io /= 0) Then
        Write (message,'(2(1x,a))') Trim(set_error), 'It appears the block has not been closed correctly. Use&
                                  & "&end_input_parameters" to close the block.&
                                  & Check if directives are set correctly.'         
        Call error_stop(message) 
      End If  
      
      Call get_word_length(word,length)
      Call capital_to_lower_case(word)
      If (Trim(word)=='&end_input_parameters') Exit
      Call check_for_rubbish(iunit, '&input_parameters')

      If (word(1:1) == '#' .Or. word(1:3) == '   ') Then
      ! Do nothing if line is a comment of we have an empty line
      Read (iunit, Fmt=*, iostat=io) word

      Else If (Trim(word)=='a0') Then
        Read (iunit, Fmt=*, iostat=io) fit_data%param%A0%ini%tag, fit_data%param%A0%ini%value
        Call set_read_status(word, io, fit_data%param%A0%ini%fread, fit_data%param%A0%ini%fail)

      Else If (Trim(word)=='a1') Then
        Read (iunit, Fmt=*, iostat=io) fit_data%param%A1%ini%tag, fit_data%param%A1%ini%value
        Call set_read_status(word, io, fit_data%param%A1%ini%fread, fit_data%param%A1%ini%fail)

      Else If (Trim(word)=='a2') Then
        Read (iunit, Fmt=*, iostat=io) fit_data%param%A2%ini%tag, fit_data%param%A2%ini%value
        Call set_read_status(word, io, fit_data%param%A2%ini%fread, fit_data%param%A2%ini%fail)

      Else If (Trim(word)=='a3') Then
        Read (iunit, Fmt=*, iostat=io) fit_data%param%A3%ini%tag, fit_data%param%A3%ini%value
        Call set_read_status(word, io, fit_data%param%A3%ini%fread, fit_data%param%A3%ini%fail)

      Else If (Trim(word)=='t1') Then
        Read (iunit, Fmt=*, iostat=io) fit_data%param%T1%ini%tag, fit_data%param%T1%ini%value
        Call set_read_status(word, io, fit_data%param%T1%ini%fread, fit_data%param%T1%ini%fail)
        
      Else If (Trim(word)=='t2') Then
        Read (iunit, Fmt=*, iostat=io) fit_data%param%T2%ini%tag, fit_data%param%T2%ini%value
        Call set_read_status(word, io, fit_data%param%T2%ini%fread, fit_data%param%T2%ini%fail)

      Else If (Trim(word)=='t3') Then
        Read (iunit, Fmt=*, iostat=io) fit_data%param%T3%ini%tag, fit_data%param%T3%ini%value
        Call set_read_status(word, io, fit_data%param%T3%ini%fread, fit_data%param%T3%ini%fail)
        
      Else
        Write (message,'(1x,5a)') Trim(set_error), ' Directive "', Trim(word),&
                                & '" is not recognised as a valid settings.',&
                                & ' Have you properly closed the block with "&end_input_parameters"?'
        Call error_stop(message)
      End If

    End Do
    
  End Subroutine read_input_parameters
  
  Subroutine read_msd(iunit, fit_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to read the settigns for mean square displacement (MSD)
    ! analysis from the &MSD block
    !
    ! author    - i.scivetti March 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Integer(Kind=wi),  Intent(In   ) :: iunit
    Type(fit_type), Intent(InOut)  :: fit_data 

    Integer(Kind=wi)   :: io, length
    Character(Len=256) :: message, word
    Character(Len=256) :: set_error
    
    set_error = '***ERROR in the &MSD block (SETTINGS file).'

    Do
      Read (iunit, Fmt=*, iostat=io) word
      If (io /= 0) Then
        Write (message,'(2(1x,a))') Trim(set_error), 'It appears the block has not been closed correctly. Use&
                                  & "&end_msd" to close the block.&
                                  & Check if directives are set correctly.'         
        Call error_stop(message) 
      End If  
      
      Call get_word_length(word,length)
      Call capital_to_lower_case(word)
      If (Trim(word)=='&end_msd') Exit
      Call check_for_rubbish(iunit, '&msd')

      If (word(1:1) == '#' .Or. word(1:3) == '   ') Then
      ! Do nothing if line is a comment of we have an empty line
      Read (iunit, Fmt=*, iostat=io) word

      Else If (Trim(word)=='select') Then
        Read (iunit, Fmt=*, iostat=io) word, fit_data%msd%select%type
        Call set_read_status(word, io, fit_data%msd%select%fread, fit_data%msd%select%fail,&
                           & fit_data%msd%select%type)

      Else If (Trim(word)=='units') Then
         Read (iunit, Fmt=*, iostat=io) word, fit_data%msd%units%type
         Call set_read_status(word, io, fit_data%msd%units%fread, fit_data%msd%units%fail, &
                             fit_data%msd%units%type)

      Else
        Write (message,'(1x,5a)') Trim(set_error), ' Directive "', Trim(word),&
                                & '" is not recognised as a valid settings.',&
                                & ' Have you properly closed the block with "&end_msd"?'
        Call error_stop(message)
      End If

    End Do
    
  End Subroutine read_msd
  
End module input_settings

