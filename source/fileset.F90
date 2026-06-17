!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
! Module for input/output files and related subroutines
!
! Author -    i.scivetti  March 2026
!!!!!!!!!!!!!!!!!!!!11!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Module fileset

  Use constants,    Only: code_name,&
                          code_VERSION, &
                          date_RELEASE
  Use numprec,      Only: wi
  use unit_output,  Only: info, &
                          set_output_unit 

  Implicit None
  Private

  ! File data
  Type, Public :: file_type
    Private
    ! Filename
    Character(Len=256), Public :: filename
    ! Fortran unit number, set with newunit=T%unit_no
    Integer(Kind=wi), Public   :: unit_no = -2
  Contains
    Procedure, Public :: init => file_type_init
    Procedure, Public :: Close => close_file
  End Type file_type

  ! SET file
  Integer(Kind=wi), Parameter, Public :: FILE_SET= 1
  ! OUT file
  Integer(Kind=wi), Parameter, Public :: FILE_OUT= 2 
  ! GNUPLOT_INPUT file
  Integer(Kind=wi), Parameter, Public :: FILE_GNUPLOT_INPUT= 3 
  ! GNUPLOT_OUTPUT file
  Integer(Kind=wi), Parameter, Public :: FILE_GNUPLOT_OUTPUT= 4
  ! FIT_PARAM file
  Integer(Kind=wi), Parameter, Public :: FILE_FIT_PARAMS= 5
  ! DUMP file
  Integer(Kind=wi), Parameter, Public :: FILE_DUMP= 6 
  ! INPUT_DATA file
  Integer(Kind=wi), Parameter, Public :: FILE_INPUT_DATA= 7 
  ! RAW_DATA file
  Integer(Kind=wi), Parameter, Public :: FILE_RAW_DATA= 8
  ! FITTED_DATA file
  Integer(Kind=wi), Parameter, Public :: FILE_FITTED_DATA= 9
  
  ! Size of filename array
  Integer(Kind=wi), Parameter, Public :: NUM_FILES = 9

  Public :: set_system_files, print_header_out, wrapping_up, refresh_out

Contains

  Subroutine refresh_out(files)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to refresh the output
    !
    ! author    - i.scivetti March 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type) :: files(NUM_FILES) 

    Call files(FILE_OUT)%close ()
    Open (Newunit=files(FILE_OUT)%unit_no, File=files(FILE_OUT)%filename, Position='Append')
    Call set_output_unit(files(FILE_OUT)%unit_no)

  End Subroutine refresh_out

  Subroutine file_type_init(T, filename)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to initialise files
    !
    ! author    - i.scivetti March 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
    Class(file_type)                :: T
    Character(Len=*), Intent(In   ) :: filename

    T%filename = Trim(filename)
  End Subroutine file_type_init


  Subroutine set_names_files(files)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Set default names for files
    !
    ! author    - i.scivetti March 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type) :: files(NUM_FILES)

    Character(Len=256), Dimension(NUM_FILES)   :: set_names
    Integer(Kind=wi)                           :: file_no

    ! Default file names array
    ! Populate default names array
    set_names(FILE_SET)             = "SETTINGS"
    set_names(FILE_OUT)             = "OUTPUT"
    set_names(FILE_GNUPLOT_INPUT)   = "GNUPLOT_INPUT"
    set_names(FILE_GNUPLOT_OUTPUT)  = "GNUPLOT_OUTPUT"
    set_names(FILE_FIT_PARAMS)      = "FIT_PARAMS"
    set_names(FILE_DUMP)            = "DUMP"
    set_names(FILE_INPUT_DATA)      = "INPUT_DATA"
    set_names(FILE_RAW_DATA)        = "RAW_DATA"
    set_names(FILE_FITTED_DATA)     = "FITTED_DATA"
    ! Set default filenames
    Do file_no = 1, NUM_FILES
      Call files(file_no)%init(set_names(file_no))
    End Do

  End Subroutine set_names_files


  Subroutine close_file(T)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to close files
    !
    ! author    - i.scivetti March 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
    Class(file_type) :: T

    Logical :: is_open

    Inquire (T%unit_no, opened=is_open)
    If (is_open) Then
      Close (T%unit_no)
      T%unit_no = -2
    End If

  End Subroutine close_file

  Subroutine set_system_files(files)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to open OUTPUT file 
    ! 
    ! author    - i.scivetti April 2020
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type) :: files(NUM_FILES)

    Call set_names_files(files)   
    Open (Newunit=files(FILE_OUT)%unit_no, File=files(FILE_OUT)%filename, Status='replace')
    Call set_output_unit(files(FILE_OUT)%unit_no)

  End Subroutine set_system_files   

  Subroutine print_header_out(files)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to print the header to OUTPUT file 
    !  
    ! author        - i.scivetti March 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type) :: files(NUM_FILES)

    Character(Len=*), Parameter :: fmt1 = '(a)'
    Character(Len=*), Parameter :: fmt2 = '(3a)'
    Character(Len=*), Parameter :: fmt3 = '(4a)'
    Character(Len=128)          :: header(14)

    Write (header(1), fmt1)   Repeat("#", 74)
    Write (header(2), fmt2)  "#                      WELCOME TO ", Trim(code_name),  Repeat(" ", 31)//"#"
    Write (header(3), fmt1)  "#  A software to analyse orientational and transfer correlations (and    #"
    Write (header(4), fmt1)  "#  mean square displacements) through automatic fitting with GNUPLOT     #"
    Write (header(5), fmt3)  "#  version:  ", Trim(code_VERSION), Repeat(' ',57),                     "#"
    Write (header(6), fmt3)  "#  release:  ", Trim(date_RELEASE), Repeat(' ',52),                     "#"
    Write (header(7), fmt1)  "#                                                                        #"
    Write (header(8), fmt1)  "#  Copyright:  2026 Ada Lovelace Centre (ALC)                            #"
    Write (header(9), fmt1)  "#              Scientific Computing Department (SCD)                     #"
    Write (header(10), fmt1) "#              Science and Technology Facilities Council (STFC)          #"
    Write (header(11), fmt1) "#                                                                        #"
    Write (header(12), fmt1) "#  Author:     Ivan Scivetti (SCD/STFC)                                  #"
    Write (header(13), fmt1) "#                                                                        #"
    Write (header(14), fmt1)  Repeat("#", 74)
    Call info(header, 14)
    
    ! Refresh OUT_EQCM
    Call refresh_out(files)

  End Subroutine print_header_out

  Subroutine wrapping_up(files)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Subroutine to print final remarks to OUT_EQCM file 
  ! and close the file 
  !  
  ! author    - i.scivetti March 2026
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type) :: files(NUM_FILES)
 
    Character(Len=*), Parameter :: fmt1 = '(1x,a)'
    Character(Len=*), Parameter :: fmt2 = '(1x,3a)'
    Character(Len=128)          :: appex(7)
     
    Write (appex(1), fmt1)   Repeat(" ", 1)
    Write (appex(2), fmt1)   Repeat("#", 38)
    Write (appex(3), fmt1)  "#                                    #" 
    Write (appex(4), fmt1)  "#  Job has finished successfully     #"
    Write (appex(5), fmt2)  "#  Thanks for using ", Trim(code_name), "!        #"
    Write (appex(6), fmt1)  "#                                    #" 
    Write (appex(7), fmt1)   Repeat("#", 38)
    Call info(appex, 7)

    Close(files(FILE_OUT)%unit_no)    

  End Subroutine wrapping_up  

End Module fileset
