!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
! Welcome to ALC_FACT: a ALC software to compute error of correlations
! obtained from MD trajectories
!
! Author:            Ivan Scivetti (i.scivetti)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Program alc_fact

  Use fileset,        Only: file_type, &
                            NUM_FILES, &
                            print_header_out, &
                            set_system_files, &
                            wrapping_up,&
                            refresh_out

  Use numprec,        Only: wi,& 
                            wp
                             
  Use input_settings, Only: read_settings

  Use gnuplot_setup,  Only: check_gnuplot_availability                         

  Use fitting_setup,  Only: fit_type
  
  Use fitting_checks, Only: check_fitting_settings   
  
  Use fitting,        Only: obtain_parameters_from_fitting

  Use unit_output,    Only: info
  
  Use models,         Only: model_type, &
                            extract_physics_from_fitting
  
  Implicit None
! Definition of variables
  Type(file_type)      :: files(NUM_FILES)
  Type(fit_type)       :: fit_data
  Type(model_type)     :: model_data

  !Time related variables
  Integer(kind=wi)   :: start,finish,rate

  ! Array to print information
  Character(Len=256) :: message
  
  ! Start of the code 
  !!!!!!!!!!!!!!!!!!!
  ! Record initial time
  Call system_clock(count_rate=rate)
  Call system_clock(start)
  ! Initialise settings for input/output files
  Call set_system_files(files)
  ! Print header of OUTPUT
  Call print_header_out(files)
  ! Check if gnuplot is available
  Call check_gnuplot_availability(files)
  ! Read settings from SETTINGS
  Call read_settings(files, fit_data)
  ! Check the specification of settings in SETTINGS
  Call check_fitting_settings(files, fit_data)
  ! Fit the data; compute the parameters and their errors
  Call obtain_parameters_from_fitting(files, fit_data)
  ! Extract physical quantities through models
  Call info(' ', 1) 
  Call info('Extract further information from parameters', 1) 
  Call info('===========================================', 1)
  Call refresh_out(files)  
  Call extract_physics_from_fitting(fit_data, model_data)

  ! Record final time
  Call system_clock(finish)

  ! Print execution time
  Call info(' ', 1)
  Call info(' ==========================================', 1)
  Write (message, '(1x,a,f9.3,a)') 'Total execution time = ',  Real(finish-start,Kind=wp)/rate,  ' seconds.' 
  Call info(message, 1)
  Call info(' ==========================================', 1)

  ! Print appendix to OUT_EQCM file
  Call wrapping_up(files)

End Program  alc_fact
