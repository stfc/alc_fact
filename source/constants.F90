!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Module containing constants and parameters for computation
!
! Author    - i.scivetti March 2026
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Module constants

  Use numprec, Only: wi, &
                     wp

  Implicit None

  ! Code reference 
  Character(Len=16), Parameter, Public  :: code_name    = "ALC_FACT" 
  Character(Len=16), Parameter, Public  :: code_VERSION = "0.1"
  Character(Len=16), Parameter, Public  :: date_RELEASE = "Jun 2026"
  
  Integer(Kind=wi),  Parameter, Public  :: min_gnuplot_version=5
  Integer(Kind=wi),  Parameter, Public  :: min_gnuplot_subversion=4
  
  Integer(Kind=wi),  Parameter, Public  :: max_std_samples=51
  Integer(Kind=wi),  Parameter, Public  :: max_iterations_fitting=2000
    
  Real(Kind=wp), Parameter, Public  :: pi    = 3.14159265358979312e0_wp
  Real(Kind=wp), Parameter, Public  :: Bohr_to_A = 0.529177249_wp
  Real(Kind=wp), Parameter, Public  :: exp_ocf_lower_bound = 0.2_wp
  
  Character(Len=16), Parameter, Public  :: fmt_mean_std  = '(1x,a,6x,2e13.4)'
  Character(Len=16), Parameter, Public  :: fmt_mean_only = '(1x,a,6x,1e13.4)'

  ! References
  Character(Len=256), Parameter, Public  :: morrone2006  = "J. A. Morrone, et al., J. Phys. Chem. B 110, 3712 (2006)"
  Character(Len=256), Parameter, Public  :: lipari1980   = "Lipari, G.; Szabo, A. Biophys. J 30, 489 (1980)"
  Character(Len=256), Parameter, Public  :: calio2021    = 'Calio, P. B. et al. JACS 143 (44), 18672-18683 (2021)'
  
End Module constants
