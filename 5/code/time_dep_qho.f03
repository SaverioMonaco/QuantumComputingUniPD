module debugmod
! ---------------------------------------------------------------------!
! --------------------------- DOCUMENTATION ---------------------------!
! ---------------------------------------------------------------------!
! TYPES: none                                                          !
! ---------------------------------------------------------------------!
! FUNCTIONS: none                                                      !
! ---------------------------------------------------------------------!
! SUBROUTINES:                                                         !
! + debugging                                                          !
!   INPUT:                                                             !
!   > condition = logic, if condition is TRUE then run the content of  !
!                        this function                                 !
!   > msg = character, optional, message to print                      !
!   > content = some variable to inspect                               !
! ---------------------------------------------------------------------!
! ---------------------------------------------------------------------!
  contains

  subroutine debugging(condition, msg, content)
    logical, intent(IN)                 :: condition
    character(*), intent(IN), optional  :: msg
    class(*), intent(IN), optional      :: content

    if(condition) then
      if (present(content)) then
        select type(content)
          type is (integer(1))
            print*, msg, " => [OK], Variable = ", content
          type is (integer(2))
            print*, msg, " => [OK], Variable = ", content
          type is (integer(4))
            print*, msg, " => [OK], Variable = ", content
          type is (integer(8))
            print*, msg, " => [OK], Variable = ", content
          type is (real(4))
            print*, msg, " => [OK], Variable = ", content
          type is (real(8))
            print*, msg, " => [OK], Variable = ", content
          type is (logical)
            print*, msg, " => [OK], Variable = ", content
        end select
      else
        print*, msg, " => [OK]"
      end if
    end if
    end subroutine
end module debugmod

module q_util
  implicit none

  contains
  function q_util_norm(psi,step)result(normpsi)
    implicit none
    double complex, dimension(:), allocatable :: psi
    double precision                          :: step, normpsi

    normpsi = (SUM(ABS(psi)**2 * step))

  end function q_util_norm

  function q_util_ExpectedValue(psi, obs, step)result(E)
    implicit none
    double complex, dimension(:), allocatable   :: psi
    double precision, dimension(:), allocatable :: obs
    double precision                            :: step, E
    integer                                     :: ii

    E = 0
    do ii = 1, size(psi)
        E = E + abs(psi(ii))**2 *obs(ii)* step
    end do
end function q_util_ExpectedValue
end module

program shroedingertimedependent
  use debugmod
  use q_util
  use, intrinsic :: iso_c_binding

  implicit none

  include 'fftw3.f03'

  double precision           :: Lx, Lp, T, dx, dt, dp
  integer                    :: Nx, Nt, xx, tt
  character(20)              :: folder
  integer                    :: iostat  ! Checking types in READ(*,*)

  double precision, dimension(:), allocatable :: xgrid, pgrid, ps
  double complex, dimension(:), allocatable   :: Uv, Ut

  ! Vectors for:
  ! tmp_k1: after Ux, Fourier transformed
  ! tmp_k2: after Uk
  ! tmp_x2: after inverse transforming
  complex(kind=8), dimension(:), allocatable   :: psi_x0, psi_x1, psi_k1, psi_k2, psi_x2, psi_x3, psi_k3

  ! Constants
  double precision :: pi

  ! For debugging
  logical :: DEBUG, DEBUG2

  ! FFTW Related
  integer*8 :: dfft_plan, idfft_plan, jdfft_plan

  ! General
  DEBUG = .FALSE.
  ! Too check if first step matches with theory
  DEBUG2 = .TRUE.

  print*, "--------------------------------------------"
  print*, "+   TIME DEPENDENT SCHROEDINGER EQUATION   +"
  print*, "+------------------------------------------+"
  print*, "+ Lx:  Lenght of x space                   +"
  print*, "+ Lp:  Lenght of p space                   +"
  print*, "+ Nx:   Number of points                   +"
  print*, "+ Nt:  Number of time-points               +"
  print*, "+ T:   Total time                          +"
  print*, "+------------------------------------------+"
  print*, "+ Type: Lx, Lp, Nx, Nt, T and folder:      +"
  write(*,"(A)",advance='no') " + "

  read (*,*, iostat=iostat) Lx, Lp, Nx, Nt, T, folder

  ! CHECK IF PARAMETERS ARE THE RIGHT TYPES
  if(iostat /= 0) then
    print*, "+ !!! INVALID PARAMETER TYPES !!!"
    print*, "+ Nx   = int"
    print*, "+ Nt  = int"
    print*, "+ Lx  = real, int"
    print*, "+ Lp  = real, int"
    print*, "+ T   = real, int"
    print*, "+ Exiting..."
    print*, "--------------------------------------------"
    stop
  end if

  ! CHECK IF PARAMETERS ARE IN THE RIGHT RANGE
  if(Lx <= 0 .OR. Lp <= 0 .OR. T<=0 .OR. Nt<2 .OR. Nx<3) then
    print*, "+ !!! INVALID PARAMETER RANGES !!!"
    print*, "+ (Lx > 0, Lp > 0, Nx > 2, Nt > 1, T > 0)"
    print*, "+ Exiting..."
    print*, "--------------------------------------------"
    stop
  end if

  print*, "+ Data will be saved in: ./"//trim(folder)
  print*, "+ Lenght of x space     (Lx): ", Lx
  print*, "+ Lenght of p space     (Lp): ", Lp
  print*, "+ Number of points       (N): ", Nx
  print*, "+ Number of time-points (Nt): ", Nt
  print*, "+ Total time             (T): ", T
  print*, "+------------------------------------------+"

  Nx = Nx + 1

  pi  = 4.d0 * datan(1.d0)
  dx  = 2*(Lx)/(Nx-1)
  dp  = 2*(Lp)/(Nx-1)
  dt  = T/Nt

  allocate(xgrid(Nx))
  allocate(pgrid(Nx))
  allocate(ps(Nx))

  ! Prepare the x-space and p-space grids
  do xx = 1, Nx, 1
    xgrid(xx) = -(Lx) + (xx-1)*dx
    pgrid(xx) = -(Lp) + (xx-1)*dp
  end do
  call debugging(DEBUG,msg="+ Allocating xgrid and pgrid")

  ! reordering necessary due to how the FFTW works
  ! positive frequencies in the first half, then negative frequencies in the second half
  do xx = 1, int(Nx/2), 1
    ps(xx) = dp * dble(xx-1)
  end do
  do xx = int(Nx/2) + 1, Nx, 1
      ps(xx) = -1d0 * dble(Nx - xx + 1) * dp
  end do

  call debugging(DEBUG,msg="+ Reordering ps")

  ! set filenames
  open(11, file = "E_x.csv")
  open(12, file = "sigma_x.csv")
  open(13, file = "E_p.csv")
  open(14, file = "sigma_p.csv")
  open(4,  file = "real_wave.csv")
  open(5,  file = "imag_wave.csv")

  ! Writing on a file the first step to see if it matches with theory
  if(DEBUG2) then
    open(40,  file = "psix0.csv")
    open(41,  file = "psix1.csv")
    open(42,  file = "psik1.csv")
    open(43,  file = "psik2.csv")
    open(44,  file = "psix2.csv")
    open(45,  file = "psix3.csv")
  end if
  call debugging(DEBUG,msg="+ Files opened")

  allocate(psi_x0(Nx), psi_x1(Nx), psi_x2(Nx), psi_x3(Nx))
  allocate(psi_k1(Nx), psi_k2(Nx), psi_k3(Nx))
  allocate(Uv(Nx),Ut(Nx))

  call debugging(DEBUG,msg="+ Allocating psi-s")

  ! DFFTW plans generation
  call dfftw_plan_dft_1d(dfft_plan,  Nx, psi_x1, psi_k1, FFTW_FORWARD, FFTW_MEASURE)
  call dfftw_plan_dft_1d(idfft_plan, Nx, psi_k2, psi_x2, FFTW_BACKWARD, FFTW_MEASURE)
  call dfftw_plan_dft_1d(jdfft_plan,  Nx, psi_x3, psi_k3, FFTW_FORWARD, FFTW_MEASURE)
  call debugging(DEBUG,msg="+ Launching dfftw_plan_dft_1d")

  ! Set the initial state of the wavefunction:
  do xx = 1, Nx
    psi_x0(xx) = pi**(-0.25d0) * EXP(-(xgrid(xx)**2d0)/2d0)
  end do
  psi_x0 = psi_x0/q_util_norm(psi_x0, dx)
  !psi_x0 = psi_x0/q_util_norm(psi_x0, dble(2*Lx/Nx) )
  call debugging(DEBUG,msg="+ Preparing first wavefunc")

  ! ###########################################################################
  ! ########################## SPLIT OPERATOR METHOD ##########################
  ! ###########################################################################
  print*, "+"
  print*, "+ Evolving:"

  do tt = 1, Nt, 1
    ! Drawing percentages
    if(mod(tt,Nt/10)==0) then
      write(*,*) '+ ', tt*100/Nt, '%'
    end if

    ! PREPARE evolution operators
    do xx = 1, Nx, 1
        ! kinetic
        Ut(xx) = EXP(dcmplx(0d0,-0.5d0) * dt * ps(xx)**2d0 )
        ! potential
        Uv(xx) = EXP(dcmplx(0d0,-0.5d0) * dt * (0.5d0 *(xgrid(xx) - tt/Nt)**2d0) )
    end do
    if(tt==1) then
      call debugging(DEBUG,msg="+ Evolution operators initialized")

      if(DEBUG2) then
        write(40, *) (/ real(psi_x0), aimag(psi_x0) /)
      end if
    end if

    ! Evolve Ux Ψx
    do xx = 1, Nx, 1
      psi_x1(xx) = Uv(xx) * psi_x0(xx)
    end do
    if(tt==1) then
      call debugging(DEBUG,msg="+ Uxpsi0")

      if(DEBUG2) then
        write(41, *) (/ real(psi_x1), aimag(psi_x1) /)
      end if
    end if

    ! Transform Ψx -> Ψk
    call dfftw_execute_dft(dfft_plan, psi_x1, psi_k1)
    if(tt==1) then
      call debugging(DEBUG,msg="+ Fpsi1x")

      if(DEBUG2) then
        write(42, *) (/ real(psi_k1), aimag(psi_k1) /)
      end if
    end if

    ! Evolve Uk Ψk
    do xx = 1, Nx, 1
      psi_k2(xx) = Ut(xx) * psi_k1(xx)/sqrt( dble(Nx) )
    end do
    if(tt==1) then
      call debugging(DEBUG,msg="+ Ukpsi1")

      if(DEBUG2) then
        write(43, *) (/ real(psi_k2), aimag(psi_k2) /)
      end if
    end if

    ! Transform Ψk -> Ψx
    call dfftw_execute_dft(idfft_plan, psi_k2, psi_x2)
    if(tt==1) then
      call debugging(DEBUG,msg="+ F^-1psi2k")
      if(DEBUG2) then
        write(44, *) (/ real(psi_x2), aimag(psi_x2) /)
      end if
    end if

    ! Evolve Ux Ψx
    do xx = 1, Nx, 1
      psi_x3(xx) = Uv(xx) * psi_x2(xx)/sqrt( dble(Nx) )
    end do
    if(tt==1) then
      call debugging(DEBUG,msg="+ Uxpsi2")

      if(DEBUG2) then
        write(45, *) (/ real(psi_x3), aimag(psi_x3) /)
      end if
    end if

    psi_x3 = psi_x3/q_util_norm(psi_x3, dble(dx) )
    ! Set initial psi as current psi for next time iteration
    psi_x0 = psi_x3
    if(tt==1) then
      call debugging(DEBUG,msg="+ Normalizing Output 1/2")
    end if

    write(4, *) real(psi_x0)
    write(5, *) aimag(psi_x0)

    call dfftw_execute_dft(jdfft_plan, psi_x3, psi_k3)
    psi_k3 = psi_k3/q_util_norm(psi_k3, dble(dx) )
    if(tt==1) then
      call debugging(DEBUG,msg="+ Normalizing Output 2/2")
    end if
  end do
  ! ###########################################################################
  ! ###################### END OF SPLIT OPERATOR METHOD #######################
  ! ###########################################################################

  ! Deallocate everything
  deallocate(xgrid)
  deallocate(pgrid)
  deallocate(psi_x0, psi_x1, psi_x2, psi_x3)
  deallocate(psi_k1, psi_k2, psi_k3)
  deallocate(ps)

  print*, "+"
  print*, "+ All done! Exiting..."
  print*, "+------------------------------------------+"
  print*, "+------------------------------------------+"
end program shroedingertimedependent

! Compile:
!  gfortran shiftingV_harmonic_oscillator.f90 -o sh_qho -llapack -Wall -lfftw3 -frecursive
