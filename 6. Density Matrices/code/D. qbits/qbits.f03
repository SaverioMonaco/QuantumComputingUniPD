module basechange
  ! ---------------------------------------------------------------------!
  ! --------------------------- DOCUMENTATION ---------------------------!
  ! ---------------------------------------------------------------------!
  ! TYPES: none                                                          !
  ! ---------------------------------------------------------------------!
  ! FUNCTIONS:                                                           !
  ! + basechange_to                                                      !
  !   INPUT:                                                             !
  !   > b_to = integer, base for the conversion 10 -> b_to               !
  !   > number = integer, number to convert number_10 -> number_b_to     !
  !   > N = integer, size of number: 4_2 is 0100 if N = 4                !
  !   OUTPUT:                                                            !
  !   > number_b_to = integer(N), number converted                       !
  !                                                                      !
  ! + basechange_from                                                    !
  !   INPUT:                                                             !
  !   > b_from = integer, base for the conversion b_to -> 10             !
  !   > number_from = integer(N), number to convert from b_from to 10    !
  !   > N = integer, size of number: 4_2 is 0100 if N = 4                !
  !   OUTPUT:                                                            !
  !   > number_b10 = integer, number converted                           !
  ! ---------------------------------------------------------------------!
  ! SUBROUTINES: none                                                    !
  ! ---------------------------------------------------------------------!
  ! ---------------------------------------------------------------------!
  contains

  function basechange_to(b_to, number, N) result(number_b_to)
    integer :: b_to, number, N, ii
    integer, dimension(N) :: number_b_to

    number_b_to = 0*number_b_to ! Allocate to 0, just to be sure
    do ii = 1, N, 1
      number_b_to(N - ii + 1) = modulo(number, b_to)
      number = number/b_to
    end do
  end function

  function basechange_from(b_from, number_from, N) result(number_b10)
    integer :: b_from, number_b10, N, ii
    integer, dimension(N) :: number_from

    do ii = 1, N, 1
      number_b10 = number_b10 + number_from(N - ii + 1)*b_from**(ii - 1)
    end do
  end function
end module

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
            print*, "+ ", msg, " => [OK], Variable = ", content
          type is (integer(2))
            print*, "+ ", msg, " => [OK], Variable = ", content
          type is (integer(4))
            print*, "+ ", msg, " => [OK], Variable = ", content
          type is (integer(8))
            print*, "+ ", msg, " => [OK], Variable = ", content
          type is (real(4))
            print*, "+ ", msg, " => [OK], Variable = ", content
          type is (real(8))
            print*, "+ ", msg, " => [OK], Variable = ", content
          type is (logical)
            print*, "+ ", msg, " => [OK], Variable = ", content
        end select
      else
        print*, "+ ", msg, " => [OK]"
      end if
    end if
    end subroutine
end module debugmod

module densmat
  use basechange
  implicit none
  ! ---------------------------------------------------------------------!
  ! --------------------------- DOCUMENTATION ---------------------------!
  ! ---------------------------------------------------------------------!
  ! TYPES:                                                               !
  ! + qsystem:                                                           !
  !           Structure for quantum multiple subsystems                  !
  !           > N           : integer, number of states                  !
  !           > d           : actual, dimension of each Hilbert space    !
  !           > separability: wether the system is separable or not      !
  !           > waves       : elements/coefficients of the systems       !
  ! ---------------------------------------------------------------------!
  ! FUNCTIONS:                                                           !
  ! + densmat_pure_init                                                  !
  !   Initializing function for qsystems                                 !
  !   INPUT:                                                             !
  !   > N: integer, number of states                                     !
  !   > d: integer, dimension of each Hilbert spaces                     !
  !   > SEP: integer if separable (SEP=1), if not (SEP=0)                !
  !   > DEBUG_: logical, prints bytes allocated (can be too many)        !
  !   OUTPUT:                                                            !
  !   > PHI: qsystem                                                     !
  ! ---------------------------------------------------------------------!
  ! SUBROUTINES:                                                         !
  ! + densmat_genstates                                                  !
  !   Generate random numbers (dxN or d^N double complex numbers)        !
  !  INPUT:                                                              !
  !  > PHI: qsystem                                                      !
  ! ---------------------------------------------------------------------!
  ! ---------------------------------------------------------------------!
  type qsystem
    integer :: N ! Number of systems
    integer :: d ! number of states

    logical :: separability ! Wether the whole system is separable
                            ! (no interactions) or not, the number of
                            ! coefficients changes:
                            ! separability == TRUE  -> #C = d X N
                            !              == FALSE -> #C = d ^ N

    double complex, dimension(:), allocatable :: waves
  end type qsystem

  contains

  function densmat_pure_separable(N, d) result(PHI)
    integer   :: N, d
    type(qsystem) :: PHI

    PHI%N = N
    PHI%d = d
    PHI%separability = .TRUE.

    allocate(PHI%waves(N*d))
  end function

  function densmat_pure_unseparable(N, d) result(PHI)
    integer   :: N, d
    type(qsystem) :: PHI

    PHI%N = N
    PHI%d = d
    PHI%separability = .FALSE.

    allocate(PHI%waves(d**N))
  end function

  function densmat_pure_init(N, d, SEP, DEBUG) result(PHI)
    integer       :: N, d, SEP
    type(qsystem) :: PHI
    logical       :: DEBUG

    if(SEP == 1) then
      PHI = densmat_pure_separable(N, d)
      if(DEBUG) then
        print*, "+ Allocating", 16*N*d, "bytes"
      end if
    else if(SEP == 0) then
      PHI = densmat_pure_unseparable(N, d)
      if(DEBUG) then
        print*, "+ Allocating", 16*(N**d), "bytes"
      end if
    else
      print*, "+ Unvalid state, exiting .... "
      stop
    end if
  end function

  subroutine densmat_genstates(PHI)
    type(qsystem) :: PHI
    double precision, dimension(:), allocatable :: psi_real, psi_imag
    double precision                            :: norm

    if(PHI%separability .eqv. .TRUE.) then
      ! If separable allocate dXN
      allocate(psi_real(PHI%d*PHI%N),psi_imag(PHI%d*PHI%N))
    else
      ! else allocate d^N
      allocate(psi_real(PHI%d**PHI%N),psi_imag(PHI%d**PHI%N))
    end if

    call random_number(psi_real) ! Generate all real
    call random_number(psi_imag) ! Generate all imag
    PHI%waves = dcmplx(psi_real, psi_imag)
    ! Normalizing
    norm = SUM(ABS(PHI%waves(:))**2)
    ! Assigning
    PHI%waves(:) =  PHI%waves(:)/(sqrt(norm))
  end subroutine
end module densmat

program density_matrices
  use debugmod
  use densmat

  implicit none

  integer :: N, d, SEP, randomstate
  real    :: randomstate_fl
  type(qsystem) :: PHI

  ! Debugging
  logical :: DEBUG
  integer :: iostat

  ! Loop variables
  integer :: ii

  DEBUG = .TRUE.
  print*, "--------------------------------------------"
  print*, "+       MULTI-SYSTEMS REPRESENTATION       +"
  print*, "+------------------------------------------+"
  print*, "+ A: This program allocates random         +"
  print*, "+    wavefunction for pure-state system    +"
  print*, "+    either separable or unseparable       +"
  print*, "+------------------------------------------+"
  print*, "+ N: Number of systems (integer)           +"
  print*, "+ d: Dimension of Hilber spaces (integer)  +"
  print*, "+ SEP: Separability (integer):             +"
  print*, "+      1 -> Separable                      +"
  print*, "+      0 -> Non-separable                  +"
  print*, "+------------------------------------------+"
  write(*,"(A)",advance='no') " + Type N, d, SEP: "
  read (*,*, iostat=iostat) N, d, SEP

  ! Check if parameters are the right types
  if(iostat /= 0) then
    print*, "+ !!! INVALID PARAMETER TYPES !!!          +"
    print*, "+ N = integer                              +"
    print*, "+ d = integer                              +"
    print*, "+ SEP = integer                            +"
    print*, "+ Exiting...                               +"
    print*, "+------------------------------------------+"
    print*, "+------------------------------------------+"
    stop
  end if

  ! Check if parameters are in the right ranges
  if((N < 2 .OR. d < 2) .OR. .NOT.(SEP == 0 .OR. SEP == 1)) then
    print*, "+ !!! INVALID PARAMETER RANGES !!!         +"
    print*, "+ (N > 1, d > 1, S = {0,1})                +"
    print*, "+ Exiting...                               +"
    print*, "+------------------------------------------+"
    print*, "+------------------------------------------+"
    stop
  end if

  ! Check on READ
  if(DEBUG) then
    write(*,"(A)",advance='no') " + N:"
    print*, N
    write(*,"(A)",advance='no') " + d:"
    print*, d
    if(SEP==1) then
      print*, "+ System is separable"
    else
      print*, "+ System is NOT separable"
    end if
  end if

  call debugging(DEBUG,"Allocating system")
  ! We initialize the qsystem, allocating either Nxd (if separable), or d^N
  PHI = densmat_pure_init(N,d,SEP,DEBUG)

  ! Generating states
  call densmat_genstates(PHI)
  call debugging(DEBUG,"Wavefunction generated")
  if(PHI%separability .eqv. .TRUE.) then
    do ii = 1, 2, 1
      print*, "+ ",ii,":", PHI%waves(ii)
    end do
    do ii = size(PHI%waves) - 2, size(PHI%waves), 1
      print*, "+ ",ii,":", PHI%waves(ii)
    end do
  end if

  if(PHI%separability .eqv. .FALSE.) then
    print*, "+ Examples of states using tensorial notation:"
    do ii = 1, 6, 1
      call random_number(randomstate_fl)
      randomstate = FLOOR((d**N)*randomstate_fl)
      print*, "+", PHI%waves(randomstate), "|", basechange_to(d,randomstate,N) ,">"
    end do
  end if



end program density_matrices
