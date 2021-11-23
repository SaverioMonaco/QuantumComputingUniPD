# Quantum Information and Computing @UniPD 2021-2022
## Weekly assignments

### Weekly exercises in `Fortran` about numerical methods for quantum computing.
#### 1. Setup, Number Precision, Test Performance. [60/100]
     I forgot to add and compare the other version of the loop method of matrix multiplication. Make sure to add both:
         
     ```Fortran          
	  function matrix_multiplication(mat_A, mat_B) result(mat_C)
	    real*4, dimension(:,:)                         :: mat_A, mat_B
	    real*4, dimension(size(mat_A,1),size(mat_B,2)) :: mat_C
	    logical                                        :: check
	    integer                                        :: ii,jj,kk

	    ! Check if multiplication is possible (shapes)
	    if (size(mat_A,2) .eq. size(mat_B,1)) then
	      check = .TRUE.
	    else
	      print*, "Input matrices cannot be multiplied"
	      check = .FALSE.
	    end if

	    ! Initialize output matrix to zeros, when allocating it's not granted every
	    ! element is 0
	    ! Begin multiplication
	    if (check .eqv. .TRUE.) then
	      do ii = 1, size(mat_A,1), 1
		do jj = 1, size(mat_B,2), 1
		    do kk = 1, size(mat_B,1), 1
		        if (kk == 1) then
		          mat_C(ii,jj) = 0
		        end if
		        mat_C(ii, jj) = mat_C(ii, jj) + mat_A(ii, kk)*mat_B(kk, jj)
		    end do
		end do
	      end do
	    end if
	  end function matrix_multiplication

	  !> Other order
	  function matrix_multiplication2(mat_A, mat_B) result(mat_C)
	    real*4, dimension(:,:)                         :: mat_A, mat_B
	    real*4, dimension(size(mat_A,1),size(mat_B,2)) :: mat_C
	    logical                                        :: check
	    integer                                        :: ii,jj,kk

	    ! Check if multiplication is possible (shapes)
	    if (size(mat_A,2) .eq. size(mat_B,1)) then
	      check = .TRUE.
	    else
	      print*, "Input matrices cannot be multiplied"
	      check = .FALSE.
	    end if

	    ! Begin multiplication
	    if (check .eqv. .TRUE.) then
	      do jj = 1, size(mat_B,2), 1
		do ii = 1, size(mat_A,1), 1
		    do kk = 1, size(mat_B,1), 1
		        if (kk == 1) then
		          mat_C(ii,jj) = 0
		        end if
		        mat_C(ii, jj) = mat_C(ii, jj) + mat_A(ii, kk)*mat_B(kk, jj)
		    end do
		end do
	      end do
	    end if
	  end function matrix_multiplication2
     ```
         
#### 2. Checkpoints, Documentation, Derived Types.
#### 3. Scaling of the matrix-matrix multiplication, Eigenproblem, Random Matrix Theory.
