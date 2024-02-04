module MOD_LinAlg
    implicit none

    ! Define the precision of the real numbers
    integer, parameter :: pv = selected_real_kind(15, 307)
    

contains
!---------------------------------------------------------------------------------------------------
!                              ** LINEAR ALGEBRA TOOLS**
subroutine LSE_Solve(A, B, N, status)
    ! Computes the solution to a real system of linear equations AX = B, where 
    ! A is a square nxn matrix, B is the coefficient matrix and X are the 
    ! unknowns. The coefficient matrix is modified to store the solutions of
    ! the system of equations. In other words, when the program ends B = X.
    ! Note that A must be a square matrix.
    !Inputs:
    !   A: The square matrix
    !   B: The coefficient matrix
    !   N: The size of the square matrix
    !Outputs:
    !   B: The solution to the system of equations

    !Checked for accuracy on 9-18-2023
    
    implicit none
    integer, intent(in) :: N                    !Size of the system
    real(pv), dimension(N,N), intent(inout) :: A !Main matrix
    real(pv), dimension(N),   intent(inout) :: B !Coefficient Matrix
    integer, intent(out) :: status !0 for sucessful solve, 1 = matrix could not be inverted.

    !Check to make sure declaried size is equal to size of matrix:
    if (SIZE(A, 1) /= N .OR. SIZE(A, 2) /= N) then
        print*, "Error: The declared size of the matrix does not match the size of the matrix."
        status = 1
        return
    else
        call Invert_Matrix(A, N, status)
        if (status == 0) then
            B = matmul(A, B)
        end if
    end if
end subroutine LSE_Solve


subroutine Invert_Matrix(A, n, status)
    !Inverts a square nxn matrix A
    !Inputs:
    !   A: The square matrix to be inverted
    !   n: The size of the square matrix
    !Outputs:
    !   A: The inverted matrix
    !   status: 0 for successful inversion, 1 for an unsuccessful inversion.


    !Checked for accuracy on 9-18-2023
    implicit none
    integer, intent(in) :: n !The size of the square nxn matrix
    integer, intent(out) :: status !0 for successful inversion, 1 for an unsuccessful inversion.
    real(pv), dimension(n, n), intent(inout) :: A
    
    !Local Variables
    real(pv), dimension(n, n) :: A_check
    real(pv), dimension(n, n) :: originalA
    real(pv), dimension(n, n) :: invA
    real(pv), dimension(n) :: B, X
    integer, dimension(n) :: P
    integer :: i, j, k, maxrow
    
    real(pv) :: sum1
    real(pv) :: norm_A, norm_invA
    real(pv) :: cond_num !cond_num = kappa

    ! Store the original matrix for later use
    originalA = A

    ! Initialize status to 0 (success)
    status = 0

    ! Initialize the permutation vector
    do i = 1, n
        P(i) = i
    end do

    ! Perform LU decomposition with partial pivoting
    do k = 1, n
        maxrow = k
        do i = k + 1, n
            if (abs(A(i, k)) > abs(A(maxrow, k))) then
                maxrow = i
            end if
        end do

        ! Swap rows if needed
        if (k /= maxrow) then
            P([k, maxrow]) = P([maxrow, k])
            A([k, maxrow], :) = A([maxrow, k], :)
        end if

        ! Check for singular matrix
        if (A(k, k) == 0.0_8) then
            status = 1 !Cannot be inverted
            return
        end if

        ! Update the matrix
        do i = k + 1, n
            A(i, k) = A(i, k) / A(k, k)
            do j = k + 1, n
                A(i, j) = A(i, j) - A(i, k) * A(k, j)
            end do
        end do
    end do

    ! Compute the inverse matrix column by column
    do j = 1, n
        ! Initialize B
        B = 0.0_8
        B(j) = 1.0_8

        ! Apply the permutation to B
        B = B(P)

        ! Forward substitution
        do i = 1, n
            sum1 = B(i)
            do k = 1, i - 1
                sum1 = sum1 - A(i, k) * B(k)
            end do
            B(i) = sum1
        end do

        ! Backward substitution
        do i = n, 1, -1
            sum1 = B(i)
            do k = i + 1, n
                sum1 = sum1 - A(i, k) * X(k)
            end do
            X(i) = sum1 / A(i, i)
        end do

        ! Store the result in invA
        invA(:, j) = X
    end do

    A = invA !Output value

    !Check for a high condition number (indicates a poorly conditioned matrix)
    ! Compute the 1-norm of the original matrix
    norm_A = maxval(sum(abs(originalA), dim=1))

    ! Compute the 1-norm of the inverse matrix
    norm_invA = maxval(sum(abs(invA), dim=1))

    ! Compute the condition number
    cond_num = norm_A * norm_invA

    ! Check the condition number
    if (cond_num > 1.0E20) then
        print*, "Warning! Poorly Conditioned Matrix! See 'invert matrix subroutine'."
    end if

end subroutine Invert_Matrix


subroutine Matrix_Det(A, det, status)
    !Computes the determinant of a square matrix.
    implicit none
    real(pv), dimension(:,:), intent(in) :: A
    real(pv), intent(out) :: det
    integer, intent(out) :: status
    real(pv), dimension(:,:), allocatable :: LU_A
    integer :: i, n, swaps

    n = size(A, 1) ! Determine the size of the square matrix dynamically
    allocate(LU_A(n, n)) ! Allocate LU_A with the determined size

    ! Copy A to LU_A
    LU_A = A

    ! Perform LU decomposition
    call LU_Decomposition(LU_A, swaps, status)

    ! Check for singular matrix
    if (status == 1) then
        det = 0.0_pv
        deallocate(LU_A) ! Ensure to deallocate before returning
        return
    end if

    ! Calculate the determinant
    det = 1.0_pv
    do i = 1, n
        det = det * LU_A(i, i)
    end do

    ! Adjust for row swaps
    if (mod(swaps, 2) == 1) then
        det = -det
    end if

    deallocate(LU_A) ! Clean up by deallocating LU_A
end subroutine Matrix_Det


subroutine LU_Decomposition(A, swaps, status)
    !Performs LU Decomp with partial pivoting on an nxn matrix. 
    !Checked for accuracy on 9-18-2023
    !Inputs:
    !   A: The square matrix to be decomposed
    !Outputs:
    !   A: The decomposed matrix
    !   swaps: The number of row swaps that were performed
    !   status: 0 for successful decomposition, 1 for an unsuccessful decomposition.

    implicit none
    real(pv), dimension(:,:), intent(inout) :: A
    integer, intent(out) :: swaps, status
    integer, dimension(:), allocatable :: P
    integer :: i, j, k, maxrow, n

    n = size(A, 1) ! Determine the size of the square matrix dynamically

    ! Initialize status and swaps
    status = 0
    swaps = 0

    ! Allocate and initialize the permutation vector
    allocate(P(n))
    do i = 1, n
        P(i) = i
    end do

    ! Perform LU decomposition with partial pivoting
    do k = 1, n
        maxrow = k
        do i = k + 1, n
            if (abs(A(i, k)) > abs(A(maxrow, k))) then
                maxrow = i
            end if
        end do

        ! Swap rows if needed
        if (k /= maxrow) then
            P([k, maxrow]) = P([maxrow, k])
            A([k, maxrow], :) = A([maxrow, k], :)
            swaps = swaps + 1 !Increment the row swap counter 
        end if

        ! Check for singular matrix
        if (A(k, k) == 0.0_8) then
            status = 1
            return
        end if

        ! Update the matrix
        do i = k + 1, n
            A(i, k) = A(i, k) / A(k, k)
            do j = k + 1, n
                A(i, j) = A(i, j) - A(i, k) * A(k, j)
            end do
        end do
    end do

    ! Clean up
    deallocate(P)
end subroutine LU_Decomposition


subroutine Transpose_Square_Matrix(A, At)
    !Performs a matrix transpose on A and returns At. Matrix must be square.
    !Inputs:
    !   A: The square matrix to be transposed
    !Outputs:
    !   At: The transposed matrix

    implicit none
    real(pv), intent(in) :: A(:,:)
    real(pv), intent(out) :: At(size(A, 1), size(A, 1))
    integer :: i, j, n

    n = size(A, 1) ! Determine the size of the square matrix

    do i = 1, n
        do j = 1, n
            At(j, i) = A(i, j)
        end do
    end do
end subroutine Transpose_Square_Matrix


subroutine Create_Identity_Matrix(n, Identity)
    !Creates an nxn identity matrix.
    !Inputs:
    !   n: The size of the square matrix that is desired.
    !Outputs:
    !   Identity: The identity matrix

    implicit none
    integer, intent(in) :: n
    real(pv), dimension(n, n), intent(out) :: Identity

    integer :: i

    ! Initialize all elements to zero
    Identity = 0.0d0

    ! Set the diagonal elements to one
    do i = 1, n
        Identity(i, i) = 1.0d0
    end do

end subroutine Create_Identity_Matrix


subroutine Create_Zero_Matrix(n,Zero_Matrix)
    !Creates an nxn zero matrix. 
    !Inputs:
    !   n: The size of the square matrix
    !Outputs:
    !   Zero_Matrix: The zero matrix

    implicit none
    integer, intent(in) :: n
    real(pv), dimension(n, n), intent(out) :: Zero_Matrix

    integer :: i, j

    ! Initialize all elements to zero
    do i = 1, n
        do j = 1, n
            Zero_Matrix(i, j) = 0.0d0
        end do
    end do

end subroutine Create_Zero_Matrix


subroutine RREF(A, B)
    ! General Description: This subroutine computes the reduced row echelon
    ! form (RREF) of a matrix A. The RREF is stored in the matrix B. The
    ! original matrix A is preserved. The RREF is computed using Gaussian
    ! elimination with partial pivoting.

    implicit none
    ! Input Variables:
    real(pv), dimension(:,:), intent(in) :: A
    ! Output Variables:
    real(pv), dimension(:,:), intent(out) :: B

    ! Local Variables:
    integer :: i, j, k, pivot_row, n, m
    real(pv) :: max, temp, factor

    ! Determine the size of the matrix A
    n = SIZE(A, 1) ! Number of rows
    m = SIZE(A, 2) ! Number of columns

    ! Copy A to B at the start to preserve A
    B = A

    ! Loop over each column to perform Gaussian elimination on B
    do k = 1, min(n,m)
        ! Step 1: Find the pivot row
        pivot_row = k
        max = abs(B(k,k))
        do i = k+1, n
            if (abs(B(i,k)) > max) then
                max = abs(B(i,k))
                pivot_row = i
            endif
        enddo

        if (B(pivot_row, k) == 0.0) cycle ! No pivot in this column, move to the next column

        ! Step 2: Swap the pivot row with the current row (if needed)
        if (pivot_row /= k) then
            do j = k, m
                temp = B(k,j)
                B(k,j) = B(pivot_row,j)
                B(pivot_row,j) = temp
            enddo
        endif

        ! Step 3: Make the pivot = 1
        temp = B(k, k)
        do j = k, m
            B(k, j) = B(k, j) / temp
        enddo

        ! Step 4: Eliminate all other entries in the current column
        do i = 1, n
            if (i /= k) then ! Skip the pivot row
                factor = B(i, k)
                do j = k, m
                    B(i, j) = B(i, j) - factor * B(k, j)
                enddo
            endif
        enddo
    enddo
end subroutine RREF


subroutine MatrixRank(A, rank)
    ! This subroutine computes the rank of a matrix A using the reduced row
    ! echelon form (RREF). The rank is returned as an integer. The RREF is
    ! computed using the subroutine RREF.
    !Inputs:
    !   A: The matrix whose rank is to be computed
    !Outputs:
    !   rank: The rank of the matrix A

    implicit none
    real(pv), dimension(:,:), intent(in) :: A
    integer, intent(out) :: rank
    real(pv), dimension(:,:), allocatable :: B
    integer :: i, j, n, m
    logical :: nonZeroRow

    ! Determine the size of the matrix A
    n = SIZE(A, 1) ! Number of rows
    m = SIZE(A, 2) ! Number of columns
    allocate(B(n, m))

    ! Compute the RREF of A and store it in B
    call RREF(A, B)

    rank = 0
    do i = 1, n
        nonZeroRow = .false.
        do j = 1, m
            if (B(i, j) /= 0.0_pv) then
                nonZeroRow = .true.
                exit ! No need to check the rest of the row
            endif
        enddo
        if (nonZeroRow) rank = rank + 1
    enddo

    ! Clean up
    deallocate(B)
end subroutine MatrixRank


subroutine Is_Symmetric(A, result)
    ! This subroutine checks if a square matrix A is symmetric. The result is
    ! returned as a logical value. The matrix A is assumed to be square.
    !Inputs:
    !   A: The square matrix to be checked
    !Outputs:
    !   result: .true. if A is symmetric, .false. otherwise


    implicit none
    real(pv), intent(in) :: A(:,:)
    logical, intent(out) :: result
    real(pv), dimension(size(A, 1), size(A, 1)) :: At
    integer :: i, j, n

    real(pv) :: A_Total_Sum, At_Total_Sum

    !Check to see if A is square:
    if (size(A, 1) /= size(A, 2)) then
        print*, "Error: The matrix is not square."
        result = .false.
        return !Exit the subroutine
    end if

    n = size(A, 1) ! Determine the size of the square matrix

    ! Transpose the matrix A and store the result in At
    call Transpose_Square_Matrix(A, At)

    call Matrix_TOtal_Sum(A, A_Total_Sum)
    call Matrix_Total_Sum(At, At_Total_Sum)

    if (A_Total_Sum == At_Total_Sum) then
        result = .true.
    else
        result = .false.
    end if
end subroutine Is_Symmetric


subroutine Matrix_Total_Sum(A, A_Total_Sum)
    ! Calculates the sum of all elements in a matrix A
    ! Inputs:
    !   A: The matrix whose elements are to be summed
    ! Outputs:
    !   sum: The sum of all elements in A

    implicit none
    real(pv), intent(in) :: A(:,:)
    real(pv), intent(out) :: A_Total_Sum
    integer :: i 

    !Local Variables
    real(pv), dimension(size(A, 1) * size(A, 2)) :: Afr
    
    call Flatten_Matrix_by_Row(A, Afr)

    A_Total_Sum = 0.0d0

    do i = 1, size(Afr)
        A_Total_Sum = A_Total_Sum + Afr(i)
    end do
end subroutine Matrix_Total_Sum


subroutine Flatten_Matrix_by_Row(A, Afr)
    !General Description: This subroutine takes a matrix A and flattens it into a
    !one-dimensional array by row. The flattened array is stored in Afr.
    !Inputs:
    !   A: The matrix to be flattened
    !Outputs:
    !   Afr: The flattened array

    !Example:
    !  
    !  A = [1 2; 3 7] --> Afr = [1 2 3 7]
    !

    implicit none
    real(pv), dimension(:,:), intent(in) :: A
    real(pv), dimension(size(A, 1) * size(A, 2)), intent(out) :: Afr

    integer :: i, j, k

    k = 1
    do i = 1, size(A, 1)
        do j = 1, size(A, 2)
            Afr(k) = A(i, j)
            k = k + 1
        end do
    end do
end subroutine Flatten_Matrix_by_Row


subroutine Flatten_Matrix_by_Column(A,Afc)
    !General Description: This subroutine takes a matrix A and flattens it into a
    !one-dimensional array by column. The flattened array is stored in Afc.
    !Inputs:
    !   A: The matrix to be flattened
    !Outputs:
    !   Afc: The flattened array

    !Example:
    !  
    !  A = [1 2; 3 7] --> Afc = [1 3 2 7]
    !

    implicit none
    real(pv), dimension(:,:), intent(in) :: A
    real(pv), dimension(size(A, 1) * size(A, 2)), intent(out) :: Afc

    integer :: i, j, k

    k = 1
    do j = 1, size(A, 2)
        do i = 1, size(A, 1)
            Afc(k) = A(i, j)
            k = k + 1
        end do
    end do
end subroutine Flatten_Matrix_by_Column



!---------------------------------------------------------------------------------------------------
!                              ** PRINTING TOOLS**

SUBROUTINE print_real_vector(A)
    ! This subroutine prints a one-dimensional array of real numbers
    IMPLICIT NONE

    ! Declare A as a one-dimensional array of unknown size
    REAL(pv), INTENT(IN) :: A(:)
    INTEGER :: i, N

    ! Determine the size of the array A
    N = SIZE(A)

    ! Loop over the array and print each element
    DO i = 1, N
        WRITE(*, '(F0.4)') A(i)
    END DO

END SUBROUTINE print_real_vector


SUBROUTINE print_integer_vector(A)
    INTEGER, INTENT(IN) :: A(:)
    INTEGER :: i
    DO i = 1, SIZE(A)
        WRITE(*, '(I0)') A(i)
    END DO
END SUBROUTINE print_integer_vector


SUBROUTINE print_real_matrix(matrix)
  REAL(pv), INTENT(IN) :: matrix(:,:)
  INTEGER :: i, j

  DO i = 1, SIZE(matrix, 1)
    DO j = 1, SIZE(matrix, 2)
      WRITE(*, '(F8.4, " ")', ADVANCE='NO') matrix(i,j)
    END DO
    WRITE(*, *) ! newline
  END DO
END SUBROUTINE print_real_matrix





end module MOD_LinAlg