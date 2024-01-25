module MOD_LinAlg
    implicit none


contains
!---------------------------------------------------------------------------------------------------
!                              ** LINEAR ALGEBRA TOOLS**
subroutine LSE_Solve(A, B, N, status)
    ! Computes the solution to a real system of linear equations AX = B, where 
    ! A is a square nxn matrix, B is the coefficient matrix and X are the 
    ! unknowns. The coefficient matrix is modified to store the solutions of
    ! the system of equations. In other words, when the program ends B = X.
    ! Note that A must be a square matrix. 

    !Checked for accuracy on 9-18-2023
    
    implicit none
    integer, intent(in) :: N                    !Size of the system
    real(pv), dimension(N,N), intent(inout) :: A !Main matrix
    real(pv), dimension(N),   intent(inout) :: B !Coefficient Matrix
    integer, intent(out) :: status !0 for sucessful solve, 1 = matrix could not be inverted.

    call Invert_Matrix(A, N, status)
    if (status == 0) then
        B = matmul(A, B)
    end if
end subroutine LSE_Solve


subroutine Invert_Matrix(A, n, status)
    !Inverts a square nxn matrix A
    !Checked for accuracy on 9-18-2023
    implicit none
    integer, intent(in) :: n !The size of the square nxn matrix
    real(pv), dimension(n, n), intent(inout) :: A
    real(pv), dimension(n, n) :: A_check
    real(pv), dimension(n, n) :: originalA
    real(pv), dimension(n, n) :: invA
    real(pv), dimension(n) :: B, X
    integer, dimension(n) :: P
    integer :: i, j, k, maxrow
    integer, intent(out) :: status !0 for successful inversion, 1 for an unsuccessful inversion.
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


subroutine Matrix_Det(A, n, det, status)
    !Computes the determinant of a square nxn matrix.
    implicit none
    integer, intent(in) :: n
    real(pv), dimension(n, n), intent(in) :: A
    real(pv), intent(out) :: det
    integer, intent(out) :: status
    real(pv), dimension(n, n) :: LU_A
    integer :: i, swaps

    ! Copy A to LU_A
    LU_A = A

    ! Perform LU decomposition
    call LU_Decomposition(LU_A, n, swaps, status)

    ! Check for singular matrix
    if (status == 1) then
        det = 0.0_8
        return
    end if

    ! Calculate the determinant
    det = 1.0_8
    do i = 1, n
        det = det * LU_A(i, i)
    end do

    ! Adjust for row swaps
    if (mod(swaps, 2) == 1) then
        det = -det
    end if
    
end subroutine Matrix_Det


subroutine LU_Decomposition(A, n, swaps, status)
    !Performs LU Decomp on an nxn matrix. Checked for accuracy on 9-18-2023

    implicit none
    integer, intent(in) :: n !Size of matrix A
    real(pv), dimension(n, n), intent(inout) :: A
    integer, intent(out) :: swaps, status
    integer, dimension(n) :: P
    integer :: i, j, k, maxrow

    ! Initialize status and swaps
    status = 0
    swaps = 0

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
end subroutine LU_Decomposition


SUBROUTINE Transpose_Square_Matrix(A, At, n)
    !Performs a matrix transpose on A and returns At. Matrix must be square.
    !Checked for accuracy on 9-19-2023
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n
    real(pv), INTENT(IN) :: A(n, n)
    real(pv), INTENT(OUT) :: At(n, n)
    INTEGER :: i, j

    DO i = 1, n
        DO j = 1, n
            At(j, i) = A(i, j)
        END DO
    END DO
END SUBROUTINE Transpose_Square_Matrix


subroutine Create_Identity_Matrix(n, Identity)
    !Creates an nxn identity matrix. Checked for accuracy on 9-19-2023
    implicit none
    integer, intent(in) :: n
    real(pv), dimension(n, n), intent(out) :: Identity

    integer :: i

    ! Initialize all elements to zero
    I = 0.0_8

    ! Set the diagonal elements to one
    do i = 1, n
        Identity(i, i) = 1.0_8
    end do

end subroutine Create_Identity_Matrix



end module MOD_LinAlg