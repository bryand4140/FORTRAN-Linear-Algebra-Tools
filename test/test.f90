program test
    use MOD_LinAlg
    implicit none

    REAL(pv), DIMENSION(3,3) :: A
    real(pv) :: Q(3,3), R(3,3)
    real(pv) :: diag(size(A,1))
    real(pv) :: eigenvalues(size(A,1))
    real(pv) :: eigenvectors(size(A,1), size(A,2))
    real(pv) :: tol
    integer :: status, maxIter

    A = RESHAPE((/1, 2, 12, 4, 5, 6, 70, 8, 9/), (/3,3/))
    tol = 1.0e-6
    maxIter = 100 

    CALL QR_Decomposition(A, Q, R, 3, status)

    if(status == 0) then
        print*, "Print Q"
        CALL print_real_matrix(Q)
        print*, "" ! newline

        print*, "Print R"
        CALL print_real_matrix(R)
        print*, "" ! newline

        print*, "Print Q*R"
        CALL print_real_matrix(matmul(Q,R))
    elseif(status == 1) then
        print*, "Matrix A must be square"
    else
        print*, "Matrix is singular or nearly singular"
    end if


    print*, " " !Spacer
    call Find_Real_Eigenvalues(A, eigenvalues, size(A,1), maxIter, tol)
    print*, "Eigenvalues"
    CALL print_real_vector(eigenvalues)



    call Find_Real_Eigenvectors(A, eigenvalues, eigenvectors, 3, status)
    if (status == 0) then
        print*, "Eigenvectors"
        CALL print_real_matrix(eigenvectors)
    else
        print*, "Failed to find eigenvectors"
    end if


    contains

SUBROUTINE Find_Real_Eigenvectors(A, eigenvalues, eigenvectors, n, status)
    !Status: Not Complete

    ! General Description: This subroutine finds eigenvectors for given 
    ! eigenvalues of matrix A. For each eigenvalue, it solves (A - lambda*I)x = 0, 
    ! storing each resultant eigenvector in the columns of the output matrix 
    ! eigenvectors. 

    !Stutus Codes:
    ! 0 - Success
    ! 1 - Failure

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n       ! Dimension of A and number of eigenvalues
    REAL(pv), DIMENSION(n,n), INTENT(IN) :: A      ! Input matrix A
    REAL(pv), DIMENSION(n), INTENT(IN) :: eigenvalues ! Array of eigenvalues
    REAL(pv), DIMENSION(n, n), INTENT(OUT) :: eigenvectors ! Matrix of eigenvectors
    INTEGER, INTENT(OUT) :: status       ! Output status (0 for success, non-zero for failure)
    
    ! Local Variables
    real(pv), DIMENSION(n,n) :: Identity ! Identity matrix
    REAL(pv), DIMENSION(n,n) :: A_prime  ! Modified matrix A - lambda*I
    REAL(pv), DIMENSION(n) :: B          ! Temporary vector for solving linear systems
    REAL(pv), DIMENSION(n) :: x          ! Solution vector for linear systems
    INTEGER :: i, j, k, solve_status

    status = 0

    if(size(A,1) /= size(A,2)) then
        print*, "ERROR: Matrix A must be square"
        status = 1
        return
    end if

    ! Initialize vector B to zeros for each eigenvalue
    do i = 1,size(B,1)
        B(i) = 0.0d0
    end do

    call Create_Identity_Matrix(size(A,1), Identity)

    DO i = 1, size(eigenvalues,1)
        A_prime = A - eigenvalues(i)*Identity

        call print_real_vector(B); print*, " "

        ! Solve the linear system (A - lambda*I)x = 0 for the current eigenvalue
        CALL LSE_Solve(A_prime, B, x, n, solve_status)

        IF (solve_status /= 0) THEN
            PRINT*, "Failed to solve for eigenvector for eigenvalue ", i, ". Status: ", solve_status
            status = 1
            RETURN
        ELSE
            ! Store the solution in the corresponding column of eigenvectors matrix
            DO j = 1, size(A,1)
                eigenvectors(j,i) = x(j)
            END DO
        END IF
    END DO
END SUBROUTINE Find_Real_Eigenvectors









end program test