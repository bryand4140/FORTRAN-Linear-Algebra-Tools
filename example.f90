program example
    USE MOD_LinAlg
    implicit none

    real(pv) :: A(3,3), B(3,3), C(3,3), det
    real(pv) :: A1(2,2), B1(2,1)
    integer :: n, status

    n = size(A,dim = 1)

    A = reshape((/1,2,12,4,5,6,70,8,9/),(/3,3/))
    B = reshape((/5,2,7,4,9,0,7,81,20/),(/3,3/))
    C = A + B

    print*, "Print A Matrix"
    call print_real_matrix(A)
    print*, " " !Spacer
    

    call Invert_Matrix(A, n, status)
    if (status /= 0) then
        write(*,*) "Matrix is singular to working precision."
    else
        print*, "Print Inverted A Matrix"
        call print_real_matrix(A)
    end if
    
    

    print*, " " !Spacer
    print*, "Print B Matrix"
    call print_real_matrix(B)


    print*, " " !Spacer
    print*, "Calculate Determinant of C"
    call Matrix_Det(C, n, det, status)
    if(status /= 0) then
        write(*,*) "Matrix Determinant cannot be determined."
    else
        write(*,"(A,F12.6)") "Determinant of C = ", det
    end if

    print*, " " !Spacer
    print*, "Solve a linear system:"
    !Solve the linear system A*x = B where 
    ! A = [2 4]      B = [16]
    !     [7 8]          [44]
    !For which the solution is x = 4, and y = 2

    A1 = reshape((/2,7,4,8/),(/2,2/))
    B1 = reshape((/16,44/),(/2,1/))

    call LSE_Solve(A1, B1, 2, status)

   
    if(status /= 0) then
        write(*,*) "Matrix is singular to working precision."
    else
        print*, "Print Solution to Linear System"
        call print_real_matrix(B1)
    end if


end program example

