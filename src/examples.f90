program examples
    !This program demonstrates the use of the MOD_LinAlg module with various
    !examples. 

    USE MOD_LinAlg
    implicit none
    integer :: i,j,k
    real(pv) :: A(3,3), B(3,3), C(3,3), det
    real(pv) :: A_sym(3,3), A_sym_check(3,3)
    real(pv),allocatable :: Afr(:), Afc(:)

    real(pv) :: A1(2,2), B1(2,1)
    real(pv) :: A2(3,4), B2(3,4)
    real(pv) :: x1(2,1)

    real(pv) :: sum
    integer :: n, status, rank
    logical :: result

    n = size(A,dim = 1)

    !Define Matricies:
    A = reshape((/1,2,12,4,5,6,70,8,9/),(/3,3/))
    A_sym = reshape((/1,2,7,2,5,6,7,6,4/),(/3,3/)) !Symmetric Matrix
    B = reshape((/5,2,7,4,9,0,7,81,20/),(/3,3/))
    C = A + B

    A1 = reshape((/2,7,4,8/),(/2,2/))
    B1 = reshape((/16,44/),(/2,1/))
    
   
    A2 = reshape((/1,2,3,2,5,5,1,1,1,9,18,24/) , (/3,4/))

    allocate(Afr(int(size(A,1)*size(A,2))), Afc(int(size(A,1)*size(A,2))))

    !***************************************************************



    !***************************************************************



    !-------------------------------------------------------------
    write(*,*) "-------------------------------------------------"
    print*, "MATRIC INVERSION EXAMPLE"
    print*, "Print A Matrix before inversion"
    call print_real_matrix(A)
    print*, " " !Spacer
    

    call Invert_Matrix(A, n, status)
    if (status /= 0) then
        write(*,*) "Matrix is singular to working precision."
    else
        print*, "Matrix A After Inversion"
        call print_real_matrix(A); print*, " " !Spacer
    end if

    !-------------------------------------------------------------
    write(*,*) "-------------------------------------------------"
    print*, "SYMMETRIC MATRIX EXAMPLE"; print*, " " !Spacer

    print*, "Print Matrix that needs to be ckecked for symmetry"
    call print_real_matrix(A_sym)

    call Is_Symmetric(A_sym, result)

    if(result) then
        write(*,*) "Result = Matrix is symmetric."
    else
        write(*,*) "Result = Matrix is not symmetric."
    end if

    print*, " " !Spacer

    !-------------------------------------------------------------
    write(*,*) "-------------------------------------------------"
    print*, "DETERMINANT EXAMPLE"
    print*, " " !Spacer
    print*, "Calculate Determinant of C"
    call Matrix_Det(C, det, status)
    if(status /= 0) then
        write(*,*) "Matrix Determinant cannot be determined."
    else
        write(*,"(A,F12.6)") "Determinant of C = ", det
    end if


    !-------------------------------------------------------------
    write(*,*) "-------------------------------------------------"
    print*, "SOLVE A LINEAR SYSTEM EXAMPLE"
    print*, " " !Spacer
    print *, "The solution to the linear system A*x = B, where"
    print *, "A = [2 4]"
    print *, "    [7 8]"
    print *, "and B = [16]"
    print *, "        [44]"
    print *, "The solution is x1 = 4, and x2 = 2"; print*, " " !Spacer


        !LSE_Solve(A, B, x, N, status)
    call LSE_Solve(A1, B1, x1, 2, status)

    if(status /= 0) then
        write(*,*) "Matrix is singular to working precision."
    else
        print*, "Print the solution matrix:"
        call print_real_matrix(x1)
    end if
    print*, " " !Spacer

    !-------------------------------------------------------------
    write(*,*) "-------------------------------------------------"
    print*, "REDUCED ROW ECHELON FORM EXAMPLE"
    print*, " " !Spacer

    print*, "Print A2 Matrix before RREF"
    call print_real_matrix(A2)

    call RREF(A2, B2)

    print*, "Print B2 = rref(A2)"
    call print_real_matrix(B2)
    print*, " " !Spacer

    !-------------------------------------------------------------
    write(*,*) "-------------------------------------------------"
    print*, "MATRIX RANK EXAMPLE"
    print*, " " !Spacer

    call MatrixRank(A2, rank)
    print*, "Rank of A2 = ", rank

    


end program examples

