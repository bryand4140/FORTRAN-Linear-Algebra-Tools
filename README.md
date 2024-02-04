# FORTRAN Linear Algebra Tools

This repository contains a collection of Modern Fortran subroutines for performing various linear algebra operations. The code is intended to offer a robust alternative to LAPACK, and other older repositories while still offers the same level of performance but in a way that is much easier to use and understand. 

## Features
Subroutines Overview
- **Linear Equation Solving**: Solve systems of linear equations using the subroutine `LSE_Solve`.
- **Matrix Inversion**: Invert square matrices with `Invert_Matrix`.
- **Determinant Calculation**: Compute the determinant of square matrices using `Matrix_Det`.
- **LU Decomposition**: Perform LU decomposition with partial pivoting via `LU_Decomposition`.
- **Matrix Transposition and Creation**: Transpose square matrices with `Transpose_Square_Matrix` and create identity or zero matrices using `Create_Identity_Matrix` and `Create_Zero_Matrix`.
- **RREF and Matrix Rank**: Calculate the reduced row echelon form and determine the rank of matrices with `RREF` and `MatrixRank`.
- **Symmetry Checking**: Check if a matrix is symmetric with `Is_Symmetric`.
- **Element Summation**: Sum all elements in a matrix using `Matrix_Total_Sum`.

## How to Use
Clone the repository to your local machine.
Include the desired subroutine(s) in your Fortran project.
Follow the input and output specifications detailed in the comments within each subroutine for integration.

## Contributing
Contributions are welcome! If you have improvements or additional functionalities you'd like to add, please fork the repository and submit a pull request.

## License
This project is open-sourced under the GNU general public license. See the LICENSE file for more details.

## Support
If you encounter any issues or have questions about implementing these subroutines, please open an issue on GitHub.

## Connect
Stay updated on this project by starring or watching this repository. Your support helps increase visibility and encourages further development.
