#include <chrono>
#include <thread>
#include <unistd.h>
#include "Vec_Tools.hpp"
#include "Math_Tools.hpp"
#include "Sys_Tools.hpp"
#include "Vector_3.hpp"
#include "Tensor4_3D.hpp"
#include "SymmTensor4_3D.hpp"
#include "IEN_FEM.hpp"
#include "NodalBC.hpp"
#include "NodalBC_3D_inflow.hpp"
#include "ElemBC_3D.hpp"
#include "ElemBC_3D_outflow.hpp"

int main(int argc, char *argv[])
{
    // ----------------
    // Constructor Test
    //-----------------
    std::cout << "------------------" << std::endl;
    std::cout << "Constructor Test: " << std::endl;
    std::cout << "------------------" << std::endl;

    MATH_T::Matrix_Dense<4> Mat1;   // Constructor

    MATH_T::Matrix_SymPos_Dense<4> Mat_SymPos1;   // Constructor 

    MATH_T::Matrix_Dense<4> Mat2(Mat1);      // Default Copy Constructor, No print

    MATH_T::Matrix_SymPos_Dense<4> Mat_SymPos2(Mat2);  // Copy Constructor 

    MATH_T::Matrix_SymPos_Dense<4> Mat_SymPos3(Mat_SymPos1);  // Copy Constructor 

    std::cout << "Mat1: " << std::endl;
    Mat1.print_info();   

    std::cout << "Mat2: " << std::endl;
    Mat2.print_info();    

    std::cout << "Mat_SymPos1: " << std::endl;
    Mat_SymPos1.print_info();  

    std::cout << "Mat_SymPos2: " << std::endl;
    Mat_SymPos2.print_info();  

    std::cout << "Mat_SymPos3: " << std::endl;
    Mat_SymPos3.print_info();   


    //----------------------------
    // Initialization Test
    //----------------------------
    std::cout << "---------------------" << std::endl;
    std::cout << "Initialization Test: " << std::endl;
    std::cout << "---------------------" << std::endl;

    std::array<double, 3*3> AA{1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0};   

    std::array<double, 3*3> BB{1.0, 2.0, 3.0, 2.0, 5.0, 6.0, 3.0, 6.0, 7.0};

    MATH_T::Matrix_Dense<3> Mat3(AA);   // Constructor

    MATH_T::Matrix_SymPos_Dense<3> Mat_SymPos4(BB);   // Constructor 

    MATH_T::Matrix_Dense<3> Mat4;      

    MATH_T::Matrix_SymPos_Dense<3> Mat_SymPos5;  

    std::cout << "Mat3: " << std::endl;
    Mat3.print_info();    

    std::cout << "Mat4: " << std::endl;
    Mat4.print_info();    

    std::cout << "Mat_SymPos4: " << std::endl;
    Mat_SymPos4.print_info();   

    std::cout << "Mat_SymPos5: " << std::endl;
    Mat_SymPos5.print_info();   

    // --------------------------------------
    // Assignment operator= Test
    // --------------------------------------
    std::cout << "---------------------------" << std::endl;
    std::cout << "Assignment operator= Test: " << std::endl; 
    std::cout << "---------------------------" << std::endl;

    Mat4 = Mat3;   // Assignment operator =

    std::cout << "Mat4 = Mat3, Mat4: " << std::endl;
    Mat4.print_info();    

    Mat4 = Mat_SymPos4;    // Assigment operator = 

    std::cout << "Mat4 = Mat_SymPos4, Mat4: " << std::endl;
    Mat4.print_info();    
   
    Mat_SymPos5 = Mat_SymPos4;   // Assigment operator =

    std::cout << "Mat_SymPos5 = Mat_SymPos4, Mat_SymPos5: " << std::endl;
    Mat_SymPos5.print_info();   

    // --------------------------
    // Function check_symm() Test
    // --------------------------
    std::cout << "----------------------------" << std::endl;
    std::cout << "Function check_symm() Test: " <<  std::endl;
    std::cout << "----------------------------" << std::endl;

    std::array<double, 3*3> CC{1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0};  

    MATH_T::Matrix_SymPos_Dense<3> Mat7{CC};  

    Mat7.check_symm();

    // ----------------------------
    // Function gen_rand() Test
    // Function gen_rand() of class Matrix_Dense gives a random regular matrix
    // Function gen_rand() of class Matrix_SymPos_Dense gives a random symmetric matrix, which is not positive definite
    //-----------------------------
    std::cout << "--------------------------" << std::endl;
    std::cout << "Function gen_rand() Test: " << std::endl;
    std::cout << "--------------------------" << std::endl;

    Mat1.gen_rand();          

    Mat2.gen_rand(); 

    Mat_SymPos1.gen_rand();

    Mat_SymPos2.gen_rand();

    Mat_SymPos3.gen_rand();

    std::cout << "Mat1: " << std::endl;
    Mat1.print_info();    

    std::cout << "Mat2: " << std::endl;
    Mat2.print_info();   

    std::cout << "Mat_SymPos1: " << std::endl;
    std::cout << "check_symm: " << std::endl;
    Mat_SymPos1.check_symm();
    Mat_SymPos1.print_info();   

    std::cout << "Mat_SymPos2: " << std::endl;
    std::cout << "check_symm: " << std::endl;
    Mat_SymPos2.check_symm();
    Mat_SymPos2.print_info();   

    std::cout << "Mat_SymPos3: " << std::endl;
    std::cout << "check_symm: " << std::endl;
    Mat_SymPos3.check_symm();
    Mat_SymPos3.print_info();   

    //-----------------------------------------
    // Assignment operator= Test via function Mult() of class Matrix_Dense
    //-----------------------------------------
    std::cout << "---------------------------------------------------------------------" << std::endl;
    std::cout << "Assignment operator= Test via function Mult() of class Matrix_Dense: " << std::endl;
    std::cout << "---------------------------------------------------------------------" << std::endl;

    MATH_T::Matrix_Dense<3> Mat5;      

    MATH_T::Matrix_SymPos_Dense<3> Mat_SymPos6; 

    MATH_T::Matrix_Dense<3> Mat6; 

    MATH_T::Matrix_Dense<3> Mat, MM, NN; 

    for (int ii = 0; ii < 1; ++ii)
    {
        Mat5.gen_rand();
        Mat.gen_rand();
        Mat_SymPos6 = Mat5;   // Copy Constructor -> Assignment operator =, check symmetricity

        std::cout << "Mat_SymPos6: " << std::endl;
        Mat_SymPos6.print_info();   

        std::cout << "Mat5: " << std::endl;
        Mat5.print_info();   

        std::cout << "Mat: " << std::endl;
        Mat.print_info(); 

        // -----------------------------------

        Mat6 = Mat_SymPos6;

        MM.Mult(Mat5, Mat);
        NN.Mult(Mat6, Mat);

        std::cout << "Mat6: " << std::endl;
        Mat6.print_info();   

        std::cout << "MM: " << std::endl;
        MM.print_info();  

        std::cout << "NN: " << std::endl;
        NN.print_info();  

        // ----------------------------------

        for (int jj = 0; jj < 3 * 3; ++jj)
        {
            if( !MATH_T::equals( MM(jj), NN(jj), 1.0e-15) )
                std::cout<<"error: Matrix MM("<<jj<<") does not match NN("<<jj<<"). \n";
        }
    }

    // ------------------------
    // Function Mult(array<N>) of class Matrix_SymPos_Dense Test
    // Ax + Bx = (A + B)x
    // A(Bx) = (AB)x, (AB): convert Matrix_SymPos_Dense to Matrix_Dense first
    // ------------------------
    std::cout << "-----------------------------------------------------------" << std::endl;
    std::cout << "Function Mult(array<N>) of class Matrix_SymPos_Dense Test: " << std::endl;
    std::cout << "-----------------------------------------------------------" << std::endl;

    MATH_T::Matrix_SymPos_Dense<3> Mat_SymPos7;

    MATH_T::Matrix_SymPos_Dense<3> Mat_SymPos8;

    std::array<double, 3> DD{1.0, 2.0, 3.0};

    std::array<double, 3> EE1, EE2, EE;

    std::array<double, 3> FF;

    Mat_SymPos7.gen_rand();

    Mat_SymPos8.gen_rand();

    // Ax + Bx = (A + B)x
    std::cout << "Ax + Bx = (A + B)x: " << std::endl;

    MATH_T::Matrix_SymPos_Dense<3> Mat_SymPos9;

    EE1 = Mat_SymPos7.Mult(DD);

    EE2 = Mat_SymPos8.Mult(DD);
    
    for (int ii = 0; ii < 3; ++ii)
      EE[ii] = EE1[ii] + EE2[ii];

    for (int ii = 0; ii < 3*3; ++ii)
      Mat_SymPos9(ii) = Mat_SymPos7(ii) + Mat_SymPos8(ii);
    
    FF = Mat_SymPos9.Mult(DD);

    for (int ii = 0; ii < 3; ++ii)
    {
      if( !MATH_T::equals( EE[ii], FF[ii], 1.0e-15))
        std::cout<<"error: EE["<<ii<<"] does not match FF["<<ii<<"]. \n";
    }

    // A(Bx) = (AB)x
    std::cout << "A(Bx) = (AB)x: " << std::endl;
    
    Mat_SymPos7.gen_rand();

    Mat_SymPos8.gen_rand();

    std::array<double, 3> GG1, GG;

    std::array<double, 3> HH;

    GG1 = Mat_SymPos8.Mult(DD);

    GG = Mat_SymPos7.Mult(GG1);

    MATH_T::Matrix_Dense<3> Mat8;

    MATH_T::Matrix_Dense<3> Mat9;

    MATH_T::Matrix_Dense<3> Mat10;

    MATH_T::Matrix_SymPos_Dense<3> Mat_SymPos10;

    Mat8 = Mat_SymPos7;
    
    Mat9 = Mat_SymPos8;

    Mat10.Mult(Mat8, Mat9);

    Mat_SymPos10 = Mat10;

    HH = Mat_SymPos10.Mult(DD);

    for (int ii = 0; ii < 3; ++ii)
    {
      if( !MATH_T::equals( GG[ii], HH[ii], 1.0e-15))
        std::cout<<"error: GG["<<ii<<"] does not match HH["<<ii<<"]. \n";
    }

    // -------------------------
    // Functions LDLT_fac() & LDLT_solve() Test
    // Ax = b  ->  A.LDLT_fac()  ->  x_a = A.LDLT_solve(b)  ->  ||x_a - x|| < tol
    // -------------------------
    std::cout << "------------------------------------------" << std::endl;
    std::cout << "Functions LDLT_fac() & LDLT_solve() Test: " << std::endl;
    std::cout << "------------------------------------------" << std::endl;

    std::array<double, 10> RHS{5.0, 1.0, 3.0, 2.0, 8.0, 6.0, 4.0, 7.0, 9.0, 10.0};

    MATH_T::Matrix_SymPos_Dense<10> Mat_SymPos11;
    MATH_T::Matrix_SymPos_Dense<10> LHS;

    LHS.gen_rand(-10, 10);
    Mat_SymPos11 = LHS;
    std::cout << "LHS: " << std::endl;
    LHS.print_info();  


    LHS.LDLt_fac();
    std::cout << "LDLT(LHS): " << std::endl;
    LHS.print_info();  

    std::array<double, 10> sol;

    sol = LHS.LDLt_solve(RHS);
    std::cout << "sol:" << std::endl;
    for (int ii = 0; ii < 3; ++ii) std::cout << sol[ii] << "\t";
    std::cout << std::endl;

    std::array<double, 10> RHS_app;

    RHS_app = Mat_SymPos11.Mult(sol);

    std::cout << "RHS:" << std::endl;
    for (int ii = 0; ii < 10; ++ii) std::cout << RHS[ii] << "\t";
    std::cout << std::endl;

    std::cout << "RHS_app:" << std::endl;
    for (int ii = 0; ii < 10; ++ii) std::cout << RHS_app[ii] << "\t";
    std::cout << std::endl;

    for (int ii = 0; ii < 10; ++ii)
    {
      if( !MATH_T::equals( RHS[ii], RHS_app[ii], 1.0e-14))
        std::cout<<"error: RHS["<<ii<<"] does not match RHS_app["<<ii<<"]. \n";
    }

  return EXIT_SUCCESS;
}

// EOF
