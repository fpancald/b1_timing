#include "b1_timing.hpp"

int main()
{
    std::chrono::time_point<std::chrono::high_resolution_clock> ct1_b1, ct2_b1, ct1_b2, ct2_b2, for_start, for_end;//initiate timing parameters
    std::chrono::duration<double> elapsed_seconds_b1, elapsed_seconds_b2, for_laps;

    FILE *fpA_b1, *fpx_b1, *fpb_b1, *fpA_b2, *fpx_b2, *fpb_b2;//initiate file pointers

		int prec=1024, size=10, iter=1;//test parameters: prec= bit precision, size=square matrix size, iter=number of repetition

    double tol=5e-324, largeChange=1e16;//parameters for some solvers in b1
    mpf_t mpftol,  mpflargeChange;

    mpf_init2(mpftol,prec);//different initialization of parameters for some solvers in b1 /*needs to be checked*/
    mpf_init2(mpflargeChange,prec);
    mpf_set_d(mpftol,tol);
    mpf_set_d(mpflargeChange,largeChange);
    double cond_num,  norm_A,  norm_A_inv;//auxiliary variables for some solvers in b1

    uint usize=(int)size;//convert size to unsigned integer /*needed for RandomOfUnits*/

    bertini::mpfr_float::default_precision(prec);//set mp precision for b2

    bertini::Vec<bertini::complex> x_b2;//initialize solution vector in b2

    mat_mp A_b1;//initialize and set precision for b1 matrix and vectors /*comment if generating multiple matrices, and uncomment in for-loop*/ /*needs better solution*/
    vec_mp x_b1, b_b1;
    init_mat_mp2(A_b1,size,size,prec);
    init_vec_mp2(x_b1,size,prec);
    init_vec_mp2(b_b1,size,prec);

    bool DontAlignCols=false;//IO parameters for Eigen IOFormat
    int sprec=prec;

    Eigen::IOFormat CommaInitFmt(sprec, DontAlignCols, ", ", ", ", "", "", " << ", ";");//different Eigen IOFormat
    Eigen::IOFormat CleanFmt(4, 0, ", ", "\n", "[", "]");
    Eigen::IOFormat OctaveFmt(sprec, 0, ", ", ";\n", "", "", "[", "]");
    Eigen::IOFormat MatlabFmtInline(sprec, 0, ", ", ";", "", "", "[", "]");//Matlab format one line
    Eigen::IOFormat HeavyFmt(sprec, 0, ", ", ";\n", "[", "]", "[", "]");

    for_start=std::chrono::system_clock::now();//start count for overall time

	for (int i=0;i<iter;i++)//start iterations
		{
		make_matrix_random_mp(A_b1,size, size, prec);//define random matrix and vectors in b1 /*Note values are in the unit circle and the Matrix is orthonormal*/
		make_vec_random_mp(b_b1,size);

    bertini::Mat<bertini::complex> A_b2=bertini::RandomOfUnits<bertini::complex>(usize,usize);//define random matrix and vectors in b2 /*Note values are in the unit circle but the Matrix is NOT orthonormal*/
    bertini::Vec<bertini::complex> b_b2=bertini::RandomOfUnits<bertini::complex>(usize);

    bertini::Mat<bertini::complex> Q=A_b2.fullPivHouseholderQr().matrixQ();

    fpA_b1 = fopen ("matrixA_b1.txt" , "a+");//write to file A and b for b1
		printMat_Matlab_mp(fpA_b1, prec, A_b1);
    fclose(fpA_b1);
    fpb_b1 = fopen ("vectorb_b1.txt" , "a+");
		printVec_Matlab_mp(fpb_b1, prec, b_b1);
    fclose(fpb_b1);

    std::stringstream buffer;//write to file A and b for b2
    // buffer << A_b2.real().format(MatlabFmtInline)<<"+i*"<<A_b2.imag().format(MatlabFmtInline) << std::endl;
    // fpA_b2 = fopen ("matrixA_b2.txt" , "a+");
    // fprintf(fpA_b2,buffer.str().c_str());
    // buffer.str( std::string() );
    // buffer.clear();
    buffer << b_b2.real().format(MatlabFmtInline)<<"+i*"<<b_b2.imag().format(MatlabFmtInline) << std::endl;
    fpb_b2 = fopen ("vectorb_b2.txt" , "a+");
    fprintf(fpb_b2,buffer.str().c_str());
    buffer.str( std::string() );
    buffer.clear();

    buffer << Q.real().format(MatlabFmtInline)<<"+i*"<<Q.imag().format(MatlabFmtInline) << std::endl;
    fpA_b2 = fopen ("matrixA_b2.txt" , "a+");
    fprintf(fpA_b2,buffer.str().c_str());
    buffer.str( std::string() );
    buffer.clear();

    ct1_b1=std::chrono::system_clock::now();//start timing of one iteration for b1

    // matrixSolve_mp( x_b1, A_b1, b_b1);//b1 solver
    // matrixSolve_cond_num_norms_mp( x_b1,  A_b1,  b_b1,  &cond_num,  &norm_A,  &norm_A_inv);
    matrixSolve_LU_mp( x_b1,  A_b1,  b_b1,  tol,  largeChange);
    // matrixSolve_svd_Least_Squares_mp( x_b1,  A_b1,  b_b1,  tol,  largeChange,  &cond_num,  &norm_A,  &norm_A_inv);  //least square svd
    // matrixSolve_svd_Least_Squares_mp2( x_b1,  A_b1,  b_b1,  mpftol,  mpflargeChange,  &cond_num,  &norm_A,  &norm_A_inv);
    // matrixSolve_Least_Squares_mp( x_b1,  A_b1,  b_b1,  tol,  largeChange,  &cond_num,  &norm_A,  &norm_A_inv);  // Least square QR
    // matrixSolve_Least_Squares_mp2( x_b1,  A_b1,  b_b1,  mpftol,  mpflargeChange,  &cond_num,  &norm_A,  &norm_A_inv);
    // matrixSolve_Hessenberg_Least_Squares_mp( x_b1,  A_b1,  b_b1, mpftol, mpflargeChange);

    ct2_b1=std::chrono::system_clock::now();//end timing of one iteration for b1
    elapsed_seconds_b1 += ct2_b1-ct1_b1;//update total solving timing for b1

    ct1_b2=std::chrono::system_clock::now();//start timing of one iteration for b2

    // x_b2 = A_b2.colPivHouseholderQr().solve(b_b2);//b2solver
    // x_b2 = A_b2.householderQr().solve(b_b2);
    // x_b2 = A_b2.partialPivLu().solve(b_b2);
    // x_b2 = A_b2.jacobiSvd().solve(b_b2);
    
    x_b2 = Q.partialPivLu().solve(b_b2);


    ct2_b2=std::chrono::system_clock::now();//end timing of one iteration for b2
    elapsed_seconds_b2 += ct2_b2-ct1_b2;//update total solving timing for b2

    fpx_b1 = fopen ("vectorx_b1.txt" , "a+");//write to file solution x for b1
		printVec_Matlab_mp(fpx_b1, prec, x_b1);
    fclose(fpx_b1);

    buffer << x_b2.real().format(MatlabFmtInline)<<"+i*"<<x_b2.imag().format(MatlabFmtInline) << std::endl;
    fpx_b2 = fopen ("vectorx_b2.txt" , "a+");
    fprintf(fpx_b2,buffer.str().c_str());
    buffer.str( std::string() );
    buffer.clear();
    }
    for_end=std::chrono::system_clock::now();//end total timing of for-loop
    for_laps=for_end-for_start;//calculate total for-loop time

    std::cout << "elapsed time for linear algebra operations in b1: " << elapsed_seconds_b1.count() << "s\n";//To screen output /*TODO change to write on file providing relevant parameters as well as timing*/
    std::cout << "elapsed time for linear algebra operations in b2: " << elapsed_seconds_b2.count() << "s\n";
    std::cout << "elapsed time during for-loop: " << for_laps.count() << "s\n";

    clear_mat_mp(A_b1);//free memory for b1 vectors and matrix /*Not needed in Eigen?*/
    clear_vec_mp(x_b1);
    clear_vec_mp(b_b1);

    mpf_clear(mpftol);//clear mpfr parameters
    mpf_clear(mpflargeChange);

    return 0;

}
