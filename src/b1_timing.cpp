#include "b1_timing.hpp"

int main()
{
    std::chrono::time_point<std::chrono::high_resolution_clock> ct1, ct2;
    std::chrono::duration<double> elapsed_seconds;
    FILE *fpA, *fpx, *fpb;
		int prec=1024, size=100, iter=10000;
    double tol=1e-16, largeChange=1e16;
    mpf_t mpftol,  mpflargeChange;
    mpf_init2(mpftol,prec);
    mpf_init2(mpflargeChange,prec);
    mpf_set_d(mpftol,tol);
    mpf_set_d(mpflargeChange,largeChange);
    double cond_num,  norm_A,  norm_A_inv;
    double diff;


	for (int i=0;i<iter;i++)
		//vec_mp vec;
		{
		mat_mp A;
    vec_mp x, b;
		init_mat_mp2(A,size,size,prec);
    init_vec_mp2(x,size,prec);
    init_vec_mp2(b,size,prec);

		make_matrix_random_mp(A,size, size, prec);
		make_vec_random_mp(b,size);

		fpA = fopen ("matrixA.txt" , "a+");
		printMat_Matlab_mp(fpA, prec, A);
    fclose(fpA);
    fpb = fopen ("vectorb.txt" , "a+");
		printVec_Matlab_mp(fpb, prec, b);
    fclose(fpb);

    ct1=std::chrono::system_clock::now();

    matrixSolve_mp( x, A, b);
    // matrixSolve_cond_num_norms_mp( x,  A,  b,  &cond_num,  &norm_A,  &norm_A_inv);
    // matrixSolve_LU_mp( x,  A,  b,  tol,  largeChange);
    // matrixSolve_svd_Least_Squares_mp( x,  A,  b,  tol,  largeChange,  &cond_num,  &norm_A,  &norm_A_inv);  //least square svd
    // matrixSolve_svd_Least_Squares_mp2( x,  A,  b,  mpftol,  mpflargeChange,  &cond_num,  &norm_A,  &norm_A_inv);
    // matrixSolve_Least_Squares_mp( x,  A,  b,  tol,  largeChange,  &cond_num,  &norm_A,  &norm_A_inv);  // Least square QR
    // matrixSolve_Least_Squares_mp2( x,  A,  b,  mpftol,  mpflargeChange,  &cond_num,  &norm_A,  &norm_A_inv);
    // matrixSolve_Hessenberg_Least_Squares_mp( x,  A,  b, mpftol, mpflargeChange);

    ct2=std::chrono::system_clock::now();
    elapsed_seconds += ct2-ct1;

    fpx = fopen ("vectorx.txt" , "a+");
		printVec_Matlab_mp(fpx, prec, x);
    fclose(fpx);
		clear_mat_mp(A);
    clear_vec_mp(x);
    clear_vec_mp(b);
		}

    std::cout << "elapsed time: " << elapsed_seconds.count() << "s\n";

    mpf_clear(mpftol);
    mpf_clear(mpflargeChange);

    return 0;

}
