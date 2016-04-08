#include "b1_timing.hpp"

int main()
{
    clock_t t1,t2;
    t1=clock();
	for (int i=0;i<1;i++)
		//vec_mp vec;
		{
		mat_mp temp;
		FILE *fp;
		int prec=16384;
		init_mat_mp2(temp,0,0,prec);
		make_matrix_random_mp(temp,10, 10, prec);
		fp = fopen ("myfile.txt" , "w");
		printMat_Matlab_mp(fp, prec, temp);
		//printMat_mp(fp, prec, temp);
		fclose(fp);
		clear_mat_mp(temp);
		}

		/*{vec_mp vec;
		init_vec_mp(10,50);
		make_vec_random_mp(vec, 10);
		}*/
		/*{
	   	 std::cout << "hi francesco\t\t";
		}*/
    t2=clock();
    double diff ((double)t2-(double)t1);
    double diffs = diff / CLOCKS_PER_SEC;
    std::cout<<diffs<<"\tseconds\n";
    
    return 0;
    
}

