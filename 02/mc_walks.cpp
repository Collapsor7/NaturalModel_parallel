#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <mpi.h>


using namespace std;

double frand(double a, double b)
{
	return a+(b-a)*(rand()/double(RAND_MAX));
}

int do_walk(int a, int b, int x, double p, double& t)
{
	int step = 0;
	while( x>a && x<b )
	{
		if( frand(0,1)<p )
			x += 1;
		else
			x -= 1;
		t += 1.0;
		step += 1;
	}
	return x;
}

void run_mc(int a, int b, int x, double p, int N, int rank, int size)
{
    double t_local = 0.0;
    double w_local = 0.0;

    // 每个进程处理一部分粒子
    int local_N = N / size;

    for( int i=0; i<local_N; i++ )
    {
        double t = 0.0;
        int out = do_walk(a, b, x, p, t);
        t_local += t;
        if( out == b )
            w_local += 1;
    }

    double t_global, w_global;
    // 汇总所有进程的结果
    MPI_Reduce(&t_local, &t_global, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&w_local, &w_global, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        ofstream f("output.dat");
        f << w_global / N << " " << t_global / N << endl;
        f.close();
    }
}

// void run_mc(int a, int b, int x, double p, int N)
// {
// 	// srand(time(0));
// 	double t = 0.0; 
// 	double w = 0.0; 

// 	for( int i=0; i<N; i++ )
// 	{
// 		int out = do_walk(a, b, x, p, t);
// 		if( out == b )
// 			w += 1;
// 	}

// 	ofstream f("output.dat");
// 	f << w/N << " " << t/N << endl;
// 	f.close();
// }

int main(int argc, char** argv)
{
	MPI_Init(&argc, &argv);
    
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);  // 获取当前进程号
    MPI_Comm_size(MPI_COMM_WORLD, &size);  // 获取进程总数

    int a = atoi(argv[1]);
    int b = atoi(argv[2]);
    int x = atoi(argv[3]);
    double p = atof(argv[4]);
    int N = atoi(argv[5]);

    run_mc(a, b, x, p, N, rank, size);

    MPI_Finalize();
    return 0;
}