#include <iostream>
#include <omp.h>
#include "conio.h"
#include <windows.h>
#include <ctime>
#include <cstdlib>

#define SIZE 4

using namespace std;

double SerialDirectMultiplication(double *MatrixA, double *MatrixB, double *MatrixC) {
	double tStart, tFinish;
	int i, j, k;
	int Size = (int)SIZE;

	tStart = omp_get_wtime();
	for (k = 0; k < Size; k++) {
		for (j = 0; j < Size; j++) {
			//	MatrixC[i * Size + j] = 0; 
			for (i = 0; i < Size; i++) {
				MatrixC[i * Size + j] += MatrixA[i * Size + k] * MatrixB[k * Size + j];
			}
			//printf("%f ", MatrixC[i * Size + j]);
		}
		//printf("\n");
	}
	tFinish = omp_get_wtime();

	return tFinish - tStart;
}

double ParallelDirectMultiplication(double *MatrixA, double *MatrixB, double *MatrixC) {
	double tStart, tFinish;
	int i, j, k;
	int Size = (int)SIZE;

	tStart = omp_get_wtime();
#pragma omp parallel for
	for (i = 0; i < Size; i++) {
		for (j = 0; j < Size; j++) {
			MatrixC[i * Size + j] = 0;
			for (k = 0; k < Size; k++) {
				MatrixC[i * Size + j] += MatrixA[i * Size + k] * MatrixB[k * Size + j];
			}
			//printf("%f ", MatrixC[i* Size + j]);
		}
		//printf("\n");
	}

	tFinish = omp_get_wtime();

	return tFinish - tStart;
}

double SerialBlockMultiplication(double *MatrixA, double *MatrixB, double *MatrixC, int N, int M, int K, int thread) {
	double tStart, tFinish;
	int m, n;

	int Size = (int)SIZE;
	int GridSize = 2; //int (sqrt((double)ThreadNum));
	int BlockSize = Size / GridSize;

	tStart = omp_get_wtime();

	for (n = 0; n < GridSize; n++)
		for (m = 0; m < GridSize; m++)
			for (int iter = 0; iter < GridSize; iter++)
				for (int i = n * BlockSize; i < (n + 1) * BlockSize; i++)
					for (int j = m * BlockSize; j < (m + 1) * BlockSize; j++)
						for (int k = iter * BlockSize; k < (iter + 1) * BlockSize; k++)
							MatrixC[i * Size + j] += MatrixA[i * Size + k] * MatrixB[k * Size + j];

	tFinish = omp_get_wtime();

	return tFinish - tStart;
}

void PrintMatrix(int *Ma, int N, int M) {
	cout << "\n";
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < M; j++) {
			//printf("%10.8f\t", Ma[i*M + j]);
			cout << Ma[i*M + j] << "\t";
		}
		cout << "\n";
	}
}

void ParallelMultiplication22(double *A, double *B, double *C, int N, int M, int K, int thread) {
	int Size = SIZE;
	int GridSize = (int)(sqrt(SIZE));
	cout << "Gridsize = " << GridSize << endl;

	int BlockSize = GridSize;
	if (pow(GridSize, 2) < SIZE) {
		GridSize++;
	}
	int * BlockSizeWight = new int[GridSize * GridSize];
	int * BlockSizeHeight = new int[GridSize * GridSize];

	if (pow(BlockSize, 2) < SIZE) {
		cout << "Gridsize = " << GridSize << endl;

		//BlockSize = GridSize; // (int)(Size / GridSize);
		int BlockSizeH = (int)(Size % BlockSize);
		cout << "Block = " << BlockSize << " blockH = " << BlockSizeH << endl;

		for (int i = 0; i < GridSize; i++) {
			for (int j = 0; j < GridSize - 1; j++) {
				BlockSizeWight[i*GridSize + j] = BlockSize;
				BlockSizeHeight[j*GridSize + i] = BlockSize;
			}
		}

		for (int t = 0; t < GridSize; t++) {
			BlockSizeWight[t*GridSize + GridSize - 1] = BlockSizeH;
			BlockSizeHeight[(GridSize - 1)*GridSize + t] = BlockSizeH;
		}
	}
	else {
		//BlockSize = GridSize;
		
		for (int i = 0; i < GridSize; i++) {
			for (int j = 0; j < GridSize; j++) {
				BlockSizeWight[i*GridSize + j] = BlockSize;
				BlockSizeHeight[j*GridSize + i] = BlockSize;
			}
		}
	}

	PrintMatrix(BlockSizeHeight, GridSize, GridSize);
	PrintMatrix(BlockSizeWight, GridSize, GridSize);

	for (int n = 0; n < GridSize; n++)
		for (int m = 0; m < GridSize; m++)
			for (int iter = 0; iter < GridSize; iter++)
				for (int i = n * BlockSizeWight[n*GridSize + m]; i < (n + 1) * BlockSizeWight[n*GridSize + m]; i++)
					for (int j = m * BlockSizeHeight[n*GridSize + m]; j < (m + 1) * BlockSizeHeight[n*GridSize + m]; j++) {
						//C[i * Size + j] = 0.0;
						for (int k = iter * BlockSizeHeight[n*GridSize + m]; k < (iter + 1) * BlockSizeHeight[n*GridSize + m]; k++)
							//C[i * BlockSizeWight[n*GridSize + m] + j] += A[i * BlockSizeWight[n*GridSize + m] + k] * B[k * BlockSizeHeight[m*GridSize + n] + j];
							C[i * Size + j] += A[i * Size + k] * B[k * Size + j];
					}
						
	delete[] BlockSizeHeight, BlockSizeWight;
}

void ParallelMultiplicationq(double *A, double *B, double *C, int N, int M, int K, int thread) {
	int d = SIZE - N;
	int i, j, k;
//	omp_set_num_threads(thread);
//#pragma omp parallel for
	for (i = 0; i < N; i++) {
		for (j = 0; j < K; j++) {
			C[i * K + j] = 0.0;
			for (k = 0; k < M; k++) {
				C[i * K + j] += A[i * M + k] * B[(k+d) * SIZE + j+d];
			}
		}
	}
}

void ParallelMultiplication(double *A, double *B, double *C, int N, int M, int K, int thread) {
	int i, j, k;
	omp_set_num_threads(thread);
//#pragma omp parallel for
	for (i = 0; i < N; i++) {
		for (j = 0; j < K; j++) {
			C[i * K + j] = 0.0;
			for (k = 0; k < M; k++) {
				C[i * K + j] += A[i * M + k] * B[k * K + j];
			}
		}
	}
}

void Sum(double *A, double *B, int N, int thread) {
	int d = SIZE - N;
	int i, j;
	omp_set_num_threads(thread);
#pragma omp parallel for
	for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
			A[(i+d) * SIZE + (j+d)] += B[i * N + j];
		}
	}
}

void MultiConst(double c, double *Ma, int N, int M, int thread) {
	omp_set_num_threads(thread);
#pragma omp parallel for
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < M; j++) {
			Ma[i * M + j] *= c;
		}
	}
}

void GenerateEye(double *Ma, int N, int thread) {
	omp_set_num_threads(thread); 
#pragma omp parallel for
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			Ma[i*N + j] = 0.0;
		}
		Ma[i*N + i] = 1.0;
	}
}

void ParallelDiff(double *A, double *B, double *C, int N, int thread) {
	int i, j;
	omp_set_num_threads(thread);
#pragma omp parallel for
	for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
			C[i * N + j] = A[i * N + j] - B[i* N + j];
		}
	}
}

void PrintMatrix(double *Ma, int N, int M) {
	cout << "\n";
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < M; j++) {
			//printf("%3.3f\t", Ma[i*M + j]);
			cout << Ma[i*M + j] << "\t";
		}
		cout << "\n";
	}
}

void Transp(double * L, double * LT, int N, int M, int thread) {
	omp_set_num_threads(thread);
#pragma omp parallel for
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < M; j++) {
			LT[j*N + i] = L[i*M + j];
		}
	}
}

double sign(double t) {
	if (abs(t) < 1e-5)
		return 0.0;
	else 
		return t < 0.0 ? -1.0 : 1.0;
}

double scal(double *x, int b) {
	double t = 0.0;
	for (int i = 0; i < b; i++) {
		t += x[i] * x[i];
	}
	return t;
}

double norm2(double *x, int b) {
	return pow(scal(x, b), 0.5);
}

void CopyMatrix(double *O, double *C, int N, int thread) {
	omp_set_num_threads(thread);
#pragma omp parallel for
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++)
			C[i*N + j] = O[i*N + j];
	}
}

void QR(double *A, double *Q, double *R, int N, int thread) {
	double *v = new double[N];
	double *vT = new double[N];

	double *P = new double[N*N];
	double *P2 = new double[N*N];
	double *repP = new double[N*N];
	double b, beta, mu;
	double *B = new double[N*N];

	int K = N;
		
	// function
	for (int j = 0; j < N; j++)
		v[j] = A[j * N];
	cout << "v" << endl; PrintMatrix(v, N, 1);

	mu = norm2(v, N);
	cout << "mu = " << mu << endl; 
	beta = v[0] + sign(v[0]) * mu;
	for (int j = 1; j < N; j++) {
		v[j] /= beta;
	}
	v[0] = 1; //end function
	cout << "v" << endl; PrintMatrix(v, N, 1);

	//---
	b = -2.0 / scal(v, N);
	cout << "b = " << b << endl;
	Transp(v, vT, 1, N, thread);
	cout << "v" << endl; PrintMatrix(v, 1, N);

	ParallelMultiplication(v, vT, B, N, 1, N, thread); // matrix
	cout << "B" << endl; PrintMatrix(B, N, N);
	MultiConst(b, B, N, N, thread);
	cout << "B" << endl; PrintMatrix(B, N, N);
	ParallelMultiplication(B, A, R, N, N, N, thread);
	cout << "R" << endl; PrintMatrix(R, N, N);
	cout << "A" << endl; PrintMatrix(A, N, N);
	Sum(R, A, N, thread);
	cout << "R" << endl; PrintMatrix(R, N, N);

	//here
	//ParallelMultiplicationq(B, Q, repP, N, N, N, thread); // A> R
		cout << "B*Q" << endl; PrintMatrix(repP, K, K);
	
	cout << "K" << K << endl;
	Sum(Q, B, N, thread); //repP
	cout << "Q" << endl; PrintMatrix(Q, N, N);
	//end here

	for (int i = 1; i < N-1; i++) { //  N > N-1 // don't touch =N-1
		K--;

		for (int j = i; j < N; j++)
			v[j - i] = R[j*N + i]; // R v A
		cout << "v" << endl; PrintMatrix(v, K, 1);
		mu = norm2(v, K);
		cout << "mu = " << mu << endl;
		beta = v[0] + sign(v[0]) * mu;
		cout << "beta = " << beta << endl;
		for (int j = 1; j < K; j++) {
			v[j] /= beta;
		}
		v[0] = 1;
		cout << "v" << endl; PrintMatrix(v, K, 1);

		b = -2.0 / scal(v, K);
		cout << "b = " << b << endl;
		Transp(v, vT, 1, K, thread);
		cout << "vT" << endl; PrintMatrix(vT, 1, K);
		ParallelMultiplication(v, vT, B, K, 1, K, thread);
		cout << "B" << endl; PrintMatrix(B, K, K);

		MultiConst(b, B, K, K, thread);
		cout << "b*B" << endl; PrintMatrix(B, K, K);

		cout << "repP" << endl; PrintMatrix(repP, N, N);
		ParallelMultiplicationq(B, R, repP, K, K, K, thread); // A> R // don't touch q
		cout << "repP" << endl; PrintMatrix(repP, K, K);
		cout << "K" << K << endl;
		
		Sum(R, repP, K, thread);
		cout << "RR" << endl; PrintMatrix(R, N, N);

		//---
		cout << "BB" << endl; PrintMatrix(B, K, K);
		//MultiConst(b, B, K, K, thread);
		//cout << "BB" << endl; PrintMatrix(B, K, K);
		ParallelMultiplicationq(B, Q, repP, K, K, K, thread);
		cout << "Q*B" << endl; PrintMatrix(repP, K, K);
		Sum(Q, repP, K, thread);
		cout << "QQ" << endl; PrintMatrix(Q, N, N);		
	}
	//here
	delete[] v, vT, P, P2, repP, B;
}

void fillMatrixA(double *A, int N, bool randomly) {
	/*double A[16] = { 5, 1, 1, 2, 1, 3, 1, 1, 1, 1, 8, 1, 2, 1, 1, 6 };
	double A[25] = {
		5, 1, 1, 3, 1,
		1, 2, 1, 1, 2,
		1, 1, 3, 1, 1,
		3, 1, 1, 8, 1,
		1, 2, 1, 1, 6
	};*/

	double t = 0.0;

	if (randomly)
		srand(clock());
	
	for (int i = 0; i < N; i++) {
		t = 2;
		for (int j = i; j < N; j++) {
			//A[i*N + j] = (double)rand() / RAND_MAX;
			A[i*N + j] = rand() % 10 + 1;
			A[j*N + i] = A[i*N + j];
			t += A[i*N + j];

		}
		A[i*N + i] += t;
	}
	//cout << "A" << endl; PrintMatrix(A, N, N);
}

double eps = 1e-5;
void Giv(double *A, double *Q, int N, int thread) {

	double c = 0.0;
	double s = 0.0;
	double a = 0.0;
	double b = 0.0;
	double t1 = 0.0;
	double t2 = 0.0;
	double q1 = 0.0;
	double q2 = 0.0;

	for (int j = 0; j < N-1; j++) {
		for (int i = N-1; i > j; i--) {
			//---------
			a = A[(i - 1)*N + j];
			b = A[i*N + j];
			//cout << "a = " << a << " b = " << b << endl;
				if ((abs(a) < eps) || (abs(b) < b)) {
					c = 1.0;
					s = 0.0;
				}
				else {
					double t = 0.0;
					if (abs(b) > abs(a)) {
						t = -a / b;
						s = 1 / sqrt(1 + t*t);
						c = s*t;
					}
					else {
						t = -b/ a;
						c = 1 / sqrt(1 + t*t);
						s = c*t;
					}
				}
				//cout << " c = " << c << " s = " << s << endl;
			//---------								
				for (int y = 0; y < N; y++) {
					t1 = A[(i-1)*N + y];
					t2 = A[i*N + y];
					A[(i - 1)*N+ y] = c * t1 - s * t2;
					A[i*N + y] = s * t1 + c * t2;
//					cout << "A: " << endl; PrintMatrix(A, N, N);

					q1 = Q[(i - 1)*N + y];
					q2 = Q[i*N + y];
					Q[(i - 1)*N + y] = c * q1 - s * q2;
					Q[i*N + y] = s * q1 + c * q2;
//					cout << "Q: " << endl; PrintMatrix(Q, N, N);
				}
								
		}
	}
}

int main(int argc, char* argv[]) {
	int N = atoi(argv[1]);
	int thread = atoi(argv[2]);
//	int N = SIZE, thread = 2;
	
	double *A = new double[N*N];
	double *R = new double[N*N];
	double *Q = new double[N*N];
	double *B = new double[N*N];
	double tStart, tFinish;

	fillMatrixA(A, N, 0);
//	cout << "expected A" << endl; PrintMatrix(A, N, N);
	GenerateEye(Q, N, thread);
//	cout << "generated Q" << endl; PrintMatrix(Q, N, N);

	//--- givens
	CopyMatrix(A, R, N, thread);
	
	tStart = omp_get_wtime();
	Giv(R, Q, N, thread);
	tFinish = omp_get_wtime();

//	cout << "R" << endl; PrintMatrix(R, N, N);
	Transp(Q, B, N, N, thread);
//	cout << "QT" << endl; PrintMatrix(B, N, N);
	//--- givens

	//---- house
////	tStart = omp_get_wtime();
////	QR(A, Q, R, N, thread);
////	tFinish = omp_get_wtime();
////
////	cout << "Q" << endl; PrintMatrix(Q, N, N);
////	Transp(Q, B, N, N, thread);
////	cout << "QT" << endl; PrintMatrix(B, N, N);
////	cout << "R" << endl; PrintMatrix(R, N, N);
	//--- house
	ParallelMultiplication(B, R, Q, N, N, N, thread);
	
//	cout << "expected A" << endl; PrintMatrix(A, N, N);
//	cout << "actual A" << endl; PrintMatrix(Q, N, N);

	ParallelDiff(A, Q, B, N, thread);
//	cout << "diff" << endl; PrintMatrix(R, N, N);

	cout << "Time:\n" << tFinish - tStart << endl;
	cout << "Norm of diff:\n" << norm2(B, N*N) << endl;
	
	/*
	double timeSerialDirect = SerialDirectMultiplication(MatrixA, MatrixB, MatrixC);
	printf("Time of serial direct algorithm:   %.12f\n", timeSerialDirect);

	double timeParallelDirect = ParallelDirectMultiplication(MatrixA, MatrixB, MatrixC);
	printf("Time of parallel direct algorithm: %.12f\n", timeParallelDirect);

	double timeSerialBlock = SerialBlockMultiplication(MatrixA, MatrixB, MatrixC);
	printf("Time of serial block algorithm:    %.12f\n", timeSerialBlock);

	double timeParallelBlock = ParallelBlockMultiplication(MatrixA, MatrixB, MatrixC);
	printf("Time of parallel block algorithm:  %.12f\n", timeParallelBlock);
	*/
	/*GenerateEye(Q, N, thread);
	cout << "Q" << endl; PrintMatrix(Q, N, N);
	cout << "A" << endl; PrintMatrix(A, N, N);
	N = 3;
	ParallelMultiplicationq(Q, A, B, N, N, N, thread);
	cout << "B" << endl; PrintMatrix(B, N, N);*/

	delete[] A, R, Q, B;
	return 0;
}