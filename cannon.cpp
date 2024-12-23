#include<stdio.h>
#include<stdlib.h>
#include<mpi.h>
#include<math.h>
#include<time.h>

int allocMatrix(int*** mat, int rows, int cols) {
	int* p = (int*)malloc(sizeof(int*) * rows * cols);
	if (!p) {
		return -1;
	}
	*mat = (int**)malloc(rows * sizeof(int*));
	if (!mat) {
		free(p);
		return -1;
	}
	for (int i = 0; i < rows; i++) {
		(*mat)[i] = &(p[i * cols]);
	}
	return 0;
}

int freeMatrix(int ***mat) {
	free(&((*mat)[0][0]));
	free(*mat);
	return 0;
}

void matrixMultiply(int **a, int **b, int rows, int cols, int ***c) {
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {
			int val = 0;
			for (int k = 0; k < rows; k++) {
				val += a[i][k] * b[k][j];
 			}
			(*c)[i][j] = val;
		}
	}
}

int main(int argc, char* argv[]) {
        clock_t t1, t2;
	MPI_Comm cartComm;
	int dim[2], period[2], reorder;
	int coord[2], id;
	int **A = NULL, **B = NULL, **C = NULL, **D = NULL;
	int **localA = NULL, **localB = NULL, **localC = NULL;
	int **localARec = NULL, **localBRec = NULL;
	int rows = 1600;
	int columns = 1600;
	int worldSize;
	int procDim;
	int blockDim;
	int left, right, up, down;
	int bCastData[4];
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &worldSize);
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (rank == 0) {
                t1 = clock();
		int n;
		char ch;
		double sqroot = sqrt(worldSize);
		int intRoot = (int)sqroot;
		procDim = intRoot;
		blockDim = columns / intRoot;
		if (allocMatrix(&A, rows, columns) != 0) {
			printf("[ERROR] Matrix alloc for A failed!\n");
			MPI_Abort(MPI_COMM_WORLD, 4);
		}
		if (allocMatrix(&B, rows, columns) != 0) {
			printf("[ERROR] Matrix alloc for B failed!\n");
			MPI_Abort(MPI_COMM_WORLD, 5);
		}
                srand(time(NULL));
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < columns; j++) {
				A[i][j] = rand() % 10;
			}
		}
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < columns; j++) {
				B[i][j] = rand() % 10;
			}
		}
                if (allocMatrix(&C, rows, columns) != 0) {
			printf("[ERROR] Matrix alloc for C failed!\n");
			MPI_Abort(MPI_COMM_WORLD, 6);
		}
                if (allocMatrix(&D, rows, columns) != 0) {
                        printf("[ERROR] Matrix alloc for C failed!\n");
                        MPI_Abort(MPI_COMM_WORLD, 6);
                }
		bCastData[0] = procDim;
		bCastData[1] = blockDim;
		bCastData[2] = rows;
		bCastData[3] = columns;
	}
	MPI_Bcast(&bCastData, 4, MPI_INT, 0, MPI_COMM_WORLD);
	procDim = bCastData[0];
	blockDim = bCastData[1];
	rows = bCastData[2];
	columns = bCastData[3];
	dim[0] = procDim; dim[1] = procDim;
	period[0] = 1; period[1] = 1;
	reorder = 1;
	MPI_Cart_create(MPI_COMM_WORLD, 2, dim, period, reorder, &cartComm);
	allocMatrix(&localA, blockDim, blockDim);
	allocMatrix(&localB, blockDim, blockDim);
	int globalSize[2] = { rows, columns };
	int localSize[2] = { blockDim, blockDim };
	int starts[2] = { 0,0 };
	MPI_Datatype type, subarrtype;
	MPI_Type_create_subarray(2, globalSize, localSize, starts, MPI_ORDER_C, MPI_INT, &type);
	MPI_Type_create_resized(type, 0, blockDim * sizeof(int), &subarrtype);
	MPI_Type_commit(&subarrtype);
	int *globalptrA = NULL;
	int *globalptrB = NULL;
	int *globalptrC = NULL;
	if (rank == 0) {
		globalptrA = &(A[0][0]);
		globalptrB = &(B[0][0]);
		globalptrC = &(C[0][0]);
	}
	int* sendCounts = (int*)malloc(sizeof(int) * worldSize);
	int* displacements = (int*)malloc(sizeof(int) * worldSize);
	if (rank == 0) {
		for (int i = 0; i < worldSize; i++) {
			sendCounts[i] = 1;
		}
		int disp = 0;
		for (int i = 0; i < procDim; i++) {
			for (int j = 0; j < procDim; j++) {
				displacements[i * procDim + j] = disp;
				disp += 1;
			}
			disp += (blockDim - 1)* procDim;
		}
	}
	MPI_Scatterv(globalptrA, sendCounts, displacements, subarrtype, &(localA[0][0]),
		rows * columns / (worldSize), MPI_INT,
		0, MPI_COMM_WORLD);
	MPI_Scatterv(globalptrB, sendCounts, displacements, subarrtype, &(localB[0][0]),
		rows * columns / (worldSize), MPI_INT,
		0, MPI_COMM_WORLD);
	if (allocMatrix(&localC, blockDim, blockDim) != 0) {
		printf("[ERROR] Matrix alloc for localC in rank %d failed!\n", rank);
		MPI_Abort(MPI_COMM_WORLD, 7);
	}
	MPI_Cart_coords(cartComm, rank, 2, coord);
	MPI_Cart_shift(cartComm, 1, coord[0], &left, &right);
	MPI_Sendrecv_replace(&(localA[0][0]), blockDim * blockDim, MPI_INT, left, 1, right, 1, cartComm, MPI_STATUS_IGNORE);
	MPI_Cart_shift(cartComm, 0, coord[1], &up, &down);
	MPI_Sendrecv_replace(&(localB[0][0]), blockDim * blockDim, MPI_INT, up, 1, down, 1, cartComm, MPI_STATUS_IGNORE);
	for (int i = 0; i < blockDim; i++) {
		for (int j = 0; j < blockDim; j++) {
			localC[i][j] = 0;
		}
	}
	int** multiplyRes = NULL;
	if (allocMatrix(&multiplyRes, blockDim, blockDim) != 0) {
		printf("[ERROR] Matrix alloc for multiplyRes in rank %d failed!\n", rank);
		MPI_Abort(MPI_COMM_WORLD, 8);
	}
	for (int k = 0; k < procDim; k++) {
		matrixMultiply(localA, localB, blockDim, blockDim, &multiplyRes);
		for (int i = 0; i < blockDim; i++) {
			for (int j = 0; j < blockDim; j++) {
				localC[i][j] += multiplyRes[i][j];
			}
		}
		MPI_Cart_shift(cartComm, 1, 1, &left, &right);
		MPI_Cart_shift(cartComm, 0, 1, &up, &down);
		MPI_Sendrecv_replace(&(localA[0][0]), blockDim * blockDim, MPI_INT, left, 1, right, 1, cartComm, MPI_STATUS_IGNORE);
		MPI_Sendrecv_replace(&(localB[0][0]), blockDim * blockDim, MPI_INT, up, 1, down, 1, cartComm, MPI_STATUS_IGNORE);
	}
	MPI_Gatherv(&(localC[0][0]), rows * columns / worldSize, MPI_INT,
		globalptrC, sendCounts, displacements, subarrtype,
		0, MPI_COMM_WORLD);
	freeMatrix(&localC);
	freeMatrix(&multiplyRes);
	if (rank == 0) {
                t2 = clock();
                double time_taken = ((double)(t2 - t1)) / CLOCKS_PER_SEC;
                printf("Row: %d, Column: %d\n", rows, columns);
                printf("Process: %d\n", worldSize);
                printf("Time: %f sec\n", time_taken);
                
                for (int a = 0; a < columns; a++)
                    for (int b = 0; b < columns; b++)
                            D[a][b] = 0;
                for (int a = 0; a < rows; a++)
	            for (int b = 0; b < rows; b++)
                        for (int c = 0; c < rows; c++)
	                    D[a][c] += A[a][b] * B[b][c];
                double diff = 0.0;
                #define abs(x) ((x) >= 0? (x) : -(x))
                for (int a = 0; a < columns; a++)
                    for (int b = 0; b < columns; b++)
	                diff += abs(D[a][b]-C[a][b]);
                if (diff > 0.000001)
                    printf ("CHECK FAIL (%lf)\n", diff);
                else
                    printf ("CHECK PASS\n");

	}
	MPI_Finalize();   
	return 0;
}
