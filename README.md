# CSE5441

   This project is to implement cannon's matrix multiplication algorithm using MPI. We expect to have linear scaling for the increasing number of nodes. The script file is run.sh, and the implementation file is cannon.cpp. Figure 1 shows the results of using single node and multiple nodes for running cannon's matrix multiplication algorithm.
   We multiply two matrices to get a new matrix.  The number of processes can be factorized into a square grid of processes. All matrices should be blocked into several blocks. The square grid simplifies things as the local products will be possible. First, figure 2 shows the initial data. Second, we shift data to the preparation. We need to shift row i of A circularly by elements to the left, and shift column j of B circularly by elements up. Figure 3 shows the data after shifting. Third, we start the multiplication. Multiply the local blocks and accumulate the results in the block of results belonging to the process. Only blocks of A and B are communicated. The blocks of results do not move. The multiplication of the local blocks of A and B is also a matrix multiplication. Figure 4 shows the code of multiplication. Last, return the blocks of A and B to the original position.

MPI_Cart_shift(cartComm, 1, coord[0], &left, &right);
MPI_Cart_shift(cartComm, 0, coord[1], &up, &down);
