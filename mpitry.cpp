/**
 * @author RookieHPC
 * @brief Original source code at https://rookiehpc.github.io/mpi/docs/mpi_accumulate/index.html
 **/

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

using namespace std;

/**
 * @brief Illustrate how to accumulate data.
 * @details This application consists of two MPI processes. MPI process 0
 * exposes a window containing an integer initialised to 0. All the other MPI
 * processes add their rank to that value. After the MPI_Accumulate is issued,
 * each MPI process calls on MPI_Win_fence to synchronise. Finally, MPI process
 * 0 prints the total value.
 **/
int main(int argc, char* argv[])
{
    MPI_Init(&argc, &argv);

    // Get my rank
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    // Create the window
    double *y = new double[2];
    if(my_rank) for(int i=0; i<2; i++) y[i]=2.0;
    else  for(int i=0; i<2; i++) y[i]=0.0;
    MPI_Win window;
    MPI_Win_create(y, sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &window);

    MPI_Win_fence(0, window);
    if(my_rank > 0)
    {
        MPI_Accumulate(y, 2, MPI_DOUBLE, 0, 0, 2, MPI_DOUBLE, MPI_SUM, window);
        //printf("[MPI process %d] I accumulate data %d in MPI process 0 window via MPI_Accumulate.\n", my_rank, my_rank);
    }
    MPI_Win_fence(0, window);

    if(my_rank == 0)
    {
        for(int i=0; i<2; i++) cout << y[i] << endl;
    }

    // Destroy the window
    MPI_Win_free(&window);

    MPI_Finalize();

    return EXIT_SUCCESS;
}
