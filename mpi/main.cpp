#include <iostream>
#include <math.h>
#include "vtkOutput.h"
#include "string"
#include "main.h"
#include <time.h>
#include "Solver.h"
#include <mpi.h>

/*
  Main program loop.
*/
int main(int argc, char *argv[]) {
    MPI_Init(&argc,&argv);
       
    int mpi_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

    int mpi_size;
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

    std::cout.precision(20);
    clock_t tStart = clock();

    Solver s;
    s.mpi_rank = mpi_rank;
    s.mpi_size = mpi_size;
    
    s.init();
    s.mpiDivideDomain();
    s.initialiseStorage();

    if(mpi_rank == 0) {
      s.setInitialState();
      s.mainLoop();
    } else {
      printf("processor %d waiting for main loop to finish! \n", mpi_rank);
    }
    

    std::cout << "Executed in " << (clock() - tStart) << " clock cycles\n";

    MPI_Finalize();
}
