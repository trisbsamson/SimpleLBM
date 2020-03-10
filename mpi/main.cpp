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
    // initialise mpi
    MPI_Init(&argc,&argv);
    // collect info on each thread
    int mpi_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

    int mpi_size;
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    // used for timing
    std::cout.precision(20);
    clock_t tStart = clock();
    // initialise the solver, pass through mpi information
    Solver s;
    s.mpi_rank = mpi_rank;
    s.mpi_size = mpi_size;
    
    // call solver methods
    s.init();
    s.mpiDivideDomain();
    s.initialiseStorage();
    s.setInitialState();
    // to be run on all threads after other methods have been updated to allow mpi
    if(mpi_rank == 0) {
      
      
      s.mainLoop();
    } else {
      printf("processor %d waiting for main loop to finish! \n", mpi_rank);
    }
    
    // report timing information
    std::cout << "Executed in " << (clock() - tStart) << " clock cycles\n";

    MPI_Finalize();
}
