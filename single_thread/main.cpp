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
       
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    std::cout.precision(20);
    clock_t tStart = clock();
    if(world_rank == 0) {
      

      Solver s;
      s.init();
      s.setInitialState();
      s.mainLoop();
    } else {
      printf("processor %d waiting for main loop to finish! \n", world_rank);
    }
    

    std::cout << "Executed in " << (clock() - tStart) << " clock cycles\n";

    MPI_Finalize();
}
