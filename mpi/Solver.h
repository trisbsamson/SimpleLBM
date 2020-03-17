#include <math.h>
#include "MPIRegion.h"
#include <mpi.h>

class Solver {
    public:

    int mpi_rank;
    int mpi_size;

    MPIRegion* regions;

    // lattice velocity directions for D2Q9 - hardcode for now
    int c[9][2] = {{0, 0}, {1, 0}, {-1, 0}, {0, 1}, {0, -1}, {1, 1}, {-1, -1}, {1, -1}, {-1, 1}};
    // pointers to opposite direction indices in c array
    int ai[9] = {0, 2, 1, 4, 3, 6, 5, 8, 7};

    int na = 9; // number of lattice velocity directions
    int D = 2; // dimension of the simulation
    int Q = na; // dDqQ

    // define weights for direction terms to be used later
    double w0;
    double w1;
    double w2;

    double* w;

    double dt;
    double dx;
    double S;

    double c1;
    double c2;
    double c3;
    double c4;

    double nu_f; // viscocity

    double tau_f; // relaxation time

    int nt; // time-steps
    int nx; // X-axis size
    int nz; // Z-axis size

    int* local_rows;
    int* displacements_2d;
    MPI_Datatype row_type_2d;

    int* domainSizes_2d;
    int* subdomainSizes_2d;
    int* startIndices_2d;

    int* domainSizes_3d_a;
    int* subdomainSizes_3d_a;
    int* startIndices_3d_a;
    
    int* domainSizes_3d_D;
    int* subdomainSizes_3d_D;
    int* startIndices_3d_D;

    // initialize storage arrays
    double*** f;
    double*** f_stream;
    double*** f_eq;
    double*** Delta_f;
    double*** u;
    double*** Pi;
    double** rho;
    double** u2;
    double** cu;

    double*** local_f;
    double*** local_f_stream;
    double*** local_f_eq;
    double*** local_Delta_f;
    double*** local_u;
    double*** local_Pi;
    double** local_rho;
    double** local_u2;
    double** local_cu;

    double rho_0; // initial density
    

    
    void init();
    void mpiDivideDomain();
    void initialiseStorage();
    void setInitialState();
    void mainLoop();
    void write_vtk_frame(char* fileName, double** data, int nx, int nz);
};