#include <math.h>
#include "MPIRegion.h"

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

    double rho_0; // initial density

    void MPIDivideDomain();
    void init();
    void setInitialState();
    void mainLoop();
    void write_vtk_frame(char* fileName, double** data, int nx, int nz);
};