#include <iostream>
#include <math.h>
#include <mpi.h>
#include "vtkOutput.h"
#include "string"
#include <time.h>
#include "Solver.h"

int allocate2DArray(double*** array, int n, int m);
int allocate3DArray(double**** array, int n, int m, int o);

/**
 * Initializes all simulation parameters.
 * 
 **/
void Solver::init() {
    // weights values for directions in D2Q9
    w0 = 4.0 / 9.0;
    w1 = 1.0 / 9.0;
    w2 = 1.0 / 36.0;
    
    // array of weights for each direction - correpsonding to c matrix
    w = new double[9] {w0, w1, w1, w1, w1, w2, w2, w2, w2};

    // grid and time increments
    dt = 1.0;
    dx = 1.0;
    S = dx / dt;

    // coefficients to compute equilibrium distribution
    c1 = 1.0;
    c2 = 3.0 / (pow(S, 2));
    c3 = 9.0 / (2.0 * pow(S, 4));
    c4 = -3.0 / (2.0 * pow(S, 2));

    // physical flow parameters
    nu_f = 0.02;
    rho_0 = 1.0;

    // calculate relaxation time
    tau_f = nu_f * 3.0 / (S * dt) + 0.5;

    // set lattice size and simulation time size
    nt = 100;
    nx = 100;
    nz = 100;

    // create MPI data type for sending to processes
    /*domainSizes_2d = new int[2] {nx, nz};
    subdomainSizes_2d = new int[2] {nx, nz / mpi_size};
    startIndices_2d = new int[2] {0, 0};

    domainSizes_3d_a = new int[3] {na, nz, nx};
    subdomainSizes_3d_a = new int[3] {na, nz / mpi_size, nx};
    startIndices_3d_a = new int[3] {0, 0, 0};

    domainSizes_3d_D = new int[3] {D, nz, nx};
    subdomainSizes_3d_D = new int[3] {D, nz / mpi_size, nx};
    startIndices_3d_D = new int[3] {0, 0, 0};

    MPI_Datatype type_2d, subType_2d, type_3d_a, subType_3d_a, type_3d_D, subType_3d_D;

    int sizes_2d[2] = {nz, nx};
    int subsizes_2d[2] = {nz / mpi_size, nx};
    int starts_2d[2] = {0, 0};
    int sizes_3d_a[3] = {na, nz, nx};
    int subsizes_3d_a[3] = {na, nz / mpi_size, nx};
    int starts_3d_a[3] = {0, 0, 0};
    int sizes_3d_D[3] = {D, nz, nx};
    int subsizes_3d_D[3] = {D, nz / mpi_size, nx};
    int starts_3d_D[3] = {0, 0, 0};
    MPI_Datatype type_2d, subarrtype_2d, type_3d_a, subarrtype_3d_a, type_3d_D, subarrtype_3d_D;
    MPI_Type_create_subarray(2, sizes_2d, subsizes_2d, starts_2d, MPI_ORDER_C, MPI_DOUBLE, &type_2d);
    MPI_Type_create_subarray(3, sizes_3d_a, subsizes_3d_a, starts_3d_a, MPI_ORDER_C, MPI_DOUBLE, &type_3d_a);
    MPI_Type_create_subarray(3, sizes_3d_D, subsizes_3d_D, starts_3d_D, MPI_ORDER_C, MPI_DOUBLE, &type_3d_D);

    MPI_Type_create_resized(type_2d, 0, )*/

    // split only over number of rows
    local_rows = new int[mpi_size];
    for(int i = 0; i < mpi_size; i++) {
        local_rows[i] = nz / mpi_size;
    }

    
    MPI_Type_contiguous(nx, MPI_DOUBLE, &row_type_2d);
    MPI_Type_commit(&row_type_2d);

    displacements_2d = new int[mpi_size];
    
    displacements_2d[0] = 0;
    for(int m = 1; m < mpi_size; m++) {
        displacements_2d[m] = displacements_2d[m-1] + local_rows[m-1];
    }
}

/**
 * Computes domain boundaries for each processor region.
 * 
 **/
void Solver::mpiDivideDomain() {
    // initialize array of region objects - one for each mpi thread
    regions = new MPIRegion[mpi_size];
    if(mpi_rank == 0) {
        for(int i = 0; i < mpi_size; i++) {
            // divide domain into sub-domains vertically
            regions[i].nz = nz / mpi_size;
            // add on any remainder to the last region
            if(i == mpi_size - 1) {
                regions[i].nz += nz % mpi_size;
            }
            // each region goes the full length of the domain
            regions[i].nx = nx;
            // define starting points of each region
            regions[i].nx0 = 0;
            regions[i].nz0 = i * (nz / mpi_size);
        }
        
    }
    // distribute this from the main thread to all other mpi threads
    MPI_Bcast(regions, mpi_size * sizeof(MPIRegion), MPI_BYTE, 0, MPI_COMM_WORLD);
}

/**
 * Sets up local and global storage arrays.
 * 
 **/
void Solver::initialiseStorage() {
    // initialise global dynamic storage arrays - we need to do this on each thread so all nodes have access to global results
    allocate3DArray(&f, na, nz, nx);
    allocate3DArray(&f_stream, na, nz, nx);
    allocate3DArray(&f_eq, na, nz, nx);
    allocate3DArray(&Delta_f, na, nz, nx);
    allocate3DArray(&u, D, nz, nx);
    allocate3DArray(&Pi, D, nz, nx);
    allocate2DArray(&rho, nz, nx);
    allocate2DArray(&u2, nz, nx);
    allocate2DArray(&cu, nz, nx);

    allocate3DArray(&local_f, na, regions[mpi_rank].nz, regions[mpi_rank].nx);
    allocate3DArray(&local_f_stream, na, regions[mpi_rank].nz, regions[mpi_rank].nx);
    allocate3DArray(&local_f_eq, na, regions[mpi_rank].nz, regions[mpi_rank].nx);
    allocate3DArray(&local_Delta_f, na, regions[mpi_rank].nz, regions[mpi_rank].nx);
    allocate3DArray(&local_u, D, regions[mpi_rank].nz, regions[mpi_rank].nx);
    allocate3DArray(&local_Pi, D, regions[mpi_rank].nz, regions[mpi_rank].nx);
    allocate2DArray(&local_rho, regions[mpi_rank].nz, regions[mpi_rank].nx);
    allocate2DArray(&local_u2, regions[mpi_rank].nz, regions[mpi_rank].nx);
    allocate2DArray(&local_cu, regions[mpi_rank].nz, regions[mpi_rank].nx);
}

/**
 * Initialises the flow state. Density is initialised to rho_0 everywhere. Source is also defined here.
 * 
 **/
void Solver::setInitialState() {
    if(mpi_rank == 0) {
        // set initial density
        for(int i = 0; i < nz; i++) {
            for(int j = 0; j < nx; j++) {
                rho[i][j] = 1.0 * rho_0;
            }
        }

        // set initial source and function values
        rho[nz / 2][3 * nx / 4] = 4.0 * rho_0;
        // compute the initial value of f based on density and weights
        for(int i = 0; i < na; i++) {
            for(int j = 0; j < nz; j++) {
                for(int k = 0; k < nx; k++) {
                    f[i][j][k] = rho[j][k] * w[i];
                }
            }
        }
    }

    // distribute to all other threads
    MPI_Bcast(&(rho[0][0]), nx * nz, MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&(f[0][0][0]), na * nx * nz, MPI_FLOAT, 0, MPI_COMM_WORLD);

    MPI_Scatterv(&(rho[0][0]), local_rows, displacements_2d, row_type_2d, &(local_rho[0][0]), local_rows[mpi_rank], row_type_2d, 0, MPI_COMM_WORLD);
    printf("Local rho at 10, 10: %f\n", local_rho[0][0]);
}

/**
 * Main calculation loop - does all iterations of LBM here.
 * 
 **/
void Solver::mainLoop() {
    

    // main simulation loop
    for(int t = 0; t <= nt; t++) {
        // (0) output current state
        if(t % 20 == 0) {
            char f_name[32];
            sprintf(f_name, "output/sim_out_%04d.vti", t);
            write_vtk_frame(f_name, rho, nx, nz);
        }

        // (1) Impose periodic boundary conditions - this is implicit from the streaming term - use of % in construction of connectivity matrix
        /*for(int i = 0; i < na; i++) {
        for(int j = 0; j < nz; j++) {
            f[i][j][0] = f[i][j][nx - 1];
            f[i][j][nx] = f[i][j][1];
        }
        }*/

        // (2) Streaming term
        for(int i = 0; i < na; i++) {
            // first put new values for f into an array f_stream. We can't do this in-place as we want to avoid double-propagating f values.
            for(int j = 0; j < nz; j++) {
                for(int k = 0; k < nx; k++) {
                    // f_stream[i][j][k] = f[i][connectivity[i][j][k][1]][connectivity[i][j][k][0]]; // Used when we use the connectivity matrix
                    f_stream[i][j][k] = f[i][(j - c[i][1] + nz) % nz][(k - c[i][0] + nx) % nx]; // Used when we calculate connectivity analytically
                }
            }
            // once all new f values are computed - copy over to f
            for(int j = 0; j < nz; j++) {
                for(int k = 0; k < nx; k++) {
                    f[i][j][k] = f_stream[i][j][k];
                }
            }
        }
        
        for(int j = 0; j < nz; j++) {
            for(int k = 0; k < nx; k++) {
                // (3) macroscopic properties: rho and u
                // sum rho from function values in all directions at each node
                double rhoSum = 0.0;
                for(int i = 0; i < na; i++) {
                    rhoSum += f[i][j][k];
                }
                rho[j][k] = rhoSum;

                // compute Pi - momentum of node in x-z direction
                double piSum[D] = {0, 0};
                for(int d = 0; d < D; d++) {
                    for(int i = 0; i < na; i++) {
                        piSum[d] += f[i][j][k] * c[i][d];
                    }
                    Pi[d][j][k] = piSum[d];
                    // calculate velocity by momentum / density
                    u[d][j][k] = Pi[d][j][k] / rho[j][k];
                }

                // (4) Equilibrium distribution
                // u^2
                u2[j][k] = u[0][j][k] * u[0][j][k] + u[1][j][k] * u[1][j][k];

                for(int i = 0; i < na; i++) {
                    cu[j][k] = c[i][0]*u[0][j][k] + c[i][1] * u[1][j][k];
                    // compute equilibrium of f at each node - used for collision step
                    f_eq[i][j][k] = rho[j][k] * w[i] * (c1 + c2*cu[j][k] + c3 * pow(cu[j][k], 2) + c4 * u2[j][k]);
                }

                // (5) Collision term
                for(int i = 0; i < na; i++) {
                    // calculate df based on equilibrium distribution and viscocity
                    Delta_f[i][j][k] = (f_eq[i][j][k] - f[i][j][k]) / tau_f;
                    // update f
                    f[i][j][k] += Delta_f[i][j][k];
                }
            }
        }
    }
}

/**
 * Generic method to output data to a vtk file.
 * 
 **/
void Solver::write_vtk_frame(char* fileName, double** data, int nx, int nz) {
    double data_flat[nx * nz];
	for (int i = 0; i < nx * nz; i++) data_flat[i] = data[i / nz][i % nx];
	vtkFileOut f;
	f.Open((char*)fileName);
	f.Init((char*)"rho", 1, 1, 1, nx, nz, 1, 1.0);
    f.WriteField((char*)"rho", (void *) data_flat, sizeof(double), (char*)"Float64", 1);
	f.Finish();
	f.Close();
}

/**
 * Standard function to allocate memory for 2-dimensional array
 * 
 * Parameters:
 *  array: pointer to root of array
 *  n: size of n-dimension
 *  m: size of m-dimension
 **/
int allocate2DArray(float*** array, int n, int m) {
    // allocate the contiguous memory
    float* p = (float*)malloc(n * m * sizeof(float));
    if(!p) {
        return -1;
    }

    // allocate pointers - use a few cases for now as I cbf figuring out how to generalize this
    (*array) = (float**)malloc(n * sizeof(char));
    if(!(*array)) {
        free(p);
        return -1;
    }

    // set up pointers to contiguous memory
    for(int i = 0; i < n; i++) {
        (*array)[i] = &(p[i * m]);
    }
    
    return 0;
}

/**
 * Standard function to allocate memory for 2-dimensional array
 * 
 * Parameters:
 *  array: pointer to root of array
 *  n: size of n-dimension
 *  m: size of m-dimension
 **/
int allocate2DArray(double*** array, int n, int m) {
    // allocate the contiguous memory
    double* p = (double*)malloc(n * m * sizeof(double));
    if(!p) {
        return -1;
    }

    // allocate row pointers
    (*array) = (double**)malloc(n * sizeof(double));
    if(!(*array)) {
        free(p);
        return -1;
    }

    // set up pointers to contiguous memory
    for(int i = 0; i < n; i++) {
        (*array)[i] = &(p[i * m]);
    }
    
    return 0;
}

/**
 * Standard function to allocate memory for 3-dimensional array
 * 
 * Parameters:
 *  array: pointer to root of array
 *  n: size of n-dimension
 *  m: size of m-dimension
 *  o: size of o-dimension
 **/
int allocate3DArray(double**** array, int n, int m, int o) {
    // allocate the contiguous memory
    double* p = (double*)malloc(n * m * o * sizeof(double));
    if(!p) {
        return -1;
    }

    // allocate pointers for first levels (planes)
    (*array) = (double***)malloc(n * sizeof(double*));
    if(!(*array)) {
        free(p);
        return -1;
    }
    double** mPointers = (double**)malloc(n * m * sizeof(double*));

    // allocate pointers for second levels (rows)
    for(int i = 0; i < n; i++) {
        (*array)[i] = &(mPointers[i * m]);
    }

    // set up pointers to contiguous memory
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < m; j++) {
            mPointers[i * m + j] = &(p[i * m * o + j * o]);
        }  
    }
    
    return 0;
}