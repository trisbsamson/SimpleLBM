#include <iostream>
#include <math.h>
#include "vtkOutput.h"
#include "string"
#include <time.h>
#include "Solver.h"
#include "MPIRegion.h"

void Solver::init() {

    w0 = 4.0 / 9.0;
    w1 = 1.0 / 9.0;
    w2 = 1.0 / 36.0;

    w = new double[9] {w0, w1, w1, w1, w1, w2, w2, w2, w2};

    dt = 1.0;
    dx = 1.0;
    S = dx / dt;

    c1 = 1.0;
    c2 = 3.0 / (pow(S, 2));
    c3 = 9.0 / (2.0 * pow(S, 4));
    c4 = -3.0 / (2.0 * pow(S, 2));

    nu_f = 0.02;
    rho_0 = 1.0;

    // calculate relaxation time
    tau_f = nu_f * 3.0 / (S * dt) + 0.5;

    nt = 200;
    nx = 201;
    nz = 201;

    // initialize dynamic storage arrays
    f = new double**[na];
    f_stream = new double**[na];
    f_eq = new double**[na];
    Delta_f = new double**[na];
    u = new double**[D];
    Pi = new double**[D];
    rho = new double*[nz];
    u2 = new double*[nz];
    cu = new double*[nz];

    for(int i = 0; i < na; i++) {
        f[i] = new double*[nz];
        f_stream[i] = new double*[nz];
        f_eq[i] = new double*[nz];
        Delta_f[i] = new double*[nz];
        for(int j = 0; j < nz; j++) {
            f[i][j] = new double[nx];
            f_stream[i][j] = new double[nx];
            f_eq[i][j] = new double[nx];
            Delta_f[i][j] = new double[nx];
        }
    }
    for(int d = 0; d < D; d++) {
        u[d] = new double*[nz];
        Pi[d] = new double*[nz];
        for(int j = 0; j < nz; j++) {
            u[d][j] = new double[nx];
            Pi[d][j] = new double[nx];
        }
    }
    for(int j = 0; j < nz; j++) {
        rho[j] = new double[nx];
        u2[j] = new double[nx];
        cu[j] = new double[nx];
    }
}

void Solver::MPIDivideDomain() {
    regions = new MPIRegion[mpi_size];
    if(mpi_rank == 0) {
        for(int i = 0; i < mpi_size; i++) {
            regions[i].nz = nz / mpi_size;
            if(i == mpi_size - 1) {
                regions[i].nz += nz % mpi_size;
            }
            regions[i].nx = nx;
            
            regions[i].nx0 = 0;
            regions[i].nz0 = i * (nz / mpi_size);
        }
    }
}

void Solver::setInitialState() {
    // set initial density
    for(int i = 0; i < nz; i++) {
        for(int j = 0; j < nx; j++) {
            rho[i][j] = 1.0 * rho_0;
        }
    }

    // set initial source and function values
    rho[nz / 2][3 * nx / 4] = 4.0 * rho_0;

    for(int i = 0; i < na; i++) {
        for(int j = 0; j < nz; j++) {
            for(int k = 0; k < nx; k++) {
                f[i][j][k] = rho[j][k] * w[i];
            }
        }
    }
}

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
            for(int j = 0; j < nz; j++) {
                for(int k = 0; k < nx; k++) {
                    // f_stream[i][j][k] = f[i][connectivity[i][j][k][1]][connectivity[i][j][k][0]]; // Used when we use the connectivity matrix
                    f_stream[i][j][k] = f[i][(j - c[i][1] + nz) % nz][(k - c[i][0] + nx) % nx]; // Used when we calculate connectivity analytically
                }
            }

            for(int j = 0; j < nz; j++) {
                for(int k = 0; k < nx; k++) {
                    f[i][j][k] = f_stream[i][j][k];
                }
            }
        }
        
        for(int j = 0; j < nz; j++) {
            for(int k = 0; k < nx; k++) {
                // (3) macroscopic properties: rho and u
                // rho = np.sum(f, axis=0)
                double rhoSum = 0.0;
                for(int i = 0; i < na; i++) {
                    rhoSum += f[i][j][k];
                }
                rho[j][k] = rhoSum;

                // Pi = np.einsum('azx,ad->dzx')
                double piSum[D] = {0, 0};
                for(int d = 0; d < D; d++) {
                    for(int i = 0; i < na; i++) {
                        piSum[d] += f[i][j][k] * c[i][d];
                    }
                    Pi[d][j][k] = piSum[d];
                    u[d][j][k] = Pi[d][j][k] / rho[j][k];
                }

                // (4) Equilibrium distribution
                u2[j][k] = u[0][j][k] * u[0][j][k] + u[1][j][k] * u[1][j][k];

                for(int i = 0; i < na; i++) {
                    cu[j][k] = c[i][0]*u[0][j][k] + c[i][1] * u[1][j][k];
                    f_eq[i][j][k] = rho[j][k] * w[i] * (c1 + c2*cu[j][k] + c3 * pow(cu[j][k], 2) + c4 * u2[j][k]);
                }

                // (5) Collision term
                for(int i = 0; i < na; i++) {
                    Delta_f[i][j][k] = (f_eq[i][j][k] - f[i][j][k]) / tau_f;
                    f[i][j][k] += Delta_f[i][j][k];
                }
            }
        }
    }
}

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