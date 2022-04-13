#ifndef lagrangian_particles_bubble_h
#define lagrangian_particles_bubble_h

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


typedef struct cartesian {

    double x, y;

} cartesian_t;

typedef struct eulerian {

    double U0, L;
    double u, v;
    double vort_z, vort_mag;
    double Du_Dt, Dv_Dt;
    double rho, mu;

} eulerian_t;

typedef struct lagrangian {

    int id;
    double x, y;
    double u, v;
    double rho, D, m;
    struct lagrangian *next;

} lagrangian_t;

typedef struct forces {

    double term_D, term_Vm;
    cartesian_t W, B, Fls;

} forces_t;



void initializing_particles(lagrangian_t *root);
void adding_particles(lagrangian_t *root, double U0, double xi, double xf, double yi, double yf, int npart, double rhop, double Dp);
lagrangian_t *searching_particle(lagrangian_t *root, int id_part);
void tracking_particles(lagrangian_t *lagrangian, double dt, double xi, double xf, double yi, double yf, double bc_x, double bc_y, double U0, double L, double rhof, double mu, cartesian_t g, int npart);
double runge_kutta4(double dt, double vp, double mp, double term_D, double vf, double Fls, double term_Vm, double Dvf_Dt);
double particle_acceleration(double dt, double vp, double mp, double term_D, double vf, double Fls, double term_Vm, double Dvf_Dt);
forces_t forces_calculation(lagrangian_t *ptr_lag, eulerian_t eul, cartesian_t g);
void saving_paraview(lagrangian_t *lagrangian, int npart, int iter, double time, char *output);
void saving_particles(lagrangian_t *lagrangian, int npart, int iter, double t, char *output);

#endif
