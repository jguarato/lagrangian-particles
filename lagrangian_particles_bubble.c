//  
//  Simulation of lagrangian particles suspended in a fluid flow
//  One-way coupling between fluid and particles is considered
//  The fluid velocity field is pre-set
//  The setup of the present case is based on a flow of air bubbles (lagrangian particles) in a two-dimensional channel
//  There are two types boundary conditions to choose
//
//  Created by: Jessica Guarato
//  Last modified on: November, 2018
//

#include "lagrangian_particles_bubble.h"

// ==========================================================================================================
// INPUT
// ==========================================================================================================
// Gravity
double gravity_x   = 0.0;                  // X-direction
double gravity_y   = -9.81;                // Y-direction

// Flow
double f_U0        = 1.0;                  // Flow velocity
double f_L         = 0.25;                 // Characteristic length
double eul_dens    = 998.02;               // Eulerian density
double eul_visc    = 1.002e-3;             // Eulerian viscosity
int num_part       = 100;                  // Number of particles
double lag_dens    = 1.2;                  // Lagrangian density
double lag_diam    = 1.0e-4;               // Lagrangian diameter

// Domain limits
double Xi          = -0.125;               // Initial position of domain in X-direction
double Xf          = 0.125;                // Final position of domain in X-direction
double Yi          = 0.0;                  // Initial position of domain in Y-direction
double Yf          = 1.0;                  // Final position of domain in Y-direction
int BC_X           = 0;                    // Boundary condition in X-direction => (0) Non-slip condition (1) Periodic condition
int BC_Y           = 1;                    // Boundary condition in Y-direction => (0) Non-slip condition (1) Periodic condition 

// Time parameters
long int iter_max  = 4000000;              // Maximum number of iterations
double dt          = 1.0e-4;               // Time step
double t_max       = 50.0;                 // Final time

// Output files
int print_step     = 100;                  // Saving frequency of results
char output[100]   = "output_np100_inlet"; // Name of output files

// ==========================================================================================================
// PROGRAM
// ==========================================================================================================
int main() {
    
    lagrangian_t *ptr_lag;
    lagrangian_t *lagrangian = (lagrangian_t *) malloc(sizeof(lagrangian_t));
    cartesian_t gravity;
    long int iter;
    double g, t, tau, taup;
    char command1[50], command2[50];
    
    
    gravity.x = gravity_x;
    gravity.y = gravity_y;

    iter = 0;
    t = 0.0;


    // Creating output directory
    sprintf(command1, "rm -rf %s",output);
    system(command1);
    sprintf(command2, "mkdir -p %s",output);
    system(command2);

    initializing_particles(lagrangian);
    
    adding_particles(lagrangian, f_U0, Xi, Xf, Yi, Yf, num_part, lag_dens, lag_diam);

    saving_paraview(lagrangian, num_part, iter, t, output); 
    //saving_particles(lagrangian, num_part, iter, t, output); 

    
    while ((iter < iter_max) && (t_max - t)/dt >= 0.0) {
    
        iter = iter + 1;
        t = t + dt;
        
        printf("\n#########################\n");
        printf("iter = %ld \n",iter);
        printf("dt   = %e  \n",dt);
        printf("t    = %e  \n",t);
        
        tracking_particles(lagrangian, dt, Xi, Xf, Yi, Yf, BC_X, BC_Y, f_U0, f_L, eul_dens, eul_visc, gravity, num_part);

        if (iter%print_step == 0 || (t_max - t)/dt == 0.0) {

            saving_paraview(lagrangian, num_part, iter, t, output); 
            //saving_particles(lagrangian, num_part, iter, t, output); 

        } 
    
    }
    
    return 0;

}


void initializing_particles(lagrangian_t *root) {

    root->next = NULL;

}


void adding_particles(lagrangian_t *root, double U0, double xi, double xf, double yi, double yf, 
                      int npart, double rhop, double Dp) {
    
    int id;
    double randx;

    for (id=0; id<npart; id++) {

        lagrangian_t *new = (lagrangian_t *) malloc(sizeof(lagrangian_t));

        randx = (double)rand()/(double)RAND_MAX;
        xi = xi - copysign(0.0001,xi);
        xf = xf - copysign(0.0001,xf);
        
        new->id = id;
        new->x = xi + randx*(xf - xi);
        new->y = 0.0;
        new->u = 0.0;
        new->v = 0.0;
        new->rho = rhop;
        new->D = Dp;
        new->m = rhop*((1.0/6.0)*M_PI*pow(Dp,3.0));
        
        new->next = NULL;
        
        if (root->next == NULL) {
            root->next = new;
        }
        else {

            lagrangian_t *tmp = root->next;
            
            while(tmp->next != NULL) {
                tmp = tmp->next;
            }
            
            tmp->next = new;

        }

    }

}


lagrangian_t *searching_particle(lagrangian_t *root, int id_part) {
    
    lagrangian_t *tmp = root->next;
    
    while(tmp != NULL) {
        
        if(tmp->id == id_part) break;
        
        tmp = tmp->next;

    }
    
    return tmp;

}


void tracking_particles(lagrangian_t *lagrangian, double dt, double xi, double xf, double yi, double yf, 
    double bc_x, double bc_y, double U0, double L, double rhof, double mu, cartesian_t g, int npart) {
    
    lagrangian_t *ptr_lag;
    eulerian_t eul;
    forces_t force;
    int id, rk4_iter, n;


    for (id=0; id<npart; id++) {

        // Searching particle
        ptr_lag = searching_particle(lagrangian,id);

            
        // Eulerian properties - Poiseuille flow
        eul.rho = rhof;
        eul.mu = mu;
        eul.u = 0.0;
        eul.v = U0*(1.0 - (ptr_lag->x*ptr_lag->x)/(0.25*L*L));
        eul.vort_z = -2.0*U0*ptr_lag->x/(0.25*L*L);
        eul.vort_mag = fabs(eul.vort_z);
        eul.Du_Dt = 0.0;
        eul.Dv_Dt = 0.0;


        // Forces calculation
        force = forces_calculation(ptr_lag, eul, g);

        // Time integration
        rk4_iter = 10;
        for (n=1; n<=rk4_iter; n++) {

            // New velocities
            ptr_lag->u = runge_kutta4(dt/rk4_iter, ptr_lag->u, ptr_lag->m, force.term_D, eul.u, force.Fls.x, force.term_Vm, eul.Du_Dt);
            ptr_lag->v = runge_kutta4(dt/rk4_iter, ptr_lag->v, ptr_lag->m, force.term_D, eul.v, force.Fls.y, force.term_Vm, eul.Dv_Dt);

        }
        
        // New positions
        ptr_lag->x = ptr_lag->x + ptr_lag->u*dt;
        ptr_lag->y = ptr_lag->y + ptr_lag->v*dt;

        // Periodic boundary condition
        // X-axis
        if (bc_x == 0) {

            if (ptr_lag->x < xi) {
                ptr_lag->x = xi - (ptr_lag->x - xi);
                ptr_lag->u = -ptr_lag->u;
            }
            else if (ptr_lag->x >= xf) {
                ptr_lag->x = xf - (ptr_lag->x - xf);
                ptr_lag->u = -ptr_lag->u;
            } 

        }
        else if (bc_x == 1) {

            if (ptr_lag->x < xi) {
                ptr_lag->x = xf + (ptr_lag->x - xi);
            }
            else if (ptr_lag->x >= xf) {
                ptr_lag->x = xi + (ptr_lag->x - xf);
            }

        }

        // Y-axis
        if (bc_y == 0) {

            if (ptr_lag->y < yi) {
                ptr_lag->y = yi - (ptr_lag->y - yi);
                ptr_lag->v = -ptr_lag->v;
            }
            else if (ptr_lag->y >= yf) {
                ptr_lag->y = yf - (ptr_lag->y - yf);
                ptr_lag->v = -ptr_lag->v;
            } 

        }
        else if (bc_y == 1) {

            if (ptr_lag->y < yi) {
                ptr_lag->y = yf + (ptr_lag->y - yi);
            }
            else if (ptr_lag->y >= yf) {
                ptr_lag->y = yi + (ptr_lag->y - yf);
            }

        }

    }
    
}


double runge_kutta4(double dt, double vp, double mp, double term_D, double vf, double Fls, 
    double term_Vm, double Dvf_Dt) {

    double k1, k2, k3, k4;
    double new_vp;

    k1 = particle_acceleration(dt, (vp),              mp, term_D, vf, Fls, term_Vm, Dvf_Dt);
    k2 = particle_acceleration(dt, (vp + k1*dt/2.0),  mp, term_D, vf, Fls, term_Vm, Dvf_Dt);
    k3 = particle_acceleration(dt, (vp + k2*dt/2.0),  mp, term_D, vf, Fls, term_Vm, Dvf_Dt);
    k4 = particle_acceleration(dt, (vp + k3*dt),      mp, term_D, vf, Fls, term_Vm, Dvf_Dt);

    new_vp = vp + (k1 + 2.0*k2 + 2.0*k3 + k4)*dt/6.0;
    
    return new_vp;

}


double particle_acceleration(double dt, double vp, double mp, double term_D, double vf, double Fls, 
    double term_Vm, double Dvf_Dt) {

    return (term_D*(vf - vp) + Fls + term_Vm*Dvf_Dt)/(mp + term_Vm);

}


forces_t forces_calculation(lagrangian_t *ptr_lag, eulerian_t eul, cartesian_t g) {

    forces_t force;
    double velr, Rep, Cd, term_D;
    double Res, beta, Cls, term_Fls;


    // ========================== Drag force ========================= //

    velr = sqrt((eul.u - ptr_lag->u)*(eul.u - ptr_lag->u) + (eul.v - ptr_lag->v)*(eul.v - ptr_lag->v));
    Rep = (eul.rho*ptr_lag->D*velr)/eul.mu;

    Cd = 0.0;

    // Solid particles
    // if (Rep > 0.0 && Rep <= 0.1) {
    //     Cd = (24.0/Rep); // Stokes regime (Stokes, 1851)
    // }    
    // else if (Rep > 0.1 && Rep <= 1000.0) {
    //     Cd = 24.0/Rep*(1.0 + 0.15*pow(Rep,0.687)); // Standard drag correlation (Schiller and Naumann, 1933)
    // }    
    // else if (Rep > 1000.0) {
    //     Cd = 0.44; // Newton regime
    // }

    // Bubbles (Lain et al., 2002)
    if (Rep > 0.0 && Rep <= 1.5) {
        Cd = 16.0/Rep;
    }    
    else if (Rep > 1.5 && Rep <= 80.0) {
        Cd = 14.9/pow(Rep,0.78);
    }
    else if (Rep > 80.0 && Rep <= 1530.0) {
        Cd = 48.0/Rep*(1 - 2.21/pow(Rep,0.5)) + 1.86e-15*pow(Rep,4.756);
    }
    else if (Rep > 1530.0) {
        Cd = 2.61;
    }

    force.term_D = 0.75*ptr_lag->m*Cd*Rep*eul.mu/(ptr_lag->rho*ptr_lag->D*ptr_lag->D); 

    // ================== Shear Lift Force - Saffman ================= // 

    Res = eul.rho*(ptr_lag->D*ptr_lag->D)*eul.vort_mag/eul.mu;
    
    beta = 0.0;
    if (Rep > 0.0) {
       beta = 0.5*Res/Rep;
    }

    Cls = 0.0; term_Fls = 0.0;
    if (beta > 0.005 && beta < 0.4) {
 
        if (Rep <= 40.0) {
            Cls = (1.0 - 0.3314*sqrt(beta))*exp(-Rep/10.0) + 0.3314*sqrt(beta);
        }
        else if (Rep > 40.0) {
            Cls = 0.0524*sqrt(beta*Rep);
        }

        term_Fls = 1.615*ptr_lag->D*eul.mu*sqrt(Res)/eul.vort_mag*Cls;

    }

    force.Fls.x = term_Fls*((eul.v - ptr_lag->v)*eul.vort_z);
    force.Fls.y = term_Fls*(-(eul.u - ptr_lag->u)*eul.vort_z);

    // ================== Virtual Mass ================= //

    force.term_Vm = (1.0/2.0)*ptr_lag->m*(eul.rho/ptr_lag->rho);

    return force;

}


void saving_paraview(lagrangian_t *lagrangian, int npart, int iter, double time, char *output) {
    
    lagrangian_t *ptr_lag;
    char file[50];
    int id;

    
    sprintf(file, "%s/paraview_iter_%010d.vtk",output,iter);
    
    FILE *pvtk;
    pvtk = fopen(file, "w");
    
    fprintf(pvtk,"# vtk DataFile Version 4.1\nPoints\nASCII\nDATASET POLYDATA\n");
    fprintf(pvtk,"FIELD FieldData 1\n");

    fprintf(pvtk,"Time 1 1 double\n%f\n",time);

    fprintf(pvtk,"POINTS	%d	double\n",npart);
    
    for (id=0; id<npart; id++) {
        ptr_lag = searching_particle(lagrangian,id);
        fprintf(pvtk,"%6.16lf	%6.16lf	%6.16lf\n",ptr_lag->x, ptr_lag->y, 0.0);
    }
    
    fprintf(pvtk,"POINT_DATA	%d\n", npart);
    
    fprintf(pvtk,"SCALARS	id double\n");
    fprintf(pvtk,"LOOKUP_TABLE default\n");
    
    for (id=0; id<npart; id++) {
        ptr_lag = searching_particle(lagrangian,id);
        fprintf(pvtk,"%d\n",ptr_lag->id);
    }
    
    fprintf(pvtk,"SCALARS	diameter double\n");
    fprintf(pvtk,"LOOKUP_TABLE default\n");
    
    for (id=0; id<npart; id++) {
        ptr_lag = searching_particle(lagrangian,id);
        fprintf(pvtk,"%6.16f\n",ptr_lag->D);
    }
    
    fprintf(pvtk,"SCALARS	u double\n");
    fprintf(pvtk,"LOOKUP_TABLE default\n");
    
    for (id=0; id<npart; id++) {
        ptr_lag = searching_particle(lagrangian,id);
        fprintf(pvtk,"%6.16f\n",ptr_lag->u);
    }
    
    fprintf(pvtk,"SCALARS	v double\n");
    fprintf(pvtk,"LOOKUP_TABLE default\n");
    
    for (id=0; id<npart; id++) {
        ptr_lag = searching_particle(lagrangian,id);
        fprintf(pvtk,"%6.16f\n",ptr_lag->v);
    }
    
    fclose(pvtk);

}



void saving_particles(lagrangian_t *lagrangian, int npart, int iter, double t, char *output) {
    
    lagrangian_t *ptr_lag;
    char file[50];
    int id;

    
    sprintf(file, "%s/iter_%010d.dat",output,iter);
    
    FILE *dat;
    dat = fopen(file, "w");
    
    fprintf(dat,"%6.16f\n",t); // Time
    
    for (id=0; id<npart; id++) {
        ptr_lag = searching_particle(lagrangian,id);
        fprintf(dat,"%6.16f    %6.16f\n",ptr_lag->x, ptr_lag->y);
    }
    
    fclose(dat);

}
