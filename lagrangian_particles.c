//  
//  Simulation of lagrangian particles suspended in a fluid flow
//  One-way coupling between fluid and particles is considered
//  The fluid velocity field is pre-set
//  The setup of the present case is based on MAXEY (1987) https://doi.org/10.1063/1.866206
//
//  Created by: Jessica Guarato
//  Last modified on: September, 2018
//


#include "lagrangian_particles.h"

// ==========================================================================================================
// SETUP
// ==========================================================================================================

// Gravity
double gravity_x  = 0.0;            // X-direction
double gravity_y  = 9.81;           // Y-direction

// Flow
double f_U0       = 1.0;            // Flow velocity
double f_L        = 1.0;            // Characteristic length
double f_W        = 0.5;            // Nondimensional settling velocity
double f_St       = 0.1;            // Stokes number
double eul_dens   = 1.2;            // Eulerian density
double eul_visc   = 1.8e-5;         // Eulerian viscosity
int num_part      = 400;            // Number of particles
double lag_dens   = 1000;           // Lagrangian density

// Domain limits
double Xi         = 0.0;           // Initial position of domain in X-direction
double Xf         = 0.0;           // Final position of domain in X-direction => Calculated in main: 2.0*M_PI*f_L;
double Yi         = 0.0;           // Initial position of domain in Y-direction
double Yf         = 0.0;           // Final position of domain in Y-direction => Calculated in main: 2.0*M_PI*f_L;

// Time parameters
long int iter_max = 1000000;       // Maximum number of iterations
double dt         = 1.0e-3;        // Time step
double t_max      = 1.0e4;         // Final time

// Output files
int print_step   = 1000;          // Saving frequency of results
char output[20]  = "output_St01"; // Name of output files


// ==========================================================================================================
// PROGRAM
// ==========================================================================================================
int main() {
    
    lagrangian_t *ptr_lag;
    lagrangian_t *lagrangian = (lagrangian_t *) malloc(sizeof(lagrangian_t));
    cartesian_t gravity;
    long int iter;
    double t, tau, taup, Dp;
    char command1[50], command2[50];
    
     
    gravity.x = gravity_x;
    gravity.y = gravity_y;
    
    Xf = 2.0*M_PI*f_L;
    Yf = 2.0*M_PI*f_L;

    iter = 0;
    t = 0.0;
    
    tau = f_L/f_U0;
    taup = f_St*tau;
    Dp = sqrt(18.0*taup*eul_visc/lag_dens); 

    // Creating output directory
    sprintf(command1, "rm -rf %s",output);
    system(command1);
    sprintf(command2, "mkdir -p %s",output);
    system(command2);

    initializing_particles(lagrangian);
    
    adding_particles(lagrangian, Xi, Xf, Yi, Yf, num_part, lag_dens, Dp);

    saving_paraview(lagrangian, num_part, iter, output); 
    //saving_particles(lagrangian, num_part, iter, t, output); 

    printf("iter_max = %ld\n",iter_max);

    
    while ((iter < iter_max) && (fabs((t_max - t)/t_max) >= 5e-13)) {
    
        iter = iter + 1;
        t = t + dt;
        

        printf("\n#########################\n");
        printf("iter = %ld / %ld\n",iter,(iter_max-iter));
        printf("dt   = %e \n",dt);
        printf("t    = %e \n",t);

        
        tracking_particles(lagrangian, dt, Xi, Xf, Yi, Yf, f_U0, f_L, eul_dens, eul_visc, gravity, num_part);

        if (iter%print_step == 0) {

            saving_paraview(lagrangian, num_part, iter, output); 
            //saving_particles(lagrangian, num_part, iter, t, output); 

        } 
    
    }
    
    return 0;

}


void initializing_particles(lagrangian_t *root) {

    root->next = NULL;

}


void adding_particles(lagrangian_t *root, double xi, double xf, double yi, double yf, 
    int npart, double rhop, double Dp) {
    
    // Defining initial positions of particles
    int id, i, j;
    double dx, dy;
    double posi_x[npart], posi_y[npart];

    id = 0;
    for (i=1; i<=sqrt(npart); i++) {
        for (j=1; j<=sqrt(npart); j++) {
            
            dx = (xf - xi)/sqrt(npart);
            dy = (yf - yi)/sqrt(npart);
            posi_x[id] = (i-0.5)*dx;
            posi_y[id] = (j-0.5)*dy;

            id = id + 1;

        }
    }
    
    for (id=0; id<npart; id++) {

        lagrangian_t *new = (lagrangian_t *) malloc(sizeof(lagrangian_t));

        new->id = id;
        new->x = posi_x[id];
        new->y = posi_y[id];
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
    double U0, double L, double rhof, double mu, cartesian_t g, int npart) {
    
    lagrangian_t *ptr_lag;
    eulerian_t eul;
    forces_t force;
    int id, rk4_iter, n;


    for (id=0; id<npart; id++) {

        // Searching particle
        ptr_lag = searching_particle(lagrangian,id);

            
        // Eulerian properties
        eul.u = U0*cos(ptr_lag->y/L)*sin(ptr_lag->x/L);
        eul.v = -U0*cos(ptr_lag->x/L)*sin(ptr_lag->y/L);
        eul.vort_z = (2.0*U0*sin(ptr_lag->x/L)*sin(ptr_lag->y/L))/L;
        eul.vort_mag = fabs(eul.vort_z);
        eul.rho = rhof;
        eul.mu = mu;


        // Forces calculation
        force = forces_calculation(ptr_lag, eul, g);


        // Time integration
        rk4_iter = 5;
        for (n=1; n<=rk4_iter; n++) {

            // New velocities
            ptr_lag->u = vel_runge_kutta4(dt/rk4_iter, eul.u, ptr_lag->u, ptr_lag->m, force.term_D, force.W.x, force.B.x, force.Fls.x);
            ptr_lag->v = vel_runge_kutta4(dt/rk4_iter, eul.v, ptr_lag->v, ptr_lag->m, force.term_D, force.W.y, force.B.y, force.Fls.y);

        }
        
        // New positions
        ptr_lag->x = ptr_lag->x + ptr_lag->u*dt;
        ptr_lag->y = ptr_lag->y + ptr_lag->v*dt;

        // Periodic boundary condition
        // X-axis
        if (ptr_lag->x < xi) {
            ptr_lag->x = xf + (ptr_lag->x - xi);
        }
        else if (ptr_lag->x >= xf) {
            ptr_lag->x = xi + (ptr_lag->x - xf);
        }
        
        // Y-axis
        if (ptr_lag->y < yi) {
            ptr_lag->y = yf + (ptr_lag->y - yi);
        }
        else if (ptr_lag->y >= yf) {
            ptr_lag->y = yi + (ptr_lag->y - yf);
        }

    }
    
}


double vel_runge_kutta4(double dt, double vf, double vp, double mp, double term_D, double W, double B, double Fls) {

    double k1, k2, k3, k4;
    double new_vp;

    k1 = (term_D*(vf - (vp            )) + W - B + Fls)/mp;
    k2 = (term_D*(vf - (vp + k1*dt/2.0)) + W - B + Fls)/mp;
    k3 = (term_D*(vf - (vp + k2*dt/2.0)) + W - B + Fls)/mp;
    k4 = (term_D*(vf - (vp + k3*dt    )) + W - B + Fls)/mp;

    new_vp = vp + (k1 + 2.0*k2 + 2.0*k3 + k4)*dt/6.0;
    
    return new_vp;

}


forces_t forces_calculation(lagrangian_t *ptr_lag, eulerian_t eul, cartesian_t g) {

    forces_t force;
    double velr, Rep, Cd, term_D;
    double Res, beta, Cls, term_Fls;


    // ========================== Drag force ========================= //

    velr = sqrt((eul.u - ptr_lag->u)*(eul.u - ptr_lag->u) + (eul.v - ptr_lag->v)*(eul.v - ptr_lag->v));
    Rep = (eul.rho*ptr_lag->D*velr)/eul.mu;

    Cd = 0.0;
    if (Rep > 0.0 && Rep <= 0.1) {
        Cd = (24.0/Rep);
    }    
    else if (Rep > 0.1 && Rep <= 1000.0) {
        Cd = 24.0*(1.0 + 0.15*pow(Rep,0.687))/Rep;
    }    
    else if (Rep > 1000.0) {
        Cd = 0.44;
    }

    force.term_D = 0.75*ptr_lag->m*Cd*Rep*eul.mu/(ptr_lag->rho*ptr_lag->D*ptr_lag->D); 

    // =========================== Weigth ============================ //

    force.W.x = ptr_lag->m*g.x;
    force.W.y = ptr_lag->m*g.y;

    // ========================== Buoyancy =========================== //    

    force.B.x = ptr_lag->m*g.x*(eul.rho/ptr_lag->rho);
    force.B.y = ptr_lag->m*g.y*(eul.rho/ptr_lag->rho);

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

    return force;

}


void saving_paraview(lagrangian_t *lagrangian, int npart, int iter, char *output) {
    
    lagrangian_t *ptr_lag;
    char file[50];
    int id;

    
    sprintf(file, "%s/paraview_iter_%010d.vtk",output,iter);
    
    FILE *pvtk;
    pvtk = fopen(file, "w");
    
    fprintf(pvtk,"# vtk DataFile Version 4.1\nPoints\nASCII\nDATASET POLYDATA\n");
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
