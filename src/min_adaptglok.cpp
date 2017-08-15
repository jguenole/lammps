/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include <math.h>
#include "min_adaptglok.h"
#include "universe.h"
#include "atom.h"
#include "force.h"
#include "update.h"
#include "output.h"
#include "timer.h"
#include "error.h"
#include "variable.h"
#include "modify.h"
#include "compute.h"
#include "domain.h"

using namespace LAMMPS_NS;

// EPS_ENERGY = minimum normalization for energy tolerance

#define EPS_ENERGY 1.0e-8

/* ---------------------------------------------------------------------- */

MinAdaptGlok::MinAdaptGlok(LAMMPS *lmp) : Min(lmp) {}

/* ---------------------------------------------------------------------- */

void MinAdaptGlok::init()
{
  Min::init();

  dt = dtinit = update->dt;
  dtmax = tmax * dt;
  dtmin = tmin * dt;
  alpha = alpha0;
  last_negative = ntimestep_fire = update->ntimestep;
  
  if (relaxbox_flag){
    int icompute = modify->find_compute("thermo_temp");
    temperature = modify->compute[icompute];
    icompute = modify->find_compute("thermo_press");
    pressure = modify->compute[icompute];

    /*
    no_pressure_flag = 0;
    xprdinit = domain->xprd;
    yprdinit = domain->yprd;
    zprdinit = domain->zprd;
    vol0 = xprdinit * yprdinit * zprdinit;
    */
  }
}

/* ---------------------------------------------------------------------- */

void MinAdaptGlok::setup_style()
{
  double **v = atom->v;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    v[i][0] = v[i][1] = v[i][2] = 0.0;
}

/* ----------------------------------------------------------------------
   set current vector lengths and pointers
   called after atoms have migrated
------------------------------------------------------------------------- */

void MinAdaptGlok::reset_vectors()
{
  // atomic dof

  nvec = 3 * atom->nlocal;
  if (nvec) xvec = atom->x[0];
  if (nvec) fvec = atom->f[0];
}

/* ----------------------------------------------------------------------
   save current box state for converting atoms to lamda coords
------------------------------------------------------------------------- */

void MinAdaptGlok::save_box_state()
{
  boxlo[0] = domain->boxlo[0];
  boxlo[1] = domain->boxlo[1];
  boxlo[2] = domain->boxlo[2];

  for (int i = 0; i < 6; i++)
    h_inv[i] = domain->h_inv[i];
}

/* ----------------------------------------------------------------------
   deform the simulation box and remap the particles
------------------------------------------------------------------------- */

void MinAdaptGlok::relax_box()
{
  int i,n;
      
  // rescale simulation box from linesearch starting point
  // scale atom coords for all atoms or only for fix group atoms

  //fprintf(screen,"boxhi %g \n", domain->boxhi[0]);

  double **x = atom->x;
  double s0,E,epsilon;
  int *mask = atom->mask;
  n = atom->nlocal + atom->nghost;
  save_box_state();

  // convert pertinent atoms and rigid bodies to lamda coords

  domain->x2lamda(n);

  // ensure the virial is tallied
  // update the corresponding flag
  pressure->addstep(update->ntimestep);
  update->vflag_global = update->ntimestep;
  
  // compute pressure and change simulation box

  //fprintf(screen,"boxhi %g \n", domain->boxhi[0]);

  pressure->compute_scalar();``
  pressure->compute_vector();
  pressure_current = pressure->vector;
  pressure_current_s = pressure->scalar;
  //fprintf(screen,"pressure_previous_s %g\n", pressure_previous_s);
  if (no_pressure_flag){
    //fprintf(screen,"-boxhip %g \n", boxhi_previous[0]);
    //fprintf(screen,"-boxhi %g \n", domain->boxhi[0]);
    E = (pressure_previous_s - pressure_current_s) / (boxhi_previous[0] - domain->boxhi[0]);
    s0 = pressure_current_s - E * domain->boxhi[0];
    //fprintf(screen,"(pressure_previous_s - pressure_current_s) / (boxhi_previous[0] - domain->boxhi[0]) %g %g %g %g\n", pressure_previous_s,pressure_current_s,boxhi_previous[0],domain->boxhi[0]);
    box0[0] = - s0 / E;
    //fprintf(screen,"box0 %g \n",box0[0]);
    epsilon = (box0[0] - domain->boxhi[0]) / domain->boxhi[0];
    epsilon = pressure_current_s / relaxbox_modulus;
    //epsilon = 0.0001;
    //epsilon = - pressure_current_s / E;
  }else{
    epsilon = 0.0001;
    no_pressure_flag = 1;
    //fprintf(screen,"TEST %g \n", pressure_current_s);
  }
  //fprintf(screen,"E s0 b0 epsilon %g %g %g %g \n", E,s0,box0[0],epsilon);

  //pressure_previous_s = pressure_current_s;
  //boxhi_previous = domain->boxhi;

  for (int i = 0; i < 3; i++) {
    pressure_previous_s = pressure_current_s;
    boxhi_previous[i] = domain->boxhi[i];
    //fprintf(screen,"--boxhip %g \n", boxhi_previous[i]);
    //fprintf(screen,"--boxhi %g \n", domain->boxhi[i]);
    //epsilon = boxhi_previous
    //if (relaxbox_flag == 2) epsilon = pressure_current[i] / relaxbox_modulus;
    domain->boxhi[i] += domain->boxhi[i] * epsilon * relaxbox_rate;
    //fprintf(screen,"boxhi %g %g \n", domain->boxhi[i],boxhi_previous[i]);
  }

  //fprintf(screen,"boxhi %g \n", domain->boxhi[0]);

  /*
  epsilon = pressure->scalar / relaxbox_modulus;
  for (int i = 0; i < 3; i++) {
    if (relaxbox_flag == 2) epsilon = pressure_current[i] / relaxbox_modulus;
    domain->boxhi[i] += domain->boxhi[i] * epsilon * relaxbox_rate;;
  }
  */

  // reset global and local box to new size/shape

  domain->set_global_box();
  domain->set_local_box();

  // convert atoms and back to box coords

  domain->lamda2x(n);
  save_box_state();

  /*
  pv2e = 1.0 / force->nktv2p;
  scale = domain->xprd/xprdinit;
  nextra_global = 1;  
  fprintf(screen,"nextra_global %i \n", nextra_global);
  fprintf(screen,"value %g \n", pv2e * (pressure_current_s)*3.0*scale*scale*vol0);
  if (nextra_global)
    for (int i = 0; i < 3; i++) {
      fextra[i] = pv2e * (pressure_current_s)*3.0*scale*scale*vol0;
    }
  */
}

/* ---------------------------------------------------------------------- */

int MinAdaptGlok::iterate(int maxiter)
{
  bigint ntimestep;
  double vmax,vdotf,vdotfall,vdotv,vdotvall,fdotf,fdotfall;
  double scale1,scale2;
  double dtvone,dtv,dtf,dtfm;
  int flag,flagall;

  alpha_final = 0.0;

  double **f = atom->f;
  double **v = atom->v;

  // Leap Frog integration initialization

  if (integrator == 2) {

    double *rmass = atom->rmass;
    double *mass = atom->mass;
    int *type = atom->type;
    int nlocal = atom->nlocal;

    energy_force(0);
    neval++;

    dtf = 0.5 * dt * force->ftm2v;
    if (rmass) {
      for (int i = 0; i < nlocal; i++) {
        dtfm = dtf / rmass[i];
        v[i][0] = dtfm * f[i][0];
        v[i][1] = dtfm * f[i][1];
        v[i][2] = dtfm * f[i][2];
      }
    } else {
      for (int i = 0; i < nlocal; i++) {
        dtfm = dtf / mass[type[i]];
        v[i][0] = dtfm * f[i][0];
        v[i][1] = dtfm * f[i][1];
        v[i][2] = dtfm * f[i][2];
      }
    }

  }

  for (int iter = 0; iter < maxiter; iter++) {

    if (timer->check_timeout(niter))
      return TIMEOUT;

    ntimestep = ++update->ntimestep;
    niter++;

    // Relax the simulation box

    //if (relaxbox_flag) relax_box();
   
   // vdotfall = v dot f

    // Euler || Leap Frog integration

    if (integrator == 0 || integrator == 2) {
      double **v = atom->v;
      double **f = atom->f;
    }
    int nlocal = atom->nlocal;
    double **x = atom->x;

    vdotf = 0.0;
    for (int i = 0; i < nlocal; i++)
      vdotf += v[i][0]*f[i][0] + v[i][1]*f[i][1] + v[i][2]*f[i][2];
    MPI_Allreduce(&vdotf,&vdotfall,1,MPI_DOUBLE,MPI_SUM,world);

    // sum vdotf over replicas, if necessary
    // this communicator would be invalid for multiprocess replicas

    if (update->multireplica == 1) {
      vdotf = vdotfall;
      MPI_Allreduce(&vdotf,&vdotfall,1,MPI_DOUBLE,MPI_SUM,universe->uworld);
    }

    // if (v dot f) > 0:
    // v = (1-alpha) v + alpha |v| Fhat
    // |v| = length of v, Fhat = unit f
    // The modificatin of v is made wihtin the Verlet integration, after v update
    // if more than delaystep since v dot f was negative:
    // increase timestep and decrease alpha

    if (vdotfall > 0.0) {
      vdotv = 0.0;
      for (int i = 0; i < nlocal; i++)
        vdotv += v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2];
      MPI_Allreduce(&vdotv,&vdotvall,1,MPI_DOUBLE,MPI_SUM,world);

      // sum vdotv over replicas, if necessary
      // this communicator would be invalid for multiprocess replicas

      if (update->multireplica == 1) {
        vdotv = vdotvall;
        MPI_Allreduce(&vdotv,&vdotvall,1,MPI_DOUBLE,MPI_SUM,universe->uworld);
      }

      fdotf = 0.0;
      for (int i = 0; i < nlocal; i++)
        fdotf += f[i][0]*f[i][0] + f[i][1]*f[i][1] + f[i][2]*f[i][2];
      MPI_Allreduce(&fdotf,&fdotfall,1,MPI_DOUBLE,MPI_SUM,world);

      // sum fdotf over replicas, if necessary
      // this communicator would be invalid for multiprocess replicas

      if (update->multireplica == 1) {
        fdotf = fdotfall;
        MPI_Allreduce(&fdotf,&fdotfall,1,MPI_DOUBLE,MPI_SUM,universe->uworld);
      }

      scale1 = 1.0 - alpha;
      if (fdotfall <= 1e-20) scale2 = 0.0;
      else scale2 = alpha * sqrt(vdotvall/fdotfall);

      if (ntimestep - last_negative > delaystep) {
        dt = MIN(dt*dt_grow,dtmax);
        update->dt = dt;
        alpha *= alpha_shrink;
      }

    // else (v dot f) <= 0
    // if more than delaystep since starting the relaxation:
    // reset alpha
    //    if dt > dtmin:
    //    decrease timestep
    // set x(t) = x(t-0.5*dt)
    // set v = 0

    } else {
      last_negative = ntimestep;
      // Limit decrease of timestep
      if (ntimestep - ntimestep_fire > delaystep) {
        alpha = alpha0;
        if (dt > dtmin) {
          dt *= dt_shrink;
          update->dt = dt;
        }
      }
      double **x = atom->x;
      if (halfstepback_flag) {
        for (int i = 0; i < nlocal; i++) {
          x[i][0] -= 0.5 * dtv * v[i][0];
          x[i][1] -= 0.5 * dtv * v[i][1];
          x[i][2] -= 0.5 * dtv * v[i][2];
        }
      }
      for (int i = 0; i < nlocal; i++)
        v[i][0] = v[i][1] = v[i][2] = 0.0;
    }

    // limit timestep so no particle moves further than dmax

    double *rmass = atom->rmass;
    double *mass = atom->mass;
    int *type = atom->type;

    dtvone = dt;

    for (int i = 0; i < nlocal; i++) {
      vmax = MAX(fabs(v[i][0]),fabs(v[i][1]));
      vmax = MAX(vmax,fabs(v[i][2]));
      if (dtvone*vmax > dmax) dtvone = dmax/vmax;
    }
    MPI_Allreduce(&dtvone,&dtv,1,MPI_DOUBLE,MPI_MIN,world);

    // min dtv over replicas, if necessary
    // this communicator would be invalid for multiprocess replicas

    if (update->multireplica == 1) {
      dtvone = dtv;
      MPI_Allreduce(&dtvone,&dtv,1,MPI_DOUBLE,MPI_MIN,universe->uworld);
    }

    // Semi-implicit Euler integration

  if (integrator == 0) {

    dtf = dtv * force->ftm2v; 

    if (rmass) {
      for (int i = 0; i < nlocal; i++) {
        dtfm = dtf / rmass[i];
        v[i][0] += dtfm * f[i][0];
        v[i][1] += dtfm * f[i][1];
        v[i][2] += dtfm * f[i][2];
        if (vdotfall > 0.0) {
          v[i][0] = scale1*v[i][0] + scale2*f[i][0];
          v[i][1] = scale1*v[i][1] + scale2*f[i][1];
          v[i][2] = scale1*v[i][2] + scale2*f[i][2];
        }
        x[i][0] += dtv * v[i][0];
        x[i][1] += dtv * v[i][1];
        x[i][2] += dtv * v[i][2];
      }
    } else {
      for (int i = 0; i < nlocal; i++) {
        dtfm = dtf / mass[type[i]];
        v[i][0] += dtfm * f[i][0];
        v[i][1] += dtfm * f[i][1];
        v[i][2] += dtfm * f[i][2];
        if (vdotfall > 0.0) {
          v[i][0] = scale1*v[i][0] + scale2*f[i][0];
          v[i][1] = scale1*v[i][1] + scale2*f[i][1];
          v[i][2] = scale1*v[i][2] + scale2*f[i][2];
        }
        x[i][0] += dtv * v[i][0];
        x[i][1] += dtv * v[i][1];
        x[i][2] += dtv * v[i][2];
      }
    }

    eprevious = ecurrent;
    ecurrent = energy_force(0);
    neval++;


  // Velocity Verlet integration

  }else if (integrator == 1) {

    dtf = 0.5 * dtv * force->ftm2v;

    if (rmass) {
      for (int i = 0; i < nlocal; i++) {
          dtfm = dtf / rmass[i];
          v[i][0] += dtfm * f[i][0];
          v[i][1] += dtfm * f[i][1];
          v[i][2] += dtfm * f[i][2];
          if (vdotfall > 0.0) {
            v[i][0] = scale1*v[i][0] + scale2*f[i][0];
            v[i][1] = scale1*v[i][1] + scale2*f[i][1];
            v[i][2] = scale1*v[i][2] + scale2*f[i][2];
          }
          x[i][0] += dtv * v[i][0];
          x[i][1] += dtv * v[i][1];
          x[i][2] += dtv * v[i][2];
      }
    } else {
      for (int i = 0; i < nlocal; i++) {
          dtfm = dtf / mass[type[i]];
          v[i][0] += dtfm * f[i][0];
          v[i][1] += dtfm * f[i][1];
          v[i][2] += dtfm * f[i][2];
          if (vdotfall > 0.0) {
            v[i][0] = scale1*v[i][0] + scale2*f[i][0];
            v[i][1] = scale1*v[i][1] + scale2*f[i][1];
            v[i][2] = scale1*v[i][2] + scale2*f[i][2];
          }
          x[i][0] += dtv * v[i][0];
          x[i][1] += dtv * v[i][1];
          x[i][2] += dtv * v[i][2];
      }
    }
    
    eprevious = ecurrent;
    ecurrent = energy_force(0);
    neval++;

    if (rmass) {
      for (int i = 0; i < nlocal; i++) {
          dtfm = dtf / rmass[i];
          v[i][0] += dtfm * f[i][0];
          v[i][1] += dtfm * f[i][1];
          v[i][2] += dtfm * f[i][2];
        }
    } else {
      for (int i = 0; i < nlocal; i++) {
        dtfm = dtf / mass[type[i]];
        v[i][0] += dtfm * f[i][0];
        v[i][1] += dtfm * f[i][1];
        v[i][2] += dtfm * f[i][2];
      }
    }

  // Leap Frog integration

  }else if (integrator == 2) {

    dtf = dtv * force->ftm2v; 

    if (rmass) {
      for (int i = 0; i < nlocal; i++) {
        dtfm = dtf / rmass[i];
        v[i][0] += dtfm * f[i][0];
        v[i][1] += dtfm * f[i][1];
        v[i][2] += dtfm * f[i][2];
        if (vdotfall > 0.0) {
          v[i][0] = scale1*v[i][0] + scale2*f[i][0];
          v[i][1] = scale1*v[i][1] + scale2*f[i][1];
          v[i][2] = scale1*v[i][2] + scale2*f[i][2];
        }
        x[i][0] += dtv * v[i][0];
        x[i][1] += dtv * v[i][1];
        x[i][2] += dtv * v[i][2];
      }
    } else {
      for (int i = 0; i < nlocal; i++) {
        dtfm = dtf / mass[type[i]];
        v[i][0] += dtfm * f[i][0];
        v[i][1] += dtfm * f[i][1];
        v[i][2] += dtfm * f[i][2];
        if (vdotfall > 0.0) {
          v[i][0] = scale1*v[i][0] + scale2*f[i][0];
          v[i][1] = scale1*v[i][1] + scale2*f[i][1];
          v[i][2] = scale1*v[i][2] + scale2*f[i][2];
        }
        x[i][0] += dtv * v[i][0];
        x[i][1] += dtv * v[i][1];
        x[i][2] += dtv * v[i][2];
      }
    }

    eprevious = ecurrent;
    ecurrent = energy_force(0);
    neval++;

  // Standard Euler integration

  }else if (integrator == 3) {

    dtf = dtv * force->ftm2v; 

    if (rmass) {
      for (int i = 0; i < nlocal; i++) {
        dtfm = dtf / rmass[i];
        if (vdotfall > 0.0) {
          v[i][0] = scale1*v[i][0] + scale2*f[i][0];
          v[i][1] = scale1*v[i][1] + scale2*f[i][1];
          v[i][2] = scale1*v[i][2] + scale2*f[i][2];
        }
        x[i][0] += dtv * v[i][0];
        x[i][1] += dtv * v[i][1];
        x[i][2] += dtv * v[i][2];
        v[i][0] += dtfm * f[i][0];
        v[i][1] += dtfm * f[i][1];
        v[i][2] += dtfm * f[i][2];
      }
    } else {
      for (int i = 0; i < nlocal; i++) {
        dtfm = dtf / mass[type[i]];
        if (vdotfall > 0.0) {
          v[i][0] = scale1*v[i][0] + scale2*f[i][0];
          v[i][1] = scale1*v[i][1] + scale2*f[i][1];
          v[i][2] = scale1*v[i][2] + scale2*f[i][2];
        }
        x[i][0] += dtv * v[i][0];
        x[i][1] += dtv * v[i][1];
        x[i][2] += dtv * v[i][2];
        v[i][0] += dtfm * f[i][0];
        v[i][1] += dtfm * f[i][1];
        v[i][2] += dtfm * f[i][2];
      }
    }

    eprevious = ecurrent;
    ecurrent = energy_force(0);
    neval++;

  }

    // energy tolerance criterion
    // only check after delaystep elapsed since velocties reset to 0
    // sync across replicas if running multi-replica minimization

    if (update->etol > 0.0 && ntimestep-last_negative > delaystep) {
      if (update->multireplica == 0) {
        if (fabs(ecurrent-eprevious) <
            update->etol * 0.5*(fabs(ecurrent) + fabs(eprevious) + EPS_ENERGY)) {
          update->dt = dtinit;
          return ETOL;
        }
      } else {
        if (fabs(ecurrent-eprevious) <
            update->etol * 0.5*(fabs(ecurrent) + fabs(eprevious) + EPS_ENERGY))
          flag = 0;
        else flag = 1;
        MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_SUM,universe->uworld);
        if (flagall == 0) {
          update->dt = dtinit;
          return ETOL;
        }
      }
    }

    // force tolerance criterion
    // sync across replicas if running multi-replica minimization

    if (update->ftol > 0.0) {
      fdotf = fnorm_sqr();
      if (update->multireplica == 0) {
        if (fdotf < update->ftol*update->ftol) {
          update->dt = dtinit;
          return FTOL;
        }
      } else {
        if (fdotf < update->ftol*update->ftol) flag = 0;
        else flag = 1;
        MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_SUM,universe->uworld);
        if (flagall == 0) {
          update->dt = dtinit;
          return FTOL;
        }
      }
    }

    // output for thermo, dump, restart files

    if (output->next == ntimestep) {
      timer->stamp();
      output->write(ntimestep);
      timer->stamp(Timer::OUTPUT);
    }
  }

  update->dt = dtinit;
  return MAXITER;
}
