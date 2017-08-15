/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef MINIMIZE_CLASS

MinimizeStyle(adaptglok,MinAdaptGlok)

#else

#ifndef LMP_MIN_ADAPTGLOK_H
#define LMP_MIN_ADAPTGLOK_H

#include "min.h"

namespace LAMMPS_NS {

class MinAdaptGlok : public Min {
 public:
  MinAdaptGlok(class LAMMPS *);
  ~MinAdaptGlok() {}
  void init();
  void setup_style();
  void reset_vectors();
  void save_box_state();
  void relax_box();
  double min_energy(double *);
  int iterate(int);

 private:
  double dt,dtmax,dtmin,dtinit;
  double alpha;
  bigint last_negative,ntimestep_fire;
  class Compute *temperature,*pressure;
  double box0[3],boxlo[3],boxhi[3],h_inv[6],boxhi_previous[3];
  double pressure_current_s,pressure_previous_s,pressure_previous[3];
  double vol0,xprdinit,yprdinit,zprdinit,scale,pv2e,eng;
  int no_pressure_flag;
  double *pressure_current;
  double *hextra;
};

}

#endif
#endif
