//////////////////////////////////////////////////////////////////
//
//  Predictor of nucleosome occupancy based on
//  Based on "Sequence-based prediction of single nucleosome
//  positioning and genome-wide nucleosome occupancy",
//  van der Heijden et al.  DOI: 10.1073/pnas.1205659109
//
//  Based on the python implementation I got from J. v. Noort
//
//////////////////////////////////////////////////////////////////

#ifndef V_NOORT_H
#define V_NOORT_H
#include "nuc_elastic.hpp"
#include <math.h>

#define PI 3.14159265

template <typename tt1, typename tt2, typename tt3, typename tt4>
Eigen::VectorXf return_dinuc_weight(tt1 w, tt2 b, tt3 p, tt4 div = 1) {
  Eigen::VectorXf vect(w);
  for (unsigned i = 0; i < w; ++i) {
    vect(i) = 0.25 + b * sin(2 * PI * i / p) / div;
  }

  return vect;
}

template <typename tt1, typename tt2> void do_vannoort(tt1 seqs, tt2 cond) {
  // There might not be an actuall reason to use Eigen here
  // but since I need it anyway for nuc_elastic, I use it

  unsigned window = 74;
  float period = 10.1;
  float amp = 0.2;
  Eigen::VectorXf AA = return_dinuc_weight(window, amp, period, 1);
  Eigen::VectorXf AC = return_dinuc_weight(window, -amp, period, 3);
  Eigen::VectorXf AG = AC;
  Eigen::VectorXf AT = AC;

  // the one below is actually 0.25 for all, decide if you want to initialize
  // to that, or compute 0*sinus to have a more compact way. I hope the
  // compiler catches that b=0
  Eigen::VectorXf CA = return_dinuc_weight(window, 0, period, 1);
  Eigen::VectorXf CC = CA;
  Eigen::VectorXf CG = CA;
  Eigen::VectorXf CT = CA;

  Eigen::VectorXf GA = return_dinuc_weight(window, amp, period, 3);
  Eigen::VectorXf GC = return_dinuc_weight(window, -amp, period, 1);
  Eigen::VectorXf GG = GA;
  Eigen::VectorXf GT = GA;

  Eigen::VectorXf TA = return_dinuc_weight(window, amp, period, 1);
  Eigen::VectorXf TC = return_dinuc_weight(window, -amp, period, 1);
  Eigen::VectorXf TG = TC;
  Eigen::VectorXf TT = TA;

  std::vector<Eigen::VectorXf> vects_As = {AA, AC, AG, AT};
  std::vector<Eigen::VectorXf> vects_Cs = {CA, CC, CG, CT};
  std::vector<Eigen::VectorXf> vects_Gs = {GA, GC, GG, GT};
  std::vector<Eigen::VectorXf> vects_Ts = {TA, TC, TG, TT};

  for (unsigned i = 0; i < window; ++i) {
    std::cout << vects_As[1](i) << std::endl;
  }
}

#endif /* end protective inclusion */
