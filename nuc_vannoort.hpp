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

//#define PI 3.1415926535897932384626L

class boxFIR {
  int numCoeffs;         // MUST be > 0
  std::vector<double> b; // Filter coefficients
  std::vector<double> m; // Filter memories

public:
  boxFIR(int _numCoeffs) : numCoeffs(_numCoeffs) {
    if (numCoeffs < 1)
      numCoeffs = 1; // Must be > 0 or bad stuff happens

    double val = 1. / numCoeffs;
    for (int ii = 0; ii < numCoeffs; ++ii) {
      b.push_back(val);
      m.push_back(0.);
    }
  }

  void filter(std::vector<double> &a) {
    double output;

    for (int nn = 0; nn < a.size(); ++nn) {
      // Apply smoothing filter to signal
      output = 0;
      m[0] = a[nn];
      for (int ii = 0; ii < numCoeffs; ++ii) {
        output += b[ii] * m[ii];
      }

      // Reshuffle memories
      for (int ii = numCoeffs - 1; ii != 0; --ii) {
        m[ii] = m[ii - 1];
      }
      a[nn] = output;
    }
  }
};

template <typename tt1, typename tt2, typename tt3, typename tt4>
Eigen::VectorXd return_dinuc_weight(tt1 w, tt2 b, tt3 p, tt4 div) {
  Eigen::VectorXd vect(w);
  for (unsigned i = 0; i < w; ++i) {
    vect(i) = 0.25 + b * sin(2 * M_PI * i / p) / div;
  }

  return vect;
}

template <typename tt1>
std::vector<std::vector<Eigen::VectorXd>> getweights(tt1 cond) {
  // There might not be an actuall reason to use Eigen here
  // but since I need it anyway for nuc_elastic, I use it

  unsigned window = 74;
  float period = 10.1;
  float amp = 0.2;
  Eigen::VectorXd AA = return_dinuc_weight(window, amp, period, 1.0);
  Eigen::VectorXd AC = return_dinuc_weight(window, -amp, period, 3.0);
  Eigen::VectorXd AG = AC;
  Eigen::VectorXd AT = AC;

  // the one below is actually 0.25 for all, decide if you want to initialize
  // to that, or compute 0*sinus to have a more compact way. I hope the
  // compiler catches that b=0
  Eigen::VectorXd CA = return_dinuc_weight(window, 0, period, 1.0);
  Eigen::VectorXd CC = CA;
  Eigen::VectorXd CG = CA;
  Eigen::VectorXd CT = CA;

  Eigen::VectorXd GA = return_dinuc_weight(window, amp, period, 3.0);
  Eigen::VectorXd GC = return_dinuc_weight(window, -amp, period, 1.0);
  Eigen::VectorXd GG = GA;
  Eigen::VectorXd GT = GA;

  Eigen::VectorXd TA = return_dinuc_weight(window, amp, period, 1.0);
  Eigen::VectorXd TC = return_dinuc_weight(window, -amp, period, 1.0);
  Eigen::VectorXd TG = TC;
  Eigen::VectorXd TT = TA;

  std::vector<Eigen::VectorXd> v_As = {AA, AC, AG, AT};
  std::vector<Eigen::VectorXd> v_Cs = {CA, CC, CG, CT};
  std::vector<Eigen::VectorXd> v_Gs = {GA, GC, GG, GT};
  std::vector<Eigen::VectorXd> v_Ts = {TA, TC, TG, TT};
  std::vector<std::vector<Eigen::VectorXd>> weights = {v_As, v_Cs, v_Gs, v_Ts};

  return weights;
}

template <typename tt1, typename tt2>
Eigen::VectorXd do_vannoort(tt1 seq, tt2 cond) {
  unsigned window = 74;
  // first one is special due to pbc and weird averaging done in original script
  // I want to avoid using anf if inside the loop, so I just make it explicit
  // notice how the two strands are slanted // offset by one ... no idea why
  std::vector<std::vector<Eigen::VectorXd>> weights = getweights(cond);
  Eigen::VectorXd pf(length(seq) - window);
  Eigen::VectorXd pr(length(seq) - window);
  double ps_f = 1.0;
  double ps_r = 1.0;
  unsigned ii = (unsigned)ordValue(seq[length(seq) - 1]);
  unsigned jj = (unsigned)ordValue(seq[0]);
  ps_f *= weights[ii][jj](0);
  unsigned ri = 3 - (unsigned)ordValue(seq[window]);
  unsigned rj = 3 - (unsigned)ordValue(seq[window - 1]);
  ps_r *= weights[ri][rj](0);
  for (unsigned s = 1; s < window; ++s) {
    unsigned ii = (unsigned)ordValue(seq[s - 1]);
    unsigned jj = (unsigned)ordValue(seq[s]);
    ps_f *= weights[ii][jj](s);
    unsigned ri = 3 - (unsigned)ordValue(seq[window - s]);
    unsigned rj = 3 - (unsigned)ordValue(seq[window - s - 1]);
    ps_r *= weights[ri][rj](s);
  }
  pf(0) = ps_f;
  pr(0) = ps_r;

  // now proceed from 1 onward
  for (unsigned i = 1; i < length(seq) - window; ++i) {
    double ps_f = 1.0;
    double ps_r = 1.0;
    for (unsigned s = 0; s < window; ++s) {
      unsigned ii = (unsigned)ordValue(seq[i + s - 1]);
      unsigned jj = (unsigned)ordValue(seq[i + s]);
      ps_f *= weights[ii][jj](s);
      unsigned ri = 3 - (unsigned)ordValue(seq[i + window - s]);
      unsigned rj = 3 - (unsigned)ordValue(seq[i + window - s - 1]);
      ps_r *= weights[ri][rj](s);
    }
    pf(i) = ps_f;
    pr(i) = ps_r;
  }
  // multiply all by 4^window
  for (unsigned i = 0; i < pf.size(); ++i) {
    pf(i) *= std::pow(4.L, window);
    pr(i) *= std::pow(4.L, window);
  }
  // circular shift pr
  for (unsigned i = 0; i < pr.size() - 1; ++i) {
    auto p = pr(i);
    pr(i) = pr(i + 1);
    pr(i + 1) = p;
  }
  // Average both strands

  Eigen::VectorXd E(length(seq) - window);
  std::vector<double> EE;
  for (unsigned i = 0; i < pr.size(); ++i) {
    E(i) = (pr(i) * log(pr(i)) + pf(i) * log(pf(i))) / (pr(i) + pf(i));
    EE.push_back((double)E(i));
  }
  boxFIR box(10);
  box.filter(EE);
  std::cout << "Unfiltered" << std::endl;
  for (unsigned i = 0; i < pr.size(); ++i) {
    std::cout << E(i) << std::endl;
  }
  std::cout << "Filtered" << std::endl;
  for (unsigned i = 0; i < pr.size(); ++i) {
    std::cout << EE[i] << std::endl;
  }
  return E;
}

template <typename tt1, typename tt2> void do_all_vannoort(tt1 seqs, tt2 cond) {

  //#pragma omp parallel for
  for (unsigned i = 0; i < length(seqs); ++i) {
    if (length(seqs[i]) >= NUC_LEN) {
      do_vannoort(seqs[i], cond);
    }
  }
}

#endif /* end protective inclusion */
