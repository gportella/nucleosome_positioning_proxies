//////////////////////////////////////////////////////////////////
//
//  Predictor of nucleosome occupancy based on
//  Based on "Sequence-based prediction of single nucleosome
//  positioning and genome-wide nucleosome occupancy",
//  van der Heijden et al.  DOI: 10.1073/pnas.1205659109
//
//  Based on the python implementation I got from J. v. Noort
//  It is amlmost a literal translation of his code.
//
//  TODO: remove Eigen vectors and change to std::vector, no
//  gain from using Eigen, and one could drop the dependency
//  it the library was to be called from a code that does not
//  require Eigen.
//
//////////////////////////////////////////////////////////////////

#ifndef V_NOORT_H
#define V_NOORT_H
#include "nuc_elastic.hpp"
#include <math.h>

// class to perform smoothing, shamelessly from StackOverflow
// http://stackoverflow.com/a/12973856
// You might want to implement another smoothing later on, e.g.
// low pass filter implemented in van Noorts python code
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

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

template <typename tt1, typename tt2>
std::vector<double> calcE_vn(tt1 seq, tt2 cond) {
  unsigned window = 74;
  // first one is special due to pbc and weird averaging done in original script
  // I want to avoid using anf if inside the loop, so I just make it explicit
  // notice how the two strands are slanted // offset by one ... no idea why
  // very verbose, but better than an if at each step of the loop
  // for clarity perhaps one could use a modulo or sth, to get the last
  // element of the sequence... would make the code nicer
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
      // ordValue gives A:0, C:1, G:2, T:3
      unsigned ii = (unsigned)ordValue(seq[i + s - 1]);
      unsigned jj = (unsigned)ordValue(seq[i + s]);
      ps_f *= weights[ii][jj](s);
      // 3 - ordValue(base) gives the ordValue of pair
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
  // overall, the fwd and rev strand will have the nuc center
  // displaced 2bp, the average will center them in the middle
  for (unsigned i = 0; i < pr.size() - 1; ++i) {
    auto p = pr(i);
    pr(i) = pr(i + 1);
    pr(i + 1) = p;
  }
  // Boltzman average both strands
  std::vector<double> EE;
  for (unsigned i = 0; i < pr.size(); ++i) {
    double E = (pr(i) * log(pr(i)) + pf(i) * log(pf(i))) / (pr(i) + pf(i));
    EE.push_back((double)E);
  }
  // smooth
  boxFIR box(10);
  box.filter(EE);
  // add padding windows/2 at begining and end to center nucleosome center
  // at the right place
  unsigned half_w = (unsigned)(ceil(window / 2.));
  std::vector<double> zeros(half_w, 0.0);
  std::vector<double> E_final;
  E_final.insert(std::end(E_final), std::begin(zeros), std::end(zeros));
  E_final.insert(std::end(E_final), std::begin(EE), std::end(EE));
  E_final.insert(std::end(E_final), std::begin(zeros), std::end(zeros));

  return E_final;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

template <typename tt1, typename tt2>
std::vector<double> vanderlick_vn(tt1 E, tt2 cond) {
  // computes Vanderlick's solution to the Percus equation
  // as said, the code is a 1-to-1 copy of van Noort's code
  double mu = -1.5;
  int footprint = NUC_LEN;
  std::vector<double> E_out;
  for (const auto &e : E) {
    E_out.push_back(e - mu);
  }
  std::vector<double> fwd(E.size(), 0);
  for (unsigned i = 0; i < E.size(); ++i) {
    double tmp = 0;
    int i_max = std::max(((int)i - footprint), 0);
    for (unsigned j = i_max; j < i; ++j) {
      tmp += fwd[j];
    }
    fwd[i] = exp(E_out[i] - tmp);
  }
  std::vector<double> bwd(E.size(), 0);
  // reverse fwd vector
  std::vector<double> r_fwd(fwd.rbegin(), fwd.rend());
  for (unsigned i = 0; i < E.size(); ++i) {
    double tmp = 0;
    int i_max = std::max(((int)i - footprint), 0);
    for (unsigned j = i_max; j < i; ++j) {
      tmp += r_fwd[j] * bwd[j];
    }
    bwd[i] = 1.0 - tmp;
  }
  std::vector<double> P(E.size(), 0);
  // bwd is iterated backwards
  for (unsigned i = 0; i < E.size(); ++i) {
    P[i] = fwd[i] * bwd[E.size() - i - 1];
  }
  return P;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

template <typename tt1, typename tt2>
std::vector<double> do_vannoort(tt1 seq, tt2 cond) {
  std::vector<double> E = calcE_vn(seq, cond);
  std::vector<double> P = vanderlick_vn(E, cond);
  return P;
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
