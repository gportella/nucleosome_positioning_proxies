///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
//
//  Predictor of nucleosome occupancy based on
//
//  "Sequence-based prediction of single nucleosome
//  positioning and genome-wide nucleosome occupancy",
//  van der Heijden et al.  DOI: 10.1073/pnas.1205659109
//
//  Based on the python implementation I got from J. v. Noort
//  It is a literal translation of his code. I checked it against
//  his results for a 2000 bp sequence, and the results match.
//
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

#ifndef V_NOORT_H
#define V_NOORT_H
#include "nuc_elastic.hpp"
#include <math.h>

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
//  Add window/2 zeros on either side of vector f
///////////////////////////////////////////////////////////////////////////////
template <typename tt1, typename tt2>
std::vector<tt1> add_zeros_padding(std::vector<tt1> const &f, tt2 win_l) {

  unsigned half_w = (unsigned)(ceil(win_l / 2.));
  std::vector<double> zeros(half_w, 0.0);
  std::vector<double> z_padded;
  z_padded.insert(std::end(z_padded), std::begin(zeros), std::end(zeros));
  z_padded.insert(std::end(z_padded), std::begin(f), std::end(f));
  z_padded.insert(std::end(z_padded), std::begin(zeros), std::end(zeros));
  return z_padded;
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// for quick&dirty debugging
////////////////////////////////////////////////////////////////////////////////
template <typename tt1> void print_debug(tt1 x) {

  for (auto &v : x) {
    std::cout << v << std::endl;
  }
  std::cout << "&" << std::endl;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// Convolution code, copied from:
// -> For 'valid' and 'full' directly from StackOverflow
// http://stackoverflow.com/a/24519913
// -> For 'same' it was addapted it from code found here
// moving from 2D to 1D (mind his code assumes vectorized 2D)
// https://github.com/jeremyfix/FFTConvolution
///////////////////////////////////////////////////////////////////////////////
template <typename T>
std::vector<T> conv_valid(std::vector<T> const &f, std::vector<T> const &g) {
  int const nf = f.size();
  int const ng = g.size();
  std::vector<T> const &min_v = (nf < ng) ? f : g;
  std::vector<T> const &max_v = (nf < ng) ? g : f;
  int const n = std::max(nf, ng) - std::min(nf, ng) + 1;
  std::vector<T> out(n, T());
  for (auto i(0); i < n; ++i) {
    for (int j(min_v.size() - 1), k(i); j >= 0; --j) {
      out[i] += min_v[j] * max_v[k];
      ++k;
    }
  }
  return out;
}
template <typename T>
std::vector<T> conv_full(std::vector<T> const &f, std::vector<T> const &g) {
  int const nf = f.size();
  int const ng = g.size();
  int const n = nf + ng - 1;
  std::vector<T> out(n, T());
  for (auto i(0); i < n; ++i) {
    int const jmn = (i >= ng - 1) ? i - (ng - 1) : 0;
    int const jmx = (i < nf - 1) ? i : nf - 1;
    for (auto j(jmn); j <= jmx; ++j) {
      out[i] += (f[j] * g[i - j]);
    }
  }
  return out;
}
// this returns a padding that is not exactly the same as numpy convolve,
// but it seems to work as well
template <typename T>
std::vector<T> conv_same(std::vector<T> const &f, std::vector<T> const &g) {

  int const n = (f.size() > g.size()) ? f.size() : g.size();
  std::vector<T> out(n, T());
  for (int i = 0; i < out.size(); ++i) {
    int low_k = std::max(0, i - (int)((g.size() - 1.0) / 2.0));
    int high_k = std::min((int)f.size(), i + (int)(g.size() / 2));
    double temp = 0.0;
    for (int k = low_k; k <= high_k; ++k) {
      temp += f[k] * g[(i - k + (int)(g.size() / 2.0))];
    }
    out[i] = temp;
  }
  return out;
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// this should do exactly the same smoothing as vanNoort
////////////////////////////////////////////////////////////////////////////////
template <typename tt1, typename tt2>
std::vector<double> smooth_vn(tt1 x, tt2 w_len) {

  // image the extremes of x : create mirrored copies of the ends
  // and paste together
  std::vector<double> b(w_len - 1);
  for (int i = 0; i < w_len - 1; ++i) {
    b[i] = x[w_len - 1 - i];
  }
  std::vector<double> e(w_len - 1);
  for (unsigned i = 0; i < w_len - 1; ++i) {
    e[i] = x[x.size() - 1 - i];
  }
  std::vector<double> s;
  s.insert(std::end(s), std::begin(b), std::end(b));
  s.insert(std::end(s), std::begin(x), std::end(x));
  s.insert(std::end(s), std::begin(e), std::end(e));

  // convolute extended curve with square window
  double norm_weight = 1.0 / (double)w_len;
  std::vector<double> w(w_len, norm_weight);
  std::vector<double> y = conv_valid(s, w);
  std::vector<double> smoothed(x.size());
  // return ignoring the added w_len - 1
  for (unsigned i = 0; i < x.size(); ++i) {
    smoothed[i] = y[w_len - 1 + i];
  }
  return smoothed;
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// Returns weight according to the prescribed periodicity functional form
// by adding div and changing the sing of amplitude b you can generalize
// the functional form for every dinucleotide (except CA, which is ctt)
////////////////////////////////////////////////////////////////////////////////
template <typename tt1, typename tt2, typename tt3, typename tt4>
std::vector<double> return_dinuc_weight(tt1 w, tt2 b, tt3 p, tt4 div) {

  std::vector<double> vect(w, 0);
  for (unsigned i = 0; i < w; ++i) {
    vect[i] = 0.25 + b * sin(2 * M_PI * i / p) / div;
  }
  return vect;
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// Probabiliis of dinucleotides based on vanNoort
////////////////////////////////////////////////////////////////////////////////
template <typename tt1>
std::vector<std::vector<std::vector<double>>> getweights(tt1 cond) {

  unsigned window = 74;
  float period = 10.2;
  float amp = 0.2;
  std::vector<double> AA = return_dinuc_weight(window, amp, period, 1.0);
  std::vector<double> AC = return_dinuc_weight(window, -amp, period, 3.0);
  std::vector<double> AG = AC;
  std::vector<double> AT = AC;

  std::vector<double> CA(window, 0.25);
  std::vector<double> CC = CA;
  std::vector<double> CG = CA;
  std::vector<double> CT = CA;

  std::vector<double> GA = return_dinuc_weight(window, amp, period, 3.0);
  std::vector<double> GC = return_dinuc_weight(window, -amp, period, 1.0);
  std::vector<double> GG = GA;
  std::vector<double> GT = GA;

  std::vector<double> TA = return_dinuc_weight(window, amp, period, 1.0);
  std::vector<double> TC = return_dinuc_weight(window, -amp, period, 1.0);
  std::vector<double> TG = TC;
  std::vector<double> TT = TA;

  // the order is really important
  std::vector<std::vector<double>> v_As = {AA, AC, AG, AT};
  std::vector<std::vector<double>> v_Cs = {CA, CC, CG, CT};
  std::vector<std::vector<double>> v_Gs = {GA, GC, GG, GT};
  std::vector<std::vector<double>> v_Ts = {TA, TC, TG, TT};
  // first index of vector[x][y][s] gives you the first dinucleotide
  // the second index the following dinucleotide, and the [s] goes
  // along the window
  std::vector<std::vector<std::vector<double>>> weights = {v_As, v_Cs, v_Gs,
                                                           v_Ts};

  return weights;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// Calculation of the free energy based on probability to find nuc dyad on a
// given base.
///////////////////////////////////////////////////////////////////////////////
template <typename tt1, typename tt2>
std::vector<double> calcE_vn(tt1 seq, tt2 cond) {
  unsigned window = 74;
  // first one is special due to pbc and averaging done in original script
  // I want to avoid using anf if inside the loop, so I just make it explicit
  // notice how the two strands are slanted // offset by one in original code
  // very verbose, but better than an if at each step of the loop
  // for clarity perhaps one could use a modulo or sth, to get the last
  // element of the sequence... would make the code nicer
  std::vector<std::vector<std::vector<double>>> weights = getweights(cond);
  std::vector<double> pf(length(seq) - window);
  std::vector<double> pr(length(seq) - window);
  double ps_f = 1.0;
  double ps_r = 1.0;
  unsigned ii = (unsigned)ordValue(seq[length(seq) - 1]);
  unsigned jj = (unsigned)ordValue(seq[0]);
  ps_f *= weights[ii][jj][0];
  unsigned ri = 3 - (unsigned)ordValue(seq[window]);
  unsigned rj = 3 - (unsigned)ordValue(seq[window - 1]);
  ps_r *= weights[ri][rj][0];
  for (unsigned s = 1; s < window; ++s) {
    unsigned ii = (unsigned)ordValue(seq[s - 1]);
    unsigned jj = (unsigned)ordValue(seq[s]);
    ps_f *= weights[ii][jj][s];
    unsigned ri = 3 - (unsigned)ordValue(seq[window - s]);
    unsigned rj = 3 - (unsigned)ordValue(seq[window - s - 1]);
    ps_r *= weights[ri][rj][s];
  }
  pf[0] = ps_f;
  pr[0] = ps_r;

  // now proceed from 1 onward
  for (unsigned i = 1; i < length(seq) - window; ++i) {
    double ps_f = 1.0;
    double ps_r = 1.0;
    for (unsigned s = 0; s < window; ++s) {
      // ordValue gives A:0, C:1, G:2, T:3, which matches order of indices
      // in weights (see getweights function)
      unsigned ii = (unsigned)ordValue(seq[i + s - 1]);
      unsigned jj = (unsigned)ordValue(seq[i + s]);
      ps_f *= weights[ii][jj][s];
      // 3 - ordValue(base) gives the ordValue of base pair
      unsigned ri = 3 - (unsigned)ordValue(seq[i + window - s]);
      unsigned rj = 3 - (unsigned)ordValue(seq[i + window - s - 1]);
      ps_r *= weights[ri][rj][s];
    }
    pf[i] = ps_f;
    pr[i] = ps_r;
  }
  // multiples all by 4^window
  for (unsigned i = 0; i < pf.size(); ++i) {
    pf[i] *= std::pow(4.L, window);
    pr[i] *= std::pow(4.L, window);
  }
  // circular shifts pr
  // overall, the fwd and rev strand will have the nuc center
  // displaced 2bp, the average will center them in the middle
  for (unsigned i = 0; i < pr.size() - 1; ++i) {
    auto p = pr[i];
    pr[i] = pr[i + 1];
    pr[i + 1] = p;
  }
  // Boltzman averages both strands
  std::vector<double> EE;
  for (unsigned i = 0; i < pr.size(); ++i) {
    double E = (pr[i] * log(pr[i]) + pf[i] * log(pf[i])) / (pr[i] + pf[i]);
    EE.push_back((double)E);
  }
  // smooths
  std::vector<double> E_smoothed = smooth_vn(EE, 10);
  return E_smoothed;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// Computes Vanderlick's solution to the Percus equation
// Again, 1-to-1 copy of van Noort's code
///////////////////////////////////////////////////////////////////////////////
template <typename tt1, typename tt2>
std::vector<double> vanderlick_vn(tt1 E, tt2 cond) {
  double mu = -1.5;
  int window = 74;
  int footprint = NUC_LEN;
  std::vector<double> E_out;
  for (const auto &e : E) {
    E_out.push_back(e - mu);
  }
  std::vector<double> fwd(E.size(), 0);
  for (unsigned i = 0; i < E.size(); ++i) {
    double tmp = 0.0;
    int i_max = std::max(((int)i - footprint), 0);
    for (unsigned j = i_max; j < i; ++j) {
      tmp += fwd[j];
    }
    fwd[i] = exp(E_out[i] - tmp);
  }
  std::vector<double> bwd(E.size(), 0.0);
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
  // bwd is reversed in place, we don't need it anymore
  std::reverse(bwd.begin(), bwd.end());
  for (unsigned i = 0; i < E.size(); ++i) {
    P[i] = fwd[i] * bwd[i];
  }
  // add padding windows/2 at begining and end to center nucleosome center
  std::vector<double> P_final = add_zeros_padding(P, window);
  return P_final;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
//  For a given sequence computes nucleosome occupancy and energy based on van
//  Noort's approach.
///////////////////////////////////////////////////////////////////////////////
template <typename tt1, typename tt2> void do_vannoort(tt1 seq, tt2 cond) {

  // change this magic number to be obtained from cond
  int window = 74;
  //  Note that the energy is smoothed before being returned.
  //  but the data is not padded (e.g. lacks window/2 on either side)
  std::vector<double> E = calcE_vn(seq, cond);
  // here P returns well padded, has the same length as bases in the orginal seq
  std::vector<double> P = vanderlick_vn(E, cond);
  // convolute the provability (?) with footprint-1 (??) to get the occupancy
  std::vector<double> zeros(NUC_LEN - 1, 1.0);
  std::vector<double> N = conv_same(P, zeros);
  // add padding to E (recall it was smoothed), to print out
  std::vector<double> e_padded = add_zeros_padding(E, window);
  print_debug(N);
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
//  For each sequece get prediction of nuc occupancy and energy
///////////////////////////////////////////////////////////////////////////////
template <typename tt1, typename tt2> void do_all_vannoort(tt1 seqs, tt2 cond) {
#pragma omp parallel for
  for (unsigned i = 0; i < length(seqs); ++i) {
    if (length(seqs[i]) >= NUC_LEN) {
      do_vannoort(seqs[i], cond);
    }
  }
}

#endif /* end protective inclusion */
