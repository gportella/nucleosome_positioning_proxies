#ifndef NUC_ELASTIC_H
#define NUC_ELASTIC_H
#include "utils_common.hpp"
#include <Eigen/Dense>
#include <fstream>
#include <iostream>
#include <map>
#include <math.h>
#include <random>
#include <seqan/bed_io.h>
#include <seqan/file.h>
#include <seqan/find.h>
#include <seqan/modifier.h>
#include <seqan/seq_io.h>
#include <seqan/stream.h>
#include <sstream>
#include <stdlib.h>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#define NUC_LEN 147
#define NUC_CORE 74
#define ELASTIC_NORM_L 1000 // interpolate/extrapolate to 1000 data points
                            // the ouput profile
using namespace seqan;

typedef Eigen::Matrix<float, 6, 6> Matrix6f;
typedef Eigen::Matrix<float, 6, 1> Vector6f;

struct Tag_NucCore {};
struct Tag_ElProf {};

typedef struct my_bpmodel {
  Matrix6f fct;
  Vector6f eq;
} NNmodel;

typedef struct results_el_nucpredict {
  std::vector<double> e;
  std::vector<double> occ;
} el_nucpred;

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// split string based on delimiter
void split(const std::string &s, char delim, std::vector<std::string> &elems) {
  std::stringstream ss(s);
  std::string item;
  while (std::getline(ss, item, delim)) {
    elems.push_back(item);
  }
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// Check if file exists
inline bool file_exists(const std::string &name) {
  struct stat buffer;
  return (stat(name.c_str(), &buffer) == 0);
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// Computes the elastic energy of a whole nucleosome considering tetranucleotide
// steps everywhere except on the first/last basepairs. Those are computed using
// the dinucleotide model.
template <typename tt1, typename tt2, typename tt3, typename tt4>
double nucElastic(tt1 tetra_model, tt2 di_model, tt3 nucref, tt4 seq) {
  std::ostringstream sstrs;
  Vector6f sup;
  unsigned l_seq = length(seq);
  double energy = 0;
  for (unsigned long i = 0; i < l_seq - 3; ++i) {
    Infix<Dna5String>::Type inf = infix(seq, i, i + 4);
    sstrs << inf;
    sup = tetra_model[sstrs.str()].eq.transpose() - nucref.row(i + 1);
    energy += 0.5 * sup.transpose() * tetra_model[sstrs.str()].fct * sup;
    sstrs.str(std::string()); // clear the contents of sstr
  }
  // add first and last bp as dinucleotides
  sstrs.str(std::string()); // clear the contents of sstr
  Infix<Dna5String>::Type inf = infix(seq, 0, 2);
  sstrs << inf;
  sup = di_model[sstrs.str()].eq.transpose() - nucref.row(0);
  energy += 0.5 * sup.transpose() * di_model[sstrs.str()].fct * sup;
  inf = infix(seq, l_seq - 2, l_seq);
  sstrs.str(std::string()); // clear the contents of sstr
  sstrs << inf;
  sup = di_model[sstrs.str()].eq.transpose() - nucref.row(l_seq - 2);
  energy += 0.5 * sup.transpose() * di_model[sstrs.str()].fct * sup;
  sstrs.str(std::string()); // clear the contents of sstr
  return energy;            /// (double)length(seq);
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// Computes the elastic energy of a nucleosome considering tetranucleotide
// and only computing the energy from the core (74 base pairs around dyad).
template <typename tt1, typename tt2, typename tt3>
double nucElastic(tt1 tetra_model, tt2 nucref, tt3 seq,
                  Tag_NucCore const & /*Tag*/) {
  std::ostringstream sstrs;
  Vector6f sup;
  // core bp is hardcoded to 74, which is the value that has been reported to
  // determine binding of nucleosomes in-vitro, according to the Van Noort
  // sequence predictor paper. 10.1073/pnas.1205659109
  unsigned int core_b = (unsigned int)std::ceil(NUC_LEN / 2 - NUC_CORE / 2);
  unsigned int core_e = (unsigned int)std::ceil(NUC_LEN / 2 + NUC_CORE / 2);
  double energy = 0;
  for (unsigned i = core_b; i < core_e; ++i) {
    Infix<Dna5String>::Type inf = infix(seq, i, i + 4);
    sstrs << inf;
    sup = tetra_model[sstrs.str()].eq.transpose() - nucref.row(i + 1);
    energy += 0.5 * sup.transpose() * tetra_model[sstrs.str()].fct * sup;
    sstrs.str(std::string()); // clear the contents of sstr
  }
  return energy;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

Eigen::MatrixXf readnucProb(CharString filename) {

  // not even checking, reading as if it was ok
  Eigen::MatrixXf mat;
  std::vector<float> v;
  std::vector<std::vector<float>> m;
  float a = 0;
  float b = 0;
  std::ifstream fileIn;

  if (!open(fileIn, toCString(filename))) {
    std::cerr << "ERROR: Cound not open nucleosome probabilities.\n";
    exit(1);
  }
  try {
    while (!fileIn.eof()) {
      fileIn >> a >> b;
      v.push_back(a);
      v.push_back(b);
      m.push_back(v);
      v.clear();
    }
  } catch (Exception const &e) {
    std::cout << "ERROR: " << e.what() << std::endl;
    exit(1);
  }
  // convert to matrix
  mat.resize(length(m), 2);
  for (unsigned j = 0; j < length(m); ++j) {
    mat(j, 0) = m[j][0];
    mat(j, 1) = m[j][1];
  }
  return mat;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

Eigen::MatrixXf loadRefNuc() {

  std::string FILE_REFNUC = "refnuc_bp.dat";
  std::ifstream ref(FILE_REFNUC);
  std::string delimiter = " ";
  std::string line;
  Eigen::MatrixXf av;
  size_t pos = 0;
  int i = 0;
  int j = 0;

  if (!file_exists(FILE_REFNUC)) {
    std::cerr << "Error: File " << FILE_REFNUC << " does not exist"
              << std::endl;
    exit(1);
  }
  // a bit dangerous, we should make sure
  // that there are no missing values, or too many
  // as it will segfault. I don't want to put an if
  av.resize(NUC_LEN - 1, 6);
  while (std::getline(ref, line)) {
    while ((pos = line.find(delimiter)) != std::string::npos) {
      av(i, j) = std::stod(line.substr(0, pos));
      line.erase(0, pos + delimiter.length());
      j++;
    }
    av(i, j) = std::stod(line);
    i++;
    j = 0;
  }
  return av;
}

// std::map<std::string, NNmodel> loadBPModel() {
template <typename tt1, typename tt2>
bool loadBPModel(tt1 &bpmodel, tt2 FILE_BP) {

  // std::string FILE_BP = "stif_bsc1_k_avg_miniabc.dat";
  std::ifstream bp(FILE_BP);
  std::string line;
  // std::map<std::string, NNmodel> bpmodel;
  int count = 0;
  int token_c = 0;

  std::vector<std::string> tokens;

  size_t pos = 0;
  NNmodel bpm;
  std::string delimiter = "  ";
  std::string token;
  std::string key;

  if (!file_exists(FILE_BP)) {
    std::cerr << "Error: File " << FILE_BP << " does not exist" << std::endl;
    return false;
  }

  while (std::getline(bp, line)) {
    // std::istringstream iss(line);
    if (line.find(">") != std::string::npos) {
      key = line.substr(1);
      count = 0;
    } else {
      count++;
      if (count <= 6) {
        // The FCT
        while ((pos = line.find(delimiter)) != std::string::npos) {
          bpm.fct(count - 1, token_c) = std::stod(line.substr(0, pos));
          line.erase(0, pos + delimiter.length());
          token_c++;
        }
        bpm.fct(count - 1, token_c) = std::stod(line);
        token_c = 0;
      }
      if (count > 6) {
        // the remaining line are the equilibrium values
        token_c = 0;
        while ((pos = line.find(delimiter)) != std::string::npos) {
          bpm.eq(token_c) = std::stod(line.substr(0, pos));
          line.erase(0, pos + delimiter.length());
          token_c++;
        }
        bpm.eq(token_c) = std::stod(line);
        token_c = 0;

        // Now add the bpmodel to the map
        bpmodel[key] = bpm;
      }
    }
  }
  return true;
  // std::cout << "Before return" << std::endl;
  // return bpmodel;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// Write the results vector

template <typename tt1, typename tt2> void dumpResults(tt1 fout, tt2 res) {
  std::ofstream file;
  file.open(toCString(fout));
  for (unsigned i = 0; i < res.size(); ++i) {
    file << i << " " << res[i] << std::endl;
  }
  file.close();
}
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// Compute the minimum elastic energy of all possible nucleosomes, whole nuc

template <typename tt1, typename tt2, typename tt3, typename tt4>
double do_min_elastic(tt1 seq, tt2 tetra_model, tt3 di_model, tt4 nucref) {
  Infix<Dna5String>::Type seq_i;
  double min_elastic = 10000000; // overkill to use std::numeric_limits::max()
  for (unsigned i = 0; i < length(seq) - NUC_LEN + 1; ++i) {
    seq_i = infix(seq, i, i + NUC_LEN);
    double E_nuc = nucElastic(tetra_model, di_model, nucref, seq_i);
    if (min_elastic > E_nuc) {
      min_elastic = E_nuc;
    }
  }
  return min_elastic;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// only compute mimimum energy of the sequence, only nuc core

template <typename tt1, typename tt2, typename tt3>
double do_min_elastic(tt1 seq, tt2 tetra_model, tt3 nucref,
                      Tag_NucCore const & /*Tag*/) {
  Infix<Dna5String>::Type seq_i;
  double min_elastic = 10000000; // overkill to use std::numeric_limits::max()
  for (unsigned i = 0; i < length(seq) - NUC_LEN + 1; ++i) {
    seq_i = infix(seq, i, i + NUC_LEN);
    double E_nuc = nucElastic(tetra_model, nucref, seq_i, Tag_NucCore());
    if (min_elastic > E_nuc) {
      min_elastic = E_nuc;
    }
  }
  return min_elastic;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// Here for the whole profile, all the nucleosome

template <typename tt1, typename tt2, typename tt3, typename tt4, typename tt5>
el_nucpred do_prof_elastic(tt1 seq, tt2 tetra_model, tt3 di_model, tt4 nucref,
                           tt5 cond) {
  Infix<Dna5String>::Type seq_i;
  std::vector<double> E_profile;
  for (unsigned i = 0; i < length(seq) - NUC_LEN + 1; ++i) {
    seq_i = infix(seq, i, i + NUC_LEN);
    double E_nuc = nucElastic(tetra_model, di_model, nucref, seq_i);
    E_profile.push_back(E_nuc);
  }
  std::vector<double> prof_smoothed = smooth_box(E_profile, cond.smooth_window);
  // normalize by the nuc_len to avoid problems of overflows in exp
  std::transform(prof_smoothed.begin(), prof_smoothed.end(),
                 prof_smoothed.begin(),
                 std::bind2nd(std::divides<double>(), NUC_LEN));
  // vanderlick could overflow due to exponential calculation, catch it
  std::vector<double> P;
  try {
    // here P returns well padded, has the same length as bases in the orginal
    // sequence
    P = vanderlick(prof_smoothed, cond);
  } catch (const std::exception &e) {
    std::cerr << "Caught " << e.what() << std::endl;
    // rethrow
    throw std::overflow_error("overflow in vanderlick");
  }
  // convolute the provability (?) with footprint-1 (??) to get the occupancy
  std::vector<double> zeros(NUC_LEN - 1, 1.0);
  std::vector<double> N = conv_same(P, zeros);
  // add padding to E (recall it was smoothed), to print out
  int window = NUC_LEN;
  std::vector<double> e_padded = add_zeros_padding(prof_smoothed, window);

  el_nucpred el_results;
  el_results.e = e_padded;
  el_results.occ = N;

  return el_results;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// Here for the whole profile, nucleosome core only

template <typename tt1, typename tt2, typename tt3, typename tt4>
el_nucpred do_prof_elastic(tt1 seq, tt2 tetra_model, tt3 nucref, tt4 cond,
                           Tag_NucCore const & /*Tag*/) {
  Infix<Dna5String>::Type seq_i;
  std::vector<double> profile;
  for (unsigned i = 0; i < length(seq) - NUC_LEN + 1; ++i) {
    seq_i = infix(seq, i, i + NUC_LEN);
    double E_nuc = nucElastic(tetra_model, nucref, seq_i, Tag_NucCore());
    profile.push_back(E_nuc);
  }
  // smooths -- I borrow the smooth setting from vn_smooth_window, rather
  // than usic 10 as the magic number...
  std::vector<double> prof_smoothed = smooth_box(profile, cond.smooth_window);
  // normalize by the nuc_len to avoid problems of overflows in exp
  std::transform(prof_smoothed.begin(), prof_smoothed.end(),
                 prof_smoothed.begin(),
                 std::bind2nd(std::divides<double>(), NUC_CORE));
  // vanderlick could overflow due to exponential calculation, catch it
  std::vector<double> P;
  try {
    // here P returns well padded, has the same length as bases in the orginal
    // sequence
    P = vanderlick(prof_smoothed, cond);
  } catch (const std::exception &e) {
    std::cerr << "Caught " << e.what() << std::endl;
    // rethrow
    throw std::overflow_error("overflow in vanderlick");
  }
  // convolute the provability (?) with footprint-1 (??) to get the occupancy
  std::vector<double> zeros(NUC_CORE - 1, 1.0);
  std::vector<double> N = conv_same(P, zeros);
  // add padding to E (recall it was smoothed), to print out
  int window = NUC_CORE;
  std::vector<double> e_padded = add_zeros_padding(prof_smoothed, window);

  el_nucpred el_results;
  el_results.e = e_padded;
  el_results.occ = N;

  return el_results;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
template <typename tt1, typename tt2, typename tt3, typename tt4, typename tt5,
          typename tt6>
void do_all_elastic(tt1 tetra_model, tt2 di_model, tt3 nucref, tt4 seqs,
                    tt5 outfile, tt6 cond) {

  std::vector<double> min_elastic_v(length(seqs), 0.0);
  if (!cond.b_nuccore) {
#pragma omp parallel for
    for (unsigned i = 0; i < length(seqs); ++i) {
      if (length(seqs[i]) >= NUC_LEN && notNInside(seqs[i])) {
        double min_E = do_min_elastic(seqs[i], tetra_model, di_model, nucref);
        min_elastic_v[i] = min_E;
      }
    }
  } else {
#pragma omp parallel for
    for (unsigned i = 0; i < length(seqs); ++i) {
      if (length(seqs[i]) >= NUC_LEN && notNInside(seqs[i])) {
        double min_E =
            do_min_elastic(seqs[i], tetra_model, nucref, Tag_NucCore());
        min_elastic_v[i] = min_E;
      }
    }
  }

  // write the result
  dumpResults(outfile, min_elastic_v);
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
template <typename tt1, typename tt2, typename tt3, typename tt4, typename tt5,
          typename tt6>
void do_all_elastic(tt1 tetra_model, tt2 di_model, tt3 nucref, tt4 seqs,
                    tt5 outfile, tt6 cond, Tag_ElProf const & /*Tag*/) {

  std::vector<double> av_e_prof(ELASTIC_NORM_L, 0.0);
  std::vector<double> av_occ_prof(ELASTIC_NORM_L, 0.0);
  // std::vector<double> Z(ELASTIC_NORM_L, 0.0); // Z = partition function
  // x value of the interpolation, from 0 to 1 with NORM_OUT_L datapoints
  std::vector<double> x_inter = linspace(0.0, 1., ELASTIC_NORM_L);
  int count_curves = 0;

  if (!cond.b_nuccore) {
    //#pragma omp parallel for
    std::cout << "Going for whole nuc" << std::endl;
    for (unsigned i = 0; i < length(seqs); ++i) {

      if (length(seqs[i]) >= NUC_LEN && notNInside(seqs[i])) {
        std::cout << "Trying" << std::endl;
        el_nucpred el_results;
        // let's protect ourselfs from exp overflows in vanderlick
        try {

          // TODO:: add try/catch
          el_results =
              do_prof_elastic(seqs[i], tetra_model, di_model, nucref, cond);
          std::vector<double> E_prof = el_results.e;
          std::vector<double> occ_prof = el_results.occ;
          // sets all predictions on a common scale to average, in the
          // case that we get sequences of different length and it makes
          // sense to do so. Otherwise, just feed me seqs of the same length
          // get x coord of fe/occ normalized from 0 to 1
          std::vector<double> el_x_inter = linspace(0.0, 1., E_prof.size());
          // interpolates
          std::vector<double> interp_el_prof =
              interp_linear(x_inter, el_x_inter, E_prof);
          std::vector<double> interp_occ_prof =
              interp_linear(x_inter, el_x_inter, occ_prof);
          // let's do Boltzman averaging, compute the exponential first
          // defined in utils_common
          // TROUBLE is that due to padding of FE with 0s, the borders have a
          // weight of 1.... which means linear averaging... Somehow I still
          // have a big peak in some places.
          // we normalize first, to avoid overflows
          /*
          std::vector<double> exp_e = exp_vect(interp_el_prof);
          std::transform(exp_e.begin(), exp_e.end(), interp_el_prof.begin(),
                         interp_el_prof.begin(), std::multiplies<double>());
            */

          av_e_prof = brave_add_vector(av_e_prof, interp_el_prof);
          av_occ_prof = brave_add_vector(av_occ_prof, interp_occ_prof);
          // and for the partition function Z
          // Z = brave_add_vector(Z, exp_e);
          count_curves++;
        } catch (const std::exception &e) {
          std::cerr << "Caught rethrowing " << e.what() << ", skipping curve."
                    << std::endl;
          // skip this iteration
          continue;
        }
      }
    }
  } else {
    // DO CALC WITH NUCCORE
    //#pragma omp parallel for
    std::cout << "Going for nuccore" << std::endl;
    for (unsigned i = 0; i < length(seqs); ++i) {
      if (length(seqs[i]) >= NUC_LEN && notNInside(seqs[i])) {
        el_nucpred el_results;
        // let's protect ourselfs from exp overflows in vanderlick
        try {

          // TODO: copy the same as above (e.g. try/catch, el_results...)
          el_results = do_prof_elastic(seqs[i], tetra_model, nucref, cond,
                                       Tag_NucCore());
          std::vector<double> E_prof = el_results.e;
          std::vector<double> occ_prof = el_results.occ;
          // sets all predictions on a common scale to average, in the
          // case that we get sequences of different length and it makes
          // sense to do so. Otherwise, just feed me seqs of the same length
          // get x coord of fe/occ normalized from 0 to 1
          std::vector<double> el_x_inter = linspace(0.0, 1., E_prof.size());
          // interpolates
          std::vector<double> interp_el_prof =
              interp_linear(x_inter, el_x_inter, E_prof);
          std::vector<double> interp_occ_prof =
              interp_linear(x_inter, el_x_inter, occ_prof);
          /* We rather use the linear avearge, actually
          // see code above, Boltzman average
          std::transform(interp_el_prof.begin(), interp_el_prof.end(),
                         interp_el_prof.begin(),
                         std::bind2nd(std::divides<double>(), NUC_CORE));
          std::vector<double> exp_e = exp_vect(interp_el_prof);
          std::transform(exp_e.begin(), exp_e.end(), interp_el_prof.begin(),
                         interp_el_prof.begin(), std::multiplies<double>());
          */
          // adds to the vect containing the sum
          av_e_prof = brave_add_vector(av_e_prof, interp_el_prof);
          av_occ_prof = brave_add_vector(av_occ_prof, interp_occ_prof);
          // Z = brave_add_vector(Z, exp_e);
          count_curves++;
        } catch (const std::exception &e) {
          std::cerr << "Caught rethrowing " << e.what() << ", skipping curve."
                    << std::endl;
          // skip this iteration
          continue;
        }
      }
    }
  }
  // only writes the files if it found something worth analysing
  if (count_curves > 0) {
    // normalize the curves

    std::transform(av_e_prof.begin(), av_e_prof.end(), av_e_prof.begin(),
                   std::bind2nd(std::divides<double>(), count_curves));
    std::transform(av_occ_prof.begin(), av_occ_prof.end(), av_occ_prof.begin(),
                   std::bind2nd(std::divides<double>(), count_curves));
    /*
    // normalize by partition function Z
    std::transform(av_prof.begin(), av_prof.end(), Z.begin(),
    av_prof.begin(),
                   std::divides<double>());
      */
    // output them, prepend a fe_ / occ_ to the specified output filename
    seqan::CharString out_e_fn = "e_elastic_";
    out_e_fn += outfile;
    write_xy(out_e_fn, x_inter, av_e_prof);
    seqan::CharString out_occ_fn = "occ_elastic_";
    out_occ_fn += outfile;
    write_xy(out_occ_fn, x_inter, av_occ_prof);
  } else {
    std::cout << "Did not find a suitable sequence to analyse, zero output"
              << std::endl;
  }
}

#endif /* end protective inclusion */
