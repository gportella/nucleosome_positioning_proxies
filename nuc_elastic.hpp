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

using namespace seqan;

typedef Eigen::Matrix<float, 6, 6> Matrix6f;
typedef Eigen::Matrix<float, 6, 1> Vector6f;

typedef struct my_bpmodel {
  Matrix6f fct;
  Vector6f eq;
} tetra;

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
template <typename tt1, typename tt2, typename tt3>
double nucElastic(tt1 bpmodel, tt2 nucref, tt3 seq) {
  std::ostringstream sstrs;
  Vector6f sup;
  double energy = 0;
  for (int i = 0; i < length(seq) - 3; ++i) {
    Infix<Dna5String>::Type inf = infix(seq, i, i + 4);
    sstrs << inf;
    sup = bpmodel[sstrs.str()].eq.transpose() - nucref.row(i + 1);
    energy += 0.5 * sup.transpose() * bpmodel[sstrs.str()].fct * sup;
    sstrs.str(std::string()); // clear the contents of sstr
  }
  return energy; /// (double)length(seq);
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

  // a bit dangerous, we should make sure
  // that there are no missing values, or too many
  // as it will segfault. I don't want to put an if
  av.resize(146, 6);
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

// std::map<std::string, tetra> loadBPModel() {
bool loadBPModel(std::map<std::string, tetra> &bpmodel) {

  std::string FILE_BP = "stif_bsc1_k_avg_miniabc.dat";
  std::ifstream bp(FILE_BP);
  std::string line;
  // std::map<std::string, tetra> bpmodel;
  int count = 0;
  int token_c = 0;

  std::vector<std::string> tokens;

  size_t pos = 0;
  tetra bpm;
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
template <typename tt1, typename tt2, typename tt3, typename tt4>
double do_all_elastic(tt1 bpmodel, tt2 nucref, tt3 seqs, tt4 cond) {}
