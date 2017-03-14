

#ifndef DO_PERIODIC_H
#define DO_PERIODIC_H
#include <algorithm> //for std::generate_n
#include <boost/regex.hpp>
#include <cfloat> // DBL_MAX
#include <fstream>
#include <functional> //for std::function
#include <iostream>
#include <math.h>
#include <omp.h>
#include <random>
#include <seqan/file.h>
#include <seqan/find.h>
#include <seqan/modifier.h>
#include <seqan/seq_io.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>
#include <sstream>
#include <stdlib.h>
#include <string>
#include <tuple>

using namespace seqan;

///////////////////////////////////////////////////////////////////////////////
// :Info:
///////////////////////////////////////////////////////////////////////////////
// TODO
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

template <typename tt1, typename tt2, typename tt3>
void dumpresults(tt1 fout, tt2 cond, tt3 prob) {
  std::ofstream file;

  file.open(toCString(fout));

  for (unsigned i = 0; i < cond.l; ++i) {
    file << i << " " << prob[i] << std::endl;
  }

  file.close();
}

/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
// borrowed from my tfo_finder code
template <typename TText1>
std::vector<unsigned int> search_by_seq(TText1 const &sequence,
                                        boost::regex const &re) {
  // keep the index of possible tfos in a tuple list
  std::vector<unsigned int> limits;

  // 1) not sure if ther is a better way to get a string from the sequence
  // record, at the moment I write it to a stream and from there to a string
  std::stringstream support_stream;
  support_stream << sequence;
  std::string str_sequence = support_stream.str();
  try {
    // notice that the pattern re is defined in the calling function
    // search, with the regular expression
    boost::sregex_iterator next(str_sequence.begin(), str_sequence.end(), re);
    boost::sregex_iterator end; // on default constructor to check against end
    while (next != end) {
      boost::smatch match = *next; // dereference match and write it as string
                                   //    cout << match.str() << endl;
      int init = match.position();
      limits.emplace_back(init);
      next++;
    }
  } catch (boost::regex_error &e) {
    // Syntax error in the regular expression
    std::cout << "Something went wrong setting up the regular_expression"
              << std::endl;
  }
  return limits;
}

/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
template <typename tt1, typename tt2, typename tt3>
void do_seq_periodic(tt1 seq, tt2 &dif_vector, tt3 cond) {
  // change 50 to the cutoff
  static const boost::regex re("(?=(AA|TT|TA))");
  std::vector<unsigned int> dinuc_pos = search_by_seq(seq, re);

  for (unsigned i = 0; i < dinuc_pos.size(); ++i) {
    for (unsigned j = i + 1;
         j < dinuc_pos.size() && (dinuc_pos[j] - dinuc_pos[i] < cond.cutoff);
         ++j) {
      unsigned dif = dinuc_pos[j] - dinuc_pos[i];
      dif_vector[dif] = dif_vector[dif] + (1.0 / (length(seq) - dif));
    }
  }
}

/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

template <typename tt1, typename tt2> void do_all_periodic(tt1 seqs, tt2 cond) {

  std::vector<double> dif_vector(cond.cutoff, 0.0);
  Iterator<StringSet<Dna5String>>::Type it = begin(seqs);
  Iterator<StringSet<Dna5String>>::Type itEnd = end(seqs);
  for (; it != itEnd; ++it) {
    do_seq_periodic(*it, dif_vector, cond);
  }
  unsigned i = 0;
  for (const auto &val : dif_vector) {
    // unsigned normalize = cond.cutoff - i;
    std::cout << i << " " << val / ((float)length(seqs)) << std::endl;
    i++;
  }
}

#endif /* end protective inclusion */
