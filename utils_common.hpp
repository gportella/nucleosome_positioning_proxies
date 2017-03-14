
#ifndef PERIOD_ELASTIC_H
#define PERIOD_ELASTIC_H

#include <functional>
#include <iostream>
#include <math.h>
#include <seqan/arg_parse.h>
#include <seqan/file.h>
#include <seqan/find.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>
#include <stdlib.h>
#include <string>
#include <vector>

using namespace seqan;

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
template <typename tt1> bool checkArguments(tt1 cond) {
  if (cond.cutoff > 147) {
    std::cerr << "A cutoff longer than the length of a nucleosome does not "
                 "make sense. "
              << std::endl;
    exit(0);
  } else {
    return true;
  }
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// check if char is not in Dna5String
template <typename tt1> bool notNInside(tt1 seq) {
  Iterator<Dna5String>::Type it = seqan::begin(seq);
  Iterator<Dna5String>::Type itEnd = seqan::end(seq);
  for (; it != itEnd; goNext(it)) {
    if (getValue(it) == 'N') {
      return false;
    }
  }
  return true;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// Functions to get random sequence

typedef std::vector<char> char_array;
char_array charset() {
  // Change this to suit
  return char_array({'A', 'T', 'C', 'G'});
};

// given a function that generates a random character,
// return a string of the requested length
std::string random_string(size_t length, std::function<char(void)> rand_char) {
  std::string str(length, 0);
  std::generate_n(str.begin(), length, rand_char);
  return str;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

Dna5String genRandSeq(int lbp) {
  // 0) create the character set.
  //   yes, you can use an array here,
  //   but a function is cleaner and more flexible
  const auto ch_set = charset();

  // 1) create a non-deterministic random number generator
  std::default_random_engine rng(std::random_device{}());
  // 2) create a random number "shaper" that will give
  //   us uniformly distributed indices into the character set
  std::uniform_int_distribution<> dist(0, ch_set.size() - 1);
  // 3) create a function that ties them together, to get:
  //   a non-deterministic uniform distribution from the
  //   character set of your choice.
  auto randchar = [ch_set, &dist, &rng]() { return ch_set[dist(rng)]; };

  Dna5String seqout = random_string(lbp, randchar);
  return seqout;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

#endif /* end protective inclusion */
