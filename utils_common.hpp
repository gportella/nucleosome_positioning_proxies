
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
// exponentiate, any faster than loops? hard to believe
// from http://stackoverflow.com/a/40642375
struct op {
  double operator()(double d) const { return std::exp(d); }
};
std::vector<double> exp_vect(const std::vector<double> &x) {
  const int n = x.size();
  std::vector<double> y(n);
  std::transform(x.begin(), x.end(), y.begin(), op());
  return y;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
//  Add window/2 zeros on either side of vector f
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
// Adds two vectors, first checks if the are of the same length
template <typename T>
std::vector<T> safe_add_vector(std::vector<T> const &v1,
                               std::vector<T> const &v2) {

  std::vector<T> result(v1.size());

  // better safe than sorry
  if (v1.size() == v2.size()) {
    // one could use boost's zip_iterator as well.
    for (unsigned i = 0; i < v1.size(); ++i) {
      result[i] = v1[i] + v2[i];
    }

  } else {
    std::cerr << "Trying to add two vectors of different length." << std::endl;
    std::cerr << "v1: " << v1.size() << " vs "
              << "v2: " << v2.size() << std::endl;
    exit(1);
  }
  return result;
}
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// Adds two vectors, ASSUMES they are both of the same size. Use only if
// you are completely sure that both have the same length or expect undefined
// behaviour
template <typename T>
std::vector<T> brave_add_vector(std::vector<T> const &v1,
                                std::vector<T> const &v2) {

  std::vector<T> result(v1.size());

  for (unsigned i = 0; i < v1.size(); ++i) {
    result[i] = v1[i] + v2[i];
  }

  return result;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// linspace  numpy-like
template <typename T> std::vector<T> linspace(T a, T b, size_t N) {
  T h = (b - a) / static_cast<T>(N - 1);
  std::vector<T> xs(N);
  typename std::vector<T>::iterator x;
  T val;
  for (x = xs.begin(), val = a; x != xs.end(); ++x, val += h)
    *x = val;
  return xs;
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
// valid returns std::max(f.size(), g.size()) - std::min(f.size(), g.size()) + 1
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
///////////////////////////////////////////////////////////////////////////////
// full returns size f.size() + g.size() - 1
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
///////////////////////////////////////////////////////////////////////////////
// same returns size std::max(f.size(), g.size() )
template <typename T>
std::vector<T> conv_same(std::vector<T> const &f, std::vector<T> const &g) {

  int const n = (f.size() > g.size()) ? f.size() : g.size();
  std::vector<T> out(n, T());
  for (unsigned i = 0; i < out.size(); ++i) {
    int low_k = std::max(0, (int)i - (int)((g.size() - 1.0) / 2.0));
    int high_k = std::min((int)f.size() - 1, (int)i + (int)(g.size() / 2.0));
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
// usese a rectangular box to average
////////////////////////////////////////////////////////////////////////////////
template <typename tt1, typename tt2>
std::vector<double> smooth_box(tt1 x, tt2 w_len) {

  // image the extremes of x : create mirrored copies of the ends
  // and paste together
  std::vector<double> b(w_len - 1);
  for (unsigned i = 0; i < w_len - 1; ++i) {
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

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// linear interpolation
// addapted from https://people.sc.fsu.edu/~jburkardt/cpp_src/interp/interp.html
// IMPORTANT: does not take into account cases in which the x_range of both
// x_ and x1 is different. If this is the case, you should change the code
// and e.g. set everything to the left/rigth of x_ to x[0]/x[x.size()-1]
///////////////////////////////////////////////////////////////////////////////
template <typename tt1, typename tt2, typename tt3>
void nearest_bracket0(tt1 x, tt2 xval, tt3 &left, tt3 &right) {
  for (unsigned i = 2; i <= x.size() - 1; ++i) {
    if (xval < x[i - 1]) {
      left = i - 2;
      right = i - 1;
      return;
    }
  }
  left = x.size() - 2;
  right = x.size() - 1;
  return;
}

template <typename T>
std::vector<T> interp_linear(std::vector<T> const &x_, std::vector<T> const &x1,
                             std::vector<T> const &y1) {
  int left;
  int right;
  std::vector<T> result(x_.size());

  for (unsigned i = 0; i < x_.size(); ++i) {
    double t = x_[i];
    //  Find the interval [ x1(LEFT), x1(RIGHT) ] that contains, or is
    //  nearest to t
    nearest_bracket0(x1, t, left, right);
    result[i] = ((x1[right] - t) * y1[left] + (t - x1[left]) * y1[right]) /
                (x1[right] - x1[left]);
  }
  return result;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// write out x y
template <typename tt1, typename tt2> void write_xy(tt1 fout, tt2 x, tt2 y) {
  std::ofstream file;
  file.open(toCString(fout));
  for (unsigned i = 0; i < x.size(); ++i) {
    file << x[i] << " " << y[i] << std::endl;
  }
  file.close();
}

#endif /* end protective inclusion */
