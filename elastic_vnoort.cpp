
#ifndef PERIOD_ELASTIC_MAIN
#define PERIOD_ELASTIC_MAIN
#include "commandline_parse.hpp"
#include "nuc_elastic.hpp"
#include "nuc_vannoort.hpp"
#include "utils_common.hpp"
#include <Eigen/Dense>
#include <functional>
#include <iostream>
#include <math.h>
#include <omp.h>
#include <seqan/arg_parse.h>
#include <seqan/file.h>
#include <seqan/find.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>
#include <stdlib.h>
#include <string>
#include <vector>

int main(int argc, char const **argv) {

  // OpenMP setting the maximum number of threads possible
  int nProcessors = omp_get_max_threads();
  omp_set_num_threads(nProcessors);

  // Parse the command line
  Options parseOptions;
  ArgumentParser::ParseResult res = parseCommandLine(parseOptions, argc, argv);

  // If parsing did not work, then exit with code 1
  // Otherwise exit with code 0
  if (res != ArgumentParser::PARSE_OK)
    return res == ArgumentParser::PARSE_ERROR;

  CharString sequenceFileName = parseOptions.needlesFileName;
  bool b_verbose = parseOptions.b_verbose;
  // declarations for fasta inputs
  StringSet<CharString> ids;
  CharString id;
  StringSet<Dna5String> dseqs;
  Dna5String dseq;
  StringSet<CharString> quals;
  SeqFileIn seqFileIn;
  bool b_fastq = false;
  std::map<std::string, NNmodel> tetra_bpmodel;
  std::map<std::string, NNmodel> dinuc_bpmodel;

  conditionPeriodic cond = {
      parseOptions.vn_window,      parseOptions.vn_mu,
      parseOptions.smooth_window,  parseOptions.b_elastic,
      parseOptions.b_elastic_prof, parseOptions.b_vnoort,
      parseOptions.b_verbose,      parseOptions.b_nuccore};

  if (b_verbose) {
    std::cout << "Set " << nProcessors << " OpenMP threads" << std::endl;
  }

  // check the extension and decide if we go for fasta or fastq
  std::string fn = toCString(sequenceFileName);
  if (fn.substr(fn.find_last_of(".") + 1) == "fq") {
    // we assume right now that the file might be valid
    b_fastq = true;
  }

  // randomly genereate or read them
  if (parseOptions.b_random) {
    for (unsigned i = 0; i < parseOptions.num_rand; ++i) {
      appendValue(dseqs, genRandSeq(NUC_LEN));
    }
  } else {
    if (!open(seqFileIn, toCString(sequenceFileName))) {
      std::cerr << "ERROR: Cound not open input file (sequences).\n";
      return 1;
    }
    try {
      if (b_fastq) {
        readRecords(ids, dseqs, quals, seqFileIn);
      } else {
        readRecords(ids, dseqs, seqFileIn);
      }
    } catch (Exception const &e) {
      std::cout << "ERROR: " << e.what() << std::endl;
      return 1;
    }
  }

  std::string fc_tetra = "stif_bsc1_k_avg_miniabc.dat";
  std::string fc_dinuc = "stif_bsc1_k_avg_miniabc_dinuc.dat";
  if (cond.b_vnoort && !cond.b_elastic && !cond.b_elastic_prof) {
    do_all_vannoort(dseqs, cond, parseOptions.outFileName);
  } else if (cond.b_elastic || cond.b_elastic_prof) {
    NNmodel tetrabp;
    NNmodel dinucp;
    if (loadBPModel(tetra_bpmodel, fc_tetra) &&
        loadBPModel(dinuc_bpmodel, fc_dinuc)) {
      Eigen::MatrixXf refnuc = loadRefNuc();
      if (cond.b_elastic) {
        do_all_elastic(tetra_bpmodel, dinuc_bpmodel, refnuc, dseqs,
                       parseOptions.outFileName, cond);
      } else if (cond.b_elastic_prof) {
        do_all_elastic(tetra_bpmodel, dinuc_bpmodel, refnuc, dseqs,
                       parseOptions.outFileName, cond, Tag_ElProf());
      } else {
        std::cerr << "Should not be here!" << std::endl;
        exit(1);
      }
    } else {
      exit(1);
    }
  } else {
    std::cerr << "You either need to set -elastic, -el_profile or -vnoort"
              << std::endl;
    exit(1);
  }
  // caca

  // caca

  // caca
}

#endif /* end protective inclusion */
