#include "commandline_parse.hpp"
#include "do_periodic.hpp"
#include "nuc_elastic.hpp"
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

int main(int argc, char const **argv) {

  // OpenMP setting the maximum number of threads possible
  int nProcessors = omp_get_max_threads();
  // for now
  nProcessors = 2;
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

  conditionPeriodic cond = {parseOptions.cutoff, parseOptions.b_elastic,
                            parseOptions.b_verbose};

  if (b_verbose) {
    std::cout << "Set " << nProcessors << " OpenMP threads" << std::endl;
  }

  // check the extension and decide if we go for fasta or fastq
  std::string fn = toCString(sequenceFileName);
  if (fn.substr(fn.find_last_of(".") + 1) == "fq") {
    // we assume right now that the file might be valid
    b_fastq = true;
  }

  checkArguments(cond);
  // open sequence file or randomly genereate one
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

  std::string fc_tetra = "stif_bsc1_k_avg_miniabc.dat";
  std::string fc_dinuc = "stif_bsc1_k_avg_miniabc_dinuc.dat";
  if (!cond.b_elastic) {
    do_all_periodic(dseqs, cond);
  } else {
    NNmodel tetrabp;
    NNmodel dinucp;
    if (loadBPModel(tetra_bpmodel, fc_dinuc) &&
        loadBPModel(dinuc_bpmodel, fc_dinuc)) {
      Eigen::MatrixXf refnuc = loadRefNuc();
      do_all_elastic(tetra_bpmodel, refnuc, dseqs, parseOptions.outFileName);
    } else {
      exit(1);
    }
  }
  // caca

  // caca

  // caca
}
