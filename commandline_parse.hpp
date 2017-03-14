
#ifndef CL_PARSER_H
#define CL_PARSER_H
#include <iostream>
#include <seqan/arg_parse.h>
#include <seqan/seq_io.h>

using namespace seqan;
// using namespace std;

// for the conditions
typedef struct my_conditions {
  unsigned int cutoff;
  unsigned int vn_window;
  double vn_mu;
  unsigned int vn_smooth_window;
  bool b_elastic;
  bool b_vnoort;
  bool b_verbose;
  bool b_nuccore;
} conditionPeriodic;

struct Options {
  bool b_verbose;
  bool b_elastic;
  bool b_vnoort;
  bool b_random;
  bool b_nuccore;
  CharString needlesFileName;
  CharString outFileName;
  unsigned int cutoff;
  unsigned int num_rand;
  unsigned int vn_window;
  double vn_mu;
  unsigned int vn_smooth_window;

  // I guess this is how to initialize
  Options()
      : b_verbose(false), b_elastic(false), b_vnoort(false), b_random(false),
        b_nuccore(false) {}
};

ArgumentParser::ParseResult parseCommandLine(Options &parseOptions, int argc,
                                             char const **argv) {
  // Setup ArgumentParser.
  ArgumentParser parser("period_elastic");
  setShortDescription(parser, "Finds dinucleotide periodicity | nucleosome "
                              "elastic energy or occupancy.");
  addDescription(parser,
                 "Finds dinucleotide periodicity "
                 "for AA/TA/AT, compute nucleosome elastic energy or"
                 "compute prediction of nucleosome occupancy/FE based on"
                 "van Noort's method.");
  addUsageLine(parser, "[\\fIOPTIONS\\fP] ");
  setVersion(parser, "0.9");
  setDate(parser, "March 2017");

  // Define Options
  addOption(parser, ArgParseOption("i", "sequence_file", "A fasta input file "
                                                         "with sequence to "
                                                         "search for.",
                                   ArgParseArgument::INPUT_FILE));
  addOption(parser, ArgParseOption("o", "result_file", "Output histogram. ",
                                   ArgParseArgument::OUTPUT_FILE));

  addOption(parser,
            ArgParseOption("c", "cutoff", "How far to search for matches .",
                           seqan::ArgParseArgument::INTEGER, "INTEGER"));
  addOption(parser, ArgParseOption("v", "be_verbose", "Be verbose in "
                                                      "what you do."));
  addOption(parser, ArgParseOption("elastic", "compute_elastic",
                                   "Compute the minimum elastic "
                                   "energy for each sequence."));
  addOption(parser,
            ArgParseOption("nuccore", "only_nuccore",
                           "Use a window of 74 bp centered  "
                           "at the dyad to compute nucleosome energy."));
  addOption(parser,
            ArgParseOption("vnoort", "do_vannoort",
                           "Use van Noort's sequence-based prediction "
                           "to compute nucleosome occupancy and free energy."));
  addOption(parser, ArgParseOption("random", "random_seq",
                                   "Generate random sequences "
                                   "instead of reading input."));
  addOption(
      parser,
      ArgParseOption("vn_window", "vannoort_window",
                     "The base window used to compute van Noort's prediction. "
                     "It should be either 74 or 147. The paper suggests 74.",
                     seqan::ArgParseArgument::INTEGER, "INTEGER"));
  addOption(parser,
            ArgParseOption("vn_smooth_window", "vannoort_smooth_window",
                           "The number or bases to use as smoothing window"
                           "The paper suggests 10.",
                           seqan::ArgParseArgument::INTEGER, "INTEGER"));
  addOption(
      parser,
      ArgParseOption(
          "vn_mu", "vannoort_mu",
          "Affinity between the histone and the DNA in van Noort's model. "
          "Units are in kT. The paper suggests -1.5 for in-vitro.",
          seqan::ArgParseArgument::DOUBLE, "DOUBLE"));
  addOption(
      parser,
      ArgParseOption("nr", "num_rand",
                     "How many random sequences of 147 base pairs to generate",
                     seqan::ArgParseArgument::INTEGER, "INTEGER"));
  setDefaultValue(parser, "result_file", "out_analysis.txt");
  setValidValues(parser, "sequence_file", "fna fa fq fna.gz fasta fa.gz fq.gz");
  setDefaultValue(parser, "cutoff", "50");
  setDefaultValue(parser, "vannoort_window", "74");
  setDefaultValue(parser, "vannoort_mu", "-1.5");
  setDefaultValue(parser, "vannoort_smooth_window", "10");
  setDefaultValue(parser, "num_rand", "1000");
  setValidValues(parser, "result_file", "txt");

  // Parse command line.
  ArgumentParser::ParseResult res = parse(parser, argc, argv);

  // Only extract options if the
  // program continues after
  // parseCommandLine()
  if (res != seqan::ArgumentParser::PARSE_OK)
    return res;

  // Extract option values.
  getOptionValue(parseOptions.needlesFileName, parser, "sequence_file");
  getOptionValue(parseOptions.outFileName, parser, "result_file");
  getOptionValue(parseOptions.cutoff, parser, "cutoff");
  getOptionValue(parseOptions.num_rand, parser, "num_rand");
  getOptionValue(parseOptions.vn_window, parser, "vannoort_window");
  getOptionValue(parseOptions.vn_mu, parser, "vannoort_mu");
  getOptionValue(parseOptions.vn_smooth_window, parser,
                 "vannoort_smooth_window");
  parseOptions.b_verbose = isSet(parser, "be_verbose");
  parseOptions.b_nuccore = isSet(parser, "only_nuccore");
  parseOptions.b_elastic = isSet(parser, "compute_elastic");
  parseOptions.b_vnoort = isSet(parser, "do_vannoort");
  parseOptions.b_random = isSet(parser, "random_seq");

  // especific for this program: check if vn_window is either 74 or 147
  // setValidValues does not let me check against numeric values, so I do it
  // here by hand

  if (parseOptions.vn_window != 74 && parseOptions.vn_window != 147) {
    std::cerr << "Please select either 74 or 147 as the VanNoort window"
              << std::endl;
    exit(1);
  }
  if (parseOptions.vn_smooth_window > 100) {
    std::cerr << "Please do not use a smoothing window larger than 100 bp."
              << std::endl;
    exit(1);
  }

  return ArgumentParser::PARSE_OK;
}
#endif /* end of protective declarations */
