#include <iostream>
#include <seqan/arg_parse.h>
#include <seqan/seq_io.h>

using namespace seqan;
// using namespace std;

// for the conditions
typedef struct my_conditions {
  unsigned int cutoff;
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

  // I guess this is how to initialize
  Options()
      : b_verbose(false), b_elastic(false), b_vnoort(false), b_random(false),
        b_nuccore(false) {}
};

ArgumentParser::ParseResult parseCommandLine(Options &parseOptions, int argc,
                                             char const **argv) {
  // Setup ArgumentParser.
  ArgumentParser parser("period_elastic");
  setShortDescription(
      parser, "Finds dinucleotide periodicity | nucleosome elastic energy.");
  addDescription(parser, "Finds dinucleotide periodicity "
                         "for AA/TA/AT, or compute nucleosome elastic energy");
  addUsageLine(parser, "[\\fIOPTIONS\\fP] ");
  setVersion(parser, "0.3");
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
      ArgParseOption("nr", "num_rand",
                     "How many random sequences of 147 base pairs to generate",
                     seqan::ArgParseArgument::INTEGER, "INTEGER"));
  setDefaultValue(parser, "result_file", "out_analysis.txt");
  setValidValues(parser, "sequence_file", "fna fa fq fna.gz fasta fa.gz fq.gz");
  setDefaultValue(parser, "cutoff", "50");
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
  parseOptions.b_verbose = isSet(parser, "be_verbose");
  parseOptions.b_nuccore = isSet(parser, "only_nuccore");
  parseOptions.b_elastic = isSet(parser, "compute_elastic");
  parseOptions.b_vnoort = isSet(parser, "do_vannoort");
  parseOptions.b_random = isSet(parser, "random_seq");

  return ArgumentParser::PARSE_OK;
}
