#include <iostream>
#include <seqan/arg_parse.h>
#include <seqan/seq_io.h>

using namespace seqan;
// using namespace std;

// for the conditions
typedef struct my_conditions {
  unsigned int cutoff;
  bool b_elastic;
  bool b_verbose;
} conditionPeriodic;

struct Options {
  bool b_verbose;
  bool b_elastic;
  CharString needlesFileName;
  CharString outFileName;
  unsigned int cutoff;

  // I guess this is how to initialize
  Options() : b_verbose(false), b_elastic(false) {}
};

ArgumentParser::ParseResult parseCommandLine(Options &parseOptions, int argc,
                                             char const **argv) {
  // Setup ArgumentParser.
  ArgumentParser parser("slider");
  setShortDescription(
      parser, "Finds dinucleotide periodicity | nucleosome elastic energy.");
  addDescription(parser, "Finds dinucleotide periodicity "
                         "for AA/TA/AT, or compute nucleosome elastic energy");
  addUsageLine(parser, "[\\fIOPTIONS\\fP] ");
  setVersion(parser, "0.1");
  setDate(parser, "February 2017");

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
                                   "Compute the maximum elastic "
                                   "energy for each sequence."));
  setDefaultValue(parser, "result_file", "out_analysis.txt");
  setValidValues(parser, "sequence_file", "fna fa fq fna.gz fasta fa.gz fq.gz");
  setDefaultValue(parser, "cutoff", "50");
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
  parseOptions.b_verbose = isSet(parser, "be_verbose");
  parseOptions.b_elastic = isSet(parser, "compute_elastic");

  return ArgumentParser::PARSE_OK;
}
