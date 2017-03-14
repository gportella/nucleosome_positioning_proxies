
# Compute dinucleotide (AA/TA/TT) periodicity, elastic energy to form a nucleosome, or prediction of nucleosome occupancy/free energy based on sequence periodicity

For the analysis of multiple short sequences (~200 bp).
Expects fasta file with one or more sequences, computes averaged periodicity
of AA/TA/TT, or the minimum elastic energy of all possible nucleosomes
that could be formed.  

## Elastic energy deformation

Based on MD-averaged helical parameters force-constants and averaged values, plus reference values of helical values in "canonical nucleosomes". Requires database files (`stif_bsc1_k_avg_miniabc_dinuc.dat`, `refnuc_bp.dat`, and `stif_bsc1_k_avg_miniabc.dat`), which were sent to me by F. Battistini, from the Orozco group.

## Computes nucleosome occupancy / free energy based on sequence using Van Noort's predictor

Based on "Sequence-based prediction of single nucleosome positioning and genome-wide nucleosome occupancy", van der Heijden et al.  DOI: [10.1073/pnas.1205659109](https://dx.doi.org/10.1073/pnas.1205659109)

The code is basically a 1-to-1 translation of a python script that  was sent to me by J. van Noort, without any mention of license. There are four parameters in the model. Two have been hardcoded as `#defines`, the other two are defaulted to the recommended values but can be changed.

