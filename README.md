
# Computes  elastic energy to form a nucleosome, or prediction of nucleosome occupancy/free energy based on sequence periodicity

This is currently a branch of `periodic_elastic`.

For the analysis of multiple short sequences (~200 bp).
Expects fasta file with one or more sequences. Does the following calculations (one or the other depending on options defined in the command line):

## Elastic energy deformation

Based on MD-averaged helical parameters force-constants and averaged values, plus reference values of helical values in "canonical nucleosomes". Requires database files (`stif_bsc1_k_avg_miniabc_dinuc.dat`, `refnuc_bp.dat`, and `stif_bsc1_k_avg_miniabc.dat`), which were sent to me by F. Battistini, from the Orozco group.

## Computes nucleosome occupancy / free energy based on sequence using Van Noort's predictor

Based on "Sequence-based prediction of single nucleosome positioning and genome-wide nucleosome occupancy", van der Heijden et al.  DOI: [10.1073/pnas.1205659109](https://dx.doi.org/10.1073/pnas.1205659109)

The code is basically a 1-to-1 translation of a python script that  was sent to me by J. van Noort, without any mention of license. There are four parameters in the model. Two have been hardcoded as `#defines`, the other two are defaulted to the recommended values but can be changed.

### Code requirements

1. Seqan 2.2.0 (other versions untested, but probably fine)
2. Eigen 3.3.0 (other versions untested, but probably fine)
3. OpenMP (but could be removed with loss of functionality)

### Install 

```
git@gitlab.com:guillemportella/periodic_elastic.git
cd period_elastic
```

Fix the CMakeLists depending on your library location. 

In clust1-headnode-1 this shoudl work

```bash
module load cmake/3.7.1 gcc/6.2.0 BOOST/1.56.0
cmake3 . -DCMAKE_C_COMPILER=/home/portel01/programs/compilers/gcc-6.3/bin/gcc -DCMAKE_CXX_COMPILER=/home/portel01/programs/compilers/gcc-6.3/bin/g++
```

