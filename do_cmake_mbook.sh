cp CMakeLists_mbook.txt CMakeLists.txt
 cmake . -DCMAKE_C_COMPILER=/opt/local/bin/clang-mp-3.9 -DCMAKE_CXX_COMPILER=/opt/local/bin/clang++-mp-3.9 -DOpenMP_C_FLAGS=/opt/local/lib/libomp -DOpenMP_CXX_FLAGS=/opt/local/lib/libomp -DSEQAN_INCLUDE_PATH=/Users/guillem/work/soft/seqan-library-2.2.0/include
