#! /bin/bash
#
g++ -c -Wall -fopenmp fft_openmp.cpp
if [ $? -ne 0 ]; then
  echo "Compile error"
  exit
fi
#
g++ -fopenmp fft_openmp.o -lm
if [ $? -ne 0 ]; then
  echo "Load error"
  exit
fi
mv a.out fft_openmp
#
#  Run with 1, 2, and 4 threads.
#
echo "Run with 1 thread."
export OMP_NUM_THREADS=1
./fft_openmp > fft_openmp.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
#
echo "Run with 2 threads."
export OMP_NUM_THREADS=2
./fft_openmp >> fft_openmp.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
#
echo "Run with 4 threads."
export OMP_NUM_THREADS=4
./fft_openmp >> fft_openmp.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm fft_openmp
#
echo "Normal end of execution."
