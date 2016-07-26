#!/bin/sh

set -v

cat header.dat >> rundata128.dat
AMR_remap/AMR_remap_openMP 128 1 4.74   18696 100 -adapt-meshgen -no-brute -plot-file
AMR_remap/AMR_remap_openMP 128 2 3.97   34017 100 -adapt-meshgen -no-brute -plot-file
AMR_remap/AMR_remap_openMP 128 3 3.11   92775 100 -adapt-meshgen -no-brute -plot-file
AMR_remap/AMR_remap_openMP 128 4 2.30  285847 100 -adapt-meshgen -no-brute -plot-file
AMR_remap/AMR_remap_openMP 128 5 1.59  789306 100 -adapt-meshgen -no-brute -plot-file
AMR_remap/AMR_remap_openMP 128 6 0.97 2071760 100 -adapt-meshgen -no-brute -plot-file
AMR_remap/AMR_remap_openMP 128 7 0.55      0  100 -adapt-meshgen -no-brute -plot-file
AMR_remap/AMR_remap_openMP 128 8 0.30      0  100 -adapt-meshgen -no-brute -plot-file

cat header.dat >> rundata200.dat
AMR_remap/AMR_remap_openMP 200 1 3.01   43574 100 -adapt-meshgen -no-brute -plot-file
AMR_remap/AMR_remap_openMP 200 2 2.54   68209 100 -adapt-meshgen -no-brute -plot-file
AMR_remap/AMR_remap_openMP 200 3 2.13  169937 100 -adapt-meshgen -no-brute -plot-file
AMR_remap/AMR_remap_openMP 200 4 1.75  543673 100 -adapt-meshgen -no-brute -plot-file
AMR_remap/AMR_remap_openMP 200 5 1.28 1626167 100 -adapt-meshgen -no-brute -plot-file
AMR_remap/AMR_remap_openMP 200 6 0.85      0  100 -adapt-meshgen -no-brute -plot-file
AMR_remap/AMR_remap_openMP 200 7 0.51      0  100 -adapt-meshgen -no-brute -plot-file
AMR_remap/AMR_remap_openMP 200 8 0.28      0  100 -adapt-meshgen -no-brute -plot-file

cat header.dat >> rundata256.dat
AMR_remap/AMR_remap_openMP 256 1 2.25   69947 100 -adapt-meshgen -no-brute -plot-file
AMR_remap/AMR_remap_openMP 256 2 2.02  103182 100 -adapt-meshgen -no-brute -plot-file
AMR_remap/AMR_remap_openMP 256 3 1.75  243663 100 -adapt-meshgen -no-brute -plot-file
AMR_remap/AMR_remap_openMP 256 4 1.46  768627 100 -adapt-meshgen -no-brute -plot-file
AMR_remap/AMR_remap_openMP 256 5 1.13 2395157 100 -adapt-meshgen -no-brute -plot-file
AMR_remap/AMR_remap_openMP 256 6 0.78      0  100 -adapt-meshgen -no-brute -plot-file
AMR_remap/AMR_remap_openMP 256 7 0.49      0  100 -adapt-meshgen -no-brute -plot-file
#../../AMR_remap/AMR_remap_openMP 255 8 0.28 474345  100 -adapt-meshgen -no-brute -plot-file

