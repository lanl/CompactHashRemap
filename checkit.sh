#!/bin/sh

set -v

AMR_remap/AMR_remap_openMP 128 1 4.74      0  1 -adapt-meshgen -no-brute
AMR_remap/AMR_remap_openMP 128 2 3.97      0  1 -adapt-meshgen -no-brute
AMR_remap/AMR_remap_openMP 128 3 3.11      0  1 -adapt-meshgen -no-brute
AMR_remap/AMR_remap_openMP 128 4 2.30      0  1 -adapt-meshgen -no-brute
AMR_remap/AMR_remap_openMP 128 5 1.59      0  1 -adapt-meshgen -no-brute
AMR_remap/AMR_remap_openMP 128 6 0.97      0  1 -adapt-meshgen -no-brute
AMR_remap/AMR_remap_openMP 128 7 0.55      0  1 -adapt-meshgen -no-brute
AMR_remap/AMR_remap_openMP 128 8 0.30      0  1 -adapt-meshgen -no-brute

AMR_remap/AMR_remap_openMP 200 1 3.01      0  1 -adapt-meshgen -no-brute
AMR_remap/AMR_remap_openMP 200 2 2.54      0  1 -adapt-meshgen -no-brute
AMR_remap/AMR_remap_openMP 200 3 2.13      0  1 -adapt-meshgen -no-brute
AMR_remap/AMR_remap_openMP 200 4 1.75      0  1 -adapt-meshgen -no-brute
AMR_remap/AMR_remap_openMP 200 5 1.28      0  1 -adapt-meshgen -no-brute
AMR_remap/AMR_remap_openMP 200 6 0.85      0  1 -adapt-meshgen -no-brute
AMR_remap/AMR_remap_openMP 200 7 0.51      0  1 -adapt-meshgen -no-brute
AMR_remap/AMR_remap_openMP 200 8 0.28      0  1 -adapt-meshgen -no-brute

AMR_remap/AMR_remap_openMP 256 1 2.25      0  1 -adapt-meshgen -no-brute
AMR_remap/AMR_remap_openMP 256 2 2.02      0  1 -adapt-meshgen -no-brute
AMR_remap/AMR_remap_openMP 256 3 1.75      0  1 -adapt-meshgen -no-brute
AMR_remap/AMR_remap_openMP 256 4 1.46      0  1 -adapt-meshgen -no-brute
AMR_remap/AMR_remap_openMP 256 5 1.13      0  1 -adapt-meshgen -no-brute
AMR_remap/AMR_remap_openMP 256 6 0.78      0  1 -adapt-meshgen -no-brute
AMR_remap/AMR_remap_openMP 256 7 0.49      0  1 -adapt-meshgen -no-brute
#../../AMR_remap/AMR_remap_openMP 255 8 0.28 474345  1 -adapt-meshgen -no-brute

