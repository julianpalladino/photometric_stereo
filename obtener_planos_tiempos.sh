#!/bin/bash

for size in `seq 50 50 1000`
do
echo > experimentacion/tiemposObtenerPlanoNormal/fact_lu_$size.txt
echo 'Running lu with' $size 
for run in {1..100}
do
  FACT_LU=1 MATRIX_SIZE=$size .././run_tests 2> /dev/null >> experimentacion/tiemposObtenerPlanoNormal/fact_lu_$size.txt
done
done
