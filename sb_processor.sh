for (( i = 1; i < 6; ++i )); do
  a=$(( i*i*i ))
  mpirun -np $a python processorMF.py $i
done