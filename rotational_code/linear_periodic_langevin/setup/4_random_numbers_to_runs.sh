cd random_input ; vim random_input.f90 ; gfortran mt95.f90 random_input.f90 -o random_input ; ./random_input ; cd .. ; for i in {k2,k3,k4}; do for sep in {3563,11768,33102,55690}; do for eps in {eps_B_lt_A,eps_B_gt_A}; do for run in {1..50} ; do head -n 1 random_input/random.txt > add1.txt ; cat add1.txt ${i}/sep_${sep}/${eps}/Run_${run}/input.txt > add2.txt ; rm add1.txt ; mv add2.txt ${i}/sep_${sep}/${eps}/Run_${run}/input.txt ; sed -i "1d" random_input/random.txt ; done ; done ; done ; done
