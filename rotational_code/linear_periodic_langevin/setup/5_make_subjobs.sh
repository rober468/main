cd subjobs_list ; python subjobs_list.py ; cd ../.. ; mv setup/k* . ; for i in {k2,k3,k4} ; do cp setup/subjobs_list/subjobs.txt setup/subjobs_list/Run_Langevin.sh ${i}/ ; done
