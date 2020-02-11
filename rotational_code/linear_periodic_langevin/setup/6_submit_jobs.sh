cd .. ; for i in {k2,k3,k4} ; do cd ${i} ; qsub Run_Langevin.sh ; cd .. ; done
