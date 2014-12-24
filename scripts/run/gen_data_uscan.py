import os
import numpy as np
import subprocess
import time
import sys
import shutil
import itertools
#from cluster_tools import *
import argparse
import collections
from pbs_pks import *

data_dir="data"
#L_array = [8, 16, 24] # array of system sizes
L_array = [48]#8, 16, 24] # array of system sizes
#L_array = [32]
dry_run = False # if true - will not make actual calculation
queue = "mid10" # queue for execution
python_path = which("python")
nprocs = 24
chebyshev = True
resume = True

T_array = [0.15]
U_array = np.arange(0,6.01, 0.25)


print "System sizes :",L_array
print "T values     :",T_array
print "U values     :",U_array

def main():
    counter = 0
    origdir = os.getcwd()
    for L in L_array:
      for T in T_array:
        for U in U_array:
            print "T=",T,"U=",U,"L=",L
            counter = counter + 1
            name = "T"+str(T)+"L"+str(L)+"U"+str(U) 
            
            call_args = [ python_path, 
                 os.getcwd()+os.path.sep+"mc_run.py",
                 "--T",str(T),
                 "--U",str(U),
                 "--L",str(L),
                 "--nprocs",str(nprocs),
                 "--data_dir", str(data_dir), 
                 "--chebyshev", str(int(chebyshev)), 
                 "--resume", str(int(resume))
               ]
            submit_mpi(
                nprocs = nprocs,
                commands = call_args,
                add_mpirun = False,
                name=name,
                dry_run=dry_run,
                output_stream_file = "stdout_"+name,
                error_stream_file = "stderr_"+name,
                cpu_time = "336:00:00",
                file_size = "24G",
                ram_size = "12G",
                prefix = [["source ~/.bash_profile"], ["module load triqs fk_mc mpich2-ch3-gnu hdf5"], ['export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/usr/lib64"']]
                )     
    
    print counter, "calcs submitted"

if __name__ == "__main__":
    main()
