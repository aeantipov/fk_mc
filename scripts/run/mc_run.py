import os
import numpy as np
import subprocess
from pbs_pks import *

def heaviside(x):
    return 0.5 * (np.sign(x) + 1)

def run(T_array, U_array,L,data_dir,dry_run,nprocs, chebyshev, resume):

    ncycles = 2**7 # number of measurements
    nwarmup = 2**0 # number of warmup steps
    cycle_len = 2**6 # how many mc steps in a measurement

    print T_array

    random_seed = True # set completely random seed
    save_plaintext = False # save data to plaintext files
    overwrite_dir = False #True # overwrite old files (useful for recalculation)
    calc_history = True # True # this gets errorbars of dos 
    calc_ipr = True # True # this gets inverse participation ratio
    cheb_prefactor = 2.5
    origdir = os.getcwd()

    print "U=",U_array
    print "T=",T_array
    print "L=",L
    print "nprocs =",nprocs

    for T in T_array:
      beta=1.0/T
      for U in U_array:
        print
        print "============="
        print "Making calc in",
        dirname=os.path.abspath(origdir+os.path.sep+"{0}/L{1}/U{2}/T{3}".format(data_dir,L,U,T))
        print dirname
        dir_exists = os.path.exists(dirname)
        make_calc = True
        if not dir_exists:
            os.makedirs(dirname)
            make_calc = True
        else:
            make_calc = not os.path.exists(dirname+os.path.sep+"output.h5") or overwrite_dir or resume
            
        print "Making calc:",make_calc

        if make_calc:
            os.chdir(dirname)

            if not os.path.exists("output.h5"):
                resume = False

            args = [
            "--U", str(float(U)),
            "--T", str(float(T)),
            "--mu", str(float(U)/2.0),
            "--L", str(L),
            "--seed" if random_seed else "",
            "-p" if save_plaintext else "",
            "--calc_history", str(int(calc_history)),
            "--calc_ipr", str(int(calc_ipr)),
            "--ncycles", str(ncycles),
            "--cyclelen", str(cycle_len),
            "--nwarmup", str(nwarmup),
            "--dos_width", str(max(6.0, 1.3*U)),
            "--chebyshev" if chebyshev else "",
            "--resume" if resume else "",
            "--cheb_prefactor", str(cheb_prefactor)
            ]

            print os.getenv("QPREFIX")
            call_args = [
            which("mpirun"), "--np", str(nprocs), which("fk_mc_cubic2d")] + args
            print call_args
            print ' '.join(call_args)
            if not dry_run:
                subprocess.call(' '.join(call_args),shell=True)
            else:
                print "skipping run"
            os.chdir(origdir)

        #    run_fkmc.submit(args)
        print "\n"
        #    counter = counter + 1
 

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='FK U loop')
    parser.add_argument('--U', '-U', nargs="+", help='values of T', type=float)
    parser.add_argument('--T', '-T', nargs="+", help='values of T', type=float)
    parser.add_argument('--data_dir', help='data directory', default = "output")
    parser.add_argument('--L', help='L',type=int, default = 12)
    parser.add_argument('--nprocs', help='n procs',type=int, default = 8)
    parser.add_argument('--dry_run', help='dry run', action="store_true", dest="dry_run", default=False)
    parser.add_argument('--chebyshev', help='chebyshev', type=int, default = 0)
    parser.add_argument('--resume', help='resume', type=int, default = 1)
    args = parser.parse_args()

    run(U_array = args.U, T_array = args.T, L = args.L, data_dir = args.data_dir, dry_run = args.dry_run, nprocs = args.nprocs, chebyshev = args.chebyshev, resume = args.resume)
