import os
import numpy as np
import scipy
import scipy.optimize
import string
import itertools
import shutil
import h5py


clean_plaintext = True

def main():

    data_dir_array = ["data"]
    origdir=os.getcwd()
    h5filename = "output.h5"

    for data_dir in data_dir_array:
        print "=================\n",data_dir, "data_dir\n","=================\n"
        os.chdir(data_dir)
        L_array = sorted(map(lambda x: eval(x.strip("L")), filter(lambda x: x.find("L")!=-1,os.listdir(os.getcwd()))),key = lambda x: float(x))
        for L in L_array:
            os.chdir('L'+str(L))
            U_array = sorted(map(lambda x: eval(x.strip("U")), filter(lambda x: x.find("U")!=-1,os.listdir(os.getcwd()))),key = lambda x: float(x))
            for U in U_array:
                os.chdir('U'+str(U))
                T_array = sorted(map(lambda x: eval(x.strip("T")), filter(lambda x: x.find("T")!=-1,os.listdir(os.getcwd()))),key = lambda x: float(x))
                for T in T_array:
                    os.chdir('T'+str(T))
                    beta=1.0/T

                    print "============="
                    print "Cleaning ",
                    print os.path.relpath(os.getcwd(),origdir)

                    if not os.path.exists(h5filename):
                        print "no data - deleting"
                        shutil.rmtree(os.getcwd())
                    else:
                        files=os.listdir(os.getcwd())
                        files_other = filter(lambda x : x.find(".dat")!=-1, np.setdiff1d(files,[h5filename]))
                        print files_other
                        for x in files_other:
                            os.remove(x)
                        
                    os.chdir("..")  #exit T
                os.chdir("..")  #exit U
            os.chdir("..")  #exit L
        os.chdir("..")  #exit data_dir

def find_nearest_index(array,value):
    idx = (np.abs(array-value)).argmin()
    return idx
def find_nearest(array,value):
    return array[find_nearest_index(array,value)]

if __name__ == "__main__":
    main()
