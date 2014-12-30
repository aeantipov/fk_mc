import os
import numpy as np
import scipy
import scipy.optimize
import string
import itertools
import h5py

def main():

    obs_binned = ["cv", "binder_pi", "fsusc_pi"]

    data_dir_array = ["pd"]
    origdir=os.getcwd()
    h5filename = "output.h5"

    for data_dir in data_dir_array:
        out_dir = data_dir + "_stats"
        out_dir = os.path.abspath(out_dir)
        print "=================\n",data_dir, "data_dir\n","=================\n"
        print "Saving output to", out_dir
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
        os.chdir(data_dir)
        L_array = sorted(map(lambda x: eval(x.strip("L")), filter(lambda x: x.find("L")!=-1,os.listdir(os.getcwd()))),key = lambda x: float(x))
        for L in L_array:
            os.chdir('L'+str(L))
            U_array = sorted(map(lambda x: eval(x.strip("U")), filter(lambda x: x.find("U")!=-1,os.listdir(os.getcwd()))),key = lambda x: float(x))
            U_array = filter(lambda x:x>0, U_array)
            for U in U_array:
                os.chdir('U'+str(U))
                streams = dict()
                for obs in obs_binned:
                    streams[obs] = open(out_dir+os.path.sep+obs+"_L"+str(L)+"U"+str(U)+".dat","w")
                T_array = sorted(map(lambda x: eval(x.strip("T")), filter(lambda x: x.find("T")!=-1,os.listdir(os.getcwd()))),key = lambda x: float(x))
                for T in T_array:
                    os.chdir('T'+str(T))
                    beta=1.0/T

                    print "============="
                    print "Statistics for",
                    print os.path.relpath(os.getcwd(),origdir)

                    if os.path.exists(h5filename):
                        h5data = h5py.File(h5filename, "r")
                        for obs in obs_binned:
                            try:
                                obs_data = h5data["stats"][obs]
                                (nbins, value, disp, error) = obs_data
                                print obs,":", value,"+/-",error
                                streams[obs].write(str(T)+"   "+str(value)+" "+str(error)+"\n")
                            except:
                                print "Couldn't get", obs
                        h5data.close()
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
