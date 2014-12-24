import os
import numpy as np
import scipy
import scipy.optimize
import string
import itertools
import h5py

import matplotlib
import matplotlib.pyplot as plt 
from matplotlib.font_manager import FontProperties
def main():

    data_obs = ["dos_err", "ipr_err"]

    font = {'family' : 'FreeSans',
        'weight' : 'light',
        }

    #print matplotlib.font_manager.findSystemFonts(fontpaths=None, fontext='ttf')

    plt.rc('font', **font) 

    data_dir_array = ["data.old"]
    outdir_array = ["data_stats"]
    origdir=os.getcwd()
    h5filename = "output.h5"

    f = plt.figure(1, (24., 12.))

    L_array = [32,48]
    U_array = [1.0, 4.0, 6.0]
    T_array = [0.15][::-1]

    shape1 = (len(T_array), len(U_array))

    from mpl_toolkits.axes_grid1 import Grid
    grid = Grid(f, 111, # similar to subplot(132)
        nrows_ncols = shape1,
        axes_pad = 0.02,
        #aspect=True,
        share_all=True,
        add_all=True,
        label_mode = "all", #"L"
        direction = "column",
    )
    
    for (data_dir,out_dir) in itertools.izip(data_dir_array,outdir_array):
        out_dir = os.path.abspath(out_dir)
        print "=================\n",data_dir, "data_dir\n","=================\n"
        print "Saving output to", out_dir
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
        os.chdir(data_dir)
        figcount = 0
        for U in U_array:
            for T in T_array:
                pls=[]
                leg=[]
                for L in L_array:
                    print "============="
                    print "Statistics for",
                    dirname=os.path.abspath(origdir+os.path.sep+"{0}/L{1}/U{2}/T{3}".format(data_dir,L,U,T))
                    os.chdir(dirname)
                    print os.path.relpath(os.getcwd(),origdir)


                    if os.path.exists(h5filename):
                        #try:
                            h5data = h5py.File(h5filename, "r")
                            for obs in data_obs:
                                obs_data = h5data["stats"][obs]
                                pls=pls+[grid[figcount].errorbar(obs_data[:,0],obs_data[:,1],obs_data[:,2], 
                                                        elinewidth=2.0, capsize=2.0, errorevery=10)]
                                grid[figcount].set_xlim([-6.,6.])
                                grid[figcount].set_ylim([0.,1.])
                                grid[figcount].text(0.0,0.9,"U = {0}, T = {1}, 7k steps".format(U,T),ha="center", family = font["family"], weight="light", size="large")
                                grid[figcount].set_xlabel("w")
                                leg = leg + [obs.split("_err")[0]+" L = {0}".format(L)]
                            h5data.close()
                        #except:
                        #    print "Couldn't get data"
                    os.chdir(origdir)
                grid[figcount].legend(pls,leg,"upper right")
                figcount = figcount + 1
        os.chdir("..")  #exit data_dir
        plt.savefig(out_dir+os.path.sep+"ipr.pdf",bbox_inches='tight', dpi=100)

def find_nearest_index(array,value):
    idx = (np.abs(array-value)).argmin()
    return idx
def find_nearest(array,value):
    return array[find_nearest_index(array,value)]

if __name__ == "__main__":
    main()
