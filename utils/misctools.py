import os
import numpy as np
from IPython.display import clear_output


def check_path(filename,overwrite=False,verbose=True,output='exist and is not overwritten'):

    if not overwrite and os.path.exists(filename):
        if verbose: print(filename+' '+output)
        skip = True
        return skip


def progress(i,index,text='Current progress',addtext=''):

    clear_output(wait=True)
    print(text+' '+addtext+":",np.round((i-min(index))/len(index)*100,2),"%")



