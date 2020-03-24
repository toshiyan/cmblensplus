import os

def check_path(filename,overwrite=False,verbose=True,output='exist and is not overwritten'):

    if not overwrite and os.path.exists(filename):
        if verbose: print(filename+' '+output)
        skip = True
        return skip


