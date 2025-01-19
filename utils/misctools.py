import os
import numpy as np
import sys
import configparser

#from IPython.display import clear_output

def check_path_core(filename,overwrite,verbose,leave,output):

    skip = False

    if not overwrite and os.path.exists(filename):
        if verbose: 
            if leave:
                print(filename+' '+output)
            else:
                print(filename+' '+output,end="\r")
        skip = True

    return skip
 


def check_path(filename,overwrite=False,verbose=True,output='exist and is not overwritten',leave=False):

    if isinstance(filename,str):
        skip = check_path_core(filename,overwrite,verbose,leave,output)

    if isinstance(filename,list):
        skips = [ check_path_core(fname,overwrite,verbose,leave,output) for fname in filename ]
        skip = all(skips)
   
    if isinstance(filename,dict):
        skips = [ check_path_core(fname,overwrite,verbose,leave,output) for key, fname in filename.items() ]
        skip = all(skips)
 
    return skip


def progress(i,index,text='Current progress',addtext=''):

    #clear_output(wait=True)
    print(text+' '+addtext+":",np.round((i-min(index))/len(index)*100,2),"%")



def load_config(section):
        
    #//// load config file ////#
    config = configparser.ConfigParser()
    
    if np.size(sys.argv) > 1 and '.ini' in sys.argv[1]:
        print('reading '+sys.argv[1])
        config.read(sys.argv[1])
    else:
        config.add_section(section)

    #//// get parameters ////#
    return config[section]



def create_directory(directory_path,verbose=True):
    try:
        # Create the directory if it doesn't exist
        os.makedirs(directory_path)
        print(f"Directory '{directory_path}' created successfully.")
    except OSError as e:
        if verbose:
            print(f"Error: {directory_path} - {e}")



