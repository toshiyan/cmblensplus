# 
# This code add f2py directives to f90 files.
#
import os
import argparse
import re
import numpy as np

# read libname and modulename
parser = argparse.ArgumentParser(description='modify f90 for f2py')
parser.add_argument('-srcname','--list',nargs='+',default='')
args = parser.parse_args()

# loop over modules
for srcname in args.list:

    # initial import
    mod = srcname.replace('.src','.f90')
    if mod==srcname:
        raise SystemExit("cannot overwrite input file")

    f = open(mod,'w')

    # open a file and read all lines
    g = open(srcname)
    lines = g.readlines()
    g.close()

    # count number of subroutines
    nsub = sum((line.count('end subroutine') for line in lines))

    # obtain line number of each subroutine
    ln = os.popen('grep -n "subroutine " '+srcname+' | cut -f1 -d:').read().strip().split("\n")

    # add header
    for line in lines[:int(ln[0])-1]:
        f.writelines(line)

    # separate lines into sublines and obtain output strings
    for ns in range(nsub):

        # lines for this subroutine
        slines = lines[int(ln[2*ns])-1:int(ln[2*ns+1])]

        modified_slines = []
        modified_slines.append(slines[0]) # subroutine name
 
        # extract I/O arguments
        args_in, args_out, args_p, args_arr = [], [], [], []
        for line in slines:
            if '::' in line:
                desc, args = line.replace('\n','').replace(' ','').split('::') # divide into description and arguments parts
                args = args.split(',')
                if 'intent(in)' not in desc and 'intent(out)' not in desc: continue
                for arg in args:
                    if 'intent(in)' in desc:
                        args_in.append(['intent(in)',arg])
                    elif 'intent(out)' in desc:
                        args_out.append(['intent(out)',arg])
                    if 'dimension(' in desc:
                        args_arr.append([desc,arg])
                    else:
                        args_p.append([desc,arg])

        args_dep = []
        # check dependence
        for desc, arg in args_arr:
            l0 = desc.find('dimension(')
            for __, p in args_in:
                matches = re.findall(r'(?<!\w)' + re.escape(p) + r'(?!\w)', desc[l0:])
                if p in matches:
                    args_dep.append(['depend('+p+')',arg])

        # Group the second items based on the first item
        args_dict = {}
        for dep, item in args_in+args_out+args_dep:
            if dep not in args_dict:
                args_dict[dep] = []
            args_dict[dep].append(item)

        #//// modify slines ////#
        # first add original lines until implicit none
        #indice = int(np.where(np.char.find(slines[1:], 'implicit none') >= 0)[0])
        indice = next((i for i, line in enumerate(slines[1:], 1) if 'implicit none' in line), None)
        for line in slines[1:indice+1]:
            modified_slines.append(line)

        # add !f2py
        for desc, args in args_dict.items():
            modified_slines.append("  !f2py " + desc + ' ' + ', '.join(args) + "\n")

        # add original lines after implicit none
        for line in slines[indice+1:]:
            modified_slines.append(line)

        # Save the modified file
        f.writelines(modified_slines)
        f.writelines("\n")

    # add end module
    for line in lines[int(ln[-1]):]:
        f.writelines(line)

    f.write('\n')
    f.close()

