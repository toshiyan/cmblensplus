# 
# This code scan f90 files and output the signature file.
#
# The code first split each subroutine, and scan lines including intent.
# The argument should be declared as 
#   [data type], [intent(...)], [dimension] (<- order insensitive) :: [arguments]
# For optional arguments, the code also look for a default value specified by !f2py ...
#
import os
import argparse
import re

# read libname and modulename
parser = argparse.ArgumentParser(description='scan f90 file and create the signature file for f2py')
parser.add_argument('-libname',default='')
parser.add_argument('-modname','--list',nargs='+',default='')
args = parser.parse_args()

libname = args.libname
modname = args.list
pyfname = libname+'.pyf'

# initial import
f = open(pyfname,'w')
f.write('python module '+libname+'\n')
f.write('    interface\n')

# loop over modules
for mod in modname:

  # read lines
  f90name = mod
  g = open(f90name)
  lines = g.readlines()
  g.close()

  # add module name
  mod = mod.replace('.f90','')
  f.write('        module '+mod+'\n')

  # add use statement
  for line in lines:
    if 'use ' in line and not '!' in line:
      f.write('          '+line)

  # count number of subroutines
  nsub = sum((line.count('end subroutine') for line in lines))

  # obtain line number of each subroutine
  ln = os.popen('grep -n "subroutine " '+f90name+' | cut -f1 -d:').read().split("\n")

  # separate lines into sublines and obtain output strings
  for ns in range(nsub):

    # lines for this subroutine
    slines = lines[int(ln[2*ns])-1:int(ln[2*ns+1])]

    # extract declaration part
    declare = []
    for line in slines:
      if '::' in line and 'intent' in line:
        dec = line.replace('\n','').split('::')
        declare.append(dec)

    # extract f2py optional part
    opt = []
    for line in slines:
      if '::' in line and '!f2py' in line:
        op = line.replace('\n','').split('::')
        opt.append(op)
        print(op)


    # extract fuction
    func = slines[0].replace('\n','')
    f.write('            '+func+'\n')

    # extract args from the function definition
    args = func[func.find('(')+1:func.find(')')].split(',')

    # extract args description from the lines
    for p in args:
      for dec in declare:
          print(dec[1])

          # split to avoid confusion of e.g. "abc" and "abctype"
          #d = re.sub("[\(].*?[\)]", "", dec[1].replace(',','')).split()
          d = dec[1].replace(',','').split()
          if p in d:

            # for optional args
            defval = ''
            if 'optional' in dec[0]:
              for op in opt:
                if p in op[1].split():
                    defval = op[1].split('=')[1]

            output = dec[0]

            # check dependence
            for q in args:
              if q in dec[0] or q in defval:
                output += ', depend('+q+')'

            output += ' :: ' + p

            if defval!='':
              output += ' =' + defval

            # write to file
            f.write('              '+output+'\n')
            break

    f.write('            '+slines[-1].replace('\n','')+'\n')

  f.write('        end module '+mod+'\n')

f.write('    end interface\n')
f.write('end python module '+libname+'\n')
f.write('\n')
f.close()

