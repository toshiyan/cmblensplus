#
# This code automatically generate python module containing docstring
# The code first split the lines into each subroutine, and divide arguments into optional/output
# This code also copies the comments starting with !*.
# 
import os
import argparse

# read libname and modulename
parser = argparse.ArgumentParser(description='scan f90 file and create the signature file for f2py')
parser.add_argument('-libname',default='')
parser.add_argument('-modname','--list',nargs='+',default='')
args = parser.parse_args()

libname = args.libname
modname = args.list

for mod in modname:
  mod = mod.replace('.f90','')
  f90name = mod+'.f90'
  pyname  = mod+'.py'

  # initial statement
  f = open(pyname,'w')
  f.write('import '+libname+'\n')
  f.write('\n')
  f.close()

  # read lines
  f = open(f90name)
  lines = f.readlines()
  f.close()

  # count number of subroutines
  nsub = sum((line.count('end subroutine') for line in lines))

  # obtain line number of each subroutine
  ln = os.popen('grep -n "subroutine " '+f90name+' | cut -f1 -d:').read().split("\n")

  # separate lines into sublines and obtain output strings
  for ns in range(nsub):

    slines = lines[int(ln[2*ns])-1:int(ln[2*ns+1])]

    # extract declaration part
    declare = []
    for line in slines:
      if '::' in line and 'intent' in line:
        dec = line.replace('\n','').split('::')
        declare.append(dec)

    # extract function
    func = slines[0].replace('\n','').replace('subroutine ','')

    # extract parameters from function
    args = func[func.find('(')+1:func.find(')')].split(',')

    # check optional or output arguments
    pout = []
    popt = []
    for p in args:
      for dec in declare:
        if p in dec[1]:
          if 'optional' in dec[0]:
            # look for def value
            for l in slines:
              if p in l and '!docstr' in l:
                defval = l[l.find('::')+3:].replace('\n','').split('=')[1]
                break
              if p in l and '!f2py' in l:
                defval = l[l.find('::')+3:].replace('\n','').split('=')[1]
            popt.append([p,defval])
            break
          elif 'intent(out)' in dec[0]:
            pout.append(p)
            break

    #print('optional args:',popt)
    #print('output args:',pout)

    # remove output args from func
    for p in pout:
      func = func.replace(','+p,'')

    # original function
    gunc = libname.replace('lib','')+'.'+mod+'.'+func

    # set value for optional args
    # set 0 at the function as sometimes the default value depends on input
    # decompose elements and rejoin to avoid confusion e.g. "abc" and "abctype"
    for p in popt:
      subf = []
      for f in func.replace(')','').split(','):
        if p[0] in f and len(p[0])==len(f): #replace to p=0, only for exact match
          subf.append(p[0]+'=0')
        else: #store default args
          subf.append(f)
      func = ','.join(subf)+')'

    #print('created func:',func)

    # add to pyname
    f = open(pyname,'a+')
    f.write('def '+func+':\n')
    f.write('  """\n')
    # extract comments
    for line in slines:
      if '!*' in line:
        f.write(line.replace('!*',''))
    f.write('  Usage:\n')
    f.write('    - e.g., '+','.join(pout)+' = '+gunc+'\n')
    f.write('  """\n')
    for p in popt:
      f.write('  if '+p[0]+'==0: '+p[0]+'='+p[1]+'\n')
    f.write('  return '+gunc+'\n')
    f.write('\n')
    f.close()

