#
# This code automatically generate python module containing docstring
# The code first split the lines into each subroutine, and divide arguments into optional/output
# This code also copies the comments starting with !*.
# 
import os
import argparse


def extract_declare(slines):
    declare = []
    for line in slines:
        if '::' in line and 'intent' in line:
            dec = line.replace('\n','').split('::')
            declare.append(dec)
    return declare


def ext_params(declare,argtype='intent(out)'):
    # extract output arguments
    pset = []
    for dec in declare:
        if argtype in dec[0]:
            p = dec[1].replace(',',' ').split()
            pset.extend(p)

    if pset: print(argtype+':',pset)

    return pset


def ext_opt4py(slines):
    # extract opt4py arguments
    pset = []
    for l in slines:
        if '!opt4py ::' in l:
            p, v = l[l.find('::')+3:].replace('\n','').split('=')
            p = p.replace(' ','')
            v = v.replace(' ','')
            pset.append([p,v])

    if pset: print('extra optional args:',pset)
    return pset


def ext_optional(pset,declare,slines):
    popt = []
    pops = [] #optional arguments whose default value depends on the input arguments
    for p in pset:
        ops = False
        # look for def value
        for l in slines:
            if p in l and '!docstr' in l:
                defval = l[l.find('::')+3:].replace('\n','').split('=')[1]
                ops = True
                break
            if p in l and '!f2py' in l:
                defval = l[l.find('::')+3:].replace('\n','').split('=')[1]
        #store optional arg + its def val
        if ops: 
            pops.append([p,defval])
        else:
            popt.append([p,defval])

    if popt: print('optional args:',popt)
    if pops: print('output args with special defval:',pops)

    return popt, pops


def funcopt(func,pset,ispops=False):
    # remove optional args
    for p in pset:
        subf = []
        # loop for argments in the function definition
        for f in func.replace(')','').split(','): 
            # decompose elements and rejoin to avoid confusion e.g. "abc" and "abctype"
            if p[0] in f and len(p[0])==len(f): #only for exact match
                if ispops:
                    subf.append(p[0]+'=None')
                else:
                    subf.append(p[0]+'='+p[1])
            else: #for not optional arguments, just pass and store default args
                subf.append(f)
        # join all arguments
        func = ','.join(subf)+')'

    return func


def ext_docstring(f,slines):
    for line in slines:
        if '!*' in line:
            # argument type to italic
            if ':' in line:
                line = line.replace('(int)','(*int*)')
                line = line.replace('(double)','(*double*)')
                line = line.replace('(dcmplx)','(*dcmplx*)')
                line = line.replace('(str)','(*str*)')
                line = line.replace('(bool)','(*bool*)')
                line = line.replace('[','[*')
                line = line.replace(']','*]')
            # remove space between ) and :
            if '*)' in line:
                line = '   '+" ".join(line.split())
                line = line.replace(') :','):')+'\n'
            line = line.replace('!*','')
            f.write(line)


# read libname and modulename
parser = argparse.ArgumentParser(description='scan f90 file and create the signature file for f2py')
parser.add_argument('-libname',default='')
parser.add_argument('-modname','--list',nargs='+',default='')
Args = parser.parse_args()

libname = Args.libname
modname = Args.list

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
      declare = extract_declare(slines)

      # check argument type
      pout = ext_params(declare)
      popt = ext_params(declare,'optional')
      pext = ext_opt4py(slines)
      popt, pops = ext_optional(popt,declare,slines)

      # extract function
      func = slines[0].replace('\n','').replace('subroutine ','')

      # remove output args from "func"
      for p in pout:
          func = func.replace(','+p+',',',')
          func = func.replace(','+p+')',')')

      # original function
      gunc = libname+'.'+mod+'.'+func
      hunc = libname.replace('lib','')+'.'+mod+'.'+func

      # set value for optional arguments in "func"
      func = funcopt(func,popt)
      func = funcopt(func,pext)
      func = funcopt(func,pops,True)
      print('create',mod+'.'+func)

      # add to pyname
      f = open(pyname,'a+')
      f.write('def '+func+':\n')
      f.write('  """\n')

      # extract comments and write to file
      ext_docstring(f,slines)

      # add example
      f.write('  Usage:\n')
      f.write('    :'+','.join(pout)+' = '+hunc+':\n')
      f.write('  """\n')
      for p in pops:
          f.write('  if '+p[0]+' is None: '+p[0]+'='+p[1]+'\n')
      f.write('  return '+gunc+'\n')
      f.write('\n')
      f.close()

