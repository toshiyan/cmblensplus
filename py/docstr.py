#
# This code automatically generate python module containing docstring
# The code first split the lines into each subroutine, and divide arguments into optional/output
# This code also copies the comments starting with !*.
# 
import os
import argparse


def extract_declare(slines):
    # extract declaration part of f90 code
    declare = []
    for line in slines:
        # extract lines which contain both "::" and "intent"
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

    return pset


def ext_charg(slines):
    # extract arguments replaced for python
    # such argument is specified by charg
    pset = []
    for l in slines:
        if '!chargs' in l and '::' in l:
            p0, p1 = l[l.find('::')+3:].replace('\n','').split('->')
            p0 = p0.replace(' ','')
            p1 = p1.replace(' ','')
            pset.append([p0,p1])

    return pset


def ext_opt4py(slines):
    # extract optional arguments for python but not for f90 
    # such argument is specified by opt4py
    pset = []
    for l in slines:
        if '!opt4py ::' in l:
            p, v = l[l.find('::')+3:].replace('\n','').split('=')
            p = p.replace(' ','')
            v = v.replace(' ','')
            pset.append([p,v])

    return pset


def ext_optional(pset,declare,slines):  # !!! This function should be removed in the future
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

    return popt, pops


def func_rm_args(func,pset,ispops=False):
    # remove args contained in pset
    # decompose elements first and rejoin to avoid confusion e.g. "abc" and "abctype"
    for p in pset:
        subf = []
        # loop for argments in the function definition
        for f in func.replace(')','').split(','): 
            if p[0] in f and len(p[0])==len(f): #only for exact match
            #if p[0] == f: #only for exact match
                if ispops:
                    subf.append(p[0]+'=None')
                else:
                    subf.append(p[0]+'='+p[1])
            else: #for not optional arguments, just pass and store default args
                subf.append(f)
        # join all arguments
        func = ','.join(subf)+')'

    return func


def func_ch_args(args,pset):
    # remove args contained in pset
    # first seperate func name and args
    chargs = args.copy()
    for p in pset:
        # loop for argments in the function definition
        if p[0] in args: #only for exact match
            chargs.insert(chargs.index(p[0]),p[1])
            chargs.remove(p[0])

    return chargs


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


# add code to function
def add_code(f,slines):
    for l in slines:
        if '!add2py' in l and '::' in l:
            code = l[l.find('::')+3:]
            f.write('  '+code)


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
        charg = ext_charg(slines)
        popt, pops = ext_optional(popt,declare,slines)

        # extract function
        func = slines[0].replace('\n','').replace('subroutine ','')

        # remove output args from "func"
        for p in pout:
            func = func.replace(','+p+',',',')
            func = func.replace(','+p+')',')')

        # call lib (should not place after chargs)
        libfunc = libname+'.'+mod+'.'+func

        # change arguments
        name = func[:func.find('(')]
        args = func[func.find('(')+1:].replace(')','').split(',')
        args = func_ch_args(args,charg)
        func = name + '(' + ','.join(args) + ')'

        # func for example usage
        hunc = libname.replace('lib','')+'.'+mod+'.'+func

        # set value for optional arguments in "func"
        func = func_rm_args(func,popt)
        func = func_rm_args(func,pext)
        func = func_rm_args(func,pops,True)
        #print('create',mod+'.'+func)

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
        # additional
        add_code(f,slines)
        f.write('  return '+libfunc+'\n')
        f.write('\n')
        f.close()

