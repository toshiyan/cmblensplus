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


def ext_params(declare,argtype='intent(out)',verbose=False):
    # extract output arguments
    pset = []
    for dec in declare:
        if argtype in dec[0]:
            p = dec[1].replace(',',' ').split()
            pset.extend(p)
    if verbose:  print(pset)

    return pset


def ext_charg(slines,verbose):
    # extract arguments replaced for python
    # such argument is specified by chargs
    pset = []
    for l in slines:
        if '!chargs' in l and '::' in l:
            p0, p1 = l[l.find('::')+3:].replace('\n','').split('->')
            p0 = p0.replace(' ','')
            p1 = p1.replace(' ','')
            pset.append([p0,p1])

    if verbose:  print(pset)

    return pset


def ext_rmarg(slines,verbose):
    # extract arguments removed for python
    # such argument is specified by rmargs
    args = []
    for l in slines:
        if '!rmargs' in l and '::' in l:
            rmargs = l[l.find('::')+3:].replace('\n','').replace(' ','').split(',')
            for p in rmargs:
                args.append([p])

    if verbose:  print(args)

    return args


def ext_opt4py(slines,verbose):
    # extract optional arguments for python but not for f90 
    # such argument is specified by opt4py
    pset = []
    for l in slines:
        if '!opt4py ::' in l:
            for oparg in l[l.find('::')+3:].replace('\n','').replace(' ','').split(','):
                p, v = oparg.split('=')
                p = p.replace(' ','')
                v = v.replace(' ','')
                pset.append([p,v])

    if verbose:  print(pset)

    return pset


def func_rm_args(Func,pset):
    """
    Remove or modify arguments in a function definition string based on pset.

    Args:
        Func (str): Function definition string (e.g., "func(a,b,c)").
        pset (list): List of tuples where each tuple contains argument name and optionally its default value.

    Returns:
        str: Modified function definition string.
    """
    
    # decompose elements first and rejoin to avoid confusion e.g. "abc" and "abctype"

    func, arg_string = Func.split('(')
    arg_list = arg_string.rstrip(')').split(',')

    arg_pset = []
    arg_pass = []
    for arg in arg_list:
        arg_exist = False # check whether arg exists in pset
        for p in pset:
            if arg == p[0]: # only for exact match
                arg_exist = True
                # add values to the args string if it has a default value, otherwise nothing to do
                if len(p)==2:
                    arg_pset.append(p[0]+'='+p[1])
        if not arg_exist:
            arg_pass.append(arg)

    # join all arguments
    arg_list_new = ','.join(arg_pass)+','+','.join(arg_pset)+')'

    return func+'('+arg_list_new


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

# restrict return index
def get_returnvals(slines):
    for l in slines:
        if '!return' in l and '::' in l:
            idx = l[l.find('::')+3:][0:1]
    return '['+idx+']'


# read libname and modulename
parser = argparse.ArgumentParser(description='scan f90 files and create docstring files')
parser.add_argument('-libname',default='')
parser.add_argument('-modname','--list',nargs='+',default='')
parser.add_argument('--verbose',action='store_true')
Args = parser.parse_args()

libname = Args.libname
modname = Args.list
verbose = Args.verbose

for mod in modname:
    
    mod = mod.replace('.f90','')
    f90name = mod+'.f90'
    pyname  = mod+'.py'

    # initial statement
    f = open(pyname,'w')
    f.write('from cmblensplus import '+libname+'\n')
    f.write('import numpy\n') # default to import numpy
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

        # extract function
        func = slines[0].replace('\n','').replace('subroutine ','')
        if verbose: print(func)

        # extract arguments and decompose them into argument types
        pout = ext_params(declare,verbose=verbose)
        #popt = ext_params(declare,'optional',verbose=verbose)
        pext = ext_opt4py(slines,verbose)
        charg = ext_charg(slines,verbose)
        rmarg = ext_rmarg(slines,verbose)
        #popt, pops = ext_optional(popt,declare,slines,verbose)

        # remove output args from "func"
        for p in pout:
            func = func.replace(','+p+',',',')
            func = func.replace(','+p+')',')')

        # extract function name
        name = func[:func.find('(')]

        # extract all arguments
        args_org = func[func.find('(')+1:].replace(')','').split(',')

        # original function derived by f2py (should not place after chargs)
        func_org = libname + '.' + mod + '.' + name.lower() + '(' + ','.join(args_org) + ')'

        # change arguments with charg
        args = func_ch_args(args_org,charg)

        # define func used for cmblensplus
        hunc = libname.replace('lib','') + '.'+mod+'.'+name + '(' + ','.join(args) + ')'

        # set value for optional arguments in "func"
        func = name + '(' + ','.join(args) + ')'
        #func = func_rm_args(func,popt)
        func = func_rm_args(func,pext)
        func = func_rm_args(func,rmarg)
        #func = func_rm_args(func,pops,True)
        func = func.replace(',,',',')
        func = func.replace(',)',')')

        # add python function to the "pyname" file
        f = open(pyname,'a+')
        f.write('def '+func+':\n')
        f.write('  """\n')

        # extract comments and write to file
        ext_docstring(f,slines)

        # add example
        f.write('  Usage:\n')
        f.write('    :'+','.join(pout)+' = '+hunc+':\n')
        f.write('  """\n')
        #for p in pops:
        #    f.write('  if '+p[0]+' is None: '+p[0]+'='+p[1]+'\n')
        # additional
        add_code(f,slines)
        f.write('  return '+func_org+'\n')
        f.write('\n')
        f.close()

