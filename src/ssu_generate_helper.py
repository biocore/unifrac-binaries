#
# Helper function that parse a header file
# and generate concrete implementations
# for all inline functions ending in _T
#

#
# ==========================
#

# replace template name with concrete type
def patch_type(s,t):
    return s.replace('TFloat',t).replace('TNum',t)

# extract argument name from the arg definition
def get_arg_name(s):
    return s.split()[-1].split('*')[-1]

# extract type from the arg definition
def get_arg_type(s):
    arg_name = get_arg_name(s)
    return s[:-len(arg_name)]

#
# ==========================
#

def print_func_direct_noargs(ftype,nmspace,fname):
    print('%s %s::%s() {'%(ftype,nmspace,fname))
    if ftype=='void':
        print('  %s_T();'%fname);
    else:
        print('  return %s_T();'%fname)
    print('}');

def print_func_direct_args(ftype,nmspace,fname,ttype,fargs):
    print('template<>')
    print('%s %s::%s('%(ftype,nmspace,fname))
    for el in fargs[:-1]:
       print('\t\t\t%s,'%patch_type(el,ttype))
    print('\t\t\t%s) {'%patch_type(fargs[-1],ttype))

    print('  %s_T('%fname)
    for el in fargs[:-1]:
       print('\t%s,'%get_arg_name(el))
    print('\t%s);'%get_arg_name(fargs[-1]))

    print('}');

#
# ==========================
#

def print_func_indirect_noargs(ftype,nmspace,fname):
    print('static %s (*dl_%s_%s)() = NULL;'%(ftype,nmspace,fname))
    print('%s %s::%s() {'%(ftype,nmspace,fname))
    if (ftype=='bool') and (fname.find('found_gpu')>=0):
        #this is the initial check, allow for the shared library to not being avaialble
        print('  if (!ssu_load_check()) return false; /* shlib not found */')
    print('  cond_ssu_load("%s_%s", (void **) &dl_%s_%s);'%(nmspace,fname,nmspace,fname))
    if ftype=='void':
        print('  (*dl_%s_%s)();'%(nmspace,fname));
    else:
        print('  return (*dl_%s_%s)();'%(nmspace,fname))
    print('}');

def print_func_indirect_args(ftype,nmspace,fname,ttype,fargs):
    print('static void (*dl_%s_%s_%s)('%(nmspace,fname,ttype))
    for el in fargs[:-1]:
       print('\t\t\t%s,'%get_arg_type(patch_type(el,ttype)))
    print('\t\t\t%s) = NULL;'%get_arg_type(patch_type(fargs[-1],ttype)))
    print('template<>')
    print('%s %s::%s('%(ftype,nmspace,fname))
    for el in fargs[:-1]:
       print('\t\t\t%s,'%patch_type(el,ttype))
    print('\t\t\t%s) {'%patch_type(fargs[-1],ttype))

    print('  cond_ssu_load("%s_%s_%s", (void **) &dl_%s_%s_%s);'%(nmspace,fname,ttype,nmspace,fname,ttype))
    print('  (*dl_%s_%s_%s)('%(nmspace,fname,ttype))
    for el in fargs[:-1]:
       print('\t%s,'%get_arg_name(el))
    print('\t%s);'%get_arg_name(fargs[-1]))

    print('}');

#
# ==========================
#

def print_func_api_h_noargs(ftype,nmspace,fname):
    print('extern "C" %s %s_%s();'%(ftype,nmspace,fname))

def print_func_api_h_args(ftype,nmspace,fname,ttype,fargs):
    print('extern "C" %s %s_%s_%s('%(ftype,nmspace,fname,ttype))
    for el in fargs[:-1]:
       print('\t\t\t%s,'%patch_type(el,ttype))
    print('\t\t\t%s);'%patch_type(fargs[-1],ttype))

#
# ==========================
#

def print_func_api_noargs(ftype,nmspace,fname):
    print('%s %s_%s() {'%(ftype,nmspace,fname))
    if ftype=='void':
        print('  %s_T();'%fname);
    else:
        print('  return %s_T();'%fname)
    print('}');

def print_func_api_args(ftype,nmspace,fname,ttype,fargs):
    print('%s %s_%s_%s('%(ftype,nmspace,fname,ttype))
    for el in fargs[:-1]:
       print('\t\t\t%s,'%patch_type(el,ttype))
    print('\t\t\t%s) {'%patch_type(fargs[-1],ttype))

    print('  %s_T('%fname)
    for el in fargs[:-1]:
       print('\t%s,'%get_arg_name(el))
    print('\t%s);'%get_arg_name(fargs[-1]))

    print('}');

#
# ==========================
#

def print_func_noargs(method,ftype,nmspace,fname):
    if method=='direct':
        print_func_direct_noargs(ftype,nmspace,fname)
    elif method=='indirect':
        print_func_indirect_noargs(ftype,nmspace,fname)
    elif method=='api':
        print_func_api_noargs(ftype,nmspace,fname)
    elif method=='api_h':
        print_func_api_h_noargs(ftype,nmspace,fname)
    else:
        raise "Unknown generation method '%s'"%method

def print_func_args(method,ftype,nmspace,fname,ttype,fargs):
    if method=='direct':
        print_func_direct_args(ftype,nmspace,fname,ttype,fargs)
    elif method=='api':
        print_func_api_args(ftype,nmspace,fname,ttype,fargs)
    elif method=='indirect':
        print_func_indirect_args(ftype,nmspace,fname,ttype,fargs)
    elif method=='api_h':
        print_func_api_h_args(ftype,nmspace,fname,ttype,fargs)
    else:
        raise "Unknown generation method '%s'"%method

#
# ==========================
#
# Print out the common header of the
# concrete implementation
# Also returns
#  namespace
#
# Args:
#  variant
#  method - direct, indirect, api, api_h
def print_header(variant,method):
    if method in ('api_h',):
        # bool and unit_t are not standard in C without these header
        print('#include <stdbool.h>')
        print('#include <stdint.h>')

    if method in ('indirect',):
        # function expected by ssu_ld
        print('static const char *ssu_get_lib_name() { return "libssu_%s.so";}'%variant)
        print('#include "ssu_ld.c"')

    print('')

    return "su_%s"%variant

#
# ==========================
#
# Print out the main body of the
# concrete implementation
#
# Args:
#  method - direct, indirect, api, api_h
#  lines - content of the header file
#  nmspace - namespace to use for the concrete implementations

def print_body(method,lines,nmspace):
    # now we can generate the concrete functions
    i=0
    while (i<len(lines)):
        line = lines[i]
        i+=1
        if not line.startswith('static inline '):
            continue # not the beginning of an intersting function
        if line.find('_T(')<0:
            continue # not the beginning of an intersting function

        larr=line.split();
        ftype = larr[2];
        fname = larr[3].split('_T(')[0]

        print('// ==================================');
        if line.find('_T()')>=0:
            # special case, no arguments
            print_func_noargs(method,ftype,nmspace,fname)
            print('');
            continue

        # assuming ftype == void simplifies the code
        if ftype!='void':
            raise "Unsupported templated void found!"

        ftypes = set()
        fargs = []
        line = lines[i]
        i+=1
        while True:
            # get all the arguments up to the optional close
            larr=line.split(') ')[0].strip().split(',')
            for el in larr:
                el = el.strip() 
                if len(el)>0:
                    fargs.append(el)
                    if el.find('TFloat')>=0:
                        ftypes |= {'float','double'}
                    elif el.find('TNum')>=0:
                        ftypes |= {'float','double','uint64_t','uint32_t','bool'}


            if line.find(') ')>=0:
                break # found end of args, exit the loop
            line = lines[i]
            i+=1

        # make sure we create always the same output, sort
        sorted_ftypes = sorted(ftypes)
        for ft in sorted_ftypes:
            print_func_args(method,ftype,nmspace,fname,ft,fargs)
            print('');


