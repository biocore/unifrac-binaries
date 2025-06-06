#
# Parse unifrac_accapi_impl.hpp
# and generate concrete implementations
# for all inline functions ending in _T
#

#
# Arguments:
#  $1 - namespace
#  $2 - method
#

from ssu_generate_helper import print_header,print_body
import sys

variant = sys.argv[1]
method =sys.argv[2]

#
# ==========================
#

with open('unifrac_accapi_impl.hpp','r') as fd:
    lines=fd.readlines()

#
# ==========================
#

# print out the header
print('// Generated from unifrac_accapi_impl.hpp (using method %s)'%method);
print('// Do not edit by hand');
print('');

if method in ('direct','indirect',):
    # we are generating unifrac_accapi.cpp
    print('#include "unifrac_accapi.hpp"');

if method in ('indirect','api',):
    # we referencing the dynamic api
    print('#include "unifrac_accapi_dyn_%s.h"'%variant);

if method in ('direct','api',):
    # we are generating the actual code
    print('#include "unifrac_accapi_impl.hpp"');

nmspace = print_header(variant,method)
print_body(method, lines, nmspace)

