#!/usr/bin/python

# this allows to run any python script as a program
#
# read the common definitions for mceep and greate a corresponding
# file to be used with root. This is the used to create a root tree
# from a mceep ntuple

import sys
import string
import re              # regular expressions


def var_ind(l):    # return the index of the variable
    p1 = re.split('\(',l)  # first split at (
    p2 = re.split('\)',p1[1])
    return int(p2[0])
#
# end of def


f=open('var_dat.cmn','r')  # var_dat.cmn is the common file
                            # used in mceep to define the variable names


lines = f.readlines()       # read the hole file and store it in lines
var=[]     # empty variable list

imax = 0
for l in lines:             # loop over each line
    if (re.search('DATA VAR_NAME',l)): # found a line of data
        

        i = var_ind(l)                 # get the variable index

        if (i > imax): imax = i        # ind largest index
        v =  re.split('/',l)           # variable name
        v[1] = re.sub("\'","\"",v[1])
        var.append ( [i,v[1]])            # add to list
    #
#
# sort the array var
var.sort()

# prepare the root output

w=file("m2root.h","w")

w.write( "// variable name definitions as used in mceep\n")
w.write( "// defined in var_dat.cmn include file\n")
w.write(" // global array of variable names \n")
w.write("TString var[" + str(imax+1) + "]; \n") # in fortran variables star at 1
w.write("\n")
w.write("Int_t init_var()\n")
w.write( "{\n")
                                                # keep it for comparison
for v in var:
    w.write( "var["+str(v[0])+ "].Append(" + str(v[1]) + ");\n")
#
w.write("return(0);\n")
w.write( "}\n")

    



