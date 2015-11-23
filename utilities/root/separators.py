#!/usr/bin/python

# this allows to run any python script as a program
#

# example: return a list with the values contained between 2 different
# separators inb a string such as (), [] or {} or any other

import sys
import string
import re              # regular expressions

line = "this is ((1) the (((th[er]e) is (2) [and] there is (text))"

print "original line : ", line

l_sep = "\("   # left separator
r_sep = "\)"   # right separator


p1 = re.split(l_sep, line)

print "after the first split : ", p1

a=[]
for i in range(1,len(p1)):
    p2 = re.split(r_sep,p1[i])
    if (p2[0]): a.append(p2[0])  # use only nonblnk values
#

print "final result : ", a


