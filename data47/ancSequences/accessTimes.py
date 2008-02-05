#! /users/ldog/muffato/python -OO

import os
import sys
import time
import random

listid = [random.randint(1,36495) for _ in xrange(int(sys.argv[1]))]

def fread(name):
	t1 = time.clock()
	s1 = time.time()
	f = open(name, "r")
	x = sum([len(l) for l in f])
	f.close()
	t2 = time.clock()
	s2 = time.time()
	return (x,t2-t1,s2-s1)

ti = 0
to1 = 0
to2 = 0
to3 = 0

si = 0
so1 = 0
so2 = 0
so3 = 0

#for i in listid:
for i in xrange(36495):
	i += 1
	#i = random.randint(1,36495)
	
	inputF = "vertFamilies/family.%d.zip" % i
	outputF1 = "test1/" + "/".join(str(i)) + "/family.zip"
	outputF2 = "test2/" + "/".join("%05d" % i) + "/family.zip"
	outputF3 = "test3/%d/family.zip" % i
	
	#shutil.copy(inputF, outputF1)

	(xi,ri,rri) = fread(inputF)
	(xo1,ro1,rro1) = fread(outputF1)
	(xo2,ro2,rro2) = fread(outputF2)
	(xo3,ro3,rro3) = fread(outputF3)
	
	if (xi != xo1) or (xi != xo2) or (xi != xo3):
		print >> sys.stderr, "PB", i
	
	ti += ri
	to1 += ro1
	to2 += ro2
	to3 += ro3

	si += rri
	so1 += rro1
	so2 += rro2
	so3 += rro3
	#print inputF, outputF1, outputF2

print ti, to1, to2, to3
print si, so1, so2, so3


