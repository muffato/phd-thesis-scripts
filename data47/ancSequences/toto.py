#! /users/ldog/muffato/python -OO

import os
import sys
import shutil

for i in xrange(36495):
	i = i+1
	inputF = "vertFamilies/family.%d.zip" % i
	dirs1 = "test1/" + "/".join(str(i))
	dirs2 = "test2/" + "/".join("%05d" % i)
	os.makedirs(dirs1)
	os.makedirs(dirs2)
	outputF1 = dirs1 + "/family.zip"
	outputF2 = dirs2 + "/family.zip"
	shutil.copy(inputF, outputF1)
	shutil.copy(inputF, outputF2)
	print inputF, outputF1, outputF2

