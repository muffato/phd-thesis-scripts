#! /users/ldog/muffato/python -OO

import sys
import time
import zipfile
import utils.myTools

(noms_fichiers, options) = utils.myTools.checkArgs( ["InputSeq"], [("IN.familyZipFile",str,"")], "!")

currID = None
txt = ""

def write():
	fz =  zipfile.ZipFile(options["IN.familyZipFile"] % currID, "a", zipfile.ZIP_DEFLATED)
	zinfo = zipfile.ZipInfo(filename="ancSequence.txt", date_time=time.localtime(time.time()))
	zinfo.compress_type =  zipfile.ZIP_DEFLATED
	zinfo.external_attr = 2175008768
	fz.writestr(zinfo, txt)
	fz.close()
	print currID

f = utils.myTools.myOpenFile(noms_fichiers["InputSeq"], "r")
for l in f:
	if l.startswith("ID"):
		newID = int(l[:-1].split()[1])
		#currF = backup
		if newID != currID:
			if currID != None:
				write()
			currID = newID
			txt = ""
			#backup = currF = utils.myTools.myOpenFile(options["outputSeq"] % currID, "w")
		# On ne veut pas la ligne "ID XXXX"
		l = ""
	#elif l.startswith("SPECIES"):
	#	x = l[:-1].split('\t')[1].split()
	#	if len(x) == 1:
	#		currF = utils.myTools.null
	#	l = ("ID\t%d\n" % currID) + l
	#elif l.startswith("PROBA"):
	#	lst = [" %g" % float(x) for x in l[:-1].split()[1:]]
	#	l = l.split()[0] + ''.join(lst) + '\n'
	#if '\t' not in l:
	#	l = l.replace(' ', '\t', 1)
	#currF.write(l)
	#print >> currF, l,
	txt += l
write()
f.close()

