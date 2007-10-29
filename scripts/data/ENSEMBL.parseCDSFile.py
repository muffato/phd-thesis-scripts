#! /users/ldog/muffato/python -OO

__doc__ = """
	Reformatte le fichier de CDS de biomart
"""


##################
# INITIALISATION #
##################

# Librairies
import os
import sys
import utils.myTools
import utils.myPhylTree

########
# MAIN #
########

# Arguments
(noms_fichiers, options) = utils.myTools.checkArgs( ["inputFile"], [], __doc__)


f = utils.myTools.myOpenFile(noms_fichiers["inputFile"], "r")
dic = {}
seqName = None
for ligne in f:
	ligne = ligne[:-1]
	if ligne.startswith('>'):
		if seqName != None:
			if (seqName not in dic) or (len(seq) > len(dic[seqName])):
				if "Sequence unavailable" not in seq:
					dic[seqName] = seq
		seqName = ligne[1:]
		seq = ""
	else:
		seq = seq + ligne

f.close()

for (gene,seq) in dic.iteritems():
	print "%s\t%s" % (gene,seq)
