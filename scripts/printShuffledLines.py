#!/usr/bin/env python2
# Copyright (c) 2010 IBENS/Dyogen : Matthieu MUFFATO, Alexandra LOUIS, Hugues ROEST CROLLIUS
# Mail : hrc@ens.fr or alouis@biologie.ens.fr
# Licences GLP v3 and CeCILL v2

__doc__ = """
	Lit un fichier et melange les lignes
"""

import random

import utils.myFile
import utils.myTools


arguments = utils.myTools.checkArgs([("fichier",file)], [("nbShuffle",int,10)], __doc__)

lst = utils.myFile.openFile(arguments["fichier"], 'r').readlines()
for _ in xrange(arguments["nbShuffle"]):
	random.shuffle(lst)

for x in lst:
	print x,

