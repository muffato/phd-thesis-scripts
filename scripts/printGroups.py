#!/usr/bin/env python2
# Copyright (c) 2010 IBENS/Dyogen : Matthieu MUFFATO, Alexandra LOUIS, Hugues ROEST CROLLIUS
# Mail : hrc@ens.fr or alouis@biologie.ens.fr
# Licences GLP v3 and CeCILL v2

__doc__ = """
	Lit des donnees dans un fichier et les regroupe
"""

import utils.myFile
import utils.myTools


arguments = utils.myTools.checkArgs([("files",utils.myTools.FileList(0))], [], __doc__, showArgs=False)

# Les fichiers a regrouper (par defaut, on lit l'entree standard)
files = arguments["files"]
if len(files) == 0:
	files.append("-")

comb = utils.myTools.myCombinator()
for f in files:
	# Ouverture du fichier
	f = utils.myFile.openFile(f, 'r')
	# Lecture & regroupement
	for l in f:
		comb.addLink(l.split())
	# Fermeture
	f.close()

# On affiche le resultat
for x in comb:
	print " ".join(x)


