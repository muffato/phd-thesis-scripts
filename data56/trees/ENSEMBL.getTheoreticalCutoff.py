#!/usr/bin/env python2
# Copyright (c) 2010 IBENS/Dyogen : Matthieu MUFFATO, Alexandra LOUIS, Hugues ROEST CROLLIUS
# Mail : hrc@ens.fr or alouis@biologie.ens.fr
# Licences GLP v3 and CeCILL v2

__doc__ = """
	Renvoie la liste des seuils theoriques en fonction du nombre d'especes a haute couverture
"""

import utils.myFile
import utils.myTools
import utils.myPhylTree

arguments = utils.myTools.checkArgs( [("phylTree.conf",file)], [("minHighCoverageGenomesProportion",float,0.5)], __doc__)

phylTree = utils.myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])

# Limites automatiques de score de duplication
for anc in phylTree.listAncestr:
	nesp = len(phylTree.species[anc])
	n2X = len(phylTree.lstEsp2X.intersection(phylTree.species[anc]))
	# La moitie des especes non 2X a vu la duplication (au minimum 1 espece)
	print utils.myFile.myTSV.printLine([anc, round(max(1., arguments["minHighCoverageGenomesProportion"]*(nesp-n2X)) / nesp, 3) - 2e-3])

