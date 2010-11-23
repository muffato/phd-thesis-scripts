#!/usr/bin/env python2
# Copyright (c) 2010 IBENS/Dyogen : Matthieu MUFFATO, Alexandra LOUIS, Hugues ROEST CROLLIUS
# Mail : hrc@ens.fr or alouis@biologie.ens.fr
# Licences GLP v3 and CeCILL v2

__doc__ = """
Prend un genome ancestral et la liste des diagonales inferees a partir des especes modernes
Renvoie le pourcentage de qualite des diagonales.
"""


##################
# INITIALISATION #
##################

# Librairies
import sys
import random
import utils.myGenomes
import utils.myTools
import utils.myPhylTree

########
# MAIN #
########

# Arguments
arguments = utils.myTools.checkArgs( [("genesFile",file), ("ancGenesFile",file)], [("geneLossRate",float,10)], "Melange un genome et enleve des genes aleatoirement" )

genome = utils.myGenomes.Genome(arguments["genesFile"])
for c in genome.lstChr:
	for i in xrange(10):
		random.shuffle(genome.lstGenes[c])
	genome.lstGenes[c] = genome.lstGenes[c][:int((1.-arguments["geneLossRate"]/100)*len(genome.lstGenes[c]))]

ancGenes = utils.myGenomes.Genome(arguments["ancGenesFile"])
for g in genome:
	(c,i) = ancGenes.dicGenes[g.names[0]]
	print g.chromosome, " ".join(ancGenes.lstGenes[c][i].names)

