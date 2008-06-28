#! /users/ldog/muffato/python -OO

__doc__ = """
Dessine un genome en coloriant ses genes a partir d'un autre genome reference.
"""

##################
# INITIALISATION #
##################

# Librairies
import sys
import utils.myGenomes
import utils.myTools
import utils.myPsOutput

########
# MAIN #
########

# Arguments
arguments = utils.myTools.checkArgs( [("lstGenomes",file)], [], __doc__)

# Chargement de la liste des genomes
lstGenomes = []
f = utils.myTools.myOpenFile(arguments["lstGenomes"], 'r')
for ligne in f:
	if ligne[0] == '.':
		gen = [int(x) for x in ligne.split()[1:]]
		lstGenomes.append(gen)

nb = len(gen)
espacementX = 19. / (len(lstGenomes))
#espacementX = 1.
espacementY = 28. / nb

# On ecrit le PostScipt
utils.myPsOutput.printPsHeader()


for j in xrange(len(lstGenomes)):
	utils.myPsOutput.drawLine(1.+espacementX*(j+1), 1., 0., espacementY*nb, None)

lastNbPos = 0
lastNbVois = 0
for i in xrange(nb):
	val = lstGenomes[0][i]
	i1 = i
	vois = set()
	if i < (nb-1):
		vois.add(lstGenomes[0][i+1])
	pos = set([i])
	for j in xrange(1,len(lstGenomes)):
		i2 = lstGenomes[j].index(val)
		if i2 < (nb-1):
			vois.add(lstGenomes[j][i2+1])
		pos.add(i2)
		if i1 != i2:
			utils.myPsOutput.drawLine(1.+espacementX*j, 1.+espacementY*i1, espacementX, espacementY*(i2-i1), None)
		i1 = i2
	#utils.myPsOutput.drawLine(1.+0.5*espacementX*(1.+float(lastNbVois)/len(lstGenomes)), 1.+espacementY*i, 0.5*espacementX*float(len(vois)-lastNbVois)/len(lstGenomes), espacementY, None)
	#utils.myPsOutput.drawLine(1.+0.5*espacementX*(1.-float(lastNbPos)/len(lstGenomes)), 1.+espacementY*i, 0.5*espacementX*float(lastNbPos-len(pos))/len(lstGenomes), espacementY, None)
	lastNbVois = len(vois)
	lastNbPos = len(pos)

utils.myPsOutput.printPsFooter()

