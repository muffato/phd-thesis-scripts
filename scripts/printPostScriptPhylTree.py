#!/usr/bin/env python2
# Copyright (c) 2010 IBENS/Dyogen : Matthieu MUFFATO, Alexandra LOUIS, Hugues ROEST CROLLIUS
# Mail : hrc@ens.fr or alouis@biologie.ens.fr
# Licences GLP v3 and CeCILL v2

# Librairies
import sys
import math

import utils.myFile
import utils.myMaths
import utils.myTools
import utils.myPhylTree
import utils.myPsOutput

# Arguments
arguments = utils.myTools.checkArgs( \
	[("phylTree.conf",file)], \
	[("landscape",bool,False), ("printSpecies",bool,True), ("printAncestors",bool,True), ("printAges",bool,False), \
	("colorFile",str,""), ("func",str,""), ("root",str,""), ("min",float,None), ("max",float,None)], \
	"Dessine l'arbre phylogenetique avec des informations dessus comme des evolutions de taux" \
)

phylTree = utils.myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])

(largeur,hauteur) = utils.myPsOutput.printPsHeader(landscape=arguments["landscape"])

if arguments["root"] in phylTree.items:
	root = arguments["root"]
else:
	root = phylTree.root
maxAge = phylTree.ages[root]


if arguments["func"] == "":
	func = lambda x: x
else:
	func = eval(arguments["func"])

colors = {}
if arguments["colorFile"] != "":
	f = utils.myFile.openFile(arguments["colorFile"], "r")
	for l in f:
		t = l.replace('\n','').split("\t")
		s1 = intern(t[0])
		s2 = intern(t[1])
		if not phylTree.isChildOf(s1, root):
			continue
		try:
			x = func(float(t[2]))
			colors[(s1,s2)] = x
			colors[(s2,s1)] = x
		except Exception, e:
			print >> sys.stderr, '"%s" on %s/%s (%s)' % (e, s1, s2, t[2])
	f.close()
	val = colors.values()
	minV = float(min(val))
	maxV = float(max(val))

if arguments["min"] != None:
	minV = arguments["min"]
if arguments["max"] != None:
	maxV = arguments["max"]


refcolors = [(0,0,127), (0,192,192), (0,192,0), (255,255,0), (242,148,0), (224,0,0)]
inter = utils.myMaths.myInterpolator.getMultDim(utils.myMaths.myInterpolator.oneDimCubic, range(len(refcolors)), refcolors)

# Le degrade 
def getColor(value):
	value = (value-minV) / (maxV-minV)
	return inter(value * len(refcolors))

	if value < 0.1:
		# bleu fonce -> bleu
		return utils.myPsOutput.alphaColor( (0,0,127), (0,0,255), 10.*value)
	elif value < 0.3:
		# bleu -> cyan
		return utils.myPsOutput.alphaColor( (0,0,255), (0,255,255), 5.*(value-.1))
	elif value < 0.5:
		# cyan -> vert
		return utils.myPsOutput.alphaColor( (0,255,255), (0,255,0), 5.*(value-.3))
	elif value < 0.7:
		# vert -> jaune
		return utils.myPsOutput.alphaColor( (0,255,0), (255,255,0), 5.*(value-.5))
	elif value < 0.9:
		# jaune -> rouge
		return utils.myPsOutput.alphaColor( (255,255,0), (255,0,0), 5.*(value-.7))
	else:
		# rouge -> rouge fonce
		return utils.myPsOutput.alphaColor( (255,0,0), (127,0,0), 10.*(value-.9))


y = 0
dy = hauteur / (len(phylTree.species[root])+1.)

margeX1 = dy
if arguments["printSpecies"]:
	margeX1 += max([len(x) for x in phylTree.species[root]]) * 0.15

margeX2 = dy
if arguments["printAncestors"]:
	margeX2 += len(root) * 0.15

dx = (largeur-margeX1-margeX2) / maxAge


def printSubTree(node):
	global y
	
	if node not in phylTree.items:
		y += dy
		if arguments["printSpecies"]:
			utils.myPsOutput.drawText(dy, y, node, "black")
		return y
	
	x = margeX1 + phylTree.ages[node]*dx
	mi = hauteur
	ma = 0
	todo = []
	for (e,a) in phylTree.items[node]:
		tmpY = printSubTree(e)
		if tmpY > ma:
			ma = tmpY
		if tmpY < mi:
			mi = tmpY

		c = None
		if (node,e) in colors:
			c = getColor( colors[(node,e)] )
		elif (e,node) in colors:
			c = getColor( colors[(e,node)] )
		todo.append( (tmpY, a, c) )

	for (tmpY,a,color) in todo:
		if color == None:
			print "0.001 cm setlinewidth"
			color = "black"
		else:
			print "0.1 cm setlinewidth"
		utils.myPsOutput.drawLine(x, tmpY, -a*dx, 0, color)
		utils.myPsOutput.drawLine(x, tmpY, 0, (mi+ma)/2.-tmpY, color)


	if arguments["printAncestors"]:
		utils.myPsOutput.drawText(x + .1, (mi+ma)/2. + .1, node, "black")
	if arguments["printAges"]:
		utils.myPsOutput.drawText(x + .1, (mi+ma)/2. - .4, "(%d)" % phylTree.ages[node], "black")

	return (mi+ma)/2.

printSubTree(root)

utils.myPsOutput.printPsFooter()

