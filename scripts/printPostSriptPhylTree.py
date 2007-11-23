#! /users/ldog/muffato/python -OO

# Librairies
import sys
import math
import utils.myTools
import utils.myPhylTree
import utils.myPsOutput

# Arguments
(noms_fichiers, options) = utils.myTools.checkArgs( \
	["phylTree.conf"], \
	[("landscape",bool,False), ("printSpecies",bool,True), ("printAncestors",bool,True), ("printAges",bool,False), \
	("colorFile",str,""), ("func",str,"")], \
	"Dessine l'arbre phylogenetique avec des informations dessus comme des evolutions de taux" \
)

phylTree = utils.myPhylTree.PhylogeneticTree(noms_fichiers["phylTree.conf"])

(largeur,hauteur) = utils.myPsOutput.printPsHeader(landscape=options["landscape"])
maxAge = phylTree.ages[phylTree.root]

if options["func"] == "":
	func = lambda x: x
else:
	func = eval(options["func"])

colors = {}
if options["colorFile"] != "":
	f = utils.myTools.myOpenFile(options["colorFile"], "r")
	for l in f:
		t = l.split("\t")
		s1 = intern(t[0])
		s2 = intern(t[1])
		if (s1 in phylTree.officialName) and (s2 in phylTree.officialName):
			x = func(float(t[2]))
			colors[(s1,s2)] = x
			colors[(s2,s1)] = x
	f.close()
	val = colors.values()
	minV = float(min(val))
	maxV = float(max(val))


def getColor(value):
	value = (value-minV) / (maxV-minV)
	if value < 0.1:
		return utils.myPsOutput.alphaColor( (0,0,127), (0,0,255), 10.*value)
	elif value < 0.3:
		return utils.myPsOutput.alphaColor( (0,0,255), (0,255,255), 5.*(value-.1))
	elif value < 0.5:
		return utils.myPsOutput.alphaColor( (0,255,255), (0,255,0), 5.*(value-.3))
	elif value < 0.7:
		return utils.myPsOutput.alphaColor( (0,255,0), (255,255,0), 5.*(value-.5))
	elif value < 0.9:
		return utils.myPsOutput.alphaColor( (255,255,0), (255,0,0), 5.*(value-.7))
	else:
		return utils.myPsOutput.alphaColor( (255,0,0), (127,0,0), 10.*(value-.9))


y = 0
dy = hauteur / (len(phylTree.listSpecies)+1.)

margeX1 = dy
if options["printSpecies"]:
	margeX1 += max([len(x) for x in phylTree.listSpecies]) * 0.15

margeX2 = dy
if options["printAncestors"]:
	margeX2 += len(phylTree.root) * 0.15

dx = (largeur-margeX1-margeX2) / phylTree.ages[phylTree.root]


def printSubTree(node):
	global y
	
	if node not in phylTree.items:
		y += dy
		if options["printSpecies"]:
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


	if options["printAncestors"]:
		utils.myPsOutput.drawText(x + .1, (mi+ma)/2. + .1, node, "black")
	if options["printAges"]:
		utils.myPsOutput.drawText(x + .1, (mi+ma)/2. - .4, "(%d)" % phylTree.ages[node], "black")

	return (mi+ma)/2.


printSubTree(phylTree.root)

utils.myPsOutput.printPsFooter()

