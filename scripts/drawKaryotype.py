#! /users/ldog/muffato/python -OO

__doc__ = """
	Dessine le karyotype d'un genome face a l'autre
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
(noms_fichiers, options) = utils.myTools.checkArgs( \
	["studiedGenome", "referenceGenome"], \
	[("orthologuesList",str,""), ("includeGaps",bool,False), ("includeScaffolds",bool,False), ("includeRandoms",bool,False), \
	("reverse",bool,False), ("dx",float,0), ("dy",float,0), ("roundedChr",bool,False), ("landscape",bool,False), \
	("defaultColor",str,"black"), ("penColor",str,"black"), ("backgroundColor",str,"")], \
	__doc__
)


# Chargement des fichiers
genome1 = utils.myGenomes.loadGenome(noms_fichiers["studiedGenome"])
genome2 = utils.myGenomes.loadGenome(noms_fichiers["referenceGenome"])
if options["reverse"]:
	x = genome1
	genome1 = genome2
	genome2 = x
if options["orthologuesList"] != "":
	genesAnc = utils.myGenomes.loadGenome(options["orthologuesList"])
else:
	genesAnc = None

# Les chromosomes a etudier
chr1 = genome1.lstChr
chr2 = genome2.lstChr
if options["includeScaffolds"]:
	chr1.extend(genome1.lstScaff)
	chr2.extend(genome2.lstScaff)
if options["includeRandoms"]:
	chr1.extend(genome1.lstRand)
	chr2.extend(genome2.lstRand)

table12 = genome1.buildOrthosTable(chr1, genome2, chr2, options["includeGaps"], genesAnc)

print >> sys.stderr, "Affichage ...",

(largeur,hauteur) = utils.myPsOutput.printPsHeader()

if len(options["backgroundColor"]) > 0:
	utils.myPsOutput.drawBox(0,0, largeur,hauteur, options["backgroundColor"], options["backgroundColor"])

# On dessine
if options["dx"] > 0:
	dx = options["dx"]
else:
	dx = ((largeur-2.)*3.)/(5.*len(chr1)-2.)
if options["dy"] > 0:
	dy = options["dy"]
else:
	dy = (hauteur-4.) / float(max([len(x) for x in table12.values()]))
xx = 1
y0 = 1.


if options["roundedChr"]:
	print "newpath"
	for c in chr1:
		if options["dy"] < 0:
			dy = (hauteur-4.) / float(len(table12[c]))
		y = y0+1
		print "%.5f cm %.5f cm %.5f cm 180 360 arc" % (xx+dx/2., y0+1+dx/2., dx/2.)
		print "%.5f %.5f 2cm rlineto" % (0,len(table12[c])*dy-dx)
		print "%.5f cm %.5f cm %.5f cm 0 180 arc" % (xx+dx/2., y0+1+len(table12[c])*dy-dx/2., dx/2.)
		print "%.5f %.5f 2cm rlineto" % (0,-len(table12[c])*dy+dx)
		xx += (5.*dx)/3.
	print "0 0.5 moveto"
	print "0 1.5 2cm rlineto"
	print "30 0 2cm rlineto"
	print "0 -1.5 2cm rlineto"
	print "-30 0 2cm rlineto"
	print "closepath"
	print "clip"

xx = 1

for c in chr1:
	
	if options["dy"] < 0:
		dy = (hauteur-4.) / float(len(table12[c]))

	utils.myPsOutput.drawText(xx, y0, str(c), options["penColor"])
	y = y0 + 1

	
	last = ""
	nb = 0
	for (i,tmp) in table12[c]:
		if len(tmp) == 0:
			col = options["defaultColor"]
		else:
			col = tmp[0][0]
		if col == last:
			nb += 1
		else:
			if nb != 0:
				utils.myPsOutput.drawBox(xx, y, dx, nb*dy, last, last)
				y += nb*dy
			last = col
			nb = 1
	if nb != 0:
		utils.myPsOutput.drawBox(xx, y, dx, nb*dy, last, last)
	xx += (5.*dx)/3.

utils.myPsOutput.printPsFooter()
print >> sys.stderr, "OK"


