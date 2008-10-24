
# Module de dessin d'un karyotype

import sys
import itertools

import myTools
import myPsOutput

myTools.addModuleOptions("karyo", \
	[("dx",float,0), ("dy",float,0), ("roundedChr",bool,False), ("landscape",bool,False), ("showText",bool,True), ("drawBorder",bool,False), \
	("defaultColor",str,""), ("penColor",str,"black"), ("backgroundColor",str,"")] \
)


def drawKaryo(data, arguments):

	(largeur,hauteur) = myPsOutput.printPsHeader()
	defaultColor = arguments["karyo:defaultColor"] if len(arguments["karyo:defaultColor"]) > 0 else None

	# Couleur de fond
	if len(arguments["karyo:backgroundColor"]) > 0:
		myPsOutput.drawBox(0,0, largeur,hauteur, arguments["karyo:backgroundColor"], arguments["karyo:backgroundColor"])

	# None correspond a la legende (tout a droite)
	lstC = [c for c in data if c is not None]
	if None in data:
		lstC.append(None)
	
	# Gestion de plusieurs pistes par chromosome
	for c in lstC:
		l = data[c]
		data[c] = l if isinstance(l, tuple) else [l]

	# On dessine
	dx = arguments["karyo:dx"] if arguments["karyo:dx"] > 0 else ((largeur-2.)*3.)/(5.*len(data)-2.)
	dy = arguments["karyo:dy"] if arguments["karyo:dy"] > 0 else (hauteur-4.) / float(max([len(x[0]) for x in data.itervalues()]))
	y0 = 2
	xx = 1

	for c in lstC:
		all = data[c]
		size = len(all[0])

		def printBorder():
			print "newpath"
			if arguments["karyo:roundedChr"]:
				print "%.5f cm %.5f cm %.5f cm 180 360 arc" % (xx+dx/2., y0+dx/2., dx/2.)
				print "%.5f %.5f 2cm rlineto" % (0,size*dy-dx)
				print "%.5f cm %.5f cm %.5f cm 0 180 arc" % (xx+dx/2., y0+size*dy-dx/2., dx/2.)
				print "%.5f %.5f 2cm rlineto" % (0,-size*dy+dx)
			else:
				print "%.5f %.5f 2cm moveto" % (xx,y0)
				print "%.5f %.5f 2cm rlineto" % (0,size*dy)
				print "%.5f %.5f 2cm rlineto" % (dx,0)
				print "%.5f %.5f 2cm rlineto" % (0,-size*dy)
			print "closepath"

		if arguments["karyo:roundedChr"]:
			printBorder()
			print "clip"

		nbl = float(len(all))
		for (i,currl) in enumerate(all):
			y = y0
			for (col,items) in itertools.groupby(currl):
				hauteur = len(list(items)) * dy
				if col is None:
					col = defaultColor
				myPsOutput.drawBox(xx+i*dx/nbl, y, dx/nbl, hauteur, col, col)
				y += hauteur
		
		if arguments["karyo:drawBorder"]:
			myPsOutput.setColor(arguments["karyo:penColor"], "color")
			printBorder()
			print "stroke"
		
		if arguments["karyo:roundedChr"]:
			print "initclip"

		if arguments["karyo:showText"] and (c is not None):
			myPsOutput.drawText(xx, y0-1, str(c), arguments["karyo:penColor"])
		
		xx += (5.*dx)/3.

	myPsOutput.printPsFooter()

