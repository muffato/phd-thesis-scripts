#! /usr/bin/python2.4

# Module d'ecriture dans un fichier PostScript

import myTools

# Les definitions des couleurs
colorTable = {}
color = {}

#
# Charge le fichier de definitions des couleurs en RGB
#
def loadColorTable():

	f = myTools.myOpenFile("~/work/scripts/utils/rgb.txt")
	for l in f:
		c = l.split()
		colorTable["".join(c[6:])] = " ".join(c[:3])

	f.close()


#
# L'en-tete PostScript
#
def printPsHeader():
	print "%!PS-Adobe-3.0"
	print "%%DocumentData: Clean8bit"
	print "%%PageOrder: Ascend"
	print "%%Pages: 1"
	print "%%DocumentFonts: Helvetica"
	print "%%EndComments"
	
	print "%%Page: 1 1"
	print "/cm {28.3464567 mul} def"
	#print "%%Bounding-Box: 0 cm 0 cm 21 cm 29.7 cm";

	print
	print "/Arial findfont"
	print "10 scalefont"
	print "setfont"
	# style de trait (coins ronds)
	print "1 setlinejoin"
	print "0.00001 cm setlinewidth"

	# definition des couleurs dns l'en-tete du Postscript
	print
	print "/color {"
	print "\taload pop setrgbcolor"
	print "} def"
	print

#
# Le pied de page PostScript
#
def printPsFooter():
	print "showpage"


#
# Definit une nouvelle couleur dans le fichier PostScript
#
def printColorDefinitionLine(C):
	print "/%s [%s] def" % (C, colorTable[C])

#
# Initialise toutes les couleurs couramment utilisees
#
def initColor():

	color["0"] = "black"
	color["1"] = "purple4"
	color["2"] = "magenta2"
	color["3"] = "grey62"
	color["4"] = "coral2"
	color["5"] = "firebrick4"
	color["6"] = "yellow4"
	color["7"] = "LightSalmon"
	color["8"] = "DarkSeaGreen4"
	color["9"] = "turquoise2"
	color["10"] = "blue2"
	color["11"] = "PeachPuff2"
	color["12"] = "DarkViolet"
	color["13"] = "OliveDrab2"
	color["14"] = "MediumAquamarine"
	color["16"] = "DarkSlateBlue"
	color["15"] = "yellow"
	color["17"] = "DarkGreen"
	color["18"] = "gold"
	color["19"] = "HotPink4"
	color["20"] = "red1"
	color["21"] = "orange"
	color["22"] = "PaleTurquoise2"
	color["23"] = "khaki1"
	color["24"] = "DarkSeaGreen1"
	color["25"] = "PaleTurquoise2"
	color["26"] = "DarkOliveGreen1"
	color["27"] = "khaki1"
	color["28"] = "lavender"
	color["29"] = "LightBlue"
	color["30"] = "salmon"
	color["31"] = "SkyBlue1"
	color["32"] = "LightGoldenrod3"
	color["33"] = "wheat1"
	color["34"] = "thistle2"
	color["35"] = "PeachPuff"

	color["41"] = "coral2"
	color["51"] = "firebrick4"
	color["61"] = "yellow4"
	color["71"] = "LightSalmon"
	color["81"] = "DarkSeaGreen4"
	color["90"] = "turquoise2"
	color["91"] = "purple4"
	color["100"] = "magenta2"
	color["101"] = "grey62"
	color["120"] = "DarkViolet"
	color["121"] = "OliveDrab2"
	color["130"] = "MediumAquamarine"
	color["131"] = "DarkSlateBlue"
	color["140"] = "yellow"
	color["141"] = "DarkGreen"

	color["A"] = "DarkSeaGreen1"
	color["B"] = "PaleTurquoise2"
	color["C"] = "DarkOliveGreen1"
	color["D"] = "khaki1"
	color["E"] = "wheat1"
	color["F"] = "thistle2"
	color["G"] = "PeachPuff"
	color["H"] = "lavender"
	color["I"] = "LightBlue"
	color["J"] = "salmon"
	color["K"] = "LightGoldenrod3"
	color["L"] = "SkyBlue1"

	color["A"] = "red1"
	color["B"] = "turquoise2"
	color["C"] = "DarkGreen"
	color["D"] = "yellow"
	color["E"] = "coral2"
	color["F"] = "OliveDrab2"
	color["G"] = "orange"
	color["H"] = "MediumAquamarine"
	color["I"] = "blue2"
	color["J"] = "firebrick4"
	color["L"] = "DarkViolet"
	color["K"] = "gold"
	color["M"] = "firebrick"
	color["N"] = "lavender"

	for c in color:
		printColorDefinitionLine(color[c])
	print


#
# Permet de rajouter une couleur creee en direct
#
def getColor(s, d):

	if s in color:
		return color[s]
		
	elif s[0] == '#':
		r = float(s[1:4])
		g = float(s[5:8])
		b = float(s[9:12])
		
		colorTable["tmp"] = "%f %f %f" % (r/255., g/255., b/255.)
		printColorDefinitionLine("tmp")
		return "tmp"
	else:
		return d



#######################
# Fonctions de dessin #
#######################

def drawLine(X, Y, L, H, C):
	print "newpath"
	if C in colorTable:
		print colorTable[C], " setrgbcolor"
		
	print X, "cm", Y, "cm", "moveto"
	print L, "cm", H, "cm", "rlineto"
	print "closepath"
	print "stroke"


def drawBox(X, Y, L, H, Cb, Cr):
	print "newpath"
	if Cb in colorTable:
		print colorTable[Cb], " setrgbcolor"
		
	print X, "cm", Y, "cm", "moveto"
	print L, "cm", 0, "cm", "rlineto"
	print 0, "cm", H, "cm", "rlineto"
	print -L, "cm", 0, "cm", "rlineto"
	print "closepath"
	
	if Cr != 0:
		print "gsave"
		print Cr, "color fill"
		print "grestore"

	print "stroke"


def drawCross(X, Y, L, H):
	print "newpath"
	print X, "cm", Y, "cm", "moveto"
	print L, "cm", H, "cm", "rlineto"
	print X+L, "cm", Y, "cm", "moveto"
	print -L, "cm", H, "cm", "rlineto"
	print "closepath"
	print "stroke"

def drawCircle(X, Y, R, A, B):
	print "newpath"
	print X, "cm", Y, "cm", R, "cm", A, "cm", B, "cm", "arc"
	print "stroke"

def drawFilledCircle(X, Y, R, A, B, C):
	print "newpath"
	print X, "cm", Y, "cm", R, "cm", A, "cm", B, "cm", "arc"
	print "gsave"
	print C, "color fill"
	print "grestore"
	print "stroke"

def drawArrowR(X, Y, L, H, P, C):
	print "newpath"
	print X, "cm", Y, "cm", "moveto"
	print L, "cm", 0, "cm", "rlineto"
	print P, "cm", H/2, "cm", "rlineto"
	print -P, "cm", H/2, "cm", "rlineto"
	print -L, "cm", 0, "cm", "rlineto"
	print "closepath"
	print "gsave"
	print C, "color fill"
	print "grestore"
	print "stroke"

def drawArrowL(X, Y, L, H, P, C):
	print "newpath";
	print X, "cm", Y+(H/2), "cm", "moveto"
	print P, "cm", H/2, "cm", "rlineto"
	print L, "cm", 0, "cm", "rlineto"
	print 0, "cm", -H, "cm", "rlineto"
	print -L, "cm", 0, "cm", "rlineto"
	print "closepath"
	print "gsave"
	print C, "color fill"
	print "grestore"
	print "stroke"

def drawText(X, Y, T, C):
    print "newpath";
    print C, "color";
    print X, "cm", Y, "cm", "moveto";
    print "(%s)" % T, "show";


