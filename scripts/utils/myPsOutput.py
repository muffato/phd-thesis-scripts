#! /usr/bin/python2.4

# Module d'ecriture dans un fichier PostScript

import myTools

# Les definitions des couleurs
colorTable = {}
color = {}


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
# Charge le fichier de definitions des couleurs en RGB
# Initialise toutes les couleurs couramment utilisees
#
def initColor():

	f = myTools.myOpenFile("~/work/scripts/utils/rgb.txt")
	for l in f:
		c = l.split()
		colorTable["".join(c[6:])] = " ".join(c[:3])

	f.close()


	lightColors = ["PaleTurquoise2", "khaki1", "DarkSeaGreen1", "PaleTurquoise2", "DarkOliveGreen1", "khaki1", "lavender", "LightBlue", "salmon", "SkyBlue1", "LightGoldenrod3", "wheat1", "thistle2", "PeachPuff"]
	
	darkColors = ["red1", "turquoise2", "DarkGreen", "yellow", "coral2", "OliveDrab2", "orange", "MediumAquamarine", "blue2", "firebrick4", "LightSalmon", "DarkViolet", "magenta2", "DarkSeaGreen4", "DarkSlateBlue", "yellow4", "grey62", "gold", "PeachPuff2", "HotPink4", "firebrick", "purple4"]
	
	for i in range(len(lightColors)):
		color[str(-(i+1))] = lightColors[i]
		color[chr(97+i)] = lightColors[i]
		color[str((i+1+len(darkColors)))] = lightColors[i]
		color[chr(65+i+len(darkColors))] = lightColors[i]
		printColorDefinitionLine(lightColors[i])
	
	for i in range(len(darkColors)):
		color[str(i+1)] = darkColors[i]
		color[chr(65+i)] = darkColors[i]
		color[str(-(i+1+len(lightColors)))] = darkColors[i]
		color[chr(97+i+len(lightColors))] = darkColors[i]
		printColorDefinitionLine(darkColors[i])
	
	printColorDefinitionLine("black")
	printColorDefinitionLine("white")

	print


#
# Permet de rajouter une couleur creee en direct
#
def getColor(s, d):

	if s in color:
		return color[s]

	elif s in colorTable:
		printColorDefinitionLine(s)
		return s

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


