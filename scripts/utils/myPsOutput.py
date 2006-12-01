#! /users/ldog/muffato/python

# Module d'ecriture dans un fichier PostScript

import myTools

# Les definitions des couleurs
color = {}
colorTable = {}

#
# L'en-tete PostScript
#
def printPsHeader(linewidth = 0.001):
	print "%!PS-Adobe-3.0"
	print "%%DocumentData: Clean7bit"
	print "%%Creator: myPsOutput"
	print "%%PageOrder: Ascend"
	print "%%Pages: 1"
	print "%%DocumentFonts: Helvetica"
	print "%%LanguageLevel: 1"
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
	print linewidth, " cm setlinewidth"

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
	(r,g,b) = colorTable[C]
	print "/%s [%f %f %f] def" % (C, r,g,b)

#
# Charge le fichier de definitions des couleurs en RGB
# Initialise toutes les couleurs couramment utilisees
#
def initColor():

	f = myTools.myOpenFile("~/work/scripts/utils/rgb.txt", 'r')
	for l in f:
		c = l.split()
		colorTable["".join(c[6:])] = tuple([float(x) for x in c[:3]])

	f.close()


	lightColors = ["salmon", "PaleTurquoise2", "DarkSeaGreen1","khaki1", "thistle2", "PeachPuff","LightBlue", "SkyBlue1","LightGoldenrod3", "wheat1",  "DarkOliveGreen1", "lavender"]
	
	darkColors = ["red1", "turquoise2", "DarkGreen", "yellow", "coral2", "OliveDrab2", "orange", "MediumAquamarine", "blue2", "firebrick4", "LightSalmon", "DarkViolet", "magenta2", "DarkSeaGreen4", "DarkSlateBlue", "yellow4", "grey62", "gold", "PeachPuff2", "HotPink4", "firebrick", "purple4"]

	greekLetters = ["ALPHA", "BETA", "DELTA", "EPSILON", "GAMMA", "PHI"]
	craniateColors = ["DarkOrange", "RoyalBlue4", "chartreuse4", "gold", "DarkOrchid4", "red3"]

	ordre = lightColors + craniateColors + darkColors
	for i in range(len(ordre)):
		printColorDefinitionLine(ordre[i])
		color[str(-(i+1))] = ordre[i]
		color[chr(i+97)] = ordre[i]

	ordre = darkColors + lightColors + craniateColors
	for i in range(len(ordre)):
		color[str(i+1)] = ordre[i]
		color[chr(i+65)] = ordre[i]
		
	for i in range(len(greekLetters)):
		color[greekLetters[i]] = craniateColors[i]

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
		
		{}["tmp"] = (r/255., g/255., b/255.)
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
		print "%f %f %f setrgbcolor" % colorTable[C]
		
	print X, "cm", Y, "cm", "moveto"
	print L, "cm", H, "cm", "rlineto"
	print "closepath"
	print "stroke"


def drawBox(X, Y, L, H, Cb, Cr):
	print "newpath"
	if Cb in colorTable:
		print "%f %f %f setrgbcolor" % colorTable[Cb]
		
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


