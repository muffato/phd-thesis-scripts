#! /users/ldog/muffato/python -OO

# Module d'ecriture dans un fichier PostScript

# Les petits noms que je donne a mes couleurs (nom -> nom UNIX)
colorTransl = {}
# Les definitions des couleurs (valeurs RGB -> nom UNIX) pour eviter de creer plusieurs fois la meme couleur
colorTableRGB2UNIX = {}
# La liste des noms de couleurs UNIX
colorListUNIX = set()


#
# L'en-tete PostScript
#
def printPsHeader(linewidth = 0.001, landscape = False):
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

	print
	print "/Arial findfont"
	print "10 scalefont"
	print "setfont"
	# style de trait (coins ronds)
	print "1 setlinejoin"
	print linewidth, " cm setlinewidth"

	if landscape:
		print "90 rotate"
		print "0 -21 cm translate"
		print
	
	# definition des couleurs dans l'en-tete du Postscript
	print
	print "/color {"
	print "\taload pop setrgbcolor"
	print "} def"
	print
	initColor()
	print

#
# Le pied de page PostScript
#
def printPsFooter():
	print "showpage"



#
# Charge le fichier de definitions des couleurs en RGB
# Initialise toutes les couleurs couramment utilisees
#
def initColor():

	# La liste des couleurs et leurs valeurs RGB
	f = open("/users/ldog/muffato/work/scripts/utils/rgb.txt", 'r')
	for l in f:
		c = l.split()
		s = "".join(c[6:])
		colorTableRGB2UNIX[tuple([int(x) for x in c[3:6]])] = s
		print "/%s [%s] def" % (s, " ".join(c[:3]))
		colorListUNIX.add(s)

	f.close()


	lightColors = ["salmon", "PaleTurquoise2", "DarkSeaGreen1","khaki1", "thistle2", "PeachPuff","LightBlue", "SkyBlue1","LightGoldenrod3", "wheat1",  "DarkOliveGreen1", "lavender"]
	
	darkColors = ["red1", "turquoise2", "DarkGreen", "yellow", "coral2", "OliveDrab2", "orange", "MediumAquamarine", "blue2", "firebrick4", "LightSalmon", "DarkViolet", "magenta2", "DarkSeaGreen4", "DarkSlateBlue", "yellow4", "grey62", "gold", "PeachPuff2", "HotPink4", "firebrick", "purple4"]

	greekLetters = ["ALPHA", "BETA", "DELTA", "EPSILON", "GAMMA", "PHI"]
	craniateColors = ["DarkOrange", "RoyalBlue4", "chartreuse4", "gold", "DarkOrchid4", "red3"]

	# Les couleurs claires sont pour les nombres negatifs et les lettres minuscules
	ordre = lightColors + craniateColors + darkColors
	for i in range(len(ordre)):
		colorTransl[str(-(i+1))] = ordre[i]
		if i < 26:
			colorTransl[chr(i+97)] = ordre[i]
		
	# Les couleurs foncees sont pour les nombres positifs et les lettres majuscules
	ordre = darkColors + lightColors + craniateColors
	for i in range(len(ordre)):
		colorTransl[str(i+1)] = ordre[i]
		if i < 26:
			colorTransl[chr(i+65)] = ordre[i]

	for i in range(len(greekLetters)):
		colorTransl[greekLetters[i]] = craniateColors[i]

#
# Permet de rajouter une couleur creee a partir de valeurs rgb
#
def getColor(s, d):

	s = str(s)
	if s in colorTransl:
		return colorTransl[s]

	elif s[0] == '#':
		(r,g,b) = [int(x) for x in s[1:].split(':')]
		if (r,g,b) in colorTableRGB2UNIX:
			return colorTableRGB2UNIX[(r,g,b)]
		else:
			colName = "tmp%d" % len(colorListUNIX)
			colorTableRGB2UNIX[(r,g,b)] = colName
			print "/%s [%f %f %f] def" % (colName, float(r)/255., float(g)/255., float(b)/255.)
			return colName
	else:
		return d



#######################
# Fonctions de dessin #
#######################

# Les coordonnes sont en cm
#   X: de gauche (0) a droite (21/29.7)
#   Y: de bas (0) en haut (29.7/21)
# Les couleurs sont facultatives, utiliser "0" ou 0 ou None pour ne pas en specifier.
#   Sinon, utiliser le nom UNIX (cf getColor pour avoir ce nom)
# Les angles sont en degres

def drawLine(X, Y, L, H, C):
	print "newpath"
	if C in colorListUNIX:
		print C, "color"
		
	print "%.3f cm %.3f cm moveto" % (X, Y)
	print "%.3f cm %.3f cm rlineto" % (L, H)
	print "closepath"
	print "stroke"


def drawBox(X, Y, L, H, Cb, Cr):
	print "newpath"
	if Cb in colorListUNIX:
		print Cb, "color"
		
	print "%.3f cm %.3f cm moveto" % (X, Y)
	print "%.3f cm 0 cm rlineto" % L
	print "0 cm %.3f cm rlineto" % H
	print "%.3f cm 0 cm rlineto" % (-L)
	print "closepath"
	
	if Cr in colorListUNIX:
		print "gsave"
		print Cr, "color fill"
		print "grestore"

	print "stroke"


def drawCross(X, Y, L, H, C):
	print "newpath"
	if C in colorListUNIX:
		print C, "color"
	print "%.3f cm %.3f cm moveto" % (X, Y)
	print "%.3f cm %.3f cm rlineto" % (L, H)
	print "%.3f cm %.3f cm moveto" % (X+L, Y)
	print "%.3f cm %.3f cm rlineto" % (-L, H)
	print "closepath"
	print "stroke"

def drawCircle(X, Y, R, A, B, Cb, Cr):
	print "newpath"
	if Cb in colorListUNIX:
		print Cb, "color"
	print "%.3f cm %.3f cm %.3f cm %.3f %.3f arc" % (X, Y, R, A, B)
	if Cr in colorListUNIX:
		print "gsave"
		print Cr, "color fill"
		print "grestore"
	print "stroke"

def drawArrowR(X, Y, L, H, P, Cb, Cr):
	print "newpath"
	if Cb in colorListUNIX:
		print Cb, "color"
	print "%.3f cm %.3f cm moveto" % (X, Y)
	print "%.3f cm 0 cm rlineto" % L
	print "%.3f cm %.3f cm rlineto" % (P, H/2)
	print "%.3f cm %.3f cm rlineto" % (-P, H/2)
	print "%.3f cm 0 cm rlineto" % (-L)
	print "closepath"
	if Cr in colorListUNIX:
		print "gsave"
		print Cr, "color fill"
		print "grestore"
	print "stroke"

def drawArrowL(X, Y, L, H, P, Cb, Cr):
	print "newpath";
	if Cb in colorListUNIX:
		print Cb, "color"
	print "%.3f cm %.3f cm moveto" % (X, Y+(H/2))
	print "%.3f cm %.3f cm rlineto" % (P, H/2)
	print "%.3f cm 0 cm rlineto" % L
	print "0 cm %.3f cm rlineto" % (-H)
	print "%.3f cm 0 cm rlineto" % (-L)
	print "closepath"
	if Cr in colorListUNIX:
		print "gsave"
		print Cr, "color fill"
		print "grestore"
	print "stroke"

def drawText(X, Y, T, C):
	print "newpath"
	if C in colorListUNIX:
		print C, "color"
	print "%.3f cm %.3f cm moveto" % (X, Y)
	print "(%s)" % T, "show"


