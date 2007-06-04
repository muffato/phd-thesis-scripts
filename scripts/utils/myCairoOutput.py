#! /users/ldog/muffato/python -OO

import cairo
import sys

# Module de sortie graphique via Cairo

# Les petits noms que je donne a mes couleurs (nom -> RGB)
colorTransl = {}



def printPsHeader(landscape = False, format = "ps", filename = '/dev/stdout'):

	format = format.strip().lower()
	
	if landscape:
		width,height = 29.7, 21
	else:
		width,height = 21, 29.7
	
	if format == 'ps':
		surface = cairo.PSSurface(filename, width*72/2.54, height*72/2.54)
	elif format == 'pdf':
		surface = cairo.PSSurface(filename, width*72/2.54, height*72/2.54)
	else:
		print >> sys.stderr, "Unknown output style '%s'" % format
		sys.exit(1)
	
	global ctx
	ctx = cairo.Context(surface)
	ctx.scale(width, height)
	ctx.set_line_width(0.001)
	#ctx.scale(width*72/2.54, height*72/2.54)

	initColor()

#
# Le pied de page PostScript
#
def printPsFooter():
	global ctx
	ctx.show_page()



#
# Charge le fichier de definitions des couleurs en RGB
# Initialise toutes les couleurs couramment utilisees
#
def initColor():

	# La liste des couleurs et leurs valeurs RGB
	f = open("/users/ldog/muffato/work/scripts/utils/rgb.txt", 'r')
	for l in f:
		c = l.split()
		s = "".join(c[6:]).lower()
		colorTransl[s] = tuple([float(x) for x in c[:3]])

	f.close()


	lightColors = ["salmon", "PaleTurquoise2", "DarkSeaGreen1","khaki1", "thistle2", "PeachPuff","LightBlue", "SkyBlue1","LightGoldenrod3", "wheat1",  "DarkOliveGreen1", "lavender"]
	
	darkColors = ["red1", "turquoise2", "DarkGreen", "yellow", "coral2", "OliveDrab2", "orange", "MediumAquamarine", "blue2", "firebrick4", "LightSalmon", "DarkViolet", "magenta2", "DarkSeaGreen4", "DarkSlateBlue", "yellow4", "grey62", "gold", "PeachPuff2", "HotPink4", "firebrick", "purple4"]
	#darkColors = ["turquoise2", "gray85", "yellow", "coral2", "OliveDrab2", "MediumAquamarine", "blue2", "LightSalmon", "magenta2", "DarkSeaGreen2", "PaleTurquoise3", "khaki2", "thistle3", "PeachPuff2", "HotPink2", "firebrick", "purple2"]

	greekLetters = ["ALPHA", "BETA", "DELTA", "EPSILON", "GAMMA", "PHI"]
	craniateColors = ["DarkOrange", "RoyalBlue4", "chartreuse4", "gold", "DarkOrchid4", "red3"]
	#craniateColors = ["DarkOrange", "RoyalBlue2", "chartreuse2", "gold", "DarkOrchid1", "red2"]

	# Les couleurs claires sont pour les nombres negatifs et les lettres minuscules
	ordre = lightColors + craniateColors + darkColors
	for i in range(len(ordre)):
		colorTransl[str(-(i+1))] = colorTransl[ordre[i].lower()]
		if i < 26:
			colorTransl[chr(i+97)] = colorTransl[ordre[i].lower()]
		
	# Les couleurs foncees sont pour les nombres positifs et les lettres majuscules
	ordre = darkColors + lightColors + craniateColors
	#ordre = craniateColors + darkColors + lightColors
	for i in range(len(ordre)):
		colorTransl[str(i+1)] = colorTransl[ordre[i].lower()]
		if i < 26:
			colorTransl[chr(i+65)] = colorTransl[ordre[i].lower()]

	for i in range(len(greekLetters)):
		colorTransl[greekLetters[i]] = colorTransl[craniateColors[i].lower()]

#
# Permet de rajouter une couleur creee a partir de valeurs rgb
#
def getColor(s, d):

	s = str(s)
	if s in colorTransl:
		return colorTransl[s]

	elif s[0] == '#':
		return tuple([float(x)/255. for x in s[1:].split(':')])
	
	else:
		#print >> sys.stderr, s
		return d



#######################
# Fonctions de dessin #
#######################

# Les coordonnes sont en cm (portrait/paysage)
#   X: de gauche (0) a droite (21/29.7)
#   Y: de bas (0) en haut (29.7/21)
# Les couleurs sont facultatives, utiliser None pour ne pas en specifier.
#   Sinon, utiliser le nom UNIX (cf getColor pour avoir ce nom)
# Les angles sont en degres

def setLineColor(C):
	global ctx

	if type(C) != tuple:
		C = getColor(C, (0,0,0))
	ctx.set_source_rgb(*C)

def drawLine(X, Y, L, H, C):

	global ctx
	
	setLineColor(C)
	ctx.move_to(X, Y)
	ctx.rel_line_to(L, H)
	ctx.stroke()


def drawBox(X, Y, L, H, Cb, Cr):

	setLineColor(Cb)
	ctx.rectangle(X, Y, L, H)

	if Cr != None:
		#print >> sys.stderr, Cr
		ctx.set_source(cairo.SolidPattern(Cr[0], Cr[1], Cr[2], 1))
		ctx.fill()


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
	if C in colorListUNIX:
		print C, "color"
	print "(%s) %.3f %.3f mytext" % (T, X,Y)


