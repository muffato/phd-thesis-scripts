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
	# En-tete Postscript
	print """
%!PS-Adobe-3.0
%%DocumentData: Clean7bit
%%Creator: myPsOutput
%%PageOrder: Ascend
%%Pages: 1
%%DocumentFonts: Helvetica
%%LanguageLevel: 1
%%EndComments

/color { aload pop setrgbcolor } def
/cm {28.3464567 mul} def
/2cm { cm exch cm exch } def
/myline { newpath 2cm moveto 2cm rlineto closepath stroke } def
/mybox { newpath 2cm moveto 2cm exch dup 0 rlineto exch 0 exch rlineto neg 0 rlineto closepath } def
/myfill { gsave color fill grestore } def
/mytext { 2cm moveto show } def

1 setlinejoin

%%BeginProlog

% Font encoding-vector for iso-8859-1 (latin-1)
/font_encoding_vector [
/.notdef /.notdef /.notdef /.notdef /.notdef /.notdef /.notdef /.notdef 
/.notdef /.notdef /.notdef /.notdef /.notdef /.notdef /.notdef /.notdef 
/.notdef /.notdef /.notdef /.notdef /.notdef /.notdef /.notdef /.notdef 
/.notdef /.notdef /.notdef /.notdef /.notdef /.notdef /.notdef /.notdef 
/space /exclam /quotedbl /numbersign /dollar /percent /ampersand /quoteright 
/parenleft /parenright /asterisk /plus /comma /hyphen /period /slash 
/zero /one /two /three /four /five /six /seven /eight /nine /colon 
/semicolon /less /equal /greater /question /at /A /B /C /D /E /F /G 
/H /I /J /K /L /M /N /O /P /Q /R /S /T /U /V /W /X /Y /Z /bracketleft 
/backslash /bracketright /asciicircum /underscore /quoteleft /a /b /c 
/d /e /f /g /h /i /j /k /l /m /n /o /p /q /r /s /t /u /v /w /x /y /z 
/braceleft /bar /braceright /asciitilde /.notdef 
/.notdef /.notdef /.notdef /.notdef /.notdef /.notdef /.notdef /.notdef 
/.notdef /.notdef /.notdef /.notdef /.notdef /.notdef /.notdef /.notdef 
/.notdef /.notdef /.notdef /.notdef /.notdef /.notdef /.notdef /.notdef 
/.notdef /.notdef /.notdef /.notdef /.notdef /.notdef /.notdef /.notdef 
/space /exclamdown /cent /sterling /currency /yen /brokenbar /section 
/dieresis /copyright /ordfeminine /guillemotleft /logicalnot /hyphen 
/registered /macron /degree /plusminus /twosuperior /threesuperior /acute 
/mu /paragraph /bullet /cedilla /dotlessi /ordmasculine /guillemotright 
/onequarter /onehalf /threequarters /questiondown /Agrave /Aacute 
/Acircumflex /Atilde /Adieresis /Aring /AE /Ccedilla /Egrave /Eacute 
/Ecircumflex /Edieresis /Igrave /Iacute /Icircumflex /Idieresis /Eth 
/Ntilde /Ograve /Oacute /Ocircumflex /Otilde /Odieresis /multiply /Oslash 
/Ugrave /Uacute /Ucircumflex /Udieresis /Yacute /Thorn /germandbls 
/agrave /aacute /acircumflex /atilde /adieresis /aring /ae /ccedilla 
/egrave /eacute /ecircumflex /edieresis /igrave /iacute /icircumflex 
/idieresis /eth /ntilde /ograve /oacute /ocircumflex /otilde /odieresis 
/divide /oslash /ugrave /uacute /ucircumflex /udieresis /yacute /thorn 
/ydieresis 
] def

/MF {   % fontname newfontname -> -     make a new encoded font
  /newfontname exch def
  /fontname exch def
  /fontdict fontname findfont def
  /newfont fontdict maxlength dict def
  fontdict { exch
    dup /FID eq { pop pop } 
      {  % copy to the new font dictionary
         exch newfont 3 1 roll put } ifelse
  } forall
  newfont /FontName newfontname put
  % insert only valid encoding vectors
  font_encoding_vector length 256 eq { newfont /Encoding font_encoding_vector put } if
  newfontname newfont definefont pop
} bind def

/Arial /font_Arial MF
%%EndProlog
%%Page: 1 1

/font_Arial findfont
10 scalefont
setfont
"""
	print
	print linewidth, "cm setlinewidth"
	print

	if landscape:
		print "90 rotate"
		print "0 -21 cm translate"
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
		s = "".join(c[6:]).lower()
		if s in colorListUNIX:
			continue
		colorTableRGB2UNIX[tuple([int(x) for x in c[3:6]])] = s
		print "/%s [%s] def" % (s, " ".join(c[:3]))
		colorListUNIX.add(s)

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
		colorTransl[str(-(i+1))] = ordre[i]
		if i < 26:
			colorTransl[chr(i+97)] = ordre[i]
		
	# Les couleurs foncees sont pour les nombres positifs et les lettres majuscules
	ordre = darkColors + lightColors + craniateColors
	#ordre = craniateColors + darkColors + lightColors
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
		return colorTransl[s].lower()

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

# Les coordonnes sont en cm (portrait/paysage)
#   X: de gauche (0) a droite (21/29.7)
#   Y: de bas (0) en haut (29.7/21)
# Les couleurs sont facultatives, utiliser None pour ne pas en specifier.
#   Sinon, utiliser le nom UNIX (cf getColor pour avoir ce nom)
# Les angles sont en degres

def setLineColor(C):
	if C in colorListUNIX:
		print C, "color"

def drawLine(X, Y, L, H, C):
	setLineColor(C)
	print "%.3f %.3f %.3f %.3f myline" % (L,H, X,Y)


def drawBox(X, Y, L, H, Cb, Cr):
	print "%.3f %.3f %.3f %.3f mybox" % (L,H, X,Y)

	if Cr in colorListUNIX:
		print Cr, "myfill"

	setLineColor(Cb)
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
	if C in colorListUNIX:
		print C, "color"
	print "(%s) %.3f %.3f mytext" % (T, X,Y)


