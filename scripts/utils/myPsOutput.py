#! /users/ldog/muffato/python -OO

# Module d'ecriture dans un fichier PostScript

# Les petits noms que je donne a mes couleurs (nom -> nom UNIX)
colorTransl = {}
# Les definitions des couleurs (valeurs RGB -> nom UNIX) pour eviter de creer plusieurs fois la meme couleur
colorTableRGB2UNIX = {}
# La table inverse
colorTableUNIX2RGB = {}


#
# L'en-tete PostScript
#
def printPsHeader(landscape = False):
	# En-tete Postscript
	print """%!PS-Adobe-3.0
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
	print "0.001 cm setlinewidth"
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
		if s in colorTableUNIX2RGB:
			continue
		rgb = tuple([int(x) for x in c[3:6]])
		print "/%s [%s] def" % (s, " ".join(c[:3]))
		colorTableRGB2UNIX[rgb] = s
		colorTableUNIX2RGB[s] = rgb

	f.close()


	c5 = []
	c5 += ["red4", "coral4", "firebrick", "red2", "coral2", "darkorange", "gold"]
	c5 += ["yellow", "khaki1","wheat1","peachpuff","lightsalmon", "hotpink2", "magenta2", "darkorchid2", "purple2", "darkorchid4"]
	c5 += ["blue2", "royalblue2", "blue4"]
	c5 += ["turquoise4", "darkseagreen4", "chartreuse4","mediumaquamarine",(121,204,61), "chartreuse2", "olivedrab2", "darkolivegreen1"]
	c5 += ["darkseagreen1", "paleturquoise2", "lightblue", "skyblue1", "turquoise2", "lavender","thistle2"]
	c5 += [(204,204,153),"lightgoldenrod3","ivory2","honeydew3","slategray"]
		
	lightColors = ["salmon", "PaleTurquoise2", "DarkSeaGreen1","khaki1", "thistle2", "PeachPuff","LightBlue", "SkyBlue1","LightGoldenrod3", "wheat1",  "DarkOliveGreen1", "lavender"]
	
	darkColors = ["red1", "turquoise2", "DarkGreen", "yellow", "coral2", "OliveDrab2", "orange", "MediumAquamarine", "blue2", "firebrick4", "LightSalmon", "DarkViolet", "magenta2", "DarkSeaGreen4", "DarkSlateBlue", "yellow4", "grey62", "gold", "PeachPuff2", "HotPink4", "firebrick", "purple4"]
	darkColors = ["red3", "chartreuse4", "turquoise2", "gold", "coral2", "OliveDrab2", "orange", "DarkSeaGreen4", "firebrick4", "RoyalBlue4", "LightSalmon", "DarkViolet", "magenta2", "MediumAquamarine"]
	#darkColors = ["turquoise2", "gray85", "yellow", "coral2", "OliveDrab2", "MediumAquamarine", "blue2", "LightSalmon", "magenta2", "DarkSeaGreen2", "PaleTurquoise3", "khaki2", "thistle3", "PeachPuff2", "HotPink2", "firebrick", "purple2"]

	greekLetters = ["ALPHA", "BETA", "DELTA", "EPSILON", "GAMMA", "PHI"]
	craniateColors = ["DarkOrange", "RoyalBlue4", "chartreuse4", "gold", "DarkOrchid4", "red3"]
	#craniateColors = ["DarkOrange", "RoyalBlue2", "chartreuse2", "gold", "DarkOrchid1", "red2"]

	# Les couleurs claires sont pour les nombres negatifs et les lettres minuscules
	ordre = lightColors + craniateColors + darkColors
	for i in xrange(len(ordre)):
		colorTransl[str(-(i+1))] = ordre[i].lower()
		if i < 26:
			colorTransl[chr(i+97)] = ordre[i].lower()
		
	# Les couleurs foncees sont pour les nombres positifs et les lettres majuscules
	ordre = darkColors + lightColors + craniateColors
	#ordre = craniateColors + darkColors + lightColors
	for i in xrange(len(ordre)):
		colorTransl[str(i+1)] = ordre[i].lower()
		if i < 26:
			colorTransl[chr(i+65)] = ordre[i].lower()

	for i in xrange(len(greekLetters)):
		colorTransl[greekLetters[i]] = craniateColors[i]


#######################
# Fonctions de dessin #
#######################

# Les coordonnes sont en cm (portrait/paysage)
#   X: de gauche (0) a droite (21/29.7)
#   Y: de bas (0) en haut (29.7/21)
# Les couleurs sont facultatives, utiliser None pour ne pas en specifier.
#   Sinon, utiliser le nom UNIX (cf getColor pour avoir ce nom)
# Les angles sont en degres

def setColor(C, txt):

	# Un nom comme on les aime
	if C in colorTableUNIX2RGB:
		s = C
	# Des raccourcis
	elif str(C) in colorTransl:
		s = colorTransl[str(C)]
	# Un triplet (r,g,b)
	else:
		if len(C) == 3:
			(r,g,b) = C
		elif C[0] == '#':
			(r,g,b) = [int(x) for x in C[1:].split(':')]
		else:
			return

		if (r,g,b) in colorTableRGB2UNIX:
			s = colorTableRGB2UNIX[(r,g,b)]
		else:
			s = "tmp%d" % len(colorTableUNIX2RGB)
			colorTableRGB2UNIX[(r,g,b)] = s
			colorTableUNIX2RGB[s] = (r,g,b)
			print "/%s [%f %f %f] def" % (s, float(r)/255., float(g)/255., float(b)/255.)
	print s, txt


def drawLine(X, Y, L, H, C):
	setColor(C, "color")
	print "%.5f %.5f %.5f %.5f myline" % (L,H, X,Y)


def drawBox(X, Y, L, H, Cb, Cr):
	print "%.5f %.5f %.5f %.5f mybox" % (L,H, X,Y)
	setColor(Cr, "myfill")
	setColor(Cb, "color")
	print "stroke"


def drawCross(X, Y, L, H, C):
	print "newpath"
	setColor(C, "color")
	print "%.5f cm %.5f cm moveto" % (X, Y)
	print "%.5f cm %.5f cm rlineto" % (L, H)
	print "%.5f cm %.5f cm moveto" % (X+L, Y)
	print "%.5f cm %.5f cm rlineto" % (-L, H)
	print "closepath"
	print "stroke"

def drawCircle(X, Y, R, A, B, Cb, Cr):
	print "newpath"
	setColor(Cb, "color")
	print "%.5f cm %.5f cm %.5f cm %.5f %.5f arc" % (X, Y, R, A, B)
	setColor(Cr, "myfill")
	print "stroke"

def drawArrowR(X, Y, L, H, P, Cb, Cr):
	print "newpath"
	setColor(Cb, "color")
	print "%.5f cm %.5f cm moveto" % (X, Y)
	print "%.5f cm 0 cm rlineto" % L
	print "%.5f cm %.5f cm rlineto" % (P, H/2)
	print "%.5f cm %.5f cm rlineto" % (-P, H/2)
	print "%.5f cm 0 cm rlineto" % (-L)
	print "closepath"
	setColor(Cr, "myfill")
	print "stroke"

def drawArrowL(X, Y, L, H, P, Cb, Cr):
	print "newpath";
	setColor(Cb, "color")
	print "%.5f cm %.5f cm moveto" % (X, Y+(H/2))
	print "%.5f cm %.5f cm rlineto" % (P, H/2)
	print "%.5f cm 0 cm rlineto" % L
	print "0 cm %.5f cm rlineto" % (-H)
	print "%.5f cm 0 cm rlineto" % (-L)
	print "closepath"
	setColor(Cr, "myfill")
	print "stroke"

def drawText(X, Y, T, C):
	setColor(C, "color")
	print "(%s) %.5f %.5f mytext" % (T, X,Y)


