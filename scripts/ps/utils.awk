function printPostScriptHeader() 
{
	while (getline < "/users/ldog/muffato/M2/scripts/ps/rgb.txt" > 0) {
		a[$7$8$9] = $1" "$2" "$3;
	}
	print "%!PS-Adobe-3.0";
	print "%%DocumentData: Clean8bit";
	print "%%PageOrder: Ascend";
	print "%%Pages: 1";
	print "%%DocumentFonts: Helvetica";
	print "%%EndComments";

	print "/cm {28.3464567 mul} def";
	print "%%Bounding-Box: 0 cm 0 cm 21 cm 29.7 cm";

	print "";
	print "/Arial findfont";
	print "10 scalefont";
	print "setfont";
	# style de trait (coins ronds)
	print "1 setlinejoin";

	# definition des couleurs dns l'en-tete du Postscript
	print "";
	print "/color {";
	print "\taload pop setrgbcolor";
	print "} def";
	print "";

	color[0] = "black"
	color[1] = "purple4"
	color[2] = "magenta2"
	color[3] = "grey62"
	color[4] = "coral2"
	color[5] = "firebrick4"
	color[6] = "yellow4"
	color[7] = "LightSalmon"
	color[8] = "DarkSeaGreen4"
	color[9] = "turquoise2"
	color[10] = "blue2"
	color[11] = "PeachPuff2"
	color[12] = "DarkViolet"
	color[13] = "OliveDrab2"
	color[14] = "MediumAquamarine"
	color[16] = "DarkSlateBlue"
	color[15] = "yellow"
	color[17] = "DarkGreen"
	color[18] = "gold"
	color[19] = "HotPink4"
	color[20] = "red1"
	color[21] = "orange"
	color[22] = "PaleTurquoise2"
	color[23] = "khaki1"
	color[24] = "DarkSeaGreen1"
	color[25] = "PaleTurquoise2"
	color[26] = "DarkOliveGreen1"
	color[27] = "khaki1"
	color[28] = "lavender"
	color[29] = "LightBlue"
	color[30] = "salmon"
	color[31] = "SkyBlue1"
	color[32] = "LightGoldenrod3"
	color[33] = "wheat1"
	color[34] = "thistle2"
	color[35] = "PeachPuff"

	color[41] = "coral2"
	color[51] = "firebrick4"
	color[61] = "yellow4"
	color[71] = "LightSalmon"
	color[81] = "DarkSeaGreen4"
	color[90] = "turquoise2"
	color[91] = "purple4"
	color[100] = "magenta2"
	color[101] = "grey62"
	color[120] = "DarkViolet"
	color[121] = "OliveDrab2"
	color[130] = "MediumAquamarine"
	color[131] = "DarkSlateBlue"
	color[140] = "yellow"
	color[141] = "DarkGreen"

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

	for (c in color) {
		couleur(color[c])
	}
	print "";

}

function getColor(s, d) {
	if (s in color)
		return color[s]
	else if (s~/#/) {
		r = substr(s, 2, 3)
		g = substr(s, 6, 3)
		b = substr(s, 10, 3)
		
		a["tmp"] = (r/255.)" "(g/255.)" "(b/255.);
		couleur("tmp")
		return "tmp"
	} else
		return d
}

function line(X, Y, L, H) {
	print "0 0 0 setrgbcolor";
	print "newpath";
	print X, "cm",Y,"cm","moveto";
	print L, "cm", H, "cm", "rlineto";
	print "closepath";
	print "stroke";
}

function box(X, Y, L, H, Cb, Cr) { 
	print "newpath";
	if (Cb in a)
		print a[Cb], " setrgbcolor";
	print X, "cm",Y,"cm","moveto";
	print L, "cm", "0", "cm", "rlineto";
	print "0", "cm", H, "cm", "rlineto";
	print -L, "cm", "0", "cm", "rlineto";
	print "closepath";
	if (Cr != 0) {
		print "gsave";
		print Cr, "color fill";
		print "grestore";
	}
	print "stroke";
}


function cross(X, Y, L, H) { 
	print "newpath";
	print X, "cm",Y,"cm","moveto";
	print L, "cm", H, "cm", "rlineto";
	print X+L, "cm",Y,"cm","moveto";
	print -L, "cm", H, "cm", "rlineto";
	print "closepath";
	print "stroke";
}

function cercle(X, Y, R, A, B) {
	print "";
	print "newpath";
	print X, "cm",Y,"cm", R, "cm", A, "cm", B, "cm", "arc";
	print "stroke";
}

function rond(X, Y, R, A, B, C) {
	print "";
	print "newpath";
	print X, "cm",Y,"cm", R, "cm", A, "cm", B, "cm", "arc";
	print "gsave";
	print C, "color fill";
	print "grestore";
	print "stroke";
}

function cmX(Z) {
	return dx*Z/nb;
}

function cmY(Z) {
	return dy*Z/nb;
}

function couleur(C) {
	print "/"C, "[", a[C], "]", "def";
}


function box(X,Y,L,H,C){
	print "newpath";
	print a[C], " setrgbcolor";
	print X, "cm",Y,"cm","moveto";
	print L, "cm", "0", "cm", "rlineto";
	print "0", "cm", H, "cm", "rlineto";
	print "-"L, "cm", "0", "cm", "rlineto";
	print "closepath";
	print "gsave";
	print C, "color fill";
	print "grestore";
	print "stroke";
}

# pour dessiner une fleche pointant a droite
# X,Y = moveto
# L,H taille du rectangle sans la pointe
# P= longueur de la pointe
# C= couleur

function arrowr(X,Y,L,H,P,C){
  print "newpath";
  print X, "cm",Y,"cm","moveto";
  print L, "cm", "0", "cm", "rlineto";
  print P, "cm", H/2, "cm", "rlineto";
  print "-"P, "cm", H/2, "cm", "rlineto";
  print "-"L, "cm", "0", "cm", "rlineto";
  print "closepath";
  print "gsave";
  print C, "color fill";
  print "grestore";
  print "stroke";
}

# pour dessiner une fleche pointant a gauche

function arrowl(X,Y,L,H,P,C){
  print "newpath";
  print X, "cm",Y+(H/2),"cm","moveto";
  print P, "cm", H/2, "cm", "rlineto";
  print L, "cm", "0", "cm", "rlineto";
  print 0, "cm", "-"H, "cm", "rlineto";
  print "-"L, "cm", "0", "cm", "rlineto";
  print "closepath";
  print "gsave";
  print C, "color fill";
  print "grestore";
  print "stroke";
}

# dessiner une ligne
# X1,Y1 coordonnées de depart
# X2,Y2 distance de deplacement (pas coordonnees d'arrivee)
# C= couleur
# W= epaisseur en cm

function rline(X1,Y1,X2,Y2,C,W){
  print "";
  print "newpath";
  print C, "color";
  print W, "cm setlinewidth";
  print X1, "cm", Y1, "cm", "moveto";
  print X2, "cm", Y2, "cm", "rlineto";
  print "stroke";
  print "";
}


function cercle(X,Y,R,A,B,C){
  print "";
  print "newpath";
  print X, "cm",Y,"cm", R, "cm", A, "cm", B, "cm", "arc";
  print "gsave";
  print C, "color fill";
  print "grestore";
  print "stroke";
}

function text(X,Y,T,C){
    print "";
    print "newpath";
    print C, "color";
    print X, "cm", Y, "cm", "moveto";
    print "(" T ")", "show";
    print "";
}


function cm(Z){
  return Z/hChr;
}


function couleur(C){
  print "/"C, "[", a[C], "]", "def";
}
