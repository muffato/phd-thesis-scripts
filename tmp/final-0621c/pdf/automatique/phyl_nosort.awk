BEGIN {
	while (getline < "/home/muff/M2/scripts/ps/rgb.txt" > 0) {
		a[$7$8$9] = $1" "$2" "$3;
	}
	print "%!PS-Adobe-3.0";
	print "%%DocumentData: Clean8bit";
	print "%%PageOrder: Ascend";
	print "%%Pages: 1";
	print "%%BoundingBox: 20 40 560 805"
	#print "%%BoundingBox: 15 40 600 805"
	print "%%DocumentFonts: Helvetica";
	print "%%EndComments";

	print "/cm {28.3464567 mul} def";
	
	print "%%Page: 1 1";
	#print "%%PageBoundingBox: 15 40 555 805"

	print "";
	print "/Arial findfont";
	print "6 scalefont";
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

	color[36] = "coral2"
	color[37] = "firebrick4"
	color[38] = "yellow4"
	color[39] = "LightSalmon"
	color[40] = "DarkSeaGreen4"
	color[41] = "turquoise2"
	color[42] = "purple4"
	color[43] = "magenta2"
	color[44] = "grey62"
	color[45] = "DarkViolet"
	color[46] = "OliveDrab2"
	color[47] = "MediumAquamarine"
	color[48] = "DarkSlateBlue"
	color[49] = "yellow"
	color[50] = "DarkGreen"
	color[51] = "DarkSeaGreen1"
	color[52] = "PaleTurquoise2"

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

	# épaisseur des lignes
	print "0.001 cm setlinewidth";

	# style de trait (coins ronds)
	print "1 setlinejoin";
 
}

$1~/>/{
	# Les parametres de la zone dans laquelle on se trouve
	x0 = $2;
	y0 = $3;
	wChr = $4;
	iChr = $5;
	hChr = $6;
	x = x0;
}


$1!~/>/{
	# on est dans une zone

	# texte
	text(x+0.1-(0.1*length($1)), y0-0.5, $1, "black")

	#barres
	y = y0
	for(i=2; i<=(NF-1); i=i+2){
		box(x, y, cm($(i+1)), wChr, color[substr($i,1,1)]);
		y = y + cm($(i+1));	  
	}
	x = x + wChr + iChr;

}

END{
	# Av -> Ahp
	rline(14.5, 3, -4.5, 0, "black", 0.05)
	rline(10, 3, 0, 3.5, "black", 0.05)
	# Av -> T
	rline(20, 3, 4.5, 0, "black", 0.05)
	rline(24.5, 3, 0, 8.5, "black", 0.05)
	arrowr(24.4, 11.5, 0.2, 0, 0.1, "black")
	# Ahp -> P
	rline(10, 6.5, 3.5, 0, "black", 0.05)
	rline(13.5, 6.5, 0, 2, "black", 0.05)
	arrowr(13.4, 8.5, 0.2, 0, 0.1, "black")
	# Ahp -> H
	rline(10, 6.5, -5.5, 0, "black", 0.05)
	rline(4.5, 6.5, 0, 4.5, "black", 0.05)
	arrowr(4.4, 11, 0.2, 0, 0.1, "black")
	
	print "/Arial-Bold-Italic findfont";
	print "12 scalefont";
	print "setfont";
	# txt T
	text(15, 26.5, "Tetraodon", "black")
	# txt T
	text(14, 14, "Poulet", "black")
	# txt T
	text(15, 7, "Homme", "black")
	# txt T
	text(1.9, 15.9, "Caryotype", "black")
	text(2, 15.4, "ancestral", "black")
	
	print "/Arial-Italic findfont";
	print "10 scalefont";
	print "setfont";
	# txt Av
	text(2.3, 14.9, "-450 Ma", "black")
	# txt Ahp
	text(5, 10.2, "-310 Ma", "black")
	print ""
	print "showpage"
}


function box(X,Y,H,L,C){
	print "newpath";
	print a[C], " setrgbcolor";
	print X, "cm",Y,"cm","moveto";
	print L, "cm", "0", "cm", "rlineto";
	print "0", "cm", H, "cm", "rlineto";
	print -L, "cm", "0", "cm", "rlineto";
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

function arrowr(Y,X,H,L,P,C){
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

function arrowl(Y,X,H,L,P,C){
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

function rline(Y1,X1,Y2,X2,C,W){
  print "";
  print "newpath";
  print C, "color";
  print W, "cm setlinewidth";
  print X1, "cm", Y1, "cm", "moveto";
  print X2, "cm", Y2, "cm", "rlineto";
  print "stroke";
  print "";
}


function cercle(Y,X,R,A,B,C){
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
