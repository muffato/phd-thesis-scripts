#! /users/ldog/muffato/python -OO

__doc__ = """
	Compare deux genomes:
		- Dessine la matrice des genes orthologues entre deux genomes.
		- Dessine le karyotype d'un genome face a l'autre
		- Renvoie les couples de chromosomes orthologues avec des stats sur l'evolution des genomes (nb genes, rearrangements)
		- Renvoie la liste des couples de genes orthologues avec les details sur leurs positions
		- Reordonne le genome 1 pour qu'il soit plus ressemblant au genome 2
		- Renvoie les diagonales entre les deux genomes
"""

##################
# INITIALISATION #
##################

# Librairies
import sys
import utils.myGenomes
import utils.myTools
import utils.myDiags
import utils.myPsOutput


#############
# FONCTIONS #
#############


def loadDCS(s):

	f = utils.myTools.myOpenFile(s, 'r')
	resL = []
	resD = {}
	curr = []
	for l in f:
		if l[0] == '-':
			if len(curr) > 0:
				resL.append(curr)
			curr = []
		else:
			(_,i,s,_,med,tet,zeb,fug,stk,anc) = l[:-1].split('\t')
			if anc == "None":
				continue
			x = (int(i),s,(med,tet,zeb,stk),anc)
			curr.append( x )
			resD[i] = x
	f.close()
	return (resL,resD)

#
# Charge le fichier qui contient les alternances que l'on doit observer
# Format de chaque ligne: "A	*Tetraodon.nigroviridis 2 3 5	*Gastero ..."
# Renvoie un dictionnaire des chromosomes ancestraux, qui contient les dictionnaires (par especes) des alternances
#
def loadChrAncIni(nom):

	chrAnc = {}
	f = utils.myTools.myOpenFile(nom, 'r')
	for ligne in f:

		c = ligne[:-1].split('\t')
		dic = phylTree.newCommonNamesMapperInstance()
		for x in c[1:]:
			(e,x) = x.split('|')
			(c1,c2) = x.split('/')

			s1 = set()
			for x in c1.split():
				try:
					x = int(x)
				except Exception:
					pass
				s1.add(x)

			s2 = set()
			for x in c2.split():
				try:
					x = int(x)
				except Exception:
					pass
				s2.add(x)

			dic[e] = (s1,s2)
		chrAnc[c[0]] = dic
	f.close()
	return chrAnc


def goto255(init, ratio):
	return int(init + ratio*(255-init))

def goto255c((r,g,b), ratio):
	return (goto255(r,ratio), goto255(g,ratio), goto255(b,ratio))


########
# MAIN #
########

# Arguments
modeMatrix = "Matrix"
modeKaryo = "Karyotype"
modeGenomeEvol = "GenomeEvolution"
modeOrthos = "OrthosGenes"
modeReindexedChr = "ReindexedChr"
modeDiags = "Diags"
modes = [modeMatrix, modeKaryo, modeOrthos, modeGenomeEvol, modeReindexedChr, modeDiags]
(noms_fichiers, options) = utils.myTools.checkArgs( \
	["phyltree", "predupGenome", "dcslist", "draftPreDupGenome.conf"], \
	[("includeGaps",bool,False), ("genesFile",str,"data44/genes/genes.%s.list.bz2"), \
	("output",str,modes), ("reverse",bool,False), ("scaleY",bool,False), ("minHomology",int,90), \
	("pointSize",float,-1), ("defaultColor",str,"black"), ("penColor",str,"black"), ("backgroundColor",str,"")], \
	__doc__
)


# Chargement des fichiers
phylTree = utils.myBioObjects.PhylogeneticTree(noms_fichiers["phyltree"])
#phylTree.loadSpeciesFromList(["Fugu"], options["genesFile"])
phylTree.loadAllSpeciesSince("Clupeocephala", options["genesFile"])
phylTree.loadSpeciesFromList(["Human"], options["genesFile"])
genomePreDup = utils.myGenomes.loadGenome(noms_fichiers["predupGenome"])
chrAnc = loadChrAncIni(noms_fichiers["draftPreDupGenome.conf"])
(dcsL,dcsD) = loadDCS(noms_fichiers["dcslist"])

	
print >> sys.stderr, "Affichage ...",

utils.myPsOutput.printPsHeader()

# On dessine
dx = 0.5
dy = 0.03
x0 = 1
y0 = 0.5
c = 'X'

#esp = [phylTree.officialName[x] for x in ["Medaka", "Tetraodon", "Zebrafish", "Fugu", "Stickleback"]]
esp = [phylTree.officialName[x] for x in ["Medaka", "Tetraodon", "Zebrafish", "Stickleback"]]

#gH = phylTree.dicGenomes["Human"]
#for i in xrange(len(gH.lstGenes["X"])):
#	s = gH.lstGenes["X"][i].names[0]
#	if s not in genomePreDup.dicGenes:
#		continue
#	(anc,_) = genomePreDup.dicGenes[s]
#	ancColor = utils.myPsOutput.getColor(str(anc), options["defaultColor"])
#	ancLightColor = utils.myPsOutput.getColor("#%d:%d:%d" % goto255c(utils.myPsOutput.colorTableUNIX2RGB[ancColor], 0.75), options["defaultColor"])
#	utils.myPsOutput.drawBox(x0, y0+i*dy, 11*dx, dy, ancLightColor, ancLightColor)


for dcs in dcsL:
	lastgood = None
	if len(dcs) == 1:
		continue
	for (i,s,poiss,anc) in dcs:
		ancColor = utils.myPsOutput.getColor(str(anc), options["defaultColor"])
		ancLightColor = utils.myPsOutput.getColor("#%d:%d:%d" % goto255c(utils.myPsOutput.colorTableUNIX2RGB[ancColor], 0.75), options["defaultColor"])

		utils.myPsOutput.drawBox(x0, y0+i*dy, 11*dx, dy, ancLightColor, ancLightColor)
		utils.myPsOutput.drawBox(x0+5*dx, y0+i*dy, dx, dy, ancColor, ancColor)

		for t in genomePreDup.getOtherNames(s):
			if t not in phylTree.dicGenes:
				continue
			(e,c,_) = phylTree.dicGenes[t]
			if e not in esp:
				continue
			#last = utils.myPsOutput.getColor(str(c), options["defaultColor"])
			j = esp.index(e)
			(s1,s2) = chrAnc[anc][e]
			if c in s1:
				off = 0
			elif c in s2:
				off = 1
			else:
				off = -1
			ppp = int((ord(anc)-65)*2+1.5*off)
			#last = utils.myPsOutput.getColor(str(ppp), options["defaultColor"])
			#if poiss[j] != "" and last != options["defaultColor"] and off != -1:
			if poiss[j] != "" and off != -1:
				#utils.myPsOutput.drawBox(xx + 3*(j+1) + off, y0+i*dy, dx/3., dy, last, last)
				#utils.myPsOutput.drawBox(xx +6*(off+1) + (j/3.), y0+i*dy, dx/3., dy, last, last)
				utils.myPsOutput.drawBox(x0 + (7*off + j)*dx, y0+i*dy, dx, dy, ancColor, ancColor)

utils.myPsOutput.printPsFooter()

print >> sys.stderr, "OK"

sys.exit(0)


table12 = buildOrthosTable(genomeH, genomeH.lstChr, genomeT, genomeT.lstChr)
#utils.myPsOutput.drawText(xx, y0, "X", options["penColor"])
y = y0 + 1

last = ""
nb = 0
res = []
for i in xrange(len(genomeH.lstGenes['X'])):
	continue
	tmp = table12[c][i]
	if len(tmp) == 0:
		col = options["defaultColor"]
	else:
		col = utils.myPsOutput.getColor(str(tmp[0][0]), options["defaultColor"])
	if col == last:
		nb += 1
	else:
		if nb > 0:
			res.append( (nb,last) )
		last = col
		nb = 1
if nb > 0:
	res.append( (nb,last) )

cote = 1
lastgood = None
for (nb,last) in res:
	continue
	if (nb > 1) and (last != options["defaultColor"]):
		if last != lastgood:
			cote *= -1
		utils.myPsOutput.drawBox(xx, y, dx*cote, nb/dy, last, last)
		lastgood = last
	y += nb/dy

dy = 0.025
for dcs in dcsL:
	lastgood = None
	for (i,s,(med,tet,zeb,fug,stk),anc) in dcs:
		#print >> sys.stderr, "*%s*%s*" % (i,s), 
		last = utils.myPsOutput.getColor(str(anc), options["defaultColor"])
		utils.myPsOutput.drawBox(xx, y0+i*dy, dx, dy, last, last)

utils.myPsOutput.printPsFooter()

print >> sys.stderr, "OK"

