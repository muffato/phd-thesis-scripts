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

def loadDiagsGenes(s):
	
	f = utils.myTools.myOpenFile(s, 'r')
	res = set()
	for l in f:
		res.update(l.split())
	f.close()

	return res


def loadXBand(s1, s2):
	
	f = utils.myTools.myOpenFile(s2, 'r')
	res1 = {}
	for l in f:
		(b,c) = l.split()
		res1[b] = int(c)
	f.close()

	f = utils.myTools.myOpenFile(s1, 'r')
	res = {}
	for l in f:
		(g,b) = l.split()
		res[g] = b
	f.close()

	return (res1,res)

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
	["phyltree", "predupGenome", "dcslist", "draftPreDupGenome.conf", "Xband", "XbandColor", "HumanInDiags"], \
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
(bandC,bands) = loadXBand(noms_fichiers["Xband"], noms_fichiers["XbandColor"])
humanInDiags = loadDiagsGenes(noms_fichiers["HumanInDiags"])



	
print >> sys.stderr, "Affichage ...",

utils.myPsOutput.printPsHeader(landscape = True)



# On dessine
dx = 0.5
dy = 0.03
x0 = 1
y0 = 0.5
c = 'X'

#esp = [phylTree.officialName[x] for x in ["Medaka", "Tetraodon", "Zebrafish", "Fugu", "Stickleback"]]
esp = [phylTree.officialName[x] for x in ["Tetraodon", "Stickleback", "Medaka", "Zebrafish"]]
bandColors = [(0,0,0), (85,85,85), (170,170,170), (255,255,255)]
centroColor = utils.myPsOutput.getColor("#170:170:255", options["defaultColor"])
centroColor = "white"
assignedGenes = []

for dcs in dcsL:
	lastgood = None
	if len(dcs) == 1:
		continue
	for (i,s,poiss,anc) in dcs:
		assignedGenes.append(i)
		ancColor = utils.myPsOutput.getColor(str(anc), options["defaultColor"])
		ancLightColor = utils.myPsOutput.getColor("#%d:%d:%d" % goto255c(utils.myPsOutput.colorTableUNIX2RGB[ancColor], 0.75), options["defaultColor"])

		utils.myPsOutput.drawBox(y0+i*dy, x0, dy, 4*dx, ancLightColor, ancLightColor)
		utils.myPsOutput.drawBox(y0+i*dy, x0+7*dx, dy, 4*dx, ancLightColor, ancLightColor)


		for t in genomePreDup.getOtherNames(s):
			if t not in phylTree.dicGenes:
				continue
			(e,c,_) = phylTree.dicGenes[t]
			if e not in esp:
				continue
			j = esp.index(e)
			(s1,s2) = chrAnc[anc][e]
			if c in s1:
				off = 0
			elif c in s2:
				off = 1
			else:
				off = -1
			if poiss[j] != "" and off != -1:
				utils.myPsOutput.drawBox(y0+i*dy, x0 + (7*off + j) * dx, dy, dx, ancColor, ancColor)

gH = phylTree.dicGenomes["Human"]
lastBand = None
lastStart = 0
nb = 0
centromere = None
for i in xrange(len(gH.lstGenes["X"])):
	s = gH.lstGenes["X"][i].names[0]
	if bands[s][0] == 'q' and lastBand[0] == 'p':
		centromere = i
	if bands[s] != lastBand and lastBand != None:
		bandColor = utils.myPsOutput.getColor("#%d:%d:%d" % bandColors[bandC[lastBand]], options["defaultColor"])
		utils.myPsOutput.drawBox(y0+lastStart*dy, x0+5*dx, nb*dy, dx, "black", bandColor)
		lastStart = i
		nb = 0
	lastBand = bands[s]
	nb += 1
	if s in humanInDiags:
		utils.myPsOutput.drawBox(y0+i*dy, x0+5*dx, dy, -dy, "black", "black")
		utils.myPsOutput.drawBox(y0+i*dy, x0+6*dx, dy, dy, "black", "black")
		
bandColor = utils.myPsOutput.getColor("#%d:%d:%d" % bandColors[bandC[lastBand]], options["defaultColor"])
utils.myPsOutput.drawBox(y0+lastStart*dy, x0+5*dx, nb*dy, dx, "black", bandColor)
#utils.myPsOutput.drawCircle(y0+centromere*dy, x0+5*dx, 0.1, 0, 360, centroColor, centroColor)
#utils.myPsOutput.drawCircle(y0+centromere*dy, x0+6*dx, 0.1, 0, 360, centroColor, centroColor)
#utils.myPsOutput.drawBox(y0+centromere*dy-dy, x0+5*dx, 2*dy, dy, "black", centroColor)
#utils.myPsOutput.drawBox(y0+centromere*dy-dy, x0+6*dx, 2*dy, -dy, "black", centroColor)
utils.myPsOutput.drawLine(y0+centromere*dy, x0+4.5*dx, 0, 2*dx, "black")
fin = i+5
utils.myPsOutput.drawText(y0+fin*dy, x0+(5+.2)*dx, "HsaX", options["penColor"])
for i in xrange(len(esp)):
	(s1,s2) = esp[i].split()
	s = s1[0].upper() + s2[:3].lower()
	utils.myPsOutput.drawText(y0+fin*dy, x0+(i+0+.2)*dx, s, options["penColor"])
	utils.myPsOutput.drawText(y0+fin*dy, x0+(i+7+.2)*dx, s, options["penColor"])


assignedGenes.sort()
assignDic = {}
for i in xrange(len(assignedGenes)):
	assignDic[assignedGenes[i]] = i

x0 = 10
dy /= float(len(assignedGenes))/float(len(gH.lstGenes["X"]))

for dcs in dcsL:
	lastgood = None
	if len(dcs) == 1:
		continue
	for (ii,s,poiss,anc) in dcs:
		i = assignDic[ii]
		assignedGenes.append(i)
		ancColor = utils.myPsOutput.getColor(str(anc), options["defaultColor"])
		ancLightColor = utils.myPsOutput.getColor("#%d:%d:%d" % goto255c(utils.myPsOutput.colorTableUNIX2RGB[ancColor], 0.75), options["defaultColor"])

		utils.myPsOutput.drawBox(y0+i*dy, x0, dy, 4*dx, ancLightColor, ancLightColor)
		utils.myPsOutput.drawBox(y0+i*dy, x0+7*dx, dy, 4*dx, ancLightColor, ancLightColor)


		for t in genomePreDup.getOtherNames(s):
			if t not in phylTree.dicGenes:
				continue
			(e,c,_) = phylTree.dicGenes[t]
			if e not in esp:
				continue
			j = esp.index(e)
			(s1,s2) = chrAnc[anc][e]
			if c in s1:
				off = 0
			elif c in s2:
				off = 1
			else:
				off = -1
			if poiss[j] != "" and off != -1:
				utils.myPsOutput.drawBox(y0+i*dy, x0 + (7*off + j) * dx, dy, dx, ancColor, ancColor)

gH = phylTree.dicGenomes["Human"]
lastBand = None
lastStart = 0
nb = 0
centromere = None
for ii in xrange(len(gH.lstGenes["X"])):
	if ii not in assignDic:
		continue
	s = gH.lstGenes["X"][ii].names[0]
	i = assignDic[ii]
	if bands[s][0] == 'q' and lastBand[0] == 'p':
		centromere = i
	if bands[s] != lastBand and lastBand != None:
		bandColor = utils.myPsOutput.getColor("#%d:%d:%d" % bandColors[bandC[lastBand]], options["defaultColor"])
		#utils.myPsOutput.drawBox(x0+5*dx, y0+lastStart*dy, dx, nb*dy, "black", bandColor)
		utils.myPsOutput.drawBox(y0+lastStart*dy, x0+5*dx, nb*dy, dx, "black", bandColor)
		lastStart = i
		nb = 0
	lastBand = bands[s]
	nb += 1
	if s in humanInDiags:
		utils.myPsOutput.drawBox(y0+i*dy, x0+5*dx, dy, -dy, "black", "black")
		utils.myPsOutput.drawBox(y0+i*dy, x0+6*dx, dy, dy, "black", "black")
		#utils.myPsOutput.drawCircle(y0+(i+.5)*dy, x0+5*dx, 0.5*dy, 180, 360, options["penColor"], centroColor)
		#utils.myPsOutput.drawCircle(y0+(i+.5)*dy, x0+6*dx, 0.5*dy, 0, 180, options["penColor"], centroColor)

bandColor = utils.myPsOutput.getColor("#%d:%d:%d" % bandColors[bandC[lastBand]], options["defaultColor"])
utils.myPsOutput.drawBox(y0+lastStart*dy, x0+5*dx, nb*dy, dx, "black", bandColor)
#utils.myPsOutput.drawCircle(y0+centromere*dy, x0+5*dx, 0.1, 0, 180, "black", centroColor)
#utils.myPsOutput.drawCircle(y0+centromere*dy, x0+6*dx, 0.1, 180, 360, "black", centroColor)
#utils.myPsOutput.drawBox(y0+centromere*dy-dy, x0+5*dx, 2*dy, dy, "black", centroColor)
#utils.myPsOutput.drawBox(y0+centromere*dy-dy, x0+6*dx, 2*dy, -dy, "black", centroColor)
utils.myPsOutput.drawLine(y0+centromere*dy, x0+4.5*dx, 0, 2*dx, "black")
fin = i+5
utils.myPsOutput.drawText(y0+fin*dy, x0+(5+.2)*dx, "HsaX", options["penColor"])
for i in xrange(len(esp)):
	(s1,s2) = esp[i].split()
	s = s1[0].upper() + s2[:3].lower()
	utils.myPsOutput.drawText(y0+fin*dy, x0+(i+0+.2)*dx, s, options["penColor"])
	utils.myPsOutput.drawText(y0+fin*dy, x0+(i+7+.2)*dx, s, options["penColor"])




utils.myPsOutput.printPsFooter()

print >> sys.stderr, "OK"


