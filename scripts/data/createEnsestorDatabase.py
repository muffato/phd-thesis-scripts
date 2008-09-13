#! /users/ldog/muffato/python

__doc__ = """
	Remplit une base de donnees a partir de fiohiers sources d'ensembl
"""

import sys
import sqlite3
import utils.myTools
import utils.myPhylTree
import utils.myProteinTree


def printIntoTable(tableName, values):

	if (tableName == "orthologues") and not arguments["printHomologs"]:
		return False

	# On modifie la table arbres pour ne garder que les sous-arbres de vertebres
	if tableName == "arbres":
		# Le noeud est trop ancien
		if not phylTree.isChildOf(values[6], arguments["cutoff"]):
			return False
		# Le pere est trop ancien, il faut renommer des champs
		if (values[-1] != None) and (not phylTree.isChildOf(values[7], arguments["cutoff"])):
			values = values[:3] + (None,None,values[5],values[6])
		else:
			values = values[:-1]
	print >> outputFiles[tableName], utils.myTools.printLine(values).replace("None", "\N")
	return True


# Enregistre l'arbre phylogenetique dans la base avec toutes les infos supplementaires
########################################################################################
def storePhylTree():

	dicEsp = {}

	# Puis les infos supplementaires d'assemblage et de genebuild
	for t in utils.myTools.readTabular(arguments["genome_db.txt"], (int,int,str,str,int,str,str)):
		dicEsp[t[2]] = (t[1],t[3],t[5])

	for esp in phylTree.allNames:
		if not phylTree.isChildOf(esp, arguments["cutoff"]):
			continue
		if esp in dicEsp:
			vers = "/".join(dicEsp[esp][1:])
		else:
			vers = ""
		printIntoTable("especes", (dicID[esp], esp, phylTree.ages[esp], vers, phylTree.fileName[esp].replace('.','_')))


# Lit chaque arbre de proteines et le stocke dans la base
# Definit la liste des genes de chaque espece (ancestrale/moderne)
###################################################################
def storeProteinTrees():

	def doProt(node, parent_name, parent_dup, parent_id, parent_distance):
		
		global gene_id, nextLeft

		dup = info[node]["Duplication"]

		# Les ancetres qui s'intercalent entre le precedent et le courant
		newAnc = info[node]['taxon_name']
		if parent_name != None:
			links = phylTree.dicLinks[parent_name][newAnc][:-1]
			if not parent_dup:
				links = links[1:]
			links.append(newAnc)
		else:
			links = [newAnc]

		# L'id du gene courant
		curr_id = gene_id
		gene_id += len(links)
		newparent_id = gene_id-1
	
		# La borne gauche evolue
		left0 = nextLeft
		nextLeft += len(links)

		ids = []
		allids = []
		if node in data:

			names = []
			# Appels recursifs sur les fils
			for (i,(fils,d)) in enumerate(data[node]):
				# Appel
				x = doProt(fils, newAnc, dup>=2, newparent_id, d)
				names.extend( x[0] )
				ids.append( x[1] )
				allids.extend( x[1] )

		else:
			# On garde le gene en memoire pour plus tard
			allGenes[newAnc].add( (newparent_id,info[node]['gene_name'],info[node].get('protein_name',None),info[node].get('transcript_name',None)) )
			names = [info[node]['gene_name']]

		# Il faut enregistrer les genes des ancetres intermediaires
		newid = []
		for (i,anc) in enumerate(links):
			
			if anc == newAnc:
				tmpdist = parent_distance
				tmpdup = dup
			else:
				tmpdist = 0
				tmpdup = 0
	
			# Il y a un gene uniquement si ce n'est pas une duplication
			if tmpdup < 2:
				# Les familles des especes modernes ne sont pas utiles
				# if anc not in phylTree.listSpecies:
				allNames[anc][frozenset(names)] = curr_id
			
			if printIntoTable("arbres", (curr_id, left0+i,nextLeft+len(links)-2*i-1, parent_id,tmpdist,tmpdup,dicID[anc],parent_name)):
				newid.append(curr_id)

			(curr_id, parent_id, parent_name) = (curr_id+1, curr_id, anc)
			nextLeft += 1
		assert parent_id == newparent_id

		# Ecriture des liens d'homologie
		if len(newid) != 0:

			# Duplication ou non ?
			if dup < 2:
				# Entre les fils
				for (l1,l2) in utils.myTools.myIterator.tupleOnStrictUpperList(ids):
					for x in l1:
						for y in l2:
							printIntoTable("orthologues", (x,y))
							printIntoTable("orthologues", (y,x))

			# Entre les peres et leurs enfants
			for x in allids:
				for y in newid:
					printIntoTable("orthologues", (x,y))
					printIntoTable("orthologues", (y,x))

			# Entre les differents peres
			for (x,y) in utils.myTools.myIterator.tupleOnStrictUpperList(newid):
				printIntoTable("orthologues", (x,y))
				printIntoTable("orthologues", (y,x))

		return (names,allids+newid)

	# Parcours du fichier d'arbres
	for (r,data,info) in utils.myProteinTree.loadTree(arguments["proteinTree"]):
		root_id = gene_id
		for g in doProt(r, None, None, None, None)[1]:
			dicRootID[g] = root_id
	print >> sys.stderr



# Enregistre un genome ancestral
#################################
def storeAncGenes(anc):
	print >> sys.stderr, "Inserting ancestral genes", anc, "...",
	dicOld = {}
	dicNames = {}
	for x in utils.myTools.readTabular(arguments["ancGenes"] % phylTree.fileName[anc], [int,str,str]):
		x = (x[0], x[1], frozenset(x[2].split()))
		# names -> old_id
		dicNames[x[2]] = x[0]
		# old_id -> (new_id,name,genes)
		dicOld[x[0]] = (allNames[anc][x[2]],x[1],x[2])
	
	if utils.myTools.fileAccess(arguments["ancGenomes"] % phylTree.fileName[anc]):
		lastChr = None
		genome = list(utils.myTools.readTabular(arguments["ancGenomes"] % phylTree.fileName[anc], [str,int,str,str]))
		genome.sort()
		for (chrom,strand,genes,strands) in genome:
			if chrom == "-":
				chrom = "Un"
			if chrom != lastChr:
				lastPos = 0
				lastChr = chrom
			lgenes = [int(x) for x in genes.split()]
			lstrands = [int(x) for x in strands.split()]
			if strand < 0:
				lgenes.reverse()
				lstrands = [-x for x in lstrands]
				lstrands.reverse()
			for (g,s) in zip(lgenes,lstrands):
				r = dicOld[g]
				printIntoTable("genes", (r[0], dicID[anc], r[1], chrom, lastPos, s, None, None, dicRootID[r[0]] ))
				lastPos += 1
				del dicNames[r[2]]
	else:
		genome = None
			
	if utils.myTools.fileAccess(arguments["blocsFile"] % phylTree.fileName[anc]):
		global synt_id
		nblock = 0
		for (_,genes,strands) in utils.myTools.readTabular(arguments["blocsFile"] % phylTree.fileName[anc], [int,str,str]):
			lgenes = [int(x) for x in genes.split()]
			if len(lgenes) == 1:
				continue
			lstrands = [int(x) for x in strands.split()]
			nblock += 1
			for (i,(g,s)) in enumerate(zip(lgenes,lstrands)):
				printIntoTable("syntenies", (synt_id,dicOld[g][0],s,i))
				if genome == None:
					printIntoTable("genes", (dicOld[g][0],dicID[anc],dicOld[g][1],"block_%d"%nblock,i,s,None,None,dicRootID[dicOld[g][0]]))
					del dicNames[dicOld[g][2]]
			synt_id += 1

	# Les genes qui ne sont pas sur des chromosomes
	for (i,(names,old_id)) in enumerate(dicNames.iteritems()):
		r = dicOld[old_id]
		printIntoTable("genes", (r[0], dicID[anc], r[1], None, i, 0, None, None, dicRootID[r[0]]))

	print >> sys.stderr, "OK"


# Enregistre un genome moderne
################################
def storeModernGenome(esp):
	print >> sys.stderr, "Inserting modern genome", esp, "...",
	global gene_id, nextLeft
	
	# Lecture du genome
	genes = utils.myTools.defaultdict(list)
	for (chrom,start,end,strand,name) in utils.myTools.readTabular(arguments["fullGenesFile"] % phylTree.fileName[esp], [str,int,int,int,str]):
		genes[chrom].append((start,end,strand,name))
	
	# Les genes deja vus dans les arbres
	dicGenesFromTree = {}
	for x in allGenes[esp]:
		dicGenesFromTree[x[1]] = (x[0],x[2],x[3])
	
	# Parcours des chromosomes
	for (chrom,l) in genes.iteritems():
		l.sort()
		i = 0
		# Et des genes
		for (start,end,strand,name) in l:

			pos = i
			i += 1
		
			# A-t-on des infos sur le transcrit
			if name in dicGenesFromTree:
				(curr_id,pr,tr) = dicGenesFromTree[name]
			else:
				# Les genes (proteines) sans arbre sont associes a des arbres singletons
				curr_id = gene_id
				gene_id += 1
				printIntoTable("arbres", (curr_id, nextLeft, nextLeft+1, None, None, 0, dicID[esp], None))
				nextLeft += 2
				dicRootID[curr_id] = curr_id
				
			printIntoTable("genes", (curr_id,dicID[esp],name,chrom,i,strand,start,end,dicRootID[curr_id]))

	print >> sys.stderr, "OK"





# Arguments
arguments = utils.myTools.checkArgs( \
	[("phylTree.conf",file), ("proteinTree",file), ("genome_db.txt",file)], \
	[("printHomologs",bool,True), \
	("fullGenesFile",str,""), ("ancGenes",str,""), ("ancGenomes",str,""), ("blocsFile",str,""), ("outputFile",str,""), ("cutoff",str,"Euteleostomi")], \
	__doc__ \
)

phylTree = utils.myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])

outputFiles = {}
for x in ["especes","genes","arbres","syntenies","orthologues"]:
	outputFiles[x] = utils.myTools.myOpenFile(arguments["outputFile"] % x, 'w')

dicID = phylTree.indNames
storePhylTree()

nextLeft = 0
gene_id = 0
allGenes = utils.myTools.defaultdict(set)
allNames = utils.myTools.defaultdict(dict)
dicRootID = {}
storeProteinTrees()

synt_id = 0
for anc in phylTree.listAncestr:
	if phylTree.isChildOf(anc, arguments["cutoff"]):
		storeAncGenes(anc)

for esp in phylTree.listSpecies:
	if phylTree.isChildOf(esp, arguments["cutoff"]):
		storeModernGenome(esp)

for f in outputFiles.itervalues():
	f.close()

