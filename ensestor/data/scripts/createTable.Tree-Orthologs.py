#!/usr/bin/env python2

__doc__ = """
	Cree la base Genomicus (tables Tree & Orthologs)
	Necessite les genomes modernes pour rajouter les genes qui ne sont pas dans des arbres
	Necessite les genes ancestraux pour creer le dictionnaire gene_name <-> gene_id
"""

import sys
import collections

import utils.myFile
import utils.myTools
import utils.myGenomes
import utils.myPhylTree
import utils.myProteinTree



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
			allGenes[newAnc].add(info[node]['gene_name'])
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
				if anc in phylTree.listSpecies:
					assert len(names) == 1
					print >> outputFiles["dicGeneID"], utils.myFile.myTSV.MySQLFileWriter([anc, names[0], curr_id, root_id])
				else:
					print >> outputFiles["dicGeneID"], utils.myFile.myTSV.MySQLFileWriter([anc, dicAncGeneNames[anc][frozenset(names)], curr_id, root_id])
			print >> outputFiles["Tree"], utils.myFile.myTSV.MySQLFileWriter((curr_id, left0+i,nextLeft+len(links)-2*i-1, parent_id,tmpdist,tmpdup,phylTree.indNames[anc],root_id))
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
							print >> outputFiles["Orthologs"], utils.myFile.myTSV.MySQLFileWriter([x,y])
							print >> outputFiles["Orthologs"], utils.myFile.myTSV.MySQLFileWriter([y,x])

			# Entre les peres et leurs enfants
			for x in allids:
				for y in newid:
					print >> outputFiles["Orthologs"], utils.myFile.myTSV.MySQLFileWriter([x,y])
					print >> outputFiles["Orthologs"], utils.myFile.myTSV.MySQLFileWriter([y,x])

			# Entre les differents peres
			for (x,y) in utils.myTools.myIterator.tupleOnStrictUpperList(newid):
				print >> outputFiles["Orthologs"], utils.myFile.myTSV.MySQLFileWriter([x,y])
				print >> outputFiles["Orthologs"], utils.myFile.myTSV.MySQLFileWriter([y,x])

		return (names,allids+newid)

	# Parcours du fichier d'arbres
	for (r,data,info) in utils.myProteinTree.loadTree(arguments["proteinTree"]):
		root_id = gene_id
		doProt(r, None, None, None, None)
	print >> sys.stderr



# Enregistre un genome moderne
################################
def storeModernGenome(esp):
	print >> sys.stderr, "Inserting modern genome", esp, "...",
	global gene_id, nextLeft
	
	# Lecture du genome
	f = utils.myFile.openFile(arguments["modernGenesFile"] % phylTree.fileName[esp], "r")
	for l in f:
		name = l[:-1].split("\t")[4]
		if name not in allGenes[esp]:
			curr_id = gene_id
			gene_id += 1
			print >> outputFiles["Tree"], utils.myFile.myTSV.MySQLFileWriter((curr_id, nextLeft, nextLeft+1, None, None, 0, phylTree.indNames[esp], curr_id))
			print >> outputFiles["dicGeneID"], utils.myFile.myTSV.MySQLFileWriter((esp, name, curr_id, curr_id))
			nextLeft += 2
	f.close()
	print >> sys.stderr, "OK"



# Arguments
arguments = utils.myTools.checkArgs( \
	[("phylTree.conf",file), ("proteinTree",file)], \
	[("modernGenesFile",str,""), ("ancGenesFile",str,""), ("outputFile",str,"")], \
	__doc__ \
)

phylTree = utils.myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])

outputFiles = {}
for x in ["Tree", "Orthologs", "dicGeneID"]:
	outputFiles[x] = utils.myFile.openFile(arguments["outputFile"] % x, 'w')

dicAncGeneNames = {}
for anc in phylTree.listAncestr:
	dic = {}
	genome = utils.myGenomes.Genome(arguments["ancGenesFile"] % phylTree.fileName[anc])
	for gene in genome:
		dic[frozenset(gene.names[1:])] = gene.names[0]
	dicAncGeneNames[anc] = dic

nextLeft = 0
gene_id = 0
allGenes = collections.defaultdict(set)
storeProteinTrees()

for esp in phylTree.listSpecies:
	storeModernGenome(esp)

for f in outputFiles.itervalues():
	f.close()

