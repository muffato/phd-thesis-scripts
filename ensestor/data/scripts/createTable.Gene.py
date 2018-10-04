#!/usr/bin/env python2

__doc__ = """
	Cree la base Genomicus (table Gene)
	Necessite la liste d'association Gene <-> gene_id definie par les arbres
	Lit les fichiers de genomes (modernes / diags+genes ancestraux)
"""

import sys
import collections

import utils.myFile
import utils.myTools
import utils.myGenomes
import utils.myPhylTree
import utils.myProteinTree



# Enregistre un genome ancestral
#################################
def storeAncGenes(anc):
	
	ancGenes = utils.myGenomes.Genome(arguments["ancGenes"] % phylTree.fileName[anc])
	genome = utils.myGenomes.Genome(arguments["diags"] % phylTree.fileName[anc], ancGenes=ancGenes)

	# Parcours des chromosomes
	print >> sys.stderr, "Inserting ancestral genome", anc, "...",
	for (chrom,l) in genome.lstGenes.iteritems():
		if len(l) == 1:
			gene = l[0]
			(gene_id,root_id) = dicGeneID[anc].pop(gene.names[0])
			res = [gene_id, phylTree.indNames[anc], gene.names[0], None, chrom, 0, None, None, root_id, None]

			print utils.myFile.myTSV.MySQLFileWriter(res)
		else:
			for (i,gene) in enumerate(l):
				(gene_id,root_id) = dicGeneID[anc].pop(gene.names[0])
				res = [gene_id, phylTree.indNames[anc], gene.names[0], arguments["ancBlockName"] % gene.chromosome, i, gene.strand, gene.beginning, gene.end, root_id, None]
				print utils.myFile.myTSV.MySQLFileWriter(res)
	
	print >> sys.stderr, "OK"


# Enregistre un genome moderne
################################
def storeModernGenome(esp):
	
	# Lecture du genome
	genome = utils.myGenomes.Genome(arguments["modernGenesFile"] % phylTree.fileName[esp])
	
	# Parcours des chromosomes
	print >> sys.stderr, "Inserting modern genome", esp, "...",
	for (chrom,l) in genome.lstGenes.iteritems():
		for (i,gene) in enumerate(l):
			(gene_id,root_id) = dicGeneID[esp].pop(gene.names[0])
			res = [gene_id, phylTree.indNames[esp], gene.names[0], gene.chromosome, i, gene.strand, gene.beginning, gene.end, root_id, None]
			print utils.myFile.myTSV.MySQLFileWriter(res)
	
	for (i,(name,(gene_id,root_id))) in enumerate(dicGeneID[esp].iteritems()):
		res = [gene_id, phylTree.indNames[esp], name, None, i, 0, None, None, root_id, None]
		print utils.myFile.myTSV.MySQLFileWriter(res)

	print >> sys.stderr, "OK (%d missing genes)" % len(dicGeneID[esp])




# Arguments
arguments = utils.myTools.checkArgs( \
	[("phylTree.conf",file), ("dicGeneID",file)], \
	[("modernGenesFile",str,""), ("ancGenes",str,""), ("diags",str,""), ("ancBlockName",str,"block_%s")], \
	__doc__ \
)

phylTree = utils.myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])

print >> sys.stderr, "Loading gene_id dictionary ...",
dicGeneID = collections.defaultdict(dict)
for (esp,name,gene_id,root_id) in utils.myFile.myTSV.readTabular(arguments["dicGeneID"], [str,str,int,int]):
	dicGeneID[esp][name] = (gene_id,root_id)
print >> sys.stderr, "OK"

for anc in phylTree.listAncestr:
	storeAncGenes(anc)
	assert len(dicGeneID[anc]) == 0, dicGeneID[anc]

for esp in phylTree.listSpecies:
	storeModernGenome(esp)
	#assert len(dicGeneID[esp]) == 0, dicGeneID[esp]

