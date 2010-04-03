#! /users/ldog/muffato/python

import sys
import collections

import utils.myFile
import utils.myTools
import utils.myGenomes
import utils.myPhylTree

arguments = utils.myTools.checkArgs( \
	[("genomesFile",file), ("dicGeneProtTrans",file), ("chromFile",file)], \
	[], \
	"Extrait les regions flanquantes des genes et calcule leur composition"
)


# Association gene_name > transcript_name
print >> sys.stderr, "Loading transcripts dict ...",
f = utils.myFile.openFile(arguments["dicGeneProtTrans"], "r")
dicProtTransc = {}
for l in f:
	t = l.split()
	dicProtTransc[t[0]] = t[2]
f.close()
print >> sys.stderr, "OK"

# Compte les nucleotides, les coordonnees sont entre 1 et n
def add(count, seq, x1, x2):
	#print "Calling", len(seq), x1, x2,
	if x1 < 1:
		x1 = 1
	if x1 > len(seq):
		#print 'None'
		return
	if x2 < 1:
		#print 'None'
		return
	if x2 > len(seq):
		x2 = len(seq)
	assert 1 <= x1 <= x2 <= len(seq)
	for x in seq[x1-1:x2]:
		count[x] += 1
	#print ">", x1, x2, "=", x2-x1+1

# Liste des genes
genome = utils.myGenomes.Genome(arguments["genomesFile"])

# Chargement du chromosome
for (chrom,seq) in utils.myGenomes.iterFastaFile(arguments["chromFile"]):
	chrom = utils.myGenomes.commonChrName(chrom.split()[0])
	print >> sys.stderr, chrom, "(len=%d)" % len(seq)
	#print chrom, len(seq), len(genome.lstGenes[chrom]), len([gene for gene in genome.lstGenes[chrom] if gene.names[0] in dicProtTransc])
	#continue

#for chrom in genome.lstGenes:
#	if not utils.myFile.hasAccess(arguments["chromFile"] % chrom):
#		continue
	
#	fasta = utils.myGenomes.loadFastaFile(arguments["chromFile"] % chrom)
#	assert len(fasta) == 1
#	seq = fasta.values()[0]

	# Parcours des genes, uniquement les transcrits dans les arbres
	for gene in genome.lstGenes[chrom]:
		if gene.names[0] not in dicProtTransc:
			continue
		#print gene

		#count = collections.defaultdict(int)
		#print gene.names[0], dicProtTransc[gene.names[0]], seq[gene.beginning-1:gene.end]
		#add(count, seq, gene.beginning, gene.beginning+2)
		#print dicProtTransc[gene.names[0]], count["A"], count["C"], count["G"], count["T"], count["N"]
		#continue
		
		# Les differentes regions avoisinantes
		count = collections.defaultdict(int)
		n = 0
		while n < 20:
			n += 1
			add(count, seq, gene.beginning-5000*n, gene.beginning-1-5000*(n-1))
			add(count, seq, gene.end+1+5000*(n-1), gene.end+5000*n)
			print dicProtTransc[gene.names[0]], "%dkb" % (5*n), count["A"], count["C"], count["G"], count["T"], count["N"], sum(count[x] for x in "ACGTN")
		
		if len(set("ACGTN").difference(count)) > 0:
			print >> sys.stderr, gene.names[0], dicProtTransc[gene.names[0]], count

