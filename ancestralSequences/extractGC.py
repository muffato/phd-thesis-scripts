#!/usr/bin/env python2
# Copyright (c) 2010 IBENS/Dyogen : Matthieu MUFFATO, Alexandra LOUIS, Hugues ROEST CROLLIUS
# Mail : hrc@ens.fr or alouis@biologie.ens.fr
# Licences GLP v3 and CeCILL v2

import sys

import utils.myFile
import utils.myMaths
import utils.myTools
import utils.myGenomes
import utils.myPhylTree

arguments = utils.myTools.checkArgs( [("fastaFile",file)], [], "Renvoie le GC des genomes modernes")

#fasta = utils.myGenomes.loadFastaFile(arguments["fastaFile"])
#for e in ["Danio.rerio_Oryzias.latipes_Gasterosteus.aculeatus_Tetraodon.nigroviridis_Takifugu.rubripes_Xenopus.tropicalis_Gallus.gallus_Ornithorhynchus.anatinus_Monodelphis.domestica_Loxodonta.africana_Echinops.telfairi_Dasypus.novemcinctus_Erinaceus.europaeus_Bos.taurus_Canis.familiaris_Macaca.mulatta_Pan.troglodytes_Homo.sapiens_Oryctolagus.cuniculus_Rattus.norvegicus_Mus.musculus"]:
#for e in fasta:
for l in sys.stdin:
	#seq = fasta[e]
	seq = l[:-1]
	seq3 = [b for (i,b) in enumerate(seq) if i%3 == 2]
	print (seq3.count("C") + seq3.count("G")) / (len(seq3)/100.)
	#print e, (seq3.count("C") + seq3.count("G")) / (len(seq3)/100.)
	#print e, (seq3.count("T")) / (len(seq3)/100.)


