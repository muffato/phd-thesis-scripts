#! /users/ldog/muffato/python -OO

__doc__ = """
	Remplit une base de donnees a partir de fiohiers sources d'ensembl
"""

import os
import sys
import time
import random
import sqlite3
import utils.myTools
import utils.myPhylTree

t0 = time.time()

# Arguments
arguments = utils.myTools.checkArgs( [("database.db",file), ("phylTree.conf",file), ("nbIter",int)], [], __doc__)

phylTree = utils.myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])
conn = sqlite3.connect(arguments["database.db"])
cursor = conn.cursor()

for i in xrange(arguments["nbIter"]):
	(e1,e2) = random.sample(phylTree.listSpecies, 2)
	a = phylTree.dicParents[e1][e2]

	cursor.execute('select count(*) from diags_proj_summary inner join chromosomes chr1 on diags_proj_summary.chrom1_id=chr1.chrom_id inner join chromosomes chr2 on diags_proj_summary.chrom2_id=chr2.chrom_id where chr1.taxon_id=%d and chr2.taxon_id=%d and diags_proj_summary.taxon_id=%d;' % (phylTree.indNames[e1],phylTree.indNames[e2],phylTree.indNames[a]))

	print e1, e2, a, cursor.fetchall()

print time.time()-t0
print os.times()
