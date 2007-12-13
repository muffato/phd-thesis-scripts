#! /users/ldog/muffato/python -OO

__doc__ = """
	Dessine le karyotype d'un genome face a l'autre
"""

# Librairies
import os
import sys
import tempfile

import itertools
import utils.myGenomes
import utils.myTools


# Arguments
(noms_fichiers, options) = utils.myTools.checkArgs( \
	["studiedGenome", "paralogsList"], \
	[("minNbParas",int,5)], \
	__doc__
)

genome = utils.myGenomes.Genome(noms_fichiers["studiedGenome"])
paralogues = utils.myGenomes.Genome(noms_fichiers["paralogsList"])


# Ecriture des paralogues
##########################
print >> sys.stderr, "Ecriture des liens ...",
(fd,fname) = tempfile.mkstemp()
nb = 0
nbPara = 0
para = utils.myTools.defaultdict(lambda : utils.myTools.defaultdict(int))
allLinks = []
for g in paralogues:
	tg = genome.getPosition(g.names)
	allLinks.append( genome.getPosition(g.names) )
	for ((c1,_),(c2,_)) in utils.myTools.myIterator.tupleOnStrictUpperList(list(tg)):
		para[c1][c2] += 1
		para[c2][c1] += 1
		nbPara += 2

lstChr = set()
for tg in allLinks:
	for ((c1,i1),(c2,i2)) in utils.myTools.myIterator.tupleOnStrictUpperList(list(tg)):
		if para[c1][c2] < options["minNbParas"]:
			continue
		lstChr.add(c1)
		lstChr.add(c2)
		#g1 = genome.lstGenes[c1][i1]
		#g2 = genome.lstGenes[c2][i2]
		#os.write(fd, "sd%d %s %d %d\n" % (nb, c1, g1.beginning, g1.end) )
		#os.write(fd, "sd%d %s %d %d\n" % (nb, c2, g2.beginning, g2.end) )
		os.write(fd, "sd%d %s %d %d\n" % (nb, c1, i1, i1) )
		os.write(fd, "sd%d %s %d %d\n" % (nb, c2, i2, i2) )
		nb += 1
os.close(fd)
options["datafile"] = fname
print >> sys.stderr, nb, nbPara/2, "OK"

# Ecriture du karyotype
########################
(fd,fname) = tempfile.mkstemp()
for c in lstChr:
	beginning = 0
	end = len(genome.lstGenes[c])
	#beginning = min([g.beginning for g in genome.lstGenes[c]])
	#end = max([g.end for g in genome.lstGenes[c]])+1
	#os.write(fd, "chr%s %d %d band black\n" % (c, beginning, end) )
	os.write(fd, "chr - %s %s %d %d black\n" % (c, c, beginning, end) )
os.close(fd)
options["karyofile"] = fname

# Ecriture du fichier de config
################################
conf = """
<colors>
 <<include /users/ldog/muffato/bin/circos-0.32/etc/colors.conf>>
</colors>

<fonts>
 <<include /users/ldog/muffato/bin/circos-0.32/etc/fonts.conf>>
</fonts>

<ideogram>
 <spacing>
  default = 50u
  break   = 1u
  axis_break_at_edge = no
  axis_break         = no
  axis_break_style   = 2
  <break_style 2>
   stroke_color     = black
   stroke_thickness = 3
   thickness        = 2
  </break>
 </spacing>
 
 # thickness (px) of chromosome ideogram
 thickness        = 5p
 stroke_thickness = 2
 # ideogram border color
 stroke_color     = black
 fill             = yes
 # the default chromosome color is set here and any value
 # defined in the karyotype file overrides it
 fill_color       = black
 
 # fractional radius position of chromosome ideogram within image
 radius         = 0.95r
 show_label     = yes
 label_with_tag = yes
 label_font     = condensedbold
 #label_radius   = dims(ideogram,radius) + 0.075r
 label_radius   = 0.975r
 label_size     = 20p
 
 # cytogenetic bands
 band_stroke_thickness = 2
 
 # show_bands determines whether the outline of cytogenetic bands
 # will be seen
 show_bands            = yes
 # in order to fill the bands with the color defined in the karyotype
 # file you must set fill_bands
 fill_bands            = yes

</ideogram>


karyotype   = %(karyofile)s

<image>
 dir = %(pngfiledir)s
 file = %(pngfilename)s

 # radius of inscribed circle in image
 radius         = 500p
 background     = white
 # by default angle=0 is at 3 o'clock position
 angle_offset   = -90
</image>

chromosomes_units = 1
chromosomes_display_default = yes

# Links (bezier curves or straight lines) are defined in <links> blocks.
# Each link data set is defined within a <link>.
# As with highlights, parameters defined
# in the root of <links> affect all data sets and are considered
# global settings. Individual parameters value can be refined by
# values defined within <link> blocks, or additionally on each
# data line within the input file.

<links>
 z      = 0
 radius = 0.95r
 bezier_radius = 0.2r
 #ribbon = yes
 <link segdup>
  show         = yes
  color        = grey
  thickness    = 1
  file         = %(datafile)s
  #record_limit = 100
 </link>
</links>

anglestep       = 0.5
minslicestep    = 10
beziersamples   = 40
debug           = no
warnings        = no
imagemap        = no

# don't touch!
units_ok        = bupr
units_nounit    = n

"""

(fd,fname) = tempfile.mkstemp()
os.close(fd)
options["pngfile"] = fname
options["pngfiledir"] = os.path.dirname(fname)
options["pngfilename"] = os.path.basename(fname)
(fd,fname) = tempfile.mkstemp()
options["configfile"] = fname
os.write(fd, conf % options)
os.close(fd)
os.system("/users/ldog/muffato/bin/circos -conf %s 1>&2" % options["configfile"])
os.system("cat %s" % options["pngfile"])
os.system("rm %(pngfile)s %(karyofile)s %(datafile)s %(configfile)s" % options )

