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
	[], \
	__doc__
)

genome = utils.myGenomes.Genome(noms_fichiers["studiedGenome"])
paralogues = utils.myGenomes.Genome(noms_fichiers["paralogsList"])

arcs = set()
for g in paralogues:
	tg = genome.getPosition(g.names)
	arcs.update( utils.myTools.myIterator.tupleOnStrictUpperList(list(tg)) )

# Ecriture des arcs
####################
(fd,fname) = tempfile.mkstemp()
for (nb,((c1,i1),(c2,i2))) in enumerate(arcs):
	g1 = genome.lstGenes[c1][i1]
	g2 = genome.lstGenes[c2][i2]
	os.write(fd, "sd%d %s %d %d\n" % (nb, c1, g1.beginning, g1.end) )
	os.write(fd, "sd%d %s %d %d\n" % (nb, c2, g2.beginning, g2.end) )
os.close(fd)
options["datafile"] = fname

# Ecriture du karyotype
########################
(fd,fname) = tempfile.mkstemp()
for c in genome.lstChr:
	beginning = min([g.beginning for g in genome.lstGenes[c]])
	end = min([g.end for g in genome.lstGenes[c]])
	#os.write(fd, "chr%s %d %d band black\n" % (c, beginning, end) )
	os.write(fd, "chr - %s %s %d %d black\n" % (c, c, beginning, end) )
os.close(fd)
options["karyofile"] = fname

# Ecriture du fichier de config
################################


conf = """
<colors>
black = 0,0,0
white = 255,255,255
# <<include etc/colors.conf>>
</colors>

#<fonts>
#<<include etc/fonts.conf>>
#</fonts>


<ideogram>

<spacing>
 default = 5u
 break   = 1u
 axis_break_at_edge = yes
 axis_break         = yes
 axis_break_style   = 2
 <break_style 1>
  stroke_color = black
  fill_color   = blue
  thickness    = 0.25
  stroke_thickness = 2
 </break>
 <break_style 2>
  stroke_color     = black
  stroke_thickness = 3
  thickness        = 2
 </break>
</spacing>

# thickness (px) of chromosome ideogram
thickness        = 100p
stroke_thickness = 2
# ideogram border color
stroke_color     = black
fill             = yes
# the default chromosome color is set here and any value
# defined in the karyotype file overrides it
fill_color       = black

# fractional radius position of chromosome ideogram within image
radius         = 0.85r
show_label     = yes
label_with_tag = yes
label_font     = condensedbold
label_radius   = dims(ideogram,radius) + 0.075r
label_size     = 60p

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

chromosomes_units = 1000000
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

<link segdup>
show         = yes
color        = black
thickness    = 2
file         = %(datafile)s
#record_limit = 1000
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

conf2 = """

<colors>
black = 0,0,0
white = 255,255,255
</colors>

karyotype   = %(karyofile)s

outputdir = %(pngfiledir)s
outputfile = %(pngfilename)s

# radius of inscribed circle in image
radius         = 1000
# distance between chromosomes
chrspacing     = 500
# thickness (px) of chromosome ideogram
chrthickness   = 7
chrstroke      = 2
# ideogram border color
chrcolor       = black

chrlabel       = yes
chrlabelradius = 0.93
chrlabelsize   = 14

# fractional radius position of chromosome ideogram within image
chrradius      = 0.90

# cytogenetic bands
showbands       = yes
fillbands       = yes

<links duplicates>
  show         = yes
  color        = black
  thickness    = 1
  offset       = 0
  bezierradius = 0
  file         = %(datafile)s
  z            = 0
</links>



# chromosomes_scaling = 14:4,17:4

anglestep       = 1
minslicestep    = 5
beziersamples   = 20



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
os.system("/users/ldog/muffato/bin/circos-0.21 -conf %s 1>&2" % options["configfile"])
os.system("cat %s" % options["pngfile"])
os.system("rm %(pngfile)s %(karyofile)s %(datafile)s %(configfile)s" % options )

