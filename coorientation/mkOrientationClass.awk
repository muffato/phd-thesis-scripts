BEGIN {
	s = ""
}

# Une ligne blanche = fin de la description de paires de genes conservees
NF == 0 {
	t = 1
	print s
	s = ""
}

# Une ligne avec des genes: $5 et $7 valent +1 ou -1 selon les orientations
NF != 0 {
	if ($5 == $7) {
		# Meme direction
		s = s" U"
	} else {
		if ($5 > 0) {
			# Convergence
			s = s" C"
		} else {
			# Divergence
			s = s" D"
		}
	}
}
