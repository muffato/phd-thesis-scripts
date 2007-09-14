BEGIN {
	t = 1
}

NF == 0 {
	t = 1
	print c[1],c[2],c[3],c[4]
}

NF != 0 {
	if ($5==$7) {
		s = "U"
	} else {
		if ($5>0) {
			s = "C"
		} else {
			s = "D"
		}
	}
	c[t] = s
	t = t+1
}
