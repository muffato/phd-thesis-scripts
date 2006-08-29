BEGIN {
	while (getline < "to_delete">0) {
		A[$2$3] = 1
	}
}

!($1$3 in A) && !($3$1 in A) {
	print
}
