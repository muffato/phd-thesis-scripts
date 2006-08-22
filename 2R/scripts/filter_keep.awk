BEGIN {
	while (getline < "to_keep">0) {
		A[$1] = 1
	}
}

$1 in A && $3 in A {
	print
}
