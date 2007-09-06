BEGIN {
	FS = "\t"
	OFS = "\t"
}

$1 == "groupI"	{$1 = 1}
$1 == "groupII"	{$1 = 2}
$1 == "groupIII"	{$1 = 3}
$1 == "groupIV"	{$1 = 4}
$1 == "groupV"	{$1 = 5}
$1 == "groupVI"	{$1 = 6}
$1 == "groupVII"	{$1 = 7}
$1 == "groupVIII"	{$1 = 8}
$1 == "groupIX"	{$1 = 9}
$1 == "groupX"	{$1 = 10}
$1 == "groupXI"	{$1 = 11}
$1 == "groupXII"	{$1 = 12}
$1 == "groupXIII"	{$1 = 13}
$1 == "groupXIV"	{$1 = 14}
$1 == "groupXV"	{$1 = 15}
$1 == "groupXVI"	{$1 = 16}
$1 == "groupXVII"	{$1 = 17}
$1 == "groupXVIII"	{$1 = 18}
$1 == "groupXIX"	{$1 = 19}
$1 == "groupXX"	{$1 = 20}
$1 == "groupXXI"	{$1 = 21}
$1 == "groupXXII"	{$1 = 22}
$1 == "groupXXIII"	{$1 = 23}
$1 == "groupXXIV"	{$1 = 24}
$1 == "groupXXV"	{$1 = 25}

{print}

END {
}
