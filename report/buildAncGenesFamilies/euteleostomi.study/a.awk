BEGIN{FS="\t"}
{
s=0
for(i=1; i<=16; i++) {
	if ($(i)>0) {s+=1}
}
if (s==1) {
print $0
}
}
