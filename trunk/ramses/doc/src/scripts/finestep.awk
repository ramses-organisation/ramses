# Strip leading line numbers and store command
(FNR==1) {
	split($0,a,"[:-]")
	print "\\newcommand{\\finestepln}{" a[1] "}" >> TEXFILE
}

{
	sub("^[0-9]+[:-]","",$0)
	print $0
}
