# Extracts from the log file:
#    the number of fine steps
#    the id of the first, last fine step
#    the id of the first, last coarse step

/Main step=/ {
	split($0,m,"[a-z]+=")
	COARSE=m[2]+0
	if(!COARSE_INIT) { COARSE_START=COARSE; COARSE_INIT=1 }
}

/Fine step=/ {
	split($0,m,"[a-z]+=")
	FINE=m[2]+0
	if(!FINE_INIT) { FINE_START=FINE; FINE_INIT=1 }
}

END {
	COARSE_END = COARSE; FINE_END = FINE
	FCNT = FINE_END-FINE_START+1; CCNT = COARSE_END-COARSE_START+1

	print "\\newcommand{\\finestepfirst}{" FINE_START "}"      >> TEXFILE
	print "\\newcommand{\\finesteplast}{" FINE_END "}"         >> TEXFILE
	print "\\newcommand{\\finestepcount}{" FCNT "}"            >> TEXFILE
	print "\\newcommand{\\coarsestepfirst}{" COARSE_START "}"  >> TEXFILE
	print "\\newcommand{\\coarsesteplast}{" COARSE_END "}"     >> TEXFILE
	print "\\newcommand{\\coarsestepcount}{" CCNT "}"          >> TEXFILE
}
