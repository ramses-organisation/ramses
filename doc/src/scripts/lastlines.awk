/Mesh structure/ { START=FNR ; OUTPUT=$0 ; NLINES=1; next }

/Run completed/ {
	print OUTPUT "\n" $0
	NLINES++
	print "\\newcommand{\\lastlinesln}{" START "}"     >> TEXFILE
	print "\\newcommand{\\lastlinescount}{" NLINES "}" >> TEXFILE
	exit
}

{ OUTPUT = OUTPUT "\n" $0 ; NLINES++; next}
