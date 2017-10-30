BEGIN { NLINES=1000 }

/Output *[0-9]+ *cells/ {
	NLINES=0
	START=FNR
	OUTPUT=$0
	sub("^ *","",$0)
	split($0,m," *"); NCELLS=m[2]
	next
}

(NLINES<7) { OUTPUT = OUTPUT "\n" $0 ; NLINES++; next}

END {
	print OUTPUT
	print "\\newcommand{\\lastoutputln}{" START "}"      >> TEXFILE
	print "\\newcommand{\\lastoutputncells}{" NCELLS "}" >> TEXFILE
}
