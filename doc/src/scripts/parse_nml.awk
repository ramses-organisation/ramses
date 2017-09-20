# Extracts from the namelist file:
#    levelmin
#    levelmax

BEGIN {
	LEVELMIN="#ERR#" 
	LEVELMAX="#ERR#" 
}

/^[ \t]*levelmin[ \t]*=/ {
	split($0,m,"=") ; LEVELMIN=m[2]+0
}


/^[ \t]*levelmax[ \t]*=/ {
	split($0,m,"=") ; LEVELMAX=m[2]+0
}

END {
	print "\\newcommand{\\nmllevelmin}{" LEVELMIN "}" >> TEXFILE
	print "\\newcommand{\\nmllevelmax}{" LEVELMAX "}" >> TEXFILE
}
