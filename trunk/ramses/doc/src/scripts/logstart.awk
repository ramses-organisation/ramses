BEGIN             { SINCE_OUTPUT=1000 }
                  { print ; SINCE_OUTPUT=SINCE_OUTPUT+1 }
(SINCE_OUTPUT==4) { exit}
/Output .* cells/ { SINCE_OUTPUT=0 }
