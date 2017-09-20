#!/bin/bash

############################################################################
# Ramses Daily Build Script
# 
# Neil Vaytet (ENS Lyon) - 07/2014 - neil.vaytet@ens-lyon.fr
#
# This script runs the RAMSES test suite at regular time intervals.
# It is currently configured to perform the test suite every day at 3am and
# upload the results of the tests to the bitbucket wiki, using git commits.
# If you just want to have the results in a local folder, without
# interacting with a wiki page, simply set the variable UPDATEWIKI to false.
############################################################################

# The source directory:
SRC="/home/ubuntu";

# Test frequency: (YY:MM:DD:hh:mm:ss)
TEST_FREQ="00:00:01:00:00:00";

# Test time offset: (YY:MM:DD:hh:mm:ss)
TEST_OFFS="00:00:00:03:00:00";

# Wiki file
WIKIFILE="AutoTests.md";

# Log file
LOGFILE="daily_test.log";

# Upload to wiki?
UPDATEWIKI=true;

# Set up variables
hline="============================================================";

RAMSESDIR="${SRC}/ramses/trunk/ramses/patch/test_suite";
WIKIDIR="${SRC}/ramses/wiki";
WIKISTOREDIR="daily_tests";
LOGFILE="${SRC}/${LOGFILE}";
COMMIT_URL="https://bitbucket.org/rteyssie/ramses/commits/";

pause=100;

############################################################################

echo "Starting RAMSES daily build script";

# Find next test time:

# Decompose frequencies and offsets into arrays
fbackupstring=$(echo $TEST_FREQ | sed 's/:/ /g');
frequencies=( $fbackupstring );
obackupstring=$(echo $TEST_OFFS | sed 's/:/ /g');
offsets=( $obackupstring );

# Find lowest time denominator
for ((i=0;i<6;i++)); do
    if [ ${frequencies[${i}]} -gt 0 ] ; then
        unit=$i;
        break;
    fi
done

# Determine reference time:
# If time now is 15:15 and a backup is performed every hour, reference time is 15:00
TIMENOW=$(date "+%Y %m %d %H %M %S");
REFTIME=( $TIMENOW );
for ((i=5;i>$unit;i--)); do
    REFTIME[${i}]="00";
done

# Remove leading zeros for arithmetic operations
for ((i=0;i<6;i++)); do
   if [ "${REFTIME[i]:0:1}" == "0" ] ; then
      REFTIME[i]="${REFTIME[i]:1:1}";
   fi
   if [ "${offsets[i]:0:1}" == "0" ] ; then
      offsets[i]="${offsets[i]:1:1}";
   fi
done

# Offset current time to find next backup time
for ((i=0;i<6;i++)); do
    OFFTIME[i]=$((${REFTIME[i]}+${offsets[i]}));
done
TNEXTBACKUP=$(date -d "${OFFTIME[0]}-${OFFTIME[1]}-${OFFTIME[2]} ${OFFTIME[3]}:${OFFTIME[4]}:${OFFTIME[5]}" +%s)

# Determine whether we are to backup in current time interval or the next
TIMENOW=$(date +%s);
if [ $TIMENOW -gt $TNEXTBACKUP ] ; then
    OFFTIME[${unit}]=$((${OFFTIME[${unit}]}+1));
    TNEXTBACKUP=$(date -d "${OFFTIME[0]}-${OFFTIME[1]}-${OFFTIME[2]} ${OFFTIME[3]}:${OFFTIME[4]}:${OFFTIME[5]}" +%s);
fi

TNEXT=$(date -d "${OFFTIME[0]}-${OFFTIME[1]}-${OFFTIME[2]} ${OFFTIME[3]}:${OFFTIME[4]}:${OFFTIME[5]}");
echo "Will perform next backup $TNEXT";

############################################################################
# Define some variables
months[0]="January";
months[1]="February";
months[2]="March";
months[3]="April";
months[4]="May";
months[5]="June";
months[6]="July";
months[7]="August";
months[8]="September";
months[9]="October";
months[10]="November";
months[11]="December";

hline="| ---:|:---:|:---:|:---:|:---:|";
############################################################################

# Now begin infinite loop which will create tests
while true ; do

    # Check whether it is time to perform test
    TIMENOW=$(date "+%s");

    if [ $TIMENOW -ge $TNEXTBACKUP ]; then

        echo "Performing test run at $(date):";

        echo "Performing test run at $(date):" >> $LOGFILE;

        # Todays date for file name:
        DATE=$(date "+%Y-%m-%d");
        YEAR=$(date "+%Y");
        YEARMONTH=$(date "+%Y-%m");
        MONTH=$(date "+%m");
        DAY=$(date "+%d");
        MONTHNOZERO=$MONTH;
        DAYNOZERO=$DAY;
        if [ "${MONTHNOZERO:0:1}" == "0" ] ; then
           MONTHNOZERO="${MONTHNOZERO:1:1}";
        fi
        if [ "${DAYNOZERO:0:1}" == "0" ] ; then
           DAYNOZERO="${DAYNOZERO:1:1}";
        fi

        # Go to RAMSES directory and make fresh start
        cd ${RAMSESDIR};
        rm test_results.pdf test_suite.log;

        # Get latest changes
        git pull >> $LOGFILE;

        # Save commit hash
        commit=$(git rev-parse HEAD);

        # Run ramses test suite
        ./run_test_suite.sh >> $LOGFILE;
        
        # Go to wiki directory
        cd "${WIKIDIR}/${WIKISTOREDIR}";
        cp ${RAMSESDIR}/test_results.pdf ${DATE}.pdf;
        cp ${RAMSESDIR}/test_suite.log ${DATE}.log;

        if $UPDATEWIKI ; then
        
            # Update wiki page ================================================
            cd ${WIKIDIR};
            
            # Pull latest wiki ================================================
	    git pull;

            # Number of lines in file
            nlines=$(wc -l ${WIKIFILE} | cut -d ' ' -f 1);
            
            # Check if year already exists
            yearcheck=$(grep "## ${YEAR}" ${WIKIFILE});
            if [ ${#yearcheck} -eq 0 ] ; then
                oldyear=$(($YEAR - 1));
                sed -i "s/## ${oldyear}/## ${YEAR}\n\n## ${oldyear}/g" ${WIKIFILE};
                # Update number of lines
                nlines=$(wc -l ${WIKIFILE} | cut -d ' ' -f 1);
            fi
            
            # Check if month already exists
            monthcheck=$(grep "${YEARMONTH}" ${WIKIFILE});
            if [ ${#monthcheck} -eq 0 ] ; then
                # Generate new month group
                m=$((($MONTHNOZERO - 1) / 4));
                monthlist="||";
                monthhead="[//]: # (";
                for ((i=0;i<4;i++)); do
                   j=$(($m * 4 + $i));
                   monthlist="${monthlist} ${months[$j]} |";
                   imonth=$(($j + 1));
                   if [ $imonth -lt 10 ] ; then
                      imonth="0${imonth}";
                   fi
                   monthhead="${monthhead} ${YEAR}-${imonth} ";
                done
                monthhead="${monthhead})";
                
                lnumber=$(grep -n "## ${YEAR}" ${WIKIFILE} | head -1 | cut -d ':' -f 1);
                fileend=$(($nlines - $lnumber));
                
                head -n $lnumber ${WIKIFILE} > tempfile;
                echo >> tempfile;
                echo "$monthhead" >> tempfile;
                echo "$monthlist" >> tempfile;
                echo "$hline" >> tempfile;
                for ((i=0;i<31;i++)); do
                    j=$(($i + 1));
                    echo "| ${j} | | | | |" >> tempfile;
                done
                echo >> tempfile;
                tail -n ${fileend} ${WIKIFILE} >> tempfile;
                cp tempfile ${WIKIFILE};
                # Update number of lines
                nlines=$(wc -l ${WIKIFILE} | cut -d ' ' -f 1);
            fi
            
            # Identify line where to insert new results
            lnumber=$(grep -n ${YEARMONTH} ${WIKIFILE} | head -1 | cut -d ':' -f 1);
            
            line1=$(($lnumber + $DAYNOZERO + 2));
            line2=$(($line1 - 1));
            line3=$(($nlines - $line1));
            
            theline=$(head -n $line1 ${WIKIFILE} | tail -1);
            
            # Extract separate columns from line
            for ((i=0;i<5;i++)); do
               j=$(($i + 2));
               segment[i]=$(echo $theline | cut -d '|' -f $j);
            done
            
            # Generate image and commit string
            failcheck=$(grep -i fail ${WIKISTOREDIR}/${DATE}.log);
            faillength=${#failcheck};
            if [ $faillength -eq 0 ] ; then
                image="![ok](ok.png)";
            else
                image="![fail](fail.png)";
            fi
            URL="${COMMIT_URL}${commit}";
            
            # Identify segment
            imonth=$(((($MONTHNOZERO - 1) % 4) + 1));
            segment[$imonth]=" [pdf](${WIKISTOREDIR}/${DATE}.pdf) [log](${WIKISTOREDIR}/${DATE}.log) [${commit:0:7}](${URL}) ${image} ";
            
            newline="|";
            for ((i=0;i<5;i++)); do
               newline="${newline}${segment[i]}|";
            done
            
            # Write to new file
            head -n $line2 ${WIKIFILE} > tempfile;
            echo "$newline" >> tempfile;
            tail -n $line3 ${WIKIFILE} >> tempfile;
            
            # Update original file
            cp tempfile ${WIKIFILE};
            rm tempfile;
            
            # Upload results to bitbucket
            git add ${WIKIFILE} ${WIKISTOREDIR}/${DATE}.log ${WIKISTOREDIR}/${DATE}.pdf >> $LOGFILE;
            git commit -m "Test results for date: ${DATE}" >> $LOGFILE;
            git push origin master >> $LOGFILE;
            
        fi

        # Return to base directory
        cd $SRC;
        
        # Update next backup time
        TNEXTBACKUP=$(date -d "1970-01-01 UTC ${TNEXTBACKUP} second ${frequencies[0]} year ${frequencies[1]} month ${frequencies[2]} day ${frequencies[3]} hour ${frequencies[4]} minute ${frequencies[5]} second" +%s);
        TNEXT=$(date -d "1970-01-01 UTC ${TNEXTBACKUP} second");
        echo "Will perform next backup $TNEXT";

    fi

    # Pause so that script is not checking constantly
    sleep $pause;

done

exit;
