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
SRC="/home/rt3504";

# Test frequency: (YY:MM:DD:hh:mm:ss)
TEST_FREQ="00:00:01:00:00:00";

# Test time offset: (YY:MM:DD:hh:mm:ss)
TEST_OFFS="00:00:00:17:04:00";

# Wiki file
WIKIFILE="AutoTests.md";

# Log file
LOGFILE="daily_test.log";

# Upload to wiki?
UPDATEWIKI=true;

# Set up variables
hline="============================================================";

RAMSESDIR="${SRC}/ramses/tests";
WIKIDIR="${SRC}/wiki";
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
        ./run_test_suite.sh -p 4 >> $LOGFILE;

        # Go to wiki directory
        cd "${WIKIDIR}/${WIKISTOREDIR}";
        cp ${RAMSESDIR}/test_results.pdf ${DATE}.pdf;
        cp ${RAMSESDIR}/test_suite.log ${DATE}.log;

        if $UPDATEWIKI ; then

            # Update wiki page ================================================
            cd ${WIKIDIR};

            # Pull latest wiki ================================================
	        git pull;

            # Call Python script to update the wikifile
            python $RAMSESDIR/update_testpage.py ${WIKIDIR} ${WIKIFILE}

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
