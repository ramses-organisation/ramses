#!/bin/bash

############################################################################
# Ramses Build and TESTScript
#
# Neil Vaytet (ENS Lyon) - 07/2014 - neil.vaytet@ens-lyon.fr
#
# This script runs the RAMSES test suite only once.
# It is currently configured to perform the test suite once and
# upload the results of the tests to the bitbucket wiki, using git commits.
# If you just want to have the results in a local folder, without
# interacting with a wiki page, simply set the variable UPDATEWIKI to false.
############################################################################

# The source directory:
SRC="/home/ubuntu";

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

############################################################################

echo "Starting RAMSES build and test script";

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

