#!/bin/bash

#######################################################################
#
# Script to run the RAMSES test suite
# 
# Neil Vaytet (ENS Lyon) - 07/2014 - neil.vaytet@ens-lyon.fr
#
# Usage:
#   ./run_test_suite.sh
#
# Options:
#   - Run the suite in parallel (on 4 cpus):
#       ./run_test_suite.sh -p 4
#   - Run using python visualization instead of gnuplot:
#       ./run_test_suite -y
#   - Run full test suite (including additional tests):
#       ./run_test_suite -f
#   - Do not delete results data:
#       ./run_test_suite -d
#   - Run in verbose mode:
#       ./run_test_suite -v
#   - Select test number (for tests 3 to 5, and 10):
#       ./run_test_suite -t 3-5,10
#
#######################################################################

#######################################################################
# Determine the parameters for running the test suite
#######################################################################
MPI=0;
NCPU=1;
USE_GNUPLOT=true;
USE_PYTHON=false;
RUN_FULL_SUITE=false;
VERBOSE=false;
DELDATA=true;
SELECTTEST=false;
while getopts "dfp:t:vy" OPTION; do
   case $OPTION in
      d)
         DELDATA=false;
      ;;
      f)
         RUN_FULL_SUITE=true;
      ;;
      p)
         MPI=1;
         NCPU=$OPTARG;
      ;;
      t)
         SELECTTEST=true;
         TESTNUMBER=$OPTARG;
      ;;
      v)
         VERBOSE=true;
      ;;
      y)
         USE_PYTHON=true;
         USE_GNUPLOT=false;
      ;;
   esac
done

#######################################################################
# Setup paths and commands
#######################################################################
TEST_DIRECTORY=$(pwd);

length=${#TEST_DIRECTORY};
icut=$(($length - 17));

BASE_DIRECTORY="${TEST_DIRECTORY:0:${icut}}";

BIN_DIRECTORY="${BASE_DIRECTORY}/bin";
VISU_DIR="${TEST_DIRECTORY}/visu";
DELETE_RESULTS="rm -rf output_* log data*.dat time.dat";
DELETE_SOURCES="rm -f units.o condinit.o";
RETURN_TO_BIN="cd ${BIN_DIRECTORY}";
EXECNAME="test_exe_";
LOGFILE="${TEST_DIRECTORY}/test_suite.log";
echo > $LOGFILE;

if [ ${MPI} -eq 1 ]; then
   RUN_TEST_BASE="mpirun -np ${NCPU} ${BIN_DIRECTORY}/${EXECNAME}";
else
   RUN_TEST_BASE="${BIN_DIRECTORY}/${EXECNAME}";
fi

STARTTIME=$(date +%s);

#######################################################################
# Welcome message
#######################################################################
echo "############################################";
echo "############################################" >> $LOGFILE;
if $RUN_FULL_SUITE ; then
   echo "Running extended RAMSES automatic test suite";
   echo "Running extended RAMSES automatic test suite" >> $LOGFILE;
else
   echo "Running standard RAMSES automatic test suite";
   echo "Running standard RAMSES automatic test suite" >> $LOGFILE;
fi
echo "############################################";
echo "############################################" >> $LOGFILE;

#######################################################################
# List of tests
#######################################################################

itest=0; # Test 1
testdir[${itest}]="sod-tube";
testname[${itest}]="sod-tube";
testpatch[${itest}]="";
testlist[${itest}]="sod-tube.nml";
ndim[${itest}]=1;
nvar[${itest}]=5;
solver[${itest}]="hydro";
flags[${itest}]="";
make_clean[${itest}]=true;
del_files[${itest}]="";

itest=$((itest + 1)); # Test 2
testdir[${itest}]="imhd-tube";
testname[${itest}]="imhd-tube";
testpatch[${itest}]="../mhd";
testlist[${itest}]="imhd-tube.nml";
ndim[${itest}]=1;
nvar[${itest}]=8;
solver[${itest}]="mhd";
flags[${itest}]="";
make_clean[${itest}]=true;
del_files[${itest}]="";

itest=$((itest + 1)); # Test 3
testdir[${itest}]="orszag-tang";
testname[${itest}]="orszag-tang";
testpatch[${itest}]="../patch/test_suite/orszag-tang";
testlist[${itest}]="orszag-tang.nml";
ndim[${itest}]=2;
nvar[${itest}]=8;
solver[${itest}]="mhd";
flags[${itest}]="";
make_clean[${itest}]=true;
del_files[${itest}]="output_*";

# Store number of standard tests
ntestsstandard=${#testname[@]};

# Additional tests: include your own tests here ==============
# ============================================================

# Store total number of tests
ntestsfull=${#testname[@]};

# Count number of tests
if ${RUN_FULL_SUITE} ; then
   ntestsall=${ntestsfull};
else
   ntestsall=${ntestsstandard};
fi

ntests=$ntestsall;

# Count number of tests
ntests=${#testname[@]};
all_tests_ok=true;

#######################################################################
# Select particular test if this was asked by user
#######################################################################
if $SELECTTEST ; then

   # Split test selection with commas
   s1=$(echo $TESTNUMBER | sed 's/,/ /'g);
   testsegs=( $s1 );
   nseg=${#testsegs[@]};
   
   # Search for dashes in individual segments
   ntests=0;
   for ((n=0;n<$nseg;n++)); do
      dashsearch=$(echo ${testsegs[n]} | grep '-');
      if [ ${#dashsearch} -gt 0 ] ; then
         istart=$(echo ${testsegs[n]} | cut -d '-' -f1);
         iend=$(echo ${testsegs[n]} | cut -d '-' -f2);
         is=$((istart - 1));
         ie=$((iend - 1));
         iep1=$(($ie + 1));
         for ((j=$is;j<$iep1;j++)); do
            if [ ${j} -ge 0 ] && [ ${j} -lt $ntestsfull ] ; then
               testnum[${ntests}]=$j;
               ntests=$((ntests + 1));
            else
               echo "Selected test ${j} does not exist! Ignoring test";
               echo "Selected test ${j} does not exist! Ignoring test" >> $LOGFILE;
            fi
         done
      else
         # No dash, just include test in list
         testnum[${ntests}]=$((${testsegs[n]} - 1));
         if [ ${testnum[${ntests}]} -gt $ntestsfull ] ; then
            echo "Selected test does not exist!";
            echo "Selected test does not exist!" >> $LOGFILE;
            exit;
         fi
         ntests=$((ntests + 1));
      fi
   done

else

   # Include all tests by default
   for ((n=0;n<$ntests;n++)); do
      testnum[n]=$n;
   done
   
fi

#######################################################################
# Write list of tests
#######################################################################
echo "Will perform the following tests:";
echo "Will perform the following tests:" >> $LOGFILE;
for ((i=0;i<$ntests;i++)); do
   n=${testnum[i]};
   j=$(($n + 1));
   if [ $ntests -gt 9 ] && [ $j -lt 10 ] ; then
      echo " [ ${j}] ${testname[n]}";
      echo " [ ${j}] ${testname[n]}" >> $LOGFILE;
   else
      echo " [${j}] ${testname[n]}";
      echo " [${j}] ${testname[n]}" >> $LOGFILE;
   fi
done
echo "--------------------------------------------";
echo "--------------------------------------------" >> $LOGFILE;

#######################################################################
# Prepare visualization software
#######################################################################
echo "Compiling visualization software";
echo "Compiling visualization software" >> $LOGFILE;
cd ${VISU_DIR};
if $VERBOSE ; then
   make clean;
   make;
else
   { make clean >> $LOGFILE; } 2>> $LOGFILE;
   { make >> $LOGFILE; } 2>> $LOGFILE;
fi
echo "--------------------------------------------";
echo "--------------------------------------------" >> $LOGFILE;

#######################################################################
# Loop through all tests
#######################################################################
itest=0;
for ((i=0;i<$ntests;i++)); do

   n=${testnum[i]};
   itest=$(($itest + 1));
   echo "Test ${itest}/${ntests}: ${testname[n]}";
   echo "Test ${itest}/${ntests}: ${testname[n]}" >> $LOGFILE;
      
   # Initial cleanup
   $RETURN_TO_BIN;
   if ${make_clean[n]}; then
      echo "Cleanup";
      echo "Cleanup" >> $LOGFILE;
      if $VERBOSE ; then
         make clean;
      else
         { make clean >> $LOGFILE; } 2>> $LOGFILE;
      fi
   fi
   rm -f ${del_files[n]};
   
   # Compile source
   echo "Compiling source";
   echo "Compiling source" >> $LOGFILE;
   if $VERBOSE ; then
      make EXEC=${EXECNAME} PATCH=${testpatch[n]} SOLVER=${solver[n]} MPI=${MPI} NDIM=${ndim[n]} NVAR=${nvar[n]} ${flags[n]};
   else
      { make EXEC=${EXECNAME} PATCH=${testpatch[n]} SOLVER=${solver[n]} MPI=${MPI} NDIM=${ndim[n]} NVAR=${nvar[n]} ${flags[n]} >> $LOGFILE; } 2>> $LOGFILE;
   fi
   
   # Run tests
   cd ${TEST_DIRECTORY}/${testdir[n]};
   $DELETE_RESULTS;
   if $VERBOSE ; then
      ./prepare-${testname[n]}.sh;
   else
      { ./prepare-${testname[n]}.sh >> $LOGFILE; } 2>> $LOGFILE;
   fi
   RUN_TEST="${RUN_TEST_BASE}${ndim[n]}d ${testlist[n]}";
   echo "Running test";
   echo "Running test" >> $LOGFILE;
   { time ${RUN_TEST} > log ; } 2> ${testname[n]}"_stats.txt";
   cat log >> $LOGFILE;

   # Plot results
   echo "Plotting results";
   echo "Plotting results" >> $LOGFILE;
   if $VERBOSE ; then
      ./plot-${testname[n]}.sh;
   else
      { ./plot-${testname[n]}.sh >> $LOGFILE; } 2>> $LOGFILE;
   fi
   if ${USE_GNUPLOT} ; then
      gnuplot plot-${testname[n]}.gp;
      ps2pdf ${testname[n]}.ps;
      rm ${testname[n]}.ps;
   fi
   if ${USE_PYTHON} ; then
      python plot-${testname[n]}.py;
   fi
   
   # Check for differences in results
   echo "Analysing results";
   echo "Analysing results" >> $LOGFILE;
   difffile="resdiff-${testname[n]}";
   diff data.dat ${testname[n]}"-ref.dat" > ${difffile};
   # Size of diff file?
   diffoutput=$(cat ${difffile});
   diffsize=${#diffoutput};
   if [ ${diffsize} -gt 0 ]; then
      diff_not_empty[n]=true;
      echo "Test failed!                          [FAIL]";
      echo "Test failed!                          [FAIL]" >> $LOGFILE;
      all_tests_ok=false;
   else
      diff_not_empty[n]=false;
      echo "Test passed                           [ OK ]";
      echo "Test passed                           [ OK ]" >> $LOGFILE;
   fi
   
#    # Check for differences in log files
#    echo "Analysing logs"; echo "Analysing logs" >> $LOGFILE;
#    # Remove grids and memory from log file as they change with MPI
#    sed -i 's/has.*$//' log;
#    sed -i 's/mem.*$//' log;
#    # Apply difference (ignoring elapsed times)
#    diff log ${testname[n]}.log -I elapsed > ${testname[n]}"_logdiff.tex";
#    # Size of diff file?
#    diffoutput=$(cat ${testname[n]}"_logdiff.tex");
#    diffsize=${#diffoutput};
#    if [ ${diffsize} -gt 0 ]; then
#       diff_not_empty[n]=true;
#       echo "Test failed!"; echo "Test failed!" >> $LOGFILE;
#       all_tests_ok=false;
#       # Format diff output
#       sed -i 's/</\$<\$/g' ${testname[n]}"_logdiff.tex";
#       sed -i 's/>/\$>\$/g' ${testname[n]}"_logdiff.tex";
#       sed -i 's/%/\\%/g' ${testname[n]}"_logdiff.tex";
#       sed -i 's/$/ \\\\/g' ${testname[n]}"_logdiff.tex";
#    else
#       diff_not_empty[n]=false;
#       echo "Test passed"; echo "Test passed" >> $LOGFILE;
#    fi
   
   echo "--------------------------------------------";
   echo "--------------------------------------------" >> $LOGFILE;

done

#######################################################################
ENDTIME=$(date +%s);
seconds=$(($ENDTIME - $STARTTIME));
hours=$((seconds / 3600));
seconds=$((seconds % 3600));
minutes=$((seconds / 60));
seconds=$((seconds % 60));
#######################################################################

#######################################################################
# Generate pdf document with test results
#######################################################################
echo "Generating pdf document with test results";
echo "Generating pdf document with test results" >> $LOGFILE;
cd ${TEST_DIRECTORY};
latexfile="test_results.tex";
echo "\documentclass[12pt]{article}" > $latexfile;
echo "\usepackage{graphicx,color}" >> $latexfile;
echo "\usepackage[colorlinks=true,linkcolor=blue]{hyperref}" >> $latexfile;
echo "\topmargin -1.3in" >> $latexfile;
echo "\textheight 10.1in" >> $latexfile;
echo "\oddsidemargin -0.7in" >> $latexfile;
echo "\evensidemargin -0.7in" >> $latexfile;
echo "\textwidth 7.7in" >> $latexfile;
echo >> $latexfile;
echo "\title{RAMSES test suite results}" >> $latexfile;
echo "\date{\today}" >> $latexfile;
echo "\author{${USER}}" >> $latexfile;
echo >> $latexfile;
echo "\begin{document}" >> $latexfile;
echo >> $latexfile;
echo "\maketitle" >> $latexfile;
echo >> $latexfile;
echo "\begin{table}[ht]" >> $latexfile;
echo "\centering" >> $latexfile;
echo "\caption{Test run summary using ${NCPU} processor(s)}" >> $latexfile;
echo "\begin{tabular}{|r|l|l|l|l|}" >> $latexfile;
echo "\hline" >> $latexfile;
echo "~ & Test name & Real time & User time & Status\\\\" >> $latexfile;
echo "\hline" >> $latexfile;
for ((i=0;i<$ntests;i++)); do
   n=${testnum[i]};
   statfile="${TEST_DIRECTORY}/${testdir[n]}/${testname[n]}_stats.txt"
   itest=$(($n + 1));
   if ${diff_not_empty[n]} ; then
      status="\hyperref[fig-${testname[n]}]{\textcolor{red}{failed}}";
   else
      status="\hyperref[fig-${testname[n]}]{\textcolor{green}{passed}}";
   fi
   echo $itest "& \hyperref[fig-${testname[n]}]{${testname[n]}} &" $(grep real ${statfile} | cut -d 'l' -f2) "&" $(grep user ${statfile} | cut -d 'r' -f2) "&" ${status} "\\\\" >> $latexfile;
done
echo "\hline" >> $latexfile;
echo "\end{tabular}" >> $latexfile;
echo "\end{table}" >> $latexfile;
echo "\begin{center}" >> $latexfile;
echo "Total run time (including compilations): ${hours}h${minutes}m${seconds}s" >> $latexfile;
echo "\end{center}" >> $latexfile;
echo "\clearpage" >> $latexfile;
echo >> $latexfile;

for ((i=0;i<$ntests;i++)); do
   n=${testnum[i]};
   echo "\begin{figure}" >> $latexfile;
   echo "\centering" >> $latexfile;
   echo "\includegraphics[scale=0.7]{${TEST_DIRECTORY}/${testdir[n]}/${testname[n]}.pdf}" >> $latexfile;
   echo "\caption{${testname[n]} test}" >> $latexfile;
   echo "\label{fig-${testname[n]}}" >> $latexfile;
   echo "\end{figure}" >> $latexfile;
   echo "\clearpage" >> $latexfile;
   echo >> $latexfile;
#    if ${diff_not_empty[n]} ; then
#       echo "{\bf Differences in log for test ${testname[n]}:}\\\\" >> $latexfile;
#       echo "\input{${TEST_DIRECTORY}/${testdir[n]}/${testname[n]}_logdiff.tex}" >> $latexfile;
#       echo "\clearpage" >> $latexfile;
#       echo >> $latexfile;
#    fi
done
echo "\end{document}" >> $latexfile;
if $VERBOSE ; then
   pdflatex $latexfile;
   pdflatex $latexfile;
else
   pdflatex $latexfile >> $LOGFILE;
   pdflatex $latexfile >> $LOGFILE;
fi
rm ${latexfile/.tex/.log};
rm ${latexfile/.tex/.aux};
rm ${latexfile/.tex/.out};
rm $latexfile;

#######################################################################
# Clean up
#######################################################################
if $all_tests_ok ; then
   echo "All tests were completed successfully";
   echo "All tests were completed successfully" >> $LOGFILE;
else
   echo "There were some failed tests";
   echo "There were some failed tests" >> $LOGFILE;
fi
if ${DELDATA} ; then
   for ((i=0;i<$ntests;i++)); do
      n=${testnum[i]};
      cd ${TEST_DIRECTORY}/${testdir[n]};
      $DELETE_RESULTS;
      rm ${testname[n]}"_stats.txt" ${testname[n]}".pdf"
   done
   if $VERBOSE ; then
      cd ${VISU_DIR}
      make clean;
      $RETURN_TO_BIN;
      make clean;
   else
      cd ${VISU_DIR}
      make clean >> $LOGFILE;
      $RETURN_TO_BIN;
      make clean >> $LOGFILE;
   fi
   rm -f ${EXECNAME}*d;
fi

exit;
