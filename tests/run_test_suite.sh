#!/bin/bash

#######################################################################
#
# Script to run the RAMSES test suite
#
# Neil Vaytet (ENS Lyon) - 07/2014 - neil.vaytet@ens-lyon.fr
# Neil Vaytet (NBI Copenhagen) - 11/2017 - neil.vaytet@nbi.ku.dk
#
# Usage:
#   ./run_test_suite.sh
#
# Options:
#   - Run the suite in parallel (on 4 cpus):
#       ./run_test_suite.sh -p 4
#   - Do not delete results data:
#       ./run_test_suite.sh -d
#   - Run in verbose mode:
#       ./run_test_suite.sh -v
#   - Select test number (for tests 3 to 5, and 10):
#       ./run_test_suite.sh -t 3-5,10
#   - Run all tests in mhd directory:
#       ./run_test_suite.sh -t mhd
#   - Run quick test suite:
#       ./run_test_suite.sh -q
#
#######################################################################

# List of directories to scan
testlist="hydro,mhd,poisson,rt,sink,turb,tracer";

#######################################################################
# Determine the parameters for running the test suite
#######################################################################
MPI=0;
NCPU=1;
VERBOSE=false;
DELDATA=true;
CLEAN_ALL=false;
SELECTTEST=false;
while getopts "cdp:qt:v" OPTION; do
   case $OPTION in
      c)
         CLEAN_ALL=true;
      ;;
      d)
         DELDATA=false;
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
   esac
done

#######################################################################
# Setup paths and commands
#######################################################################
TEST_DIRECTORY=$(pwd);                    # The test suite directory
BASE_DIRECTORY="${TEST_DIRECTORY}/.."; # The main RAMSES directory
BIN_DIRECTORY="${BASE_DIRECTORY}/bin";    # The bin directory
VISU_DIR="${TEST_DIRECTORY}/visu";        # The visualization directory

export PYTHONPATH=${VISU_DIR}:$PYTHONPATH;
DELETE_RESULTS="rm -rf output_* *.tex data*.dat *.pdf *.pyc";
RETURN_TO_BIN="cd ${BIN_DIRECTORY}";
EXECNAME="test_exe_";
LOGFILE="${TEST_DIRECTORY}/test_suite.log";
GIT_URL=$(git config --get remote.origin.url | sed 's/git@bitbucket.org:/https:\/\/bitbucket.org\//g');
GIT_URL=${GIT_URL:0:$((${#GIT_URL}-4))};
THIS_COMMIT=$(git rev-parse HEAD);
echo > $LOGFILE;
if [ ${MPI} -eq 1 ]; then
   RUN_TEST_BASE="mpirun -np ${NCPU} ${BIN_DIRECTORY}/${EXECNAME}";
else
   RUN_TEST_BASE="${BIN_DIRECTORY}/${EXECNAME}";
fi
line="--------------------------------------------";
blankline="                         ";
BEFORETEST="before-test.sh";
AFTERTEST="after-test.sh";

STARTTIME=$(python3 -c 'import time; print(int(time.time()*1000))');

#######################################################################
# Welcome message
#######################################################################
echo "############################################" | tee -a $LOGFILE;
echo "#   Running RAMSES automatic test suite    #" | tee -a $LOGFILE;
echo "############################################" | tee -a $LOGFILE;
if $VERBOSE ; then
   echo "Repository url: ${GIT_URL}" | tee -a $LOGFILE;
   echo "Commit hash: ${THIS_COMMIT}" | tee -a $LOGFILE;
   echo $line | tee -a $LOGFILE;
else
   echo "Repository url: ${GIT_URL}" >> $LOGFILE;
   echo "Commit hash: ${THIS_COMMIT}" >> $LOGFILE;
   echo $line >> $LOGFILE;
fi

#######################################################################
# Generate list of tests from scanning directories
#######################################################################

# Split test list with commas
s1=$(echo $testlist | sed 's/,/ /g');
testsegs_all=( $s1 );
nseg_all=${#testsegs_all[@]};
testlist="";
for ((m=0;m<$nseg_all;m++)); do
   testlist="${testlist} ${testsegs_all[m]}/*";
done

# Count number of tests
testname=( $testlist );
ntestsall=${#testname[@]};
ntests=$ntestsall;
all_tests_ok=true;

#######################################################################
# Clean all directories and exit
#######################################################################
if $CLEAN_ALL ; then
   for ((i=0;i<$ntests;i++)); do
      cd ${TEST_DIRECTORY}/${testname[i]};
      $DELETE_RESULTS;
      if [ -f after_test.sh ]; then
         # rm_list=$(cat to_be_removed);
         # rm -f $rm_list;
         # rm to_be_removed;
         $SHELL after_test.sh;
      fi
   done
   $RETURN_TO_BIN;
   make clean;
   rm -f ${EXECNAME}*d;
   exit;
fi

#######################################################################
# Select particular test if this was asked by user
#######################################################################
if $SELECTTEST ; then

   # Split test selection with commas
   s1=$(echo $TESTNUMBER | sed 's/,/ /g');
   testsegs=( $s1 );
   nseg=${#testsegs[@]};

   # Check if entire directory is submitted
   dir_list="";
   for ((n=0;n<$nseg;n++)); do
      for ((m=0;m<$nseg_all;m++)); do
         if [ ${testsegs[n]} == ${testsegs_all[m]} ] ; then
            dir_list="${dir_list} ${testsegs[n]}/*";
         fi
      done
   done

   # Split list of directories into array
   s1=$(echo $dir_list);
   submit_dirs=( $s1 );
   nsubs=${#submit_dirs[@]};
   ntests=0;
   if [ ${nsubs} -gt 0 ] ; then
      for ((n=0;n<$nsubs;n++)); do
         for ((m=0;m<$ntestsall;m++)); do
            # If directory requested is found in global test list,
            # add it to the current test list
            if [ ${submit_dirs[n]} == ${testname[m]} ] ; then
               testnum[${ntests}]=$m;
               ntests=$((ntests + 1));
            fi
         done
      done

   else

      # Search for dashes in individual segments
      for ((n=0;n<$nseg;n++)); do
         dashsearch=$(echo ${testsegs[n]} | grep '-');
         if [ ${#dashsearch} -gt 0 ] ; then
            istart=$(echo ${testsegs[n]} | cut -d '-' -f1);
            iend=$(echo ${testsegs[n]} | cut -d '-' -f2);
            is=$((istart - 1));
            ie=$((iend - 1));
            iep1=$(($ie + 1));
            for ((j=$is;j<$iep1;j++)); do
               if [ ${j} -ge 0 ] && [ ${j} -lt $ntestsall ] ; then
                  testnum[${ntests}]=$j;
                  ntests=$((ntests + 1));
               else
                  echo "Selected test ${j} does not exist! Ignoring test" | tee -a $LOGFILE;
               fi
            done
         else
            # No dash, just include test in list
            if [ ${testsegs[n]} -gt 0 ] && [ ${testsegs[n]} -le $ntestsall ] ; then
               testnum[${ntests}]=$((${testsegs[n]} - 1));
               ntests=$((ntests + 1));
            else
               echo "Selected test ${testsegs[n]} does not exist! Ignoring test" | tee -a $LOGFILE;
            fi

         fi
      done
   fi

else

   # Include all tests by default
   for ((n=0;n<$ntests;n++)); do
      testnum[n]=$n;
   done

fi

#######################################################################
# Write list of tests
#######################################################################
if [ $ntests -eq 0 ] ; then
   echo "The test list is empty." | tee -a $LOGFILE;
   exit;
fi
echo "Will perform the following tests:" | tee -a $LOGFILE;
for ((i=0;i<$ntests;i++)); do
   n=${testnum[i]};
   j=$(($n + 1));
   if [ $j -lt 10 ] ; then
      echo " [ ${j}] ${testname[n]}" | tee -a $LOGFILE;
   else
      echo " [${j}] ${testname[n]}" | tee -a $LOGFILE;
   fi
done
echo $line | tee -a $LOGFILE;

#######################################################################
# Loop through all tests
#######################################################################
for ((i=0;i<$ntests;i++)); do

   # Start timer for test, including compilations
   STARTTIME_GLOB=$(python3 -c 'import time; print(int(time.time()*1000))');

   # Get test number
   n=${testnum[i]};
   ip1=$(($i + 1));
   echo "Test ${ip1}/${ntests}: ${testname[n]}" | tee -a $LOGFILE;

   # Get raw test name for namelist, pdf and tex files
   nslash=$(grep -o "/" <<< "${testname[n]}" | wc -l);
   if [ $nslash -gt 0 ] ; then
      np1=$(($nslash + 1));
      rawname[i]=$(echo ${testname[n]} | cut -d '/' -f$np1);
   else
      rawname[i]=${testname[n]};
   fi

   # Read test configuration file and extract NDIM (needed for executable)
   FLAGS=$(grep FLAGS ${TEST_DIRECTORY}/${testname[n]}/config.txt | cut -d ':' -f2);
   flag_split=( $FLAGS );
   nflags=${#flag_split[@]};
   for ((k=0;k<$nflags;k++)); do
      if [ ${flag_split[$k]:0:4} = "NDIM" ] ; then
         ndim=$(echo ${flag_split[$k]} | cut -d '=' -f2);
      fi
   done

   # Initial cleanup
   $RETURN_TO_BIN;
   if ${make_clean[n]}; then
      echo "Cleanup" | tee -a $LOGFILE;
      if $VERBOSE ; then
         make clean 2>&1 | tee -a $LOGFILE;
      else
         make clean >> $LOGFILE 2>&1;
      fi
   fi

   # Compile source
   echo "Compiling source" | tee -a $LOGFILE;
   MAKESTRING="make EXEC=${EXECNAME} MPI=${MPI} ${FLAGS}";
   # if [ ${MPI} -eq 1 ]; then
   #    MAKESTRING="${MAKESTRING} -j ${NCPU}";
   # fi
   if $VERBOSE ; then
      $MAKESTRING 2>&1 | tee -a $LOGFILE;
   else
      $MAKESTRING >> $LOGFILE 2>&1;
   fi

   # Run test
   cd ${TEST_DIRECTORY}/${testname[n]};
   $DELETE_RESULTS;
   RUN_TEST="${RUN_TEST_BASE}${ndim}d ${rawname[i]}.nml";
   echo -n "Running test:" | tee -a $LOGFILE;
   STARTTIME_TEST=$(python3 -c 'import time; print(int(time.time()*1000))');
   # prepname="prepare-${rawname[i]}.sh";
   if $VERBOSE ; then
      if [ -f ${BEFORETEST} ]; then
         ${SHELL} ${BEFORETEST} 2>&1 | tee -a $LOGFILE;
      fi
      ${RUN_TEST} 2>&1 | tee -a $LOGFILE;
   else
      if [ -f ${BEFORETEST} ]; then
         ${SHELL} ${BEFORETEST} >> $LOGFILE 2>&1;
      fi
      ${RUN_TEST} >> $LOGFILE 2>&1;
   fi
   # Record test time
   ENDTIME_TEST=$(python3 -c 'import time; print(int(time.time()*1000))');
   milliseconds=$(($ENDTIME_TEST - $STARTTIME_TEST));
   seconds=$(($milliseconds / 1000));
   hours=$(($seconds / 3600));
   seconds=$(($seconds % 3600));
   minutes=$(($seconds / 60));
   seconds=$(($seconds % 60));
   hours_test[${i}]=$hours;
   minutes_test[${i}]=$minutes;
   seconds_test[${i}]=$seconds;
   echo " ${hours_test[i]}h${minutes_test[i]}m${seconds_test[i]}s" | tee -a $LOGFILE;

   # Plot and analyse results
   echo "Plotting and analysing results" | tee -a $LOGFILE;
   status=$(python3 plot-${rawname[i]}.py 2>&1);
   if $VERBOSE ; then
      echo $status;
   fi
   echo $status >> $LOGFILE;

   # Print message on test status
   ispassed=$(echo $status | grep PASSED);
   length=${#testname[n]};
   if [ ${#ispassed} -gt 0 ]; then
      test_failed[n]=false;
      echo "Test ${testname[n]} passed ${blankline:$length}[ OK ]" | tee -a $LOGFILE;
   else
      test_failed[n]=true;
      echo "Test ${testname[n]} failed!${blankline:$length}[FAIL]" | tee -a $LOGFILE;
      all_tests_ok=false;
   fi

   echo $line | tee -a $LOGFILE;

   # Record global time including compilations
   ENDTIME_GLOB=$(python3 -c 'import time; print(int(time.time()*1000))');
   milliseconds=$(($ENDTIME_GLOB - $STARTTIME_GLOB));
   seconds=$(($milliseconds / 1000));
   hours=$(($seconds / 3600));
   seconds=$(($seconds % 3600));
   minutes=$(($seconds / 60));
   seconds=$(($seconds % 60));
   hours_glob[${i}]=$hours;
   minutes_glob[${i}]=$minutes;
   seconds_glob[${i}]=$seconds;

done

# Total time ##########################################################
ENDTIME=$(python3 -c 'import time; print(int(time.time()*1000))');
milliseconds=$(($ENDTIME - $STARTTIME));
seconds=$(($milliseconds / 1000));
hours=$(($seconds / 3600));
seconds=$(($seconds % 3600));
minutes=$(($seconds / 60));
seconds=$(($seconds % 60));
#######################################################################

#######################################################################
# Generate pdf document with test results
#######################################################################
echo "Generating pdf document with test results" | tee -a $LOGFILE;
cd ${TEST_DIRECTORY};
latexfile="test_results.tex";
echo "\documentclass[12pt]{article}" > $latexfile;
echo "\usepackage{graphicx,color,caption}" >> $latexfile;
echo "\usepackage[colorlinks=true,linkcolor=blue]{hyperref}" >> $latexfile;
echo "\topmargin -1.3in" >> $latexfile;
echo "\textheight 10.1in" >> $latexfile;
echo "\oddsidemargin -0.7in" >> $latexfile;
echo "\evensidemargin -0.7in" >> $latexfile;
echo "\textwidth 7.7in" >> $latexfile;
echo "\title{RAMSES test suite results}" >> $latexfile;
echo "\date{\today}" >> $latexfile;
echo "\author{${USER}}" >> $latexfile;
echo "\nonstopmode" >> $latexfile;
echo "\begin{document}" >> $latexfile;
echo "\maketitle" >> $latexfile;
echo "\begin{center}" >> $latexfile;
SAFE_URL=$(echo ${GIT_URL})
echo "Commit hash: \href{${GIT_URL}/commits/${THIS_COMMIT}}{${THIS_COMMIT:0:6}}" >> $latexfile;
echo "\end{center}" >> $latexfile;
echo "\begin{table}[ht]" >> $latexfile;
echo "\centering" >> $latexfile;
echo "\caption*{Test run summary using ${NCPU} processor(s)}" >> $latexfile;
echo "\begin{tabular}{|r|l|l|l|l|}" >> $latexfile;
echo "\hline" >> $latexfile;
echo "~ & Test name & Run time & Total time & Status\\\\" >> $latexfile;
echo "\hline" >> $latexfile;
for ((i=0;i<$ntests;i++)); do
   n=${testnum[i]};
   itest=$(($n + 1));
   if ${test_failed[n]} ; then
      status="\hyperref[fig-${testname[n]}]{\textcolor{red}{failed}}";
   else
      status="\hyperref[fig-${testname[n]}]{\textcolor{green}{passed}}";
   fi
   echo "$itest & \hyperref[fig-${testname[n]}]{${testname[n]}} & ${hours_test[i]}h${minutes_test[i]}m${seconds_test[i]}s & ${hours_glob[i]}h${minutes_glob[i]}m${seconds_glob[i]}s & ${status} \\\\" >> $latexfile;
done
echo "\hline" >> $latexfile;
echo "\end{tabular}" >> $latexfile;
echo "\end{table}" >> $latexfile;
echo "\begin{center}" >> $latexfile;
echo "Total run time (including compilations): ${hours}h${minutes}m${seconds}s" >> $latexfile;
echo "\end{center}" >> $latexfile;
echo "\clearpage" >> $latexfile;

for ((i=0;i<$ntests;i++)); do
   n=${testnum[i]};
   echo "\begin{figure}" >> $latexfile;
   echo "\centering" >> $latexfile;
   pdfname=${TEST_DIRECTORY}/${testname[n]}/${rawname[i]}.pdf;
   if [ ! -f $pdfname ]; then
      echo "\begin{tabular}{|c|}" >> $latexfile;
      echo "\hline" >> $latexfile;
      echo "~\\\\" >> $latexfile;
      echo "{\LARGE MISSING:}\\\\" >> $latexfile;
      echo "{\LARGE PDF FILE}\\\\" >> $latexfile;
      echo "~\\\\" >> $latexfile;
      echo "\hline" >> $latexfile;
      echo "\end{tabular}" >> $latexfile;
   else
      echo "\includegraphics[height=0.5\textheight,width=\textwidth,keepaspectratio]{$pdfname}" >> $latexfile;
   fi
   echo "\caption{${testname[n]} test}" >> $latexfile;
   echo "\label{fig-${testname[n]}}" >> $latexfile;
   echo "\end{figure}" >> $latexfile;
   texname=${TEST_DIRECTORY}/${testname[n]}/${rawname[i]}.tex;
   if [ ! -f $texname ]; then
      echo "\begin{table}[ht]" >> $latexfile;
      echo "\centering" >> $latexfile;
      echo "\begin{tabular}{|c|}" >> $latexfile;
      echo "\hline" >> $latexfile;
      echo "~\\\\" >> $latexfile;
      echo "{\LARGE MISSING:}\\\\" >> $latexfile;
      echo "{\LARGE STATS FILE}\\\\" >> $latexfile;
      echo "~\\\\" >> $latexfile;
      echo "\hline" >> $latexfile;
      echo "\end{tabular}" >> $latexfile;
      echo "\end{table}" >> $latexfile;
   else
      echo "\input{$texname}" >> $latexfile;
   fi
   echo "\clearpage" >> $latexfile;
done
echo "\end{document}" >> $latexfile;
if $VERBOSE ; then
   pdflatex $latexfile 2>&1 | tee -a $LOGFILE;
   pdflatex $latexfile 2>&1 | tee -a $LOGFILE;
else
   pdflatex $latexfile >> $LOGFILE 2>&1;
   pdflatex $latexfile >> $LOGFILE 2>&1;
fi
rm ${latexfile/.tex/.log};
rm ${latexfile/.tex/.aux};
rm ${latexfile/.tex/.out};
rm $latexfile;

#######################################################################
# Clean up
#######################################################################
if $all_tests_ok ; then
   echo "All tests were completed successfully" | tee -a $LOGFILE;
else
   echo "There were some failed tests" | tee -a $LOGFILE;
fi
if ${DELDATA} ; then
   for ((i=0;i<$ntests;i++)); do
      n=${testnum[i]};
      cd ${TEST_DIRECTORY}/${testname[n]};
      $DELETE_RESULTS;
      if [ -f ${AFTERTEST} ]; then
         # rm_list=$(cat to_be_removed);
         # rm -f $rm_list;
         # rm to_be_removed;
         ${SHELL} ${AFTERTEST};
      fi
   done
   $RETURN_TO_BIN;
   if $VERBOSE ; then
      make clean 2>&1 | tee -a $LOGFILE;
   else
      make clean >> $LOGFILE 2>&1;
   fi
   rm -f ${EXECNAME}*d;
fi

exit;
