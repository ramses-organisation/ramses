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
#   - Run the suite in parallel with MPI (on 4 cpus):
#       ./run_test_suite.sh -p 4
#       An individual test can override this and pin itself to a fixed number
#       of MPI processes with a line "NPROC: <n>" in its config.txt. This is
#       for tests whose reference solution depends on the decomposition.
#   - Run the suite in parallel with OPENMP (on 4 cpus):
#       ./run_test_suite.sh -m 4
#   - Run the suite in parallel with MPI+OPENMP (on 4 cpus):
#       ./run_test_suite.sh -p 2 -m 2
#   - OpenMP behaviour:
#       Giving -m runs every test with OpenMP. Without it, only the tests that
#       opt in via their config.txt (a line "OPENMP: true") are built and run
#       with OpenMP, the others are built without it.
#       ./run_test_suite.sh            # the opted-in tests only, 1 thread
#       ./run_test_suite.sh -o -m 2    # same tests, 2 threads
#       ./run_test_suite.sh -o all     # force OpenMP on ALL tests
#       ./run_test_suite.sh -o none    # disable OpenMP for ALL tests
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
#   - Run test suite with coverage:
#       ./run_test_suite.sh -s
#   - Restart behaviour:
#       By default, only tests that opt in via their config.txt
#       (a line "RESTART: true") are run through the restart mechanism.
#       ./run_test_suite.sh            # restart the opted-in tests only
#       ./run_test_suite.sh -r         # force restart on ALL tests
#       ./run_test_suite.sh -r all     # same as -r
#       ./run_test_suite.sh -r none    # disable restart for ALL tests
#
#######################################################################

# List of directories to scan
testlist="hydro,mhd,poisson,rt,sink,star,turb,tracer";

#######################################################################
# Determine the parameters for running the test suite
#######################################################################
MPI=0;
GCOV=0;
OPENMP=0;
NPROC=1;
NTHREADS=1;
VERBOSE=false;
DELDATA=true;
COVERAGE=false;
CLEAN_ALL=false;
SELECTTEST=false;
# Restart mode: "default" (only tests opted-in via config.txt), "all" or "none".
RESTART_MODE="default";
# OpenMP mode: "default" (only tests opted-in via config.txt), "all" or "none".
# Giving -m implies "all", so that -m N keeps its old meaning of running the
# whole suite with N threads, but an explicit -o always wins, whatever the
# order of the options.
OPENMP_MODE="default";
OPENMP_MODE_SET=false;
while getopts "cdsp:qm:t:vro" OPTION; do
   case $OPTION in
      c)
         CLEAN_ALL=true;
      ;;
      d)
         DELDATA=false;
      ;;
      s)
         GCOV=1;
         COVERAGE=true;
      ;;
      p)
         MPI=1;
         NPROC=$OPTARG;
      ;;
      m)
         OPENMP=1;
         NTHREADS=$OPTARG;
      ;;
      t)
         SELECTTEST=true;
         TESTNUMBER=$OPTARG;
      ;;
      v)
         VERBOSE=true;
      ;;
      r)
         # -r takes an optional argument ("all" or "none"). getopts treats -r as
         # a plain flag, so peek at the next token and consume it only if it is
         # one of the recognised keywords; otherwise -r means "all".
         eval nextarg=\${$OPTIND};
         if [ "$nextarg" = "none" ] || [ "$nextarg" = "all" ] ; then
            RESTART_MODE=$nextarg;
            OPTIND=$((OPTIND + 1));
         else
            RESTART_MODE="all";
         fi
      ;;
      o)
         # -o takes an optional argument ("all" or "none"), same trick as -r.
         eval nextarg=\${$OPTIND};
         if [ "$nextarg" = "none" ] || [ "$nextarg" = "all" ] ; then
            OPENMP_MODE=$nextarg;
            OPTIND=$((OPTIND + 1));
         else
            OPENMP_MODE="all";
         fi
         OPENMP_MODE_SET=true;
      ;;
   esac
done

# -m on its own means "run everything with OpenMP", unless -o said otherwise.
if ! ${OPENMP_MODE_SET} && [ ${OPENMP} -eq 1 ] ; then
   OPENMP_MODE="all";
fi
# Conversely, tests opting in via config.txt need the OpenMP settings exported
# even when -m was not given, in which case they run on a single thread.
if [ "${OPENMP_MODE}" = "none" ] ; then
   OPENMP=0;
else
   OPENMP=1;
fi

#######################################################################
# Setup paths and commands
#######################################################################
TEST_DIRECTORY=$(pwd);                    # The test suite directory
BASE_DIRECTORY="${TEST_DIRECTORY}/..";    # The main RAMSES directory
BIN_DIRECTORY="${BASE_DIRECTORY}/bin";    # The bin directory
VISU_DIR="${TEST_DIRECTORY}/visu";        # The visualization directory

export PYTHONPATH=${VISU_DIR}:$PYTHONPATH;
DELETE_RESULTS="rm -rf output_* *.tex data*.dat *.pdf *.pyc *.gc* coverage_stats.txt movie1";
RETURN_TO_BIN="cd ${BIN_DIRECTORY}";
EXECNAME="test_exe_";
LOGFILE="${TEST_DIRECTORY}/test_suite.log";
GIT_URL=$(git config --get remote.origin.url | sed 's/git@github.com:/https:\/\/github.com\//g');
GIT_URL=${GIT_URL:0:$((${#GIT_URL}-4))};
THIS_COMMIT=$(git rev-parse HEAD);
echo > $LOGFILE;
# The thread placement settings are the same for every test, only the number of
# threads varies, and that is set per test inside the loop below.
if [ ${OPENMP} -eq 1 ]; then
   export OMP_PLACES=cores
   export OMP_PROC_BIND=true
   export OMP_STACKSIZE=2048M
fi
if [ ${MPI} -eq 1 ]; then
   RUN_TEST_BASE="mpirun --map-by slot:pe=${NTHREADS} --np ${NPROC} ${BIN_DIRECTORY}/${EXECNAME}";
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

   # Decide whether this test is run through the restart mechanism.
   #   -r all  -> every test ; -r none -> no test ;
   #   default -> only tests whose config.txt contains "RESTART: true".
   DO_RESTART=false;
   case $RESTART_MODE in
      all)
         DO_RESTART=true;
      ;;
      none)
         DO_RESTART=false;
      ;;
      *)
         testrestart=$(grep -i '^[[:space:]]*RESTART[[:space:]]*:' ${TEST_DIRECTORY}/${testname[n]}/config.txt | cut -d ':' -f2 | tr -d '[:space:]' | tr '[:upper:]' '[:lower:]');
         if [ "$testrestart" = "true" ] ; then
            DO_RESTART=true;
         fi
      ;;
   esac

   # Decide whether this test is compiled and run with OpenMP.
   #   -o all  -> every test ; -o none -> no test ;
   #   default -> only tests whose config.txt contains "OPENMP: true".
   DO_OPENMP=false;
   case $OPENMP_MODE in
      all)
         DO_OPENMP=true;
      ;;
      none)
         DO_OPENMP=false;
      ;;
      *)
         testopenmp=$(grep -i '^[[:space:]]*OPENMP[[:space:]]*:' ${TEST_DIRECTORY}/${testname[n]}/config.txt | cut -d ':' -f2 | tr -d '[:space:]' | tr '[:upper:]' '[:lower:]');
         if [ "$testopenmp" = "true" ] ; then
            DO_OPENMP=true;
         fi
      ;;
   esac

   # Set the number of threads to use for this test. A test that does not use
   # OpenMP still has to be built and run, just with a single thread, so that
   # the reference solutions stay comparable.
   TEST_OPENMP=0;
   TEST_NTHREADS=1;
   if ${DO_OPENMP} ; then
      TEST_OPENMP=1;
      TEST_NTHREADS=${NTHREADS};
      echo "Test uses OpenMP with ${TEST_NTHREADS} thread(s)" | tee -a $LOGFILE;
   fi
   export OMP_NUM_THREADS=${TEST_NTHREADS};

   # Set the number of MPI processes to use for this test
   TEST_NPROC=${NPROC};
   if [ ${MPI} -eq 1 ] ; then
      # Check of NPROC was specified in config
      testnproc=$(grep -i '^[[:space:]]*NPROC[[:space:]]*:' ${TEST_DIRECTORY}/${testname[n]}/config.txt | cut -d ':' -f2 | tr -d '[:space:]');
      if [ -n "$testnproc" ] ; then
         TEST_NPROC=$testnproc;
         echo "Test fixed to ${TEST_NPROC} MPI process(es)" | tee -a $LOGFILE;
      fi
      RUN_TEST_BASE="mpirun --map-by slot:pe=${TEST_NTHREADS} --np ${TEST_NPROC} ${BIN_DIRECTORY}/${EXECNAME}";
   fi

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
   MAKESTRING="make EXEC=${EXECNAME} MPI=${MPI} GCOV=${GCOV} ${FLAGS} OPENMP=${TEST_OPENMP}";
   # if [ ${MPI} -eq 1 ]; then
   #    MAKESTRING="${MAKESTRING} -j ${NPROC}";
   # fi
   if $VERBOSE ; then
      $MAKESTRING 2>&1 | tee -a $LOGFILE;
   else
      $MAKESTRING >> $LOGFILE 2>&1;
   fi

   # Run test
   cd ${TEST_DIRECTORY}/${testname[n]};
   $DELETE_RESULTS;

   if $VERBOSE ; then
      function run_before_test
      {
         (${SHELL} ${BEFORETEST} 2>&1 | tee -a ${LOGFILE})
      }
      function run_test
      {
         (${RUN_TEST_BASE}${ndim}d ${rawname[i]}.nml 2>&1 | tee -a ${LOGFILE})
      }
   else
      function run_before_test
      {
         (${SHELL} ${BEFORETEST} >> ${LOGFILE} 2>&1)
      }
      function run_test
      {
         (${RUN_TEST_BASE}${ndim}d ${rawname[i]}.nml >> ${LOGFILE} 2>&1)
      }
   fi
   echo "Running test:" | tee -a $LOGFILE;
   STARTTIME_TEST=$(python3 -c 'import time; print(int(time.time()*1000))');

   if [ -f ${BEFORETEST} ]; then
         run_before_test;
   fi

   if $DO_RESTART ; then
      echo  "Restart: step 1 ..." | tee -a $LOGFILE;
      python3 ../../run_with_restart.py -s 1 -t ${rawname[i]}  | tee -a $LOGFILE;
      run_test;
      echo  "Restart: step 2 ..." | tee -a $LOGFILE;
      python3 ../../run_with_restart.py -s 2 -t ${rawname[i]}  | tee -a $LOGFILE;
      run_test;
      echo  "Restart: step 3 ..." | tee -a $LOGFILE;
      python3 ../../run_with_restart.py -s 3 -t ${rawname[i]}  | tee -a $LOGFILE;
   else
      run_test;
   fi

   # Record test time
   ENDTIME_TEST=$(python3 -c 'import time; print(int(time.time()*1000))');
   milliseconds=$((${ENDTIME_TEST} - ${STARTTIME_TEST}));
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
      echo "$status";
   fi
   echo "$status" >> $LOGFILE;

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

   # move coverage files to test dir
   if ${COVERAGE} ; then
      $RETURN_TO_BIN;
      gcov *.gcno > coverage_stats.txt
      cd -
      mv ${BIN_DIRECTORY}/*.gc* .
   fi
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
echo "\caption*{Test run summary using ${NPROC} processes with up to ${NTHREADS} threads,\\\\except where a test fixes its own NPROC or does not use OpenMP}" >> $latexfile;
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

# Compile latex file
# Use pdflatex if available, otherwise use tectonic
if command -v pdflatex &> /dev/null; then
   echo "Using pdflatex to compile the pdf document" | tee -a $LOGFILE;
   LATEX_COMPILER="pdflatex";
else
   echo "Using tectonic to compile the pdf document" | tee -a $LOGFILE;
   LATEX_COMPILER="tectonic";
fi
if $VERBOSE ; then
   $LATEX_COMPILER $latexfile 2>&1 | tee -a $LOGFILE;
   $LATEX_COMPILER $latexfile 2>&1 | tee -a $LOGFILE;
else
   $LATEX_COMPILER $latexfile >> $LOGFILE 2>&1;
   $LATEX_COMPILER $latexfile >> $LOGFILE 2>&1;
fi
rm ${latexfile/.tex/.log};
rm ${latexfile/.tex/.aux};
rm ${latexfile/.tex/.out};
rm $latexfile;

#######################################################################
# Generate total coverage data
#######################################################################
if ${COVERAGE} ; then
   rm -r coverage
   ALL_TEST_DIRS=""
   for ((i=0;i<$ntests;i++)); do
      n=${testnum[i]};
      test_dir_name=${TEST_DIRECTORY}/${testname[n]};
      ALL_TEST_DIRS="${ALL_TEST_DIRS} ${test_dir_name}"
   done
   mkdir coverage
   python3 multi_gcov_aggregator.py ${ALL_TEST_DIRS} coverage
fi

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

if $all_tests_ok ; then
   exit;
else
   exit 1;
fi
