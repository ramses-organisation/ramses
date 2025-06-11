#!/bin/bash

#TEST_NAME=${rawname[i]}

TEST_NAME=hydro-sod-tube

# Define possible values for each parameter
RIEMANN_METHODS=("'hllc'" "'hll'" "'acoustic'" "'exact'")
#SLOPE_TYPES=(1 2 3 4 5 6 7 8)
#INTERPOL_VARS=(0 1 2)
#INTERPOL_TYPES=(0 1 2 3 4)
SLOPE_TYPES=(1)
INTERPOL_VARS=(0)
INTERPOL_TYPES=(0)

# Path to the test executable and the Python script
TEST_EXECUTABLE="./run_test"
COMPARE_SCRIPT="./compare_results.py"

# Loop over all valid combinations
for RIEMANN in "${RIEMANN_METHODS[@]}"; do
    for SLOPE_TYPE in "${SLOPE_TYPES[@]}"; do
        for INTERPOL_VAR in "${INTERPOL_VARS[@]}"; do
            for INTERPOL_TYPE in "${INTERPOL_TYPES[@]}"; do
                # Ensure INTERPOL_TYPE=4 is only used when INTERPOL_VAR=2
                if [[ "$INTERPOL_TYPE" -eq 4 && "$INTERPOL_VAR" -ne 2 ]]; then
                    continue
                fi

                # Create a unique name for the test case
                #NAMELIST_FILE="${TEST_NAME}_${RIEMANN}_${SLOPE_TYPE}_${INTERPOL_VAR}_${INTERPOL_TYPE}.nml"
                NAMELIST_FILE="${TEST_NAME}.nml"

                # Generate the namelist file
                source namelist.sh
                                
                # Run the test with the generated namelist
                echo "parameters: ${RIEMANN} ${SLOPE_TYPE} ${INTERPOL_VAR} ${INTERPOL_TYPE}"
                ${RUN_TEST_BASE}${ndim}d $NAMELIST_FILE >> $LOGFILE 2>&1
                
                # Call the Python script for comparison
                #python3 "$COMPARE_SCRIPT" "$TEST_NAME"

                # Plot and analyse results
                #status=$(python3 check_solution.py 2>&1);
                #if $VERBOSE ; then
                #    echo $status;
                #fi
                #echo $status >> $LOGFILE;

            done
        done
    done
done

#python3 plot-$TEST_NAME.py