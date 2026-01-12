if [ ! -d "ic" ]; then
   # go to the root directory to avoid python import problems
   cd ../../..

   # Generate initial conditions
   echo "Generating ICs..."
   python3 tests/hydro/decaying-turbulence/initial_conditions.py

   # return to the test directory
   cd tests/hydro/decaying-turbulence
fi
