if [ ! -d "ICs_lowres" ]; then
    echo "ICs_lowres directory does not exist. Generating initial conditions..."
    mkdir ICs_lowres
    cd ../../../utils/ic/galaxy_ic/Galic/
    make clean
    make
    ./galic_MakeDiskGalaxy ../../../../tests/poisson/galaxy/ICs_lowres/
    cd -
    echo "ICs generated."
else
    echo "ICs_lowres directory already exists."
fi