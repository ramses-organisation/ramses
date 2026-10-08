if [ -d "cosmo_64" ]; then
    echo "cosmo_64 directory already exists."
else
    # Reuse the 128^3 ICs from cosmo test if present
    ICDIR="../cosmo/cosmo_128"
    if [ ! -d "${ICDIR}" ]; then
        echo "cosmo_128 directory does not exist, downloading it..."
        wget --timeout=10 --tries=3 --no-check-certificate https://ramses.cnrs.fr/wp-content/uploads/2025/01/cosmo_128.zip
        unzip cosmo_128.zip
        rm cosmo_128.zip
        ICDIR="cosmo_128"
        echo "cosmo_128 directory downloaded."
    fi

    # Degrade to level 6
    echo "Degrading ${ICDIR} to cosmo_64..."
    ${FC:-gfortran} -O2 -o degrade_grafic ../../../utils/f90/degrade_grafic.f90
    mkdir -p cosmo_64
    ./degrade_grafic ${ICDIR} cosmo_64
    rm -f degrade_grafic
fi
