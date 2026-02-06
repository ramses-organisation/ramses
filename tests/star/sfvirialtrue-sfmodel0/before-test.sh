if [ ! -d "dwarf_galaxy" ]; then
    echo "dwarf_galaxy ICs directory does not exist, dowloading it..."
    wget --timeout=10 --tries=3 --content-disposition https://seafile.unistra.fr/f/b294343dac9c4ce7b22e/?dl=1 --no-check-certificate
    tar -xvf dwarf_galaxy.tar
    rm dwarf_galaxy.tar

    echo "dwarf_galaxy ICs downloaded."
else
    echo "dwarf_galaxy ICs  directory already exists."
fi
