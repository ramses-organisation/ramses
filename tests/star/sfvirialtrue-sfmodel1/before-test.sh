if [ ! -d "dwarf_galaxy" ]; then
    echo "dwarf_galaxy ICs directory does not exist, dowloading it..."
    wget --timeout=10 --tries=3 --content-disposition https://ramses.cnrs.fr/wp-content/uploads/2026/02/dwarf_galaxy.tar  --no-check-certificate
    tar -xvf dwarf_galaxy.tar
    rm dwarf_galaxy.tar

    echo "dwarf_galaxy ICs downloaded."
else
    echo "dwarf_galaxy ICs  directory already exists."
fi
