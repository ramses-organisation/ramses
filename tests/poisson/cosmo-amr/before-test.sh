if [ ! -d "cosmo_128" ]; then
    if [ ! -d "../cosmo/cosmo_128" ]; then
        mv ../cosmo/cosmo_128 ./
        echo "cosmo_128 directory moved from ../cosmo/cosmo_128."
    else
        echo "cosmo_128 directory does not exist, dowloading it..."
        wget --timeout=10 --tries=3 --no-check-certificate https://ramses.cnrs.fr/wp-content/uploads/2025/01/cosmo_128.zip
        unzip cosmo_128.zip
        rm cosmo_128.zip

        echo "cosmo_128 directory downloaded."
    fi
else
    echo "cosmo_128 directory already exists."
fi
