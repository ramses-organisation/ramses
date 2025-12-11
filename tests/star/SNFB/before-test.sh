if [ ! -d "cosmo_128" ]; then
    echo "cosmo_128 directory does not exist, dowloading it..."
    wget --timeout=10 --tries=3 --no-check-certificate https://ramses.cnrs.fr/wp-content/uploads/2025/01/cosmo_128.zip
    unzip cosmo_128.zip
    rm cosmo_128.zip

    echo "cosmo_128 directory downloaded."
else
    echo "cosmo_128 directory already exists."
fi
