#!/bin/bash
#MSUB -r tstranex1           # Nom du job                
#MSUB -n 16                  # Reservation de 4 processus au total
#MSUB -N 8                  # Les processus sont répartis sur 2 noeuds
#MSUB -T 36000                # Limite de temps elapsed du job ici 600s      
#MSUB -o stdout_test7ramses          # Sortie standard
#MSUB -e stderr_test7ramses          # Sortie d'erreur       
#MSUB -p gen2191      # Allocation
#MSUB -q gpu                # sélection de la queue GPU (groupe genci ou challenge uniquement)

set -x

cd $SCRATCHDIR/tstranex/test7/ramses

mpirun $BRIDGE_MSUB_PWD/ramses3d_otsa $BRIDGE_MSUB_PWD/test7_ramses.nml