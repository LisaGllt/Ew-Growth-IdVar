#!/bin/bash

cd $BRIDGE_MSUB_PWD		# $BRIDGE_MSUB_PWD est une variable d'environnement representant le repertoire de soumission
set -x   

module load r

#MSUB -r VDT_Part2	            # Nom du job
#MSUB -o Bootstrap_VI_Part2.out			# Sortie standard
#MSUB -e Bootstrap_VI_Part2.err			# Sortie d'erreur
#MSUB -@ remy.beaudouin@ineris.fr:end	# envoie un mail a l'adresse indiquee en fin de job 

#MSUB -T 86400  # maximum walltime of the batch job in seconds
#MSUB -n 1  # 1 processus
#MSUB -c 28 # 28 coeurs par processus
#MSUB -N 1		# number of NODES (CPU) per task (default=1)
#MSUB -q milan	# type of node : choosing standard nodes (see ccc_mpinfo)
#MSUB -m scratch,work

cd ${BRIDGE_MSUB_PWD}

ccc_mprun Rscript --vanilla Bootstrap_VI_Part2.R

