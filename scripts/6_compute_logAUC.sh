#!/bin/bash

#here is the plan: 


#making directories
#run this in the compute_logAUC trajectory
#create the indentity lists
#create properly formated scores and name lists

#run the logAUC script

#setting variables
SCRIPT_DIR=/home/limcaoco/turbo/limcaoco/diffdock_v_screen/scripts/
PROJ_DIR=/home/limcaoco/turbo/limcaoco/diffdock_v_screen

for RECEPTOR in ampcfep;
do
    mkdir -p ${RECEPTOR}/results
    #mkdir -p ${RECEPTOR}/results/smina_relax
    
    cd ${RECEPTOR}
    echo "$(pwd)"
#    for SCORE_TYPE in DOCK_score l1_norm 12_norm l3_norm l4_norm l5_norm l6_norm l7_norm l8_norm confidence_score_mean confidence_score_max average_rank smina_affinity; do
    for SCORE_TYPE in fep; do
	mkdir -p results/${SCORE_TYPE}
	cd results/${SCORE_TYPE}
	echo "$(pwd)"
	echo "now working on ${SCORE_TYPE} of ${RECEPTOR}"
	python ${SCRIPT_DIR}/bootstrap_tldr_full_metrics_fixed.py -l ../../data/${RECEPTOR}_ligands.txt -d ../../data/${RECEPTOR}_decoys.txt -s1 ../../${SCORE_TYPE}.txt -b -m all > ${RECEPTOR}_${SCORE_TYPE}_tldr.log
	cd ${PROJ_DIR}/analysis/${RECEPTOR} 
    done
  #  cd ${PROJ_DIR}/analysis
done


