#!/bin/bash
# Bash script to run the data generation on the server


while getopts ":c:G:r:S:f:N:T:" opt;
do
	case $opt in 
		c) COND="$OPTARG"
			;;
		G) GENES="$OPTARG"
			;;
		T) TF="$OPTARG"
			;;
		r) REP="$OPTARG"
			;;
		S) SIGMA="$OPTARG"
			;;
		f) FRAC="$OPTARG"
			;;
		N) FRACNOISE="$OPTARG"
			;;
		\?) echo "Invalid option -$OPTARG" >&2
			;;
	esac
done


R2METHOD=mean_V_mean_M

## add data and data/stats/ folders if not existent
if [! -d data/stats/]; then
  	mkdir -p data/stats/
fi


FILENAME=data/stats/output_dataGeneration_C_${COND}_G_${GENES}_TF_${TF}_rep_${REP}_Sigma_${SIGMA}_FracNoise_${FRACNOISE}_Rsquared_${R2METHOD}_FRAC_${FRAC}.txt
echo Fit generated data with 
echo 	${COND} conditions, 
echo 	${GENES} genes, 
echo 	${TF} TFs, 
echo 	${REP} repititions and 
echo 	Sigma of the form ${SIGMA} with 
echo 	${FRAC}% variation explained.
echo 	${FRACNOISE}% variation in noise.

# python command:
time -p python generate_data.py -c ${COND} -G ${GENES} -T ${TF} -r ${REP} -S ${SIGMA} -R ${R2METHOD} -f ${FRAC} --fN ${FRACNOISE} |&tee -a ${FILENAME}
