## regNet: a python analysis for motif influence on gene expression

To determine the influence of motif activity on gene expression, we compare usual Ridge Regression with a Bayesian Linear Mixed Model (https://github.com/limix/limix)

# Technical details
    - code is written for python2.7
    - entire analysis can be run via command line with bash scripts
    - conda environment provided

## CONDA-environment	
You can create a conda environment with the provided YAML file:
```bash
conda env create --name NAME --file cmapPy_conda.yml
```
you can then use the environment with 
```bash 
source activate NAME
```
and when you're done using it, 
```bash 
source deactivate NAME
```
For more information, please go to 	[conda project](https://docs.conda.io/projects/conda/en/latest/user-guide/index.html).
	

## SIMULATION OF DATA
We generate gene expression data where the gene expression is determined as linear relationship between motif scores and noise:
    $$\mathbf{Y}_{GENES, COND} = \mathbf{M}_{GENES, TF} \beta_{TF, COND} + \text{noise}$$
    
### Features:
    - play with degree of signal between motifs $\mathbf{M}\beta$ and noise
    - play with degree of structure in the noise

### USAGE
*  $REP:  integer   
    generate $REP randomly drawn sets of gene expression data.

* $COND:    integer   
    number of conditions/samples
    
* $GENES:   integer   
    number of genes/peaks
    
*  Motif file for chosen genes/peaks must be provided in data/motif/ and can be set in 
            code/01_simulation/datageneration/motif.py \\\
   and can be performed with gimmemotifs (https://github.com/vanheeringen-lab/gimmemotifs)

* $TF:  integer   
    number of Transcription Factors

*  $FRACTION: value between 0 and 1   
    control degree of signal in gene expression from motifs or noise (FRACTION) .


* $SIGMA: string of 
   "random", "randomV", "identity", "identityV"    
    Specify shape of generated noise: 
            where "randomV" and "identityV" exhibit structure similar to covariance 
            between motif scores V. 
*  $NOISEFRACTION: floar between 0 and 1   
intensity of structure in noise.
   
Generated data is then stored in './data/simuliaton/' and the printed output of the console saved in './data/stats/' 
        
 ```bash
 ./scripts/dataGeneration.sh -r ${REP} -c ${COND} -G ${GENES} -T ${TF}  -S ${SIGMA} -f ${FRACTION} -N ${NOISEFRACTION}
 ```

The script automatically generates datasets for a covariance of shape:
- identity (independence assumption)
- random  (no-pre-defined) correlation assumption
- lowrank_2 (block matrix with 2 block matrices along the diagonal)
- lowrank_(0.5\*$COND) (block matrix with 0.5\*$COND blocks along the diagonal)

## RUN LIMIX

Compute Limix and Ridge Regression on simulated data


### FEATURES
    -run limix on generated dataset
    -set noise-structure to be fit in limix
 

### USAGE
* $COND:    integer   
    number of conditions
    
* $GENES:   integer   
    number of genes   
    Motif file for chosen genes/peaks must be provided in data/motif/ and can be set in 
            code/01_simulation/datageneration/motif.py \\\
   and can be performed with [gimmemotif]{https://
github.com/vanheeringen-lab/gimmemotifs}

* $TF:  integer   
    number of Transcription Factors
    
* $SIGMA: string of 
   "random", "randomV", "identity", "identityV"    
    Specify shape of generated noise: 
            where "randomV" and "identityV" exhibit structure similar to covariance 
    shape of Sigma (of generated data)

* $NOISEFRACTION:   float between 0 and 1 
  intensity of structure in noise in generated data
   
* $NOISE:      string ('random', 'id', 'diag')   
    noise structure to be fit to data
    
* $GEN: string ('random', 'id', 'lowrank_2', 'lowrank_'$(0.5\*$COND))   
    covariance structure between samples (COND) used to generate data
    
* $ESTIM: string ('freeform', 'lowrank_$RANK', 'diagonal', 'block',...) (see limix for complete overview)

* $FRAC:    float (value between 0 and 1)   
    fraction of signal and noise ratio of generated data
    
* $INIT:    boolean   
    initialize limix with real data
  
* $PERTURN: boolean 
    perturbation in limix (see limix for more detail)
   
* $NUM:     integer   
    number of cores to use per replicate to compute limix
   
```bash
./scripts/compLimix_withParam.sh -c $COND -G $GENES -T $TF -r $REP -S $SIGMA	-N $NOISEFRACTION -E $NOISE	-g $GEN -e $ESTIM -f $FRAC
		-i $INIT -p $PERTURB -n $NUM 
```


### Exemplary use of how to use data
The generated data is stored as a cPickle-object. Those objects are generally loaded by
```python
import cPickle as pickle
import pandas as pd

with open(FILENAME, 'rb') as f:
    df = pickle.load(f)
```

The data itself is stored in a Ybetaparams class, with slots for 
-Y		the gene expression data (pandas dataframe)
-beta 		the motif influential weights (pandas dataframe)
-params		parameter set that was used to generate data (class params_dict)

To access and use data further, I suggest the following:
```python
### load packages
import Ybetaparam as Ybp
import cPickle as pickle

### FILENAMES
FILE_PATH = "data/simulation/"
'''#exemplary file name with 
-"C_4":		4 conditions/samples, 
-"G_978":	978 genes, 
-"TF_623":	623 motifs, 
-"V_lowrank_2":	assumend correlation between conditions a lowrank matrix of rank 2, 
-"rep_2":	2 repititions, 
-"Sigma_random":	a random noise matrix Sigma, 
-"R2method_mean_V_mean_M":	initialization strategy to control for the fraction in signal between motif influence and noise, and 
-"frac_20":	20% of the signal being explained by the motifs, the rest being noise
'''
FILE_NAME = "Ybetaparams_generated_C_4_G_978_TF_623_V_lowrank_2_rep_2_Sigma_random_R2method_mean_V_mean_M_frac_20.pkl"

# load Ybp object
with open(FILE_PATH + FILE_NAME, 'rb') as f:
    Ybp = pickle.load(f)

# repetitions are stored in columns
Y = Ybp.Y

### get parameters used for data generation
beta = Ybp.beta
params = Ybp.parameter

# motifs
motif = params.motif

# covariance structure used to generate data
covarV = params["V"]


# information about shape of data
Genes = params["G"]
Cond = params["C"]

### reshape Y
## get ith repetition of Y back into shape:
# set i
i = 0
# subset Y and reshape
Y_i = Y.iloc[:,i].reshape(Genes, Cond, order='F')




#### Your analysis comes here
.....
``` 
