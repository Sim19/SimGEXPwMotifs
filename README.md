## regNet: a python analysis for motif influence on gene expression

To determine the influence of motif activity on gene expression, we compare usual Ridge Regression with a Bayesian Linear Mixed Model (https://github.com/limix/limix)

# Technical details
    - code is written for python2.7
    - entire analysis can be run via command line with bash scripts

# CONDA-environment	
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

