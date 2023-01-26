# rec-mast
Recombination detection with mixtures across sites and trees

## Dependencies  
nextflow 22.04.5  
[iqtree 2.2.0.7.mix](https://github.com/iqtree/iqtree2/releases/tag/v2.2.0.7.mx)  
R 4.1.2  
devtools 2.4.3  
[MixtureModelHMM](https://github.com/fredjaya/MixtureModelHMM)  
[AMAS](https://github.com/marekborowiec/AMAS)  

## Installation  

```
conda env create -f env.yml
conda activate rec-mast
Rscript install_pacakges.R
```

- Install IQ-TREE manually and add symlink to /bin  

## Usage
```
# call python3 to ensure conda version is used
python3 run.py -h  

usage: run.py [-h] -a ALN -f {fasta,phylip} [-d OUTPUT_DIR] [-T NUM_THREADS] [-e {local,slurm}]

options:
  -h, --help            show this help message and exit
  -a ALN, --aln ALN     Input sequence alignment
  -f {fasta,phylip}, --aln_format {fasta,phylip}
                        Alignment format
  -d OUTPUT_DIR, --output_dir OUTPUT_DIR
                        Path to output directory
  -T NUM_THREADS, --num_threads NUM_THREADS
                        [IQ-TREE2] No. cores/threads or AUTO-detect (default: 1)
  -e {local,slurm}, --executor {local,slurm}
                        [NEXTFLOW] Where to run nextflow processes (default: local)

```

Example:
```
python3 run.py -a data/test1.phy -f fasta
```
