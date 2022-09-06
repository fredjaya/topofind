# Prototyping and development using data from Lytras et al. (2022)  

## Methods  
Construct a single unconstrained tree estimating the best (model/) number of rates to delimit trees. Output site RHAS likelihoods and posterior probabilities:  
```
iqtree-2.2.0.7.mix-Linux/bin/iqtree2 -s ../2207_lytras/78sarbeco_wgal_fn_gisaid.fas -pre rates_i0 -alninfo -wslr -wspr -nt AUTO
```

Run `hmm.Rmd` to assign sites to a rate class, and generate a partition file of class boundaries. Then, run loci trees according to partitions.  
```
iqtree-2.2.0.7.mix-Linux/bin/iqtree2 -s ../2207_lytras/78sarbeco_wgal_fn_gisaid.fas -S rates_i0.nex -pre rates_i0_hmmclasses -alninfo -wslr -wspr -nt AUTO
```

Run MAST using +T and +TR with HMM trees:
```
iqtree-2.2.0.7.mix-Linux/bin/iqtree2 -s ../2207_lytras/78sarbeco_wgal_fn_gisaid.fas -te rates_i0_hmmclasses.treefile -pre rates_i0_mast -m"TMIX{GTR+FO+G,GTR+FO+G,GTR+FO+G,GTR+FO+G}+T" -nt AUTO
iqtree-2.2.0.7.mix-Linux/bin/iqtree2 -s ../2207_lytras/78sarbeco_wgal_fn_gisaid.fas -te rates_i0_hmmclasses.treefile -pre rates_i0_mast -m"TMIX{GTR+FO+G,GTR+FO+G,GTR+FO+G,GTR+FO+G}+TR" -nt AUTO
```

Run best partition model using HMM class boundaries to compare with MAST:  
```
iqtree-2.2.0.7.mix-Linux/bin/iqtree2 -s ../2207_lytras/78sarbeco_wgal_fn_gisaid.fas -p rates_i0_rates.nex -m MFP+MERGE -pre rates_i0_partition -alninfo -wslmr -wspmr -nt AUTO
```

## Results  

Best fitting whole-alignment single-tree model is GTR+I+R4.  


### BIC scores  

| Step                | BIC         |
| --------------      | ----------- |
| Single tree         | 699483.2846 |
| Post-HMM RHAS	trees | 670486.8827 |
| Partition           | 693574.6871 |
| MAST +T             | 669540.3526 |
| MAST +TR            | 672496.2097 |
