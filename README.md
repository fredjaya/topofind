# rec-mast
Recombination detection with mixtures across sites and trees

## Dependencies  
iqtree 2.2.0.7.mix  
R 4.1.2  
devtools 2.4.3  
[MixtureModelHMM](https://github.com/roblanf/MixtureModelHMM)  

## Workflow 

Currently, use RHAS categories to delimit sites for input MAST topologies:  

```mermaid  
flowchart TD
	0[MSA] --> 1["1. Estimate the best-fitting (RHAS) model for a single unconstrained tree; output site lnLs"] --> 2[2. Run HMM to assign sites according to RHAS classes and save partitioning scheme] --> 3a[3a. Construct individual trees for each HMM class partition] & 3b[3b. Run merged partition model for comparison]
	3a --> 4a[4a. Check for and remove identical topologies] & 4b[4b. What happens when all topologies are identical? hehe]
	4a --> 5a[5a. Run MAST using each unique topology with +T] & 5b[5b. And with +TR]
	1 & 3a & 3b & 5a & 5b -.- 6[6. Assess and compare BICs] --> 7[7. Do either MAST models show an improvement in BIC?] --> |No| STOP
	7 --> |Yes| 8b[8b. Repeat process on individual HMM class partitions] --> 2
```

## Methods to delimit trees  
- Rate heterogeneity classes  
- Site likelihoods  
- sCF  

## Testing and development  

Currently focusing on building a working MVP, utilising the 78 sarbecovirus alignment from Lytras et al. (2022)
