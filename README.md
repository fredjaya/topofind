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
git clone https://github.com/fredjaya/rec-mast.git
cd rec-mast
conda env create -f env.yml
conda activate rec-mast
Rscript install_pacakges.R
```

- Install IQ-TREE manually and add symlink to /bin  

## Usage  
Call python3 to ensure conda version is used.  

- Infer the (best number of) topologies from a sequence alignment (`test1.fa`):  
```
python3 run.py -a data/test1.fa -f fasta
```  

- Show all available options:  
```
python3 run.py -h  
```

- Run unit tests (for my own reference):  
```
python3 -m run.py -h  
```

## Workflow  
```mermaid
flowchart TD
	INPUT ==> split_aln
	split_aln ==> |Update| PartitionedTrees
	PartitionedTrees ==> mast
	mast ==> |Update| MastResults
	MastResults ==> bic_improving
	bic_improving ==> |Y| get_new_trees
	bic_improving --> |N| Stop
	PartitionedTrees -.-> get_new_trees
	get_new_trees --> |Update| PartitionedTrees
	get_new_trees --> split_aln

	%% Subgraphs
	subgraph split_aln
	t1_iqtree_per_split --> |test2| 2a[Error: Only one state is observed in alignment]
	2a --> 2b["Outcome: tree_files = []"]
	end
	
	subgraph PartitionedTrees
	tree_names
	tree_files
	end
	
	subgraph mast
	concat_trees --> run_nf
	end
	
	subgraph MastResults
	run_name
	bic
	input_trees
	aln
	end
```

## Outputs  

## Issues and contributing  

