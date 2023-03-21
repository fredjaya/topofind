# TopoFind
Finding tree topologies from sequence alignments.

## Dependencies  
iqtree 2.2.3.hmmster  
R 4.1.2  
devtools 2.4.3  
[MixtureModelHMM](https://github.com/fredjaya/MixtureModelHMM)  
[AMAS](https://github.com/marekborowiec/AMAS)  

## Installation  

```
git clone https://github.com/fredjaya/topofind.git
cd topofind
conda env create -f env.yml
conda activate topofind
Rscript install_pacakges.R
```

- Install IQ-TREE manually and add symlink to /bin  

## Usage  

- Infer the (best number of) topologies from a sequence alignment (`test1.fa`):  
```
TopoFind.py -a data/test1.fa
```  

- Show all available options:  
```
TopoFind.py -h  
```

- Run unit tests (for my own reference):  
```
python3 -m run.py -h  
```
