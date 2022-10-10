# rec-mast
Recombination detection with mixtures across sites and trees

## Dependencies  
nextflow 22.04.5  
[iqtree 2.2.0.7.mix](https://github.com/iqtree/iqtree2/releases/tag/v2.2.0.7.mx)  
R 4.1.2  
devtools 2.4.3  
[MixtureModelHMM](https://github.com/fredjaya/MixtureModelHMM)  
[AMAS](https://github.com/marekborowiec/AMAS)  

## Simulations  

|       | Best $R_n$ | $BIC_{t1;r2}$ | $BIC_{MAST+TR}$ | $T_n$ | Partitions        | Simulation parameters                      | Outcome               |
| ----- | ---------- | ------------- | --------------- | ----- | ----------------- | ------------------------------------------ | --------------------- |
| test1 | JC+R2      | 83582.392     | **81891.4474**  | 2     | 1-3500; 3501-5000 | JC model, 3500 tree1, 1500 tree2           | PERFECT               |
| test2 | JC+I       | 77944.946     |                 | 1     |                   | JC model, 500 tree1, 4500 tree2            | **2T missed**         |
| test3 | JC         | 66541.139     |                 | 1     |                   | 37 taxa, GTR+I+G, 2000 tree1, 3000 tree2   | PERFECT               |
| test4 | F81+F+R5   | 276162.269    | **236542.8737** | 2     | 1-1999; 2000-5000 | 37 taxa, GTR+I+G, 2000 tree1, 3000 tree2   | Correct, but MF wrong |
| test5 | F81+F+R4   | 216576.703    |                 | 1     |                   | 37 taxa, GTR+I+G, 5000 tree1               | Correct, but MF wrong |
| test6 | F81+F+R4   | 216576.703    |                 | 1     |                   | 10 taxa, GTR+I+R12, 5000 tree1             | Correct, but MF wrong |
| test7 | F81+F+R3   | 68415.617     | **67845.8588**  | 2     | 1-3995; 3995-5000 | 10 taxa, GTR+I+R12, 4000 tree1, 1000 tree2 | **Off by 5 bases!**   |
| test8 | JC+R2      | 75814.794     | **70325.4333**  | 2     | 1-2500; 2501-5000 | 10 taxa, JC, 2500 tree1, 2500 tree2        | PERFECT               |


