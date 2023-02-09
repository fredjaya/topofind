# Simulations  
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

**test1**  
[2_mast_A_B]    BIC: 81891.4428
[3_mast_B_AA_AB]        No new splits, MAST not run. (All sites assigned to AA)
[3_mast_A_BA_BB]        BIC: 81972.3383

```  
OrderedDict([
	('A', '/home/frederickjaya/Dropbox/treemix_rc/04_testing/python/test1/2_split_A_B/test1_class_1-out.treefile'), 
	('B', '/home/frederickjaya/Dropbox/treemix_rc/04_testing/python/test1/2_split_A_B/test1_class_2-out.treefile'), 
	('AA', '/home/frederickjaya/Dropbox/treemix_rc/04_testing/python/test1/3_split_AA_AB/test1_class_1-out_class_1-out.treefile'), 
	('AB', None), 
	('BA', '/home/frederickjaya/Dropbox/treemix_rc/04_testing/python/test1/3_split_BA_BB/test1_class_2-out_class_1-out.treefile'),
	('BB', '/home/frederickjaya/Dropbox/treemix_rc/04_testing/python/test1/3_split_BA_BB/test1_class_2-out_class_2-out.treefile')
])

OrderedDict([
	('2_mast_A_B', {'bic': 81891.4428, 'input_trees': ['A', 'B'], 'aln': {'A′': '/home/frederickjaya/Dropbox/treemix_rc/04_testing/python/test1/2_mast_A_B/test1_class_1-out.fas', 'B′': '/home/frederickjaya/Dropbox/treemix_rc/04_testing/python/test1/2_mast_A_B/test1_class_2-out.fas'}}),
	('3_mast_B_AA_AB', None),
	('3_mast_A_BA_BB', {'bic': 81972.3383, 'input_trees': ['A', 'BA', 'BB'], 'aln': {'A′': '/home/frederickjaya/Dropbox/treemix_rc/04_testing/python/test1/3_mast_A_BA_BB/test1_class_1-out.fas', 'BA′': '/home/frederickjaya/Dropbox/treemix_rc/04_testing/python/test1/3_mast_A_BA_BB/test1_class_2-out.fas', 'BB′': '/home/frederickjaya/Dropbox/treemix_rc/04_testing/python/test1/3_mast_A_BA_BB/test1_class_3-out.fas'}})])
```  

**test2**  
Only one state is observed in alignment  
[2_mast_A_B]    No new partitions, MAST not run.  

**test3**  
All sites assigned to the one tree.   
[2_mast_A_B]    No new partitions, MAST not run.  


**concat**  


**orchid set2**  
|                           |              | ID     | Parent |
| ------------------------- | ------------ | ------ | ------ |
| 2_mast_A_B                | 2004337.8210 | T2     |        |
| 3_mast_A_BA_BB            | 2002991.7827 | T3Best |        |
| 3_mast_B_AA_AB            | 2003296.2824 |        |        |
| 4_mast_AA_AB_BA_BB        | 2002380.9617 |        |        |
| 4_mast_A_BA_BBA_BBB       | 2002293.2420 | T4Best | T3Best |
| 4_mast_A_BB_BAA_BAB       | 2002660.2608 |        |        |
| 4_mast_B_AA_ABA_ABB       | 2002606.5423 |        |        |
| 4_mast_B_AB_AAA_AAB       | 2002797.4753 |        |        |
| 4_mast_BA_BB_AA_AB        | 2002384.0475 |        |        |
| 5_mast_AA_ABA_ABB_BA_BB   | 2001934.3887 |        |        |
| 5_mast_AA_AB_BA_BBA_BBB   | 2002223.2472 |        |        |
| 5_mast_AA_AB_BB_BAA_BAB   | 2001855.2690 |        |        |
| 5_mast_AA_BA_BB_ABA_ABB   | 2001936.8025 |        |        |
| 5_mast_AB_AAA_AAB_BA_BB   | 2002019.3600 |        |        |
| 5_mast_A_BAA_BAB_BBA_BBB  | 2002045.6453 |        |        |
| 5_mast_AB_BA_BB_AAA_AAB   | 2002023.5144 |        |        |
| 5_mast_A_BB_BAA_BABA_BABB | 2001865.4787 |        |        |
| 5_mast_A_BB_BAB_BAAA_BAAB | 2001579.0009 |        |        |
| 5_mast_B_AAA_AAB_ABA_ABB  | 2001989.9314 |        |        |
| 5_mast_BA_AA_AB_BBA_BBB   | 2001938.0071 |        |        |
| 5_mast_B_AA_ABA_ABBA_ABBB | 2001954.9565 |        |        |
| 5_mast_B_AA_ABB_ABAA_ABAB | 2001807.2188 |        |        |
| 5_mast_B_AB_AAA_AABA_AABB | 2001949.4601 |        |        |
| 5_mast_B_AB_AAB_AAAA_AAAB | 2001728.9198 |        |        |
| 5_mast_B_ABA_ABB_AAA_AAB  | 2002304.2706 |        |        |
| 5_mast_BA_BB_AA_ABA_ABB   | 2001937.7858 |        |        |
| 5_mast_BA_BB_AB_AAA_AAB   | 2001779.0554 |        |        |
| 5_mast_BA_BBA_BBB_AA_AB   | 2001984.5023 |        |        |
| 5_mast_BB_AA_AB_BAA_BAB   | 2001867.1141 |        |        |
| 5_mast_BB_BAA_BAB_AA_AB   | 2001837.2068 |        |        |


