# STO-MS
quantitative profiling of STOichiometry by Mass Shift

## tutorial
```
cd examples/
python2 ../scripts/pcp_data_paraser.py ../database/uniprot_sprot_HUMAN_160504.fasta > results.txt
python2 ../scripts/processing_data.py results.txt ../database/molecular_weight.txt > data.txt
python3 ../scripts/find_the_stoichiometry_DMSO.py data.txt ../database/molecular_weight.txt DMSO > stoichiometry_DMSO.txt
cd ..; sort */ref_list_DMSO.txt | sort | uniq > database/ref_list.txt
python3 ../scripts/find_the_stoichiometry_HNE.py data.txt ../database/molecular_weight.txt HNE ../database/ref_list.txt > stoichiometry_HNE.txt
