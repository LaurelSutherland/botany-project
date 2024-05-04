# botany-project

|Software|Description|Strengths|Weaknesses|Assumptions|User Choices|
|--------|-----------|---------|----------|-----------|------------|
|EVOLVER|models gene evolution|simulates full-sized multichromosome genome evolution|only haploid genome and lacks a population model|produces reference-like genomes||
|--------|-----------|---------|----------|-----------|------------|
|PSAR|alignment tool that compares simulation results to an original alignment|calculates a PSAR pair score similar to probability of accuracy|only allows substitutions, insertions, deletions|PSAR pair score is a good predictor of accuracy||
|--------|-----------|---------|----------|-----------|------------|
|Clustalw (MSA)|builds a sequence alignment by adding one sequence at a time|able to change weights so it is sensitive to difference in diverse sequences, fast, works with whole genome or protein sequences|very sensitive to any errors during the alignment process and a single misalignment will alter the final tree, does not use all available data for each sequence|all alignments already made are correct||
|--------|-----------|---------|----------|-----------|------------|
|T-Coffee (MSA)|uses a progressive alignment strategy to consider alignments between all pairs from a library of alignments|builds on ClustalW and other programs to create a more accurate tree|slower than ClustalW, cannot remove preexisting gaps, overweights small segments|this weighting scheme tolerates noise from small similar segments arising due to chance better|This is my preferred method. Also, able to run mcoffee to compare multiple methods at once.|
|--------|-----------|---------|----------|-----------|------------|
|MUSCLE (MSA)|calculates two distance measures for a pair of sequences and cluster distance matrices to create phylogeny|fastest and most accurate|not consistently the most accurate, have to look at results carefully|the program will find the global optimal tree each time||

