#Reproducing "One thousand plant transcriptomes and the phylogenomics of green plants" on a Small Gene Subset

##Dataset
The data for the study will include a subset of the genes from the original “One thousand plant transcriptomes and the phylogenomics of green plants” study (2019). The paper of the original study can be found at https://doi.org/10.1038/s41586-019-1693-2. The data used in the study was freely available through CyVerse Data Commons and could be downloaded from https://doi.org/10.25739/8m7t-4e85 in /oneKP_capstone_2019/alignments_and_trees/alignments/alignments-FAA-unmasked.tar.gz. The purpose of the study is to create two unique species trees using two different alignment methods with a smaller subset of the available genes and compare these trees to the tree provided in the original paper to examine the predicted evolutionary history of major clades within Viridiplantae. 

The original study utilized genetic information from 1124 species across 410 genes. To make my substudy more manageable, I utilized a smaller set of 8 genes: 4471, 4519, 4527, 4603, 4691, 4724, 4744, and 4757.
The downloaded genetic data is organized into folders assigned to each of the 410 genes, each of which contains a .fasta file of the pre-aligned transcriptomes for each species at that gene, both masked (meaning sequences with a high proportion of gaps were removed) and unmasked. I will be using the unmasked data for this study to allow my own phylogenetic pipelinen to manipulate the sequences and eliminate the chance that any similarities in results were due to filtering done before I obtained the data. All commands were run on Ubuntu v22.04.

#Extract the Data
After downloading tar.gz file from CyVerse Data Commons with the transcriptome data, it is saved to ~/botany-project/data where it needs to be unzipped and extracted. 
```
tar -xvzf alignments-FAA-unmasked.tar.gz
```
This produces a folder named genes_FAA which contains 410 subfolders for each of the genes included in the original study. This is the raw data we will be manipulating.



##Multiple-Sequence Alignment Methods
=====================================
To compare the efficiency of different alignment methods and their impact on the product of the phylogenetic pipeline, I will utilize two different MSA methods on the dataset. I selected MUSCLE and KALIGN.

#MUSCLE
MUSCLE uses a progressive alignment method guided by an initial pairwise alignment to match pairs that are most similar and build upon these pairs to create a full alignment. The algorithm calculates a Sum-of-Pairs (SP) scores which sums the scores of alignments over sequences and performs iterations to find the best score, which it concludes as the optimal tree. Some of the strengths of MUSCLE as a method for MSA include its ability to manipulate relatively large datasets rather fast and its high accuracy from the UPGMA tree construction and progressive alignment method it utilizes which make it adept at aligning homologous sequences. One of the weaknesses of MUSCLE is that early errors have a continuing impact on later alignments because of the progressive alignment method, which means MUSCLE may struggle if used on distantly related species. One of the assumptions made by the MUSCLE algorithm is that all taxa share a common ancestor. In a study such as this which examines such a wide array of species, it is very difficult to say whether that is true and may cause errors in the output. MUSCLE does not allow the user to decide the value of gap penalties, instead using the provided data to assign this value itself, which makes it user-friendly but may lead to less-than-optimal alignments.

Below is an example of the code used to download MUSCLE and align one of the genes.
1. Open the terminal and install MUSCLE using apt-get.
```
sudo apt-get update
sudo apt-install muscle
```
    It may ask you to approve the use of the amount of storage the installation of MUSCLE requires. If this happens, press "Y" on your keyboard. Once the download is complete, you will be able to execute MUSCLE by typing "muscle" in the command line of the terminal.
2. Run MUSCLE on gene 4471.
Below is the code used to produce an alignment with MUSCLE from the original .FASTA file provided in the open-source data from the original study. The same process will be repeated on all 8 genes utilized in this substudy.
```
cd ~/botany-project/data/genes_FAA/4471
muscle -in 4471.fasta -out 4471-muscle.afa

MUSCLE v3.8.1551 by Robert C. Edgar

http://www.drive5.com/muscle
This software is donated to the public domain.
Please cite: Edgar, R.C. Nucleic Acids Res 32(5), 1792-97.

4471 976 seqs, lengths min 326, max 984, avg 769
00:00:09    28 MB(-8%)  Iter   1  100.00%  K-mer dist pass 1
00:00:09    28 MB(-8%)  Iter   1  100.00%  K-mer dist pass 2
00:00:38  520 MB(-157%)  Iter   1  100.00%  Align node      
00:00:38  520 MB(-157%)  Iter   1  100.00%  Root alignment
00:01:14  523 MB(-158%)  Iter   2  100.00%  Refine tree   
00:01:14  523 MB(-158%)  Iter   2  100.00%  Root alignment
00:01:15  523 MB(-158%)  Iter   2  100.00%  Root alignment
00:19:35  523 MB(-158%)  Iter   3  100.00%  Refine biparts
00:33:42  523 MB(-158%)  Iter   4  100.00%  Refine biparts
00:45:45  523 MB(-158%)  Iter   5  100.00%  Refine biparts
00:52:07  523 MB(-158%)  Iter   6  100.00%  Refine biparts
00:52:07  523 MB(-158%)  Iter   6  100.00%  Refine biparts
```
3. Without changing directories, save the output file to a different folder with all of the alignments.
```
mv 4471-muscle.afa ~/botany-project/data/muscle/alignment/
```
This moves the output file from the folder for the individual gene where it was produced to a separate folder. Keeping all alignments produced using MUSCLE here will be helpful later.


#KALIGN
KALIGN utilizes the Wu-Manber string-matching algorithm which divides patterns into subpatterns and searches for those subpatterns within sequences in addition to a progressive alignment method similar than that used by MUSCLE. The strength of KALIGN is primarily that the Wu-Manber string-matching algorithm makes it significantly faster at aligning a large number of sequences, which is beneficial in a study such as this where other methods may be slowed by the scale of taxa and genes. Its high accuracy from the UPGMA tree construction and progressive alignment method it utilizes which make it adept at aligning homologous sequences. KALIGN is also capable of removing gaps from previously-aligned sequences, as is the case of the dataset used in this study, and realigning it using the Wu-Manber string-matching algorithm, which makes it an ideal choice when reproducing the work of others. Weaknesses of KALIGN include its previous difficulties associated with building guide trees using identical sequences. Many of the species are similar and therefore may cause conflict in tree-building. One assumption made by KALIGN is that pattern-based matching provides an acceptable accuracy and produces a correct guide tree, which has required a few version revisions to ensure. There is very little user choice required, and KALIGN is user-friendly with the algorithm determining the values used in progressive alignment.

Below is an example of the code used to download KALIGN and align one of the genes.
1. Open the terminal and install MUSCLE using apt-get.
```
sudo apt-get update
sudo apt -y install kalign
```
    It may ask you to approve the use of the amount of storage the installation of KALIGN requires. If this happens, press "Y" on your keyboard. Once the download is complete, you will be able to execute KALIGN by typing "kalign" in the command line of the terminal.
2. Run KALIGN on gene 4471.
Below is the code used to produce an alignment with KALIGN from the original .FASTA file provided in the open-source data from the original study. The same process will be repeated on all 8 genes utilized in this substudy.
```
cd ~/botany-project/data/genes_FAA/4471
kalign -in 4471.fasta -out 4471-kalign.afa

Kalign (3.3.1)

Copyright (C) 2006,2019,2020,2021 Timo Lassmann

This program comes with ABSOLUTELY NO WARRANTY; for details type:
`kalign -showw'.
This is free software, and you are welcome to redistribute it
under certain conditions; consult the COPYING file for details.

Please cite:
  Lassmann, Timo.
  "Kalign 3: multiple sequence alignment of large data sets."
  Bioinformatics (2019) 
  https://doi.org/10.1093/bioinformatics/btz795


WARNING: AVX2 instruction set not found!
         Kalign will not run optimally.

[2024-05-04 19:16:53] :     LOG : Detected DNA sequences.
[2024-05-04 19:16:53] :     LOG : CPU Time: 0.03u 00:00:00.02 Elapsed: 00:00:00.00
[2024-05-04 19:16:53] :     LOG : Detected: 976 sequences.
[2024-05-04 19:16:53] :     LOG : Calculating pairwise distances
[2024-05-04 19:16:54] :     LOG : CPU Time: 1.03u 00:00:01.02 Elapsed: 00:00:01.00
[2024-05-04 19:16:54] :     LOG : 32 anchors 
[2024-05-04 19:16:54] :     LOG : Building guide tree.
[2024-05-04 19:16:56] :     LOG : CPU Time: 2.17u 00:00:02.16 Elapsed: 00:00:02.00
[2024-05-04 19:16:56] :     LOG : Aligning
[2024-05-04 19:17:03] :     LOG : CPU Time: 11.93u 00:00:11.93 Elapsed: 00:00:07.00
```
3. Without changing directories, save the output file to a different folder with all of the alignments.
```
mv 4471-kalign.afa ~/botany-project/data/kalign/alignment/
```
This moves the output file from the folder for the individual gene where it was produced to a separate folder. Keeping all alignments produced using KALIGN here will be helpful later.



##Maximum Likelihood with IQ-TREE
================================
For the purpose of conducting the maximum likelihood comparison, I have selected IQ-TREE. IQ-Tree uses aligned gene sequences to assemble an optimal gene tree for each gene. Each gene tree will include a node for each of the species that contains that gene. Some of the strengths of IQ-TREE are that it has shown to return higher likelihood scores of the final gene tree than other maximum likelihood methods as a result of the hill-climbing approach and stochastic perturbation methods employed within the algorithm. Additionally, IQ-TREE automatically uses ModelFinder to select the best model of substitution, as well as providing bootstrap approximation for branch support, both of which help increase the accuracy. A weakness of IQ-TREE is that it can be very slow with large data sets, like the one used in this study. The IQ-TREE model assumes homogeneity, a constant rate of evolution that allows for reversibility, and constant nucleotide frequencies over time. User choices within the algorithm include the model that is used and the number of bootstrap iterations.

Below is an example of the code used to download IQ-TREE and find the maximum likelihood of one aligned gene. The example will be demonstrated on the MUSCLE alignment of gene 4471, but the same method will be used on the KALIGN alignment of the same gene and for the other genes.
1. Open the terminal and install IQ-TREE using apt-get.
```
sudo apt-get update
sudo apt -y install iqtree2
```
    It may ask you to approve the use of the amount of storage the installation of IQ-TREE requires. If this happens, press "Y" on your keyboard. Once the download is complete, you will be able to execute IQ-TREE by typing "iqtree2" in the command line of the terminal.
2. Run IQ-TREE on the MUSCLE alignment of gene 4471.
Below is the code used to find the maximum likelihood tree of gene 4471 using the MUSCLE-aligned sequences.
```
cd ~/botany-project/data/muscle/alignment/
iqtree2 -s 4471-muscle.afa 

IQ-TREE multicore version 2.0.7 for Linux 64-bit built Jan 21 2022
Developed by Bui Quang Minh, Nguyen Lam Tung, Olga Chernomor,
Heiko Schmidt, Dominik Schrempf, Michael Woodhams.

Host:    laurel-virtual-machine (AVX512, FMA3, 3 GB RAM)
Command: iqtree2 -s /home/laurel/Downloads/genes_FAA/4471/4471-muscle.fa
Seed:    561050 (Using SPRNG - Scalable Parallel Random Number Generator)
Time:    Sat Mar 30 10:30:49 2024
Kernel:  AVX+FMA - 1 threads (2 CPU cores detected)

HINT: Use -nt option to specify number of threads because your CPU has 2 cores!
HINT: -nt AUTO will automatically determine the best number of threads to use.

Reading alignment file /home/laurel/Downloads/genes_FAA/4471/4471-muscle.fa ... Fasta format detected
Alignment most likely contains DNA/RNA sequences
WARNING: 1 sites contain only gaps or ambiguous characters.
Alignment has 976 sequences with 1611 columns, 1484 distinct patterns
1199 parsimony-informative, 189 singleton sites, 223 constant sites
                     Gap/Ambiguity  Composition  p-value
   1  VRGZ                  76.91%    failed      0.14%
   2  YRMA                  70.70%    failed      0.43%
   3  RWXW                  70.83%    failed      2.07%
   4  WCLV                  43.14%    failed      0.01%
   5  NQYP                  44.01%    failed      0.00%
   6  XAXW                  78.52%    failed      0.85%
   7  Cyame_v1.0            40.04%    failed      0.00%
   8  PVGP                  44.75%    failed      2.02%
   9  OBUY                  68.22%    failed      1.70%
  10  NPRL                  60.77%    passed     39.91%
  11  OQWW                  70.33%    passed     28.80%
  12  YSBD                  56.05%    failed      0.00%
  13  OQON                  42.89%    failed      4.32%
  14  CQQP                  44.51%    failed      0.22%
  15  OAEZ                  43.02%    failed      0.08%
  16  Ostlu_v2.0            40.16%    failed      0.00%
  17  LSHT                  43.02%    failed      0.11%
  18  BOGT                  59.53%    passed     20.12%
  19  JQFK                  42.15%    passed     53.66%
  20  RTLC                  42.02%    failed      0.01%
  21  JMTE                  42.40%    failed      0.00%
  22  LLXJ                  42.27%    passed     10.42%
  23  RRSV                  42.89%    failed      2.80%
  24  PYDB                  42.89%    failed      1.36%
  25  URSB                  41.65%    failed      5.00%
  26  IKIZ                  42.64%    failed      1.76%
  27  ZJOJ                  43.39%    failed      1.95%
  28  DVYE                  70.20%    failed      0.67%
  29  JJZR                  42.02%    failed      2.15%
  30  KSFK                  43.76%    failed      0.00%
  31  HYHN                  42.40%    failed      0.00%
  32  KADG                  41.90%    passed      7.00%
  33  DZPJ                  41.40%    failed      0.00%
  34  UYFR                  44.88%    passed     10.18%
  35  NMAK                  41.96%    failed      0.43%
  36  LLEN                  42.64%    failed      0.93%
  37  RFAD                  42.64%    failed      0.93%
  38  VJED                  42.83%    failed      1.23%
  39  ALZF                  42.77%    failed      0.00%
  40  JIWJ                  43.64%    passed     25.63%
  41  GJIY                  41.28%    failed      0.02%
  42  BILC                  41.53%    failed      0.00%
  43  AJAU                  45.13%    failed      0.21%
  44  ZULJ                  43.02%    failed      0.00%
  45  Chabr_v0.1            39.04%    passed     17.89%
  46  Micpu_v3.0            40.29%    failed      0.00%
  47  PQED                  42.02%    passed     54.65%
  48  POOW                  58.16%    failed      4.02%
  49  QFND                  42.15%    failed      4.29%
  50  JKHA                  65.11%    passed     15.60%
  51  FXHG                  41.40%    passed     57.86%
  52  XOAL                  61.64%    failed      0.00%
  53  PRIQ                  41.15%    failed      0.00%
  54  QRTH                  41.03%    failed      0.00%
  55  ZFXU                  42.02%    failed      4.58%
  56  SYJM                  42.27%    passed     14.46%
  57  UKUC                  69.58%    failed      1.26%
  58  WDWX                  61.51%    passed      5.39%
  59  XKWQ                  41.15%    failed      4.05%
  60  VJDZ                  40.91%    failed      0.05%
  61  MWAN                  57.91%    failed      0.00%
  62  PZIF                  41.03%    failed      0.01%
  63  XDLL                  74.05%    failed      1.22%
  64  AOUJ                  40.91%    failed      0.13%
  65  AJUW                  41.03%    failed      0.00%
  66  LBRP                  41.90%    failed      0.00%
  67  JKKI                  41.03%    failed      0.00%
  68  ISPU                  41.15%    failed      0.00%
  69  Chlre_v5.5            40.78%    failed      0.00%
  70  Volca_v2.0            40.41%    failed      0.00%
  71  WRSL                  62.01%    failed      0.07%
  72  JWGT                  40.66%    failed      0.00%
  73  BAJW                  43.76%    failed      0.80%
  74  AZZW                  40.41%    failed      0.02%
  75  BAZF                  69.46%    failed      0.04%
  76  MFYC                  53.45%    failed      0.00%
  77  IJMT                  41.78%    passed      8.18%
  78  BZSH                  41.15%    failed      2.35%
  79  ZLQE                  42.15%    failed      0.00%
  80  XRTZ                  40.29%    failed      0.00%
  81  NSTT                  41.15%    failed      0.43%
  82  XMCL                  41.65%    failed      2.78%
  83  ZRMT                  42.52%    passed     39.52%
  84  MMKU                  41.15%    failed      0.05%
  85  AYPS                  41.65%    failed      1.65%
  86  KYIO                  41.40%    failed      0.08%
  87  DUMA                  75.42%    passed     16.95%
  88  QPDY                  62.51%    failed      0.04%
  89  VQBJ                  42.02%    failed      0.11%
  90  NKXU                  41.40%    failed      2.80%
  91  AKCR                  46.12%    failed      0.22%
  92  EEJO                  41.78%    failed      0.00%
  93  Chlva_v1.0            40.53%    failed      0.00%
  94  MNCB                  42.02%    failed      1.28%
  95  XIVI                  41.65%    passed      5.96%
  96  SNOX                  42.77%    passed     43.87%
  97  FFGR                  68.96%    failed      2.39%
  98  TSBQ                  41.03%    failed      0.18%
  99  VIAU                  71.94%    failed      0.99%
 100  MNNM                  46.55%    failed      0.01%
 101  QWFV                  41.53%    failed      1.35%
 102  GYRP                  57.91%    failed      0.24%
 103  ZZOL                  40.78%    passed     62.56%
 104  Selmo_v1.0            41.96%    passed     71.38%
 105  ABIJ                  59.78%    passed     44.90%
 106  LGDQ                  57.79%    passed     71.37%
 107  GGWH                  41.78%    failed      0.93%
 108  RPRU                  41.90%    passed      8.16%
 109  WCQU                  41.40%    failed      0.07%
 110  WCZB                  57.54%    passed     80.39%
 111  VFIV                  42.02%    failed      0.00%
 112  RAWF                  41.15%    failed      0.01%
 113  ZIVZ                  41.90%    failed      0.00%
 114  BNCU                  41.78%    passed     33.68%
 115  NRWZ                  58.16%    passed     56.60%
 116  JHFI                  42.02%    passed     64.71%
 117  PIUF                  62.38%    passed     97.33%
 118  QZZU                  40.04%    passed      8.84%
 119  YFGP                  62.38%    passed     23.41%
 120  KRUQ                  41.78%    passed     48.56%
 121  UUHD                  66.60%    passed     21.52%
 122  TGKW                  41.15%    passed     94.50%
 123  BHBK                  47.11%    failed      1.95%
 124  FOYQ                  72.81%    passed      6.25%
 125  SILJ                  39.79%    passed     46.70%
 126  LNIL                  74.55%    failed      0.04%
 127  KVAY                  47.73%    passed     65.47%
 128  NNGU                  39.42%    passed     55.59%
 129  XJGM                  72.32%    failed      0.00%
 130  CWYJ                  40.29%    passed     53.50%
 131  TQKZ                  42.02%    passed     58.13%
 132  HERT                  62.14%    passed      7.09%
 133  DCCI                  39.04%    passed     31.26%
 134  HDWF                  54.44%    passed     49.51%
 135  YGCX                  78.77%    passed     67.60%
 136  RKLL                  39.29%    passed     15.22%
 137  RKFX                  52.58%    passed     56.33%
 138  JETM                  39.54%    passed     40.94%
 139  Phavu_v1.0            40.47%    passed     19.72%
 140  KEGA                  40.29%    passed     26.53%
 141  TVSH                  48.73%    passed     40.62%
 142  SUAK                  39.79%    passed     33.94%
 143  SLYR                  43.39%    passed     22.80%
 144  ZSSR                  39.29%    passed     29.42%
 145  VLNB                  39.29%    passed     13.95%
 146  KNMB                  49.35%    passed     80.66%
 147  MYMP                  65.74%    passed     73.87%
 148  JTQQ                  70.08%    passed     31.13%
 149  IWMW                  60.27%    failed      4.43%
 150  RXEN                  39.54%    passed     66.73%
 151  RNBN                  40.29%    passed     41.39%
 152  IZLO                  40.16%    passed     76.97%
 153  IUSR                  78.27%    passed     88.92%
 154  FFFY                  62.76%    passed     13.61%
 155  ILBQ                  66.36%    passed     90.47%
 156  ZCDJ                  39.54%    passed     31.97%
 157  PJSX                  39.54%    passed     27.28%
 158  YZRI                  39.54%    passed     77.93%
 159  PSJT                  78.90%    passed     37.85%
 160  Vitvi_Genoscope.12X   46.00%    passed     82.99%
 161  BGZG                  61.02%    passed     40.43%
 162  OHAE                  74.05%    passed     31.92%
 163  SZPD                  75.79%    passed     46.76%
 164  SMMC                  39.54%    passed     17.38%
 165  AAXJ                  39.54%    passed     29.89%
 166  CBJR                  40.04%    passed     24.04%
 167  ONLQ                  39.54%    passed     32.20%
 168  TXVB                  41.15%    passed     51.67%
 169  ZHMB                  58.91%    passed     87.89%
 170  JTRM                  39.17%    passed     70.66%
 171  CAQZ                  40.66%    passed     70.62%
 172  GSZA                  61.51%    passed     17.61%
 173  XSZI                  48.98%    passed     56.70%
 174  VYGG                  72.44%    passed     31.08%
 175  CCID                  40.29%    passed     46.36%
 176  WEQK                  40.16%    passed     61.84%
 177  TKEK                  74.30%    passed     36.57%
 178  IHPC                  39.66%    passed     61.57%
 179  RTMU                  57.42%    passed     35.39%
 180  XSSD                  39.79%    passed     14.02%
 181  WMLW                  42.02%    passed     17.82%
 182  XPAF                  39.91%    passed     39.93%
 183  CLNU                  57.05%    passed     43.05%
 184  CLMX                  39.29%    passed     58.12%
 185  BERS                  40.16%    passed     39.10%
 186  GJNX                  41.28%    passed     38.39%
 187  HZTS                  40.78%    passed     22.75%
 188  EDIT                  40.29%    passed     33.73%
 189  HQRJ                  39.17%    passed     40.01%
 190  TIUZ                  39.79%    passed     53.12%
 191  ZSGF                  58.41%    passed     49.70%
 192  TPEM                  40.04%    passed     43.39%
 193  AUIP                  62.63%    passed     49.41%
 194  QUTB                  56.55%    passed     32.75%
 195  HAEU                  68.96%    passed      6.88%
 196  WJLO                  42.40%    passed      5.58%
 197  ZJRC                  46.37%    passed     80.03%
 198  EDBB                  39.54%    passed     27.99%
 199  AZBL                  40.91%    passed     34.08%
 200  NUZN                  38.92%    passed     42.31%
 201  ZRIN                  79.14%    passed     27.11%
 202  HLJG                  39.42%    passed     50.50%
 203  LKKX                  39.42%    passed      7.38%
 204  BJKT                  52.20%    passed     43.29%
 205  PKMO                  68.09%    passed     54.86%
 206  WBIB                  39.91%    passed     59.24%
 207  YBQN                  41.78%    passed     70.27%
 208  SFKQ                  61.39%    passed     56.79%
 209  CPKP                  40.53%    passed     30.16%
 210  IRBN                  41.78%    passed     71.21%
 211  UQCB                  49.47%    passed     27.65%
 212  EZGR                  39.42%    passed     19.65%
 213  QEHE                  52.20%    passed     80.21%
 214  BYNZ                  62.76%    passed     30.52%
 215  KDCH                  63.87%    passed     44.73%
 216  HMHL                  43.89%    passed     17.34%
 217  HDSY                  40.04%    passed     36.56%
 218  MVSE                  69.46%    failed      4.89%
 219  CUTE                  41.40%    passed     25.03%
 220  LLQV                  39.79%    passed     43.89%
 221  OHKC                  44.51%    passed     48.51%
 222  PDQH                  39.79%    passed     71.47%
 223  KTWL                  40.91%    passed     81.15%
 224  OLXF                  74.30%    passed     68.39%
 225  COCP                  74.92%    passed     76.02%
 226  KIIX                  66.48%    passed     69.42%
 227  THEW                  74.55%    passed     79.17%
 228  EYRD                  44.63%    passed     33.32%
 229  ZBPY                  39.79%    passed     33.07%
 230  TOKV                  74.55%    passed     19.04%
 231  QURC                  77.16%    passed     40.61%
 232  ZETY                  74.55%    passed      8.35%
 233  JMUI                  72.81%    failed      4.88%
 234  FEDW                  61.39%    passed     69.13%
 235  SALZ                  61.64%    passed     29.63%
 236  VHZV                  60.77%    passed     12.53%
 237  IWIS                  47.24%    passed     38.27%
 238  TNAW                  63.00%    passed     11.25%
 239  UMUL                  78.52%    passed     42.29%
 240  FDMM                  74.30%    passed     46.27%
 241  GGJD                  77.78%    passed     71.24%
 242  UOYN                  61.02%    passed     10.47%
 243  OMYK                  69.09%    passed     53.37%
 244  CPLT                  68.72%    passed     28.04%
 245  LDEL                  68.84%    passed     36.63%
 246  BWRK                  68.59%    passed     59.74%
 247  NMGG                  65.11%    passed     56.74%
 248  GRRW                  69.58%    passed     70.10%
 249  MUNP                  68.84%    passed     29.91%
 250  FYTP                  77.90%    passed     68.94%
 251  ADHK                  65.98%    passed     48.53%
 252  THHD                  68.59%    passed     37.63%
 253  NHAG                  68.96%    passed     21.35%
 254  PTLU                  74.30%    passed     70.62%
 255  QIKZ                  75.67%    failed      4.87%
 256  SBZH                  64.99%    passed     20.95%
 257  NTEO                  74.55%    passed     22.80%
 258  SUVN                  67.97%    passed     60.54%
 259  OINM                  67.97%    passed     60.52%
 260  PCNH                  75.67%    passed     25.78%
 261  UWFU                  70.70%    passed     17.88%
 262  GEHT                  72.19%    passed      6.71%
 263  QZXQ                  74.55%    passed     16.90%
 264  HJMP                  61.64%    passed     44.96%
 265  CMFF                  74.92%    passed     15.61%
 266  KZED                  62.14%    passed      8.46%
 267  JKNQ                  70.95%    passed     97.98%
 268  PWSG                  62.01%    passed     88.90%
 269  QOXT                  75.79%    passed     27.11%
 270  ZQYU                  72.56%    passed     49.66%
 271  CVEG                  75.79%    passed     57.70%
 272  OFTV                  73.81%    passed     73.40%
 273  TFYI                  69.09%    passed      6.98%
 274  LGOW                  73.31%    passed     15.78%
 275  DFDS                  79.64%    passed      5.48%
 276  JPYU                  72.44%    passed     55.36%
 277  KFEB                  40.91%    failed      0.00%
 278  EATP                  41.78%    failed      0.00%
 279  ZNUM                  41.28%    passed     14.03%
 280  GXBM                  42.02%    failed      1.80%
 281  SBLT                  41.53%    failed      0.18%
 282  WEJN                  41.90%    failed      0.34%
 283  UGPM                  41.90%    failed      0.67%
 284  WDCW                  40.91%    failed      0.38%
 285  VZWX                  40.91%    failed      0.01%
 286  BFIK                  42.02%    failed      0.00%
 287  RQFE                  40.91%    failed      1.98%
 288  WSJO                  41.65%    failed      0.00%
 289  JOJQ                  40.91%    passed      7.26%
 290  MFZO                  40.91%    passed      7.41%
 291  NBYP                  40.91%    failed      0.06%
 292  ISIM                  60.02%    failed      0.00%
 293  APTP                  42.27%    failed      1.15%
 294  IEHF                  41.90%    passed      5.23%
 295  DRFX                  71.69%    failed      0.73%
 296  RPQV                  42.89%    failed      0.00%
 297  WDGV                  42.02%    failed      1.63%
 298  QWRA                  69.71%    failed      0.00%
 299  GUBD                  57.91%    passed      9.80%
 300  VHIJ                  76.78%    passed     44.89%
 301  BSNI                  40.78%    passed     27.01%
 302  TWUW                  40.78%    passed     27.01%
 303  FAJB                  58.54%    passed     77.46%
 304  WEEQ                  41.90%    passed     57.62%
 305  DXOU                  53.45%    passed     51.58%
 306  TCBC                  43.64%    passed     23.03%
 307  UCRN                  42.15%    passed     32.74%
 308  AKXB                  41.78%    passed     34.86%
 309  GYBH                  72.19%    passed     70.73%
 310  AEKF                  49.72%    failed      0.00%
 311  GBGT                  73.06%    failed      0.05%
 312  MOYY                  41.78%    failed      0.01%
 313  PKOX                  40.66%    passed     75.93%
 314  YOXI                  41.03%    passed      8.45%
 315  RPGL                  57.91%    failed      1.37%
 316  VAZE                  40.66%    failed      0.08%
 317  JKAA                  40.91%    passed     77.70%
 318  ZYCD                  40.91%    passed     85.59%
 319  ROZZ                  76.41%    passed     64.43%
 320  VIBO                  61.51%    passed     81.46%
 321  RFMZ                  65.49%    passed     58.11%
 322  UOMY                  41.40%    passed     87.18%
 323  PIVW                  41.65%    passed     62.09%
 324  PBUU                  57.17%    passed     95.52%
 325  IKWM                  41.15%    failed      0.29%
 326  NOKI                  41.53%    passed     69.65%
 327  YIXP                  41.53%    passed     64.15%
 328  MEKP                  46.87%    passed     95.83%
 329  NDUV                  41.40%    passed     99.74%
 330  DCDT                  65.61%    passed     66.40%
 331  XDDT                  56.92%    passed     83.49%
 332  WQML                  72.07%    passed     28.29%
 333  LJPN                  41.03%    passed     10.46%
 334  VNAL                  41.03%    passed     10.46%
 335  ORJE                  41.78%    passed     73.26%
 336  UJWU                  43.89%    passed     79.39%
 337  CJNT                  41.78%    passed     76.01%
 338  YLJA                  49.10%    passed     95.37%
 339  GYFU                  42.64%    passed     72.89%
 340  GETL                  53.45%    passed     85.53%
 341  EAAA                  39.91%    passed     29.11%
 342  SNNC                  70.95%    passed     14.51%
 343  CSSK                  40.66%    passed     93.58%
 344  DMLT                  39.04%    passed     86.62%
 345  ATYL                  39.17%    passed     36.40%
 346  BMIF                  58.16%    passed     96.81%
 347  EWXK                  48.36%    passed     52.98%
 348  UCNM                  55.31%    passed     70.24%
 349  UWOD                  41.53%    passed     69.70%
 350  YOWV                  45.13%    passed     41.64%
 351  XXHP                  48.98%    passed     78.18%
 352  PWKQ                  74.05%    passed      8.89%
 353  PAWA                  76.54%    passed     61.55%
 354  HTFH                  59.90%    passed     85.68%
 355  UFJN                  41.40%    passed     73.81%
 356  PNZO                  74.80%    passed     82.62%
 357  NNHQ                  79.27%    passed     33.62%
 358  KJZG                  67.10%    passed     77.20%
 359  GNPX                  39.04%    passed     42.63%
 360  RTNA                  69.58%    passed     24.29%
 361  OPDF                  78.90%    passed     43.83%
 362  ZJUL                  69.83%    passed     56.03%
 363  QIAD                  79.39%    passed     69.20%
 364  VVRN                  72.81%    passed     98.54%
 365  OCZL                  71.69%    passed     45.80%
 366  UJTT                  41.40%    passed     91.47%
 367  POPJ                  41.65%    passed     54.20%
 368  FLTD                  61.76%    passed     26.35%
 369  ZXJO                  74.55%    passed     97.63%
 370  YCKE                  77.41%    passed     57.59%
 371  BMJR                  43.14%    passed     64.67%
 372  WCLG                  74.05%    passed     75.04%
 373  FOMH                  75.67%    failed      4.18%
 374  VYER                  75.29%    failed      0.03%
 375  HFIK                  61.64%    failed      0.01%
 376  JGGD                  65.61%    failed      0.01%
 377  QLMZ                  75.42%    failed      0.05%
 378  ULXR                  69.09%    failed      0.00%
 379  ASZK                  79.64%    failed      0.01%
 380  FSQE                  45.75%    failed      0.00%
 381  ZFGK                  63.50%    passed     70.37%
 382  RCBT                  42.89%    passed     33.56%
 383  GOWD                  45.25%    passed     47.31%
 384  UHLI                  49.10%    passed     46.22%
 385  SKQD                  50.34%    passed     11.61%
 386  WOGB                  41.90%    passed     64.43%
 387  HRWG                  41.78%    passed     53.58%
 388  AWOI                  40.91%    passed     33.30%
 389  Phypa_v3.0            39.79%    passed     31.47%
 390  YEPO                  47.98%    passed     31.28%
 391  KEFD                  55.43%    passed     52.53%
 392  XWHK                  69.34%    passed     94.90%
 393  JMXW                  41.03%    passed     50.07%
 394  CMEQ                  41.03%    passed     55.88%
 395  BGXB                  55.18%    passed     70.21%
 396  WNGH                  41.03%    passed     56.39%
 397  DHWX                  65.86%    passed     82.02%
 398  YWNF                  57.54%    passed     91.12%
 399  ORKS                  41.03%    passed     53.61%
 400  EEMJ                  53.82%    passed     80.59%
 401  JADL                  40.41%    passed     83.24%
 402  QKQO                  41.03%    passed     69.51%
 403  ZACW                  41.03%    passed     77.43%
 404  IGUH                  41.03%    passed     76.08%
 405  VBMM                  50.09%    passed     74.32%
 406  TMAJ                  41.03%    passed     81.86%
 407  QMWB                  54.69%    passed     38.42%
 408  WSPM                  75.42%    failed      1.67%
 409  WVEF                  74.30%    passed     37.49%
 410  TAVP                  56.05%    passed     86.64%
 411  LNSF                  41.03%    passed     76.32%
 412  KYNE                  40.53%    passed     36.53%
 413  VMXJ                  41.03%    passed     60.74%
 414  RGKI                  67.23%    passed     89.21%
 415  NGTD                  57.54%    passed     95.28%
 416  ABCD                  42.15%    passed     36.60%
 417  RDOO                  41.28%    passed     44.78%
 418  GRKU                  55.18%    passed     90.18%
 419  BPSG                  41.40%    passed     65.47%
 420  FFPD                  66.36%    passed     17.52%
 421  RDYY                  79.14%    passed     58.41%
 422  HVBQ                  41.15%    passed     28.30%
 423  SZYG                  51.96%    passed     54.04%
 424  ZTHV                  41.53%    passed     69.42%
 425  FPCO                  41.03%    passed      5.45%
 426  HVNO                  41.15%    failed      0.13%
 427  MULF                  63.13%    failed      0.01%
 428  Spipo_v2              39.54%    passed     55.50%
 429  YMES                  40.29%    passed     53.61%
 430  MFIN                  40.29%    passed     12.10%
 431  Elagu_v2.0            62.51%    passed     56.31%
 432  MTHW                  39.91%    passed     80.10%
 433  XZME                  41.28%    passed     55.04%
 434  LELS                  40.16%    passed     57.49%
 435  VTUS                  39.54%    passed     15.88%
 436  BSTR                  49.84%    passed     55.18%
 437  BYQM                  39.79%    passed     52.60%
 438  SVTS                  66.73%    passed     26.31%
 439  YJUG                  68.47%    passed      6.74%
 440  EMJJ                  39.79%    passed     51.21%
 441  JVBR                  40.04%    passed     24.27%
 442  IXEM                  63.00%    passed     33.94%
 443  MUMD                  39.91%    passed     18.57%
 444  LTZF                  56.67%    passed     72.06%
 445  ROEI                  41.15%    passed     83.34%
 446  Sorbi_v2.1            39.29%    passed     74.79%
 447  SOHV                  39.42%    passed     95.04%
 448  ZMGN                  39.66%    passed     91.81%
 449  BPKH                  39.66%    passed     87.67%
 450  XUAB                  39.29%    passed     92.17%
 451  WCOR                  39.29%    passed     98.11%
 452  XBKS                  39.29%    passed     97.49%
 453  NNOK                  39.29%    passed     97.49%
 454  BXAY                  49.47%    passed     49.68%
 455  VQYB                  40.41%    passed     93.32%
 456  ZENX                  59.78%    passed     85.77%
 457  RCAH                  62.01%    passed     86.01%
 458  YXNR                  39.79%    passed     92.71%
 459  EFCZ                  39.42%    passed     89.11%
 460  Orysa_v7.0            39.29%    passed     60.27%
 461  RMVB                  73.56%    passed     86.54%
 462  IADP                  69.58%    passed     81.41%
 463  WTDE                  39.79%    passed     15.35%
 464  FGRF                  39.17%    passed     18.15%
 465  PLBZ                  42.15%    passed     10.08%
 466  CMCY                  60.40%    passed     36.49%
 467  KXSK                  41.78%    passed     16.04%
 468  ICNN                  39.42%    passed      7.12%
 469  SART                  39.66%    passed     11.37%
 470  ONBE                  56.80%    passed     45.75%
 471  HATH                  39.66%    passed     76.89%
 472  FCEL                  43.76%    passed     46.47%
 473  MVRF                  39.04%    passed     16.85%
 474  XFJG                  39.79%    passed      7.86%
 475  RCUX                  39.29%    passed      7.13%
 476  BLAJ                  74.67%    passed     77.19%
 477  JHUL                  40.04%    passed     43.59%
 478  UZXL                  51.46%    passed     32.08%
 479  HOKG                  69.96%    passed     44.02%
 480  EQDA                  71.57%    passed     36.14%
 481  RQZP                  73.43%    passed     36.91%
 482  PRFO                  72.56%    passed     40.43%
 483  YPIC                  75.17%    passed     90.63%
 484  LSJW                  67.10%    passed     21.03%
 485  MTII                  39.91%    passed     31.74%
 486  Ambtr_v1.0.27         39.17%    passed     23.73%
 487  URDJ                  40.04%    passed     19.15%
 488  SJEV                  48.60%    passed     71.47%
 489  QDVW                  39.29%    passed     32.60%
 490  VYLQ                  39.54%    passed     71.02%
 491  OBPL                  39.29%    passed     39.16%
 492  FALI                  39.42%    passed     54.08%
 493  WPHN                  40.66%    passed     37.74%
 494  WZFE                  39.42%    passed     48.44%
 495  WKSU                  39.17%    passed     30.72%
 496  DHPO                  39.17%    passed     59.38%
 497  WBOD                  39.29%    passed     75.81%
 498  XQWC                  57.67%    passed     74.85%
 499  BSVG                  61.14%    passed     60.03%
 500  KRJP                  39.42%    passed     72.38%
 501  MAQO                  39.17%    passed     40.46%
 502  ABSS                  59.28%    passed     71.64%
 503  WAIL                  70.33%    passed     53.06%
 504  OOSO                  39.54%    passed     30.31%
 505  AFLV                  44.13%    passed     36.52%
 506  WOBD                  40.29%    passed     88.13%
 507  THDM                  57.91%    passed     86.05%
 508  QNPH                  39.79%    passed     19.68%
 509  GDKK                  59.16%    passed     67.18%
 510  FNEN                  48.36%    passed     90.70%
 511  PXYR                  39.79%    passed     47.30%
 512  RHAU                  40.29%    passed     43.22%
 513  OCWZ                  45.00%    passed     34.63%
 514  JHCN                  40.84%    passed     74.47%
 515  VNMY                  40.16%    passed     60.61%
 516  FZJL                  39.66%    passed     35.45%
 517  VVVV                  39.42%    passed     23.27%
 518  ZUHO                  39.79%    passed     73.22%
 519  UPOG                  40.16%    passed     59.32%
 520  GBVZ                  53.45%    passed     61.07%
 521  Aquco_v1.1            39.04%    passed     30.58%
 522  AALA                  39.29%    passed     78.35%
 523  VGHH                  51.71%    passed     66.14%
 524  EVOD                  39.54%    passed     76.79%
 525  DGXS                  49.97%    passed     48.04%
 526  IRAF                  39.42%    passed     10.20%
 527  RQNK                  38.92%    passed     30.09%
 528  QCOU                  39.42%    passed     25.79%
 529  SSDU                  39.42%    passed     44.72%
 530  BEKN                  38.92%    passed     48.86%
 531  XMVD                  39.42%    passed     29.22%
 532  XHKT                  73.81%    passed     43.39%
 533  VZCI                  62.14%    passed     25.54%
 534  MWYQ                  57.17%    passed     70.32%
 535  OBTI                  57.79%    passed     69.15%
 536  SIIK                  38.92%    passed     62.17%
 537  RQUG                  38.92%    passed     63.83%
 538  CIEA                  52.82%    passed      6.23%
 539  WFBF                  40.04%    passed     78.93%
 540  FAKD                  75.42%    passed     54.20%
 541  YGAT                  62.88%    passed     46.33%
 542  FYSJ                  67.47%    passed     30.32%
 543  PPQR                  61.64%    passed     60.51%
 544  YHFG                  58.29%    passed     72.61%
 545  YQIJ                  63.50%    passed     41.21%
 546  BAHE                  65.61%    passed     49.43%
 547  VTLJ                  65.11%    passed     11.87%
 548  PAZJ                  55.68%    passed     63.50%
 549  QTJY                  60.15%    passed     53.50%
 550  HWUP                  40.29%    passed     43.37%
 551  UYED                  70.58%    passed     36.50%
 552  YZVJ                  78.03%    passed     53.66%
 553  UDHA                  40.04%    passed     73.71%
 554  Manes_v4.1            57.91%    passed     72.52%
 555  AUGV                  39.42%    passed     65.90%
 556  FMVB                  79.02%    failed      0.88%
 557  EILE                  68.09%    passed     28.92%
 558  ZGQD                  67.60%    passed     23.38%
 559  OSMU                  39.42%    passed     70.38%
 560  LWCK                  65.11%    passed     21.07%
 561  XNLP                  74.55%    passed     32.79%
 562  VVPY                  74.05%    passed     34.50%
 563  HKMQ                  68.72%    passed     47.01%
 564  HXJE                  71.32%    passed     59.50%
 565  FXGI                  68.96%    passed     35.03%
 566  TPHT                  72.69%    failed      1.85%
 567  TOXE                  39.17%    passed     32.28%
 568  QVMR                  68.59%    passed     50.83%
 569  ALVQ                  57.29%    passed     75.65%
 570  DZTK                  53.82%    passed     71.24%
 571  YADI                  39.91%    passed     41.96%
 572  DSUV                  45.38%    passed     78.74%
 573  JGYZ                  39.29%    passed     36.04%
 574  JCLQ                  39.79%    passed     86.34%
 575  YFQX                  39.79%    passed     86.34%
 576  EDEQ                  39.91%    passed     52.26%
 577  KPUM                  39.91%    passed     70.33%
 578  ECTD                  54.07%    passed     89.72%
 579  Eucgr_v1.1            40.91%    passed     80.75%
 580  FGDU                  39.66%    passed     83.71%
 581  NEBM                  41.03%    passed     96.18%
 582  UPMJ                  41.78%    passed     78.92%
 583  SWPE                  51.33%    passed     58.27%
 584  LAPO                  58.16%    passed     43.35%
 585  TZWR                  79.39%    passed      9.92%
 586  JSAG                  79.02%    passed     72.38%
 587  AJFN                  60.52%    passed     33.35%
 588  VMNH                  50.96%    passed     21.61%
 589  CSUV                  57.91%    passed     38.18%
 590  Arath_TAIR10          40.04%    passed     15.91%
 591  UAXP                  40.53%    passed     27.12%
 592  RTTY                  40.78%    passed     30.43%
 593  QSKP                  42.52%    passed     44.55%
 594  VDKG                  40.04%    passed     50.43%
 595  UPZX                  52.08%    passed     29.58%
 596  LVUS                  40.16%    passed     46.46%
 597  YNUE                  39.79%    passed     74.51%
 598  RJNQ                  72.81%    passed     23.12%
 599  RPPC                  62.01%    passed     29.79%
 600  GKAG                  40.66%    passed     79.04%
 601  YHZW                  40.66%    passed     79.04%
 602  CBAE                  78.77%    passed     65.32%
 603  PUDI                  40.91%    passed     55.95%
 604  AWJM                  64.00%    passed     61.44%
 605  ZTLR                  72.44%    passed     50.42%
 606  LPGY                  39.91%    passed     76.54%
 607  NJLF                  40.04%    passed     82.36%
 608  COAQ                  40.04%    passed     86.97%
 609  EZZT                  41.40%    passed     93.80%
 610  SIZE                  40.29%    passed     89.60%
 611  HBUQ                  77.65%    passed     96.44%
 612  ZBVT                  58.66%    passed     97.72%
 613  KWGC                  55.56%    passed     89.82%
 614  GJPF                  64.37%    passed     16.12%
 615  Poptr_v3.0            39.79%    passed     73.75%
 616  IEPQ                  53.32%    passed     76.95%
 617  INQX                  56.55%    passed     72.52%
 618  TDTF                  40.16%    passed     61.46%
 619  GLVK                  39.79%    passed     61.23%
 620  LFOG                  57.79%    passed     93.48%
 621  KKDQ                  39.79%    passed     64.02%
 622  TXMP                  48.23%    passed     74.23%
 623  VXOD                  39.91%    passed     79.40%
 624  XVRU                  40.53%    passed     82.18%
 625  KCPT                  73.18%    passed     14.92%
 626  JWEY                  64.37%    passed     22.93%
 627  MZOB                  39.66%    passed     97.43%
 628  ABEH                  57.17%    passed     44.00%
 629  DIHD                  39.91%    passed     79.14%
 630  OEKO                  55.43%    passed     70.32%
 631  NIGS                  39.66%    passed     87.98%
 632  IDGE                  39.66%    passed     58.31%
 633  OUER                  56.18%    passed     45.56%
 634  MDJK                  62.63%    passed     77.32%
 635  HNCF                  39.91%    passed     56.65%
 636  AEPI                  39.91%    passed     61.27%
 637  BHYC                  40.16%    passed     72.27%
 638  BVOF                  39.91%    passed     63.89%
 639  MYVH                  39.91%    passed     39.65%
 640  POZS                  40.16%    passed     57.30%
 641  OODC                  63.38%    passed     60.04%
 642  HRUR                  39.66%    passed     59.43%
 643  MKZR                  69.96%    passed     38.16%
 644  AIOU                  40.66%    passed     64.60%
 645  BOLZ                  39.42%    passed     56.17%
 646  Solly_iTAGv2.3        39.04%    passed     72.11%
 647  GHLP                  39.42%    passed     64.92%
 648  LQJY                  39.42%    passed     66.42%
 649  NMDZ                  67.72%    passed     40.20%
 650  WKCY                  40.29%    passed     33.50%
 651  FCCA                  39.91%    passed     65.37%
 652  BCAA                  39.91%    passed     45.21%
 653  YUOM                  39.91%    passed     36.93%
 654  VFFP                  40.29%    passed     34.44%
 655  HBHB                  71.45%    passed     53.05%
 656  WAXR                  39.91%    passed     72.41%
 657  JCMU                  39.04%    passed     10.46%
 658  MXFG                  39.42%    passed      8.47%
 659  EDHN                  40.04%    passed     61.92%
 660  XVJB                  57.79%    passed     55.95%
 661  KYAD                  40.04%    passed     56.21%
 662  BJSW                  40.04%    passed     50.70%
 663  AQGE                  74.05%    passed     21.91%
 664  SZUO                  79.76%    passed     54.02%
 665  RBYC                  70.95%    passed      9.18%
 666  IANR                  39.79%    passed     87.43%
 667  Prupe_v1.0            43.39%    passed     44.97%
 668  VCIN                  66.36%    passed     36.90%
 669  EAVM                  39.91%    passed     67.14%
 670  SQCF                  59.03%    passed     93.78%
 671  XFFT                  56.42%    passed     72.47%
 672  Mimgu_v2.0            39.17%    passed     63.47%
 673  FZQN                  64.74%    passed     42.48%
 674  OQHZ                  76.54%    passed     46.92%
 675  BFJL                  39.42%    passed     44.67%
 676  PQTO                  41.78%    passed     98.10%
 677  XPBC                  76.29%    passed     30.02%
 678  SJAN                  79.52%    passed     42.36%
 679  AXNH                  40.04%    passed     94.89%
 680  PMTB                  40.04%    passed     94.89%
 681  GVCB                  39.91%    passed     99.07%
 682  QFAE                  66.85%    passed     71.15%
 683  HUSX                  68.96%    passed     41.83%
 684  KBRW                  64.25%    passed     61.20%
 685  EQYT                  64.00%    passed     74.23%
 686  ZZEI                  65.61%    passed     54.23%
 687  XNXF                  68.84%    passed     89.48%
 688  JRGZ                  66.60%    failed      0.55%
 689  YHLF                  72.94%    passed     52.68%
 690  BZDF                  40.04%    passed     97.01%
 691  IMZV                  40.04%    passed     97.18%
 692  ARYD                  40.04%    passed     97.18%
 693  ZINQ                  40.04%    passed     97.18%
 694  IDAU                  40.04%    passed     97.18%
 695  ROLB                  49.97%    passed     69.13%
 696  UJGI                  40.04%    passed     95.55%
 697  TLCA                  40.29%    passed     98.24%
 698  GDZS                  39.04%    passed     68.70%
 699  OAGK                  39.66%    passed     74.96%
 700  GUMF                  39.66%    passed     62.89%
 701  JEBK                  41.03%    passed     16.16%
 702  WGTU                  41.40%    passed     84.00%
 703  BWVJ                  41.53%    passed     16.83%
 704  BMSE                  56.67%    passed     72.82%
 705  GCFE                  39.66%    passed     67.96%
 706  WAFT                  42.27%    passed     79.63%
 707  FQGQ                  79.27%    passed     86.97%
 708  MQIV                  51.46%    passed     66.74%
 709  DUNJ                  50.34%    passed     99.02%
 710  WGET                  65.36%    passed     54.39%
 711  ENQF                  65.61%    passed     44.10%
 712  YSQT                  65.74%    passed     12.82%
 713  DUQG                  76.29%    passed     41.65%
 714  Klefl_v1.0            39.79%    failed      0.00%
 715  EMAL                  76.29%    passed     47.98%
 716  FWBF                  70.45%    passed     44.84%
 717  XMQO                  69.71%    passed     45.69%
 718  WQRD                  72.69%    passed     50.11%
 719  FQLP                  75.29%    failed      0.32%
 720  RSOF                  67.85%    passed     28.57%
 721  PSHB                  71.20%    passed     48.53%
 722  MGVU                  66.73%    passed     30.99%
 723  DZLN                  66.98%    passed     67.90%
 724  HXCD                  75.67%    passed     38.98%
 725  NWWI                  67.60%    passed     79.75%
 726  IHJY                  65.74%    passed     27.83%
 727  FTRP                  66.73%    passed     51.26%
 728  CKXF                  74.67%    passed     38.65%
 729  STKJ                  79.76%    failed      1.75%
 730  WQUF                  59.34%    passed     34.02%
 731  EEAQ                  40.41%    passed     39.73%
 732  DJSE                  64.37%    passed     96.50%
 733  WWQZ                  39.91%    passed     68.49%
 734  SWGX                  41.28%    passed     73.08%
 735  DRIL                  61.27%    passed     10.06%
 736  BNDE                  41.40%    passed     64.54%
 737  HKZW                  76.16%    passed     21.20%
 738  CRNC                  40.41%    passed     73.30%
 739  Carpa_ASGPBv0.4       43.76%    passed     77.13%
 740  CZPV                  40.04%    passed     80.81%
 741  FWCQ                  40.91%    passed     83.83%
 742  Theca_v1.1            40.04%    passed     28.25%
 743  EHNF                  43.89%    passed     69.61%
 744  JNUB                  39.54%    passed     15.05%
 745  KJAA                  39.29%    passed     26.26%
 746  OLES                  43.89%    passed     64.79%
 747  DTNC                  39.29%    passed     63.91%
 748  NXTS                  39.29%    passed     35.56%
 749  Musac_v1.0            39.17%    passed     34.11%
 750  YSRZ                  39.79%    passed     69.63%
 751  PUCW                  42.27%    passed     46.79%
 752  FUMQ                  48.36%    passed     70.48%
 753  TORX                  40.41%    passed     60.00%
 754  KTAR                  39.66%    passed     81.27%
 755  MZLD                  39.17%    passed     68.04%
 756  LEMW                  39.79%    passed     25.16%
 757  BDJQ                  47.36%    passed     61.12%
 758  DFYF                  55.80%    passed     51.75%
 759  SXML                  50.22%    passed     88.97%
 760  ASMV                  39.04%    passed     83.36%
 761  AHRN                  39.54%    passed     31.63%
 762  VQFW                  76.04%    passed     97.84%
 763  TZNS                  39.54%    passed     18.84%
 764  YNFJ                  67.35%    passed     37.95%
 765  GIWN                  40.29%    passed     20.85%
 766  SKNL                  47.11%    passed     46.15%
 767  FCBJ                  78.03%    passed     58.71%
 768  VJPU                  40.41%    passed     25.24%
 769  ZBTA                  39.42%    passed     38.51%
 770  JGAB                  39.66%    passed     31.15%
 771  DVXD                  39.66%    passed     39.57%
 772  EGOS                  39.66%    passed     37.93%
 773  DKFZ                  40.04%    passed     89.51%
 774  IKFD                  40.04%    passed     53.39%
 775  CLRW                  79.27%    passed     75.42%
 776  GIPR                  76.54%    passed     62.26%
 777  YKFU                  73.93%    passed     51.76%
 778  OOVX                  76.16%    passed     49.35%
 779  FYUH                  62.38%    passed     63.51%
 780  JAFJ                  39.54%    passed     36.13%
 781  QACK                  72.69%    passed     47.81%
 782  DAYQ                  54.56%    passed     69.28%
 783  UHBY                  39.79%    passed     51.17%
 784  IYDF                  69.21%    passed     56.42%
 785  XHHU                  48.85%    passed     45.65%
 786  SMUR                  59.65%    passed     74.51%
 787  CKKR                  67.97%    passed     42.15%
 788  XTZO                  64.74%    passed     71.40%
 789  RJIM                  70.45%    passed     28.79%
 790  EYKJ                  39.79%    passed     69.00%
 791  LSKK                  40.66%    passed     41.54%
 792  RSCE                  39.91%    passed     99.97%
 793  CPOC                  39.42%    passed     73.92%
 794  DMIN                  39.42%    passed     13.91%
 795  ACWS                  51.96%    passed     88.38%
 796  MIXZ                  73.56%    passed     38.52%
 797  LDME                  39.79%    passed     19.33%
 798  EMBR                  57.91%    passed     48.96%
 799  IZNU                  74.05%    passed     18.19%
 800  VXKB                  39.54%    passed     48.34%
 801  ALUC                  68.96%    passed     14.43%
 802  QSLH                  40.91%    passed     42.73%
 803  OQBM                  39.42%    passed     50.60%
 804  XQRV                  39.42%    passed     52.69%
 805  JDTY                  39.42%    passed     22.96%
 806  NAUM                  79.14%    passed     57.05%
 807  PHCE                  73.31%    passed     21.04%
 808  QICX                  79.02%    passed     57.77%
 809  UVDC                  79.02%    passed     94.95%
 810  IHCQ                  79.14%    passed     59.04%
 811  HANM                  78.03%    passed     68.87%
 812  MRKX                  40.29%    passed     30.74%
 813  WHNV                  39.66%    passed     62.07%
 814  FROP                  77.90%    passed     54.81%
 815  XMBA                  79.64%    passed     38.02%
 816  BKQU                  65.24%    passed     29.31%
 817  FVXD                  78.03%    passed     32.12%
 818  SHEZ                  76.91%    passed      6.01%
 819  CFRN                  68.96%    passed     68.56%
 820  DPFW                  76.91%    passed     59.42%
 821  UOEL                  77.28%    passed     19.63%
 822  MIRS                  76.54%    passed     86.20%
 823  BYPY                  75.54%    passed     30.98%
 824  QKMG                  62.88%    passed     41.31%
 825  XGFU                  39.42%    passed     84.68%
 826  RSPO                  74.43%    passed     52.30%
 827  BSEY                  62.88%    passed     17.18%
 828  JEXA                  52.95%    passed     57.89%
 829  BNTL                  41.65%    passed     58.90%
 830  PZRT                  70.08%    passed     33.98%
 831  ATFX                  53.94%    passed     53.22%
 832  MYZV                  53.45%    passed     47.85%
 833  CVDF                  44.51%    passed     58.86%
 834  RUUB                  61.64%    passed     42.70%
 835  UFHF                  73.68%    passed     80.69%
 836  AVJK                  68.59%    passed     27.07%
 837  LNER                  44.01%    passed     42.75%
 838  PTFA                  39.73%    passed     66.76%
 839  AWQB                  50.47%    passed     82.42%
 840  IIOL                  39.66%    passed     64.10%
 841  DZQM                  53.69%    passed     79.13%
 842  MFTM                  40.04%    passed     65.97%
 843  JBND                  46.49%    passed     74.79%
 844  GGEA                  56.80%    passed     58.25%
 845  IOVS                  39.79%    passed     81.92%
 846  WVWN                  39.79%    passed     86.59%
 847  VSRH                  54.56%    passed     51.84%
 848  GAMH                  40.04%    passed     93.53%
 849  AQFM                  71.69%    passed     43.68%
 850  CKDK                  40.29%    passed     97.08%
 851  WMUK                  44.75%    passed     54.73%
 852  OCTM                  78.27%    passed     33.40%
 853  KPTE                  40.41%    passed     32.01%
 854  ODDO                  38.92%    passed     70.43%
 855  DAAD                  39.04%    passed     50.76%
 856  HTIP                  39.17%    passed     54.62%
 857  SVVG                  39.79%    passed     22.36%
 858  NHUA                  50.34%    passed     65.87%
 859  HENI                  73.43%    passed     21.14%
 860  UZWG                  51.09%    passed     45.16%
 861  NVGZ                  42.02%    passed     80.22%
 862  LWDA                  59.78%    passed     42.55%
 863  IAJW                  49.35%    passed     95.81%
 864  EFMS                  52.58%    passed     89.12%
 865  YFZK                  46.24%    passed     95.47%
 866  TJLC                  51.96%    passed     25.57%
 867  QSNJ                  40.04%    passed     90.95%
 868  ZYAX                  40.91%    passed     88.13%
 869  WWSS                  49.10%    passed     94.03%
 870  CDFR                  53.45%    passed     78.78%
 871  FMWZ                  61.39%    passed     73.73%
 872  IZGN                  55.68%    passed     96.80%
 873  ROWR                  40.29%    passed     98.69%
 874  QHBI                  40.04%    passed     98.47%
 875  OWFC                  49.35%    passed     91.66%
 876  JZVE                  40.04%    passed     76.08%
 877  ZQWM                  40.04%    passed     97.94%
 878  BBDD                  49.84%    passed     66.01%
 879  MHGD                  40.04%    passed     90.40%
 880  VGSX                  40.04%    passed     98.97%
 881  UUJS                  40.04%    passed     98.99%
 882  HILW                  52.82%    passed     87.86%
 883  XLGK                  40.04%    passed     99.61%
 884  SCEB                  40.04%    passed     98.60%
 885  UDUT                  75.54%    passed     73.62%
 886  QCGM                  77.41%    passed     35.34%
 887  INSP                  78.03%    passed     91.08%
 888  JVSZ                  74.05%    passed     54.28%
 889  PUAN                  71.07%    passed     89.11%
 890  JRNA                  50.59%    passed     90.96%
 891  NRXL                  40.04%    passed     85.52%
 892  XIRK                  40.04%    passed     82.06%
 893  ZQVF                  40.04%    passed     92.03%
 894  CGDN                  40.91%    passed     91.97%
 895  UEVI                  40.04%    passed     90.55%
 896  FRPM                  40.04%    passed     81.83%
 897  AIGO                  40.04%    passed     85.25%
 898  BUWV                  40.91%    passed     85.48%
 899  XQSG                  56.92%    passed     94.89%
 900  GMHZ                  41.40%    passed     90.60%
 901  QNGJ                  40.04%    passed     78.40%
 902  VFYZ                  40.04%    passed     88.31%
 903  HQOM                  70.83%    passed     51.74%
 904  YLPM                  49.35%    passed     90.52%
 905  BTTS                  73.56%    passed     42.10%
 906  EGLZ                  75.67%    passed     10.09%
 907  KLGF                  61.89%    passed     89.52%
 908  NKIN                  72.32%    passed     45.89%
 909  XMGP                  66.36%    passed     70.67%
 910  FHST                  67.85%    passed     38.56%
 911  KUXM                  67.60%    passed     95.94%
 912  OXGJ                  72.32%    passed     45.75%
 913  YYPE                  40.04%    passed     59.85%
 914  AUDE                  40.04%    passed     67.94%
 915  ETCJ                  40.04%    passed     55.65%
 916  GKCZ                  40.29%    passed     65.62%
 917  RMMV                  39.54%    passed     96.86%
 918  JDQB                  39.54%    passed     94.53%
 919  IFLI                  74.55%    passed     15.10%
 920  KVFU                  39.17%    passed     67.47%
 921  RWKR                  74.30%    passed     39.21%
 922  JPDJ                  67.47%    passed     38.10%
 923  VKJD                  56.67%    passed     52.14%
 924  QAUE                  39.17%    passed     98.01%
 925  LRTN                  59.78%    passed     58.94%
 926  SERM                  74.18%    passed     75.09%
 927  PPPZ                  65.98%    passed     74.93%
 928  WXVX                  44.38%    passed     92.99%
 929  EJBY                  39.17%    passed     67.17%
 930  VUSY                  39.42%    passed     59.43%
 931  GRFT                  39.29%    passed     61.32%
 932  XRLM                  39.29%    passed     69.92%
 933  NGRR                  62.14%    passed     32.94%
 934  WRPP                  56.67%    passed     60.61%
 935  BEFC                  40.04%    passed     54.23%
 936  XXYA                  39.29%    passed     83.39%
 937  OXYP                  76.29%    passed     76.12%
 938  SIBR                  66.36%    passed     26.32%
 939  GNRI                  41.40%    passed     46.47%
 940  PTBJ                  52.58%    passed     97.39%
 941  YKZB                  41.65%    passed     49.42%
 942  PCGJ                  39.42%    passed     61.39%
 943  KJYC                  78.09%    passed     62.20%
 944  YRHD                  39.29%    passed     59.97%
 945  NBMW                  54.07%    passed     89.29%
 946  AYIY                  62.14%    passed     69.17%
 947  UTQR                  67.72%    passed     62.23%
 948  COBX                  59.03%    passed     67.99%
 949  WOHL                  62.01%    passed     39.03%
 950  EDXZ                  39.04%    passed     56.85%
 951  VYDM                  39.79%    passed     65.09%
 952  FAMO                  64.25%    passed     30.56%
 953  DDRL                  39.54%    passed     58.32%
 954  FUPX                  39.91%    passed     57.38%
 955  KGJF                  39.54%    passed     77.29%
 956  KFZY                  59.03%    passed     90.13%
 957  EJCM                  48.36%    passed     72.35%
 958  JNKW                  39.54%    passed     21.14%
 959  CTYH                  65.24%    passed     32.87%
 960  TAGM                  57.67%    passed     66.84%
 961  UBLN                  62.63%    passed     70.82%
 962  AQZD                  53.13%    passed     75.71%
 963  OWAS                  63.63%    passed     39.04%
 964  TEZA                  56.30%    passed     63.26%
 965  HUQC                  77.41%    passed     57.50%
 966  XRCX                  79.76%    passed     41.72%
 967  BIDT                  59.53%    passed     48.61%
 968  MHYG                  39.66%    passed     66.04%
 969  DESP                  71.32%    passed     74.52%
 970  NVSO                  72.94%    passed     13.24%
 971  JYMN                  44.57%    passed     67.24%
 972  LYPZ                  56.30%    passed     50.86%
 973  YRBQ                  39.91%    passed     66.58%
 974  ZCUA                  39.91%    passed     71.24%
 975  QXWF                  39.66%    passed     68.73%
 976  JEPE                  39.66%    passed     54.38%
WARNING: 428 sequences contain more than 50% gaps/ambiguity
****  TOTAL                 52.27%  136 sequences failed composition chi2 test (p-value<5%; df=3)
NOTE: RFAD is identical to LLEN but kept for subsequent analysis
NOTE: TWUW is identical to BSNI but kept for subsequent analysis
NOTE: VNAL is identical to LJPN but kept for subsequent analysis
NOTE: YFQX is identical to JCLQ but kept for subsequent analysis
NOTE: YHZW is identical to GKAG but kept for subsequent analysis
NOTE: PMTB is identical to AXNH but kept for subsequent analysis
NOTE: ARYD is identical to IMZV but kept for subsequent analysis
NOTE: 2 identical sequences (see below) will be ignored for subsequent analysis
NOTE: ZINQ (identical to IMZV) is ignored but added at the end
NOTE: IDAU (identical to IMZV) is ignored but added at the end
Alignment was printed to /home/laurel/botany-project/data/muscle/alignment/4471-muscle.fa.uniqueseq.phy

For your convenience alignment with unique sequences printed to /home/laurel/botany-project/data/muscle/alignment/4471-muscle.fa.uniqueseq.phy


Create initial parsimony tree by phylogenetic likelihood library (PLL)... 0.657 seconds
Perform fast likelihood tree search using GTR+I+G model...
Estimate model parameters (epsilon = 5.000)
Perform nearest neighbor interchange...
Estimate model parameters (epsilon = 1.000)
1. Initial log-likelihood: -198446.277
2. Current log-likelihood: -198444.309
Optimal log-likelihood: -198444.143
Rate parameters:  A-C: 2.26251  A-G: 3.21812  A-T: 1.22477  C-G: 2.09511  C-T: 3.76028  G-T: 1.00000
Base frequencies:  A: 0.288  C: 0.198  G: 0.264  T: 0.251
Proportion of invariable sites: 0.012
Gamma shape alpha: 0.930
Parameters optimization took 2 rounds (6.900 sec)
Time for fast ML tree search: 34.475 seconds

NOTE: ModelFinder requires 469 MB RAM!
ModelFinder will test up to 286 DNA models (sample size: 1611) ...
 No. Model         -LnL         df  AIC          AICc         BIC
  1  GTR+F         215803.133   1953 435512.265   8067836.265  446028.409
  2  GTR+F+I       215076.005   1954 434060.010   8074200.010  444581.539
  3  GTR+F+G4      198534.265   1954 400976.530   8041116.530  411498.059
  4  GTR+F+I+G4    198444.088   1955 400798.175   8048758.175  411325.089
  5  GTR+F+R2      201918.508   1955 407747.015   8055707.015  418273.929
  6  GTR+F+R3      199410.565   1957 402735.130   8066347.130  413272.812
  7  GTR+F+R4      198395.799   1959 400709.598   8079989.598  411258.050
  8  GTR+F+R5      198037.244   1961 399996.488   8094960.488  410555.709
  9  GTR+F+R6      197765.754   1963 399457.509   8110121.509  410027.499
 10  GTR+F+R7      197625.850   1965 399181.701   8125561.701  409762.460
 11  GTR+F+R8      197571.787   1967 399077.573   8141189.573  409669.102
 12  GTR+F+R9      197526.607   1969 398991.213   8156851.213  409593.511
 13  GTR+F+R10     197514.255   1971 398970.510   8172594.510  409583.577
 25  SYM+R9        197310.796   1966 398553.593   8132797.593  409139.737
 26  SYM+R10       197302.940   1968 398541.879   8148525.879  409138.793
 38  TVM+F+R9      197548.969   1968 399033.938   8149017.938  409630.851
 39  TVM+F+R10     197537.713   1970 399015.425   8164755.425  409623.108
 51  TVMe+R9       197350.002   1965 398630.005   8125010.005  409210.764
 52  TVMe+R10      197333.549   1967 398601.098   8140713.098  409192.627
 64  TIM3+F+R9     197547.263   1967 399028.526   8141140.526  409620.054
 65  TIM3+F+R10    197534.746   1969 399007.492   8156867.492  409609.790
 77  TIM3e+R9      197362.669   1964 398653.338   8117173.338  409228.713
 78  TIM3e+R10     197356.244   1966 398644.488   8132888.488  409230.632
 90  TIM2+F+R9     198191.364   1967 400316.728   8142428.728  410908.257
 91  TIM2+F+R10    198174.529   1969 400287.059   8158147.059  410889.357
103  TIM2e+R9      197658.691   1964 399245.383   8117765.383  409820.757
104  TIM2e+R10     197648.642   1966 399229.284   8133473.284  409815.428
116  TIM+F+R9      198204.867   1967 400343.735   8142455.735  410935.263
117  TIM+F+R10     198197.057   1969 400332.113   8158192.113  410934.411
129  TIMe+R9       197701.537   1964 399331.074   8117851.074  409906.449
130  TIMe+R10      197688.131   1966 399308.263   8133552.263  409894.407
142  TPM3u+F+R9    197571.489   1966 399074.979   8133318.979  409661.123
143  TPM3u+F+R10   197560.163   1968 399056.326   8149040.326  409653.240
155  TPM3+F+R9     197571.489   1966 399074.979   8133318.979  409661.123
156  TPM3+F+R10    197560.163   1968 399056.326   8149040.326  409653.240
168  TPM2u+F+R9    198214.939   1966 400361.877   8134605.877  410948.021
169  TPM2u+F+R10   198210.008   1968 400356.015   8150340.015  410952.929
181  TPM2+F+R9     198214.939   1966 400361.877   8134605.877  410948.021
182  TPM2+F+R10    198210.008   1968 400356.015   8150340.015  410952.929
194  K3Pu+F+R9     198228.629   1966 400389.259   8134633.259  410975.403
195  K3Pu+F+R10    198219.984   1968 400375.969   8150359.969  410972.882
207  K3P+R9        197728.992   1963 399383.984   8110047.984  409953.974
208  K3P+R10       197722.421   1965 399374.843   8125754.843  409955.602
220  TN+F+R9       198205.171   1966 400342.342   8134586.342  410928.486
221  TN+F+R10      198197.400   1968 400330.801   8150314.801  410927.714
233  TNe+R9        197701.802   1963 399329.605   8109993.605  409899.595
234  TNe+R10       197688.474   1965 399306.948   8125686.948  409887.707
246  HKY+F+R9      198228.950   1965 400387.901   8126767.901  410968.660
247  HKY+F+R10     198220.030   1967 400374.061   8142486.061  410965.589
259  K2P+R9        197729.245   1962 399382.490   8102194.490  409947.095
260  K2P+R10       197722.659   1964 399373.318   8117893.318  409948.693
272  F81+F+R9      200601.509   1964 405131.018   8123651.018  415706.392
273  F81+F+R10     200598.780   1966 405129.561   8139373.561  415715.705
285  JC+R9         200381.237   1961 404684.474   8099648.474  415243.695
286  JC+R10        200366.651   1963 404659.302   8115323.302  415229.292
Akaike Information Criterion:           SYM+R10
Corrected Akaike Information Criterion: GTR+F+G4
Bayesian Information Criterion:         SYM+R10
Best-fit model: SYM+R10 chosen according to BIC

All model information printed to /home/laurel/botany-project/data/muscle/alignment/4471-muscle.fa.model.gz
CPU time for ModelFinder: 12524.236 seconds (3h:28m:44s)
Wall-clock time for ModelFinder: 12561.644 seconds (3h:29m:21s)

NOTE: 469 MB RAM (0 GB) is required!
Estimate model parameters (epsilon = 0.100)
1. Initial log-likelihood: -197302.940
2. Current log-likelihood: -197302.835
Optimal log-likelihood: -197302.700
Rate parameters:  A-C: 1.92376  A-G: 3.52733  A-T: 1.27149  C-G: 1.67358  C-T: 3.00529  G-T: 1.00000
Base frequencies:  A: 0.250  C: 0.250  G: 0.250  T: 0.250
Site proportion and rates:  (0.099,0.038) (0.132,0.123) (0.095,0.244) (0.105,0.423) (0.107,0.670) (0.135,0.996) (0.098,1.330) (0.105,1.751) (0.090,2.525) (0.033,4.972)
Parameters optimization took 2 rounds (51.658 sec)
Computing ML distances based on estimated model parameters... 24.748 sec
Computing BIONJ tree...
1.704 seconds
Log-likelihood of BIONJ tree: -202045.354
--------------------------------------------------------------------
|             INITIALIZING CANDIDATE TREE SET                      |
--------------------------------------------------------------------
Generating 98 parsimony trees... 92.595 second
Computing log-likelihood of 98 initial trees ... 138.930 seconds
Current best score: -197302.700

Do NNI search on 20 best initial trees
Estimate model parameters (epsilon = 0.100)
BETTER TREE FOUND at iteration 1: -197277.663
Estimate model parameters (epsilon = 0.100)
BETTER TREE FOUND at iteration 2: -196632.143
Iteration 10 / LogL: -196755.290 / Time: 0h:16m:8s
Estimate model parameters (epsilon = 0.100)
BETTER TREE FOUND at iteration 11: -196624.323
Iteration 20 / LogL: -196742.899 / Time: 0h:26m:16s
Finish initializing candidate tree set (20)
Current best tree score: -196624.323 / CPU time: 1448.447
Number of iterations: 20
--------------------------------------------------------------------
|               OPTIMIZING CANDIDATE TREE SET                      |
--------------------------------------------------------------------
Estimate model parameters (epsilon = 0.100)
BETTER TREE FOUND at iteration 27: -196612.620
Iteration 30 / LogL: -196683.603 / Time: 0h:32m:20s (1h:48m:11s left)
Estimate model parameters (epsilon = 0.100)
BETTER TREE FOUND at iteration 37: -196582.234
Iteration 40 / LogL: -196607.348 / Time: 0h:39m:1s (1h:37m:3s left)
Estimate model parameters (epsilon = 0.100)
BETTER TREE FOUND at iteration 42: -196575.455
Iteration 50 / LogL: -196587.606 / Time: 0h:45m:5s (1h:24m:39s left)
Estimate model parameters (epsilon = 0.100)
BETTER TREE FOUND at iteration 51: -196562.542
Estimate model parameters (epsilon = 0.100)
BETTER TREE FOUND at iteration 58: -196547.118
Iteration 60 / LogL: -196637.077 / Time: 0h:51m:48s (1h:26m:2s left)
Iteration 70 / LogL: -196662.681 / Time: 0h:57m:0s (1h:12m:42s left)
Estimate model parameters (epsilon = 0.100)
BETTER TREE FOUND at iteration 77: -196523.456
Iteration 80 / LogL: -196566.244 / Time: 1h:2m:45s (1h:17m:3s left)
Iteration 90 / LogL: -196563.107 / Time: 1h:7m:40s (1h:6m:9s left)
Iteration 100 / LogL: -196593.439 / Time: 1h:12m:47s (0h:56m:37s left)
Iteration 110 / LogL: -196583.579 / Time: 1h:17m:37s (0h:47m:43s left)
Iteration 120 / LogL: -196523.476 / Time: 1h:22m:45s (0h:39m:38s left)
Iteration 130 / LogL: -196574.554 / Time: 1h:27m:42s (0h:31m:57s left)
Estimate model parameters (epsilon = 0.100)
BETTER TREE FOUND at iteration 137: -196509.040
Iteration 140 / LogL: -196560.053 / Time: 1h:33m:42s (1h:5m:23s left)
Iteration 150 / LogL: -196553.587 / Time: 1h:38m:35s (0h:57m:34s left)
Iteration 160 / LogL: -196521.369 / Time: 1h:43m:43s (0h:50m:13s left)
Iteration 170 / LogL: -196528.211 / Time: 1h:48m:51s (0h:43m:9s left)
Iteration 180 / LogL: -196590.634 / Time: 1h:54m:18s (0h:36m:23s left)
Iteration 190 / LogL: -196571.712 / Time: 1h:59m:35s (0h:29m:44s left)
Iteration 200 / LogL: -196513.225 / Time: 2h:4m:37s (0h:23m:10s left)
Estimate model parameters (epsilon = 0.100)
BETTER TREE FOUND at iteration 210: -196505.151
Iteration 210 / LogL: -196505.151 / Time: 2h:10m:23s (1h:2m:23s left)
Estimate model parameters (epsilon = 0.100)
BETTER TREE FOUND at iteration 215: -196504.504
Iteration 220 / LogL: -196542.974 / Time: 2h:15m:59s (0h:58m:59s left)
Estimate model parameters (epsilon = 0.100)
BETTER TREE FOUND at iteration 226: -196504.085
Iteration 230 / LogL: -196545.894 / Time: 2h:21m:51s (0h:59m:28s left)
Iteration 240 / LogL: -196523.971 / Time: 2h:27m:7s (0h:52m:56s left)
Iteration 250 / LogL: -196525.089 / Time: 2h:32m:13s (0h:46m:27s left)
Iteration 260 / LogL: -196554.335 / Time: 2h:37m:21s (0h:40m:5s left)
Iteration 270 / LogL: -196562.460 / Time: 2h:42m:28s (0h:33m:49s left)
Iteration 280 / LogL: -196615.401 / Time: 2h:47m:44s (0h:27m:39s left)
Estimate model parameters (epsilon = 0.100)
BETTER TREE FOUND at iteration 288: -196498.134
Iteration 290 / LogL: -196508.201 / Time: 2h:53m:45s (0h:58m:55s left)
Iteration 300 / LogL: -196560.349 / Time: 2h:59m:21s (0h:52m:47s left)
Iteration 310 / LogL: -196529.577 / Time: 3h:5m:27s (0h:46m:48s left)
Iteration 320 / LogL: -196568.388 / Time: 3h:11m:11s (0h:40m:45s left)
Iteration 330 / LogL: -196530.156 / Time: 3h:19m:58s (0h:35m:15s left)
Iteration 340 / LogL: -196538.808 / Time: 3h:35m:32s (0h:30m:31s left)
Iteration 350 / LogL: -196513.769 / Time: 3h:49m:15s (0h:24m:57s left)
Iteration 360 / LogL: -196584.942 / Time: 4h:2m:6s (0h:18m:52s left)
Estimate model parameters (epsilon = 0.100)
BETTER TREE FOUND at iteration 366: -196488.467
Iteration 370 / LogL: -196608.513 / Time: 4h:15m:28s (1h:6m:28s left)
Iteration 380 / LogL: -196531.332 / Time: 4h:26m:58s (1h:0m:34s left)
Estimate model parameters (epsilon = 0.100)
BETTER TREE FOUND at iteration 388: -196478.204
Iteration 390 / LogL: -196535.568 / Time: 4h:34m:20s (1h:9m:6s left)
Iteration 400 / LogL: -196519.667 / Time: 4h:39m:46s (1h:1m:42s left)
Estimate model parameters (epsilon = 0.100)
BETTER TREE FOUND at iteration 403: -196476.428
Iteration 410 / LogL: -196543.289 / Time: 4h:46m:16s (1h:5m:5s left)
Iteration 420 / LogL: -196587.010 / Time: 4h:51m:34s (0h:57m:45s left)
Estimate model parameters (epsilon = 0.100)
BETTER TREE FOUND at iteration 425: -196475.478
Iteration 430 / LogL: -196530.481 / Time: 4h:57m:43s (1h:5m:55s left)
Iteration 440 / LogL: -196507.340 / Time: 5h:3m:26s (0h:58m:45s left)
Iteration 450 / LogL: -196522.269 / Time: 5h:9m:7s (0h:51m:38s left)
Iteration 460 / LogL: -196515.342 / Time: 5h:14m:50s (0h:44m:35s left)
Iteration 470 / LogL: -196549.699 / Time: 5h:20m:13s (0h:37m:33s left)
Iteration 480 / LogL: -196523.925 / Time: 5h:28m:53s (0h:30m:53s left)
Iteration 490 / LogL: -196702.169 / Time: 5h:36m:19s (0h:24m:4s left)
Iteration 500 / LogL: -196540.819 / Time: 5h:42m:58s (0h:17m:10s left)
Iteration 510 / LogL: -196494.492 / Time: 5h:50m:10s (0h:10m:19s left)
Iteration 520 / LogL: -196549.527 / Time: 5h:57m:20s (0h:3m:26s left)
TREE SEARCH COMPLETED AFTER 526 ITERATIONS / Time: 6h:1m:26s

--------------------------------------------------------------------
|                    FINALIZING TREE SEARCH                        |
--------------------------------------------------------------------
Performs final model parameters optimization
Estimate model parameters (epsilon = 0.010)
1. Initial log-likelihood: -196475.478
2. Current log-likelihood: -196475.420
3. Current log-likelihood: -196475.355
4. Current log-likelihood: -196475.314
5. Current log-likelihood: -196475.272
6. Current log-likelihood: -196475.225
7. Current log-likelihood: -196475.188
8. Current log-likelihood: -196475.163
9. Current log-likelihood: -196475.128
10. Current log-likelihood: -196475.099
11. Current log-likelihood: -196475.064
12. Current log-likelihood: -196475.042
13. Current log-likelihood: -196474.994
14. Current log-likelihood: -196474.957
15. Current log-likelihood: -196474.924
16. Current log-likelihood: -196474.902
17. Current log-likelihood: -196474.871
18. Current log-likelihood: -196474.855
19. Current log-likelihood: -196474.840
20. Current log-likelihood: -196474.819
21. Current log-likelihood: -196474.800
Optimal log-likelihood: -196474.784
Rate parameters:  A-C: 1.91518  A-G: 3.50216  A-T: 1.26132  C-G: 1.64351  C-T: 3.00585  G-T: 1.00000
Base frequencies:  A: 0.250  C: 0.250  G: 0.250  T: 0.250
Site proportion and rates:  (0.091,0.033) (0.125,0.114) (0.100,0.222) (0.104,0.395) (0.111,0.642) (0.128,0.965) (0.123,1.324) (0.103,1.800) (0.086,2.607) (0.030,5.149)
Parameters optimization took 21 rounds (409.057 sec)
BEST SCORE FOUND : -196474.784
Total tree length: 104.600

Total number of iterations: 526
CPU time used for tree search: 21433.304 sec (5h:57m:13s)
Wall-clock time used for tree search: 21558.525 sec (5h:59m:18s)
Total CPU time used: 21968.976 sec (6h:6m:8s)
Total wall-clock time used: 22095.682 sec (6h:8m:15s)

Analysis results written to: 
  IQ-TREE report:                /home/laurel/botany-project/data/muscle/alignment/4471-muscle.fa.iqtree
  Maximum-likelihood tree:       /home/laurel/botany-project/data/muscle/alignment/4471-muscle.fa.treefile
  Likelihood distances:          /home/laurel/botany-project/data/muscle/alignment/4471-muscle.fa.mldist
  Screen log file:               /home/laurel/botany-project/data/muscle/alignment/4471-muscle.fa.log

Date and Time: Sun Mar 31 16:57:18 2024
```
The same process is then repeated on the KALIGN-aligned gene 4471 sequence and the other 7 genes used in the substudy.



##Coestimation with ASTRAL
I chose to perform coestimation using ASTRAL of a subset of genes from the genes examined in the original study. ASTRAL is a software that is used to inference species trees from a set of individual gene trees while accounting for gene tree discordance. ASTRAL assumes there will be some gene tree discordance from the true species tree and newer versions may remove branches with lower branch support from consideration as it assumes this improves accuracy. Some strengths of ASTRAL include statistical consistency so that as the number of genes increases it continues to converge to the true species tree, scalability to handle large data sets with many species and genes, and high relative accuracy. ASTRAL is limited by the quality of gene tree input, reliance on modern organisms so that incomplete data sets can alter results, and lower accuracy in areas of low branch support. Users can choose to exclude low-support branches using bootstrap.

The input data for ASTRAL is a .txt file containing the Newick formatted data from each gene tree. This information can be found in the IQ-TREE output file with the ending .treefile. Create a new .txt file in a separate folder, and manually copy and paste the Newick data for all 8 MUSCLE-aligned genes into the .txt file. Create a new file and do the same for the KALIGN-aligned data. The output from 4471-muscle.fa.treefil looks like this:
```
(VRGZ:0.4064704192,(YRMA:0.0151024863,RWXW:0.0734826841):0.1630104964,(((((((((WCLV:0.7691868146,ISIM:0.2017596127):0.0657708962,(((((((((((((NQYP:0.5484215663,LSHT:0.4064285397):0.0729534353,KSFK:0.2581717112):0.0265075424,(JIWJ:0.2989506795,GJIY:0.1877847617):0.0960434741):0.0497751041,(OQON:0.3102138392,(CQQP:0.3234176308,OAEZ:0.1633281414):0.1034023197):0.1378217103):0.0344908294,AJAU:0.3592753208):0.1220177937,Cyame_v1.0:0.8924113290):0.1132058176,(VFIV:0.1491083969,(RAWF:0.1524362840,JMUI:0.1917234101):0.0200802175):0.2032306249):0.0151944562,BAZF:0.4587312719):0.0320787617,(DZPJ:0.4560229327,(ALZF:0.3933979235,NSTT:0.2376738700):0.0851021282):0.0251981833):0.0298762185,(((((FXHG:0.3451191620,BZSH:0.2922397906):0.0528881627,((XKWQ:0.2417903485,VJDZ:0.1593310855):0.1099835950,(MWAN:0.2960957221,PZIF:0.1357699603):0.0678589698):0.0494493316):0.0356000620,(((((((PRIQ:0.0504831472,QRTH:0.0404807705):0.1562335055,ZLQE:0.3748091669):0.0459539372,(KFEB:0.2439672304,GUBD:0.2864722353):0.1039347118):0.0498883215,(ZFXU:0.3069148249,(SYJM:0.2141442568,(UKUC:0.0310053365,WDWX:0.1464703786):0.3924869767):0.1281812507):0.1444780086):0.0428307013,(((TSBQ:0.2160685755,JRGZ:0.3151639969):0.0526818881,(ZIVZ:0.1737815287,LNIL:0.2393794608):0.0676209525):0.0148499064,MULF:0.3641836840):0.0297370610):0.0267723755,(JKKI:0.1807965513,(ISPU:0.1188610880,(Chlre_v5.5:0.1017108009,(Volca_v2.0:0.0944391713,(WRSL:0.0155402807,JWGT:0.0055636503):0.0859249192):0.0712423409):0.0326425106):0.1094719057):0.0952311218):0.0259444046,(AOUJ:0.1485206116,(AJUW:0.0456315847,LBRP:0.0443814954):0.0733491133):0.1372402477):0.0232840237):0.0382534347,FOYQ:0.3059230309):0.0279554315,(XDLL:0.3862903980,VIAU:0.5018696173):0.1405509921):0.0328290866):0.0379483008,(((((BILC:0.3922262064,(EEJO:0.1216108403,Chlva_v1.0:0.0890926605):0.0910358589):0.0664688557,NKXU:0.2829636903):0.0272737007,((MFYC:0.4074615598,MNCB:0.2303816926):0.0475443939,AKCR:0.2088271486):0.0184265166):0.0556184669,(AYPS:0.4096055993,((DUMA:0.3120321396,HVNO:0.2385853598):0.0264372104,FMVB:0.2394047288):0.1254309415):0.0564403814):0.0190771115,((EATP:0.3044526824,GXBM:0.1654439227):0.0515970853,ZNUM:0.4230598382):0.0526488761):0.0137974324):0.0175227421,(XIVI:0.2323449829,TNAW:0.4349008364):0.0865150228):0.0238189374,(((((((((XAXW:1.2433181099,SNOX:0.4159090210):0.1267430826,((((((((ZZOL:0.0023295877,Selmo_v1.0:0.0154189854):0.0263123896,(ABIJ:0.0708606148,(LGDQ:0.0546175206,KJYC:0.3032589244):0.0385866909):0.0699184991):0.1223585090,((JKAA:0.0354187776,ZYCD:0.0246393290):0.0956814478,ZFGK:0.1411645473):0.0431358826):0.0474073905,KUXM:0.2048448074):0.0131663996,((UPMJ:0.0997063085,(((GKAG:0.0000010944,YHZW:0.0000010944):0.0136530179,CBAE:0.0414397566):0.0089990389,ZZEI:0.0985519304):0.0407273416):0.0135590841,((PQTO:0.0410890749,(XNXF:0.0661794873,ENQF:0.0066817725):0.0023294696):0.0040848416,WAFT:0.0802315437):0.0188561088):0.1429934306):0.0233819615,((((BNCU:0.1500321956,((NRWZ:0.1337661000,YFGP:0.1544350827):0.0272431362,((JHFI:0.0350405836,PIUF:0.0092251803):0.0943397930,(((KRUQ:0.0311684757,UUHD:0.0311548261):0.0816546025,((RTMU:0.1075178629,(OFTV:0.0176051981,LGOW:0.2697106398):0.0193915282):0.0168529666,(YBQN:0.0300226568,IRBN:0.0332021442):0.0376289103):0.0503640989):0.0206333533,TGKW:0.1380934119):0.0126120775):0.0113290901):0.0097138587):0.0495855176,((((HERT:0.0784017085,(HMHL:0.0501160311,TFYI:0.0518328275):0.0000023143):0.0115707948,TXVB:0.0770806070):0.0087598258,WJLO:0.0667738915):0.0044355490,(ILBQ:0.0864724829,JPYU:0.1102826085):0.0485371111):0.0453106621):0.0690305018,((BSNI:0.0000010944,TWUW:0.0000010944):0.0564258593,((FAJB:0.0332398996,WEEQ:0.0257634479):0.0405883412,(((DXOU:0.0266619961,TCBC:0.0060423908):0.0452360672,AKXB:0.0046901519):0.0107042647,UCRN:0.0342609748):0.0693204817):0.0664597025):0.1233552093):0.0253991048,(((RCBT:0.0135638044,(GOWD:0.0037112754,UHLI:0.0096785146):0.0030813151):0.0997872313,(WOGB:0.0883998131,(((HRWG:0.0867086990,HVBQ:0.0610927650):0.0111854528,(AWOI:0.0514465515,((Phypa_v3.0:0.0025340068,YEPO:0.0080756266):0.0682960879,(KEFD:0.0344444750,(((XWHK:0.0000024676,JMXW:0.0138634886):0.0388478367,(((((CMEQ:0.0352930152,(((((((DHWX:0.0283825722,(ZACW:0.0000030032,IGUH:0.0013925680):0.0082545283):0.0015986877,(VBMM:0.0353102406,(TMAJ:0.0041847899,(QMWB:0.0234136306,WSPM:0.0250574278):0.0000021865):0.0000025047):0.0014279099):0.0027614899,QKQO:0.0084472254):0.0027974717,LNSF:0.0028063456):0.0014037757,(EEMJ:0.0133687681,MIRS:0.1408492844):0.0000021308):0.0014064106,TAVP:0.0058947526):0.0028113044,JADL:0.0059087559):0.0154670933):0.0020509680,WNGH:0.0091993761):0.0030669222,YWNF:0.0289859997):0.0041195818,ORKS:0.0122835664):0.0020145659,BGXB:0.0136167350):0.0088728168):0.0164113202,(((((VMXJ:0.0000029375,RGKI:0.0000010944):0.0313903385,NGTD:0.0316472969):0.0043161800,(ABCD:0.0029639905,RDOO:0.0031436445):0.0185991770):0.0017948726,GRKU:0.0694361520):0.0070065572,(BPSG:0.0220293509,FFPD:0.0476544187):0.0165959489):0.0145050179):0.0311292060):0.0076397693):0.0527993843):0.0318561546):0.0112754463,(SZYG:0.0082079554,ZTHV:0.0416747977):0.0568356056):0.0439479445):0.0320609886):0.0470782266,SKQD:0.0686102371):0.0479792719):0.0301077692):0.0181667244,(((((((((((((QZZU:0.1818392105,((((RKLL:0.0513923274,(RKFX:0.0146546267,JETM:0.0377549983):0.0162173862):0.0143726346,(((((Phavu_v1.0:0.0214404981,((KEGA:0.0071251864,TVSH:0.0368038991):0.0163779325,SUAK:0.0266144083):0.0079653875):0.0229189047,((KNMB:0.0612817345,(MYMP:0.0220196165,HJMP:0.0645770210):0.0023826078):0.0266625408,JTQQ:0.0194316919):0.0128393868):0.0040728570,VLNB:0.0349601913):0.0114766695,(SLYR:0.0205655524,ZSSR:0.0208816903):0.0111892382):0.0194943512,(((ZCDJ:0.0025953936,PJSX:0.0103340833):0.0506936070,KZED:0.0733614355):0.0161401337,(VHZV:0.0148488603,(GEHT:0.0000030636,(QZXQ:0.0000030663,CMFF:0.1147516604):0.0460998645):0.0637007156):0.0425123507):0.0122320223):0.0077452190):0.0450549624,OQHZ:0.0524985946):0.0122169070,OHAE:0.1081339800):0.0167118665):0.0025957932,PTLU:0.0756660166):0.0039912781,(KVAY:0.0764579970,ZHMB:0.1035808590):0.0406542000):0.0167148987,(((WVEF:0.0702892763,EILE:0.0101310088):0.0995782941,((WKCY:0.0857340089,(EDHN:0.0209009822,XVJB:0.0161318285):0.0114004286):0.0075220177,(KYAD:0.0185358257,(BJSW:0.0182499575,AQGE:0.0247743007):0.0308944980):0.0353281372):0.0053543707):0.0196594936,(RBYC:0.1309770724,(IANR:0.0444365728,((Prupe_v1.0:0.0323442542,(VCIN:0.0437127530,EAVM:0.0058597412):0.0230493055):0.0102924128,(SQCF:0.0133637421,XFFT:0.0019816344):0.0112140708):0.0076787286):0.0308608432):0.0107182169):0.0206702324):0.0057824032,((((((HDWF:0.0884179447,OBTI:0.0745936236):0.0064213999,((((FCCA:0.0352346342,BCAA:0.0327342034):0.0043488700,YUOM:0.0309029322):0.0043169083,((VFFP:0.0309905425,HBHB:0.0266834099):0.0052570776,WAXR:0.0287837189):0.0233640505):0.0026634021,IKFD:0.0634212150):0.0180923036):0.0095013311,(((TIUZ:0.0297141325,THHD:0.0792124870):0.0055120371,YZVJ:0.0504027687):0.0077117790,JHCN:0.0772740621):0.0222886624):0.0047464323,((((PXYR:0.0261756053,RHAU:0.0328981787):0.0356617095,(PAZJ:0.0373497762,((Manes_v4.1:0.0162609762,XNLP:0.0149094018):0.0206237250,VVPY:0.0307552373):0.0073298041):0.0093742116):0.0178197732,(YGAT:0.0842478649,(Poptr_v3.0:0.0065990882,(IEPQ:0.0023513659,((INQX:0.0019925905,LFOG:0.0000010944):0.0000025665,((TDTF:0.0013332311,GLVK:0.0013313917):0.0000026244,KKDQ:0.0000010944):0.0000025000):0.0011807781):0.0164853657):0.0383849441):0.0017671470):0.0033800396,((VNMY:0.0807356664,(((COAQ:0.0562067144,(EZZT:0.0095762333,SIZE:0.0050625056):0.0492239647):0.0224907253,((((TXMP:0.0067319303,VXOD:0.0071386484):0.0073013645,KCPT:0.0306004902):0.0164904373,(HNCF:0.0287779652,(((AEPI:0.0013468740,(BHYC:0.0000022343,BVOF:0.0013310588):0.0026467616):0.0155445943,MYVH:0.0143230587):0.0014500364,(POZS:0.0000024113,OODC:0.0000024140):0.0130996181):0.0074912742):0.0157299644):0.0809680014,CKDK:0.0954931354):0.0060998722):0.0049917464,(BNDE:0.0761420846,FWCQ:0.0788817666):0.0406536375):0.0039874248):0.0045652918,((RPPC:0.1247845797,(ZTLR:0.0450741992,(LPGY:0.0176264068,NJLF:0.0048525792):0.0289875197):0.0081764268):0.0185847252,((HBUQ:0.0118681651,ZBVT:0.0072800245):0.0555162721,(XPBC:0.0892460425,IHCQ:0.1121351880):0.0253694281):0.0097338959):0.0157420809):0.0023849389):0.0252610486):0.0114416981,((((PKMO:0.1227204504,WMUK:0.0220362448):0.0290265898,KPTE:0.0475004105):0.0053302060,ATFX:0.0781642502):0.0085411609,((((ZJUL:0.0932162509,DRIL:0.1029500385):0.0162634734,KWGC:0.0389825735):0.0928438158,(PUDI:0.0463257823,AWJM:0.0103329043):0.0850306091):0.0124693542,Theca_v1.1:0.0525849521):0.0088484799):0.0134308297):0.0030908728,((YGCX:0.0764076022,(((DZTK:0.1472039759,(SWPE:0.0745737544,(((((LAPO:0.0659497067,((VMNH:0.0189574275,CSUV:0.0163807830):0.0050114338,Arath_TAIR10:0.0149279214):0.0124259574):0.0068589373,TZWR:0.0503784621):0.0567012149,((QSKP:0.0248998142,LVUS:0.0161131142):0.0061049691,(VDKG:0.0268338622,UPZX:0.0173889402):0.0056904253):0.0152532844):0.0049632504,RTTY:0.0428950792):0.0245645253,UAXP:0.0407007240):0.0103654215):0.0176427746):0.0131290916,CRNC:0.0749158257):0.0152555397,(Carpa_ASGPBv0.4:0.0528066016,CZPV:0.0445013565):0.0078908838):0.0047723397):0.0082742845,MYZV:0.0937588793):0.0150074260):0.0082853385):0.0106496885,((((LNER:0.0734509971,LWDA:0.0333393661):0.0130418529,(SVVG:0.0272021706,((NHUA:0.0136980140,HENI:0.0489269466):0.0344667223,UZWG:0.0017268862):0.0300674420):0.0146281958):0.0026536248,INSP:0.0260965648):0.0088285640,TJLC:0.0415206561):0.0210700565):0.0046996970,(((((((FEDW:0.0631136358,((KBRW:0.0044096994,DZLN:0.0148831520):0.0094812246,EQYT:0.0035511573):0.0586801888):0.0293200635,((JKNQ:0.0305243879,(UJGI:0.0013489303,TLCA:0.0013747227):0.0020491379):0.0019743254,((((GVCB:0.0096085392,ARYD:0.0000010944):0.0000025686,((IMZV:0.0000000000,IDAU:0.0000000000):0.0000000000,ZINQ:0.0000000000):0.0000010944):0.0000022348,(YHLF:0.0739639163,BZDF:0.0013486256):0.0000012948):0.0000021365,ROLB:0.0000010944):0.0014092719):0.0030599466):0.0108705792,HKMQ:0.0656706044):0.0031441090,(AXNH:0.0000029597,PMTB:0.0000026943):0.0047925780):0.0163814320,SJAN:0.0000021629):0.0624914410,(YNUE:0.0260990991,RJNQ:0.0432070528):0.0282106349):0.0155575194,((Eucgr_v1.1:0.0162736524,(FGDU:0.0183689615,NEBM:0.0039693429):0.0120024228):0.0488462608,(WWQZ:0.0143262304,SWGX:0.0208604695):0.0650559498):0.0206120385):0.0420814806):0.0095923521,((((((((((((((SILJ:0.0521505244,QOXT:0.0595036414):0.1066311096,(VVVV:0.0470854774,DGXS:0.0578688746):0.0077027653):0.0102206479,OCWZ:0.0680419000):0.0263434567,(((((((((((XPAF:0.0381774313,WBIB:0.0291352090):0.0028348764,PWSG:0.0548342517):0.0453456411,CIEA:0.0868688551):0.1018918344,(BSTR:0.1170325239,((((((((ROEI:0.0066014988,Sorbi_v2.1:0.0093724823):0.0070190236,YPIC:0.0173172495):0.0118015586,(SOHV:0.0146311984,(((((ZMGN:0.0051061451,BPKH:0.0063795231):0.0025643901,XUAB:0.0076311203):0.0012609956,NNOK:0.0076310065):0.0012390833,((BXAY:0.0221657336,ZENX:0.0025403156):0.0018176277,VQYB:0.0078699414):0.0049964823):0.0045077397,(WCOR:0.0012675348,XBKS:0.0012529345):0.0049496100):0.0037800673):0.0061870593):0.0132304710,HATH:0.0267423407):0.0068132124,(YXNR:0.0199761192,EFCZ:0.0187877574):0.0093635504):0.0033004751,RCAH:0.0723737189):0.0192727606,(RMVB:0.0159080663,IADP:0.0026201665):0.0133633535):0.0110831920,Orysa_v7.0:0.0347990712):0.0613858520):0.0186748844):0.0385671100,PPQR:0.0652373778):0.0263348246,(JSAG:0.1659265201,BYPY:0.0540942130):0.0153405594):0.0114112846,(((Elagu_v2.0:0.1043129089,HXJE:0.0422544519):0.0000029663,HWUP:0.0298532625):0.0250409568,((((JNUB:0.0585481013,TZNS:0.0331893227):0.0192144415,Musac_v1.0:0.0555821587):0.0030257951,(XHHU:0.0330269315,(LSKK:0.0364405497,UOEL:0.0524668562):0.0084888861):0.0030902960):0.0076052320,(LEMW:0.0205944624,BDJQ:0.0097358513):0.0534113395):0.0800949594):0.0043413113):0.0225708610,((((KYNE:0.0354480984,RDYY:0.1225813660):0.0190082369,YJUG:0.0924299195):0.0255716747,((((SVTS:0.0437700850,MUMD:0.0595698146):0.0005383170,((((((IXEM:0.0609736809,FGRF:0.0399171438):0.0088264657,ONBE:0.0468467868):0.0068288902,(((MVRF:0.0296866624,((UZXL:0.0139423261,HOKG:0.0192480189):0.0064841376,RQZP:0.0030104854):0.0023230200):0.0024053392,((XFJG:0.0038318546,RCUX:0.0038548540):0.0039001905,LSJW:0.0500732325):0.0010948827):0.0155247155,PRFO:0.0211708245):0.0107079479):0.0060527656,(((PLBZ:0.0157186419,ICNN:0.0059529279):0.0043183279,CMCY:0.0092472882):0.0000022743,KXSK:0.0042488143):0.0298762127):0.0014928301,(DMIN:0.0060397576,(LDME:0.0049997032,(JDTY:0.0051538822,DPFW:0.0582441858):0.0227377161):0.0063717446):0.0351381401):0.0017471080,GJPF:0.1116658365):0.0022612662):0.0082049040,SART:0.0397577266):0.0046485087,((EMJJ:0.0513047657,(JVBR:0.0558265667,(WTDE:0.0580103714,(FCEL:0.0270668051,(BLAJ:0.0000026774,JHUL:0.0000030114):0.0227018410):0.0027125611):0.0260726568):0.0112603643):0.0083357770,LTZF:0.0649235378):0.0034774805):0.0070762636):0.0065813951,(((MTHW:0.0505048477,(LELS:0.0302947432,VTUS:0.0226946676):0.0236761627):0.0049662868,XZME:0.0458325299):0.0327768941,THDM:0.0913654134):0.0909291544):0.0133266577):0.0135600006,AFLV:0.0434942914):0.0042644167,OOSO:0.0500513560):0.0027600362,((THEW:0.1676446116,MWYQ:0.0472947530):0.0068921914,(QNPH:0.0159044094,GDKK:0.0138137638):0.0527706605):0.0159350601):0.0218295629):0.0263175227,((COCP:0.1631945478,BYQM:0.0719689848):0.0401213080,(Spipo_v2:0.0784472200,(YMES:0.0258089672,MFIN:0.0445094429):0.0342829716):0.0718760895):0.0191849606):0.0236886880,(((MTII:0.1410797855,(SJEV:0.0949586763,WKSU:0.0339981968):0.0121967430):0.0076279441,(QDVW:0.0898411942,WZFE:0.0408196371):0.0102886232):0.0025509880,((Ambtr_v1.0.27:0.0000021439,URDJ:0.0000028064):0.1145148552,(FZJL:0.0496307717,VZCI:0.0701245292):0.0135295012):0.0239083470):0.0064442215):0.0030040679,(((YZRI:0.1176942880,DHPO:0.0287796041):0.0083311481,(OBPL:0.0587584495,(WBOD:0.0044419951,XQWC:0.0075787507):0.0171543977):0.0054619825):0.0246076359,(((((((XSZI:0.0367634184,MUNP:0.0621698346):0.0906112115,(CSSK:0.0415298548,OPDF:0.0423452169):0.0289514099):0.0743264678,(BSVG:0.0942442269,KRJP:0.0290192672):0.0087021109):0.0058650852,((MAQO:0.0318083954,WAIL:0.0074967748):0.0096490013,ABSS:0.0247860253):0.0080428088):0.0095946624,VYLQ:0.0873773388):0.0073528418,(PAWA:0.1962405349,(FALI:0.0159126236,WPHN:0.0312705122):0.0298248702):0.0104250877):0.0042040773,PZRT:0.2087017503):0.0018981932):0.0080138275):0.0215253973,(((CCID:0.0616849520,((((ZUHO:0.0193741030,UPOG:0.0122204349):0.0654653172,(GBVZ:0.0175653735,Aquco_v1.1:0.0111867715):0.0471879443):0.0132053598,VGHH:0.0548569214):0.0210332265,(WFBF:0.0477437823,YHFG:0.0474774212):0.0075603028):0.0109368359):0.0030875075,((NMGG:0.0870240641,(((EVOD:0.0550789280,(IRAF:0.0306611998,(((RQNK:0.0000010944,QCOU:0.0012546904):0.0113427620,BEKN:0.0061771617):0.0068924706,SSDU:0.0038171148):0.0291046344):0.0065638741):0.0072413846,XMVD:0.0219612373):0.0018048790,XHKT:0.0161178871):0.0124650589):0.0056734126,((UDHA:0.0335723913,ZGQD:0.0138840107):0.0045174480,AUGV:0.0133284632):0.0341729301):0.0118713196):0.0059596159,QTJY:0.0495077595):0.0051721894):0.0036288152,(IWMW:0.0811162817,(((GRRW:0.0696538945,(SIIK:0.0024154161,RQUG:0.0013084040):0.0006469458):0.0450724868,VQFW:0.0839287643):0.0373181812,(AALA:0.0809577504,FAKD:0.0556836545):0.0154134447):0.0040210523):0.0052363926):0.0257242712,(((((((((IUSR:0.1082465282,TOKV:0.0330552146):0.0522732495,(OLXF:0.0961166705,((RJIM:0.0000021371,UVDC:0.0623384761):0.0138556508,QICX:0.0670080607):0.0442834325):0.0267362815):0.0068534093,YKFU:0.0801858830):0.0000022963,(SZPD:0.1160959253,FWBF:0.0647121010):0.0097760219):0.0070858250,XMQO:0.0826550499):0.0055728618,FYTP:0.0358379934):0.0047963440,(PSJT:0.0799044073,VYGG:0.0737111321):0.0000022596):0.0052606424,HAEU:0.0349846142):0.0505503361,QUTB:0.0974125355):0.0000023567):0.0076755675,(((((RXEN:0.0755050544,((FZQN:0.0557618256,OLES:0.0541815190):0.0036916063,(SKNL:0.0253920290,SHEZ:0.1392681686):0.0260813302):0.0096474199):0.0120868857,((((((SMMC:0.0210113928,((AAXJ:0.0061721996,CBJR:0.0146399191):0.0089009508,ONLQ:0.0069585524):0.0204859679):0.0411664055,((XSSD:0.0142553289,WMLW:0.0000021626):0.0342025532,((HDSY:0.0356308406,(((OHKC:0.0344135442,(ZBPY:0.0091349040,BWRK:0.0519737826):0.0048631612):0.0033734309,EYRD:0.0144006132):0.0082144024,PDQH:0.0162908599):0.0082084753):0.0010250952,CUTE:0.0257868856):0.0218455456):0.0068102642):0.0222849056,(WGET:0.0313931672,FVXD:0.0972339125):0.0074945870):0.0284248330,YNFJ:0.1083306304):0.0044207251,((((((BERS:0.0089234518,(GJNX:0.0112212477,(HZTS:0.0014014634,EDIT:0.0041656573):0.0027745365):0.0114766224):0.0000026823,OMYK:0.0808548131):0.0225521305,(MRKX:0.0000021573,BKQU:0.0227523710):0.0178549856):0.0053056037,((LKKX:0.0204035145,(CPKP:0.0268627072,(((UQCB:0.0178554518,EZGR:0.0070809406):0.0110523558,(BYNZ:0.0036510764,KDCH:0.0067379788):0.0715123432):0.0042765280,((LLQV:0.0232903410,LDEL:0.0338454375):0.0000025268,(IWIS:0.0071522395,CPLT:0.0459916774):0.0080593980):0.0029928796):0.0122027794):0.0096194848):0.0101102668,CTYH:0.0844171736):0.0233714885):0.0049769390,BJKT:0.0552769988):0.0069341588,(((AZBL:0.0368408505,SFKQ:0.0615017094):0.0185285607,GIWN:0.0509447588):0.0070123998,(((VJPU:0.0039503991,ZBTA:0.0042642649):0.0072684467,(JGAB:0.0161407186,(DVXD:0.0013174230,EGOS:0.0000025396):0.0110375521):0.0142606321):0.0403834947,JAFJ:0.0185405843):0.0191952597):0.0110261486):0.0159583396):0.0039442763,(KJAA:0.0548715774,NXTS:0.0543717888):0.0168258995):0.0082181584):0.0101863550,RNBN:0.0641843459):0.0311991020,(WQUF:0.0563051876,RUUB:0.0907291057):0.0176694547):0.0064544611,((WOBD:0.0953099727,FYSJ:0.0938230965):0.0153471956,CVDF:0.0631990805):0.0042380977):0.0428051560):0.0066030429,((HQRJ:0.0495030349,NUZN:0.0146248800):0.0236417769,(((((UWFU:0.1034269011,OOVX:0.0240853200):0.0516257580,FCBJ:0.0327260076):0.0018456967,(DAYQ:0.0428686174,UHBY:0.0166004154):0.0133295549):0.0093228211,CKKR:0.0270238258):0.0271626324,HTIP:0.0454207707):0.0072335419):0.0067459521):0.0035388333,(NNGU:0.1006374904,(Vitvi_Genoscope.12X:0.0134233957,BGZG:0.0512646019):0.0401806961):0.0144247160):0.0032102011,((((((((((((CWYJ:0.0141008621,TQKZ:0.0202220688):0.0687034855,WEQK:0.0624638119):0.0182516920,AJFN:0.0296805681):0.0027932058,(((EDBB:0.0399384227,SUVN:0.0393648983):0.0042989201,(SALZ:0.0028070150,NTEO:0.0476597573):0.0674358568):0.0064986129,OINM:0.0623302271):0.0035436975):0.0127230862,((MVSE:0.0085416244,QIKZ:0.0221951431):0.0129765979,((QURC:0.0411945553,ZETY:0.0284794626):0.0136626415,(SZUO:0.0922822191,UDUT:0.1171450213):0.0732757873):0.0735235326):0.0256115146):0.0208567345,(((IZLO:0.0312420835,IHPC:0.0420365760):0.0155410507,(TPEM:0.0237743905,(AUIP:0.0673805981,((((((UYED:0.6916839696,(((OAGK:0.0114848024,GUMF:0.0057139934):0.0120166689,DUQG:0.0327110447):0.0343392476,BMSE:0.0609534441):0.0087395785):0.0000029902,(EYKJ:0.0179767436,((TEZA:0.0508869633,XRCX:0.0000010944):0.0000028962,(MHYG:0.0088932607,DESP:0.0000028423):0.0090047448):0.0396606604):0.0107956142):0.0071386989,(DUNJ:0.0545408736,((((((HXCD:0.0070175457,BIDT:0.0248732916):0.0000021139,QXWF:0.0021434118):0.0022041914,((NVSO:0.0168521916,ZCUA:0.0040289578):0.0040341573,JEPE:0.0278861872):0.0014756705):0.0061674077,JYMN:0.0094033967):0.0047318497,((AQZD:0.0339735198,YRBQ:0.0000010944):0.0000026640,LYPZ:0.0039048768):0.0063671632):0.0171574657,UBLN:0.0544271202):0.0078045053):0.0173285107):0.0035124188,((CFRN:0.0462509574,JNKW:0.0296227822):0.0318965439,((DDRL:0.0105396533,(KGJF:0.0000020939,KFZY:0.0023537457):0.0078790442):0.0012967344,FUPX:0.0014442305):0.0339601850):0.0137497484):0.0335551801,HUQC:0.1035077367):0.0090888792,FXGI:0.1264868266):0.0280278357):0.0090100862):0.0085509322):0.0185186070,SBZH:0.0525276765):0.0000024855):0.0093879143,(((FFFY:0.1220630211,(JTRM:0.0360449223,(CAQZ:0.0173514143,GSZA:0.0475394083):0.0137744051):0.0040733992):0.0210058350,(ZJRC:0.0391453302,HLJG:0.0288722433):0.0101280996):0.0016545732,QACK:0.0619755750):0.0020260562):0.0092590304,(CLNU:0.0130171809,CLMX:0.0052698785):0.0583774265):0.0237856340,VUSY:0.0398157574):0.0137406835,((VTLJ:0.0579738991,BFJL:0.0458923798):0.0065725406,GIPR:0.0327002870):0.0069872428):0.0011717666,((((((((ADHK:0.0896852158,((WRPP:0.0203519523,BEFC:0.0026555687):0.0081622986,OXYP:0.0266729613):0.0271505034):0.0055500094,NGRR:0.0458370512):0.0016775092,UFHF:0.0523891558):0.0052773110,JPDJ:0.0568306132):0.0022026043,((AVJK:0.0653636176,(LRTN:0.0538738546,SERM:0.0565269396):0.0119573734):0.0022007460,(PPPZ:0.0360572705,WXVX:0.0023179038):0.0121187923):0.0174033771):0.0032767786,YSRZ:0.0534203676):0.0000021805,((HUSX:0.1169564435,QAUE:0.0437686673):0.0008354946,(((PTFA:0.0597271136,(ODDO:0.0123606558,DAAD:0.0183289397):0.0466501105):0.0156717715,VKJD:0.0294323346):0.0177526946,KVFU:0.0756570681):0.0035098758):0.0075318066):0.0133057615,(FNEN:0.1175710622,(JEXA:0.0637852572,BNTL:0.0342607727):0.0310377741):0.0207120181):0.0129618399):0.0113489929,((((((((DCCI:0.0600474976,(DTNC:0.0418055205,RWKR:0.0560743033):0.0061477922):0.0116033795,CLRW:0.0457093290):0.0057249073,(((((((((((((((TKEK:0.0375068096,FROP:0.0646608916):0.0153278062,FDMM:0.0640376506):0.0157008732,SNNC:0.0213545456):0.0025198365,RTNA:0.0253940957):0.0505837706,ATYL:0.0297304123):0.0000022709,DMLT:0.0190907309):0.0034029852,GNPX:0.0275224838):0.0166138422,(GETL:0.0443948098,EAAA:0.0364614179):0.0071234436):0.0114887955,(UCNM:0.0439968553,((EQDA:0.0283063250,BAHE:0.0514040699):0.0104958376,((((PUCW:0.0095243705,FUMQ:0.0438085507):0.0132719007,FYUH:0.0449033199):0.0079266117,TAGM:0.0373298928):0.0040337290,((IYDF:0.0254900384,(WHNV:0.0079952285,XMBA:0.0576284787):0.0230264592):0.0000021696,PHCE:0.0211509217):0.0157041331):0.0000028901):0.0313380652):0.0268638326):0.0107119967,FAMO:0.0555166229):0.0067421874,(EJBY:0.0507974534,((GRFT:0.0051370633,XRLM:0.0025789088):0.0141099729,(XXYA:0.0000030304,SIBR:0.0056761601):0.0365919995):0.0154054477):0.0216709129):0.0026034460,UTQR:0.0302105440):0.0099723361,(GCFE:0.0244932285,(MQIV:0.0168873008,PSHB:0.0920980245):0.0204872348):0.0148790448):0.0213561849,((((((ZRIN:0.0509334256,UMUL:0.0136632730):0.0914206598,EJCM:0.0360079434):0.0160670535,VYDM:0.0255102381):0.0069746209,OWAS:0.0514513570):0.0054811489,Mimgu_v2.0:0.0603836750):0.0079244266,((((HRUR:0.0970905251,(JCMU:0.0030095223,MXFG:0.0072010687):0.0393732174):0.0278918350,GDZS:0.0424484269):0.0095497458,(PCGJ:0.0558026783,(NBMW:0.0160060203,AYIY:0.0625790469):0.0372509966):0.0357862638):0.0059129378,(WOHL:0.0164143127,EDXZ:0.0071622688):0.0147051745):0.0066145844):0.0039575737):0.0108022249,((GNRI:0.0284861576,(PTBJ:0.0341332754,YKZB:0.0189743454):0.0198303685):0.0316636351,YRHD:0.0378106506):0.0099256359):0.0182992363):0.0212071583,(((TORX:0.0090208033,KTAR:0.0110652899):0.0118306693,MZLD:0.0106023936):0.0297634647,COBX:0.0864169882):0.0162167400):0.0326483545,((((ZSGF:0.0279085624,(KTWL:0.0133782976,PCNH:0.0836795682):0.0478103385):0.0000029045,WQRD:0.1041989422):0.0678091048,(((((QEHE:0.0249710358,UOYN:0.0296433424):0.0190841972,MGVU:0.1010077358):0.0198995826,((((YADI:0.0047792042,DSUV:0.0140568974):0.0524415322,JGYZ:0.0321058772):0.0018430938,(JCLQ:0.0000023681,YFQX:0.0000022544):0.0253609563):0.0073995320,EDEQ:0.0247397527):0.0098387260):0.0060699481,GGJD:0.0712767140):0.0064672402,(KPUM:0.0550669280,ECTD:0.0458540078):0.0408228644):0.0231969371):0.0165120143,(((((((NHAG:0.0041924019,IZNU:0.0065959165):0.0425953855,NAUM:0.1080780794):0.0524955702,(QSLH:0.0014354013,(OQBM:0.0012755661,XQRV:0.0000023436):0.0000021348):0.0012734914):0.0000016478,VXKB:0.0000010944):0.0000022831,(EMBR:0.0043384711,ALUC:0.0000010944):0.0000024016):0.0054611387,(AHRN:0.0629211631,CPOC:0.0076981017):0.0079536911):0.0507986828,(((((OSMU:0.0203062193,LWCK:0.0787524806):0.0075130231,AIOU:0.0177219998):0.0000027184,(Solly_iTAGv2.3:0.0075553938,(GHLP:0.0079410717,(LQJY:0.0088679176,NMDZ:0.0027407414):0.0000028333):0.0025597902):0.0014706367):0.0038748999,BOLZ:0.0129814822):0.0090943608,MKZR:0.0220390935):0.0460170910):0.0246784444):0.0142323112):0.0138265244,(YQIJ:0.0597415586,(((((XVRU:0.0120579127,DIHD:0.0175982068):0.0038634406,(JWEY:0.0272546785,((((MZOB:0.0132012215,(IDGE:0.0045950213,OUER:0.0000022834):0.0046734251):0.0012846262,NIGS:0.0090105978):0.0000029625,(OEKO:0.0041850004,MDJK:0.0025838719):0.0020815599):0.0000020819,ABEH:0.0130392555):0.0000026778):0.0053333100):0.0335376207,(SMUR:0.0335313642,HANM:0.0131439673):0.0342795476):0.0073144521,EMAL:0.0393550414):0.0178590758,DKFZ:0.0836578526):0.0037647612):0.0146051128):0.0381060643,(DFYF:0.0072957282,(SXML:0.0000010944,ASMV:0.0086384656):0.0024247491):0.0320938251):0.0192177935,OCTM:0.1565560756):0.0061837120):0.0101790569):0.0090699518):0.0000028780,(QKMG:0.0722942122,(XGFU:0.0376479607,(RSPO:0.0227886357,BSEY:0.0417580507):0.0127388802):0.0106841831):0.0450084843):0.0303460644,EHNF:0.0772306625):0.0821529105,ROZZ:0.5460373231):0.0466970721,(TOXE:0.1564044084,(((((((((QFAE:0.0709969480,NRXL:0.0099395230):0.0028230158,XIRK:0.0119015076):0.0094459303,(((((((CGDN:0.0093880053,(FRPM:0.0052847871,(BUWV:0.0073251027,XQSG:0.0000023707):0.0055030096):0.0018141359):0.0012943884,(QNGJ:0.0056762685,XMGP:0.0377806379):0.0053588714):0.0033820940,(UEVI:0.0093018045,AIGO:0.0092851012):0.0025984407):0.0000026213,(VFYZ:0.0065700286,NKIN:0.0256518739):0.0000025814):0.0038601129,(YYPE:0.0039597975,(AUDE:0.0079046425,((ETCJ:0.0092319696,(RMMV:0.0079610510,(JDQB:0.0010657842,IFLI:0.0193927747):0.0000976307):0.0320396917):0.0000026116,GKCZ:0.0079075158):0.0000027823):0.0013194799):0.0110056891):0.0028881686,FHST:0.0130865082):0.0027589251,(GMHZ:0.0059241498,OXGJ:0.0313549399):0.0000021800):0.0067993121):0.0081639414,ZQVF:0.0170474315):0.0037481835,QSNJ:0.0161241104):0.0295876960,(NVGZ:0.0316662658,((IAJW:0.0102519487,(EFMS:0.0000030677,HQOM:0.0281383378):0.0118757284):0.0156860732,((ZYAX:0.0034350635,WWSS:0.0037588401):0.0070502995,(YLPM:0.0000022779,BTTS:0.0676103932):0.0166370348):0.0139341900):0.0020036056):0.0020354176):0.0112019725,(((XTZO:0.0303396299,(ACWS:0.0000030229,MIXZ:0.0042025737):0.0101231182):0.0026895818,RSCE:0.0128829315):0.0344307020,(((((CDFR:0.0267762082,ZQWM:0.0090300744):0.0053945333,(((((FMWZ:0.0405081079,(IZGN:0.0020985191,(ROWR:0.0000010944,QHBI:0.0000010944):0.0039931351):0.0018337650):0.0101912475,(((VGSX:0.0079579253,XLGK:0.0053038006):0.0013139687,UUJS:0.0052835948):0.0013700658,SCEB:0.0167420261):0.0013154237):0.0000028012,QCGM:0.0263786526):0.0044583833,HILW:0.0075831892):0.0062278479,MHGD:0.0149492841):0.0054597794):0.0015226846,JRNA:0.0173911894):0.0014501389,BBDD:0.0214401051):0.0014821027,((OWFC:0.0229832220,KLGF:0.0650208016):0.0033370676,(JZVE:0.0209877632,EGLZ:0.0690040977):0.0039632544):0.0037643309):0.0359678174):0.0078793045):0.0038218808,YFZK:0.0390058023):0.0354105223,(((AWQB:0.0238439882,(IIOL:0.0145483661,(DZQM:0.0057847689,(MFTM:0.0000023164,JBND:0.0035340558):0.0020491461):0.0103387559):0.0143021498):0.0136031944,(IOVS:0.0093304889,WVWN:0.0081499545):0.0162444080):0.0050902966,(((GGEA:0.0440264696,VSRH:0.0296597963):0.0070703950,AQFM:0.0879678593):0.0035391821,GAMH:0.0197083098):0.0065115600):0.0261805942):0.0416904871):0.0400945489):0.0958769754,((((((KIIX:0.1490151051,(((((((((ZQYU:0.2159675654,(ORJE:0.0202880697,(UJWU:0.0254659687,((CJNT:0.0060055120,GYFU:0.0015281066):0.0016357404,YLJA:0.0038178021):0.0029392271):0.0046066356):0.0021313634):0.0283326947,NWWI:0.0679358731):0.0104821678,(((BMIF:0.0479185201,KJZG:0.0594661186):0.0325467591,(WGTU:0.0373828389,FQGQ:0.0254559834):0.0200749738):0.0027878992,(HTFH:0.0334317784,(UFJN:0.0156721016,OCZL:0.0215692857):0.0004748898):0.0163390740):0.0226896964):0.0013004741,(YOWV:0.0020334551,XXHP:0.0020147727):0.0254159809):0.0157556531,((PIVW:0.1172398597,(UJTT:0.0679733293,(POPJ:0.0198873018,FLTD:0.0515078709):0.0032319192):0.0228555861):0.0123342045,(((NDUV:0.0935689534,WCLG:0.0268384920):0.0152721547,(((DCDT:0.0545308787,(ZXJO:0.0108860879,YCKE:0.1183393095):0.0308565108):0.0042672852,XDDT:0.0464062942):0.0370588396,BMJR:0.0629073245):0.0178596242):0.0185571449,WQML:0.0531669445):0.0015841691):0.0466681116):0.0200907241,VVRN:0.1146623720):0.0123976281,(NOKI:0.0013114846,YIXP:0.0047544232):0.0842401479):0.0207294323,(((CVEG:0.1320043112,QIAD:0.1920826002):0.0493267937,PNZO:0.0490900680):0.0525359584,UWOD:0.0230071973):0.0041296249):0.0266905126,EWXK:0.0333021906):0.0122507810):0.0645419027,MEKP:0.0936999914):0.0198620881,PBUU:0.1001222082):0.0567895645,(VIBO:0.0018601155,(RFMZ:0.0081962680,UOMY:0.0047540276):0.0103420867):0.0797849271):0.0465073172,((QVMR:0.0457014641,ALVQ:0.0161413169):0.0791513958,(EEAQ:0.0486403568,DJSE:0.1238482602):0.0562217301):0.0340205479):0.0381608364,JVSZ:0.3002769080):0.0186775389):0.0327868724):0.0277152252,PKOX:0.2118518610):0.0626138054):0.0435535561,((ZRMT:0.3774733625,HKZW:0.3001578574):0.0643176241,(JOJQ:0.1365724485,(MFZO:0.0813883463,NBYP:0.0599315185):0.0227439200):0.0657440984):0.0456235647):0.0119759710,((BFIK:0.2100941501,(FPCO:0.1333089389,(Klefl_v1.0:0.0310784778,FQLP:0.0422572345):0.1626245724):0.1710693588):0.0342621968,TPHT:0.2288616366):0.0142509003):0.0199622445,(WDCW:0.2101325929,(YOXI:0.0607228187,(RPGL:0.1115422160,VAZE:0.0201561724):0.0454207630):0.1849285992):0.0369806327):0.0316848543,(FFGR:0.3459331258,WSJO:0.2180792583):0.0464869522):0.0232768954,((Chabr_v0.1:0.4055559104,(QPDY:0.2279213167,VQBJ:0.1681901465):0.0766756818):0.0436622496,((((((((MNNM:0.1814216588,QWFV:0.1683515320):0.0336214932,(GGWH:0.1452785380,DFDS:0.1840311506):0.0135251925):0.0335792635,((RPRU:0.1656591785,WCQU:0.1344242168):0.0299847503,BHBK:0.1620538529):0.0313685610):0.0284451787,(GYRP:0.2059983131,RQFE:0.2554023039):0.0363436802):0.0202454731,((WDGV:0.1684027926,(MOYY:0.1233009948,STKJ:0.3312737280):0.0265189720):0.0152551794,GBGT:0.1342865748):0.0267237512):0.0291280041,(RPQV:0.1806944343,YSQT:0.1587899150):0.0274231102):0.0519798087,(DRFX:0.3624850201,AEKF:0.2073875594):0.0164160110):0.0372483973,NNHQ:0.5703039270):0.0201793232):0.0402951786):0.0303980668,(AZZW:0.3073405856,XRTZ:0.4562022712):0.0516484591):0.0398415034,(((((NPRL:0.3377467037,OQWW:0.3956423596):0.2378088721,(JQFK:0.5097816518,((NMAK:0.0382743998,VJED:0.0102498419):0.0200997959,(LLEN:0.0000021304,RFAD:0.0000010944):0.0296075809):0.2865083464):0.0695120271):0.0362859388,(IJMT:0.2924609667,WCZB:0.7374284408):0.1186404591):0.0326371523,((PQED:0.2442691055,(POOW:0.2541256074,QFND:0.1377627509):0.0557877800):0.0283249408,JKHA:0.1858212732):0.0841685185):0.0670380644,(((Ostlu_v2.0:0.4423898620,Micpu_v3.0:0.3641455255):0.1828388103,XMCL:0.2797524795):0.0531200996,(((DVYE:0.6096568718,VHIJ:0.4273180937):0.1407177767,KYIO:0.2505090861):0.0303793715,UYFR:0.5755965343):0.0541792763):0.0200345792):0.0160731706):0.0349926251):0.0307977800):0.0044568053,MMKU:0.3026290607):0.0302737932,(((JMTE:0.3808229053,XJGM:0.4000194269):0.0720972655,HYHN:0.4061050979):0.1178378218,QWRA:0.4355941101):0.0594569385):0.0118013319,((((((PVGP:0.0677976566,OBUY:0.2429606263):0.3991893399,((YSBD:0.6326683384,(((PYDB:0.0769383292,(URSB:0.0655837952,(IKIZ:0.0280938063,ZJOJ:0.0311063836):0.0263359465):0.0550677110):0.2930986962,((((((SBLT:0.0616047601,(WEJN:0.0147555641,UGPM:0.0111778231):0.0490069624):0.0679019299,CKXF:0.1718314943):0.1170755704,IEHF:0.2106609160):0.0176010800,((JEBK:0.0395063704,BWVJ:0.0324335108):0.1554468776,IHJY:0.5088514165):0.0911938307):0.0589845223,VZWX:0.3137769276):0.0305595481,(((IKWM:0.0980228412,(LJPN:0.0000010944,VNAL:0.0000010944):0.0843322705):0.0591699176,FTRP:0.1323592784):0.0958517393,PWKQ:0.2805579591):0.2107010873):0.0677376214):0.0999696315,ZULJ:0.4209898502):0.0624230904):0.0364825484,RTLC:0.4552142861):0.0234166915):0.0404691142,(LLXJ:0.3190776565,(JJZR:0.2331241756,RSOF:0.3224196100):0.0469912262):0.0595697336):0.0691984759,((RRSV:0.2162849121,PUAN:0.4817189856):0.3063209405,BAJW:0.3821323581):0.0522086677):0.0425136725,KADG:0.4582094979):0.0170205484,GYBH:0.5580368061):0.0153705692):0.0630320113,(BOGT:0.9269854475,XOAL:0.4457187623):0.0578842522):0.2553185927,(((FOMH:0.0243205985,HFIK:0.0250198470):0.0343282054,(VYER:0.1089412736,JGGD:0.0319633952):0.0696995444):0.1136719860,((QLMZ:0.0542199620,ULXR:0.0860545180):0.1101892141,ASZK:0.1072817500):0.0791171037):0.0254987742):0.0487590770,FSQE:0.1339569021):0.0221567704,APTP:0.3389292774):0.3419266433);
```
Create a new .txt file that contains this information from all of the MUSCLE-aligned genes named muscle-ml-8.gene.tre and another with all the Newick data for the KALIGN-aligned genes named kalign-ml-8.gene.tre. Save these in ~/botany-project/data/muscle/coestimation/ and ~/botany-project/data/kalign/coestimation/ respectively.

1. Install ASTRAL
Install ASTRAL by downloading the zip file and extracting the contents. I downloaded the software from GitHub at https://github.com/smirarab/ASTRAL/blob/master/README.md#installation and extracted the file to ~/Downloads/Programs using the file manager. ASTRAL can be run from any directory with:
```
java -jar ~/Downloads/Programs/Astral.5.7.8/Astral/astral.5.7.8.jar
```
2. Perform coestimation of MUSCLE-aligned genes.
```
java -jar ~/Downloads/Programs/Astral.5.7.8/Astral/astral.5.7.8.jar -i ~/botany-project/data/muscle/coestimation/muscle-ml-8.gene.tre -o muscle-species.tre 2>muscle-species.tre
```
    The output will be recorded in muscle-species.tre and include Newick data for the final species tree.
```

================== ASTRAL ===================== 

This is ASTRAL version 5.7.8
Gene trees are treated as unrooted
8 trees read from /home/laurel/botany-project/muscle-ml-8.gene.tre
index0
All output trees will be *arbitrarily* rooted at VRGZ

======== Running the main analysis
Number of taxa: 1176 (1176 species)
Taxa: [VRGZ, YRMA, RWXW, WCLV, ISIM, NQYP, LSHT, KSFK, JIWJ, GJIY, OQON, CQQP, OAEZ, AJAU, Cyame_v1.0, VFIV, RAWF, JMUI, BAZF, DZPJ, ALZF, NSTT, FXHG, BZSH, XKWQ, VJDZ, MWAN, PZIF, PRIQ, QRTH, ZLQE, KFEB, GUBD, ZFXU, SYJM, UKUC, WDWX, TSBQ, JRGZ, ZIVZ, LNIL, MULF, JKKI, ISPU, Chlre_v5.5, Volca_v2.0, WRSL, JWGT, AOUJ, AJUW, LBRP, FOYQ, XDLL, VIAU, BILC, EEJO, Chlva_v1.0, NKXU, MFYC, MNCB, AKCR, AYPS, DUMA, HVNO, FMVB, EATP, GXBM, ZNUM, XIVI, TNAW, XAXW, SNOX, ZZOL, Selmo_v1.0, ABIJ, LGDQ, KJYC, JKAA, ZYCD, ZFGK, KUXM, UPMJ, GKAG, YHZW, CBAE, ZZEI, PQTO, XNXF, ENQF, WAFT, BNCU, NRWZ, YFGP, JHFI, PIUF, KRUQ, UUHD, RTMU, OFTV, LGOW, YBQN, IRBN, TGKW, HERT, HMHL, TFYI, TXVB, WJLO, ILBQ, JPYU, BSNI, TWUW, FAJB, WEEQ, DXOU, TCBC, AKXB, UCRN, RCBT, GOWD, UHLI, WOGB, HRWG, HVBQ, AWOI, Phypa_v3.0, YEPO, KEFD, XWHK, JMXW, CMEQ, DHWX, ZACW, IGUH, VBMM, TMAJ, QMWB, WSPM, QKQO, LNSF, EEMJ, MIRS, TAVP, JADL, WNGH, YWNF, ORKS, BGXB, VMXJ, RGKI, NGTD, ABCD, RDOO, GRKU, BPSG, FFPD, SZYG, ZTHV, SKQD, QZZU, RKLL, RKFX, JETM, Phavu_v1.0, KEGA, TVSH, SUAK, KNMB, MYMP, HJMP, JTQQ, VLNB, SLYR, ZSSR, ZCDJ, PJSX, KZED, VHZV, GEHT, QZXQ, CMFF, OQHZ, OHAE, PTLU, KVAY, ZHMB, WVEF, EILE, WKCY, EDHN, XVJB, KYAD, BJSW, AQGE, RBYC, IANR, Prupe_v1.0, VCIN, EAVM, SQCF, XFFT, HDWF, OBTI, FCCA, BCAA, YUOM, VFFP, HBHB, WAXR, IKFD, TIUZ, THHD, YZVJ, JHCN, PXYR, RHAU, PAZJ, Manes_v4.1, XNLP, VVPY, YGAT, Poptr_v3.0, IEPQ, INQX, LFOG, TDTF, GLVK, KKDQ, VNMY, COAQ, EZZT, SIZE, TXMP, VXOD, KCPT, HNCF, AEPI, BHYC, BVOF, MYVH, POZS, OODC, CKDK, BNDE, FWCQ, RPPC, ZTLR, LPGY, NJLF, HBUQ, ZBVT, XPBC, IHCQ, PKMO, WMUK, KPTE, ATFX, ZJUL, DRIL, KWGC, PUDI, AWJM, Theca_v1.1, YGCX, DZTK, SWPE, LAPO, VMNH, CSUV, Arath_TAIR10, TZWR, QSKP, LVUS, VDKG, UPZX, RTTY, UAXP, CRNC, Carpa_ASGPBv0.4, CZPV, MYZV, LNER, LWDA, SVVG, NHUA, HENI, UZWG, INSP, TJLC, FEDW, KBRW, DZLN, EQYT, JKNQ, UJGI, TLCA, GVCB, ARYD, IMZV, IDAU, ZINQ, YHLF, BZDF, ROLB, HKMQ, AXNH, PMTB, SJAN, YNUE, RJNQ, Eucgr_v1.1, FGDU, NEBM, WWQZ, SWGX, SILJ, QOXT, VVVV, DGXS, OCWZ, XPAF, WBIB, PWSG, CIEA, BSTR, ROEI, Sorbi_v2.1, YPIC, SOHV, ZMGN, BPKH, XUAB, NNOK, BXAY, ZENX, VQYB, WCOR, XBKS, HATH, YXNR, EFCZ, RCAH, RMVB, IADP, Orysa_v7.0, PPQR, JSAG, BYPY, Elagu_v2.0, HXJE, HWUP, JNUB, TZNS, Musac_v1.0, XHHU, LSKK, UOEL, LEMW, BDJQ, KYNE, RDYY, YJUG, SVTS, MUMD, IXEM, FGRF, ONBE, MVRF, UZXL, HOKG, RQZP, XFJG, RCUX, LSJW, PRFO, PLBZ, ICNN, CMCY, KXSK, DMIN, LDME, JDTY, DPFW, GJPF, SART, EMJJ, JVBR, WTDE, FCEL, BLAJ, JHUL, LTZF, MTHW, LELS, VTUS, XZME, THDM, AFLV, OOSO, THEW, MWYQ, QNPH, GDKK, COCP, BYQM, Spipo_v2, YMES, MFIN, MTII, SJEV, WKSU, QDVW, WZFE, Ambtr_v1.0.27, URDJ, FZJL, VZCI, YZRI, DHPO, OBPL, WBOD, XQWC, XSZI, MUNP, CSSK, OPDF, BSVG, KRJP, MAQO, WAIL, ABSS, VYLQ, PAWA, FALI, WPHN, PZRT, CCID, ZUHO, UPOG, GBVZ, Aquco_v1.1, VGHH, WFBF, YHFG, NMGG, EVOD, IRAF, RQNK, QCOU, BEKN, SSDU, XMVD, XHKT, UDHA, ZGQD, AUGV, QTJY, IWMW, GRRW, SIIK, RQUG, VQFW, AALA, FAKD, IUSR, TOKV, OLXF, RJIM, UVDC, QICX, YKFU, SZPD, FWBF, XMQO, FYTP, PSJT, VYGG, HAEU, QUTB, RXEN, FZQN, OLES, SKNL, SHEZ, SMMC, AAXJ, CBJR, ONLQ, XSSD, WMLW, HDSY, OHKC, ZBPY, BWRK, EYRD, PDQH, CUTE, WGET, FVXD, YNFJ, BERS, GJNX, HZTS, EDIT, OMYK, MRKX, BKQU, LKKX, CPKP, UQCB, EZGR, BYNZ, KDCH, LLQV, LDEL, IWIS, CPLT, CTYH, BJKT, AZBL, SFKQ, GIWN, VJPU, ZBTA, JGAB, DVXD, EGOS, JAFJ, KJAA, NXTS, RNBN, WQUF, RUUB, WOBD, FYSJ, CVDF, HQRJ, NUZN, UWFU, OOVX, FCBJ, DAYQ, UHBY, CKKR, HTIP, NNGU, Vitvi_Genoscope.12X, BGZG, CWYJ, TQKZ, WEQK, AJFN, EDBB, SUVN, SALZ, NTEO, OINM, MVSE, QIKZ, QURC, ZETY, SZUO, UDUT, IZLO, IHPC, TPEM, AUIP, UYED, OAGK, GUMF, DUQG, BMSE, EYKJ, TEZA, XRCX, MHYG, DESP, DUNJ, HXCD, BIDT, QXWF, NVSO, ZCUA, JEPE, JYMN, AQZD, YRBQ, LYPZ, UBLN, CFRN, JNKW, DDRL, KGJF, KFZY, FUPX, HUQC, FXGI, SBZH, FFFY, JTRM, CAQZ, GSZA, ZJRC, HLJG, QACK, CLNU, CLMX, VUSY, VTLJ, BFJL, GIPR, ADHK, WRPP, BEFC, OXYP, NGRR, UFHF, JPDJ, AVJK, LRTN, SERM, PPPZ, WXVX, YSRZ, HUSX, QAUE, PTFA, ODDO, DAAD, VKJD, KVFU, FNEN, JEXA, BNTL, DCCI, DTNC, RWKR, CLRW, TKEK, FROP, FDMM, SNNC, RTNA, ATYL, DMLT, GNPX, GETL, EAAA, UCNM, EQDA, BAHE, PUCW, FUMQ, FYUH, TAGM, IYDF, WHNV, XMBA, PHCE, FAMO, EJBY, GRFT, XRLM, XXYA, SIBR, UTQR, GCFE, MQIV, PSHB, ZRIN, UMUL, EJCM, VYDM, OWAS, Mimgu_v2.0, HRUR, JCMU, MXFG, GDZS, PCGJ, NBMW, AYIY, WOHL, EDXZ, GNRI, PTBJ, YKZB, YRHD, TORX, KTAR, MZLD, COBX, ZSGF, KTWL, PCNH, WQRD, QEHE, UOYN, MGVU, YADI, DSUV, JGYZ, JCLQ, YFQX, EDEQ, GGJD, KPUM, ECTD, NHAG, IZNU, NAUM, QSLH, OQBM, XQRV, VXKB, EMBR, ALUC, AHRN, CPOC, OSMU, LWCK, AIOU, Solly_iTAGv2.3, GHLP, LQJY, NMDZ, BOLZ, MKZR, YQIJ, XVRU, DIHD, JWEY, MZOB, IDGE, OUER, NIGS, OEKO, MDJK, ABEH, SMUR, HANM, EMAL, DKFZ, DFYF, SXML, ASMV, OCTM, QKMG, XGFU, RSPO, BSEY, EHNF, ROZZ, TOXE, QFAE, NRXL, XIRK, CGDN, FRPM, BUWV, XQSG, QNGJ, XMGP, UEVI, AIGO, VFYZ, NKIN, YYPE, AUDE, ETCJ, RMMV, JDQB, IFLI, GKCZ, FHST, GMHZ, OXGJ, ZQVF, QSNJ, NVGZ, IAJW, EFMS, HQOM, ZYAX, WWSS, YLPM, BTTS, XTZO, ACWS, MIXZ, RSCE, CDFR, ZQWM, FMWZ, IZGN, ROWR, QHBI, VGSX, XLGK, UUJS, SCEB, QCGM, HILW, MHGD, JRNA, BBDD, OWFC, KLGF, JZVE, EGLZ, YFZK, AWQB, IIOL, DZQM, MFTM, JBND, IOVS, WVWN, GGEA, VSRH, AQFM, GAMH, KIIX, ZQYU, ORJE, UJWU, CJNT, GYFU, YLJA, NWWI, BMIF, KJZG, WGTU, FQGQ, HTFH, UFJN, OCZL, YOWV, XXHP, PIVW, UJTT, POPJ, FLTD, NDUV, WCLG, DCDT, ZXJO, YCKE, XDDT, BMJR, WQML, VVRN, NOKI, YIXP, CVEG, QIAD, PNZO, UWOD, EWXK, MEKP, PBUU, VIBO, RFMZ, UOMY, QVMR, ALVQ, EEAQ, DJSE, JVSZ, PKOX, ZRMT, HKZW, JOJQ, MFZO, NBYP, BFIK, FPCO, Klefl_v1.0, FQLP, TPHT, WDCW, YOXI, RPGL, VAZE, FFGR, WSJO, Chabr_v0.1, QPDY, VQBJ, MNNM, QWFV, GGWH, DFDS, RPRU, WCQU, BHBK, GYRP, RQFE, WDGV, MOYY, STKJ, GBGT, RPQV, YSQT, DRFX, AEKF, NNHQ, AZZW, XRTZ, NPRL, OQWW, JQFK, NMAK, VJED, LLEN, RFAD, IJMT, WCZB, PQED, POOW, QFND, JKHA, Ostlu_v2.0, Micpu_v3.0, XMCL, DVYE, VHIJ, KYIO, UYFR, MMKU, JMTE, XJGM, HYHN, QWRA, PVGP, OBUY, YSBD, PYDB, URSB, IKIZ, ZJOJ, SBLT, WEJN, UGPM, CKXF, IEHF, JEBK, BWVJ, IHJY, VZWX, IKWM, LJPN, VNAL, FTRP, PWKQ, ZULJ, RTLC, LLXJ, JJZR, RSOF, RRSV, PUAN, BAJW, KADG, GYBH, BOGT, XOAL, FOMH, HFIK, VYER, JGGD, QLMZ, ULXR, ASZK, FSQE, APTP, MCPK, TGNL, ZDIZ, GFUR, SDPC, EBWI, JTIG, HAOX, PYHZ, TRRQ, XAYK, LRRR, WEAC, HGSM, JNVS, RFSD, IXVJ, NFXV, RVGH, HYZL, HABV, CJGZ, TTRG, ACFP, TJQY, ZPKK, DXQW, CWZU, VGVI, SCAO, Betvu_v1.1, Nelnu_v1.0, XZUY, GNQG, KAWQ, GTHK, UGNK, DFHO, MTGC, SKYV, GSXD, YQEC, IXLH, HNDZ, LHLE, HEGQ, VITX, YJJY, CQPW, WZYK, HPXA, ANON, ZQRI, ISHC, LIRF, FIKG, FIDQ, JCXF, YDCQ, NATT, BTFM, EGNB, DRGY, ULKT, YLWW, BRUD, NHIX, VBHQ, TCYS, KOFB, WIGA, BCGB, OSHQ, UHJR, JSZD, AYMT, WYIG, ZSAB, OSIP, RZTJ, ZIWB, NPND, NXOH, PEZP, DOVJ, AFQQ, DLJZ, DLAI, UGJI, ERWT, OMDH, YKQR, XKPS, TJES, GCYL, QNOC, SXCE, QIEH, IPWB, ERIA, CIAC, XISJ, AXAF, WWKL, SWOH, WTKZ, NWMY, VDAO, Pinta_v2.0, JUWL, AREG, OVIJ, SGTW, MROH, RICC, FCHS, AFPO, URCP, PSKY, JBLI, RFRB, GANB, BEGM, CAPN, KEYW, MCHJ, HIDG, RAPY, VKVG, IRZA, VBLH, PFUD, QYXY, LETF, ENAU, OTQG, XTON, RHVC, OFUE, VALZ, RNAT, RYJX, KUJU, IHOI, UTRE, ISGT, IRYH, MXDS, GAON, WXNT, YBML, NSPR, BXBF, TVCU, AJJE, CQMG, ZHEE, CTSS, SLOI, SYHW, AXBO, WLIC, ACRY, USIX, BCYF, WXRI, ETGN, BAKF, IAYV, SRSQ, HHXJ, CBNG, LDRY, QXSZ, MWXT, HBGV, YZGX, JLOV, RMWJ, XOOE, AXPJ, JOIS, KBXS, ZKPF, DWZT, KAYP, GIOY, BLVL, DHAW, NCVK, BBBA, GTSV, LJQF, MFEA, LXRN, NHCM, KMNX, WVMY, DBYD, RUIF]
Taxon occupancy: {POOW=6, SVTS=5, KSFK=6, EBWI=3, AIOU=8, XFFT=7, KUJU=4, JRGZ=5, PQTO=8, JAFJ=8, RQNK=8, POPJ=8, UZWG=8, WKSU=7, RWXW=6, OCZL=7, RSPO=4, SIZE=8, AVJK=6, UGPM=4, YGCX=7, CVDF=8, GKAG=8, SERM=6, AXNH=8, WTDE=8, RUUB=6, PBUU=6, HNDZ=7, NQYP=5, OWAS=7, EZGR=6, TWUW=6, TMAJ=7, DAYQ=8, CVEG=8, XNXF=8, WXNT=5, JEPE=7, LIRF=5, ROLB=8, XDDT=7, HANM=7, GMHZ=8, MFIN=7, ATFX=8, JCLQ=8, YOWV=8, BBBA=4, AEKF=4, ZHEE=5, NEBM=8, XWHK=5, MYMP=7, VQBJ=6, JPDJ=8, AXPJ=4, JIWJ=6, XHKT=8, UZXL=7, Manes_v4.1=8, ACFP=6, QRTH=5, SVVG=8, BFIK=7, SALZ=7, GKCZ=7, BYNZ=6, LGOW=7, OWFC=7, PZIF=5, XFJG=7, Sorbi_v2.1=6, JCMU=7, HAOX=6, QAUE=8, KJYC=5, AKXB=7, TQKZ=7, ISPU=5, AREG=5, DFDS=6, JVSZ=7, QNOC=6, BFJL=8, YKQR=6, VBHQ=4, NVGZ=8, BLVL=2, VDKG=8, CZPV=8, IUSR=7, BYQM=7, YOXI=8, BJSW=7, WEJN=5, QHBI=8, BYPY=8, LEMW=8, TDTF=8, JTQQ=7, OHKC=6, NTEO=6, JETM=8, JRNA=7, BSEY=8, SRSQ=2, QFAE=7, GXBM=8, FWCQ=8, BBDD=8, SJAN=8, LELS=8, WVMY=2, VQFW=7, ZNUM=8, EQYT=7, UCNM=8, LXRN=2, ZLQE=4, RKFX=7, KJZG=8, EEAQ=7, WXRI=2, RICC=6, QNPH=6, RDYY=6, FWBF=7, GDZS=7, PKMO=5, GSXD=7, GBVZ=8, VBLH=4, TOKV=8, YRBQ=7, PKOX=8, RMMV=8, AEPI=8, KOFB=6, Aquco_v1.1=7, ODDO=8, Solly_iTAGv2.3=8, MULF=6, HRUR=8, NXOH=6, YXNR=7, HATH=6, PMTB=8, BDJQ=7, DFHO=7, MUMD=7, RBYC=7, WZYK=7, ZJOJ=5, WTKZ=5, TFYI=7, Carpa_ASGPBv0.4=6, JTRM=7, UTQR=7, AZZW=7, QLMZ=6, WAFT=6, ZSAB=6, YZRI=7, GSZA=8, UTRE=6, SART=7, OUER=7, RQUG=8, ZFGK=7, Volca_v2.0=8, JGYZ=7, MDJK=7, EILE=7, WEQK=7, IKFD=8, XDLL=4, WAIL=8, VHZV=8, JEXA=6, RXEN=8, CKXF=7, PVGP=5, ZJRC=8, WXVX=5, XOAL=8, VUSY=8, MUNP=7, QWFV=6, BWRK=7, LPGY=7, POZS=7, GGEA=7, NKXU=7, UPMJ=8, RKLL=7, HLJG=8, VBMM=6, CCID=8, UEVI=7, ZHMB=8, PTBJ=7, Cyame_v1.0=6, UCRN=7, SJEV=8, WCLG=8, SNNC=8, ETCJ=7, DJSE=7, WGTU=8, MYVH=6, TKEK=7, HRWG=8, WCLV=6, Nelnu_v1.0=7, WCOR=8, OQBM=7, LETF=3, LVUS=8, AALA=6, FFGR=4, URSB=6, LNER=7, EGLZ=7, JNKW=7, RQZP=6, YKZB=8, GOWD=7, UPOG=8, FLTD=6, PTFA=8, VSRH=8, EAAA=7, NGRR=7, SNOX=7, DHPO=8, YADI=8, ZUHO=7, Pinta_v2.0=6, OLXF=6, GVCB=7, FFFY=8, WPHN=8, YZVJ=6, EEJO=6, DUMA=7, JCXF=5, OSHQ=6, MYZV=7, NGTD=8, MHYG=8, YCKE=7, KKDQ=7, BHYC=8, WNGH=7, HPXA=7, SYHW=5, BSNI=6, Phypa_v3.0=5, ZULJ=2, SLOI=4, BWVJ=5, WCQU=8, XMBA=5, WVWN=8, BOGT=5, PZRT=8, SHEZ=4, ZJUL=7, NXTS=7, IKIZ=7, ZSGF=8, ZDIZ=6, MQIV=8, HYHN=4, MFTM=8, EGNB=3, ETGN=2, ACRY=3, KFZY=7, VFYZ=7, IZGN=8, XBKS=8, KUXM=6, YRHD=8, ROWR=8, KBRW=8, NVSO=7, RMVB=6, VZCI=8, MWXT=4, AJAU=5, RVGH=7, HJMP=7, LNIL=5, LLEN=4, JHCN=7, UHBY=7, SYJM=6, VIAU=7, UJGI=8, LAPO=5, ROZZ=5, EEMJ=7, YIXP=7, YEPO=7, EGOS=8, TORX=8, SUAK=8, DUNJ=7, OSIP=6, WRPP=8, PCGJ=6, IXEM=7, XMCL=6, GGJD=6, RCAH=8, IOVS=7, YRMA=5, DHWX=6, GIPR=6, MFYC=4, DUQG=8, GEHT=7, FYSJ=6, Prupe_v1.0=7, OQHZ=7, TVCU=4, IZLO=7, UYFR=2, YPIC=7, WRSL=6, PIVW=8, ACWS=8, UYED=8, RCBT=8, LRRR=7, MWYQ=6, EVOD=8, ULKT=3, VIBO=8, RMWJ=1, KZED=5, BUWV=7, GIOY=1, PRFO=7, PIUF=8, SWGX=6, RGKI=7, OODC=5, CRNC=7, HFIK=5, FHST=8, OFUE=3, BOLZ=8, VKJD=8, FYUH=8, WLIC=4, NNGU=7, JLOV=2, CPKP=7, CTSS=4, PRIQ=5, EZZT=8, TOXE=8, KMNX=1, BSTR=7, GTHK=7, FUMQ=8, FFPD=5, LTZF=8, ATYL=7, RZTJ=6, VMNH=6, IZNU=7, JHFI=7, UAXP=8, LRTN=6, XSSD=8, DDRL=7, ERIA=5, PTLU=8, OSMU=8, CPLT=8, MFZO=8, FYTP=8, QFND=5, AJFN=7, YNFJ=6, KBXS=3, XMGP=7, BMIF=7, GCFE=7, OFTV=7, ENAU=3, JPYU=8, KVAY=7, CAPN=5, BZDF=7, OMDH=6, PCNH=8, IXLH=7, BMJR=8, ANON=4, XQRV=8, FQGQ=8, UWFU=8, TCBC=6, SSDU=8, JWGT=6, GRFT=8, WCZB=5, NPND=2, BSVG=8, SFKQ=8, IIOL=8, ARYD=8, JWEY=6, TZNS=8, NNHQ=7, JNUB=7, QWRA=5, APTP=6, Chlre_v5.5=7, FDMM=6, ZFXU=6, TXMP=8, PXYR=7, BXBF=5, BIDT=7, AAXJ=7, JSAG=6, NCVK=3, QYXY=5, LWCK=8, DFYF=8, BXAY=7, JYMN=6, SLYR=7, IRBN=8, ZBPY=7, GIWN=8, CAQZ=7, SWOH=3, CPOC=7, RTLC=7, QSKP=8, JNVS=6, XQSG=8, XOOE=1, OBPL=8, VQYB=8, IRAF=8, DWZT=3, KIIX=7, MMKU=8, Chlva_v1.0=5, FUPX=7, TGKW=5, QSLH=8, IEHF=6, AJJE=5, UHJR=5, UPZX=7, IMZV=8, XSZI=8, NHAG=7, HBGV=2, RNAT=4, RTMU=7, DMIN=7, DSUV=7, RTNA=8, HBHB=7, TGNL=5, NATT=1, NPRL=6, VGHH=8, KVFU=8, IKWM=5, MXDS=4, CTYH=7, OQON=5, SQCF=7, UFHF=6, LWDA=8, WAXR=8, ZMGN=8, PLBZ=7, ZSSR=8, QURC=7, LNSF=7, AYIY=7, UHLI=7, QSNJ=8, KTAR=6, WFBF=7, SWPE=8, FSQE=5, NHCM=4, NRWZ=7, ZQRI=7, RPGL=8, CJGZ=7, IADP=5, LLQV=7, MXFG=7, Ambtr_v1.0.27=8, GGWH=7, OBUY=6, KEGA=7, EJBY=7, NRXL=8, YLJA=8, WYIG=6, AUDE=7, OBTI=8, XQWC=7, RNBN=5, ECTD=7, FQLP=6, NAUM=7, VXKB=6, UFJN=8, GRKU=5, QUTB=8, NNOK=8, UUHD=7, KEFD=7, KGJF=7, DKFZ=8, ZBTA=8, WJLO=8, VMXJ=7, MZLD=8, AYMT=6, BMSE=7, TZWR=8, HWUP=7, CLMX=7, ZZEI=8, UUJS=8, PEZP=6, ZBVT=8, UWOD=7, HHXJ=4, GAMH=7, LHLE=7, XGFU=8, UQCB=7, Phavu_v1.0=7, LJPN=6, BILC=5, AWJM=7, ALUC=7, UDHA=5, CLNU=7, EJCM=7, DMLT=8, UJTT=6, ULXR=4, GETL=8, XMQO=8, TIUZ=8, HUQC=7, IXVJ=6, TPEM=7, PPQR=8, FOMH=5, YJJY=7, ZQVF=8, VKVG=2, HUSX=7, NJLF=8, IGUH=7, UJWU=6, DBYD=2, TNAW=6, TXVB=8, GAON=4, OOSO=7, QIAD=3, XKPS=4, XXHP=8, HYZL=7, SDPC=4, AYPS=5, WDCW=7, PPPZ=7, GANB=5, Arath_TAIR10=6, LJQF=4, QDVW=8, SBLT=7, GTSV=4, BKQU=8, XZME=8, TAGM=7, EATP=7, ALVQ=8, YHFG=8, LDEL=8, AUGV=7, BCAA=8, IEPQ=8, VXOD=7, HDSY=7, AUIP=8, IAJW=7, XMVD=7, YNUE=7, SOHV=6, Musac_v1.0=6, NHIX=5, GRRW=5, LLXJ=6, WDGV=7, JJZR=7, ABCD=6, Mimgu_v2.0=7, ENQF=8, PAWA=8, RAPY=4, BEFC=7, VITX=7, DOVJ=5, RTTY=7, DIHD=6, PYDB=6, MZOB=8, MOYY=7, ZQWM=8, OQWW=3, WWKL=4, ERWT=5, USIX=3, JQFK=5, TVSH=8, AFLV=7, CLRW=6, WHNV=8, CJNT=7, EAVM=7, LYPZ=8, BEGM=5, JHUL=7, ICNN=7, ADHK=7, YLPM=7, KADG=7, PHCE=8, HOKG=8, RRSV=4, HQOM=6, BTFM=5, VRGZ=4, OOVX=7, ALZF=6, VZWX=5, AWOI=8, ILBQ=8, YDCQ=1, TTRG=7, VTLJ=7, AJUW=5, RPPC=8, JBLI=4, GYBH=6, ABEH=6, Eucgr_v1.1=7, QIEH=5, YFGP=8, XVJB=8, KXSK=8, QICX=8, OXGJ=8, TPHT=6, FIDQ=4, BZSH=8, RCUX=7, DZLN=7, ZQYU=7, VCIN=7, SUVN=7, AHRN=8, XIRK=8, BCGB=6, PAZJ=7, RPQV=8, XRCX=7, VNAL=6, RHAU=7, JBND=7, ZKPF=3, KRJP=7, GCYL=5, MTHW=8, VVRN=8, PUAN=6, AWQB=8, MEKP=8, CBAE=7, YHLF=8, VVPY=7, UOEL=7, MTGC=7, HMHL=8, UBLN=8, JDQB=7, AQFM=7, HDWF=8, Ostlu_v2.0=7, LSHT=6, LFOG=8, TLCA=7, CDFR=8, Theca_v1.1=8, YSBD=5, EFCZ=7, MIRS=6, BEKN=8, AFPO=6, IANR=7, VGSX=8, PYHZ=7, SIBR=8, JKAA=6, LSKK=6, WSJO=6, WBIB=8, VGVI=6, JOIS=3, GNPX=7, XZUY=5, AQGE=8, XPBC=8, QKMG=7, NWMY=5, DZQM=8, FXGI=7, RAWF=7, ZINQ=8, TRRQ=6, YLWW=6, AFQQ=5, ABIJ=6, PUCW=7, GNQG=6, JDTY=7, HQRJ=8, WWQZ=8, XISJ=5, MTII=7, XKWQ=4, RPRU=7, ZXJO=6, DZPJ=5, FGDU=7, WOBD=8, LDME=5, EYKJ=8, BVOF=8, XPAF=6, ZZOL=6, HBUQ=7, LSJW=7, FZJL=8, DCCI=8, OINM=7, YYPE=8, MCHJ=4, EDBB=8, JUWL=6, ASMV=8, ZCDJ=8, UDUT=6, YJUG=7, XIVI=6, DCDT=8, VALZ=5, HXCD=8, OMYK=7, SMMC=8, DXOU=6, GYFU=8, BAHE=8, HIDG=4, WWSS=7, OVIJ=5, JOJQ=8, PUDI=8, FXHG=7, RFAD=4, GJIY=6, MKZR=7, RJIM=8, YWNF=6, KCPT=8, TEZA=7, DAAD=5, GNRI=7, FIKG=4, FOYQ=6, COAQ=8, TJES=5, PWKQ=7, IYDF=8, GUBD=5, FMVB=6, PNZO=7, DXQW=7, WOGB=7, SBZH=7, KJAA=6, LDRY=3, QKQO=7, OEKO=8, MVSE=6, YQEC=7, Poptr_v3.0=5, SXCE=6, BNCU=7, QIKZ=7, VVVV=7, EDEQ=8, DZTK=8, FEDW=7, LQJY=8, MVRF=8, MIXZ=8, EHNF=8, QOXT=7, TAVP=7, CFRN=7, DRFX=3, QMWB=7, BAKF=1, EDHN=8, HKMQ=7, COCP=5, EYRD=7, WBOD=8, FMWZ=7, GJNX=7, FZQN=7, WQML=8, BNDE=7, XRLM=8, BAJW=4, YUOM=7, JKHA=6, THDM=6, IRYH=5, BERS=6, KEYW=4, RJNQ=8, COBX=7, KTWL=6, YFQX=8, SIIK=8, UOMY=8, ZTHV=7, NHUA=8, MNCB=6, FCBJ=5, XVRU=7, BLAJ=5, TCYS=6, IHCQ=7, BGXB=5, WSPM=6, ZACW=7, XTON=3, RYJX=3, VTUS=8, WOHL=5, MRKX=7, CBJR=8, PJSX=8, FCCA=8, MGVU=7, HVBQ=7, THEW=4, FCEL=8, GHLP=8, OPDF=7, AKCR=7, EQDA=6, NMAK=4, YQIJ=6, ZIWB=4, IPWB=6, DRIL=7, JSZD=4, OAGK=8, LBRP=5, IAYV=2, PSHB=6, ZIVZ=5, MROH=4, KNMB=8, ZENX=8, DRGY=5, XXYA=8, OAEZ=4, KLGF=8, EFMS=7, ZGQD=3, EDIT=8, BGZG=7, INQX=6, Selmo_v1.0=8, SKNL=5, KYAD=8, IRZA=3, BPKH=7, QEHE=7, MCPK=5, GJPF=4, JGAB=8, YBML=4, Betvu_v1.1=7, RSCE=8, PQED=8, WDWX=6, ONBE=7, CBNG=2, TJLC=8, RWKR=6, VNMY=8, SKQD=8, BTTS=8, SMUR=6, Vitvi_Genoscope.12X=8, THHD=8, SCAO=7, AZBL=7, QACK=8, QCGM=7, ABSS=7, KRUQ=7, INSP=7, VJDZ=6, WKCY=8, Klefl_v1.0=7, DTNC=8, HXJE=6, ZTLR=8, HILW=8, SILJ=8, ZRIN=7, NBMW=8, GLVK=8, ORJE=5, VJED=4, CUTE=6, JKKI=6, DPFW=5, CMCY=7, NWWI=6, DESP=8, AXAF=3, Elagu_v2.0=6, HEGQ=6, RUIF=1, CMEQ=7, IHJY=6, DLAI=6, MAQO=7, GYRP=5, NFXV=5, KAWQ=7, QTJY=7, YSQT=7, PSKY=4, FGRF=8, FPCO=8, CQMG=5, VYER=5, CWYJ=7, KPUM=7, OTQG=4, DGXS=8, NSPR=5, VYDM=6, WQRD=8, YHZW=8, PSJT=7, FCHS=4, IJMT=7, QVMR=7, IDAU=8, ORKS=6, CWZU=7, NMDZ=7, CMFF=6, KPTE=8, SZPD=7, PWSG=8, JVBR=7, BRUD=6, VLNB=8, WQUF=8, QZXQ=8, UVDC=6, RFMZ=8, WMLW=7, WZFE=7, VYGG=5, YSRZ=8, SCEB=8, KAYP=2, YBQN=7, QPDY=6, XNLP=8, JMTE=5, CSSK=7, YFZK=7, EMBR=7, ZETY=8, AXBO=3, GDKK=7, DVXD=7, SXML=8, OXYP=8, JKNQ=8, JEBK=7, QXSZ=3, CQPW=7, YMES=8, HTFH=5, CKDK=8, XLGK=8, EMAL=8, NDUV=8, ZRMT=8, Spipo_v2=8, JGGD=5, IWIS=6, NMGG=8, NOKI=8, Micpu_v3.0=7, ASZK=5, QZZU=8, MNNM=4, WIGA=6, EWXK=5, UMUL=8, UOYN=8, CSUV=6, AIGO=7, LKKX=8, QXWF=7, BHBK=6, BCYF=3, CQQP=3, HTIP=7, FAJB=7, CIAC=6, XTZO=8, RQFE=6, KFEB=4, IDGE=7, IHOI=4, NUZN=8, PFUD=5, NSTT=5, Chabr_v0.1=8, GUMF=8, TSBQ=6, GBGT=5, XRTZ=4, HABV=3, PDQH=8, HZTS=8, JMUI=6, KYIO=5, TJQY=7, ZPKK=5, BPSG=8, VAZE=7, DVYE=7, RFRB=5, AQZD=8, RHVC=5, YKFU=7, HAEU=8, Orysa_v7.0=7, FNEN=6, VFFP=8, SKYV=6, WGET=8, FALI=8, WEAC=7, HENI=8, ZCUA=7, AOUJ=6, QCOU=7, KDCH=6, SZUO=6, OCTM=7, LGDQ=6, KWGC=7, HVNO=6, OLES=7, IHPC=8, XJGM=4, MWAN=5, FAKD=7, VDAO=3, UGJI=6, IFLI=7, IWMW=6, VHIJ=4, XAXW=5, ONLQ=8, UKUC=6, KYNE=7, NKIN=7, ROEI=8, ISHC=5, FROP=6, VJPU=7, SZYG=7, OCWZ=6, XAYK=6, OHAE=6, ZYAX=8, MHGD=8, CIEA=7, FVXD=5, GFUR=7, JMXW=5, VYLQ=8, HKZW=7, RFSD=7, DHAW=3, STKJ=6, FTRP=7, BJKT=7, RDOO=7, BAZF=5, ISGT=4, HGSM=6, FAMO=7, EDXZ=8, NIGS=8, WMUK=7, URDJ=8, HNCF=7, EMJJ=8, ZYCD=8, JTIG=5, CKKR=6, MFEA=2, WEEQ=6, XHHU=8, BNTL=7, JADL=6, JZVE=8, HERT=8, QNGJ=7, WVEF=7, YGAT=5, FRPM=7, UGNK=7, ISIM=6, URCP=6, NBYP=8, CGDN=7, SGTW=3, VFIV=6, XUAB=8, DLJZ=6, YZGX=4, RSOF=5}
Number of gene trees: 8
8 trees have missing taxa
Calculating quartet distance matrix (for completion of X)
Species tree distances calculated ...
Will attempt to complete bipartitions from X before adding using a distance matrix.
Building set of clusters (X) from gene trees 
------------------------------
gradient0: 14169
Number of Clusters after addition by distance: 14169
calculating extra bipartitions to be added at level 1 ...
Adding to X using resolutions of greedy consensus ...
Limit for sigma of degrees:29450
polytomy size limit : 94
discarded polytomies:  [3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 9, 9, 9, 10, 10, 10, 11, 11, 11, 11, 11, 12, 13, 13, 13, 13, 13, 13, 13, 13, 17, 17, 17, 17, 17, 19, 20, 21, 24, 49, 85, 94]
Threshold 0.0:
polytomy of size 3; rounds with additions with at least 5 support: 0; clusters: 14175
polytomy of size 3; rounds with additions with at least 5 support: 0; clusters: 14181
polytomy of size 13; rounds with additions with at least 5 support: 0; clusters: 15519
polytomy of size 8; rounds with additions with at least 5 support: 0; clusters: 15815
polytomy of size 4; rounds with additions with at least 5 support: 0; clusters: 15833
polytomy of size 3; rounds with additions with at least 5 support: 1; clusters: 15839
polytomy of size 5; rounds with additions with at least 5 support: 0; clusters: 15883
polytomy of size 3; rounds with additions with at least 5 support: 0; clusters: 15887
polytomy of size 3; rounds with additions with at least 5 support: 0; clusters: 15891
polytomy of size 3; rounds with additions with at least 5 support: 1; clusters: 15897
polytomy of size 3; rounds with additions with at least 5 support: 0; clusters: 15903
polytomy of size 5; rounds with additions with at least 5 support: 0; clusters: 15953
polytomy of size 3; rounds with additions with at least 5 support: 0; clusters: 15957
polytomy of size 3; rounds with additions with at least 5 support: 0; clusters: 15963
polytomy of size 8; rounds with additions with at least 5 support: 0; clusters: 16237
polytomy of size 17; rounds with additions with at least 5 support: 1; clusters: 18697
polytomy of size 4; rounds with additions with at least 5 support: 0; clusters: 18717
polytomy of size 3; rounds with additions with at least 5 support: 0; clusters: 18723
polytomy of size 7; rounds with additions with at least 5 support: 0; clusters: 18895
polytomy of size 3; rounds with additions with at least 5 support: 1; clusters: 18901
polytomy of size 3; rounds with additions with at least 5 support: 0; clusters: 18907
polytomy of size 3; rounds with additions with at least 5 support: 0; clusters: 18913
polytomy of size 4; rounds with additions with at least 5 support: 0; clusters: 18933
polytomy of size 3; rounds with additions with at least 5 support: 0; clusters: 18937
polytomy of size 4; rounds with additions with at least 5 support: 0; clusters: 18957
Threshold 0.01:
polytomy of size 3; rounds with additions with at least 5 support: 0; clusters: 18957
polytomy of size 3; rounds with additions with at least 5 support: 0; clusters: 18957
polytomy of size 13; rounds with additions with at least 5 support: 0; clusters: 19931
polytomy of size 8; rounds with additions with at least 5 support: 0; clusters: 20017
polytomy of size 4; rounds with additions with at least 5 support: 0; clusters: 20019
polytomy of size 3; rounds with additions with at least 5 support: 0; clusters: 20019
polytomy of size 5; rounds with additions with at least 5 support: 0; clusters: 20023
polytomy of size 3; rounds with additions with at least 5 support: 0; clusters: 20023
polytomy of size 3; rounds with additions with at least 5 support: 0; clusters: 20023
polytomy of size 3; rounds with additions with at least 5 support: 0; clusters: 20023
polytomy of size 3; rounds with additions with at least 5 support: 0; clusters: 20023
polytomy of size 5; rounds with additions with at least 5 support: 0; clusters: 20023
polytomy of size 3; rounds with additions with at least 5 support: 0; clusters: 20023
polytomy of size 3; rounds with additions with at least 5 support: 0; clusters: 20023
polytomy of size 8; rounds with additions with at least 5 support: 0; clusters: 20075
polytomy of size 17; rounds with additions with at least 5 support: 0; clusters: 21507
polytomy of size 4; rounds with additions with at least 5 support: 0; clusters: 21507
polytomy of size 3; rounds with additions with at least 5 support: 0; clusters: 21507
polytomy of size 7; rounds with additions with at least 5 support: 0; clusters: 21537
polytomy of size 3; rounds with additions with at least 5 support: 0; clusters: 21537
polytomy of size 3; rounds with additions with at least 5 support: 0; clusters: 21537
polytomy of size 3; rounds with additions with at least 5 support: 0; clusters: 21537
polytomy of size 4; rounds with additions with at least 5 support: 0; clusters: 21537
polytomy of size 3; rounds with additions with at least 5 support: 0; clusters: 21537
polytomy of size 4; rounds with additions with at least 5 support: 0; clusters: 21537
Threshold 0.02:
polytomy of size 3; rounds with additions with at least 5 support: 0; clusters: 21537
polytomy of size 3; rounds with additions with at least 5 support: 0; clusters: 21537
polytomy of size 13; rounds with additions with at least 5 support: 0; clusters: 22297
polytomy of size 8; rounds with additions with at least 5 support: 0; clusters: 22355
polytomy of size 4; rounds with additions with at least 5 support: 0; clusters: 22355
polytomy of size 3; rounds with additions with at least 5 support: 0; clusters: 22355
polytomy of size 5; rounds with additions with at least 5 support: 0; clusters: 22355
polytomy of size 3; rounds with additions with at least 5 support: 0; clusters: 22355
polytomy of size 3; rounds with additions with at least 5 support: 0; clusters: 22357
polytomy of size 3; rounds with additions with at least 5 support: 0; clusters: 22357
polytomy of size 3; rounds with additions with at least 5 support: 0; clusters: 22357
polytomy of size 5; rounds with additions with at least 5 support: 0; clusters: 22357
polytomy of size 3; rounds with additions with at least 5 support: 0; clusters: 22357
polytomy of size 3; rounds with additions with at least 5 support: 0; clusters: 22357
polytomy of size 8; rounds with additions with at least 5 support: 0; clusters: 22397
polytomy of size 17; rounds with additions with at least 5 support: 0; clusters: 23915
polytomy of size 4; rounds with additions with at least 5 support: 0; clusters: 23915
polytomy of size 3; rounds with additions with at least 5 support: 0; clusters: 23915
polytomy of size 7; rounds with additions with at least 5 support: 0; clusters: 23931
polytomy of size 3; rounds with additions with at least 5 support: 0; clusters: 23931
polytomy of size 3; rounds with additions with at least 5 support: 0; clusters: 23931
polytomy of size 3; rounds with additions with at least 5 support: 0; clusters: 23931
polytomy of size 4; rounds with additions with at least 5 support: 0; clusters: 23931
polytomy of size 3; rounds with additions with at least 5 support: 0; clusters: 23931
polytomy of size 4; rounds with additions with at least 5 support: 0; clusters: 23931
Threshold 0.05:
polytomy of size 3; rounds with additions with at least 5 support: 0; clusters: 23931
polytomy of size 3; rounds with additions with at least 5 support: 0; clusters: 23931
polytomy of size 13; rounds with additions with at least 5 support: 0; clusters: 23933
polytomy of size 8; rounds with additions with at least 5 support: 0; clusters: 23933
polytomy of size 4; rounds with additions with at least 5 support: 0; clusters: 23933
polytomy of size 3; rounds with additions with at least 5 support: 0; clusters: 23933
polytomy of size 5; rounds with additions with at least 5 support: 0; clusters: 23933
polytomy of size 3; rounds with additions with at least 5 support: 0; clusters: 23933
polytomy of size 3; rounds with additions with at least 5 support: 0; clusters: 23933
polytomy of size 3; rounds with additions with at least 5 support: 0; clusters: 23933
polytomy of size 3; rounds with additions with at least 5 support: 0; clusters: 23933
polytomy of size 5; rounds with additions with at least 5 support: 0; clusters: 23933
polytomy of size 3; rounds with additions with at least 5 support: 0; clusters: 23933
polytomy of size 3; rounds with additions with at least 5 support: 0; clusters: 23933
polytomy of size 8; rounds with additions with at least 5 support: 0; clusters: 23933
polytomy of size 17; rounds with additions with at least 5 support: 0; clusters: 23947
polytomy of size 4; rounds with additions with at least 5 support: 0; clusters: 23947
polytomy of size 3; rounds with additions with at least 5 support: 0; clusters: 23947
polytomy of size 7; rounds with additions with at least 5 support: 0; clusters: 23947
polytomy of size 3; rounds with additions with at least 5 support: 0; clusters: 23947
polytomy of size 3; rounds with additions with at least 5 support: 0; clusters: 23947
polytomy of size 3; rounds with additions with at least 5 support: 0; clusters: 23947
polytomy of size 4; rounds with additions with at least 5 support: 0; clusters: 23947
polytomy of size 3; rounds with additions with at least 5 support: 0; clusters: 23947
polytomy of size 4; rounds with additions with at least 5 support: 0; clusters: 23947
Threshold 0.1:
polytomy of size 3; rounds with additions with at least 5 support: 0; clusters: 23947
polytomy of size 3; rounds with additions with at least 5 support: 0; clusters: 23947
polytomy of size 13; rounds with additions with at least 5 support: 0; clusters: 23951
polytomy of size 8; rounds with additions with at least 5 support: 0; clusters: 23951
polytomy of size 4; rounds with additions with at least 5 support: 0; clusters: 23951
polytomy of size 3; rounds with additions with at least 5 support: 0; clusters: 23951
polytomy of size 5; rounds with additions with at least 5 support: 0; clusters: 23951
polytomy of size 3; rounds with additions with at least 5 support: 0; clusters: 23951
polytomy of size 3; rounds with additions with at least 5 support: 0; clusters: 23951
polytomy of size 3; rounds with additions with at least 5 support: 0; clusters: 23951
polytomy of size 3; rounds with additions with at least 5 support: 0; clusters: 23951
polytomy of size 5; rounds with additions with at least 5 support: 0; clusters: 23951
polytomy of size 3; rounds with additions with at least 5 support: 0; clusters: 23951
polytomy of size 3; rounds with additions with at least 5 support: 0; clusters: 23951
polytomy of size 8; rounds with additions with at least 5 support: 0; clusters: 23955
polytomy of size 17; rounds with additions with at least 5 support: 0; clusters: 23957
polytomy of size 4; rounds with additions with at least 5 support: 0; clusters: 23957
polytomy of size 3; rounds with additions with at least 5 support: 0; clusters: 23957
polytomy of size 7; rounds with additions with at least 5 support: 0; clusters: 23957
polytomy of size 3; rounds with additions with at least 5 support: 0; clusters: 23957
polytomy of size 3; rounds with additions with at least 5 support: 0; clusters: 23957
polytomy of size 3; rounds with additions with at least 5 support: 0; clusters: 23957
polytomy of size 4; rounds with additions with at least 5 support: 0; clusters: 23957
polytomy of size 3; rounds with additions with at least 5 support: 0; clusters: 23957
polytomy of size 4; rounds with additions with at least 5 support: 0; clusters: 23957
Threshold 0.2:
polytomy of size 9; rounds with additions with at least 5 support: 0; clusters: 23983
polytomy of size 4; rounds with additions with at least 5 support: 0; clusters: 23983
polytomy of size 4; rounds with additions with at least 5 support: 0; clusters: 23985
polytomy of size 3; rounds with additions with at least 5 support: 0; clusters: 23987
polytomy of size 4; rounds with additions with at least 5 support: 0; clusters: 23995
polytomy of size 3; rounds with additions with at least 5 support: 0; clusters: 23999
polytomy of size 4; rounds with additions with at least 5 support: 0; clusters: 24001
polytomy of size 3; rounds with additions with at least 5 support: 0; clusters: 24001
polytomy of size 3; rounds with additions with at least 5 support: 0; clusters: 24001
polytomy of size 3; rounds with additions with at least 5 support: 0; clusters: 24003
polytomy of size 4; rounds with additions with at least 5 support: 0; clusters: 24009
polytomy of size 9; rounds with additions with at least 5 support: 0; clusters: 24041
polytomy of size 3; rounds with additions with at least 5 support: 0; clusters: 24043
polytomy of size 3; rounds with additions with at least 5 support: 0; clusters: 24043
polytomy of size 6; rounds with additions with at least 5 support: 0; clusters: 24069
polytomy of size 3; rounds with additions with at least 5 support: 0; clusters: 24069
polytomy of size 3; rounds with additions with at least 5 support: 0; clusters: 24073
polytomy of size 5; rounds with additions with at least 5 support: 0; clusters: 24073
polytomy of size 3; rounds with additions with at least 5 support: 0; clusters: 24075
polytomy of size 4; rounds with additions with at least 5 support: 0; clusters: 24075
polytomy of size 5; rounds with additions with at least 5 support: 0; clusters: 24081
polytomy of size 11; rounds with additions with at least 5 support: 0; clusters: 24151
polytomy of size 3; rounds with additions with at least 5 support: 0; clusters: 24153
polytomy of size 4; rounds with additions with at least 5 support: 0; clusters: 24159
polytomy of size 3; rounds with additions with at least 5 support: 0; clusters: 24161
polytomy of size 3; rounds with additions with at least 5 support: 0; clusters: 24161
polytomy of size 3; rounds with additions with at least 5 support: 0; clusters: 24161
polytomy of size 11; rounds with additions with at least 5 support: 0; clusters: 24207
polytomy of size 3; rounds with additions with at least 5 support: 0; clusters: 24209
polytomy of size 3; rounds with additions with at least 5 support: 0; clusters: 24211
polytomy of size 3; rounds with additions with at least 5 support: 0; clusters: 24211
polytomy of size 3; rounds with additions with at least 5 support: 0; clusters: 24211
polytomy of size 3; rounds with additions with at least 5 support: 0; clusters: 24215
polytomy of size 3; rounds with additions with at least 5 support: 0; clusters: 24217
polytomy of size 13; rounds with additions with at least 5 support: 0; clusters: 24281
polytomy of size 3; rounds with additions with at least 5 support: 0; clusters: 24283
polytomy of size 3; rounds with additions with at least 5 support: 0; clusters: 24283
polytomy of size 3; rounds with additions with at least 5 support: 0; clusters: 24283
polytomy of size 6; rounds with additions with at least 5 support: 0; clusters: 24303
polytomy of size 3; rounds with additions with at least 5 support: 0; clusters: 24305
polytomy of size 4; rounds with additions with at least 5 support: 0; clusters: 24311
polytomy of size 5; rounds with additions with at least 5 support: 0; clusters: 24319
polytomy of size 3; rounds with additions with at least 5 support: 0; clusters: 24321
polytomy of size 3; rounds with additions with at least 5 support: 0; clusters: 24321
polytomy of size 12; rounds with additions with at least 5 support: 0; clusters: 24405
polytomy of size 24; rounds with additions with at least 5 support: 1; clusters: 24705
polytomy of size 4; rounds with additions with at least 5 support: 0; clusters: 24705
polytomy of size 3; rounds with additions with at least 5 support: 0; clusters: 24705
polytomy of size 20; rounds with additions with at least 5 support: 0; clusters: 24837
polytomy of size 5; rounds with additions with at least 5 support: 0; clusters: 24847
polytomy of size 9; rounds with additions with at least 5 support: 0; clusters: 24913
polytomy of size 4; rounds with additions with at least 5 support: 0; clusters: 24919
polytomy of size 4; rounds with additions with at least 5 support: 0; clusters: 24923
polytomy of size 3; rounds with additions with at least 5 support: 0; clusters: 24925
polytomy of size 4; rounds with additions with at least 5 support: 0; clusters: 24925
polytomy of size 13; rounds with additions with at least 5 support: 0; clusters: 24959
polytomy of size 4; rounds with additions with at least 5 support: 0; clusters: 24967
polytomy of size 5; rounds with additions with at least 5 support: 0; clusters: 24979
polytomy of size 5; rounds with additions with at least 5 support: 0; clusters: 24987
polytomy of size 49; rounds with additions with at least 5 support: 0; clusters: 25677
polytomy of size 3; rounds with additions with at least 5 support: 0; clusters: 25679
polytomy of size 3; rounds with additions with at least 5 support: 0; clusters: 25681
Threshold 0.3333333333333333:
polytomy of size 4; rounds with additions with at least 5 support: 0; clusters: 25685
polytomy of size 7; rounds with additions with at least 5 support: 0; clusters: 25691
polytomy of size 8; rounds with additions with at least 5 support: 0; clusters: 25705
polytomy of size 3; rounds with additions with at least 5 support: 0; clusters: 25705
polytomy of size 3; rounds with additions with at least 5 support: 0; clusters: 25705
polytomy of size 11; rounds with additions with at least 5 support: 0; clusters: 25713
polytomy of size 4; rounds with additions with at least 5 support: 0; clusters: 25717
polytomy of size 6; rounds with additions with at least 5 support: 0; clusters: 25717
polytomy of size 5; rounds with additions with at least 5 support: 0; clusters: 25719
polytomy of size 3; rounds with additions with at least 5 support: 0; clusters: 25719
polytomy of size 3; rounds with additions with at least 5 support: 0; clusters: 25719
polytomy of size 3; rounds with additions with at least 5 support: 0; clusters: 25721
polytomy of size 3; rounds with additions with at least 5 support: 0; clusters: 25721
polytomy of size 3; rounds with additions with at least 5 support: 0; clusters: 25721
polytomy of size 7; rounds with additions with at least 5 support: 0; clusters: 25723
polytomy of size 10; rounds with additions with at least 5 support: 0; clusters: 25731
polytomy of size 6; rounds with additions with at least 5 support: 0; clusters: 25737
polytomy of size 4; rounds with additions with at least 5 support: 0; clusters: 25739
polytomy of size 5; rounds with additions with at least 5 support: 0; clusters: 25741
polytomy of size 3; rounds with additions with at least 5 support: 0; clusters: 25741
polytomy of size 4; rounds with additions with at least 5 support: 0; clusters: 25743
polytomy of size 11; rounds with additions with at least 5 support: 0; clusters: 25779
polytomy of size 4; rounds with additions with at least 5 support: 0; clusters: 25779
polytomy of size 4; rounds with additions with at least 5 support: 0; clusters: 25783
polytomy of size 3; rounds with additions with at least 5 support: 0; clusters: 25783
polytomy of size 4; rounds with additions with at least 5 support: 0; clusters: 25787
polytomy of size 5; rounds with additions with at least 5 support: 0; clusters: 25787
polytomy of size 3; rounds with additions with at least 5 support: 0; clusters: 25787
polytomy of size 10; rounds with additions with at least 5 support: 0; clusters: 25821
polytomy of size 6; rounds with additions with at least 5 support: 0; clusters: 25821
polytomy of size 4; rounds with additions with at least 5 support: 0; clusters: 25823
polytomy of size 11; rounds with additions with at least 5 support: 0; clusters: 25843
polytomy of size 8; rounds with additions with at least 5 support: 0; clusters: 25849
polytomy of size 3; rounds with additions with at least 5 support: 0; clusters: 25849
polytomy of size 6; rounds with additions with at least 5 support: 0; clusters: 25851
polytomy of size 4; rounds with additions with at least 5 support: 0; clusters: 25851
polytomy of size 3; rounds with additions with at least 5 support: 0; clusters: 25851
polytomy of size 4; rounds with additions with at least 5 support: 0; clusters: 25851
polytomy of size 3; rounds with additions with at least 5 support: 0; clusters: 25851
polytomy of size 3; rounds with additions with at least 5 support: 0; clusters: 25851
polytomy of size 4; rounds with additions with at least 5 support: 0; clusters: 25855
polytomy of size 6; rounds with additions with at least 5 support: 0; clusters: 25859
polytomy of size 3; rounds with additions with at least 5 support: 0; clusters: 25859
polytomy of size 19; rounds with additions with at least 5 support: 0; clusters: 25941
polytomy of size 4; rounds with additions with at least 5 support: 0; clusters: 25941
polytomy of size 5; rounds with additions with at least 5 support: 0; clusters: 25947
polytomy of size 13; rounds with additions with at least 5 support: 0; clusters: 25975
polytomy of size 6; rounds with additions with at least 5 support: 0; clusters: 25981
polytomy of size 3; rounds with additions with at least 5 support: 0; clusters: 25981
polytomy of size 6; rounds with additions with at least 5 support: 0; clusters: 25985
polytomy of size 3; rounds with additions with at least 5 support: 0; clusters: 25985
polytomy of size 7; rounds with additions with at least 5 support: 0; clusters: 25989
polytomy of size 3; rounds with additions with at least 5 support: 0; clusters: 25989
polytomy of size 94; rounds with additions with at least 5 support: 0; clusters: 26901
polytomy of size 3; rounds with additions with at least 5 support: 0; clusters: 26901
polytomy of size 4; rounds with additions with at least 5 support: 0; clusters: 26901
polytomy of size 4; rounds with additions with at least 5 support: 0; clusters: 26909
polytomy of size 10; rounds with additions with at least 5 support: 0; clusters: 26937
polytomy of size 5; rounds with additions with at least 5 support: 0; clusters: 26941
polytomy of size 5; rounds with additions with at least 5 support: 0; clusters: 26945
polytomy of size 3; rounds with additions with at least 5 support: 0; clusters: 26945
polytomy of size 3; rounds with additions with at least 5 support: 0; clusters: 26945
polytomy of size 4; rounds with additions with at least 5 support: 0; clusters: 26947
polytomy of size 6; rounds with additions with at least 5 support: 0; clusters: 26947
polytomy of size 5; rounds with additions with at least 5 support: 0; clusters: 26953
polytomy of size 21; rounds with additions with at least 5 support: 0; clusters: 26977
polytomy of size 3; rounds with additions with at least 5 support: 0; clusters: 26977
polytomy of size 85; rounds with additions with at least 5 support: 0; clusters: 27703
polytomy of size 5; rounds with additions with at least 5 support: 0; clusters: 27715
max k is :2
Number of Clusters after addition by greedy: 27715
gradient0 in heuristiic: 27715
partitions formed in 104.011 secs
Dynamic Programming starting after 104.011 secs
Using tree-based weight calculation.
Using polytree-based weight calculation.
Polytree max score: 317135137206
Polytree building time: 0.118 seconds.
Number of quartet trees in the gene trees: 317135137206
Size of largest cluster: 1176
Greedy score: 251749701553
estimationFactor: 1.2597239847739585
Sub-optimal score: 262732553206
Total Number of elements weighted: 88296
Normalized score (portion of input quartet trees satisfied before correcting for multiple individuals): 0.8411102114766025
Optimization score: 266745602322
Optimal tree inferred in 140.584 secs.
(DBYD,(BOGT,(VKVG,((BAJW,(EBWI,(NMAK,RFAD))),((ISIM,AYPS),(((Cyame_v1.0,((RTLC,RSOF),(LLXJ,(JJZR,(PVGP,OBUY))))),(ZULJ,(VZWX,((XAXW,YSBD),((PYDB,(URSB,(IKIZ,ZJOJ))),((IKWM,((LJPN,VNAL),(FTRP,PWKQ))),((BWVJ,(JEBK,IHJY)),(IEHF,(CKXF,(SBLT,(WEJN,UGPM))))))))))),(((ROZZ,(BAKF,(IRZA,IAYV))),((FSQE,APTP),(LIRF,((ASZK,(QLMZ,(JCXF,(VRGZ,ULXR)))),((FIDQ,(RAPY,SRSQ)),(FIKG,(HFIK,(FOMH,(LDRY,(VYER,(JGGD,(YRMA,RWXW)))))))))))),((XOAL,(VBLH,(XJGM,HYHN))),(((XMCL,(XIVI,TNAW)),((JMTE,(MCPK,(Ostlu_v2.0,(Micpu_v3.0,QXSZ)))),((MMKU,EGNB),((RRSV,PUAN),((((GYBH,JTIG),(YDCQ,NATT)),(MNCB,((QYXY,ETGN),((GXBM,ZNUM),(MFYC,(WCLV,WXRI)))))),((NKXU,KADG),((EATP,((BILC,JQFK),(AKCR,(RUIF,(EEJO,Chlva_v1.0))))),((((ALZF,NSTT),(FMVB,(DUMA,(HVNO,HHXJ)))),((JIWJ,(LETF,(AJAU,CBNG))),(GJIY,((OQON,(CQQP,OAEZ)),(NQYP,(LSHT,KSFK)))))),(IRYH,((PFUD,(DVYE,SDPC)),((DZPJ,(FOYQ,((TGNL,((MWAN,PZIF),(VJDZ,(XKWQ,XTON)))),(((VFIV,JMUI),(IJMT,(RAWF,ISGT))),(VIAU,(BAZF,(ENAU,(OTQG,MXDS)))))))),((FXHG,(BZSH,(UTRE,(AOUJ,(AJUW,LBRP))))),(((QWRA,USIX),(JKKI,(BCYF,((KUJU,IHOI),(RYJX,((Volca_v2.0,(WRSL,JWGT)),(Chlre_v5.5,(ISPU,RNAT)))))))),((XDLL,OFUE),(((ZIVZ,(LNIL,ACRY)),((TSBQ,GFUR),(JRGZ,(MULF,VALZ)))),(((PRIQ,QRTH),(KFEB,GUBD)),(ZLQE,(ZFXU,(SYJM,(UKUC,(RHVC,(WDWX,ZDIZ)))))))))))))))))))))),(((VHIJ,KYIO),(JKHA,(QFND,(PQED,POOW)))),((LXRN,(VJED,LLEN)),((OQWW,(UYFR,BTFM)),((NNHQ,AZZW),(((DRGY,(QPDY,VQBJ)),(TPHT,(BFIK,(FPCO,(Klefl_v1.0,FQLP))))),((KMNX,(KEYW,(DRFX,(AEKF,(RPQV,(YSQT,(HIDG,((WDGV,(MOYY,(MCHJ,(STKJ,GBGT)))),(GYRP,((BHBK,(RPRU,WCQU)),(RQFE,((MNNM,ISHC),(QWFV,(GGWH,DFDS)))))))))))))),(((HAOX,(SNOX,XRTZ)),((FFGR,(WDCW,WSJO)),(ZRMT,((YOXI,(RPGL,VAZE)),(HKZW,(JOJQ,(MFZO,NBYP))))))),((Chabr_v0.1,MWXT),(((((PKOX,PYHZ),(KUXM,((ZFGK,(JKAA,ZYCD)),(KJYC,((ABIJ,LGDQ),(ZZOL,Selmo_v1.0)))))),((UPMJ,ULKT),((WAFT,(ENQF,(PQTO,XNXF))),(ZZEI,((CBAE,GAON),(GKAG,YHZW)))))),(((HERT,((WJLO,ILBQ),((TFYI,TXVB),(HMHL,JPYU)))),(YFGP,((HPXA,(JHFI,PIUF)),(BNCU,(NRWZ,(TGKW,((KRUQ,UUHD),(LGOW,((RTMU,WZYK),(YBQN,(OFTV,IRBN))))))))))),(((BSNI,TWUW),(ANON,((FAJB,(WEEQ,WCZB)),(AKXB,(UCRN,(DXOU,TCBC)))))),(SKQD,((RCBT,(GOWD,UHLI)),(WOGB,((HVBQ,(SZYG,ZTHV)),(HRWG,(AWOI,(((Phypa_v3.0,YEPO),(KEFD,(BPSG,(FFPD,(GRKU,((ABCD,RDOO),(NGTD,(VMXJ,RGKI)))))))),(BGXB,(ZQRI,(CMEQ,((YWNF,(XWHK,JMXW)),(ORKS,(WNGH,((LNSF,MIRS),(TAVP,((ZACW,IGUH),((DHWX,QKQO),((EEMJ,(TMAJ,WSPM)),(QMWB,(VBMM,JADL))))))))))))))))))))))),(((JVSZ,CAPN),(((QVMR,ALVQ),(BEGM,(EEAQ,DJSE))),((DFHO,(UGNK,NHCM)),((RFMZ,(VIBO,UOMY)),(MEKP,((PBUU,CQPW),((GANB,(EWXK,(PNZO,UWOD))),((KIIX,(CVEG,QIAD)),(VVRN,((NOKI,YIXP),((MTGC,(PIVW,((UJTT,(POPJ,FLTD)),(WQML,(((NDUV,SKYV),(WCLG,BMJR)),((XDDT,GSXD),(ZXJO,(DCDT,YCKE)))))))),((WGTU,((JBLI,RFRB),(FQGQ,(NWWI,((ZQYU,ORJE),(UJWU,((CJNT,YLJA),(GYFU,IXLH)))))))),((KJZG,(BMIF,PSKY)),((HEGQ,(HNDZ,((YOWV,XXHP),(LHLE,RICC)))),(YQEC,(MROH,((OCZL,(HTFH,YJJY)),((VITX,FCHS),(UFJN,(AFPO,URCP)))))))))))))))))))),(((SGTW,(XZUY,(WLIC,(GNQG,KAWQ)))),(((VDAO,(TOXE,GTHK)),((GGEA,((GAMH,AREG),(JUWL,(VSRH,AQFM)))),(NPRL,((IOVS,WVWN),(AWQB,(IIOL,((DZQM,Pinta_v2.0),(MFTM,JBND)))))))),(((XTZO,(RSCE,(ACWS,MIXZ))),((BBDD,(JRNA,((OWFC,(KLGF,EGLZ)),(JZVE,(CDFR,ZQWM))))),(QCGM,(MHGD,(HILW,(FMWZ,((IZGN,(ROWR,QHBI)),((XLGK,SCEB),(VGSX,UUJS))))))))),(YFZK,((NVGZ,((IAJW,(EFMS,HQOM)),((ZYAX,WWSS),(YLPM,BTTS)))),(ZQVF,(QSNJ,(XIRK,((QFAE,(NRXL,HBGV)),((GMHZ,(FHST,OXGJ)),((OVIJ,(YYPE,(AUDE,(ETCJ,(GKCZ,(JDQB,(RMMV,IFLI))))))),((VFYZ,NKIN),((UEVI,AIGO),((QNGJ,XMGP),(CGDN,(FRPM,(BUWV,XQSG))))))))))))))))),(((((Ambtr_v1.0.27,URDJ),(PZRT,WTKZ)),((QDVW,PAWA),((XSZI,MUNP),(CSSK,OPDF)))),((((WZFE,OSHQ),(SJEV,WKSU)),(((DWZT,(FZJL,(VZCI,NWMY))),(OBPL,((WBOD,XQWC),(DHPO,(YZRI,PSJT))))),((FALI,WPHN),(VYLQ,((MAQO,WAIL),(BSVG,(KRJP,((WIGA,KAYP),(ABSS,BCGB))))))))),(MTII,(((COCP,BYQM),(Spipo_v2,(YMES,MFIN))),((OCWZ,((SILJ,QOXT),(VBHQ,(VVVV,DGXS)))),((((Elagu_v2.0,(NSPR,(HXJE,HWUP))),((XHHU,LSKK),(UOEL,((JNUB,TZNS),(Musac_v1.0,(LEMW,BDJQ)))))),(BYPY,((PPQR,BRUD),((CIEA,(XPAF,(WBIB,PWSG))),(BSTR,(WXNT,(Orysa_v7.0,((RMVB,IADP),(((YXNR,EFCZ),(RCAH,YLWW)),(HATH,((YPIC,(ROEI,Sorbi_v2.1)),(SOHV,((WCOR,XBKS),((ZMGN,(BPKH,XUAB)),((NNOK,BXAY),(ZENX,VQYB)))))))))))))))),((NPND,((AFLV,OOSO),(MWYQ,(GDKK,(QNPH,NHIX))))),((THDM,(JSAG,((MTHW,XZME),(LELS,VTUS)))),(THEW,((YJUG,EMJJ),(SWOH,((KYNE,RDYY),((SART,LTZF),((JVBR,(FCEL,(WTDE,(BLAJ,JHUL)))),((PRFO,((GJPF,KBXS),(TRRQ,((DMIN,ZKPF),(JDTY,(LDME,DPFW)))))),(IXEM,(FGRF,((MUMD,(KXSK,(CMCY,(PLBZ,(ICNN,YBML))))),((SVTS,KOFB),(ONBE,((MVRF,LSJW),(RCUX,((UZXL,TCYS),(XFJG,(HOKG,RQZP))))))))))))))))))))))))),(((IWMW,(AALA,((FAKD,Nelnu_v1.0),(VQFW,(GRRW,(SIIK,RQUG)))))),(QTJY,((CCID,((WFBF,YHFG),(VGHH,((GBVZ,Aquco_v1.1),(ZUHO,UPOG))))),((NMGG,(AUGV,(UDHA,ZGQD))),(EVOD,((XMVD,XHKT),(IRAF,(SSDU,(BEKN,(RQNK,QCOU)))))))))),(EHNF,((XMQO,(((BGZG,(Vitvi_Genoscope.12X,(SZPD,BBBA))),((VGVI,XKPS),((XGFU,RSPO),(QKMG,BSEY)))),((((IUSR,TOKV),(KWGC,(ZJUL,DRIL))),((HTIP,(FYTP,(NUZN,(YKQR,(HQRJ,OMDH))))),(UWFU,(FCBJ,(GIOY,(SYHW,((ERIA,(DAYQ,(CTSS,SLOI))),((UHBY,CKKR),(CIAC,(YKFU,OOVX)))))))))),((ZHMB,(HDWF,(AJJE,(((WWQZ,SWGX),((FGDU,NEBM),(Eucgr_v1.1,AYMT))),((YNUE,RJNQ),(FEDW,(((UJGI,TLCA),(JKNQ,(AXNH,PMTB))),((EQYT,SJAN),((ZINQ,((ARYD,IMZV),(ROLB,(IDAU,BZDF)))),((KBRW,YHLF),(GVCB,(HKMQ,(DZLN,CJGZ))))))))))))),((PTLU,VYGG),(((KVAY,UDUT),((((RBYC,(ZHEE,(WVEF,EILE))),(CQMG,((KYAD,(BJSW,AQGE)),((EDHN,XVJB),(WKCY,ACFP))))),(ZPKK,((DHAW,(IANR,QNOC)),((SQCF,XFFT),((TJQY,SXCE),((Prupe_v1.0,NCVK),(EAVM,(QIEH,(VCIN,BLVL))))))))),((TJLC,((SVVG,(HENI,(NHUA,UZWG))),((INSP,DXQW),(LNER,(LWDA,CWZU))))),(OHAE,(OQHZ,(((QZXQ,(VHZV,GEHT)),(KZED,(XOOE,(ZCDJ,PJSX)))),((RKLL,(RKFX,JETM)),((SLYR,ZSSR),((CMFF,TTRG),(VLNB,((NXOH,((KEGA,TVSH),(Phavu_v1.0,SUAK))),((KNMB,(MYMP,HJMP)),(RMWJ,(JTQQ,PEZP)))))))))))))),(YGCX,((NNGU,OCTM),(((ATFX,((KPTE,(PUDI,AWJM)),((PKMO,WMUK),((OLXF,ZSAB),(Theca_v1.1,WYIG))))),((CRNC,((Carpa_ASGPBv0.4,CZPV),(MYZV,HYZL))),(JOIS,(DZTK,((SWPE,UAXP),(RTTY,((QSKP,(MFEA,(LVUS,UPZX))),(AXAF,(VDKG,(((CSUV,Arath_TAIR10),(VMNH,IPWB)),(TZWR,(LAPO,((BXBF,LJQF),(HABV,GTSV)))))))))))))),(IHCQ,(((JHCN,(THHD,(TIUZ,YZVJ))),(OBTI,((WAXR,(VFFP,HBHB)),((YUOM,(FCCA,JSZD)),(BCAA,((IKFD,QICX),(UHJR,(RJIM,UVDC)))))))),(((CKDK,TVCU),(YGAT,(HBUQ,ZBVT))),(WWKL,((XPBC,((RPPC,ZTLR),(((BNDE,NFXV),(FWCQ,OSIP)),((AXPJ,(KCPT,(TXMP,VXOD))),(HNCF,((BVOF,(AEPI,BHYC)),(MYVH,(POZS,OODC)))))))),((LPGY,NJLF),((COAQ,(EZZT,(SIZE,ZIWB))),(VNMY,(RVGH,(((PXYR,RHAU),(PAZJ,(VVPY,(Manes_v4.1,XNLP)))),(Poptr_v3.0,(INQX,((LFOG,GLVK),(RZTJ,(IEPQ,(TDTF,KKDQ))))))))))))))))))))))))),(((WOBD,FYSJ),(WQUF,(CVDF,(RUUB,(YNFJ,(((RXEN,(TJES,((SKNL,SHEZ),(FZQN,OLES)))),((WGET,(FVXD,Betvu_v1.1)),(((XSSD,WMLW),(SMMC,(CBJR,(AAXJ,ONLQ)))),((HDSY,PDQH),(CUTE,((OHKC,EYRD),(ZBPY,BWRK))))))),((RNBN,(KJAA,(NXTS,SCAO))),(((LKKX,CTYH),((CPKP,JLOV),((UQCB,EZGR),(LLQV,((BYNZ,KDCH),(CPLT,(GCYL,(LDEL,IWIS)))))))),((BJKT,(OMYK,(BERS,(GJNX,(HZTS,EDIT))))),(GIWN,((AZBL,SFKQ),((MRKX,BKQU),(JAFJ,(JGAB,((DVXD,EGOS),(VJPU,ZBTA)))))))))))))))),((HAEU,QUTB),(((((FWBF,BFJL),(VUSY,(VTLJ,(QURC,ZETY)))),((FNEN,WVMY),((JEXA,BNTL),((YSRZ,(KVFU,(VKJD,(PTFA,(ODDO,DAAD))))),(((JPDJ,AXBO),(NGRR,(ADHK,UFHF))),((OXYP,(WRPP,BEFC)),((HUSX,QAUE),((LRTN,SERM),(AVJK,(PPPZ,WXVX)))))))))),(((SZUO,GIPR),(QZZU,(QACK,(DFYF,(SXML,ASMV))))),(((ZJRC,HLJG),((FFFY,JTRM),(CAQZ,GSZA))),(((CLNU,CLMX),(SBZH,((MVSE,QIKZ),(((WEQK,AJFN),(CWYJ,TQKZ)),((SALZ,NTEO),(OINM,(EDBB,SUVN))))))),((IZLO,IHPC),(FXGI,(TPEM,(AUIP,(IXVJ,(HUQC,(((EYKJ,CFRN),((JNKW,RFSD),(DDRL,(FUPX,(KGJF,KFZY))))),((AFQQ,((BMSE,(DUQG,(OAGK,GUMF))),(DOVJ,(XRCX,(TEZA,(MHYG,DESP)))))),(DUNJ,(UBLN,(BIDT,((LYPZ,(AQZD,YRBQ)),((UYED,JYMN),((HXCD,QXWF),(JEPE,(NVSO,ZCUA)))))))))))))))))))),(((AHRN,(CPOC,((ALUC,(IZNU,ERWT)),((OQBM,(NHAG,QSLH)),(EMBR,(NAUM,(XQRV,VXKB))))))),((AIOU,JNVS),((MKZR,(BOLZ,(OSMU,LWCK))),(GHLP,((LQJY,(NMDZ,UGJI)),(DLJZ,(Solly_iTAGv2.3,DLAI))))))),((YQIJ,(DKFZ,(EMAL,((SMUR,HANM),(DIHD,(XVRU,((NIGS,(OEKO,MDJK)),((IDGE,OUER),(JWEY,(MZOB,ABEH)))))))))),(((WQRD,(ZSGF,(KTWL,PCNH))),((KPUM,ECTD),(HGSM,(GGJD,((MGVU,(QEHE,UOYN)),(EDEQ,((YADI,DSUV),(JGYZ,(JCLQ,YFQX))))))))),((COBX,((TORX,MZLD),(KTAR,YZGX))),((DCCI,(DTNC,RWKR)),((CLRW,((YRHD,XISJ),(GNRI,(PTBJ,YKZB)))),((EJBY,((XXYA,SIBR),(GRFT,XRLM))),(GDZS,(((HRUR,(JCMU,MXFG)),(((WOHL,EDXZ),(PCGJ,(AYIY,(NBMW,WEAC)))),((XAYK,(TKEK,UTQR)),(GCFE,(MQIV,PSHB))))),(ZRIN,(UMUL,(Mimgu_v2.0,((OWAS,(EJCM,(VYDM,(FROP,FAMO)))),(((DMLT,(RTNA,GNPX)),(ATYL,((UCNM,LRRR),(GETL,(SNNC,EAAA))))),(((BAHE,FYUH),(TAGM,(FDMM,EQDA))),(PHCE,((PUCW,FUMQ),(IYDF,(WHNV,XMBA))))))))))))))))))))))))))))))))))))))))))))));
Final quartet score is: 266745602322
Final normalized quartet score is: 0.8411102114766025
```
2. Repeat this coestimation with the KALIGN-aligned gene trees.



##Tree Visualization
The species tree files obtained using ASTRAL can be input into many tree visualization programs to make them more readable and for further analysis. Due to the high number of taxa, I chose to upload the .tre file from each MSA method into the online iTOL database and perform visual edits there, which was significantly easier and produced two species trees that were much more readable than if they had been entered into FigTree.
