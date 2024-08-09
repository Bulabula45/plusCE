# ANALYSIS--PART I

### Prepare dependencies

**software**

This project requires Python version >= 3.9, uses cd-hit 4.8.1 and RNAfold 2.6.4, bedtools 2.31.0.



**python library**

`pip install -r requirements.txt`

<br>

### Source data and preprocess

| Database          | url                                                          |
| ----------------- | ------------------------------------------------------------ |
| Ensembl ncRNA     | https://ftp.ensembl.org/pub/release-112/fasta/homo_sapiens/ncrna/ |
| Smprot            | http://bigdata.ibp.ac.cn/SmProt/download.htm                 |
| SPENCER           | https://spencer.renlab.org/#/download                        |
| mRNA MANE         | https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.3/ |
| FANTOM            | https://fantom.gsc.riken.jp/5/suppl/Hon_et_al_2016/data/assembly/lv1_raw/ |
| GENCODE           | https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.transcripts.fa.gz |
| AHIGH(annotation) | https://ars.els-cdn.com/content/image/1-s2.0-S1097276522006062-mmc3.xls |

<br>

### Codes

1. Run  `data_preprocess.ipynb` to generate transcripts and ORF annotation files for MANE mRNA and ncRNA_ribo from Smprot, SPENCER, and predicted ORFs of ncRNA.

   We also provide the above files associated with ncRNA_ribo in `ribo_rawfiles.zip`.

2. Run `data_produce.ipynb` to generate 50nt downstream start codon, 5' UTR, ORF sequence files of mRNA, ncRNA_ribo, ncRNA.

3. Calculate features and plot.

   3.1 RNA folding energy: `RNAfold < mRNA/ribo/ncRNA_d50.fasta > res_mRNA/ribo/ncRNA_d50`.

   3.2 Other features include edit distance, uStart, uORF, last uORF distance, pair kozak distances, fickett score and hexamer score are in `Cal_features.ipynb`.

<br>

# ANALYSIS--PART II

### Source data

| Data              | url                                                          |
| ----------------- | ------------------------------------------------------------ |
| AHIGH(expression) | https://ars.els-cdn.com/content/image/1-s2.0-S1097276522006062-mmc4.xlsx |
| Phast score       | https://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/phastConsElements30way.txt.gz |

<br>

### Codes

1. Run `plot.ipynb` to analysis translational features of mRNA with long 5' UTR grouped by translational efficiency.

2. Run `ConservedEle.ipynb` to intersect 5' UTR with conserved elements within genome, and compare the distribution of mRNA with long 5' UTR grouped by TE.

<br>

# MODEL

### Prepare dependencies

**software**

This project requires Python version >= 3.9, uses RNAfold 2.6.4.

**python library**

`pip install -r requirements.txt`

<br>

### Source data

| Data                                     | url                                                        |
| ---------------------------------------- | ---------------------------------------------------------- |
| 5'UTR sequences and annotations from Cao | https://github.com/zzz2010/5UTR_Optimizer/tree/master/data |

<br>

### Codes

1. Download the sequence file and annotation files.

   ```bash
   mkdir data output
   
   wget https://github.com/zzz2010/5UTR_Optimizer/blob/master/data/gencode_v17_5utr_15bpcds.fa ./data
   
   wget https://github.com/zzz2010/5UTR_Optimizer/blob/master/data/df_counts_and_len.TE_sorted.HEK_Andrev2015.with_annot.txt ./data
   
   wget https://github.com/zzz2010/5UTR_Optimizer/blob/master/data/df_counts_and_len.TE_sorted.Muscle.with_annot.txt ./data
   
   wget https://github.com/zzz2010/5UTR_Optimizer/blob/master/data/df_counts_and_len.TE_sorted.Muscle.with_annot.txt ./data
   ```

   

2. Prepare features of all sequences(except conserved elements) and database for querying conserved elements.

   ```bash
   python FeatureExtraction_final.py data/gencode_v17_5utr_15bpcds.fa output/
   ```

   Run plusCE_database.ipynb.

3. Train model and output spearman $\rho$.

   ```bash
   python plusCE.py --data data --annot df_counts_and_len.TE_sorted.Muscle.with_annot.txt --feature output --rna 5 --ribo 0.1 --querydb DB
   ```



<br>

# Citations

Li Y, Zhou H, Chen X, Zheng Y, Kang Q, Hao D, Zhang L, Song T, Luo H, Hao Y, et al: SmProt: A Reliable Repository with Comprehensive Annotation of Small Proteins Identified from Ribosome Profiling**.** *Genomics Proteomics Bioinformatics* 2021, 19**:**602-610. http://dx.doi.org/10.1016/j.gpb.2021.09.002.

Luo X, Huang Y, Li H, Luo Y, Zuo Z, Ren J, Xie Y: SPENCER: a comprehensive database for small peptides encoded by noncoding RNAs in cancer patients**.** *Nucleic Acids Res* 2022, 50**:**D1373-D1381. http://dx.doi.org/10.1093/nar/gkab822.

Chothani SP, Adami E, Widjaja AA, Langley SR, Viswanathan S, Pua CJ, Zhihao NT, Harmston N, D'Agostino G, Whiffin N, et al: A high-resolution map of human RNA translation**.** *Mol Cell* 2022, 82**:**2885-2899 e2888. http://dx.doi.org/10.1016/j.molcel.2022.06.023.

Cao J, Novoa EM, Zhang Z, Chen WCW, Liu D, Choi GCG, Wong ASL, Wehrspaun C, Kellis M, Lu TK: High-throughput 5' UTR engineering for enhanced protein production in non-viral gene therapies**.** *Nat Commun* 2021, 12**:**4138. http://dx.doi.org/10.1038/s41467-021-24436-7.
