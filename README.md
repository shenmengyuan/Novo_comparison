# Comparative genomics of degradative Novosphingobium analysis pipeline

>The pipeline is developed by  by Juaning Wang (juanping_wang@163.com). For questions and comments, please contact Juanping or submit an issue on github.

**The pipeline can be broadly separated into six main sections：**

1. The THN1 genome assembly, and annotation
2. Identification of orthologous proteins
3. Phylogenetic analysis
4. Insertion sequence elements, genomic islands, and CRISPR detection
5. Functional metabolic analysis and comparison
6. Comparison analysis of the microcystin-degrading mlr gene cluster


# 1 The THN1 genome assembly, and annotation

### 1.1 Assembly

*De novo* assembly was performed using the [SMRT Analysis pipeline v2.3.0](https://www.pacb.com/products-and-services/analytical-software/smrt-analysis/analysis-applications/de-novo-assembly/) in conjunction with the **HGAP assembler**  (Chien et al., 2016). Additional assemblies were performed using **minimus2** (Treangen et al., 2011). The assembly was validated by aligning the raw reads onto the finished contigs using the **BWA version 0.7.9a** (Li and Durbin, 2009).

### 1.2 Annotation

Genes in the assembled genome were predicted using **Prodigal**(Hyatt et al., 2010). tRNA, rRNA, and ncRNA were identified by **Infernal** (Nawrocki et al., 2009) and **RNAmmer** (Lagesen et al., 2007). Protein-coding sequences were annotated based on **BLASTP** searches against the downloaded **Rfam** (Griffiths-Jones et al., 2005), **NCBI NR, COG, KEGG** (Kanehisa et al., 2008), and **SwissProt** (Watanabe and Harayama, 2001) databases with an E-value cut-off of 1e−20 and subsequent filtering for the best hits.

```shell
### Gene Prediction
~/software/prodigal.v2_60/prodigal  -i  THN1.fa -f gbk -o THN1 -a THN1_pro
### NCBI NR annotation
~/software/ncbi-blast-2.2.31+/binastp -num_threads 12 -query THN1_pro.fa -out THN1.csv -db nr -outfmt 7 -evalue 1e-10 -max_target_seqs 20
```

# 2 Identification of orthologous proteins

Groups of orthologous sequences were classified by clustering with **OrthoMCL version 2.0.9** (Li et al., 2003) using a Markov cluster algorithm. Protein families were constructed using rule of **60%identity** and **60% coverage** in the alignments (Snipen and Ussery, 2010). 

```shell
### InstallSchema
cp ~/doc/OrthoMCLEngine/Main/orthomcl.config.template ./ 
orthomclInstallSchema orthomcl.config.template orthomcl.config.log 
### AdjustFasta
orthomclAdjustFasta Aro ./proteins/Novosphingobium_aromaticivorans_DSM12444.faa 1
orthomclAdjustFasta Len ./proteins/Novosphingobium_lentum_NBRC107847.faa 1
                       ······
                       ......
                       ......
### FilterFasta
orthomclFilterFasta compliantFasta/ 10 20 
### blastp
makeblastdb -in goodProteins.fasta -dbtype prot -out good_proteins.fasta
blastp -db good_proteins.fasta -query goodProteins.fasta -outfmt 7 -out goodProteins_blastp.out
grep -v -P "^#" goodProteins_blastp.out > goodProteins_v1_blastp.out
### BlastParser,LoadBlast,Pairs,DumpPairsFiles
orthomclBlastParser goodProteins_v1_blastp.out ./compliantFasta >similarSequences.txt
orthomclLoadBlast orthomcl.config.template similarSequences.txt
orthomclPairs orthomcl.config.template pairs.log cleanup=no 
orthomclDumpPairsFiles orthomcl.config.template 
### mcl,MclToGroups
mcl mclInput --abc -I 1.5 -o mclOutput
orthomclMclToGroups No_ 1 <mclOutput >groups.txt
```

# 3 Phylogenetic analysis

### 3.1 Build 16S rRNA sequences phylogenetic trees

We performed a multiple sequence alignment of 16S rRNA sequences using **Muscle**(Edgar, 2004). Unaligned sequences were trimmed from the edges. Then,phylogenetic relationships based on the 16S rRNA sequences of the *Novosphingobium* strains were analyzed using **MEGA v7.0.26** (Kumar et al., 2016) with the neighbor-joining method. The robustness of clustering was evaluated by 1000 bootstrap replicates.

### 3.2 Build core genes phylogenetic trees

The core genes were computed by **OrthoMCLversion 2.0.9** (Li et al., 2003) and retrieved using inhouse scripts. We performed a multiple sequence alignment of the core genes using **Clustal Omega**(Sievers et al., 2011). The edge-trimmed core sequences were joined and the phylogenetic tree was constructed using **MEGA v7.0.26**. 

```shell
Count_core.py
ExtractCoreSequence.pl
ClustalOmegaAlignment.sh
CoreEdgeTrim.sh
CoreJoin.pl
```

### 3.3  Calculate ANI values

Average nucleotide identity (ANI) values were calculated using the MUMmer algorithm of **JSpecies v 1.2.1** (Kurtz et al., 2004; Richter andRossello-Mora, 2009; Chan et al., 2012) and visualized as a heatmap using **Morpheus software** (https://software.broadinstitute.org/morpheus/).

```shell
java -jar -Xms512m -Xmx512m ~/software/jspecies/jspecies1.2.1.jar
```

### 3.4 Build whole-genome-based phylogeny trees

 We employed **CVTree (version 3.9.6)** (Qi et al., 2004) with a K value of 6 to compute the whole-genome composition vector of the 22 *Novosphingobium* strains. The robustness of this whole-genome-based phylogeny was validated by performing statistical resampling tests by 1000 bootstrap replicates using **Seqkit** (Shen et al., 2016). **PHYLIP** (Retief, 2000) was used to construct an eighbor-joining phylogenetic tree, which was visualized using **MEGA v7.0.26**(Kumar et al., 2016). 

```shell
cvtree -i species.list -p data -o CVTree_k6.txt -k 6
neighbor
```

# 4 Insertion sequence elements, genomic islands, and CRISPR detection

Insertion sequences (ISs) were detected by **BLAST comparisons** (E-value ≤1e−5) against the ISFinder database (Siguier et al., 2006). **IslandViewer 4** (Bertelli et al., 2017) was used to predict genomic islands (GIs), as described previously (Zhang et al., 2016). CRISPR arrays were detected using the **CRISPRFinder** (Grissa et al., 2007) online server to perform BLAST searches against dbCRISPR (CRISPR database). 

# 5 Functional metabolic analysis and comparison

### 5.1 COG

We assigned the protein-coding genes to COG categories by searches against eggNOG version 4.5.1 (Huerta-Cepas et al., 2016) using the **eggNOG-mapper tool**(Huerta-Cepas et al., 2017). 

```shell
python eggnog-mapper-1.0.3/emapper.py --cpu 30 -i total.fa  -d bact -o total
```

### 5.2 KEGG

Predicted proteins from the *Novosphingobium* genomes were also annotated using the **KEGG KAAS server** (Moriya et al., 2007). Metabolic pathways were deduced by manual inspection of KEGG Orthology (KO) that was predicted based on comparisons to the KEGG database (Kanehisa and Goto, 2000; Kanehisa et al., 2014).

# 6 Comparison analysis of the microcystin-degrading mlr gene cluster

We compared the *mlr* genes of each pair of strains by **BLAST alignments**. Codon usage of selected core genes and *mlr* genes was computed and analyzed by graphical codon usage analyser (<http://gcua.schoedl.de>) platform.

# Datasets

```txt
Datasets/
├── Genomes
│   ├── Novosphingobium_aromaticivorans_DSM12444.fna
│   ├── Novosphingobium_barchaimii_LL02.fna
│   ├── Novosphingobium_lentum_NBRC_107847.fna
│   ├── Novosphingobium_lindaniclasticum_LE124.fna
│   ├── Novosphingobium_naphthalenivorans_NBRC_102051.fna
│   ├── Novosphingobium_panipatense_P5.fna
│   ├── Novosphingobium_pentaromativorans_US6-1.fna
│   ├── Novosphingobium_resinovorum_SA1.fna
│   ├── Novosphingobium_sp_63-713.fna
│   ├── Novosphingobium_sp_B-7.fna
│   ├── Novosphingobium_sp_Chol11.fna
│   ├── Novosphingobium_sp_Fuku2-ISO-50.fna
│   ├── Novosphingobium_sp_KN65.fna
│   ├── Novosphingobium_sp_MBES04.fna
│   ├── Novosphingobium_sp_P6W.fna
│   ├── Novosphingobium_sp_PC22D.fna
│   ├── Novosphingobium_sp_PP1Y.fna
│   ├── Novosphingobium_sp_SCN63-17.fna
│   ├── Novosphingobium_sp_SCN64-18.fna
│   ├── Novosphingobium_sp_ST904.fna
│   ├── Novosphingobium_sp_THN1.fna
│   └── Novosphingobium_subterraneum_NBRC16086.fna
└── Proteins
    ├── Novosphingobium_aromaticivorans_DSM12444.faa
    ├── Novosphingobium_barchaimii_LL02.faa
    ├── Novosphingobium_lentum_NBRC107847.faa
    ├── Novosphingobium_lindaniclasticum_LE124.faa
    ├── Novosphingobium_naphthalenivorans_NBRC_102051.faa
    ├── Novosphingobium_panipatense_P5.faa
    ├── Novosphingobium_pentaromativorans_US6-1.faa
    ├── Novosphingobium_resinovorum_SA1.faa
    ├── Novosphingobium_sp_63-713.faa
    ├── Novosphingobium_sp_B-7.faa
    ├── Novosphingobium_sp_Chol11.faa
    ├── Novosphingobium_sp_Fuku2-ISO-50.faa
    ├── Novosphingobium_sp_KN65.faa
    ├── Novosphingobium_sp_MBES04.faa
    ├── Novosphingobium_sp_P6W.faa
    ├── Novosphingobium_sp_PC22D.faa
    ├── Novosphingobium_sp_PP1Y.faa
    ├── Novosphingobium_sp_SCN63-17.faa
    ├── Novosphingobium_sp_SCN64-18.faa
    ├── Novosphingobium_sp_ST904.faa
    ├── Novosphingobium_sp_THN1.faa
    └── Novosphingobium_subterraneum_NBRC16086.faa
```

# Reference
- Bertelli, C., Laird, M.R., Williams, K.P., Lau, B.Y., Hoad, G., Winsor, G.L., et al. (2017). IslandViewer 4: expanded prediction of genomic islands for larger-scale datasets. Nucleic Acids Research 45(W1), W30-W35. doi: 10.1093/nar/gkx343.
- Chan, J.Z., Halachev, M.R., Loman, N.J., Constantinidou, C., and Pallen, M.J. (2012). Defining bacterial species in the genomic era: insights from the genus Acinetobacter. BMC Microbiol 12, 302. doi: 10.1186/1471-2180-12-302.
- Chien, J.T., Pakala, S.B., Geraldo, J.A., Lapp, S.A., Humphrey, J.C., Barnwell, J.W., et al. (2016). High-quality genome assembly and annotation for Plasmodium coatneyi, generated using single-molecule real-time PacBio technology. Genome announcements 4(5). doi: 10.1128/genomeA.00883-16.
- Edgar, R.C. (2004). MUSCLE: multiple sequence alignment with high accuracy and high throughput. Nucleic Acids Research 32(5), 1792-1797. doi: 10.1093/nar/gkh340.
- Griffiths-Jones, S., Moxon, S., Marshall, M., Khanna, A., Eddy, S.R., and Bateman, A. (2005). Rfam: annotating non-coding RNAs in complete genomes. Nucleic Acids Research 33, D121-D124. doi: 10.1093/nar/gri081.
- Grissa, I., Vergnaud, G., and Pourcel, C. (2007). CRISPRFinder: a web tool to identify clustered regularly interspaced short palindromic repeats. Nucleic Acids Research 35, W52-W57. doi: 10.1093/nar/gkm360.
- Huerta-Cepas, J., Forslund, K., Coelho, L.P., Szklarczyk, D., Jensen, L.J., von Mering, C., et al. (2017). Fast genome-wide functional annotation through orthology assignment by eggNOG-Mapper. Molecular Biology and Evolution 34(8), 2115-2122. doi: 10.1093/molbev/msx148.
- Huerta-Cepas, J., Szklarczyk, D., Forslund, K., Cook, H., Heller, D., Walter, M.C., et al. (2016). eggNOG 4.5: a hierarchical orthology framework with improved functional annotations for eukaryotic, prokaryotic and viral sequences. Nucleic Acids Research 44(D1), D286-D293. doi: 10.1093/nar/gkv1248.
- Hyatt, D., Chen, G.L., Locascio, P.F., Land, M.L., Larimer, F.W., and Hauser, L.J. (2010). Prodigal: prokaryotic gene recognition and translation initiation site identification. BMC Bioinformatics 11, 119. doi: 10.1186/1471-2105-11-119.
- Kanehisa, M., and Goto, S. (2000). KEGG: kyoto encyclopedia of genes and genomes. Nucleic Acids Res 28(1), 27-30.
- Kanehisa, M., Araki, M., Goto, S., Hattori, M., Hirakawa, M., Itoh, M., et al. (2008). KEGG for linking genomes to life and the environment. Nucleic Acids Research 36, D480-D484. doi: 10.1093/nar/gkm882.
- Kanehisa, M., Goto, S., Sato, Y., Kawashima, M., Furumichi, M., and Tanabe, M. (2014). Data, information, knowledge and principle: back to metabolism in KEGG. Nucleic Acids Research 42(D1), D199-D205. doi: 10.1093/nar/gkt1076.
- Kumar, S., Stecher, G., and Tamura, K. (2016). MEGA7: molecular evolutionary genetics analysis version 7.0 for bigger datasets. Molecular Biology and Evolution 33(7), 1870-1874. doi: 10.1093/molbev/msw054.
- Kurtz, S., Phillippy, A., Delcher, A.L., Smoot, M., Shumway, M., Antonescu, C., et al. (2004). Versatile and open software for comparing large genomes. Genome Biol 5(2), R12. doi: 10.1186/gb-2004-5-2-r12.
- Lagesen, K., Hallin, P., Rodland, E.A., Staerfeldt, H.H., Rognes, T., and Ussery, D.W. (2007). RNAmmer: consistent and rapid annotation of ribosomal RNA genes. Nucleic Acids Research 35(9), 3100-3108. doi: 10.1093/nar/gkm160.
- Li, H., and Durbin, R. (2009). Fast and accurate short read alignment with Burrows-Wheeler Transform. Bioinformatics. 25, 1754–1760. doi: 10.1093/bioinformaticsp324
- Li, L., Stoeckert, C.J., and Roos, D.S. (2003). OrthoMCL: Identification of ortholog groups for eukaryotic genomes. Genome Research 13(9), 2178-2189. doi: 10.1101/gr.1224503.
- Moriya, Y., Itoh, M., Okuda, S., Yoshizawa, A.C., and Kanehisa, M. (2007). KAAS: an automatic genome annotation and pathway reconstruction server. Nucleic Acids Research 35, W182-W185. doi: 10.1093/nar/gkm321.
- Nawrocki, E.P., Kolbe, D.L., and Eddy, S.R. (2009). Infernal 1.0: inference of RNA alignments. Bioinformatics 25(10), 1335-1337. doi: 10.1093/bioinformatics/btp157.
- Qi, J., Wang, B., and Hao, B.I. (2004). Whole proteome prokaryote phylogeny without sequence alignment: A K-string composition approach. Journal of Molecular Evolution 58(1), 1-11. doi: 10.1007/s00239-003-2493-7.
- Richter, M., and Rossello-Mora, R. (2009). Shifting the genomic gold standard for the prokaryotic species definition. Proc Natl Acad Sci U S A 106(45), 19126-19131. doi: 10.1073/pnas.0906412106.
- Shen, W., Le, S., Li, Y., and Hu, F. (2016). SeqKit: a cross-platform and ultrafast toolkit for FASTA/Q file manipulation. PLoS One 11(10), e0163962. doi: 10.1371/journal.pone.0163962.
- Sievers, F., Wilm, A., Dineen, D., Gibson, T.J., Karplus, K., Li, W.Z., et al. (2011). Fast, scalable generation of high-quality protein multiple sequence alignments using Clustal Omega. Molecular Systems Biology 7. doi: 10.1038/msb.2011.75.
- Snipen, L., and Ussery, D.W. (2010). Standard operating procedure for computing pangenome trees. Standards in Genomic Sciences 2(1), 135-141. doi: 10.4056/sigs.38923.
- Treangen, T.J., Sommer, D.D., Angly, F.E., Koren, S., and Pop, M. (2011). Next generation sequence assembly with AMOS. Curr Protoc Bioinformatics Chapter 11, Unit 11 18. doi: 10.1002/0471250953.bi1108s33.
- Watanabe, K., and Harayama, S. (2001). SWISS-PROT: the curated protein sequence database on Internet. Tanpakushitsu Kakusan Koso 46(1), 80-86.





Copyright (C) 2018 Juanping Wang.
