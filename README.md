# Molecular evolution and depth-related adaptations in rhodopsin in the adaptive radiation of cichlid fishes in Lake Tanganyika
Virginie Ricci, Fabrizia Ronco, Zuzana Musilova & Walter Salzburger (2021)

## Scripts

* 01_Reads_mapping.sh: Illumina raw reads mapping of Tanganyikan cichlid species to Nile tilapia RefSeq (RefSeq accession GCF_001858045.2)
  * WA and SAMtools

* 02_Extract_BAM_region.sh: Extraction of reads mapped to the RefSeq opsin mRNA
  * SAMtools

* In-between steps on Geneious (see Materials and Methods)

* 03_MAFFT.sh: Multiple alignment of all consensus opsin mRNA including RefSeq mRNA, extraction of all consensus opsin CDS including RefSeq CDS
  * Python, Biopython, Extract_CDS.py and MAFFT

* 04_Filter_CDS.sh: Check for the presence of start/stop codons, early stop codons and sequence lengths dividable by 3, opsin CDS translation to AA
  * Python, Biopython, and Filter_CDS.py

* 05_BLASTN.sh: BLASTN of opsin CDS to Tanganyikan cichlid genome assemblies
  * Python, Biopython, BLAST, and Extract_individual_CDS.py

* 06_Reads_coverage.sh: Extraction of mean/median reads coverage in the entired RefSeq genome, extraction of mean/median reads coverage in the RefSeq CDS
  * SAMtools, R, and Reads_coverage.R

* 07_SubstitutionModel_tests.sh: Identification of appropriate nucleotide/amino acid substitution models to build phylogenetic trees
  * jModelTest, PAUP, and ProtTest

* 08_PhylogeneticTrees.sh: Build nucleotide/amino acid phylogenetic trees
  * MrBayes, Beagle, MrBayes_CDS.txt, MrBayes_AA.txt, and IQ-TREE

* 09_TopologyTest.sh: Phylogenetic topology tests
  * PAUP

* 10_AAsubstitutions_on_BestTree.sh: Mapping of AA substitutions on the best phylogenetic trees
  * PAUP, R, and PAUP_R_trees.R

* 11_CodeML.sh: Positive selection test using CodeML
  * PAML, CodeML_M1a.ctl, CodeML_M2a.ctl, CodeML_M7.ctl, and CodeML_M8.ctl

* 12_HyPhy.sh: Positive selection test using HyPhy
  * HyPhy, HyPhy_FEL.txt, HyPhy_FUBAR.txt, HyPhy_aBSREL_deep.txt, and HyPhy_aBSREL_shallow.txt

