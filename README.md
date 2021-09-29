# Repository for scripts used in Molecular evolution and depth-related adaptations in rhodopsin in the adaptive radiation of cichlid fishes in Lake Tanganyika
Virginie Ricci, Fabrizia Ronco, Zuzana Musilova & Walter Salzburger (2021)

## Scripts

* 01_Reads_mapping.sh: Illumina raw reads mapping of Tanganyikan cichlid species to Nile tilapia RefSeq (RefSeq accession GCF_001858045.2)
* 02_Extract_BAM_region.sh: Extraction of reads mapped to the RefSeq opsin mRNA
* In-between steps on Geneious (see Materials and Methods)
* 03_MAFFT.sh: Multiple alignment of all consensus opsin mRNA including RefSeq mRNA, extraction of all consensus opsin CDS including RefSeq CDS
* 04_Filter_CDS.sh: Check for the presence of start/stop codons, early stop codons and sequence lengths dividable by 3, opsin CDS translation to AA
* 05_BLASTN.sh: BLASTN of opsin CDS to Tanganyikan cichlid genome assemblies
* 06_Reads_coverage.sh: Extraction of mean/median reads coverage in the entired RefSeq genome, extraction of mean/median reads coverage in the RefSeq CDS
* 07_SubstitutionModel_tests.sh: Identification of appropriate nucleotide/amino acid substitution models to build phylogenetic trees
* 08_PhylogeneticTrees.sh: Build nucleotide/amino acid phylogenetic trees
* 09_TopologyTest.sh: Phylogenetic topology tests
* 10_AAsubstitutions_on_BestTree.sh: Mapping of AA substitutions on the best phylogenetic trees
* 11_CodeML.sh: Positive selection test using CodeML
* 12_HyPhy.sh: Positive selection test using HyPhy 

