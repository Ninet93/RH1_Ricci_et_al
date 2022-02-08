# Molecular evolution and depth-related adaptations in rhodopsin in the adaptive radiation of cichlid fishes in Lake Tanganyika
Virginie Ricci, Fabrizia Ronco, Zuzana Musilova & Walter Salzburger (2021)

## Scripts

* [`01_Reads_mapping.sh`](Scripts/01_Reads_mapping.sh): Illumina raw reads mapping of Tanganyikan cichlid species to Nile tilapia RefSeq (RefSeq accession GCF_001858045.2)
  * BWA and SAMtools

* [`02_Extract_BAM_region.sh`](Scripts/02_Extract_BAM_region.sh): Extraction of reads mapped to the RefSeq opsin mRNA
  * SAMtools

* In-between steps on Geneious (https://www.geneious.com, see Materials and Methods)

* [`03_MAFFT.sh`](Scripts/03_MAFFT.sh): Multiple alignment of all consensus opsin mRNA including RefSeq mRNA, extraction of all consensus opsin CDS including RefSeq CDS
  * Python, Biopython, Extract_CDS.py and MAFFT

* [`04_Filter_CDS.sh`](Scripts/04_Filter_CDS.sh): Check for the presence of start/stop codons, early stop codons and sequence lengths dividable by 3, opsin CDS translation to AA
  * Python, Biopython, and Filter_CDS.py

* [`05_BLASTN.sh`](Scripts/05_BLASTN.sh): BLASTN of opsin CDS to Tanganyikan cichlid genome assemblies
  * Python, Biopython, BLAST, and Extract_individual_CDS.py

* [`06_Reads_coverage.sh`](Scripts/06_Reads_coverage.sh): Extraction of mean/median reads coverage in the entire RefSeq genome, extraction of mean/median reads coverage in the RefSeq CDS
  * SAMtools, R, and Reads_coverage.R

* [`07_SubstitutionModel_tests.sh`](Scripts/07_SubstitutionModel_tests.sh): Identification of appropriate nucleotide/amino acid substitution models to build phylogenetic trees
  * jModelTest, PAUP, and ProtTest

* [`08_PhylogeneticTrees.sh`](Scripts/08_PhylogeneticTrees.sh): Build nucleotide/amino acid phylogenetic trees
  * MrBayes, Beagle, MrBayes_CDS.txt, MrBayes_AA.txt, and IQ-TREE

* [`09_TopologyTest.sh`](Scripts/09_TopologyTest.sh): Phylogenetic topology tests
  * PAUP

* Haplotype Viewer (http://www.cibiv.at/~greg/haploviewer, see Materials and Methods)

* [`10_AAsubstitutions_on_BestTree.sh`](Scripts/10_AAsubstitutions_on_BestTree.sh): Mapping of AA substitutions on the best phylogenetic trees
  * PAUP, R, and PAUP_R_trees.R

* [`11_CodeML.sh`](Scripts/11_CodeML.sh): Positive selection test using CodeML
  * PAML, CodeML_M1a.ctl, CodeML_M2a.ctl, CodeML_M7.ctl, CodeML_M8.ctl,  CodeML_branchsite_H0.ctl, and CodeML_branchsite_HA.ctl

* [`12_HyPhy.sh`](Scripts/12_HyPhy.sh): Positive selection test using HyPhy
  * HyPhy, HyPhy_FEL.txt, HyPhy_FUBAR.txt, HyPhy_aBSREL_deep.txt, and HyPhy_aBSREL_shallow.txt

* [`13_BayesTraits.sh`](Scripts/13_BayesTraits.sh): Depth-related opsin substitutions analysis
  * R, BayesTraits, BayesTraits_inputs.R, parfile_dependent_MC.txt, parfile_independent_MC.txt, and BayesTraits.R

## Data

### Inputs

* [`Dataset.txt`](Data/Dataset.txt): Table of the dataset

* [`Color_tribes.txt`](Data/Color_tribes.txt): Color of tribes

* [`Nucleotide_table.txt`](Data/Nucleotide_table.txt): Table of nucleotides

* [`AA_table.txt`](Data/AA_table.txt): Table of amino acids

* [`b1.tre`](Data/b1.tre) species tree (see Ronco et al. 2021)

* [`Coo_mRNA_Exon_CDS_RefSeq.txt`](Data/Coo_mRNA_Exon_CDS_RefSeq.txt): coordinates of opsin in the RefSeq genome

### Outputs

* [`RH1_CDS_w_haplotypes.fasta`](Data/RH1_CDS_w_haplotypes.fasta): RH1 CDS multiple alignment (haplo0: RefSeq CDS)

* [`RH1_codons_manualcorrection.txt`](Data/RH1_codons_manualcorrection.txt): Table of manual corrections after visual inspection of reads mapping on Geneious

* [`RH1_ALL.blastn`](Data/RH1_ALL.blastn): RH1 BLASTN results

* [`exoRH1_CDS_w_indiv.fasta`](Data/exoRH1_CDS_w_indiv.fasta): exoRH1 CDS multiple alignment (XM_003438995_5: RefSeq CDS)

* [`exoRH1_ALL.blastn`](Data/exoRH1_ALL.blastn): exoRH1 BLASTN results

* IQ-Tree phylogenetic trees:
  * 1) CDS: [`RH1_CDS_w_haplotypes_outgroup_bootstrap_GTR_I_G_1.treefile`](Data/RH1_CDS_w_haplotypes_outgroup_bootstrap_GTR_I_G_1.treefile) (BestTree)
  * 2) CDS: RH1_CDS_w_haplotypes_outgroup_bootstrap_GTR_I_G_2.treefile
  * 3) CDS: RH1_CDS_w_haplotypes_outgroup_bootstrap_GTR_I_G_allnni_1.treefile
  * 4) CDS: RH1_CDS_w_haplotypes_outgroup_bootstrap_GTR_I_G_allnni_2.treefile

  * 5) AA: RH1_AA_w_haplotypes_outgroup_bootstrap_JTT_I_G_F_1.treefile
  * 6) AA: RH1_AA_w_haplotypes_outgroup_bootstrap_JTT_I_G_F_2.treefile
  * 7) AA: [`RH1_AA_w_haplotypes_outgroup_bootstrap_JTT_I_G_F_allnni_1.treefile`](Data/RH1_AA_w_haplotypes_outgroup_bootstrap_JTT_I_G_F_allnni_1.treefile) (BestTree)
  * 8) AA: RH1_AA_w_haplotypes_outgroup_bootstrap_JTT_I_G_F_allnni_2.treefile

* MrBayes phylogenetic trees:
  * 1) CDS: RH1_CDS_w_haplotypes_1.fasta.nexus.con.tre
  * 2) CDS: RH1_CDS_w_haplotypes_2.fasta.nexus.con.tre
  * 3) CDS: RH1_CDS_w_haplotypes_3.fasta.nexus.con.tre
  * 4) CDS: RH1_CDS_w_haplotypes_4.fasta.nexus.con.tre

  * 5) AA: RH1_AA_w_haplotypes_1.fasta_wo_ambiguouschar.nexus.con.tre
  * 6) AA: RH1_AA_w_haplotypes_2.fasta_wo_ambiguouschar.nexus.con.tre
  * 7) AA: RH1_AA_w_haplotypes_3.fasta_wo_ambiguouschar.nexus.con.tre
  * 8) AA: RH1_AA_w_haplotypes_4.fasta_wo_ambiguouschar.nexus.con.tre

* [`Extra figures`](Data/Extra_figures):
  * AA changes at variable RH1 sites mapped on the species tree (see Ronco et al., 2021)
