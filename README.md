<h1>First steps - directory structure and setting up conda</h1>

All analyses described below were performed on the Lund University Bioinformatics server. The software conda (v4.11.0) was already installed; detailed installation instructions can be found at https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html. To facilitate the installation of the software used in the analyses described below, specific channels should be configured before setting up the working environment.

```bash
#After installing conda, configure the channels
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
#Then set up a new environment dedicated to the project
conda create -n MasterThesis
conda activate MasterThesis
#This environment runs python v3.6.15
```

For this project, all files are stored within a master directory. Directories used for storing files associated with specific steps of the project are described below as each part of the analysis is performed. All given absolute paths consider the main project directory as "root" and should be adapted accordingly when recreating the analyses.

```bash
#Main directory that holds all files associated with the project
mkdir MasterThesis
#Directory used for storing all in-house scripts
mkdir MasterThesis/scripts
```

<h1>Identification and phylogenetic analysis of RQUA sequences</h1>

```bash
#New directory to hold the files associated with this step of the analysis
mkdir MasterThesis/01_phylogenetic_analysis
```

<h2>Identification of RQUA homologues on EukProt</h2>

Two versions of the EukProt database are used: v2, which was the only released version at the time of running these analyses, and a subset of v3 that was received through personal communication before publication.

First, the databases are downloaded, indexed using an in-house python script, concatenated into one FASTA file for each database, and converted to a single-line format.

<h3>Downloading and processing EukProt</h3>

```bash
#Create directories to store the two versions of EukProt
mkdir MasterThesis/01_phylogenetic_analysis/EukProt2
mkdir MasterThesis/01_phylogenetic_analysis/EukProt3
#Download a compressed version of EukProt v2 into the respective directory
cd MasterThesis/01_phylogenetic_analysis/EukProt2
wget https://figshare.com/ndownloader/articles/12417881/versions/2
#Uncompress the main file (called "2" in this case) and resulting datasets
unzip 2
tar -xvzf EukProt_assembled_transcriptomes.v02.2020_06_30.tgz
tar -xvzf EukProt_proteins.v02.2020_06_30.tgz
tar -xvzf EukProt_unannotated_genomes.v02.2020_06_30.tgz
#EukProt v3 was received through personal communication and placed in the corresponding directory
cd MasterThesis/01_phylogenetic_analysis/EukProt3
#Uncompress all datasets
tar -xvzf EukProt_assembled_transcriptomes.v03.2021_11_22.tgz
tar -xvzf EukProt_genome_annotations.v03.2021_11_22.tgz
tar -xvzf EukProt_proteins.v03.2021_11_22.tgz
```

The data that will be used in the subsequent analyses is found inside the _proteins_ directories.

The predicted proteomes from EukProt v2 contain headers that are inconsistent and difficult to parse. An in-house python script (called `fasta_index.py`) is used to rename these headers that adds the names of the organisms that the predicted proteomes belong to to the headers. This step is not neccessary for EukProt v3.

```bash
cd MasterThesis/01_phylogenetic_analysis/EukProt2
#Make directory to store the indexed fasta files and indexing files produced by the python script
mkdir MasterThesis/01_phylogenetic_analysis/EukProt2/Indexing
cd MasterThesis/01_phylogenetic_analysis/EukProt2/Indexing
#fasta_index.py is located inside MasterThesis/scripts
ls ../proteins/*.fasta | while read file; do python3 ../../../scripts/fasta_index.py -f $file -o Indexed_$(basename $file); done
```

The resulting files are fasta files beginning with "Indexed", which are the actual indexed fasta files, and .tsv files beginning with "Index", which contain information needed for unindexing.

Then, both versions of the EukProt databases are concatenated to facilitate subsequent analyses.

```bash
cd MasterThesis/01_phylogenetic_analysis/EukProt2/Indexing
#Concatenate EukProt v2
cat Indexed_* > all_indexed_eukprot.fasta
#Remove the individual files to save space
rm Indexed_*
#Concatenate EukProt v3
cd MasterThesis/01_phylogenetic_analysis/EukProt3/proteins
cat EP* > all_eukprot3.fasta
```

Next, the two EukProt databases are converted to single-line format.

```bash
#EukProt v2
cd MasterThesis/01_phylogenetic_analysis/EukProt2/Indexing
cat all_indexed_eukprot.fasta | awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' | tail -n +2 > all_indexed_eukprot_oneline.fasta
#EukProt v3
cd MasterThesis/01_phylogenetic_analysis/EukProt3/proteins
cat all_eukprot3.fasta | awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' | tail -n +2 > all_eukprot3_oneline.fasta
```

<h3>Running BLAST for RQUA on EukProt</h3>

All files related to this step of the analysis are stored in a separate directory.

```bash
mkdir MasterThesis/01_phylogenetic_analysis/BLAST_RQUA_EukProt
cd MasterThesis/01_phylogenetic_analysis/BLAST_RQUA_EukProt
```

In the same directory, a fasta file containing the query RQUA sequence from _Pygsuia biforma_ (GenBank: AKA62179.1) is placed and named `pygsuya_rqua.fasta`. The amino acid sequence can be obtained from https://www.uniprot.org/uniprot/A0A0K0MD65.

To perform the actual BLASTP analyses, the EukProt datasets are first converted into databases. Then, the actual BLASTP analysis is performed.

```bash
#First, install BLAST
conda install blast=2.12.0
#Make database of EukProt v2
makeblastdb -in ../EukProt2/Indexing/all_indexed_eukprot.fasta -out EukProt -dbtype prot
#Make database of EukProt v3
makeblastdb -in ../EukProt3/proteins/all_eukprot3.fasta -out EukProt3 -dbtype prot
#Run the two BLASTP analyses
blastp -query pygsuia_rqua.fasta -db EukProt -evalue 1e-10 -out eukprot_rqua.blastp -num_threads 2
blastp -query pygsuia_rqua.fasta -db EukProt3 -evalue 1e-10 -out eukprot3_rqua.blastp -num_threads 2
```

Then, hits are retrieved from the BLASTP results files and placed into FASTA files.

```bash
grep "^EP00" eukprot_rqua.blastp | cut -d " " -f 1 | while read line; do grep -A 1 $line$ ../EukProt2/Indexing/all_indexed_eukprot_oneline.fasta; done > eukprot_rqua.fasta
grep "^EP00" eukprot3_rqua.blastp | cut -d " " -f 1 | while read line; do grep -A 1 $line ../EukProt3/proteins/all_eukprot3_oneline.fasta; done > eukprot3_rqua.fasta
```

<h2>Retrieval of RQUA sequence from other eukaryotic species</h2>


<h3>Other Porifera species</h3>

Datasets belonging to diverse representatives of phylum Porifera were downloaded from the following links:

 - _Ephydatia muelleri_ genome : https://spaces.facsci.ualberta.ca/ephybase/resources/, BioProject: PRJNA579531
 - _Ephydatia muelleri_ transcriptome : https://era.library.ualberta.ca/items/d643c77d-aefd-45f3-8ddb-985189eebcf2
 - _Spongilla lacustris_ transcriptome 1 : https://era.library.ualberta.ca/items/e5895ef8-cc1d-4cb9-bfb9-9516d3bdfa9e
 - _Spongilla lacustris_ transcriptome 2 : https://git.embl.de/musser/profiling-cellular-diversity-in-sponges-informs-animal-cell-type-and-nervous-system-evolution/-/tree/master/transcriptome_proteome_and_phylome
 - _Eunapius fragilis_ transcriptome : https://era.library.ualberta.ca/items/6139a88f-895d-44a7-bd0d-e22d455d2785
 - _Aphrocallistes vastus_ transcriptome : https://era.library.ualberta.ca/items/3df7eb24-1f79-43d9-bcb4-49ef043ea226
 - _Sycon coactum_ transcriptome : https://era.library.ualberta.ca/items/0ecd7382-3dc1-43f1-b2dd-050b4e9e159a

All datasets were placed inside a new directory:

```bash
mkdir MasterThesis/02_sponge_datasets
```

All datasets are converted to BLAST databases, then the RQUA sequence of _E. muelleri_ (file `E_muelleri_RQUA.fasta`) is used as query to identify potential RQUA homologues. For the _Ephydatia fluviatilis_ transcriptome (BioProject: PRJNA244851) the TBLASTN analysis was performed separately, using the web version of BLAST, since no assembly is available.

```bash
#Directory to store all files associated with the BLAST analysis
mkdir MasterThesis/02_sponge_datasets/BLAST_RQUA
cd MasterThesis/02_sponge_datasets/BLAST_RQUA
#Create the BLAST databases
makeblastdb -in ../Emu_genome_v1.fa -dbtype nucl -out E_muelleri_genome
makeblastdb -in ../Ephydatia_muelleri_Trinity.fa -dbtype nucl -out E_muelleri_transcriptome
makeblastdb -in ../Eufr_Trinity.fa -dbtype nucl -out E_fragilis
makeblastdb -in ../Spongilla_lacustris_ar.fa -dbtype nucl -out S_lacustris_1
makeblastdb -in ../spongilla_transcriptome_annotated_200501.fasta -dbtype nucl -out S_lacustris_2
makeblastdb -in ../AV4_1_Trinity.fa -dbtype nucl -out A_vastus
makeblastdb -in ../Sycon_coactum_ar.fa -dbtype nucl -out S_coactum
#Run TBLASTN using the E. muelleri RQUA amino acid sequence retrieved from EukProt
ls *.ndb | while read file; do database=$(echo $file | cut -d "." -f 1); tblastn -query E_muelleri_RQUA.fasta -db $database -out $database\_RQUA.tblastn; done
```

<h3>Other Euglenozoa species</h3>

RQUA sequences belonging to other euglenid species that were received through personal communication are stored inside the file `other_euglenids_rqua_oneline.fasta`, that is located in MasterThesis/01_phylogenetic_analysis/Phylogenetic_tree (see below).

<h2>Phylogenetic tree of RQUA</h2>

```bash
#Directory used to store all files associated with generating the phylogenetic tree
mkdir MasterThesis/01_phylogenetic_analysis/Phylogenetic_tree
cd MasterThesis/01_phylogenetic_analysis/Phylogenetic_tree
```

<h3>Input file of RQUA homologues</h3>

The files that will be used to construct the phylogenetic tree are as follows:

- rqua.seqs.fasta (used in the phylogenetic reconstruction of Stairs et al., 2018)
- rqua_outgroups.fasta (used in the phylogenetic reconstruction of Stairs et al., 2018)
- eukprot_rqua.fasta
- eukprot3_rqua.fasta
- all_sponge_RQUA.fasta (created manually using the sponge RQUA sequences of _E. muelleri_, _S. lacustris_, and _E. fragilis_)
- other_euglenids_rqua_oneline.fasta
- Pelomyxa_schiedti_RQUA.fasta (dditional sequence belonging to the protist _Pelomyxa schiedti_: https://www.ncbi.nlm.nih.gov/protein/KAH3742582.1?report=genbank&log$=protalign&blast_rank=8&RID=3PUASV4K01N)

All files are placed inside the MasterThesis/01_phylogenetic_analysis/Phylogenetic_tree directory and concatenated.

```bash
cat *.fasta > all_rqua_final_iqtree.fasta
```

<h3>Multiple Sequence Alignment and trimming</h3>

In order to build a phylogenetic tree, the sequences are first aligned using MAFFT v7.490. Then the alignment is trimmed using BMGE v1.12.

```bash
#Install MAFFT
conda install mafft=7.490
#Align all RQUA sequences
mafft all_rqua_final_iqtree.fasta > all_rqua_final_iqtree_aligned.fasta
#Install BMGE
conda install bmge=1.12
#Remove ambiguously aligned residues
bmge -i all_rqua_final_iqtree_aligned.fasta -of all_rqua_final_iqtree_trimmed.fasta -t AA -m BLOSUM30 -h 0.7
```

<h3>Phylogenetic reconstruction</h3>

The phylogenetic reconstruction is performed using IQ-Tree v2.1.4_beta, under the LG+C60+F+G model of evolution, with 1000 ultrafast bootstrap.

```bash
#Install IQ-Tree
conda install iqtree=2.1.4_beta
#Create the phylogenetic tree
iqtree -quiet -s all_rqua_final_iqtree_trimmed.fasta -pre all_rqua_final -bb 1000 -nt 4 -alrt 1000 -m LG+C60+F+G
```

<h1>Proteome prediction from the transcriptomes of <i>Ephydatia muelleri</i>, <i>Eunapius fragilis</i>, and <i>Spongilla lacustris</i></h1>

Protein sequences are predicted using TransDecoder (v 5.5.0) (documentation can be found at https://github.com/TransDecoder/TransDecoder/wiki) from the transcriptomes of the three RQUA-encoding freshwater sponges.

First, the coding sequences are extracted from the transcriptomes. Then, protein sequences are predicted.

```bash
#Install TransDecoder
conda install transdecoder=5.5.0
#Run on E. muelleri transcriptome
TransDecoder.LongOrfs -t Ephydatia_muelleri_Trinity.fa #CDS
TransDecoder.Predict -t Ephydatia_muelleri_Trinity.fa #Proteins
#Run on E. fragilis transcriptome
TransDecoder.LongOrfs -t Eufr_Trinity.fa
TransDecoder.Predict -t Eufr_Trinity.fa
#Run on S. lacustris transcriptome 1
TransDecoder.LongOrfs -t Spongilla_lacustris_ar.fa
TransDecoder.Predict -t Spongilla_lacustris_ar.fa
#Run on S. lacustris transcriptome 2
TransDecoder.LongOrfs -t spongilla_transcriptome_annotated_200501.fasta
TransDecoder.Predict -t spongilla_transcriptome_annotated_200501.fasta
```

After running TransDecoder, the output files are organized into separate directories.

```bash
#Directory for storing the predicted coding sequence files
mkdir MasterThesis/02_sponge_datasets/Predicted_cds
cd MasterThesis/02_sponge_datasets/Predicted_cds
mv ../*.cds . #Move all .cds files here
#Directory for storing the predicted proteomes
mkdir MasterThesis/02_sponge_datasets/Predicted_proteins
cd MasterThesis/02_sponge_datasets/Predicted_proteins
mv ../*.pep . #Move all .pep files here
```

<h1>Investigation of potential contamination of the datasets of RQUA-encoding sponges</h1>

```bash
#New directory to store files related to this step of the analysis
mkdir MasterThesis/03_contamination_checks
cd MasterThesis/03_contamination_checks
```

<h2>Analysis of sponge predicted proteomes</h2>

The sequences in the sponge predicted proteomes are used as query in a BLASTX analysis against a modified version of EukProt. This database contains EukProt v2 (without the _Ephydatia muelleri_ sequences)and the predicted proteomes of the three sponges. Self-hits for each freshwater sponge species are excluded using the `-negarive_seqidlist` option of BLAST. In order to do so, the BLAST database is created using the `-parse_seqids` option, making the database easier to parse and manipulate.

Commands to create the database:

```bash
#Concatenate the files used to create the modified EukProt database into one dataset
cat ../01_phylogenetic_analysis/EukProt2/Indexing/all_indexed_eukprot_oneline.fasta ../02_sponge_datasets/Predicted_proteins/*.pep > eukprot_and_sponges.fasta
#Remove the E. muelleri sequences
awk -v A=-1 '/Ephydatia/ { A = 1 } A-- >= 0 {next } 1' eukprot_and_sponges.fasta > eukprot_and_sponges_no_ephy.fasta
#In order to create a BLAST database using the -parse_seqids flag, the headers need to be shorter than 50 characters (up to the first space)
#Run an in-house python script to solve this issue
python3 ../scripts/fix_local_ids.py eukprot_and_sponges.fasta eukprot_and_sponges_fixed.fasta
#Create the database using the -parse_seqids flag
makeblastdb -in eukprot_and_sponges_fixed.fasta -dbtype prot -parse_seqids -out eukprot_and_sponges
```

Commands to perform the BLASTX analyses:

```bash
#Create files of seqids to exclude for each species of sponge
cat ../02_sponge_datasets/Predicted_proteins/Eufr_Trinity.fa.transdecoder.pep | grep "^>" | cut -d " " -f 1 | tr -d ">" > E_fragilis_exclude_seqid.txt #Eunapius fragilis
cat ../02_sponge_datasets/Predicted_proteins/Ephydatia_muelleri_Trinity.fa.transdecoder.pep | grep "^>" | cut -d " " -f 1 | tr -d ">" > E_muelleri_exclude_seqid.txt #Ephydatia muelleri
#For Spongilla lacustris there are two datasets to exclude self-hits from. One of these datasets has had all of the headers indexed, therefore pull out the indexed headers and add them to the exclusion list
#To check that all of the headers are indexed in the database:
#All headers coming from the original file contain the string "Spongilla"
#The output of the command below is empty
cat eukprot_and_sponges_fixed.fasta | grep "Spongilla"
#Extract the S. lacustris headers that were indexed
cat ../02_sponge_datasets/Predicted_proteins/Spongilla_lacustris_ar.fa.transdecoder.pep | grep "^>" | while read line; do grep "$line" renamed_headers.tsv; done | cut -f 1 | tr -d ">" > S_lacustris_exclude_seqid.txt
#The headers in the other S. lacustris dataset did not get converted
#Command below does not return anything, proving that no sequence headers were converted
cat ../02_sponge_datasets/Predicted_proteins/spongilla_transcriptome_annotated_200501.fasta.transdecoder.pep | grep "^>" | while read line; do grep "$line" renamed_headers.tsv; done
#Extract the seqids and add them to the S. lacustris exclusion list
cat ../02_sponge_datasets/Predicted_proteins/spongilla_transcriptome_annotated_200501.fasta.transdecoder.pep | grep "^>" | cut -d " " -f 1 | tr -d ">" >> S_lacustris_exclude_seqid.txt
#Run the BLAST analyses
blastx -query ../02_sponge_datasets/Ephydatia_muelleri_Trinity.fa -db eukprot_and_sponges -out E_muelleri_transc_vs_EukProt.blastx -evalue 1e-5 -outfmt 6 -max_target_seqs 100 -negative_seqidlist E_muelleri_exclude_seqid.txt -num_threads 8
blastx -query ../02_sponge_datasets/Eufr_Trinity.fa -db eukprot_and_sponges -out E_fragilis_transc_vs_EukProt.blastx -evalue 1e-5 -outfmt 6 -max_target_seqs 100 -negative_seqidlist E_fragilis_exclude_seqid.txt -num_threads 8
blastx -query ../02_sponge_datasets/Spongilla_lacustris_ar.fa -db eukprot_and_sponges -out S_lacustris_transc_1_vs_EukProt.blastx -evalue 1e-5 -outfmt 6 -max_target_seqs 100 -negative_seqidlist S_lacustris_exclude_seqid.txt -num_threads 4
blastx -query ../02_sponge_datasets/spongilla_transcriptome_annotated_200501.fasta -db eukprot_and_sponges -out S_lacustris_transc_2_vs_EukProt.blastx -evalue 1e-5 -outfmt 6 -max_target_seqs 100 -negative_seqidlist S_lacustris_exclude_seqid.txt -num_threads 4
```

Then, the results of the BLASTX analyses are summarized by creating pie charts that show the taxonomic distribution of the top 10 hits of each query sequence. These hits are classified at phylum level. All phyla accounting for less than 1% of the top 10 hits is counted under "Other" in the final pie charts. Commands to create the pie charts using in-house created python scripts:

```bash
#Reformat the EukProt_included_datasets file for EukProt v2 so that it can be used by the scripts
cat ../01_phylogenetic_analysis/EukProt2/EukProt_included_data_sets.v02.2020_06_30.txt | grep "^EP" | cut -f 1,10 > EukProt_taxonomy.txt
#Parse the BLAST result files flagging all instances where hits belonging to Euglenozoa are found among the top 10 hits
python3 ../scripts/parse_blast.py -i E_muelleri_transc_vs_EukProt.blastx -o E_muelleri_transc_vs_EukProt_parsed.txt -e EukProt_taxonomy.txt -c Euglenozoa -n 10
python3 ../scripts/parse_blast.py -i E_fragilis_transc_vs_EukProt.blastx -o E_fragilis_transc_vs_EukProt_parsed.txt -e EukProt_taxonomy.txt -c Euglenozoa -n 10
python3 ../scripts/parse_blast.py -i S_lacustris_transc_1_vs_EukProt.blastx -o S_lacustris_transc_1_vs_EukProt_parsed.txt -e EukProt_taxonomy.txt -c Euglenozoa -n 10
python3 ../scripts/parse_blast.py -i S_lacustris_transc_2_vs_EukProt.blastx -o S_lacustris_transc_2_vs_EukProt_parsed.txt -e EukProt_taxonomy.txt -c Euglenozoa -n 10
#Create pie charts
#It is possible that some modules, such as matplotlib, seaborn, and ete3 need to be installed. This can be done using e.g. pip
pip install matplotlib
pip install seaborn
pip install ete3
#Create the actual pie charts grouping the top 10 hits at phylum level
#Note that this script will download the NCBI taxonomy database if it is not already present in the working directory if the -n flag is added
python3 ../scripts/plot_blast_results.py -i E_muelleri_transc_vs_EukProt_parsed.txt -o E_muelleri.pdf -r phylum -t 10
python3 ../scripts/plot_blast_results.py -i E_fragilis_transc_vs_EukProt_parsed.txt -o E_fragilis.pdf -r phylum -t 10
python3 ../scripts/plot_blast_results.py -i S_lacustris_transc_1_vs_EukProt_parsed.txt -o S_lacustris_1.pdf -r phylum -t 10
python3 ../scripts/plot_blast_results.py -i S_lacustris_transc_2_vs_EukProt_parsed.txt -o S_lacustris_2.pdf -r phylum -t 10
```

<h2>Check for contamination using <i>E. gracilis</i> 18S rRNA gene</h2>

The gene for 18S rRNA of _Euglena gracilis_ (GenBank: AF283309.1) is used as query in a BLAST analysis against the genome and transcriptomes of RQUA-encoding sponges.

```bash
#Set up new directory
mkdir MasterThesis/03_contamination_checks/18S
cd MasterThesis/03_contamination_checks/18S
```

Databaes of the datasets of interest were created already and can be found at MasterThesis/02_sponge_datasets_BLAST_RQUA. The query sequence is saved inside a file called `E_gracilis_18S.fasta`.

```bash
#Run the BLASTN analyses
blastn -task blastn -query E_gracilis_18S.fasta -db ../../02_sponge_datasets/BLAST_RQUA/E_muelleri_genome -out E_muelleri_genome_18S.blastn
blastn -task blastn -query E_gracilis_18S.fasta -db ../../02_sponge_datasets/BLAST_RQUA/E_muelleri_transcriptome -out E_muelleri_transcriptome_18S.blastn
blastn -task blastn -query E_gracilis_18S.fasta -db ../../02_sponge_datasets/BLAST_RQUA/E_fragilis -out E_fragilis_18S.blastn
blastn -task blastn -query E_gracilis_18S.fasta -db ../../02_sponge_datasets/BLAST_RQUA/S_lacustris_1 -out S_lacustris_1_18S.blastn
blastn -task blastn -query E_gracilis_18S.fasta -db ../../02_sponge_datasets/BLAST_RQUA/S_lacustris_2 -out S_lacustris_2_18S.blastn
```

<h2>Check for contamination using <i>E. gracilis</i> plastidial 16S rRNA gene</h2>

```bash
#Set up new directory
mkdir MasterThesis/03_contamination_checks/16S
cd MasterThesis/03_contamination_checks/16S
```

The gene for plastidial 16S rRNA of _Euglena gracilis_ (GenBank: V00159.1) is used as query in a BLAST analysis against the genome and transcriptomes of RQUA-encoding sponges.

Databaes of the datasets of interest were created already and can be found at MasterThesis/02_sponge_datasets_BLAST_RQUA. The query sequence is saved inside a file called `E_gracilis_16S.fasta`.

```bash
#Run the BLASTN analyses
blastn -task blastn -query E_gracilis_16S.fasta -db ../../02_sponge_datasets/BLAST_RQUA/E_muelleri_genome -out E_muelleri_genome_16S.blastn
blastn -task blastn -query E_gracilis_16S.fasta -db ../../02_sponge_datasets/BLAST_RQUA/E_muelleri_transcriptome -out E_muelleri_transcriptome_16S.blastn
blastn -task blastn -query E_gracilis_16S.fasta -db ../../02_sponge_datasets/BLAST_RQUA/E_fragilis -out E_fragilis_16S.blastn
blastn -task blastn -query E_gracilis_16S.fasta -db ../../02_sponge_datasets/BLAST_RQUA/S_lacustris_1 -out S_lacustris_1_16S.blastn
blastn -task blastn -query E_gracilis_16S.fasta -db ../../02_sponge_datasets/BLAST_RQUA/S_lacustris_2 -out S_lacustris_2_16S.blastn
```

<h2>Presence of spliced leaders</h2>

The presence of spliced leaders is investigated in _rqua_ sequences present in the transcriptomes of euglenid species _Apiculatamorpha spiralis_, _Euglena gracilis_, _Eutreptiella gymnastica_, _Euglena longa_, and _Euglena mutabilis_. The RQUA amino acid sequences can be obtained from the fasta file used to build the phylogenetic tree. The _rqua_ transcripts are retrieved from the transcriptomes used to generate the proteomes on EukProt (links are available in the EukProt_included_datasets file).

```bash
#Set up new directory for this analysis
mkdir MasterThesis/03_contamination_checks/Spliced_leaders
cd MasterThesis/03_contamination_checks/Spliced_leaders
#Transcriptomes of the euglenids are downloaded from the links indicated by EukProt
#Make the databases of the euglenid transcriptomes
makeblastdb -in A_spiralis_transcriptome.fasta -dbtype nucl -out A_spiralis
makeblastdb -in E_gracilis_transcriptome.fasta -dbtype nucl -out E_gracilis
makeblastdb -in E_gymnastica_transcriptome.fasta -dbtype nucl -out E_gymnastica
makeblastdb -in E_longa_transcriptome.fasta -dbtype nucl -out E_longa
makeblastdb -in E_mutabilis_transcriptome.fasta -dbtype nucl -out E_mutabilis
#Run TBLASTN using the corresponding RQUA sequences
tblastn -query A_spiralis_rqua.fasta -db A_spiralis -out A_spiralis.tblastn
tblastn -query E_gracilis_rqua.fasta -db E_gracilis -out E_gracilis.tblastn
tblastn -query E_gymnastica_rqua.fasta -db E_gymnastica -out E_gymnastica.tblastn
tblastn -query E_longa_rqua.fasta -db E_longa -out E_longa.tblastn
tblastn -query E_mutabilis_rqua.fasta -db E_mutabilis -out E_mutabilis.tblastn
```

<h1>Analysis of <i>Ephydatia muelleri</i> genomic <i>rqua</i> sequence</h1>

<h2>Analysis of neighboring genes of <i>rqua</i></h2>

```bash
#Set up new directory
mkdir MasterThesis/04_genomic_rqua/Neighboring_genes
cd MasterThesis/04_genomic_rqua/Neighboring_genes
```

Nucleotide sequences of neighboring genes of _rqua_ on the _E. muelleri_ chromosome are retrieved manually from the CDS file available at https://bitbucket.org/EphydatiaGenome/ephydatiagenome/downloads/Emu_genome_v1.codingseq.gz. The FASTA file containing these sequences (`rqua_and_neighboring_genes.fasta`) is placed inside the current working directory.

The aforementioned sequences are used as queries gainst a combined dataset consisting of EukProt v2, a pre-release subset of EukProt v3, and predicted proteomes of _S. lacustris_ and _E. fragilis_.

<h3>Creating the dataset to use as BLAST database</h3>

```bash
#Concatenate files to generate a combined dataset
cat ../../01_phylogenetic_analysis/EukProt2/Indexing/all_indexed_eukprot_oneline.fasta ../../01_phylogenetic_analysis/EukProt3/proteins/all_eukprot3_oneline.fasta ../../02_sponge_datasets/Predicted_proteins/Eufr_Trinity.fa.transdecoder.pep ../../02_sponge_datasets/Predicted_proteins/Spongilla_lacustris_ar.fa.transdecoder.pep ../../02_sponge_datasets/Predicted_proteins/spongilla_transcriptome_annotated_200501.fasta.transdecoder.pep > eukprot_and_sponges.fasta
#Remove the E. muelleri sequences
awk -v A=-1 '/Ephydatia/ { A = 1 } A-- >= 0 {next } 1' eukprot_and_sponges.fasta > eukprot_and_sponges_no_ephy.fasta
```

<h3>Run BLAST on the combined dataset</h3>

```bash
#Create the blast database
makeblastdb -in eukprot_and_sponges_no_ephy.fasta -dbtype prot -out eukprot_and_sponges
#Run BLASTX
blastx -query rqua_and_neighboring_genes.fasta -db eukprot_and_sponges -evalue 1e-5 -out rqua_and_neighboring_genes.blastx
```

<h2>Investigation of PacBio reads from <i>E. muelleri</i></h2>

```bash
#Set up new directory
mkdir MasterThesis/04_genomic_rqua/PacBio
cd MasterThesis/04_genomic_rqua/PacBio
```

<h3>Download of SRA datasets and BLASTN</h3>

The PacBio sequencing datasets from _E. muelleri_ are first downloaded using SRA Toolkit (v 2.8.0) (https://github.com/ncbi/sra-tools), specifically, the `fastq-dump` program that is included in this toolkit. Running this program requires a list of SRA accession numbers. These can be found at https://www.ncbi.nlm.nih.gov/sra?linkname=bioproject_sra_all&from_uid=579531. The accession numbers of interest are:

- SRR10983238
- SRR10983239
- SRR10983240
- SRR10983241
- SRR10983242

These are added to a file called `accessions_ephydatia.txt`.

```bash
#First install SRA Toolkit
conda install sra-tools=2.8.0
#Then iterate through the list of accession numbers and download gzip FASTA files of each
dos2unix accessions_ephydatia.txt #Run if there is an encoding error with the file
for line in $(cat accessions_ephydatia.txt); do fastq-dump --fasta 0 --gzip $line; done
```

To search for long sequencing reads containing the _rqua_ nucleotide sequence, the SRA datasets are first converted into BLAST databases. Then BLASTN is run. The query sequence is the _Ephydatia muelleri rqua_ sequence.

```bash
#Create BLAST databases of all SRA datasets
ls *.gz | while read file; do database=$(echo $file | cut -d "." -f 1); gunzip -c $file | makeblastdb -in - -dbtype nucl -out $database -title $database; done
#Run BLASTN
ls *.gz | while read file; do database=$(echo $file | cut -d "." -f 1); blastn -task blastn -query Ephydatia_RQUA_nt.fasta -db $database -evalue 1e-5 -out $database.blastn; done
```

Finally, the reads matching _rqua_ are extracted.

```bash
cat *.blastn | grep "^>" | while read line; do gunzip -c *.gz | grep -A 1 "$line" >> all_rqua_hits_blastn.fasta; done
```

<h3>Mapping of <i>rqua</i> hits to the reference genome</h3>

Two software are used to map the long sequencing reads matching _rqua_ to the _E. muelleri_ reference genome: MashMap v2.0 (https://github.com/marbl/MashMap) and minimap2 v2.24 (https://github.com/lh3/minimap2). Both programs can be installed easily using conda.

```bash
conda install mashmap=2.0
conda install minimap2=2=24
```

The programs are then run using the FASTA files containing the long sequencing reads matching _rqua_. First MashMap:

```bash
mashmap -r Ephydatia_muelleri_genome.fa -q all_rqua_hits_blastn.fasta
mv mashmap.out mashmap_blastn.out #Rename the output file
#This output file can be used to extract information about the mapping, but cannot be opened by most genome browser software
#Instead, MashMap comes with a script that can generate plots of the result files
git clone https://github.com/marbl/MashMap #Download the GitHub repository
#The script is found inside the MashMap/scripts directory and requires the gnuplot dependency to run
conda install gnuplot #Install gnuplot
#Run the plotting script
perl generateDotPlot png small ../../mashmap_blastn.out
mv out.png blastn_out.png #Rename the output file
```

Then minimap2, to generate a SAM file than can be converted to the BAM format and easily visualized using genome browser software such as GenomeView.

```bash
#Generate a PAF file that contains easy to parse information about the alignment
minimap2 -o all_rqua_hits_minimap.paf Ephydatia_muelleri_genome.fa all_rqua_hits_blastn.fasta
#Generate a SAM file that can then be converted into BAM format and visualized in GenomeView
minimap2 -a Ephydatia_muelleri_genome.fa all_rqua_hits_blastn.fasta > all_rqua_hits_blastn.sam
#Convert the SAM format output to BAM format using samtools
conda install samtools #Install samtools
#Convert SAM format to BAM format
samtools view -S -b all_rqua_hits_blastn.sam > all_rqua_hits_blastn.bam
#Sort resulting BAM file
samtools sort all_rqua_hits_blastn.bam > all_rqua_hits_blastn.sorted.bam
#Index the sorted BAM file
samtools index all_rqua_hits_blastn.sorted.bam
```

<h1>Identification of CI and CII subunits in RQUA-encoding sponges</h1>

Quinone-interacting components of CI and CII are identified in the RQUA-encoding sponges using BLAST and HMM profiles. First, a query file is created that contains amino acid sequences of the various subunits (SDHA, SDHB, SHDC, SDHD, NDUFS2, NDUFS7, ND1); the sequences used as queries belong to _Homo sapiens_, _Saccharomyces cerevisiae_, and _Amphimedon queenslandica_ and are retrieved from Uniprot, NCBI, or by using the web version of BLAST to look for _A. queenslandica_ sequences (human sequences were used as queries). HMM profiles are retrieved from NCBI. All query sequences are placed in the file `ETC_components.fasta`, which is used as query file.

```bash
#Directory to store all files associated with this step of the analysis
mkdir MasterThesis/05_ETC_components
cd MasterThesis/05_ETC_components
```

<h2>Using BLAST to search for the ETC components</h2>

First, BLAST databases of all the predicted protein files, as well as the coding sequence files are created.

```bash
#Ephydatia muelleri cds
makeblastdb -in ../02_sponge_datasets_Predicted_cds/Ephydatia_muelleri_Trinity.fa.transdecoder.cds -dbtype nucl -out E_muelleri_cds
#Eunapius fragilis cds
makeblastdb -in ../02_sponge_datasets_Predicted_cds/Eufr_Trinity.fa.transdecoder.cds -dbtype nucl -out E_fragilis_cds
#Spongilla lacustris cds 1
makeblastdb -in ../02_sponge_datasets_Predicted_cds/Spongilla_lacustris_ar.fa.transdecoder.cds -dbtype nucl -out S_lacustris_cds_1
#Spongilla lacustris cds 2
makeblastdb -in ../02_sponge_datasets_Predicted_cds/spongilla_transcriptome_annotated_200501.fasta.transdecoder.cds -dbtype nucl -out S_lacustris_cds_2
#Ephydatia muelleri proteins
makeblastdb -in ../02_sponge_datasets_Predicted_proteins/Ephydatia_muelleri_Trinity.fa.transdecoder.pep -dbtype prot -out E_muelleri_pep
#Eunapius fragilis proteins
makeblastdb -in ../02_sponge_datasets_Predicted_proteins/Eufr_Trinity.fa.transdecoder.pep -dbtype prot -out E_fragilis_pep
#Spongilla lacustris proteins 1
makeblastdb -in ../02_sponge_datasets_Predicted_proteins/Spongilla_lacustris_ar.fa.transdecoder.pep -dbtype prot -out S_lacustris_pep_1
#Spongilla lacustris proteins 2
makeblastdb -in ../02_sponge_datasets_Predicted_proteins/spongilla_transcriptome_annotated_200501.fasta.transdecoder.pep -dbtype prot -out S_lacustris_pep_2
```

Then, TBLASTN is run on the .cds databases and BLASTP on the .pep databases using the `ETC_components.fasta` file as query.

```bash
#For the .cds databases
ls *.ndb | while read line; do database=$(echo $line | cut -d . -f 1); tblastn -query ETC_components.fasta -db $database -out $database.tblastn; done
#For the .pep databases
ls *.pdb | while read line; do database=$(echo $line | cut -d . -f 1); blastp -query ETC_components.fasta -db $database -out $database.blastp; done
```

Subunit ND1 is encoded on the mitochondrial genome. Therefore, available mitochondrial genomes are used in a TBLASTN analysis; the query file, `ND1_query.fasta` contains ND1 sequences from _Homo sapiens_ and _Amphimedon queenslandica_. The sponge mitochondrial genomes are available from:

- _Ephydatia muelleri_: LT158504.1 (https://www.ncbi.nlm.nih.gov/nuccore/LT158504.1)
- _Spongilla lacustris_: LT158503.1 (https://www.ncbi.nlm.nih.gov/nuccore/LT158503.1)

First, the mitochondrial genomes are converted into BLAST databases. Then, the TBLASTN analysis is run.

```bash
#Create databases of the downloaded mitochondrial genomes
makeblastdb -in E_muelleri_mt_genome.fasta -dbtype nucl -out E_muelleri_mt
makeblastdb -in S_lacustris_mt_genome.fasta -dbtype nucl -out S_lacustris_mt
#Run TBLASTN using the ND1 sequences as query
tblastn -query ND1_query.fasta -db E_muelleri_mt -out E_muelleri_mt.tblastn
tblastn -query ND1_query.fasta -db S_lacustris_mt -out S_lacustris_mt.tblastn
```

<h2>Using HMMER and HMM profiles to look for ETC components</h2>

To identify the subunits of Complex II, HMM profiles are also used. These are downloaded from NCBI and queried using the program HMMER (v3.3.2).

```bash
#Install HMMER
conda install hmmer=3.3.2
#Make a new directory to store all files associated with hmmsearch
mkdir MasterThesis/05_ETC_components/HMM
cd MasterThesis/05_ETC_components/HMM
#Use wget to download the HMM profiles into this directory
#Use hmmsearch to look for CII subunits in the predicted proteomes of the freshwater sponges
ls *.HMM | while read line; do profile=${line::-4}; ls ../../02_sponge_datasets/Predicted_proteins/*.pep | while read line2; do proteome=$(basename $line2); hmmsearch $line $line2 > $profile"_"$proteome.txt; done; done
```

<h2>Using BUSCO to assess the level of completion of the sponge transcriptomes</h2>

The degree of completion of the sponge transcriptomes is assessed using BUSCO v5.0.0 (https://busco.ezlab.org/; user guide at https://busco.ezlab.org/busco_userguide.html#running-busco) to determine how many orthologs can be found in the datasets out of the total number of core proteins expected to be found in metazoans.

BUSCO is run on all the transcriptomes belonging to RQUA-encoding sponges, as well as the predicted proteome of _Amphimedon queenslandica_ (available as part of EukProt) as reference, in order to assess the level of completion that is expected of sponges when using BUSCO.

```bash
#New directory to store BUSCO results
mkdir MasterThesis/05_ETC_components/BUSCO
cd MasterThesis/05_ETC_components/BUSCO
#BUSCO can be installed using conda
conda install busco=5.0.0
#Check the lineages available for the BUSCO analysis
busco --list-datasets
#Metazoa will be used as reference lineage
#Run BUSCO for E. muelleri
busco -m transcriptome -i ../../02_sponge_datasets/Ephydatia_muelleri_Trinity.fa -o
BUSCO_E_muelleri -l metazoa
#Run BUSCO for E. fragilis
busco -m transcriptome -i ../../02_sponge_datasets/Eufr_Trinity.fa -o BUSCO_E_fragilis -l metazoa
#Run BUSCO for the two S. lacustris transcriptomes
busco -m transcriptome -i ../../02_sponge_datasets/Spongilla_lacustris_ar.fa -o BUSCO_S_lacustris_1 -l metazoa
busco -m transcriptome -i ../../02_sponge_datasets/spongilla_transcriptome_annotated_200501.fasta -o BUSCO_S_lacustris_2 -l metazoa #Returns an error
#Run BUSCO for the A. queenslandica predicted proteome
busco -m proteins -i ../../01_phylogenetic_analysis/EukProt2/proteins/EP00119_Amphimedon_queenslandica.fasta -o BUSCO_A_queenslandica -l metazoa
```
