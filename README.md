Peng et al., 2018 doi:
======================

Below, details are given with respect to steps of sequence data processing that involved custom _python_ scripts. The individual subsections are named to match the corresponding methods sections of the publication and are complementary to those.

All used custom _python_ scripts were written using _python_ (v.2.7.12). Required packages include: 
* biopython
* textwrap
* sqlite3

All of these can be installed via the python package installer `pip`.

**1. Bioinformatic processing of total RNA-derived sequences**
----------------------------------------------------------

OTU clustering of rRNA-derived sequences was done using _usearch_ (v.7) (REF) applying the implemented _uclust_. A custom _python_ script was subsequently used to convert the resulting OTU table (in fact a usearch cluster format file) into a legacy OTU table matching the specifications of _qiime_ (v.1.9.1) (REF), which was in turn converted into a _biom_ table.

The original _python_ script (uc2otutab.py) by Robert Edgar is publicly available and its usage is in detail explained here: https://drive5.com/python/. We modified the script to match our internal naming conventions for datasets, and the script (uc2otutab_mod.py) is availabe from this repository.

The usage is as follows:
```
python uc2otutab_mod.py input_uc.file > output_otutab.file
```

**2. Assembly of full-length 16S rRNA sequences**
---------------------------------------------

We used _emirge_ (REF) for the reconstruction of full-length 16S rRNA sequences. The subsampling of our 16S rRNA-derived sequence data for sequences of interest was a necessary prerequisite to reduce the computational load of 16S rRNA sequence reconstruction. For the subsampling of the data we wrote a _python_ script that processes various _qiime_ (v.1.9.1) output files. The script is available from this repository (extract_seqs_based_on_taxonomy.py).

In its current version the paths of the necessary input and the to be generated output files have to be modified in the code of the script. A more convenient version, where paths and parameters will be passed via commandline arguments will be available shortly.

After modifying the path variables, as well as defining the taxonomic group of interest the script can be called as follows:
```
python extract_seqs_based_on_taxonomy.py
```

**3. Analysis of CAZyme-related sequences**
---------------------------------------

CAZyme-affiliated mRNA reads were identified by querying mRNA-derived sequences against a local copy of the dbCAN database (Ref) that has been indexed for being used with _diamond_ (Ref). _diamond_ was used as follows:

```
diamond blastx -d dbCAN_072017.dmd -q mRNA_derived_seqs.fna -o dbCAN_hits.out -f 6
```

Functional CAZyme modules were defined as outlined the manuscript main text. Annotations summaries for the mRNA-derived sequences queried against dbCAN were generated using the custom _python_ script dbCAN_annotator.py.

dbCAN_annotator makes use of a mapping file provided by the dbCAN consortium (http://csbl.bmb.uga.edu/dbCAN/download.php, CAZyDB-ec-info.txt.07-20-2017). This mapping file contains information about the CAZyme-family affiliation of each deposited sequence, as well as information about assigned enzyme commission numbers. The latter were used for defining the aforementioned CAZyme functional modules.

The mapping file was used for generating an indexed sqlite3 database object to facilitate fast querying, and thus a time-efficient processing of generated _diamond_ output.

In order to use dbCAN_annotator, the script (dbcan_annotator.py) and the sqlite3 database (dbcan.db) have to be downloaded. If the database is stored at a different place than the script, the path variable referring to the database has to be modified in dbCAN_annotator.py.

dbCAN_annotator can be used as follows.

#3 (i) Generating annotation summaries
```
python dbCAN_annotator.py --mode annot --diamond diamond_output.tab --seq queried_sequence_data.fasta 
```
#3 (ii) Extracting sequences linked to functional modules
```
python dbCAN_annotator.py --mode filter --diamond diamond_output.tab --seq queried_sequence_data.fasta 
```
#3 (iii) Extracting sequences linked to particular CAZyme functions of interest
```
python dbCAN_annotator.py --mode func_filter --diamond diamond_output.tab --seq queried_sequence_data.fasta --func list_of_CAZyme_functions_of_interest.txt
```

(c) Peng et al., 2018
