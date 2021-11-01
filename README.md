Peng et al., 2018 https://doi.org/10.1186/s40168-018-0546-9
============================================================

Below, details are given with respect to steps of sequence data processing that involved custom _python_ scripts. The individual subsections are named to match the corresponding methods sections of the publication and are complementary to those.

All used custom _python_ scripts were written using _python_ (v.2.7.12). Required packages include: 
* biopython
* textwrap
* sqlite3

All of these can be installed via the python package installer `pip`.

**1. Bioinformatic processing of total RNA-derived sequences**
----------------------------------------------------------

OTU clustering of rRNA-derived sequences was done using _usearch_ (v.7) (Edgar, 2010) applying the implemented _uclust_ algorithm. A custom _python_ script was subsequently used to convert the resulting OTU table (in fact a usearch cluster format file) into a legacy OTU table matching the specifications of _qiime_ (v.1.9.1) (Caporaso et al., 2010), which was in turn converted into a _biom_ table.

The original _python_ script (uc2otutab.py) by Robert Edgar is publicly available and its usage is in detail explained here: https://drive5.com/python/. We modified the script to match our internal naming conventions for datasets, and the script (uc2otutab_mod.py) is availabe from this repository.

The usage is as follows:
```
python uc2otutab_mod.py input_uc.file > output_otutab.file
```

**2. Assembly of full-length 16S rRNA sequences**
---------------------------------------------

We used _emirge_ (Miller et al., 2011) for the reconstruction of full-length 16S rRNA sequences. The subsampling of our 16S rRNA-derived sequence data for sequences of interest was a necessary prerequisite to reduce the computational load of 16S rRNA sequence reconstruction. For the subsampling of the data we wrote a _python_ script that processes various _qiime_ (v.1.9.1) output files. The script is available from this repository (extract_seqs_based_on_taxonomy.py).

The paths of the necessary input and the to be generated output files have to be modified in the code of the script. 

After modifying the path variables, as well as defining the taxonomic group of interest the script can be called as follows:
```
python extract_seqs_based_on_taxonomy.py
```

**3. Analysis of CAZyme-related sequences**
---------------------------------------

CAZyme-affiliated mRNA reads were identified by querying mRNA-derived sequences against a local copy of the dbCAN database (Yin et al., 2012) that has been indexed for being used with _diamond_ (Buchfink et al., 2015). _diamond_ was used as follows:

```
diamond blastx -d dbCAN_072017.dmd -q mRNA_derived_seqs.fna -o dbCAN_hits.out -f 6
```

Functional CAZyme modules were defined as outlined the manuscript main text. Annotation summaries for the mRNA-derived sequences queried against dbCAN were generated using the custom _python_ script dbCAN_annotator.py.

dbCAN_annotator makes use of a mapping file provided by the dbCAN consortium (http://csbl.bmb.uga.edu/dbCAN/download.php, CAZyDB-ec-info.txt.07-20-2017). This mapping file contains information about the CAZyme-family affiliation of each deposited sequence, as well as information about assigned enzyme commission numbers. The latter were used for defining the aforementioned CAZyme functional modules.

The mapping file was used for generating an indexed sqlite3 database object to facilitate fast querying, and thus a time-efficient processing of generated _diamond_ output.

In order to use dbCAN_annotator, the script (dbcan_annotator.py) and the sqlite3 database (dbcan.db) have to be downloaded. If the database is stored at a different place than the script, the path variable referring to the database has to be modified in dbCAN_annotator.py.

dbCAN_annotator can be used as follows.

### (i) Generating annotation summaries
```
python dbCAN_annotator.py --mode annot --diamond diamond_output.tab --seq queried_sequence_data.fasta 
```
### (ii) Extracting sequences linked to functional modules
```
python dbCAN_annotator.py --mode filter --diamond diamond_output.tab --seq queried_sequence_data.fasta 
```
### (iii) Extracting sequences linked to particular CAZyme functions of interest
```
python dbCAN_annotator.py --mode func_filter --diamond diamond_output.tab --seq queried_sequence_data.fasta --func list_of_CAZyme_functions_of_interest.txt
```

(c) Peng et al., 2018

**References**
--------------

* Edgar,RC Search and clustering orders of magnitude faster than BLAST, Bioinformatics 2010;26(19):2460-2461.

* Caporaso JG, Kuczynski J, Stombaugh J, Bittinger K, Bushman FD, Costello EK, Fierer N, Pena AG, Goodrich JK, Gordon JI, Huttley GA, Kelley ST, Knights D, Koenig JE, Ley RE, Lozupone CA, McDonald D, Muegge BD, Pirrung M, Reeder J, Sevinsky JR,
Tumbaugh PJ, Walters WA, Widmann J, Yatsunenko T, Zaneveld J, Knight R. QIIME allows analysis of high-throughput community sequencing data. Nat Methods. 2010;7(5):335-336.

* Miller CS, Baker BJ, Thomas BC, Singer SW, Banfield JF. EMIRGE: reconstruction of full-length ribosomal genes from microbial community short read sequencing data.Genome Biol. 2011;12(5).

* Yin YB, Mao XZ, Yang JC, Chen X, Mao FL, Xu Y. dbCAN: a web resource for automated carbohydrate-active enzyme annotation. Nucleic Acids Res. 2012;40(W1):789.

* Buchfink B, Xie C, Huson DH. Fast and sensitive protein alignment using DIAMOND.Nat Methods. 2015;12(1):59-60.
