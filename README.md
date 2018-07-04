Peng et al., 2018 doi:
======================

Below, details are given with respect to steps of sequence data processing that involved custom _python_ scripts. The individual subsections are named to match the corresponding methods sections of the publication and are complementary to those.

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

We used _emirge_ (REF) for the reconstruction of full-length 16S rRNA sequences. The subsampling of our 16S rRNA-derived sequence data for sequences of interest was a necessary prerequisite to reduce the computational load of 16S rRNA sequence reconstruction. For the subsampling of the data we wrote a _python_ script that processes various _qiime_ (v.1.9.1) output files. 

**3. Analysis of CAZyme-related sequences**
---------------------------------------

UVW
