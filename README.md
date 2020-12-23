# SARS-CoV-2 Immunoanalytics Platform
### Website Resource: http://genomics.lshtm.ac.uk/immuno/
### Preprint Paper: https://www.biorxiv.org/content/10.1101/2020.05.11.089409v2

In this repo you can find the most up-to-date raw data files for the data and analyses shown genomics.lshtm.ac.uk/immuno website. The code underpinning the website, including visualisations, static data, and tables in JSON format can be found at [Jody Phelan's COVID-profiler repo](https://github.com/jodyphelan/covid-profiler).

  
### Raw data

All raw data files can be downloaded using this [link](https://lshtm-my.sharepoint.com/:f:/g/personal/lsh1603403_lshtm_ac_uk/Es5YMHN19nlGkBB0zR6Y1o8BrEVDMXwywSqZouQClT9cyg?e=V1wPT9). Because these files are too large to be hosted through GitHub we will be using Microsoft OneDrive for the foreseeable future.

##### If you would like a specific dataset that isn't easily accessed here or on the website, email daniel.ward1@lshtm.ac.uk

  
  
  
  
### Immunoanalytics Plot
#### as noted above, all raw data files can be found on OneDrive with the following [link](https://lshtm-my.sharepoint.com/:f:/g/personal/lsh1603403_lshtm_ac_uk/Es5YMHN19nlGkBB0zR6Y1o8BrEVDMXwywSqZouQClT9cyg?e=V1wPT9).

**Uniprot annotations** - these are feature annotations scraped form the [UniProt database](https://covid-19.uniprot.org/uniprotkb?query=*) . These have been parsed in to a tabular format and mapped to the relevant proteome position.

**Non-synonymous mutation data** - This dataset was generated using the [dataset_builder.sh](scripts/database_builder.sh) script. We used an in-house script designed by [Matt Higgins](https://github.com/MatthewHiggins2017) to count non-synonymous mutations in aligned FASTA formatted amino-acid sequence files.

**Epitope mapping** - B/T cell epitope inference was mapped to the main data table using positional metadata. The IEDB database was scraped using a looped *wget* command, each epitope sequence was mapped using BLASTp and merged with the primary data table.

**Homology analysis** - We obtained a canonical amino-acid sequence of each orthologous protein, separated it in to kmers and mapped it to the SARS-CoV-2 reference sequence using BLASTp.

  
  
  
### Mutation Tracker

The source for this page can be found here. This page was build using Google's Charts JS libraries. 
