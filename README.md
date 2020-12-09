# SARS-CoV-2 Immunoanalytics Platform
## http://genomics.lshtm.ac.uk/immuno/

In this repo you can find the most up-to-date raw data files for the data and analyses shown genomics.lshtm.ac.uk/immuno website. The code underpinning the website, including visualisations, static data, and tables in JSON format can be found at [Jody Phelan's COVID-profiler repo](https://github.com/jodyphelan/covid-profiler).

### Immunoanalytics Plot

**Uniprot annotations** - these are feature annotations scraped form the [UniProt database](https://covid-19.uniprot.org/uniprotkb?query=*) . These have been parsed in to a tabular format and mapped to the relevant proteome position.

**Non-synonymous mutation data** - This dataset was generated using the [dataset_builder.sh](./scripts/database_builder.sh) script. We used an in-house script designed by [Matt Higgins](https://github.com/MatthewHiggins2017) to count non-synonymous mutations in aligned FASTA formatted amino-acid sequence files.

**Epitope mapping** - B/T cell epitope inference was mapped to the main data table using positional metadata. The IEDB database was scraped using a looped *wget* command, each epitope sequence was mapped using BLASTp and merged with the primary data table.

**Homology analysis** - We obtained a canonical amino-acid sequence of each orthologous protein, separated it in to kmers and mapped it to the SARS-CoV-2 reference sequence using BLASTp.
