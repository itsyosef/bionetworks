# Project Notes

**Should be able to access these files on Rivanna at /project/biocomplexity/jho5ze

## BGCs

**Gist**:
* Working with [my fork of Andrew's pangenome_graph pipeline](https://github.com/itsyosef/pangenome_graphs),  the scripts in scripts/ generally try to achieve the task of annotating these graphs with BGCs from sources like antismash, get features like the associated pgfams, pfams and BGC flag for intervals on the graph, and build data sets for training a classifier on these splits of the graph.
* There is also code for identifying superbubbles in scripts/superbubble.py

## Coevolution

**Gist:**
* Mainly only has the days till detection maps for VDH
* Also has the example notebook for MSAFroMetadata used to help Alan subsample the GISAID data
* The script in scripts/ can be used to subsample the MSA from GISAID to only include those from the USA and generate a screed database for easy access and for use in the entropy calculation

## Covid

**Gist:**
* Has all the notebooks for the entropy calculation and some of the initial forays into the nextstrain vizualizations. 
* scripts/ has all the scripts used for the msa entropy graph generation and calculation, especially run_pipeline.sh
* scripts/process_msa_headers.sh can be used to subsample the MSA from GISAID to only include those from the USA and generate a screed database for easy access and for use in the entropy calculation

## Epihiper

**Gist:**
* Has all the notebooks and scripts for dealing with the vulnerability calculations (computed in scripts/node_profiles.py), the heatmaps for V2 and V1 prevalence by different cuts in the population, and map plotting code. 

## Operons

**Gist:**
* scripts/coli.sh has the pipeline very roughly written out as I imagined it, but there will need to be a lot of cleaning of the input to make sure we are getting the data we need and we should verify the output is as expected as well. 
* SRA/scripts/pipeline.sh has a more fleshed out/automated version of the pipeline that attempts to do the data pull/prep given an input file which is a CSV from pysradb pull. It attempts to grab the relevant SRA runs and råun the operon pipeline from annogesic on it, but is not in a stable state
* There are also some scripts in compare_tools/ that were trying to compare the output/process of different tools that could get to operon annotation, but most of this isn't useful
