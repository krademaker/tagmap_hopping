# Mapping of (putative) transposon insertion sites.
#
# This script maps (putative) insertion sites with the tagmapAnalyseR package
# for R, which can be used in various experimental contexts to investigate the
# functional impact of transposon hopping.
#
# This script is incorporated in the Snakemake pipeline for TagMap-based mapping
# of transposon hopping by the Van Steensel lab at the Netherlands Cancer
# Institute (NKI) at https://github.com/vansteensellab/tagmap_hopping/.
#
# Rule 'find_putative_insertions' precedes rule 'map_insertions' in which this
# script is called. Rule 'find_putative_insertions' stores putative insertions
# in the output folder 'putative_insertions', with detailed paired-end R2 read
# data.
# The output of 'find_putative_insertions', together with the raw BAM files are
# used as input for this script in rule 'map_insertions', which will map:
# insertion sites, their genomic (region) location, overhang sequence, strand,
# read count(s) and mapping quality (median MAPQ).
#
# Dependencies: data.table,
#               argparse,
#               tagmapAnalyseR (https://github.com/krademaker/tagmapAnalyseR/)
#
# Author:       Koen Rademaker (k.rademaker@nki.nl)
# Date:         6 December 2021
# -------------------------------------------------------------------------



# Initialise script -------------------------------------------------------
library(data.table)
library(argparse)
library(tagmapAnalyseR)

# Parse command line arguments --------------------------------------------
parser <- argparse::ArgumentParser(description = 'Map insertion sites using the tagmapAnalyseR package.')
parser$add_argument('--input', help = 'Putative integration sites reported by the pipeline script "find_insertions.sh".')
parser$add_argument('--bam', help = 'Sorted BAM files with both forward and reverse TagMap reads.')
parser$add_argument('--overhang', help = 'Expected overhang of the integration site (e.g. TTAA, TA).')
parser$add_argument('--depth', help = 'Minimum read depth of integration sites.')
parser$add_argument('--samtools', help = 'Path to samtools executable.')
parser$add_argument('--out', help = 'File for result output.')
argv <- parser$parse_args()

# Split BAM path into separate items for forward and/or reverse BAM files.
argv$bam <- unlist(strsplit(argv$bam, split = ' '))


# Map insertion sites -----------------------------------------------------
dt <- tagmapAnalyseR::readPutativeInsertions(argv$input)
ambiguous <- tagmapAnalyseR::findAmbiguousInsertionSites(dt, padding = 2)
mapped <- tagmapAnalyseR::mapInsertionSites(dt = dt,
                                            bam = argv$bam,
                                            overhang = argv$overhang,
                                            samtoolsPath = argv$samtools,
                                            depth = argv$depth,
                                            ambiguousInsertions = ambiguous)
# TODO: Store 'mapped' as .RDS for potential detailed investigation of the
# individual reads underlying insertion sites.
mapped <- within(mapped, rm(read_names, width))


# Write output file -------------------------------------------------------
data.table::fwrite(mapped, file = argv$out, sep = '\t', na = '.', quote = FALSE)
