#!/bin/bash
#
# 04_Claident_Lib2018i.sh
# A shell script to assign taxonomy names to ASVs using the program Claident
#  (Tanabe and Toju 2013)
# Reference: https://github.com/astanabe/MetabarcodingTextbook/
# claident 0.9.2025.04.13

# Paths of input fasta file and output directory
	INFILE=02_DADA2_Lib2018i/ASV.fa
	OUTDIR=04_Claident_Lib2018i
	
# Claident command options
## Reference database
	DB=animals_mt_genus # Animal (Metazoa) mitochondrial sequences which have species or lower level taxonomic information
## Number of CPU cores
  NCPU=$(sysctl -n hw.logicalcpu_max) # Use all available NCPU cores
## Two taxonomy assignment methods used
  METHOD1=QC # QCauto algorithm (Tanabe and Toju 2013)
  METHOD2=1,99% # Assign the name of top-1 sequence(s) with ≥99% identity
## Minimum required number of neighborhood sequences, necessary in METHOD2
  MINNNEIGHBORHOODSEQ=1
## Maximum permissible proportion of opposer sequences
  MAXPORROSER=0.10 # Allows ≤10% inconsistency (misidentification and/or synonyms) 
## Minimum permissible ratio of supporter:opposer sequences, the option to be
##  set together with 'MAXPORROSER'
  MINSORATIO=9 # Means that supporter:opposer = 9:1
## Minimum required number of supporter sequences, necessary in METHOD2
  MINNSUPPORTER=1

## Suffixes for taxonomy assignment files
	SUFFIX1=QC
	SUFFIX2=99NN
	SUFFIX3=merger

# Prepare an output directory
mkdir -p ${OUTDIR}

# Create a cached database for faster sequence search
clmakecachedb \
  --numthreads=${NCPU} \
  --blastdb=${DB} \
  ${INFILE} \
  "${OUTDIR}/blastdb_cache"

# METHOD1  
## Identify sequences that match the query
clidentseq \
  --numthreads=${NCPU} \
  --blastdb="${OUTDIR}/blastdb_cache" \
  --method=${METHOD1} \
  ${INFILE} \
  "${OUTDIR}/reference_${SUFFIX1}.dat"
## Assign taxonomy names using lowest common ancestor (LCA) algorithm
classigntax \
  --taxdb=${DB} \
  --maxpopposer=${MAXPORROSER} \
  --minsoratio=${MINSORATIO} \
  "${OUTDIR}/reference_${SUFFIX1}.dat" \
  "${OUTDIR}/taxonomy_${SUFFIX1}.tsv"

# METHOD2  
## Identify sequences that match the query
clidentseq \
  --numthreads=${NCPU} \
  --blastdb="${OUTDIR}/blastdb_cache" \
  --method=${METHOD2} \
  --minnneighborhoodseq=${MINNNEIGHBORHOODSEQ} \
  ${INFILE} \
  "${OUTDIR}/reference_${SUFFIX2}.dat"
## Assign taxonomy namesusing lowest common ancestor (LCA) algorithm
classigntax \
  --taxdb=${DB} \
  --maxpopposer=${MAXPORROSER} \
  --minsoratio=${MINSORATIO} \
  --minnsupporter=${MINNSUPPORTER} \
  "${OUTDIR}/reference_${SUFFIX2}.dat" \
  "${OUTDIR}/taxonomy_${SUFFIX2}.tsv"

# Merge the results of METHOD1 and METHOD2, prefering the former result
clmergeassign \
	--priority=descend \
	"${OUTDIR}/taxonomy_${SUFFIX1}.tsv" \
	"${OUTDIR}/taxonomy_${SUFFIX2}.tsv" \
	"${OUTDIR}/taxonomy_${SUFFIX3}.tsv"
