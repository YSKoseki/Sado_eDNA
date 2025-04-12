#!/bin/bash
#
# 01_Cutadapt_Lib2018i.sh
# A shell script to remove a pair of MiFish PCR primer sequences from MiSeq paired-end sequencing reads
# Reference: https://cutadapt.readthedocs.io/en/stable/guide.html
# cutadapt 5.0

# Parameters
INDIR='00_FASTQ_Lib2018i'
OUTDIR='01_Cutadapt_Lib2018i/FASTQ_trim'
LOGFILE='01_Cutadapt_Lib2018i/process.log'
STATSFILE='01_Cutadapt_Lib2018i/summary.tsv'

# Prepare directories and files
mkdir -p "${OUTDIR}"
: > "${LOGFILE}"
: > "${STATSFILE}"

# Write a header in the stats file, "${STATSFILE}"
echo 'sample\ttotalRP\tfiltdRP\tpctRP\ttotBP\tfiltdBP\tpctBP' > "${STATSFILE}"

# Find input files (paired-end Fastq)
FASTQR1=`find "${INDIR}" -type f -name '*_R1_001.fastq.gz' | sed 's!^.*/!!' | sort`

# Remove primers and filter reads by cutadapt
echo "${FASTQR1}" | while read line
do
  READ1=`echo "${line}"`
  READ2=`echo "${READ1}" | sed -e 's/R1/R2/g'`
  READ1TRIM=`echo "${READ1}" | sed -e 's/001.fastq/001_trim.fastq/g'`
  READ2TRIM=`echo "${READ2}" | sed -e 's/001.fastq/001_trim.fastq/g'`
  LABEL=`echo "${line}" | cut -f 1 -d '_'`
  # Run cutadapt
  OUTPUT=`cutadapt --cores=0 \
    -e 0.1 \
    -g ^NNNNNNGTCGGTAAAACTCGTGCCAGC \
    -G ^NNNNNNCATAGTGGGGTATCTAATCCCAGTTTG \
    --discard-untrimmed \
    -o "${OUTDIR}/${READ1TRIM}" -p "${OUTDIR}/${READ2TRIM}" \
    "${INDIR}/${READ1}" "${INDIR}/${READ2}"`
  # Save logs
  echo "${LABEL}\n\n${OUTPUT}\n\n\n" >> "${LOGFILE}" 2>&1
  # Save the selected stats
  if echo "${OUTPUT}" | grep 'No reads processed!' >/dev/null; then
    echo "${LABEL}\t0\t0\t\t0\t0\t" >> "${STATSFILE}" 2>&1
  else
    TOTRPS=`echo "${OUTPUT}" | grep 'Total read pairs' | sed -e 's/[^0-9]//g'`
    FILTRPS=`echo "${OUTPUT}" | grep 'Pairs written' | cut -f 2 -d '(' | sed -e 's/[^0-9]//g'`
    PCT_FILTRPS=`echo "${OUTPUT}" | grep 'Pairs written' | cut -f 3 -d '(' | tr -d ')'`
    TOTBPS=`echo "${OUTPUT}" | grep 'Total basepairs' | sed -e 's/[^0-9]//g'`
    FILTBPS=`echo "${OUTPUT}" | grep 'Total written' | cut -f 2 -d '(' | sed -e 's/[^0-9]//g'`
    PCT_FILTBPS=`echo "${OUTPUT}" | grep 'Total written' | cut -f 3 -d '(' | tr -d ')'`
    echo "${LABEL}\t${TOTRPS}\t${FILTRPS}\t${PCT_FILTRPS}\t${TOTBPS}\t${FILTBPS}\t${PCT_FILTBPS}" >> "${STATSFILE}" 2>&1
  fi
done

