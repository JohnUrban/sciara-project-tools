#!/usr/bin/env python
import gzip, sys, argparse
import numpy as np
from Bio import SeqIO
##from fqEncoding import detectEncoding, rosettaStone, baseQC, baseQCmeans, baseQCquantile
##import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description="""

DESCRIPTION
    
    Takes in fastq assembly file. Reports shit.

    q = -10 * log_{10} p_{error}
    p{error} = 10^{ q/-10 }
    p{accurate} = 1 - p{error}
    accuracy = 100*p{accuracte}

    Quiver notes:
    - https://github.com/PacificBiosciences/GenomicConsensus/blob/master/doc/HowToQuiver.rst
    - https://github.com/PacificBiosciences/GenomicConsensus/blob/master/doc/QuiverFAQ.rst

    - If there are errors remaining applying Quiver, they will almost invariably be homopolymer run-length errors (insertions or deletions).

    - At present, the Quiver model is trained per-chemistry, so it is very important that Quiver knows the sequencing chemistries used.

    - If SMRT Analysis software was used to build the cmp.h5 file, the cmp.h5 will be loaded with information about the sequencing chemistry used for each SMRT Cell, and Quiver will automatically identify the right parameters to use.
    - If custom software was used to build the cmp.h5, or an override of Quiver's autodetection is desired, then the chemistry or model must be explicity entered.
    - For example:
        $ quiver -p P4-C2
        $ quiver -p P4-C2.AllQVsMergingByChannelModel

    - The most likely cause for true errors made by Quiver is that the coverage in the region was low.
    - If there is 5x coverage over a 1000-base region, then 10 errors in that region can be expected.
    - It is important to understand that the effective coverage available to Quiver is not the full coverage apparent in plots
        - Quiver filters out ambiguously mapped reads by default.
        - The remaining coverage after filtering is called the "effective coverage".
        - If you have verified that there is high effective coverage in the region in question, it is highly possible---given the high accuracy Quiver can achieve---that the apparent errors actually reflect true sequence variants.
        - Inspect the FASTQ output file to ensure that the region was called at high confidence; if an erroneous sequence variant is being called at high confidence, please report a bug to us.

    - For regions with no effective coverage, no variants are outputted, and the FASTQ confidence is 0.
    - The output in the FASTA and FASTQ consensus sequence tracks is dependent on the setting of the --noEvidenceConsensusCall flag.
    - Assuming the reference in the window is "ACGT", the options and outputs are:
        --noEvidenceConsensusCall={option}
            options = "nocall" {outputs NNNN}, "reference" {outputs ACGT}, "lowercasereference {outputs acgt}

    - MapQV is a single scalar Phred-scaled QV per aligned read that reflects the mapper's degree of certainty that the read aligned to this part of the reference and not some other.
        - Unambigously mapped reads will have a high MapQV (typically 255), while a read that was equally likely to have come from two parts of the reference would have a MapQV of 3.
    - MapQV is pretty important when you want highly accurate variant calls.
        - Quiver and Plurality both filter out aligned reads with a MapQV below 20 (by default), so as not to call a variant using data of uncertain genomic origin.
    - This can be problematic if using Quiver to get a consensus sequence.
        - If the genome of interest contains long (relative to the library insert size) highly-similar repeats, the effective coverage (after MapQV filtering) may be reduced in the repeat regions---this is termed these MapQV dropouts.
        - If the coverage is sufficiently reduced in these regions, Quiver will not call consensus in these regions.
    - If you want to use ambiguously mapped reads in computing a consensus for a denovo assembly, the MapQV filter can be turned off entirely.
        - In this case, the consensus for each instance of a genomic repeat will be calculated using reads that may actually be from other instances of the repeat, so the exact trustworthiness of the consensus in that region may be suspect.
    - The MapQV filter can be disabled using the flag --mapQvThreshold=0 (shorthand: -m=0).
        - Consider this in de novo assembly projects, but it is not recommended for variant calling applications.

    - Quiver limits read coverage, filters reads by MapQV, and filters variants by quality and coverage.
        - The overall read coverage used to call consensus in every window is 100x by default, but can be changed using -X=value.
        - The MapQV filter, by default, removes reads with MapQV < 20.
            - This is configured using --mapQvThreshold=value / -m=value
        - Variants are only called if the read coverage of the site exceeds 5x, by default---this is configurable using -x=value.
            - Further, they will not be called if the confidence (Phred-scaled) does not exceed 40---configurable using -q=value.

    - iterating Quiver
    - Some customers have found it useful for repetitive genomes. It is recommended for repetitive genomes, but there are warnings. See links above.

    --> Also email PacBio with the question!
        support@pacb.com
        1.877.920.PACB (7222)
    """, formatter_class= argparse.RawTextHelpFormatter)

parser_input = parser.add_mutually_exclusive_group(required=True)
parser_input.add_argument("-f", "--fastx",
                   type= str, default=False,
                   help='''Input fastq file.''')
parser.add_argument("--case",
                    action="store_true", default=False,
                    help=''' For each contig, also count number of lower-case and upper-case letters. Reports percents as well. This is useful, since quiver makes bases without support lower-case.''')
parser.add_argument("--caseqv",
                    action = "store_true", default=False,
                    help=''' For each contig, partition into upper and lower-case components and compute QV stats on bases of each.''')

parser.add_argument("-H", "--header",
                    action = "store_true", default=False,
                    help=''' Print header.''')

parser.add_argument("-g", "--genomestats",
                    action = "store_true", default=False,
                    help=''' Also include genome stats at end...''')

parser.add_argument("-b", "--getlowercasebed", action='store_true', default=False,
                    help='''This option over-rides default output behavior and outputs BED coordinates of stretches of lower-case sequence.
                            In the context of a quiver fastq of an assembly, these are regions that had low effective coverage,
                            and may be candidates for mis-assemblies.''')
parser.add_argument("-i", "--internalonly", action='store_true', default=False,
                    help='''This option is used with "-b"/"--getlowercasebed". Use it to only return internal stretches of lower-case letters,
                            i.e. not the stretches that usually occur at the beginning and end of a contig.''')
args = parser.parse_args()


contig_qvs = []
contig_lengths = []
allqv = []

def case_counter(seq):
    return sum(1 for b in seq if b.isupper())

def lowercase_regions(seqname, seq, internalonly=False):
##    regions = []
    lenseq = len(seq)
    s = None
    e = None
    i = -1
    for b in seq:
        i += 1
        if b.islower():
            if s is None:
                s = i
        else: #upper
            if s is not None:
                e = i
##                regions.append((s,e))
                if internalonly and s == 0:
                    s = None
                    e = None
                    continue
                print ("\t").join([seqname, str(s), str(e)])
                s = None
                e = None
    ## process end of contig
    if s is not None and not internalonly:
        print ("\t").join([seqname, str(s), str(lenseq)])
        

def get_length(contig):
    return len(contig)

def get_qv_info(contig):
    qv_mean = np.mean(contig.letter_annotations['phred_quality'])
    qv_median = np.median(contig.letter_annotations['phred_quality'])
    qv_min = np.min(contig.letter_annotations['phred_quality'])
    qv_max = np.max(contig.letter_annotations['phred_quality'])
    return [qv_mean, qv_median, qv_min, qv_max]
    

if args.header and not args.getlowercasebed:
    if args.case:
        sys.stdout.write(("\t").join(["contig", "length", "mean_QV", "median_QV", "min_QV", "max_QV", "n_uppercase", "pct_uppercase"]) +"\n")
    else:
        sys.stdout.write(("\t").join(["contig", "length", "mean_QV", "median_QV", "min_QV", "max_QV"]) +"\n")

if args.getlowercasebed:
    for contig in SeqIO.parse(args.fastx, "fastq"):
        lowercase_regions(contig.description, str(contig.seq), args.internalonly)
        
else:
    for contig in SeqIO.parse(args.fastx, "fastq"):
        length = get_length(contig)
        contig_lengths.append(length)
        qv_mean = np.mean(contig.letter_annotations['phred_quality'])
        qv_median = str(np.median(contig.letter_annotations['phred_quality']))
        qv_min = str(np.min(contig.letter_annotations['phred_quality']))
        qv_max = str(np.max(contig.letter_annotations['phred_quality']))
        contig_qvs.append(qv_mean)
        if args.case:
            nlow = case_counter(str(contig.seq))
            sys.stdout.write(("\t").join([contig.name, str(length), str(qv_mean), qv_median, qv_min, qv_max, str(nlow), str(100.0*nlow/length)]) +"\n")
        else:
            sys.stdout.write(("\t").join([contig.name, str(length), str(qv_mean), qv_median, qv_min, qv_max]) +"\n")

    if args.genomestats:
        G=sum(contig_lengths)
        sys.stdout.write(("\t").join(["contig_mean_qv\t" + str(G), str(np.mean(contig_qvs)), str(np.median(contig_qvs)),  str(np.min(contig_qvs)), str(np.max(contig_qvs))]) + "\n")

        ##genome_qv = np.mean(allqv)
        genome_qv = 0
        for i in range(len(contig_qvs)):
            genome_qv += contig_lengths[i]*contig_qvs[i]
        genome_qv = float(genome_qv)/G
        sys.stdout.write("genome_qv\t" + str(genome_qv) + "\n")
