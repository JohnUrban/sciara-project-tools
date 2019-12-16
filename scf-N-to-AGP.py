#!/usr/bin/env python2.7
import argparse, sys
import re
from Bio import SeqIO

# Parse command-line arguments
parser = argparse.ArgumentParser(description='''
Description: Takes in Fasta.
Outputs AGP given N-gap locations.
Also outputs contigs and scaffolds fasta file.

For now, it is hard-coded to assume you're breaking apart a BioNano Map scaffolded assembly.

This script is similar to scf-N-to-BED.py (copied from there for starters).

Use '-' or 'stdin' if coming from stdin.
Use -l/--length to control minimum gap length (default 25).
Use -r/--recognition to allow regex to expand gaps through recognition sequence islands dispersed in N-gaps (e.g. from optical maps).
When using -r provide only the recognition sequence. For now it is case-sensitive.
By default N-gaps are simply where Ns start to the first non-N.

Example recognition sequences:
BssSI: CACGAG
BspQI: GCTCTTC

AGP coordinates are 1-based inclusive.

AGP entry for contig:
name start end part_number component_type name_id contig_start contig_end orientation
component_type = W
orientation = +


AGP entry for gap:
name start end part_number component_type gapLen gapType linkage evidence
component_type = N
gapType = scaffold
linkage = yes
evidence = map


Example of scaffold with 2 contigs separated by one gap from map-evidence:
trial_scaffold_1        1       9       1       W       trial_scaffold_1_1      1       9       +
trial_scaffold_1        10      33      2       N       24      scaffold        yes     map
trial_scaffold_1        34      39      3       W       trial_scaffold_1_3      1       6       +


I started using a similar tool (fasta2agp.pl) found here:
    wget --no-check-certificate --content-disposition https://gist.githubusercontent.com/IdoBar/8a324fff97c6a866ea7e0cbcec41d570/raw/7faa9e342d62355b11c79f114ec6367e1de9525e/fasta2agp.pl


- I found that fasta2agp.pl can only filter by scaffold size, but
        - can't discriminate gap sizes (breaks at 1 bp)
        - can't mow over recseqs (so treats as multiple gaps)
            - BioNano scaffolding has huge gaps that can have multiple short (6 bp) recognition sequences within them
            - This can be useful in some contexts
            - It may not be useful to treat them all as contigs though...
            
- I created my own tool to do the same thing as fasta2agp.pl with the additional features I needed.
- It also default outputs what I needed.
        - It is called:
            - scf-N-to-AGP.py (this script)
        - It follows fasta2agp.pl behavior on contigs that do not start with N
                - If contig starts with N
                        - fasta2agp.pl ignores it as if it were not there and all coordinates are shifted
                        - scf-N-to-AGP.py gives a warning, but includes it as a "part" of the object
        - scf-N-to-AGP.py is able to mow over the recseqs and only break contigs at gap sizes exceeding some minlength.

        

''', formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("fasta")
parser.add_argument("-l", "--length", type=int, default=25, help='''Do not report as gap if less than this length. Default=25.''')
parser.add_argument("-r", "--recognition", type=str, default=False, help='''Allow N-gaps to encompass occurences of given recognition sequence. Default=False.''')
parser.add_argument("-G", "--gapchar", type=str, default='N', help='''The single letter gap character. Default =  N.''')
parser.add_argument("-P", "--prefix", type=str, default='', help='''Add this prefix in front of all sequence names.''')
parser.add_argument("-o", "--outprefix", type=str, default='', help='''Add this prefix in front of all output filenames. Can be path to directory or path + prefix as well.''')
#parser.add_argument("-D", "--do_not_allow_start_end_with_N", type=str, default='', help='''Add this prefix in front of all sequence names.''')

args = parser.parse_args()



def get_contig_info(record, gap, scfstart, part, prefix=''):
    scfend = gap.start() # No need to subtract 1 as this is 0-based (need to add 1 for real gap start)
    ctgend = scfend - scfstart + 1
    return (str(e) for e in (prefix+str(record.id), scfstart, scfend, part, 'W', prefix+str(record.id)+'_'+str(part), 1, ctgend, '+'))


def get_last_contig_info(record, gap, scflength, scfstart, part, prefix=''):
    scfend = scflength
    ctgend = scfend - scfstart + 1
    return (str(e) for e in (prefix+str(record.id), scfstart, scfend, part, 'W', prefix+str(record.id)+'_'+str(part), 1, ctgend, '+'))


def get_gap_info(record, gap, part, gapType='scaffold', linkage='yes', evidence='map', prefix=''):
    gapLen = gap.end()-gap.start()
    return (str(e) for e in (prefix+str(record.id), gap.start()+1, gap.end(), part, 'N', gapLen, gapType, linkage, evidence))

def get_file_name(x,pre, sfx):
    if not pre:
        pre = x
    else:
        pre = pre + '.' + x
    name = pre + '.' + sfx
    sys.stderr.write('Will write to ' + name + ' ....\n')
    return name

def append_sequence_to_fasta(fh, record, prefix='', sep='', part='', start=-1, end=-1):
    if start == -1 or end == -1:
        # Return whole record (probably to scaffolds fasta)
        fh.write('>'+str(prefix)+str(record.id)+str(sep)+str(part)+'\n')
        fh.write(str(record.seq)+'\n')
    else:
        # Return slice of record (probably to scaffolds fasta)
        fh.write('>'+str(prefix)+str(record.id)+str(sep)+str(part)+'\n')
        fh.write(str(record.seq)[start:end]+'\n')

    
# OPEN CONNECTION TO FILE/STDIN
handle = sys.stdin if args.fasta in ('-','stdin') else open(args.fasta)

# DEFINE REGEX, 
N = args.gapchar + '+' ## 'N+'
pattern = N
if args.recognition:
    pattern = N + '((' + args.recognition + '){0,1}' + N + ')*' ## 1+ Ns and optionally 0 or more (recseq, Ns) 

# Shorten some variable names
P = args.prefix


# Open contig, fasta, and AGP files for writing
scaffolds = open(get_file_name('scaffolds', args.outprefix, 'fasta'), 'w')
contigs = open(get_file_name('contigs', args.outprefix, 'fasta'), 'w')
AGP = open(get_file_name('AGP', args.outprefix, 'agp'), 'w')

# LOOP OVER FASTA
for record in SeqIO.parse(handle, "fasta"):
    # Write out scaffold
    append_sequence_to_fasta(scaffolds, record, prefix=P)
    
    # Length of scaffold
    scflength = len(record.seq)

    # Gap locations
    allgaps = [e for e in re.finditer(pattern, str(record.seq))]

    # Number of gaps found
    numgaps = len(allgaps)

    # Process scaffold accordingly
    if numgaps == 0:
        # Single contig - write AGP 
        outline = ("\t").join([str(e) for e in (P+str(record.id), 1, scflength, 1, 'W', P+str(record.id)+'_1', 1, scflength, '+')])
        AGP.write( outline + '\n')
        # Write contig 
        append_sequence_to_fasta(contigs, record, prefix=P, sep='_', part=1)

    elif numgaps > 0: # 1 or more gaps
        
        # Describe initial contig
        scfstart = 1
        part = 1
        ended_on_gap = False
        for gap in allgaps:
            # If current gap is not minimum gap length, do not report it, do not update current scfstart variable, just continue to next gap
            if gap.end()-gap.start()  < args.length:
                continue
            
            # Warn about leading Ns if present, report gap, and move on
            if gap.start() == 0:
                sys.stderr.write("Warning: "+str(record.id)+" started with N. The leading N-gap will be reported as the first part of this scaffold. NCBI does not accept this though.\n")

                # Report gap info
                AGP.write( ("\t").join(get_gap_info(record, gap, part, prefix=P)) + '\n' )

                # Update scfstart s.t. it starts after this gap
                scfstart = gap.end()+1

                #Update part number for next item
                part += 1
                continue

            

            # Process previous contig (and update part number):
            #   Since iterating over gaps, contigs precede current gap starting at most recent scfstart variable and ending at gap
            # Write AGP
            AGP.write( ("\t").join( get_contig_info(record, gap, scfstart, part, prefix=P) ) + '\n' )
            # Write contig
            append_sequence_to_fasta(contigs, record, prefix=P, sep='_', part=part, start=scfstart-1, end=gap.start())
            # Update part
            part += 1


            # Process current gap, update scfstart for next contig and part number
            AGP.write( ("\t").join(get_gap_info(record, gap, part, prefix=P)) + '\n' )
            scfstart = gap.end()+1
            part += 1
            
            # Warn about trailing Ns if present, report gap (above), and break
            if gap.end() == scflength:
                ended_on_gap = True
                sys.stderr.write("Warning: "+str(record.id)+" ended with N. The trailing N-gap will be reported as the last part of this scaffold. NCBI does not accept this though.\n")
                break



        # Process last contig (if it did not end on a gap)
        if not ended_on_gap:
            # Write AGP
            AGP.write( ("\t").join( get_last_contig_info(record, gap, scflength, scfstart, part, prefix=P) ) + '\n' )
            # Write contig
            append_sequence_to_fasta(contigs, record, prefix=P, sep='_', part=part, start=scfstart-1, end=scflength)
                


# Close input file            
if args.fasta not in ('-','stdin'):
    handle.close()


# Close output files
scaffolds.close()
contigs.close()
AGP.close()



# END


