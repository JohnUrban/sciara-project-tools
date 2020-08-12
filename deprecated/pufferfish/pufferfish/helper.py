

def run(parser, args):
    help_msg = '''


Was designed to be used on Sciara genome assemblies.
Best to use a quiver-polished, pilon-polished assembly.

Example:

Map reads from samples 1-N:
python pufferfish_main.py mapreads --bt2 sciara 1.fq.gz 2.fq.gz ... N.fq.gz




Get coverage:
python pufferfish_main.py getcov -g polished_assembly.genome -w 500 -s 500 -Q 30 bams/*bam



Make a "stage file" that has 2 columns (tab-sep): stagenumber (e.g. 1,2,3,4,5) and cov.bedGraph file
Each line represents 1 of the files generated above.



Run "find puffs" for various correlation/slope related approaches:
python ../pufferfish/pufferfish_main.py findpuffs -i stagefile -op name.pk.gz




After generating the pickled MultiCovBed object, get desired bedGraphs:
python pufferfish/pufferfish_main.py dump -ip data.fp.pk.gz -p prefix -c -fc -vc -cc -scc -sc



Use various awk commands to extract peaks from those bedGraphs:
From slope:
$ awk '$4 >1.5' q30.w500.s500.default.slope_collapsed.bedGraph | mergeBed -d 10000 -i - -c 4 -o median,max | awk '$3-$2 > 50e3'


From smoothed correlations:
$ awk '$4 >0.7' q30.w500.s500.default.smoothed_corscore_collapsed.bedGraph | mergeBed -d 20000 -i - -c 4 -o median,max | awk '$3-$2 > 50e3' | less


From state path through correlations:
$ awk '$4 == 3' q30.w500.s500.default.viterbi_collapsed.bedGraph | awk '$3-$2 >50e3'



Can also look at difference and fold enrichment from counts.bedGraph (-c) - example latest stage vs earliest:

Diff:
awk 'OFS="\\t" {print $1,$2,$3,$8-$4}' q30.w500.s500.default.counts.bedGraph | awk '$4 > 2' | mergeBed -i - -d 10000 -c 4 -o median,max | awk '$3-$2 > 50e3' 


FE:
awk 'OFS="\\t" {print $1,$2,$3,($8+0.1)/($4+0.1)}' q30.w500.s500.default.counts.bedGraph | awk '$4 > 2.5' | mergeBed -i - -d 10000 -c 4 -o median,max | awk '$3-$2 > 50e3' 


log2FE:
awk 'OFS="\\t" {print $1,$2,$3,log(($8+0.1)/($4+0.1))/log(2)}' q30.w500.s500.default.counts.bedGraph | awk '$4 > 1.5' | mergeBed -i - -d 10000 -c 4 -o median,max | awk '$3-$2 > 50e3'




For the 7-state HMM approach, run puffcn. There are 4 protocols that can be run with or without some time of normalization control (e.g. earliest stage).
After looking at the outcomes of all approaches, I'd recommend protocol 1 with early stage normalization. Protocol 4 with early normalization is good too.
The smaller the smoothing window, the more similar to P1 it becomes.
python ../pufferfish/pufferfish_main.py puffcn -l stage2.q30.w500.s500.bedGraph -e stage1.q30.w500.s500.bedGraph -1 -c > cn.stage2.p1.withearly.bedGraph
python ../pufferfish/pufferfish_main.py puffcn -l stage3.q30.w500.s500.bedGraph -e stage1.q30.w500.s500.bedGraph -1 -c > cn.stage3.p1.withearly.bedGraph
python ../pufferfish/pufferfish_main.py puffcn -l stage4.q30.w500.s500.bedGraph -e stage1.q30.w500.s500.bedGraph -1 -c > cn.stage4.p1.withearly.bedGraph
python ../pufferfish/pufferfish_main.py puffcn -l stage5.q30.w500.s500.bedGraph -e stage1.q30.w500.s500.bedGraph -1 -c > cn.stage5.p1.withearly.bedGraph

One is then able to extract peak coordinates from the state path bedGraph - for example:
awk '$4>1' cn.stage2.withearly.bedGraph | mergeBed -d 10000 -i - -c 4 -o median | awk '$3-$2 > 50e3' > cn.sag43.d10k.width50k.bed
awk '$4>1' cn.stage3.withearly.bedGraph | mergeBed -d 10000 -i - -c 4 -o median | awk '$3-$2 > 50e3' > cn.sag44.d10k.width50k.bed
awk '$4>1' cn.stage4.withearly.bedGraph | mergeBed -d 10000 -i - -c 4 -o median | awk '$3-$2 > 50e3' > cn.sag46.d10k.width50k.bed
awk '$4>1' cn.stage5.withearly.bedGraph | mergeBed -d 10000 -i - -c 4 -o median | awk '$3-$2 > 50e3' > cn.sag48.d10k.width50k.bed





    '''
    print help_msg
