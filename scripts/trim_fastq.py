##  S Martínez-Puchol,  Itarte M ,  Rusiol M , et al. Exploring the diversity of coronavirus in sewage during COVID-19 pandemic: Don't miss the forest for the trees[J]. Science of The Total Environment, 2021, 800:149562.
##  Fontenele R S, Kraberger S, Hadfield J, et al. High-throughput sequencing of SARS-CoV-2 in wastewater provides insights into circulating variants[J]. Water Research, 2021, 205: 117710.
##  Baaijens J A, Zulli A, Ott I M, et al. Lineage abundance estimation for SARS-CoV-2 in wastewater using transcriptome quantification techniques[J]. Genome biology, 2022, 23(1): 236.
##  Email:yucai.fan@illumina.com
##  2023.03.31 version:1.0
##  2023.10.10 version:2.0 bug fix:
    # fastp reads longer than length_limit will be discarded, default 0 means no limitation. (35bp)
    # plot coverage (sample id)

import os
import subprocess
import argparse
import time
import json
import matplotlib.pyplot as plt
import numpy as np
import statistics
import pandas as pd

parser=argparse.ArgumentParser("A pipeline for lineage abundance estimation from wastewater sequencing data.\nEmail:yucai.fan@illumina.com\n\n")
parser.add_argument("-p1","--pe1",help="R1 fastq",required=True)
parser.add_argument("-p2","--pe2",help="R2 fastq",required=True)
parser.add_argument("-a","--adapter",help="fasta file adapter sequence",required=True)
parser.add_argument("-i","--index",help="prefix of covidRefSequences wuhan bowtie2 index",required=True)
parser.add_argument("-r","--reference",help="reference fasta sequence:NC_045512.2")
parser.add_argument("-g","--gff",help="gff file:NC_045512.2")
parser.add_argument("-t","--thread",help="number threads",type=int,default=24)
parser.add_argument("-b","--bed",help="primer bed file",required=True)
parser.add_argument('--min_base_qual', required=False, type=int, help="Minimum Base Quality for consensus sequence",default=20)
parser.add_argument("-p","--prefix",help="prefix of output",default=time.strftime("%Y-%m-%d"))
parser.add_argument("-o","--outdir",help="output directory",default=os.getcwd())
args=parser.parse_args()

start = time.time()
args.outdir=os.path.abspath(args.outdir)
args.adapter=os.path.abspath(args.adapter)
args.pe1=os.path.abspath(args.pe1)
args.pe2=os.path.abspath(args.pe2)
args.bed=os.path.abspath(args.bed)
args.reference=os.path.abspath(args.reference)
args.gff=os.path.abspath(args.gff)

out=args.outdir+"/"+args.prefix
if not os.path.exists(args.outdir):
    subprocess.check_call("mkdir -p %s"%(args.outdir),shell=True)

# http://www.usadellab.org/cms/?page=trimmomatic
# https://github.com/usadellab/Trimmomatic
# trim adapters with trimmomatic
# https://github.com/baymlab/wastewater_analysis/blob/main/pipeline/trim_reads.sh
# https://github.com/CFSAN-Biostatistics/C-WAP/blob/main/startWorkflow.nf
# https://github.com/niemasd/SD-COVID-Sequencing/blob/main/scripts/pipeline.sh
# https://github.com/niemasd/ViReflow/blob/main/ViReflow.py

#输出最短序列30 or 35 bp
cmd="fastp --in1 %s --in2 %s --out1 %s_1.fastq.gz --out2 %s_2.fastq.gz " \
    "--html %s.html --json %s.json --report_title %s " \
    "--thread 8 --adapter_fasta %s --length_required 35 --qualified_quality_phred 20"%(args.pe1,args.pe2,out,out,out,out,out,args.adapter)
print(cmd)
subprocess.check_call(cmd,shell=True)

# align reads with bowtie2 and sort bam with samtools
cmd="bowtie2 --threads %s -x %s -1 %s_1.fastq.gz -2 %s_2.fastq.gz |samtools view -bh |samtools sort > %s.bam && samtools index %s.bam && rm -rf %s_1.fastq.gz %s_2.fastq.gz"%(args.thread,args.index,out,out,out,out,out,out)
print(cmd)
subprocess.check_call(cmd,shell=True)

# trim primers with ivar (soft clipping)
# https://andersen-lab.github.io/ivar/html/manualpage.html
# -e    Include reads with no primers
cmd="ivar trim -e -i %s.bam -b %s -p %s.soft.clipped | tee %s.ivar.stdout && rm -rf %s.bam %s.bam.bai"%(out,args.bed,out,out,out,out)
print(cmd)
subprocess.check_call(cmd,shell=True)

## remove soft-clipped primers
#https://jvarkit.readthedocs.io/en/latest/Biostar84452/
# source activate && conda deactivate
cmd ="/software/samtools-v1.17/bin/samtools sort %s.soft.clipped.bam -o %s.soft.clipped.sort.bam"%(out,out)
print(cmd)
subprocess.check_call(cmd,shell=True)

cmd="java -jar /software/jvarkit.jar biostar84452 --samoutputformat BAM %s.soft.clipped.sort.bam |samtools sort -n >%s.trimmed.bam" %(out,out)
print(cmd)
subprocess.check_call(cmd,shell=True)

# extract fastqs
# https://www.htslib.org/doc/samtools-fasta.html
cmd="/software/samtools-v1.17/bin/samtools fastq -1 %s.R1.fq -2 %s.R2.fq -s %s.singleton.fastq %s.trimmed.bam &>%s.bam2fastq.stdout"%(out,out,out,out,out)
print(cmd)
subprocess.check_call(cmd,shell=True)

# Generate Pile-Up and variantCalling
# https://github.com/CFSAN-Biostatistics/C-WAP/blob/main/startWorkflow.nf
# -m    Minimum read depth to call variants (Default: 10)
cmd="/software/samtools-v1.17/bin/samtools mpileup -A -aa -d 0 -Q 0 -o %s.pile.up --reference %s %s.soft.clipped.sort.bam"%(out,args.reference,out)
print(cmd)
subprocess.check_call(cmd,shell=True)
cmd="cat %s.pile.up | ivar variants -p %s.rawVarCalls -g %s -r %s -m 10"%(out,out,args.gff,args.reference)
print(cmd)
subprocess.check_call(cmd,shell=True)

# Calculation of the consensus sequence using bcftools
# https://github.com/niemasd/ViReflow/blob/main/ViReflow.py
# https://github.com/CFSAN-Biostatistics/C-WAP/blob/main/startWorkflow.nf
cmd="bcftools mpileup -Ou -f %s %s.soft.clipped.sort.bam | bcftools call --ploidy 1 -mv -Oz -o %s.calls.vcf.gz"%(args.reference,out,out)
cmd+=" && bcftools index %s.calls.vcf.gz"%(out)
cmd+=" && cat %s | bcftools consensus %s.calls.vcf.gz -p %s > %s.consensus.fa"%(args.reference,out,args.prefix,out)
print(cmd)
subprocess.check_call(cmd,shell=True)
subprocess.check_call('sed -i \'1s/.*/>%s/\' %s.consensus.fa'%(args.prefix,out),shell=True)

# Calculation of the consensus sequence is used to determine the predominant lineage
cmd="pangolin --alignment %s.consensus.fa --threads 2 --outdir %s --outfile %s.lineage_report.csv --alignment-file %s.alignment.fasta"%(out,args.outdir,args.prefix,args.prefix)
print(cmd)
subprocess.check_call(cmd,shell=True)

# Calculate depth from the trimmed mapped reads
# http://www.htslib.org/doc/samtools-depth.html
# "-J Include reads with deletions in depth computation."
# "-q only count reads with base quality greater than or equal to INT"
# https://github.com/niemasd/ViReflow/blob/main/ViReflow.py
cmd="/software/samtools-v1.17/bin/samtools depth -J -d 8000 -Q 0 -q %s -aa %s.soft.clipped.sort.bam >%s.depth.txt"%(args.min_base_qual,out,out)
if not os.path.exists("%s.depth.tx"%(out)):
    subprocess.check_call(cmd,shell=True)

##############statistics result###########################
outfile = open("%s.stat.tsv" % (out), "w")
outfile.write("SampleID\tRaw_reads\tQ20_rate(%)\tQ30_rate(%)\tClean_reads\tReads aligned(Trimmed primer)\tGenomic coordinates 0X(bp)\tGenomic coordinates <10X(bp)\n")
with open("%s.json" % (out), "r") as load_f:
    load_dict = json.load(load_f)
outfile.write("%s\t%s\t%s\t%s\t"
              % (args.prefix,int(load_dict['summary']['before_filtering']['total_reads']/2),
                 format(float(load_dict['summary']['before_filtering']['q20_rate'])*100,".2f"),
                 format(float(load_dict['summary']['before_filtering']['q30_rate']) * 100, ".2f")))
outfile.write("%s\t"
              % (int(load_dict['summary']['after_filtering']['total_reads']/2)))

infile=open("%s.bam2fastq.stdout"%out,"r")
effect=[]
for line in infile:
    line=line.strip()
    array=line.split(" ")
    effect.append(int(array[2]))
outfile.write("%s\t"%(int((effect[1]-effect[0])/2)))
# A maximum of 1 000 000 reads are kept to limit the computation time of variant calling processes.
if(int((effect[1]-effect[0])/2)>1000000):
    subprocess.check_call("seqtk sample -s100 %s.R1.fq 1000000 > %s.sub.R1.fq && "
                          "seqtk sample -s100 %s.R2.fq 1000000 > %s.sub.R2.fq"%(out,out,out,out),shell=True)

infile.close()
infile=open("%s.depth.txt"%out,"r")
cov=[0,0]
for line in infile:
    line=line.strip()
    array=line.split("\t")
    if int(array[2])==0:
        cov[0]+=1
    if int(array[2])<10:
        cov[1]+=1
infile.close()
outfile.write("%s\t%s"%(cov[0],cov[1]))
outfile.close()
end=time.time()
##############plot coverage###########################
df = pd.read_csv("%s.depth.txt"%(out),
                 sep='\t',
                 # engine='python',
                 names=["ref", "pos", "depth"]
                 )

median_depth = statistics.median(df["depth"])
plt.figure(figsize=[10, 4])
plt.axhline(median_depth, linestyle='--', color='red', linewidth=1, label="median: %.0f" % median_depth)
plt.axhline(10, linestyle='--', color='grey', linewidth=1,label="<10X(%s bp)"%(cov[1]))

max = np.max(10000)
maxlog10 = np.ceil(np.log10(max))
plt.ylim(top=10 ** maxlog10)

plt.title("Sample: %s\nAccession: NC_045512.2"%(args.prefix), fontsize=10, wrap=True)
plt.xlabel("Position along genome [bp]")
plt.ylabel("Coverage depth")
plt.yscale("log")
plt.margins(x=0.01)
plt.legend()
plt.ylim(bottom=1)
plt.yscale("log")
plt.plot(df["pos"],df["depth"])
plt.savefig("%s.coverage.png"%(out), dpi=300)
print("\nPre-process Done.\n")
print("Elapse time is %g seconds" % (end - start))