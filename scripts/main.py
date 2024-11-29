import argparse
import os
import subprocess
import json
import matplotlib.pyplot as plt
import numpy as np
import statistics
import pandas as pd
import time

docker="micro2amplicon:latest"
def run(pe1,pe2,index,bed,outdir,prefix,ref,gff):
    start = time.time()
    # 输出最短序列30 or 35 bp
    pe1 = os.path.abspath(pe1)
    in_dir = os.path.dirname(pe1)
    outdir = os.path.abspath(outdir)
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    a = pe1.split("/")[-1]
    cmd = "docker run -v %s:/raw_data/ -v %s:/outdir/ %s " % (in_dir, outdir, docker)
    cmd += ("sh -c \'export PATH=/opt/conda/bin/:$PATH && fastp -i /raw_data/%s -o /outdir/%s.clean_R1.fastq "
            "--length_required 35 --dedup --thread 16 --low_complexity_filter --qualified_quality_phred 20 "
            "--html /outdir/%s.fastp.html --json /outdir/%s.fastp.json ") % (a, prefix, prefix, prefix)
    if pe2 is not None:
        pe2 = os.path.abspath(pe2)
        if in_dir != os.path.dirname(pe2):
            print("read1 and reads2 must be in the same directory.")
            exit()
        b = pe2.split("/")[-1]
        cmd += ("-I /raw_data/%s -O /outdir/%s.clean_R2.fastq\'") % (b, prefix)
    else:
        cmd += '\''
    print(cmd)
    subprocess.check_call(cmd, shell=True)

    # align reads with bowtie2 and sort bam with samtools
    index=os.path.abspath(index)
    host_index = ""
    for i in os.listdir(index):
        if i.endswith(".rev.2.bt2"):
            host_index = i.split(".rev.2.bt2")[0]
    cmd=f'docker run -v {index}:/ref/ -v {outdir}:/outdir/ {docker} sh -c \'export PATH=/opt/conda/bin/:\$PATH && bowtie2 '
    if pe2!=None:
        cmd+= f"--threads 48 -x /ref/{host_index} -1 /outdir/{prefix}.clean_R1.fastq -2 /outdir/{prefix}.clean_R2.fastq|samtools view -bh |samtools sort > /outdir/{prefix}.bam && samtools index /outdir/{prefix}.bam\'"
    else:
        cmd += f"--threads 48 -x /ref/{host_index} -U /outdir/{prefix}.clean_R1.fastq|samtools view -bh |samtools sort > /outdir/{prefix}.bam && samtools index /outdir/{prefix}.bam\'"
    subprocess.check_call(cmd, shell=True)

    # trim primers with ivar (soft clipping)
    if bed!=None:
        bed=os.path.abspath(bed)
        in_dir = os.path.dirname(bed)
        in_file=bed.split("/")[-1]
        cmd = f'docker run -v {in_dir}:/raw_data/ -v {outdir}:/outdir/ {docker} sh -c \'export PATH=/opt/conda/bin/:\$PATH && ivar '
        cmd += f"trim -e -i /outdir/{prefix}.bam -b /raw_data/{in_file} -p /outdir/{prefix}.soft.clipped | tee /outdir/{prefix}.ivar.stdout\'"
        print(cmd)
        subprocess.check_call(cmd, shell=True)

        ## remove soft-clipped primers
        cmd = f'docker run -v {outdir}:/outdir/ {docker} sh -c \'export PATH=/opt/conda/bin/:\$PATH && '
        cmd += f"samtools sort /outdir/{prefix}.soft.clipped.bam -o /outdir/{prefix}.soft.clipped.sort.bam\'"
        print(cmd)
        subprocess.check_call(cmd, shell=True)

        cmd = (f'docker run -v {outdir}:/outdir/ {docker} sh -c '
               f'\'java -jar /software/jvarkit.jar biostar84452 --samoutputformat BAM /outdir/{prefix}.soft.clipped.sort.bam |samtools sort -n >/outdir/{prefix}.trimmed.bam\'')
        print(cmd)
        subprocess.check_call(cmd, shell=True)

        # extract fastqs
        # https://www.htslib.org/doc/samtools-fasta.html
        if pe2 != None:
            cmd = (f'docker run -v {outdir}:/outdir/ {docker} sh -c '
                   f'\'export PATH=/opt/conda/bin/:\$PATH && samtools fastq -1 /outdir/{prefix}_trim_primer.R1.fq -2 /outdir/{prefix}_trim_primer.R2.fq -s /outdir/{prefix}_trim_primer.singleton.fastq /outdir/{prefix}.trimmed.bam &>/outdir/{prefix}.bam2fastq.stdout\'')
        else:
            cmd = (f'docker run -v {outdir}:/outdir/ {docker} sh -c '
                   f'\'export PATH=/opt/conda/bin/:\$PATH && samtools fastq -0 /outdir/{prefix}_trim_primer.R1.fq /outdir/{prefix}.trimmed.bam &>/outdir/{prefix}.bam2fastq.stdout\'')
        print(cmd)
        subprocess.check_call(cmd, shell=True)

    # Generate Pile-Up and variantCalling
    # https://github.com/CFSAN-Biostatistics/C-WAP/blob/main/startWorkflow.nf
    # -m    Minimum read depth to call variants (Default: 10)
    ref=os.path.abspath(ref)
    cmd = f'docker run -v {outdir}:/outdir/ -v {os.path.dirname(ref)}:/ref/ {docker} sh -c \'export PATH=/opt/conda/bin/:\$PATH && samtools mpileup -A -aa -d 0 -Q 0 -o /outdir/{prefix}.pile.up --reference /ref/{ref.split("/")[-1]}  /outdir/{prefix}.soft.clipped.sort.bam\''
    print(cmd)
    subprocess.check_call(cmd, shell=True)
    if gff!=None:
        gff=os.path.abspath(gff)
        cmd = (f'docker run -v {outdir}:/outdir/ -v {os.path.dirname(ref)}:/ref/ -v {os.path.dirname(gff)}:/raw_data/ {docker} sh -c \'export PATH=/opt/conda/bin/:\$PATH && cat /outdir/{prefix}.pile.up | '
               f'ivar variants -p /outdir/{prefix}.rawVarCalls -g /raw_data/{gff.split("/")[-1]} -r /ref/{ref.split("/")[-1]} -m 10\'')
    else:
        cmd = (
            f'docker run -v {outdir}:/outdir/ -v {os.path.dirname(ref)}:/ref/ -v {os.path.dirname(gff)}:/raw_data/ {docker} sh -c \'export PATH=/opt/conda/bin/:\$PATH && cat /outdir/{prefix}.pile.up | '
            f'ivar variants -p /outdir/{prefix}.rawVarCalls -r /ref/{ref.split("/")[-1]} -m 10\'')
    print(cmd)
    subprocess.check_call(cmd, shell=True)

    # Calculation of the consensus sequence using bcftools
    # https://github.com/niemasd/ViReflow/blob/main/ViReflow.py
    # https://github.com/CFSAN-Biostatistics/C-WAP/blob/main/startWorkflow.nf
    cmd = (f'docker run -v {outdir}:/outdir/ -v {os.path.dirname(ref)}:/ref/ {docker} sh -c \'export PATH=/opt/conda/bin/:\$PATH && '
           f'bcftools mpileup -Ou -f /ref/{ref.split("/")[-1]} /outdir/{prefix}.soft.clipped.sort.bam | '
           f'bcftools call --ploidy 1 -mv -Oz -o /outdir/{prefix}.calls.vcf.gz && '
           f'bcftools index /outdir/{prefix}.calls.vcf.gz && cat /ref/{ref.split("/")[-1]} '
           f'| bcftools consensus /outdir/{prefix}.calls.vcf.gz -p {prefix} > /outdir/{prefix}.consensus.fa\'')
    print(cmd)
    subprocess.check_call(cmd, shell=True)
    subprocess.check_call('sed -i \'1s/.*/>%s/\' %s/%s.consensus.fa' % (prefix, outdir,prefix), shell=True)

    # Calculation of the consensus sequence is used to determine the predominant lineage
    cmd = (f'docker run -v {outdir}:/outdir/ {docker} sh -c \'export PATH=/opt/conda/envs/pangolin/bin/:\$PATH && '
           f'pangolin --alignment /outdir/{prefix}.consensus.fa --threads 2 --outdir /outdir/ --outfile /outdir/{prefix}.lineage_report.csv --alignment-file /outdir/{prefix}.alignment.fasta\'')
    print(cmd)
    subprocess.check_call(cmd, shell=True)

    # Calculate depth from the trimmed mapped reads
    # http://www.htslib.org/doc/samtools-depth.html
    # "-J Include reads with deletions in depth computation."
    # "-q only count reads with base quality greater than or equal to INT"
    # https://github.com/niemasd/ViReflow/blob/main/ViReflow.py
    cmd = (f'docker run -v {outdir}:/outdir/ {docker} sh -c '
           f'\'export PATH=/opt/conda/bin/:\$PATH && samtools depth -J -d 8000 -Q 0 -q 20 -aa /outdir/{prefix}.soft.clipped.sort.bam >/outdir/{prefix}.depth.txt\'')

    subprocess.check_call(f'rm -rf {outdir}/{prefix}.depth.txt', shell=True)
    subprocess.check_call(cmd, shell=True)

    out=outdir+'/'+prefix
    ##############statistics result###########################
    outfile = open("%s.stat.tsv" % (out), "w")
    outfile.write(
        "SampleID\tRaw_reads\tQ20_rate(%)\tQ30_rate(%)\tClean_reads\tReads aligned(Trimmed primer)\tGenomic coordinates 0X(bp)\tGenomic coordinates <10X(bp)\n")
    with open("%s.json" % (out), "r") as load_f:
        load_dict = json.load(load_f)
    outfile.write("%s\t%s\t%s\t%s\t"
                  % (args.prefix, int(load_dict['summary']['before_filtering']['total_reads'] / 2),
                     format(float(load_dict['summary']['before_filtering']['q20_rate']) * 100, ".2f"),
                     format(float(load_dict['summary']['before_filtering']['q30_rate']) * 100, ".2f")))
    outfile.write("%s\t"
                  % (int(load_dict['summary']['after_filtering']['total_reads'] / 2)))

    infile = open("%s.bam2fastq.stdout" % out, "r")
    effect = []
    for line in infile:
        line = line.strip()
        array = line.split(" ")
        effect.append(int(array[2]))
    outfile.write("%s\t" % (int((effect[1] - effect[0]) / 2)))
    # A maximum of 1 000 000 reads are kept to limit the computation time of variant calling processes.
    if (int((effect[1] - effect[0]) / 2) > 1000000):
        subprocess.check_call("seqtk sample -s100 %s.R1.fq 1000000 > %s.sub.R1.fq && "
                              "seqtk sample -s100 %s.R2.fq 1000000 > %s.sub.R2.fq" % (out, out, out, out), shell=True)

    infile.close()
    infile = open("%s.depth.txt" % out, "r")
    cov = [0, 0]
    for line in infile:
        line = line.strip()
        array = line.split("\t")
        if int(array[2]) == 0:
            cov[0] += 1
        if int(array[2]) < 10:
            cov[1] += 1
    infile.close()
    outfile.write("%s\t%s" % (cov[0], cov[1]))
    outfile.close()
    end = time.time()
    ##############plot coverage###########################
    df = pd.read_csv("%s.depth.txt" % (out),
                     sep='\t',
                     # engine='python',
                     names=["ref", "pos", "depth"]
                     )

    median_depth = statistics.median(df["depth"])
    plt.figure(figsize=[10, 4])
    plt.axhline(median_depth, linestyle='--', color='red', linewidth=1, label="median: %.0f" % median_depth)
    plt.axhline(10, linestyle='--', color='grey', linewidth=1, label="<10X(%s bp)" % (cov[1]))

    max = np.max(10000)
    maxlog10 = np.ceil(np.log10(max))
    plt.ylim(top=10 ** maxlog10)

    plt.title("Sample: %s\n" % (args.prefix), fontsize=10, wrap=True)
    plt.xlabel("Position along genome [bp]")
    plt.ylabel("Coverage depth")
    plt.yscale("log")
    plt.margins(x=0.01)
    plt.legend()
    plt.ylim(bottom=1)
    plt.yscale("log")
    plt.plot(df["pos"], df["depth"])
    plt.savefig("%s.coverage.png" % (out), dpi=300)
    print("\nPre-process Done.\n")
    print("Elapse time is %g seconds" % (end - start))

if __name__ == '__main__':
    parser = argparse.ArgumentParser("A pipeline for .\nEmail:yucai.fan@illumina.com\n\n")
    parser.add_argument("-p1","--pe1",help="R1 fastq file",required=True)
    parser.add_argument("-p2","--pe2",help="R2 fastq file",default=None)
    parser.add_argument("-o","--outdir",help="output directory",default=os.getcwd())
    parser.add_argument("-p","--prefix",help="prefix of output",required=True)
    parser.add_argument("-i","--index",help="directory contains reference bowtie2 index",required=True)
    parser.add_argument("-b","--bed",help="bed file",default=None)
    parser.add_argument("-r", "--ref", help="reference fasta sequence",required=True)
    parser.add_argument("-g", "--gff", help="gff file",default=None)
    args = parser.parse_args()
    run(args.pe1,args.pe2,args.index,args.bed,args.outdir,args.prefix,args.ref,args.gff)