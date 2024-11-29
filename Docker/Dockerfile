FROM alpine AS nextclade
# glibc+conda+nextclade v3.8.2
RUN apk update && \
    apk add --no-cache bash openjdk21 git && mkdir -p /lib64/ /ref/ /script/ /raw_data/ /outdir/ && \
    wget -q -O /etc/apk/keys/sgerrand.rsa.pub https://alpine-pkgs.sgerrand.com/sgerrand.rsa.pub && \
    wget https://github.com/sgerrand/alpine-pkg-glibc/releases/download/2.35-r1/glibc-2.35-r1.apk && \
    apk add --no-cache --force-overwrite glibc-2.35-r1.apk && \
    rm glibc-2.35-r1.apk && ln -s /usr/glibc-compat/lib/* /lib64/ && \
    wget -q -O /opt/Miniforge3-Linux-x86_64.sh https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh && \
    bash /opt/Miniforge3-Linux-x86_64.sh -f -b -p /opt/conda/ && rm -rf /opt/Miniforge3-Linux-x86_64.sh /var/cache/apk/* && \
    cd /bin/ && wget -O nextclade https://github.com/nextstrain/nextclade/releases/latest/download/nextclade-x86_64-unknown-linux-musl && chmod u+x ./nextclade

FROM nextclade AS pangolin
RUN git clone https://github.com/cov-lineages/pangolin.git && cd pangolin/ &&  \
    /opt/conda/bin/mamba env create -f environment.yml --name pangolin &&  \
    /opt/conda/envs/pangolin/bin/pip install . &&  rm -rf /opt/pangolin/ && \
    /opt/conda/bin/conda clean -a -y

FROM pangolin AS micro2amplicon
COPY jvarkit.jar /software/
RUN /opt/conda/bin/mamba install --channel conda-forge --channel bioconda --channel defaults ivar=1.3 trimmomatic  \
    bowtie2 bbmap fastp seqtk samtools bedtools bcftools bwa prinseq cutadapt drep
RUN /opt/conda/bin/pip3 install seaborn matplotlib numpy pysam pandas
