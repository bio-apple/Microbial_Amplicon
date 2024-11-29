# This pipeline could analysis the microbial amplicon.

## step1: Build your docker

    cd Docker/
    docker build -t micro2amplicon ./

## step2: Download the nextclade database

    cd ref/
    python3 download_nextclade_db.py -d nextclade_db/

## step3: prepare reference fasta(required)，bed file(required)



    

 