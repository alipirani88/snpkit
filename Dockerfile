FROM ubuntu:18.04

RUN apt-get update && apt-get upgrade -y

RUN apt-get install wget -y && apt-get install unzip -y

RUN wget https://repo.anaconda.com/archive/Anaconda3-2020.02-Linux-x86_64.sh

Run bash Anaconda3-2020.02-Linux-x86_64.sh -b -p

ENV PATH "/root/anaconda3/bin/:$PATH"

RUN wget https://github.com/alipirani88/variant_calling_pipeline/archive/master.zip

RUN unzip master.zip

#RUN conda create -f code/variant_calling_pipeline-master/docker_env -n varcall

#RUN conda activate varcall

RUN conda install -c bioconda trimmomatic

RUN conda install -c bioconda bcftools
RUN conda install -c bioconda bedtools
RUN conda install -c bioconda bioawk
RUN conda install -c bioconda bowtie2
RUN conda install -c bioconda bwa
RUN conda install -c bioconda bzip2
RUN conda install -c bioconda gatk
RUN conda install -c bioconda gatk4
RUN conda install -c bioconda vcftools
RUN conda install -c bioconda mummer
RUN conda install -c bioconda mash
RUN conda install -c bioconda picard
RUN conda install -c bioconda pilon
RUN conda install -c bioconda qualimap
RUN conda install -c bioconda raxml
RUN conda install -c bioconda samtools
RUN conda install -c bioconda snpeff
#RUN conda install -c conda-forge subprocess32
RUN pip install subprocess32
RUN conda install -c bioconda tabix
RUN conda install -c bioconda pyfasta
