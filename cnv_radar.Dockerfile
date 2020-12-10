#################################################################
# Dockerfile
#
# This can be used to run the CNV Radar tool
#
# Software: R
# Software Version: 3.6.3 Holding the windsock
#
# Software: SNPEFF
# Software Version: 4_3t
#
# Software: Bedtools
# Software Version: 2.29.2
#
# Software: CNV Radar
# Software Version: 1.2.1
#
# Description: Create the ROI summary from a bam file
# 	       Create the reference panel
#	       Call CNV on target sample
#	       Annotate VCF files with snpEff
#
# Provides:
# Base Image: Ubuntu 16.04
# Build Cmd: docker build --rm -t eagenomics/cnvradar:v1.2.1 -f cnv_radar.Dockerfile .
# Pull Cmd: docker pull eagenomics/cnvradar:v1.2.1
# Run Cmd: docker run --rm -v ${PWD}:/data -w /data eagenomics/cnvradar:v1.2.1 /bin/bash -c ""
#################################################################

FROM ubuntu:16.04

#----------------------------------------------------------------
# Container Metadata
#----------------------------------------------------------------
LABEL base.image="ubuntu:16.04"
LABEL version="2"
LABEL software="CNVRadar, bedtools, snpEff, snpSift"
LABEL software.version="1.2.1"
LABEL about.summary="a copy number variant detection algorithm from Expression Analysis"
LABEL license="https://github.com/ExpressionAnalysis/CNV_Radar/blob/master/LICENSE.txt"
LABEL about.tags="Genomics"

#----------------------------------------------------------------
# Install command line tools and packages
#----------------------------------------------------------------
ENV DEBIAN_FRONTEND noninteractive

RUN apt-get -qq update && apt-get -y upgrade && \
	apt-get install -y --no-install-recommends \
    build-essential \
    pkg-config \
    python3-pip \
    python3-dev \
    ca-certificates \
    wget \
    bzip2 \
    unzip \
    tabix \
    libz-dev \
    libbz2-dev \
    liblzma-dev \
    openjdk-8-jdk \
    openjdk-8-jre \
    cmake && \
    apt-get clean && \
    apt-get autoremove && \
	ln -s /usr/bin/python3 /usr/bin/python

#----------------------------------------------------------------
# Setup Java
#----------------------------------------------------------------
ENV JAVA_HOME /usr/lib/jvm/java-8-openjdk-amd64/jre/bin/java
ENV PATH ${PATH}:${JAVA_HOME}
RUN echo "JAVA_HOME is set to ${JAVA_HOME}"

#----------------------------------------------------------------
# Install Bedtools
#----------------------------------------------------------------
ENV BEDTOOLS_VERSION 2.29.2
RUN wget --no-check-certificate https://github.com/arq5x/bedtools2/releases/download/v${BEDTOOLS_VERSION}/bedtools-${BEDTOOLS_VERSION}.tar.gz -O /tmp/bedtools.tar.gz && \
	tar zxvf /tmp/bedtools.tar.gz -C /usr/local/bin/ && \
    rm /tmp/bedtools.tar.gz && \
	cd /usr/local/bin/bedtools2 && \
	make

ENV PATH ${PATH}:/usr/local/bin/bedtools2/bin/
RUN bedtools --version

#----------------------------------------------------------------
# Install snpEff and annotation databases
#----------------------------------------------------------------
ENV SNPEFF_VERSION 4_3t

RUN wget --no-check-certificate http://sourceforge.net/projects/snpeff/files/snpEff_v${SNPEFF_VERSION}_core.zip -O /usr/local/bin/snpeff.zip && \
    unzip /usr/local/bin/snpeff.zip -d /usr/local/bin/ && \
    rm /usr/local/bin/snpeff.zip && \
    #bash -c 'echo -e "java -jar /usr/local/bin/snpEff/snpEff.jar  \$@" > /usr/local/bin/snpeff' && \
    #chmod +x /usr/local/bin/snpeff && \
    java -jar /usr/local/bin/snpEff/snpEff.jar -version 
    #bash -c 'echo -e "#!/bin/bash\njava -jar /usr/local/bin/snpEff/SnpSift.jar \$@" > /usr/local/bin/snpsift' && \
    #chmod +x /usr/local/bin/snpsift && \
    
ENV PATH ${PATH}:/usr/local/bin/snpEff/scripts:/usr/local/bin/snpEff/scripts/gsa

#----------------------------------------------------------------
# Install R 
#----------------------------------------------------------------

ENV R_VERSION 3.6.3
ENV OS_IDENTIFIER "ubuntu-1604"

RUN wget --no-check-certificate https://cdn.rstudio.com/r/${OS_IDENTIFIER}/pkgs/r-${R_VERSION}_1_amd64.deb && \
    apt-get update -qq && \
    DEBIAN_FRONTEND=noninteractive apt-get install -f -y ./r-${R_VERSION}_1_amd64.deb && \
    ln -s /opt/R/${R_VERSION}/bin/R /usr/bin/R && \
    ln -s /opt/R/${R_VERSION}/bin/Rscript /usr/bin/Rscript && \
    ln -s /opt/R/${R_VERSION}/lib/R /usr/lib/R && \
    rm r-${R_VERSION}_1_amd64.deb && \
    rm -rf /var/lib/apt/lists/*

#----------------------------------------------------------------
# Install R Packages
#----------------------------------------------------------------
RUN echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" > ~/.Rprofile
RUN Rscript -e "install.packages('getopt')"
RUN Rscript -e "install.packages('data.table')"
RUN Rscript -e "install.packages('R.utils', dependencies = T)"
RUN Rscript -e "install.packages('yaml', dependencies = T)"
RUN Rscript -e "library('getopt');##### R SESSION INFORMATION #####; sessionInfo()"

#----------------------------------------------------------------
# Install Third Party Annotation Databases
#----------------------------------------------------------------
ENV GRCH38_VERSION 86
ENV GRCH37_VERSION 75

#RUN java -jar /usr/local/bin/snpEff/snpEff.jar download GRCh38.${GRCH38_VERSION} && \
#    java -jar /usr/local/bin/snpEff/snpEff.jar download GRCh37.${GRCH37_VERSION} 

#RUN mkdir -p /usr/local/bin/snpEff/db/GRCh37/dbSnp && \
#    cd /usr/local/bin/snpEff/db/GRCh37/dbSnp && \
#    wget --no-check-certificate ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/All_20180423.vcf.gz && \
#    wget --no-check-certificate ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/All_20180423.vcf.gz.tbi && \
#    wget --no-check-certificate ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/All_20180423.vcf.gz.md5 && \
#    md5sum --check All_20180423.vcf.gz.md5 && \
#    rename 's/All_20180423/dbSnp/g' * && \
#    mkdir -p /usr/local/bin/snpEff/db/GRCh38/dbSnp && \
#    cd /usr/local/bin/snpEff/db/GRCh38/dbSnp && \
#    wget --no-check-certificate ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/All_20180418.vcf.gz && \
#    wget --no-check-certificate ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/All_20180418.vcf.gz.tbi && \
#    wget --no-check-certificate ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/All_20180418.vcf.gz.md5 && \
#    md5sum --check All_20180418.vcf.gz.md5 && \
#    rename 's/All_20180418/dbSnp/g' *


#----------------------------------------------------------------
# Copy over analysis scripts
#----------------------------------------------------------------
RUN wget --no-check-certificate https://github.com/ExpressionAnalysis/CNV_Radar/archive/master.zip -O /opt/CNVRadar.zip && \
    unzip /opt/CNVRadar.zip -d /opt/CNVRadar && \
    mv /opt/CNVRadar/CNV_Radar-master/* /opt/CNVRadar/ && \
    rm -r /opt/CNVRadar/CNV_Radar-master/ && \
    rm /opt/CNVRadar.zip

#RUN Rscript /opt/CNVRadar/CNV_Radar.r -h

#----------------------------------------------------------------
# Set working dir
#----------------------------------------------------------------
WORKDIR /data/
