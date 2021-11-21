FROM bioconductor/bioconductor_docker:devel
MAINTAINER Sagnik Banerjee <sagnikbanerjee15@gmail.com>

ENV TZ=America/New_York
ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get -y update
RUN apt-get -y install libcurl4-doc libidn11-dev libkrb5-dev libldap2-dev librtmp-dev libssh2-1-dev git python3 curl libboost-all-dev apt-transport-https less vim wget time zlib1g zlib1g-dev lzma-dev libssl-dev libncurses5-dev libxml2-dev libxml2 liblzma-dev libncursesw5-dev make unzip zip build-essential gcc g++ cmake ca-certificates libbz2-dev xz-utils htop autoconf automake binutils bison flex gettext libtool make patch pkg-config dirmngr gnupg apt-transport-https ca-certificates software-properties-common r-base texlive-latex-base texlive-latex-extra python3-distutils trimmomatic
RUN wget https://bootstrap.pypa.io/get-pip.py && python3 get-pip.py
RUN pip install ruffus

RUN apt-get clean all

# Install Trimmomatic
RUN chmod a+x /usr/share/java/trimmomatic.jar
ENV PATH ${PATH}:/usr/share/java

# Install DESeq2
RUN R -e 'BiocManager::install("DESeq2", dependencies=TRUE)'

# Make directory for installation
RUN mkdir /software

# Install STAR
ARG STAR_VERSION=2.7.9a
RUN mkdir -p /software/STAR_${STAR_VERSION}
RUN cd /software/STAR_${STAR_VERSION} && \
	wget --no-check-certificate https://github.com/alexdobin/STAR/archive/${STAR_VERSION}.zip && \
	unzip ${STAR_VERSION}.zip && \
	cd STAR-${STAR_VERSION}/source && \
	make STAR STARlong

ENV PATH ${PATH}:/software/STAR_${STAR_VERSION}/STAR-${STAR_VERSION}/bin/Linux_x86_64

# Install Samtools
ARG SAMTOOLS_VERSION=1.14
RUN mkdir -p /software/samtools_${SAMTOOLS_VERSION}
RUN cd /software/samtools_${SAMTOOLS_VERSION} && \
	wget https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
	tar jxf samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
	cd samtools-${SAMTOOLS_VERSION} && make && make install



# Install salmon
ARG SALMON_VERSION=1.5.2
RUN mkdir -p /software/salmon_${SALMON_VERSION}
RUN cd /software/salmon_${SALMON_VERSION} && \
	wget https://github.com/COMBINE-lab/salmon/releases/download/v${SALMON_VERSION}/salmon-${SALMON_VERSION}_linux_x86_64.tar.gz && \
	tar xzf salmon-${SALMON_VERSION}_linux_x86_64.tar.gz && \
	cd salmon-${SALMON_VERSION}_linux_x86_64


ENV PATH ${PATH}:/software/salmon_${SALMON_VERSION}/salmon-${SALMON_VERSION}_linux_x86_64/bin

# Download and install NGPINT
ARG NGPINT_VERSION=1.0.0
RUN mkdir -p /software/ngpint_${NGPINT_VERSION} 
RUN cd /software/ngpint_${NGPINT_VERSION} && \
	git clone https://github.com/sagnikbanerjee15/NGPINT.git
RUN chmod a+x /software/ngpint_${NGPINT_VERSION}/NGPINT/ngpint

ENV PATH ${PATH}:/software/ngpint_${NGPINT_VERSION}/NGPINT

# Download gffread
ARG GFFREAD_VERSION=0.12.7
RUN mkdir -p /software/gffread_${GFFREAD_VERSION} # to force download
RUN cd /software/gffread_${GFFREAD_VERSION} && \
	wget --no-check-certificate https://github.com/gpertea/gffread/releases/download/v${GFFREAD_VERSION}/gffread-${GFFREAD_VERSION}.Linux_x86_64.tar.gz && \
	tar xzf gffread-${GFFREAD_VERSION}.Linux_x86_64.tar.gz && \
	chmod a+x gffread-${GFFREAD_VERSION}.Linux_x86_64/gffread

ENV PATH ${PATH}:/software/gffread_${GFFREAD_VERSION}/gffread-${GFFREAD_VERSION}.Linux_x86_64/
