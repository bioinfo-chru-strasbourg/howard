
##############################################################
# Dockerfile Version:   1.3
# Software:             HOWARD
# Software Version:     0.9.14b
# Software Website:     https://gitlab.bioinfo-diag.fr/Strasbourg/HOWARD
# Licence:              GNU Affero General Public License (AGPL)
# Description:          HOWARD
# Usage:                docker run -ti [-v [DATA FOLDER]:/data -v [DATABASE_FOLDER]:/databases] howard:version
##############################################################

##########
# README #
##########

# Config parameters
#    identify yum packages for installation
#    identify yum packages to remove
#
# Dependecies installation
#    identify tools dependences
#    config each tool
#    write installation procedure for each tools
#
# Tool
#    configure tool
#    write isntallation procedure for the tool
#    add link to current and root tool folder
#
# Workdir / Entrypoint / Cmd
#    configure workdir, endpoint and command
#    /!\ no variables in endpoint



########
# FROM #
########

FROM centos:7
MAINTAINER Antony Le Bechec <antony.lebechec@gmail.com>
LABEL Software="HOWARD" \
	Version="0.9.14b" \
	Website="https://gitlab.bioinfo-diag.fr/Strasbourg/HOWARD" \
	Description="HOWARD" \
	License="GNU Affero General Public License (AGPL)" \
	Usage="docker run -ti [-v [DATA FOLDER]:/data -v [DATABASE_FOLDER]:/databases] howard:version"



##############
# PARAMETERS #
##############

ENV TOOLS=/home/TOOLS/tools
ENV DATA=/data
ENV TOOL=/tool
ENV DATABASES=/databases
ENV YUM_INSTALL="gcc bc make wget perl-Switch perl-Digest-MD5 perl-Data-Dumper which zlib-devel zlib 	zlib2-devel zlib2 	bzip2-devel bzip2 	lzma-devel lzma 	xz-devel xz 	ncurses-devel 	unzip"
ENV YUM_REMOVE="zlib-devel zlib2-devel bzip2-devel lzma-devel xz-devel ncurses-devel unzip gcc"



###############
# YUM INSTALL #
###############

RUN yum install -y $YUM_INSTALL ;



################
# DEPENDENCIES #
################


##########
# HTSLIB #
##########

ENV TOOL_NAME=htslib
ENV TOOL_VERSION=1.8
ENV TARBALL_LOCATION=https://github.com/samtools/$TOOL_NAME/releases/download/$TOOL_VERSION/
ENV TARBALL=$TOOL_NAME-$TOOL_VERSION.tar.bz2
ENV DEST=$TOOLS/$TOOL_NAME/$TOOL_VERSION
ENV PATH=$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin:$PATH

# INSTALL
RUN wget $TARBALL_LOCATION/$TARBALL ; \
    tar xf $TARBALL ; \
    rm -rf $TARBALL ; \
    cd $TOOL_NAME-$TOOL_VERSION ; \
    make prefix=$TOOLS/$TOOL_NAME/$TOOL_VERSION install ; \
    cd ../ ; \
    rm -rf $TOOL_NAME-$TOOL_VERSION



############
# BCFTOOLS #
############

ENV TOOL_NAME=bcftools
ENV TOOL_VERSION=1.8
ENV TARBALL_LOCATION=https://github.com/samtools/$TOOL_NAME/releases/download/$TOOL_VERSION/
ENV TARBALL=$TOOL_NAME-$TOOL_VERSION.tar.bz2
ENV DEST=$TOOLS/$TOOL_NAME/$TOOL_VERSION
ENV PATH=$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin:$PATH

# INSTALL
RUN wget $TARBALL_LOCATION/$TARBALL ; \
    tar xf $TARBALL ; \
    rm -rf $TARBALL ; \
    cd $TOOL_NAME-$TOOL_VERSION ; \
    make prefix=$TOOLS/$TOOL_NAME/$TOOL_VERSION install ; \
    cd ../ ; \
    rm -rf $TOOL_NAME-$TOOL_VERSION


########
# JAVA #
########

ENV TOOL_NAME=java
ENV TOOL_VERSION=1.8.0
RUN yum install -y java-$TOOL_VERSION && \
	mkdir -p $TOOLS/$TOOL_NAME/$TOOL_VERSION/bin && \
	ln -s /usr/bin/java $TOOLS/$TOOL_NAME/$TOOL_VERSION/bin/java && \
	ln -s $TOOLS/$TOOL_NAME/$TOOL_VERSION $TOOLS/$TOOL_NAME/current ;



##########
# SNPEFF #
##########

ENV DATABASES=/databases
ENV TOOL_NAME=snpeff
ENV TOOL_VERSION=4.3t
ENV TOOL_VERSION_FOR_FILE=4_3t
ENV TARBALL_LOCATION=https://sourceforge.net/projects/snpeff/files
ENV TARBALL="snpEff_v"$TOOL_VERSION_FOR_FILE"_core.zip"
ENV TARBALL_FOLDER=snpeff_folder
ENV TOOL_DATABASE_FOLDER=$DATABASES/snpeff/$TOOL_VERSION
ENV DEST=$TOOLS/$TOOL_NAME/$TOOL_VERSION
ENV PATH=$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin:$PATH


# INSTALL
RUN wget $TARBALL_LOCATION/$TARBALL ; \
    unzip $TARBALL -d $TARBALL_FOLDER ; \
    mkdir -p $TOOLS/$TOOL_NAME/$TOOL_VERSION/bin/ ; \
    cp $TARBALL_FOLDER/*/*jar $TOOLS/$TOOL_NAME/$TOOL_VERSION/bin/ -R ; \
    cp $TARBALL_FOLDER/*/*.config $TOOLS/$TOOL_NAME/$TOOL_VERSION/bin/ -R ; \
    ln -s $TOOLS/$TOOL_NAME/$TOOL_VERSION $TOOLS/$TOOL_NAME/current ; \
    mkdir -p $TOOL_DATABASE_FOLDER ; \
    ln -s $TOOL_DATABASE_FOLDER $TOOLS/$TOOL_NAME/$TOOL_VERSION/bin/data ;



###########
# ANNOVAR #
###########

ENV DATABASES=/databases
ENV TOOL_NAME=annovar
ENV TOOL_VERSION=2018Apr16
ENV TARBALL_LOCATION=http://www.openbioinformatics.org/annovar/download/0wgxR2rIVP/
ENV TARBALL=annovar.latest.tar.gz
ENV TARBALL_FOLDER=$TOOL_NAME
ENV TOOL_DATABASE_FOLDER=$DATABASES/annovar
ENV DEST=$TOOLS/$TOOL_NAME/$TOOL_VERSION
ENV PATH=$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin:$PATH
# http://www.openbioinformatics.org/annovar/download/0wgxR2rIVP/annovar.latest.tar.gz


# INSTALL
RUN wget $TARBALL_LOCATION/$TARBALL ; \
    tar xf $TARBALL ; \
    rm -rf $TARBALL ; \
    cd $TARBALL_FOLDER ; \
    mkdir -p $TOOLS/$TOOL_NAME/$TOOL_VERSION/bin ; \
    cp *.pl $TOOLS/$TOOL_NAME/$TOOL_VERSION/bin/ -R ; \
    cd ../ ; \
    rm -rf $TARBALL_FOLDER ; \
    ln -s $TOOLS/$TOOL_NAME/$TOOL_VERSION $TOOLS/$TOOL_NAME/current ; \
    mkdir -p $TOOL_DATABASE_FOLDER ; \
    ln -s $TOOL_DATABASE_FOLDER $TOOLS/$TOOL_NAME/$TOOL_VERSION/bin/databases ;



###########
# HOWARD #
###########

ENV DATABASES=/databases
ENV TOOL_NAME=howard
ENV TOOL_VERSION=0.9.14b
ENV TARBALL_LOCATION=https://gitlab.bioinfo-diag.fr/Strasbourg/HOWARD/repository/$TOOL_VERSION
ENV TARBALL=archive.tar.gz
ENV TARBALL_FOLDER=archive
ENV TOOL_DATABASE_FOLDER=/home/TOOLS/databases
ENV DEST=$TOOLS/$TOOL_NAME/$TOOL_VERSION
ENV PATH=$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin:$PATH


RUN wget $TARBALL_LOCATION/$TARBALL ; \
    tar xf $TARBALL ; \
    rm -rf $TARBALL ; \
    mkdir -p $TOOLS/$TOOL_NAME/$TOOL_VERSION/ ; \
    cp $(ls ${TOOL_NAME^^}-$TOOL_VERSION* -d)/* $TOOLS/$TOOL_NAME/$TOOL_VERSION/ -R ; \
    rm -rf $(ls ${TOOL_NAME^^}-$TOOL_VERSION* -d) ; \
    ln -s $TOOLS/$TOOL_NAME/$TOOL_VERSION $TOOLS/$TOOL_NAME/current ; \
    chmod 0775 $TOOLS/$TOOL_NAME/$TOOL_VERSION $TOOLS/$TOOL_NAME/current -R ; \
	mkdir -p $DATABASES ; \
	ln -s $DATABASES $TOOL_DATABASE_FOLDER ; \
	ln -s $TOOLS/$TOOL_NAME/$TOOL_VERSION $TOOLS/$TOOL_NAME/current ; \
	ln -s $TOOLS/$TOOL_NAME/$TOOL_VERSION/ /tool ;



######################
# YUM REMOVE & CLEAR #
######################

RUN yum erase -y $YUM_REMOVE ; yum clean all ;



##############################
# WORKDIR / ENTRYPOINT / CMD #
##############################


WORKDIR "/data"

ENTRYPOINT [ "/tool/bin/HOWARD.sh" ]
