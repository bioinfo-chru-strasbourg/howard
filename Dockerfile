
##############################################################
# Dockerfile Version:   1.5
# Software:             HOWARD
# Software Version:     0.9.15.5
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

FROM almalinux:8
LABEL Software="HOWARD" \
	Version="1.0.0" \
	Website="https://github.com/bioinfo-chru-strasbourg/howard" \
	Description="HOWARD" \
	License="GNU Affero General Public License (AGPL)" \
    maintener="Antony Le Bechec <antony.lebechec@gmail.com>" \
	Usage='docker run -v $DATA FOLDER:/data -v $DATABASE_FOLDER:/databases -ti howard:1.0.0'



########
# ARGS #
########

ARG THREADS="8"



##############
# PARAMETERS #
##############

ENV TOOLS=/tools
ENV DATA=/data
ENV TOOL=/tool
ENV DATABASES=/databases
ENV YUM_INSTALL="gcc bc make wget perl-devel which zlib-devel zlib bzip2-devel bzip2 xz-devel xz ncurses-devel unzip curl-devel python39 java-17 htop"
ENV YUM_REMOVE="zlib-devel bzip2-devel xz-devel ncurses-devel gcc"



###############
# YUM INSTALL #
###############

RUN yum install -y epel-release;
RUN yum install -y $YUM_INSTALL;



################
# DEPENDENCIES #
################




###########
# PYTHON3 #
###########

ENV TOOL_NAME="python"
ENV TOOL_VERSION="3"
ENV PATH=$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin:$PATH
RUN echo "#[INFO] TOOL installation '$TOOL_NAME:$TOOL_VERSION'" && \
	ln -s /usr/bin/python$TOOL_VERSION /usr/bin/python; 





########
# PERL #
########

# PERL installation
RUN	echo "#[INFO] System Perl packages installation"; \
	yum -y --enablerepo=powertools install perl perl-Switch perl-Time-HiRes perl-Data-Dumper perl-Digest-MD5 perl-Tk perl-devel




##########
# HTSLIB #
##########

ENV TOOL_NAME=htslib
ENV TOOL_VERSION=1.15.1
ENV TARBALL_LOCATION=https://github.com/samtools/$TOOL_NAME/releases/download/$TOOL_VERSION/
ENV TARBALL=$TOOL_NAME-$TOOL_VERSION.tar.bz2
ENV DEST=$TOOLS/$TOOL_NAME/$TOOL_VERSION
ENV PATH=$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin:$PATH

# INSTALL
RUN echo "#[INFO] TOOL installation '$TOOL_NAME:$TOOL_VERSION'" && \
    wget $TARBALL_LOCATION/$TARBALL && \
    tar xf $TARBALL && \
    rm -rf $TARBALL && \
    cd $TOOL_NAME-$TOOL_VERSION && \
    make -j $THREADS prefix=$TOOLS/$TOOL_NAME/$TOOL_VERSION install && \
    cd ../ && \
    ln -s $TOOL_VERSION $TOOLS/$TOOL_NAME/current && \
    rm -rf $TOOL_NAME-$TOOL_VERSION



############
# BCFTOOLS #
############

ENV TOOL_NAME=bcftools
ENV TOOL_VERSION=1.15.1
ENV TARBALL_LOCATION=https://github.com/samtools/$TOOL_NAME/releases/download/$TOOL_VERSION/
ENV TARBALL=$TOOL_NAME-$TOOL_VERSION.tar.bz2
ENV DEST=$TOOLS/$TOOL_NAME/$TOOL_VERSION
ENV PATH=$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin:$PATH

# INSTALL
RUN echo "#[INFO] TOOL installation '$TOOL_NAME:$TOOL_VERSION'" && \
    wget $TARBALL_LOCATION/$TARBALL && \
    tar xf $TARBALL && \
    rm -rf $TARBALL && \
    cd $TOOL_NAME-$TOOL_VERSION && \
    make -j $THREADS prefix=$TOOLS/$TOOL_NAME/$TOOL_VERSION install && \
    cd ../ && \
    ln -s $TOOL_VERSION $TOOLS/$TOOL_NAME/current && \
    rm -rf $TOOL_NAME-$TOOL_VERSION


########
# JAVA #
########

# ENV TOOL_NAME=java
# ENV TOOL_VERSION=1.8.0
# RUN yum install -y java-$TOOL_VERSION && \
# 	mkdir -p $TOOLS/$TOOL_NAME/$TOOL_VERSION/bin && \
# 	ln -s /usr/bin/java $TOOLS/$TOOL_NAME/$TOOL_VERSION/bin/java && \
# 	ln -s $TOOL_VERSION $TOOLS/$TOOL_NAME/current ;



##########
# SNPEFF #
##########

ENV DATABASES=/databases
ENV TOOL_NAME=snpeff
ENV TOOL_VERSION=5.1d
ENV TARBALL="snpEff_latest_core.zip"
ENV TARBALL_LOCATION=https://snpeff.blob.core.windows.net/versions
ENV TARBALL_FOLDER=snpeff_folder
ENV TOOL_DATABASE_FOLDER=$DATABASES/snpeff/current
ENV DEST=$TOOLS/$TOOL_NAME/$TOOL_VERSION
ENV PATH=$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin:$PATH
ENV SNPEFF_DATABASES=$TOOL_DATABASE_FOLDER


# INSTALL
RUN echo "#[INFO] TOOL installation '$TOOL_NAME:$TOOL_VERSION'" && \
    wget $TARBALL_LOCATION/$TARBALL && \
    unzip $TARBALL -d $TARBALL_FOLDER && \
    mkdir -p $TOOLS/$TOOL_NAME/$TOOL_VERSION/bin/ && \
    cp $TARBALL_FOLDER/*/*jar $TOOLS/$TOOL_NAME/$TOOL_VERSION/bin/ -R && \
    cp $TARBALL_FOLDER/*/*.config $TOOLS/$TOOL_NAME/$TOOL_VERSION/bin/ -R && \
    ln -s $TOOL_VERSION $TOOLS/$TOOL_NAME/current && \
    mkdir -p $TOOL_DATABASE_FOLDER && \
    ln -s $TOOL_DATABASE_FOLDER $TOOLS/$TOOL_NAME/$TOOL_VERSION/bin/data ;



###########
# ANNOVAR #
###########

ENV DATABASES=/databases
ENV TOOL_NAME=annovar
ENV TOOL_VERSION=2020Jun08
ENV TARBALL_LOCATION=http://www.openbioinformatics.org/annovar/download/0wgxR2rIVP
ENV TARBALL=annovar.latest.tar.gz
ENV TARBALL_FOLDER=$TOOL_NAME
ENV TOOL_DATABASE_FOLDER=$DATABASES/annovar/current
ENV DEST=$TOOLS/$TOOL_NAME/$TOOL_VERSION
ENV PATH=$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin:$PATH
ENV ANNOVAR_DATABASES=$TOOL_DATABASE_FOLDER
# http://www.openbioinformatics.org/annovar/download/0wgxR2rIVP/annovar.latest.tar.gz


# INSTALL
RUN echo "#[INFO] TOOL installation '$TOOL_NAME:$TOOL_VERSION'" && \
    wget $TARBALL_LOCATION/$TARBALL && \
    tar xf $TARBALL && \
    rm -rf $TARBALL && \
    cd $TARBALL_FOLDER && \
    mkdir -p $TOOLS/$TOOL_NAME/$TOOL_VERSION/bin && \
    cp *.pl $TOOLS/$TOOL_NAME/$TOOL_VERSION/bin/ -R && \
    cd ../ && \
    rm -rf $TARBALL_FOLDER && \
    ln -s $TOOL_VERSION $TOOLS/$TOOL_NAME/current && \
    mkdir -p $TOOL_DATABASE_FOLDER && \
    ln -s $TOOL_DATABASE_FOLDER $TOOLS/$TOOL_NAME/$TOOL_VERSION/bin/databases ;



###########
# HOWARD #
###########

ENV DATABASES=/databases
ENV TOOL_NAME=howard
ENV TOOL_VERSION=devel
#ENV TARBALL_LOCATION=https://github.com/bioinfo-chru-strasbourg/howard/repository/$TOOL_VERSION
ENV TARBALL=archive.tar.gz
ENV TARBALL_FOLDER=archive
ENV DEST=$TOOLS/$TOOL_NAME/$TOOL_VERSION
ENV PATH=$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin:$PATH


ADD . $TOOLS/$TOOL_NAME/$TOOL_VERSION

RUN echo "#[INFO] TOOL installation '$TOOL_NAME:$TOOL_VERSION'" && \
    ln -s $TOOLS/$TOOL_NAME/$TOOL_VERSION $TOOLS/$TOOL_NAME/current && \
    chmod 0775 $TOOLS/$TOOL_NAME/$TOOL_VERSION $TOOLS/$TOOL_NAME/current -R && \
	mkdir -p $DATABASES && \
	ln -s $TOOL_VERSION $TOOLS/$TOOL_NAME/current && \
	ln -s $TOOLS/$TOOL_NAME/$TOOL_VERSION/ /tool && \
    (cd /tool && python -m pip install -e .) ;



######################
# YUM REMOVE & CLEAR #
######################

RUN yum erase -y $YUM_REMOVE ; yum clean all ;



##############################
# WORKDIR / ENTRYPOINT / CMD #
##############################


WORKDIR "/data"

ENTRYPOINT [ "howard" ]
