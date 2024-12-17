
##############################################################
# Dockerfile Version:   1.6
# Software:             HOWARD
# Software Version:     0.12.1
# Software Website:     https://github.com/bioinfo-chru-strasbourg/howard
# Licence:              GNU Affero General Public License (AGPL)
# Description:          HOWARD - Highly Open Workflow for Annotation & Ranking toward genomic variant
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
	Version="0.12.1" \
	Website="https://github.com/bioinfo-chru-strasbourg/howard" \
	Description="HOWARD" \
	License="GNU Affero General Public License (AGPL)" \
    maintener="Antony Le Bechec <antony.lebechec@gmail.com>" \
	Usage='docker run -v $DATA FOLDER:/data -v $DATABASE_FOLDER:/databases -ti howard:0.12.1'



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

# YUM and Perl
ENV YUM_INSTALL="gcc bc make wget perl-devel which zlib-devel zlib bzip2-devel bzip2 xz-devel xz ncurses-devel unzip curl-devel java-17 htop libgomp aria2 docker-ce"
ENV PERL_INSTALL="perl-Switch perl-Time-HiRes perl-Data-Dumper perl-Digest-MD5 perl-Tk perl-devel"


###############
# YUM INSTALL #
###############

# AlmaLinux GPG key and YUM install and Perl install and bashrc
RUN echo "#[INFO] System YUM packages installation" && \
    rpm --import https://repo.almalinux.org/almalinux/RPM-GPG-KEY-AlmaLinux && \
    yum install -y yum-utils && \
    yum config-manager --add-repo=https://download.docker.com/linux/centos/docker-ce.repo && \
    yum install -y epel-release && \
    yum install -y $YUM_INSTALL && \
    yum install -y --enablerepo=powertools $PERL_INSTALL && \
    yum clean all && \
    echo 'alias ll="ls -lah"' >> ~/.bashrc;



################
# DEPENDENCIES #
################


##########
# HTSLIB #
##########

ENV TOOL_NAME=htslib
ENV TOOL_VERSION=1.19.1
ENV TARBALL_LOCATION=https://github.com/samtools/$TOOL_NAME/releases/download/$TOOL_VERSION/
ENV TARBALL=$TOOL_NAME-$TOOL_VERSION.tar.bz2
ENV DEST=$TOOLS/$TOOL_NAME/$TOOL_VERSION
ENV PATH=$DEST/bin:$PATH

# INSTALL
RUN echo "#[INFO] TOOL installation '$TOOL_NAME:$TOOL_VERSION'" && \
    wget $TARBALL_LOCATION/$TARBALL && \
    tar xf $TARBALL && \
    rm -rf $TARBALL && \
    cd $TOOL_NAME-$TOOL_VERSION && \
    make -j $THREADS prefix=$DEST install && \
    cd ../ && \
    ln -s $TOOL_VERSION $TOOLS/$TOOL_NAME/current && \
    rm -rf $TOOL_NAME-$TOOL_VERSION



############
# BCFTOOLS #
############

ENV TOOL_NAME=bcftools
ENV TOOL_VERSION=1.19
ENV TARBALL_LOCATION=https://github.com/samtools/$TOOL_NAME/releases/download/$TOOL_VERSION/
ENV TARBALL=$TOOL_NAME-$TOOL_VERSION.tar.bz2
ENV DEST=$TOOLS/$TOOL_NAME/$TOOL_VERSION
ENV PATH=$DEST/bin:$PATH

# INSTALL
RUN echo "#[INFO] TOOL installation '$TOOL_NAME:$TOOL_VERSION'" && \
    wget $TARBALL_LOCATION/$TARBALL && \
    tar xf $TARBALL && \
    rm -rf $TARBALL && \
    cd $TOOL_NAME-$TOOL_VERSION && \
    make -j $THREADS prefix=$DEST install && \
    cd ../ && \
    ln -s $TOOL_VERSION $TOOLS/$TOOL_NAME/current && \
    rm -rf $TOOL_NAME-$TOOL_VERSION



##########
# SNPEFF #
##########

ENV TOOL_NAME=snpeff
ENV TOOL_VERSION=5.2a
ENV TARBALL="snpEff_v5_2a_core.zip"
ENV TARBALL_LOCATION=https://snpeff.blob.core.windows.net/versions
ENV TARBALL_FOLDER=snpeff_folder
ENV TOOL_DATABASE_FOLDER=$DATABASES/snpeff/current
ENV DEST=$TOOLS/$TOOL_NAME/$TOOL_VERSION
ENV PATH=$DEST/bin:$PATH
ENV SNPEFF_DATABASES=$TOOL_DATABASE_FOLDER

# INSTALL
RUN echo "#[INFO] TOOL installation '$TOOL_NAME:$TOOL_VERSION'" && \
    wget $TARBALL_LOCATION/$TARBALL && \
    unzip $TARBALL -d $TARBALL_FOLDER && \
    mkdir -p $DEST/bin/ && \
    cp $TARBALL_FOLDER/*/*jar $DEST/bin/ -R && \
    cp $TARBALL_FOLDER/*/*.config $DEST/bin/ -R && \
    ln -s $TOOL_VERSION $TOOLS/$TOOL_NAME/current && \
    mkdir -p $TOOL_DATABASE_FOLDER && \
    ln -s $TOOL_DATABASE_FOLDER $DEST/bin/data && \
    rm -rf $TARBALL $TARBALL_FOLDER;



###########
# ANNOVAR #
###########

ENV TOOL_NAME=annovar
ENV TOOL_VERSION=2020Jun08
ENV TARBALL_LOCATION=http://www.openbioinformatics.org/annovar/download/0wgxR2rIVP
ENV TARBALL=annovar.latest.tar.gz
ENV TARBALL_FOLDER=$TOOL_NAME
ENV TOOL_DATABASE_FOLDER=$DATABASES/annovar/current
ENV DEST=$TOOLS/$TOOL_NAME/$TOOL_VERSION
ENV PATH=$DEST/bin:$PATH
ENV ANNOVAR_DATABASES=$TOOL_DATABASE_FOLDER

# INSTALL
RUN echo "#[INFO] TOOL installation '$TOOL_NAME:$TOOL_VERSION'" && \
    wget $TARBALL_LOCATION/$TARBALL && \
    tar xf $TARBALL && \
    rm -rf $TARBALL && \
    cd $TARBALL_FOLDER && \
    mkdir -p $DEST/bin && \
    cp *.pl $DEST/bin/ -R && \
    cd ../ && \
    rm -rf $TARBALL_FOLDER && \
    ln -s $TOOL_VERSION $TOOLS/$TOOL_NAME/current && \
    mkdir -p $TOOL_DATABASE_FOLDER && \
    ln -s $TOOL_DATABASE_FOLDER $DEST/bin/databases ;


##############
# MICROMAMBA #
##############

ENV TOOL_NAME=micromamba
ENV TOOL_VERSION=2.0.5-0
ENV DEST=$TOOLS/$TOOL_NAME/$TOOL_VERSION
ENV MICROMAMBA=$DEST/bin/micromamba
ENV PATH=$TOOLS/$TOOL_NAME/current/bin:$PATH

RUN echo "#[INFO] TOOL installation '$TOOL_NAME:$TOOL_VERSION'" && \
    mkdir -p $DEST/bin && \
    wget "micro.mamba.pm/install.sh" && \
    chmod 0755 "install.sh" && \
    VERSION=$TOOL_VERSION && yes | bash install.sh && \
    ln -s ~/.local/bin/micromamba $MICROMAMBA



##########
# PYTHON #
##########

ENV TOOL_NAME=python
ENV PATH=$TOOLS/$TOOL_NAME/current/bin:$PATH

# PYTHON 3.10 - current
ENV TOOL_NAME=python
ENV TOOL_VERSION=3.10
ENV DEST=$TOOLS/$TOOL_NAME/$TOOL_VERSION

# INSTALL
RUN echo "#[INFO] TOOL installation '$TOOL_NAME:$TOOL_VERSION'" && \
    $MICROMAMBA create python=$TOOL_VERSION -p $TOOLS/$TOOL_NAME/$TOOL_VERSION && \
    $MICROMAMBA clean --all && \
    cd $TOOLS/$TOOL_NAME && \
    ln -s $TOOL_VERSION current



############
# EXOMISER #
############

ENV TOOL_NAME=exomiser
ENV TOOL_VERSION=14.0.0
ENV TARBALL_LOCATION=https://data.monarchinitiative.org/exomiser/$TOOL_VERSION
ENV TARBALL=exomiser-cli-$TOOL_VERSION-distribution.zip
ENV TARBALL_FOLDER=exomiser-cli-$TOOL_VERSION
ENV DEST=$TOOLS/$TOOL_NAME/$TOOL_VERSION
ENV PATH=$DEST/bin:$PATH

# INSTALL
RUN echo "#[INFO] TOOL installation '$TOOL_NAME:$TOOL_VERSION'" && \
    wget $TARBALL_LOCATION/$TARBALL && \
    mkdir -p $DEST && \
    unzip $TARBALL -d $DEST && \
    mv $DEST/$TARBALL_FOLDER $DEST/bin && \
    rm -rf $TARBALL && \
    ln -s $TOOL_VERSION $TOOLS/$TOOL_NAME/current



###########
# HOWARD #
###########

ENV TOOL_NAME=howard
ENV TOOL_VERSION=0.12.1
ENV DEST=$TOOLS/$TOOL_NAME/$TOOL_VERSION
ENV PATH=$DEST/bin:$PATH
ENV USER_HOME=/root
ENV HOWARD_HOME=$USER_HOME/howard

# INSTALL
ADD . $DEST
RUN echo "#[INFO] TOOL installation '$TOOL_NAME:$TOOL_VERSION'" && \
    ln -s $DEST $TOOLS/$TOOL_NAME/current && \
    chmod 0775 $DEST $TOOLS/$TOOL_NAME/current -R && \
	mkdir -p $DATABASES && \
	ln -s $TOOL_VERSION $TOOLS/$TOOL_NAME/current && \
	ln -s $DEST/ /tool && \
    (cd /tool && python -m pip install -e .) && \
    mkdir -p $DEST/bin && \
    cp $(whereis howard | cut -d" " -f2) $DEST/bin/ && \
    mkdir -p $HOWARD_HOME && \
    ln -s $TOOLS $HOWARD_HOME$TOOLS && \
    ln -s $DATABASES $HOWARD_HOME$DATABASES && \
    ln -s $DATA $HOWARD_HOME$DATA && \
    python -m pip cache purge && \
    rm -rf $DEST/.git* $DEST/howard.*


##############################
# WORKDIR / ENTRYPOINT / CMD #
##############################


WORKDIR "/data"

ENTRYPOINT [ "howard" ]
