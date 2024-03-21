
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

# YUM
ENV YUM_INSTALL="gcc bc make wget perl-devel which zlib-devel zlib bzip2-devel bzip2 xz-devel xz ncurses-devel unzip curl-devel java-17 htop libgomp aria2 docker"
#ENV YUM_REMOVE="zlib-devel bzip2-devel xz-devel ncurses-devel gcc"



##########
# BASHRC #
##########

RUN echo 'alias ll="ls -lah"' >> ~/.bashrc



###############
# YUM INSTALL #
###############

# AlmaLinux GPG key and YUM install
RUN rpm --import https://repo.almalinux.org/almalinux/RPM-GPG-KEY-AlmaLinux && \
    yum install -y epel-release && \
    yum install -y $YUM_INSTALL;



################
# DEPENDENCIES #
################



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
#ENV TOOL_VERSION=1.15.1
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
#ENV TOOL_VERSION=1.15.1
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
#ENV TOOL_VERSION=5.1d
ENV TOOL_VERSION=5.2a
#ENV TARBALL="snpEff_latest_core.zip"
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
    ln -s $TOOL_DATABASE_FOLDER $DEST/bin/data ;



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



#########
# MAMBA #
#########

ENV TOOL_NAME=mamba
ENV TOOL_VERSION=23.3.1-1
ENV TARBALL_LOCATION=https://github.com/conda-forge/miniforge/releases/download/$TOOL_VERSION
ENV DEST=$TOOLS/$TOOL_NAME/$TOOL_VERSION
ENV MAMBA=$DEST/bin/mamba
ENV PIP=$DEST/bin/pip

# INSTALL
RUN echo "#[INFO] TOOL installation '$TOOL_NAME:$TOOL_VERSION'" && \
    wget "$TARBALL_LOCATION/Mambaforge-$TOOL_VERSION-$(uname)-$(uname -m).sh" && \
    bash "Mambaforge-$TOOL_VERSION-$(uname)-$(uname -m).sh" -b -p $DEST && \
    rm -f "Mambaforge-$TOOL_VERSION-$(uname)-$(uname -m).sh"
    


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
    $MAMBA create python=$TOOL_VERSION -p $TOOLS/$TOOL_NAME/$TOOL_VERSION && \
    cd $TOOLS/$TOOL_NAME && \
    ln -s $TOOL_VERSION current



############
# EXOMISER #
############

ENV TOOL_NAME=exomiser
#ENV TOOL_VERSION=13.2.0
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
ENV TOOL_VERSION=devel
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
    ln -s $DATA $HOWARD_HOME$DATA



######################
# YUM REMOVE & CLEAR #
######################

RUN [ "${YUM_REMOVE}" != "" ] && yum erase -y $YUM_REMOVE || echo "Nothing to clean"
RUN yum clean all



##############################
# WORKDIR / ENTRYPOINT / CMD #
##############################


WORKDIR "/data"

ENTRYPOINT [ "howard" ]
