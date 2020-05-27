
##############################################################
# Dockerfile Version:   1.4.1
# Software:             HOWARD
# Software Version:     0.9.15.1b
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

FROM base_howard:latest
MAINTAINER Antony Le Bechec <antony.lebechec@gmail.com>
LABEL Software="HOWARD" \
	Version="0.9.15b" \
	Website="https://gitlab.bioinfo-diag.fr/Strasbourg/HOWARD" \
	Description="HOWARD" \
	License="GNU Affero General Public License (AGPL)" \
	Usage="docker run -ti [-v [DATA FOLDER]:/data -v [DATABASE_FOLDER]:/databases] howard:version"



##############
# PARAMETERS #
##############

ENV TOOLS=/tools
ENV DATA=/data
ENV TOOL=/tool
ENV DATABASES=/databases
ENV YUM_INSTALL="gcc bc make wget perl-Switch perl-Digest-MD5 perl-Data-Dumper which zlib-devel zlib zlib2-devel zlib2 bzip2-devel bzip2 lzma-devel lzma xz-devel xz ncurses-devel unzip"
ENV YUM_REMOVE="zlib-devel bzip2-devel xz-devel ncurses-devel unzip gcc"



###############
# YUM INSTALL #
###############

RUN yum install -y $YUM_INSTALL ;



################
# DEPENDENCIES #
################




###########
# HOWARD #
###########

ENV DATABASES=/databases
ENV TOOL_NAME=howard
ENV TOOL_VERSION=0.9.15.1b
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

ENTRYPOINT [ "/tool/bin/HOWARD" ]
