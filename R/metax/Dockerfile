FROM centos/systemd
MAINTAINER Peter Li <peter@gigasciencejournal.com>

# Minimal installation
RUN yum install -y epel-release
RUN yum update -y
RUN yum upgrade -y
RUN yum clean all
EXPOSE 80

# R devtools pre-requisites:
RUN yum install -y git xml2 libxml2-devel curl curl-devel openssl-devel pandoc boost-devel ed
# Required by RCurl R package
RUN yum install -y libcurl-devel netcdf-devel

# Required by rgl R package
RUN	yum install -y libpng-devel libX11-devel mesa-libGL mesa-libGL-devel mesa-libGLU mesa-libGLU-devel

# Required to download metaX
RUN	yum install -y wget unzip ImageMagick graphviz cairo-devel libXt-devel

WORKDIR /tmp
RUN yum install -y R

RUN R -e "install.packages('Nozzle.R1', repos='https://mirror.lzu.edu.cn/CRAN/')"
RUN R -e "install.packages('ggplot2', repos='https://mirror.lzu.edu.cn/CRAN/')"

# metaX says it requires parallel but problem with its installation
# 1: package 'parallel' is not available (for R version 3.4.1)
# 2: package 'parallel' is a base package, and should not be updated
# RUN R -e "install.packages('parallel', repos='https://mirror.lzu.edu.cn/CRAN/')"

RUN R -e "install.packages('reshape2', repos='https://mirror.lzu.edu.cn/CRAN/')"
RUN R -e "install.packages('plyr', repos='https://mirror.lzu.edu.cn/CRAN/')"
RUN R -e "install.packages('BBmisc', repos='https://mirror.lzu.edu.cn/CRAN/')"
RUN R -e "source('http://bioconductor.org/biocLite.R'); biocLite('rgl')"
RUN R -e "install.packages('devtools', repos='https://mirror.lzu.edu.cn/CRAN/')"
RUN R -e "library(devtools); install_github('gaborcsardi/pkgconfig'); install_github('igraph/rigraph')"
RUN R -e "install.packages('mixOmics', repos='https://mirror.lzu.edu.cn/CRAN/')"
RUN R -e "install.packages('missForest', repos='https://mirror.lzu.edu.cn/CRAN/')"
RUN R -e "install.packages('doParallel', repos='https://mirror.lzu.edu.cn/CRAN/')"
RUN R -e "install.packages('DiscriMiner', repos='https://mirror.lzu.edu.cn/CRAN/')"
RUN R -e "install.packages('ape', repos='https://mirror.lzu.edu.cn/CRAN/')"
RUN R -e "install.packages('scatterplot3d', repos='https://mirror.lzu.edu.cn/CRAN/')"
RUN R -e "install.packages('pheatmap', repos='https://mirror.lzu.edu.cn/CRAN/')"
RUN R -e "install.packages('bootstrap', repos='https://mirror.lzu.edu.cn/CRAN/')"
RUN R -e "install.packages('boot', repos='https://mirror.lzu.edu.cn/CRAN/')"
RUN R -e "install.packages('caret', repos='https://mirror.lzu.edu.cn/CRAN/')"
RUN R -e "install.packages('dplyr', repos='https://mirror.lzu.edu.cn/CRAN/')"
RUN R -e "install.packages('stringr', repos='https://mirror.lzu.edu.cn/CRAN/')"
RUN R -e "install.packages('RColorBrewer', repos='https://mirror.lzu.edu.cn/CRAN/')"
RUN R -e "install.packages('RCurl', repos='https://mirror.lzu.edu.cn/CRAN/')"
RUN R -e "install.packages('lattice', repos='https://mirror.lzu.edu.cn/CRAN/')"
RUN R -e "install.packages('data.table', repos='https://mirror.lzu.edu.cn/CRAN/')"
RUN R -e "install.packages('tidyr', repos='https://mirror.lzu.edu.cn/CRAN/')"
RUN R -e "install.packages('scales', repos='https://mirror.lzu.edu.cn/CRAN/')"
RUN R -e "install.packages('VennDiagram', repos='https://mirror.lzu.edu.cn/CRAN/')"
RUN R -e "install.packages('pROC', repos='https://mirror.lzu.edu.cn/CRAN/')"
RUN R -e "install.packages('readr', repos='https://mirror.lzu.edu.cn/CRAN/')"
RUN R -e "install.packages('e1071', repos='https://mirror.lzu.edu.cn/CRAN/')"
RUN R -e "install.packages('randomForest', repos='https://mirror.lzu.edu.cn/CRAN/')"
RUN R -e "install.packages('coop', repos='https://mirror.lzu.edu.cn/CRAN/')"
RUN R -e "install.packages('fpc', repos='https://mirror.lzu.edu.cn/CRAN/')"
RUN R -e "install.packages('ROCR', repos='https://mirror.lzu.edu.cn/CRAN/')"
RUN R -e "source('http://bioconductor.org/biocLite.R'); biocLite(c('pcaMethods', 'multtest'))"
RUN R -e "install.packages('DiffCorr', repos='https://mirror.lzu.edu.cn/CRAN/')"
RUN R -e "source('http://bioconductor.org/biocLite.R'); biocLite('preprocessCore')"
RUN R -e "source('http://bioconductor.org/biocLite.R'); biocLite('vsn')"
RUN R -e "source('http://bioconductor.org/biocLite.R'); biocLite('pls')"
RUN R -e "install.packages('XML', repos='https://mirror.lzu.edu.cn/CRAN/')"
RUN R -e "source('http://bioconductor.org/biocLite.R'); biocLite('mzID')"
RUN R -e "source('http://bioconductor.org/biocLite.R'); biocLite('mzR')"
RUN R -e "source('http://bioconductor.org/biocLite.R'); biocLite('MSnbase')"
RUN R -e "source('http://bioconductor.org/biocLite.R'); biocLite('faahKO')"
RUN R -e "source('http://bioconductor.org/biocLite.R'); biocLite('xcms')"
RUN R -e "source('http://bioconductor.org/biocLite.R'); biocLite('impute')"
RUN R -e "source('http://bioconductor.org/biocLite.R'); biocLite('CAMERA')"
RUN R -e "source('http://bioconductor.org/biocLite.R'); biocLite('ropls')"
RUN R -e "install.packages('RSQLite', repos='https://mirror.lzu.edu.cn/CRAN/')"
RUN R -e "source('http://bioconductor.org/biocLite.R'); biocLite('sva')"
RUN R -e "source('http://bioconductor.org/biocLite.R'); biocLite('SSPA')"

RUN R -e "install.packages('matrixStats', repos='https://mirror.lzu.edu.cn/CRAN/')"
RUN R -e "install.packages('naturalsort', repos='https://mirror.lzu.edu.cn/CRAN/')"

RUN R -e "install.packages('mice', repos='https://mirror.lzu.edu.cn/CRAN/')"
RUN R -e "install.packages('ggfortify', repos='https://mirror.lzu.edu.cn/CRAN/')"

# Download metaX
RUN wget --directory-prefix=/tmp http://metax.genomics.cn/files/metaX_1.4.16.tar.gz


RUN export R_MAX_NUM_DLLS=150

# Install metaX
RUN R -e "install.packages('/tmp/metaX_1.4.16.tar.gz', repos=NULL, type='source')"


