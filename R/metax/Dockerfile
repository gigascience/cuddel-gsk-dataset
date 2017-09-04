FROM centos/systemd
MAINTAINER Peter Li <peter@gigasciencejournal.com>

# Minimal installation
RUN yum install -y epel-release
RUN yum update -y
RUN yum upgrade -y
RUN yum clean all
EXPOSE 80

# R devtools pre-requisites:
RUN yum install -y git xml2 libxml2-devel curl curl-devel openssl-devel pandoc
# Required by RCurl R package
RUN yum install -y libcurl-devel netcdf-devel

# Required by rgl R package
RUN	yum install -y libpng-devel libX11-devel mesa-libGL mesa-libGL-devel mesa-libGLU mesa-libGLU-devel

# Required to download metaX
RUN	yum install -y wget unzip ImageMagick graphviz cairo-devel libXt-devel

WORKDIR /tmp
RUN yum install -y R

RUN R -e "install.packages('Nozzle.R1', repos='http://cran.rstudio.com/', INSTALL_opts='--no-html')"
RUN R -e "install.packages('ggplot2', repos='http://cran.rstudio.com/', INSTALL_opts='--no-html')"
RUN R -e "install.packages('parallel', repos='http://cran.rstudio.com/', INSTALL_opts='--no-html')"
RUN R -e "install.packages('reshape2', repos='http://cran.rstudio.com/', INSTALL_opts='--no-html')"
RUN R -e "install.packages('plyr', repos='http://cran.rstudio.com/', INSTALL_opts='--no-html')"
RUN R -e "install.packages('BBmisc', repos='http://cran.rstudio.com/', INSTALL_opts='--no-html')"
RUN R -e "source('http://bioconductor.org/biocLite.R'); biocLite('rgl')"
RUN R -e "install.packages('mixOmics', repos='http://cran.rstudio.com/', INSTALL_opts='--no-html')"
RUN R -e "install.packages('missForest', repos='http://cran.rstudio.com/', INSTALL_opts='--no-html')"
RUN R -e "install.packages('doParallel', repos='http://cran.rstudio.com/', INSTALL_opts='--no-html')"
RUN R -e "install.packages('DiscriMiner', repos='http://cran.rstudio.com/', INSTALL_opts='--no-html')"
RUN R -e "install.packages('ape', repos='http://cran.rstudio.com/', INSTALL_opts='--no-html')"
RUN R -e "install.packages('scatterplot3d', repos='http://cran.rstudio.com/', INSTALL_opts='--no-html')"
RUN R -e "install.packages('pheatmap', repos='http://cran.rstudio.com/', INSTALL_opts='--no-html')"
RUN R -e "install.packages('bootstrap', repos='http://cran.rstudio.com/', INSTALL_opts='--no-html')"
RUN R -e "install.packages('boot', repos='http://cran.rstudio.com/', INSTALL_opts='--no-html')"
RUN R -e "install.packages('caret', repos='http://cran.rstudio.com/', INSTALL_opts='--no-html')"
RUN R -e "install.packages('dplyr', repos='http://cran.rstudio.com/', INSTALL_opts='--no-html')"
RUN R -e "install.packages('stringr', repos='http://cran.rstudio.com/', INSTALL_opts='--no-html')"
RUN R -e "install.packages('RColorBrewer', repos='http://cran.rstudio.com/', INSTALL_opts='--no-html')"
RUN R -e "install.packages('RCurl', repos='http://cran.rstudio.com/', INSTALL_opts='--no-html')"
RUN R -e "install.packages('lattice', repos='http://cran.rstudio.com/', INSTALL_opts='--no-html')"
RUN R -e "install.packages('data.table', repos='http://cran.rstudio.com/', INSTALL_opts='--no-html')"
RUN R -e "install.packages('tidyr', repos='http://cran.rstudio.com/', INSTALL_opts='--no-html')"
RUN R -e "install.packages('scales', repos='http://cran.rstudio.com/', INSTALL_opts='--no-html')"
RUN R -e "install.packages('VennDiagram', repos='http://cran.rstudio.com/', INSTALL_opts='--no-html')"
RUN R -e "install.packages('pROC', repos='http://cran.rstudio.com/', INSTALL_opts='--no-html')"
RUN R -e "install.packages('readr', repos='http://cran.rstudio.com/', INSTALL_opts='--no-html')"
RUN R -e "install.packages('e1071', repos='http://cran.rstudio.com/', INSTALL_opts='--no-html')"
RUN R -e "install.packages('randomForest', repos='http://cran.rstudio.com/', INSTALL_opts='--no-html')"
RUN R -e "install.packages('coop', repos='http://cran.rstudio.com/', INSTALL_opts='--no-html')"
RUN R -e "install.packages('fpc', repos='http://cran.rstudio.com/', INSTALL_opts='--no-html')"
RUN R -e "install.packages('ROCR', repos='http://cran.rstudio.com/', INSTALL_opts='--no-html')"
RUN R -e "install.packages('devtools', repos='http://cran.rstudio.com/', INSTALL_opts='--no-html')"
RUN R -e "library(devtools); install_github('gaborcsardi/pkgconfig'); install_github('igraph/rigraph')"
RUN R -e "source('http://bioconductor.org/biocLite.R'); biocLite(c('pcaMethods', 'multtest'))"
RUN R -e "install.packages('DiffCorr', repos='http://cran.rstudio.com/', INSTALL_opts='--no-html')"
RUN R -e "source('http://bioconductor.org/biocLite.R'); biocLite('vsn')"
RUN R -e "source('http://bioconductor.org/biocLite.R'); biocLite('pls')"
RUN R -e "install.packages('XML', repos='http://cran.rstudio.com/', INSTALL_opts='--no-html')"
RUN R -e "install.packages('mzID', repos='http://cran.rstudio.com/', INSTALL_opts='--no-html')"
RUN R -e "source('http://bioconductor.org/biocLite.R'); biocLite('MSnbase')"
RUN R -e "source('http://bioconductor.org/biocLite.R'); biocLite('xcms')"
RUN R -e "source('http://bioconductor.org/biocLite.R'); biocLite('faahKO')"
RUN R -e "source('http://bioconductor.org/biocLite.R'); biocLite('mzR')"
RUN R -e "source('http://bioconductor.org/biocLite.R'); biocLite('preprocessCore')"
RUN R -e "source('http://bioconductor.org/biocLite.R'); biocLite('impute')"
RUN R -e "source('http://bioconductor.org/biocLite.R'); biocLite('CAMERA')"
RUN R -e "source('http://bioconductor.org/biocLite.R'); biocLite('ropls')"
RUN R -e "source('http://bioconductor.org/biocLite.R'); biocLite('sva')"
RUN R -e "source('http://bioconductor.org/biocLite.R'); biocLite('SSPA')"

# Download metaX
RUN wget http://metax.genomics.cn/files/metaX_1.4.16.tar.gz

# q()
RUN export R_MAX_NUM_DLLS=150

# Install metaX
RUN R -e "install.packages('/tmp/metaX_1.4.16.tar.gz', repos=NULL, type='source')"