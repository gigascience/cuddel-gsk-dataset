#
# Cookbook Name:: cuddel
# Recipe:: default
#
# Copyright 2016, GigaScience
#
# All rights reserved - Do Not Redistribute
#

#####################
### Install metaX ###
#####################

bash 'Install R' do
  user 'root'
  cwd '/tmp'
  code <<-EOH
    # yum groupinstall "Development Tools"
    yum install epel-release -y
    yum install R -y
    # Required by RCurl R package
    yum install libcurl-devel netcdf-devel libxml2-devel -y
    # For rgl R package
    yum install libpng-devel libX11-devel mesa-libGL mesa-libGL-devel mesa-libGLU mesa-libGLU-devel -y
  EOH
end

bash 'Install metaX dependencies in R' do
  user 'root'
  cwd '/tmp'
  code <<-EOH
    R -e \"install.packages('Nozzle.R1', dependencies=TRUE, repos='http://cran.rstudio.com/')\"
    R -e \"install.packages('ggplot2', dependencies=TRUE, repos='http://cran.rstudio.com/')\"
    # Warning messages by install parallel package:
    # 1: package ‘parallel’ is not available (for R version 3.5.0)
    # 2: package ‘parallel’ is a base package, and should not be updated
    R -e \"install.packages('parallel', dependencies=TRUE, repos='http://cran.rstudio.com/')\"
    R -e \"install.packages('reshape2', dependencies=TRUE, repos='http://cran.rstudio.com/')\"
    R -e \"install.packages('plyr', dependencies=TRUE, repos='http://cran.rstudio.com/')\"
    R -e \"install.packages('BBmisc', dependencies=TRUE, repos='http://cran.rstudio.com/')\"
    R -e \"source('http://bioconductor.org/biocLite.R'); biocLite('rgl')\"
    R -e \"install.packages('mixOmics', dependencies=TRUE, repos='http://cran.rstudio.com/')\"
    R -e \"install.packages('missForest', dependencies=TRUE, repos='http://cran.rstudio.com/')\"
    R -e \"install.packages('doParallel', dependencies=TRUE, repos='http://cran.rstudio.com/')\"
    R -e \"install.packages('DiscriMiner', dependencies=TRUE, repos='http://cran.rstudio.com/')\"
    R -e \"install.packages('ape', dependencies=TRUE, repos='http://cran.rstudio.com/')\"
    R -e \"install.packages('scatterplot3d', dependencies=TRUE, repos='http://cran.rstudio.com/')\"
    R -e \"install.packages('pheatmap', dependencies=TRUE, repos='http://cran.rstudio.com/')\"
    R -e \"install.packages('bootstrap', dependencies=TRUE, repos='http://cran.rstudio.com/')\"
    R -e \"install.packages('boot', dependencies=TRUE, repos='http://cran.rstudio.com/')\"
    R -e \"install.packages('caret', dependencies=TRUE, repos='http://cran.rstudio.com/')\"
    R -e \"install.packages('dplyr', dependencies=TRUE, repos='http://cran.rstudio.com/')\"
    R -e \"install.packages('stringr', dependencies=TRUE, repos='http://cran.rstudio.com/')\"
    R -e \"install.packages('RColorBrewer', dependencies=TRUE, repos='http://cran.rstudio.com/')\"
    R -e \"install.packages('RCurl', dependencies=TRUE, repos='http://cran.rstudio.com/')\"
    R -e \"install.packages('lattice', dependencies=TRUE, repos='http://cran.rstudio.com/')\"
    R -e \"install.packages('data.table', dependencies=TRUE, repos='http://cran.rstudio.com/')\"
    R -e \"install.packages('tidyr', dependencies=TRUE, repos='http://cran.rstudio.com/')\"
    R -e \"install.packages('scales', dependencies=TRUE, repos='http://cran.rstudio.com/')\"
    R -e \"install.packages('VennDiagram', dependencies=TRUE, repos='http://cran.rstudio.com/')\"
    R -e \"install.packages('pROC', dependencies=TRUE, repos='http://cran.rstudio.com/')\"
    R -e \"install.packages('readr', dependencies=TRUE, repos='http://cran.rstudio.com/')\"
    R -e \"install.packages('e1071', dependencies=TRUE, repos='http://cran.rstudio.com/')\"
    R -e \"install.packages('randomForest', dependencies=TRUE, repos='http://cran.rstudio.com/')\"
    R -e \"install.packages('coop', dependencies=TRUE, repos='http://cran.rstudio.com/')\"
    R -e \"install.packages('fpc', dependencies=TRUE, repos='http://cran.rstudio.com/')\"
    R -e \"install.packages('ROCR', dependencies=TRUE, repos='http://cran.rstudio.com/')\"
    R -e \"install.packages('devtools', dependencies=TRUE, repos='http://cran.rstudio.com/')\"
    R -e \"library(devtools); install_github('gaborcsardi/pkgconfig'); install_github('igraph/rigraph')\"
    R -e \"source('http://bioconductor.org/biocLite.R'); biocLite(c('pcaMethods', 'multtest'))\"
    R -e \"install.packages('DiffCorr', dependencies=TRUE, repos='http://cran.rstudio.com/')\"
    R -e \"source('http://bioconductor.org/biocLite.R'); biocLite('vsn')\"
    R -e \"source('http://bioconductor.org/biocLite.R'); biocLite('pls')\"
    R -e \"install.packages('XML', dependencies=TRUE, repos='http://cran.rstudio.com/')\"
    R -e \"source('http://bioconductor.org/biocLite.R'); biocLite('mzID')\"
    R -e \"source('http://bioconductor.org/biocLite.R'); biocLite('MSnbase')\"
    R -e \"source('http://bioconductor.org/biocLite.R'); biocLite('xcms')\"
    R -e \"source('http://bioconductor.org/biocLite.R'); biocLite('faahKO')\"
    R -e \"source('http://bioconductor.org/biocLite.R'); biocLite('mzR')\"
    R -e \"source('http://bioconductor.org/biocLite.R'); biocLite('preprocessCore')\"
    R -e \"source('http://bioconductor.org/biocLite.R'); biocLite('impute')\"
    R -e \"source('http://bioconductor.org/biocLite.R'); biocLite('CAMERA')\"
    R -e \"source('http://bioconductor.org/biocLite.R'); biocLite('ropls')\"
    R -e \"source('http://bioconductor.org/biocLite.R'); biocLite('sva')\"
    R -e \"source('http://bioconductor.org/biocLite.R'); biocLite('SSPA')\"
  EOH
end

# Download metaX
remote_file "/tmp/metaX_1.4.16.tar.gz" do
  source "http://metax.genomics.cn/files/metaX_1.4.16.tar.gz"
end

# Fix maximal number of DLLs reached error message
ENV['R_MAX_NUM_DLLS'] = '150'

bash 'Install metaX' do
  user 'root'
  code <<-EOH
    R -e \"install.packages('/tmp/metaX_1.4.16.tar.gz', repos=NULL, type='source')\"
  EOH
end