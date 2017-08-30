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
	  yum install libcurl-devel -y
	  yum install netcdf-devel -y
	  yum install libxml2-devel -y
	  # For rgl R package
		yum install libpng-devel libX11-devel mesa-libGL mesa-libGL-devel mesa-libGLU mesa-libGLU-devel -y
	  EOH
end

bash 'Install metaX dependencies in R' do
  user 'root'
  cwd '/tmp'
  code <<-EOH
	  R -e \"install.packages('Nozzle.R1', repos='http://cran.rstudio.com/')\"
	  R -e \"install.packages('ggplot2', repos='http://cran.rstudio.com/')\"
	  R -e \"install.packages('parallel', repos='http://cran.rstudio.com/')\"
	  R -e \"install.packages('reshape2', repos='http://cran.rstudio.com/')\"
	  R -e \"install.packages('plyr', repos='http://cran.rstudio.com/')\"
	  R -e \"install.packages('BBmisc', repos='http://cran.rstudio.com/')\"
		R -e \"source('http://bioconductor.org/biocLite.R'); biocLite('rgl')\"
	  R -e \"install.packages('mixOmics', repos='http://cran.rstudio.com/')\"
	  R -e \"install.packages('missForest', repos='http://cran.rstudio.com/')\"
	  R -e \"install.packages('doParallel', repos='http://cran.rstudio.com/')\"
	  R -e \"install.packages('DiscriMiner', repos='http://cran.rstudio.com/')\"
	  R -e \"install.packages('ape', repos='http://cran.rstudio.com/')\"
	  R -e \"install.packages('scatterplot3d', repos='http://cran.rstudio.com/')\"
	  R -e \"install.packages('pheatmap', repos='http://cran.rstudio.com/')\"
	  R -e \"install.packages('bootstrap', repos='http://cran.rstudio.com/')\"
		R -e \"install.packages('boot', repos='http://cran.rstudio.com/')\"
    R -e \"install.packages('caret', repos='http://cran.rstudio.com/')\"
    R -e \"install.packages('dplyr', repos='http://cran.rstudio.com/')\"
    R -e \"install.packages('stringr', repos='http://cran.rstudio.com/')\"
 	  R -e \"install.packages('RColorBrewer', repos='http://cran.rstudio.com/')\"
    R -e \"install.packages('RCurl', repos='http://cran.rstudio.com/')\"
	  R -e \"install.packages('lattice', repos='http://cran.rstudio.com/')\"
	  R -e \"install.packages('data.table', repos='http://cran.rstudio.com/')\"
    R -e \"install.packages('tidyr', repos='http://cran.rstudio.com/')\"
    R -e \"install.packages('scales', repos='http://cran.rstudio.com/')\"
    R -e \"install.packages('VennDiagram', repos='http://cran.rstudio.com/')\"
	  R -e \"install.packages('pROC', repos='http://cran.rstudio.com/')\"
		R -e \"install.packages('readr', repos='http://cran.rstudio.com/')\"
    R -e \"install.packages('e1071', repos='http://cran.rstudio.com/')\"
    R -e \"install.packages('randomForest', repos='http://cran.rstudio.com/')\"
		R -e \"install.packages('coop', repos='http://cran.rstudio.com/')\"
    R -e \"install.packages('fpc', repos='http://cran.rstudio.com/')\"
    R -e \"install.packages('ROCR', repos='http://cran.rstudio.com/')\"
    R -e \"install.packages('devtools', repos='http://cran.rstudio.com/')\"
	  R -e \"library(devtools); install_github('gaborcsardi/pkgconfig'); install_github('igraph/rigraph')\"
    R -e \"source('http://bioconductor.org/biocLite.R'); biocLite(c('pcaMethods', 'multtest'))\"
    R -e \"install.packages('DiffCorr', repos='http://cran.rstudio.com/')\"
	  R -e \"source('http://bioconductor.org/biocLite.R'); biocLite('vsn')\"
	  R -e \"source('http://bioconductor.org/biocLite.R'); biocLite('pls')\"
	  R -e \"install.packages('XML', repos='http://cran.rstudio.com/')\"
		R -e \"install.packages('mzID', repos='http://cran.rstudio.com/')\"
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

# For stopping the maximal number of DLLs reached error message
ENV['R_MAX_NUM_DLLS'] = '150'

bash 'Install metaX' do
  user 'root'
  code <<-EOH
		R -e \"install.packages('/tmp/metaX_1.4.16.tar.gz', repos=NULL, type='source')\"
  EOH
end