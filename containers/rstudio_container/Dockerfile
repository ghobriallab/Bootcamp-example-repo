FROM bioconductor/tidyverse:3.18

RUN Rscript -e "install.packages('Cairo', repos='http://cran.us.r-project.org')"

RUN Rscript -e 'BiocManager::install(c("pheatmap", "viridis", \
                       "scran", "scuttle", "scater", \
                       "bluster", "scDblFinder", "scds", "Rtsne",\ 
                       "UCell","DropletUtils","fgsea", \
                       "ggridges", "cowplot", "gridGraphics", \
                       "scales", "batchelor","zellkonverter"))'
                       
RUN Rscript -e "install.packages(c('SeuratObject','Seurat','harmony','hdf5r', \
                      'ggpubr','RColorBrewer','vegan','ComplexHeatmap','circlize', \
                      'survival','survminer','msigdbr'), \
                      repos='http://cran.us.r-project.org')"
                      
RUN Rscript -e "install.packages(c('hdf5r'), repos='http://cran.us.r-project.org')"
