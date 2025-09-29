FROM jupyter/datascience-notebook

USER root
RUN apt update

RUN pip install anndata scanpy umap-learn leidenalg igraph 
RUN pip install scrublet harmonypy
RUN pip install seaborn --upgrade
RUN pip install -U pymde

#RUN pip install -U scvi-tools
#RUN pip install cellbender

USER foster
WORKDIR /workspace
