FROM nfcore/base:1.14
LABEL authors="chris.cheshire@crick.ac.uk"

# Install conda packages
COPY ./environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/compbio-variant-calling/bin:$PATH
