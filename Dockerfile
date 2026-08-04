FROM rocker/r-ver:4.3

# ── System dependencies ───────────────────────────────────────────────────
# curl/https support for R package installers + git for devtools/remotes
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libgit2-dev \
    && rm -rf /var/lib/apt/lists/*

# bzip2/wget needed by reticulate::install_miniconda()
RUN apt-get update && apt-get install -y \
    wget \
    bzip2 \
    && rm -rf /var/lib/apt/lists/*

# ── R package installs (kept small/separate so a failure is cheap to retry) ─
RUN Rscript -e "install.packages('BiocManager', repos='https://cloud.r-project.org')"

RUN Rscript -e "BiocManager::install('org.Hs.eg.db', update = FALSE, ask = FALSE)"

# fontconfig/freetype2 etc. needed to compile ragg/textshaping/systemfonts
# (transitive deps of devtools -> pkgdown), kept as its own layer for caching
RUN apt-get update && apt-get install -y \
    libfontconfig1-dev \
    libfreetype6-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libwebp-dev \
    && rm -rf /var/lib/apt/lists/*

RUN apt-get update && apt-get install -y libuv1-dev && rm -rf /var/lib/apt/lists/*

RUN Rscript -e "install.packages('data.table', repos='https://cloud.r-project.org')"
RUN Rscript -e "install.packages('devtools', repos='https://cloud.r-project.org')"
RUN Rscript -e "install.packages('remotes', repos='https://cloud.r-project.org')"

RUN Rscript -e "install.packages('reticulate', repos='https://cloud.r-project.org')"

# ── Python / conda env used by reticulate::use_condaenv('r-tf23') ──────────
RUN Rscript -e "reticulate::install_miniconda()"

RUN Rscript -e "reticulate::conda_create('r-tf23', python_version = '3.7')"

# ── R wrappers around the python env + DeepCC itself ────────────────────────
RUN Rscript -e "install.packages('tensorflow', repos='https://cloud.r-project.org')"

RUN Rscript -e "install.packages('keras', repos='https://cloud.r-project.org')"

RUN Rscript -e "remotes::install_github('CityUHK-CompBio/DeepCC')"

# h5py needed by keras::load_model_hdf5() to read .hdf5 model files
RUN Rscript -e "reticulate::conda_install('r-tf23', packages = 'h5py')"

# tensorflow==2.3.0 isn't published on conda-forge (only PyPI / the default
# Anaconda channel, which isn't enabled here), so install it via pip instead.
RUN Rscript -e "reticulate::conda_install('r-tf23', packages = 'tensorflow==2.3.0', pip = TRUE)"

WORKDIR /app/214-NANOSTRING-18-ANALYSIS

CMD ["R"]
