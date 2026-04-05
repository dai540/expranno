FROM rocker/r-ver:4.5.1

RUN apt-get update && apt-get install -y --no-install-recommends \
    git \
    libcurl4-openssl-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    libfribidi-dev \
    libgit2-dev \
    libharfbuzz-dev \
    libjpeg-dev \
    libpng-dev \
    libssl-dev \
    libtiff5-dev \
    libxml2-dev \
    qpdf \
    zlib1g-dev \
  && rm -rf /var/lib/apt/lists/*

RUN R -q -e "install.packages(c('pak', 'BiocManager', 'remotes'), repos = 'https://cloud.r-project.org')"
RUN R -q -e "pak::pak('dai540/expranno')"
RUN R -q -e "source(system.file('scripts', 'install_optional_backends.R', package = 'expranno')); expranno_install_optional_backends()"

WORKDIR /work
