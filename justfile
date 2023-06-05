# works for debian releases
# Run just to deply the env.

all: dev blast
    echo "Check the msg to if the dev env has been deployed."

dev:
    sudo apt update && sudo apt install libxml2 pandoc libglpk-dev -y
    # https://github.com/r-lib/rig
    curl -Ls https://github.com/r-lib/rig/releases/download/latest/rig-linux-latest.tar.gz | sudo tar xz -C /usr/local
    rig add release

    pip3 install -U radian
    R -e 'pak::pkg_install(c("devtools", "languageserver", "pkgdown"))'
    R -e 'if (file.exists("DESCRIPTION")) remotes::install_local(dependencies = TRUE, force = TRUE)'

blast:
    cd /workspaces &&\
    wget -c https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.14.0+-x64-linux.tar.gz &&\
    tar zxvf ncbi-blast-2.14.0+-x64-linux.tar.gz &&\
    mv ncbi-blast-2.14.0+ blast

[no-cd]
test:
    -mkdir /workspaces/phylotaR-test
    Rscript test.R