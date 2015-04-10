#!/bin/bash

###########################################################
# Quick ncbi database fetch script
# Author: Tyghe Vallard
# Date: 7/18/2014
# Purpose:
#   Fetch/update ncbi databases
# Usage:
#  ./get_blast_dbs.sh /path/to/put/databases db1 db2 ...
#    Where db1, db2... is a space separated list of database names like
#    nt taxdb nr gene
###########################################################

function get_db_list {
    db=$1
    echo "Fetching file list from ftp site"
    _LIST=( $(curl -s ftp.ncbi.nih.gov/blast/db/ | grep -Eo "[[:space:]]${db}(.[0-9]+){0,1}.tar.gz(.md5){0,1}" | tr -d ' ') )
}

function download_file {
    file=$1
    echo "Downloading $file"
    wget ftp.ncbi.nih.gov/blast/db/${file}
}

function download_db {
    db=$1
    echo "Downloading $db"
    get_db_list $db
    files=( ${_LIST[@]} )

    i=0
    for file in ${files[@]}
    do
        download_file $file
    done
    cat *.md5 > md5sums && md5sum -c md5sums
    rm *.md5
    _FILES=${files[@]}
    return $?
}

function update_db {
    db=$1
    mkdir -p $db
    pushd $db

    download_db $db
    if [ $? -ne 0 ]
    then
        return $?
    fi

    for file in ${_FILES[@]}
    do
        if [[ $file =~ .tar.gz$ ]]
        then
            tar xzvf $file && rm $file
        fi
    done
    popd
}

function main {
    dst=$1
    shift

    mkdir -p $dst && cd $dst

    for db in $@
    do
        echo "Fetching/Updating $db from ncbi"
        if ! update_db $db
        then
            echo "Failed to update $db"
        fi
    done
}

main $@
