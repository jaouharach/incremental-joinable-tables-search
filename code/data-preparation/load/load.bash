TAR_DIR=$1       # directory where tar files will be stored
OUT_DIR=$2       # directory where tar files will be extracted

for f in "$TAR_DIR/*.tar.gz"; do tar -xvf $f -C $OUT_DIR ; done # extract all tar files
