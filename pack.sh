#!/bin/bash

source ./version.in
VER="${version_yy}.${version_month}"

SYS=`uname -m`
COD=`lsb_release -c | cut -f2`

if [ ! -d ./release ] ; then
  mkdir ./release
fi

NAME="xmv_${VER}_${COD}_${SYS}"
DIR=release/$NAME

mkdir ./$DIR
cp src/xmolview $DIR/$NAME

date > $DIR/README
uname -a >> $DIR/README
lsb_release -a >> $DIR/README


cd release
  zip $NAME".zip" ./$NAME/*
cd -
