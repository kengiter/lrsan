#!/bin/bash

# binary choices
BASE="http://releases.llvm.org"
VERN="7.0.0"

# platform specifics
case "$OSTYPE" in
  linux*)
    DIST="x86_64-linux-gnu-ubuntu-16.04"
    ;;
  darwin*)
    DIST="x86_64-apple-darwin"
    ;;
  bsd*)
    DIST="amd64-unknown-freebsd10"
    ;;
  *)
    echo "Unknown OS ($OSTYPE)"; exit 
    ;;
esac

# download pre-built package 
PKGN=clang+llvm-$VERN-$DIST
TARF=$PKGN.tar.xz

LINK=$BASE/$VERN/$TARF
wget $LINK

# extract to location
tar xvf $TARF
mv $PKGN build

# remove file
rm $TARF
