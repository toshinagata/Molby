#!/bin/sh
#  Usage: sh ../record_build_date.sh [--with-git-status] > buildInfo.c
LAST_BUILD=`date '+%Y-%m-%d %H:%M:%S %Z'`
echo "last_build = \"$LAST_BUILD\"" > buildInfo.txt
echo "volatile char *gLastBuildString = \"$LAST_BUILD\";"
REVISION_INFO=""
if [ "$1" = "--with-git-status" ]; then
  REVISION_INFO=`git rev-list HEAD --count`
#  REVISION_INFO=`(cd ..; svn status -v . --depth=empty | awk '{print $1}')`
  if [ "$REVISION_INFO" != "" ]; then
    echo $REVISION_INFO > ../revisionInfo.txt
  fi
fi
if [ "$REVISION_INFO" = "" ]; then
  REVISION_INFO=`cat ../revisionInfo.txt 2>/dev/null`
  if [ "$REVISION_INFO" = "" ]; then
    REVISION_INFO=0
  fi
fi
echo "volatile int gRevisionNumber = $REVISION_INFO;"
