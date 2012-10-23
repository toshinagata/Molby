#!/bin/sh
#  Usage: sh ../record_build_date.sh [--with-svn-status]
LAST_BUILD=`date '+%Y-%m-%d %H:%M:%S %Z'`
echo "last_build = \"$LAST_BUILD\"" > buildInfo.txt
echo "char *gLastBuildString = \"$LAST_BUILD\";" > build/buildInfo.c
if [ "$1" = "--with-svn-status" ]; then
  REVISION_INFO=`(cd ..; svn status -v . --depth=empty | awk '{print $1}')`
  echo $REVISION_INFO > ../revisionInfo.txt
else
  REVISION_INFO=`cat ../revisionInfo.txt 2>/dev/null`
  if [ "$REVISION_INFO" = "" ]; then
    REVISION_INFO=0
  fi
fi
echo "int gRevisionNumber = $REVISION_INFO;" >> build/buildInfo.c
