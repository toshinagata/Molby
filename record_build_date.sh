#!/bin/sh
#  Usage: sh ../record_build_date.sh [--with-svn-status]
echo "last_build = \"`date '+%Y-%m-%d %H:%M:%S %Z'`\"" > buildInfo.txt
if [ "$1" = "--with-svn-status" ]; then
  (cd ..; svn status -v . --depth=empty | awk '{print $1}' > revisionInfo.txt)
fi
