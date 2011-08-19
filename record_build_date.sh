#!/bin/sh
#  Usage: sh ../record_build_date.sh
echo "last_build = \"`date '+%Y-%m-%d %H:%M:%S %Z'`\"" > buildInfo.txt
