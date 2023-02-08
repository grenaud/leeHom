#!/bin/bash

version=$(git describe --tags $(git rev-list --tags --max-count=1))

if [ $? -eq 0 ]; then
    echo "#define VERSION \"$version\"" > version.h
else
    echo "#define VERSION "NA"" > version.h
fi
