#!/bin/bash
SAMPLE="$( cat $1 | cut -d '|' -f1 )"
URL="$( cat $1 | cut -d'|' -f2)"
EXTENSION="$(basename $URL | cut -d'.' -f2,3)"

wget -c -O "$SAMPLE"."$EXTENSION" "$URL"
