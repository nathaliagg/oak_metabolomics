#!/bin/bash

# $1 == output dir name
# $2 == input mgf file

$SIRIUS --output $1 \
        --input $2 \
        formula --database='ALL' --tree-timeout=300 --compound-timeout=500 \
        --candidates=5 --elements-considered='SBrCl' --elements-enforced='CHON[3]P[1]' \
        structure \
        canopus
