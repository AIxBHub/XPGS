#!/bin/bash

input_vcf=$1
output_vcf=$2

awk 'BEGIN { OFS="\t" }
{
    if(/^##/ || (NR==1 && $0 ~ /^#/)) {
    print $0;
    } else {
    print $1, $2, $3, $4, $5, $6, $7, $8;
    }
}' "$input_vcf" > "$output_vcf"

