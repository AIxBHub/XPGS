#!/bin/bash

#1 loop through ontology ids
#2 use pgs catalog rest api to search for all pgs records available for each ontology
#3 filter pgs by ancestry and by r2 performance metric

ids=(OBA_VT0000217)
#    EFO_0004305
#    EFO_0004509
#    EFO_0004348
#    OBA_0003460
#    EFO_0004527
#    EFO_0004528
#    EFO_0009188
#    EFO_0004309
#    OBA_0003277
#    EFO_0007985
#    EFO_0007984)

# get records for each ontology
for id in ${ids[@]}; do
    trait_json=$(curl -X GET "https://www.pgscatalog.org/rest/score/search?trait_id=$id" -H  "accept: application/json")
    python parseJson.py $trait_json
done

