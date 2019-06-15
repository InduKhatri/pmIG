#!/bin/bash

for f in ~/vcf/*.vcf

do
	python VCF2Ped.py $f
done