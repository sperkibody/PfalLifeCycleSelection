#!/bin/sh
#use .bcftools-1.9
use Tabix 
file=$1
bgzip $file
