#!/bin/bash

#goal is to use synapse client to upload files.

#first get parentId
parentid='syn5061128'

#then upload pdfs
for file in `ls abs_res/*pdf`
do
    synapse store $file --parentId $parentid
done


#then upload rdata files
for file in `ls abs_res/*RData`
do
    synapse store $file --parentId $parentid
done
