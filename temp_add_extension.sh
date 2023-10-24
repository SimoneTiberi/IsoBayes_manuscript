#!/bin/bash


# Specify the source directory and name pattern

src_dir="./"

pattern="Merged_validation_res*"


# Specify the destination directory and the new extension

dest_dir="./"

new_extension=".RData"

find "$src_dir" -type f -name "$pattern" | while read -r file

do

    # Create new file name
	mv $file $file.RData
done
