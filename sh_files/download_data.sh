#!/bin/bash

echo "Downloading Containers.zip"
mkdir tmp
wget -O ./tmp/Containers.zip -c https://figshare.com/ndownloader/files/43243176
cd tmp
unzip Containers.zip
rsync -avhz Containers/* ../Containers/
cd ..
rm -r tmp

echo "Downloading Data.zip from https://figshare.com/ndownloader/files/42929275"
wget -O ./Data.zip -c https://figshare.com/ndownloader/files/42929275

echo "Unzipping Data.zip"
unzip Data.zip

