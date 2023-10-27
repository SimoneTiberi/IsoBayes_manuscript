#!/bin/bash

echo "Downloading Data.zip from https://figshare.com/ndownloader/files/42929275"
wget -O ./Data.zip -c https://figshare.com/ndownloader/files/42929275

echo "Unzipping Data.zip"
unzip Data.zip

