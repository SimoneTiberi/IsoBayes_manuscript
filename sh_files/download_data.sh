#!/bin/bash

echo "Downloading Data.zip from https://figshare.com/ndownloader/files/42894736"
wget -O ./Data.zip -c https://figshare.com/ndownloader/files/42894736

echo "Unzipping Data.zip"
unzip Data.zip

