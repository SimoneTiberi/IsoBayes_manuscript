#!/bin/bash

sh sh_files/download_data.sh
sh sh_files/pre_process_data.sh

##################################################################################################
# Once data have been downloaded and preprocessed, to run the entire benchmark analysis
# you must move the IsoBayes_paper folder on a HPC infrastructure endowed with OpenPBS software
# and then execute the next two commented lines:

#sh sh_files/run_benchmark.sh
#sh sh_files/internal_main_benchmark.sh

##################################################################################################

sh sh_files/get_models_and_results.sh

