#!/bin/bash

sh sh_files/download_data.sh
sh sh_files/pre_process_data.sh
sh sh_files/run_benchmark.sh
sh sh_files/internal_main_benchmark.sh
sh sh_files/get_models_and_results.sh

