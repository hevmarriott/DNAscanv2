#!/bin/bash
#Usage: bash add_read_group.sh $RG_ID $RG_LB $RG_PL $RG_PU $RG_SM
RG_ID=$1
RG_LB=$2
RG_PL=$3
RG_PU=$4
RG_SM=$5
sed "s|RG_ID = \"\"|RG_ID = \"$RG_ID\"|" scripts/paths_configs.py > scripts/paths_configs.py_temp
sed "s|RG_LB = \"\"|RG_LB = \"$RG_LB\"|" scripts/paths_configs.py_temp > scripts/paths_configs.py
sed "s|RG_SM = \"\"|RG_SM = \"$RG_SM\"|" scripts/paths_configs.py > scripts/paths_configs.py_temp
sed "s|RG_PU = \"\"|RG_PU = \"$RG_PU\"|" scripts/paths_configs.py_temp > scripts/paths_configs.py
sed "s|RG_PL = \"\"|RG_PL = \"$RG_PL\"|" scripts/paths_configs.py > scripts/paths_configs.py_temp
mv scripts/paths_configs.py_temp scripts/paths_configs.py
chmod +x scripts/*
