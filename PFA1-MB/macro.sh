#!/usr/bin/env bash

# Lista opcji do uruchomienia w makrze ROOT
OPTIONS=(
  "ElectronPEndCapOptions"
  "ElectronNEndCapOptions"
  "ElectronBarrelOptions"
  "PionLFOptions"
  "PionPEndCapOptions"
  "PionNEndCapOptions"
  "PionBarrelOptions"
)

echo "=================================================="
echo " Starting benchmark processing for all options... "
echo "=================================================="

# Tworzymy katalog na wykresy, jeśli nie istnieje
mkdir -p Plots

for OPT in "${OPTIONS[@]}"; do
    echo ""
    echo "--------------------------------------------------"
    echo " Running benchmark for: ${OPT}"
    echo "--------------------------------------------------"
    
    # Kompilacja i uruchomienie makra z wybraną opcją
    root -b -l -q "BenchmarkTrackCluster.cxx++(${OPT})"
    
    if [ $? -eq 0 ]; then
        echo "--> SUCCESS: ${OPT}"
    else
        echo "--> ERROR: Failed to run for ${OPT}"
    fi
done
echo ""
echo "=================================================="
echo " All benchmark tasks finished!                    "
echo "=================================================="