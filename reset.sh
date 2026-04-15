#!/bin/bash

# Navigate to the pipeline directory
cd $HOME/Downloads/ont-wgs-pipeline || { echo "Failed to change directory"; exit 1; }

# Remove test_results directory
rm -rf test_results/
echo "Removed test_results/ directory"

# Loop to clean nextflow until no output
while true; do
    # Run nextflow clean command and capture output
    output=$(nextflow clean -f 2>&1)
    
    # Check if the output indicates nothing to clean
    if echo "$output" | grep -q "no pipeline was executed"; then
        # Nothing to clean - we're done
        echo "Nextflow clean completed - no pipeline history found"
        break
    elif [ -z "$output" ]; then
        # No output at all - also done
        echo "Nextflow clean completed - nothing to clean"
        break
    else
        # There was actual cleaning output
        echo "DEBUG: Nextflow cleaned something - continuing..."
        echo "Output: $output"
    fi
done

# Remove work directory contents (with sudo)
sudo rm -rf work/*
echo "Removed contents of work/ directory"

echo ""
read -p "Do you want to run nextflow again? (y/n): " answer

if [[ "$answer" == "y" || "$answer" == "Y" ]]; then
    echo "Running nextflow main..."
    nextflow run main.nf -profile test,docker -resume --region_bed_file data/region_hg002_90.bed
else
    echo "Skipping nextflow run. Script completed."
fi
