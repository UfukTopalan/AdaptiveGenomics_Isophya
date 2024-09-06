#!/bin/bash

# Read the list of population names
pop_list="pop11_list"

# Initialize a variable to hold the intersection result
intersection_result="nonparalog_sites"

# Read each population name from the list
while IFS= read -r pop_name; do
    # Construct the file name from the population name
    sites_file="${pop_name}_sites"
    
    # If this is the first file, just copy it to the result
    if [ ! -f "$intersection_result" ]; then
        cp "$sites_file" "$intersection_result"
    else
        # Use awk to extract common sites across populations
        awk 'NR==FNR{a[$0];next} $0 in a' "$sites_file" "$intersection_result" > temp_sites
        mv temp_sites "$intersection_result"
    fi
done < "$pop_list"

echo "Sites present across all populations are saved in $intersection_result"
