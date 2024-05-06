#!/bin/bash

# Path to the text file
file_path="mlc"

# Use awk to process the file
awk '
BEGIN {process=0}  # Initialize a flag to false
/dN\/dS \(w\) for site classes \(K=3\)/ {process=1; next}  # Set the flag when header is found
process==1 && /^p:/ {p1=$2; p2=$3; p3=$4; next}  # Process p line if flag is true
process==1 && /^w:/ {  # Process w line if flag is true
    w1=$2; w2=$3; w3=$4;
    # Calculate the result and print it
    result = (p1 * w1) + (p2 * w2) + (p3 * w3);
    print result;  # Print just the result
    process=0;  # Reset the flag
}' "$file_path"
