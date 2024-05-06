#!/bin/bash

# Initialize an array to store the results for each value of x
results1=()
results2=()
mins1=()
mins2=()
maxs1=()
maxs2=()

# Loop through values of x from 1 to 20
for ((x=1; x<=20; x++))
do
    # Initialize an array to store the results for this value of x
    x_results1=()
    x_results2=()

    # Perform the operation 10 times
    for ((i=1; i<=10; i++))
    do
        # Run the python script with x as input and get the output
        output=$(python3.10 phylogeny_generator.py $x)
	cp normal.ctl codeml.ctl 
        output=$(./bin/codeml)
        double=$(./testing_calc.sh)
	echo $double
        # Add the double to the array of results for this value of x
        x_results1+=($double)
        cp shadow.ctl codeml.ctl
        output=$(./bin/codeml)
        # Run the executable on the output to get a double
	double=$(./testing_calc.sh)
        echo $double
        # Add the double to the array of results for this value of x
        x_results2+=($double)
    done

    # Calculate the average for this value of x
    sum=0
    for num in "${x_results1[@]}"
    do
        sum=$(echo $sum + $num | bc)
    done
    average=$(echo "scale=2; $sum / ${#x_results1[@]}" | bc)
    echo "Hello"
    echo $average
    # Add the average to the array of overall results
    results1+=("$x: $average")

    min_val=${x_results1[0]}
    max_val=${x_results1[0]}
    for num in "${x_results1[@]}"
    do
        min_val=$(echo "if ($num < $min_val) $num else $min_val" | bc)
        max_val=$(echo "if ($num > $max_val) $num else $max_val" | bc)
    done    
    mins1+=("$x: $min_val")
    maxs1+=("$x: $max_val")


    sum=0
    for num in "${x_results2[@]}"
    do
        sum=$(echo $sum + $num | bc)
    done
    average=$(echo "scale=2; $sum / ${#x_results2[@]}" | bc)
    echo "Hello"
    echo $average
    # Add the average to the array of overall results
    results2+=("$x: $average")

    sum_of_squares=0
    for num in "${x_results2[@]}"
    do
        sum_of_squares=$(echo "$sum_of_squares + (($num - $average)^2)" | bc)
    done
    std_dev=$(echo "scale=2; sqrt($sum_of_squares / ${#x_results2[@]})" | bc -l)
    std_dev2+=("$x: $std_dev")

    min_val=${x_results2[0]}
    max_val=${x_results2[0]}
    for num in "${x_results2[@]}"
    do
        min_val=$(echo "if ($num < $min_val) $num else $min_val" | bc)
        max_val=$(echo "if ($num > $max_val) $num else $max_val" | bc)
    done    
    mins2+=("$x: $min_val")
    maxs2+=("$x: $max_val")

done

# Print the results for each value of x
printf '%s\n' "${results1[@]}" > r1.txt
printf '%s\n' "${results2[@]}" > r2.txt
printf '%s\n' "${mins1[@]}" > min1.txt
printf '%s\n' "${mins2[@]}" > min2.txt
printf '%s\n' "${maxs1[@]}" > max1.txt
printf '%s\n' "${maxs2[@]}" > max2.txt
for result in "${results1[@]}"
do
    echo $result
done
