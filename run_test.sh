#!/bin/bash
ver=`echo $3 | awk -F- '{print $2}'`
num=$2
bench=$1
shift 2
for j in $(seq 1 1 $num); do
  for i in $(cat ../batch/$bench); do
    echo "Running $@ --bench $i $ver >> ../results/$bench-$ver"
    $@ --bench $i $ver >> ../results/$bench-$ver
  done
done
