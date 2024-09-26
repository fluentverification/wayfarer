#!/bin/sh

# This test script is really hacky. Please don't judge me.

rm -rf results/isr
rm -rf results/isr_piped
# rm -rf results/single_order
# rm -rf results/single_order_piped

mkdir -p results/isr
mkdir -p results/isr_piped
# mkdir -p results/single_order
# mkdir -p results/single_order_piped

# Test ISR
python3 test.py
mv results/*.csv results/isr

# Test ISR with Piped space scaling
python3 test.py --piped
mv results/*.csv results/isr_piped

echo "WARNING! Not testing SOP"

# Test SOP
# python3 test.py --sop
# mv results/*.csv results/single_order
#
# # Test SOP with Piped space scaling
# python3 test.py --sop --piped
# mv results/*.csv results/single_order_piped
