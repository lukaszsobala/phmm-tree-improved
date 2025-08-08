#!/bin/bash

# Simple test to demonstrate the improved time formatting
echo "=== Time Formatting Test ==="

# Test 1: Quick execution (should show seconds with decimal places)
echo "Test 1: Quick execution (< 1 minute)"
echo "Expected: X.XXX seconds"
echo

# Test 2: Medium execution (simulated - should show minutes and seconds)
echo "Test 2: Medium execution (1-60 minutes)"
echo "Expected: X minutes, Y seconds"
echo

# Test 3: Long execution (simulated - should show hours, minutes, seconds)  
echo "Test 3: Long execution (> 1 hour)"
echo "Expected: X hours, Y minutes, Z seconds"
echo

# Test 4: Very long execution (simulated - should show days, hours, minutes, seconds)
echo "Test 4: Very long execution (> 1 day)"
echo "Expected: X days, Y hours, Z minutes, W seconds"
echo

echo "The phmm-tree program will now show runtime in human-readable format:"
echo "- Short runs (< 60s): 'X.XXX seconds'"
echo "- Medium runs: 'X minutes, Y seconds'"  
echo "- Long runs: 'X hours, Y minutes, Z seconds'"
echo "- Very long runs: 'X days, Y hours, Z minutes, W seconds'"
echo

echo "Previous format was always 'XXXX ms' which was not useful for long analyses."
echo "New format automatically adapts to the duration for better readability."
