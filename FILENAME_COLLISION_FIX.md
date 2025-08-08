# Filename Collision Bug Fix

## Problem Description

The PHMM-Tree program had a critical bug in its name shortening algorithm that caused collisions when creating temporary filenames. This occurred when similar initial names were processed, resulting in identical shortened names that would overwrite each other.

## Root Cause

The bug was located in `matrix_deal.cpp` in the name shortening logic (lines ~500-540). The original algorithm:

1. Shortened long names using: `first_3_chars + last_7_chars`
2. Attempted collision resolution by replacing the end of the string with numbers
3. **Bug**: The collision resolution used `str.replace()` incorrectly, often creating new collisions or invalid names

### Original Problematic Code
```cpp
// Original buggy collision resolution
while(shorted_names_map.find(str_hmms_names) != shorted_names_map.end() ){
    num++;
    std::string str_shorted_num = int_2_string(num+1);
    str_hmms_names=str_hmms_names.replace(str_hmms_names.length() - str_shorted_num.length(), str_shorted_num.length(),str_shorted_num.c_str());
}
```

## Solution

Implemented a robust collision resolution algorithm that:

1. Creates a base shortened name: `first_3_chars + last_7_chars`
2. Detects collisions using the existing hashmap
3. Resolves collisions by appending unique suffixes (`_1`, `_2`, etc.)
4. Handles length constraints by truncating the base name when necessary
5. Falls back to numeric-only suffixes for very long collision counters

### New Fixed Code
```cpp
// Create base shortened name: first 3 + last 7 characters
std::string base_short_name = str_hmms_names.substr(0,3) + str_hmms_names.substr(str_hmms_names.length()-7,7);
str_hmms_names = base_short_name;

// Handle collisions by appending a unique suffix
int collision_counter = 1;
while(shorted_names_map.find(str_hmms_names) != shorted_names_map.end() ){
    std::string suffix = "_" + int_2_string(collision_counter);
    
    // If the suffix would make the name too long, truncate the base name
    if(base_short_name.length() + suffix.length() > 10) {
        int available_length = 10 - suffix.length();
        if(available_length < 3) {
            // If suffix too long, use shorter suffix format
            suffix = int_2_string(collision_counter);
            available_length = 10 - suffix.length();
        }
        str_hmms_names = base_short_name.substr(0, available_length) + suffix;
    } else {
        str_hmms_names = base_short_name + suffix;
    }
    
    collision_counter++;
    std::cout<<"Collision resolved: " << str_hmms_names << std::endl;
}
```

## Fix Details

### Files Modified
- `matrix_deal.cpp` (lines ~500-540): Fixed both NAME and ACC collision resolution

### Key Improvements
1. **Proper Suffix Handling**: Appends `_1`, `_2`, etc. instead of replacing characters
2. **Length Management**: Dynamically truncates base name to fit within 10-character limit
3. **Guaranteed Uniqueness**: Incremental counter ensures no duplicate names
4. **Fallback Strategy**: Uses numeric-only suffixes when underscore format is too long
5. **Debug Output**: Added logging for collision resolution process

## Testing Results

Tested with collision-prone names:
```
ABC_test_1234567 -> ABC1234567   (original)
ABC_test_1234567 -> ABC12345_1   (collision resolved)
ABC_diff_1234567 -> ABC12345_2   (collision resolved) 
DEF_other_7654321 -> DEF7654321  (no collision)
ABC_another_4567 -> ABCer_4567   (no collision)
```

✓ All collisions successfully resolved  
✓ All names stay within 10-character limit  
✓ All names are unique  

## Impact

This fix prevents:
- Temporary file overwrites that caused data loss
- Program crashes due to file conflicts  
- Incorrect phylogenetic analysis results
- Workflow interruptions in high-throughput scenarios

The fix maintains backward compatibility while ensuring robust collision-free filename generation.
