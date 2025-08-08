#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

// Simple version of int_2_string for testing
std::string int_2_string(int num) {
    return std::to_string(num);
}

void test_collision_resolution() {
    // Test the collision resolution algorithm
    std::unordered_map<std::string, int> shorted_names_map;
    std::vector<std::string> vec_unshorted_names;
    std::vector<std::string> vec_shorted_names;
    
    // Test names that would create collisions with the old algorithm
    std::vector<std::string> test_names = {
        "ABC_test_1234567",  // Would shorten to: ABC1234567
        "ABC_test_1234567",  // Same - should cause collision
        "ABC_diff_1234567",  // Would shorten to: ABC1234567 (collision!)
        "DEF_other_7654321", // Would shorten to: DEF7654321
        "ABC_another_4567"   // Would shorten to: ABC_another_4567 (collision with first!)
    };
    
    std::cout << "Testing collision resolution fix:\n";
    std::cout << "================================\n\n";
    
    int shorted_num = 0;
    int all_names_num = 0;
    
    for(const auto& original_name : test_names) {
        std::cout << "Processing: " << original_name << std::endl;
        
        std::string str_hmms_names = original_name;
        
        if(str_hmms_names.length() > 10){
            vec_unshorted_names.push_back(str_hmms_names);
            std::cout << "  Original name: " << str_hmms_names << std::endl;
            
            // Create base shortened name: first 3 + last 7 characters
            std::string base_short_name = str_hmms_names.substr(0,3) + str_hmms_names.substr(str_hmms_names.length()-7,7);
            str_hmms_names = base_short_name;
            
            std::cout << "  Base shortened: " << base_short_name << std::endl;
            
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
                std::cout << "  Collision detected, trying: " << str_hmms_names << std::endl;
            }
            
            vec_shorted_names.push_back(str_hmms_names);
            shorted_num++;
        }
        
        shorted_names_map[str_hmms_names] = all_names_num;
        all_names_num++;
        
        std::cout << "  Final shortened: " << str_hmms_names << " (length: " << str_hmms_names.length() << ")" << std::endl;
        std::cout << std::endl;
    }
    
    std::cout << "Results Summary:\n";
    std::cout << "================\n";
    for(size_t i = 0; i < vec_unshorted_names.size(); i++) {
        std::cout << vec_unshorted_names[i] << " -> " << vec_shorted_names[i] << std::endl;
    }
    
    // Verify no collisions
    std::unordered_map<std::string, int> collision_check;
    bool has_collisions = false;
    for(const auto& short_name : vec_shorted_names) {
        if(collision_check.find(short_name) != collision_check.end()) {
            std::cout << "ERROR: Collision detected for: " << short_name << std::endl;
            has_collisions = true;
        }
        collision_check[short_name] = 1;
    }
    
    if(!has_collisions) {
        std::cout << "\n✓ SUCCESS: No collisions detected!" << std::endl;
    } else {
        std::cout << "\n✗ FAILURE: Collisions still exist!" << std::endl;
    }
}

int main() {
    test_collision_resolution();
    return 0;
}
