#include "HMMTree.h"
#include <cstdio>

// Create folder named after the input file or path
void HMMTree::create_files_folder(std::string input_file_or_folder_path, int prc_hhsuite, int int_run_mode){
    std::string str_new_folder = "";
    std::string str_input_folder_or_file_name = "";
    std::string str_input_temp = input_file_or_folder_path;

    if(prc_hhsuite == 0){
        switch(int_run_mode){
        case 0:
            str_input_folder_or_file_name =str_input_temp.substr(str_input_temp.find_last_of("/")+1,str_input_temp.find_last_of(".")-str_input_temp.find_last_of("/")-1);
            str_new_folder = "prc_uals_mode_"+str_input_folder_or_file_name;
            break;
        case 1:
            if(str_input_temp[str_input_temp.length() -1] == '/'){
                str_input_temp = str_input_temp.substr(0,str_input_temp.length()-1);
            }
            str_input_folder_or_file_name =str_input_temp.substr(str_input_temp.find_last_of("/")+1,str_input_temp.length()-str_input_temp.find_last_of("/")-1);
            str_new_folder = "prc_als_mode_"+str_input_folder_or_file_name;
            break;
        case 2:
            if(str_input_temp[str_input_temp.length() -1] == '/'){
                str_input_temp = str_input_temp.substr(0,str_input_temp.length()-1);
            }
            str_input_folder_or_file_name =str_input_temp.substr(str_input_temp.find_last_of("/")+1,str_input_temp.length()-str_input_temp.find_last_of("/")-1);
            str_new_folder = "prc_hmms_mode_"+str_input_folder_or_file_name;
            break;
        case 3:
            if(str_input_temp[str_input_temp.length() -1] == '/'){
                str_input_temp = str_input_temp.substr(0,str_input_temp.length()-1);
            }
            str_input_folder_or_file_name =str_input_temp.substr(str_input_temp.find_last_of("/")+1,str_input_temp.length()-str_input_temp.find_last_of("/")-1);
            str_new_folder = "prc_als_phmms_mode_"+str_input_folder_or_file_name;
            break;
        default:
            break;

        }
    }else{
        switch(int_run_mode){
        case 0:
            str_input_folder_or_file_name =str_input_temp.substr(str_input_temp.find_last_of("/")+1,str_input_temp.find_last_of(".")-str_input_temp.find_last_of("/")-1);
            str_new_folder = "hhsuite_uals_mode_"+str_input_folder_or_file_name;
            break;
        case 1:
            if(str_input_temp[str_input_temp.length() -1] == '/'){
                str_input_temp = str_input_temp.substr(0,str_input_temp.length()-1);
            }
            str_input_folder_or_file_name =str_input_temp.substr(str_input_temp.find_last_of("/")+1,str_input_temp.length()-str_input_temp.find_last_of("/")-1);
            str_new_folder = "hhsuite_als_mode_"+str_input_folder_or_file_name;
            break;
        case 2:
            if(str_input_temp[str_input_temp.length() -1] == '/'){
                str_input_temp = str_input_temp.substr(0,str_input_temp.length()-1);
            }
            str_input_folder_or_file_name =str_input_temp.substr(str_input_temp.find_last_of("/")+1,str_input_temp.length()-str_input_temp.find_last_of("/")-1);
            str_new_folder = "hhsuite_hhms_mode_"+str_input_folder_or_file_name;
            break;
        case 3:
            if(str_input_temp[str_input_temp.length() -1] == '/'){
                str_input_temp = str_input_temp.substr(0,str_input_temp.length()-1);
            }
            str_input_folder_or_file_name =str_input_temp.substr(str_input_temp.find_last_of("/")+1,str_input_temp.length()-str_input_temp.find_last_of("/")-1);
            str_new_folder = "hhsuite_als_phhms_mode_"+str_input_folder_or_file_name;
            break;
        default:
            break;

        }
    }


    unsigned int folder_num = 0;
    std::string str_folder_temp =str_new_folder;
    while(dir_exist_opendir(str_new_folder)){
        std::string str_num_temp = uint2str(folder_num);
        str_new_folder=str_folder_temp+"_"+str_num_temp;
        folder_num++;
    }

    system_return(system(("mkdir -m 777 "+str_new_folder).c_str()));

    files_folder = "./"+str_new_folder+"/";
    folder_hmms =files_folder + "hmms"+"/";
    folder_hmms_from_als =files_folder + "hmms_from_als"+"/";
    folder_prcfiles = files_folder + "prcfiles"+"/";
    folder_tree_files =files_folder + "tree_files"+"/";
    folder_matrices =files_folder + "matrices"+"/";
    folder_clusters =files_folder + "clusters"+"/";
    folder_unalign_seqs =files_folder + "unalign_seqs"+"/";
    folder_invalid_clusters =files_folder + "invalid_clusters"+"/";
    folder_aligned =files_folder + "aligned"+"/";
    folder_hmmer2 =files_folder + "hmmer2"+"/";
    folder_hmmer3 =files_folder + "hmmer3"+"/";
    folder_hhms = files_folder + "hhms"+"/";
    folder_hhms_from_als = files_folder + "hhms_from_als"+"/";

    return;
}


// Test if the required program dependencies exist
void HMMTree::test_depend_programs(int mode_num){
    bool all_find = true;
    int prc_ = PRC_exist();
    if( prc_ != 0){
        if(prc_ > 0){
            bool_PRC_in_folder = true;
        }
    }else{
        std::cout<<"Can not find 'prc' program in the PATH variable or folder !"<<std::endl;
        all_find = false;
    }
    // prcX is optional; detect availability
    int prcx_ = PRCX_exist();
    if( prcx_ != 0){
        if(prcx_ > 0){
            bool_PRCX_in_folder = true;
        }
    }
    int usearch_ = 0;
    int mafft_ = 0;
    int hmmbuild_ = 0;
    switch(mode_num){
        case 0:
            //usearch
            usearch_ = USEARCH_exist();
            if( usearch_ != 0){
                if(usearch_ > 0){
                    bool_userach_in_folder = true;
                }
            }else{
                std::cout<<"Can not find 'usearch' program in the PATH variable or folder !"<<std::endl;
                all_find = false;
            }

            //mafft
            mafft_ = MAFFT_exist();
            if( mafft_ != 0){
                if(mafft_ > 0){
                    bool_MAFFT_in_folder = true;
                }
            }else{
                std::cout<<"Can not find 'mafft' program in the PATH variable or folder !"<<std::endl;
                all_find = false;
            }

            //hmmbuild
            hmmbuild_ = HMMER_hmmbuild_exist();
            if( hmmbuild_ != 0){
                if(hmmbuild_ > 0){
                    bool_HMMER_hmmbuild_in_folder = true;
                }
            }else{
                std::cout<<"Can not find 'hmmbuild' program in the PATH variable or folder !"<<std::endl;
                all_find = false;
            }
            break;
        case 1:
            //hmmbuild
            hmmbuild_ = HMMER_hmmbuild_exist();
            if( hmmbuild_ != 0){
                if(hmmbuild_ > 0){
                    bool_HMMER_hmmbuild_in_folder = true;
                }
            }else{
                std::cout<<"Can not find 'hmmbuild' program in the PATH variable or folder !"<<std::endl;
                all_find = false;
            }
            break;
        default:
            break;
    }

    if(!all_find){
        exit(1);
    }
    return ;
}


// Test if the required program dependencies exist
void HMMTree::test_depend_programs_2(int prc_hhsuite, int mode_num){

    bool all_find = true;
    if(prc_hhsuite == 0){
        int PRC_ = PRC_exist();
        if( PRC_ != 0){
            if(PRC_ > 0){
                bool_PRC_in_folder = true;
            }
        }else{
            std::cout<<"Can not find 'PRC' program in the PATH variable or folder !"<<std::endl;
            all_find = false;
        }
        // Detect optional prcX
        int PRCX_ = PRCX_exist();
        if( PRCX_ != 0){
            if(PRCX_ > 0){
                bool_PRCX_in_folder = true;
            }
        }
    }else{
        int hhsuite_hhalign_ = hhsuite_hhalign_exist();
        if( hhsuite_hhalign_ != 0){
            if(hhsuite_hhalign_ > 0){
                bool_hhsuite_hhalign_in_folder = true;
            }
        }else{
            std::cout<<"Can not find 'hhalign' program in the PATH variable or folder !"<<std::endl;
            all_find = false;
        }

        int hhsuite_hhamke_ = hhsuite_hhmake_exist();
        if( hhsuite_hhamke_ != 0){
            if(hhsuite_hhamke_ > 0){
                bool_hhsuite_hhmake_in_folder = true;
            }
        }else{
            std::cout<<"Can not find 'hhmake' program in the PATH variable or folder !"<<std::endl;
            all_find = false;
        }
    }


    int usearch_ = 0;
    int mafft_ = 0;
    int hmmbuild_ = 0;
    switch(mode_num){
        case 0:
            //usearch
            usearch_ = USEARCH_exist();
            if( usearch_ != 0){
                if(usearch_ > 0){
                    bool_userach_in_folder = true;
                }
            }else{
                std::cout<<"Can not find 'usearch' program in the PATH variable or folder !"<<std::endl;
                all_find = false;
            }

            //mafft
            mafft_ = MAFFT_exist();
            if( mafft_ != 0){
                if(mafft_ > 0){
                    bool_MAFFT_in_folder = true;
                }
            }else{
                std::cout<<"Can not find 'mafft' program in the PATH variable or folder !"<<std::endl;
                all_find = false;
            }

            //hmmbuild
            hmmbuild_ = HMMER_hmmbuild_exist();
            if( hmmbuild_ != 0){
                if(hmmbuild_ > 0){
                    bool_HMMER_hmmbuild_in_folder = true;
                }
            }else{
                std::cout<<"Can not find 'hmmbuild' program in the PATH variable or folder !"<<std::endl;
                all_find = false;
            }
            break;
        case 1:
            //hmmbuild
            hmmbuild_ = HMMER_hmmbuild_exist();
            if( hmmbuild_ != 0){
                if(hmmbuild_ > 0){
                    bool_HMMER_hmmbuild_in_folder = true;
                }
            }else{
                std::cout<<"Can not find 'hmmbuild' program in the PATH variable or folder !"<<std::endl;
                all_find = false;
            }
            break;
        default:
            break;
    }

    if(!all_find){
        exit(1);
    }
    return ;
}

// Check if directories exist, otherwise create them
void HMMTree::dir_exist_or_create(int prc_hhsuite, int mode_num){

    if(prc_hhsuite == 0){
        //test the hmms folder
        if(!dir_exist_opendir(folder_hmms))
        {
            system_return(system(("mkdir  -m 777 "+folder_hmms).c_str()));
        }

        //test the hmms folder
        if(!dir_exist_opendir(folder_prcfiles))
        {
            system_return(system(("mkdir  -m 777 "+folder_prcfiles).c_str()));
        }


    }else{
         //test the hmms folder
        if(!dir_exist_opendir(folder_hhms))
        {
            system_return(system(("mkdir  -m 777 "+folder_hhms).c_str()));
        }
    }
    //test the tree_files folder
    if(!dir_exist_opendir(folder_tree_files))
    {
        system_return(system(("mkdir  -m 777 "+folder_tree_files).c_str()));
    }

    //test the matrices folder
    if(!dir_exist_opendir(folder_matrices))
    {
        system_return(system(("mkdir  -m 777 "+folder_matrices).c_str()));
    }
    switch(mode_num){
        case 0:
            //test the clusters folder
            if(!dir_exist_opendir(folder_clusters))
            {
                system_return(system(("mkdir  -m 777 "+folder_clusters).c_str()));
            }
            //test the unalign_seqs folder
            if(!dir_exist_opendir(folder_unalign_seqs))
            {
                system_return(system(("mkdir  -m 777 "+folder_unalign_seqs).c_str()));
            }

            //test the invalid_clusters folder
            if(!dir_exist_opendir(folder_invalid_clusters))
            {
                system_return(system(("mkdir  -m 777 "+folder_invalid_clusters).c_str()));
            }

            //test the aligned folder
            if(!dir_exist_opendir(folder_aligned))
            {
                system_return(system(("mkdir  -m 777 "+folder_aligned).c_str()));
            }

            break;
        case 1:
            //test the aligned folder
            if(!dir_exist_opendir(folder_aligned))
            {
                system_return(system(("mkdir  -m 777 "+folder_aligned).c_str()));
            }
            break;
        default :
            break;
    }

}

// Process input FASTA sequences using PRC
void HMMTree::process_fasta_sequences(std::string file_path_name, double identity){
   // Test if the required program dependencies exist
    test_depend_programs_2(prc_hhsuite,0);

    // Use USEARCH for clustering
    if(!usearch_cluster(file_path_name, identity)){
        output_error_("USEARCH cluster ");
        return;
    }


    // Align all sequences in the current folder
    if(!align_do_mafft_all_from_file()){
        output_error_("MAFFT align ");
        return;
    }

    // Build HMM profiles for all FASTA files in the folder
    if(!hmm_do_hmmbuild_all_from_file()){
        output_error_("HMMER hmmbuild ");
        return;
    }

    // Resolve backend and convert only if needed
    {
        int fmt = prc_check_profile_HMM_format();
        bool has_hmmer3 = (fmt == 2);
        std::cout << "HMM profile format detected: " << (has_hmmer3 ? "HMMER3" : "HMMER2") << std::endl;
        // Choose effective backend
        if (prc_backend == PRC_BACKEND_LEGACY) {
            effective_prc_backend = PRC_BACKEND_LEGACY;
        } else if (prc_backend == PRC_BACKEND_PRCX) {
            // Forced prcX: use prcX if available regardless of input format; else fallback to legacy
            if (bool_PRCX_in_folder || PRCX_exist() != 0) effective_prc_backend = PRC_BACKEND_PRCX; else effective_prc_backend = PRC_BACKEND_LEGACY;
        } else {
            // auto
            if (has_hmmer3 && (bool_PRCX_in_folder || PRCX_exist() != 0)) effective_prc_backend = PRC_BACKEND_PRCX; else effective_prc_backend = PRC_BACKEND_LEGACY;
        }
        std::cout << "PRC backend selected: " << (effective_prc_backend == PRC_BACKEND_PRCX ? "prcX" : "legacy (convert+prc)") << std::endl;
        if (effective_prc_backend == PRC_BACKEND_LEGACY && has_hmmer3){
            if(!hmm_hmmconvert()){
                std::cout<<"Failed to convert the profile HMM files to HMMER2.0 format !"<<std::endl;
                return;
            }
        }
        // If prcX is used, rename the output base folder prefix from prc_ to prcx_
        if (effective_prc_backend == PRC_BACKEND_PRCX) {
            if (files_folder.rfind("./prc_", 0) == 0) {
                std::string old_dir = files_folder;
                if (!old_dir.empty() && old_dir.back() == '/') old_dir.pop_back();
                std::string new_dir = std::string("./prcx_") + old_dir.substr(6);
                if (rename(old_dir.c_str(), new_dir.c_str()) == 0) {
                    files_folder = new_dir + "/";
                    // Recompute subfolder paths
                    folder_hmms =files_folder + "hmms"+"/";
                    folder_hmms_from_als =files_folder + "hmms_from_als"+"/";
                    folder_prcfiles = files_folder + "prcfiles"+"/";
                    folder_tree_files =files_folder + "tree_files"+"/";
                    folder_matrices =files_folder + "matrices"+"/";
                    folder_clusters =files_folder + "clusters"+"/";
                    folder_unalign_seqs =files_folder + "unalign_seqs"+"/";
                    folder_invalid_clusters =files_folder + "invalid_clusters"+"/";
                    folder_aligned =files_folder + "aligned"+"/";
                    folder_hmmer2 =files_folder + "hmmer2"+"/";
                    folder_hmmer3 =files_folder + "hmmer3"+"/";
                    folder_hhms = files_folder + "hhms"+"/";
                    folder_hhms_from_als = files_folder + "hhms_from_als"+"/";
                }
            }
        }
    }

    // Time the PRC analysis phase
    struct timeb prcStartTime, prcEndTime;
    ftime(&prcStartTime);
    std::cout << "PRC processing: distance calculations..." << std::endl;

    // Compute distances between HMM pairs
    if(!pairwise_mode){

        if(!prc_library()){
            output_error_("PRC processing 'prc_library()' ");
            return;
        }

    }else{
        if(!prc_each2()){
                output_error_("PRC processing 'prc_each2()' ");
                return;
            }
    }

    ftime(&prcEndTime);
    std::cout << "PRC analysis completed in: " << format_time_duration((prcEndTime.time-prcStartTime.time)*1000 + (prcEndTime.millitm - prcStartTime.millitm)) << std::endl;

/*
    // Test function to compare each pair of HMMs
    if(!matrix_get_each2_hmms_result()){
         output_error_("'matrix_get_each2_hmms_result()' ");
        return;
    }
*/
    // Output the distance matrix to a MEGA format file
    if(!matrix_mega_output_dist_matrix_to_file()){
        output_error_("'matrix_mega_output_dist_matrix_to_file()' ");
        return;
    }

    // Output the distance matrix to a PHYLIP format file
    if(!matrix_phylip_output_dist_matrix_to_file()){
        output_error_("'matrix_phylip_output_dist_matrix_to_file()'  ");
        return;
    }

    // Time the phylogenetic tree building phase
    struct timeb treeStartTime, treeEndTime;
    ftime(&treeStartTime);
    std::cout << "Building phylogenetic trees..." << std::endl;
    
    draw_tree_test();
    
    ftime(&treeEndTime);
    std::cout << "Phylogenetic tree building completed in: " << format_time_duration((treeEndTime.time-treeStartTime.time)*1000 + (treeEndTime.millitm - treeStartTime.millitm)) << std::endl;
    trees_replace_shorted_names(folder_tree_files);
    system_return(system(("rm -rf "+ folder_clusters).c_str()));
    system_return(system(("rm -rf "+ folder_prcfiles).c_str()));
    return ;
}


// Process input FASTA sequences using HH-suite
void HMMTree::process_hhsuite_fasta_sequences(std::string file_path_name, double identity)
{
    // Test if the required program dependencies exist
    test_depend_programs_2(prc_hhsuite,0);

    // Use USEARCH for clustering
    if(!usearch_cluster(file_path_name, identity)){
        output_error_("USEARCH cluster ");
        return;
    }


    // Align all sequences in the current folder
    if(!align_do_mafft_all_from_file()){
        output_error_("MAFFT align ");
        return;
    }

    // Build HHM profiles for all alignments
    if(!hhsuite_hhmake_all()){
        output_error_("hhsuite hhmake");
        return;
    }

    // Calculate distances between HHM pairs
    if(!hhsuite_hhalign_each2()){
        output_error_("hhsuite hhalign");
        return;
    }

/*
    // Test function to compare each pair of HMMs
    if(!matrix_get_each2_hmms_result()){
         output_error_("'matrix_get_each2_hmms_result()' ");
        return;
    }
*/
    // Output the distance matrix to a MEGA format file
    if(!matrix_mega_output_dist_matrix_to_file()){
        output_error_("'matrix_mega_output_dist_matrix_to_file()' ");
        return;
    }

    // Output the distance matrix to a PHYLIP format file
    if(!matrix_phylip_output_dist_matrix_to_file()){
        output_error_("'matrix_phylip_output_dist_matrix_to_file()'  ");
        return;
    }

    draw_tree_test();

    system_return(system(("rm -rf "+ folder_clusters).c_str()));
    trees_replace_shorted_names(folder_tree_files);
    return ;
}

// Process input alignments using PRC
void HMMTree::process_prc_alignments(std::string file_path_name){
    // Test if the required program dependencies exist
    test_depend_programs_2(prc_hhsuite,1);

    // Move the alignments to the alignment folder
    if(!align_move_aligned(file_path_name)){
        output_error_("'align_move_aligned()'  ");
        return;
    }

    // Build HMM profiles for all alignments
    if(!hmm_do_hmmbuild_all_from_file()){
        output_error_("HMMER hmmbuild  ");
        return;
    }

    {
        int fmt = prc_check_profile_HMM_format();
        bool has_hmmer3 = (fmt == 2);
        std::cout << "HMM profile format detected: " << (has_hmmer3 ? "HMMER3" : "HMMER2") << std::endl;
        if (prc_backend == PRC_BACKEND_LEGACY) {
            effective_prc_backend = PRC_BACKEND_LEGACY;
        } else if (prc_backend == PRC_BACKEND_PRCX) {
            if (bool_PRCX_in_folder || PRCX_exist() != 0) effective_prc_backend = PRC_BACKEND_PRCX; else effective_prc_backend = PRC_BACKEND_LEGACY;
        } else {
            if (has_hmmer3 && (bool_PRCX_in_folder || PRCX_exist() != 0)) effective_prc_backend = PRC_BACKEND_PRCX; else effective_prc_backend = PRC_BACKEND_LEGACY;
        }
        std::cout << "PRC backend selected: " << (effective_prc_backend == PRC_BACKEND_PRCX ? "prcX" : "legacy (convert+prc)") << std::endl;
        if (effective_prc_backend == PRC_BACKEND_LEGACY && has_hmmer3){
            if(!hmm_hmmconvert()){
                std::cout<<"Failed to convert the profile HMM files to HMMER2.0 format !"<<std::endl;
                return;
            }
        }
        if (effective_prc_backend == PRC_BACKEND_PRCX) {
            if (files_folder.rfind("./prc_", 0) == 0) {
                std::string old_dir = files_folder;
                if (!old_dir.empty() && old_dir.back() == '/') old_dir.pop_back();
                std::string new_dir = std::string("./prcx_") + old_dir.substr(6);
                if (rename(old_dir.c_str(), new_dir.c_str()) == 0) {
                    files_folder = new_dir + "/";
                    folder_hmms =files_folder + "hmms"+"/";
                    folder_hmms_from_als =files_folder + "hmms_from_als"+"/";
                    folder_prcfiles = files_folder + "prcfiles"+"/";
                    folder_tree_files =files_folder + "tree_files"+"/";
                    folder_matrices =files_folder + "matrices"+"/";
                    folder_clusters =files_folder + "clusters"+"/";
                    folder_unalign_seqs =files_folder + "unalign_seqs"+"/";
                    folder_invalid_clusters =files_folder + "invalid_clusters"+"/";
                    folder_aligned =files_folder + "aligned"+"/";
                    folder_hmmer2 =files_folder + "hmmer2"+"/";
                    folder_hmmer3 =files_folder + "hmmer3"+"/";
                    folder_hhms = files_folder + "hhms"+"/";
                    folder_hhms_from_als = files_folder + "hhms_from_als"+"/";
                }
            }
        }
    }

    // Time the PRC analysis phase
    struct timeb prcStartTime, prcEndTime;
    ftime(&prcStartTime);
    std::cout << "PRC processing: distance calculations..." << std::endl;

    // Compute distances between HMM pairs
    if(!pairwise_mode){

        if(!prc_library()){
            output_error_("PRC processing 'prc_library()' ");
            return;
        }

    }else{
        if(!prc_each2()){
                output_error_("PRC processing 'prc_each2()' ");
                return;
            }
    }

    ftime(&prcEndTime);
    std::cout << "PRC analysis completed in: " << format_time_duration((prcEndTime.time-prcStartTime.time)*1000 + (prcEndTime.millitm - prcStartTime.millitm)) << std::endl;

/*
    // Test function to compare each pair of HMMs
    if(!matrix_get_each2_hmms_result()){
         output_error_("'matrix_get_each2_hmms_result()' ");
        return;
    }
*/
    // Output the distance matrix to a MEGA format file
    if(!matrix_mega_output_dist_matrix_to_file()){
        output_error_("'matrix_mega_output_dist_matrix_to_file()' ");
        return;
    }

    // Output the distance matrix to a PHYLIP format file
    if(!matrix_phylip_output_dist_matrix_to_file()){
        output_error_("'matrix_phylip_output_dist_matrix_to_file()'  ");
        return;
    }

    // Time the phylogenetic tree building phase
    struct timeb treeStartTime, treeEndTime;
    ftime(&treeStartTime);
    std::cout << "Building phylogenetic trees..." << std::endl;
    
    draw_tree_test();
    
    ftime(&treeEndTime);
    std::cout << "Phylogenetic tree building completed in: " << format_time_duration((treeEndTime.time-treeStartTime.time)*1000 + (treeEndTime.millitm - treeStartTime.millitm)) << std::endl;
    
    system_return(system(("rm -rf "+ folder_prcfiles).c_str()));
    trees_replace_shorted_names(folder_tree_files);
    return ;
}


// Process input alignments using HH-suite
void HMMTree::process_hhsuite_alignments(std::string file_path_name)
{
    // Test if the required program dependencies exist
    test_depend_programs_2(prc_hhsuite,1);

    // Move the alignments to the alignment folder
    if(!align_move_aligned(file_path_name)){
        output_error_("'align_move_aligned()'  ");
        return;
    }

    // Build HHM profiles for all alignments
    if(!hhsuite_hhmake_all()){
        output_error_("hhsuite hhmake");
        return;
    }

    // Calculate distances between HHM pairs
    if(!hhsuite_hhalign_each2()){
        output_error_("hhsuite hhalign");
        return;
    }

/*
    // Test function to compare each pair of HMMs
    if(!matrix_get_each2_hmms_result()){
         output_error_("'matrix_get_each2_hmms_result()' ");
        return;
    }
*/
    // Output the distance matrix to a MEGA format file
    if(!matrix_mega_output_dist_matrix_to_file()){
        output_error_("'matrix_mega_output_dist_matrix_to_file()' ");
        return;
    }

    // Output the distance matrix to a PHYLIP format file
    if(!matrix_phylip_output_dist_matrix_to_file()){
        output_error_("'matrix_phylip_output_dist_matrix_to_file()'  ");
        return;
    }

    draw_tree_test();
    trees_replace_shorted_names(folder_tree_files);
    return ;
}

// Process input HMM profiles using PRC
void HMMTree::process_prc_HMMs(std::string infile_path_and_name){

    // Test if the required program dependencies exist
    test_depend_programs_2(prc_hhsuite,2);
    // Copy the input HMMs to the hmms folder
    if(!hmm_copy_hmmfiles(infile_path_and_name)){
        output_error_("hmm_copy_hmmfiles():  ");
        return ;
    }

    // Decide effective PRC backend based on user preference, availability, and input format
    // 0: HMMER2 ok; 2: HMMER3 detected
    int fmt = prc_check_profile_HMM_format();
    bool has_hmmer3 = (fmt == 2);
    std::cout << "HMM profile format detected: " << (has_hmmer3 ? "HMMER3" : "HMMER2") << std::endl;
        if (prc_backend == PRC_BACKEND_LEGACY && has_hmmer3) {
            effective_prc_backend = PRC_BACKEND_LEGACY;
        } else if (prc_backend == PRC_BACKEND_PRCX) {
            // Forced prcX regardless of input format; fallback to legacy if prcX unavailable
            if (bool_PRCX_in_folder || PRCX_exist() != 0) effective_prc_backend = PRC_BACKEND_PRCX; 
            else effective_prc_backend = PRC_BACKEND_LEGACY;
        } else {
            // auto: prefer prcX if hmmer3 and prcX available; else legacy
            if (has_hmmer3 && (bool_PRCX_in_folder || PRCX_exist() != 0)) effective_prc_backend = PRC_BACKEND_PRCX; 
            else effective_prc_backend = PRC_BACKEND_LEGACY;
    }
    std::cout << "PRC backend selected: " << (effective_prc_backend == PRC_BACKEND_PRCX ? "prcX" : "legacy (convert+prc)") << std::endl;

    if (effective_prc_backend == PRC_BACKEND_LEGACY && has_hmmer3){
        if(!hmm_hmmconvert()){
            std::cout<<"Failed to convert the profile HMM files to HMMER2.0 format !"<<std::endl;
            return;
        }
    }
    if (effective_prc_backend == PRC_BACKEND_PRCX) {
        if (files_folder.rfind("./prc_", 0) == 0) {
            std::string old_dir = files_folder;
            if (!old_dir.empty() && old_dir.back() == '/') old_dir.pop_back();
            std::string new_dir = std::string("./prcx_") + old_dir.substr(6);
            if (rename(old_dir.c_str(), new_dir.c_str()) == 0) {
                files_folder = new_dir + "/";
                folder_hmms =files_folder + "hmms"+"/";
                folder_hmms_from_als =files_folder + "hmms_from_als"+"/";
                folder_prcfiles = files_folder + "prcfiles"+"/";
                folder_tree_files =files_folder + "tree_files"+"/";
                folder_matrices =files_folder + "matrices"+"/";
                folder_clusters =files_folder + "clusters"+"/";
                folder_unalign_seqs =files_folder + "unalign_seqs"+"/";
                folder_invalid_clusters =files_folder + "invalid_clusters"+"/";
                folder_aligned =files_folder + "aligned"+"/";
                folder_hmmer2 =files_folder + "hmmer2"+"/";
                folder_hmmer3 =files_folder + "hmmer3"+"/";
                folder_hhms = files_folder + "hhms"+"/";
                folder_hhms_from_als = files_folder + "hhms_from_als"+"/";
            }
        }
    }

    // Time the PRC analysis phase
    struct timeb prcStartTime, prcEndTime;
    ftime(&prcStartTime);
    std::cout << "PRC processing: Processing distance calculations..." << std::endl;

    // Compute distances between HMM pairs
    if(!pairwise_mode){

        if(!prc_library()){
            output_error_("PRC processing 'prc_library()' ");
            return;
        }

    }else{
        if(!prc_each2()){
                output_error_("PRC processing 'prc_each2()' ");
                return;
            }
    }

    ftime(&prcEndTime);
    std::cout << "PRC analysis completed in: " << format_time_duration((prcEndTime.time-prcStartTime.time)*1000 + (prcEndTime.millitm - prcStartTime.millitm)) << std::endl;

/*
    // Test function to compare each pair of HMMs
    if(!matrix_get_each2_hmms_result()){
         output_error_("'matrix_get_each2_hmms_result()' ");
        return;
    }
*/
    // Output the distance matrix to a MEGA format file
    if(!matrix_mega_output_dist_matrix_to_file()){
        output_error_("'matrix_mega_output_dist_matrix_to_file()' ");
        return;
    }

    // Output the distance matrix to a PHYLIP format file
    if(!matrix_phylip_output_dist_matrix_to_file()){
        output_error_("'matrix_phylip_output_dist_matrix_to_file()'  ");
        return;
    }

    // Time the phylogenetic tree building phase
    struct timeb treeStartTime, treeEndTime;
    ftime(&treeStartTime);
    std::cout << "Building phylogenetic trees..." << std::endl;
    
    draw_tree_test();
    
    ftime(&treeEndTime);
    std::cout << "Phylogenetic tree building completed in: " << format_time_duration((treeEndTime.time-treeStartTime.time)*1000 + (treeEndTime.millitm - treeStartTime.millitm)) << std::endl;
    
    system_return(system(("rm -rf "+ folder_prcfiles).c_str()));
    trees_replace_shorted_names(folder_tree_files);
    return ;
}


// Process input HHM profiles using HH-suite
void HMMTree::process_hhsuite_HHMs(std::string infile_path_and_name)
{
    // Test if the required program dependencies exist
    test_depend_programs_2(prc_hhsuite,2);
    // Copy the input HHMs to the hhms folder
    if(!hhsuite_copy_hhmfiles(infile_path_and_name)){
        output_error_("hhsuite_copy_hhmfiles():  ");
        return ;
    }

    // Calculate distances between HHM pairs
    if(!hhsuite_hhalign_each2()){
        output_error_("hhsuite hhalign");
        return;
    }


    // Output the distance matrix to a MEGA format file
    if(!matrix_mega_output_dist_matrix_to_file()){
        output_error_("'matrix_mega_output_dist_matrix_to_file()' ");
        return;
    }

    // Output the distance matrix to a PHYLIP format file
    if(!matrix_phylip_output_dist_matrix_to_file()){
        output_error_("'matrix_phylip_output_dist_matrix_to_file()'  ");
        return;
    }
    draw_tree_test();
    trees_replace_shorted_names(folder_tree_files);
    return;
}


// Process both input alignments and profile HMMs using PRC
void HMMTree::process_prc_alignments_phmms(std::string file_path_name, std::string file_path_name_2){
    // Test if the required program dependencies exist
    test_depend_programs_2(prc_hhsuite,1);


    if (!dir_exist_opendir(file_path_name))
	{
        std::cout<<"Failed to open the file directory '"<<file_path_name <<"' !"<<std::endl;
		return;
	}

    if (!dir_exist_opendir(file_path_name_2))
	{
        std::cout<<"Failed to open the file directory '"<<file_path_name_2 <<"' !"<<std::endl;
		return;
	}


    // Move the alignments to the alignment folder
    if(!als_phmms_phhms_move_aligned(file_path_name)){
        output_error_("'als_phmms_phhms_move_aligned()'  ");
        return;
    }

    // Build HMM profiles for all alignments
    if(!hmm_do_hmmbuild_all_from_file()){
        output_error_("HMMER hmmbuild  ");
        return;
    }

    // Copy the HMM files to the hmms folder
    if(!hmm_als_phmms_copy_hmmfiles(file_path_name_2)){
        output_error_("'hmm_als_phmms_copy_hmmfiles()' ");
        return;
    }

    // Resolve backend and convert only if needed
    {
        int fmt = prc_check_profile_HMM_format();
        bool has_hmmer3 = (fmt == 2);
        std::cout << "HMM profile format detected: " << (has_hmmer3 ? "HMMER3" : "HMMER2") << std::endl;
        if (prc_backend == PRC_BACKEND_LEGACY) {
            effective_prc_backend = PRC_BACKEND_LEGACY;
        } else if (prc_backend == PRC_BACKEND_PRCX) {
            if (bool_PRCX_in_folder || PRCX_exist() != 0) effective_prc_backend = PRC_BACKEND_PRCX; else effective_prc_backend = PRC_BACKEND_LEGACY;
        } else {
            if (has_hmmer3 && (bool_PRCX_in_folder || PRCX_exist() != 0)) effective_prc_backend = PRC_BACKEND_PRCX; else effective_prc_backend = PRC_BACKEND_LEGACY;
        }
        std::cout << "PRC backend selected: " << (effective_prc_backend == PRC_BACKEND_PRCX ? "prcX" : "legacy (convert+prc)") << std::endl;
        if (effective_prc_backend == PRC_BACKEND_LEGACY && has_hmmer3){
            if(!hmm_hmmconvert()){
                std::cout<<"Failed to convert the profile HMM files to HMMER2.0 format !"<<std::endl;
                return;
            }
        }
        if (effective_prc_backend == PRC_BACKEND_PRCX) {
            if (files_folder.rfind("./prc_", 0) == 0) {
                std::string old_dir = files_folder;
                if (!old_dir.empty() && old_dir.back() == '/') old_dir.pop_back();
                std::string new_dir = std::string("./prcx_") + old_dir.substr(6);
                if (rename(old_dir.c_str(), new_dir.c_str()) == 0) {
                    files_folder = new_dir + "/";
                    folder_hmms =files_folder + "hmms"+"/";
                    folder_hmms_from_als =files_folder + "hmms_from_als"+"/";
                    folder_prcfiles = files_folder + "prcfiles"+"/";
                    folder_tree_files =files_folder + "tree_files"+"/";
                    folder_matrices =files_folder + "matrices"+"/";
                    folder_clusters =files_folder + "clusters"+"/";
                    folder_unalign_seqs =files_folder + "unalign_seqs"+"/";
                    folder_invalid_clusters =files_folder + "invalid_clusters"+"/";
                    folder_aligned =files_folder + "aligned"+"/";
                    folder_hmmer2 =files_folder + "hmmer2"+"/";
                    folder_hmmer3 =files_folder + "hmmer3"+"/";
                    folder_hhms = files_folder + "hhms"+"/";
                    folder_hhms_from_als = files_folder + "hhms_from_als"+"/";
                }
            }
        }
    }

    // Time the PRC analysis phase
    struct timeb prcStartTime, prcEndTime;
    ftime(&prcStartTime);
    std::cout << "PRC processing: distance calculations..." << std::endl;

    // Compute distances between HMM pairs
    if(!pairwise_mode){

        if(!prc_library()){
            output_error_("PRC processing 'prc_library()' ");
            return;
        }

    }else{
        if(!prc_each2()){
                output_error_("PRC processing 'prc_each2()' ");
                return;
            }
    }

    ftime(&prcEndTime);
    std::cout << "PRC analysis completed in: " << format_time_duration((prcEndTime.time-prcStartTime.time)*1000 + (prcEndTime.millitm - prcStartTime.millitm)) << std::endl;

/*
    // Test function to compare each pair of HMMs
    if(!matrix_get_each2_hmms_result()){
         output_error_("'matrix_get_each2_hmms_result()' ");
        return;
    }
*/
    // Output the distance matrix to a MEGA format file
    if(!matrix_mega_output_dist_matrix_to_file()){
        output_error_("'matrix_mega_output_dist_matrix_to_file()' ");
        return;
    }

    // Output the distance matrix to a PHYLIP format file
    if(!matrix_phylip_output_dist_matrix_to_file()){
        output_error_("'matrix_phylip_output_dist_matrix_to_file()'  ");
        return;
    }

    // Time the phylogenetic tree building phase
    struct timeb treeStartTime, treeEndTime;
    ftime(&treeStartTime);
    std::cout << "Building phylogenetic trees..." << std::endl;
    
    draw_tree_test();
    
    ftime(&treeEndTime);
    std::cout << "Phylogenetic tree building completed in: " << format_time_duration((treeEndTime.time-treeStartTime.time)*1000 + (treeEndTime.millitm - treeStartTime.millitm)) << std::endl;
    system_return(system(("rm -rf "+ folder_prcfiles).c_str()));
    trees_replace_shorted_names(folder_tree_files);
    return ;
}


// Process both input alignments and profile HHMs using HH-suite
void HMMTree::process_hhsuite_alignments_phhms(std::string file_path_name, std::string file_path_name_2){
    // Test if the required program dependencies exist
    test_depend_programs_2(prc_hhsuite,1);

    if (!dir_exist_opendir(file_path_name))
	{
        std::cout<<"Failed to open the file directory '"<<file_path_name <<"' !"<<std::endl;
		return;
	}

    if (!dir_exist_opendir(file_path_name_2))
	{
        std::cout<<"Failed to open the file directory '"<<file_path_name_2 <<"' !"<<std::endl;
		return;
	}

    // Move the alignments to the alignment folder
    if(!als_phmms_phhms_move_aligned(file_path_name)){
        output_error_("'als_phmms_phhms_move_aligned()'  ");
        return;
    }

    // Build HHM profiles for all alignments
    if(!hhsuite_hhmake_all()){
        output_error_("hhsuite hhmake");
        return;
    }

    // Copy the HHM files to the hhms folder
    if(!hhsuite_als_phhms_copy_hhmfiles(file_path_name_2)){
        output_error_("'hhsuite_als_phhms_copy_hhmfiles()'  ");
        return;
    }

    // Calculate distances between HHM pairs
    if(!hhsuite_hhalign_each2()){
        output_error_("hhsuite hhalign");
        return;
    }

/*
    // Test function to compare each pair of HMMs
    if(!matrix_get_each2_hmms_result()){
         output_error_("'matrix_get_each2_hmms_result()' ");
        return;
    }
*/
    // Output the distance matrix to a MEGA format file
    if(!matrix_mega_output_dist_matrix_to_file()){
        output_error_("'matrix_mega_output_dist_matrix_to_file()' ");
        return;
    }

    // Output the distance matrix to a PHYLIP format file
    if(!matrix_phylip_output_dist_matrix_to_file()){
        output_error_("'matrix_phylip_output_dist_matrix_to_file()'  ");
        return;
    }

    draw_tree_test();
    trees_replace_shorted_names(folder_tree_files);
    return ;
}

