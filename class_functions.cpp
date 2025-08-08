#include "HMMTree.h"

//create fils folder named after the input file or path
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
    folder_matrixs =files_folder + "matrixs"+"/";
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


//test the depend programs exist or not
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


//test the depend programs exist or not
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

//test if the dir exists, or create them
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

    //test the matrixs folder
    if(!dir_exist_opendir(folder_matrixs))
    {
        system_return(system(("mkdir  -m 777 "+folder_matrixs).c_str()));
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

//function deal a fasta sequences input
void HMMTree::prc_fasta_seqs_deal(std::string file_path_name, double identity){
   //test the depend programs exist or not
    test_depend_programs_2(prc_hhsuite,0);

    //use usearch do clusters
    if(!usearch_cluster(file_path_name, identity)){
        output_error_("USEARCH cluster ");
        return;
    }


    //aligen all the squences in the current folder
    if(!align_do_mafft_all_from_file()){
        output_error_("MAFFT align ");
        return;
    }

    //hmmbuild all the fasta files in the folder
    if(!hmm_do_hmmbuild_all_from_file()){
        output_error_("HMMER hmmbuild ");
        return;
    }

        //compute the distance of two hmms
   if(2 == prc_check_profile_HMM_format()){
        if(!hmm_hmmconvert()){
            std::cout<<"Failed to convert the profile HMM files to HMMER2.0 format !"<<std::endl;
            return;
        }
     }
        //compute the distance of two hmms
    if(!pairwise_mode){

        if(!prc_library()){
            output_error_("PRC deal 'prc_library()' ");
            return;
        }

    }else{
        if(!prc_each2()){
                output_error_("PRC deal 'prc_each2()' ");
                return;
            }
    }

/*
    //temple test function to compare each 2 hmms
    if(!matrix_get_each2_hmms_result()){
         output_error_("'matrix_get_each2_hmms_result()' ");
        return;
    }
*/
    //out put the distance matrix to a file
    if(!matrix_mega_out_put_dist_matrix_to_file()){
        output_error_("'matrix_mega_out_put_dist_matrix_to_file()' ");
        return;
    }

    //out put the distance matrix to a file in phylip format
    if(!matrix_phylip_out_put_dist_matrix_to_file()){
        output_error_("'matrix_phylip_out_put_dist_matrix_to_file()'  ");
        return;
    }

    draw_tree_test();
    trees_replace_shorted_names(folder_tree_files);
    system_return(system(("rm -rf "+ folder_clusters).c_str()));
    system_return(system(("rm -rf "+ folder_prcfiles).c_str()));
    return ;
}


//function deal a fasta sequences input
void HMMTree::hhsuite_fasta_seqs_deal(std::string file_path_name, double identity)
{
    //test the depend programs exist or not
    test_depend_programs_2(prc_hhsuite,0);

    //use usearch do clusters
    if(!usearch_cluster(file_path_name, identity)){
        output_error_("USEARCH cluster ");
        return;
    }


    //aligen all the squences in the current folder
    if(!align_do_mafft_all_from_file()){
        output_error_("MAFFT align ");
        return;
    }

    //hmmbuild all the fasta files in the folder
    if(!hhsuite_hhmake_all()){
        output_error_("hhsuite hhmake");
        return;
    }

    if(!hhsuite_hhalign_each2()){
        output_error_("hhsuite hhalign");
        return;
    }

/*
    //temple test function to compare each 2 hmms
    if(!matrix_get_each2_hmms_result()){
         output_error_("'matrix_get_each2_hmms_result()' ");
        return;
    }
*/
    //out put the distance matrix to a file
    if(!matrix_mega_out_put_dist_matrix_to_file()){
        output_error_("'matrix_mega_out_put_dist_matrix_to_file()' ");
        return;
    }

    //out put the distance matrix to a file in phylip format
    if(!matrix_phylip_out_put_dist_matrix_to_file()){
        output_error_("'matrix_phylip_out_put_dist_matrix_to_file()'  ");
        return;
    }

    draw_tree_test();

    system_return(system(("rm -rf "+ folder_clusters).c_str()));
    trees_replace_shorted_names(folder_tree_files);
    return ;
}

//function deal alignments input
void HMMTree::prc_alignments_deal(std::string file_path_name){
     //test the depend programs exist or not
    test_depend_programs_2(prc_hhsuite,1);

    //move the alignments to the aligment folder
    if(!align_move_aligned(file_path_name)){
        output_error_("'align_move_aligned()'  ");
        return;
    }

    //hmmbuild all the fasta files in the folder
    if(!hmm_do_hmmbuild_all_from_file()){
        output_error_("HMMER hmmbuild  ");
        return;
    }

    if(2 == prc_check_profile_HMM_format()){
        if(!hmm_hmmconvert()){
            std::cout<<"Failed to convert the profile HMM files to HMMER2.0 format !"<<std::endl;
            return;
        }
     }
        //compute the distance of two hmms
    if(!pairwise_mode){

        if(!prc_library()){
            output_error_("PRC deal 'prc_library()' ");
            return;
        }

    }else{
        if(!prc_each2()){
                output_error_("PRC deal 'prc_each2()' ");
                return;
            }
    }

/*
    //temple test function to compare each 2 hmms
    if(!matrix_get_each2_hmms_result()){
         output_error_("'matrix_get_each2_hmms_result()' ");
        return;
    }
*/
    //out put the distance matrix to a file
    if(!matrix_mega_out_put_dist_matrix_to_file()){
        output_error_("'matrix_mega_out_put_dist_matrix_to_file()' ");
        return;
    }

    //out put the distance matrix to a file in phylip format
    if(!matrix_phylip_out_put_dist_matrix_to_file()){
        output_error_("'matrix_phylip_out_put_dist_matrix_to_file()'  ");
        return;
    }

    draw_tree_test();
    system_return(system(("rm -rf "+ folder_prcfiles).c_str()));
    trees_replace_shorted_names(folder_tree_files);
    return ;

}



//function deal alignments input
void HMMTree::hhsuite_alignments_deal(std::string file_path_name)
{
    //test the depend programs exist or not
    test_depend_programs_2(prc_hhsuite,1);

    //move the alignments to the aligment folder
    if(!align_move_aligned(file_path_name)){
        output_error_("'align_move_aligned()'  ");
        return;
    }

    //hmmbuild all the fasta files in the folder
    if(!hhsuite_hhmake_all()){
        output_error_("hhsuite hhmake");
        return;
    }

    if(!hhsuite_hhalign_each2()){
        output_error_("hhsuite hhalign");
        return;
    }

/*
    //temple test function to compare each 2 hmms
    if(!matrix_get_each2_hmms_result()){
         output_error_("'matrix_get_each2_hmms_result()' ");
        return;
    }
*/
    //out put the distance matrix to a file
    if(!matrix_mega_out_put_dist_matrix_to_file()){
        output_error_("'matrix_mega_out_put_dist_matrix_to_file()' ");
        return;
    }

    //out put the distance matrix to a file in phylip format
    if(!matrix_phylip_out_put_dist_matrix_to_file()){
        output_error_("'matrix_phylip_out_put_dist_matrix_to_file()'  ");
        return;
    }

    draw_tree_test();
    trees_replace_shorted_names(folder_tree_files);
    return ;
}

//function deal HMMs input
void HMMTree::prc_HMMs_deal(std::string infile_path_and_name){

    //test the depend programs exist or not
    test_depend_programs_2(prc_hhsuite,2);
    //copy the input hmms to the hmms folder
    if(!hmm_copy_hmmfiles(infile_path_and_name)){
        output_error_("hmm_copy_hmmfiles():  ");
        return ;
    }

    if(2 == prc_check_profile_HMM_format()){
        if(!hmm_hmmconvert()){
            std::cout<<"Failed to convert the profile HMM files to HMMER2.0 format !"<<std::endl;
            return;
        }
     }
        //compute the distance of two hmms
    if(!pairwise_mode){

        if(!prc_library()){
            output_error_("PRC deal 'prc_library()' ");
            return;
        }

    }else{
        if(!prc_each2()){
                output_error_("PRC deal 'prc_each2()' ");
                return;
            }
    }
/*
    //temple test function to compare each 2 hmms
    if(!matrix_get_each2_hmms_result()){
         output_error_("'matrix_get_each2_hmms_result()' ");
        return;
    }
*/
    //out put the distance matrix to a file
    if(!matrix_mega_out_put_dist_matrix_to_file()){
        output_error_("'matrix_mega_out_put_dist_matrix_to_file()' ");
        return;
    }

    //out put the distance matrix to a file in phylip format
    if(!matrix_phylip_out_put_dist_matrix_to_file()){
        output_error_("'matrix_phylip_out_put_dist_matrix_to_file()'  ");
        return;
    }
    draw_tree_test();
    system_return(system(("rm -rf "+ folder_prcfiles).c_str()));
    trees_replace_shorted_names(folder_tree_files);
    return ;


}


//function deal HMMs input
void HMMTree::hhsuite_HHMs_deal(std::string infile_path_and_name)
{
    //test the depend programs exist or not
    test_depend_programs_2(prc_hhsuite,2);
    //copy the input hmms to the hmms folder
    if(!hhsuite_copy_hhmfiles(infile_path_and_name)){
        output_error_("hhsuite_copy_hhmfiles():  ");
        return ;
    }

    if(!hhsuite_hhalign_each2()){
        output_error_("hhsuite hhalign");
        return;
    }


    //out put the distance matrix to a file
    if(!matrix_mega_out_put_dist_matrix_to_file()){
        output_error_("'matrix_mega_out_put_dist_matrix_to_file()' ");
        return;
    }

    //out put the distance matrix to a file in phylip format
    if(!matrix_phylip_out_put_dist_matrix_to_file()){
        output_error_("'matrix_phylip_out_put_dist_matrix_to_file()'  ");
        return;
    }
    draw_tree_test();
    trees_replace_shorted_names(folder_tree_files);
    return;
}


//function deal alignments and phmms input
void HMMTree::prc_alignments_phmms_deal(std::string file_path_name, std::string file_path_name_2){
    //test the depend programs exist or not
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


    //move the alignments to the aligment folder
    if(!als_phmms_phhms_move_aligned(file_path_name)){
        output_error_("'als_phmms_phhms_move_aligned()'  ");
        return;
    }

    //hmmbuild all the fasta files in the folder
    if(!hmm_do_hmmbuild_all_from_file()){
        output_error_("HMMER hmmbuild  ");
        return;
    }

    //copy the hmms files to the hmms folder
    if(!hmm_als_phmms_copy_hmmfiles(file_path_name_2)){
        output_error_("'hmm_als_phmms_copy_hmmfiles()' ");
        return;
    }

        //compute the distance of two hmms
    if(2 == prc_check_profile_HMM_format()){
        if(!hmm_hmmconvert()){
            std::cout<<"Failed to convert the profile HMM files to HMMER2.0 format !"<<std::endl;
            return;
        }
     }
        //compute the distance of two hmms
    if(!pairwise_mode){

        if(!prc_library()){
            output_error_("PRC deal 'prc_library()' ");
            return;
        }

    }else{
        if(!prc_each2()){
                output_error_("PRC deal 'prc_each2()' ");
                return;
            }
    }

/*
    //temple test function to compare each 2 hmms
    if(!matrix_get_each2_hmms_result()){
         output_error_("'matrix_get_each2_hmms_result()' ");
        return;
    }
*/
    //out put the distance matrix to a file
    if(!matrix_mega_out_put_dist_matrix_to_file()){
        output_error_("'matrix_mega_out_put_dist_matrix_to_file()' ");
        return;
    }

    //out put the distance matrix to a file in phylip format
    if(!matrix_phylip_out_put_dist_matrix_to_file()){
        output_error_("'matrix_phylip_out_put_dist_matrix_to_file()'  ");
        return;
    }

    draw_tree_test();
    system_return(system(("rm -rf "+ folder_prcfiles).c_str()));
    trees_replace_shorted_names(folder_tree_files);
    return ;
}



//function deal alignments and phhms input
void HMMTree::hhsuite_alignments_phhms_deal(std::string file_path_name, std::string file_path_name_2){
    //test the depend programs exist or not
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

    //move the alignments to the aligment folder
    if(!als_phmms_phhms_move_aligned(file_path_name)){
        output_error_("'als_phmms_phhms_move_aligned()'  ");
        return;
    }

    //hmmbuild all the fasta files in the folder
    if(!hhsuite_hhmake_all()){
        output_error_("hhsuite hhmake");
        return;
    }

    if(!hhsuite_als_phhms_copy_hhmfiles(file_path_name_2)){
        output_error_("'hhsuite_als_phhms_copy_hhmfiles()'  ");
        return;
    }


    if(!hhsuite_hhalign_each2()){
        output_error_("hhsuite hhalign");
        return;
    }

/*
    //temple test function to compare each 2 hmms
    if(!matrix_get_each2_hmms_result()){
         output_error_("'matrix_get_each2_hmms_result()' ");
        return;
    }
*/
    //out put the distance matrix to a file
    if(!matrix_mega_out_put_dist_matrix_to_file()){
        output_error_("'matrix_mega_out_put_dist_matrix_to_file()' ");
        return;
    }

    //out put the distance matrix to a file in phylip format
    if(!matrix_phylip_out_put_dist_matrix_to_file()){
        output_error_("'matrix_phylip_out_put_dist_matrix_to_file()'  ");
        return;
    }

    draw_tree_test();
    trees_replace_shorted_names(folder_tree_files);
    return ;
}
