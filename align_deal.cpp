#include "HMMTree.h"

//seqs align by mafft LINS_i approch
int HMMTree::align_mafft_LINS_i_one(std::string family_name){
	std::string str_cmd = "";
	std::string str_pickname_outof_ext = "";
	int begin = family_name.find_last_of('/') + 1;
	int end = family_name.find_last_of('.');
	str_pickname_outof_ext = family_name.substr(begin, end - begin);
	std::cout<< family_name <<std::endl;
	if(!bool_MAFFT_in_folder){
        str_cmd = "mafft --localpair --maxiterate 1000 " + family_name + " > "+folder_aligned + str_pickname_outof_ext + "_aligned.fasta";
	}else{
        str_cmd = "./mafft --localpair --maxiterate 1000 " + family_name + " > "+folder_aligned + str_pickname_outof_ext + "_aligned.fasta";
    }
	return system(str_cmd.c_str());
}


//aligen all the squences in the current folder
int HMMTree::align_do_mafft_all_from_file()
{
    //if the aligned folder is not empty then empty it
    if(dir_noempty_opendir_readir(folder_aligned))
	{
		//system_return(system("rm ./aligned/*"));
		if(!delete_files(folder_aligned)){
            return 0;
		}
	}

    std::vector<std::string> vec_unalign_filenames;
	if(!get_file_names(folder_unalign_seqs,vec_unalign_filenames,"")){
        return 0;
	}
	if(vec_unalign_filenames.size() < 3){
        return 0;
	}
	size_t unint_unalign = 0;
	bool error_flag = false;
	while (unint_unalign < vec_unalign_filenames.size())
	{
		int return_flag = 0;
		return_flag = align_mafft_LINS_i_one(folder_unalign_seqs+vec_unalign_filenames[unint_unalign]);
		if(system_return(return_flag)){
            error_flag = true;
		}
		if (error_flag)
		{
			vec_unalign_filenames.clear();
			break;
		}
		unint_unalign++;
	}

	if (!error_flag)
	{
		vec_unalign_filenames.clear();
	}
	else {
		return 0;
	}
	return 1;
}


//fasta file format and exists check
bool HMMTree::align_fasta_file_exist_format_check(std::string infile_path_and_name)
{
	//check the file exists or not by test to open it
	std::ifstream ifstream_test_fasta_format;
	ifstream_test_fasta_format.open(infile_path_and_name.c_str());
	if (!ifstream_test_fasta_format.is_open())
	{
		std::cout << "File '" << infile_path_and_name << "' doesn't exist !" << std::endl;
		return false;
	}

	std::string str_temp_one_line = "";
	std::getline(ifstream_test_fasta_format, str_temp_one_line);
	ifstream_test_fasta_format.close();
	//check the '>' word to test the file is fasta format or not
	if (str_temp_one_line.find('>') == -1)
	{
		std::cout << "The content in file '" << infile_path_and_name << "' is not in fasta format !" << std::endl;
		return false;
	}

	std::string str_fasta_ext = infile_path_and_name.substr( infile_path_and_name.find_last_of("."),6);

	if (str_fasta_ext == ".fasta")
	{
		return true;
	}

	std::cout << "'" << infile_path_and_name.c_str() << "' is not in fasta format, the file should be ended by '.fasta' exitemsion !" << std::endl;
	return false;
}

//move the alignments to the aligned folder
int HMMTree::align_move_aligned(std::string file_path_name){

    if (!dir_exist_opendir(file_path_name))
	{
        std::cout<<"Failed to open the file directory !"<<std::endl;
		return 0;
	}
    create_files_folder(file_path_name,prc_hhsuite,1);
    //test if the dir exists, or create them
    dir_exist_or_create(prc_hhsuite,1);

    //if the aligned folder is not empty then empty it
    if(dir_noempty_opendir_readir(folder_aligned))
	{
		//system_return(system("rm ./aligned/*"));
		if(!delete_files(folder_aligned)){
            return 0;
		}
	}

    std::string  path_temp=file_path_name;
    if(path_temp[path_temp.length()-1] != '/'){
        path_temp = path_temp+"/";
    }

/*
    std::string str_move_alenments_cmd= "cp "+file_path_name+"/*.fasta ./aligned/";
    int return_flag = 0;
    return_flag = system(str_move_alenments_cmd.c_str());
    if(system_return(return_flag)){
        return 0;
    }
*/
    copy_files(path_temp,folder_aligned,".fasta");


	return 1;

}

//move the alignments to the aligned folder
int HMMTree::als_phmms_phhms_move_aligned(std::string file_path_name){
    if (!dir_exist_opendir(file_path_name))
	{
        std::cout<<"Failed to open the file directory !"<<std::endl;
		return 0;
	}
    create_files_folder(file_path_name,prc_hhsuite,3);
    //test if the dir exists, or create them
    dir_exist_or_create(prc_hhsuite,1);

    //if the aligned folder is not empty then empty it
    if(dir_noempty_opendir_readir(folder_aligned))
	{
		//system_return(system("rm ./aligned/*"));
		if(!delete_files(folder_aligned)){
            return 0;
		}
	}

    std::string  path_temp=file_path_name;
    if(path_temp[path_temp.length()-1] != '/'){
        path_temp = path_temp+"/";
    }

/*
    std::string str_move_alenments_cmd= "cp "+file_path_name+"/*.fasta ./aligned/";
    int return_flag = 0;
    return_flag = system(str_move_alenments_cmd.c_str());
    if(system_return(return_flag)){
        return 0;
    }
*/
    copy_files(path_temp,folder_aligned,".fasta");


	return 1;
}
