#include "HMMTree.h"
//hmmbuild a hmm model
int HMMTree::hmm_hmmbuild_one(std::string family_name)
{
	std::string str_cmd = "";
	std::string str_pickname_outof_ext = "";
	int begin = family_name.find_last_of('/') + 1;
	int end = family_name.find_last_of('_');
	str_pickname_outof_ext = family_name.substr(begin, end - begin);
	//str_pickname_outof_ext = family_name.substr(family_name.find_last_of('/') + 1, family_name.find_last_of('_'));
	if(!bool_HMMER_hmmbuild_in_folder){
        str_cmd = " hmmbuild  --informat afa -n " + str_pickname_outof_ext + " " +folder_hmms+ str_pickname_outof_ext + ".hmm " + family_name;
	}else{
        str_cmd = "./hmmbuild  --informat afa -n " + str_pickname_outof_ext + " " +folder_hmms+ str_pickname_outof_ext + ".hmm " + family_name;
	}

	return system(str_cmd.c_str());
}

//hmmbuild all the fasta files in the folder
int HMMTree::hmm_do_hmmbuild_all_from_file()
{
    if(dir_noempty_opendir_readir(folder_hmms))
	{
		//system_return(system("rm ./hmms/*.hmm"));
		if(!delete_files(folder_hmms)){
            return 0;
		}
	}
	std::vector<std::string> vec_aligned_filenames;
	if(!get_file_names(folder_aligned,vec_aligned_filenames,"")){
        return 0;
	}
	if((vec_aligned_filenames.size() < 3) && (!HMMTree::als_phmms_phhms) ){
        return 0;
	}


	size_t unint_aligned = 0;
	bool error_flag = false;
	while (unint_aligned < vec_aligned_filenames.size())
	{
		if(!align_fasta_file_exist_format_check(folder_aligned+vec_aligned_filenames[unint_aligned]))
		{
			vec_aligned_filenames.clear();
			error_flag = true;
			break;
		}
		int return_flag = 0;
		return_flag = hmm_hmmbuild_one(folder_aligned+vec_aligned_filenames[unint_aligned]);
		if(system_return(return_flag)){
            error_flag = true;
		}
		if (error_flag)
		{
			vec_aligned_filenames.clear();
			break;
		}
		unint_aligned++;
	}
	if (!error_flag)
	{
		vec_aligned_filenames.clear();
	}
	else {
		return 0;
	}
	return 1;
}


//divide hmms into single hmm files and the hmms in the file must be ended by "//"
int HMMTree::hmm_divide_hmms_to_single_hmm(std::string infile_path_and_name, std::string outfiles_path, std::string str_key)
{
	if (!file_exists_and_empty_check(infile_path_and_name))
	{
		return 0;
	}

	if(dir_noempty_opendir_readir(folder_hmms))
	{
		//system_return(system("rm ./hmms/*.hmm"));
		if(!delete_files(folder_hmms)){
            return 0;
		}
	}

	std::string str_temp_one_line_file = "";
	std::string str_temp_hmmbegin_name_line = "";

	std::string str_temp_hmm_outfiles_path_and_nmae = "";

	bool bool_flag_hmm_begin = true;
	bool bool_flag_ofstream_opened = false;
	bool boll_flag_all_divide_true_or_only_one_false = true;
	bool error_flag = false;
	//check the divide model : all or only one kind
	if (str_key.length() > 0)
	{
		boll_flag_all_divide_true_or_only_one_false = false;
	}

	std::ifstream ifstream_file_hmms;			//the file steam variable to process the input hmms file
	ifstream_file_hmms.open(infile_path_and_name.c_str());

	//read the input file and divide the hmms into single hmm files
	std::ofstream ofstream_file_hmm;
	while (std::getline(ifstream_file_hmms, str_temp_one_line_file))
	{
		if (bool_flag_hmm_begin)
		{

			if(str_temp_one_line_file.find("HMMER") == -1)
			{
				ifstream_file_hmms.close();
				std::cout<<"There is an format error in input file of hmm models!"<<std::endl;
				error_flag = true;
				break;
			}
			bool_flag_hmm_begin = false;
			std::getline(ifstream_file_hmms, str_temp_hmmbegin_name_line);
			int begin = str_temp_hmmbegin_name_line.find_first_of(" ");
			int end = str_temp_hmmbegin_name_line.find_last_not_of(" ");
			std::string str_hmm_name = str_temp_hmmbegin_name_line.substr(begin, end - begin + 1);
			begin = str_hmm_name.find_first_not_of(" ");
			str_hmm_name = str_hmm_name.substr(begin, end - begin + 1);


			str_temp_hmm_outfiles_path_and_nmae = outfiles_path + str_hmm_name;
			//std::cout << str_temp_hmm_outfiles_path_and_nmae << std::endl;


			if (boll_flag_all_divide_true_or_only_one_false || (str_hmm_name.find(str_key) != -1))
			{
				ofstream_file_hmm.open(str_temp_hmm_outfiles_path_and_nmae.c_str(),std::ios_base::out|std::ios_base::trunc);
                ofstream_file_hmm.close();
                chmod(str_temp_hmm_outfiles_path_and_nmae.c_str(),S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP|S_IWOTH|S_IROTH);
                ofstream_file_hmm.open(str_temp_hmm_outfiles_path_and_nmae.c_str(),std::ios_base::out|std::ios_base::trunc);
				if (!ofstream_file_hmm.is_open())
				{
					std::cout << "Failed to open the output hmm file!" << std::endl;
					ifstream_file_hmms.close();
					return 0;
				}

				ofstream_file_hmm << str_temp_one_line_file << std::endl;
				ofstream_file_hmm << str_temp_hmmbegin_name_line << std::endl;
				bool_flag_ofstream_opened = true;
			}
			continue;
		}

		if (bool_flag_ofstream_opened)
		{
			ofstream_file_hmm << str_temp_one_line_file << std::endl;
		}

		if (str_temp_one_line_file.find("//") != -1)
		{
			ofstream_file_hmm.close();
			bool_flag_hmm_begin = true;
			bool_flag_ofstream_opened = false;
		}
	}
	if(error_flag)
	{
		return 0;
	}
	ifstream_file_hmms.close();
	return 1;
}

//HMM file format and exists check
bool HMMTree::hmm_file_exist_format_check(std::string infile_path_and_name)
{
	//check the file exists or not by test to open it
	std::ifstream ifstream_test_fasta_format;
	ifstream_test_fasta_format.open(infile_path_and_name.c_str());
	if (!ifstream_test_fasta_format.is_open())
	{
		std::cout << "File '" << infile_path_and_name << "' doesn't exist!" << std::endl;
		return false;
	}
	ifstream_test_fasta_format.close();

	//put the show numbers of 'HMMER' and end mark '//' into different files
	std::string str_cmd_check_HMM_fromat="";
	str_cmd_check_HMM_fromat = "grep -o 'HMMER' " + infile_path_and_name + " |wc -l > hmm_HMMER_num.txt";
	system_return(system(str_cmd_check_HMM_fromat.c_str()));
	str_cmd_check_HMM_fromat = "grep -o '//' " + infile_path_and_name + " |wc -l > hmm_ends_num.txt";
	system_return(system(str_cmd_check_HMM_fromat.c_str()));

	//read in the  show number of 'HMMER'
	std::ifstream ifstream_hmm_HMMER_num;
	ifstream_hmm_HMMER_num.open("./hmm_HMMER_num.txt");
	if (!ifstream_hmm_HMMER_num.is_open())
	{
		std::cout << "Failed to open 'hmm_HMMER_num.txt' when check the HMM format!" << std::endl;
		return false;
	}
	std::string str_hmm_HMMER_num = "";
	std::getline(ifstream_hmm_HMMER_num, str_hmm_HMMER_num);
	ifstream_hmm_HMMER_num.close();

	//read in the  show number of  end mark '//'
	std::ifstream ifstream_hmm_ends_num;
	ifstream_hmm_ends_num.open("./hmm_ends_num.txt");
	if (!ifstream_hmm_ends_num.is_open())
	{
		std::cout << "Failed to open 'hmm_ends_num.txt' when check the HMM format!" << std::endl;
		return false;
	}

	std::string str_hmm_ends_num = "";
	std::getline(ifstream_hmm_ends_num, str_hmm_ends_num);
	ifstream_hmm_ends_num.close();
	//system_return(system("rm hmm_HMMER_num.txt"));
	if(unlink("./hmm_HMMER_num.txt")){
        return false;
	}
	//system_return(system("rm hmm_ends_num.txt"));
	if(unlink("./hmm_ends_num.txt")){
        return false;
	}
	unsigned int int_hmm_HMMER_num = atoi(str_hmm_HMMER_num.c_str());
	unsigned int int_hmm_ends_num = atoi(str_hmm_ends_num.c_str());

	if ((0 == int_hmm_HMMER_num) || (0 == int_hmm_ends_num) || (int_hmm_HMMER_num != int_hmm_ends_num))
	{
		std::cout << "The content in file '" << infile_path_and_name << "' is not in HMM format!" << std::endl;
		return false;
	}

	std::string str_fasta_ext = infile_path_and_name.substr(infile_path_and_name.find_last_of("."),4);

	if (str_fasta_ext == ".hmm")
	{
		return true;
	}

	std::cout << "'" << infile_path_and_name.c_str() << "' is not in HMM format, the file should be ended by '.hmm' exitemsion!" << std::endl;
	return false;
}



//hmmconvert one hmmer3 hmm model to hmmer2 hmm model
int HMMTree::hmm_hmmconvert_3_to_2_one(std::string hmmer3_name)
{
	std::string str_cmd = "";
	std::string str_pickname_outof_ext = "";
	unsigned int begin = hmmer3_name.find_last_of('/') + 1;
	unsigned int end = hmmer3_name.length();
	str_pickname_outof_ext = hmmer3_name.substr(begin, end - begin);
    if(!bool_HMMER_hmmconvert_in_folder){
        str_cmd = "hmmconvert -2 " + hmmer3_name + " > "+folder_hmmer2 + str_pickname_outof_ext;
    }else{
        str_cmd = "./hmmconvert -2 " + hmmer3_name + " > "+folder_hmmer2 + str_pickname_outof_ext;
    }


	//std::cout << str_cmd << std::endl;
	return system(str_cmd.c_str());
}

//hmmconvert all hmmer3 hmm models to hmmer2 hmm models
int HMMTree::hmm_hmmconvert_3_to_2()
{
	std::vector<std::string> vec_hmm3_2_filenames;
	if(!get_file_names(folder_hmmer3,vec_hmm3_2_filenames,"")){
        return 0;
	}
	if(vec_hmm3_2_filenames.size() < 3 && (!HMMTree::als_phmms_phhms)){
        return 0;
	}
	size_t unint_hmm3_2 = 0;
	bool error_flag = false;
	while (unint_hmm3_2 < vec_hmm3_2_filenames.size())
	{
		int return_flag = 0;
		return_flag = hmm_hmmconvert_3_to_2_one(folder_hmmer3 + vec_hmm3_2_filenames[unint_hmm3_2]);
		if(system_return(return_flag)){
            error_flag = true;
		}

		if (error_flag)
		{
			vec_hmm3_2_filenames.clear();
			break;
		}
		unint_hmm3_2++;
	}
	if (!error_flag)
	{
		vec_hmm3_2_filenames.clear();
	}
	else {
		return 0;
	}
	return 1;
}


/*
//set the pfam clan vector and the hmms data vector and map
int HMMTree::hmm_set_pfamclan_vector_hmmdat_vector_map() {


    //init the vector and the unordered map of hmms of pfam families
	std::ifstream ifstream_hmmdat;
	ifstream_hmmdat.open("./pfam.hmm/Pfam-A.hmm.dat");

	if (!ifstream_hmmdat.is_open()) {
		std::cout << "hmm_set_hmmclan_vector_hmmdat_vector_map(): Failed to open file Pfam-A.hmm.dat" << std::endl;
		return 0;
	}

	std::string str_oneline_stream="";
	while (std::getline(ifstream_hmmdat, str_oneline_stream)) {
		if (str_oneline_stream.find("STOCKHOLM") != -1) {
			HMM_DAT hmm_dat_temp;
			vector_hmm_dat.push_back(hmm_dat_temp);

		}

		if (str_oneline_stream.find("#=GF ID") != -1) {
			vector_hmm_dat.back().ID= str_oneline_stream.substr(10, str_oneline_stream.length()-10);
		}

		if (str_oneline_stream.find("#=GF AC") != -1) {
			vector_hmm_dat.back().AC = str_oneline_stream.substr(10);
			vector_hmm_dat.back().nopoint_AC = str_oneline_stream.substr(10, (str_oneline_stream.find_last_of('.') - 10));
		}
	}

	for (unsigned int uint_i_vec_hmmdat = 0; uint_i_vec_hmmdat < vector_hmm_dat.size(); uint_i_vec_hmmdat++) {
		hmm_dat_nopointAC_datVectorIndex_list_unordered_map[vector_hmm_dat[uint_i_vec_hmmdat].nopoint_AC] = uint_i_vec_hmmdat;
	}

	//init the vector of pfam clans
	std::ifstream ifstream_hmmclan;
	ifstream_hmmclan.open("./pfam.hmm/Pfam-C");

	if (!ifstream_hmmclan.is_open()) {
		std::cout << "hmm_set_hmmclan_vector_hmmdat_vector_map(): Failed to open file Pfam-C" << std::endl;
		return 0;
	}

	str_oneline_stream = "";
	while (std::getline(ifstream_hmmclan, str_oneline_stream)) {
		if (str_oneline_stream.find("STOCKHOLM") != -1) {
			PFAM_C pfam_clan_temp;
			vector_pfam_clans.push_back(pfam_clan_temp);
		}

		if (str_oneline_stream.find("#=GF ID") != -1) {
			vector_pfam_clans.back().ID = str_oneline_stream.substr(10, str_oneline_stream.length() - 10);
		}

		if (str_oneline_stream.find("#=GF AC") != -1) {
			vector_pfam_clans.back().AC = str_oneline_stream.substr(10);
			vector_pfam_clans.back().nopoint_AC = str_oneline_stream.substr(10, str_oneline_stream.find_last_of('.') - 10);
		}

		if (str_oneline_stream.find("#=GF MB") != -1) {
			std::string str_pfam_family_name=str_oneline_stream.substr(10, str_oneline_stream.length() - 11);
			if(hmm_dat_nopointAC_datVectorIndex_list_unordered_map.find(str_pfam_family_name) != hmm_dat_nopointAC_datVectorIndex_list_unordered_map.end()){
                vector_pfam_clans.back().families.push_back(str_pfam_family_name);
			}
		}
	}
	return 1;
}

*/

//copy the hmm files from user dir to hmms
int HMMTree::hmm_copy_hmmfiles(std::string path){

    if (!dir_exist_opendir(path))
	{
        std::cout<<"Failed to open the file directory!"<<std::endl;
		return 0;
	}
    create_files_folder(path,prc_hhsuite,2);

    //test if the dir exists, or create them
    dir_exist_or_create(prc_hhsuite, 2);

    if(dir_noempty_opendir_readir(folder_hmms))
	{
		//system_return(system("rm ./hmms/*.hmm"));
		if(!delete_files(folder_hmms)){
            return 0;
		}
	}


	std::string  path_temp=path;
    if(path_temp[path_temp.length()-1] != '/'){
        path_temp = path_temp+"/";
    }
    copy_files(path_temp,folder_hmms,".hmm");

    return 1;
}

//copy the hmm files from user dir to hmms
int HMMTree::hmm_als_phmms_copy_hmmfiles(std::string path){

    if (!dir_exist_opendir(path))
	{
        std::cout<<"Failed to open the file directory!"<<std::endl;
		return 0;
	}
   // create_files_folder(path,prc_hhsuite,2);

    //test if the dir exists, or create them
    //dir_exist_or_create(prc_hhsuite, 2);
/*
    if(dir_noempty_opendir_readir(folder_hmms))
	{
		//system_return(system("rm ./hmms/*.hmm"));
		if(!delete_files(folder_hmms)){
            return 0;
		}
	}
*/
    //test the hmms folder
    if(!dir_exist_opendir(folder_hmms_from_als))
    {
        system_return(system(("mkdir  -m 777 "+folder_hmms_from_als).c_str()));
    }

    copy_files(folder_hmms,folder_hmms_from_als,".hmm");

	std::string  path_temp=path;
    if(path_temp[path_temp.length()-1] != '/'){
        path_temp = path_temp+"/";
    }
    copy_files(path_temp,folder_hmms,".hmm");

    return 1;
}


//hmmconvert
int HMMTree::hmm_hmmconvert(){
    //usearch
    int hmmconvert_ = HMMER_hmmconvert_exist();
    if( hmmconvert_ != 0){
        if(hmmconvert_ > 0){
            bool_HMMER_hmmconvert_in_folder = true;
        }
    }else{
        std::cout<<"Can not find 'hmmconvert' program in the PATH variable or folder!"<<std::endl;
        exit(1);
    }
    std::cout<<"Convert HMM file version to hmmer2 by hmmconvert command!"<<std::endl;
    //test the hmmer2 folder
    if(!dir_exist_opendir(folder_hmmer2))
    {
        system_return(system(("mkdir  -m 777 "+folder_hmmer2).c_str()));
    }
    //test the hmmer3 folder
    if(!dir_exist_opendir(folder_hmmer3))
    {
        system_return(system(("mkdir  -m 777 "+folder_hmmer3).c_str()));
    }

    if(dir_noempty_opendir_readir(folder_hmmer3))
    {
        //system_return(system("rm ./hmmer3/*.hmm"));
        if(!delete_files(folder_hmmer3)){
            return 0;
        }
    }

    //system_return(system("mv ./hmms/*.hmm ./hmmer3/"));
    if(!mv_files(folder_hmms,folder_hmmer3,"")){
        std::cout<<"alignments_deal(): hmmconvert: mv_files()1: error!"<<std::endl;
        return 0;
    }

    if(!hmm_hmmconvert_3_to_2()){
        std::cout<<"HMMER hmmconvert ERROR!"<<std::endl<<"Please be sure your hmm file version be right!"<<std::endl;
        return 0;
    }

    if(dir_noempty_opendir_readir(folder_hmmer2))
    {
        //system_return(system("mv ./hmmer2/*.hmm ./hmms/"));
        if(!mv_files(folder_hmmer2,folder_hmms,"")){
            std::cout<<"alignments_deal(): hmmconvert: mv_files()2: error!"<<std::endl;
            return 0;
        }
    }

    return 1;
}
