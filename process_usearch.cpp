#include "HMMTree.h"

int HMMTree::usearch_rename_result()
{
	std::vector<std::string> vec_usearch_filenames;
	if(!get_file_names(folder_clusters,vec_usearch_filenames,"")){
        return 0;
	}
	if(vec_usearch_filenames.size() < 3){
        return 0;
	}
	size_t unint_usearch = 0;
	std::vector<std::string> vec_usearch_cluster_names;
	while (unint_usearch < vec_usearch_filenames.size())
	{
		vec_usearch_cluster_names.push_back(folder_clusters + vec_usearch_filenames[unint_usearch]);
		unint_usearch++;
	}
	vec_usearch_filenames.clear();

	//std::ifstream ifstream_usearch_cluster_shownum_file;
	std::vector<std::string> vec_usearch_unclustered_seq_names;  //to record the invalid clusters' file names
	std::string str_shownum = "";

	unsigned int valid_seqs_num = 0;                            //the valid sequences number
	unsigned int invalid_seqs_num = 0;
	unsigned int valid_clusters_num = 0;							//the invalid sequences number
																//std::unordered_map<std::string, unsigned int> unordered_map_usearch_cluster_seq_names_notenum;
	for (unsigned int vec_clusternames_i = 0; vec_clusternames_i < vec_usearch_cluster_names.size(); vec_clusternames_i++)
	{
		//put the sequences' number into the file
        std::ifstream if_seqnum;
        if_seqnum.open(vec_usearch_cluster_names[vec_clusternames_i]);
        if(!if_seqnum.is_open()){
            std::cout<<"usearch_rename_result(): Open file error !"<<std::endl;
            return -1;
        }
        size_t uint_seqnum = 0;
        std::string str_one_line="";
        while(std::getline(if_seqnum,str_one_line)){
            if(str_one_line.find(">") != -1){
                uint_seqnum++;
            }
        }
        if_seqnum.close();


		//if the number of the sequences in the cluster is less than 3 , means that this cluster is invalid
		if (uint_seqnum < 3)
		{
			//record the invalid sequences number
			invalid_seqs_num += uint_seqnum;
			//vec_usearch_unclustered_seq_names.push_back(vec_usearch_cluster_names[vec_clusternames_i]);
			continue;
		}

		//record the valid sequences number
		valid_seqs_num += uint_seqnum;
        valid_clusters_num++;
		//rename the valid clusters' files
		//std::string str_temp_cmd_rename = "";
		//str_temp_cmd_rename = "mv " + vec_usearch_cluster_names[vec_clusternames_i] + " " + vec_usearch_cluster_names[vec_clusternames_i] + ".fasta";
		if(rename(vec_usearch_cluster_names[vec_clusternames_i].c_str(),(vec_usearch_cluster_names[vec_clusternames_i] + ".fasta").c_str())){
            std::cout<<"usearch_rename_result(): rename(): error !"<<std::endl;
            return -1;
		}
		//system_return(system(str_temp_cmd_rename.c_str()));
	}
    //system_return(system("rm usearch_cluster_names.txt"));

    /*
    if(unlink("./usearch_cluster_names.txt")){
        std::cout<<"usearch_rename_result(): faild to remove file 'usearch_cluster_names.txt' !"<<std::endl;
        return -1;
    }
    */


	//if the cluster result reach the threshold mvoe the valid clusters to the unlign folder
	if (valid_clusters_num < 3){
        std::cout<<"USEARCH ucluster failed, 'id = "<< CLUSTER_OK_THRESHOLD <<"'"<<std::endl;
		if(dir_noempty_opendir_readir(folder_clusters)){
            //system_return(system("rm ./clusters/*"));
            if(!delete_files(folder_clusters)){
                std::cout<<"usearch_cluster_unalign_seqs(): Clear "<<folder_clusters<<" error!"<<std::endl;
                return -1;
            }
        }
        return 1;
	}
	//clear the clusters folder
    if(dir_noempty_opendir_readir(folder_unalign_seqs)){
        //system_return(system("rm ./unalign_seqs/*"));
        if(!delete_files(folder_unalign_seqs)){
            std::cout<<"usearch_cluster_unalign_seqs(): Clear "<<folder_unalign_seqs<<" error!"<<std::endl;
            return -1;
        }
    }
    //clear the clusters folder
    if(dir_noempty_opendir_readir(folder_invalid_clusters))
    {
        //system_return(system("rm ./invalid_clusters/*"));
        if(!delete_files(folder_invalid_clusters)){
            std::cout<<"usearch_cluster_unalign_seqs(): Clear "<<folder_invalid_clusters<<" error!"<<std::endl;
            return -1;
        }
    }
	std::cout<<"validly clustered sequences: "<<valid_seqs_num<<std::endl;
    std::cout<<"invalidly clustered sequences: "<<valid_seqs_num<<std::endl;
    //system_return(system("mv  ./clusters/*.fasta ./unalign_seqs/"));
    if(!mv_files(folder_clusters,folder_unalign_seqs,".fasta")){
        std::cout<<"usearch_rename_result(): mv files to "<<folder_unalign_seqs  <<" error!"<<std::endl;
        return -1;
    }
    //system_return(system("mv ./clusters/* ./invalid_clusters/"));
    if(dir_noempty_opendir_readir(folder_clusters)){
        if(!mv_files(folder_clusters,folder_invalid_clusters,"")){
            std::cout<<"usearch_rename_result(): mv files to "<<folder_invalid_clusters  <<" error!"<<std::endl;
            return -1;
        }
    }

	return 0;
}

//cluster the unalign sequences by usearch
int HMMTree::usearch_cluster_unalign_seqs(std::string file_path_name, double identity)
{

	//clear the clusters folder
	if(dir_noempty_opendir_readir(folder_clusters))
	{
		//system_return(system("rm ./clusters/*"));
		if(!delete_files(folder_clusters)){
            std::cout<<"usearch_cluster_unalign_seqs(): Clear "<<folder_clusters<<" error!"<<std::endl;
            return 0;
        }
	}

	if(!align_fasta_file_exist_format_check(file_path_name))
	{
		return 0;
	}
	std::string str_double_to_str = double_2_string(identity);

	std::string str_usearch_cmd = "";
	if(bool_userach_in_folder){
        str_usearch_cmd = "./usearch -cluster_fast  "+ file_path_name + " -id " + str_double_to_str + " -clusters "+folder_clusters+"cluster_";
	}else{
        str_usearch_cmd = "usearch -cluster_fast  "+ file_path_name + " -id " + str_double_to_str + " -clusters "+folder_clusters+"cluster_";
	}
    std::cout << str_usearch_cmd << std::endl;
    system_return(system(str_usearch_cmd.c_str()));

    //mark if successfully clustered or not
    if(usearch_rename_result() != 0){
        exit(1);
    }

    //clear the clusters folder
	if(dir_noempty_opendir_readir(folder_clusters))
	{
		//system_return(system("rm ./clusters/*"));
		if(!delete_files(folder_clusters)){
            std::cout<<"usearch_cluster_unalign_seqs(): Clear "<<folder_clusters<<" error!"<<std::endl;
            return 0;
        }
	}
	return 1;
}


//cluster the unalign sequences by usearch
int HMMTree::usearch_cluster(std::string file_path_name, double identity){

    if (identity <= 0.0 || identity >= 1.0)
	{
		std::cout << "Identity error !" << std::endl;
		return 0;
	}

	if (!file_exists_and_empty_check(file_path_name))
	{
		return 0;
	}

    create_files_folder(file_path_name,prc_hhsuite,0);

    //test if the dir exists, or create them
    dir_exist_or_create(prc_hhsuite,0);

    if(!usearch_cluster_unalign_seqs(file_path_name, identity)){
        return 0;
    }
    return 1;
}

