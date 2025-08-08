#include "HMMTree.h"

int HMMTree::Phylip_draw_tree2(){
     kitsch_build_tree((folder_matrixs+"file_dist_matrix_out_phylip.txt").c_str(),(folder_tree_files+"f-m").c_str(), 0);
     kitsch_build_tree((folder_matrixs+"file_dist_matrix_out_phylip.txt").c_str(),(folder_tree_files+"min").c_str(), 1);
     fitch_build_tree((folder_matrixs+"file_dist_matrix_out_phylip.txt").c_str(),(folder_tree_files+"f-m").c_str(), 0);
     fitch_build_tree((folder_matrixs+"file_dist_matrix_out_phylip.txt").c_str(),(folder_tree_files+"min").c_str(), 1);
     neighbor_build_tree((folder_matrixs+"file_dist_matrix_out_phylip.txt").c_str(),(folder_tree_files+"neighbor").c_str());
     upgma_build_tree((folder_matrixs+"file_dist_matrix_out_phylip.txt").c_str(),(folder_tree_files+"upgma").c_str());
     //UPGMA_build_tree("ssssss");
return 1;
}

//function to replace the shortednames to the before-replaced names
void HMMTree::trees_replace_shorted_names(std::string str_treefiles_path){
    if(!dir_exist_opendir(str_treefiles_path)){
        return;
    }

    std::vector<std::string>  vec_treefiles_names;
    if(!get_file_names(str_treefiles_path,vec_treefiles_names,"")){
        std::cout<<"Failed to get the tree files names in shorted names replace step !"<<std::endl;
        return;
	}

	for(unsigned int int_shortednames_i = 0; int_shortednames_i<vec_treefiles_names.size();int_shortednames_i++){
        std::ifstream input_treefile;
        std::string infile_path_name = folder_tree_files+vec_treefiles_names[int_shortednames_i];
        input_treefile.open(infile_path_name.c_str());
        if(!input_treefile.is_open()){
            std::cout<<"Failed to open the file: '"<<infile_path_name<<"' in the shorted names replace step !"<<std::endl;
            return;
        }



        std::string str_treefile_name_after_replacing_shorted_names = vec_treefiles_names[int_shortednames_i]+"_renamed";
        std::ofstream output_treefile;
        std::string outfile_path_name = folder_tree_files +str_treefile_name_after_replacing_shorted_names;
        output_treefile.open(outfile_path_name.c_str());
        if(!output_treefile.is_open()){
            std::cout<<"Failed to open the file: '"<<outfile_path_name<<"' in the shorted names replace step !"<<std::endl;
            return;
        }

        std::string str_one_line="";
        while(std::getline(input_treefile,str_one_line)){
            for(int uint_shorted_names_i=0; uint_shorted_names_i < vec_unshorted_names.size(); uint_shorted_names_i++){
                std::string str_unshorted_name = vec_unshorted_names[uint_shorted_names_i];
                std::string str_shorted_name = vec_shorted_names[uint_shorted_names_i];

                string_replace(str_one_line,str_shorted_name,str_unshorted_name,uint_shorted_names_i,vec_unshorted_names);
            }
            output_treefile<<str_one_line<<std::endl;
        }
        output_treefile.close();
        input_treefile.close();
        unlink(vec_treefiles_names[int_shortednames_i].c_str());
        if(rename(outfile_path_name.c_str(),infile_path_name.c_str())){
            std::cout<<"Failed to rename the file '"<<outfile_path_name <<"' to '"<< infile_path_name<<"' in the step of replacing the shorted names in tree files !"<<std::endl;
            return;
        }
	}
	return;
 }



