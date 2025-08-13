#include "HMMTree.h"
#include <omp.h>
#include <sys/timeb.h>

int HMMTree::Phylip_draw_tree2(){
     kitsch_build_tree((folder_matrices+"file_dist_matrix_out_phylip.txt").c_str(),(folder_tree_files+"f-m").c_str(), 0, phylo_threads_count);
     kitsch_build_tree((folder_matrices+"file_dist_matrix_out_phylip.txt").c_str(),(folder_tree_files+"min").c_str(), 1, phylo_threads_count);
     fitch_build_tree((folder_matrices+"file_dist_matrix_out_phylip.txt").c_str(),(folder_tree_files+"f-m").c_str(), 0, phylo_threads_count);
     fitch_build_tree((folder_matrices+"file_dist_matrix_out_phylip.txt").c_str(),(folder_tree_files+"min").c_str(), 1, phylo_threads_count);
     neighbor_build_tree((folder_matrices+"file_dist_matrix_out_phylip.txt").c_str(),(folder_tree_files+"neighbor").c_str(), phylo_threads_count);
     upgma_build_tree((folder_matrices+"file_dist_matrix_out_phylip.txt").c_str(),(folder_tree_files+"upgma").c_str(), phylo_threads_count);
     //UPGMA_build_tree("ssssss");
return 1;
}

//function to run selected phylogenetic analyses
void HMMTree::draw_tree_selective(bool run_fitch, bool run_kitsch, bool run_upgma, bool run_nj, bool fm_only, bool min_only) {
    std::string matrix_file = folder_matrices + "file_dist_matrix_out_phylip.txt";
    struct timeb startTime, endTime;
    
    if (run_kitsch) {
        std::cout << "Running Kitsch analysis (contemporary tips method)..." << std::endl;
        ftime(&startTime);
        
        // Determine which variants to run
        bool run_fm = !min_only;    // Run f-m unless min_only is specified
        bool run_min = !fm_only;    // Run min unless fm_only is specified
        
        if (run_fm) {
            kitsch_build_tree(matrix_file.c_str(), (folder_tree_files + "kitsch_f-m").c_str(), 0, phylo_threads_count);
        }
        if (run_min) {
            kitsch_build_tree(matrix_file.c_str(), (folder_tree_files + "kitsch_min").c_str(), 1, phylo_threads_count);
        }
        
        ftime(&endTime);
        std::cout << "Kitsch analysis completed in: " << format_time_duration((endTime.time-startTime.time)*1000 + (endTime.millitm - startTime.millitm)) << std::endl;
    }
    
    if (run_fitch) {
        std::cout << "Running Fitch-Margoliash analysis..." << std::endl;
        ftime(&startTime);
        
        // Determine which variants to run
        bool run_fm = !min_only;    // Run f-m unless min_only is specified
        bool run_min = !fm_only;    // Run min unless fm_only is specified
        
        if (run_fm) {
            fitch_build_tree(matrix_file.c_str(), (folder_tree_files + "fitch_f-m").c_str(), 0, phylo_threads_count);
        }
        if (run_min) {
            fitch_build_tree(matrix_file.c_str(), (folder_tree_files + "fitch_min").c_str(), 1, phylo_threads_count);
        }
        
        ftime(&endTime);
        std::cout << "Fitch-Margoliash analysis completed in: " << format_time_duration((endTime.time-startTime.time)*1000 + (endTime.millitm - startTime.millitm)) << std::endl;
    }
    
    if (run_nj) {
        std::cout << "Running Neighbor-Joining analysis..." << std::endl;
        ftime(&startTime);
        neighbor_build_tree(matrix_file.c_str(), (folder_tree_files + "neighbor").c_str(), phylo_threads_count);
        ftime(&endTime);
        std::cout << "Neighbor-Joining analysis completed in: " << format_time_duration((endTime.time-startTime.time)*1000 + (endTime.millitm - startTime.millitm)) << std::endl;
    }
    
    if (run_upgma) {
        std::cout << "Running UPGMA analysis..." << std::endl;
        ftime(&startTime);
        upgma_build_tree(matrix_file.c_str(), (folder_tree_files + "upgma").c_str(), phylo_threads_count);
        ftime(&endTime);
        std::cout << "UPGMA analysis completed in: " << format_time_duration((endTime.time-startTime.time)*1000 + (endTime.millitm - startTime.millitm)) << std::endl;
    }
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



