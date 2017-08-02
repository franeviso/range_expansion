////////////////// RANGE EXPANSION MODEL //////////////////////////////
////////////////// @ FRANCISCO ENCINAS VISO 2017 //////////////////////

/* Changes to make to previous model:
 * 1) Selection acting on parent's selection (choosing parents with high fitness)
 * 2) Modify output files to facilitate interpretation of results and parsing
 * 3) No selfing rates differences between SC/SI heterozygotes and SC/SC homozygotes (only one conditional for selfing)
 * 4) Be explicit about recombination rates.
 * 5) Add warning signs (STOP or break) to calculation of inbreeding depression
 * 6) Change carrying capacity according to average population fitness sensu Peischl et al 2015
 * 7) Try fixed selection coefficient and dominance h = 0
 * 8) Calculate polllen/seed discounting
 * 9) Potentially make lattice smaller to get faster results (5 X 50)
 * 10) Heterozygotes and homozygotes with SC allele treated equally sensu Gervais et al 2014
 * 11) Consider adding mutation rate of new S alleles
 * 12) Perhaps add STOP sign if no compatible mates are found after X iterations
 */


#include <iostream>
#include <cstdlib>
#include <vector>
#include <time.h>
#include <math.h>  
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <fstream>
#include <iomanip> 
#include <algorithm>
#include <boost/format.hpp>
#include <time.h>
#include <sstream>




using namespace std;
using boost::format;


struct indiv{ //properties of individual plants
        vector<int> genome;
        int patch_x;
        int patch_y;                       // Patch (or site) in the landscape
        int chosen;                   //for when an adult survives to the next year; to make sure the same plant doesn't survive twice. 10(chosen, -10(not chosen) / Also to stablish which individuals survive in a local extinction
        double fitness;
        int id;
        int S_allele1(void){return genome[(genome.size()/4)-1];}
        int S_allele2(void){return genome[3*(genome.size()/4)-1];}
};



int random_gen_selected(gsl_rng *r){
	int x;
	if(gsl_rng_uniform(r)>0.5) x = 21;
	else  x = 23;
	return x;
}

int random_gen_neutral(gsl_rng *r, int neutral_alleles){
	int x = gsl_rng_uniform_int(r, neutral_alleles) + 30000;
	return x;
}


int random_S_locus(gsl_rng *r, int S_alleles){
	int x = gsl_rng_uniform_int(r, S_alleles) + 100;
	return x;
}



void mutate_poisson_neutral(gsl_rng *r,vector<int> &adn, double lambda, vector<int> &Ntype) {
	int x = gsl_ran_poisson(r, lambda); // number of selected sites
	int selected, diffs;

	if(x>0){
		for(int i=0; i<x; i++){
		    selected = gsl_rng_uniform_int(r,adn.size()-1); //Selected sites
			if(adn[selected] >10000){ // Neutral loci - step-wise mutation model
			    
			    //cout << " JACKPOT !!! " << endl; 
			    //cout << "Selected " << selected << " Site " << adn[selected] << endl;
			    int s=0;
				if(gsl_rng_uniform(r) > 0.5) s=1;
				else s=-1;
				int temp = adn[selected];
				adn[selected] = temp + s;
				//Update neutral alleles list
				diffs=0;
				for(vector<int>::iterator iter=Ntype.begin(); iter!=Ntype.end(); iter++){
				     if(adn[selected] == *iter)  diffs++;
				}
				if(diffs == 0)     Ntype.push_back(adn[selected]);
				
			}

		}    
	}
}	


void mutate_poisson_selected(gsl_rng *r,vector<int> &adn, double lambda) {
	int x = gsl_ran_poisson(r, lambda); // number of selected sites
	//cout << "Number of mutations " << x << " lambda "  << lambda << endl;
	int selected;
	double muback = 0.01;
	if(x>0){
		for(int i=0; i<x; i++){
		    selected = gsl_rng_uniform_int(r,adn.size()-1); //Selected sites
			if(adn[selected]==21) adn[selected] = 23;
			else if(adn[selected]==23){
				if(muback > gsl_rng_uniform(r))  adn[selected] = 21;
			}

		}    
	}
}	



void mutate_S_locus(gsl_rng *r,vector<int> &adn, double mu_slocus, vector<int> &Stype, int flag) {
	int L = adn.size(), diffs;
	int S_allele = adn[(L/2)-1];
	if(mu_slocus > gsl_rng_uniform(r) && S_allele>0){
		if(flag == 1){
		     adn[(L/2)-1] = -100;//S_allele;
		     //cout << " ------------------------------------------------------------------------------------> MUTATION S LOCUS " << adn[(L/2)-1] << endl;
		}
		//Update S alleles list
		diffs=0;
		for(vector<int>::iterator iter=Stype.begin(); iter!=Stype.end(); iter++){
			if(adn[(L/2)-1] == *iter)  diffs++;
		}
		if(diffs == 0){
			Stype.push_back(adn[(L/2)-1]);
			//cout << "Mutation S locus " << adn[(L/2)-1] << endl;
		}
	}
}	


void mutate_S_locus2(gsl_rng *r,vector<int> &adn, double mu_slocus, vector<int> &Stype, int flag) {
	int L = adn.size(), diffs, S_alleles=1000;
	int S_allele = adn[(L/2)-1];
	if(mu_slocus > gsl_rng_uniform(r) && S_allele>0){
		if(flag == 1){
		     if(0.1 > gsl_rng_uniform(r))  adn[(L/2)-1] = -S_allele;
		     else adn[(L/2)-1] = random_S_locus(r, S_alleles);
		}else adn[(L/2)-1] = random_S_locus(r, S_alleles);
		//cout << " ------------------------------------------------------------------------------------> MUTATION S LOCUS " << adn[(L/2)-1] << endl;
		//Update S alleles list
		diffs=0;
		for(vector<int>::iterator iter=Stype.begin(); iter!=Stype.end(); iter++){
			if(adn[(L/2)-1] == *iter)  diffs++;
		}
		if(diffs == 0){
			Stype.push_back(adn[(L/2)-1]);
			//cout << "Mutation S locus " << adn[(L/2)-1] << endl;
		}
	}
}	


void mutations(gsl_rng *r, vector<int> &adn, double lambda_neutral, double lambda_sel, double mu_slocus, vector<int> &Ntype, vector<int> &Stype, int flag){
	mutate_poisson_neutral(r, adn, lambda_neutral, Ntype);
	mutate_poisson_selected(r, adn, lambda_sel);
	mutate_S_locus(r,adn, mu_slocus, Stype, flag);
}

vector<int> genome_fun(gsl_rng *r,int selgenes,int neutral_alleles, int S_alleles){
	vector<int> genome_vec;
	vector<int> genome_vec2;
	int L = selgenes;

	for(int gen=0; gen<L; gen++){
		
		if(gen%10 == 0)
		     genome_vec.push_back(random_gen_neutral(r,neutral_alleles));
		else if(gen==(L/2)-1)
		     genome_vec.push_back(random_S_locus(r,S_alleles));
		else
		     genome_vec.push_back(21); //random_gen_selected(r)
		
	}
	
	for(int gen=0; gen<L; gen++){
		
		if(gen%10 == 0)
		     genome_vec2.push_back(random_gen_neutral(r,neutral_alleles));
		else if(gen==(L/2)-1)
		     genome_vec2.push_back(random_S_locus(r,S_alleles));
		else
		     genome_vec2.push_back(21); //random_gen_selected(r)
	}

	genome_vec.insert(genome_vec.end(), genome_vec2.begin(), genome_vec2.end());
	return genome_vec;
}

void show_vector(vector<int> &vec){
	cout << "Vector of length: " << vec.size() << endl;
	cout << "Vector contains: ";
    for (vector<int>::iterator it=vec.begin(); it!=vec.end(); ++it)
         cout << ' ' << *it;
    cout << endl;
}


void show_sc(vector<indiv> *vec){
	int sc=0,n=0;
	//cout << "Self-compatibles: ";
    for (vector<indiv>::iterator it=vec->begin(); it!=vec->end(); it++){
         if(it->S_allele1() < 0){
			  sc++;
			 // cout << "S allele 1 " << it->S_allele1() << endl;
			  //cout << "ind " << n << endl;
		 }
         if(it->S_allele2() < 0){
			  sc++;
			   //cout << "S allele 2 " << it->S_allele2() << endl;
		 }
		 n++;    
	 }
	
    //cout << sc << endl;
}

void number_genes(vector<int> &vec){
	int sel=0, neutral=0, S=0;
	//cout << "Genes: ";
    for (vector<int>::iterator it=vec.begin(); it!=vec.end(); ++it){
		if(*it==21 || *it==23)   sel++;
		else if(*it>=30000)         neutral++;
		else if(*it>=100 && *it<500)         S++;
	}
    //cout <<  "Sel " << sel << " N " << neutral << " S " << S << endl;
}


void refresh_chosen(vector<indiv> *adults){
	for(vector<indiv>::iterator iter=adults->begin(); iter!=adults->end(); iter++){
		iter->chosen=-10;
	}
}

/// New gamete production and mating function
bool sorted(int i, int j){ return (i<j); }


vector<int> choose_dna(gsl_rng *r, vector<int> adn1, vector<int> adn2){
	vector<int> adn;
	if(gsl_rng_uniform(r)>0.5){
		adn = adn1;
	}else{
		adn = adn2;
	}
	return adn;
}



vector<int> create_haplotype(gsl_rng *r, vector<int> genome, double lambda, int genes){
	int i=0,x;
	vector<int> haplotype;
	vector<int> adn;
	vector<int> selected;
	vector<int> pool;
	/// First chromosome:
	vector<int>::const_iterator first_chrom1 = genome.begin();
    vector<int>::const_iterator last_chrom1 = genome.begin() + genome.size()/2;
    vector<int> adn1(first_chrom1,last_chrom1); // First chromosome - pair one
    /// Second chromosome:
    vector<int>::const_iterator first_chrom2 = genome.begin() + genome.size()/2;
    vector<int>::const_iterator last_chrom2 = genome.end();
    vector<int> adn2(first_chrom2,last_chrom2); // Second chromosome - pair one
	// Pick cross-over points:
    // number of selected crossover points
	do{
		x = gsl_ran_poisson(r, lambda); 
	}while(x>genes);
	if(x>0){
		 // Fisher-Yates algorithm
		 for(int j=0; j<genes; j++)
	          pool.push_back(j);
	     random_shuffle(pool.begin(), pool.end()); // Create shuffled array 
         selected.insert(selected.begin(), pool.begin(), pool.begin()+x); // Select x sites
	     sort(selected.begin(), selected.end(), sorted);
         for(vector<int>::iterator iter=selected.begin(); iter!=selected.end(); iter++){ // Recombination process
		        if(i%2==0) adn = adn1;
		        else       adn = adn2;
		        if(i==0){ 
					 haplotype.insert(haplotype.end(), adn.begin(), adn.begin() + *iter);
				}
				else if(i>0){
		             haplotype.insert(haplotype.end(), adn.begin() + *(iter-1), adn.begin() + *iter);
				}
		        i++;
		  }
		  if(i%2==0) adn = adn1;
		  else       adn = adn2;
		  haplotype.insert(haplotype.end(), adn.begin() + selected[x-1], adn.end());
	}else{
		  haplotype = choose_dna(r,adn1,adn2);
	}
	if(haplotype.size() != 500)  show_vector(haplotype);
	return haplotype;
}



double calculatefitness(gsl_rng *r,vector<int> &gen, int genes, double mu, double h){
	int i;
	int L=genes, allele1, allele2;
	double x1, x2,x3,s = mu;
	vector<double> fitness;
	int count=0, countna=0, homoB=0, homob=0, hetero=0;
	for(i=0; i<genes; i++){
        allele1 = gen[i];
        allele2 = gen[i+L];
  
		if (allele1==21 || allele1 ==23 || allele2==21 || allele2==23){
			//cout << "Allele 1 " << allele1  <<  " Allele 2 " << allele2 << endl;
			if(allele1+allele2 == 42){
			    x1 = 1.0; 
			    homoB++;
			    fitness.push_back(x1);
			    //cout << "Non-del " << fitness[countna] << " " << x1 << endl;
		    }else if(allele2+allele1==46){
			    x2 = 1.0 - s;//*(mat_b + (countna)*3 + j);
			    homob++;
			    fitness.push_back(x2);
			    //cout << "Deleterious " << fitness[countna] << " " << allele1 << " " << allele2 << endl;
		    }else if(allele1+allele2==44){
			    x3 =  1.0 - h*s;//*(mat_b + (countna)*3 + j);
			    hetero++;
			    fitness.push_back(x3);
			    //cout << "Hetero "  << fitness[countna] << " " << allele1 << " " << allele2 << endl;
			}
			countna++;
		}
		count++;
			
	}
	//cout << "Fitness size vector " << fitness.size() << endl;
	// Product of fitness effects:
	double fitness_product = fitness[0];
	for(unsigned int j=1;j<fitness.size();j++){
		fitness_product *= fitness[j];
	}
	//cout << "Homo non-del: " << homoB  << " Hetero: "  << hetero << " Homo del: " << homob << " Fitness: " << fitness_product << endl;
	return fitness_product;
}




// Population initialization for fully populated space (metapopulation) or range expansion case 
void initializeplants(gsl_rng *r, vector<indiv> *p, int selgenes, int nalleles, int S_alleles, int subpop, int patchesx, int patchesy, double mu, double dominance, int **space)
{ 
	    int ii=0;
	    for(int i=0; i<patchesx; i++){
			for(int j=0; j<patchesy; j++){
				for(int inds=0; inds<subpop; inds++){
					(*p)[inds + ii].genome = genome_fun(r, selgenes, nalleles, S_alleles); // ngenes ---> number of neutral alleles   
                    (*p)[inds + ii].patch_x = i;
                    (*p)[inds + ii].patch_y = j;
                    (*p)[inds + ii].chosen = -10;
                    (*p)[inds + ii].fitness = calculatefitness(r, (*p)[inds+ii].genome , (*p)[inds+ii].genome.size()/2, mu, dominance); 
                    //(*p)[inds + ii].S_allele = s_alleles_fun( (*p)[inds+ii].genome );
                    //cout << "S alleles: " << (*p)[inds + ii].S_allele1() << " " << (*p)[inds + ii].S_allele2() << endl;
				}
				space[i][j]=1;
				ii=ii+subpop;
			}
		}
							
}

void show_vector_pointer(vector<indiv> *vec){

    for (vector<indiv>::iterator it=vec->begin(); it!=vec->end(); ++it){
	    cout << "Vector contains: ";
		for(vector<int>::iterator gen=it->genome.begin(); gen!=it->genome.end(); gen++)
            cout << " " << *gen;
        cout << endl;
    }
}

void glass(vector<int> &vec){

    for (vector<int>::iterator it=vec.begin(); it!=vec.end(); ++it){
		if(*it == 0){
			show_vector(vec);
			cout << "FLAG 1" << endl;
			exit(EXIT_FAILURE);
		}else if(*it < -100){
			show_vector(vec);
			cout << "FLAG 2" << endl;
			exit(EXIT_FAILURE);
		}else if(*it > 100000){
			show_vector(vec);
			cout << "FLAG 3" << endl;
			exit(EXIT_FAILURE);
		}
	}

}


vector<int> fertilization(vector<int> *ovule, vector<int> *pollen){
	vector<int> seed;
    seed.reserve(ovule->size() + pollen->size());
	seed.insert(seed.end(), ovule->begin(), ovule->end());
	seed.insert(seed.end(), pollen->begin(), pollen->end());
    return seed;
}


vector<int> fertilization2(vector<int> ovule, vector<int> pollen){
	vector<int> seed;
    seed.reserve(ovule.size() + pollen.size());
	seed.insert(seed.end(), ovule.begin(), ovule.end());
	seed.insert(seed.end(), pollen.begin(), pollen.end());
    return seed;
}


void build_subpop_vector_pointer2D(int patchx, int patchy, vector<indiv> *vectors, vector<indiv> *subpop){
	int popsize = (*vectors).size();
	for(int i=0; i<popsize; i++){
			if((*vectors)[i].patch_x == patchx && (*vectors)[i].patch_y==patchy){
				subpop->push_back((*vectors)[i]);
			}
	}
}

void build_subpop_vector_pointer2D_indexes(int patchx, int patchy, vector<indiv> *vectors, vector<int> *subpop){
	int popsize = (*vectors).size();
	for(int i=0; i<popsize; i++){
			if((*vectors)[i].patch_x == patchx && (*vectors)[i].patch_y==patchy){
				subpop->push_back(i);
			}
	}
}




vector<int> choose_migrants2(gsl_rng *r, vector<int> *subpop_indexes, int migrants){ 
	vector<int> chosenmigrants;
    random_shuffle(subpop_indexes->begin(), subpop_indexes->end()); // Create shuffled array 
    chosenmigrants.insert(chosenmigrants.begin(), subpop_indexes->begin(), subpop_indexes->begin()+migrants); // Select migrants
    return chosenmigrants;
}


void choose_migrants(gsl_rng *r, vector<int> *subpop_indexes, vector<int> *chosenmigrants, int migrants){ 
    random_shuffle(subpop_indexes->begin(), subpop_indexes->end()); // Create shuffled array 
    chosenmigrants->insert(chosenmigrants->begin(), subpop_indexes->begin(), subpop_indexes->begin()+migrants); // Select migrants
}

void find_immigrant(int patchx, int patchy, vector<indiv> *seeds, int *index){
	int i=0;
	for(vector<indiv>::iterator it=seeds->begin(); it!=seeds->end(); it++){
		if(it->patch_x == patchx && it->patch_y == patchy && it->chosen ==-10){
			  *index = i;
			  break;
		}
		i++;
	}
}

void build_subpop_indexes(int patchx, int patchy, vector<indiv> *seeds, vector<int> *indexes){
	int i=0;
	for(vector<indiv>::iterator it=seeds->begin(); it!=seeds->end(); it++){
		if(it->patch_x == patchx && it->patch_y == patchy && it->chosen==-10){
			indexes->push_back(i);
		}
		i++;
	}
}


// This function picks a neighbor patch to and provide the coordinates of the patch selected to make list colonizers
int * new_patches_neighbor(gsl_rng *r, int patchx, int patchy, int patchesx, int patchesy, int **spatial_occupancy){
	int x,y,dirx,diry;
	vector<int> list;
	static int patchvec[2];
	//int *patch_vector_pt = patchvec;
	for (int i=-1; i<2; ++i) list.push_back(i); // -2 -1 0 1 2 -> 24 neighbor patches 
	//cout << "central Patch " << patchx << " " << patchy << endl;
	do{ // Busca solo dentro del lattice
	    random_shuffle(list.begin(), list.end());
        dirx = list[0];
        random_shuffle(list.begin(), list.end());
        diry = list[0];
		x = patchx + dirx;
		y = patchy + diry;
		//cout << "Patch: " << x << " Patch: " << y << endl;
	}while(x >= patchesx || y >= patchesy || x < 0 || y < 0 || spatial_occupancy[x][y] ==0 || (dirx ==0 && diry ==0) );
	
	
	//patch_vector_pt[0] = x;
	//patch_vector_pt[1] = y;
	patchvec[0] = x;
	patchvec[1] = y;
	//cout << "Selected patch " << patchvec[0] << " " << patchvec[1] << endl;
	return patchvec;
}

// Choose colonists at random from the whole population (migrant pool)
void choose_colonists_propagule(gsl_rng *r, vector<indiv> *population, vector<int> *subpop ,vector<int> *chosencol, int colonists){
	int k=0;
	random_shuffle(subpop->begin(),subpop->end());
	
	while(k < colonists){
		if( (*population)[(*subpop)[k]].chosen== -10){
			chosencol->push_back(k);
		}
		k++;
	}
}


// CHANGE IT TO STEPPING-STONE MIGRATION
void migration_stepping_stone(gsl_rng *r, vector<indiv> *seeds, double mig, int patchesx, int patchesy, int **spatial_occupancy){
	// Pick subpop*migrants
	int index, migrants, migrantnumber=0, subpopsize;
	int *pt;
	int *pt_index = &index;
	for(int i=0; i<patchesx; i++){
		for(int j=0; j<patchesy; j++){
			if(spatial_occupancy[i][j] == 1){
				vector<int> *chosen_migrants = new vector<int>;
				vector<int> *subpop_indexes = new vector<int>;
				vector<indiv> *subpop = new vector<indiv>;
				build_subpop_vector_pointer2D(i,j,seeds,subpop);
				subpopsize = subpop->size();
				do{
					migrants = gsl_ran_poisson(r,mig*subpop->size());
				}while(migrants > subpopsize);
				//Leaving individuals from patch (i,j)
				build_subpop_indexes(i,j,seeds,subpop_indexes);	
				choose_migrants(r, subpop_indexes, chosen_migrants, migrants);
			    for(vector<int>::iterator iter=chosen_migrants->begin(); iter!=chosen_migrants->end(); iter++){
					pt = new_patches_neighbor(r, i, j, patchesx, patchesy, spatial_occupancy);
					find_immigrant(*pt,*(pt+1),seeds,pt_index);
					(*seeds)[*iter].patch_x = *pt;
				    (*seeds)[*iter].patch_y = *(pt+1);
				    (*seeds)[*iter].chosen=10;
                    (*seeds)[*pt_index].patch_x = i;
                    (*seeds)[*pt_index].patch_y = j;
                    (*seeds)[*pt_index].chosen = 10;
					migrantnumber++;
					
					//cout << "Individual emigrant: " << *iter << " from patch " << i << " " << j << " to patch " << *pt << "  " << *(pt+1) << endl;
					//cout << "Individual immigrant: " << *pt_index << " from patch " << *pt << " " << *(pt+1) << " to patch " << i << "  " << j << endl;
				}
				// Choose individuals without replacement to come to patch (i,j)
				delete subpop_indexes;
				delete chosen_migrants;
				delete subpop;
			}
		}
	}
	//cout << "Number of migrants: " << migrantnumber << endl;
    //cout << "Popsize after migration: " << seeds->size() << endl;		
		
}



void neighbor_list(vector<vector<int> > *pt_neilist1, int patchx, int patchy, int L){
	int neighbors [] = {-1,0,1};
	vector<int> list(2);
	if(patchx>0 && patchx<L && patchy>0 && patchy<L){
		for(int i=0; i<3; i++){
			for(int j=0; j<3; j++){
			  if(patchx + neighbors[i] == patchx &&  patchy + neighbors[j] ==patchy){
				  continue;
			  }else{
					list[0] = patchx + neighbors[i]; list[1] = patchy + neighbors[j];
					(*pt_neilist1).push_back(list);
			  }
			}
		}		
	}else if(patchx==0 && patchy >0 && patchy<L){
		for(int i=1; i<3; i++){
			for(int j=0; j<3; j++){
			  if(patchx + neighbors[i] == patchx &&  patchy + neighbors[j] ==patchy){
				  continue;
			  }else{
					list[0] = patchx + neighbors[i]; list[1] = patchy + neighbors[j];
					(*pt_neilist1).push_back(list);
			  }
			}
		}		
	}else if(patchx==L && patchy >0 && patchy<L){
		for(int i=0; i<2; i++){
			for(int j=0; j<3; j++){
			  if(patchx + neighbors[i] == patchx &&  patchy + neighbors[j] ==patchy){
				  continue;
			  }else{
					list[0] = patchx + neighbors[i]; list[1] = patchy + neighbors[j];
					(*pt_neilist1).push_back(list);
			  }
			}
		}		
	}else if(patchx>0 && patchx<L && patchy==0){
		for(int i=0; i<3; i++){
			for(int j=1; j<3; j++){
			  if(patchx + neighbors[i] == patchx &&  patchy + neighbors[j] ==patchy){
				  continue;
			  }else{
					list[0] = patchx + neighbors[i]; list[1] = patchy + neighbors[j];
					(*pt_neilist1).push_back(list);
			  }
			}
		}		
	}else if(patchx>0 && patchx<L && patchy==0){
		for(int i=0; i<3; i++){
			for(int j=1; j<3; j++){
			  if(patchx + neighbors[i] == patchx &&  patchy + neighbors[j] ==patchy){
				  continue;
			  }else{
					list[0] = patchx + neighbors[i]; list[1] = patchy + neighbors[j];
					(*pt_neilist1).push_back(list);
			  }
			}
		}		
	}else if(patchx>0 && patchx<L && patchy==L){
		for(int i=0; i<3; i++){
			for(int j=0; j<2; j++){
			  if(patchx + neighbors[i] == patchx &&  patchy + neighbors[j] ==patchy){
				  continue;
			  }else{
					list[0] = patchx + neighbors[i]; list[1] = patchy + neighbors[j];
					(*pt_neilist1).push_back(list);
			  }
			}
		}		
	}else if(patchx==0 && patchy==0){
		for(int i=1; i<3; i++){
			for(int j=1; j<3; j++){
			  if(patchx + neighbors[i] == patchx &&  patchy + neighbors[j] ==patchy){
				  continue;
			  }else{
					list[0] = patchx + neighbors[i]; list[1] = patchy + neighbors[j];
					(*pt_neilist1).push_back(list);
			  }
			}
		}		
	}else if(patchx==L && patchy==0){
		for(int i=0; i<2; i++){
			for(int j=1; j<3; j++){
			  if(patchx + neighbors[i] == patchx &&  patchy + neighbors[j] ==patchy){
				  continue;
			  }else{
					list[0] = patchx + neighbors[i]; list[1] = patchy + neighbors[j];
					(*pt_neilist1).push_back(list);
			  }
			}
		}		
	}else if(patchx==0 && patchy==L){
		for(int i=1; i<3; i++){
			for(int j=0; j<2; j++){
			  if(patchx + neighbors[i] == patchx &&  patchy + neighbors[j] ==patchy){
				  continue;
			  }else{
					list[0] = patchx + neighbors[i]; list[1] = patchy + neighbors[j];
					(*pt_neilist1).push_back(list);
			  }
			}
		}		
	}else if(patchx==L && patchy==L){
		for(int i=0; i<2; i++){
			for(int j=0; j<2; j++){
			  if(patchx + neighbors[i] == patchx &&  patchy + neighbors[j] ==patchy){
				  continue;
			  }else{
					list[0] = patchx + neighbors[i]; list[1] = patchy + neighbors[j];
					(*pt_neilist1).push_back(list);
			  }
			}
		}		
	}
}

void list_migrants(gsl_rng *r, vector<vector<int> > *pt_neilist, vector<int> *chosenmigrantsindex, vector<indiv> *population, double mig ){
	int subpopsize, migrants;
	for(vector<vector<int> >::iterator iter= pt_neilist->begin(); iter!=pt_neilist->end(); iter++){
		vector<int> *subpop_indexes = new vector<int>;
		vector<int> *chosenmigrants = new vector<int>;
		build_subpop_vector_pointer2D_indexes((*iter)[0], (*iter)[1], population, subpop_indexes);
		subpopsize = subpop_indexes->size();
		migrants = gsl_ran_poisson(r, mig*subpopsize); // Number of migrants depending on subpop size
		//double product = mig*subpopsize;
		//cout << " Migrants: " << migrants << " Subpop size " << subpopsize << " migration " << mig << " product " << product << endl;
		*chosenmigrants =  choose_migrants2(r, subpop_indexes, migrants);
		chosenmigrantsindex->insert(chosenmigrantsindex->begin(), chosenmigrants->begin(), chosenmigrants->end());
		delete subpop_indexes;
		delete chosenmigrants;
	}
}

bool rankfunc (indiv i, indiv j) { return (i.fitness>j.fitness); }

double maxfit_indiv(vector<indiv> *subpop){
	sort (subpop->begin(), subpop->end(), rankfunc);
	return (*subpop)[0].fitness;
}

vector<int> find_mate_GSI(gsl_rng *r,vector<indiv> *subpop, vector<int> &ovule_S_allele, double rec, double lambda_neutral, double lambda_sel, double mu_slocus, int genes, vector<int> &Ntype, vector<int> &Stype, int flag){
	int k, L, fail = 0;
	double maxfit;
	vector<int> pollen_haplotype, temp;
	int pop_size = subpop->size();
	k = gsl_rng_uniform_int(r,pop_size);
	pollen_haplotype = create_haplotype(r, (*subpop)[k].genome, rec, genes);
	L = pollen_haplotype.size(); 
	//cout << " First Ovule S: " << ovule_S_allele[0] << " " << ovule_S_allele[1] << "  Pollen S: " << pollen_haplotype[(L/2)-1]  << endl;
	maxfit = maxfit_indiv(subpop);
	for(int counter = 0; counter < pop_size; counter++){
     /// Functional S allele
			//cout << "Ovule S: " << ovule_S_allele[0] << " " << ovule_S_allele[1] << "  Pollen S: " << pollen_haplotype[(L/2)-1] << endl;
		    if(pollen_haplotype[(L/2)-1] != ovule_S_allele[0] &&  pollen_haplotype[(L/2)-1] != ovule_S_allele[1] && gsl_rng_uniform(r) < (*subpop)[k].fitness/maxfit){
				   temp = pollen_haplotype;
				   fail = 1;
			       break;
		    }
		k = gsl_rng_uniform_int(r,pop_size);
		pollen_haplotype = create_haplotype(r, (*subpop)[k].genome, rec, genes);
	}
	
	if(fail == 0){
		temp.push_back(-1);
		//cout << "NO compatible MATES !!" << endl;
	}else mutations(r, temp, lambda_neutral, lambda_sel, mu_slocus, Ntype, Stype,flag);
	
	return temp;
}


vector<int> find_mate_GSI_SC(gsl_rng *r,vector<indiv> *subpop, vector<int> &ovule_S_allele, double rec, double lambda_neutral, double lambda_sel, double mu_slocus, int genes, vector<int> &Ntype, vector<int> &Stype, int flag){
	int counter=0,k, L;
	double maxfit;
	vector<int> pollen_haplotype, temp;
	int pop_size = subpop->size();
	k = gsl_rng_uniform_int(r,pop_size);
	pollen_haplotype = create_haplotype(r, (*subpop)[k].genome, rec, genes);
	L = pollen_haplotype.size(); 
	//cout << " First Ovule S: " << ovule_S_allele[0] << " " << ovule_S_allele[1] << "  Pollen S: " << pollen_haplotype[(L/2)-1]  << endl;
	maxfit = maxfit_indiv(subpop);	
	while(counter<pop_size){
     /// Functional S allele
		//cout << "Ovule S: " << ovule_S_allele[0] << " " << ovule_S_allele[1] << "  Pollen S: " << pollen_haplotype[(L/2)-1] << endl;
		if(pollen_haplotype[(L/2)-1] > 0){	
		    if(pollen_haplotype[(L/2)-1] != ovule_S_allele[0] &&  pollen_haplotype[(L/2)-1] != ovule_S_allele[1] && gsl_rng_uniform(r) < (*subpop)[k].fitness/maxfit){
				   temp = pollen_haplotype;
			       break;
		    }
        }else{ // SC pollen haplotype - functional S allele of ovule has to be different than pollen haplotype
			if(ovule_S_allele[0] > 0 &&  pollen_haplotype[(L/2)-1] != ovule_S_allele[0] && gsl_rng_uniform(r) < (*subpop)[k].fitness/maxfit){ 
				   temp = pollen_haplotype;
			       break;
		    }else if(ovule_S_allele[1] > 0 &&  pollen_haplotype[(L/2)-1] != ovule_S_allele[1] && gsl_rng_uniform(r) < (*subpop)[k].fitness/maxfit){
				   temp = pollen_haplotype;
			       break;
			}
		}
		k = gsl_rng_uniform_int(r,pop_size);
		pollen_haplotype = create_haplotype(r, (*subpop)[k].genome, rec, genes);
		counter++;
	}
	
	//cout << "Counter " << counter << endl;
	if(counter >= pop_size)   temp.push_back(-1);
	else mutations(r, temp, lambda_neutral, lambda_sel, mu_slocus, Ntype, Stype,flag);
	
	return temp;
}

vector<int> find_pollen(gsl_rng *r,vector<indiv> *subpop, vector<int> &ovule_S_allele, double rec, double lambda_neutral, double lambda_sel, double mu_slocus, int genes, vector<int> &Ntype, vector<int> &Stype, int flag){
	int chosen_father;
	double maxfit;
	vector<int> pollen_haplotype;
	int pop_size = subpop->size();
    maxfit = maxfit_indiv(subpop);
	chosen_father = gsl_rng_uniform_int(r, pop_size);
	do{
	   chosen_father = gsl_rng_uniform_int(r, pop_size);
    }while( gsl_rng_uniform(r) > (*subpop)[chosen_father].fitness/maxfit);
	pollen_haplotype = create_haplotype(r, (*subpop)[chosen_father].genome, rec, genes);
	mutations(r, pollen_haplotype, lambda_neutral, lambda_sel, mu_slocus, Ntype, Stype, flag);
	return pollen_haplotype;
}




// Check which patches are next to occupied patches and make list

void list_expansion(int patchesx, int patchesy, int **spatial_occupancy, vector< vector<int> > *listexp){
	vector<int> coords(2);
	for(int i=0; i<patchesx; i++){
		for(int j=0; j<patchesy-1; j++){
			if(spatial_occupancy[i][j] == 1 && spatial_occupancy[i][j+1] == 0){
				coords[0] = i;
				coords[1] = j+1;
				//cout << coords[0] << "  " << coords[1] << endl;
				listexp->push_back(coords);		
			}
		}
	}
}


int subpopsize(int patchx, int patchy, vector<indiv> *pop){
	int inds=0;
	for(vector<indiv>::iterator iter=pop->begin(); iter!=pop->end(); iter++){
		if(iter->patch_x == patchx && iter->patch_y == patchy)   inds++;
	}
	return inds;
}





bool rank (indiv i, indiv j) { return (i.fitness>j.fitness); }

void colonization(gsl_rng *r, vector<indiv> *seeds, vector<vector<int> > *listexp, int N, int genes, double mig, int patchesx, int patchesy, int **spatial_occupancy, int flag, int expa){
    int i,j,l;
    list_expansion(patchesx, patchesy, spatial_occupancy, listexp);
    for(vector< vector<int> >::iterator patch = listexp->begin(); patch!=listexp->end(); patch++){
		for(l=0;l<2;l++){
			i = (*patch)[0]; 
			j = (*patch)[1];
		}
		int colonists=0;
		vector<vector<int> > *pt_neilist = new vector<vector<int> >;
		vector<int> *pt_chosen = new vector<int>;  
		// Loop through list of colonists:
	    neighbor_list(pt_neilist,i,j, patchesy); //generates list of neighbor patches
		list_migrants(r,pt_neilist, pt_chosen, seeds, mig); // generates list of colonists (migrants) to colonize patch (i,j)
		for(vector<int>::iterator iter=pt_chosen->begin(); iter!=pt_chosen->end(); iter++){
			(*seeds)[*iter].patch_x = i;
			(*seeds)[*iter].patch_y = j;
			(*seeds)[*iter].chosen = 10; 
			colonists++;
		}
				cout << "---------------------------------------------------------------> Size after colonization " << subpopsize(i,j,seeds) << endl;
		if(subpopsize(i,j,seeds)>0) spatial_occupancy[i][j]=1;
		else spatial_occupancy[i][j]=0;
		delete pt_neilist;
		delete pt_chosen;
	}
	
}




void checkpopsize(int patchesx, int patchesy, vector<indiv> *seeds, int **spatial_occupancy){
	for(int i=0; i<patchesx; i++){
		for(int j=0; j<patchesy; j++){
			if(spatial_occupancy[i][j] == 1){
				vector<indiv> *subpop = new vector<indiv>;
				build_subpop_vector_pointer2D(i,j,seeds,subpop);
				//cout << "Subpop size patch: " << i << " " << j << " : " << subpop->size() << endl;
				delete subpop;
			}
		}
	} 
}

double mean_fitness(vector<indiv> *subpop){
	double sumfit=0;
	for(vector<indiv>::iterator iter=subpop->begin(); iter!=subpop->end(); iter++)
		 sumfit += iter->fitness;
    return sumfit/subpop->size();
}




void gameteprod_mating(gsl_rng *r, vector<indiv> *parents, vector<indiv> *seeds, int genes, int patchesx, int patchesy, int N, double growth, double lambda_neutral, double lambda_sel, double mu_slocus, double selection_exp_mu, double dominance, double rec, double selfing_prob, int **space, vector<int> &Ntype, vector<int> &Stype, int flag, double ***p2outputSelfrate, int expa){
	int chosen_mother, Nj, popsize, self, L, stop;
	double sumfit, meanfit, maxfit;
	vector<int> ovule_haplotype, pollen_haplotype;
	pair <vector<int>, int > pollen_pair;
	for(int i=0; i<patchesx; i++){
		for(int j=0; j<patchesy; j++){

			if(space[i][j]==1){

				vector<indiv> *subpop = new  vector<indiv>;
				build_subpop_vector_pointer2D(i,j,parents,subpop);
				meanfit = mean_fitness(subpop);
				//if(meanfit > 0.0001){
	               stop = 0, self = 0;				
				   vector<int> S_alleles_mother;
				   //cout << "subpop size: " << subpop->size() << endl;			
				   cout << "mean fitness: " << meanfit << endl;	
				   maxfit = maxfit_indiv(subpop);	            
                   popsize = subpop->size();
				   if(expa > 0 && popsize < N)   Nj = min(N, (int) ceil(growth*popsize));
				   else Nj = N;
				   cout << "Nj " << Nj << endl;
				   sumfit = 0;
			       for(vector<indiv>::iterator iter=subpop->begin(); iter!=subpop->end(); iter++)
		              sumfit += iter->fitness;
				   for(int n=0; n < Nj; n++){
					  //cout << "Counter for loop: " << n << endl;
					  chosen_mother = gsl_rng_uniform_int(r, popsize);
					  do{
					    	chosen_mother = gsl_rng_uniform_int(r, popsize);
					  }while( gsl_rng_uniform(r) > ((*subpop)[chosen_mother].fitness)/maxfit);
					  //cout << "Fitness " << ((*subpop)[chosen_mother].fitness) << endl;
					
					  S_alleles_mother.push_back((*subpop)[chosen_mother].S_allele1());
					  S_alleles_mother.push_back((*subpop)[chosen_mother].S_allele2());
					  ovule_haplotype = create_haplotype(r, (*subpop)[chosen_mother].genome, rec, genes);
					  glass(ovule_haplotype);
					  L = ovule_haplotype.size();
					  if (S_alleles_mother[0] < 0 || S_alleles_mother[1] < 0){
						
						//selfing
					    	if(selfing_prob > gsl_rng_uniform(r)){
						    	pollen_haplotype = create_haplotype(r, (*subpop)[chosen_mother].genome, rec, genes);
							    glass(pollen_haplotype);
							    mutations(r, pollen_haplotype, lambda_neutral, lambda_sel, mu_slocus, Ntype, Stype, flag);
							    self++;
						    }else // outcrossing
						        if((S_alleles_mother[0] < 0 && S_alleles_mother[1] > 0)  || (S_alleles_mother[0] > 0 && S_alleles_mother[1] < 0)){
							       pollen_haplotype = find_mate_GSI_SC(r, subpop, S_alleles_mother, rec, lambda_neutral, lambda_sel, mu_slocus, genes, Ntype, Stype, flag);
						        }else if(S_alleles_mother[0] < 0 && S_alleles_mother[1] < 0){
								   pollen_haplotype = find_pollen(r, subpop, S_alleles_mother, rec, lambda_neutral, lambda_sel, mu_slocus, genes, Ntype, Stype, flag);
							    }
						
				 	  }else // GSI - outcrossing
					  {
						    pollen_haplotype = find_mate_GSI(r, subpop, S_alleles_mother, rec, lambda_neutral, lambda_sel, mu_slocus, genes, Ntype, Stype, flag);
					  }
						
					
					  vector<int>().swap(S_alleles_mother);
					  // If there aren't compatible mates (Allee effect) or all with fitness ~= 0
					  if(pollen_haplotype[0] == -1)   stop++;	

					  if(stop >= Nj)   continue;
					  if(pollen_haplotype[0] != -1){
					     indiv seed;
					     seed.genome = fertilization2(ovule_haplotype, pollen_haplotype);
					     seed.patch_x = i;
					     seed.patch_y = j;
					     seed.chosen = -10;
					     seed.fitness = calculatefitness(r, seed.genome, seed.genome.size()/2, selection_exp_mu, dominance);
					     seeds->push_back(seed);
				      } 
                      vector<int>().swap(pollen_haplotype);
                      vector<int>().swap(ovule_haplotype);
				  //}
				}

				if(self>0)  cout << "-------------------- >>> SELFING RATE: " << (double) self/Nj << endl;
				if(expa>0)   p2outputSelfrate[expa][i][j] = (double) self/Nj;
				
			    if(subpopsize(i,j,seeds) > 0) space[i][j] = 1;
			    else space[i][j] = 0;
				
			    delete subpop;
			   
			}
		}
	}
	
}



void show_vector_pointer(vector<double> *vec){
	cout << "Vector of length: " << vec->size() << endl;
	cout << " contains: ";
    for (vector<double>::iterator it=vec->begin(); it!=vec->end(); ++it)
         cout << ' ' << *it;
    cout << endl;
}


double variance(double mean, double data[], int size){
	double sqsum=0;
	int n=0;
	while(n<size){
		sqsum += (data[n] - mean)*(data[n]-mean);
		n++;
	}
	return sqsum/size;
	
}




double inbreeding_depression_within_deme(gsl_rng *r, vector<indiv> *subpop, double selection_exp_mu, double rec, double dominance){
	vector<int> pollen_haplotype, ovule_haplotype;
	int chosen_father, chosen_mother, chosen_ind, N = subpop->size(), n=0, L = (*subpop)[0].genome.size()/2;
	double sumrand=0.0, suminb=0.0; 
	vector<double> offspring(N);
    vector<double> inbred(N);
    vector<int> seed_genomer, seed_genomei;
	while(n<N){
		// Random mating
		chosen_mother = gsl_rng_uniform_int(r,subpop->size());
		ovule_haplotype = create_haplotype(r, (*subpop)[chosen_mother].genome, rec, L);
		chosen_father = gsl_rng_uniform_int(r,subpop->size());
		pollen_haplotype = create_haplotype(r, (*subpop)[chosen_father].genome, rec, L);
		seed_genomer = fertilization2(ovule_haplotype, pollen_haplotype);
        offspring[n] = calculatefitness(r,seed_genomer, L, selection_exp_mu, dominance);
		// Inbreeding
		chosen_ind = gsl_rng_uniform_int(r,subpop->size());
		ovule_haplotype = create_haplotype(r, (*subpop)[chosen_ind].genome, rec, L);
		pollen_haplotype = create_haplotype(r, (*subpop)[chosen_ind].genome, rec, L);
		seed_genomei = fertilization2(ovule_haplotype, pollen_haplotype);
		inbred[n] = calculatefitness(r, seed_genomei, L, selection_exp_mu, dominance);
		n++;
	}
	//Calculate mean fitness of offspring and inbred
	for(vector<double>::iterator iter=offspring.begin(); iter!= offspring.end(); iter++){
	      sumrand += *iter;
	}
    for(vector<double>::iterator iter=inbred.begin(); iter!= inbred.end(); iter++){
	      suminb += *iter;
	} 
	return 1.0 - ( (suminb/N)/(sumrand/N) );
}

double inbreeding_depression_between_deme(gsl_rng *r, vector<indiv> *seeds, vector<indiv> *subpop, int patchesx, int patchx, int patchy, double selection_exp_mu, double rec, double dominance, int **space){
	///vector<int> pollen_haplotype, ovule_haplotype;
	int chosen_father, chosen_mother, chosen_ind, N = subpop->size(), n=0, L = (*subpop)[0].genome.size()/2, chosen_row;
	double sumrand=0.0, suminb=0.0; 
	vector<double> offspring(N);
    vector<double> inbred(N);
    vector<int> seed_genomer, seed_genomei;
    
	while(n<N){
		vector<int> *pollen_haplotype = new vector<int>;
		vector<int> *ovule_haplotype = new vector<int>;
		vector<int> *pollen_haplotype_inb = new vector<int>;
		vector<int> *ovule_haplotype_inb = new vector<int>;
		// Random mating
		vector<indiv> *subpop_temp = new vector<indiv>;
		// First pick first a mother from patch (X,Y)
		chosen_mother = gsl_rng_uniform_int(r,subpop->size());
		(*ovule_haplotype) = create_haplotype(r, (*subpop)[chosen_mother].genome, rec, L);
		// Second pick a father from another patch in the column (e.g. patch (X,j) j != Y
		do{
			chosen_row = gsl_rng_uniform_int(r, patchesx);
			do{
				chosen_row = gsl_rng_uniform_int(r, patchesx);
			}while(chosen_row==patchx);
		}while(space[chosen_row][patchy]==0);
        build_subpop_vector_pointer2D(chosen_row,patchy,seeds,subpop_temp);
        //cout << "Subpop size " << subpop_temp->size() << " " << space[chosen_row][patchy] << " " << chosen_row << " " << patchy << endl;
        chosen_father = gsl_rng_uniform_int(r,subpop_temp->size());
		(*pollen_haplotype) = create_haplotype(r, (*subpop_temp)[chosen_father].genome, rec, L);
		seed_genomer = fertilization(ovule_haplotype, pollen_haplotype);
        offspring[n] = calculatefitness(r,seed_genomer, L, selection_exp_mu, dominance);
		// Inbreeding
		chosen_ind = gsl_rng_uniform_int(r,subpop->size());
		(*ovule_haplotype_inb) = create_haplotype(r, (*subpop)[chosen_ind].genome, rec, L);
		(*pollen_haplotype_inb) = create_haplotype(r, (*subpop)[chosen_ind].genome, rec, L);
		seed_genomei = fertilization(ovule_haplotype_inb, pollen_haplotype_inb);
		inbred[n] = calculatefitness(r, seed_genomei, L, selection_exp_mu, dominance);
		n++;
		delete subpop_temp;
		delete pollen_haplotype;
		delete ovule_haplotype;
		delete pollen_haplotype_inb;
		delete ovule_haplotype_inb;
	}
	//Calculate mean fitness of offspring and inbred
	for(vector<double>::iterator iter=offspring.begin(); iter!= offspring.end(); iter++)
	      sumrand += *iter;
    for(vector<double>::iterator iter=inbred.begin(); iter!= inbred.end(); iter++)
	      suminb += *iter;

	return 1.0 - ( (suminb/N)/(sumrand/N) );
}


//---------> Create function to Calculate mate availability based on counting the number of compatible pollen donors for each plant available in the sample of the same population (Glemin et al 2008)
double mate_availability(gsl_rng *r, vector<indiv> *subpop)
{
	// Loop through all ovules in the subpop:
	int  N = subpop->size(), x, n=0, chosen_mother, chosen_father;
	double compat=0;
	vector<int> S_alleles_mother, S_alleles_father;
	while(n<N){
		chosen_mother = gsl_rng_uniform_int(r,subpop->size());
		S_alleles_mother.push_back((*subpop)[chosen_mother].S_allele1());
		S_alleles_mother.push_back((*subpop)[chosen_mother].S_allele2());
		chosen_father = gsl_rng_uniform_int(r,subpop->size());
	    S_alleles_father.push_back((*subpop)[chosen_father].S_allele1());
	    S_alleles_father.push_back((*subpop)[chosen_father].S_allele2());
	    if(0.5 > gsl_rng_uniform(r)) x = 0;
	    else x = 1;
	    if(S_alleles_mother[0] != S_alleles_father[x] && S_alleles_mother[1] != S_alleles_father[x])   compat++;
	    //if(S_alleles_mother[0] != S_alleles_father[1] && S_alleles_mother[1] != S_alleles_father[1])   comp++;
	    n++;
	}
	//cout << " ------------------------------------------------------------->>>>>> MATE AVAIL: " << compat << endl;
	return compat/N;
}

bool check_expansion(vector<indiv> *population, int patchesx, int **space){
	int hit=0;
	bool temp;
	for(int i=0; i<patchesx; i++){
		if(space[i][99]==1)    hit++;
	}
	if(hit==patchesx)  temp=true;
	else     temp=false;
	return temp;
}

/// Function to record S alleles 

void s_alleles_record(int gen, vector<indiv> *population, int patchesx, int patchesy, int **space, double ***** p2_S_alleles, vector<int> &Stype){
	int L = (*population)[0].genome.size();
	unsigned int s_alleles =200;

	for(int k=0; k<patchesx; k++){
		for(int l=0; l<patchesy; l++){
			if(space[k][l]==1){
				int total_S_allele = 0;
				vector<indiv> *subpop = new vector<indiv>;
				vector<int> *s_alleles_count = new vector<int>;
				int *counts_S = new int[200];
				double *freqs_S = new double[200];
				build_subpop_vector_pointer2D(k,l,population,subpop);
				
				
				for(unsigned int s_al=0; s_al < s_alleles; s_al++){
					int delet = 0;
					for(std::vector<indiv>::iterator it = subpop->begin(); it != subpop->end(); it++){
						s_alleles_count->push_back(it->genome[(L/4) - 1]);
			            s_alleles_count->push_back(it->genome[(3*L/4) - 1]);
						if(Stype[s_al] ==  it->genome[(L/4) - 1]){
							// Counting deleterious alleles associated to one S allele
							for(int kk=0; kk<50; kk++){ 
							   if(it->genome[(L/4)+kk] == 23) delet++;
							   if(it->genome[(L/4)-kk] == 23) delet++;
						    }

					    }
					    if(Stype[s_al] ==  it->genome[(3*L/4) - 1]){
							// Counting deleterious alleles associated to one S allele
							for(int kk=0; kk<50; kk++){ 
							   if(it->genome[(3*L/4)+kk] == 23) delet++;
							   if(it->genome[(3*L/4)-kk] == 23) delet++;
						    }
					    }					    	            
					}
					//cout << "S allele " << Stype[s_al] << " index S-type " << s_al << " patch " << k << " " << l << " generation " << gen << endl;
					p2_S_alleles[gen][k][l][s_al][1] = (double) delet/(2*subpop->size());  // average number of deleterious mutations
					//cout << "S allele " << Stype[s_al] << " delet " << p2_S_alleles[gen][k][l][s_al][1] << endl;
				}
				for(unsigned int s_al=0; s_al < s_alleles; s_al++){
					counts_S[s_al] = std::count(s_alleles_count->begin(), s_alleles_count->end(), Stype[s_al]);
					total_S_allele += counts_S[s_al];
				}
				for(unsigned int s_al=0; s_al < s_alleles; s_al++){
					freqs_S[s_al] =  (double) counts_S[s_al]/total_S_allele;
					//cout << "Frequency of S allele: " << Stype[s_al] << "  "  << freqs_S[s_al] << endl;
					// Putting output to pointer array 
					p2_S_alleles[gen][k][l][s_al][0] = freqs_S[s_al];
				}
				
			    delete subpop;
			    delete s_alleles_count;
			    delete[] counts_S;
		        delete[] freqs_S;
			}
		}
	}
}


/// Selection coefficient distribution function 
void selection_coef_deme2(gsl_rng *r, vector<indiv> *subpop, double *p2output, int genes, double mu){
	int L=genes,allele1, allele2, count=0;
	double s;
	for(vector<indiv>::iterator it=subpop->begin(); it!= subpop->end(); it++){
		for(int i=0; i<L; i++){
		    do{
			  s = gsl_ran_exponential(r, mu);
		    }while(s > 1.0);
            allele1 = it->genome[i];
            allele2 = it->genome[i+L];
            cout.precision(15);
            //cout << " Sel coef " << fixed << s << endl;
		    if (allele1==21 || allele1 ==23 || allele2==21 || allele2==23){
			//cout << "Allele 1 " << allele1  <<  " Allele 2 " << allele2 << endl;
		        if(allele2+allele1==46){		       
			       p2output[count]=s;
			       //cout << "Deleterious " << fitness[count] << " " << allele1 << " " << allele2 << endl;
		        }else if(allele1+allele2==44){
			       p2output[count]=s;
			       //cout << "Hetero "  << (double) fitness[count] << " " << allele1 << " " << allele2 << endl;
			    }
			    count++;
		   }
			
	    }
		
	}
	
}

// Calculate allele freqs for all patches
void allele_frequencies_2(gsl_rng *r, vector<indiv> *population, double **stat_year_vec, int interval, int patchesx, int patchesy, int selgenes, int ngenes, int S_alleles, vector<int> &Stype, vector<int> &Ntype, int year, int **space, double selection_exp_mu, double rec, double dominance, double ****p2output, int gen, double mu_selection){
    double sumfit_total=0, sum_inbreeding=0, sum_Ae=0, sum_Hs=0, sum_Hs2=0, sum_genAA=0.0, sum_Hneutral=0, sum_quad_freq_Npop=0, sum_quad_freq_Npop_far=0, sum_Hexp_nfar=0, sum_p=0,sum_p2=0, sum_p_onelocus=0, sum_p2_onelocus=0,sum_q=0, sum_freq_quad_S_total=0, sum_Hs_S=0, Ht, Hs, Hs2, Hs_S, Ht_S, Hs_N, Ht_N, Hs_N_far, Ht_N_far, freq_p_total, freq_q_total, freq_p_onelocus ,average_genA, genAA, genaa, genAa, variance_pop, sum_between_inbdep=0; //p_ave_vec[patchesx*patchesy], 
    int L, number_of_demes, demes=0, nalleles = Ntype.size(), totselfcomp=0;
    vector<double> *sum_freq_S = new vector<double>;
    for(unsigned int s_al=0; s_al < Stype.size(); s_al++)  sum_freq_S->push_back(0.0);
    vector<double> *sum_freq_N = new vector<double>;
    vector<double> *sum_freq_N_far = new vector<double>;
    vector<double> *mean_subpop_allele_size_vec = new vector<double>;
    vector<double> *sum_sd_allele_size_vec = new vector<double>;
    vector<double> *neutral_loci_total = new vector<double>;
    for(unsigned int n_al=0; n_al <Ntype.size(); n_al++){
		  sum_freq_N->push_back(0.0);
		  sum_freq_N_far->push_back(0.0);
	}
	vector<int> selfcomps;
	for(int k=0; k<patchesx; k++){
		for(int l=0; l<patchesy; l++){
			if(space[k][l]==1){
				vector<indiv> *subpop = new vector<indiv>;
				build_subpop_vector_pointer2D(k,l,population,subpop);
				double allele_b, allele_B, hetero_b=0, heteron_far=0, deleterious_alleles=0, Hobs_nfar=0, Hobs_nlinked=0, Hexp_nlinked=0, Hexp_nfar=0, heteron_linked=0, homo_b =0, homo_B=0, H_neutral, sumb=0, sumB=0, p_average=0, q_average=0, Hs_S_d=0 ,Hexp_b, Hobs_b, p_onelocus=0, sum_quad_freq_S=0, sum_quad_freq_neutral=0, sum_quad_freq_all_loci_neutral=0, betdeminbdep, inbdep;
				int total_S_allele=0, total_neutral=0, nalleles = Ntype.size(), cc, allele1, allele2;
				vector<double> *locus_b = new vector<double>;
				vector<double> *locus_B = new vector<double>;
				vector<int> *s_alleles_count = new vector<int>;
				vector<double> *freqs_S = new vector<double>;
				vector<double> *mean_freq_nloci = new vector<double>;
				vector<double> *sum_freq_N = new vector<double>(nalleles);
				vector<int> *counts_S = new vector<int>;
				vector< vector<int> > *neutral_loci = new vector< vector<int> >;
				double sumfitness=0.0;
                //cout<< "subpop size: "<< subpop->size() << " patch: " << k << " " << l << endl;
                //if(k==0 && l==10){
					//cout << "Contains :";
					//for(vector<indiv>::iterator it=subpop->begin(); it!=subpop->end(); it++)
					  //   cout << it->S_allele1() << "-" << it->S_allele2() << " ";
					//cout << endl; 
				//}
			    //Calculate allele frequencies for each selected locus and neutral locus
			    L = (*subpop)[0].genome.size();
			    int delet = 0;
			    for(int gene=0; gene<L; gene++){
					allele_b=0;
					allele_B=0;
					// Selected and neutral loci
					for(vector<indiv>::iterator it=subpop->begin(); it!=subpop->end(); it++){
						if(it->genome[gene] == 21)  allele_b++;
						if(it->genome[gene] == 23)  allele_B++;
				    }
				    locus_b->push_back(allele_b/subpop->size());
				    locus_B->push_back(allele_B/subpop->size());
				    
				    delet += allele_B;
				    //cout << "Deleterious alleles per locus: " << delet << endl;
			    }
			    

			     // Neutral loci
			    for(int gene=0; gene<L/2; gene++){
					if(gene%10==0){
						vector<int> *neutral_locus_subpop = new vector<int>;
						for(vector<indiv>::iterator it=subpop->begin(); it!=subpop->end(); it++){
							neutral_locus_subpop->push_back(it->genome[gene]);
							neutral_locus_subpop->push_back(it->genome[gene + L/2]);
						}
					    neutral_loci->push_back(*neutral_locus_subpop);
					    delete neutral_locus_subpop;
					}
				}
			    
			    //cout << "Number of neutral loci: " << neutral_loci->size() << "Selected: " << (L/2) - neutral_loci->size() - 1 << endl;
			    number_genes((*subpop)[0].genome);
			    int selected_genes = (L/2) - neutral_loci->size() - 1;
			    //cout << "Size genome " << (*subpop)[0].genome.size() << endl;

			    // Calculate number of heterozygotes
			    int selfcomp=0;		
			    for(std::vector<indiv>::iterator it = subpop->begin(); it != subpop->end(); it++){
			        for (int j=0; j < L/2; j++){
					   allele1 = it->genome[j];
					   allele2 = it->genome[j+L/2];
					   
					   if(allele1+allele2 == 44){
						   hetero_b++;
					   }else if(allele2+allele1 == 42){
						   homo_b++;
					   }else if(allele2+allele1 == 46){
						   homo_B++;
					   }
					   if(j==0){
						   if(it->genome[j] != it->genome[j+L/2]) heteron_far++;
					   }else if(j==240){
						   if(it->genome[j] != it->genome[j+L/2])  heteron_linked++;
					   }
			        }
			        sumfitness += it->fitness;
			        s_alleles_count->push_back(it->genome[(L/4) - 1]);
			        s_alleles_count->push_back(it->genome[(3*L/4) - 1]);
			        if(it->genome[(L/4) - 1] < 0) selfcomp++;
			        if(it->genome[(3*L/4) - 1] < 0) selfcomp++;
			        
			        
			    }
			    selfcomps.push_back(selfcomp);
			    cout << "hetero " << hetero_b << " homo_b " << homo_b << " homo_B " << homo_B << endl;
			    cout << hetero_b/(selected_genes*subpop->size()) << " " << homo_b/(selected_genes*subpop->size()) << " " << homo_B/(selected_genes*subpop->size()) << " " << endl;
			    
			    /// Calculate neutral allele frequencies:
			    double sum_neutral_allele_sizes = 0, mean_subpop_allele_size;
			    //for(vector<vector<int> >::iterator iter=neutral_loci->begin(); iter!=neutral_loci->end(); iter++){
			    for(unsigned int nloc=0; nloc<neutral_loci->size(); nloc++){
					cc=0, total_neutral=0;
					vector<int> *counts_neutral = new vector<int>;
					for(vector<int>::iterator neu=Ntype.begin(); neu!=Ntype.end(); neu++){
						counts_neutral->push_back(std::count((*neutral_loci)[nloc].begin(), (*neutral_loci)[nloc].end(), *neu));
						total_neutral += (*counts_neutral)[cc];
						cc++;
					}
					sum_quad_freq_neutral=0.0;
					vector<double> *freqs_neutral = new vector<double>;
                    for(int ne=0; ne<nalleles; ne++){
						freqs_neutral->push_back( (double) (*counts_neutral)[ne]/total_neutral);
						//cout << "Frequency of neutral allele: " << Ntype[ne] << "  "  << (*freqs_neutral)[ne] << endl;
						(*sum_freq_N)[ne] += (*freqs_neutral)[ne]; // Sum of neutral allele frequencies
						sum_quad_freq_neutral += (*freqs_neutral)[ne]*(*freqs_neutral)[ne]; // ---> to have Hs
						if(nloc == 0)  (*sum_freq_N_far)[ne] += (*freqs_neutral)[ne];
					}
					if(nloc == 0) Hexp_nfar = 1.0 - (sum_quad_freq_neutral);
					else if(nloc == 24) Hexp_nlinked = 1.0 - (sum_quad_freq_neutral);
					sum_quad_freq_all_loci_neutral += 1.0 - (sum_quad_freq_neutral);
					delete freqs_neutral;
					delete counts_neutral;
					// Sum of Subpopulation of allele sizes:
                    if(nloc == 10){
					   for(vector<int>::iterator al=(*neutral_loci)[nloc].begin(); al!=(*neutral_loci)[nloc].end(); al++){
						  sum_neutral_allele_sizes += *al;
						  neutral_loci_total->push_back(*al);						
					   }
				   }
				}
				mean_subpop_allele_size = (double) sum_neutral_allele_sizes/(2*subpop->size());
				// Mean Subpopulation of allele sizes:
				mean_subpop_allele_size_vec->push_back( (double) sum_neutral_allele_sizes/(2*subpop->size())); // divided by number of neutral loci * subpop size
				// Sum of squared difference (Rst)
				double sum_sd_allele_size = 0;
				
				for(vector<int>::iterator al=(*neutral_loci)[10].begin(); al!=(*neutral_loci)[10].end(); al++){
						sum_sd_allele_size += (*al - mean_subpop_allele_size)*(*al - mean_subpop_allele_size);
				}

				//double sum_sd_allele_size2 = (double) sum_sd_allele_size/(2*subpop->size());
				sum_sd_allele_size_vec->push_back( (double) sum_sd_allele_size/(2*subpop->size()));
				//---------------------------------------------------
				//cout << "loci size " << neutral_loci->size() << endl;
				H_neutral = sum_quad_freq_all_loci_neutral/(neutral_loci->size());
				sum_Hneutral += H_neutral;
				//cout << "Average expected heterozygocity " << H_neutral << endl; //<<------------------------
				// Mean frequency value for each neutral allele
                for(int ne=0; ne<nalleles; ne++){
					mean_freq_nloci->push_back((double) (*sum_freq_N)[ne]/neutral_loci->size());
					//cout << "Frequency of neutral allele: " << Ntype[ne] << "  "  << (*mean_freq_nloci)[ne] << endl;
				}
					
			    /// Calculate S allele frequencies:
			    for(unsigned int s_al=0; s_al < Stype.size(); s_al++){
					counts_S->push_back(std::count(s_alleles_count->begin(), s_alleles_count->end(), Stype[s_al]));
					total_S_allele += (*counts_S)[s_al];
				} 
				for(unsigned int s_al=0; s_al < Stype.size(); s_al++){
					freqs_S->push_back( (double) (*counts_S)[s_al]/total_S_allele);
					//cout << "Frequency of S allele: " << Stype[s_al] << "  "  << (*freqs_S)[s_al] << endl;
					sum_quad_freq_S += (*freqs_S)[s_al]*(*freqs_S)[s_al];
					(*sum_freq_S)[s_al] += (*freqs_S)[s_al];
				}
               // Calculate average allele frequecies:
                for(vector<double>::iterator iter=locus_b->begin(); iter!=locus_b->end(); iter++){
				   sumb += *iter;
	            }
	            p_average = ((homo_b/selected_genes)+0.5*(hetero_b/selected_genes))/subpop->size();
	            for(vector<double>::iterator iter2=locus_B->begin(); iter2!=locus_B->end(); iter2++){
				   sumB += *iter2;
			    }
			    q_average = ((homo_B/selected_genes)+0.5*(hetero_b/selected_genes))/subpop->size();  // <<--- Deleterious allele frequency
			    deleterious_alleles = ((homo_B)+0.5*(hetero_b))/subpop->size();
			    p_onelocus = (*locus_B)[2];
                Hexp_b = 1.0 - q_average*q_average - p_average*p_average;
                Hs_S_d = 1.0 - sum_quad_freq_S;
                Hobs_b = hetero_b/(selgenes*subpop->size());
                Hobs_nfar = heteron_far/subpop->size();
                Hobs_nlinked = heteron_linked/subpop->size();
                inbdep = inbreeding_depression_within_deme(r, subpop, selection_exp_mu, rec, dominance);
                //cout << "Selected locus freq: " << k << " " << l << " " << p_onelocus << endl;
		        //cout << "Allele freq A in patch: " << k << " " << l << " " << p_average << endl;
		        //cout << "Allele freq a in patch: " << k << " " << l << " " << q_average << endl;
		        //cout << "Observed HET-A freq in patch: " << k << " " << l << " " << Hobs_b << endl;
		        //cout << "Expected HET-A freq in patch: " << k << " " << l << " " << Hexp_b << endl;
	            //cout << "Inbreeding coef (Fis): " << k << " " << l << " " << 1-(Hobs_b/Hexp_b) << endl;
	            //cout << "Mean fitness: " << sumfitness/subpop->size() << endl;
	            //cout << "Genetic load: " << 1 - sumfitness/subpop->size() << endl;
	            cout << "Inbreeding depression:  " << inbdep << endl;
	            //cout << "Effective number of S alleles: " << 1/sum_quad_freq_S << endl;
	            //cout << "Mate availability: " << mate_availability(r,subpop) << endl;
	            //cout << "Heterozygotes observed linked: " << Hobs_nlinked << endl;
	            //cout << "Inbreeding coef linked: " << 1.0 - (Hobs_nlinked/Hexp_nlinked) << endl;
	            //cout << "Heterozygotes observed far: " << Hobs_nfar << endl;
	            //cout << "Inbreeding coef far: " << 1.0 - (Hobs_nfar/Hexp_nfar) << endl;
	            //cout << "Deleterious allele number: " << deleterious_alleles << " allele B "  << delet/subpop->size() << endl;
	            betdeminbdep = inbreeding_depression_between_deme(r, population, subpop, patchesx, k, l, selection_exp_mu, rec, dominance,space);
	            sum_between_inbdep += betdeminbdep;
	            sum_Hs+= Hobs_b;
	            sum_Hs2 += Hexp_b;
                sum_Hs_S += Hs_S_d;
                sum_p_onelocus += p_onelocus;
                sum_p2_onelocus += p_onelocus*p_onelocus;
	            sum_p+= p_average;
	            sum_p2+= p_average*p_average;
	            sum_q+= q_average;
	            sum_genAA += homo_b/(selgenes*subpop->size());
	            sumfit_total += sumfitness/subpop->size();
	            sum_inbreeding += inbreeding_depression_within_deme(r, subpop, selection_exp_mu, rec, dominance);
	            sum_Ae += 1/sum_quad_freq_S;
	            sum_Hexp_nfar += Hexp_nfar;
	            /// OUTPUT DEME STATS
	            if(gen == 2000){
					p2output[0][k][l][0] = inbdep;
					p2output[0][k][l][1] = 1 - sumfitness/subpop->size();
					p2output[0][k][l][2] = 1/sum_quad_freq_S;    // Fis mean/var Fis G (glemin expected) Expected fitness Whitlock 
					p2output[0][k][l][3] = 1.0 - (Hobs_nfar/Hexp_nfar);
					p2output[0][k][l][4] = 1.0 - (Hobs_nlinked/Hexp_nlinked);
					p2output[0][k][l][5] = -1.0/(S_alleles + 1 + 4*rec*subpop->size());
					p2output[0][k][l][6] = mate_availability(r,subpop);
					p2output[0][k][l][7] = deleterious_alleles;
					p2output[0][k][l][8] = selfcomp;
				}else if(gen == 2100){
					p2output[1][k][l][0] = inbdep;
					p2output[1][k][l][1] = 1 - sumfitness/subpop->size();
					p2output[1][k][l][2] = 1/sum_quad_freq_S;    // Fis mean/var Fis G (glemin expected) Expected fitness Whitlock 
					p2output[1][k][l][3] = 1.0 - (Hobs_nfar/Hexp_nfar);
					p2output[1][k][l][4] = 1.0 - (Hobs_nlinked/Hexp_nlinked);
					p2output[1][k][l][5] = -1.0/(S_alleles + 1 + 4*rec*subpop->size());
					p2output[1][k][l][6] = mate_availability(r,subpop);
					p2output[1][k][l][7] = deleterious_alleles;
					p2output[1][k][l][8] = selfcomp;
				
				}else if(gen == 2200){
					p2output[2][k][l][0] = inbdep;
					p2output[2][k][l][1] = 1 - sumfitness/subpop->size();
					p2output[2][k][l][2] = 1/sum_quad_freq_S;    // Fis mean/var Fis G (glemin expected) Expected fitness Whitlock 
					p2output[2][k][l][3] = 1.0 - (Hobs_nfar/Hexp_nfar);
					p2output[2][k][l][4] = 1.0 - (Hobs_nlinked/Hexp_nlinked);
					p2output[2][k][l][5] = -1.0/(S_alleles + 1 + 4*rec*subpop->size());
					p2output[2][k][l][6] = mate_availability(r,subpop);
					p2output[2][k][l][7] = deleterious_alleles;
					p2output[2][k][l][8] = selfcomp;
					
				}else if(gen == 2300){
					p2output[3][k][l][0] = inbdep;
					p2output[3][k][l][1] = 1 - sumfitness/subpop->size();
					p2output[3][k][l][2] = 1/sum_quad_freq_S;    // Fis mean/var Fis G (glemin expected) Expected fitness Whitlock 
					p2output[3][k][l][3] = 1.0 - (Hobs_nfar/Hexp_nfar);
					p2output[3][k][l][4] = 1.0 - (Hobs_nlinked/Hexp_nlinked);
					p2output[3][k][l][5] = -1.0/(S_alleles + 1 + 4*rec*subpop->size());
					p2output[3][k][l][6] = mate_availability(r,subpop);
					p2output[3][k][l][7] = deleterious_alleles;
					p2output[3][k][l][8] = selfcomp;
				
				}else if(gen == 3000){
					p2output[4][k][l][0] = inbdep;
					p2output[4][k][l][1] = 1 - sumfitness/subpop->size();
					p2output[4][k][l][2] = 1/sum_quad_freq_S;    // Fis mean/var Fis G (glemin expected) Expected fitness Whitlock 
					p2output[4][k][l][3] = 1.0 - (Hobs_nfar/Hexp_nfar);
					p2output[4][k][l][4] = 1.0 - (Hobs_nlinked/Hexp_nlinked);
					p2output[4][k][l][5] = -1.0/(S_alleles + 1 + 4*rec*subpop->size());
					p2output[4][k][l][6] = mate_availability(r,subpop);
					p2output[4][k][l][7] = deleterious_alleles;
					p2output[4][k][l][8] = selfcomp;
					
				}
	            //vector<double>().swap(ib_block);
 		        delete subpop;
		        delete locus_b;
		        delete locus_B;
		        delete s_alleles_count;
		        delete freqs_S;
		        delete counts_S;
		        delete neutral_loci;
		        delete mean_freq_nloci;
		        delete sum_freq_N;
		        demes++;
			}
		}
	}
	/// METAPOP AVERAGE
	number_of_demes = demes;
	cout << "number of demes occupied: " << number_of_demes << endl;
	// Calculate Fst at each interval
	Hs = sum_Hs/(number_of_demes);
	Hs2 = sum_Hs2/(number_of_demes);
	cout << "Hs: " << Hs << endl;
	cout << "Hs2:    " << Hs2 << endl;
	freq_p_total = sum_p/(number_of_demes);
    freq_q_total = sum_q/(number_of_demes);
    freq_p_onelocus = sum_p_onelocus/number_of_demes;
    genAA = freq_p_total*freq_p_total; //Pooled HW
    genAa = 2*freq_p_total*freq_q_total; //Pooled HW
    genaa = freq_q_total*freq_q_total; //Pooled HW
    average_genA = sum_p2_onelocus/number_of_demes;//sum_genAA/number_of_demes;
    cout << "Homozygote A total freq: " << genAA << endl;
    cout << "Homozygote a total freq: " << genaa << endl;
    cout << "Heterozygote a total freq: " << genAa << endl;
    cout << "Average genA " << average_genA << endl;
    variance_pop = average_genA - (freq_p_onelocus*freq_p_onelocus);
    //cout << "Variance pop: " << variance_pop << endl;
    // Calculations for neutral loci
    for(int ne=0; ne<nalleles; ne++){
		sum_quad_freq_Npop += (*sum_freq_N)[ne]/(number_of_demes*ngenes)*(*sum_freq_N)[ne]/(number_of_demes*ngenes);
		sum_quad_freq_Npop_far += (*sum_freq_N_far)[ne]/(number_of_demes)*(*sum_freq_N_far)[ne]/(number_of_demes);
	}
    Ht_N = 1.0 - sum_quad_freq_Npop;
    Ht_N_far = 1.0 - sum_quad_freq_Npop_far;
    Hs_N = sum_Hneutral/number_of_demes;
    Hs_N_far = sum_Hexp_nfar/number_of_demes;
    // Calculate Rst
    double sum_mean_subpop_allele_size_vec = 0, sum_sum_sd_allele_size = 0, Sw, Sb, Rst, sum_square_dif_total_allele_size = 0, mean_sum_square_dif_total_allele_size, mean_total_allele_size, sum_all_allele_sizes = 0;
    /// Calculate mean of mean allele size across populations
    for(int d=0; d < number_of_demes; d++){
		sum_mean_subpop_allele_size_vec += (*mean_subpop_allele_size_vec)[d];
		//cout << " Sum sd allele size " << (*sum_sd_allele_size_vec)[d] << endl;
		sum_sum_sd_allele_size += (*sum_sd_allele_size_vec)[d];
	}
	
	for(vector<double>::iterator al=neutral_loci_total->begin(); al!=neutral_loci_total->end(); al++)  sum_all_allele_sizes += *al;
	mean_total_allele_size = (double) sum_all_allele_sizes/(2*population->size());
	for(vector<double>::iterator al=neutral_loci_total->begin(); al!=neutral_loci_total->end(); al++){
		sum_square_dif_total_allele_size += (*al - mean_total_allele_size)*(*al - mean_total_allele_size);
	}
	mean_sum_square_dif_total_allele_size = sum_square_dif_total_allele_size/(2*population->size());
	Sb = mean_sum_square_dif_total_allele_size ;
	Sw = sum_sum_sd_allele_size/number_of_demes;
	Rst = (Sb - Sw)/Sb;
	//cout << " Rst - Genetic structure of microsatellite: " << Rst << " Sb " << Sb << " Sw " << Sw << endl;
    /// Calculate sum of squared difference between populations - Sb
    // Calculate Gst for S-locus:
    for(int s_al=0; s_al < S_alleles; s_al++)    (*sum_freq_S)[s_al] = (*sum_freq_S)[s_al]/number_of_demes;
    for(int s_al=0; s_al < S_alleles; s_al++)    sum_freq_quad_S_total += (*sum_freq_S)[s_al]*(*sum_freq_S)[s_al];
    
    for(int patch=0; patch<patchesy; patch++){
		cout << " SC " << selfcomps[patch];
		totselfcomp += selfcomps[patch];
	}
    cout << endl;
    Hs_S = sum_Hs_S/number_of_demes;
    Ht_S = 1.0 - sum_freq_quad_S_total;
	Ht = 2*freq_p_total*freq_q_total; //- 2*variance_pop;
	//cout << "Ht: " << Ht << endl;
	//cout << "Wright Heterozygocity: " << 2*freq_p_total*freq_q_total*( 1 - (variance_pop/(freq_p_total*(1-freq_p_total)) )) << endl;
    //cout << " ********************   Fst selected loci: " << variance_pop/(freq_p_onelocus*(1-freq_p_onelocus)) << endl;
    //cout << " ********************   Gst neutral loci: " << 1.0 - Hs_N/Ht_N  << endl;
    //cout << " ********************   Gst S locus: " << 1.0 - Hs_S/Ht_S << endl;
	stat_year_vec[interval][0] = (double) year;
	stat_year_vec[interval][1] = freq_p_total;  /// selected locus
	stat_year_vec[interval][2] = Hs; 
	stat_year_vec[interval][3] = Ht;
	stat_year_vec[interval][4] = variance_pop/(freq_p_onelocus*(1-freq_p_onelocus));
	stat_year_vec[interval][5] = Hs_S; /// S locus
	stat_year_vec[interval][6] = Ht_S;//
	stat_year_vec[interval][7] = 1.0 - Hs_S/Ht_S;//(1.0 - (Hs/Ht));
	stat_year_vec[interval][8] = Hs_N; /// S locus
	stat_year_vec[interval][9] = Ht_N;//
	stat_year_vec[interval][10] = 1.0 - Hs_N/Ht_N;
    stat_year_vec[interval][11] = sumfit_total/number_of_demes; /// mean fitness
    stat_year_vec[interval][12] = 1.0 - sumfit_total/number_of_demes; /// mean genetic load
    stat_year_vec[interval][13] = sum_inbreeding/number_of_demes; /// Inbreeding depression mean
    stat_year_vec[interval][14] = sum_between_inbdep/number_of_demes; /// Inbreeding depression mean between demes
    stat_year_vec[interval][15] = sum_Ae/number_of_demes; /// Mean effective allele number
    stat_year_vec[interval][16] = 1.0 - Hs_N_far/Ht_N_far; /// Gst far linked neutral locus
    stat_year_vec[interval][17] = totselfcomp; /// Self-compatibles total metapop
    stat_year_vec[interval][18] = Rst; /// Self-compatibles total metapop
	delete sum_freq_S;
	delete sum_freq_N;
	delete sum_freq_N_far;
	delete mean_subpop_allele_size_vec;
	delete sum_sd_allele_size_vec;
	delete neutral_loci_total;
}

void final_year_stat(double **stat_year_vec, double **p2output_sims, int sim, int interval, int generations)
{
	int final_year = (generations/interval) - 1;
	p2output_sims[sim][0] = (double) sim;
	p2output_sims[sim][1] = stat_year_vec[final_year][1]; // freq p
	p2output_sims[sim][2] = stat_year_vec[final_year][2]; // Hs
	p2output_sims[sim][3] = stat_year_vec[final_year][3]; // Ht
	p2output_sims[sim][4] = stat_year_vec[final_year][4]; // Fst sel
	p2output_sims[sim][5] = stat_year_vec[final_year][5]; // Hs S
	p2output_sims[sim][6] = stat_year_vec[final_year][6]; // Ht S
	p2output_sims[sim][7] = stat_year_vec[final_year][7]; // Gst S
	p2output_sims[sim][8] = stat_year_vec[final_year][8]; // Hs N
	p2output_sims[sim][9] = stat_year_vec[final_year][9];  // Ht N
	p2output_sims[sim][10] = stat_year_vec[final_year][10]; // Gst N
	p2output_sims[sim][11] = stat_year_vec[final_year][11]; // Mean W
	p2output_sims[sim][12] = stat_year_vec[final_year][12]; // Load
	p2output_sims[sim][13] = stat_year_vec[final_year][13]; // Inb Dep
	p2output_sims[sim][13] = stat_year_vec[final_year][14]; // Inb Dep between
	p2output_sims[sim][14] = stat_year_vec[final_year][15]; // AE number S
	p2output_sims[sim][15] = stat_year_vec[final_year][16]; // Gst N far (nloc=0)
	p2output_sims[sim][16] = stat_year_vec[final_year][17]; // Total SC metapop
	p2output_sims[sim][17] = stat_year_vec[final_year][18]; // Rst
}


/// Output function to quantify the number of SC alleles in the metapopulation

void self_compatibles_output(gsl_rng *r, int generation, int sim, int patchesx, int patchesy, double ****p2outputSC, double ***p2outputSelfrate, vector<indiv> *seeds, int **space, double selection_exp_mu, double rec, double dominance){
	int subpopsize;
	double sumfit=0, sc, s100;
	for(int i=0; i<patchesx; i++){
		for(int j=0; j<patchesy; j++){
			if(space[i][j] == 1){
			  vector<indiv> *subpop = new vector<indiv>;
			  build_subpop_vector_pointer2D(i,j,seeds,subpop);
			  //cout << "Subpop size " << subpop->size() << " " << i << " " << j << endl;
			  sc=0,s100=0,sumfit=0;
			  for(vector<indiv>::iterator iter=subpop->begin(); iter != subpop->end(); iter++){
			  	if(iter->S_allele1() < 0) sc++;
				if(iter->S_allele2() < 0) sc++;
				if(iter->S_allele1() == 100) s100++;
				if(iter->S_allele2() == 100) s100++;
				sumfit += iter->fitness;
			  }
			  subpopsize = subpop->size();
			  p2outputSC[generation][i][j][0] = sc/(2*subpopsize);
		      p2outputSC[generation][i][j][1] = sumfit/subpopsize;
		      p2outputSC[generation][i][j][2] = s100/(2*subpopsize);
		      p2outputSC[generation][i][j][3] = inbreeding_depression_within_deme(r, subpop, selection_exp_mu, rec, dominance);
		      p2outputSC[generation][i][j][4] = p2outputSelfrate[generation][i][j];
		      delete subpop;
		  }
		}
	}
	
}



/// Create output of one simulation -  statistics for every  X generations (intervals) - to show spatial stats

void makeoutput_between_years_spatial(double ****stat, int stat_length, int patchesx, int patchesy, int intervals, ofstream &outfile)
{
	    outfile << "Deme\tInbreedingDep\tLoad\tEf-S\tFis-far\tFis-linked\tFis-glemin\tMate-avail\tDeleterious\tSelfcomp\n";
	    
	    // output write
	    for(int gen=0; gen<intervals; gen++){
			outfile << "\t" << gen << endl;
			for(int i=0; i<patchesx; i++){
				for(int j=0; j <patchesy; j++){
					outfile << i << "\t" << j << "\t";
			        for(int k=0; k<stat_length; k++){
				       outfile << format("%+05.4f") % stat[gen][i][j][k] << "\t";//setprecision(4) << stat[gen][i][j][k] << setfill(' ') << setw(10);
				    }
				    outfile << endl;
				}
			}
			outfile << endl;
		}
        
}


/// Create output of one simulation -  statistics for every  X generations - to show spatial stats of SC during expansion

void makeoutput_between_years_spatial_sc(double ****stat, int stat_length, int patchesx, int patchesy, int generations, ofstream &outfile)
{
	    outfile << "Deme\tInbreedingDep\tLoad\tEf-Allele S\tFis-far\tFis-linked\tFis-glemin\tMate-avail\tDeleterious\tSelfcomp\tID-2\tID-1\tID0\tID1\tID2";
	    
	    // output write
	    for(int gen=0; gen<generations; gen++){
			outfile << "\t" << gen << endl;
			for(int i=0; i<patchesx; i++){
				for(int j=0; j <patchesy; j++){
					outfile << i << "\t" << j << "\t";
			        for(int k=0; k<stat_length; k++){
				       outfile << format("%+05.4f") % stat[gen][i][j][k] << "\t";//setprecision(4) << stat[gen][i][j][k] << setfill(' ') << setw(10);
				    }
				    outfile << endl;
				}
			}
			outfile << endl;
		}
        
}


/// Create output of one simulation -  statistics for every  generations- to show spatial stats of S alleles before and during expansion

void makeoutput_between_years_spatial_s_alleles(double *****stat, int stat_length, int patchesx, int patchesy, int generations, vector<int> &Stype, int **space, ofstream &outfile)
{
	    outfile << "Deme\tS allele\tFreq S allele\tDelet mutations";
	    // add for loop s alleles
	    // output write
	    for(int gen=0; gen<generations; gen++){
			outfile << "Generation: " << gen << endl;
			for(int i=0; i<patchesx; i++){
				for(int j=0; j <patchesy; j++){
					if(space[i][j] == 1){
					   outfile << "Deme " << i << " " << j << endl;
					   for(unsigned int s_al=1; s_al<200; s_al++){
						  outfile << Stype[s_al] << "\t";
						  for(int k=0; k<stat_length; k++){
				              outfile << format("%+05.4f") % stat[gen][i][j][s_al][k] << "\t";
				          }
				          outfile << endl;
					   }
				       outfile << endl;
				    }
				}
			}
			outfile << endl;
		}
        
}


/// Create output of one simulation -  statistics for every  X generations (intervals) - to show the dynamics

void makeoutput_onesim(double **stat, int stat_length, int intervals, ofstream &outfile, double mig, double growth, double mu_selection, double mu_slocus, int generations, double rec, double selfing_prob, double h, double mutrates)
{	    
        //first line: column titles
        outfile << " ** Migration (m): " << mig << " " << "Growth: " << growth << " " << "Selection mean : " << mu_selection << "  " << "Mutation S locus: " << mu_slocus << "  " << "Generations: " << generations << " " << "Recombination " << rec << " " << "Selfing prob " << selfing_prob << "  " << " Dominance " << h << "  " << " Mutation rate sel" << mutrates << endl;
	    outfile << "Year\tp-freq\tHs\tHt\tFst\tHs_S\tHt_S\tGst_S\tHs_N\tHt_N\tMean W\tLoad\tInb-Dep\tInb-dBet\tAE-S\tGst_N_far"; // Freq of p, mean Fis, Hs, Ht, Fst, Sdiv, Mate availability
	    outfile << endl;
	    
	    // output write
	    for(int gen=0; gen<intervals; gen++){
			for(int k=0; k<stat_length; k++){
				if(k==0){
					outfile << (int) stat[gen][k] << "\t";
				}else{
					outfile << format("%+05.4f") % stat[gen][k] << "\t";
				}
			 }
			 outfile << endl;
		}

}


/// Create output of several simulations -  statistics for every simulation run at equilibrium

void makeoutput(double **stat, int stat_length, int simulations, ofstream &outfile, double mig, double growth, double mu_selection, double mu_slocus, int generations, double rec, double selfing_prob, double h, double mutrates)
{	    
        //first line: column titles
                //first line: column titles
        outfile << " Parameters:  Migration (m): " << mig << " " << "Growth: " << growth << " " << "Selection mean : " << mu_selection << "  " << "Mutation rate sel" << mutrates << "  Mutation S locus: " << mu_slocus << "  " << "Generations: " << generations << " " << "Recombination " << rec << "  " << "Selfing prob " << selfing_prob << "  " << "Dominance " << h << endl;
	    outfile << "Run\tp-freq\tHs\tHt\tFst\tHs_s\tHt_s\tGst_S\tHs_N\tHt_N\tMean W\tLoad\tID-within\tID-bet\tAE-S\tGst-far\tSC"; 
	    outfile << endl;
	    
	    // output write
	    for(int sim=0; sim<simulations; sim++){
			for(int k=0; k<stat_length; k++){
				if(k==0){
					outfile << (int) stat[sim][k] << "\t";
				}else{
					outfile << setprecision(4) << stat[sim][k] << "\t";
				}
			 }
			 outfile << endl;
		}

}




//////////////////////////////////////////////////////////////////

#define RUNS 1
#define GENS 401
#define EXPANSION 2001
#define INTERS 122
#define INTERVALS 5
#define STATS 19
#define STATSD 9
#define STATSC 5
#define STATS_A 2
#define DIMX 5
#define DIMY 50
#define S_AL 200

//////////////////////////////////////////////////////////////////

int main(int argc, char* argv[]) {
	srand (time(NULL));
	
	const gsl_rng_type * T;
    gsl_rng * r;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);
    gsl_rng_set(r, rand()%100);
    
    //Make output file
    string outfilename_year = argv[1];            
	//ofstream outfile_year (outfilename_year.c_str());
	string outfilename = argv[2];            
	ofstream outfile (outfilename.c_str());
    string outfilename_space = argv[3];            
	//ofstream outfilespace (outfilename_space.c_str());
	string outfilename_sc = argv[4];
    
	
	double mig = atof(argv[5]);
	double mu_selection = atof(argv[6]);
	double dominance = atof(argv[7]);
	double mutrates = atof(argv[8]);
	double mutraten = atof(argv[9]);
	double selfing_prob = atof(argv[10]);
	double rec = atof(argv[11]);
	double growth = 8.0;              
    
    int subpop= 100, patchesx=DIMX, patchesy=DIMY;
	int simulations = RUNS;    // number of simulations
	int generations = GENS;
	int interval = 100;
	int popsize_init = subpop*patchesx*patchesx;
	int nalleles =20;
	int genes = 500;
	int ngenes = genes/10;
	int S_alleles = S_AL;

	

	double lambda_neutral =  genes*mutraten;
	double lambda_sel = genes*mutrates;
	double mu_slocus = 0.00001;
	int stat_length=STATS;
	int stat_lengthD = STATSD;
	int stat_lengthSC = STATSC;
	int intervalsD = INTERVALS;
    int expansion = EXPANSION;

    clock_t t;
	cout << "Flag" << endl;

	// Vectors:
    //Vector of S alleles
    vector<int> Stype(S_alleles);  //vector of S alleles
    for( int s_al=0; s_al<S_alleles; s_al++){
		Stype[s_al] = 100 + s_al;
        //cout << "S type: " << Stype[s_al] << endl;
	}
	// Vector of neutral alleles
	vector<int> Ntype(nalleles);  //vector of neutral alleles
    for( int n_al=0; n_al<nalleles; n_al++){
		Ntype[n_al] = 30 + n_al;
        //cout << "N initial type: " << Ntype[n_al] << endl;
	}
    // Create dyn mem pointers to mult arrays for output
    double **p2stat_year_vec; /// pointers to stats per interval ---> single simulation
    double **p2output_sim;  /// pointers to stats per simulation
    int **p2space_occup;
    double ****p2output;  /// pointers to stats per interval and deme(x,y) ---> single simulation
    double ****p2outputSC;  /// pointers to SC freq and ID per generation during expansion and deme(x,y)
    double ***p2outputSelfrate;  /// pointers to Selfing rate per generation during expansion and deme(x,y) 
    

    
    // Allocating memeory for simulation statistics pointer
 	p2output_sim = new double*[RUNS];
	for(int run=0;  run<RUNS; ++run){
		p2output_sim[run] = new double[STATS];
	}
    
    
    /// ------------ %Simulation tuns % -------------- //
    for( int sim=0; sim<simulations; sim++){
		vector<indiv> *parents = new vector<indiv>(popsize_init); 
        
        
        //Allocate memory deme-interval pointers
        p2output = new double***[INTERVALS];
        for(int i=0;  i<INTERVALS; ++i){
			p2output[i] = new double**[DIMX];
			for(int j=0; j<DIMX; ++j){
				p2output[i][j] = new double*[DIMY];
				for(int l=0; l<DIMY; ++l)
				    p2output[i][j][l] = new double[STATSD];
			}
		}
        
        //Allocate memory deme-expansion SC pointers
        p2outputSC = new double***[EXPANSION];
        for(int i=0;  i<EXPANSION; ++i){
			p2outputSC[i] = new double**[DIMX];
			for(int j=0; j<DIMX; ++j){
				p2outputSC[i][j] = new double*[DIMY];
				for(int l=0; l<DIMY; ++l)
				    p2outputSC[i][j][l] = new double[STATSC];
			}
		}

		
		//Allocate memory deme-expansion SC pointers
        p2outputSelfrate = new double**[EXPANSION];
        for(int i=0;  i<EXPANSION; ++i){
			p2outputSelfrate[i] = new double*[DIMX];
			for(int j=0; j<DIMX; ++j){
				p2outputSelfrate[i][j] = new double[DIMY];
			}
		}
		
		//Allocate memory interval pointers
		p2stat_year_vec = new double*[INTERS];
		for(int i=0;  i<INTERS; ++i){
			p2stat_year_vec[i] = new double[STATS];
		}
 

    // Allocate memory space occupancy  --- Fill spatial occupancy matrix with zeroes
         p2space_occup = new int*[DIMX];
	     for(int i=0; i<DIMX; ++i)
	          p2space_occup[i] = new int[DIMY];
	     for(int i=0; i<patchesx; ++i){
			 for(int j=0; j<patchesy; ++j){
				 p2space_occup[i][j]=0;
		     }
	     }



		// Initialize individuals
        initializeplants(r, parents, genes, nalleles, S_alleles, subpop,  DIMX, DIMX, mu_selection, dominance, p2space_occup);
        show_vector((*parents)[0].genome);
        int inter=0, expa=0, t1, flag=0;
        t = clock();
	// ------------  % Generation loop % ----------- //
	    for(int gen=0; gen<generations ; gen++){
			cout << "Generation: " << gen << " of " << generations << endl;
			// Allocate memory
		    vector<indiv> *seeds = new vector<indiv>;
            vector<vector<int> > *pt_neilist = new vector<vector<int> >;
            vector<vector<int> > *pt_listexp = new vector<vector<int> >;  
            vector<int> *pt_chosen = new vector<int>; 

			int popsize = parents->size();
			cout << "Popsize : " << popsize << endl;
			
			/// 1) Gamete production (recombination, mutation) and local mating
			t1 = clock();
			gameteprod_mating(r, parents, seeds, genes, patchesx, patchesy, subpop, growth, lambda_neutral, lambda_sel, mu_slocus, mu_selection, dominance, rec, selfing_prob, p2space_occup, Ntype, Stype, flag, p2outputSelfrate, expa);
			cout << "Time gametoprod: " << float(clock() - t1) << endl;
			
			if(gen>0){
			/// 3) Seed dispersal and colonization:
			    refresh_chosen(seeds);
			    t1 = clock();
			    migration_stepping_stone(r,seeds, mig, patchesx, patchesy,  p2space_occup);
			    cout << "Time migration: " << float(clock() - t1) << endl;
			    //cout << "migration: " << mig << endl;
			    refresh_chosen(seeds);
			}
			/// Introduction of null mutation at the S locus (introduction of Self-Compatibility)
			if (gen > 1000){
				flag = 1; // --> 0 means no SC invasion
				//cout << "Invasion SC starts: " << mu_slocus << endl;
			}
			cout << "S locus mutation rate " << mu_slocus << endl;
			
			
			/// Range expansion
			if (gen > 2000){
				t1 = clock();
				colonization(r, seeds, pt_listexp, subpop, genes, mig, patchesx, patchesy, p2space_occup, flag, expa);
				//cout << "Time colonization: " << float(clock() - t1) << endl;
				self_compatibles_output(r,expa, sim, patchesx, patchesy, p2outputSC, p2outputSelfrate, seeds, p2space_occup,mu_selection, rec, dominance);
				//cout << "FLAG" << endl;
			    expa++;
			}
			show_sc(seeds);

		    /// Calculate statistics:
	        if(gen%interval == 0){
				 allele_frequencies_2(r,seeds, p2stat_year_vec, inter, patchesx, patchesy, genes, ngenes, S_alleles, Stype, Ntype, gen, p2space_occup, mu_selection, rec, dominance, p2output, gen, mu_selection);
				 //cout << "Generation: " << gen << endl;
				 //cout << "Size: " << seeds->size() << " Selection " << mu_selection << " Migration: " << mig << " Dominance " << dominance << endl;
				 inter++;
		    }


            /// Break when expansion reaches the boundary
			//if(check_expansion(seeds,patchesx,p2space_occup) == true)  break;
			
			parents->swap(*seeds);
			
			// Release memory
	        delete seeds;
	        delete pt_neilist;
            delete pt_listexp;
            delete pt_chosen;
           
		} /// END OF GENERATION LOOP
		cout << "It took one generation: " << clock() - t << endl; 
		delete parents;
		
       /// Output spatial stats 
		stringstream extension_file;//create a stringstream
		extension_file << sim;//add number to the stream
		string outfilespacefull1 = outfilename_space + extension_file.str();
		string outfilespacefull1_year = outfilename_year + extension_file.str();
		string outfilespacefull1_sc= outfilename_sc + extension_file.str();
		
		ofstream outfilespacefull2 (outfilespacefull1.c_str());
		ofstream outfilespacefull2_year (outfilespacefull1_year.c_str());
		ofstream outfilespacefull2_sc (outfilespacefull1_sc.c_str());

		makeoutput_between_years_spatial(p2output, stat_lengthD, patchesx, patchesy, intervalsD, outfilespacefull2);
		makeoutput_onesim(p2stat_year_vec, stat_length, inter, outfilespacefull2_year, mig, growth, mu_selection, mu_slocus, generations, rec, selfing_prob, dominance);  
		makeoutput_between_years_spatial_sc(p2outputSC, stat_lengthSC, patchesx, patchesy, expansion, outfilespacefull2_sc);

		
		// Put final year statistics (equilibrium) at double pointer p2output_sim
		final_year_stat(p2stat_year_vec, p2output_sim, sim, interval, generations);
		
		/// Deallocating memory:
        // Deallocating memory pointers stats-interval
        for(int i=0; i<INTERS; i++)
			delete [] p2stat_year_vec[i];
		delete [] p2stat_year_vec;
		
		//Deallocate memory pointers space-occupancy (DIMX-DIMY)
        for(int i=0; i<DIMX; i++)
            delete [] p2space_occup[i];
        delete [] p2space_occup;
        

        
        // Deallocating memory pointers stats-deme-interval
		for (int i = 0; i < INTERVALS; ++i){
			for (int j = 0; j < DIMX; ++j){
				for(int l=0; l <DIMY; ++l)
				     delete [] p2output[i][j][l];
				delete [] p2output[i][j];
		     }
			 delete [] p2output[i];
		}	    
        delete [] p2output;
        
        // Deallocating memory pointers generations(expansion)-demes
		for (int i=0; i < EXPANSION; ++i){
			for (int j=0; j < DIMX; ++j){
				for(int l=0; l < DIMY; ++l)
				     delete [] p2outputSC[i][j][l];
				delete [] p2outputSC[i][j];
		     }
			 delete [] p2outputSC[i];
		}	    
        delete [] p2outputSC;
        
        // Deallocating memory pointers (p2outputSelfrate) generations(expansion)-demes
		for (int i=0; i < EXPANSION; ++i){
			for (int j=0; j < DIMX; ++j){
				delete [] p2outputSelfrate[i][j];
		     }
			 delete [] p2outputSelfrate[i];
		}	    
        delete [] p2outputSelfrate;
        
      
	}/// END OF SIMULATION LOOP

    // Write output for all simulations
    makeoutput(p2output_sim, stat_length, simulations, outfile, mig, growth, mu_selection, mu_slocus, generations, rec, selfing_prob, dominance, mutrates);
    
    // Deallocate memory pointers stats-simulation
    for(int run=0; run<RUNS; run++)
		delete [] p2output_sim[run];
	delete [] p2output_sim;



	gsl_rng_free (r);
	return 0;

}
