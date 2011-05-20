#include <math.h>
#include <stdlib.h>

#include "parition-dangle.h"
#include "energy.h"

double RT=(0.00198721 * 310.15)/100.00; 

double randdouble(){
	    return rand()/(double(RAND_MAX)+1);
} 

//P_0
double single_stranded_prob(int i, int j, dangle_struct d_struct){
	return 1 / d_struct.u[i][j];
}

//P_{ij}
double ends_pair_prob(int i, int j, dangle_struct d_struct){
	return d_struct.up[i][j] * exp(-auPenalty(i,j)/RT)/d_struct.u[i][j];
}

//P_{hj}
double prob_h_pairs_with_5_end
		(int i, int h, int j, dangle_struct d_struct){

	return d_struct.up[h][j] * 
			exp(-(ed5(h,j,h - 1) + 
			auPenalty(h,j))/RT) /
			d_struct.u[i][j];
}

//P_{il}
double prob_l_pairs_with_3_end
	(int i, int l, int j dangle_struct d_struct){

	return d_struct.up[i][l] * exp(-auPenalty(i,l)/RT) *
			(exp(-Ed3(i,l,l+1)/RT) * d_struct.u[l + 2][j] +
			d_struct.u[l + 1][j] -
			d_struct.u[l + 2][j]) /
			d_struct.u[i][j];
}

//P_{s1h}
double prob_h_pairs_in_ij(int i, int h, int j, dangle_struct d_struct){
	return d_struct.s1[h][j] / d_struct.u[i][j];
}

//P_{hl}
//This is the probability h pairs with l 
//**Given that h pairs before j**
double prob_h_pairs_with_l(int h, int l, int j, dangle_struct d_struct){
	return d_struct.up[h][l] * exp(-(Ed5(h,l,h - 1) + auPenalty(h,l))/RT) *
			(exp(-Ed3(h,l,l + 1)/RT) * d_struct.u[l + 2][j] +
			d_struct.u[l + 1][j] -
			d_struct.u[l + 2][j]) /
			d_struct.s1[h][j];
}

//Q_{ijH}
double prob_ij_pair_is_hairpin(int i, int j, dangle_struct d_struct){
	return exp(-eH(i,j)/RT)/d_struct.up[i][j];
}

//Q_{ijS}
double prob_ij_pair_is_stack(int i, int j, dangle_struct d_struct){
	return exp(-eS(i,j)/RT) * d_struct.up[i + 1][j - 1]/d_struct.u[i][j];
}
//To find Q_{ijBl} take the rest of the probabilities and subtract
//them from 1.

//Q_{ijM}
double prob_ij_pair_is_multiloop(int i, int j, dangle_struct d_struct){
	return d_struct.upm[i][j]/d_struct.up[i][j];
}

//This the probability hl closes a multiloop
//**GIVEN IJ CLOSES A MULTILOOP**
double prob_ijhl_is_internal_loop
	(int i, int j, int h, int l, double partition_sum, dangle_struct d_struct){
	return exp(-eL(i,j,h,l)/RT) * d_struct.up[h][l] / partition_sum;
}

//P_{ij(1+1)l}
double first_multi_prob_pair_next_base
	(int i, int l, int j, dangle_struct d_struct){
	return d_struct.up[i + 1][l] * 
			(exp( - (Ea + 2 * Ec + auPenalty(i + 1, l))/RT) *
			(exp(-Ed3(i + 1,l,l +1))/RT) * d_struct.u1[l + 2][j - 1] +
			d_struct.u1[l + 1][j - 1] -
			d_struct.u1[l + 1][j - 1]) /
			d_struct.upm[i][j];
}

//P_{ij(1+2)l}
double first_multi_prob_pair_after_base
	(int i, int l, int j, dangle_struct d_struct){
	
	return d_struct[i + 2][l] * 
			exp(-(Ea + 2 * Ec + Eb + Ed3(j, i, i + 1) + auPenalty(i + 2,l))/RT) *
			(exp(-Ed3(i + 2,l,l+1)/RT) * d_struct.u1[l + 2][j - 1] +
			d_struct.u1[l + 1][j - 1] -
			d_struct.u1[l + 2][j - 1]) /
			d_struct.upm[i][j];
}

//P_{ijs2h)
double first_multi_prob_pair_rest_base(int i, int h, int j, dangle_struct d_struct){
	return exp(-(Ea + 2 * Ec + (h - i - 1) * b + Ed3(j,i,i+1))/RT) *
			d_struct.s2[h][j] / d_struct.upm[i][j];
}
//P_{ijhl}
double first_multi_prob_base_rest
	(int i, int h, int l, int j, dangle_struct d_struct){
	
	return d_struct.up[h][l] * exp(-(Ed5(h,l,h-1) + auPenalty(h,l))/RT) * 
			(exp( -Ed3(h,l,l+1)/RT) * d_struct.u1[l + 2][j - 1] +
			d_struct.u1[l + 1][j - 1] -
			d_struct.u1[l + 2][j - 1])/d_struct.s2[h][j];
}

//Q_{(l1 + 1)(j -1)(l1 +1)}
double rest_multi_prob_next_base_pairs
	(int l1, int l, int j, dangle_struct d_struct){
	
	return d_struct.up[l1 +1][l] * exp( - (Ec + auPenalty(l1 + 1, l))/RT)
			(cond_dangle(j, l1 + 1, l) * exp( -(j - 1 - l) * Eb /RT) +
			exp(-(j - 1 - l) * Eb) / RT + 
			exp(-Ed3(l1 + 1, l, l + 1)/RT) * d_struct.u1[l + 2][j - 1] +
			d_struct.u1[l + 1][j - 1] -
			d_struct.u1[l + 1][j - 1]) /
			d_struct.u1[l1 + 1][j - 1];
}

//Q_{(l1 + 1)(j - 1)(l1 + 2)l}
double rest_multi_prob_after_base_pairs
	(int l1, int l, int j, dangle_struct d_struct){
	
	return d_struct.up[l1 + 2][l] * exp(-(Ec + Eb + auPenalty(l1 + 2,l))/RT) *
			(cond_dangle(j, l1 + 2, l) * exp( -(j - 1 - l) * Eb / RT) +
			exp(-Ed3(l1 + 1,l ,l + 1)/RT) * d_struct.u1[l + 2][j - 1] +
			d_struct.u1[l + 1][j - 1] - 
			d_struct.u1[l + 2][j - 1]) /
			d_struct.u1[l1 + 1][j - 1];
			
}

//Q_{(l1 + 1)(j - 1)s3h}
double rest_multi_prob_rest_base_pairs
	(int l1, int h, int j, dangle_struct d_struct){
	
	return exp( -(Ec + (h - l1 - 1) * Eb)/RT) * d_struct.s3[h][j - 1] /
			d_struct.u1[l1 + 1][j - 1];
}

//Q_{(j-1)hl}
//GIVEN THAT H PAIRS BEFORE J
double rest_multi_prob_h_pairs_with_l
	(int l1, int h, int l, int j, dangle_struct d_struct){
	
	return d_struct.up[h][l] * exp(-(Ed5(h,l,h - 1) + auPenalty(h,l))/RT) *
			(cond_dangle(j,h,l) * exp(-(j - 1 -l) * Eb/RT) + 
			exp(-Ed3(h,l,l + 1)/RT) * d_struct.u1[l + 2][j - 1] +
			d_struct.u1[l + 1][j - 1] -
			d_struct.u1[l + 2][j - 1]) /
			d_struct.s3[h][j - 1];
}

int there_is_another_helix(int h2, int l2, int j, dangle_struct d_struct){
	double counter = randdouble();
	double coin-flip = cond_dangle(j,h2,l2) * exp(-(j - 1 - l2) * Eb / RT) /
					  (cond_dangle(j,h2,l2) * exp(-(j - 1 - l2) * Eb / RT) + 
					exp(Ed3(h2,l2,l2 + 1)/RT) * d_struct.u1[l2 + 2][j - 1] +
					d_struct.u1[l2 + 1][j - 1] -
					d_struct.u1[l2 + 2][j - 1]);
	
	//If I was slick and mean I would have put
		
	//return counter > coin-flip

	if(counter <= coin-flip){
		return 0;
	}
	else{
		return 1;
	}
}
