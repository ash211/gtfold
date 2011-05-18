//The paper (Dirks, Pierce) offhandedly mentions a G-U terminal pair penalty
//I don't think we do that...

/**s1 -> partial_external
*  s2 -> partial_multi
*  s3 -> partial_multi2
*  u1 -> u_multi
*  f -> cond_dangle (conditional)
**/

#include energy.h


int cond_dangle(int j, int h, int l){
	return j - 1 == l ? 1 : exp(-Ed3(h,l,l+1)/RT);
}

double RT = (0.00198721 * 310.15)/100.00; 
//So that we don't have the /100 in every exponent
//Seeing as we missed some last time

int seg_length;
int i;
for(seg_length = MIN_TURN; seg_length < len; len ++){
	for(i = 1; i < seg_length - MIN_TURN; i++){
		j = i + seg_length;
		int l;
		/**
		*  For each spot in partial_external at index h,j we'll need:
		**/
		for(l = i+2; l < j; l++){
			upm[i][j] += up[i+1][l] * exp(-(Ea + 2 * Ec + auPenalty(i + 1,l)/RT)) *
						(exp(-Ed3(i + 1,l,l+1)/RT) * u_multi[l + 2][j - 1] + 
						u_multi[l + 1][j - 1] -
						u_multi[l + 2][j - 1]);
			if(l != i+2){ //goes from i + 2 < l < j
				upm[i][j] += up[i + 2][l] * 
						exp(-(Ea + 2 * Ec + Eb + Ed3(j,i,i + 1) + auPenalty(i + 2, l))/RT) *
						(exp(-Ed3(i + 2,l,l + 1)/RT)*u_multi[l + 2][j - 1] + 
						u_multi[l + 1][j - 1]  - 
						u_multi[l + 2][j - 1]);

				if(l != j - 1){
					upm[i][j] += exp(-Ed3(j,i,i + 1)/RT)*exp(-(Ea + 2 * Ec + (l - i - 1) * Eb)/RT) * 
								partial_multi[l][j]; 
					//Changed h in the paper to l it makes no difference
					//other than unifying loops
				}
			}
		}


		up[i][j] += exp(-eh(i,j)/RT) + exp(-es(i,j)/RT) * up[i + 1][j - 1] + upm[i][j];
		for(l = j - 1; l > i + 1; l--){
			if((j - l) - 2> MAX_LOOP){
				break;
			}
			for(h = i + 1 ; h < l; h++){
				if(!(i == h - 1 && j == l - 1)){ //If this is true we have a stack
					up[i][j] += exp(-ebi(i,j,h,l)/RT);
				}
				if((j - l) + (h - i) - 2> MAX_LOOP){
					break;
				}
			}
		}

		for(l = i + 1; l <= j ; l++){
			u_multi[i][j] += up[i][l] * exp(-(Ec + auPenalty(i,l))/RT) *
					(cond_dangle(j + 1, i, l) * exp(-(j - l) * Eb/RT) + 
					exp(-Ed3(i,l,l + 1)/RT) * u_multi[l + 2][j] + 
					u_multi[l + 1][j] -
					u_multi[l + 2][j]);
			if(l != i+1){	
				u_multi[i][j] += up[i + 1][j]*exp(-(Ec + Eb + auPenalty(i + 1,l))/RT)*
						(cond_dangle(j + 1, i + 1, l) * exp(-(j - l)*Eb/RT) + 
						exp(-Ed3(i + 1,l,l + 1)/RT) * u_multi[l + 2][j] + 
						u_multi[l + 1][j] - 
						u_multi[l + 2][j]);
				if(l != j){
					u_multi[i][j] += exp(-(Ec + (h - i) * Eb)/RT) * partial_multi2[l][j];

				}
			}
		}


		int h = i; //To stay consistant with the notation in the paper.

		for(l = h+1; l < j; l++){
			partial_external[h][j] += up[h][l] * exp( -(Ed5(h,l,h - 1) + auPenalty(h,l)) / RT) *
						(exp(-Ed3(h,l,l + 1)/RT)*u[l + 2,j - 1] + 
						u[l + 1][j] -
						u[l + 2][j]);

			partial_multi[h][j] += up[h][l] * exp( -(Ed5(h,l,h - 1) + auPenalty(h,l)) / RT) *
						(exp(-Ed3(h,l,l + 1)/RT)*u_multi[l + 2][j - 1] + 
						u_multi[l + 1][j - 1] -
						u_multi[l + 2][j - 1]);

			partial_multi2[h][j] += up[h][l] * exp( -(Ed5(h,l,h - 1) + auPenalty(h,l)) / RT) *
						(cond_dangle(j + 1, h, l)*exp(-(j - l) * Eb / RT) + 
						exp(-Ed3(h,l,l + 1)/RT)*u_multi[l + 2][j] + 
						u_multi[l + 1][j] - 
						u_multi[l + 2][j]);
						
		}
		/** partial_multi2 goes up to j **/
		partial_multi2[h][j] += up[h][l] * exp( -(Ed5(h,l,h - 1) + auPenalty(h,l)) / RT) *
					(cond_dangle(j + 1, h, l)*exp(-(j-l) * Eb / RT) + 
					exp(-Ed3(h,l,l + 1)/RT) * u_multi[l + 2][j] + 
					u_multi[l + 1][j] - 
					u_multi[l + 2][j]);

		//Finally we do the u matrix
		u[i][j] = 1 + up[i][j]*exp(-auPenalty(i,j) /RT);

		for(l = i + 1; l < j; l++){
			//Replaced l with h in the next line to unify loops
			u[i][j] += up[l][j] * exp(-(Ed5(l, j , l - 1) + auPenalty(h,j))/RT)

			u[i][j] += exp(-(auPenalty(i,l))/RT)*(exp(-Ed3(i,l,l + 1)/RT) * u[l + 2][j] +
				   u[l + 1][j] - u[l + 2][j]);
			if(l != j - 1) {
				u[i][j] += partial_external[l][j];
			}
				   
		}
	}
}
