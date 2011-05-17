/**s1 -> partial_external
*  s2 -> partial_multi
*  s3 -> partial_multi2
*  u1 -> u_helper
*  f -> cond_dangle (conditional)
**/

/**TODO cond_dangle, ed3, ed5,**/

/**
*  For each spot in partial_external at index h,j we'll need:
**/

double RT = 0.00198721 * 310.15;

int l;

for(l = i+2; l < j; l++){
	upm[i][j] += up[i+1][l] * exp(-(a + 2 * c + /**TODO terminal pair i+1,l**/)/RT) *
				(exp(-ed3(i+1,l,l+1)/RT) * 
				u_helper[l+2][j-1] + 
				u_helper[l+1][j-1] -
				u_helper[l+2][j-1]);
	if(l != i+2){ //goes from i + 2 < l < j
		upm[i][j] += up[i + 2][l] * 
					exp(-(a + 2c + b + ed3(j,i,i + 1) + /**TODO terminal pair i+2, l**/)/RT) *
					(exp(-ed3(i+2,l,l+1)/RT)*u_helper[l+2][j-1] + u_helper[l+1][j-1]  - 
					u_helper[l+2][j-1]);

		if(l != j - 1){
			upm[i][j] += exp(-ed3(j,i,i+1)/RT)*exp(-(a + 2 * c + (h - i -1) * b)/RT) * 
						partial_multi[l][j]; 
			//Changed h in the paper to l it makes no difference
			//other than unifying loops
		}
	}
}


for(l = h+1; l < j; l++){
	partial_external[h][j] += up[h][l] * exp( -(ed5(h,l,h-1) + /**TODO terminal pair hl**/) / RT) *
				(exp(-ed3(h,l,l+1)/RT)*u_helper[l+2,j-1] + u_helper[l+1][j] -
				u_helper[l+2][j]);

	partial_multi[h][j] += up[h][l] * exp( -(ed5(h,l,h-1) + /**TODO terminal pair hl**/) / RT) *
				(exp(-ed3(h,l,l+1)/RT)*u_helper[l+2,j-1] + u_helper[l+1][j-1] -
				u_helper[l+2][j-1]);

	partial_multi2[h][j] += up[h][l] * exp( -(ed5(h,l,h-1) + /**TODO terminal pair hl**/) / RT) *
				(cond_dangle(j + 1, h, l)*exp(-(j-l) * b / RT) + exp(-ed3(h,l,l + 1)/RT)*u_helper(l+2,j)
				+u_helper(l + 1,j) - u1(l + 2,j));
				
}
/** partial_multi2 goes up to j **/
partial_multi2[h][j] += up[h][l] * exp( -(ed5(h,l,h-1) + /**TODO terminal pair hl**/) / RT) *
				(cond_dangle(j + 1, h, l)*exp(-(j-l) * b / RT) + exp(-ed3(h,l,l + 1)/RT)*u_helper(l+2,j)
				+u_helper(l + 1,j) - u1(l + 3,j));
