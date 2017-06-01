#include <stdio.h>
#include <stdint.h>
typedef struct {
	uint16_t rnd1 : 11;
	uint8_t sign1 : 1;
	uint16_t rnd2 : 11;
	uint8_t sign2 : 1;
} __attribute__((__packed__)) three_bytes_packed;
const size_t CDF_LENGTH_D3 = 6;
const uint16_t CDF_D3[6] = {602, 1521, 1927, 2031, 2046, 2047}; // out of [0, 2047]
const uint16_t b_bar=11;
//752*8
//n=752*8=6016
void lwe_sample_n_inverse_12_alice(uint16_t s[752][8],uint16_t rnd[9026]) {
	uint16_t n =6016;
	size_t rndlen =3 * ((n + 1) / 2);  // 12 bits of unif randomness per output element


	size_t i;
	//n=752*8=6016
	
	
	for (i = 0; i < n; i += 2) {  // two output elements at a time
		three_bytes_packed *ptr_packed = (three_bytes_packed *)(rnd + 3 * i / 2);
		
		uint16_t rnd1 = ptr_packed->rnd1;
		uint16_t rnd2 = ptr_packed->rnd2;

		uint8_t sample1 = 0;
		uint8_t sample2 = 0;

		size_t j;
		// No need to compare with the last value.
		for (j = 0; j < CDF_LENGTH_D3 - 1; j++) {
			// Constant time comparison: 1 if LWE_CDF_TABLE[j] < rnd1, 0 otherwise.
			// Critically uses the fact that LWE_CDF_TABLE[j] and rnd1 fit in 15 bits.
			sample1 += (uint16_t)(CDF_D3[j] - rnd1) >> 15;
			sample2 += (uint16_t)(CDF_D3[j] - rnd2) >> 15;
		}

		uint8_t sign1 = ptr_packed->sign1;
		uint8_t sign2 = ptr_packed->sign2;

		// Assuming that sign1 is either 0 or 1, flips sample1 iff sign1 = 1

		s[i/8][i%8] = ((-sign1) ^ sample1) + sign1;

		if (i + 1 < n) {
			s[(i+1)/8][(i+1)%8] = ((-sign2) ^ sample2) + sign2;
		}
	}
}


void lwe_sample_n_inverse_12_bob(uint16_t s[8][752],uint16_t rnd[9026]) {
	uint16_t n =6016;
	size_t rndlen =3 * ((n + 1) / 2);  // 12 bits of unif randomness per output element


	size_t i;
	//n=752*8=6016
	
	
	for (i = 0; i < n; i += 2) {  // two output elements at a time
		three_bytes_packed *ptr_packed = (three_bytes_packed *)(rnd + 3 * i / 2);
		
		uint16_t rnd1 = ptr_packed->rnd1;
		uint16_t rnd2 = ptr_packed->rnd2;

		uint8_t sample1 = 0;
		uint8_t sample2 = 0;

		size_t j;
		// No need to compare with the last value.
		for (j = 0; j < CDF_LENGTH_D3 - 1; j++) {
			// Constant time comparison: 1 if LWE_CDF_TABLE[j] < rnd1, 0 otherwise.
			// Critically uses the fact that LWE_CDF_TABLE[j] and rnd1 fit in 15 bits.
			sample1 += (uint16_t)(CDF_D3[j] - rnd1) >> 15;
			sample2 += (uint16_t)(CDF_D3[j] - rnd2) >> 15;
		}

		uint8_t sign1 = ptr_packed->sign1;
		uint8_t sign2 = ptr_packed->sign2;

		// Assuming that sign1 is either 0 or 1, flips sample1 iff sign1 = 1

		s[i%8][i/8] = ((-sign1) ^ sample1) + sign1;

		if (i + 1 < n) {
			s[(i+1)%8] [(i+1)/8]= ((-sign2) ^ sample2) + sign2;
		}
	}
}

void lwe_sample_n_inverse_12_64(uint16_t s[8][8],uint16_t rnd[98]) {
	uint16_t n =64;
	size_t rndlen =3 * ((n + 1) / 2);  // 12 bits of unif randomness per output element


	size_t i;
	//n=752*8=6016
	
	
	for (i = 0; i < n; i += 2) {  // two output elements at a time
		three_bytes_packed *ptr_packed = (three_bytes_packed *)(rnd + 3 * i / 2);
		
		uint16_t rnd1 = ptr_packed->rnd1;
		uint16_t rnd2 = ptr_packed->rnd2;

		uint8_t sample1 = 0;
		uint8_t sample2 = 0;

		size_t j;
		// No need to compare with the last value.
		for (j = 0; j < CDF_LENGTH_D3 - 1; j++) {
			// Constant time comparison: 1 if LWE_CDF_TABLE[j] < rnd1, 0 otherwise.
			// Critically uses the fact that LWE_CDF_TABLE[j] and rnd1 fit in 15 bits.
			sample1 += (uint16_t)(CDF_D3[j] - rnd1) >> 15;
			sample2 += (uint16_t)(CDF_D3[j] - rnd2) >> 15;
		}

		uint8_t sign1 = ptr_packed->sign1;
		uint8_t sign2 = ptr_packed->sign2;

		// Assuming that sign1 is either 0 or 1, flips sample1 iff sign1 = 1

		s[i%8][i/8] = ((-sign1) ^ sample1) + sign1;

		if (i + 1 < n) {
			s[(i+1)%8] [(i+1)/8]= ((-sign2) ^ sample2) + sign2;
		}
	}
}

uint16_t rounding (uint16_t input){
	return ((input+1024)%32768)>>b_bar;
}

uint16_t cross (uint16_t input){
	return ((input%32768)>>(b_bar+1))%2;
}

//outputs 1 if a>b, 0 else
uint16_t greater_than(uint16_t a, uint16_t b){
	uint16_t changed = 0;
	uint16_t output = 0;

	int i;
	uint16_t pow2[16]={32768, 16384, 8192, 4096, 2048, 1024, 512,  256, 128, 64, 32, 16, 8, 4, 2, 1};
	for (i=0; i<16; i++){
		uint16_t bita = (a&pow2[i])>>(15-i);//(a<<(i))>>15;
		uint16_t bitb = (b&pow2[i])>>(15-i);//(b<<(i))>>15;
		//printf("should be getting binary %u %u \n",bita, bitb);
		output = (changed ^ 0)*output+ (changed^1)*bita*(bita ^ bitb);
		changed = (changed ^ 0) * changed + (changed ^ 1) * ( changed ^ (bita ^ bitb));
	}
	return output;

}

uint16_t rec(uint16_t w, uint16_t b){
	uint16_t equal = (cross(w) & b) + ((1^cross(w)) & (b^1));
	uint16_t output = equal*w;

	uint16_t want_zero=greater_than(cross(w),b);
	uint16_t output1=(want_zero ^ 1) * (greater_than(w,2047)*4096 +(greater_than(w,2047) ^1)*32767 );
	output1+=(want_zero ) * (greater_than(w,18433)*0 +(greater_than(w,18433) ^1)*4095 );
	output=(equal ^1)*output1;
	return output;

}

//#include <math.h>
void main3(){
	int i,j;
	int counter =0;
	/*
	for (i=0; i<32768; i++){
		for (j=0; j<32768; j++){
			uint16_t result=0;
			if (i>j) result=1;
			if (greater_than(i,j) != result) {
				counter ++;
				break;
			}
		}
		//if (rounding(i) != ((int)(floor(i/2048+0.5)))) counter ++;
	}
	*/
	//printf ("we have that %i\n",counter);
	//printf ("we have that for greater %u \n",greater_than(1,5));
	//printf ("we have that %u %u %u %u %u %u \n",cross(2047),cross(2048),cross(2049), cross(4095), cross(4096), cross(4097));

}
void main2(){
	int i;
	int counter =0;
	for (i=0; i<32768; i++){
		//if (rounding(i) != ((int)(floor(i/2048+0.5)))) counter ++;
		//if (cross(i) != ((int)(floor(i/4096)))) counter ++;
	}
	printf ("we have that %i\n",counter);

	printf ("we have that %u %u %u %u %u %u \n",cross(2047),cross(2048),cross(2049), cross(4095), cross(4096), cross(4097));

}

void main(){
	uint16_t A[752][752];
	uint16_t q=32768;
	int i,j;

	//This is what Alice does
	////////////////////////////////////////////////////////////

	for ( i = 0; i < 752; i++ ){
		for ( j = 0; j < 752; j++ ){
			A[i][j]=(i+j)%32768;

		}

	}
	uint16_t rnd1[9026];
	uint16_t rnd2[9026];
	for ( i = 0; i < 9026; i++ ){
		rnd1[i]=i;
		rnd2[i]=(i*5303)%32768;
	}

	uint16_t s[752][8];
	uint16_t e[752][8];

	lwe_sample_n_inverse_12_alice(s,rnd1);
	lwe_sample_n_inverse_12_alice(e,rnd2);


	uint16_t b[752][8];
	for ( i = 0; i < 752; i++ ){
		for ( j = 0; j < 8; j++ ){
			b[i][j]=e[i][j];

		}

	}

	int i1,j1;
	for ( i = 0; i < 752; i++ ){
		for ( j = 0; j < 8; j++ ){


			for ( i1 = 0; i1 < 752; i1++ ){
				for ( j1 = 0; j1 < 752; j1++ ){
					b[i][j]=b[i][j]+A[i][i1]*s[j1][j];
				}
			}


		}

	}
	////////////////////////////////////////////////////////////
	//This is what Bob does
	uint16_t Ab[752][752];
	for ( i = 0; i < 752; i++ ){
		for ( j = 0; j < 752; j++ ){
			Ab[i][j]=(i+j)%32768;

		}

	}
	uint16_t rnd3[9026];
	uint16_t rnd4[9026];
	for ( i = 0; i < 9026; i++ ){
		rnd3[i]=(i*1103	)%32768;
		rnd4[i]=(i*7717)%32768;
	}

	uint16_t sp[8][752];
	uint16_t ep[8][752];

	lwe_sample_n_inverse_12_bob(sp,rnd3);
	lwe_sample_n_inverse_12_bob(ep,rnd4);


	uint16_t bp[8][752];
	for ( i = 0; i < 8; i++ ){
		for ( j = 0; j < 752; j++ ){
			bp[i][j]=ep[i][j];

		}

	}

	for ( i = 0; i < 8; i++ ){
		for ( j = 0; j < 8; j++ ){


			for ( i1 = 0; i1 < 752; i1++ ){
				for ( j1 = 0; j1 < 752; j1++ ){
					bp[i][j]=bp[i][i1]+sp[i][i1]*b[j1][j];
				}
			}


		}

	}
	uint16_t epp[8][8];
	uint16_t rnd5[98];
	for ( i = 0; i < 98; i++ ){
		rnd5[i]=(i*5503)%32768;
	}
	lwe_sample_n_inverse_12_bob(epp,rnd5);
	uint16_t v[8][8];
	for ( i = 0; i < 8; i++ ){
		for ( j = 0; j < 8; j++ ){
			v[i][j]=epp[i][j];
		}
	}
	for ( i = 0; i < 8; i++ ){
		for ( j = 0; j < 8; j++ ){


			for ( i1 = 0; i1 < 752; i1++ ){
				for ( j1 = 0; j1 < 752; j1++ ){
					v[i][j]=v[i][j]+sp[i][i1]*Ab[j1][j];
				}
			}


		}

	}


	uint16_t c[8][8];
	for ( i = 0; i < 8; i++ ){
		for ( j = 0; j < 8; j++ ){
			c[i][j]=cross(v[i][j]);
		}
	}




	uint16_t k_bob[8][8];
	for ( i = 0; i < 8; i++ ){
		for ( j = 0; j < 8; j++ ){
			k_bob[i][j]=rounding(v[i][j]);
		}
	}

	//This is waht alice does
	/////////////////////////////
	uint16_t bp_s[8][8];
	for ( i = 0; i < 8; i++ ){
		for ( j = 0; j < 8; j++ ){
			bp_s[i][j]=0;
		}
	}
	for ( i = 0; i < 8; i++ ){
		for ( j = 0; j < 8; j++ ){


			for ( i1 = 0; i1 < 752; i1++ ){
				for ( j1 = 0; j1 < 752; j1++ ){
					bp_s[i][j]=bp_s[i][j]+bp[i][i1]*s[j1][j];
				}
			}


		}

	}


	

	uint16_t k_alice[8][8];
	for ( i = 0; i < 8; i++ ){
		for ( j = 0; j < 8; j++ ){
			bp_s[i][j]=bp_s[i][j]%q;
			k_alice[i][j]=rec(bp_s[i][j],c[i][j]);
		}
	}
	printf("Done...................\n");

}
