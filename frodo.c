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
	return (((input+1024)%32768)>>b_bar);
	//return (((input+1024)%32768)/2048);
}

uint16_t cross (uint16_t input){
	return ((input%32768)/4096)%2;
	//return ((input%32768)>>(b_bar+1))%2;
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
	/*
	if (cross(w)==b) return rounding(w);
	if (w%4096 <2048){
		return rounding(  ((w-w%4096)-1)%32768 );
	}
	return rounding(  ((w-w%4096+4096))%32768 );
	*/

	uint16_t equal = (cross(w) & b) + ((1^cross(w)) & (b^1));
	uint16_t output = equal*w;

	uint16_t want_high=greater_than(w%4096,2047);
	uint16_t output1=(want_high ^ 1)*(  ((w-w%4096)-1)%32768 )+(want_high ^ 0)*(  ((w-w%4096+4096))%32768 );
	output=output+(equal ^1)*output1;
	//uint16_t output1=(want_zero ^ 1) * (greater_than(w,2047)*4096 +(greater_than(w,2047) ^1)*32767 );
	//output1+=(want_zero ) * (greater_than(w,18433)*0 +(greater_than(w,18433) ^1)*4095 );
	//output=(equal ^1)*output1;
	return rounding(output);


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

void main4(){
	uint16_t A[2][2];
	A[0][0]=5;
	A[0][1]=1;
	A[1][0]=3;
	A[1][1]=7;

	uint16_t x[2][1];
	x[0][0]=10;
	x[0][1]=100;

	uint16_t b[2][1];
	b[0][0]=0;
	b[0][1]=0;

	int i,j,i1,j1;
	for ( i = 0; i < 2; i++ ){
		for ( j = 0; j < 1; j++ ){


			for ( i1 = 0; i1 < 2; i1++ ){
				//for ( j1 = 0; j1 < 2; j1++ ){
					b[i][j]=b[i][j]+A[i][i1]*x[i1][j];
				//}
			}


		}

	}
	printf("we have %u %u \n",b[0][0],b[0][1]);

}

void main5(){
	unsigned short i;
	printf("We have %u %u %u %u %u\n",cross(32766),cross(32767),cross(32768), cross(0), cross(1));
	unsigned short last=0;
	for (i=0; i<32768; i++){
		if (cross(i)!=last){
			printf("%u %u \n",last, i);
			last=cross(i);
		}

	}

}
void print_mat(uint16_t mat[][8]){
	int i,j;
	printf("\n MATRIX\n");
	for ( i = 0; i < 8; i++ ){
		for ( j = 0; j < 8; j++ ){
			
			printf("%u, ",mat[i][j]);
		}
		printf("\n");
	}
	printf("\n \n");
}

void print_matb(uint16_t mat[][752]){
	printf("\n MATRIX\n");
	int i,j;
	for ( i = 0; i < 8; i++ ){
		for ( j = 0; j < 8; j++ ){
			
			printf("%u, ",mat[i][j]);
		}
		printf("\n");
	}
	printf("\n \n");
}
void main(){
	uint16_t A[752][752];
	uint16_t q=32768;
	int i,j;

	//This is what Alice does
	////////////////////////////////////////////////////////////

	for ( i = 0; i < 752; i++ ){
		for ( j = 0; j < 752; j++ ){
			//A[i][j]=(1675*j+1513*i*j+689*i+j*j*500+9850)%q;
			A[i][j]=i+j;

		}

	}
	uint16_t rnd1[9026];
	uint16_t rnd2[9026];
	for ( i = 0; i < 9026; i++ ){
		rnd1[i]=(i*i*i+2000*i-9999*i)%q;
		rnd2[i]=(i*i*5303-200)%q;
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
				//for ( j1 = 0; j1 < 752; j1++ ){
					b[i][j]=b[i][j]+A[i][i1]*s[i1][j];
				//}
			}


		}

	}
	
	//print_matb(A);
	//print_mat(b);
	////////////////////////////////////////////////////////////
	//This is what Bob does
	uint16_t Ab[752][752];
	for ( i = 0; i < 752; i++ ){
		for ( j = 0; j < 752; j++ ){
			Ab[i][j]=(i+j)%q;

		}

	}
	uint16_t rnd3[9026];
	uint16_t rnd4[9026];
	for ( i = 0; i < 9026; i++ ){
		rnd3[i]=(i*i*i*i -i*i*i*1103	+10000*i)%q;
		rnd4[i]=(i*i*i*i*959+i*i*i*105+i*i*7717)%q;
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
		for ( j = 0; j < 752; j++ ){


			for ( i1 = 0; i1 < 752; i1++ ){
				//for ( j1 = 0; j1 < 752; j1++ ){
					bp[i][j]=bp[i][j]+sp[i][i1]*Ab[i1][j];
				//}
			}


		}

	}


	uint16_t epp[8][8];
	uint16_t rnd5[98];
	for ( i = 0; i < 98; i++ ){
		rnd5[i]=(i*i*i*i*8587- i*i*i*505+ i*5503)%32768;
	}
	lwe_sample_n_inverse_12_64(epp,rnd5);

	uint16_t v[8][8];
	for ( i = 0; i < 8; i++ ){
		for ( j = 0; j < 8; j++ ){
			v[i][j]=epp[i][j];
		}
	}
	for ( i = 0; i < 8; i++ ){
		for ( j = 0; j < 8; j++ ){


			for ( i1 = 0; i1 < 752; i1++ ){
				//for ( j1 = 0; j1 < 752; j1++ ){
					v[i][j]=v[i][j]+sp[i][i1]*b[i1][j];
				//}
			}


		}

	}


	uint16_t c[8][8];
	for ( i = 0; i < 8; i++ ){
		for ( j = 0; j < 8; j++ ){
			c[i][j]=cross((v[i][j]%q));
		}
	}


/*
	printf("This is v..................\n");
	for ( i = 0; i < 8; i++ ){
		for ( j = 0; j < 8; j++ ){
			
			printf("%u, ",v[i][j]);
		}
	}
	printf("\n");
*/

	uint16_t k_bob[8][8];
	for ( i = 0; i < 8; i++ ){
		for ( j = 0; j < 8; j++ ){
			k_bob[i][j]=rounding((v[i][j]%q));
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
				//for ( j1 = 0; j1 < 752; j1++ ){
					bp_s[i][j]=bp_s[i][j]+bp[i][i1]*s[i1][j];
				//}
			}


		}

	}
	print_matb(bp);
	printf("s..............\n");
	print_mat(s);
	print_mat(bp_s);

	
 
	uint16_t k_alice[8][8];
	for ( i = 0; i < 8; i++ ){
		for ( j = 0; j < 8; j++ ){
			bp_s[i][j]=bp_s[i][j]%q;
			k_alice[i][j]=rec(bp_s[i][j],c[i][j]);
		}
	}


///////////////////////////////////////////////////////////////////////////////////////
	printf("\n v MATRIX\n");
	for ( i = 0; i < 8; i++ ){
		for ( j = 0; j < 8; j++ ){
			
			printf("%u, ",v[i][j]);
		}
		printf("\n");
	}
	printf("\n \n");


	printf("\nbp_s MATRIX\n");
	for ( i = 0; i < 8; i++ ){
		for ( j = 0; j < 8; j++ ){
			
			printf("%u, ",bp_s[i][j]);
		}
		printf("\n");
	}
	printf("\n c MATRIX\n");		

/*

	printf("\n");
	for ( i = 0; i < 8; i++ ){
		for ( j = 0; j < 8; j++ ){
			
			printf("[%u, %u], ",bp_s[i][j],c[i][j]);
		}
		printf("\n");
	}
	printf("\n");


	printf("\n v MATRIX\n");
	for ( i = 0; i < 8; i++ ){
		for ( j = 0; j < 8; j++ ){
			
			printf("%u, ",v[i][j]);
		}
		printf("\n");
	}
	printf("\n \n");


	printf("\nbp_s MATRIX\n");
	for ( i = 0; i < 8; i++ ){
		for ( j = 0; j < 8; j++ ){
			
			printf("%u, ",bp_s[i][j]);
		}
		printf("\n");
	}
	printf("\n c MATRIX\n");

	printf("\n");
	for ( i = 0; i < 8; i++ ){
		for ( j = 0; j < 8; j++ ){
			
			printf("%u, ",c[i][j]);
		}
		printf("\n");
	}
	printf("\n");
*/

	printf("\n K_ALICE\n");
	for ( i = 0; i < 8; i++ ){
		for ( j = 0; j < 8; j++ ){
			
			printf("%u, ",k_alice[i][j]);
		}
		printf("\n");
	}
	printf("\n");

	printf("\n K_BOB\n");
	for ( i = 0; i < 8; i++ ){
		for ( j = 0; j < 8; j++ ){
			
			printf("%u, ",k_bob[i][j]);
		}
		printf("\n");
	}
	printf("\n");


	printf("\n \n");
	for ( i = 0; i < 8; i++ ){
		for ( j = 0; j < 8; j++ ){
			
			if (epp[i][j]>5 && (32768- epp[i][j])>6) printf("%u,",epp[i][j]);
		}
		
	}
	printf("\n");


	printf("Done...................\n");

}
