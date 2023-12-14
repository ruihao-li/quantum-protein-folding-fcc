
#define _CRT_SECURE_NO_WARNINGS

#include <stdlib.h>
#include <stdio.h>
#include <cstdint> // include this header for uint64_t and int64_t


// # of amino acids
#define	AMINOACIDS	17


// NOTE: LIMITATION OF THIS VERSION:
// MAX QBITS = 128 (using two 64-bit unsigned integers for bitstring)
// AND, conformation bitstring and interaction bitstring lengths each must not exceed 64

//#define QBITS	22
#define QBITS	69
#define MAXTERMS	16384
#define TOP	100

//#define THREADS		8
#define MAXTHREADS		64


int main()
{
	int64_t index, permindex;
	uint64_t perm, perm_int, ilong, temp, temp2;
	int a, b, i, j, k, kk, threadnum, threads, bitval, p, q, ctr, bitshift;
	char filename[4096];
	//char x[QBITS], xbest[MAXTHREADS][QBITS], xtop[MAXTHREADS][TOP][QBITS], xtop_all[TOP][QBITS];
	uint64_t x, xbest[MAXTHREADS], xtop[MAXTHREADS][TOP], xtop_all[TOP];
	uint64_t x_2, xbest_2[MAXTHREADS], xtop_2[MAXTHREADS][TOP], xtop_all_2[TOP];
	double obj, objbest[MAXTHREADS], objtop[MAXTHREADS][TOP], objtop_all[TOP];
	float *qubo;
	double dotproduct[QBITS];
	
	int n0, n1, n2, n3, d2, qbits_int, qbits_config;
	// turnseq_abs[] uses the lattice directions {0,1,2,3} or {0-bar,1-bar,2-bar,3-bar}
	unsigned char turnseq_abs[AMINOACIDS - 1];
	// turnseq_rel[] is the 3-turn sequence in Boulebnane et al, each turn has value {0,1,2}
	unsigned char turnseq_rel[AMINOACIDS - 1];


	char str1[64], str2[64], str3[2048], str4[2048];
	int terms, ising[QBITS], spinproduct;
	double h, val, offset, coef[MAXTERMS];
	uint64_t y, y1, y2, mask[QBITS], temp_uint64;
	uint64_t bits[MAXTERMS], hbits, bitstring, bitstring_config, bitstring_inter;
	uint64_t bits_2[MAXTERMS], hbits_2, bitstring_2;
	FILE *f;


	
	// constants for Hamming weight calculation
	// see Wikipedia:  https://en.wikipedia.org/wiki/Hamming_weight
	// uses C-standard types (defined in C99 version of C language)
	const uint64_t m1 = 0x5555555555555555; //binary: 0101...
	const uint64_t m2 = 0x3333333333333333; //binary: 00110011..
	const uint64_t m4 = 0x0f0f0f0f0f0f0f0f; //binary:  4 zeros,  4 ones ...
	const uint64_t m8 = 0x00ff00ff00ff00ff; //binary:  8 zeros,  8 ones ...
	const uint64_t m16 = 0x0000ffff0000ffff; //binary: 16 zeros, 16 ones ...
	const uint64_t m32 = 0x00000000ffffffff; //binary: 32 zeros, 32 ones
	const uint64_t h01 = 0x0101010101010101; //the sum of 256 to the power of 0,1,2,3...

	const uint64_t one64u = 1;


	// read in Hamiltonian text file
	f = fopen("hamiltonian.txt", "r");

	// first line - constant term
	k = fscanf(f, "%s %s %s\n", str1, str2, str3);
	sscanf(str1, "%lf", &offset);


	// parse each term (Pauli Z's and I's),
	// store the bitmasks as a uint64 array
	// (note: the text strings are in reverse order)
	terms = 0;
	while (terms >= 0)
	{
		k = fscanf(f, "%s %s %s %s\n", str1, str2, str3, str4);
		if (k != 4)
			break;

		sscanf(str2, "%lf", &val);
		if (str1[0] == '-')
		{
			val *= -1.;
		}
		coef[terms] = val;

		bits[terms] = 0;
		bits_2[terms] = 0;
		for (i = 0; i < QBITS; i++)
		{
			switch (str4[i])
			{
			case 'Z':
				bitshift = (QBITS - 1 - i);
				if (bitshift <= 64)
					bits[terms] += ((uint64_t)1) << bitshift;
				else
					bits_2[terms] += ((uint64_t)1) << (bitshift - 64);
				break;
			case 'I':
				break;
			}
		}

		terms++;
	}

	fclose(f);

	printf("Number of terms in Hamiltonian = %d\n", terms);


	threads = 1;
	printf("max number of threads = %d\n", threads);


	for (j = 0; j < threads; j++)
	{
		objbest[j] = 1e12;
		for (i = 0; i < TOP; i++)
		{
			objtop[j][i] = 1e12;
		}
	}


	// load in the amino acid sequence and interaction energies from a file
	// 
	// for now, assume no side chains!...
	// 
	// aa_parent[] lists which amino acid is its parent on the relative turn sequence (note, first three always {0,0,1} )
	//unsigned char aa_parent[AMINOACIDS] = { 0,0,1,2,3,4,5,6,7,8,9,10 };
	// aa_main[] lists the order of the amino acid on the main backbone
	//unsigned char aa_main[AMINOACIDS] = { 0,1,2,3,4,5,6,7,8,9,10,11 };
	// aa_side[] lists the order of the amino acid on a side chain (or, = 0 on the main chain)
	//unsigned char aa_side[AMINOACIDS] = { 0,0,0,0,0,0,0,0,0,0,0,0 };
	unsigned char aa_parent[AMINOACIDS], aa_main[AMINOACIDS], aa_side[AMINOACIDS];
	aa_parent[0] = aa_main[0] = aa_side[0] = 0;
	for (i = 1; i < AMINOACIDS; i++)
	{
		aa_parent[i] = i - 1;
		aa_main[i] = i;
		aa_side[i] = 0;
	}

	// encode all 8 absolute directions as a single 64-bit integer
	// lowest-to-highest 8-bit chunks:  
	//     set = 1 to represent absolute directions: 0, 1, 2, 3, 0bar, 1bar, 2bar, 3bar
	// then, just sum these up to easily compute n_a(i,j) and n_abar(i,j)
	//     no need to define and loop over turn indicator functions f_a(i)
	//     having 8 bits per direction limits problem to 256 amino acids, not a problem
	//
	uint64_t turn_index[AMINOACIDS-1], turn_accum, turn_accum_i, turn_accum_j;


	// loop over all 3-turn permutations
	//
	// Note: turns #0 and #1 are specified.
	// Also, if AA#1 (second bead) has no side chain, then turn #2 has 2 possibilities, so total permutations = 2 * 3^(AMINOACIDS-4) 
	// For simplicity (for now), use 3^(AMINOACIDS-3) permutations and reject.
	// Get this working first, then optimize
	perm = 1;
	for (i = 2; i < (AMINOACIDS-1); i++)
	{
		perm *= 3;
	}	// perm = 3^(N-3)
		// 
	if (aa_side[1] == 0)	// no side chain on second AA, fix q6=1 (see Robert et al paper)
	{
		qbits_config = 2 * (AMINOACIDS - 3) - 1;
	}
	else
	{
		qbits_config = 2 * (AMINOACIDS - 3);
	}
	qbits_int = QBITS - qbits_config;
	perm_int = (uint64_t)1 << qbits_int;


	//
	for (permindex = 0; permindex < perm; permindex++)
		//for (permindex = (perm-1); permindex >=0; permindex--)	// reverse order, for debugging
	{

		// generate the turn sequence array from the permutation index
		// see the Boulebnane et al paper about converting from relative to absolute turn directions
		//
		// first two turns are fixed
		turnseq_rel[0] = 3;		// doesn't matter, no prior AA...
		turnseq_abs[0] = 1;
		turn_index[0] = ((uint64_t)1 << (8 * (turnseq_abs[0] + 4 * ((0 + 1) % 2))));
		turnseq_rel[1] = 2;
		turnseq_abs[1] = 0;		// (turnseq_rel[1] + turnseq_abs[aa_parent[1]] + 1) % 4;
		turn_index[1] = ((uint64_t)1 << (8 * (turnseq_abs[1] + 4 * ((1 + 1) % 2))));
		// initialize:
		temp = permindex;
		// second turn
		temp2 = temp / 3;
		turnseq_rel[2] = temp - temp2 * 3;
		turnseq_abs[2] = (turnseq_rel[2] + turnseq_abs[aa_parent[2]] + 1) % 4;
		turn_index[2] = ((uint64_t)1 << (8 * (turnseq_abs[2] + 4 * ((2 + 1) % 2))));
		temp = temp2;
		// check if second amino acid does not have a side chain
		if (aa_side[1] == 0)
			if ((turnseq_abs[2] % 2) == 0)	// q6=1, reject if not
				continue;			// so only proceed if turn#2 (absolute) = 1 or 3
		// other turns
		for (i = 3; i < (AMINOACIDS - 1); i++)
		{
			temp2 = temp / 3;
			turnseq_rel[i] = temp - temp2 * 3;
			turnseq_abs[i] = (turnseq_rel[i] + turnseq_abs[aa_parent[i]] + 1) % 4;
			turn_index[i] = ((uint64_t)1 << (8 * (turnseq_abs[i] + 4 * ((i + 1) % 2))));
			temp = temp2;
		}

		// generate qubit configuration bitstring from turn sequence
		bitstring_config = 0;
		if (aa_side[1] == 0)	// no side chain, start with just one bit 
		{
			bitstring_config += turnseq_abs[2] / 2;		// turn direction 1 and 3 are q5=0 and q5=1, respectively
			for (i = 3; i < (AMINOACIDS - 1); i++)
			{
				// Least significant bit is the higher qubit (why qiskit, why?)
				temp_uint64 = turnseq_abs[i] % 2;
				bitstring_config += (temp_uint64 << (2 * (i - 2)));
				// and, the most significant bit is the lower qubit
				temp_uint64 = turnseq_abs[i] / 2;
				bitstring_config += (temp_uint64 << (2 * (i - 2)) - 1);
			}
		}
		else	// side chain, start with two bits
		{
			for (i = 2; i < (AMINOACIDS - 1); i++)
			{
				// LSB
				temp_uint64 = turnseq_abs[i] % 2;
				bitstring_config += (temp_uint64 << (2 * (i - 2)) + 1);
				// MSB
				temp_uint64 = turnseq_abs[i] / 2;
				bitstring_config += (temp_uint64 << (2 * (i - 2)));
			}
		}

		// initialize interaction bitstring
		bitstring_inter = 0;

		// for debugging...
		//if (bitstring_config == 0xdc72 || bitstring_config == 0xdc73)
		//	p = 0;


		// calculate Hamiltonian

		// for each AA pair under consideration, compute the lattice distance d(i,j)
		// for now, consider only first-NN interactions:  j >= (i+5)
		// for second-NN interactions, need to consider j >= (i+4)
		// 
		// Not really sure if this works perfectly, need to trace and debug!
		//
		// loop over all AA pairs, regardless of main chain or side chain
		ctr = -1;
		for (a = 0; a < (AMINOACIDS - 5); a++)
		{
			for (b = (a + 5); b < AMINOACIDS; b++)
			{
				if (((b - a)) % 2 == 1)
				{
					ctr++;
				}
				else
				{
					continue;
				}

				// determine the main chain starting point i and main chain ending point j
				i = a;	//i = aa_parent[a];
				j = b;	//j = aa_parent[b];

				// initialize
				turn_accum = turn_accum_i = turn_accum_j = 0;
				//
				// main chain:
				for (k = i; k < j; k++)
				{
					//if (aa_side[k+1] == 0)	
					turn_accum += turn_index[k];
				}
				// side chain i
				//if (aa_side[a] != 0)
				//{
				//	for (k = i; k < a; k++)
				//		turn_accum_i += turn_index[k];
				//}
				// side chain j		
				//if (aa_side[b] != 0)
				//{
				//	for (k = j; k < b; k++)
				//		turn_accum_j += turn_index[k];
				//}

				// calculate distance
				//
				// main chain
				n0 = (int)(turn_accum & 0x00000000000000ff) - (int)((turn_accum & 0x000000ff00000000) >> 32);
				n1 = (int)((turn_accum & 0x000000000000ff00) >> 8) - (int)((turn_accum & 0x0000ff0000000000) >> 40);
				n2 = (int)((turn_accum & 0x0000000000ff0000) >> 16) - (int)((turn_accum & 0x00ff000000000000) >> 48);
				n3 = (int)((turn_accum & 0x00000000ff000000) >> 24) - (int)((turn_accum & 0xff00000000000000) >> 56);
				// side chain i
				//n0 -= ((int)(turn_accum_i & 0x00000000000000ff) - (int)((turn_accum_i & 0x000000ff00000000) >> 32));
				//n1 -= ((int)((turn_accum_i & 0x000000000000ff00) >> 8) - (int)((turn_accum_i & 0x0000ff0000000000) >> 40));
				//n2 -= ((int)((turn_accum_i & 0x0000000000ff0000) >> 16) - (int)((turn_accum_i & 0x00ff000000000000) >> 48));
				//n3 -= ((int)((turn_accum_i & 0x00000000ff000000) >> 24) - (int)((turn_accum_i & 0xff00000000000000) >> 56));
				// side chain j
				//n0 += ((int)(turn_accum_j & 0x00000000000000ff) - (int)((turn_accum_j & 0x000000ff00000000) >> 32));
				//n1 += ((int)((turn_accum_j & 0x000000000000ff00) >> 8) - (int)((turn_accum_j & 0x0000ff0000000000) >> 40));
				//n2 += ((int)((turn_accum_j & 0x0000000000ff0000) >> 16) - (int)((turn_accum_j & 0x00ff000000000000) >> 48));
				//n3 += ((int)((turn_accum_j & 0x00000000ff000000) >> 24) - (int)((turn_accum_j & 0xff00000000000000) >> 56));

				d2 = n0 * n0 + n1 * n1 + n2 * n2 + n3 * n3;

				if (d2 == 1)
				{
					// add the 1NN Hamiltonian term

					bitstring_inter += one64u << ctr;
				}

			}
		}

		// stuff the bits in the bitstrings
		// this version assumes config bits <= 64 and interaction bits <= 64
		bitstring = bitstring_config + (bitstring_inter << qbits_config);
		if (QBITS > 64)
			bitstring_2 = bitstring_inter >> (64 - qbits_config);
		else
			bitstring_2 = 0;


		threadnum = 0;
		obj = 0.;


		for (i = 0; i < terms; i++)
		{
			hbits = bits[i];	// uint64 representing spin state

			y = hbits & bitstring;	// uint64 matching bits

			// calculate Hamming weight
			// see Wikipedia:  https://en.wikipedia.org/wiki/Hamming_weight
			//
			//
			y -= (y >> 1) & m1;             //put count of each 2 bits into those 2 bits
			y = (y & m2) + ((y >> 2) & m2); //put count of each 4 bits into those 4 bits 
			y = (y + (y >> 4)) & m4;        //put count of each 8 bits into those 8 bits 
			y += y >> 8;  //put count of each 16 bits into their lowest 8 bits
			y += y >> 16;  //put count of each 32 bits into their lowest 8 bits
			y += y >> 32;  //put count of each 64 bits into their lowest 8 bits
			y = (y & 0x7f);
			//

			if (QBITS > 64)		// repeat for second 64-bit chunk
			{
				hbits_2 = bits_2[i];	// uint64 representing spin state;

				y2 = hbits_2 & bitstring_2;	// uint64 matching bits

				y2 -= (y2 >> 1) & m1;             //put count of each 2 bits into those 2 bits
				y2 = (y2 & m2) + ((y2 >> 2) & m2); //put count of each 4 bits into those 4 bits 
				y2 = (y2 + (y2 >> 4)) & m4;        //put count of each 8 bits into those 8 bits 
				y2 += y2 >> 8;  //put count of each 16 bits into their lowest 8 bits
				y2 += y2 >> 16;  //put count of each 32 bits into their lowest 8 bits
				y2 += y2 >> 32;  //put count of each 64 bits into their lowest 8 bits
				y2 = (y2 & 0x7f);

				y += y2;
			}

			// if y is odd, then the ising product is negative, so subtract the coefficient
			// if y is even, then the ising product is positive, so add the coefficient
			if (y & one64u)
			{
				obj -= coef[i];
			}
			else
			{
				obj += coef[i];
			}

		}

		if (obj < objbest[threadnum])
		{
			objbest[threadnum] = obj;
			xbest[threadnum] = (uint64_t)bitstring;
			xbest_2[threadnum] = (uint64_t)bitstring_2;

			printf("new best objective = %lf", objbest[threadnum]);
			printf(", bitstring = %llu", xbest[threadnum]);
			if (QBITS > 64)
				printf(", bitstring_2 = %llu", xbest_2[threadnum]);
			printf("\n");

		}

		if (obj < objtop[threadnum][TOP - 1])	// cracked the top ten
		{
			// determine rank on list
			for (j = 0; j < TOP; j++)
			{
				if (obj < objtop[threadnum][j])
					break;
			}
			for (i = (TOP - 1); i > j; i--)
			{
				objtop[threadnum][i] = objtop[threadnum][i - 1];
				xtop[threadnum][i] = xtop[threadnum][i - 1];
				xtop_2[threadnum][i] = xtop_2[threadnum][i - 1];
			}
			objtop[threadnum][j] = obj;
			xtop[threadnum][j] = (uint64_t)bitstring;
			xtop_2[threadnum][j] = (uint64_t)bitstring_2;

		}

	}

	// compile top list from all threads
	for (i = 0; i < TOP; i++)
	{
		objtop_all[i] = 1e12;
	}
	for (threadnum = 0; threadnum < threads; threadnum++)
	{
		for (k = 0; k < TOP; k++)
		{
			if (objtop[threadnum][k] < objtop_all[TOP - 1])	// cracked the top group
			{
				// determine rank on list
				for (j = 0; j < TOP; j++)
				{
					if (objtop[threadnum][k] < objtop_all[j])
						break;
				}
				for (i = (TOP - 1); i > j; i--)
				{
					objtop_all[i] = objtop_all[i - 1];
					xtop_all[i] = xtop_all[i - 1];
					xtop_all_2[i] = xtop_all_2[i - 1];
				}
				objtop_all[j] = objtop[threadnum][k];
				xtop_all[j] = xtop[threadnum][k];
				xtop_all_2[j] = xtop_2[threadnum][k];
			}
		}
	}
	f = fopen("topobj.txt", "w");
	for (i = 0; i < TOP; i++)
	{
		fprintf(f, "%lf\t", (objtop_all[i] + offset));
		for (k = 0; k < QBITS; k++)		// right to left (LSB to MSB)
		{
			bitshift = (QBITS - 1 - k);
			if (bitshift <= 64)
			{
				bitval = (xtop_all[i] & (one64u << bitshift)) ? 1 : 0;
			}
			else
			{
				bitval = (xtop_all_2[i] & (one64u << (bitshift - 64))) ? 1 : 0;
			}
			fprintf(f, "%d", bitval);
		}
		fprintf(f, "\n");
	}
	fclose(f);


	return 0;
}
