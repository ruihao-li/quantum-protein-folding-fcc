
#define _CRT_SECURE_NO_WARNINGS

#include <stdlib.h>
#include <stdio.h>
#include <cstdint> // include this header for uint64_t and int64_t
#include <math.h>


#ifdef _OPENMP
#include <omp.h>
#define OMP 1
#endif


// define peptide sequence:
// 
//#define AMINOACIDS	7
//#define PEPTIDE "APRLRFY"	// from IBM paper
#define AMINOACIDS	10
//#define PEPTIDE		"YYDPETGTWY"	// chignolin
#define PEPTIDE			"QYQFWKNFQT"	// 2MZX
//#define PEPTIDE		"DRVYVHPFHL"	// incorrect angiotensin sequence from Robert et al paper
//#define PEPTIDE		"DRVYIHPFHL"	// correct angiotensin sequence
//#define AMINOACIDS	16
//#define PEPTIDE		"RHYYKFNSTGRHYHYY"	// 8T61	designer peptide BH33
//#define PEPTIDE		"WHMWNTVPNAKQVIAA"	// 8T63	designer peptide PH1
//#define AMINOACIDS	20
//#define PEPTIDE		"DAYAQWLKDGGPSSGRPPPS"	// modified trp cage
//#define PEPTIDE		"NLYIQWLKDGGPSSGRPPPS"	// original Trp Cage
// #define AMINOACIDS	25
// #define PEPTIDE		"KKPGASLAALQALQALQAAQAAKKY"	// 25AA 8B1X helical peptide
//#define AMINOACIDS	26
//#define PEPTIDE		"YYHFWHRGVTKRSLSPHRPRHSRLQR"	// 26AA 6A8Y helix + turn
//#define PEPTIDE		"GNDYEDRYYRENMYRYPNQVYYRPVC"	// 26AA 1G04 sheep prion


#define TOP	100

//#define THREADS		8
#define MAXTHREADS		64


int main()
{
	int64_t index, permindex;
	uint64_t perm, perm_int, ilong, temp, temp2;
	int a, b, i, j, k, kk, threadnum, threads, bitval, p, q, ctr, bitshift, collision;
	char filename[4096];
	//char x[QBITS], xbest[MAXTHREADS][QBITS], xtop[MAXTHREADS][TOP][QBITS], xtop_all[TOP][QBITS];
	uint64_t x, xbest[MAXTHREADS], xtop[MAXTHREADS][TOP], xtop_all[TOP];
	uint64_t permbest[MAXTHREADS], permtop[MAXTHREADS][TOP], permtop_all[TOP];
	uint64_t x_2, xbest_2[MAXTHREADS], xtop_2[MAXTHREADS][TOP], xtop_all_2[TOP];
	double obj, objbest[MAXTHREADS], objtop[MAXTHREADS][TOP], objtop_all[TOP];
	
	int n0, n1, n2, n3, d2, qbits_int, qbits_config;
	// turnseq_abs[] uses the lattice directions {0,1,2,3} or {0-bar,1-bar,2-bar,3-bar}
	unsigned char turnseq_abs[AMINOACIDS - 1], turnseq_abs_best[AMINOACIDS - 1][MAXTHREADS];
	unsigned char turnseq_abs_top[AMINOACIDS - 1][MAXTHREADS][TOP], turnseq_abs_top_all[AMINOACIDS - 1][TOP];
	// turnseq_rel[] is the 3-turn sequence in Boulebnane et al, each turn has value {0,1,2}
	unsigned char turnseq_rel[AMINOACIDS - 1];


	char str1[64], str2[64], str3[2048], str4[2048];
	uint64_t y, y1, y2, temp_uint64;
	uint64_t bitstring, bitstring_config, bitstring_inter;
	FILE *f, *g;

	char aa_string[] = PEPTIDE;
	char aa_list[20];
	int aa_index[AMINOACIDS];
	double eij[20][20], aa_eij[AMINOACIDS][AMINOACIDS], xx, yy, zz;
	double aa_scale = 3.8;

	double d1_d2_ratio = sqrt(3./8.);

	const uint64_t one64u = 1;



	// read in the interaction energies from .csv file
	f = fopen("mj_matrix.csv", "r");
	//f = fopen("mj_matrix_1985.csv", "r");
	//f = fopen("contactenergies_2003.csv", "r");
	// first line is the amino acid list
	for (i=0; i<20; i++)
	{
		fscanf(f, "%c,", &aa_list[i]);
	}
	// the next 20 lines is the energy matrix
	for (j = 0; j < 20; j++)
	{
		for (i = 0; i < 20; i++)
		{
			fscanf(f, "%lf,", &eij[i][j]);
		}
	}
	fclose(f);

	// for each AA in the peptide, figure out its index
	for (i = 0; i < AMINOACIDS; i++)
	{
		for (k = 0; k < 20; k++)
		{
			if (aa_string[i] == aa_list[k])
			{
				aa_index[i] = k;
			}
		}
	}
	// load up the interaction energies into the look-up matrix
	for (j = 0; j < AMINOACIDS; j++)
	{
		for (i = 0; i < AMINOACIDS; i++)
		{
			aa_eij[i][j] = eij[aa_index[i]][aa_index[j]];
		}
	}


#ifdef OMP
	threads = omp_get_max_threads();
#else
	threads = 1;
#endif
	printf("max number of threads = %d\n", threads);


	// initialize the top objective energies with large positive number
	for (j = 0; j < threads; j++)
	{
		objbest[j] = 1e12;
		for (i = 0; i < TOP; i++)
		{
			objtop[j][i] = 1e12;
		}
	}


	// define the amino acid sequence information (main vs side chain, and parent index)
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
	//
	uint64_t turn_index[AMINOACIDS], turn_accum, turn_accum_prior, turn_accum_next, turn_accum_i, turn_accum_j;


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
	

	// for qubit formulation...
	if (aa_side[1] == 0)	// no side chain on second AA, fix q6=1 (see Robert et al paper)
	{
		qbits_config = 2 * (AMINOACIDS - 3) - 1;
	}
	else
	{
		qbits_config = 2 * (AMINOACIDS - 3);
	}
	//qbits_int = QBITS - qbits_config;
	//perm_int = (uint64_t)1 << qbits_int;




	//
#pragma omp parallel for schedule(dynamic,1), default(shared) \
	private(permindex,turnseq_rel,turnseq_abs,turn_index,temp,temp2,i,j,k,ctr,collision,obj,a,b) \
	private(turn_accum,turn_accum_i,turn_accum_j,n0,n1,n2,n3,d2,turn_accum_prior,turn_accum_next,threadnum) \
	shared(perm,aa_side,aa_parent,aa_eij,objbest,objtop,permbest,permtop,turnseq_abs_best,turnseq_abs_top) \
	num_threads(threads)
	//
	for (permindex = 0; permindex < perm; permindex++)
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
		// turn_index[AMINOACIDS] only represents (AMINOACIDS-1) turns, but set the last value to zero (no turn)
		turn_index[(AMINOACIDS - 1)] = 0;


		// calculate Hamiltonian on-the-fly

		//
		// non-backtracking is enforced by relative 3-turn sampling
		// self-avoiding walk is enforced by collision check (d2 == 0)
		// to do:  check chirality constraints, if considering side chains
		//

		// for each AA pair under consideration, compute the lattice distance d(i,j)
		// consider only 1-NN interactions:  j >= (i+5)
		// 
		// if a 1-NN is found, then check the prior and next AA to see if d2 == 0
		// 
		// Not 100% sure if this works perfectly, need to trace and debug!
		//
		// loop over all AA pairs, regardless of main chain or side chain
		// 
		// initialize:
		ctr = -1;
		collision = 0;
		obj = 0.;
		//
		for (a = 0; a < (AMINOACIDS - 4); a++)
		{
			// d2 = 0 had been found, break out of loop, and go to next turn sequence permutation
			if (collision)
			{
				obj = 10000.;
				break;
			}

			for (b = (a + 4); b < AMINOACIDS; b++)
			{
				// can only have a 1-NN interaction if an odd difference between a and b
				//if (((b - a)) % 2 == 0)	// both on same sub-lattice, cannot be 1-NN
				//{
				//	continue;
				//}

				//ctr++;	// currently not used...

				// determine the main chain starting point i and main chain ending point j
				i = a;	//i = aa_parent[a];		// if considering side chains...
				j = b;	//j = aa_parent[b];

				// initialize the accumulation of turns
				turn_accum = 0;
				//turn_accum_i = turn_accum_j = 0;
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
					// check if the prior AA collides
					turn_accum_prior = turn_accum - turn_index[j - 2];
					if ((turn_accum_prior & 0x00000000ffffffff) == (turn_accum_prior >> 32))
					{
						collision = 1;
					}
					// check if the next AA collides
					turn_accum_next = turn_accum + turn_index[j];
					if ((turn_accum_next & 0x00000000ffffffff) == (turn_accum_next >> 32))
					{
						collision = 1;
					}

					if (collision)
					{
						obj = 10000.;
						break;
					}
					else
					{
						// add the 1-NN interaction energy
						obj += aa_eij[i][j];
					}

				}

				if (d2 == 2) //for 2NN

				{
					// check if the prior AA collides
					turn_accum_prior = turn_accum - turn_index[j - 2];
					if ((turn_accum_prior & 0x00000000ffffffff) == (turn_accum_prior >> 32))
					{
						collision = 1;
					}
					// check if the next AA collides
					turn_accum_next = turn_accum + turn_index[j];
					if ((turn_accum_next & 0x00000000ffffffff) == (turn_accum_next >> 32))
					{
						collision = 1;
					}

					if (collision)
					{
						obj = 10000.;
						break;
					}
					else
					{
					// add the 2-NN interaction energy
						obj += aa_eij[i][j] * d1_d2_ratio;
					}

				}

			}
		}



#ifdef OMP
		threadnum = omp_get_thread_num();
#else
		threadnum = 0;
#endif
		
		

		if (obj < objbest[threadnum])
		{
			objbest[threadnum] = obj;
			permbest[threadnum] = (uint64_t)permindex;
			for (k = 0; k < (AMINOACIDS - 1); k++)
			{
				turnseq_abs_best[k][threadnum] = turnseq_abs[k];
			}

			printf("new best objective = %lf", objbest[threadnum]);
			printf(", perm = %llu", permbest[threadnum]);
			printf(", threadnum = %d\n", threadnum);

		}

		if (obj < objtop[threadnum][TOP - 1])	// cracked the top list
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
				permtop[threadnum][i] = permtop[threadnum][i - 1];
				for (k = 0; k < (AMINOACIDS - 1); k++)
				{
					turnseq_abs_top[k][threadnum][i] = turnseq_abs_top[k][threadnum][i-1];
				}
			}
			objtop[threadnum][j] = obj;
			permtop[threadnum][j] = (uint64_t)permindex;
			for (k = 0; k < (AMINOACIDS - 1); k++)
			{
				turnseq_abs_top[k][threadnum][j] = turnseq_abs[k];
			}
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
					permtop_all[i] = permtop_all[i - 1];
					for (a = 0; a < (AMINOACIDS - 1); a++)
					{
						turnseq_abs_top_all[a][i] = turnseq_abs_top_all[a][i - 1];
					}
				}
				objtop_all[j] = objtop[threadnum][k];
				permtop_all[j] = permtop[threadnum][k];
				for (a = 0; a < (AMINOACIDS - 1); a++)
				{
					turnseq_abs_top_all[a][j] = turnseq_abs_top[a][threadnum][k];
				}
			}
		}
	}
	f = fopen("topobj.txt", "w");
	for (i = 0; i < TOP; i++)
	{
		fprintf(f, "%lf\t", objtop_all[i]);
		//fprintf(f, "%llu\t", permtop_all[i]);
		for (k = 0; k < (AMINOACIDS - 1); k++)
		{
			fprintf(f, "%d", turnseq_abs_top_all[AMINOACIDS - 2 - k][i]);
		}
		fprintf(f, "\t");

		// generate qubit configuration bitstring from turn sequence
		bitstring_config = 0;
		if (aa_side[1] == 0)	// no side chain, start with just one bit 
		{
			bitstring_config += turnseq_abs_top_all[2][i] / 2;		// turn direction 1 and 3 are q5=0 and q5=1, respectively
			for (k = 3; k < (AMINOACIDS - 1); k++)
			{
				// Least significant bit is the higher qubit (why qiskit, why?)
				temp_uint64 = turnseq_abs_top_all[k][i] % 2;
				bitstring_config += (temp_uint64 << (2 * (k - 2)));
				// and, the most significant bit is the lower qubit
				temp_uint64 = turnseq_abs_top_all[k][i] / 2;
				bitstring_config += (temp_uint64 << (2 * (k - 2)) - 1);
			}
		}
		else	// side chain, start with two bits
		{
			for (k = 2; k < (AMINOACIDS - 1); k++)
			{
				// LSB
				temp_uint64 = turnseq_abs_top_all[k][i] % 2;
				bitstring_config += (temp_uint64 << (2 * (k - 2)) + 1);
				// MSB
				temp_uint64 = turnseq_abs_top_all[k][i] / 2;
				bitstring_config += (temp_uint64 << (2 * (k - 2)));
			}
		}
		// print out configuration bitstring
		for (k = 0; k < qbits_config; k++)		// right to left (LSB to MSB)
		{
			bitshift = (qbits_config - 1 - k);
			//if (bitshift < 64)
			//{
				bitval = (bitstring_config & (one64u << bitshift)) ? 1 : 0;
			//}
			//else
			//{
			//	bitval = (bitstring_config & (one64u << (bitshift - 64))) ? 1 : 0;
			//}
			fprintf(f, "%d", bitval);
		}

		fprintf(f, "\n");

		// output the coordinates to a file
		sprintf(filename, "top_%d_e%.3f_coord.xyz", i, objtop_all[i]);
		g = fopen(filename, "w");
		fprintf(f, "%d\n\n", AMINOACIDS);
		fprintf(g, "%d\n\n", AMINOACIDS);

		// print out the coordinates
		xx = yy = zz = 0.;
		fprintf(f, "%c %f %f %f\n", aa_string[0], xx, yy, zz);
		fprintf(g, "%c %f %f %f\n", aa_string[0], xx, yy, zz);
		for (k = 0; k < (AMINOACIDS-1); k++)
		{
			/*
			    # Coordinates of the 4 edges of a tetrahedron centered at 0. The vectors are normalized.
				COORDINATES = (1.0 / np.sqrt(3)) * np.array(
					[[-1, 1, 1], [1, 1, -1], [-1, -1, -1], [1, -1, 1]]
				)
			*/
			switch (turnseq_abs_top_all[k][i])
			{
			case 0:
				if ((k % 2) == 0)	// even:  B lattice
				{
					xx -= 1. / sqrt(3.);
					yy += 1. / sqrt(3.);
					zz += 1. / sqrt(3.);
				}
				else	// odd:  A lattice
				{
					xx += 1. / sqrt(3.);
					yy -= 1. / sqrt(3.);
					zz -= 1. / sqrt(3.);
				}
				break;
			case 1:
				if ((k % 2) == 0)	// even:  B lattice
				{
					xx += 1. / sqrt(3.);
					yy += 1. / sqrt(3.);
					zz -= 1. / sqrt(3.);
				}
				else	// odd:  A lattice
				{
					xx -= 1. / sqrt(3.);
					yy -= 1. / sqrt(3.);
					zz += 1. / sqrt(3.);
				}
				break;
			case 2:
				if ((k % 2) == 0)	// even:  B lattice
				{
					xx -= 1. / sqrt(3.);
					yy -= 1. / sqrt(3.);
					zz -= 1. / sqrt(3.);
				}
				else	// odd:  A lattice
				{
					xx += 1. / sqrt(3.);
					yy += 1. / sqrt(3.);
					zz += 1. / sqrt(3.);
				}
				break;
			case 3:
				if ((k % 2) == 0)	// even:  B lattice
				{
					xx += 1. / sqrt(3.);
					yy -= 1. / sqrt(3.);
					zz += 1. / sqrt(3.);
				}
				else	// odd:  A lattice
				{
					xx -= 1. / sqrt(3.);
					yy += 1. / sqrt(3.);
					zz -= 1. / sqrt(3.);
				}
				break;
			}
			fprintf(f, "%c %f %f %f\n", aa_string[k + 1], xx* aa_scale, yy* aa_scale, zz* aa_scale);
			fprintf(g, "%c %f %f %f\n", aa_string[k + 1], xx* aa_scale, yy* aa_scale, zz* aa_scale);
		}
		fprintf(f, "\n\n");

		fclose(g);

	}
	fclose(f);




	return 0;
}
