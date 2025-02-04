
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
// #define AMINOACIDS	5
// #define PEPTIDE "GNLVS"	// 4QXX
// #define AMINOACIDS	6
// #define PEPTIDE "SNQNNF" // 2OL9
// #define AMINOACIDS	7
// #define PEPTIDE "GSNQNNF"	// linear 7 AA peptide
//#define AMINOACIDS	7
//#define PEPTIDE "APRLRFY"	// from IBM paper
//#define AMINOACIDS	8
//#define PEPTIDE "GCVLYPWC"	// 2M6C from PDB
#define AMINOACIDS	8
#define PEPTIDE "YQFWKNFQ"   // middle 8 AAs from 2MZX
// #define AMINOACIDS	9
// #define PEPTIDE "DIRALKTLV"	// 1LVR from PDB
//#define AMINOACIDS	10
//#define PEPTIDE		"YYDPETGTWY"	// chignolin
//#define PEPTIDE		"DRVYVHPFHL"	// incorrect angiotensin sequence from Robert et al paper
//#define PEPTIDE		"DRVYIHPFHL"	// correct angiotensin sequence
// #define AMINOACIDS	12
// #define PEPTIDE		"FATMRYPSDSDE"	// 1IXU protease inhibitor - loop
//#define PEPTIDE		"VRRFDLLKRILK"	// 2N5R calmodulin binding peptide - helix + loop
// #define AMINOACIDS	13
// #define PEPTIDE		"IFGAIAGFIKNIW"	// 2L24 antimicrobial peptide - helix
// #define AMINOACIDS	14
//#define PEPTIDE		"RGKWTYNGITYEGR"	// mbh12 ph=5, T=283, beta sheet
// #define PEPTIDE		"IFGAIAGFIKNIWX"	// 2L24 antimicrobial peptide - helix
//#define AMINOACIDS	15
//#define PEPTIDE		"INWLKLGKKIIASLX"	// synoeca peptide
//#define PEPTIDE		"VLAMWKVGFFKRNRP"	// Jun's integrin fragment
//#define AMINOACIDS	16
//#define PEPTIDE		"RHYYKFNSTGRHYHYY"	// 8T61	designer peptide BH33
//#define PEPTIDE		"WHMWNTVPNAKQVIAA"	// 8T63	designer peptide PH1
//#define AMINOACIDS	20
//#define PEPTIDE		"DAYAQWLKDGGPSSGRPPPS"	// modified trp cage
//#define PEPTIDE		"NLYIQWLKDGGPSSGRPPPS"	// original Trp Cage
//#define AMINOACIDS	25
//#define PEPTIDE		"KKPGASLAALQALQALQAAQAAKKY"	// 25AA 8B1X helical peptide
//#define AMINOACIDS	26
//#define PEPTIDE		"YYHFWHRGVTKRSLSPHRPRHSRLQR"	// 26AA 6A8Y helix + turn
//#define PEPTIDE		"GNDYEDRYYRENMYRYPNQVYYRPVC"	// 26AA 1G04 sheep prion


#define TOP	50

// #define THREADS		8
#define MAXTHREADS		64


int main()
{
	int64_t index, permindex;
	uint64_t perm, perm_int, ilong, temp, temp2;
	int a, b, i, j, k, kk, threadnum, threads, bitval, p, q, ctr, bitshift, collision, backtrack;
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

	int xstep[AMINOACIDS - 1], ystep[AMINOACIDS - 1], zstep[AMINOACIDS - 1];
	int xstep_sum, ystep_sum, zstep_sum, direction_index;
	int step_lookup[36];

	char str1[64], str2[64], str3[2048], str4[2048];
	uint64_t y, y1, y2, temp_uint64;
	uint64_t bitstring, bitstring_config, bitstring_inter;
	FILE *f, *g;

	char aa_string[] = PEPTIDE;
	char aa_list[20];
	int aa_index[AMINOACIDS];
	double eij[20][20], aa_eij[AMINOACIDS][AMINOACIDS], xx, yy, zz;
	double aa_scale = 3.8;

	const uint64_t one64u = 1;



	// FCC lattice step lookup
	// This might not the most optimal ordering,
	// But importantly, the second half is the inverse of the first half
	step_lookup[0] = 1;		step_lookup[1] = 1;		step_lookup[2] = 0;
	step_lookup[3] = -1;	step_lookup[4] = 1;		step_lookup[5] = 0;
	step_lookup[6] = 0;		step_lookup[7] = 1;		step_lookup[8] = 1;
	step_lookup[9] = 0;		step_lookup[10] = 1;	step_lookup[11] = -1;
	step_lookup[12] = 1;	step_lookup[13] = 0;	step_lookup[14] = 1;
	step_lookup[15] = 1;	step_lookup[16] = 0;	step_lookup[17] = -1;
	step_lookup[18] = -1;	step_lookup[19] = -1;	step_lookup[20] = 0;
	step_lookup[21] = 1;	step_lookup[22] = -1;	step_lookup[23] = 0;
	step_lookup[24] = 0;	step_lookup[25] = -1;	step_lookup[26] = -1;
	step_lookup[27] = 0;	step_lookup[28] = -1;	step_lookup[29] = 1;
	step_lookup[30] = -1;	step_lookup[31] = 0;	step_lookup[32] = -1;
	step_lookup[33] = -1;	step_lookup[34] = 0;	step_lookup[35] = 1;



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
		aa_index[i] = -1;
		for (k = 0; k < 20; k++)
		{
			if (aa_string[i] == aa_list[k])
			{
				aa_index[i] = k;
			}
		}
		if (aa_index[i] < 0)
		{
			printf("amino acid %c not found, exiting...\n\n", aa_string[i]);
			exit(2);
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
	// threads = 48;
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
	// First step is defined to be 0
	// Second step has 4 directions, considering rotational symmetry
	// The remaining steps have 11 permutations (no backtracking)
	perm = 1;
	perm *= 4;
	for (i = 2; i < (AMINOACIDS-1); i++)
	{
		perm *= 11;
	}	// perm = 4*11^(N-3)
	




	//
#pragma omp parallel for schedule(dynamic,1), default(shared) \
	private(permindex,turnseq_abs,xstep,ystep,zstep,temp,backtrack,temp2,i,j,k,collision,obj,a,b) \
	private(turn_accum,d2,turn_accum_prior,turn_accum_next,threadnum) \
	private(xstep_sum,ystep_sum,zstep_sum, direction_index) \
	shared(perm,step_lookup,aa_side,aa_parent,aa_eij,objbest,objtop,permbest,permtop,turnseq_abs_best,turnseq_abs_top) \
	num_threads(threads)
	//
	for (permindex = 0; permindex < perm; permindex++)
	{

		//
		// first step is fixed
		turnseq_abs[0] = 0;
		xstep[0] = step_lookup[turnseq_abs[0] * 3 + 0];
		ystep[0] = step_lookup[turnseq_abs[0] * 3 + 1];
		zstep[0] = step_lookup[turnseq_abs[0] * 3 + 2];
		// initialize:
		temp = permindex;
		backtrack = 0;
		// other turns
		for (i = 1; i < (AMINOACIDS - 1); i++)
		{
			if (i == 1)		// second step has 4 unique permutations
			{
				temp2 = temp / 4;
				direction_index = temp - temp2 * 4;
				switch (direction_index)
				{
				case 0:
					turnseq_abs[i] = 0;
					break;
				case 1:
					turnseq_abs[i] = 1;
					break;
				case 2:
					turnseq_abs[i] = 4;
					break;
				case 3:
					turnseq_abs[i] = 10;
					break;
				}
			}
			else
			{
				temp2 = temp / 11;
				direction_index = temp - temp2 * 11;
				turnseq_abs[i] = (turnseq_abs[i-1] + direction_index + 7) % 12;
			}
			// !!! Check if this step is in the opposite direction as the prior
			// !!! If so, break out to next permutation
			//if (((turnseq_abs[i] % 6) == (turnseq_abs[i - 1] % 6)) && (turnseq_abs[i] != turnseq_abs[i - 1]))
			//{
			//	backtrack = 1;
			//	break;
			//}
			// look up the x, y, z steps for this step index
			xstep[i] = step_lookup[turnseq_abs[i] * 3 + 0];
			ystep[i] = step_lookup[turnseq_abs[i] * 3 + 1];
			zstep[i] = step_lookup[turnseq_abs[i] * 3 + 2];
			temp = temp2;
		}
		// turn_index[AMINOACIDS] only represents (AMINOACIDS-1) turns, but set the last value to zero (no turn)
		//turn_index[(AMINOACIDS - 1)] = 0;

		//if (backtrack)
		//	continue;

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
		collision = 0;
		obj = 0.;
		//
		// on the FCC lattice, 
		// a 1-NN interaction has to be an odd number of steps
		// a minimum of 2 steps are needed to be 1-NN apart
		// a collision has to be an even number of steps
		// a minimum of 3 steps are needed for a collision
		for (a = 0; a < (AMINOACIDS - 2); a++)
		{
			// if a collision had been found, break out of loop, and go to next turn sequence permutation
			if (collision)
			{
				obj = 10000.;
				break;
			}

			for (b = (a + 2); b < AMINOACIDS; b++)
			{

				// determine the main chain starting point i and main chain ending point j
				i = a;	//i = aa_parent[a];		// if considering side chains...
				j = b;	//j = aa_parent[b];

				// initialize the lattice step sums
				xstep_sum = ystep_sum = zstep_sum = 0;
				//
				// main chain:
				for (k = i; k < j; k++)
				{
					xstep_sum += xstep[k];
					ystep_sum += ystep[k];
					zstep_sum += zstep[k];
				}
				// side chain i:  N/A
				// side chain j:  N/A		


				// calculate distance
				//

				d2 = xstep_sum * xstep_sum + ystep_sum * ystep_sum + zstep_sum * zstep_sum;

				if (d2 == 0)	// collision
				{
					collision = 1;
					obj = 10000.;
					break;
				}

				if (d2 == 2)	// First nearest neighbor
				{
					// add the 1-NN interaction energy
					obj += aa_eij[i][j];
				}

				if (d2 == 4)    // Second nearest neighbor
				{

                	// add the 2-NN interaction energy, scale by 1/r

                	obj += 0.707106781 * aa_eij[i][j];

				}

			}
		}


		if (collision)
			continue;



#ifdef OMP
		threadnum = omp_get_thread_num();
#else
		threadnum = 0;
		// threadnum = 48;
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
			fprintf(f, "%x", turnseq_abs_top_all[AMINOACIDS - 2 - k][i]);
		}

		fprintf(f, "\n");

		// output the coordinates to a file
		// sprintf(filename, "top_%d_e%.3f_coord.xyz", i, objtop_all[i]);
		// g = fopen(filename, "w");
		fprintf(f, "%d\n\n", AMINOACIDS);
		// fprintf(g, "%d\n\n", AMINOACIDS);

		// print out the coordinates
		xx = yy = zz = 0.;
		fprintf(f, "%c %f %f %f\n", aa_string[0], xx, yy, zz);
		// fprintf(g, "%c %f %f %f\n", aa_string[0], xx, yy, zz);
		for (k = 0; k < (AMINOACIDS-1); k++)
		{
			xx += (double)step_lookup[turnseq_abs_top_all[k][i] * 3 + 0] / sqrt(2.);
			yy += (double)step_lookup[turnseq_abs_top_all[k][i] * 3 + 1] / sqrt(2.);
			zz += (double)step_lookup[turnseq_abs_top_all[k][i] * 3 + 2] / sqrt(2.);

			fprintf(f, "%c %f %f %f\n", aa_string[k + 1], xx* aa_scale, yy* aa_scale, zz* aa_scale);
			// fprintf(g, "%c %f %f %f\n", aa_string[k + 1], xx* aa_scale, yy* aa_scale, zz* aa_scale);
		}
		fprintf(f, "\n\n");

		// fclose(g);

	}
	fclose(f);




	return 0;
}
