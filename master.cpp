// ==================================================================================
// Name        : master.cpp
// Author      : Felix Baessler
// Version     : 11.04.2023
// Copyright   : Felix Baessler, felix.baessler@gmail.com
// SEE TLDR; VERSION OF THE LICENSE: https://creativecommons.org/licenses/by-nc/4.0/legalcode
// SEE FULL LICENSE DETAILS HERE   : https://creativecommons.org/licenses/by-nc/4.0/
//
// Description :
// Master generates on disk a long enough IID distributed byte sequence of minimum N= ns+NS bytes.
// The master disk file comprises the
// -	reference data set (ns bytes) concatenated to the
// -	test data set (NS bytes)
// Seamless shingling of the master file, requires that the reference data overlaps
// with the first L-1 bytes of the test data, so that
// -	the last reference shingle starts at position: ns - 1
// -	the first test shingle starts at position    : ns
// Because of this overlap, the test string appears lengthened by L-1 bytes in the present implementation.
// Run on an ordinary laptop, the execution of master takes about 1 minute for N= 2GB
// ==================================================================================

#include <time.h>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <random>
#include <cstring>
#include <chrono>
using namespace std;
using Time = std::chrono::time_point<std::chrono::high_resolution_clock>;

// maximum length of reference string s and test string S (ns + NS bytes)
#define N   2000000000ULL	//   2 GB
#define BUFFER_SIZE  32768ULL // 2 ^ 15

int
main()
{
	// master file name
	const string master_string_file_name= "C:\\cr\\master.txt";

	// attach output stream to master file
	ofstream output_stream(master_string_file_name, ios::binary);
	// check stream status
	if (!output_stream) {
		cerr << "Can't open master file!";
		exit(10);
	}

	printf("\n");
	printf("master_v0 11.04.23 \n");
	printf("========= \n");
	printf("maximum total string length: %llu  (reference + test string)\n", N);
	printf("buffer size: %llu  \n", BUFFER_SIZE);
	fflush(stdout);

	printf("\n");
	printf("storage allocation \n");
	// random byte string
	uint8_t *buffer;
	buffer= (uint8_t *)malloc(BUFFER_SIZE);
	if (buffer == NULL) exit(11);

	printf("\n");
	printf("generate and write master string \n");
	fflush(stdout);
	time_t  cur_time = time(NULL);		// current time
	mt19937 mt_rand(time(&cur_time));	// random number initialization

	// generate master string
	// ----------------------
    std::uniform_int_distribution<uint8_t> dist(0, 255);
    uint64_t i= 0;
	// N : total length
	for (uint64_t j= 0; j < N; j++ ) {
		buffer[i]= dist(mt_rand);
		if (++i == BUFFER_SIZE) {
			output_stream.write((char *)buffer, BUFFER_SIZE);
			i= 0;
			// printf("%x ... %x \n", buffer[0], buffer[BUFFER_SIZE - 1]);
		}
	}
	if (i > 0) {
		output_stream.write((char *)buffer, i);
		// printf("%x ... %x \n", buffer[0], buffer[i - 1]);
	}
	output_stream.close();
	printf("\n");
	printf("end \n");
}
