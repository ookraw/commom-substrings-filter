// ================================================================================
// Name        : scatter_v1.cpp  (v1: 8 cofilters)
// Author      : Felix Baessler
// Version     : 25.04.2023
// Copyright   : Felix Baessler, felix.baessler@gmail.com
// SEE TLDR; VERSION OF THE LICENSE: https://creativecommons.org/licenses/by-nc/4.0/legalcode
// SEE FULL LICENSE DETAILS HERE   : https://creativecommons.org/licenses/by-nc/4.0/
//
// Description : Scatter
//  - reads the reference data from the master file (n= ns shingles),
//  - creates the fingerprint map and
//  - writes the result to the map file.
// The map can be viewed as a minimalistic hash table reduced to m one-bit slots.
// Scatter will mark those slots that correspond to the hash value modulo m
// of the reference shingles (reference fingerprint).
// Run on an ordinary laptop, scatter takes about 70 seconds for ns= 1GB.
//
// Demo-String:
// For demonstration, scatter smuggles an artificial 20 byte long demo-string
// into the reference data. Gather doing the same with the test data, this
// string will be detected as a common substring (search for "Demo-String" in
// the code).
//
// Compilation flags:
// -O3 -g3 -Wall         : optimization
// -Wl,--stack,0xFFFFFF  : long arrays
// if required: upgrade minGW to x86_64 : 64 bit executable
//
// Include files :
// - #include "mingw.thread.h"
// - #include "mingw.mutex.h"
// - #include "mingw.condition_variable.h"
// available from:
// Standard threads implementation currently still missing on MinGW GCC on Windows
// Author: Alexander Vassilev
// https://github.com/meganz/mingw-std-threads
// ================================================================================

#include <iostream>
#include <fstream>
#include <random>
#include "mingw.thread.h"
#include "mingw.mutex.h"
#include "mingw.condition_variable.h"
using namespace std;
// time measurement
using Time = std::chrono::time_point<std::chrono::high_resolution_clock>;
Time start_timer() {
    return std::chrono::high_resolution_clock::now();
}
double get_elapsed_time(Time start) {
    Time end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> d = end - start;
    std::chrono::microseconds us = std::chrono::duration_cast<std::chrono::microseconds>(d);
    return us.count() / 1000.0f;
}
// random cyclic permutations
void rcp_generator(std::mt19937& mt_rand, uint8_t p[]);

// file names
// ----------
// input:  read S from the master file
const string master_string_file_name= "C:\\cr\\master.txt";
// output: write the hash map
const string map_file_name_prefix=    "C:\\cr\\v1_map_";
// batch size of input buffers and hash values (cf. CONTAINERS)
#define BATCH_SIZE   (8*1024)

// GLOBAL PARAMETERS: same values in scatter and gather!
// =================

#define DV   8		// number of diversified hashes == 8 for this implementation (byte packed!!)
#define L    5		// shingle length L
#define LC  (L-1)	// shingle carry length : LC == L - 1
#define ns      1000000000ULL 		// length of the reference string s [bytes]
// n is lengthened by the first L-1 test bytes (overlap with reference/test string)
#define n		ns 					// number of reference shingles
#define M_COM   1000000007ULL		// modulus of the common hashes
#define B_COM   257ULL				// base of the common hashes (first prime > 256)
#define M_DIV   67ULL

// scatter_v1: diversified fingerprint bases
// -----------------------------------------
// 1 cache-line	 ( 64 bytes)
const uint64_t B_DIV[DV]= {257, 263, 269, 271, 277, 281, 283, 293};

// C_COM= (B_COM ^ L) % M_COM
auto c_com= [](){
		uint64_t result = 1ULL;
		for (int k = 0; k < L; k++) {
			result *= B_COM;
			result %= M_COM;
		}
		return result;
};
const uint64_t C_COM= c_com();

// C_DIV= (B_DIV[id] ^ L) % M_DIV
auto c_div= [](int id){
		uint64_t result = 1ULL;
		for (int k = 0; k < L; k++) {
			result *= B_DIV[id];
			result %= M_DIV;
		}
		return result;
};
const uint64_t C_DIV[DV]= {c_div(0), c_div(1), c_div(2), c_div(3),
						   c_div(4), c_div(5), c_div(6), c_div(7)};

// THREADS (1,2,3) working on three containers (A,B,C)
// =======
// WORKER 1 : read input
void worker1_thread();
// fill the batch from the reference file
mutex mx1;
condition_variable cv1;
bool cv1_scheduler_enabled= false;
bool cv1_worker1_enabled= false;
double worker1_process_time;
double worker1_waiting_time;
// WORKER 2 : produce hash
void worker2_thread();
// hash the shingles of the batch
void hash_batch(uint8_t s[], uint32_t hash_count, uint64_t com_hash[], uint8_t div_hash[]);
mutex mx2;
condition_variable cv2;
bool cv2_scheduler_enabled= false;
bool cv2_worker2_enabled= false;
double worker2_process_time;
double worker2_waiting_time;
// WORKER 3 : consume hash
void worker3_thread();
// record in the hash map the hash values of the batch
mutex mx3;
condition_variable cv3;
bool cv3_scheduler_enabled= false;
bool cv3_worker3_enabled= false;
double worker3_process_time;
double worker3_waiting_time;

// CONTAINERS (A,B,C)
// ==================
// CONTAINER A
bool ctr_A_busy= false;					// true: container A busy (a worker is processing its contents)
uint8_t  buf_ctr_A[BATCH_SIZE + LC];	// string buffer: length= buffer size + carry length
uint64_t com_ctr_A[BATCH_SIZE];			// hash buffer  : common hashes
uint8_t  div_ctr_A[BATCH_SIZE * DV];	// hash buffer  : diversified hashes (DV number of cofilters)
// CONTAINER B
bool ctr_B_busy= false;					// true: container B busy
uint8_t  buf_ctr_B[BATCH_SIZE + LC];
uint64_t com_ctr_B[BATCH_SIZE];
uint8_t  div_ctr_B[BATCH_SIZE * DV];
// CONTAINER C
bool ctr_C_busy= false;					// true: container C busy
uint8_t  buf_ctr_C[BATCH_SIZE + LC];
uint64_t com_ctr_C[BATCH_SIZE];
uint8_t  div_ctr_C[BATCH_SIZE * DV];

// THREAD INTERFACE
// ================
// hash map
uint8_t *map;
// cyclic permutation vector
uint8_t shuffle[256];

// ****************************************************************************************************************************

int main() {
	uint32_t stage_id;			// current stage
	uint32_t batch_count;		// total number of required batches (>= 3)
	// compatible map file: same M_COM, M_DIV and L in scatter and gather!!
	string map_file_name= map_file_name_prefix
			+ to_string(M_DIV) + "_"
			+ to_string(L) + ".txt";

	// total number of batches
	// -----------------------
	// note that the last batch is not necessarily completely full
	batch_count= n / BATCH_SIZE;
	if ((n - (batch_count * BATCH_SIZE)) > 0) batch_count++;
	if (batch_count < 3) exit(10);

	printf("\n");
	printf("scatter_v1 \n");
	printf("========== \n");
	printf("master file           : %s \n", master_string_file_name.c_str() );
	printf("map    file           : %s \n", map_file_name.c_str());
	printf("string  s length ns   : %llu \t(reference string) \n", ns);
	printf("shingle length L      : %d \n", L);
	printf("carry   length LC     : %d \n", LC);
	printf("batch count           : %d \n", batch_count);
	printf("batch size            : %d \n", BATCH_SIZE);
	printf("common modulus        : %llu \n", M_COM);
	printf("diversity modulus     : %llu \n", M_DIV);
	printf("\n");
	fflush(stdout);

	// reset time expenditures
    worker1_waiting_time= 0;
    worker1_process_time= 0;
	worker2_waiting_time= 0;
    worker2_process_time= 0;
    worker3_waiting_time= 0;
    worker3_process_time= 0;

	// hash map allocation / reset
	// ---------------------------
	map= (uint8_t *)malloc(M_COM + M_DIV); if (map == NULL) exit(11);
	// reset hash map
	memset(map, 0b11111111, M_COM + M_DIV);

	// random number initialization with current time
	// ==============================================
	//mt19937 mt_rand(123456789);
	time_t  cur_time;
	mt19937 mt_rand(time(&cur_time));
	//printf("map setup_time :  %s \n", ctime(&cur_time));

	// generate random cyclic permutations: shuffle
	// -----------------------------------
    rcp_generator(mt_rand, shuffle);
	// printf("first 20 shuffle values: \n");
	// for (uint32_t i= 0; i<20; i++) printf(" %d ", shuffle[i]);
	// printf("\n");

	// elapsed time
	double elapsed_time= 0;
	Time start_elapsed_time= start_timer();
	// work time
	double work_time= 0;
	Time start_work_time;
	// schedule time
	double schedule_time= 0;
	Time start_schedule_time;
	// overhead time
	double overhead_time= 0;
	Time start_overhead_time;

	// start threads (they will wait for the first start-signal)
	start_overhead_time= start_timer();
    thread worker1(worker1_thread);
    thread worker2(worker2_thread);
    thread worker3(worker3_thread);
    overhead_time+= get_elapsed_time(start_overhead_time);

    // first stage: 1A (worker1 processes container A)
	{{lock_guard<mutex> lk(mx1); cv1_worker1_enabled= true;} cv1.notify_one();}
	{unique_lock<mutex> lk(mx1); cv1.wait(lk, []{return cv1_scheduler_enabled;}); cv1_scheduler_enabled= false;}
	// second stage: 1B, 2A
	{{lock_guard<mutex> lk(mx1); cv1_worker1_enabled= true;} cv1.notify_one();}
	{{lock_guard<mutex> lk(mx2); cv2_worker2_enabled= true;} cv2.notify_one();}
	{unique_lock<mutex> lk(mx1); cv1.wait(lk, []{return cv1_scheduler_enabled;}); cv1_scheduler_enabled= false;}
	{unique_lock<mutex> lk(mx2); cv2.wait(lk, []{return cv2_scheduler_enabled;}); cv2_scheduler_enabled= false;}
	// third stage: 1C, 2B, 3A
	// etc...
	for (stage_id= 3; stage_id <= batch_count; stage_id++) {
		//p cout << "send start-signal to workers \n"; fflush(stdout);
		// <== send start-signal to workers
		start_schedule_time= start_timer();
		{{lock_guard<mutex> lk(mx1); cv1_worker1_enabled= true;} cv1.notify_one();}
		{{lock_guard<mutex> lk(mx2); cv2_worker2_enabled= true;} cv2.notify_one();}
		{{lock_guard<mutex> lk(mx3); cv3_worker3_enabled= true;} cv3.notify_one();}
		schedule_time+= get_elapsed_time(start_schedule_time);
		//p cout << "wait for end-signals from workers \n"; fflush(stdout);
		// <== wait for end-signal from workers
		start_work_time= start_timer();
		{unique_lock<mutex> lk(mx1); cv1.wait(lk, []{return cv1_scheduler_enabled;}); cv1_scheduler_enabled= false;}
		{unique_lock<mutex> lk(mx2); cv2.wait(lk, []{return cv2_scheduler_enabled;}); cv2_scheduler_enabled= false;}
		{unique_lock<mutex> lk(mx3); cv3.wait(lk, []{return cv3_scheduler_enabled;}); cv3_scheduler_enabled= false;}
		work_time+= get_elapsed_time(start_work_time);
	}
	// second last stage (batch_count + 1): 2x, 3x
	{{lock_guard<mutex> lk(mx2); cv2_worker2_enabled= true;} cv2.notify_one();}
	{{lock_guard<mutex> lk(mx3); cv3_worker3_enabled= true;} cv3.notify_one();}
	{unique_lock<mutex> lk(mx2); cv2.wait(lk, []{return cv2_scheduler_enabled;}); cv2_scheduler_enabled= false;}
	{unique_lock<mutex> lk(mx3); cv3.wait(lk, []{return cv3_scheduler_enabled;}); cv3_scheduler_enabled= false;}
	// last stage (batch_count + 2): 3y
	{{lock_guard<mutex> lk(mx3); cv3_worker3_enabled= true;} cv3.notify_one();}
	{unique_lock<mutex> lk(mx3); cv3.wait(lk, []{return cv3_scheduler_enabled;}); cv3_scheduler_enabled= false;}

	// end threads
	start_overhead_time= start_timer();
    worker1.join();
    worker2.join();
    worker3.join();
    overhead_time+= get_elapsed_time(start_overhead_time);
	elapsed_time= get_elapsed_time(start_elapsed_time);

	// result
	// ======
	// write hash map to disk
	// ----------------------
	ofstream map_output_stream(map_file_name, ios::binary);
	if (!map_output_stream) cerr << "Can't open map output file!";
	map_output_stream.write((char *)&cur_time, (int)sizeof(time_t));
	map_output_stream.write((char *)map, M_COM + M_DIV);
	map_output_stream.close();
	printf("\nmap setup_time :  %s \n", ctime(&cur_time));
	// printf("first 20 map values: \n");
	// for (uint32_t i= 0; i<20; i++) printf(" %d ", map[i]);
	// printf("\n");

	// time expenditure
	// ================
	printf("\n");
	printf("time expenditure [milliseconds]\n");
	printf("---------------- \n");
	printf("elapsed     : %9.0f  \n", elapsed_time);
	printf("work        : %9.0f  \n", work_time);
	printf("schedule    : %9.0f  \n", schedule_time);
	printf("overhead    : %9.0f  \n", overhead_time);
	printf("worker1     : %9.0f  \n", worker1_waiting_time + worker1_process_time);
    printf(" - wait     : %9.0f  \n", worker1_waiting_time);
    printf(" - process  : %9.0f  \n", worker1_process_time);
	printf("worker2     : %9.0f  \n", worker2_waiting_time + worker2_process_time);
    printf(" - wait     : %9.0f  \n", worker2_waiting_time);
    printf(" - process  : %9.0f  \n", worker2_process_time);
	printf("worker3     : %9.0f  \n", worker3_waiting_time + worker3_process_time);
    printf(" - wait     : %9.0f  \n", worker3_waiting_time);
    printf(" - process  : %9.0f  \n", worker3_process_time);
	fflush(stdout);

}

// ****************************************************************************************************************************

void worker1_thread()
{
	SetThreadAffinityMask(GetCurrentThread(), 7ULL);
	Time start_time;			// start of time measurement
	uint8_t *buffer;			// byte string buffer
	uint8_t *input_buffer;		// start at: buffer + LC (carry length)
	uint8_t *carry_buffer;		// start at: buffer + BATCH_SIZE (buffer size)

	uint32_t batch_id;			// current batch
	uint32_t batch_size;		// current batch size
	uint32_t last_batch_size;	// batch size of the last batch
	uint32_t batch_count;		// total number of required batches
	uint32_t demo_batch_id;		// "Demo-String"

	batch_count= n / BATCH_SIZE;
	demo_batch_id= batch_count / 2;
	last_batch_size= n - (batch_count * BATCH_SIZE);
	if (last_batch_size > 0) batch_count++;
	else last_batch_size= BATCH_SIZE;
	if (batch_count == 1) {
		// the first batch is also the last one
		batch_size= last_batch_size;
	} else {
		// for full batches:
		batch_size= BATCH_SIZE;
	}

	// attach input stream to master file
	ifstream string_input_stream(master_string_file_name, ios::in|ios::binary|ios::ate);
	// check stream status
	if (!string_input_stream) cerr << "Can't open master file!";

	// get/check length of master file
	string_input_stream.seekg (0, string_input_stream.end);
	// check file size:
	if ((uint64_t)string_input_stream.tellg() < ns) {
		printf("master file length < ns : %llu \n", ns);
		fflush(stdout);
		exit(12);
	}
	fflush(stdout);

	// set first artificial carry (skipped by worker3 -> first_batch_j0= LC)
	// ---------------------------------------------------------------------
	// set position of the input stream at the begin of S0
	string_input_stream.seekg (0, string_input_stream.beg);
	// locate current carry buffer
	carry_buffer= buf_ctr_C + BATCH_SIZE;
	// fill carry with dummy values
	for (uint32_t j= 0; j < LC; j++) {
		carry_buffer[j]= 0;
	}

	batch_id= 1;
	while (true) {	// as long as batch_id <= batch_count

	  // worker1 container A
	  // ===================
	  //p cout << "worker1 waits for start-signal from scheduler \n"; fflush(stdout);
	  // <== worker1 waits for start-signal from scheduler
	  start_time= start_timer();
	  {unique_lock<mutex> lk(mx1); cv1.wait(lk, []{return cv1_worker1_enabled;}); cv1_worker1_enabled= false;}
	  worker1_waiting_time+= get_elapsed_time(start_time);
	  start_time= start_timer();

	  //p cout << "worker1 processes container A \n"; fflush(stdout);
	  //p cout << "2a "; fflush(stdout);
	  if (ctr_A_busy) exit(13); else ctr_A_busy= true;
	  // ***********************************************************
	  // worker1 produces/processes the current batch in container A
	  // ***********************************************************
	  // define current buffer
	  buffer= buf_ctr_A;
	  // move carry from the previous buffer (buf_ctr_C) to the beginning of the current buffer
	  for (uint32_t i= 0; i < LC; i++) buffer[i]= carry_buffer[i];
	  // define current carry buffer
	  carry_buffer= buffer + BATCH_SIZE;
	  // define current input buffer
	  input_buffer= buffer + LC;
	  // fill input buffer
	  string_input_stream.read((char *)input_buffer, batch_size);
	  // check number of bytes read
	  if (batch_size != string_input_stream.gcount()) exit(14);
	  // shuffle
	  for (uint32_t j= 0; j < batch_size; j++) {
		  input_buffer[j]= shuffle[input_buffer[j]];
	  }

///*
	  // "Demo-String"
	  if (batch_id == demo_batch_id) {
		  for (uint32_t j= 0; j < 20; j++) {
			  input_buffer[j]= 0;
		  }
	  }
//*/

	  // this_thread::sleep_for(chrono::milliseconds(200));
	  //shuffle cout << "worker1 completed container A \n"; fflush(stdout);
	  //p cout << "2A "; fflush(stdout);
	  ctr_A_busy= false;

	  //p cout << "worker1 sends end-signal to scheduler \n"; fflush(stdout);
	  // ==> worker1 sends end-signal to scheduler
	  {lock_guard<mutex> lk(mx1); cv1_scheduler_enabled= true;} cv1.notify_one();

	  worker1_process_time+= get_elapsed_time(start_time);

	  if (batch_id == batch_count) {
		  // the current batch was the last one
		  cout << "worker1 terminates on A \n"; fflush(stdout);
		  return;
	  }
	  if (++batch_id == batch_count) {
		  // the next batch is the last one
		  // the last container is not necessarily completely full
		  //p cout << "worker1 last batch on B \n"; fflush(stdout);
		  batch_size= last_batch_size;
	  }

	  // worker1 container B
	  // ===================
	  //p cout << "worker1 waits for start-signal from scheduler \n"; fflush(stdout);
	  // <== worker1 waits for start-signal from scheduler
	  start_time= start_timer();
	  {unique_lock<mutex> lk(mx1); cv1.wait(lk, []{return cv1_worker1_enabled;}); cv1_worker1_enabled= false;}
	  worker1_waiting_time+= get_elapsed_time(start_time);
	  start_time= start_timer();

	  //p cout << "worker1 processes container B \n"; fflush(stdout);
	  //p cout << "2b "; fflush(stdout);
	  if (ctr_B_busy) exit(15); else ctr_B_busy= true;
	  // ***********************************************************
	  // worker1 produces/processes the current batch in container B
	  // ***********************************************************
	  // define current buffer
	  buffer= buf_ctr_B;
	  // move carry from the previous buffer (buf_ctr_A) to the beginning of the current buffer
	  for (uint32_t i= 0; i < LC; i++) buffer[i]= carry_buffer[i];
	  // define current carry buffer
	  carry_buffer= buffer + BATCH_SIZE;
	  // define current input buffer
	  input_buffer= buffer + LC;
	  // fill input buffer
	  string_input_stream.read((char *)input_buffer, batch_size);
	  // check number of bytes read
	  if (batch_size != string_input_stream.gcount()) exit(16);
	  // shuffle
	  for (uint32_t j= 0; j < batch_size; j++) {
		  input_buffer[j]= shuffle[input_buffer[j]];
	  }

///*
	  // "Demo-String"
	  if (batch_id == demo_batch_id) {
		  for (uint32_t j= 0; j < 20; j++) {
			  input_buffer[j]= 0;
		  }
	  }
//*/

	  // this_thread::sleep_for(chrono::milliseconds(200));
	  //shuffle cout << "worker1 completed container B \n"; fflush(stdout);
	  //p cout << "2B "; fflush(stdout);
	  ctr_B_busy= false;

	  //p cout << "worker1 sends end-signal to scheduler \n"; fflush(stdout);
	  // ==> worker1 sends end-signal to scheduler
	  {lock_guard<mutex> lk(mx1); cv1_scheduler_enabled= true;} cv1.notify_one();

	  worker1_process_time+= get_elapsed_time(start_time);

	  if (batch_id == batch_count) {
		  // the current batch was the last one
		  cout << "worker1 terminates on B \n"; fflush(stdout);
		  return;
	  }
	  if (++batch_id == batch_count) {
		  // the next batch is the last one
		  // the last container is not necessarily completely full
		  //p cout << "worker1 last batch on C \n"; fflush(stdout);
		  batch_size= last_batch_size;
	  }

	  // worker1 container C
	  // ===================
	  //p cout << "worker1 waits for start-signal from scheduler \n"; fflush(stdout);
	  // <== worker1 waits for start-signal from scheduler
	  start_time= start_timer();
	  {unique_lock<mutex> lk(mx1); cv1.wait(lk, []{return cv1_worker1_enabled;}); cv1_worker1_enabled= false;}
	  worker1_waiting_time+= get_elapsed_time(start_time);
	  start_time= start_timer();

	  //p cout << "worker1 processes container C \n"; fflush(stdout);
	  //p cout << "2c "; fflush(stdout);
	  if (ctr_C_busy) exit(17); else ctr_C_busy= true;
	  // ***********************************************************
	  // worker1 produces/processes the current batch in container C
	  // ***********************************************************
	  // define current buffer
	  buffer= buf_ctr_C;
	  // move carry from the previous buffer (buf_ctr_B) to the beginning of the current buffer
	  for (uint32_t i= 0; i < LC; i++) buffer[i]= carry_buffer[i];
	  // define current carry buffer
	  carry_buffer= buffer + BATCH_SIZE;
	  // define current input buffer
	  input_buffer= buffer + LC;
	  // fill input buffer
	  string_input_stream.read((char *)input_buffer, batch_size);
	  // check number of bytes read
	  if (batch_size != string_input_stream.gcount()) exit(18);
	  // shuffle
	  for (uint32_t j= 0; j < batch_size; j++) {
		  input_buffer[j]= shuffle[input_buffer[j]];
	  }

///*
	  // "Demo-String"
	  if (batch_id == demo_batch_id) {
		  for (uint32_t j= 0; j < 20; j++) {
			  input_buffer[j]= 0;
		  }
	  }
//*/
	  // this_thread::sleep_for(chrono::milliseconds(200));
	  //p cout << "worker1 completed container C \n"; fflush(stdout);
	  //p cout << "2C "; fflush(stdout);
	  ctr_C_busy= false;

	  //p cout << "worker1 sends end-signal to scheduler \n"; fflush(stdout);
	  // ==> worker1 sends end-signal to scheduler
	  {lock_guard<mutex> lk(mx1); cv1_scheduler_enabled= true;} cv1.notify_one();

	  worker1_process_time+= get_elapsed_time(start_time);

	  if (batch_id == batch_count) {
		  // the current batch was the last one
		  cout << "worker1 terminates on C \n"; fflush(stdout);
		  return;
	  }
	  if (++batch_id == batch_count) {
		  // the next batch is the last one
		  // the last container is not necessarily completely full
		  //p cout << "worker1 last batch on A \n"; fflush(stdout);
		  batch_size= last_batch_size;
	  }
	}
}

void worker2_thread()
{
	SetThreadAffinityMask(GetCurrentThread(), 7ULL);
	Time start_time;			// start of time measurement
	uint32_t batch_id;			// current batch
	uint32_t batch_size;		// current batch size
	uint32_t last_batch_size;	// batch size of the last batch
	uint32_t batch_count;		// total number of required batches

	batch_count= n / BATCH_SIZE;
	last_batch_size= n - (batch_count * BATCH_SIZE);
	if (last_batch_size > 0) batch_count++;
	else last_batch_size= BATCH_SIZE;
	if (batch_count == 1) {
		// the first batch is also the last one
		batch_size= last_batch_size;
	} else {
		// for full batches:
		batch_size= BATCH_SIZE;
	}

  batch_id= 1;
  while (true) {	// as long as batch_id <= batch_count

    // worker2 container A
    // ===================
	//p cout << "worker2 waits for start-signal from scheduler \n"; fflush(stdout);
	// <== worker2 waits for start-signal from scheduler
	start_time= start_timer();
	{unique_lock<mutex> lk(mx2); cv2.wait(lk, []{return cv2_worker2_enabled;}); cv2_worker2_enabled= false;}
	worker2_waiting_time+= get_elapsed_time(start_time);
	start_time= start_timer();

		//p cout << "worker2 processes container A \n"; fflush(stdout);
		//p cout << "1a "; fflush(stdout);
		if (ctr_A_busy) exit(19); else ctr_A_busy= true;
		// ***********************************************************
		// worker2 produces/processes the current batch in container A
		// ***********************************************************
		hash_batch(buf_ctr_A, batch_size, com_ctr_A, div_ctr_A);
		// this_thread::sleep_for(chrono::milliseconds(100));
		//p cout << "worker2 completed container A \n"; fflush(stdout);
		//p cout << "1A "; fflush(stdout);
		ctr_A_busy= false;

		//p cout << "worker2 sends end-signal to scheduler \n"; fflush(stdout);
		// ==> worker2 sends end-signal to scheduler
		{lock_guard<mutex> lk(mx2); cv2_scheduler_enabled= true;} cv2.notify_one();

	worker2_process_time+= get_elapsed_time(start_time);
	if (batch_id == batch_count) {
		// the current batch was the last one
		cout << "worker2 terminates on A \n"; fflush(stdout);
		return;
	}
	if (++batch_id == batch_count) {
		// the next batch is the last one
		// the last container is not necessarily completely full
		//p cout << "worker2 last batch on B \n"; fflush(stdout);
		batch_size= last_batch_size;
	}

    // worker2 container B
    // ===================
	//p cout << "worker2 waits for start-signal from scheduler \n"; fflush(stdout);
	// <== worker2 waits for start-signal from scheduler
	start_time= start_timer();
	{unique_lock<mutex> lk(mx2); cv2.wait(lk, []{return cv2_worker2_enabled;}); cv2_worker2_enabled= false;}
	worker2_waiting_time+= get_elapsed_time(start_time);
	start_time= start_timer();

		//p cout << "worker2 processes container B \n"; fflush(stdout);
		//p cout << "1b "; fflush(stdout);
		if (ctr_B_busy) exit(20); else ctr_B_busy= true;
		// ***********************************************************
		// worker2 produces/processes the current batch in container B
		// ***********************************************************
		hash_batch(buf_ctr_B, batch_size, com_ctr_B, div_ctr_B);
		// this_thread::sleep_for(chrono::milliseconds(100));
		//p cout << "worker2 completed container B \n"; fflush(stdout);
		//p cout << "1B "; fflush(stdout);
		ctr_B_busy= false;

		//p cout << "worker2 sends end-signal to scheduler \n"; fflush(stdout);
		// ==> worker2 sends end-signal to scheduler
		{lock_guard<mutex> lk(mx2); cv2_scheduler_enabled= true;} cv2.notify_one();

	worker2_process_time+= get_elapsed_time(start_time);
	if (batch_id == batch_count) {
		// the current batch was the last one
		cout << "worker2 terminates on B \n"; fflush(stdout);
		return;
	}
	if (++batch_id == batch_count) {
		// the next batch is the last one
		// the last container is not necessarily completely full
		//p cout << "worker2 last batch on C \n"; fflush(stdout);
		batch_size= last_batch_size;
	}

    // worker2 container C
    // ===================
	//p cout << "worker2 waits for start-signal from scheduler \n"; fflush(stdout);
	// <== worker2 waits for start-signal from scheduler
	start_time= start_timer();
	{unique_lock<mutex> lk(mx2); cv2.wait(lk, []{return cv2_worker2_enabled;}); cv2_worker2_enabled= false;}
	worker2_waiting_time+= get_elapsed_time(start_time);
	start_time= start_timer();

		//p cout << "worker2 processes container C \n"; fflush(stdout);
		//p cout << "1c "; fflush(stdout);
		if (ctr_C_busy) exit(21); else ctr_C_busy= true;
		// ***********************************************************
		// worker2 produces/processes the current batch in container C
		// ***********************************************************
		hash_batch(buf_ctr_C, batch_size, com_ctr_C, div_ctr_C);
		// this_thread::sleep_for(chrono::milliseconds(100));
		//p cout << "worker2 completed container C \n"; fflush(stdout);
		//p cout << "1C "; fflush(stdout);
		ctr_C_busy= false;

		//p cout << "worker2 sends end-signal to scheduler \n"; fflush(stdout);
		// ==> worker2 sends end-signal to scheduler
		{lock_guard<mutex> lk(mx2); cv2_scheduler_enabled= true;} cv2.notify_one();

	worker2_process_time+= get_elapsed_time(start_time);
  if (batch_id == batch_count) {
	  // the current batch was the last one
	  cout << "worker2 terminates on C \n"; fflush(stdout);
	  return;
  }
  if (++batch_id == batch_count) {
	  // the next batch is the last one
	  // the last container is not necessarily completely full
	  //p cout << "worker2 last batch on A \n"; fflush(stdout);
	  batch_size= last_batch_size;
  }
  }
  //p cout << "end thread 1 \n"; fflush(stdout);
}

void worker3_thread()
{
	SetThreadAffinityMask(GetCurrentThread(), 8ULL);
	Time start_time;			// start of time measurement
	uint32_t batch_id;			// current batch
	uint32_t first_batch_j0= LC;// skip first (artificial) carry  <-------------
	uint32_t batch_size;		// current batch size
	uint32_t last_batch_size;	// batch size of the last batch
	uint32_t batch_count;		// total number of required batches

	batch_count= n / BATCH_SIZE;
	last_batch_size= n - (batch_count * BATCH_SIZE);
	if (last_batch_size > 0) batch_count++;
	else last_batch_size= BATCH_SIZE;
	if (batch_count == 1) {
		// the first batch is also the last one
		batch_size= last_batch_size;
	} else {
		// for full batches:
		batch_size= BATCH_SIZE;
	}

  batch_id= 1;
  while (true) {	// as long as batch_id <= batch_count

    // worker3 container A
    // ===================
	//p cout << "worker3 waits for start-signal from scheduler \n"; fflush(stdout);
	// <== worker3 waits for start-signal from scheduler
	start_time= start_timer();
	{unique_lock<mutex> lk(mx3); cv3.wait(lk, []{return cv3_worker3_enabled;}); cv3_worker3_enabled= false;}
	worker3_waiting_time+= get_elapsed_time(start_time);
	start_time= start_timer();

		//p cout << "worker3 processes container A \n"; fflush(stdout);
		//p cout << "3a "; fflush(stdout);
		if (ctr_A_busy) exit(22); else ctr_A_busy= true;
		// ***********************************************************
		// worker3 produces/processes the current batch in container A
		// ***********************************************************
		// record in the hash map the hash values of the shingle j
		for (uint32_t j= first_batch_j0; j < batch_size; j++) {
			// if it doesn’t help, it doesn’t hurt?
			for (uint64_t i= 0; i < ((M_DIV + 32) / 64); i++) {
				__builtin_prefetch (&com_ctr_A[j] + i * 64, 1, 3);
			}
			// keep track of the hash occurrence (TIME CRITICAL)
			for (uint8_t id= 0; id < DV; id++) {
				// ClearBit
				map[com_ctr_A[j] + div_ctr_A[j*DV+id]] &= ~(1<<id);
			}
		}
		// skip only the first (artificial) carry
		first_batch_j0= 0;

		// this_thread::sleep_for(chrono::milliseconds(300));
		//p cout << "worker3 completed container A \n"; fflush(stdout);
		//p cout << "3A "; fflush(stdout);
		ctr_A_busy= false;

		//p cout << "worker3 sends end-signal to scheduler \n"; fflush(stdout);
		// ==> worker3 sends end-signal to scheduler
		{lock_guard<mutex> lk(mx3); cv3_scheduler_enabled= true;} cv3.notify_one();

	worker3_process_time+= get_elapsed_time(start_time);
	if (batch_id == batch_count) {
		// the current batch was the last one
		cout << "worker3 terminates on A \n"; fflush(stdout);
		return;
	}
	if (++batch_id == batch_count) {
		// the next batch is the last one
		// the last container is not necessarily completely full
		//p cout << "worker3 last batch on B \n"; fflush(stdout);
		batch_size= last_batch_size;
	}

    // worker3 container B
    // ===================
	//p cout << "worker3 waits for start-signal from scheduler \n"; fflush(stdout);
	// <== worker3 waits for start-signal from scheduler
	start_time= start_timer();
	{unique_lock<mutex> lk(mx3); cv3.wait(lk, []{return cv3_worker3_enabled;}); cv3_worker3_enabled= false;}
	worker3_waiting_time+= get_elapsed_time(start_time);
	start_time= start_timer();

		//p cout << "worker3 processes container B \n"; fflush(stdout);
		//p cout << "3b "; fflush(stdout);
		if (ctr_B_busy) exit(23); else ctr_B_busy= true;
		// ***********************************************************
		// worker3 produces/processes the current batch in container B
		// ***********************************************************
		// record in the hash map the hash values of the shingle j
		for (uint32_t j= 0; j < batch_size; j++) {
			// if it doesn’t help, it doesn’t hurt?
			for (uint64_t i= 0; i < ((M_DIV + 32) / 64); i++) {
				__builtin_prefetch (&com_ctr_B[j] + i * 64, 1, 3);
			}
			// keep track of the hash occurrence (TIME CRITICAL)
			for (uint8_t id= 0; id < DV; id++) {
				// ClearBit
				map[com_ctr_B[j] + div_ctr_B[j*DV+id]] &= ~(1<<id);
			}
		}

		// this_thread::sleep_for(chrono::milliseconds(300));
		//p cout << "worker3 completed container B \n"; fflush(stdout);
		//p cout << "3B "; fflush(stdout);
		ctr_B_busy= false;

		//p cout << "worker3 sends end-signal to scheduler \n"; fflush(stdout);
		// ==> worker3 sends end-signal to scheduler
		{lock_guard<mutex> lk(mx3); cv3_scheduler_enabled= true;} cv3.notify_one();

	worker3_process_time+= get_elapsed_time(start_time);
	if (batch_id == batch_count) {
		// the current batch was the last one
		cout << "worker3 terminates on B \n"; fflush(stdout);
		return;
	}
	if (++batch_id == batch_count) {
		// the next batch is the last one
		// the last container is not necessarily completely full
		//p cout << "worker3 last batch on C \n"; fflush(stdout);
		batch_size= last_batch_size;
	}

    // worker3 container C
    // ===================
	//p cout << "worker3 waits for start-signal from scheduler \n"; fflush(stdout);
	// <== worker3 waits for start-signal from scheduler
	start_time= start_timer();
	{unique_lock<mutex> lk(mx3); cv3.wait(lk, []{return cv3_worker3_enabled;}); cv3_worker3_enabled= false;}
	worker3_waiting_time+= get_elapsed_time(start_time);
	start_time= start_timer();

		//p cout << "worker3 processes container C \n"; fflush(stdout);
		//p cout << "3c "; fflush(stdout);
		if (ctr_C_busy) exit(24); else ctr_C_busy= true;
		// ***********************************************************
		// worker3 produces/processes the current batch in container C
		// ***********************************************************
		// record in the hash map the hash values of the shingle j
		for (uint32_t j= 0; j < batch_size; j++) {
			// if it doesn’t help, it doesn’t hurt?
			for (uint64_t i= 0; i < ((M_DIV + 32) / 64); i++) {
				__builtin_prefetch (&com_ctr_C[j] + i * 64, 1, 3);
			}
			// keep track of the hash occurrence (TIME CRITICAL)
			for (uint8_t id= 0; id < DV; id++) {
				// ClearBit
				map[com_ctr_C[j] + div_ctr_C[j*DV+id]] &= ~(1<<id);
			}
		}

		// this_thread::sleep_for(chrono::milliseconds(300));
		//p cout << "worker3 completed container C \n"; fflush(stdout);
		//p cout << "3C "; fflush(stdout);
		ctr_C_busy= false;

		//p cout << "worker3 sends end-signal to scheduler \n"; fflush(stdout);
		// ==> worker3 sends end-signal to scheduler
		{lock_guard<mutex> lk(mx3); cv3_scheduler_enabled= true;} cv3.notify_one();

	worker3_process_time+= get_elapsed_time(start_time);
	if (batch_id == batch_count) {
		// the current batch was the last one
		cout << "worker3 terminates on C \n"; fflush(stdout);
		return;
	}
	if (++batch_id == batch_count) {
		// the next batch is the last one
		// the last container is not necessarily completely full
		//p cout << "worker3 last batch on A \n"; fflush(stdout);
		batch_size= last_batch_size;
	}
  }
  //p cout << "end thread 3 \n"; fflush(stdout);
}

// ****************************************************************************************************************************

inline void update_div_hashes(uint8_t hash[], uint8_t hash1[], const uint8_t s[]) {
	for (uint8_t id= 0; id < DV; id++) {
		hash1[id]= hash[id]= (256*M_DIV + s[L]  +  hash[id] * B_DIV[id]  -  C_DIV[id] * s[0]) % M_DIV;
	}
}

void hash_batch(
	uint8_t  s[], 			// input : current string buffer
	uint32_t hash_count, 	// input : number of hashes
	uint64_t com_hash[], 	// output: batch of (hash_count)      common hashes
	uint8_t  div_hash[]) 	// output: batch of (hash_count * DV) diversified hashes
{
	// produce batch of hashes (common & diversified)
	// for the shingles in the current input buffer s

	// compute diversified hashes
	// --------------------------
	// compute the hashes of the first, leftmost shingles
	for (uint8_t id= 0; id < DV; id++) {
		div_hash[id]= 0;
		for (uint32_t j= 0; j < L; j++) {
			div_hash[id]= (div_hash[id] * B_DIV[id] + s[j]) % M_DIV;
		}
		div_hash[DV+id]= div_hash[id];
	}
	// compute the hashes of the following shingles (rolling forward)
	// note: j= 1!
	for (uint32_t j= 1; j < hash_count; j++) {
		update_div_hashes(&div_hash[j*DV], &div_hash[(j+1)*DV], s+j-1);
	}

	// compute common hashes
	// ---------------------
	// compute the hashes of the first, leftmost shingle
	com_hash[0]= 0;
	for (uint32_t j= 0; j < L; j++) {
		com_hash[0]= (com_hash[0] * B_COM + s[j]) % M_COM;
	}
	// compute the hashes of the following shingles in the buffer
	for (uint32_t j= 0; j < hash_count; j++) {
		com_hash[j+1]= ((com_hash[j] + M_COM) * B_COM   -  C_COM * s[j]   +   s[j+L]) % M_COM;
	}
}

//	*********************************************************************************************************************************************

void rcp_generator(std::mt19937& mt_rand, uint8_t p[]) {
    std::uniform_int_distribution<uint8_t> dist(0, 255);
	// p: random cyclic permutations (see Sattolo / Fisher–Yates)
	int i, j;
	uint8_t rand;	// random number
	bool a[256];	// assigned random number
	for (j= 0; j < 256; j++) a[j]= false;
	for (j= 0; j < 256; j++) {
		rand= dist(mt_rand);
		while (a[rand]) {rand++;}
		a[rand]= true;
		p[j]= rand;
	}
	// test only
	for (j= 0; j < 256; j++ ) {
		for (i= j+1; i < 256; i++ ) {
			if (p[j] == p[i]) {
				printf("RCP ERROR\n");
				exit(25);
			}
		}
	}
}

//	*********************************************************************************************************************************************


