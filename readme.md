
Common Substrings Filter
------------------------

Designed and coded by Felix Baessler, felix.baessler@gmail.com

### Summary

This project is closely related to the file based version of the well-known **common substring problem**. <br/> 
Given on auxiliary storage:
-	a test data set containing a very big string S, and
-	a reference data set of a much smaller string s,<br/>

find all cross-repetitions defined by those strings of lengths >= LP, <br/>
which are substrings of both the reference and the test string.<br/>
A variant is LCS, the longest common substring problem.

Reference and test string are composed of IID distributed byte sequences. 
Their lengths underlie certain restrictions: 
-	ns [bytes]: the length of the reference string s, <br/>
is limited by the size of the available system memory (usually several giga bytes)
-	NS [bytes]: the length of the test string S, is orders of magnitudes bigger; <br/>
its limitation is implied by the filtration ratio.

**The objective is to reduce big &nbsp;‘common substring problems'&nbsp; to a more manageable size.**

Based on the concepts: 
-	Shingle: the smallest information unit consists of a substring of length L that overlaps on L-1 bytes with the neighbors, like roof shingles.
-	Fingerprint: the hash(b,m) value of a shingle, base b and modulus m are the parameters of the Rabin-Karp Rolling Hash (aka signature algorithm).
-	Diversification: to alleviate the **locality problem** of fingerprinting, the hashes are diversified as presented in [On_Finding_Common_Substrings_between_two_Large_Files](https://www.researchgate.net/publication/370411448_On_Finding_Common_Substrings_between_two_Large_Files_by_Diversified_Hashing_and_Prefix_Shingling).<br/>

the general idea is to **eliminate those test shingles that have no matching fingerprint among the reference shingles.**<br/> 

### Implementation 

The project consists of three C++ programs.

**A) master** <br/>
generates on disk a long enough IID distributed byte sequence of minimum ns + NS bytes. The master disk file comprises the
-	reference data set (ns bytes) concatenated with the
-	test data set (NS bytes)

Seamless shingling of the master file, suggests that the reference data overlaps with the first L-1 bytes of the test data, so that
-	the last reference shingle starts at position : ns - 1
-	the first test shingle starts at position  : ns <br/>

Because of this overlap, the reference string appears lengthened by L-1 bytes in the present implementation.

**B) scatter** <br/>
reads the reference data (n= ns shingles), creates the fingerprint map and writes the result to the map file.<br/>
The map can be viewed as a minimalistic hash table reduced to m one-bit slots.
Scatter will mark those slots that correspond to the hash value modulo m (fingerprint) of the reference shingles.

**C) gather** <br/>
loads the map file into RAM, reads the big test data set (N= NS-L+1 shingles) and filters the test shingles by means of the map.<br/>
It turns out that the most time consuming operation consists in reading the map, when the fingerprints of a large amount of test shingles are checked via map against the fingerprints of the reference shingles.<br/>
Run on an ordinary laptop, the throughput is of the order of 20 MB/s.<br/>
An output example is given in the Appendix of the long write-up:  &nbsp;
[On_Finding_Common_Substrings_between_two_Large_Files](https://www.researchgate.net/publication/370411448_On_Finding_Common_Substrings_between_two_Large_Files_by_Diversified_Hashing_and_Prefix_Shingling).<br/>

**Batchwise Processing** <br/>
both scatter and gather distribute their workload on three threads:
-	thread 1: reads a batch of shingles into memory (RAM)
-	thread 2: from the shingle batch the thread produces a hash batch
-	thread 3: the hash batch is mapped:<br/>
  &nbsp; -	scatter (write to map) &nbsp;&nbsp;: slots are marked &nbsp;free -> occupied <br/>
  &nbsp; -	gather  (read from map): slots are checked free / occupied <br/>
  
In the present implementation logical processor 3 is reserved for thread 3, which guaranties that mapping takes place within the same thread. <br/>

Using three containers (a,b,c), each containing a shingle and the corresponding hash batch, scatter and gather advances stage by stage as illustrated by the example: <br/>

workload processed in 5 batches on 7 stages <br/>
thread 1 2 3 4 5 6 7       -> stage (time) <br/>
&nbsp;&nbsp; 1 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;         a b c a b <br/>
&nbsp;&nbsp; 2 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;            a b c a b <br/>
&nbsp;&nbsp; 3 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;               a b c a b <br/>
For example at stage 4 the threads (1,2,3) are simultaneously busy with the batches in containers (a,c,b). <br/>

### Description
For a more detailed write-up see: &nbsp;
[On_Finding_Common_Substrings_between_two_Large_Files](https://www.researchgate.net/publication/370411448_On_Finding_Common_Substrings_between_two_Large_Files_by_Diversified_Hashing_and_Prefix_Shingling).<br/>


### LICENSE
This project is released under [CC-BY-NC 4.0](https://creativecommons.org/licenses/by-nc/4.0/).<br/>
The licensing TLDR is: You are free to use, copy, distribute and transmit this Software for personal, non-commercial purposes, as long as you give attribution and share any modifications under the same license. Commercial or for-profit use requires a license. <br/>
For more details see the [LICENSE](https://github.com/ookraw/OOK-Raw-Data-Receiver/blob/master/LICENSE)
