
Common Substrings Filter
------------------------

Designed and coded by Felix Baessler, felix.baessler@gmail.com

### Summary

This project is closely related to the file based version of the well-known 'common substring problem'. <br/> 
Given on auxiliary storage:
-	a test data set containing a big string S, and
-	a reference data set of a much smaller string s,<br/>

find all cross-repetitions defined by those strings of lengths >= LP,<br/>
which are substrings of both the reference and the test string.<br/>
A variant is LCS, the longest common substring problem.

Reference and test string are composed of IID distributed byte sequences. 
Their lengths underlie certain restrictions: 
-	ns [bytes]: the length of the reference string s, <br/>
is limited by the size of the available system memory (RAM)
-	NS [bytes]: the length of the test string S, is orders of magnitudes bigger; <br/>
its limitation is implied by the filtration ratio f, defined below.

**The objective is to reduce big ‘common substring problems' to a more manageable size.**

Based on the concepts: 
-	Shingle: the smallest information unit consists of a substring of length L that overlaps on L-1 bytes with the neighbors, like roof shingles.
-	Fingerprint: the hash(b,m) value of a shingle, base b and modulus m are the parameters of the Rabin fingerprint.
-	Diversification: to alleviate the ’locality’ drawback of fingerprinting, the hashes are diversified as presented in section 3 of the report.<br/>

the general idea is to eliminate those test shingles that have no matching fingerprint among the reference shingles. Accordingly we define:
-	residue r 	: the set of remaining test shingles with matching fingerprints
-	filtration ratio	: f= r/N,  N= NS - L + 1, the original number of test shingles.

### Implementation

The project consists of three C++ programs.

**A) master** <br/>
generates on disk a long enough IID distributed byte sequence of minimum ns+NS bytes. The master disk file comprises the
-	reference data set (ns ~ 1GB) concatenated to the
-	test data set (NS ~ 1TB)

Seamless shingling of the master file, requires that the reference data overlaps with the first L-1 bytes of the test data, so that
-	the last reference shingle starts at position: ns - 1
-	the first test shingle starts at position  : ns <br/>

Because of this overlap, the test string appears lengthened by L-1 bytes in the present implementation.

**B) scatter** <br/>
reads the reference data set (n= ns shingles), creates the fingerprint map and writes the result to the map file.
The map can be viewed as a minimalistic hash table reduced to m one-bit slots.
Scatter will mark those slots that correspond to the hash value modulo m (fingerprint) of the reference shingles.

**C) gather** <br/>
loads the map file into RAM, reads the big test data set (N= NS-L+1 shingles) and filters the test shingles by means of the map.
It turns out that the most time consuming operation consists in reading the map, when the fingerprints of a large amount of test shingles are checked against the fingerprints of the reference shingles marked by the map.
Run on an ordinary laptop, the throughput is of the order of 20 MB/s.

**Batchwise Processing** <br/>
Both scatter and gather distribute their workload on three threads:
-	thread 1: reads a batch of shingles into memory (RAM)
-	thread 2: from the shingle batch the thread produces a hash batch
-	thread 3: the hash batch is mapped:<br/>
  &nbsp; -	scatter (write): map slots are marked free -> occupied <br/>
  &nbsp; -	gather  (read) : map slots are checked free / occupied

With three containers, each containing a shingle- and the corresponding hash- batch, both scatter and gather, advance synchronously from stage to stage.

### Project Presentation
A more detailed report is available on https://sites.google.com/view/repsieve

### LICENSE
This project is released under [CC-BY-NC 4.0](https://creativecommons.org/licenses/by-nc/4.0/).<br/>
The licensing TLDR is: You are free to use, copy, distribute and transmit this Software for personal, non-commercial purposes, as long as you give attribution and share any modifications under the same license. Commercial or for-profit use requires a license. <br/>
For more details see the [LICENSE](https://github.com/ookraw/OOK-Raw-Data-Receiver/blob/master/LICENSE)
