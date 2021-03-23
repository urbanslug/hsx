#include <iostream>
#include <string>
#include <string.h>
#include <fstream>
#include <vector>
#include <valarray>
#include <math.h>

#define MAX_SEQUENCE 1000000
#define MAX_HEADER 255
#define HEADER_LENGTH 0x0000001C


class Hsx {

/*
"Hashed sequence index" (hsx) file reader (for a fasta file)
-------------------------------------------------------------------

offset 0x00: D2 52 70 95        big endian magic number
																.. (95 70 52 D2 => little endian)
offset 0x04: 00 00 01 xx        version 1.0 (see note 1)
offset 0x08: 00 00 00 1C        header length (in bytes, including this
									 							.. field)
offset 0x0C: xx xx xx xx        FN, number of files (see note 2)
offset 0x10: xx xx xx xx        FO, offset to file table
offset 0x14: xx xx xx xx        HN, number of hash buckets (see notes 3 and 4)
offset 0x18: xx xx xx xx        HO, offset to hash table
offset 0x1C: xx xx xx xx        SN, number of sequences
offset 0x20: xx xx xx xx        SO, offset to sequence index table (see
								 								.. note 5)

offset FO:   xx xx xx xx        FIO0, offset to file info for file 0
			 			  ...               (FN-1 more entries, at 4 bytes per)

offset FIOn: LL xx ..           type of file (ascii "fa", "2bit", etc., see
                                note 6)
						 LL xx ..           name of file (see note 7)
			 			  ...               (FN-1 more entries, variable length)

offset HO:   xx xx xx xx xx     SIOn, offset into sequence index table (see
								 								.. notes 8, 9 and 10)
		 	 			  ...               (HN-1 more entries, at 5 bytes per)
				 	 	 xx xx xx xx xx     offset past end of sequence index table

offset SO:   xx xx xx xx xx     length of the sequence (see note 11)
			 xx                 file number (index into file table)
			 xx xx xx xx xx xx  offset to the sequence data (see note 12)
			 LL xx ..           name of sequence (see note 13)
			  ...               (SN-1 more entries, variable length)

Notes:

	(1)  The least significant byte of the version is the "sub version".
	     For version 1, this is 00 (secondary hashes are not in use) or 01
	     (secondary hashes are in use).
	(2)  The number of files is limited to 255.
	(3)  It is assumed that the number of buckets is set so that the average
	     number of sequences per bucket (SN/HN) is reasonably small (e.g. 10).
	(4)  The hash table actually includes HN+1 buckets.  The extra bucket has
	     size zero and gives the offset to just past the end of the sequence
	     index table.
	(5)  Entries in the sequence index table are necessarily stored in hash
	     order.  Entries with the same hash are stored in alphabetical order;
	     actually, in lexicographic order over the bytes of their names.
	(6)  Strings are stored as a length byte followed by ascii text.
	(7)  If a file info record contains an empty name, the name of the file is
	     the same as the index file itself, with the file type used as the
	     extension (e.g. "reads.hsx" becomes "reads.fa").  This allows files to
	     be renamed without rebuilding the index.
	(8)  SIOn is the file offset for the nth entry in the sequence index table.
	     When this is in a hash table entry, it is the index for the first
	     sequence in that hash's bucket.
	(9)  The most significant bit in a bucket's SIOn value is used to indicate
	     whether the bucket is empty or not.  If a bucket is empty, this bit is
	     set (1), otherwise it is clear.
	(10) The end of a bucket can be determined from the SIOn entry for the
	     start of the next bucket.
	(11) A sequence may be empty, so zero is a legitimate value for the
	     sequence length.
	(12) The offset to the sequence data is an offset into the sequence file.
	     For fasta it can point to the ">" at the start of the sequence's
	     header, or directly to the sequence data.
	(13) When secondary hashes are in use, the sequence name (including the
	     terminating zero) is replaced by the four-byte secondary hash.

 */

public:
  const uint32_t magic_number  = 0x957052D2; // little endian
  const uint32_t version       = 0x00000100;
  const uint32_t header_length = HEADER_LENGTH;
  uint32_t number_of_files = 0;                 // FLEN
  uint32_t file_table_offset = 0;               // FOFF
  uint32_t number_of_buckets_in_hash_table = 0; // HLEN
  uint32_t hash_table_offset = 0;               // HOFF
  uint32_t number_of_sequences = 0;
  uint32_t sequence_table_offset = 0;



  Hsx(uint32_t num_of_files, uint32_t num_of_seqs, uint32_t num_buckets) :
    number_of_sequences(num_of_seqs),
    number_of_files(num_of_files),
    number_of_buckets_in_hash_table(num_buckets)
  {}

  void write_file(std::string output_filename) {
    std::ofstream f(output_filename, std::ios::binary | std::ios::app);
    f.write((char*)this, sizeof(*this));
    f.close();
  }

  void info() {
    std::cout << "Number of files:\t"       << number_of_files
              << "\nNumber of sequences:\t" << number_of_sequences
              << std::endl;
  }

  void read_file(std::string input_filename) {
    std::ifstream f(input_filename, std::ios::binary | std::ios::in);
    while (f.read((char*)this, sizeof(this))) {
      this->info();
    }
    f.close();
  }

};

class Fasta {
public:
  char header[MAX_HEADER];
  char sequence[MAX_SEQUENCE];
  int header_length;
  int sequence_length;

  Fasta() {}
  Fasta (std::string h, int h_len, std::string s, int s_len) :
    header_length(h_len), sequence_length(s_len) {
    strcpy(header, h.c_str());
    strcpy(sequence, s.c_str());
  }

  void print() {
    std::cout << '>' << header << std::endl;
    for (int i = 0; i < sequence_length; i += 80 ) {
      for (int j = i; (j < i+80 && j < sequence_length) ; j++) {
        std::cout << sequence[j];
      }
      std::cout << std::endl;
    }
  }
};

uint32_t hassock_hash (const void* key, uint32_t len) {
  const uint32_t seed = 0x5C3FC4D3;
  const uint32_t m    = 0x87C10417;
  const uint8_t* data = ((const uint8_t*) key) + len;
  const uint8_t* stop = ((const uint8_t*) key) + 4;
  uint32_t       h, k;

  h = seed ^ len;
  while (data >= stop)
    {
      k  = *(--data);
      k |= *(--data) << 8;
      k |= *(--data) << 16;
      k |= *(--data) << 24;
      k *= m;
      k ^= k >> 24;
      k *= m;
      h *= m;
      h ^= k;
      len -= 4;
    }
  switch (len)
    {
    case 3: h ^= *(--data) << 16;
    case 2: h ^= *(--data) << 8;
    case 1: h ^= *(--data);
      h *= m;
    }
  h ^= h >> 13;
  h *= m;
  h ^= h >> 15;
  return h;
}

std::vector<Fasta> read_fasta_file(std::string fasta_file_path) {
  std::vector<Fasta> sequences;

  std::string header;
  std::string seq;

  std::string line;

  std::ifstream myfile;
  myfile.open(fasta_file_path);

  if (myfile.is_open()) {
    while (getline(myfile, line)) {
      if (line[0] == '>') {
        if (header != "") {
          Fasta fasta(header, header.size(), seq, seq.size());
          sequences.push_back(fasta);
          seq.clear();
          header.clear();
        }
        header = line.substr(1, line.size());
      } else {
        seq += line;
      }
    }

    // add the last one
    Fasta fasta(header, header.size(), seq, seq.size());
    sequences.push_back(fasta);

    myfile.close();
  } else {
    std::cerr << "Unable to open " << fasta_file_path << std::endl;
    exit(EXIT_FAILURE);
  }

  return sequences;
}

int main() {
  std::string filename;
  std::vector<std::string> filenames;

  while (std::cin >> filename) {
    filenames.push_back(filename);
  }

  int file_count = filenames.size();
  if (file_count > 255) {
    std::cerr << "Maximum number of files (255) exceeded" << std::endl;
    exit(EXIT_FAILURE);
  }

  std::vector<Fasta> all_sequences;

  // collate sequences from all the files into a single vector
  for (std::string f : filenames) {
    std::vector<Fasta> f_sequences = read_fasta_file(f);

    // preallocate memory, faster, to avoid constant re-allocation.
    std::vector<Fasta> temp;
    temp.reserve( all_sequences.size() + f_sequences.size() );

    // merge all_sequences and sequences from this iteration, f_sequences
    temp.insert( temp.end(), all_sequences.begin(), all_sequences.end() );
    temp.insert( temp.end(), f_sequences.begin(), f_sequences.end() );
    all_sequences = temp;
  }

  int num_sequences = all_sequences.size();

  // std::vector<Fasta> index(HEADER_LENGTH, Fasta());
  // index is a 2D array
  int avg_bucket = 10;
  int num_buckets = floor((num_sequences+avg_bucket-1)/avg_bucket);
  std::vector<std::vector<Fasta> > index(HEADER_LENGTH, std::vector<Fasta>());
  for (Fasta i : all_sequences ) {
    uint32_t bucket = hassock_hash(i.header,i.header_length) % num_buckets;
    index[bucket].push_back(i);
  }

  for(auto i = 0; i < HEADER_LENGTH; i++) {
    std::vector<Fasta> bucket = index[i];
    if (!bucket.empty()) {
      std::cout << i << std::endl;
      for (auto e : bucket) {
        std::cout << e.header << " " << std::endl;
      }
    }
  }


  Hsx obj(file_count, num_sequences, num_buckets);
  obj.write_file("/home/sluggie/src/learning/C++/misc/out.hsx");
  obj.read_file("/home/sluggie/src/learning/C++/misc/out.hsx");


  return 0;
}
