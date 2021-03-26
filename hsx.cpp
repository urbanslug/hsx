/*
 HSX Hashed Sequence Index
 */
#include <iostream>
#include <string>
#include <string.h>
#include <fstream>
#include <vector>
#include <valarray>
#include <utility>
#include <math.h>
#include <algorithm>
#include <set>

#include <filesystem> // TODO: remove C++ 17 features

#define MAX_FILES 255
#define MAX_SEQUENCE 1000000
#define MAX_HEADER 255
#define HEADER_LENGTH 0x0000001C
#define MAX_STRING 255

// Structs and classes
// -------------------
struct Seq {
  uint32_t hash;
  char     name[MAX_HEADER];
  uint32_t length;
  char     file_name[MAX_HEADER];
  uint32_t file_num;
  uint32_t offset;
};

struct Short_Seq {
  uint32_t length;
  uint32_t file_num;
  uint32_t offset;
  char     name[MAX_HEADER];
};

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
  const uint32_t magic_number    = 0x957052D2;    // Little endian. This implementation only supports little endian.
  const uint32_t version         = 0x00000100;    // HSX Version number

  /*
    Header
   */
  const uint32_t header_length   = HEADER_LENGTH;
  uint32_t number_of_files       = 0;             // FLEN
  uint32_t file_table_offset     = 0x00000030;    // FOFF
  uint32_t number_of_buckets     = 0;             // HLEN
  uint32_t hash_table_offset     = 0;             // HOFF
  uint32_t number_of_sequences   = 0;             // SLEN
  uint32_t sequence_table_offset = 0;             // SOFF

  uint32_t header_padding[4]     = {0};           // Padding

  /*
    FINFO File table
    It holds offsets into the FTYPE Table
   */
  std::vector<uint32_t> file_table;               // FINFO_0 .. number of files
  uint32_t file_table_padding = 0;                // Padding

  /*
    FTYPE & FNAME file information records
    Contigous records of the the file extension (FTYPE) and base name (FNAME)
    from 0 .. number of files
   */
  std::vector<std::string > file_information_records;
  uint32_t info_padding = 0;                      // Padding

  /*
    SOFF Hash table/sequence offsets

    Buckets hold offsets from the file start to the sequence index array
    The hash table provides direct access to the sequence index array
   */
  std::vector<size_t> hash_table;
  uint32_t hash_table_padding = 0;

  /*
    IXLEN Length (in nucleotides) of the sequence
    IXFILE index into the file table (FINFO) for the file containing the sequence
    IXOFF  offset (from the start of the appropriate sequence file) pointing to the sequence.
    IXNAME name of the sequence ...

    Sequence index
   */
  std::vector<Short_Seq> sequence_index_array;         // SOFF table
  uint32_t soff_padding = 0;                      // Padding

  Hsx(uint32_t num_of_files, uint32_t num_of_seqs, uint32_t num_buckets, std::vector<std::string> file_info_records, std::vector<Seq> index) :
    number_of_files(num_of_files),
    number_of_buckets(num_buckets),
    number_of_sequences(num_of_seqs),
    file_information_records(file_info_records)
  {

    std::vector<Short_Seq> short_index;
    for (auto p: index) {
      Short_Seq foo;
      foo.length = p.length;
      foo.file_num = p.file_num;
      foo.offset = p.offset;
      strcpy(foo.name, p.name);

      short_index.push_back(foo);
    }

    sequence_index_array = short_index;

    // Compute offsets
    size_t file_start = 0;
    size_t constants_end = file_start + sizeof(magic_number) + sizeof(version);

    size_t header_end = constants_end + header_length + sizeof(header_padding);

    file_table_offset = header_end; // FOFF

    // TODO: check that this number is divisible by two
    uint32_t file_table_len = file_information_records.size()/2;
    file_table.reserve(file_table_len);

    size_t file_table_end = sizeof(file_table) + sizeof(file_table_padding);

    // populate the file table
    uint32_t bar = 2; // the length of a file info record
    int y = file_table_end;
    for (size_t i = 0; i+bar < file_information_records.size(); i += bar) {
      file_table.push_back(sizeof(std::string)*bar + y);
      y += sizeof(std::string)*bar;
    }

    size_t info_records_end = file_table_end + sizeof(file_information_records);

    hash_table_offset = info_records_end + info_padding; // HOFF

    std::vector<size_t> bucket_offsets; // offsets within the sequence index array
    int current_bucket, previous_bucket = 0;
    for (size_t i = 0; i < number_of_sequences; i++) {
      current_bucket = index[i].hash;
      if (current_bucket != previous_bucket) {
        bucket_offsets.push_back(i+1);
        previous_bucket = current_bucket;
      }
    }

    hash_table.reserve(number_of_buckets);
    size_t hash_table_end = hash_table_offset + sizeof(hash_table_padding);

    for (auto bucket_offset : bucket_offsets) {
      hash_table.push_back(sizeof(Seq)*bucket_offset + hash_table_end);
    }

    sequence_table_offset = hash_table_end; // SOFF


  }

  void write_file(std::string output_filename) {
    std::ofstream f(output_filename, std::ios::binary | std::ios::out);
    f.write((char*)this, sizeof(*this));
    f.close();
  }

  void read_file(std::string input_filename) {
    std::ifstream f(input_filename, std::ios::binary | std::ios::in);
    while (f.read((char*)this, sizeof(*this))) {
      this->info();
    }
    f.close();
  }


  void info() {
    // Header
    std::cout << "Number of files:\t"       << number_of_files
              << "\nNumber of sequences:\t" << number_of_sequences
              << std::endl;
    std::cout << "-----------------------------------------\n";
    // file information
    std::cout << "File information records" << std::endl;
    for (auto r : file_information_records) {
      std::cout << r << std::endl;
    }
    std::cout << "-----------------------------------------\n";

    // Hash table

    // Sequence index
    for (auto s: sequence_index_array) {
      std::cout << s.length   << "\t"
                << s.file_num << "\t"
                << s.offset   << "\t"
                << s.name     << std::endl;
    }
  }
};

class Fasta {
public:
  char header[MAX_HEADER];
  char sequence[MAX_SEQUENCE];
  int header_length;
  int sequence_length;
  int header_offset;
  int sequence_offset;
  int file_number;
  uint32_t short_hash = 98;
  uint32_t long_hash;


  Fasta() {}
  Fasta (std::string h, int h_len, std::string s, int s_len, int h_offset, int s_offset, int f_num) :
    header_length(h_len),
    sequence_length(s_len),
    header_offset(h_offset),
    sequence_offset(s_offset),
    file_number(f_num)
  {
    long_hash = hassock_hash(h.c_str(), h_len);
    strcpy(header, h.c_str());
    strcpy(sequence, s.c_str());
  }


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

  void set_hash(int num_buckets) {
    short_hash = long_hash % num_buckets;
  }

  uint32_t get_hash() {
    return short_hash;
  }

  void print() {
    std::cout << '>' << header << std::endl;
    for (int i = 0; i < sequence_length; i += 80 ) {
      for (int j = i; (j < i+80 && j < sequence_length) ; j++) {
        std::cout << sequence[j];
      }
      std::cout << std::endl;
    }

    std::cout << "Header length "     << header_length
              << "\thash "            << short_hash
              << "\tsequence length " << sequence_length
              << "\theader offset "   << header_offset
              << "\tsequence offset " << sequence_offset
              << "\tfile number "     << file_number
              << "\tTotal "           << header_length + sequence_length
              << std::endl;
  }
};

/*
 Parse a fasta file...
 */
std::vector<Fasta> read_fasta_file(std::string fasta_file_path, size_t file_num) {
  std::vector<Fasta> sequences;

  int line_offset = 0;
  int header_offset = 0;
  int sequence_offset = 0;
  int counter = 0;

  std::string header;
  std::string seq;

  std::string line;

  std::ifstream myfile;
  myfile.open(fasta_file_path);

  if (myfile.is_open()) {
    while (getline(myfile, line)) {
      line_offset = counter;
      counter += line.size() + 1;

      // Process a header
      if (line[0] == '>') {
        if (header != "") {
          // if there was a header here previously
          // save the current header and seq
          Fasta fasta(header, header.size(), seq, seq.size(), header_offset, sequence_offset, file_num);

          // clear the header and seq
          sequences.push_back(fasta);
          seq.clear();
          header.clear();
        }

        sequence_offset = 0; // None in the python
        header_offset = line_offset;
        // set header to the current line
        header = line.substr(1, line.size());
      } else {
        if (sequence_offset == 0) sequence_offset = line_offset;
        seq += line;
      }
    }

    // add the last one
    if (sequence_offset == 0) sequence_offset = line_offset;
    Fasta fasta(header, header.size(), seq, seq.size(), header_offset, sequence_offset, file_num);
    sequences.push_back(fasta);

    myfile.close();
  } else {
    std::cerr << "Unable to open " << fasta_file_path << std::endl;
    exit(EXIT_FAILURE);
  }

  return sequences;
}

/*
  Sort Seqs in ascending order of the buckets they are in and
  also sort them lexicographically within those buckets.

  Sort based on the hash.
  If hashes are equal, sort based on the lexicographic order of thier names.
 */
bool sorter (Seq i, Seq j) {
  if (i.hash == j.hash) {
    return std::string(i.name) < std::string(j.name);
  } else {
    return i.hash < j.hash;
  }
}


int main() {
  std::string filename;
  std::vector<std::string> filenames;

  while (std::cin >> filename) {
    filenames.push_back(filename);
  }

  // Check for max files
  int file_count = filenames.size();
  if (file_count > MAX_FILES) {
    std::cerr << "Maximum number of files (255) exceeded" << std::endl;
    exit(EXIT_FAILURE);
  }

  // Generate the file information records
  std::set<std::string> seen_filenames;
  std::vector<std::string> file_information_records;
  for (auto filename : filenames) {

    // Ensure there are no duplicate filenames
    if (seen_filenames.find(filename) != seen_filenames.end()) {
      std::cerr << "Hsx does not support two files with the same filename but duplicate filename found "
                << filename << std::endl;
      exit(EXIT_FAILURE);
    }
    seen_filenames.insert(filename);

    // TODO: make backwards compatible stop using std::filesystem
    std::filesystem::path p (filename);
    std::string dot_ext = p.extension();
    std::string ext = dot_ext.substr(1, dot_ext.size());
    std::string base_name = p.stem();

    // TODO: check that extension is fasta or fa
    if (false) {
      std::cerr << "Hsx only support fasta or fa extensions" << std::endl;
      exit(EXIT_FAILURE);
    }

    file_information_records.push_back(ext);
    file_information_records.push_back(base_name);
  }

  std::vector<Fasta> all_sequences;

  // collate sequences from all the files into a single vector
  for (size_t i = 0; i < filenames.size(); i++) {
    std::string f = filenames[i];
    std::vector<Fasta> f_sequences = read_fasta_file(f, i);

    // preallocate memory, faster, to avoid constant re-allocation.
    std::vector<Fasta> temp;
    temp.reserve( all_sequences.size() + f_sequences.size() );

    // merge all_sequences and sequences from this iteration, f_sequences
    temp.insert( temp.end(), all_sequences.begin(), all_sequences.end() );
    temp.insert( temp.end(), f_sequences.begin(), f_sequences.end() );
    all_sequences = temp;
  }

  int num_sequences = all_sequences.size();
  int avg_bucket = 10; // TODO: move this upwards?
  int num_buckets = floor((num_sequences+avg_bucket-1)/avg_bucket);

  // Generate the sequence index from the Fasta objects (all_sequences)
  // The sequence index is an array of Seqs
  std::vector<Seq> index;
  index.reserve(num_sequences);
  for (auto i =  all_sequences.begin(); i !=  all_sequences.end(); i++) {
    i->set_hash(num_buckets);

    Seq s;
    s.hash = i->get_hash();
    strcpy(s.name, &i->header[0]);
    s.length = i->sequence_length;
    s.file_num = i->file_number;
    s.offset = i->header_offset;
    index.push_back(s);
  }

  // sort index
  std::sort(index.begin(), index.end(), sorter);

  Hsx obj(file_count, num_sequences, num_buckets, file_information_records, index);

  obj.write_file("/home/sluggie/src/bio/hsx/out.hsx");
  obj.read_file("/home/sluggie/src/bio/hsx/out.hsx");

  return 0;
}
