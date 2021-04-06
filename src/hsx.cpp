/*
 HSX Hashed Sequence Index
 Spec: http://www.bx.psu.edu/miller_lab/dist/README.lastz-1.02.00/hsx_format.html
 Lastz documentation http://www.bx.psu.edu/~rsharris/lastz/README.lastz-1.04.03.html#fmt_hsx
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

#define MAGIC_BIG 0xD2527095
#define MAGIC_LITTLE 0x957052D2
#define HSX_VERSION 0x00000100
#define HEADER_LENGTH 0x0000001C
#define MAX_FILES 255

bool DEBUG = false;


// Endianness switching templates from https://github.com/alipha/cpp/blob/master/endian/endian.hpp
template<typename Char>
Char *uint64_5_to_le(std::uint64_t src, Char *dest) {
  static_assert(sizeof(Char) == 1, "Char must be a byte-sized type");
  //dest[7] = static_cast<Char>(static_cast<std::uint8_t>(src >> 56));
  //dest[6] = static_cast<Char>(static_cast<std::uint8_t>(src >> 48));
  //dest[5] = static_cast<Char>(static_cast<std::uint8_t>(src >> 40));
  dest[4] = static_cast<Char>(static_cast<std::uint8_t>(src >> 32));
  dest[3] = static_cast<Char>(static_cast<std::uint8_t>(src >> 24));
  dest[2] = static_cast<Char>(static_cast<std::uint8_t>(src >> 16));
  dest[1] = static_cast<Char>(static_cast<std::uint8_t>(src >> 8));
  dest[0] = static_cast<Char>(static_cast<std::uint8_t>(src));
  return dest;
}

template<typename Char>
Char *uint64_6_to_le(std::uint64_t src, Char *dest) {
  static_assert(sizeof(Char) == 1, "Char must be a byte-sized type");
  // dest[7] = static_cast<Char>(static_cast<std::uint8_t>(src >> 56));
  // dest[6] = static_cast<Char>(static_cast<std::uint8_t>(src >> 48));
  dest[5] = static_cast<Char>(static_cast<std::uint8_t>(src >> 40));
  dest[4] = static_cast<Char>(static_cast<std::uint8_t>(src >> 32));
  dest[3] = static_cast<Char>(static_cast<std::uint8_t>(src >> 24));
  dest[2] = static_cast<Char>(static_cast<std::uint8_t>(src >> 16));
  dest[1] = static_cast<Char>(static_cast<std::uint8_t>(src >> 8));
  dest[0] = static_cast<Char>(static_cast<std::uint8_t>(src));
  return dest;
}

template<typename Char>
Char *uint64_6_to_be(std::uint64_t src, Char *dest) {
  static_assert(sizeof(Char) == 1, "Char must be a byte-sized type");
  dest[0] = static_cast<Char>(static_cast<std::uint8_t>(src >> 56));
  dest[1] = static_cast<Char>(static_cast<std::uint8_t>(src >> 48));
  dest[2] = static_cast<Char>(static_cast<std::uint8_t>(src >> 40));
  dest[3] = static_cast<Char>(static_cast<std::uint8_t>(src >> 32));
  dest[4] = static_cast<Char>(static_cast<std::uint8_t>(src >> 24));
  dest[5] = static_cast<Char>(static_cast<std::uint8_t>(src >> 16));
  //dest[6] = static_cast<Char>(static_cast<std::uint8_t>(src >> 8));
  //dest[7] = static_cast<Char>(static_cast<std::uint8_t>(src));
  return dest;
}

template<typename Char>
Char *uint64_5_to_be(std::uint64_t src, Char *dest) {
  static_assert(sizeof(Char) == 1, "Char must be a byte-sized type");
  dest[0] = static_cast<Char>(static_cast<std::uint8_t>(src >> 56));
  dest[1] = static_cast<Char>(static_cast<std::uint8_t>(src >> 48));
  dest[2] = static_cast<Char>(static_cast<std::uint8_t>(src >> 40));
  dest[3] = static_cast<Char>(static_cast<std::uint8_t>(src >> 32));
  dest[4] = static_cast<Char>(static_cast<std::uint8_t>(src >> 24));
  //dest[5] = static_cast<Char>(static_cast<std::uint8_t>(src >> 16));
  //dest[6] = static_cast<Char>(static_cast<std::uint8_t>(src >> 8));
  //dest[7] = static_cast<Char>(static_cast<std::uint8_t>(src));
  return dest;
}

template<typename Char>
Char *uint32_to_le(std::uint32_t src, Char *dest) {
  static_assert(sizeof(Char) == 1, "Char must be a byte-sized type");
  dest[3] = static_cast<Char>(static_cast<std::uint8_t>(src >> 24));
  dest[2] = static_cast<Char>(static_cast<std::uint8_t>(src >> 16));
  dest[1] = static_cast<Char>(static_cast<std::uint8_t>(src >> 8));
  dest[0] = static_cast<Char>(static_cast<std::uint8_t>(src));
  return dest;
}

template<typename Char>
Char *uint32_to_be(std::uint32_t src, Char *dest) {
  static_assert(sizeof(Char) == 1, "Char must be a byte-sized type");
  dest[0] = static_cast<Char>(static_cast<std::uint8_t>(src >> 24));
  dest[1] = static_cast<Char>(static_cast<std::uint8_t>(src >> 16));
  dest[2] = static_cast<Char>(static_cast<std::uint8_t>(src >> 8));
  dest[3] = static_cast<Char>(static_cast<std::uint8_t>(src));
  return dest;
}

template<typename Char>
std::uint32_t be_to_uint32(const Char *src) {
  static_assert(sizeof(Char) == 1, "Char must be a byte-sized type");
  return static_cast<std::uint32_t>(
    static_cast<std::uint32_t>(static_cast<std::uint8_t>(src[0])) << 24
    | static_cast<std::uint32_t>(static_cast<std::uint8_t>(src[1])) << 16
    | static_cast<std::uint32_t>(static_cast<std::uint8_t>(src[2])) << 8
    | static_cast<std::uint32_t>(static_cast<std::uint8_t>(src[3]))
    );
}

template<typename Char>
std::uint32_t le_to_uint32(const Char *src) {
  static_assert(sizeof(Char) == 1, "Char must be a byte-sized type");
  return static_cast<std::uint32_t>(
    static_cast<std::uint32_t>(static_cast<std::uint8_t>(src[3])) << 24
    | static_cast<std::uint32_t>(static_cast<std::uint8_t>(src[2])) << 16
    | static_cast<std::uint32_t>(static_cast<std::uint8_t>(src[1])) << 8
    | static_cast<std::uint32_t>(static_cast<std::uint8_t>(src[0]))
  );
}



/*
  Structs and classes
  -------------------
*/
struct Seq {
  uint32_t    hash;
  std::string name;
  uint32_t    length;
  uint32_t    file_num;
  uint32_t    offset;
};

struct Short_Seq {
  uint64_t    length;
  uint8_t     file_num;
  uint64_t    offset;
  std::string name;
};

struct Path {
  std::string directory;
  std::string stub;
  std::string extension;
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
  /*
    Header
  */
  uint32_t magic_number    = MAGIC_LITTLE;     // Big endian
  uint32_t version         = HSX_VERSION;   // HSX version number

  uint32_t header_length   = HEADER_LENGTH;
  uint32_t number_of_files       = 0;             // FLEN
  uint32_t file_table_offset     = 0;             // FOFF
  uint32_t number_of_buckets     = 0;             // HLEN
  uint32_t hash_table_offset     = 0;             // HOFF
  uint32_t number_of_sequences   = 0;             // SLEN
  uint32_t sequence_table_offset = 0;             // SOFF

  std::vector<uint8_t> header_padding;            // Padding

  /*
    FINFO File table
    It holds offsets into the FTYPE Table
   */
  std::vector<uint32_t> file_table;               // FINFO_0 .. number of files
  std::vector<uint8_t> file_table_padding;        // Padding

  /*
    FTYPE & FNAME file information records
    Contigous records of the the file extension (FTYPE) and base name (FNAME)
    from 0 .. number of files
   */
  std::vector<std::string> file_information_records;
  std::vector<uint8_t> file_information_records_padding;     // Padding

  /*
    SOFF Hash table/sequence offsets

    Buckets hold offsets from the file start to the sequence index array
    The hash table provides direct access to the sequence index array
   */
  std::vector<uint64_t> hash_table;                   // This is HOFF
  std::vector<uint8_t> hash_table_padding;     // Padding

  /*
    IXLEN  Length (in nucleotides) of the sequence
    IXFILE index into the file table (FINFO) for the file containing the sequence
    IXOFF  offset (from the start of the appropriate sequence file) pointing to the sequence.
    IXNAME name of the sequence ...

    Sequence index
   */
  std::vector<Short_Seq> sequence_index_array;               // SOFF table
  uint32_t soff_padding = 0;                                 // Padding

  uint32_t pad_for_16(uint32_t n) {
    return ((16 - (n % 16)) % 16);
  }

  Hsx() {}

  Hsx(uint32_t num_of_files, uint32_t num_of_seqs, uint32_t num_buckets,
      std::vector<std::string> file_info_records, std::vector<Seq> index) :
    number_of_files(num_of_files),
    number_of_buckets(num_buckets),
    number_of_sequences(num_of_seqs),
    file_information_records(file_info_records)
  {

    // ----------------------------------------------
    //     Header
    // ----------------------------------------------

    uint32_t file_start = 0;
    uint32_t constants_end = file_start + sizeof(magic_number) + sizeof(version);

    uint32_t th = constants_end + sizeof(header_length) + sizeof(number_of_files)
      + sizeof(file_table_offset) + sizeof(number_of_buckets)
      + sizeof (hash_table_offset) + sizeof(number_of_sequences)
      + sizeof(sequence_table_offset);

    uint32_t header_padding_bytes = pad_for_16(8+header_length);
    for (uint32_t i = 0; i < header_padding_bytes; i++){
      header_padding.push_back(0);
    }
    uint32_t header_padding_size = header_padding.size() * sizeof(char);

    uint32_t header_end = th + header_padding_size;

    if (DEBUG)
      printf("header end \t %#010x\n", header_end);

    file_table_offset = header_end; // set FOFF

    // ----------------------------------------------
    //     FINFO
    // ----------------------------------------------

    file_table.reserve(number_of_files);
    uint32_t file_table_size = sizeof(uint32_t)*number_of_files;

    uint32_t file_table_length = number_of_files * 4;
    uint32_t file_table_padding_bytes = pad_for_16(file_table_length);
    for (uint32_t i = 0; i < file_table_padding_bytes; i++){
      file_table_padding.push_back(0);
    }
    uint32_t file_table_padding_size = file_table_padding.size() * sizeof(char);

    uint32_t file_table_end = file_table_size + file_table_padding_size + header_end;

    // populate the file table (FINFO)
    uint32_t info_record_offsets = file_table_end;
    uint32_t step = 2;
    uint32_t inc = 1;
    for (auto i = 0; i+inc < file_information_records.size(); i += step) {
      file_table.push_back(info_record_offsets);
      info_record_offsets += file_information_records[i].size() + file_information_records[i+inc].size();
    }

    if (DEBUG)
      printf("file table end \t %#010x\n", file_table_end);

    // ----------------------------------------------
    //     FYPE and FNAME
    // ----------------------------------------------

    // the length of the data held by information records
    // each char is a byte
    int total_string_length = 0;
    int fileInfoLength = 0;
    for (auto &s: file_information_records) {
      total_string_length += s.size();
      fileInfoLength += s.size()+1;
    }

    uint32_t file_information_records_padding_bytes = pad_for_16(fileInfoLength);
    for (uint32_t i = 0; i < file_information_records_padding_bytes; i++){
      file_information_records_padding.push_back(0);
    }
    uint32_t file_information_records_padding_size = file_information_records_padding.size() * sizeof(char);

    if (DEBUG) {
        printf("file info len \t %#010x file info pad  %#010x\n", fileInfoLength, file_information_records_padding_bytes);
      }

    uint32_t info_records_end = file_table_end + fileInfoLength + file_information_records_padding_size;


    // ----------------------------------------------
    //     SOFF
    // ----------------------------------------------
    hash_table_offset = info_records_end; // set HOFF
    if (DEBUG) {
      printf("hash table offset \t %#010x\n", hash_table_offset);
    }

    hash_table.reserve(number_of_buckets+1);

    int hashTableLength = (number_of_buckets+1) * 5;
    uint32_t  hashTablePad = pad_for_16(hashTableLength);

    for (uint32_t i = 0; i < hashTablePad; i++){
      hash_table_padding.push_back(0);
    }

    uint32_t hash_table_padding_size = hash_table_padding.size() * sizeof(char);

    // the number of elements in the hash table * size of each element
    //uint32_t hash_table_length = (number_of_buckets+1)*sizeof(uint32_t);
    // the hash table entries are 5 wide
    uint32_t hash_table_length = (number_of_buckets+1)*5;

    // TODO: check that number_of_buckets == hash_table.size()

    uint32_t hash_table_end = hash_table_offset
      + hash_table_length
      + hash_table_padding_size;

    if (DEBUG) {
      printf("hash table len \t %#010x hash table pad  %#010x\n", hash_table_length, hash_table_padding_size);
    }

    // ----------------------------------------------
    //     IXLEN, IXFILE, IXOFF, & IXNAME
    // ----------------------------------------------
    sequence_table_offset = hash_table_end; // set SOFF

    if (DEBUG) {
      printf("sequence Table offset \t %#010x\n", sequence_table_offset);
    }

    std::vector<std::pair<uint32_t, uint32_t>> member_sizes; // bucket and size
    uint32_t acc = 0;
    std::vector<Short_Seq> short_index;
    for (auto p: index) {
      Short_Seq foo;
      foo.length   = p.length;
      foo.file_num = p.file_num;
      foo.offset   = p.offset;
      foo.name     = std::string(p.name);

      short_index.push_back(foo);

     uint32_t sz = sizeof(foo.length) + sizeof(foo.file_num) + sizeof(foo.offset) + foo.name.size();
     acc += sz;
     member_sizes.push_back(std::make_pair(p.hash, sz));
    }

    sequence_index_array = short_index;

    // ---------------------------

    // Populate the hash table that provides direct access to the sequence index

    std::vector<uint32_t> l (number_of_buckets, 0);
    for (auto &m : member_sizes) {
      l[m.first] += m.second;
    }

    uint32_t start = sequence_table_offset;
    for (auto &sz : l) {
      hash_table.push_back(start);
      start += sz;
      if (sz == 0) {
        // set the most significant bit to 1
      }
    }

    // set the sentinel bucket
    hash_table.push_back(sequence_table_offset + acc);

    /*
    // TODO: remove this debug print statement
    std::cout << "number of buckets " << number_of_buckets
              << " hash table size " << hash_table.size()
              << " sentinel " << sequence_table_offset + acc
              << std::endl;
    */
  }

  void write_file(std::string output_filename) {
    char buffer[8];
    std::ofstream f(output_filename, std::ios::binary | std::ios::out);


    // TODO: Why do I have to write this magic number in big endian?
    f.write(uint32_to_be(magic_number, buffer), sizeof(magic_number));
    f.write(uint32_to_le(version, buffer), sizeof(version));

    // Write header
    f.write(uint32_to_le(header_length, buffer),     sizeof(header_length));
    f.write(uint32_to_le(number_of_files, buffer), sizeof(number_of_files));
    f.write(uint32_to_le(file_table_offset, buffer), sizeof(file_table_offset));
    f.write(uint32_to_le(number_of_buckets, buffer), sizeof(number_of_buckets));
    f.write(uint32_to_le(hash_table_offset, buffer), sizeof(hash_table_offset));
    f.write(uint32_to_le(number_of_sequences, buffer), sizeof(number_of_sequences));
    f.write(uint32_to_le(sequence_table_offset, buffer), sizeof(sequence_table_offset));

    for (auto &p : header_padding)
      f.write((char*)&p, sizeof(p));

    // -------
    for (auto &i : file_table) {
      f.write(uint32_to_le(i, buffer), sizeof(i));
    }

    for (auto &p : file_table_padding)
      f.write((char*)&p, sizeof(p));

    // -------
    for (auto &k : file_information_records) {
      // TODO: Why this +1 ?
      char len = static_cast<char>(k.size());
      f.write(&len, sizeof(len));
      f.write(k.data(), k.size());
    }

    for (auto &p : file_information_records_padding)
      f.write((char*)&p, sizeof(p));

    // --- hash table ---
    for (auto &i : hash_table) {
      f.write(uint64_5_to_le(i, buffer), 5);
    }

    for (auto &p : hash_table_padding)
      f.write((char*)&p, sizeof(p));

    // --- sequence index -----
    for (auto &i : sequence_index_array) {

      f.write(uint64_5_to_le(i.length, buffer), 5);
      f.write((char*)&i.file_num, sizeof(i.file_num));
      f.write(uint64_6_to_le(i.offset, buffer), 6);

      char len = static_cast<char>(i.name.size());
      f.write(&len, sizeof(len));
      f.write(i.name.data(), i.name.size());
    }

    f.close();
  }

  void read_file(std::string input_filename) {
    std::ifstream f(input_filename, std::ios::binary | std::ios::in);
    char buffer[8];

    f.read(buffer, sizeof(this->magic_number));
    this->magic_number = be_to_uint32(buffer);

    f.read(buffer, sizeof(version));
    this->version = le_to_uint32(buffer);

    // Header
    f.read((char*)&this->header_length, sizeof(this->header_length));
    f.read((char*)&this->number_of_files, sizeof(this->number_of_files));
    f.read((char*)&this->file_table_offset, sizeof(this->file_table_offset));
    f.read((char*)&this->number_of_buckets, sizeof(this->number_of_buckets));
    f.read((char*)&this->hash_table_offset, sizeof(this->hash_table_offset));
    f.read((char*)&this->number_of_sequences, sizeof(this->number_of_sequences));
    f.read((char*)&this->sequence_table_offset, sizeof(this->sequence_table_offset));

    // f.read((char*)&this->header_padding, sizeof(this->header_padding));

    f.close();
  }

  void info() {
    printf("magic number \t %#010x\n", magic_number);
    printf("version \t %#010x\n", version);

    // Header
    printf("header length \t %#010x\n", header_length);
    printf("number of files \t %#010x\n", number_of_files);
    printf("file table offset \t %#010x\n", file_table_offset);
    printf("number of buckets \t %#010x\n", number_of_buckets);
    printf("hash table offset \t %#010x\n", hash_table_offset);
    printf("number of sequences \t %#010x\n", number_of_sequences);
    printf("sequence table offset \t %#010x\n", sequence_table_offset);

    for (auto &i : header_padding) {
      printf("header padding \t %#010x\n", i);
    }
  }
};

class Fasta {
public:
  std::string header;
  std::string sequence;
  int header_length;
  int sequence_length;
  int header_offset;
  int sequence_offset;
  int file_number;
  uint32_t short_hash = 98;
  uint32_t long_hash;


  Fasta() {}
  Fasta (std::string h, int h_len, std::string s, int s_len, int h_offset, int s_offset, int f_num) :
    header(h),
    sequence(s),
    header_length(h_len),
    sequence_length(s_len),
    header_offset(h_offset),
    sequence_offset(s_offset),
    file_number(f_num)
  {
    long_hash = hassock_hash(h.c_str(), h_len);
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
  Helper functions
  ----------------
*/

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

std::string lowercase_string(std::string str) {
  auto lowercase_char = [](unsigned char c){ return std::tolower(c); };
  std::transform(str.begin(), str.end(), str.begin(), lowercase_char);
  return str;
}

Path parse_path (std::string filename) {
  std::string::size_type idx_ext, idx_base;
  Path p;

  idx_ext = filename.rfind('.');

  if(idx_ext != std::string::npos) {
      p.extension = filename.substr(idx_ext+1);
  } else {
    std::cerr << "Could not find the extension on file. HSX expects the extension as .fa or .fasta"
              << filename << std::endl;
    exit(EXIT_FAILURE);
  }

  idx_base = filename.rfind('/');
  if(idx_base != std::string::npos) {
    p.stub = filename.substr(idx_base+1, idx_ext-(idx_base+1));
    p.directory = filename.substr(0, idx_base);
  } else {
    p.stub = filename.substr(0, idx_ext);
  }

  return p;
}

int generate_hsx(std::vector<std::string>& filenames, std::string hsx_file) {
  // Check that we don't have more than 255 files
  int file_count = filenames.size();
  if (file_count > MAX_FILES) {
    std::cerr << "Maximum number of files (255) exceeded" << std::endl;
    exit(EXIT_FAILURE);
  }

  // Generate the file information records
  std::set<std::string> seen_filenames;
  std::vector<std::string> file_information_records;
  for (auto &filename : filenames) {

    // Ensure there are no duplicate filenames
    if (seen_filenames.find(filename) != seen_filenames.end()) {
      std::cerr << "Hsx does not support duplicates, but duplicate found "
                << filename << std::endl;
      exit(EXIT_FAILURE);
    }
    seen_filenames.insert(filename);

    Path p = parse_path(filename);
    std::string ext = p.extension;
    std::string base_name = p.stub;

    // check that extension is fasta or fa
    std::string lowercase_extension = lowercase_string(ext);
    if (lowercase_extension != "fa" && lowercase_extension != "fasta") {
      std::cerr << "Hsx only support fasta or fa extensions. Found "
                << ext << " from " << filename << "." << std::endl;
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
    s.name = std::string(i->header);
    s.length = i->sequence_length;
    s.file_num = i->file_number;
    s.offset = i->header_offset;
    index.push_back(s);
  }

  // sort index
  std::sort(index.begin(), index.end(), sorter);

  // Create a hsx object and write a hsx file
  Hsx obj(file_count, num_sequences, num_buckets, file_information_records, index);
  obj.write_file(hsx_file);


  /*
  Hsx obj_read;
  obj_read.read_file(hsx_file);
  obj_read.info();
  */


  return 0;
}
