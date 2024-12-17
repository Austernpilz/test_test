
//markov Chain almost finished
//build class object tiles next
//SeqAn3 Alphabets, hashing, vectors and strings, sequence_object???

#pragma once

#include "kmc/include/kmc_runner.h"
#include "kmc/kmc_api/kmc_file.h"

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/views/to_char.hpp>
#include <seqan3/io/sequence_file/output.hpp>
#include <seqan3/io/sequence_file/record.hpp>
#include <seqan3/alphabet/views/char_to.hpp>
#include <seqan3/io/sequence_file/format_fasta.hpp>


#include <filesystem>
#include <algorithm>
#include <iostream>
#include <string>
#include <map>
#include <random>
#include <vector>

#include <fstream>
#include <sstream>

// #include <typeinfo>


// void print_info(void);




struct configuration
{
    std::filesystem::path fastq_input{};
    std::filesystem::path fasta_output{};
    bool verbose{}; // Default is false.
};

std::map<char, uint32_t> translate = 
{
  {'a', 0}, {'A', 0}, {'c', 1}, {'C', 1}, {'g', 2}, {'G', 2}, {'t', 3}, {'T', 3} 
};

uint32_t hash2(const std::string& kmer)
{
  uint32_t returnvalue = translate[kmer[0]];

  for (uint8_t i = 1; i<kmer.size(); ++i)
  {
    returnvalue <<= 2;
    returnvalue |= translate[kmer[i]];
  }

  return returnvalue;
}


class Markov_Chain_3
{
  private:
    std::map<uint16_t, std::string> kmer_map;

  public:
        
    void build_kmer_db(std::string input_file_name, auto input_file_type, std::string output_file_name, int threads)
    {
      try
      {
        KMC::Runner runner;

        KMC::Stage1Params stage1Params;

        stage1Params
            .SetKmerLen(4)
            .SetInputFiles({input_file_name})
            .SetNThreads(threads)
            .SetInputFileType(input_file_type) 
            //KMC::InputFileType::FASTA, KMC::InputFileType::MULTILINE_FASTA, KMC::InputFileType::FASTQ
            .SetCanonicalKmers(false)
            ;

        auto stage1Result = runner.RunStage1(stage1Params);

        std::cout << "Stage 1 finished \n";


        KMC::Stage2Params stage2Params;

        stage2Params
            .SetOutputFileName(output_file_name);

        auto stage2Result = runner.RunStage2(stage2Params);


        std::cout << "Stage 2 finished \n";
      }

      catch(const std::exception& e)
      {
        std::cerr << e.what() << '\n';
      }
    }

          
    std::map<uint8_t, std::map<char, double> > build_markov_chain(std::string output_file_name, int threads)
    {
      CKMCFile KMC_Database; 

      if (!KMC_Database.OpenForListing(output_file_name))
      {
        std::cerr << "Couldn't build KMC_Database" << '\n';
      }

      else
      {
        std::map<uint8_t, std::map<char, double> > markov_chain;
        std::vector<uint32_t> kmer_overall;
        kmer_overall.resize(64, 4);

        for (uint8_t i = 0; i < 64; ++i)
        {
          markov_chain[i] = {{'A', 1}, {'C', 1}, {'G', 1}, {'T', 1} };
        }

        CKmerAPI kmer_object(4);
        std::string str;
        uint32_t counter;

        while (KMC_Database.ReadNextKmer(kmer_object, counter))
        {
          kmer_object.to_string(str);

          std::transform(str.begin(), str.end(), str.begin(), ::toupper);

          uint32_t hashval = hash2(str);
          kmer_map[hashval] = counter;

          hashval = hash2(str.substr(0,3));
          kmer_overall[hashval] += counter;

          markov_chain[hashval][str[3]] += counter;
        }

        for (uint8_t i = 0; i < 64; ++i)
        {
          markov_chain[i]['A'] /= kmer_overall[i];
          markov_chain[i]['C'] /= kmer_overall[i];
          markov_chain[i]['G'] /= kmer_overall[i];
          markov_chain[i]['T'] /= kmer_overall[i];
        }

        std::cout << "Markov_Chain build \n"; 

        uint32_t i = 0;
        for(const auto& elem : kmer_overall)
        {
          std::cout << i << " | " << elem << " " << "\n";
          ++i;
        }


        return markov_chain;
      }
    }

};


void generate_sequence(std::map<uint8_t, std::map<char, double> > & markovchain, seqan3::dna4_vector & testsequence, size_t length)
{
  

    for (size_t i = 0; i<length; ++i)
    {   

      std::vector<double> probal;
      std::vector<char> picks;

      using namespace seqan3::literals;
      // seqan3::dna4_vector s{"ACTTTGATAA"_dna4};
      using iterator = seqan3::dna4_vector::iterator;
      auto v1 = std::ranges::subrange<iterator, iterator>{std::ranges::end(testsequence) - 3, std::ranges::end(testsequence)}
              | seqan3::views::to_char; // == "TTTGATAA"

      // std::vector<char> testsequence = {'A', 'C', 'G', 'T'};
      std::string str << v1;
      uint32_t hashval = hash2(str);
      
      for(const auto& elem : markovchain[hashval])
      {
        picks.push_back(elem.first);
        probal.push_back(elem.second);
      }
      
      //probably slow but good
      std::mt19937 gen(std::random_device{}());
      std::discrete_distribution<std::size_t> d{probal.begin(), probal.end()};
  
      auto s = picks[d(gen)];
      // std::cout<< sampled_value << '\n';

      auto s_as_dna = s | seqan3::views::char_to<seqan3::dna4>;

      testsequence.push_back(s_as_dna); 
    }
}


void print(seqan3::dna4_vector & testsequence)
{
  
  auto fasta_file = std::filesystem::current_path() / "my.fasta"; 

  using namespace seqan3::literals;

  seqan3::sequence_file_output fout{fasta_file, seqan3::format_fasta{}};
 
  using types = seqan3::type_list<std::vector<seqan3::dna4>, std::string>;
  using fields = seqan3::fields<seqan3::field::seq, seqan3::field::id>;
  using sequence_record_type = seqan3::sequence_record<types, fields>;


  std::string id{"test_id"};
  sequence_record_type record{std::move(testsequence), std::move(id)};
 
  fout.push_back(record);
}

int main(){
  Markov_Chain_3 testchain;

  std::string input_file_name = "/group/ag_abi/manuel/testdata/GCF_917046035.1_P2_DC1_genomic.fna";
  KMC::InputFileType input_file_type = KMC::InputFileType::MULTILINE_FASTA;
  std::string output_file_name = "4mers";
  int threads = 16;

  testchain.build_kmer_db(input_file_name, input_file_type, output_file_name, threads);

  std::map<uint8_t, std::map<char, double> > markovchain;
  markovchain = testchain.build_markov_chain(output_file_name, threads);

  uint32_t i = 0;
  for(const auto& elem : markov_chain)
  {
    std::cout << elem.first << " | " << i << "\n";
    ++i;
    for(const auto& elem2 : elem.second)
    {
      std::cout << elem2.first << " | " << elem2.second << " " << "\n";
    }
  }

  // std::vector<char> testsequence = {'A', 'C', 'G', 'T'};

  seqan3::dna4_vector testsequence = "ACG"_dna4;
  size_t length = 10000000;


  // std::cout << "Testsequence before " << '\n';
  // for (auto& elem : testsequence)
  // {
  //     std::cout << elem;
  // }

  generate_sequence(markovchain, testsequence, length);


  print(testsequence);



  Markov_Chain_3 testchain_new;

  auto input_file_name_new = std::filesystem::current_path() / "my.fasta"; 

  KMC::InputFileType input_file_type = KMC::InputFileType::MULTILINE_FASTA;
  std::string output_file_name_new = "4mers_new";

  testchain_new.build_kmer_db(input_file_name_new, input_file_type, output_file_name_new, threads);

  std::map<uint8_t, std::map<char, double> > markovchain_new;
  markovchain_new = testchain_new.build_markov_chain(output_file_name_new, threads);

  uint32_t j = 0;
  for(const auto& e : markovchain_new)
  {
    std::cout << e.first << " | " << i << "\n";
    ++j;
    for(const auto& e2 : e.second)
    {
      std::cout << e2.first << " | " << e2.second << " " << "\n";
    }
  }
}
  // std::cout << "Testsequence after " << '\n';
  // for (auto& elem : testsequence)
  // {
  //     std::cout << elem;
  // }

  // std::cout << '\n' <<"Complete" << '\n';


  // std::string input_SNP_file = "/group/ag_abi/manuel/testdata/gruppe6chr1.vcf";

  // std::vector<SNP> testsnp = read_SNPs(input_SNP_file);

  // for (auto& elem : testsnp)
  // {
  //     std::cout << elem.chr << '\t' << elem.pos << '\t' << elem.ref << '\t' << elem.alt << '\n';
  // }

  // std::cout << '\n' <<"Complete2" << '\n';

  // std::vector<std::vector<char>> test_kmer_subset; = 
  // {
  //   {'C', 'C', 'A', 'G', 'C', 'G', 'G', 'A', 'G', 'T', 'A', 'A', 'T'},
  //   {'T', 'G', 'G', 'T'},
  //   {'A', 'A', 'C', 'T', 'A', 'G', 'C'},
  //   {'C', 'A', 'G', 'A'},
  //   {'T', 'C', 'G', 'A', 'C', 'G', 'T', 'C', 'G', 'T'},
  //   {'A', 'T', 'C', 'T'},
  //   {'T', 'A', 'G', 'T', 'T', 'G', 'T', 'A', 'G', 'T', 'A', 'G', 'T', 'A', 'G', 'T'},
  //   {'T', 'C', 'G', 'T', 'A', 'C', 'C', 'C', 'A', 'T'},
  //   {'A', 'A', 'C', 'G', 'C', 'C', 'T', 'G', 'G', 'A', 'C', 'G', 'T'},
  //   {'A', 'C', 'A', 'T', 'C', 'C', 'T'}
  // };




  // class Tile
  // {
  //   std::string name;
  //   size_t size;
  //   // border;
  //   std::vector<std::vector<char>>kmer_subset;
  //   std::string path_reference_sequence;
  //   std::map<uint8_t, std::map<char, double> > markovchain;
  //   // Polymorphisms SNPs, INDELs, CNVs;

    
  // }



  // std::map<char, std::map<char, double> > testmarkovchain
  // {
  //     {'A', {{'A', 0.1},{'C', 0.4}, {'G', 0.2}, {'T', 0.3} } },
  //     {'C', {{'A', 0.3},{'C', 0.2}, {'G', 0.3}, {'T', 0.2} } },
  //     {'G', {{'A', 0.7},{'C', 0.1}, {'G', 0.1}, {'T', 0.1} } },
  //     {'T', {{'A', 0.1},{'C', 0.2}, {'G', 0.6}, {'T', 0.1} } }
  // };






  


// // struct SNP
// // {
// //   std::string chr;
// //   std::string pos;
// //   std::string ref;
// //   std::string alt;
// // };


// std::vector<SNP> read_SNPs(std::string path_to_SNP)
// {
//   std::vector<SNP> SNP_table;
//   std::ifstream fin(path_to_SNP);
//   std::string line;

//   while (std::getline(fin, line))
//   {
//     if (line[0] == '#')
//     {
//       continue;
//     }

//     std::vector<std::string> line_vec;
//     std::istringstream iss(line);
//     std::string field;

//     while(std::getline(iss, field, '\t'))
//         line_vec.push_back(field);

//     // for (auto& elem : testsequence)
//     // {
//     //   std::cout << line_vec[0] << '\t' << line_vec[1] << '\t' << line_vec[3] << '\t' << line_vec[4] << '\n';
//     // }

//     SNP SNP_field;
//     SNP_field.chr = line_vec[0];
//     SNP_field.pos = line_vec[1];
//     SNP_field.ref = line_vec[3];
//     SNP_field.alt = line_vec[4];

//     SNP_table.push_back(SNP_field);
//   }
  
//   fin.close();

//   return SNP_table;
// }


  
  // for (auto& elem : testsequence)
  // {
  //     std::cout << elem;
  // }


  // for(const auto& elem : kmer_map)
  // {
  //   std::cout << elem.first << " | " << elem.second << " " << "\n";
  // }

  // uint32_t i = 0;
  // for(const auto& elem : kmer_overall)
  // {
  //   std::cout << i << " | " << elem << " " << "\n";
  //   ++i;
  // }

  // uint32_t i = 0;
  // for(const auto& elem : markov_chain)
  // {
  //   std::cout << elem.first << " | " << i << "\n";
  //   ++i;
  //   for(const auto& elem2 : elem.second)
  //   {
  //     std::cout << elem2.first << " | " << elem2.second << " " << "\n";
  //   }
  // }
    
    // bool OpenForRA(const std::string &file_name)
    // CKMCFile() 
    //bool GetCountersForRead(const std::string& read, std::vector<uint32_t>& counters)
  // }

//https://stackoverflow.com/questions/167735/fast-pseudo-random-number-generator-for-procedural-content
// fast random number gen???
// v = 36969*(v & 65535) + (v >> 16);
// u = 18000*(u & 65535) + (u >> 16);
// return (v << 16) + (u & 65535);