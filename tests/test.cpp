#include <UnitTest++.h>
#include "../source/local_alignment.hpp"
#include "../source/global_alignment.hpp"
#include "../source/local_global_alignment.hpp"
#include "../source/uncertain_local_global_alignment.hpp"
#include "../source/align_lib.hpp"

using namespace std;

TEST(NeedlemanWunsch)
{
  cout << "Testing Global Alignment...\n";

  TSequence seq1 = "GGGCCCCTGCGTCCGACCG";
  TSequence seq2 = "TCTA";
  Score<int, Simple> score(5, -5, -2);
  TAlign test_align;
  TAlign corr_align;
  resize(rows(corr_align), 2);
  assignSource(row(corr_align,0),seq1);
  assignSource(row(corr_align,1),seq2);
  

  int my_score = GlobalAlignment::needle(test_align, seq1, seq2);
  int official_score = globalAlignment(corr_align, score);
  float pm = AlignLib::percent_match(corr_align, test_align);

  cout << test_align << endl;
  //cout << "Percent Match = " << pm << endl;

  CHECK(pm == 1.0);
}

TEST(SmithWaterman)
{
  cout << "Testing Local Alignment...\n";

  TSequence seq1 = "AGCGTAGCTCGTCATTGCT";
  TSequence seq2 = "TTGCT";
  Score<int, Simple> score(5, -5, -2);
  TAlign test_align;
  TAlign corr_align;
  resize(rows(corr_align), 2);
  assignSource(row(corr_align,0),seq1);
  assignSource(row(corr_align,1),seq2);
  

  int my_score = LocalAlignment::smith_waterman(test_align, seq1, seq2);
  int official_score = localAlignment(corr_align, score);
  float pm = AlignLib::percent_match(corr_align, test_align);

  cout << test_align << endl;
  //cout << "Percent Match = " << pm << endl;

  CHECK(pm == 1.0);
}

TEST(LocalGlobalAlignment)
{
  cout << "Testing Local Global Alignment ...\n";

  cout << "LGA equals basic local alignment\n";
  TSequence seq1 = "AGCGTAGCTCGTCATTGCT";
  TSequence seq2 = "CTCGTC";
  Score<int, Simple> score(5, -5, -2);
  TAlign test_align_eq;
  TAlign corr_align_eq;
  resize(rows(corr_align_eq), 2);
  assignSource(row(corr_align_eq,0),seq1);
  assignSource(row(corr_align_eq,1),seq2);
  

  int my_score = LocalGlobalAlignment::lga(test_align_eq, seq1, seq2);
  int official_score = localAlignment(corr_align_eq, score);
  float pm = AlignLib::percent_match(corr_align_eq, test_align_eq);

  cout << test_align_eq << endl;
  //cout << "Percent Match = " << pm << endl;

  CHECK(pm == 1.0);

  cout << "LGA not equal local alignment\n";
  seq1 = "GCCCCCACCCCCCT";
  seq2 = "GGACT";
  TAlign test_align_neq;
  TAlign corr_align_neq;
  resize(rows(corr_align_neq), 2);
  assignSource(row(corr_align_neq,0),seq1);
  assignSource(row(corr_align_neq,1),seq2);
  

  my_score = LocalGlobalAlignment::lga(test_align_neq, seq1, seq2);
  official_score = localAlignment(corr_align_neq, score);
  pm = AlignLib::percent_match(corr_align_neq, test_align_neq);

  cout << "Local Alignment: \n" << corr_align_neq << endl;
  cout << "LGA Alignment: \n" << test_align_neq << endl;
  //cout << "Percent Match = " << pm << endl;

  CHECK(pm != 1.0);

}

TEST(UncertainLocalGlobalAlignment)
{
  cout << "Testing Uncertain Local Global Alignment ...\n";

  cout << "ULGA equals ULGA\n";
  TSequence seq1 = "AGCGTAGCTCGTCATTGCT";
  TSequence seq2 = "CTCGTC";
  CharString qual1 = "135231!";
  TAlign ulga_align;
  TAlign lga_align;

  int ulga_score = UncertainLocalGlobalAlignment::ulga(ulga_align, seq1, seq2, qual1);
  int lga_score = LocalGlobalAlignment::lga(lga_align, seq1, seq2);
  float pm = AlignLib::percent_match(lga_align, ulga_align);

  cout << ulga_align << endl;
  //cout << "Percent Match = " << pm << endl;

  CHECK(pm == 1.0);

  cout << "ULGA not equal LGA\n";
  seq1 = "GCCCACCCCT";
  seq2 = "GGAAAACT";
  qual1 = "91111919";
  qual1[1] = (char) ( (int)qual1[1] - 40 );
  qual1[2] = (char) ( (int)qual1[2] - 40 );
  qual1[3] = (char) ( (int)qual1[3] - 40 );
  qual1[4] = (char) ( (int)qual1[4] - 40 );
  qual1[6] = (char) ( (int)qual1[6] - 40 );

  TAlign ulga_align_neq;
  TAlign lga_align_neq;
 
  ulga_score = UncertainLocalGlobalAlignment::ulga(ulga_align_neq, seq1, seq2, qual1);
  lga_score = LocalGlobalAlignment::lga(lga_align_neq, seq1, seq2);

  cout << "LGA Alignment: \n" << lga_align_neq << endl;
  cout << "ULGA Alignment: \n" << ulga_align_neq << endl;
  //cout << "Percent Match = " << pm << endl;

  CHECK(ulga_align_neq != lga_align_neq);
}

int main()
{
  return UnitTest::RunAllTests();
}
