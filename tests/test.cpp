#include <UnitTest++.h>
#include "../source/local_alignment.hpp"
using namespace std;
TEST(SmithWaterman)
{
  TSequence seq1 = "AGCGTAG";
  TSequence seq2 = "CTCGTC";
  Score<int, Simple> score(10, -5, -7);
  int my_score = LocalAlignment::smith_waterman( "AGCGTAG", "CTCGTC", score);
  TAlign align;
  resize(rows(align), 2);
  assignSource(row(align,0),seq1);
  assignSource(row(align,1),seq2);
  int official_score = localAlignment(align, score);
  CHECK(my_score == official_score);
}

int main()
{
  return UnitTest::RunAllTests();
}
