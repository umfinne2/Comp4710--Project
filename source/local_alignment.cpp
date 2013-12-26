#include "local_alignment.hpp"
#include "align_lib.hpp"

using namespace std;

float LocalAlignment::max(float * array, int size) {
  int i = 0;
  float max = 0;
  for (i = 0; i < size; i++){
    if(array[i] > max) {
      max = array[i];
    }
  }
  return max;
}

int LocalAlignment::smith_waterman(TAlign &align, TSequence seq1, TSequence seq2){
  int score = 0;
  resize(rows(align), 2);
  assignSource(row(align,0),seq1);
  assignSource(row(align,1),seq2);
  TRow &row1 = row(align,0);
  TRow &row2 = row(align,1);

  int seq1Len = length(seq1);
  int seq2Len = length(seq2);
  //int matrix[seq1Len + 1][seq2Len + 1];
  float currentCosts[3];
  int i = 0;
  int j= 0;

  int max_i = 0;
  int max_j = 0;

  float **matrix = new float * [seq1Len + 1];
  for (i = 0; i <= seq1Len; i++)
  {
    matrix[i] = new float[seq2Len + 1];
  }

  for(i = 0; i < seq1Len + 1; i++) {
    matrix[i][0];
  }

  for(j = 0; j < seq2Len + 1; j++) {
    matrix[0][j] = 0;
  }


  for(i = 1; i < seq1Len + 1; i++) {
    for(j = 1; j < seq2Len + 1; j++) {
      currentCosts[DIAGONAL] = matrix[i-1][j-1] + AlignLib::get_score(seq1[i-1], seq2[j-1]);
      currentCosts[TOP] = matrix[i-1][j] + AlignLib::gapcost;
      currentCosts[LEFT] = matrix[i][j-1] + AlignLib::gapcost;
      matrix[i][j] = LocalAlignment::max(currentCosts, 3);
      if(matrix[i][j] > matrix[max_i][max_j]) {
        max_i = i;
        max_j = j;
      }
    }
  }

  score = matrix[max_i][max_j];

  i--;
  j--;

  //since we know that seq2 will be shorter than seq1 we can just figure out how many gaps to insert
  //before we start the local alignment given the max_i + max_j
  int num_gaps = (i - max_i) - (j - max_j);
  insertGaps(row2, j, num_gaps);

  i = max_i;
  j = max_j;
  //need to refactor this backtrace as it's pretty much a copy paste of the other one
  float pos, dmap, vgap, hgap;
  while ( (i > 0 && j > 0) && matrix[i][j] > 0)
  {
    pos = matrix[i][j];
    dmap = matrix[i-1][j-1] + AlignLib::get_score(seq1[i-1], seq2[j-1]);
    vgap = matrix[i][j-1] + AlignLib::gapcost;
    hgap = matrix[i-1][j] + AlignLib::gapcost;
    //printf("Pos:%d, Match:%d, Current:%d\n", pos, match, matrix[i][j]);
    if (i > 0 && j > 0 && pos == dmap)
    {
      i--;
      j--;
    }

    else if (j > 0 && pos == vgap)
    {
      insertGap(row1, i);
      j--;
    }

    else if (i > 0 && pos == hgap)
    {
      insertGap(row2, j);
      i--;
    }

    else
    {
      cout << "You distroyed the universe! You're on your own now :( \n";
      exit(1);
    }
  }

  //insert some more gaps after we're done
  insertGaps(row2, 0, (i - j));

  //delete the dp matrix
  for (i = 0; i <= seq1Len; i++)
  {
    delete matrix[i];
  }

  delete[] matrix;

  //cout << align << "\n";
  return score;
}
