#include "local_alignment.hpp"
using namespace std;

int LocalAlignment::match_cost(TSequence seq1, TSequence seq2, int i, int j, Score<int, Simple> scores) {
  if (seq1[i-1] == seq2[j-1]) return scoreMatch(scores);
  else return scoreMismatch(scores);
}

int LocalAlignment::max(int * array, int size) {
  int i = 0;
  int max = 0;
  for (i = 0; i < size; i++){
    if(array[i] > max) {
      max = array[i];
    }
  }
  return max;
}

int LocalAlignment::smith_waterman(TSequence seq1, TSequence seq2, Score<int, Simple> scores){
  TAlign align;
  resize(rows(align), 2);
  assignSource(row(align,0),seq1);
  assignSource(row(align,1),seq2);
  TRow &row1 = row(align,0);
  TRow &row2 = row(align,1);

  int seq1Len = length(seq1);
  int seq2Len = length(seq2);
  int matrix[seq1Len + 1][seq2Len + 1];
  int currentCosts[3];
  int i = 0;
  int j= 0;

  int max_i = 0;
  int max_j = 0;
  for(i =0; i < seq1Len + 1; i++) {
    for(j = 0; j < seq2Len + 1; j++) {
      matrix[i][j] = 0;
    }
  }

  /*for(i =0; i < seq1Len + 1; i++) {
    for(j = 0; j < seq2Len + 1; j++) {
      cout << matrix[i][j] << " ";
    }
    cout << "\n";
  }*/

  for(i = 1; i < seq1Len + 1; i++) {
    for(j = 1; j < seq2Len + 1; j++) {
      currentCosts[DIAGONAL] = matrix[i-1][j-1] + LocalAlignment::match_cost(seq1, seq2, i, j, scores);
      currentCosts[TOP] = matrix[i-1][j] + scoreGap(scores);
      currentCosts[LEFT] = matrix[i][j-1] + scoreGap(scores);
      matrix[i][j] = LocalAlignment::max(currentCosts, 3);
      if(matrix[i][j] > matrix[max_i][max_j]) {
        max_i = i;
        max_j = j;
      }
    }
  }

/*
  for(i =0; i < seq1Len + 1; i++) {
    for(j = 0; j < seq2Len + 1; j++) {
      cout << matrix[i][j] << " ";
    }
    cout << "\n";
  }
  */

  i = max_i;
  j = max_j;
  //need to refactor this backtrace as it's pretty much a copy paste of the other one
  int pos, match, mismatch, vgap, hgap;
  while ( (i > 0 || j > 0) && matrix[i][j] > 0)
  {
    pos = matrix[i][j];
    match = matrix[i-1][j-1] + scoreMatch(scores);
    mismatch = matrix[i-1][j-1] + scoreMismatch(scores);
    vgap = matrix[i][j-1] + scoreGap(scores);
    hgap = matrix[i-1][j] + scoreGap(scores);
    //printf("Pos:%d, Match:%d, Current:%d\n", pos, match, matrix[i][j]);
    if (i > 0 && j > 0 && (pos == match || pos == mismatch))
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

  cout << align << "\n";
  return matrix[max_i][max_j];
}
