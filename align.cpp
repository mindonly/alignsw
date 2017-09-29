#include <iostream>
#include <iterator>
#include <fstream>
#include <string>
#include <vector>
#include <thread>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

#define GAP_PENALTY 2
#define MATCH_BONUS 1 

using namespace std;
using namespace boost::numeric::ublas;

/* 
 * import a nucleotide sequence file
 * returns: [vector<char>]
 */
std::vector<char> importSeqFile(const string &filename) {
    // ifstream inFile(filename, ios::binary);
    ifstream inFile(filename, ios::in);
    std::vector<char> fileContents( (istreambuf_iterator<char>(inFile)),
                                     istreambuf_iterator<char>() );

    return fileContents;
}

/*
 * print a nucleotide sequence
 */
void printSeq(const std::vector<char> &seq) {
    cout << "\ncontents = ";
    for (auto &ch : seq)
        cout << ch;
    cout << endl;
}

/*
 * print a uBLAS matrix<int> (similarity)
 */ 
void printSimMatrix(const matrix<int> &mat) {
    for (int i = 0; i < mat.size1(); i++) {
        for (int j = 0; j < mat.size2(); j++) {
                cout << mat(i, j) << " ";
        }
        cout << endl;
    }
}

/*
 * print a uBLAS matrix<tup> (tuple)
 */
void printTupMatrix(const matrix<tuple<int, int>> &mat) {
    for (int i = 0; i < mat.size1(); i++) {
        for (int j = 0; j < mat.size2(); j++) {
            auto tup = mat(i, j);
            cout << "[" << get<0>(tup) << ", " << get<1>(tup) << "] ";
        }
        cout << endl;
    }
}

/*
 * return a Smith-Waterman score for a North cell
 * returns: [int]
 */
int North(const matrix<int> &mat, int row, int col) { 
    if (row == 0 || col == 0) {
        cerr << "\nerror: nucleotide coordinates cannot be zero.\n";
        exit(-1);
    }

    // North-adjusted row
    int adjRow = row - 1;

    return mat(adjRow, col) - 2;
}

/*
 * return a Smith-Waterman score for a West cell
 * returns: [int]
 */
int West(const matrix<int> &mat, int row, int col) {
    if (row == 0 || col == 0) {
        cerr << "\nerror: nucleotide coordinates cannot be zero.\n";
        exit(-1); 
    }

    // West-adjusted col
    int adjCol = col - 1;

    return mat(row, adjCol) - 2;
}

/*
 * determine if two nucleotides match
 * returns: [int]
 */
int match(const std::vector<char> &s, 
          const std::vector<char> &t, 
          int row, int col) {
    
    if (row == 0 || col == 0) {
        cerr << "\nerror: nucleotide coordinates cannot be zero.\n";
        exit(-1); 
    }

    // s and t are 0-indexed
    if ( s[row-1] == t[col-1] || s[row-1] == '?' || t[col-1] == '?' )
        return 1;
    else
        return -1;
}

/*
 * return a Smith-Waterman score for a NorthWest cell
 * returns: [int]
 */
int NorthWest(const matrix<int> &mat,
              const std::vector<char> &s,
              const std::vector<char> &t,
              int row, int col) {

    if (row == 0 || col == 0) {
        cerr << "\nerror: nucleotide coordinates cannot be zero.\n";
        exit(-1); 
    }

    // North-West adjustment
    int adjRow = row - 1;
    int adjCol = col - 1;

    return mat(adjRow, adjCol) + match(s, t, row, col);
}

/*
 * compute the source cell for a Smith-Waterman score
 * return [tuple<int, int>]
 */
tuple<int, int> source(int idx, int row, int col) {
    if (row == 0 || col == 0) {
        cerr << "\nerror: nucleotide coordinates cannot be zero.\n";
        exit(-1); 
    }

    switch (idx) {
        case 0:     // North
            row = row - 1;
            break;
        case 1:     // NW
            row = row - 1;
            col = col - 1;
            break;
        case 2:     // West
            col = col - 1;
            break;
        default:
            cerr << "\nerror: invalid input.\n";
            break;
    }

    return make_tuple(row, col);
}

/*
 * update Smith-Waterman score for similarity matrix cell
 */
void SmithWaterman(matrix<int> &mat, matrix<tuple<int, int>> &trc,
                 const std::vector<char> &s,
                 const std::vector<char> &t,
                 int row, int col) {

    if (row == 0 || col == 0) {
        cerr << "\nerror: nucleotide coordinates cannot be zero.\n";
        exit(-1); 
    }
    
    std::vector<int> scores;
    scores.push_back(North(mat, row, col));
    scores.push_back(NorthWest(mat, s, t, row, col));
    scores.push_back(West(mat, row, col));

    auto top_score = max_element(scores.begin(), scores.end());
    auto top_index = distance(scores.begin(), top_score);
    
    if (*top_score < 0)
        *top_score = 0;

    mat(row, col) = *top_score;

    trc(row, col) = source(top_index, row, col);
}

/*
 * find the largest, rightmost, lowermost score in a 
 * uBLAS S-W similarity matrix and its [x, y] coordinates
 * returns: [tuple<int, int, int>]
 */
tuple<int, int, int> maxScore(const matrix<int> &mat) {
    int x_sz = mat.size1() - 1;
    int y_sz = mat.size2() - 1;
    
    int cur_max = mat(x_sz, y_sz);
    int row = 0;
    int col = 0;

    for (int i = x_sz; i >= 1; i--)
        for (int j = y_sz; j >= 1; j--) {
            if (mat(i, j) > cur_max) {
                cur_max = mat(i, j);
                row = i;
                col = j;
            }
        }

    return make_tuple(cur_max, row, col);
}


/*
 * main program
 */
int main(int argc, char* argv[]) {
    if (argc != 3) {
        cerr << "usage: align sequence_file unknown_file\n";
        exit(-1);
    }
    string seqFilNam = argv[1];
    string unkFilNam = argv[2];

    std::vector<char> s = importSeqFile(seqFilNam);
    s.shrink_to_fit();
    s.pop_back();
    cout << "\nSEQUENCE(S): " << seqFilNam << " size: " << s.size();
    printSeq(s);

    std::vector<char> t = importSeqFile(unkFilNam);
    t.shrink_to_fit();
    t.pop_back();
    cout << "\nUNKNOWN(T): " << unkFilNam << " size: " << t.size();
    printSeq(t);

        // create and zero-out similarity matrix
    matrix<int> sim_mat(s.size() + 1, t.size() + 1);
    sim_mat.clear();

        // mark all sim. matrix cells except first row 
        // and first column as not ready (-999)
    for (int i = 1; i <= s.size(); i++)
        for (int j = 1; j <= t.size(); j++)
            sim_mat(i, j) = -999;

    // printSimMatrix(sim_mat);

        // create and zero-out traceback() tuple matrix
    matrix<tuple<int, int>> tup_mat(s.size() + 1, t.size() + 1);
    tup_mat.clear();

    // printTupMatrix(tup_mat);

        // for each cell in the similarity matrix, 
        // compute and update the S-W score
    for (int i = 1; i <= s.size(); i++) 
        for (int j = 1; j <= t.size(); j++)
            SmithWaterman(sim_mat, tup_mat, s, t, i, j);

    cout << endl;
    // printSimMatrix(sim_mat);
    cout << endl;
   
    auto tup = maxScore(sim_mat);
    cout << "(" << get<0>(tup) << ", [" << get<1>(tup) << ", " << get<2>(tup) << "])\n";
    cout << "similarity matrix dims: (" << sim_mat.size1() << "x" << sim_mat.size2() << ")" << endl;
   
    /*
    cout << "[1751, 48] " << sim_mat(1751, 48) << endl;
    cout << "[1685, 48] " << sim_mat(1685, 48) << endl; 
    */
}