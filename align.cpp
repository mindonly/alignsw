#include <iostream>
#include <iterator>
#include <fstream>
#include <string>
#include <vector>
#include <thread>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

#define RESET "\e[m"
#define GREEN "\e[32m"
#define RED   "\e[31m"

// #define DEBUG

using namespace std;
using namespace boost::numeric::ublas;

std::vector<char> importSeqFile(const std::string &filename) {
    ifstream inFile(filename, ios::binary);
    std::vector<char> fileContents((istreambuf_iterator<char>(inFile)),
                                    istreambuf_iterator<char>());

    return fileContents;
}

void printSeq(const std::vector<char> &seq) {
    cout << "\ncontents = ";
    for (auto& ch : seq)
        cout << ch;
    cout << endl;
}

int gapNorth(matrix<int> &mat, int row, int col) { 
    if (row == 0 || col == 0) {
        cerr << "\nerror: nucleotide coordinates cannot be zero.\n";
        exit(-1);
    }

    #ifdef DEBUG
        cout << "\ngapNorth():\n";
        cout << "(" << row << ", " << col << ")\n";
    #endif
    
    // North-adjusted row
    int adjRow = row - 1;

    #ifdef DEBUG
        cout << "(" << adjRow << ", " << col << ")\n";
        cout << mat(adjRow, col) - 2;
    #endif

    return mat(adjRow, col) - 2;
}

int gapWest(matrix<int> &mat, int row, int col) {
    if (row == 0 || col == 0) {
        cerr << "\nerror: nucleotide coordinates cannot be zero.\n";
        exit(-1); 
    }

    #ifdef DEBUG
        cout << "\ngapWest():\n";
        cout << "(" << row << ", " << col << ")\n";
    #endif

    // West-adjusted col
    int adjCol = col - 1;
    
    #ifdef DEBUG
        cout << "(" << row << ", " << adjCol << ")\n";
        cout << mat(row, adjCol) - 2;
    #endif

    return mat(row, adjCol) - 2;
}

int match(std::vector<char> &s, 
          std::vector<char> &t, 
          int row, int col) {
    
    if (row == 0 || col == 0) {
        cerr << "\nerror: nucleotide coordinates cannot be zero.\n";
        exit(-1); 
    }

    #ifdef DEBUG
        cout << s[row-1] << " " << t[col-1] << endl;
    #endif
    
    if (s[row-1] == t[col-1])
        return 1;
    else
        return -1;
}

int northWest(matrix<int> &mat,
              std::vector<char> &s,
              std::vector<char> &t,
              int row, int col) {

    if (row == 0 || col == 0) {
        cerr << "\nerror: nucleotide coordinates cannot be zero.\n";
        exit(-1); 
    }

    #ifdef DEBUG
        cout << "\nnorthWest():\n";
        cout << "(" << row << ", " << col << ")\n";
    #endif

    int adjRow = row - 1;
    int adjCol = col - 1;

    #ifdef DEBUG
        cout << "(" << adjRow << ", " << adjCol << ")\n";
        cout << mat(adjRow, adjCol) + match(s, t, row, col) << endl;
    #endif

    return mat(adjRow, adjCol) + match(s, t, row, col);
}

void updateScore(matrix<int> &mat,
                 std::vector<char> &s,
                 std::vector<char> &t,
                 int row, int col) {

    if (row == 0 || col == 0) {
        cerr << "\nerror: nucleotide coordinates cannot be zero.\n";
        exit(-1); 
    }
    
    #ifdef DEBUG
        cout << "\nupdateScore():\n";
        cout << "(" << row << ", " << col << ")\n";
    #endif
    
    std::vector<int> scores;
    scores.push_back(gapNorth(mat, row, col));
    scores.push_back(gapWest(mat, row, col));
    scores.push_back(northWest(mat, s, t, row, col));

    #ifdef DEBUG
        cout << "\n\nscores:\n";
        for (int i = 0; i < scores.size(); i++)
            cout << scores[i] << " " << endl;

        cout << "\nmaximum score:\n";
    #endif
    
    auto top_score = max_element(scores.begin(), scores.end());
    
    if (*top_score < 0)
        mat(row, col) = 0;
    else
        mat(row, col) = *top_score;
}

int main(int argc, char* argv[]) {
    if (argc != 3) {
        cerr << "usage: align sequence_file unknown_file\n";
        exit(-1);
    }
    string seqFilenm = argv[1];
    string unkFilenm = argv[2];

    std::vector<char> s = importSeqFile(seqFilenm);
    s.shrink_to_fit();
    s.pop_back();
    cout << "\nSEQUENCE: " << seqFilenm << " size: " << s.size();
    printSeq(s);

    std::vector<char> t = importSeqFile(unkFilenm);
    t.shrink_to_fit();
    t.pop_back();
    cout << "\nUNKNOWN: " << unkFilenm << " size: " << t.size();
    printSeq(t);

    // create and zero-out similarity matrix
    matrix<int> sim_mat(s.size() + 1, t.size() + 1);
    sim_mat.clear();

    for (int i = 1; i <= s.size(); i++) 
        for (int j = 1; j <= t.size(); j++) {
            updateScore(sim_mat, s, t, i, j);
            // cout << endl;
        }

    cout << endl;
    cout << sim_mat << endl;
}