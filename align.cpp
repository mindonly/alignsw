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

    cout << "\ngapNorth():\n";
    cout << "(" << row << ", " << col << ")\n";
    // North-adjusted row
    int adjRow = row - 1;
    cout << "(" << adjRow << ", " << col << ")\n";

    cout << mat(adjRow, col) - 2;
    return mat(adjRow, col) - 2;
}

int gapWest(matrix<int> &mat, int row, int col) {
    if (row == 0 || col == 0) {
        cerr << "\nerror: nucleotide coordinates cannot be zero.\n";
        exit(-1); 
    }

    cout << "\ngapWest():\n";
    cout << "(" << row << ", " << col << ")\n";
    // West-adjusted col
    int adjCol = col - 1;
    cout << "(" << row << ", " << adjCol << ")\n";

    cout << mat(row, adjCol) - 2;
    return mat(row, adjCol) - 2;
}

int match(std::vector<char> &s, 
          std::vector<char> &t, 
          int row, int col) {
    
    if (row == 0 || col == 0) {
        cerr << "\nerror: nucleotide coordinates cannot be zero.\n";
        exit(-1); 
    }
            
    // cout << s[row-1] << " " << t[col-1] << endl;
    
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

    cout << "\nnorthWest():\n";
    cout << "(" << row << ", " << col << ")\n";
    int adjRow = row - 1;
    int adjCol = col - 1;
    cout << "(" << adjRow << ", " << adjCol << ")\n";

    cout << mat(adjRow, adjCol) + match(s, t, row, col);
    return mat(adjRow, adjCol) + match(s, t, row, col);
}

/* void computePos() {

} */

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
    matrix<int> sim_mat(s.size(), t.size());
    sim_mat.clear();

    for (int i = 1; i <= 2; i++) 
        for (int j = 1; j <= 2; j++) {
            gapNorth(sim_mat, i, j);
            gapWest(sim_mat, i, j);
            northWest(sim_mat, s, t, i, j);
            cout << endl;
        }
}