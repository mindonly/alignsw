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

/* int gapWest() {

}

int gapNorth() {

} */

int match(std::vector<char> &s, 
           std::vector<char> &t, 
           int row, int col) {
    
    cout << s[row-1] << " " << t[col-1] << endl;
    
    if (s[row-1] == t[col-1])
        return 1;
    else
        return -1;
}

/* int northWest() {

} */

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
}