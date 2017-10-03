/*
 * align2.cpp
 * --
 * Rob Sanchez
 * Parallelized Local Sequence Alignment
 * programming assignment #2
 * CIS 677, F2017
 * Wolffe
 * --
 * multi-threaded readyqueue version
*/


#include <iostream>
#include <fstream>
#include <vector>
#include <thread>
#include <mutex>
#include <queue>
#include <chrono>
#include <utility>
#include <algorithm>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

#define GAP_PENALTY 2
#define MATCH_BONUS 1

using namespace std;
using namespace boost::numeric::ublas;


std::queue<pair<int, int>> rq;     // readyqueue
std::mutex rq_mutex;

/*
 * Timer class pilfered from https://gist.github.com/gongzhitaao/7062087
 * Timer() constructs the timer.
 * .reset() resets the timer.
 * .elapsed() returns elapsed seconds (double) since last reset.
 */
class Timer
{
public:
    Timer() : beg_(clock_::now()) {}
    void reset() { beg_ = clock_::now(); }
    double elapsed() const { 
        return std::chrono::duration_cast<second_>
            (clock_::now() - beg_).count(); }

private:
    typedef std::chrono::high_resolution_clock clock_;
    typedef std::chrono::duration<double, std::ratio<1> > second_;
    std::chrono::time_point<clock_> beg_;
};

/* 
 * import a nucleotide sequence file
 * returns: [vector<char>]
 */
std::vector<char> importSeqFile(const string &filename) {
    ifstream inFile(filename, ios::binary);
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
 * print a uBLAS matrix<tup<int, int>> (tuple)
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
 * retrieve S-W score for a North cell, minus GAP_PENALTY
 * returns: [int]
 */
int North(const matrix<int> &smat, int row, int col) { 
    if (row == 0 || col == 0) {
        cerr << "\nNorth() error: nucleotide coordinates cannot be zero.\n";
        exit(-1);
    }

    // North-adjusted row
    int adjRow = row - 1;

    return smat(adjRow, col) - GAP_PENALTY;
}

/*
 * retrieve S-W score for a West cell, minus GAP_PENALTY
 * returns: [int]
 */
int West(const matrix<int> &smat, int row, int col) {
    if (row == 0 || col == 0) {
        cerr << "\nWest() error: nucleotide coordinates cannot be zero.\n";
        exit(-1); 
    }

    // West-adjusted col
    int adjCol = col - 1;

    return smat(row, adjCol) - GAP_PENALTY;
}

/*
 * determine if two nucleotides match, see MATCH_BONUS
 * returns: [int]
 */
int similarity(const std::vector<char> &s, 
               const std::vector<char> &t, 
               int row, int col) {
    
    if (row == 0 || col == 0) {
        cerr << "\nsimilarity() error: nucleotide coordinates cannot be zero.\n";
        exit(-1); 
    }

        // s and t are 0-indexed
    if ( (s[row-1] == t[col-1]) || s[row-1] == '?' || t[col-1] == '?' )
        return MATCH_BONUS;
    else
        return MATCH_BONUS * -1;
}

/*
 * return a Smith-Waterman score for a NorthWest cell
 * returns: [int]
 */
int NorthWest(const matrix<int> &smat,
              const std::vector<char> &s,
              const std::vector<char> &t,
              int row, int col) {

    if (row == 0 || col == 0) {
        cerr << "\nNorthWest() error: nucleotide coordinates cannot be zero.\n";
        exit(-1); 
    }

    // North-West adjustment
    int adjRow = row - 1;
    int adjCol = col - 1;

    return smat(adjRow, adjCol) + similarity(s, t, row, col);
}

/*
 * compute the source cell for a Smith-Waterman score
 * return [tuple<int, int>]
 */
tuple<int, int> source(int idx, int row, int col) {
    if (row == 0 || col == 0) {
        cerr << "\nsource() error: nucleotide coordinates cannot be zero.\n";
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
            cerr << "\nsource() error: invalid index.\n";
            break;
    }

    return make_tuple(row, col);
}

/*
 * update Smith-Waterman score for each sim. matrix (smat) cell;
 * also update source for each score in tuple matrix (tmat)
 * threads: use mutex to lock global readyqueue
 */
void SmithWaterman(matrix<int> &smat, matrix<tuple<int, int>> &tmat,
                   const std::vector<char> &s,
                   const std::vector<char> &t,
                //    std::queue<tuple<int, int>> &rq,
                //    std::vector<tuple<int, int>> &rq,
                   int row, int col) {

    if (row == 0 || col == 0) {
        cerr << "\nSmithWaterman() error: nucleotide coordinates cannot be zero.\n";
        exit(-1); 
    }
   
        // set up vector<int> of neighbor scores
    std::vector<int> scores;
    scores.push_back(North(smat, row, col));
    scores.push_back(NorthWest(smat, s, t, row, col));
    scores.push_back(West(smat, row, col));

        // get largest neighbor score
    auto top_score = max_element(scores.begin(), scores.end());
    
    if (*top_score < 0)
        *top_score = 0;

        // update similarity matrix
    smat(row, col) = *top_score;
        
        // get top score index; update tuple matrix with source() 
    int  top_index = distance(scores.begin(), top_score);
    tmat(row, col) = source(top_index, row, col);

    rq_mutex.lock();
            // push neighbors onto ready queue
        if ( (row == s.size() && col != t.size()) || (row != s.size() && col != t.size()) ) {
            auto east_n  = make_tuple(row, col+1);
            rq.push(east_n); 
        }
        if (col == 1 && row != s.size()) {
            auto south_n = make_tuple(row+1, col);
            auto east_n  = make_tuple(row, col+1);
            rq.push(south_n);
            rq.push(east_n);
        } 
    rq_mutex.unlock();    
}

/*
 * find the largest, rightmost, lowermost score in a 
 * uBLAS S-W similarity matrix and its [x, y] coordinates
 * returns: [tuple<int, int, int>]
 */
tuple<int, int, int> maxScore(const matrix<int> &smat) {
    int s_sz = smat.size1() - 1;
    int t_sz = smat.size2() - 1;
    
    int cur_max = 0;
    int row = s_sz;
    int col = t_sz;

    for (int i = s_sz; i >= 1; i--)
        for (int j = t_sz; j >= 1; j--) {
            if (smat(i, j) > cur_max) {
                cur_max = smat(i, j);
                row = i;
                col = j;
            }
        }

    return make_tuple(cur_max, row, col);
}

/*
 * wrap SmithWaterman() in a thread-join. 
 */ 
void threadedSW(matrix<int> &smat, matrix<tuple<int, int>> &tmat,
                const std::vector<char> &s,
                const std::vector<char> &t,
                //    std::queue<tuple<int, int>> &rq,
                //    std::vector<tuple<int, int>> &rq,
                int row, int col) {

    std::thread th(SmithWaterman, std::ref(smat), std::ref(tmat),
                 std::ref(s), std::ref(t), row, col);
    th.join();
}

/*
 *
 */
void traceback(const matrix<tuple<int, int>> &tmat, const tuple<int, int> &p) {
    std::vector<tuple<int, int>> route;
    route.push_back(p);

    auto cell = tmat(get<0>(p), get<1>(p));

    while (cell != make_tuple(0, 0)) {
       route.push_back(cell);
       cell = tmat(get<0>(cell), get<1>(cell)); 
    }

    std::reverse(route.begin(), route.end());

    for (auto it=route.begin(); it!= route.end(); it++)
        cout << "[" << get<0>(*it) << ", " << get<1>(*it) << "] ";
    cout << endl;
}

/*
 * main program
 */
int main(int argc, char* argv[]) {
    if (argc != 3) {
        cerr << "usage: align sequence_file unknown_file\n";
        exit(-1);
    }

        // start the timer
    Timer tmr;

        // command-line args
    string seqFilNam = argv[1];
    string unkFilNam = argv[2];

        // import sequences
    std::vector<char> s = importSeqFile(seqFilNam);
    s.shrink_to_fit();
    s.pop_back();
    cout << "\nSEQUENCE(S): " << seqFilNam << " size: " << s.size();
    // printSeq(s);

    std::vector<char> t = importSeqFile(unkFilNam);
    t.shrink_to_fit();
    t.pop_back();
    cout << "\n UNKNOWN(T): " << unkFilNam << " size: " << t.size();
    // printSeq(t);

        // create and zero-out similarity matrix
    matrix<int> sim_mat(s.size() + 1, t.size() + 1);
    sim_mat.clear();

        // create and zero-out traceback() tuple matrix
    matrix<tuple<int, int>> tup_mat(s.size() + 1, t.size() + 1);
    tup_mat.clear();

        // mark sim. matrix cells as not ready (except row, col = 0) 
    for (int i = 1; i <= s.size(); i++)
        for (int j = 1; j <= t.size(); j++)
            sim_mat(i, j) = -999;

        // main task:
        // kick off multi-threaded sequence
        // main() pop()s the readyqueue while
        // SmithWaterman() push()es onto it
    int row = 1;
    int col = 1;
    auto seed = make_pair(row, col);

    rq.push(seed);
    while (! rq.empty()) {
        
        rq_mutex.lock();
            auto p = rq.front();
            row = p.first;
            col = p.second;
            rq.pop();
        rq_mutex.unlock();
    
        threadedSW(sim_mat, tup_mat, s, t, row, col);
    }

    // cout << endl;
    // printSimMatrix(sim_mat);
    // cout << endl;
    // printTupMatrix(tup_mat);

        // retrieve max score and output its location
    auto tup = maxScore(sim_mat);
    cout << "\n\n(" << get<0>(tup) << ", [" << get<1>(tup) << ", " << get<2>(tup) << "])\n";
    cout << "similarity matrix dims: (" << sim_mat.size1() << "x" << sim_mat.size2() << ")" << endl;
   
        // stop the timer
    double elapsed = tmr.elapsed();
    cout << "\nelapsed time: " << elapsed << " seconds." << endl;

        // print the traceback path
    auto maxop = make_tuple(get<1>(tup), get<2>(tup));
    traceback(tup_mat, maxop);
}