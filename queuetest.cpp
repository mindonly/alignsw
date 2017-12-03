#include <iostream>
#include <queue>

using namespace std;

int main() {
    queue<int> iq;
    iq.push(1);
    iq.push(2);
    iq.push(3);

    cout << "queue (iq):" << endl;
    cout << "front element: " << iq.front() << "\nback element: " << iq.back() << endl;
   
    iq.pop();
    cout << "queue (iq):" << endl;
    cout << "front element: " << iq.front() << "\nback element: " << iq.back() << endl;


}
