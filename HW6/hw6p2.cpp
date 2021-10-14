#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

int main()
{
    vector<int> vec;
    ofstream out("size.csv");
    out << "size\tcapacity\n";

    for(int i = 0; i < 999; i++)
    {
        vec.push_back(i);
        out << vec.size() << "\t" << vec.capacity() << "\n";
    }
}