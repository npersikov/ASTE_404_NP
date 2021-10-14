#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

struct Container
{
  int i;
  int num_copies;
  Container(int i) : i{i}, num_copies{0} {}  // use initializer list to set i
  Container(const Container &o) : i{o.i}, num_copies{o.num_copies + 1} {}   // copy constructor
  void operator=(const Container &o) {i=o.i;} // copy assignment
};

int main() 
{
  vector<Container> vec;   // array of container objects

  /*
    When this line is used, the number of copies outputted is zero,
    as opposed to the plot in part 2c.
  */
  vec.reserve(20); 
  
  for (int i=0;i<20;i++) 
  {
    vec.emplace_back(i);   // emplace_back calls the constructor
  }
  
  int* copyTotals = new int[20]; // Make an array to store numbers of copy operations

  int x = 0;
  for (const Container& cont:vec) 
  {
    cout << cont.i << " ";
    // copyTotals[x] = cont.num_copies;
  }
  cout << endl;
  for(int index = 0; index < 20; index++)
    cout << vec.at(index).num_copies << " ";

  cout << endl;  // new line
  return 0;
}