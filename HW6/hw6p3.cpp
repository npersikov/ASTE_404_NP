#include <iostream>
#include <thread>
using namespace std;
struct MyObject{       // some custom object
  MyObject(int i) {
    size_t N=1024*1024*20;
    data = new double[N];  // some large buffer
    for (int n=0;n<N;n++) data[n]=i;
 }
 ~MyObject()
 {
     delete[] data;
 }

 double val() {return data[0];}
 protected:
  double *data;
};
int main() {
  for (int i=0;i<20;i++) {
    MyObject object(i);
    cout<<object.val()<<endl;  // placeholder for some action done
    this_thread::sleep_for(3s);// wait for 3 second before continuing
  }   // object goes out of scope here
  return 0;
}