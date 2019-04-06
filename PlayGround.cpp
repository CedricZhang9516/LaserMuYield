//#include <boost/numeric/odeint.hpp>
#include <iostream>
#include <sstream>
#include <string>

using namespace std;

void PlayGround(){

	const char* filename = "test.root";
	//char* 
	stringstream filename2;filename2<<"filename2";
	
	string filename3;// = "filename3";
	
	//filename3.str(filename);
	//const char *s = "Hello, World!";
	std::string str(filename);

	filename3 = filename;

	//filename = filename3.c_str();
	//filename2<<filename;
	//string filename2 = "_" + ss;
	//string filename3 = "test2.root";
	
	//filename2.str(filename3);
	//cout<<filename3.c_str()<<endl;
	cout<<filename3<<endl;
	cout<<"odeint"<<endl;
}

int main(){
	
	return 0;
}
