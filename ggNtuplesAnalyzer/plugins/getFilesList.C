#include <iostream>
#include <fstream>
using namespace std;

vector<string>  getFilesList (string filesList) {

  vector<string> filesVector;

  string line;
  ifstream filesData (filesList);
  if (filesData.is_open()) {
    while ( getline (filesData,line) ) {
      filesVector.push_back(line);
    }
    filesData.close();
  } else {
    cout << "Unable to open file" << endl;
  } 

  return filesVector;
}