#pragma once
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
//#include <iostream>
#include <unordered_map>
#include <set>
#include <vector>
#include <map>
#include <iterator>
#include <algorithm>
#include <cassert>
using namespace std;
/*template <size_t N>
  void splitString(string (&arr)[N], string str)
{
    int n = 0;
    istringstream iss(str);
    for (auto it = istream_iterator<string>(iss); it != istream_iterator<string>() && n < N; ++it, ++n)
        arr[n] = *it;
}*/
class HyperGraphLoad {
public:
    string dataname;
    string inputpath;
    int number_of_hedges;
    int number_of_nodes;
    vector< vector<int> > node2hyperedge; 
	vector< vector<int> > hyperedge2node;
   // unordered_map<int, int> nodename2index;
    //unordered_map<int, int> index2nodename;
    unordered_map<string, int> nodename2index;
    unordered_map<int, string> index2nodename;
    unordered_map<int, string> edgename;
    bool exist_edgename;
    
  
   // HyperGraph(string inputpath, string dataname);
    HyperGraphLoad(string inputpath, string dataname){
    string path;
    this->inputpath = inputpath;
    this->dataname = dataname;
      path = inputpath + dataname + ".txt";
        this->exist_edgename = false;

	ifstream graphFile(path.c_str());
	string line;
	int num_hyperedge = 0;
    int num_nodes = 0;
    cout << "Start Read Dataset" << endl; 
    
    //std::string delimiter = ",";
    std::string delimiter = " ";
    size_t pos = 0;
    vector<int> tokens;
	while (getline(graphFile, line)){
        string ename;
		vector<int> nodes;
        tokens.clear();
        std::string token;
        pos = 0;
		while (pos != line.npos) {
			pos = line.find(delimiter);
			token = line.substr(0, pos);
			nodes.push_back( atoi(token.c_str()));
			line.erase(0, pos + delimiter.length());
		}
		for (int i = 0; i < nodes.size(); i++){
			// tokens.push_back(stoi(nodes[i]));
            if (nodename2index.find(to_string(nodes[i])) == nodename2index.end()){
                int index = num_nodes++;
                nodename2index[to_string(nodes[i])] = index;
                index2nodename[index] = to_string(nodes[i]);
                this->node2hyperedge.push_back(vector<int>());
            }
            int node_index = nodename2index[to_string(nodes[i])];
            tokens.push_back(node_index);
            this->node2hyperedge[node_index].push_back(num_hyperedge);
        }
        sort(tokens.begin(), tokens.end());
        this->hyperedge2node.push_back(tokens);
        if (this->exist_edgename){
            edgename[num_hyperedge] = ename;
        }
      
        num_hyperedge++;
	}
    this->number_of_hedges = num_hyperedge;
    this->number_of_nodes = (int)this->node2hyperedge.size();
    cout << "Load " << number_of_hedges << " hyperedges" << endl;

   unordered_map<string, long long> tmp;
   // unordered_map<string, pair<int,int>> tmp_com_vector_size;
   unordered_map<string, pair<double,double>> tmp_com_vector_size;
   unordered_map<int, long long> intersect;
   //	time_t start = 0,end = 0;  
   // time(&start);	
    // intersect
    cout << "graphInfo" << endl;
   // string its_name = "graphInfo/" + dataname + "_graph";
   string its_name = "graphInfo/" + dataname;
    cout << "Load " << number_of_hedges << " hyperedges++++++" << endl;
        for(int h = 0 ; h < number_of_hedges ; h++){
            int hsize = (int)hyperedge2node[h].size();
            for(int i = 0 ; i < hsize-1 ; i++){
                for(int j = (i + 1) ; j < hsize ; j++){
                    int va = hyperedge2node[h][i];
                    int vb = hyperedge2node[h][j];
                    string key;
                    if (va > vb){
                        key = to_string(vb) + " " + to_string(va);
                    }
                    else{
                        key = to_string(va) + " " + to_string(vb);
                    }
                    tmp[key]++;
                  //  tmp_com_vector[key].push_back(h);
                  //  tmp_com_vector[key]++;
                }
            }
        }
        std::string token;
        ofstream resultFile(its_name.c_str());
       // std::string delimiter = " ";
       // size_t pos = 0;
        // string line;

       vector<int> arr;
    cout << "Load " << number_of_hedges << " hyperedges********************" << endl;
    resultFile <<this->number_of_nodes<<" "<<tmp.size()*2<<endl;
        for(auto elem : tmp)
        {
           // stringstream ssin(elem.first);
           // ssin >> arr[0];
           // ssin >> arr[1];
           istringstream ss(elem.first);
           int i;
         while(ss>>i)
        {
            arr.push_back(i);
        }
        
          // cout << elem.first<<endl;
          // splitString(arr, elem.first);
           
            int deg_1 = (int)node2hyperedge[arr[0]].size();
            int deg_2= (int)node2hyperedge[arr[1]].size();
            resultFile << arr[0]<<" "<<arr[1]<<" "<<deg_1<<endl;
            resultFile << arr[1]<<" "<<arr[0]<<" "<<deg_2<<endl;
           // cout << arr[0]<<" "<<arr[1]<<" "<<deg_1<<endl;
            //cout << arr[1]<<" "<<arr[0]<<" "<<deg_2<<endl;
            arr.clear();
              // << elem.first<<" "<<elem.second<<" "<node2hyperedge[elem.second]<<endl;
          //  cout << elem.second<<" "<<elem.first<<" "<<node2hyperedge[elem.first].size()<<endl;
         } 
    }
void analyzeSeed(string inputpath)
{
    std::map<int,int> *tran=new std::map<int,int>;
    tran->clear();
    uint32_t srcId;
    string  path = "result/seed/seed_" + inputpath;
    ifstream seedFile(path .c_str());
    if (!seedFile.is_open())
        {
            std::cout << "The file \"" + path + "\" can NOT be opened\n";
            return;
        }
        int totalVCnt=0;
        int totalVcnt_size=0;
        int totalHe=0;
    for (auto i = 100; i--;)
      {
        seedFile >> srcId;
        uint32_t cntHyper=this->node2hyperedge[srcId].size();
        totalHe=totalHe+(int)(this->node2hyperedge[srcId].size());
        totalVCnt=0;
        for(int j=0;j<cntHyper;j++)
        {
            int he=this->node2hyperedge[srcId][j];
            totalVCnt=totalVCnt+(int)hyperedge2node[he].size();
            totalVcnt_size=totalVcnt_size+(int)hyperedge2node[he].size();

        }
       // (*tran)[cntHyper]=(*tran)[cntHyper]+1;
        (*tran)[totalVCnt]=(*tran)[totalVCnt]+1;
      }
    float ratio=0.0;
    string its_name_r = "result/" + inputpath+"_dis";
    cout<<"totalHe: "<<totalHe<<" totalVcnt_size: "<<totalVcnt_size<<endl;
    ofstream resultFile_r(its_name_r.c_str());
 for(map<int,int>::iterator iter=(*tran).begin();iter!=(*tran).end();++iter)
	 {
	   ratio=(float)iter->second/100;
       resultFile_r<< iter->first<<" "<<iter->second<<" "<<ratio<<endl;
	 }
    delete tran;   

}
    
};