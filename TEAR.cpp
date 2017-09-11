#include <iostream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <map>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <random>
#include <time.h>
#include <math.h>
#include <algorithm>


using namespace std;

#define MAX 99999
string dataDirectory, demandDirectory;
string resultdir;
int nFile;
int linkCap, nodeCap;
vector<vector<double> > uvStructPairs;
vector<vector<vector<double> > > uvDemandPairs;
vector<vector<vector<double> > > uvDemandPairsS;
vector<vector<double> > demandIndex(0, vector<double>(0));
int N, M;
int unusedw;
ofstream outputFile;
double beta;
double gama;
double link_p;



vector<vector<double> > loadStructure() {
	vector<vector<double> > uvStructurePair;
	string fileName = resultdir+string("structure.txt");
	ifstream inputFile(fileName);
	string line;

	//get #node and #link in first line
	getline(inputFile, line);
	istringstream ss(line);
	ss >> N >> M;

	//get graph structure
	while (getline(inputFile, line)) {
		istringstream ss(line);
		vector<double> temp;
		double u, v, c;
		ss >> u >> v;
		temp.push_back(u);
		temp.push_back(v);
		uvStructurePair.push_back(temp);
	}
	return uvStructurePair;
}

void loadDemand() {
	for (int i = 1; i <= nFile; i++) {
		string fileName = demandDirectory + "data_" + to_string(i) + ".txt";
		cout << fileName << endl;
		ifstream inputFile(fileName);
		string line;
		vector<vector<double> > uvDemandPair;
		while (getline(inputFile, line)) {
			istringstream ss(line);
			vector<double> temp; //source, target, traffic
			double u, v, d,tag;
			ss >> u >> v >> d >> tag;
			temp.push_back(u);
			temp.push_back(v);
			temp.push_back(d);
			temp.push_back(tag);
			uvDemandPair.push_back(temp);
		}
		uvDemandPairs.push_back(uvDemandPair);
	}
}

void getDemandIndex() {
	string s = resultdir+"demandIndex.txt";
	ifstream inputFile(s);
	string line;
	while (getline(inputFile, line)) {
		istringstream ss(line);
		vector<double> temp; //source, target, traffic
		double u, v;
		ss >> u >> v;
		temp.push_back(u);
		temp.push_back(v);
		demandIndex.push_back(temp);
	}
	/*demandIndex.clear();
	vector<double> temp;
	temp.push_back(5);
	temp.push_back(9);
	temp.push_back(5);
	demandIndex.push_back(temp);

	temp.clear();
	temp.push_back(4);
	temp.push_back(9);
	temp.push_back(4);
	demandIndex.push_back(temp);

	temp.clear();
	temp.push_back(1);
	temp.push_back(9);
	temp.push_back(6);
	demandIndex.push_back(temp);

	temp.clear();
	temp.push_back(1);
	temp.push_back(10);
	temp.push_back(1);
	demandIndex.push_back(temp);*/
}

vector<int> SP(vector<vector<double> > adj, int A, int B) {
	vector<int> dist(N, 1e9);
	vector<int> parent(N, -1);
	vector<bool> isVisit(N, false);

	dist[A] = 0; //source dist is 0
	parent[A] = A; //source parent is source
	for (int i = 0; i < N; i++) {
		int csrn = -1; //current shortest reachable node

		double max = 1e9;
		//search cur node
		for (int j = 0; j < N; j++) {
			if (!isVisit[j] && dist[j] < max) {
				max = dist[j];
				csrn = j;
			}
		}
		//if (csrn == B) break;
		isVisit[csrn] = true;
		//update dist
		for (int j = 0; j < N; j++) {
			int path = dist[csrn] + adj[csrn][j];
			if (path < dist[j]) {
				dist[j] = path;
				parent[j] = csrn;
			}
		}
	}
	int x = B;
	vector<int> P;
	while (parent[x] != x) {
		P.push_back(x);
		x = parent[x];
	}
	P.push_back(A);
	//get shortest path from source
	reverse(P.begin(), P.end());

	//system("pause");
	return P;

}

bool pathIsReachable(vector<vector<int> > G, vector<int> P) {
	for (int i = 0; i < P.size() - 1; i++) {
		if (G[P[i]][P[i + 1]] == 0){
			//outputFile << "demand: " << P[0] <<" "<< P[P.size()-1] << ":";
			//outputFile << "link: "<< P[i] << " "<< P[i+1] << " ";
			return false;
		}	
			
	}
	return true;
}

int findDemandIndex(int s, int t) {
	for (int i = 0; i < demandIndex.size(); i++) {
		if (demandIndex[i][0] == s && demandIndex[i][1] == t)
			return i;
	}
	return -1;
}

int findMostFrequentPort(vector<tuple<int, int, int> > table) {

	map<int, int> occur;
	for (int i = 0; i < table.size(); i++) {
		occur[get<2>(table[i])]++;
	}
	int max = -1;
	int port = -1;
	map<int, int>::iterator it;
	for (it = occur.begin(); it != occur.end(); it++) {
		if (it->second > max) {
			max = it->second;
			port = it->first;
		}
	}
	return port;
}

int OurFindMostFrequentPort(vector<vector<int> > table) {

	map<int, int> occur;
	for (int i = 0; i < table.size(); i++) {
		occur[table[i][2]]++;
	}
	int max = -1;
	int port = -1;
	map<int, int>::iterator it;
	for (it = occur.begin(); it != occur.end(); it++) {
		if (it->second > max) {
			max = it->second;
			port = it->first;
		}
	}
	return port;
}

void shrinkTable(vector<tuple <int, int, int> >& Fu, vector<tuple <int, int, int> >& Gu) {
	int port = findMostFrequentPort(Fu);

	//logFile << "before shrink table:" << endl;
	for (int i = 0; i < Fu.size(); i++) {
		//logFile << get<0>(Fu[i]) << " " << get<1>(Fu[i]) << " " << get<2>(Fu[i]) << endl;
	}

	//logFile << "most frequent port: " << port << endl;

	for (int k = 0; k < Fu.size(); k++) {
		if (get<2>(Fu[k]) == port) {
			Gu.push_back(Fu[k]);
			Fu.erase(Fu.begin() + k);
			k--;
		}
	}

	//logFile << "after shrink table:" << endl;
	//logFile << "--Fu" << endl;
	for (int i = 0; i < Fu.size(); i++) {
		//logFile << get<0>(Fu[i]) << " " << get<1>(Fu[i]) << " " << get<2>(Fu[i]) << endl;
	}
	//logFile << "--Gu" << endl;
	for (int i = 0; i < Gu.size(); i++) {
		//logFile << get<0>(Gu[i]) << " " << get<1>(Gu[i]) << " " << get<2>(Gu[i]) << endl;
	}
}
void OurShrinkTable(vector<vector<int> >& ruleSpace, vector<vector<int> >& DefaultruleSpace){
	int port = OurFindMostFrequentPort(ruleSpace);

	for (int k = 0; k < ruleSpace.size(); k++) {
		if (ruleSpace[k][2] == port) {
			DefaultruleSpace.push_back(ruleSpace[k]);
			ruleSpace.erase(ruleSpace.begin() + k);
			k--;
		}
	}

}

//input variable exist & ==> call by reference
void updatingFuAndGu(vector<int>& path, vector<vector<tuple <int, int, int> > >& Fu, vector<vector<tuple<int, int, int> > >& Gu) {
	for (int i = 0; i < path.size() - 1; i++) {
		tuple<int, int, int> t;
		t = tuple<int, int, int>{ path[0], path[path.size() - 1], path[i + 1] };
		//node u already assigned default port
		if (Gu[path[i]].size() != 0) {
			if (get<2>(t) == get<2>(Gu[path[i]][0])) {
				Gu[path[i]].push_back(t);
			}
			else {
				Fu[path[i]].push_back(t);
			}
		}
		else {
			Fu[path[i]].push_back(t);
		}
	}
}
int NewFindingFeasibleRouting(vector<vector<int> > G, vector<vector<double> > &Re, vector<vector<double> > D, vector<vector<int> > &p, vector<vector<double> > &WTemp,vector<vector<tuple<int, int, int> > > &FuTemp,vector<vector<tuple<int, int, int> > > &GuTemp,vector<bool> &isShrinkedTemp) {
	vector<vector<double> > flowOnLink(uvStructPairs.size(), vector<double>());
	vector<vector<double> > W(N, vector<double>(N, MAX));//weight
	vector<vector<tuple<int, int, int> > > Fu(N, vector<tuple<int, int, int> >(0, tuple<int, int, int>()));
	vector<vector<tuple<int, int, int> > > Gu(N, vector<tuple<int, int, int> >(0, tuple<int, int, int>()));
	vector<bool> isShrinked(N, false); //node is already shrinked or not
	vector<vector<int> > Path; //demands's path

	W=WTemp;
	Fu=FuTemp;
	Gu=GuTemp;
	isShrinked=isShrinkedTemp;

	for (int i = 0; i < D.size(); i++) {
		//logFile << "demand: " << D[i][0] << " " << D[i][1] << " " << D[i][2] << endl;

		

		int deindex=findDemandIndex(D[i][0],D[i][1]);

		//cout << "or: ";
		/*for(int j=0;j<p[deindex].size();j++){
			cout  << p[deindex][j] << " ";
		}
		cout << endl;*/

		
		for(int j=0;j<p[deindex].size()-1;j++) {
			

			Re[p[deindex][j]][p[deindex][j + 1]] += D[i][2];
			Re[p[deindex][j + 1]][p[deindex][j]] += D[i][2];

			bool tag=false;
			for(int k=0;k<Fu[p[deindex][j]].size();k++){
				if( get<0>(Fu[p[deindex][j]][k])==p[deindex][0] && get<1>(Fu[p[deindex][j]][k])==p[deindex][p[deindex].size()-1] ){
					Fu[p[deindex][j]].erase(Fu[p[deindex][j]].begin() + k);
					k--;
					tag=true;
					break ;
				}
			}
			if(Gu[p[deindex][j]].size()!=0 && tag==false){
				for(int k=0;k<Gu[p[deindex][j]].size();k++){
					if( get<0>(Gu[p[deindex][j]][k])==p[deindex][0] && get<1>(Gu[p[deindex][j]][k])==p[deindex][p[deindex].size()-1] ){
						Gu[p[deindex][j]].erase(Gu[p[deindex][j]].begin() + k);
						k--;
						break ;
					}
				}
				if(Gu[p[deindex][j]].size()==0){
					isShrinked[p[deindex][j]]=false;
				}
			}

		}

		for (int u = 0; u < N; u++) {
			for (int v = 0; v < N; v++) {
				if (W[u][v] != MAX ) {
					double Uv = linkCap * Fu[v].size() / linkCap;
					if (Uv > 1) {
						W[u][v] = Uv;
						W[v][u] = Uv;
					}
					else {
						W[u][v] = 1;
						W[v][u] = 1;
					}
				}
			}
		}
		
		


		vector<vector<double> > Wp = W;
		//check link capacity
		for (int j = 0; j < N; j++) {
			for (int k = 0; k < N; k++) {
				if (Re[j][k] < D[i][2] && j != k && G[j][k] == 1) {
					Wp[j][k] = MAX;
					Wp[k][j] = MAX;
					//logFile << "remove link " << j << " " << k << endl;
				}
				if(G[j][k]==0 && G[k][j]==0){
					Wp[j][k] = MAX;
					Wp[k][j] = MAX;
				}
			}
		}
		//logFile << "Wp:" << endl;
		for (int j = 0; j < uvStructPairs.size(); j++) {
			//logFile << Wp[uvStructPairs[j][0]][uvStructPairs[j][1]] << " ";
		}
		//logFile << endl;
		vector<int> path = SP(Wp, D[i][0], D[i][1]);
		////////////////////cout << D[i][0] << " " << D[i][1] << " "<<D[i][2]<< endl;
		//logFile << "path: ";
		//for (int j = 0; j < path.size(); j++) {
		//	cout<< path[j] << " "  ;
		//}
		//cout << endl;

		if (pathIsReachable(G, path)) {
			//check node capactiy
			for (int j = 0; j < path.size(); j++) {

				if (Fu[path[j]].size() == nodeCap) {
					//logFile << "table " << path[j] << "is full, try to shrink" << endl;
					//logFile << "shrink " << path[j]; system("pause");
					if (!isShrinked[path[j]]) {
						isShrinked[path[j]] = true;
						shrinkTable(Fu[path[j]], Gu[path[j]]);
					}
					else {
						//logFile << "table already shrink, no more space" << endl;
						vector<vector<int> > p;
						return -1;
					}
				}
			}
			//logFile << "update Fu and Gu" << endl;
			updatingFuAndGu(path, Fu, Gu);
			//logFile << "table size:" << endl;
			//logFile << "--Fu" << endl;
			for (int j = 0; j < Fu.size(); j++) {
				//logFile << Fu[j].size() << " ";
			}
			//logFile << endl;
			//logFile << "--Gu" << endl;
			for (int j = 0; j < Gu.size(); j++) {
				//logFile << Gu[j].size() << " ";
			}
			//logFile << endl;

			for (int j = 0; j < path.size() - 1; j++) {
				Re[path[j]][path[j + 1]] -= D[i][2];
				Re[path[j + 1]][path[j]] -= D[i][2];
				/*
				for (int k = 0; k < M; k++) {
					if ((path[j] == uvStructPairs[k][0] && path[j + 1] == uvStructPairs[k][1]) || (path[j] == uvStructPairs[k][1] && path[j + 1] == uvStructPairs[k][0]))
						flowOnLink[k].push_back(D[i][2]);
				}*/

			}

			//updating link weights
			for (int u = 0; u < N; u++) {
				for (int v = 0; v < N; v++) {
					if (W[u][v] != MAX) {
						double Uv = linkCap * Fu[v].size() / linkCap;
						if (Uv > 1) {
							W[u][v] = Uv;
							W[v][u] = Uv;
						}
						else {
							W[u][v] = 1;
							W[v][u] = 1;
						}
					}
				}
			}
			p[deindex] = path;
			//Path.push_back(path);
		}//if (pathIsReachable(W, path))
		else {
			//path for D[i] is not reachable
			//logFile << "path not reachable" << endl;
			//cout << "here" << endl;
			return -1;
		}
		/*
		for (int j = 0; j < flowOnLink.size(); j++) {
			//logFile << uvStructPairs[j][0] << " " << uvStructPairs[j][1] << " : ";
			for (int k = 0; k < flowOnLink[j].size(); k++) {
				//logFile << flowOnLink[j][k] << " ";
			}
			//logFile << endl;
		}
		*/
	}
	WTemp=W;
	FuTemp=Fu;
	GuTemp=Gu;
	isShrinkedTemp=isShrinked;
	
	//for (int i = 0; i < demands.size(); i++)
	 //system("pause");
	return 1;
}

vector<vector<double> > GetDemandpair(int index,vector<vector<int> > paths,vector<vector<double> > uvDemandPair){
	vector<vector<double> > result;
	for(int i=0;i<paths.size();i++){
		vector <double> temp;
		if(paths[i].size()!=0){
			int tag=1;
			double src;
			double dst; 
			for(int j=0;j<paths[i].size()-1;j++){
				
				if(paths[i][j]==uvStructPairs[index][0] && paths[i][j+1]==uvStructPairs[index][1] || paths[i][j+1]==uvStructPairs[index][0] && paths[i][j]==uvStructPairs[index][1]){
					tag=0;
					src= paths[i][0];
					dst= paths[i][paths[i].size()-1];
					temp.push_back(src);
					temp.push_back(dst);
					break ;
				}
			}

			if(tag==0){
				for(int j=0;j<uvDemandPair.size();j++){
					if(uvDemandPair[j][0]==src && uvDemandPair[j][1] == dst){
						temp.push_back(uvDemandPair[j][2]);
						break;
					}
				}
				result.push_back(temp);
			}
		}
		//result.push_back(temp);
	}
	return result;

}

int FindingFeasibleRouting(vector<vector<int> > G, vector<vector<double> > &Re, vector<vector<double> > D, vector<vector<int> > &p, vector<vector<double> > &WTemp,vector<vector<tuple<int, int, int> > > &FuTemp,vector<vector<tuple<int, int, int> > > &GuTemp,vector<bool> &isShrinkedTemp) {
	vector<vector<double> > flowOnLink(uvStructPairs.size(), vector<double>());
	vector<vector<double> > W(N, vector<double>(N, MAX)); //weight
	vector<vector<tuple<int, int, int> > > Fu(N, vector<tuple<int, int, int> >(0, tuple<int, int, int>()));
	vector<vector<tuple<int, int, int> > > Gu(N, vector<tuple<int, int, int> >(0, tuple<int, int, int>()));
	vector<bool> isShrinked(N, false); //node is already shrinked or not
	vector<vector<int> > Path; //demands's path

	//construct
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			if (G[i][j] == 1) {
				Re[i][j] = linkCap;
				Re[j][i] = linkCap;
				W[i][j] = 1;
				W[j][i] = 1;
			}
		}
	}

	for (int i = 0; i < D.size(); i++) {

		//logFile << "demand: " << D[i][0] << " " << D[i][1] << " " << D[i][2] << endl;

		vector<vector<double> > Wp = W;

		//check link capacity
		for (int j = 0; j < N; j++) {
			for (int k = 0; k < N; k++) {
				if (Re[j][k] < D[i][2] && j != k && G[j][k] == 1) {
					Wp[j][k] = MAX;
					Wp[k][j] = MAX;
					//logFile << "remove link " << j << " " << k << endl;
				}
			}
		}
		/*
		//logFile << "Wp:" << endl;
		for (int j = 0; j < uvStructPairs.size(); j++) {
			//logFile << Wp[uvStructPairs[j][0]][uvStructPairs[j][1]] << " ";
		}
		//logFile << endl;
		*/
		vector<int> path = SP(Wp, D[i][0], D[i][1]);
		//logFile << "path: ";
		/*
		for (int j = 0; j < path.size(); j++) {
			//logFile << path[j] << " ";
		}
		//logFile << endl;
		*/
		if (pathIsReachable(G, path)) {
			//check node capactiy
			for (int j = 0; j < path.size(); j++) {


				if (Fu[path[j]].size() == nodeCap) {
					//logFile << "table " << path[j] << "is full, try to shrink" << endl;
					//logFile << "shrink " << path[j]; system("pause");
					if (!isShrinked[path[j]]) {
						isShrinked[path[j]] = true;
						shrinkTable(Fu[path[j]], Gu[path[j]]);
					}
					else {
						//logFile << "table already shrink, no more space" << endl;
						vector<vector<int> > p;
						return -1;
					}
				}
			}
			//logFile << "update Fu and Gu" << endl;
			updatingFuAndGu(path, Fu, Gu);
			//logFile << "table size:" << endl;
			//logFile << "--Fu" << endl;
			/*
			for (int j = 0; j < Fu.size(); j++) {
				//logFile << Fu[j].size() << " ";
			}
			//logFile << endl;
			//logFile << "--Gu" << endl;
			for (int j = 0; j < Gu.size(); j++) {
				//logFile << Gu[j].size() << " ";
			}
			//logFile << endl;
			*/

			for (int j = 0; j < path.size() - 1; j++) {

				Re[path[j]][path[j + 1]] -= D[i][2];
				Re[path[j + 1]][path[j]] -= D[i][2];
				/*
				for (int k = 0; k < M; k++) {
					if ((path[j] == uvStructPairs[k][0] && path[j + 1] == uvStructPairs[k][1]) || (path[j] == uvStructPairs[k][1] && path[j + 1] == uvStructPairs[k][0]))
						flowOnLink[k].push_back(D[i][2]);
				}
				*/

			}

			//updating link weights
			for (int u = 0; u < N; u++) {
				for (int v = 0; v < N; v++) {
					if (W[u][v] != MAX) {
						double Uv = linkCap * Fu[v].size() / linkCap;
						if (Uv > 1) {
							W[u][v] = Uv;
							W[v][u] = Uv;
						}
						else {
							W[u][v] = 1;
							W[v][u] = 1;
						}
					}
				}
			}
			p[findDemandIndex(D[i][0], D[i][1])] = path;
			//Path.push_back(path);
		}//if (pathIsReachable(W, path))
		else {
			//path for D[i] is not reachable
			//logFile << "path not reachable" << endl;
			return -1;
		}
		/*
		for (int j = 0; j < flowOnLink.size(); j++) {
			//logFile << uvStructPairs[j][0] << " " << uvStructPairs[j][1] << " : ";
			for (int k = 0; k < flowOnLink[j].size(); k++) {
				//logFile << flowOnLink[j][k] << " ";
			}
			//logFile << endl;
		}*/
		
	}//for (int i = 0; i < demands.size(); i++)
	 //system("pause");
	WTemp=W;
	FuTemp=Fu;
	GuTemp=Gu;
	isShrinkedTemp=isShrinked;
	return 1;
}

void demandleaving(int demandIndex,vector<double> demandinfo, vector<vector<double> > &R,vector<vector<int> > &flowOnLink,vector<vector<int> > &paths,vector<vector<vector<int> > > &ruleSpace,vector<vector<vector<int> > > &DefaultruleSpace){
	if(paths[demandIndex].size()!=0){
		//remove rule from table for each node on path
		for (int j = 0; j < paths[demandIndex].size()-1; j++) {
						
						
			if(DefaultruleSpace[paths[demandIndex][j]].size()!=0 && DefaultruleSpace[paths[demandIndex][j]][0][2]==paths[demandIndex][j+1]){
			////already shrinked, if the demand match the default port , remove it in the default table
				for (int k = 0; k < DefaultruleSpace[paths[demandIndex][j]].size(); k++) {	
							
					if (DefaultruleSpace[paths[demandIndex][j]][k][0] == demandinfo[0] && DefaultruleSpace[paths[demandIndex][j]][k][1] == demandinfo[1]) {
							
						DefaultruleSpace[paths[demandIndex][j]].erase(DefaultruleSpace[paths[demandIndex][j]].begin() + k);
						k--;
					}
				}
			}
			else{
			//not shrink yet or demand doesn't match the default port, directly remove it in normal table
				for (int k = 0; k < ruleSpace[paths[demandIndex][j]].size(); k++) {	
			
					if (ruleSpace[paths[demandIndex][j]][k][0] == demandinfo[0] && ruleSpace[paths[demandIndex][j]][k][1] == demandinfo[1]) {
						
						ruleSpace[paths[demandIndex][j]].erase(ruleSpace[paths[demandIndex][j]].begin() + k);
						k--;
					}
				}

			}
						
		}
		
		//add back residual space, reduce # flow on link for each link on path
		for (int j = 0; j < paths[demandIndex].size() - 1; j++) {
				
			R[paths[demandIndex][j]][paths[demandIndex][j + 1]] += demandinfo[2];
			R[paths[demandIndex][j + 1]][paths[demandIndex][j]] += demandinfo[2];
			flowOnLink[paths[demandIndex][j]][paths[demandIndex][j + 1]]--;
			flowOnLink[paths[demandIndex][j + 1]][paths[demandIndex][j]]--;
				
		}
		paths[demandIndex].clear();
	}
}
///assign weight to the link accroding to the utilization of node and link 
void updatingLp(vector<vector<int> > &G,vector<vector<double> > &Lp,vector<vector<int> > &flowOnLink,vector<vector<double> > &R,vector<vector<vector<int> > > &ruleSpace){
	for (int u = 0; u < N; u++) {
		for (int v = 0; v < N; v++) {
			if (flowOnLink[u][v] == 0 && G[u][v] == 1){
				Lp[v][u] = M;
				Lp[u][v] = M;//uv has link and # flow on uv is 0
			}

			
			double node_p=1-link_p;
			double linkutil=(linkCap-R[u][v])/linkCap;
			double nodeutil=(ruleSpace[u].size()+ruleSpace[v].size())*0.5/nodeCap;

			

			//double weight= link_p*(pow(M,linkutil)) + node_p*(pow(M,nodeutil));
			double weight= link_p*(exp(gama*linkutil)) + node_p*(exp(beta*nodeutil));

			/*
			if(ruleSpace[u].size()>=nodeCap*0.97 || ruleSpace[v].size()>=nodeCap*0.97){
				weight=1;
			}
			*/

			if (flowOnLink[u][v] != 0 && G[u][v] == 1){
				Lp[v][u] = weight;
				Lp[u][v] = weight;//uv has link and # flow on uv is not 0
			}
		}	
	}
}

void demandarriving(int demandIndex,vector<double> demandinfo, vector<vector<double> > &R,vector<vector<int> > &flowOnLink,vector<int> path,vector<vector<vector<int> > > &ruleSpace,vector<vector<vector<int> > > &DefaultruleSpace){
	//substract residual space, increase flowOnlink
	for (int j = 0; j < path.size() - 1; j++) {
		//cout << "before R: " << path[j] << " " << path[j + 1] << ": " << R[path[j]][path[j + 1]] << endl;
		//cout << "before R: " << path[j] << " " << path[j + 1] << ": " << R[path[j]][path[j + 1]] << endl;
		//////cout << "before flowOnLink: " << path[j] << " " << path[j + 1] << ": " << flowOnLink[path[j]][path[j + 1]] << endl;
		//cout << "before flowOnLink: " << path[j] << " " << path[j + 1] << ": " << flowOnLink[path[j]][path[j + 1]] << endl;
		////////coutFile << path[j] << " " << path[j + 1] <<endl;
		R[path[j]][path[j + 1]] -= demandinfo[2];
		flowOnLink[path[j]][path[j + 1]] ++;

		R[path[j + 1]][path[j]] -= demandinfo[2];
		flowOnLink[path[j + 1]][path[j]] ++;

		//cout << "after R: " << path[j] << " " << path[j + 1] << ": " << R[path[j]][path[j + 1]] << endl;
		//cout << "after R: " << path[j] << " " << path[j + 1] << ": " << R[path[j]][path[j + 1]] << endl;
		//cout << "after flowOnLink: " << path[j] << " " << path[j + 1] << ": " << flowOnLink[path[j]][path[j + 1]] << endl;
		//cout << "after flowOnLink: " << path[j] << " " << path[j + 1] << ": " << flowOnLink[path[j]][path[j + 1]] << endl;
	}
	//cout << "add rule for each node on path" << endl;
	//cout << "add rule for each node on path" << endl;
	//add rule for each node on path
	for (int j = 0; j < path.size() - 1; j++) {

		vector<int> temp;
		temp.push_back(demandinfo[0]);
		temp.push_back(demandinfo[1]);
		temp.push_back(path[j + 1]);
		//ruleSpace[path[j]].push_back(temp);
				
		//already shrinked, if the demand match the default port , put into the default table
		if(DefaultruleSpace[path[j]].size()!=0 && DefaultruleSpace[path[j]][0][2]==path[j+1]){
			
			DefaultruleSpace[path[j]].push_back(temp);
		}

		//not shrink yet or demand dosen't match the deafult port, directly put into normal table
		else{
			ruleSpace[path[j]].push_back(temp);
		}
					

		//ruleSpace[path[j]].push_back(temp);	
	}
}

int compare(pair<int,int> elem1 , pair<int,int> elem2){
	return elem1.second>elem2.second;
}
//sort the demands in order to accept more demand because some demand may can not be routed due to table size.
vector<vector<double> > Sortingdemand(vector<vector<double> > uvDemandPair,vector<vector<vector<int> > > ruleSpace){

	vector<vector<double> > uvDemandPairnew;
	vector<pair<int,int> >tablesizelist;
	for(int i=0;i<N;i++){
		tablesizelist.push_back(make_pair(i,ruleSpace[i].size()));
	}
	sort(tablesizelist.begin(),tablesizelist.end(),compare);
	for(int i=0;i<tablesizelist.size();i++){
		if(tablesizelist[i].second==0){
			break;
		}
		for(int j=0;j<uvDemandPair.size();j++){
			if(uvDemandPair[j][0]==tablesizelist[i].first || uvDemandPair[j][1]==tablesizelist[i].first){
				uvDemandPairnew.push_back(uvDemandPair[j]);
				uvDemandPair.erase(uvDemandPair.begin()+j);
				j--;
			}
		}
	}
	if(uvDemandPair.size()!=0){
		for(int i=0;i<uvDemandPair.size();i++){
			uvDemandPairnew.push_back(uvDemandPair[i]);
		}
	}

	return uvDemandPairnew;

}


void oursWithTime() {
	vector<vector<int> > G(N, vector<int>(N, 0)); //adjacency matrix, connect:1 notConnect:0
	vector<vector<double> > R(N, vector<double>(N, 0)); //link recidual space
	vector<vector<int> > flowOnLink(N, vector<int>(N, 0)); //#flows on link
	vector<vector<int> > paths(demandIndex.size(), vector<int>()); //flow path
	
	vector<vector<vector<int> > > ruleSpace(N, vector<vector<int> >());//[source,dest,port]
	vector<vector<vector<int> > > DefaultruleSpace(N, vector<vector<int> >());//shrink the same port to the default port [source,dest,port]

	vector<vector<vector<int> > > pathAll; //flow path
	vector<vector<int> > previousFlowOnLink (N, vector<int>(N, 0));
	int totalOpenCost=0;
	int totalCloseCost=0;
	double totalcount=0;
	double demand_no_path=0;
	
	//construct graph
	for (int i = 0; i < M; i++) {
		G[uvStructPairs[i][0]][uvStructPairs[i][1]] = 1;
		G[uvStructPairs[i][1]][uvStructPairs[i][0]] = 1;

		R[uvStructPairs[i][0]][uvStructPairs[i][1]] = linkCap;
		R[uvStructPairs[i][1]][uvStructPairs[i][0]] = linkCap;
	}
	//cout << uvDemandPairs.size() << endl;
	//algorithm begin
	for (int t = 0; t < uvDemandPairs.size(); t++) {

		clock_t start, finish;
		start = clock();
		cout << "t:" << t + 1 << ": ";
		//outputFile << "t:" << t + 1 << ": ";
		//cout << endl;
		vector<vector<double> > uvDemandPair = uvDemandPairs[t];
		//vector<bool> isInPrevious(uvDemandPair.size(), false);
		int openCost=0;
		int closeCost=0;
		//vector<vector<double> > Outoftable;//to check the demands which are out of table on source or destination
		
		//sorting algorithm

		//uvDemandPair=Sortingdemand(uvDemandPair,ruleSpace);

		/*
		for(int i=0;i<uvDemandPair.size();i++){

			cout << uvDemandPair[i][0] << " "<< uvDemandPair[i][1] << " "<< uvDemandPair[i][2] << " "<< uvDemandPair[i][3]<< " " << endl; 
			//outputFile << uvDemandPair[i][0] << " "<< uvDemandPair[i][1] << " "<< uvDemandPair[i][2]<< " "<< uvDemandPair[i][3] << " " << endl;
		}*/

		
		for (int i = 0; i < uvDemandPair.size(); i++) {

			vector<double> demandinfo;

			demandinfo.push_back(uvDemandPair[i][0]);
			demandinfo.push_back(uvDemandPair[i][1]);
			demandinfo.push_back(uvDemandPair[i][2]);
		    
			
			int demandIndex= findDemandIndex(uvDemandPair[i][0],uvDemandPair[i][1]);

			
			//demand is leaving
			if(uvDemandPair[i][3]==1){

				/*
				//to prevent demands in the ouoftable from leaving the topology because they didn't go into topology before
				if(Outoftable.size()!=0){
					int tag=0;
					for(int j=0;j<Outoftable.size();j++){
						//cout << Outoftable.size()<<endl;
						//cout << Outoftable[j][0] << " "<< Outoftable[j][1]<< " "<< Outoftable[j][2]<<endl;
						if(Outoftable[j][0]==uvDemandPair[i][0] && Outoftable[j][1]==uvDemandPair[i][1]){
							tag=1;
							Outoftable[j].push_back(1);
							break;
						}
					}
					if(tag==1){
						continue;
					}
				}
				*/

				demandleaving(demandIndex,demandinfo,R,flowOnLink,paths,ruleSpace,DefaultruleSpace);
				/*
				//reroute demands in the outoftable if source and distination of the demand are available,
				//becacuse some demands have been left the topology,so some source and distination of the demand may have space.
				if(Outoftable.size()!=0){
					for(int j=0;j<Outoftable.size();j++){
						
						if(ruleSpace[Outoftable[j][0]].size()<nodeCap && ruleSpace[Outoftable[j][1]].size()<nodeCap){

							demandinfo.clear();
							demandinfo.push_back(Outoftable[j][0]);
							demandinfo.push_back(Outoftable[j][1]);
							demandinfo.push_back(Outoftable[j][2]);

							demandIndex= findDemandIndex(Outoftable[j][0],Outoftable[j][1]);

							vector<vector<int> > Gp = G;
							vector<vector<double> > Lp(N, vector<double>(N, MAX));

							updatingLp(G,Lp,flowOnLink,R,ruleSpace);

							for (int l = 0; l < N; l++) {
								if (ruleSpace[l].size() >= nodeCap) {//out of space
									
									if(DefaultruleSpace[l].size()==0){
										//not shrink yet, so we shrink it
										OurShrinkTable(ruleSpace[l],DefaultruleSpace[l]);

									}
									else{
										//has been shrinked, so we need to remove it
										for (int k = 0; k < N; k++) {
											if (G[l][k] == 1) {
												Lp[l][k] = MAX;
												Lp[k][l] = MAX;
												Gp[l][k] = 0;
												Gp[k][l] = 0;
											}
										}
									}
								}
							}

							for (int k = 0; k < M; k++) {
								if (R[uvStructPairs[k][0]][uvStructPairs[k][1]] < Outoftable[j][2]) {
									//cout << "remove link " << uvStructPairs[k][0] << " " << uvStructPairs[k][1] << endl;
									//remove link
									Lp[uvStructPairs[k][0]][uvStructPairs[k][1]] = MAX;
									Lp[uvStructPairs[k][1]][uvStructPairs[k][0]] = MAX;

									Gp[uvStructPairs[k][0]][uvStructPairs[k][1]] = 0;
									Gp[uvStructPairs[k][1]][uvStructPairs[k][0]] = 0;
								}
							}

							vector<int> path = SP(Lp, Outoftable[j][0], Outoftable[j][1]);

							if (!pathIsReachable(Gp, path)) { 

								outputFile <<" demand: " <<Outoftable[j][0]<<" "<<Outoftable[j][1]<<" "<<Outoftable[j][2]<< " : ";
									
								//outputFile << endl;
								for (int k = 0; k < M; k++) {
									if (R[uvStructPairs[k][0]][uvStructPairs[k][1]] < Outoftable[j][2]) {
										outputFile << uvStructPairs[k][0] << uvStructPairs[k][1] << R[uvStructPairs[k][0]][uvStructPairs[k][1]] << endl;
									}
								}
								outputFile << "nodes out of space: ";
								for (int k = 0; k < N; k++) {
									if (ruleSpace[k].size() >= nodeCap) {//out of space
										//outputFile<<ruleSpace[k].size()<<" ";
										outputFile << k << " ";
											
									}
								}			
								return;
							}

							paths[demandIndex] = path ;
							
							demandarriving(demandIndex,demandinfo,R,flowOnLink,path,ruleSpace,DefaultruleSpace);

							if(Outoftable[j].size()==4){
								demandleaving(demandIndex,demandinfo,R,flowOnLink,paths,ruleSpace,DefaultruleSpace);
							}

							Outoftable.erase(Outoftable.begin() + j);
							j--;
						}
					}
				}
				*/
			}
			//demand is arriving
			else{
				vector<vector<int> > Gp = G;
				vector<vector<double> > Lp(N, vector<double>(N, MAX));
				
				updatingLp(G,Lp,flowOnLink,R,ruleSpace);

				bool Isoutoftable=false;

				//check rule space capacity constrain
				for (int j = 0; j < N; j++) {
					if (ruleSpace[j].size() >= nodeCap) {//out of space
						//not shrink yet, so we shrink it
						if(DefaultruleSpace[j].size()==0){
							
							OurShrinkTable(ruleSpace[j],DefaultruleSpace[j]);
						}
						//has been shrinked, so we need to remove it
						else{
							/*
							//if source and distination of demand are full, put into the outoftable to wait
							if(j==uvDemandPair[i][0] || j==uvDemandPair[i][1]){
								Isoutoftable=true;
							}
							*/
							for (int k = 0; k < N; k++) {
								if (G[j][k] == 1) {
									Lp[j][k] = MAX;
									Lp[k][j] = MAX;
									Gp[j][k] = 0;
									Gp[k][j] = 0;
								}
							}		
						}	
					}
				}

				/*
				if(Isoutoftable){
					vector<double> temp;
					temp.push_back(uvDemandPair[i][0]);
					temp.push_back(uvDemandPair[i][1]);
					temp.push_back(uvDemandPair[i][2]);
					Outoftable.push_back(temp);
					continue;
				}
				*/

				//check link capacity constrain
				for (int j = 0; j < M; j++) {
					if (R[uvStructPairs[j][0]][uvStructPairs[j][1]] < uvDemandPair[i][2]) {
						//cout << "remove link " << uvStructPairs[j][0] << " " << uvStructPairs[j][1] << endl;
						//remove link
						Lp[uvStructPairs[j][0]][uvStructPairs[j][1]] = MAX;
						Lp[uvStructPairs[j][1]][uvStructPairs[j][0]] = MAX;

						Gp[uvStructPairs[j][0]][uvStructPairs[j][1]] = 0;
						Gp[uvStructPairs[j][1]][uvStructPairs[j][0]] = 0;
					}
				}
				

				//use weight matrix Lp to find shortest path
				vector<int> path = SP(Lp, uvDemandPair[i][0], uvDemandPair[i][1]);

				//check path exist
				if (!pathIsReachable(Gp, path)) { 
					//outputFile <<" demand: " <<uvDemandPair[i][0]<<" "<<uvDemandPair[i][1]<<" "<<uvDemandPair[i][2]<< " : ";
					/*
					//outputFile << endl;
					for (int j = 0; j < M; j++) {
						if (R[uvStructPairs[j][0]][uvStructPairs[j][1]] < uvDemandPair[i][2]) {
							outputFile << uvStructPairs[j][0] <<" "<< uvStructPairs[j][1] <<" "<< R[uvStructPairs[j][0]][uvStructPairs[j][1]] << endl;
						}
					}
					outputFile << "nodes out of space: ";
					for (int j = 0; j < N; j++) {
						if (ruleSpace[j].size() >= nodeCap) {//out of space
							//outputFile<<ruleSpace[j].size()<<endl;
							outputFile << j << " ";

							
						}
					}
					*/
					demand_no_path++;
					
					//cout << "nrbl" << endl;
					//cout << "nrbl" << endl;
					//return;
				}
				else{
					paths[demandIndex] = path ;
				
					demandarriving(demandIndex,demandinfo,R,flowOnLink,path,ruleSpace,DefaultruleSpace);

				}

				

			}
			
		}
		/*
		if(Outoftable.size()!=0){
				
			//outputFile << endl;
			for (int j = 0; j < Outoftable.size(); j++) {
				outputFile <<" demand: " <<Outoftable[j][0]<<" "<<Outoftable[j][1]<<" "<<Outoftable[j][2]<< " : ";
			}
			outputFile << "TAG"<<" ";
			outputFile << "nodes out of space: ";
			for (int j = 0; j < N; j++) {
				if (ruleSpace[j].size() >= nodeCap) {//out of space
					//outputFile<<ruleSpace[j].size()<<endl;
					outputFile << j << " ";
						
				}
			}
			
			return;
		}
		*/
		

		
		int count = 0;
		for (int i = 0; i < M; i++) {
			if (flowOnLink[uvStructPairs[i][0]][uvStructPairs[i][1]] == 0) count++;
		}
		totalcount=totalcount+count;
		cout << count;
		//outputFile << count;
		
		//cout << count << endl;
		pathAll.push_back(paths);
		
		finish = clock();
		double pT = (double(finish - start) / CLOCKS_PER_SEC);
		cout << " pt: " << pT << " ";
		//outputFile << " pt: " << pT << " ";
		
		//calculate the open cost and close cost in every demand 
		for(int i=0;i<N;i++){
			for(int j=0;j<N;j++){
				//unused-->used
				if(previousFlowOnLink[i][j]==0 && flowOnLink[i][j]!=0){
					openCost++;
					totalOpenCost++;

				}
				//uesd-->unused
				else if(previousFlowOnLink[i][j]!=0 && flowOnLink[i][j]==0){
					closeCost++;
					totalCloseCost++;
				}
			}
		}
		openCost=openCost/2;
		closeCost=closeCost/2;
		

		previousFlowOnLink=flowOnLink;
		cout <<"OpenLinks: " << openCost << " ";
		//outputFile <<"OpenLinks: " << openCost << " ";
		cout <<"CloseLinks: " << closeCost << endl;
		//outputFile <<"CloseLinks: " << closeCost << endl;
	} //for different demand t

	
	totalOpenCost=totalOpenCost/2;
	totalCloseCost=totalCloseCost/2;
	cout <<"TotalOpenLinks: " << totalOpenCost << endl;
	//outputFile <<"TotalOpenLinks: " << totalOpenCost << endl;
	cout <<"TotalCloseLinks: " << totalCloseCost << endl;
	//outputFile <<"TotalCloseLinks: " << totalCloseCost << endl;

	outputFile << (totalcount/288) << " : "; 
	outputFile << (totalOpenCost + totalCloseCost)<< " : ";
	outputFile << demand_no_path << " ";
	

	

}

void spWithTime() {
	vector<vector<int> > G(N, vector<int>(N, 0)); //adjacency matrix, connect:1 notConnect:0
	vector<vector<double> > R(N, vector<double>(N, 0)); //link recidual space
	vector<vector<int> > flowOnLink(N, vector<int>(N, 0)); //#flows on link
	vector<vector<int> > paths (demandIndex.size(), vector<int>()); //flow path
	vector<vector<vector<int> > > ruleSpace(N, vector<vector<int> >());//[source,dest,port]
	vector<vector<vector<int> > > pathAll; //flow path
	vector<vector<int> > previousFlowOnLink (N, vector<int>(N, 0));
	int totalOpenCost=0;
	int totalCloseCost=0;

	int totalcount=0;
										   //construct graph
	for (int i = 0; i < M; i++) {
		G[uvStructPairs[i][0]][uvStructPairs[i][1]] = 1;
		G[uvStructPairs[i][1]][uvStructPairs[i][0]] = 1;

		R[uvStructPairs[i][0]][uvStructPairs[i][1]] = linkCap;
		R[uvStructPairs[i][1]][uvStructPairs[i][0]] = linkCap;
	}
	//algorithm begin
	for (int t = 0; t < uvDemandPairs.size(); t++) {
		clock_t start, finish;
		start = clock();
		cout << "t:" << t + 1 << ": ";
		//outputFile << "t:" << t + 1 << ": ";
		//cout << endl;
		vector<vector<double> > uvDemandPair = uvDemandPairs[t];
		//vector<bool> isInPrevious(uvDemandPair.size(), false);
		int openCost=0;
		int closeCost=0;
		
		
		
		for (int i = 0; i < uvDemandPair.size(); i++) { 
		    
			
			int demandIndex= findDemandIndex(uvDemandPair[i][0],uvDemandPair[i][1]);
			
			//demand is leaving
			if(uvDemandPair[i][3]==1){
				
				//remove rule from table for each node on path
				for (int j = 0; j < paths[demandIndex].size(); j++) {
					for (int k = 0; k < ruleSpace[paths[demandIndex][j]].size(); k++) {
						if (ruleSpace[paths[demandIndex][j]][k][0] == uvDemandPair[i][0] && ruleSpace[paths[demandIndex][j]][k][1] == uvDemandPair[i][1]) {
							
							ruleSpace[paths[demandIndex][j]].erase(ruleSpace[paths[demandIndex][j]].begin() + k);
							k--;
						}
					}
				}

				//add back residual space, reduce # flow on link for each link on path
				for (int j = 0; j < paths[demandIndex].size() - 1; j++) {
						
					R[paths[demandIndex][j]][paths[demandIndex][j + 1]] += uvDemandPair[i][2];
					R[paths[demandIndex][j + 1]][paths[demandIndex][j]] += uvDemandPair[i][2];
					flowOnLink[paths[demandIndex][j]][paths[demandIndex][j + 1]]--;
					flowOnLink[paths[demandIndex][j + 1]][paths[demandIndex][j]]--;
						
				}
				paths[demandIndex].clear();
			        
			}


			//demand is arriving
			else{
				vector<vector<int> > Gp = G;
				vector<vector<double> > Lp(N, vector<double>(N, MAX));
				for (int u = 0; u < N; u++) {
					for (int v = 0; v < N; v++) {
						if (G[u][v] == 1) Lp[u][v] = 1;//uv has link and # flow on uv is 0
					}
				}

				//check rule space capacity constrain
				for (int j = 0; j < N; j++) {
					if (ruleSpace[j].size() >= nodeCap) {
						//remove link connected to node j
						//////cout << "remove node " << j << endl;
						//cout << "remove node " << j << endl;
						
						for (int k = 0; k < N; k++) {
							if (G[j][k] == 1) {
								Lp[j][k] = MAX;
								Lp[k][j] = MAX;

								Gp[j][k] = 0;
								Gp[k][j] = 0;
							}
						}
						
						
					}
				}

				//check link capacity constrain
				for (int j = 0; j < M; j++) {
					if (R[uvStructPairs[j][0]][uvStructPairs[j][1]] < uvDemandPair[i][2]) {
						//cout << "remove link " << uvStructPairs[j][0] << " " << uvStructPairs[j][1] << endl;
						//remove link
						Lp[uvStructPairs[j][0]][uvStructPairs[j][1]] = MAX;
						Lp[uvStructPairs[j][0]][uvStructPairs[j][1]] = MAX;

						Gp[uvStructPairs[j][0]][uvStructPairs[j][1]] = 0;
						Gp[uvStructPairs[j][0]][uvStructPairs[j][1]] = 0;
					}
				}
				//for (int j = 0; j < N; j++) {
				//	for (int k = 0; k < N; k++) {
				//		////////cout << "check " << j << " " << k << ", " << R[j][k] << " " << uvDemandPair[i][2] << ", " << G[j][k] << endl;
				//		//not the same node,link capacity < demand,uv has link
				//		if (j != k && R[j][k] < uvDemandPair[i][2] && G[j][k] == 1) {
				//			//////cout << "remove link " << j << " " << k << endl;
				//			//cout << "remove link " << j << " " << k << endl;
				//			//remove link
				//			Lp[j][k] = MAX;
				//			Lp[k][j] = MAX;
				//			Gp[j][k] = 0;
				//			Gp[k][j] = 0;
				//		}
				//	}
				//}

				//use weight matrix Lp to find shortest path
				vector<int> path = SP(Lp, uvDemandPair[i][0], uvDemandPair[i][1]);

				//check path exist
				if (!pathIsReachable(Gp, path)) {
					outputFile <<" demand: " <<uvDemandPair[i][0]<<" "<<uvDemandPair[i][1]<<" "<<uvDemandPair[i][2]<< " : ";
					
					
					//outputFile << endl;
					for (int j = 0; j < M; j++) {
						if (R[uvStructPairs[j][0]][uvStructPairs[j][1]] < uvDemandPair[i][2]) {
							outputFile << uvStructPairs[j][0] << uvStructPairs[j][1] << R[uvStructPairs[j][0]][uvStructPairs[j][1]] << endl;
						}
					}
					outputFile << "nodes out of space: ";
					for (int j = 0; j < N; j++) {
						if (ruleSpace[j].size() >= nodeCap) {//out of space
							//outputFile<<ruleSpace[j].size()<<endl;
							outputFile << j << " ";

							
						}
					}
					
					//cout << "nrbl" << endl;
					//cout << "nrbl" << endl;
					return;
				}


				paths[demandIndex] = path;
				//cout << "path: ";
				//cout << "path: ";
				for (int i = 0; i < path.size(); i++) {
					//////cout << path[i] << " ";
					//cout << path[i] << " ";
				}
				//cout << endl;
				//cout << endl;

				//cout << "substract residual space, increase flowOnlink" << endl;
				//cout << "substract residual space, increase flowOnlink" << endl;
				//substract residual space, increase flowOnlink
				for (int j = 0; j < path.size() - 1; j++) {
					//cout << "before R: " << path[j] << " " << path[j + 1] << ": " << R[path[j]][path[j + 1]] << endl;
					//cout << "before R: " << path[j] << " " << path[j + 1] << ": " << R[path[j]][path[j + 1]] << endl;
					//////cout << "before flowOnLink: " << path[j] << " " << path[j + 1] << ": " << flowOnLink[path[j]][path[j + 1]] << endl;
					//cout << "before flowOnLink: " << path[j] << " " << path[j + 1] << ": " << flowOnLink[path[j]][path[j + 1]] << endl;
					////////coutFile << path[j] << " " << path[j + 1] <<endl;
					R[path[j]][path[j + 1]] -= uvDemandPair[i][2];
					flowOnLink[path[j]][path[j + 1]] ++;

					R[path[j + 1]][path[j]] -= uvDemandPair[i][2];
					flowOnLink[path[j + 1]][path[j]] ++;

					//cout << "after R: " << path[j] << " " << path[j + 1] << ": " << R[path[j]][path[j + 1]] << endl;
					//cout << "after R: " << path[j] << " " << path[j + 1] << ": " << R[path[j]][path[j + 1]] << endl;
					//cout << "after flowOnLink: " << path[j] << " " << path[j + 1] << ": " << flowOnLink[path[j]][path[j + 1]] << endl;
					//cout << "after flowOnLink: " << path[j] << " " << path[j + 1] << ": " << flowOnLink[path[j]][path[j + 1]] << endl;
				}
				//cout << "add rule for each node on path" << endl;
				//cout << "add rule for each node on path" << endl;
				//add rule for each node on path
				for (int j = 0; j < path.size() - 1; j++) {
					//cout << "before ruleSpace" << path[j] << ": " << ruleSpace[path[j]].size() << endl;
					//cout << "before ruleSpace" << path[j] << ": " << ruleSpace[path[j]].size() << endl;
					vector<int> temp;
					temp.push_back(uvDemandPair[i][0]);
					temp.push_back(uvDemandPair[i][1]);
					temp.push_back(path[j + 1]);
					ruleSpace[path[j]].push_back(temp);
					//////cout << "after ruleSpace" << path[j] << ": " << ruleSpace[path[j]].size() << endl;
					//cout << "after ruleSpace" << path[j] << ": " << ruleSpace[path[j]].size() << endl;
				}
			}
			
		}


		int count = 0;
		for (int i = 0; i < M; i++) {
			if (flowOnLink[uvStructPairs[i][0]][uvStructPairs[i][1]] == 0) count++;
		}
		cout << count;
		//outputFile << count;
		//cout << count << endl;
		pathAll.push_back(paths);
		/*for (int pp = 0; pp < paths.size(); pp++) {
		cout << "demand " << pp << ": ";
		for (int pp1 = 0; pp1 < paths[pp].size(); pp1++) {
		cout << paths[pp][pp1] << " ";
		}
		cout << endl;
		}*/
		//system("pause");
		//system("cls");
		finish = clock();
		double pT = (double(finish - start) / CLOCKS_PER_SEC);
		cout << " pt: " << pT << " ";
		//outputFile << " pt: " << pT << " ";
		
		//calculate the open cost and close cost in every demand 
		for(int i=0;i<N;i++){
			for(int j=0;j<N;j++){
				//unused-->used
				if(previousFlowOnLink[i][j]==0 && flowOnLink[i][j]!=0){
					openCost++;
					totalOpenCost++;

				}
				//uesd-->unused
				else if(previousFlowOnLink[i][j]!=0 && flowOnLink[i][j]==0){
					closeCost++;
					totalCloseCost++;
				}
			}
		}
		openCost=openCost/2;
		closeCost=closeCost/2;
		

		previousFlowOnLink=flowOnLink;
		cout <<"OpenLinks: " << openCost << " ";
		//outputFile <<"OpenLinks: " << openCost << " ";
		cout <<"CloseLinks: " << closeCost << endl;
		//outputFile <<"CloseLinks: " << closeCost << endl;
		
	} //for different demand t

	
	totalOpenCost=totalOpenCost/2;
	totalCloseCost=totalCloseCost/2;
	cout <<"TotalOpenLinks: " << totalOpenCost << endl;
	//outputFile <<"TotalOpenLinks: " << totalOpenCost << endl;
	cout <<"TotalCloseLinks: " << totalCloseCost << endl;
	//outputFile <<"TotalCloseLinks: " << totalCloseCost << endl;
	outputFile << (totalcount/288) << " : "; 
	outputFile << (totalOpenCost + totalCloseCost)<< " "; 

	for (int i = 0; i < pathAll.size(); i++) {
		//cout << "file====" << i << endl;
		//ofstream outPath;
		//string s = "data/ABILENE-META/_path" + to_string(i) + ".txt";
		//outPath.open(s);
		for (int j = 0; j < pathAll[i].size(); j++) {
			for (int k = 0; k < pathAll[i][j].size(); k++) {
				//outPath  << pathAll[i][j][k] << " ";
				//cout  << pathAll[i][j][k] << " ";
			}
			//outPath << endl;
			//cout << endl;
		}
		//outPath.close();
	}



	int totalChange = 0;
	for (int i = 0; i < nFile; i++) {
		int change = 0;

		if (i == 0) {
			//cout << "i = 0," << endl;
			for (int j = 0; j < pathAll[i].size(); j++) {
				if (pathAll[i][j].size() != 0) {
					//cout << "add: " << pathAll[i][j].size() - 1 << endl;
					change += (pathAll[i][j].size() - 1);
				}
			}
			/*
			cout << "change: " << change << endl;
			outputFile << "change: " << change << endl;
			totalChange += change;*/
			
			continue;
		}
		for (int j = 0; j < pathAll[i].size(); j++) {
			//cout << "path now: ";
			for (int k = 0; k < pathAll[i][j].size(); k++) {
				//cout << pathAll[i][j][k] << " ";
			}
			//cout << "compare to path previous: ";
			for (int k = 0; k < pathAll[i - 1][j].size(); k++) {
				//cout << pathAll[i - 1][j][k] << " ";
			}
			//cout << ":" << endl;

			//pathAll[i][j] compare to pathAll[i-1][j]
			if (pathAll[i][j].size() == 0 && pathAll[i - 1][j].size() == 0) {
				//cout << "all 0" << endl;
			}
			else if (pathAll[i][j].size() == 0 && pathAll[i - 1][j].size() != 0) {
				//cout << ": add " << pathAll[i - 1][j].size() - 1 << endl;
				change += pathAll[i - 1][j].size() - 1;
			}
			else if (pathAll[i][j].size() != 0 && pathAll[i - 1][j].size() == 0) {
				//cout << ": add " << pathAll[i][j].size() - 1 << endl;
				change += pathAll[i][j].size() - 1;
			}
			else {
				int repeat = 0;
				for (int k = 0; k < pathAll[i][j].size() - 1; k++) {
					//cout << "k: " << k << endl;
					//pathAll[i][j]'s link (k,k+1) compare to pathAll[i-1][j]
					//cout << "link" << pathAll[i][j][k] << ", " << pathAll[i][j][k + 1];
					for (int m = 0; m < pathAll[i - 1][j].size() - 1; m++) {
						//pathAll[i][j]'s link (k,k+1) compare to pathAll[i-1][j]'s link (m,m+1)
						if (pathAll[i][j][k] == pathAll[i - 1][j][m] && pathAll[i][j][k + 1] == pathAll[i - 1][j][m + 1]) {
							repeat++;
							break;
						}
					}
				}
				change += (pathAll[i][j].size() - 1) + (pathAll[i - 1][j].size() - 1) - 2 * repeat;
			}
		}
		/*
		cout << "change: " << change << endl;
		outputFile << "change: " << change << endl;*/
		totalChange += change;
	}
	/*
	cout << "total change: " << totalChange << endl;
	outputFile << "total change: " << totalChange << endl;*/
}

void G2014WithTime() {
	vector<vector<int> > paths(demandIndex.size(), vector<int>()); //flow path
	vector<vector<vector<int> > > pathAll; //flow path
	vector<vector<int> > G(N, vector<int>(N, 0)); //adjacency matrix, connect:1 notConnect:0
	vector<vector<double> > Re(N, vector<double>(N, 0)); //link recidual space
	vector<vector<int> > previousFlowOnLink(N,vector<int>(N,0));


	vector<vector<double> > WTemp(N, vector<double>(N, MAX));//weight
	vector<vector<tuple<int, int, int> > > FuTemp(N, vector<tuple<int, int, int> >(0, tuple<int, int, int>()));
	vector<vector<tuple<int, int, int> > > GuTemp(N, vector<tuple<int, int, int> >(0, tuple<int, int, int>()));
	vector<bool> isShrinkedTemp(N, false); //node is already shrinked or not

	for (int i = 0; i < M; i++) {
		G[uvStructPairs[i][0]][uvStructPairs[i][1]] = 1;
		G[uvStructPairs[i][1]][uvStructPairs[i][0]] = 1;
	}
	int removeEdge;
	int totalCloseCost=0;
	int totalOpenCost=0;

	 
	for (int t = 0; t < uvDemandPairs.size(); t++) {
		//cout << "here" ;
		
		clock_t start, finish;
		start = clock();
		cout << "t: " << t + 1 << ": ";
		outputFile << "t: " << t + 1 << ": ";
		//cout << "here" ;
		vector<vector<int> > flowOnLink(N,vector<int>(N,0));
		 
		
		int openCost=0;
		int closeCost=0;

		//check is it nessary to route again
		if (t != 0) {
			if (uvDemandPairs[t].size() == uvDemandPairs[t - 1].size()) {
				//cout << "size the same" << endl;
				//system("pause");
				//check is the set the same
				bool isTheSame = true;
				for (int i = 0; i < uvDemandPairs[t].size(); i++) {
					bool isFound = false;
					for (int j = 0; j < uvDemandPairs[t - 1].size(); j++) {
						if (uvDemandPairs[t][i][0] == uvDemandPairs[t - 1][j][0] && uvDemandPairs[t][i][1] == uvDemandPairs[t - 1][j][1]) {
							isFound = true;
							break;
						}
					}
					if (!isFound) {
						//cout << uvDemandPairs[t][i][0] << "," << uvDemandPairs[t][i][1] << "not found in previous" << endl;
						isTheSame = false;
						//system("pause");
						break;
					}
				}
				//if the set is the same, check is the link overflow
				if (isTheSame) {
					//cout << "if is the same" << endl;
					//system("pause");
					bool isExceed = false;
					for (int i = 0; i < uvDemandPairs[t].size(); i++) {
						int ind = findDemandIndex(uvDemandPairs[t][i][0], uvDemandPairs[t][i][1]);
						double addTraffic = uvDemandPairs[t][i][2] - uvDemandPairs[t - 1][i][2];

						for (int j = 0; j < paths[ind].size() - 1; j++) {
							//cout << "R " << Re[paths[ind][j]][paths[ind][j + 1]] << endl;
							//cout << addTraffic << endl;
							//cout << linkCap << endl;
							if (Re[paths[ind][j]][paths[ind][j + 1]] + addTraffic > linkCap) {
								//cout << "link exceed" << endl;
								isExceed = true;
								break;
							}
						}
						if (isExceed) {
							break;
						}
					}
					if (!isExceed) {
						cout << "=====no reroute=====";
						pathAll.push_back(paths);
						finish = clock();
						double pT = (double(finish - start) / CLOCKS_PER_SEC);
						cout << removeEdge;
						outputFile << removeEdge;
						cout << " pt: " << pT << endl;
						outputFile << " pt: " << pT << endl;
						continue;
					}
				}
				else {
					//cout << "the demand set is not the same" << endl;
					isTheSame = false;
					//system("pause");
				}
			}
			else {
				//cout << "size not the same" << endl;
			}

		}

		//cout << "reroute" << endl;
		paths = vector<vector<int> >(demandIndex.size(), vector<int>());

		vector<vector<double> > uvDemandPair = uvDemandPairs[t];
		
		int c = FindingFeasibleRouting(G, Re, uvDemandPair, paths,WTemp,FuTemp,GuTemp,isShrinkedTemp);
		
		
		//cout << c << endl; 
		vector<vector<double> > preComputeRe = Re;
		//logFile << "precompute Re" << endl;
		//cout << "precompute Re" << endl;
		for (int i = 0; i < uvStructPairs.size(); i++) {
			//logFile << i << " " << Re[uvStructPairs[i][0]][uvStructPairs[i][1]] << endl;
		}
		
		if (c == -1) {
			//logFile << "some demand has not route" << endl;
			cout << "some demand has no route" << endl;
			return;
		}
		else {
			//cout << "try to remove edges" << endl;
			//logFile << "try to remove edges" << endl;
			//try to remove edges
			
			vector<bool> edgeHasBeenChosen(M, false);
			removeEdge = 0;
			vector<vector<int> > Gp = G;
			vector<vector<int> > Pr = paths;

			//clock_t t1, t2;
			
			//t1 = clock();
			for (int i = 0; i < M; i++) {
				
				//find the index of least load link
				int minIndex = -1;
				double max = -1;
				for (int j = 0; j < M; j++) {
					if (preComputeRe[uvStructPairs[j][0]][uvStructPairs[j][1]] > max && !edgeHasBeenChosen[j]) {
						minIndex = j;
						max = preComputeRe[uvStructPairs[j][0]][uvStructPairs[j][1]];
					}
				}
				
				//can find index
				if (max != -1) {

					//cout << "min index: " << minIndex << endl;
					edgeHasBeenChosen[minIndex] = true;

					//remove link and find feasible routing again
					Gp[uvStructPairs[minIndex][0]][uvStructPairs[minIndex][1]] = 0;
					Gp[uvStructPairs[minIndex][1]][uvStructPairs[minIndex][0]] = 0;

					if(max==linkCap){
						Pr = paths;
						removeEdge++;
					}
					else{

						vector<vector<double> > uvDemandPairtemp=GetDemandpair(minIndex,paths,uvDemandPair);

						int c = NewFindingFeasibleRouting(Gp, Re, uvDemandPairtemp , paths,WTemp,FuTemp,GuTemp,isShrinkedTemp);
					
						if (c != -1) {
							//cout << "remove" << endl;
							Pr = paths;
							removeEdge++;
						}
						else {
							//cout << "cant remove" << endl;
							Gp[uvStructPairs[minIndex][0]][uvStructPairs[minIndex][1]] = 1;
							Gp[uvStructPairs[minIndex][1]][uvStructPairs[minIndex][0]] = 1;
						}

					}
				}	
			}
			//t2 = clock();
			//double spendt = (double(t2 - t1) / CLOCKS_PER_SEC);
			//cout << " Time: " << spendt;

			//logFile << removeEdge << endl;
			
			cout << removeEdge;
			outputFile << removeEdge;
			
			/*
			for (int i = 0; i < Pr.size(); i++) {
				//hopCountFile << Pr[i].size() - 1 << ": ";
				for (int j = 0; j < Pr[i].size(); j++) {
					//hopCountFile << Pr[i][j] << " ";
				}
				//hopCountFile << endl;
			}*/

		}
		//cout << t << endl; 
		pathAll.push_back(paths);
		/*ofstream outPath;
		string s = demandDirectory + to_string(t) + "g2014path.txt";
		outPath.open(s);
		for (int pp = 0; pp < paths.size(); pp++) {
		outPath << "demand " << pp << ": ";
		for (int pp1 = 0; pp1 < paths[pp].size(); pp1++) {
		outPath << paths[pp][pp1] << " ";
		}
		outPath << endl;
		}*/
		finish = clock();
		double pT = (double(finish - start) / CLOCKS_PER_SEC);
		cout << " pt: " << pT << endl;
		outputFile << " pt: " << pT << endl;
		//cout <<"here: " << endl;
		//cout <<paths.size() << endl;
		
		for(int i=0;i<paths.size();i++){
			if(paths[i].size()!=0){
				for(int j=0;j<paths[i].size()-1;j++){
					flowOnLink[paths[i][j]][paths[i][j+1]]=1;
					flowOnLink[paths[i][j+1]][paths[i][j]]=1;
				}
			}
		}
		//cout <<"here: " << endl;
		for(int i=0;i<N;i++){
			for(int j=0;j<N;j++){
				//unused-->used
				if(previousFlowOnLink[i][j]==0 && flowOnLink[i][j]!=0){
					openCost++;
					totalOpenCost++;

				}
				//uesd-->unused
				else if(previousFlowOnLink[i][j]!=0 && flowOnLink[i][j]==0){
					closeCost++;
					totalCloseCost++;
				}
			}
		}
		openCost=openCost/2;
		closeCost=closeCost/2;
		

		previousFlowOnLink=flowOnLink;
		cout <<"OpenLinks: " << openCost << " ";
		outputFile << "OpenLinks: " << openCost << " ";
		cout <<"CloseLinks: " << closeCost << endl;
		outputFile << "CloseLinks: " << closeCost << endl;
		
	}

	
	totalOpenCost=totalOpenCost/2;
	totalCloseCost=totalCloseCost/2;
	cout <<"TotalOpenLinks: " << totalOpenCost << endl;
	outputFile <<"TotalOpenLinks: " << totalOpenCost << endl;
	cout <<"TotalCloseLinks: " << totalCloseCost << endl;
	outputFile<<"TotalCloseLinks: " << totalCloseCost << endl;
	
	/*
	int totalChange = 0;
	for (int i = 0; i < nFile; i++) {
		int change = 0;
		if (i == 0) {
			//cout << "i = 0," << endl;
			for (int j = 0; j < pathAll[i].size(); j++) {
				if (pathAll[i][j].size() != 0) {
					//cout << "add: " << pathAll[i][j].size() - 1 << endl;
					change += (pathAll[i][j].size() - 1);
				}
			}
			
			cout << "change: " << change << endl;
			outputFile << "change: " << change << endl;
			totalChange += change;
			
			continue;
		}
		for (int j = 0; j < pathAll[i].size(); j++) {
			//cout << "path now: ";
			for (int k = 0; k < pathAll[i][j].size(); k++) {
				//cout << pathAll[i][j][k] << " ";
			}
			//cout << "compare to path previous: ";
			for (int k = 0; k < pathAll[i - 1][j].size(); k++) {
				//cout << pathAll[i - 1][j][k] << " ";
			}
			//cout << ":" << endl;

			//pathAll[i][j] compare to pathAll[i-1][j]
			if (pathAll[i][j].size() == 0 && pathAll[i - 1][j].size() == 0) {
				//cout << "all 0" << endl;
			}
			else if (pathAll[i][j].size() == 0 && pathAll[i - 1][j].size() != 0) {
				//cout << ": add " << pathAll[i - 1][j].size() - 1 << endl;
				change += pathAll[i - 1][j].size() - 1;
			}
			else if (pathAll[i][j].size() != 0 && pathAll[i - 1][j].size() == 0) {
				//cout << ": add " << pathAll[i][j].size() - 1 << endl;
				change += pathAll[i][j].size() - 1;
			}
			else {
				int repeat = 0;
				for (int k = 0; k < pathAll[i][j].size() - 1; k++) {
					//cout << "k: " << k << endl;
					//pathAll[i][j]'s link (k,k+1) compare to pathAll[i-1][j]
					//cout << "link" << pathAll[i][j][k] << ", " << pathAll[i][j][k + 1];
					for (int m = 0; m < pathAll[i - 1][j].size() - 1; m++) {
						//pathAll[i][j]'s link (k,k+1) compare to pathAll[i-1][j]'s link (m,m+1)
						if (pathAll[i][j][k] == pathAll[i - 1][j][m] && pathAll[i][j][k + 1] == pathAll[i - 1][j][m + 1]) {
							repeat++;
							break;
						}
					}
				}
				change += (pathAll[i][j].size() - 1) + (pathAll[i - 1][j].size() - 1) - 2 * repeat;
			}
		}
		
		cout << "change: " << change << endl;
		outputFile << "change: " << change << endl;
		totalChange += change;
	}
	
	cout << "total change: " << totalChange << endl;
	outputFile << "total change: " << totalChange << endl;
	*/

}

int main()
{
	srand(time(NULL));
	//initial parameter
	gama;
	beta=1;
	link_p=0.5;
	unusedw;
	resultdir="tradeoff_big/demand/";

	
	for(gama=6;gama<=7;gama++){
		string s = resultdir + "gama"+ to_string(gama) +".txt";
		outputFile.open(s);

		
		for(int i=11;i<=15;i++){

			outputFile << "demandset_ID : " << to_string(i) << " : ";
			{
				dataDirectory = "tradeoff_big/demand/"+to_string(i)+"/";
				demandDirectory = dataDirectory;
				nFile = 288;
				uvStructPairs = loadStructure();
				loadDemand();
			}

			getDemandIndex();
			nodeCap = 10240;
			linkCap = 50;
			
			
			//dynamic energy aware routing
			//string s = dataDirectory + "ourcostnode250.txt";
			//outputFile.open(s);
			clock_t start, finish;
			start = clock();
			

			oursWithTime();
			//spWithTime();
			//G2014WithTime();
			finish = clock();
			double pT = (double(finish - start) / CLOCKS_PER_SEC);
			outputFile << "process time: " << pT << endl;
			cout << "process time: " << pT << endl;
			//outputFile.close();

			dataDirectory="";
			demandDirectory="";
			uvStructPairs.clear();
			uvDemandPairs.clear();
			uvDemandPairsS.clear();
			demandIndex.clear();

		}

		outputFile.close();
	}


	system("PAUSE");
	return 0;
}
