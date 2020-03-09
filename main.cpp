#include<iostream>
#include<vector>
#include<cmath>
#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include <fstream>

using namespace std;

#define ANT_NUM 32
#define MAX_GEN 100
//#define V_NUM 48
#define INF 0x3f3f3f
#define alpha 1
#define beta 5
#define Q 100
#define rho 0.6
#define c 0.00001
#define random(x) (rand()%x)
int GEN = 1;

struct POS {
	double x;
	double y;
};

struct Next_POS {
	POS n_pos;
	double n_distance;
	double probability;
};

struct Ants {
	POS current_pos;
	vector<Next_POS> next_pos;
	POS best_next_pos;
	vector<POS> shortest_path;
	double len_path = 0;
	vector<POS> visited_table;
};

void Init();
bool operator == (const POS pos_a, const POS pos_b);
bool IsVertexVisited(Ants ant, POS pos);
int GetVertexIndex(POS pos);
void GetNextPos(Ants& ant);
void Move(Ants& ant);
void UpdatePheromon();
void Clear();
void PrintSolution();

vector<Ants> ant;
vector<POS> vertex;
vector<vector<double>> _distance;
vector<vector<double>> pheromon;

double GetDistance(POS a, POS b) {
	return sqrt(pow(a.x - b.x, 2) + pow(a.y - b.y, 2));
}

int main() {
	srand((int)time(0));

	Init();
	while (GEN <= MAX_GEN) {
		for (unsigned int i = 0; i < vertex.size() - 1; ++i) {
			for (int j = 0; j < ANT_NUM; ++j) {
				//calculate next position
				GetNextPos(ant[j]);
				//move to next postion, update pheromon
				Move(ant[j]);
			}
		}
		//update pheromon data
		UpdatePheromon();	//output pheromon matrix to check
		cout << GEN << ": ";
		PrintSolution();

		/*double min_len = INF;
		int flag = 0;
		for (int i = 0; i < ANT_NUM; ++i) {
			if (ant[i].len_path < min_len) {
				min_len = ant[i].len_path;
				flag = i;
			}
		}
		cout << ant[flag].len_path << endl;*/


		if (GEN != MAX_GEN)
			Clear();
		GEN += 1;
	}

	string filename = "result.txt";
	fstream fout;
	fout.open(filename.c_str(), ios::out);
	double min_len = INF;
	int flag = 0;
	for (int i = 0; i < ANT_NUM; ++i) {
		if (ant[i].len_path < min_len) {
			min_len = ant[i].len_path;
			flag = i;
		}
	}
	cout << ant[flag].len_path << endl;
	for (auto i : ant[flag].shortest_path) {
		fout << i.x << " " << i.y << endl;
	}
	fout.close();
}

void Init() {
	//fill the vertex
	/*int pos_x, pos_y;
	while (true) {
		cin >> pos_x >> pos_y;
		POS temp_v;
		temp_v.x = pos_x;
		temp_v.y = pos_y;
		vertex.push_back(temp_v);
	}*/

	string filename = "data.txt";		//input txt file
	fstream fin;
	fin.open(filename.c_str(), ios::in);
	int pos_x, pos_y;
	while (fin >> pos_x >> pos_y) {
		POS temp_v;
		temp_v.x = pos_x;
		temp_v.y = pos_y;
		vertex.push_back(temp_v);
	}
	fin.close();

	//construct the distance & pheromon matrix
	unsigned int size = vertex.size();
	for (unsigned int i = 0; i < size; ++i) {
		vector<double> row;
		for (unsigned int j = 0; j < size; ++j) {
			if (i == j) row.push_back(INF);
			else row.push_back(GetDistance(vertex[i], vertex[j]));
		}
		_distance.push_back(row);
	}

	for (unsigned int i = 0; i < size; ++i) {
		vector<double> row;
		for (unsigned int j = 0; j < size; ++j) {
			row.push_back(c);
		}
		pheromon.push_back(row);
	}

	//init ants
	for (int i = 0; i < ANT_NUM; ++i) {
		Ants temp_ant;
		int rand_pos = random(vertex.size());
		temp_ant.current_pos = vertex[rand_pos];
		temp_ant.shortest_path.push_back(temp_ant.current_pos);
		temp_ant.visited_table.push_back(temp_ant.current_pos);
		ant.push_back(temp_ant);
	}
}

bool operator == (const POS pos_a, const POS pos_b) {
	if (pos_a.x == pos_b.x && pos_a.y == pos_b.y) return true;
	else return false;
}

bool IsVertexVisited(Ants ant, POS pos) {
	for (auto j : ant.visited_table) {
		if (pos == j) return true;
	}
	return false;
}

int GetVertexIndex(POS pos) {
	for (unsigned int i = 0; i < vertex.size(); ++i) {
		if (pos == vertex[i]) return i;
	}
	return -1;
}

void GetNextPos(Ants &ant) {
	//construct next position table according to visited table
	for (auto i : vertex) {
		if (!IsVertexVisited(ant, i)) {
			Next_POS temp_next_pos;
			temp_next_pos.n_pos = i;
			ant.next_pos.push_back(temp_next_pos);
		}
	}

	//get next position's distance
	for (auto& i : ant.next_pos) {
		i.n_distance = GetDistance(ant.current_pos, i.n_pos);
	}

	//calculate probability of each next position, get the best next postion, update shortest path & visited table
	for (auto& i : ant.next_pos) {
		double ita = 1 / i.n_distance;
		int row = GetVertexIndex(ant.current_pos);
		int column = GetVertexIndex(i.n_pos);
		double tao = pheromon[row][column];
		i.probability = pow(tao, alpha) * pow(ita, beta);
	}
	double sum = 0.0;
	for (auto i : ant.next_pos) {
		sum += i.probability;
	}
	for (auto i : ant.next_pos) {
		//if (sum == 0.0) i.probability = 0.0;
		//else 
		i.probability /= sum;
	}
}

void Move(Ants &ant) {
	//select best next position
	double max = ant.next_pos[0].probability;
	int max_index = 0;
	for (unsigned int i = 0; i < ant.next_pos.size(); ++i) {
		if (ant.next_pos[i].probability > max) {
			max = ant.next_pos[i].probability;
			max_index = i;
		}
	}
	//if (max == 0.0) {
	//	max_index = random(ant.next_pos.size());
	//	ant.best_next_pos = ant.next_pos[max_index].n_pos;
	//}
	//else {
		ant.best_next_pos = ant.next_pos[max_index].n_pos;
	//}

	//move, update postion
	ant.current_pos = ant.best_next_pos;
	ant.shortest_path.push_back(ant.best_next_pos);
	ant.len_path += ant.next_pos[max_index].n_distance;
	ant.visited_table.push_back(ant.best_next_pos);
	ant.next_pos.clear();
}

void UpdatePheromon() {
	unsigned int size = vertex.size();
	for (unsigned int i = 0; i < size; ++i) {
		for (unsigned int j = 0; j < size; ++j) {
			if (i != j) pheromon[i][j] *= rho;
		}
	}
	for (auto i : ant) {
		for (unsigned int j = 0; j < size - 1; ++j) {
			int row = GetVertexIndex(i.shortest_path[j]);
			int column = GetVertexIndex(i.shortest_path[j + 1]);
			pheromon[row][column] += Q / i.len_path;
			pheromon[column][row] = pheromon[row][column];
		}
	}

	/*for (unsigned int i = 0; i < size; ++i) {
		for (unsigned int j = 0; j < size; ++j) {
			cout << pheromon[i][j] << " ";
		}
		cout << endl;
	}*/
}

void Clear() {
	for (auto& i : ant) {
		i.len_path = 0;
		i.next_pos.clear();
		i.shortest_path.clear();
		i.visited_table.clear();
		i.shortest_path.push_back(i.current_pos);
		i.visited_table.push_back(i.current_pos);
	}
}

void PrintSolution() {
	for (auto i : ant) {
		cout << i.len_path << " ";
		/*for (auto j : i.shortest_path) {
			cout << j.x << " " << j.y << endl;
		}*/
	}
	cout << endl;
}