#include <iostream>
#include <vector>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <fstream>

using namespace std;

/* global paremeter definition */
#define ANT_NUM 50			//ants number
#define MAX_GEN 500			//iteration number
#define alpha 1				//exponent of pheromone
#define beta 2				//exponent of visibility
#define rho_global 0.1		//global pheromone
#define rho_local 0.1		//local pheromone
#define qzero 0.9			//judgment condition
#define INF 0x3f3f3f		//infinity
int GEN = 0;				//current generation
double Lnn;					//length of nearest neighbor
unsigned int vertex_num;	//problem vertex number
string inputfilename = "datafile/data51.txt";
string outputfilename = "datafile/result.txt";

/* objective struct definition */
struct POS {
	int index;	//index of vertex array
	double x;
	double y;
};

struct Ants {
	POS pos;
	vector<int> tourpath;
	double path_len;
	vector<int> visited_table;
};

/* function statement */
void Init();
void ReInitInfo();
double GetNearestNeighborPath(int);
void ConstructRouters();
void Search(Ants&);
double CalTransitionProb(int, int);
void UpdateLocalPathpheromone(int, int);
void UpdateGlobalPathpheromone();
void PrintSolution();
double GetDistance(POS, POS);
double power(double, int);

/* global parameter statement */
vector<Ants> ant;						//ants array
vector<POS> vertex;						//vertex array
vector<vector<double>> Distance;		//vertex adjacency matrix
vector<vector<double>> Pheromone;		//pheromone matrix
vector<int> global_shortest_path;
double global_shortest_path_len = INF;



int main() {
	srand((int)time(NULL));
	Init();
	while (GEN++ < MAX_GEN) {
		ReInitInfo();
		ConstructRouters();
		UpdateGlobalPathpheromone();
		cout << GEN << ": " << global_shortest_path_len << endl;
	}
	PrintSolution();
}

/* initialize global parameters */
void Init() {
	//input data, construct vertex array
	fstream fin;
	fin.open(inputfilename.c_str(), ios::in);
	int index;
	double pos_x, pos_y;
	while (fin >> index >> pos_x >> pos_y) {
		POS temp_v;
		temp_v.index = index - 1;
		temp_v.x = pos_x;
		temp_v.y = pos_y;
		vertex.push_back(temp_v);
	}
	fin.close();

	//construct the distance & pheromone matrix
	vertex_num = vertex.size();
	for (unsigned int i = 0; i < vertex_num; ++i) {
		vector<double> row;
		for (unsigned int j = 0; j < vertex_num; ++j) {
			if (i == j) row.push_back(INF);
			else row.push_back(GetDistance(vertex[i], vertex[j]));
		}
		Distance.push_back(row);
	}

	Lnn = GetNearestNeighborPath(rand() % vertex_num);
	for (unsigned int i = 0; i < vertex_num; ++i) {
		vector<double> row;
		for (unsigned int j = 0; j < vertex_num; ++j) {
			row.push_back(1 / (vertex_num * Lnn));
		}
		Pheromone.push_back(row);
	}

	//init ants
	for (int i = 0; i < ANT_NUM; ++i) {
		Ants temp_ant;
		temp_ant.pos = vertex[i % vertex_num];
		ant.push_back(temp_ant);
	}
}

/* re-initialize ants infomation */
void ReInitInfo() {
	for (int i = 0; i < ANT_NUM; ++i) {
		ant[i].path_len = 0.0;
		ant[i].tourpath.clear();
		ant[i].tourpath.push_back(ant[i].pos.index);
		ant[i].visited_table.clear();
		ant[i].visited_table.resize(vertex_num);
		ant[i].visited_table[ant[i].pos.index] = 1;
	}
}

/* construct a tourpath with greedy algorithm from start vertex(ver_idx) */
double GetNearestNeighborPath(int ver_idx) {
	double sum = 0.0;

	vector<int> visited_vertex(vertex_num);
	for (unsigned int i = 0; i < vertex_num; ++i) {
		visited_vertex[i] = 0;
	}
	visited_vertex[ver_idx] = 1;

	int current = ver_idx;
	int next = -1;

	for (unsigned int i = 0; i < vertex_num - 1; ++i) {
		double min_dis = DBL_MAX;
		for (unsigned int j = 0; j < vertex_num; ++j) {
			if (visited_vertex[j] == 0) {
				if (Distance[current][j] < min_dis) {
					min_dis = Distance[current][j];
					next = j;
				}
			}
		}
		sum += Distance[current][next];
		current = next;
		visited_vertex[current] = 1;
	}
	sum += Distance[current][ver_idx];
	return sum;
}

/* construct tourpath & update global best tourpath */
void ConstructRouters() {
	vector<int> local_shortest_path;
	double local_shortest_path_len = INF;

	for (int i = 0; i < ANT_NUM; ++i) {
		Search(ant[i]);

		if (ant[i].path_len < local_shortest_path_len) {
			local_shortest_path_len = ant[i].path_len;
			local_shortest_path = ant[i].tourpath;
		}
	}

	if (local_shortest_path_len < global_shortest_path_len) {
		global_shortest_path_len = local_shortest_path_len;
		global_shortest_path = local_shortest_path;
	}
}

/* construct a complete tourpath */
void Search(Ants &ant) {
	int next;
	do {
		next = -1;
		double q = rand() / (double)RAND_MAX;

		if (q <= qzero) {
			double probability = DBL_MIN;
			for (unsigned int i = 0; i < vertex_num; i++) {
				if (ant.visited_table[i] != 1) {
					double prob = CalTransitionProb(ant.pos.index, vertex[i].index);
					if (prob > probability) {		
						probability = prob;
						next = i;
					}
				}
			}
		}
		else
		{
			double p = rand() / (double)RAND_MAX;
			double sum = 0.0;
			double probability = 0.0;
			for (unsigned int i = 0; i < vertex_num; ++i) {
				if (ant.visited_table[i] != 1) {
					sum += CalTransitionProb(ant.pos.index, vertex[i].index);
				}
			}
			for (unsigned int i = 0; i < vertex_num; ++i) {
				if (ant.visited_table[i] != 1 && sum > 0.0) {
					probability += CalTransitionProb(ant.pos.index, vertex[i].index) / sum;
					if (probability > p || (p > 0.9999 && probability > 0.9999)) {
						next = i;
						break;
					}
				}
			}
		}
		if (next >= 0) {
			UpdateLocalPathpheromone(ant.pos.index, next);
			ant.visited_table[next] = 1;
			ant.tourpath.push_back(next);
			ant.path_len += Distance[ant.pos.index][next];
			ant.pos = vertex[next];
		}
	} while (next >= 0);

	next = ant.tourpath[0];		//if next == -1, means all the vertices are visited, then comeback to start vertex

	UpdateLocalPathpheromone(ant.pos.index, next);
	ant.tourpath.push_back(next);
	ant.path_len += Distance[ant.pos.index][next];
	ant.pos = vertex[next];

}

/* update the pheromone of local path after moving a step */
void UpdateLocalPathpheromone(int row, int column) {
	Pheromone[row][column] = (1.0 - rho_local) * Pheromone[row][column] + rho_local * (1.0 / (vertex_num * Lnn));
	Pheromone[column][row] = Pheromone[row][column];
}

/* update the pheromone of global shortest path */
void UpdateGlobalPathpheromone() {
	unsigned int global_tour_size = global_shortest_path.size();
	for (unsigned int i = 1; i < global_tour_size; i++) {
		int row = global_shortest_path[i - 1];
		int column = global_shortest_path[i];
		Pheromone[row][column] = (1.0 - rho_global) * Pheromone[row][column] + rho_global * (1.0 / global_shortest_path_len);
		Pheromone[column][row] = Pheromone[row][column];
	}
}

/* calculate transition probability */
double CalTransitionProb(int row, int column) {
	double tao = Pheromone[row][column];
	double ita = 1.0 / Distance[row][column];
	return power(tao, alpha) * power(ita, beta);
}

/* print solution into outputfile */
void PrintSolution() {
	fstream fout;
	fout.open(outputfilename.c_str(), ios::out);
	for (auto i : global_shortest_path) {
		fout << vertex[i].x << " " << vertex[i].y << endl;
	}
	fout.close();
}

/* calculate Euclidean distance */
double GetDistance(POS a, POS b) {
	return sqrt(power(a.x - b.x, 2) + power(a.y - b.y, 2));
}

/* quick power opration */
double power(double x, int y) {
	double ans = 1;
	while (y) {
		if (y & 1) ans *= x;
		x *= x;
		y >>= 1;
	}
	return ans;
}