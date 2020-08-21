#include <cstdio>
#include <queue>
#include <utility>
#include <string>
#include <tuple>
#include <ctime>
#include <cstdlib>

using std::priority_queue;
using std::pair;
using std::tuple;

const int maxn = 30;
const int maxm = 30;

FILE *fin;
int n;   //number of ions
double cost[maxn][maxn];
double avg_cost[maxn];
double min_cost = 1.0;
int ions[maxn];

int m;  //number of qubits
int p;  //number of edges
int ngates[maxm][maxm];
int ngates_total = 0;
int degree[maxm];
int qubits[maxm];

double time_limit;
double time_start;

struct node {
	int v;
	node *left, *right;
};

node *head;
node nodes[maxn];

int mapping[maxm];
int optimal_mapping[maxm];
double optimal_cost = 1e9;
int counter = 0;
double heur_mat[maxm][maxn];

tuple<double, int, int>  edges[maxn*maxn];

struct cmp {
	int i;
	cmp(int _i): i(_i) {}
	bool operator() (int u, int v) { return heur_mat[i][u] < heur_mat[i][v]; }
};

struct cmp2 {
	bool operator() (int u, int v) { return avg_cost[u] < avg_cost[v]; }
};

void build_ion_list() {
	for (int i = 0; i < n; i++) {
		avg_cost[i] = 0;
		for (int j = 0; j < n; j++) {
			if (i == j)  continue;
			cost[i][j] -= min_cost;
			avg_cost[i] += cost[i][j];
		}
		avg_cost[i] /= n-1;
	}

	head = new node;
	head->left = head->right = head;
	for (int i = 0; i < n; i++)  ions[i] = i;

	std::sort(ions, ions+n, cmp2());
	for (int i = 0; i < n; i++) {
		node *v = nodes + ions[i];
		v->v = ions[i];
		node *pt = head->left;
		v->right = head;
		v->left = pt;
		pt->right = head->left = v;
	}
}


void reorder_qubits() {
	int start = 0;
	for (int i = 1; i < m; i++) {
		if (degree[i] > degree[start])  start = i;
	}

	int tail = 0;
	bool visited[maxn];
	for (int i = 0; i < m; i++)  visited[i] = false;
	visited[start] = true;
	priority_queue< pair<int, int> >  H;
	H.push( std::make_pair(degree[start], start) );
	while (!H.empty()) {
		int u = H.top().second;
		H.pop();
		qubits[tail++] = u;
		for (int v = 0; v < m; v++) {
			if (ngates[u][v] > 0 && !visited[v]) {
				visited[v] = true;
				H.push(std::make_pair(degree[v], v));
			}
		}
	}
}

void search(int depth, double cost_temp) {
	double current_time = clock();

	if ((current_time - time_start) / CLOCKS_PER_SEC > time_limit) {
		printf("%.6lf\n", optimal_cost + min_cost * ngates_total);
		for (int i = 0; i < m; i++)  printf("%d ", optimal_mapping[i]);
		puts("");
		exit(0);
	}
	double heuristic = 0.0;

	//the lower bound of weights between {qubits[0], ..., qubits[depth-1]} and {qubits[depth], ..., qubits[m-1]}
	for (int i = depth; i < m; i++) {
		int &qi = qubits[i];
		double min = 1e9;
		for (node *pt = head->right; pt != head; pt = pt->right) {
			int &v = pt->v;
			if (min > heur_mat[qi][v])  min = heur_mat[qi][v];
		}
		heuristic += min;
	}

	if (cost_temp + heuristic >= optimal_cost)  return ;

	++counter;
	if (depth == m) {
		optimal_cost = cost_temp;
		memcpy(optimal_mapping, mapping, sizeof(int) * m);
		return ;
	}


	int &q = qubits[depth];


	for (node *pt = head->right; pt != head; pt = pt->right) {
		int &v = pt->v;

		mapping[q] = v;
		double increment = 0.0;

		pt->left->right = pt->right;
		pt->right->left = pt->left;

		for (int i = depth+1; i < m; i++) {
			int &qi = qubits[i];
			for (node *pt2 = head->right; pt2 != head; pt2 = pt2->right) {
				heur_mat[qi][pt2->v] += ngates[q][qi] * cost[v][pt2->v];
			}	
		}
		search(depth+1, cost_temp + heur_mat[q][v]);

		for (int i = depth+1; i < m; i++) {
			int &qi = qubits[i];
			for (node *pt2 = head->right; pt2 != head; pt2 = pt2->right) {
				heur_mat[qi][pt2->v] -= ngates[q][qi] * cost[v][pt2->v];
			}	
		}
		pt->left->right = pt->right->left = pt;
	}
}

int main(int argc, char *argv[]) {
	fin = fopen(argv[1], "r");
	sscanf(argv[2], "%lf", &time_limit);
	fscanf(fin, "%d", &n);
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			fscanf(fin, "%lf", &cost[i][j]);
			if (i == j)  continue;
			if (min_cost > cost[i][j])  min_cost = cost[i][j];
		}
	}


	fscanf(fin, "%d%d", &m, &p);
	for (int i = 0; i < p; i++) {
		int u, v, w;
		fscanf(fin, "%d%d%d", &u, &v, &w);
		ngates[u][v] = ngates[v][u] = w;
		degree[u] ++;
		degree[v] ++;
		ngates_total += w;
	}

	reorder_qubits();
	build_ion_list();

	time_start = clock();
	search(0, 0.0);

	printf("%.6lf\n", optimal_cost + min_cost * ngates_total);
	for (int i = 0; i < m; i++)  printf("%d ", optimal_mapping[i]);
	puts("");

	return 0;
}
