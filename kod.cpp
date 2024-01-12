/*
Author: Péter Gyimesi

To compile:
g++ kod.cpp -o srlg-path -std=c++11
g++ kod.cpp -o srlg-path2 -std=c++11 -D TWOPATH
g++ kod.cpp -o srlg-path3 -std=c++11 -D BINSEARCHCUT -D FAST
g++ kod.cpp -o srlg-path4 -std=c++11 -D GREEDY -D FAST
*/

#include <cstdio>
#include <cstring>
#include <cmath>
#include <iostream>
#include <vector>
#include <algorithm>
#include <queue>
#include <map>
#include <set>
#include <chrono>

using namespace std;

const int Max_Node_Count=1000005, Max_Edge_Count=4*Max_Node_Count, Max_Dual_Node_Count=Max_Edge_Count;

int Node_Count, X_Coord[Max_Node_Count], Y_Coord[Max_Node_Count], Node_Id[Max_Node_Count], Starting_Node, Finishing_Node;
map<int, int> Node_Id_Inverse;

int Edge_Count, Primal_Edge_Left_Endpoint[Max_Edge_Count], Primal_Edge_Right_Endpoint[Max_Edge_Count], Primal_Edge_Id[Max_Edge_Count];
int Next_Node_Left_Fixed[Max_Edge_Count], Next_Node_Right_Fixed[Max_Edge_Count];
int Prev_Node_Left_Fixed[Max_Edge_Count], Prev_Node_Right_Fixed[Max_Edge_Count]; // hiányzik, order neighboursnál kéne
map<int, int> Edge_Id_Inverse;
map<pair<int, int>, int> Edge_Inverse;

int Primal_Edge_Left_Side[Max_Edge_Count], Primal_Edge_Right_Side[Max_Edge_Count], Primal_Area_Count, Primal_Edge_Weight[Max_Edge_Count];
bool is_Special_Edge[Max_Edge_Count];

vector<int> Primal_Adj_List[Max_Dual_Node_Count];
vector<int> SRLG_list[Max_Edge_Count];
vector<int> Special_st_Path;
int max_path;
int max_srlg_size;

#ifdef FAST
bool check_input=0;
bool check_output=0;
bool eliminate_cutting_srlg=1;
#else
bool check_input=1;
bool check_output=1;
bool eliminate_cutting_srlg=1;
#endif
bool Edge_Disjoint_Paths=1;

struct weighted_edge {
    int dest, weight, edge_id;
};


int Dual_Node_Count;
vector<weighted_edge> Dual_Adj_List[Max_Dual_Node_Count];
int Min_Cut_Approx;
int SRLG_Count;

int Number_Of_SRLG_Disjoint_Paths;
vector<vector<int> > SRLG_Disjoint_Paths;
vector<int> Min_Ratio_Cut_Dual_Node, Min_Ratio_Cut_SRLG_Id;
vector<int> Min_Cut_SRLG;

vector<vector<int> > Greedy_Paths;

void print_vector(vector<int> v) {
    for (auto x:v) {
        cout << x << " ";
    }
    cout << "\n";
}

void check_the_answer() {
    set<int> edges;
    vector<int> Path_id(SRLG_Count+1, 0);
    for (int i=0; i<Number_Of_SRLG_Disjoint_Paths; i++) {
        int si=SRLG_Disjoint_Paths[i].size();
        for (int j=1; j<si; j++) {
            int a=SRLG_Disjoint_Paths[i][j-1], b=SRLG_Disjoint_Paths[i][j];
            int edge_id=Edge_Inverse[{a, b}];
            if (Edge_Disjoint_Paths) {
                if (edges.count(edge_id)) {
                    cout << "Wrong answer: not edge disjoint paths\n";
                    exit(0);
                }
            }
            edges.insert(edge_id);
            for (auto x:SRLG_list[edge_id]) {
                if (Path_id[x] && Path_id[x]!=i+1) {
                    cout << "Wrong answer: path " << Path_id[x] << " and path " << i+1 << " are not SRLG disjoint\n";
                    exit(0);
                }
                Path_id[x]=i+1;
            }
        }
    }
    cout << "<valid>1</valid>\n";
}
void print_the_answer() {
    cout << "<minCutVec>";
    for (auto x:Min_Cut_SRLG) {
        cout << x << " ";
    }
    cout << "</minCutVec>\n";

    Number_Of_SRLG_Disjoint_Paths=SRLG_Disjoint_Paths.size();
    //cout << "Number of SRLG disjoint paths: " << Number_Of_SRLG_Disjoint_Paths << "\n";
    cout << "<pathnum>" << Number_Of_SRLG_Disjoint_Paths << "</pathnum>\n";
    //cout << "Paths\n";
    for (int i=0; i<Number_Of_SRLG_Disjoint_Paths; i++) {
        cout << "<path>\n<pathnodes>";
        int length=0;
        for (auto x:SRLG_Disjoint_Paths[i]) {
            cout << Node_Id[x] << " ";
            length++;
        }
        cout << "</pathnodes>\n<pathlength>"<<length<<"</pathlength></path>\n";
    }
    if (check_output) {
        check_the_answer();
    }
    //exit(0);
}

void check_nodes() {
    for (int i=1; i<=Node_Count; i++) {
        for (int j=i+1; j<=Node_Count; j++) {
            if (Node_Id[i]==Node_Id[j]) {
                cout << "Wrong input: same id for node " << i << " and node " << j << "\n";
                exit(0);
            }
            if (X_Coord[i]==X_Coord[j] && Y_Coord[i]==Y_Coord[j]) {
                cout << "Wrong input: same coordinates for node " << i << " and node " << j << "\n";
                exit(0);
            }
        }
    }
}


void read_nodes() {
    cin >> Node_Count;
    for (int i=1; i<=Node_Count; i++) {
        cin >> X_Coord[i] >> Y_Coord[i] >> Node_Id[i];
        Node_Id_Inverse[Node_Id[i]]=i;
    }
    cin >> Starting_Node >> Finishing_Node;
    Starting_Node=Node_Id_Inverse[Starting_Node], Finishing_Node=Node_Id_Inverse[Finishing_Node];
    //cout<<endl<< Starting_Node <<" and "<< Finishing_Node;
    assert(1<=Starting_Node && Starting_Node<=Node_Count && 1<=Finishing_Node && Finishing_Node<=Node_Count);
    if (check_input) {
        check_nodes();
    }
}


int direction(int a, int b, int c) {
    assert(1<=min({a, b, c}) && max({a, b, c})<=Node_Count);
    assert(a!=b && a!=c && b!=c);

    long long val=1ll*(X_Coord[b]-X_Coord[a])*(Y_Coord[c]-Y_Coord[a])-1ll*(X_Coord[c]-X_Coord[a])*(Y_Coord[b]-Y_Coord[a]);
    if (val>0) return 1;
    if (val==0) return 0;
    if (val<0) return -1;
    return -1;
}

bool is_between(int a, int b, int c) {
    assert(direction(a, b, c)==0);
    if ((X_Coord[a]<X_Coord[b])+(X_Coord[a]<X_Coord[c])==1) return true;
    if ((Y_Coord[a]<Y_Coord[b])+(Y_Coord[a]<Y_Coord[c])==1) return true;
    return false;
}
void check_planarity() {
    for (int i=1; i<=Edge_Count; i++) {
        for (int j=i+1; j<=Edge_Count; j++) {
            int a=Primal_Edge_Left_Endpoint[i], b=Primal_Edge_Right_Endpoint[i];
            int c=Primal_Edge_Left_Endpoint[j], d=Primal_Edge_Right_Endpoint[j];
            int Equal_Count=(a==c)+(a==d)+(b==c)+(b==d);
            assert(Equal_Count<=2);
            if (Equal_Count==2) {
                cout << "Wrong input: intersection between edge " << i << " and edge " << j << "\n";
                exit(0);
            }

            if (Equal_Count==1) {
                if (b==c || b==d) swap(a, b);
                if (a==d) swap(c, d);

                assert(a==c);

                if (direction(a, b, d)==0 && !is_between(a, b, d)) {
                    cout << "Wrong input: intersection between edge " << i << " and edge " << j << "\n";
                    exit(0);
                }
            }
            if (Equal_Count==0) {
                int dir_abc=direction(a, b, c), dir_abd=direction(a, b, d);

                if (dir_abc==0 && dir_abd==0) {
                    assert(direction(c, d, a)==0 && direction(c, d, b)==0);
                    if ((is_between(a, c, d) || is_between(b, c, d)) && (is_between(c, a, b) || is_between(d, a, b))) {
                        cout << "Wrong input: intersection between edge " << i << " and edge " << j << "\n";
                        exit(0);
                    }
                } else if (dir_abc+dir_abd==0) {
                    int dir_cda=direction(c, d, a), dir_cdb=direction(c, d, b);
                    if (dir_cda+dir_cdb==0) {
                        cout << "Wrong input: intersection between edge " << i << " and edge " << j << "\n";
                        exit(0);
                    }
                }
            }
        }
    }
}


void read_edges() {
    cin >> Edge_Count;
    for (int i=1; i<=Edge_Count; i++) {
        int a, b, c;
        cin >> a >> b >> c;
        if (a==b) {
            cout << "Wrong input: edge " << i << " " << "must have different endpoints\n";
            exit(0);
        }
        a=Node_Id_Inverse[a], b=Node_Id_Inverse[b];
        if (!a || !b) {
            cout << "Wrong input: there is a problem with edge " << i << "\n";
            exit(0);
        }

        Primal_Edge_Left_Endpoint[i]=a, Primal_Edge_Right_Endpoint[i]=b;
        Primal_Adj_List[a].push_back(b), Primal_Adj_List[b].push_back(a);
        Primal_Edge_Id[i]=c, Edge_Id_Inverse[c]=i;
        Edge_Inverse[{a, b}]=i, Edge_Inverse[{b, a}]=i;
    }
    if (check_input) {
        check_planarity();
    }
}

void check_primal_connectivity() {
    vector<bool> Is_Visited(Node_Count+1, 0);
    queue<int> q;
    q.push(1);
    Is_Visited[1]=1;
    while (q.size()>0) {
        int Cur_Node=q.front();
        q.pop();
        for (auto Next_Node:Primal_Adj_List[Cur_Node]) {
            if (!Is_Visited[Next_Node]) {
                q.push(Next_Node);
                Is_Visited[Next_Node]=1;
            }
        }
    }

    for (int i=1; i<=Node_Count; i++) {
        if (!Is_Visited[i]) {
            cout << "Wrong input: the primal graph is not connected\n";
            exit(0);
        }
    }
}


void order_neighbours() {
    for (int i=1; i<=Node_Count; i++) {
        vector<pair<long double, int> > New_Order;
        for (auto adj:Primal_Adj_List[i]) {
            New_Order.push_back({atan2(X_Coord[adj]-X_Coord[i], Y_Coord[adj]-Y_Coord[i]), adj});
        }
        sort(New_Order.begin(), New_Order.end());

        Primal_Adj_List[i].clear();
        for (auto adj:New_Order) {
            Primal_Adj_List[i].push_back(adj.second);
        }

        int si=Primal_Adj_List[i].size();

        for (int j=0; j<si; j++) {
            int a=i, b=Primal_Adj_List[i][j], c=Primal_Adj_List[i][(j+1)%si];
            int edge_id=Edge_Inverse[{a, b}];
            if (Primal_Edge_Left_Endpoint[edge_id]==a) {
                Next_Node_Left_Fixed[edge_id]=c;
            } else {
                assert(Primal_Edge_Right_Endpoint[edge_id]==a);
                Next_Node_Right_Fixed[edge_id]=c;
            }
            edge_id=Edge_Inverse[{a, c}];
            if (Primal_Edge_Left_Endpoint[edge_id]==a) {
                Prev_Node_Left_Fixed[edge_id]=b;
            } else {
                assert(Primal_Edge_Right_Endpoint[edge_id]==a);
                Prev_Node_Right_Fixed[edge_id]=b;
            }
        }
    }
}
int find_next_in_order(int x, int y) {
    int edge_id=Edge_Inverse[{x, y}];
    assert(1<=edge_id && edge_id<=Edge_Count);

    if (Primal_Edge_Left_Endpoint[edge_id]==y) {
        return Next_Node_Left_Fixed[edge_id];
    } else {
        return Next_Node_Right_Fixed[edge_id];
    }
}
int find_prev_in_order(int x, int y) {
    int edge_id=Edge_Inverse[{x, y}];
    assert(1<=edge_id && edge_id<=Edge_Count);
    if (Primal_Edge_Left_Endpoint[edge_id]==y) {
        return Prev_Node_Left_Fixed[edge_id];
    } else {
        return Prev_Node_Right_Fixed[edge_id];
    }
}
void find_primal_area(int x, int y, int id) {
    assert(1<=x && x<=Edge_Count && 1<=y && y<=Edge_Count);

    int edge=Edge_Inverse[{x, y}];
    assert(edge!=0);
    bool side=(Primal_Edge_Left_Endpoint[edge]==x ? 0 : 1); // 0 - left, 1 - right side
    if (side==0) {
        if (Primal_Edge_Left_Side[edge]) {
            return;
        }
        Primal_Edge_Left_Side[edge]=id;
    }
    if (side==1) {
        if (Primal_Edge_Right_Side[edge]) {
            return;
        }
        Primal_Edge_Right_Side[edge]=id;
    }

    int z=find_next_in_order(x, y);
    find_primal_area(y, z, id);

}
void find_areas() {
    for (int i=1; i<=Edge_Count; i++) {
        if (!Primal_Edge_Left_Side[i]) {
            Primal_Area_Count++;
            int x=Primal_Edge_Left_Endpoint[i], y=Primal_Edge_Right_Endpoint[i];
            find_primal_area(x, y, Primal_Area_Count);
        }
        if (!Primal_Edge_Right_Side[i]) {
            Primal_Area_Count++;
            int x=Primal_Edge_Right_Endpoint[i], y=Primal_Edge_Left_Endpoint[i];
            find_primal_area(x, y, Primal_Area_Count);
        }
    }
    Dual_Node_Count=Primal_Area_Count;
}

void find_special_st_path() {
    vector<bool> Is_Visited(Node_Count+1, 0);
    vector<int> Previous_Node(Node_Count+1, 0);
    queue<int> q;
    q.push(Starting_Node);
    Is_Visited[Starting_Node]=1;

    while (q.size()>0) {
        int Cur_Node=q.front();
        q.pop();
        for (auto Next_Node:Primal_Adj_List[Cur_Node]) {
            if (!Is_Visited[Next_Node]) {
                q.push(Next_Node);
                Is_Visited[Next_Node]=1;
                Previous_Node[Next_Node]=Cur_Node;
            }
        }
    }
    vector<int> st_path;
    int pos=Finishing_Node;
    while (pos!=Starting_Node) {
        st_path.push_back(pos);
        pos=Previous_Node[pos];
    }
    st_path.push_back(Starting_Node);
    reverse(st_path.begin(), st_path.end());

    Special_st_Path=st_path;
}

void process_the_st_path() {
    int path_len=Special_st_Path.size();
    for (int i=1; i<path_len; i++) {
        int a=Special_st_Path[i-1], b=Special_st_Path[i], edge_id=Edge_Inverse[{a, b}];
        assert(1<=edge_id && edge_id<=Edge_Count);
        if (a==Primal_Edge_Left_Endpoint[edge_id]) {
            Primal_Edge_Weight[edge_id]=1;
        } else {
            Primal_Edge_Weight[edge_id]=-1;
        }
    }
}



void read_SRLGs() {
    cin >> SRLG_Count;

    vector<pair<int, int> > Area_Adj_List[Primal_Area_Count+1];
    vector<bool> Is_Visited(Primal_Area_Count+1, 0);
    vector<int> Dist(Primal_Area_Count+1, 0);
    set<int> Important_Areas;
    int special_SRLG_id=0, special_SRLG_start=0;
    queue<int> q;
    vector<int> edges;

    for (int i=1; i<=SRLG_Count; i++) {
        bool cutting_srlg=0;
        int SRLG_Size;
        cin >> SRLG_Size;
        if (SRLG_Size==0) {
            continue;
        }
        max_srlg_size=max(max_srlg_size,SRLG_Size);
        int st=0;
        for (int j=1; j<=SRLG_Size; j++) {
            int Edge_id;
            cin >> Edge_id;
            Edge_id=Edge_Id_Inverse[Edge_id];
            if (Edge_id==0) {
                cout << "Wrong input: there is a problem with SRLG " << j << "\n";
                exit(0);
            }
            edges.push_back(Edge_id);
            int a=Primal_Edge_Left_Side[Edge_id], b=Primal_Edge_Right_Side[Edge_id], c=Primal_Edge_Weight[Edge_id];
            Important_Areas.insert(a), Important_Areas.insert(b);
            Area_Adj_List[a].push_back({b, c}), Area_Adj_List[b].push_back({a, -c});
            if (!st) {
                st=a;
            }
        }
        int vis_count=0, important_count=Important_Areas.size();
        q.push(st);
        Is_Visited[st]=1;
        while (q.size()>0) {
            int Cur_Area=q.front();
            vis_count++;
            q.pop();
            for (auto x:Area_Adj_List[Cur_Area]) {
                int Next_Area=x.first, Next_Dist=Dist[Cur_Area]+x.second;
                if (!Is_Visited[Next_Area]) {
                    Is_Visited[Next_Area]=1;
                    Dist[Next_Area]=Next_Dist;
                    q.push(Next_Area);
                } else {
                    if (Dist[Next_Area]!=Next_Dist) {
                        special_SRLG_id=i;
                        special_SRLG_start=Cur_Area;
                        cutting_srlg=1;
                    }
                }
            }
        }

        if (vis_count!=important_count) {
            cout << "Wrong input: SRLG " << i << " " << "is not connected\n";
            exit(0);
        }
        if (!cutting_srlg || !eliminate_cutting_srlg) {
            Dual_Node_Count++;
            for (auto x:Important_Areas) {
                Dual_Adj_List[x].push_back({Dual_Node_Count, -Dist[x], 0});
                Dual_Adj_List[Dual_Node_Count].push_back({x, Dist[x], i});
            }
            for (auto x:edges) {
                SRLG_list[x].push_back(i);
            }
        }

        for (auto x:Important_Areas) {
            Area_Adj_List[x].clear();
            Is_Visited[x]=0;
            Dist[x]=0;
        }
        Important_Areas.clear();
        edges.clear();
    }

    if (!eliminate_cutting_srlg && special_SRLG_id) {
        Min_Ratio_Cut_Dual_Node.push_back(special_SRLG_start);
        Min_Ratio_Cut_SRLG_Id.push_back(special_SRLG_id);
        SRLG_Disjoint_Paths.push_back(Special_st_Path);
        print_the_answer();
    }
    return;
}
void add_dual_edges() {
    for (int i=1; i<=Edge_Count; i++) {
        int a=Primal_Edge_Left_Side[i], b=Primal_Edge_Right_Side[i], w=Primal_Edge_Weight[i];
        Dual_Adj_List[a].push_back({b, w, -i});
        Dual_Adj_List[b].push_back({a, -w, -i});
    }
    return;
}
int min_cut_approximation() {
    //vector<int> Odd_Dist(Dual_Node_Count+1, 0), Even_Dist(Dual_Node_Count+1, 0);
    vector<int> Dist[2], Prev_SRLG[2];
    vector<pair<int, int> > Prev_Node[2];
    for (int i=0; i<2; i++) {
        Dist[i].resize(Dual_Node_Count+1);
        Prev_SRLG[i].resize(Dual_Node_Count+1);
        Prev_Node[i].resize(Dual_Node_Count+1);
    }

    int Inf=2*Dual_Node_Count+5;
    int Min_Cut=Inf;
    for (int First_Node=1; First_Node<=Primal_Area_Count; First_Node++) {
        for (int i=1; i<=Dual_Node_Count; i++) {
            Dist[0][i]=Inf;
            Dist[1][i]=Inf;
        }
        //deque<pair<pair<int, bool>, int> > q; // 0-1 bfs
        deque<pair<int, int> > q;
        Dist[0][First_Node]=0;
        q.push_front({First_Node, 0});
        while (q.size()>0) {
            int Current_Node=q.front().first, Parity=q.front().second;
            int Current_Dist=Dist[Parity][Current_Node];
            q.pop_front();

            for (auto x:Dual_Adj_List[Current_Node]) {
                int Next_Node=x.dest;
                int Next_Parity=abs((Parity+x.weight)%2), Next_Dist=Current_Dist+(x.edge_id!=0 ? 1 : 0);
                if (Next_Dist<Dist[Next_Parity][Next_Node]) {
                    Dist[Next_Parity][Next_Node]=Next_Dist;
                    Prev_SRLG[Next_Parity][Next_Node]=x.edge_id;
                    Prev_Node[Next_Parity][Next_Node]={Parity, Current_Node};
                    if (Current_Dist==Next_Dist) {
                        q.push_front({Next_Node, Next_Parity});
                    } else {
                        q.push_back({Next_Node, Next_Parity});
                    }
                }

            }
        }
        int res=Dist[1][First_Node];
        if (res<Min_Cut) {
            Min_Cut=res;
            Min_Cut_SRLG.clear();
            int Current_Node=First_Node, Parity=1;
            while (Current_Node!=First_Node || Parity!=0) {
                auto Prev=Prev_Node[Parity][Current_Node];
                int SRLG=Prev_SRLG[Parity][Current_Node];
                if (SRLG) {
                    Min_Cut_SRLG.push_back(SRLG);
                }
                Parity=Prev.first, Current_Node=Prev.second;
            }
        }
    }



    if (Min_Cut==Inf) {
        cout << "There is a path with safe edges, so the answer is infinity\n";
        exit(0);
    }

    return Min_Cut;
}

vector<int> find_potential_SPFA(int n, vector<vector<pair<int, int> > > &Adj_List) {
    vector<int> Potential(n+1, 0);
    vector<int> Opt_Dist(n+1, 0);
    vector<bool> in_queue(n+1, 0);
    queue<int> q;
    int First_Node=1;
    for (int i=1; i<=n; i++) {
        if (Adj_List[i].size()>Adj_List[First_Node].size()) {
            First_Node=i;
        }
    }
    int inf=1e9;
    for (int i=1; i<=n; i++) {
        Potential[i]=inf;
    }
    Potential[First_Node]=0;
    q.push(First_Node);
    in_queue[First_Node]=1;

    while (q.size()>0) {
        int Cur_Node=q.front();
        if (Cur_Node==First_Node && Potential[First_Node]<0) {
            Potential[0]=-1;
            return Potential;
        }
        q.pop();
        in_queue[Cur_Node]=0;
        if (Opt_Dist[Cur_Node]>n) {
            Potential[0]=-1;
            return Potential;
        }

        for (auto x:Adj_List[Cur_Node]) {
            int Next_Node=x.first, Next_Dist=Potential[Cur_Node]+x.second;
            if (Potential[Next_Node]>Next_Dist) {
                Potential[Next_Node]=Next_Dist;
                Opt_Dist[Next_Node]=Opt_Dist[Cur_Node]+1;
                if (!in_queue[Next_Node]) {
                    in_queue[Next_Node]=1;
                    q.push(Next_Node);
                }
            }
        }
    }
    return Potential;
}




bool find_srlg_disjoint_path(int Val) {
    assert(Val>0);
    vector<vector<pair<int, int> > > Adj_List(Dual_Node_Count+1);
    for (int Cur_Node=1; Cur_Node<=Dual_Node_Count; Cur_Node++) {
        for (auto x:Dual_Adj_List[Cur_Node]) {
            int Next_Node=x.dest, edge_weight=x.weight, id=x.edge_id;
            int Next_Weight=edge_weight*Val+(id!=0 ? 1 : 0);
            Adj_List[Cur_Node].push_back({Next_Node, Next_Weight});
        }
    }
    vector<int> Potential;
    Potential=find_potential_SPFA(Dual_Node_Count, Adj_List);

    if (Potential[0]==-1) {
        return false;
    }
    SRLG_Disjoint_Paths.clear();
    if (Val==1) {
        SRLG_Disjoint_Paths.push_back(Special_st_Path);
        return true;
    }
    /*cout << "potencial:\n";
    for (int i=1; i<=Dual_Node_Count; i++) {
        cout << Potential[i] << " ";
    }
    cout << "\n";
    cout << "ends " << Val << "\n";
    cout << "start " << Starting_Node << " " << "finishing node " << Finishing_Node << "\n";*/
    for (int path_id=0; path_id<Val; path_id++) {
        set<int> used_edges;
        vector<int> general_st_path, single_st_path;

        int Prev_Node=Primal_Adj_List[Starting_Node][0], Cur_Node=Starting_Node; // itt prev_node nem is szamit

        while (Cur_Node!=Finishing_Node) {
            general_st_path.push_back(Cur_Node);
            /*if (general_st_path.size()>100) {
                for (auto x:general_st_path) {
                    cout << x << " ";
                }
                cout << "\n";
                exit(0);
            }*/
            assert(general_st_path.size()<=10*Node_Count);
            int Next_Node=find_next_in_order(Prev_Node, Cur_Node);
            while (true) {
                int edge_id=Edge_Inverse[{Cur_Node, Next_Node}];
                //cout << "proba " << Prev_Node << " " << Cur_Node << " " << Next_Node << "\n";
                if (used_edges.count(edge_id)) {
                    Next_Node=find_next_in_order(Next_Node, Cur_Node);
                    continue;
                }
                int left_side=Potential[Primal_Edge_Left_Side[edge_id]], right_side=Potential[Primal_Edge_Right_Side[edge_id]];
                int weight=Primal_Edge_Weight[edge_id];
                left_side+=weight*Val;

                //cout << "fontos " << Cur_Node << " " << Next_Node << " " << edge_id << " " << Primal_Edge_Left_Endpoint[edge_id] << " " << Primal_Edge_Right_Endpoint[edge_id] << "\n";
                if (Primal_Edge_Left_Endpoint[edge_id]!=Cur_Node) {
                    assert(Primal_Edge_Right_Endpoint[edge_id]==Cur_Node);
                    swap(left_side, right_side);
                }
                int left_mod=(left_side%Val+Val)%Val, right_mod=(right_side%Val+Val)%Val;

                if (right_mod==path_id && left_mod!=path_id) {
                    /*if ((right_mod+1)%Val!=left_mod) {
                        cout << "baj " << Cur_Node << " "  << Next_Node << "\n";
                        cout << "ut: " << path_id << "\n";
                        cout << "jobb: " << right_side << " " << right_mod << "\n";
                        cout << "bal: " << left_side << " " <<left_mod << "\n";
                        exit(0);
                    }
                    assert((right_mod+1)%Val==left_mod);*/
                    if (right_side+1==left_side) {
                        used_edges.insert(edge_id);
                        break;
                    }
                }
                Next_Node=find_next_in_order(Next_Node, Cur_Node);
            }
            Prev_Node=Cur_Node;
            Cur_Node=Next_Node;
        }
        general_st_path.push_back(Finishing_Node);
        vector<int> Node_in_Path(Node_Count+1, 0);
        for (auto x:general_st_path) {
            Node_in_Path[x]++;
        }
        int si=general_st_path.size(), pos=0;
        while (pos<si) {
            int Cur_Node=general_st_path[pos];
            single_st_path.push_back(Cur_Node);
            while (Node_in_Path[Cur_Node]!=0) {
                Node_in_Path[general_st_path[pos]]--;
                pos++;
            }
        }
        SRLG_Disjoint_Paths.push_back(single_st_path);
    }
    return true;
}

bool is_legal_edge(int x, int y, set<int> &Blocked_SRLG) {
    //cout << "kerdes " << x << " " << y << "\n";
    int edge_id=Edge_Inverse[{x, y}];
    for (auto p:SRLG_list[edge_id]) {
        if (Blocked_SRLG.count(p)) {
            return false;
        }
    }
    if (Blocked_SRLG.count(-edge_id)) {
        return false;
    }
    return true;
}

bool check_two_paths(vector<int> &Path_One, vector<int> &Path_Two) {
    int Size_One=Path_One.size(), Size_Two=Path_Two.size();
    set<int> Blocked_SRLG;
    for (int i=0; i<Size_One-1; i++) {
        int edge_id=Edge_Inverse[{Path_One[i], Path_One[i+1]}];
        for (auto x:SRLG_list[edge_id]) {
            Blocked_SRLG.insert(x);
        }
        if (Edge_Disjoint_Paths) {
            Blocked_SRLG.insert(-edge_id);
        }
    }

    for (int i=0; i<Size_Two-1; i++) {
        int x=Path_Two[i], y=Path_Two[i+1];
        if (!is_legal_edge(x, y, Blocked_SRLG)) {
            return false;
        }
    }
    return true;
}

bool is_common_edge(vector<int> &Path_One, vector<int> &Path_Two) {
    int Size_One=Path_One.size(), Size_Two=Path_Two.size();
    set<int> edges;
    for (int i=0; i<Size_One-1; i++) {
        int edge_id=Edge_Inverse[{Path_One[i], Path_One[i+1]}];
        edges.insert(edge_id);
    }

    for (int i=0; i<Size_Two-1; i++) {
        int edge_id=Edge_Inverse[{Path_Two[i], Path_Two[i+1]}];
        if (edges.count(edge_id)) {
            return true;
        }
    }
    return false;
}

bool find_greedy_path_dfs(set<int> &Blocked_SRLG, vector<bool> &visited, vector<int> &answer, int spec_node) {
    //cout << "greedy_path_dfs: ";
    //print_vector(answer);
    //return true;
    if (answer.back()==Finishing_Node) {
        return true;
    }

    int si=answer.size();
    int Prev_Node=(si==1 ? spec_node : answer[si-2]), Cur_Node=answer[si-1], Next_Node=find_prev_in_order(Prev_Node, Cur_Node);
    bool crossing=0;
    while (true) {
        int edge_id=Edge_Inverse[{Cur_Node, Next_Node}];
        if (Blocked_SRLG.count(-edge_id)) {
            crossing=1-crossing;
        }
        /*if (Cur_Node==3 && Next_Node==8) {
            //cout << "elromlik........................\n";
            cout << crossing << " " << visited[Next_Node] << " " << is_legal_edge(Cur_Node, Next_Node, Blocked_SRLG) << "\n";
        }*/
        if (!crossing && !visited[Next_Node] && is_legal_edge(Cur_Node, Next_Node, Blocked_SRLG)) {
            answer.push_back(Next_Node);
            visited[Next_Node]=1;
            if (find_greedy_path_dfs(Blocked_SRLG, visited, answer, 0)) {
                return true;
            }
            //visited[Next_Node]=0;
            answer.pop_back();
        }
        Next_Node=find_prev_in_order(Next_Node, Cur_Node);
        if (Next_Node==Prev_Node) {
            //if (Cur_Node==Starting_Node) {
            //    cout << "nem talal uj utat\n";
            //    exit(0);
            //}
            assert(Cur_Node!=Starting_Node);
            return false;
        }
    }

    //cout << "nem kene itt lennie\n";
    assert(false);
    exit(0);

}

vector<int> find_next_greedy_path(vector<int> &Last_Path) {
    //cout << "alap\n";
    //for (auto x:Last_Path) {
    //    cout << x << " ";
    //}
    cout << "\n";
    set<int> Blocked_SRLG;
    int Path_Len=Last_Path.size();
    for (int i=0; i<Path_Len-1; i++) {
        int edge_id=Edge_Inverse[{Last_Path[i], Last_Path[i+1]}];
        for (auto x:SRLG_list[edge_id]) {
            //cout << "specialis SRLG: " << x << "\n";
            Blocked_SRLG.insert(x);
        }
        if (Edge_Disjoint_Paths) {
            Blocked_SRLG.insert(-edge_id);
        }
    }
    vector<bool> visited(Node_Count+1, 0);

    vector<int> answer;
    answer.push_back(Starting_Node);
    visited[Starting_Node]=1;
    bool x=find_greedy_path_dfs(Blocked_SRLG, visited, answer, Last_Path[1]);
    assert(x);
    //cout << "vege " << Finishing_Node << " " << answer.back() << "\n";

    /*cout << "find_next_greedy_path:\n";
    for (auto x:Last_Path) {
        cout << x << " ";
    }
    cout << "\n";
    for (auto x:answer) {
        cout << x << " ";
    }
    cout << "\n";*/
    //exit(0);

    return answer;
}

vector<int> clockwise_maximum(vector<int> &Ref_Path, vector<int> &Path_One, vector<int> &Path_Two) {
    vector<int> answer;
    answer.push_back(Starting_Node);
    int Prev_Node=Starting_Node, Cur_Node=Ref_Path[1];
    do {
        Cur_Node=find_next_in_order(Cur_Node, Prev_Node);
    } while (Cur_Node!=Path_One[1] && Cur_Node!=Path_Two[1]);

    vector<int> Node_Inverse_One(Node_Count+1, -1), Node_Inverse_Two(Node_Count+1, -1);
    for (int i=0; i<Path_One.size(); i++) {
        Node_Inverse_One[Path_One[i]]=i;
    }
    for (int i=0; i<Path_Two.size(); i++) {
        Node_Inverse_Two[Path_Two[i]]=i;
    }

    /*cout << "Ref_path: ";
    print_vector(Ref_Path);
    cout << "Path_one: ";
    print_vector(Path_One);
    cout << "Path_two: ";
    print_vector(Path_Two);*/

    answer.push_back(Cur_Node);
    while (Cur_Node!=Finishing_Node) {
        int x=Node_Inverse_One[Cur_Node], y=Node_Inverse_Two[Cur_Node];
        int Next_Node=Prev_Node;
        assert(x!=-1 || y!=-1);
        if (x==-1) {
            Next_Node=Path_Two[y+1];
        } else if (y==-1) {
            Next_Node=Path_One[x+1];
        } else {
            //cout << "fontos " << Path_One[x+1] << " " << Path_Two[y+1] << "\n";
            while (Next_Node!=Path_One[x+1] && Next_Node!=Path_Two[y+1]) {
                //cout << "baj " << Next_Node << " " << Cur_Node << "\n";
                Next_Node=find_next_in_order(Next_Node, Cur_Node);
            }
        }
        answer.push_back(Next_Node);
        Prev_Node=Cur_Node;
        Cur_Node=answer.back();
    }


    return answer;
}

int main()
{
   chrono::steady_clock::time_point begin = chrono::steady_clock::now();
    max_srlg_size=0;
    read_nodes();
    read_edges();
    check_primal_connectivity();
    cout<<"<graph>\n";
    cout<<"<V>"<< Node_Count<<"</V>\n";
    cout<<"<E>"<< Edge_Count<<"</E>\n";
    order_neighbours();
    find_areas();
    find_special_st_path();
    process_the_st_path();
    read_SRLGs();
    chrono::steady_clock::time_point begin2 = chrono::steady_clock::now();
    cout<<"<srlg>\n";
    cout<<"<srlgnum>"<< SRLG_Count <<"</srlgnum>\n";
    cout<<"<srlgMaxSize>"<< max_srlg_size <<"</srlgMaxSize>\n";
    cout<<"</srlg>\n";
    cout<<"<run>\n";
    cout<<"<s>"<< Starting_Node <<"</s>\n";
    cout<<"<t>"<< Finishing_Node <<"</t>\n";
    int start_degree=0;
    for (auto Next_Node:Primal_Adj_List[Starting_Node]) start_degree++;
    int finish_degree=0;
    for (auto Next_Node:Primal_Adj_List[Finishing_Node]) finish_degree++;
    max_path=std::min(start_degree,finish_degree);
    if (Edge_Disjoint_Paths) {
        add_dual_edges();
    }

    string alg_name="regular";
    #ifdef GREEDY
      bool greedy_like_algorithm=1;
      alg_name="greedy";
    #else
      bool greedy_like_algorithm=0;
    #endif
    #ifdef TWOPATH
      bool try_to_find_two_disjoint_paths=1; //!!!!
      alg_name+="-two-path";
    #else
      bool try_to_find_two_disjoint_paths=0; //!!!!
    #endif
    #ifdef BINSEARCHCUT
      bool replace_min_cut_approx_with_binsearch=1; //!!!!
      alg_name+="-binsearch";
    #else
      bool replace_min_cut_approx_with_binsearch=0; //!!!!
    #endif

    if (greedy_like_algorithm) {
          bool start=find_srlg_disjoint_path(2);
          if (!start) {
              find_srlg_disjoint_path(1);
              print_the_answer();
          } else {
              Greedy_Paths.push_back(SRLG_Disjoint_Paths[0]);
              int last=0;
              int cnt=0, len=1;
              while (cnt<=Primal_Area_Count+10) {
                  Greedy_Paths.push_back(find_next_greedy_path(Greedy_Paths[len-1]));
                  if (Greedy_Paths.size()>2) {
                      vector<int> Next_Path=clockwise_maximum(Greedy_Paths[len-1], Greedy_Paths[len], Greedy_Paths[last]);
                      if (Greedy_Paths[last]!=Next_Path || Greedy_Paths[len]==Greedy_Paths[last]) {
                          last++, cnt++;
                          Greedy_Paths[len]=Next_Path;
                      } else {
                          cnt=0;
                      }
                  }
                  len++;
              }
              int si=Greedy_Paths.size();
              int first=si-1;
              if (!check_two_paths(Greedy_Paths[last], Greedy_Paths[first])) {
                  last++;
              }
              SRLG_Disjoint_Paths.clear();
              for (int i=last; i<=first; i++) {
                  SRLG_Disjoint_Paths.push_back(Greedy_Paths[i]);
              }
              print_the_answer();
          }
      } else if (try_to_find_two_disjoint_paths) {
          if (find_srlg_disjoint_path(2)) {
              cout << "There are at least two disjoint paths:\n";
              print_the_answer();
          } else {
              cout << "No two disjoint paths\n";
          }
      } else if (replace_min_cut_approx_with_binsearch) {
          int lo=0, hi=Edge_Count+1, mid;
          //cout<<"max_path: "<<max_path;
          hi=max_path+1;
          while (hi-lo>1) {
              mid=(hi+lo)/2;
              if (find_srlg_disjoint_path(mid)) {
                  lo=mid;
              } else {
                  hi=mid;
              }
          }
          print_the_answer();
      } else {
          Min_Cut_Approx=min_cut_approximation();
        cout<<"<mincut>"<< Min_Cut_Approx <<"</mincut>\n";
        for (int i=Min_Cut_Approx+1; i>=Min_Cut_Approx-3; i--) {
            assert(i!=Min_Cut_Approx-3);
            if (find_srlg_disjoint_path(i)) {
                print_the_answer();
                break;
            }
        }
    }
    chrono::steady_clock::time_point end = chrono::steady_clock::now();
    cout << "<runtime>" << std::chrono::duration_cast<std::chrono::microseconds> (end - begin).count() << "</runtime>\n";
    cout << "<runtime2>" << std::chrono::duration_cast<std::chrono::microseconds> (end - begin2).count() << "</runtime2>\n";
    //std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::nanoseconds> (end - begin).count() << "[ns]" << std::endl;
    cout << "<algorithm>" <<alg_name<< "</algorithm>";
    cout<<"</run>\n";
    cout<<"</graph>\n";
    return 0;
}
// THE EXISTING CODE SUPPLEMENTED WITH SOME UNWRITTEN FUNCTIONS REQUIRED FOR THE OLD ALGORITHM
/*
Greedy algorithm
1.
Starting from an S-T path that has an independent SRLG
(one of the two paths must be chosen)
2.
After a path, the first SRLG-independent to its left needs to be selected
note the SRLGs intersected by the path (edges need special attention)
if the path starts with S-X then start from the X-S edge and S vertex
always search for the first vertex after the current edge that is not part of the selected SRLG set
this will work with a dfs-like algorithm; if stuck or returning to the old path, simply return false
3.
For a simple greedy algorithm (slow and not necessarily correct), repeating step 2 with the last active path to examine
and if not SRLG-independent, then discard the last path
this way, I think the correct answer will be obtained in most cases
4.
For the Dervish algorithm a bit more is needed
need to determine the maximum of two paths, for which a third is needed (for comparison)
besides the next one on the left, it would be good to have the next one on the right as well
*/
/*
easy test, with many independent paths
7
0 0 1
1 0 2
1 1 3
1 2 4
1 3 5
1 4 6
2 0 7
1 7
10
1 2 1
1 3 2
1 4 3
1 5 4
1 6 5
2 7 6
3 7 7
4 7 8
5 7 9
6 7 10
0

*/
