#ifndef FUNCTION_HPP
#define FUNCTION_HPP

#include <vector>
#include <Eigen/Sparse>
#include <cassert>

using namespace std;

typedef Eigen::Triplet<double> T;
typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::Matrix<double,6,6> RigMat;
typedef Eigen::Matrix<double,6,3> DivMat;
typedef Eigen::Matrix<double,3,3> PenMat;
typedef Eigen::Matrix<double,6,2> Mat62;
typedef Eigen::VectorXd VEigen;

typedef vector<double> vecd;
typedef vector<int> veci;

void Dim(int& Nodes, int& Element, int& Edge, int& Boundary, const char* filename);

void Load(VEigen* Nodes, veci* Element, veci* Boundary, const char* filename);

void Trian(int Nbv, int Nbt, veci* Element, veci* trian);

//void Edge(int Nbv, int Nbt, int Nbe, VEigen* Nodes, veci* Element, veci* edge);

void Edge(int Nbv, int Nbt, int Nbe, veci* Element, veci* edge, veci* Numedge);

veci IdNodesBound(int Nbv, int Nbe, int Nbeb, veci* edge, veci* Boundary);

veci IdNodesLbl(veci& idnodesbound, VEigen* MidPnt, int lbl);

void ndofP2(int& Nbv, int& Nbe, VEigen* Nodes, veci* edge, VEigen* MidPnt);

void dofP2(int Nbt, int Nbv, veci* Element, veci* Numedge, veci* ddl);

bool IsSymetric(const SpMat& M);

double g(VEigen pnt);

VEigen Gradg(VEigen pnt);

double lambda(const VEigen& pnt, int i);

VEigen Gradlambda(int i);

double divlambda(int i);

double phi(const VEigen& pnt, int i);

VEigen Gradphi(const VEigen& pnt, int i);

double dxdyphi(const VEigen& pnt, int i, int j);

RigMat MassEle(const VEigen* Nodes);

SpMat Mass(VEigen* MidPnt, veci* ddl, veci* Element, int Nbv, int Nbt, int Nbe);

RigMat RigEle(const VEigen* Nodes);

SpMat Rig(VEigen* MidPnt, veci* ddl, veci* Element, int Nbv, int Nbt, int Nbe);

DivMat DivEle(const VEigen* Nodes, int l);

SpMat Div(VEigen* MidPnt, veci* ddl, veci* Element, int Nbv, int Nbt, int Nbe);

PenMat PenEle(const VEigen* Nodes);

SpMat Pen(VEigen* MidPnt, veci* ddl, veci* Element, int Nbv, int Nbt, int Nbe);

SpMat AssemSys(VEigen* MidPnt, veci* ddl, veci* Element, const veci& IdNodesBound, const veci& IdNodesLbl, int Nbv, int Nbt, int Nbe, double eps, double nu, double alp, VEigen& RightHandSide);

//SpMat AssemSysNS(veci* trian, double dt, const VEigen& u_n, VEigen* MidPnt, veci* ddl, veci* Element, const veci& IdNodesBound, const veci& IdNodesLbl, int Nbv, int Nbt, int Nbe, double eps, double nu, VEigen& RightHandSide);

int StructMatGlob(const SpMat& M, bool coef);

void SavetoPlot(int Nbe, int Nbv, int Nbt, const veci* Element, const VEigen& x, const char* filename);

//void SavetoPlot(int Nbv, int Nbe, int Nbt, const VEigen* MidPnt, const veci* ddl, const VEigen& x);

//int WhereWasI(int Nbv, int Nbe, VEigen* Nodes, VEigen& a_l, int n_t, veci* Element, veci* trian, VEigen& u_n, double& dt);

//Mat62 RightHandSideEle(double dt, const VEigen& u_n, const VEigen* Nodes, int Nbv, int Nbe, const VEigen* NodesP2, int n_t, const veci* Element, veci* trian);










#endif
