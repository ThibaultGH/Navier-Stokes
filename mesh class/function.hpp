#ifndef FUNCTION_HPP
#define FUNCTION_HPP

#include <Eigen/Sparse>

#include "mesh.hpp"

using namespace std;
using namespace Eigen;

typedef Vector2d R2;
typedef Matrix<double,2,2> R22;
typedef VectorXi Tri;
typedef Vector4i Edge;

typedef Triplet<double> T;
typedef SparseMatrix<double> SpMat;
typedef VectorXd Vec;
typedef Matrix<double,6,6> R66;
typedef Matrix<double,6,3> R63;
typedef Matrix<double,3,3> R33;

double g(R2 pnt); // FONCTION IMPOSER SUR LE BORD D'ENTREE

double psi(const R2& pnt, int i); // FONCTION DE FORME P1-LAGRANGE SUR LE TRIANGLE DE REF

R2 Gradpsi(int i); // GRADIENT FONCTION DE FORME P1-LAGRANGE SUR LE TRIANGLE DE REF

double phi(const R2& pnt, int i); // FONCTION DE FORME P2-LAGRANGE SUR LE TRIANGLE DE REF

R2 Gradphi(const R2& pnt, int i); // GRADIENT FONCTION DE FORME P2-LAGRANGE SUR LE TRIANGLE DE REF

R2 F_T(const R2* T, const R2& u); // APPLICATION LINEAIRE DU TRIANGLE DE REF VERS UN TRIANGLE QCQ

R2 invF_T(const R2* T, const R2& u); // APPLICATION LINEAIRE D'UN TRIANGLE QCQ VERS LE TRIANGLE DE REF

R66 MassEle(const R2* T); // CALCUL DE LA MATRICE DE MASSE ELEMENTAIRE ASSOCIÉE A UN TRIANGLE QCQ

SpMat Mass(int Nt, int Nv, int Ne, map<int, Tri> Element, map<int, R2> Nodes); // ASSEMBLAGE DE MATRICE DE MASSE DE TOUT LE MAILLAGE POUR TEST

R66 RigEle(const R2* T); // CALCUL DE LA MATRICE DE RIGIDITÉ ELEMENTAIRE ASSOCIÉE A UN TRIANGLE QCQ

SpMat Rig(int Nt, int Nv, int Ne, map<int, Tri> Element, map<int, R2> Nodes); // ASSEMBLAGE DE MATRICE DE RIGIDITÉ DE TOUT LE MAILLAGE POUR TEST

R63 DivEle(const R2* T, int l); // CALCUL DE LA MATRICE DE DIVERGENCE ELEMENTAIRE ASSOCIÉE A UN TRIANGLE QCQ

SpMat Div(int Nt, int Nv, int Ne, map<int, Tri> Element, map<int, R2> Nodes, int l); // ASSEMBLAGE DE MATRICE DE DIVERGENCE DE TOUT LE MAILLAGE POUR TEST

R33 PenEle(const R2* T); // CALCUL DE LA MATRICE DE PÉNALISATION ELEMENTAIRE ASSOCIÉE A UN TRIANGLE QCQ

SpMat Pen(int Nt, int Nv, int Ne, map<int, Tri> Element, map<int, R2> Nodes, int l); // ASSEMBLAGE DE MATRICE DE PÉNALISATION DE TOUT LE MAILLAGE POUR TEST

SpMat AssemSys(mesh& Th, double nu, double alp, double eps, double tgv, Vec& v); // FONCTION ASSEMBLAGE DU SYSTÈME LINÉAIRE

Vec AssemSnd(const mesh& Th, const Vec& u, double dt); // FONCTION ASSEMBLAGE DU SECOND MEMBRE POUR NAVIER-STOKES

void SndBorderCond(mesh& Th, Vec& v, double tgv); // FONCTION PÉNALISATION DU SECOND MEMBRE POUR NAVIER-STOKES

void SavetoPlot(const mesh Th, const Vec& x, const char* filename); // FONCTION DE SAUVEGARDE DE LA SOLUTION DANS UN FORMAT PROMPT POUR FF++

int StructMatGlob(const SpMat& M, bool coef); // FONCTION QUI REVELE LA STRUCTURE CREUSE DE LA MATRICE DU SYSTEME LINEAIRE POUR TEST

#endif
