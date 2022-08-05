#ifndef MESH_HPP
#define MESH_HPP

#include <iostream>
#include <map>
#include <Eigen/Sparse>

using namespace std;
using namespace Eigen;

typedef Vector2d R2;
typedef Matrix<double,2,2> R22;
typedef VectorXi Tri;
typedef Vector4i Edge;

typedef Triplet<double> T;
typedef SparseMatrix<double> SpMat;
typedef Matrix<double,6,6> R66;
typedef Matrix<double,6,3> R63;
typedef Matrix<double,3,3> R33;

class mesh{

private :

  int m_Nv; // NOMBRE DE SOMMET DU MAILLAGE

  int m_Nt; // NOMBRE DE TRIANGLE DU MAILLAGE

  int m_Ne; // NOMBRE D'ARETE DU MAILLAGE

  int m_Neb; // NOMBRE D'ARETE SUR LE BORD DU MAILLAGE

  int m_Nlbl; // NOMBRE DE LABEL SUR LE BORD DU MAILLAGE

  map<int, R2> m_Nodes; // MAPAGE DES COORDONNÉES DES NOEUDS DU MAILLAGE

  map<int, vector<int> > m_lbl_Nodes; // MAPAGE DES NOEUDS QUI APPARTIENNENT A UN BORD LABELISÉ

  map<int, Tri> m_Element; // MAPAGE DES TRIANGLES DU MAILLAGE

  map<int, vector<int> > m_lbl_Element; // MAPAGE DES TRIANGLES QUI ONT UNE ARETE SUR UN BORD LABELISÉ
  
  map<int, Edge> m_Boundary; // MAPAGE DES ARETES QUI APPARTIENNENT AU BORD DU  MAILLAGE
  
  map<int, vector<int> > m_lbl_Boundary; // MAPAGE DES ARETES QUI APPARTIENNENT A UN BORD LABELISÉ

  map<int, Edge> m_Edge; // MAPAGE DE TOUTES LES ARETES

  map<int, vector<Edge> > m_Lbl_Edge; // MAPAGE DE TOUTES LES ARETES PAR LABEL. P.S. : LE LABEL 0 CORRESPOND AUX ARETES QUI SONT A L'INTÉRIEUR DU MAILLAGE

  map<int, int> m_Lbl_Size; // MAPAGE DU NOMBRE D'ARETES QUI APPARTIENNENT A UN LABEL

  map<int, vector<int> > m_Neighbors; // MAPAGE DES TRIANGLES QUI PARTAGENT UN MÊME SOMMET

  map<int, int> m_NbofNeighbors; // MAPAGE DU NOMBRE DE TRIANGLE QUI PARTAGE UN MÊME SOMMET
  
  char* m_filename; // NOM DU FICHIER DE MAILLAGE

public :

  mesh(char* filename); // CONSTRUCTEUR DE MON OBJET DE LA CLASSE MESH

  int Nv() const ; // GET DU MEMBRE Nv

  int Nt() const ; // GET DU MEMBRE Nt

  int Ne() const ; // GET DU MEMBRE Ne

  int Neb() const ; // GET DU MEMBRE Neb

  int Nlbl() const; // GET DU MEMBRE Nlbl

  map<int, Edge> GetEdge(); // GET DU MEMBRE Edge

  map<int, int> Lbl_Size() const; // GET DU MEMBRE Lbl_Size
  
  map<int, R2> Nodes() const; // GET DU MEMBRE Nodes

  map<int, vector<int> > Lbl_Nodes() const; // GET DU MEMBRE Lbl_Nodes

  map<int, Tri> Element() const; // GET DU MEMBRE Element
  
  map<int, vector<int> > Lbl_Element() const; // GET DU MEMBRE Lbl_Element

  map<int, Edge> Boundary() const; // GET DU MEMBRE Boundary
  
  map<int, vector<int> > Lbl_Boundary() const; // GET DU MEMBRE Lbl_Boundary

  map<int, Edge> Edg(); // CONSTRUCTEUR ET GET DU MEMBRE Edge

  map<int, vector<Edge> > Lbl_Edge(); // CONSTRUCTEUR ET GET DU MEMBRE Lbl_Edge

  map<int, vector<int> > Neighbors() const; // GET DU MEMBRE Neighbors

  map<int, int> NbofNeighbors() const; // GET DU MEMBRE NbofNeighbors

};

#endif
