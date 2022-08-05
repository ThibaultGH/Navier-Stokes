#include <iostream>
#include <fstream>

#include "mesh.hpp"

using namespace std;
using namespace Eigen;

typedef VectorXd VEigen;
typedef Vector2d R2;
typedef VectorXi Tri;
typedef Vector4i Edge;

typedef vector<int> veci;

mesh::mesh(char* filename)
{
  m_filename = filename;
  
  ifstream myflux(m_filename, ios::in);
  
  if(myflux){

    string word;
    
    myflux >> word;
    m_Nv = atoi(word.c_str());
    myflux >> word;
    m_Nt = atoi(word.c_str());
    myflux >> word;
    m_Neb = atoi(word.c_str());

    m_Ne = (3*m_Nt-m_Neb)/2+m_Neb;

    // CONSTRUCTION DE LA MAP DES SOMMETS

    for(int i = 0; i < m_Nv; ++i){

      m_NbofNeighbors[i] = 0;

      R2 temp;

      myflux >> word;
      temp[0] = atof(word.c_str());

      myflux >> word;
      temp[1] = atof(word.c_str());

      m_Nodes[i] = temp;

      myflux >> word;
      m_lbl_Nodes[atoi(word.c_str())].push_back(i);
      
    }

    // FIN CONSTRUCTION DE LA MAP DES SOMMETS

    // CONSTRUCTION DE LA MAP DES TRIANGLES
    
    for(int i = 0; i < m_Nt; ++i){
      Tri temp(6);

      myflux >> word;
      temp[0] = atoi(word.c_str())-1;

      myflux >> word;
      temp[1] = atoi(word.c_str())-1;

      myflux >> word;
      temp[2] = atof(word.c_str())-1;

      m_Neighbors[temp[0]].push_back(i);
      m_NbofNeighbors[temp[0]] += 1;
      m_Neighbors[temp[1]].push_back(i);
      m_NbofNeighbors[temp[1]] += 1;
      m_Neighbors[temp[2]].push_back(i);
      m_NbofNeighbors[temp[2]] += 1;
      
      m_Element[i] = temp;

      myflux >> word;
      m_lbl_Element[atoi(word.c_str())].push_back(i);
      
    }

    // FIN CONSTRUCTION DE LA MAP DES TRIANGLES

    cout << "ENTER YOUR NUMBER OF LABEL ON YOUR BOUNDARY PLEASE : " << endl;
    cin >> m_Nlbl;
    cout << endl;
    for(int i = 0; i < m_Nlbl; ++i){m_Lbl_Size[i] = 0;}

    // CONSTRUCTION DE LA MAP DES ARETES SUR LE BORD
    
    for(int i = 0; i < m_Neb; ++i){

      Edge temp;

      myflux >> word;
      temp[0] = atoi(word.c_str())-1;

      myflux >> word;
      temp[1] = atoi(word.c_str())-1;

      myflux >> word;
      temp[3] = atoi(word.c_str());
      m_lbl_Boundary[atoi(word.c_str())].push_back(i);

      m_Lbl_Size[atoi(word.c_str())] += 1;
      
      m_Boundary[i] = temp;
      
    }
    
    // FIN CONSTRUCTION DE LA MAP DES ARETES SUR LE BORD

    m_Lbl_Size[0] = m_Ne-m_Neb;

  }
}

int mesh::Nv() const {return m_Nv;}

int mesh::Nt() const {return m_Nt;}

int mesh::Ne() const {return m_Ne;}

int mesh::Neb() const {return m_Neb;}

int mesh::Nlbl() const {return m_Nlbl;}

map<int, Edge> mesh::GetEdge() {return m_Edge;}

map<int, int> mesh::Lbl_Size() const {return m_Lbl_Size;}

map<int, R2> mesh::Nodes() const {return m_Nodes;}

map<int, vector<int> > mesh::Lbl_Nodes() const {return m_lbl_Nodes;}

map<int, Tri> mesh::Element() const {return m_Element;}

map<int, vector<int> > mesh::Lbl_Element() const {return m_lbl_Element;}

map<int, Edge> mesh::Boundary() const {return m_Boundary;}

map<int, vector<int> > mesh::Lbl_Boundary() const {return m_lbl_Boundary;}

map<int, Edge> mesh::Edg() {

  /*
    La fonction Edge construit à partir des noeuds P1, ie les sommets des triangles, la suite de la numérotation pour les noeuds P2, cad on a rajouté les milieux des aretes des triangles.
    Elle fait cela en parcourant les triangles et en rangeant les aretes dans des classes définient par : "Deux arètes à appartiennent à la même classe si le plus petit numéro des deux sommets dont eelles sont composées est le même." Par exemple les arêtes (3,85) et (42,3) appartiennent à la même classe, mais les arêtes (3,85) et (2,3) ou bien (3,85) et (16,75) n'appartiennent pas aux mêmes classes.
    On peut ainsi par la même occasion détecté les aretes qui ne sont pas sur le bord, et donc indirectement celles qui le sont, car celles qui ne sont pas sur le bord seront celles qui sont trouvés une deuxième fois par l'algorithme car : "Une arête appartient au bord du maillage ssi elle n'appartient qu'à un seul triangle."
    A noter que ici une arete est modéliser par un vecteur de Eigen de 4 int de la manière suivante :
    _ Premier int : plus petit numéro des deux sommets qui compose l'arete.
    _ Deuxième int : l'autre numéro, cad le plus grand
    _ Troisième int : le numéro de l'arete décaler du nombre de sommet du maillage, cad le numéro du sommet P2 présent au milieu de l'arete.
    _ Quatrième int : Un 0 si l'arete appartient au bord et un 1 si elle appartient à l'intérieur.

    A la fin, on change ce dernier int pour les aretes qui sont sur le bord, ie si il vaut 0 par le label attribué dans le fichier de maillage généré par ff++.
   */
  
  veci first(m_Nv,-1);
  veci next;
  
  int v0,v1;
  int I;
  bool found;
  int Na = 0;
  
  for(int t = 0; t < m_Nt; ++t){
    for(int i = 0; i < 3; ++i){

      v0 = min(m_Element[t][i], m_Element[t][(i+1)%3]);
      v1 = max(m_Element[t][i], m_Element[t][(i+1)%3]);

      I = first[v0];
      found = false;

      while( (found == false) && (I != -1) ){
	
	if ( (v0 == m_Edge[I][0]) && (v1 == m_Edge[I][1]) ){
	  found = true;
	  m_Element[t][i+3] = I+m_Nv;
	  m_Edge[I][2] = I+m_Nv;
	  m_Edge[I][3] = 1;
	}
	
	I = next[I];
      }
      
      if(found == false){
	
	m_Edge[Na][0] = v0;
	m_Edge[Na][1] = v1;
	m_Edge[Na][2] = m_Nv+Na;
	m_Edge[Na][3] = 0;

	R2 x0 = m_Nodes[v0];
	R2 x1 = m_Nodes[v1];
	R2 x = (x0+x1)/2;
	m_Nodes[m_Nv+Na] = x;
	
	m_Element[t][i+3] = m_Nv+Na;

	next.push_back(first[v0]);
	
	first[v0] = Na;
	
	Na += 1;
      }
    }    
  }

  for(int i = 0; i < m_Neb; ++i){
    int v0 = min( m_Boundary[i][0], m_Boundary[i][1] );
    int v1 = max( m_Boundary[i][1], m_Boundary[i][0] );
    int lbl = m_Boundary[i][2];

    for(int j = 0; j < m_Neb; ++j){
      if( (v0 == m_Edge[j][0]) && (v1 == m_Edge[j][1]) ){
	// cout << lbl << endl;
	m_Edge[j][3] = lbl;
	// cout << m_Edge[i] << endl << endl;
      }
    }
    
  }
  
  return m_Edge;

}

map<int, vector<Edge> > mesh::Lbl_Edge() {
  /*
    Cette fonction utilise la construction précedente de Edge pour crée une map qui associe à un numéro de label toutes les aretes, et donc tout les noeuds P1 et P2, qui appartiennent à ce label.
    On en profite aussi pour passer le label des aretes intérieur au maillage à 0.
   */
  
  for(int i = 0; i < m_Ne; ++i){
    int v0 = m_Edge[i][0];
    int v1 = m_Edge[i][1];
    m_Edge[i][3] = 1 - m_Edge[i][3];
    int lbl = m_Edge[i][3];

    for(int j = 0; j < m_Neb; ++j){
      if( (v0 == min(m_Boundary[j][0],m_Boundary[j][1])) && (v1 == max(m_Boundary[j][0],m_Boundary[j][1])) ){
	lbl = m_Boundary[j][3];
	m_Edge[i][3] = lbl;
      }
    }
    
    m_Lbl_Edge[lbl].push_back(m_Edge[i]);

  }
    
  return m_Lbl_Edge;
}

map<int, vector<int> > mesh::Neighbors() const {return m_Neighbors;}

map<int, int> mesh::NbofNeighbors() const {return m_NbofNeighbors;}
