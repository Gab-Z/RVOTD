#ifndef DEF_ETP_Triangle
#define DEF_ETP_Triangle

#include "./Edge.h"
#include "./Component.h"
#include <vector>
#include <memory>

namespace ETP {



template <typename P, typename D>

class Triangle {

  std::shared_ptr<ETP::Edge<P,D>> edge1 = nullptr;
  std::shared_ptr<ETP::Edge<P,D>> edge2 = nullptr;
  std::shared_ptr<ETP::Edge<P,D>> edge3 = nullptr;
  int level = -1;
  int component = -1;
  std::shared_ptr<ETP::Triangle<P,D>> adjacents[ 3 ] = { nullptr, nullptr, nullptr };
  std::shared_ptr<ETP::Triangle<P,D>> connectedNodes[ 3 ] = { nullptr, nullptr, nullptr };
  D chokes[ 3 ]  = { (D) 0.0f, (D) 0.0f, (D) 0.0f };
  D widths[ 3 ]  = { (D) 0.0f, (D) 0.0f, (D) 0.0f };
  D lowerBoundDistances[ 3 ]  = { (D) 0.0f, (D) 0.0f, (D) 0.0f };
  D angles[ 3 ]  = { (D) 0.0f, (D) 0.0f, (D) 0.0f };
  unsigned int id ;

public:
  Triangle(){};
  Triangle( std::shared_ptr<ETP::Edge<P,D>> _edge1, std::shared_ptr<ETP::Edge<P,D>> _edge2, std::shared_ptr<ETP::Edge<P,D>> _edge3 ){
    orderEdges( _edge1, _edge2, _edge3 );
  };
  Triangle( std::shared_ptr<ETP::Edge<P,D>> _edge1, std::shared_ptr<ETP::Edge<P,D>> _edge2, std::shared_ptr<ETP::Edge<P,D>> _edge3, unsigned int _id ): id( _id ){
    orderEdges( _edge1, _edge2, _edge3 );
  };
  ~Triangle(){};
  void orderEdges( std::shared_ptr<ETP::Edge<P,D>> _edge1, std::shared_ptr<ETP::Edge<P,D>> _edge2, std::shared_ptr<ETP::Edge<P,D>> _edge3 ){
    edge1 = _edge1;
    std::shared_ptr<ETP::Point<P,D>> e1pt2 = edge1->getPointVal( 1 );
    if( _edge2->getPointVal( 0 )->equals( e1pt2 ) || _edge2->getPointVal( 1 )->equals( e1pt2 ) ){
      edge2 = _edge2;
      edge3 = _edge3;
    }else{
      edge3 = _edge2;
      edge2 = _edge3;
    }
  }
  void setComponent( int _c ){
    component = _c;
  }
  int getComponent(){
    return component;
  }

  void setWidth( int _i, D _width ){
    widths[ _i ] = _width;
  }
  D getWidth( int _i ){
    return widths[ _i ];
  }
  D getWidthbetweenEdges( std::shared_ptr<ETP::Edge<P,D>> _edge1, std::shared_ptr<ETP::Edge<P,D>> _edge2 ){
    if( _edge1 == nullptr || _edge2 == nullptr ){ return 0.0f; }
    if( (_edge1 == edge1 && _edge2 == edge2) || (_edge1 == edge2 && _edge2 == edge1) ){
      return widths[ 0 ];
    }else if( (_edge1 == edge2 && _edge2 == edge3) || (_edge1 == edge3 && _edge2 == edge2) ){
      return widths[ 1 ];
    }else{
      return widths[ 2 ];
    }
  }
  void setAngle( size_t _i, D _angle ){
    angles[ _i ] = _angle;
  }
  D getAngle( int _i ){
    return angles[ _i ];
  }
  void setChoke( size_t _i, P _choke ){
    chokes[ _i ] = _choke;
  }
  D getChoke( int _i ){
    return chokes[ _i ];
  }
  void setLowerBound( size_t _i, D _angle ){
    lowerBoundDistances[ _i ] = _angle;
  }
  D getLowerBound( size_t _i ){
    return lowerBoundDistances[ _i ];
  }
  void setLevel( int _l ){
    level = _l;
  }
  int getLevel(){
    return level;
  }
  std::shared_ptr<ETP::Edge<P,D>>& getEdgeRef( int _i ){
    if( _i == 0 ){
      return edge1;
    }else if( _i == 1 ){
      return edge2;
    }else{
      return edge3;
    }
  }
  std::shared_ptr<ETP::Edge<P,D>> getEdgeVal( int _i ){
    if( _i == 0 ){
      return edge1;
    }else if( _i == 1 ){
      return edge2;
    }else{
      return edge3;
    }
  }
  void setAdjacent( size_t _i, std::shared_ptr<ETP::Triangle<P,D>> adjacentTri ){
    adjacents[ _i ] = adjacentTri;
  }
  std::shared_ptr<ETP::Triangle<P,D>>& getAdjacentRef( size_t _i ){
    return adjacents[ _i ];
  }
  std::shared_ptr<ETP::Triangle<P,D>> getAdjacentVal( size_t _i ){
    return adjacents[ _i ];
  }
  void setConnectedNode( size_t _i, std::shared_ptr<ETP::Triangle<P,D>> _tri ){
    connectedNodes[ _i ] = _tri;
  }
  std::shared_ptr<ETP::Triangle<P,D>> getConnectedNode( size_t _i ){
    return connectedNodes[ _i ];
  }
  unsigned int getId(){
    return id;
  }
  std::vector<P> getPositions(){
    std::vector<P> ret;
    std::shared_ptr<ETP::Point<P,D>> pt1 = edge1->getPointVal( 0 );
    std::shared_ptr<ETP::Point<P,D>> pt2 = edge1->getPointVal( 1 );
    std::shared_ptr<ETP::Point<P,D>> e2pt1 = edge2->getPointVal( 0 );
    std::shared_ptr<ETP::Point<P,D>> pt3;
    if( !e2pt1->equals( pt1 ) && !e2pt1->equals( pt2 ) ){
      pt3 = e2pt1;
    }else{
      pt3 = edge2->getPointVal( 1 );
    }
    ret.push_back( pt1->x );
    ret.push_back( pt1->y );
    ret.push_back( pt2->x );
    ret.push_back( pt2->y );
    ret.push_back( pt3->x );
    ret.push_back( pt3->y );
    return ret;
  }
  std::shared_ptr<ETP::Triangle<P,D>> triangleOpposite( std::shared_ptr<ETP::Edge<P,D>> _edge ){
    if( _edge->hasSamePoints( edge1 ) ){
      return adjacents[ 0 ];
    }else if( _edge->hasSamePoints( edge2 ) ){
      return adjacents[ 1 ];
    }else if( _edge->hasSamePoints( edge3 ) ){
      return adjacents[ 2 ];
    }
  }
  std::vector<std::shared_ptr<ETP::Edge<P,D>>> otherEdges( std::shared_ptr<ETP::Edge<P,D>> _edge ){
    std::vector<std::shared_ptr<ETP::Edge<P,D>>> ret( 2 );
    if( _edge == edge1 ){
      ret[ 0 ] = edge2;
      ret[ 1 ] = edge3;
    }else if( _edge == edge2 ){
      ret[ 0 ] = edge3;
      ret[ 1 ] = edge1;
    }else if( _edge == edge3 ){
      ret[ 0 ] = edge1;
      ret[ 1 ] = edge2;
    }
    return ret;
  }
  std::shared_ptr<ETP::Edge<P,D>> edgeOpposite( std::shared_ptr<ETP::Point<P,D>> _pt ){
    if( !edge1->hasPoint( _pt ) ){
      return edge1;
    }else if( !edge2->hasPoint( _pt ) ){
      return edge2;
    }else{
      return edge3;
    }
  }
  std::shared_ptr<ETP::Point<P,D>> vertexOpposite( std::shared_ptr<ETP::Edge<P,D>> _edge ){
    if( !_edge->hasPoint( edge1->getPointVal( 0 ) ) ){
      return edge1->getPointVal( 0 );
    }else if( !_edge->hasPoint( edge1->getPointVal( 1 ) ) ){
      return edge1->getPointVal( 1 );
    }else if( !_edge->hasPoint( edge2->getPointVal( 0 ) ) ){
      return edge2->getPointVal( 0 );
    }else{
      return edge2->getPointVal( 1 );
    }
  }
  int getNumConstrainedEdges(){
    int n = 0;
    /*
    if(  edge1->isConstrained() == true ){ n++; }
    if(  edge2->isConstrained() == true ){ n++; }
    if(  edge3->isConstrained() == true ){ n++; }
    */
    if(  adjacents[ 0 ] == nullptr ){ n++; }
    if(  adjacents[ 1 ] == nullptr ){ n++; }
    if(  adjacents[ 2 ] == nullptr ){ n++; }
    return n;
  }
  int getNumAdjacentLevel( int _lvl ){
    int n = 0;
    if(  adjacents[ 0 ] != nullptr && adjacents[ 0 ]->getLevel() == _lvl ){ n++; }
    if(  adjacents[ 1 ] != nullptr && adjacents[ 1 ]->getLevel() == _lvl ){ n++; }
    if(  adjacents[ 2 ] != nullptr && adjacents[ 2 ]->getLevel() == _lvl ){ n++; }
    return n;
  }
  bool hasSamePoints( std::shared_ptr<ETP::Triangle<P,D>> _tri ){
    std::shared_ptr<ETP::Edge<P,D>> tEdge1 = _tri->getEdgeVal( 0 );
    std::shared_ptr<ETP::Edge<P,D>> tEdge2 = _tri->getEdgeVal( 1 );
    std::shared_ptr<ETP::Edge<P,D>> tEdge3 = _tri->getEdgeVal( 2 );
    if(     ( edge1->hasSamePoints( tEdge1 ) && edge2->hasSamePoints( tEdge2 ) && edge3->hasSamePoints( tEdge3 ) )
        ||  ( edge1->hasSamePoints( tEdge2 ) && edge2->hasSamePoints( tEdge3 ) && edge3->hasSamePoints( tEdge1 ) )
        ||  ( edge1->hasSamePoints( tEdge3 ) && edge2->hasSamePoints( tEdge1 ) && edge3->hasSamePoints( tEdge2 ) )
        ||  ( edge1->hasSamePoints( tEdge1 ) && edge2->hasSamePoints( tEdge3 ) && edge3->hasSamePoints( tEdge2 ) )
        ||  ( edge1->hasSamePoints( tEdge2 ) && edge2->hasSamePoints( tEdge1 ) && edge3->hasSamePoints( tEdge3 ) )
        ||  ( edge1->hasSamePoints( tEdge3 ) && edge2->hasSamePoints( tEdge2 ) && edge3->hasSamePoints( tEdge1 ) )
    ){ return true; }
  }
  std::shared_ptr<ETP::Edge<P,D>> getSharedEdgeWith( std::shared_ptr<ETP::Triangle<P,D>> _tri ){
    for( int i = 0; i < 3; i++ ){
      std::shared_ptr<ETP::Edge<P,D>> edge = getEdgeVal( i );
      for( int j = 0; j < 3; j++ ){
        if( _tri->getEdgeVal( j ) == edge ){
          return edge;
        }
      }
    }
    return nullptr;
  }
  int getAdjacentIndex( std::shared_ptr<ETP::Triangle<P,D>> _tri ){
    if( _tri == adjacents[ 0 ] ){
      return 0;
    }else if( _tri == adjacents[ 1 ] ){
      return 1;
    }else if( _tri == adjacents[ 2 ] ){
      return 2;
    }
    return -1;
  }
};

}

#endif
