#ifndef DEF_ETP_Triangle
#define DEF_ETP_Triangle

#include "./Edge.h"
#include "./Component.h"
#include <vector>
#include <memory>

namespace ETP {

template <typename T>

class Triangle {

  std::shared_ptr<ETP::Edge<T>> edge1 = nullptr;
  std::shared_ptr<ETP::Edge<T>> edge2 = nullptr;
  std::shared_ptr<ETP::Edge<T>> edge3 = nullptr;
  enum level { ZERO, ONE, TWO, THREE };
  std::shared_ptr<ETP::Component> component = nullptr;
  std::shared_ptr<ETP::Triangle<T>> adjacents[ 3 ] = { nullptr, nullptr, nullptr };
  T chokes[ 3 ]  = { (T) 0.0f, (T) 0.0f, (T) 0.0f };
  T widths[ 3 ]  = { (T) 0.0f, (T) 0.0f, (T) 0.0f };
  T lowerBoundDistances[ 3 ]  = { (T) 0.0f, (T) 0.0f, (T) 0.0f };
  T angles[ 3 ]  = { (T) 0.0f, (T) 0.0f, (T) 0.0f };
  unsigned int id;

public:
  Triangle(){};
  Triangle( std::shared_ptr<ETP::Edge<T>> _edge1, std::shared_ptr<ETP::Edge<T>> _edge2, std::shared_ptr<ETP::Edge<T>> _edge3 ){
    orderEdges( _edge1, _edge2, _edge3 );
  };
  Triangle( std::shared_ptr<ETP::Edge<T>> _edge1, std::shared_ptr<ETP::Edge<T>> _edge2, std::shared_ptr<ETP::Edge<T>> _edge3, unsigned int _id ): id( _id ){
    orderEdges( _edge1, _edge2, _edge3 );
  };
  ~Triangle(){};

  void orderEdges( std::shared_ptr<ETP::Edge<T>> _edge1, std::shared_ptr<ETP::Edge<T>> _edge2, std::shared_ptr<ETP::Edge<T>> _edge3 ){
    edge1 = _edge1;
    std::shared_ptr<ETP::Point<T>> e1pt2 = edge1->getPointVal( 1 );
    if( _edge2->getPointVal( 0 )->equals( e1pt2 ) || _edge2->getPointVal( 1 )->equals( e1pt2 ) ){
      edge2 = _edge2;
      edge3 = _edge3;
    }else{
      edge3 = _edge2;
      edge2 = _edge3;
    }
  }
  void setComponent( std::shared_ptr<Component> _c ){
    component = _c;
  }

  void setWidth( int _i, T _width ){
    widths[ _i ] = _width;
  }
  T getWidth( int _i ){
    return widths[ _i ];
  }
  void setAngle( size_t _i, T _angle ){
    angles[ _i ] = _angle;
  }
  T getAngle( int _i ){
    return angles[ _i ];
  }
  void setChoke( size_t _i, T _choke ){
    chokes[ _i ] = _choke;
  }
  std::shared_ptr<ETP::Edge<T>>& getEdgeRef( int _i ){
    if( _i == 0 ){
      return edge1;
    }else if( _i == 1 ){
      return edge2;
    }else{
      return edge3;
    }
  }
  std::shared_ptr<ETP::Edge<T>> getEdgeVal( int _i ){
    if( _i == 0 ){
      return edge1;
    }else if( _i == 1 ){
      return edge2;
    }else{
      return edge3;
    }
  }
  void setAdjacent( size_t _i, std::shared_ptr<ETP::Triangle<T>> adjacentTri ){
    adjacents[ _i ] = adjacentTri;
  }
  std::shared_ptr<ETP::Triangle<T>>& getAdjacentRef( size_t _i ){
    return adjacents[ _i ];
  }
  unsigned int getId(){
    return id;
  }
  std::vector<T> getPositions(){
    std::vector<T> ret;
    std::shared_ptr<ETP::Point<T>> pt1 = edge1->getPointVal( 0 );
    std::shared_ptr<ETP::Point<T>> pt2 = edge1->getPointVal( 1 );
    std::shared_ptr<ETP::Point<T>> e2pt1 = edge2->getPointVal( 0 );
    std::shared_ptr<ETP::Point<T>> pt3;
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
  std::shared_ptr<ETP::Triangle<T>> triangleOpposite( std::shared_ptr<ETP::Edge<T>> _edge ){
    if( _edge == edge1 ){
      return adjacents[ 0 ];
    }else if( _edge == edge2 ){
      return adjacents[ 1 ];
    }else if( _edge == edge3 ){
      return adjacents[ 2 ];
    }
  }
  std::vector<std::shared_ptr<ETP::Edge<T>>> otherEdges( std::shared_ptr<ETP::Edge<T>> _edge ){
    std::vector<std::shared_ptr<ETP::Edge<T>>> ret( 2 );
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
  std::shared_ptr<ETP::Edge<T>> edgeOpposite( std::shared_ptr<ETP::Point<T>> _pt ){
    if( !edge1->hasPoint( _pt ) ){
      return edge1;
    }else if( !edge2->hasPoint( _pt ) ){
      return edge2;
    }else{
      return edge3;
    }
  }
  std::shared_ptr<ETP::Point<T>> vertexOpposite( std::shared_ptr<ETP::Edge<T>> _edge ){
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

};

}

#endif
