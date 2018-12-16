#ifndef DEF_ETP_Triangle
#define DEF_ETP_Triangle

#include "./Edge.h"
#include "./Point.h"
#include "./SearchNode.h"
#include <vector>
#include <memory>

namespace ETP {

template <typename sP = int, typename sD = float>
class SearchNode;

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
  ETP::Point<D,D> centroid;
  std::shared_ptr<ETP::SearchNode<P,D>> searchNode = nullptr;

public:
  Triangle(){};
  Triangle( std::shared_ptr<ETP::Edge<P,D>> _edge1, std::shared_ptr<ETP::Edge<P,D>> _edge2, std::shared_ptr<ETP::Edge<P,D>> _edge3 ){
    init( std::move( _edge1 ), std::move( _edge2 ), std::move( _edge3 ) );
  };
  Triangle( std::shared_ptr<ETP::Edge<P,D>> _edge1, std::shared_ptr<ETP::Edge<P,D>> _edge2, std::shared_ptr<ETP::Edge<P,D>> _edge3, unsigned int _id ): id( _id ){
    init( std::move( _edge1 ), std::move( _edge2 ), std::move( _edge3 ) );
  };
  ~Triangle(){};
  void init( std::shared_ptr<ETP::Edge<P,D>> _edge1, std::shared_ptr<ETP::Edge<P,D>> _edge2, std::shared_ptr<ETP::Edge<P,D>> _edge3 ){
    orderEdges( std::move( _edge1 ), std::move( _edge2 ), std::move( _edge3 ) );
    centroid = calculateCentroid();
  }
  void orderEdges( std::shared_ptr<ETP::Edge<P,D>> _edge1, std::shared_ptr<ETP::Edge<P,D>> _edge2, std::shared_ptr<ETP::Edge<P,D>> _edge3 ){
    edge1 = std::move( _edge1 );
    std::shared_ptr<ETP::Point<P,D>> e1pt2 = edge1->getPointVal( 1 );
    if( _edge2->getPointVal( 0 )->equals( e1pt2 ) || _edge2->getPointVal( 1 )->equals( e1pt2 ) ){
      edge2 = std::move( _edge2 );
      edge3 = std::move( _edge3 );
    }else{
      edge3 = std::move( _edge2 );
      edge2 = std::move( _edge3 );
    }
  }
  void setComponent( int _c ){
    component = _c;
  }
  int getComponent(){
    return component;
  }
  std::shared_ptr<ETP::SearchNode<P,D>> getSearchNode(){
    return searchNode;
  }
  void setSearchNode( ETP::SearchNode<P,D> _searchNode ){
    searchNode = std::make_shared<ETP::SearchNode<P,D>>( _searchNode );
  }
  void clearSearchNode(){
    searchNode = nullptr;
  }
  void setSearchHValue( D _hValue ){
    if( searchNode != nullptr ){ searchNode->hValue = _hValue; }
  }
  void setSearchGValue( D _gValue ){
    if( searchNode != nullptr ){ searchNode->gValue = _gValue; }
  }
  void setSearchComeFrom( std::shared_ptr<ETP::Triangle<P,D>> _comeFrom ){
    if( searchNode != nullptr ){ searchNode->comeFrom = _comeFrom; }
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
  D getAngle( std::shared_ptr<ETP::Edge<P,D>> _edge1, std::shared_ptr<ETP::Edge<P,D>> _edge2 ){
    if( ( _edge1 == edge1 && _edge2 == edge2 ) || ( _edge1 == edge2 && _edge2 == edge1 ) ){
      return angles[ 0 ];
    }else if( ( _edge1 == edge2 && _edge2 == edge3 ) || ( _edge1 == edge3 && _edge2 == edge2 ) ){
      return angles[ 1 ];
    }else if( ( _edge1 == edge3 && _edge2 == edge1 ) || ( _edge1 == edge1 && _edge2 == edge3 ) ){
      return angles[2];
    }
    //return 0.0f;
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
  ETP::Point<D,D> getCentroid(){
    return centroid;
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
  int getEdgeIndex( std::shared_ptr<ETP::Edge<P,D>> _edge ){
    if( _edge == edge1 ){ return 0; }
    else if( _edge == edge2 ){ return 1; }
    else if( _edge == edge3 ){ return 2; }
    else{ return -1; }
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
  bool isPointInside( std::shared_ptr<ETP::Point<P,D>> _pt ){
    std::vector<std::shared_ptr<ETP::Point<P,D>>> points = getPoints();
    std::shared_ptr<ETP::Point<P,D>> p0 = points[ 0 ];
    std::shared_ptr<ETP::Point<P,D>> p1 = points[ 1 ];
    std::shared_ptr<ETP::Point<P,D>> p2 = points[ 2 ];
    D A = 0.5 *  (-p1->y * p2->x + p0->y * (-p1->x + p2->x) + p0->x * (p1->y - p2->y) + p1->x * p2->y);
    P sign = A < 0.0f ? (P) -1 : (P) 1;
    P s = (p0->y * p2->x - p0->x * p2->y + (p2->y - p0->y) * _pt->x + (p0->x - p2->x) * _pt->y) * sign;
    P t = (p0->x * p1->y - p0->y * p1->x + (p0->y - p1->y) * _pt->x + (p1->x - p0->x) * _pt->y) * sign;
    return  s > 0 && t > 0 &&  (D) ( s + t ) <  (D) 2.0f * A *  (D) sign;
  }
  bool isPointInside( ETP::Point<P,D> _pt ){
    std::vector<std::shared_ptr<ETP::Point<P,D>>> points = getPoints();
    std::shared_ptr<ETP::Point<P,D>> p0 = points[ 0 ];
    std::shared_ptr<ETP::Point<P,D>> p1 = points[ 1 ];
    std::shared_ptr<ETP::Point<P,D>> p2 = points[ 2 ];
    D A = 0.5 *  (-p1->y * p2->x + p0->y * (-p1->x + p2->x) + p0->x * (p1->y - p2->y) + p1->x * p2->y);
    P sign = A < 0.0f ? (P) -1 : (P) 1;
    P s = (p0->y * p2->x - p0->x * p2->y + (p2->y - p0->y) * _pt.x + (p0->x - p2->x) * _pt.y) * sign;
    P t = (p0->x * p1->y - p0->y * p1->x + (p0->y - p1->y) * _pt.x + (p1->x - p0->x) * _pt.y) * sign;
    return  s > 0 && t > 0 &&  (D) ( s + t ) <  (D) 2.0f * A *  (D) sign;
  }
  bool isPointInside( ETP::Point<P,D> p_pt, D _scale ){
    std::vector<std::shared_ptr<ETP::Point<P,D>>> points = getPoints();
    ETP::Point<P,D> bp0 = *points[ 0 ].get();
    ETP::Point<D,D> p0 = ETP::Point<D,D>( (D) bp0.x * _scale, (D) bp0.y * _scale );
    ETP::Point<P,D> bp1 = *points[ 1 ].get();
    ETP::Point<D,D> p1 = ETP::Point<D,D>( (D) bp1.x * _scale, (D) bp1.y * _scale );
    ETP::Point<P,D> bp2 = *points[ 2 ].get();
    ETP::Point<D,D> p2 = ETP::Point<D,D>( (D) bp2.x * _scale, (D) bp2.y * _scale );
    ETP::Point<D,D> _pt = ETP::Point<D,D>( (D) p_pt.x, (D) p_pt.y );
    D A = 0.5 *  (-p1.y * p2.x + p0.y * (-p1.x + p2.x) + p0.x * (p1.y - p2.y) + p1.x * p2.y );
    D sign = A < 0.0f ? (D) -1.0 : (D) 1.0;
    D s = (p0.y * p2.x - p0.x * p2.y + (p2.y - p0.y) * _pt.x + (p0.x - p2.x) * _pt.y) * sign;
    D t = (p0.x * p1.y - p0.y * p1.x + (p0.y - p1.y) * _pt.x + (p1.x - p0.x) * _pt.y) * sign;
    return  s > 0.0 && t > 0.0 && ( s + t ) < 2.0 * A * sign;
  }
  std::vector<ETP::Point<P,D>> getBoundingBox(){
    std::vector<std::shared_ptr<ETP::Point<P,D>>> pts = getPoints();
    std::shared_ptr<ETP::Point<P,D>> p1 = pts[ 0 ];
    std::shared_ptr<ETP::Point<P,D>> p2 = pts[ 1 ];
    std::shared_ptr<ETP::Point<P,D>> p3 = pts[ 2 ];
    ETP::Point<P,D> min( p1->x, p1->y );
    min.x = std::min( min.x, p2->x );
    min.x = std::min( min.x, p3->x );
    min.y = std::min( min.y, p2->y );
    min.y = std::min( min.y, p3->y );
    ETP::Point<P,D> max( p1->x, p1->y );
    max.x = std::max( max.x, p2->x );
    max.x = std::max( max.x, p3->x );
    max.y = std::max( max.y, p2->y );
    max.y = std::max( max.y, p3->y );
    std::vector<ETP::Point<P,D>> ret = std::vector<ETP::Point<P,D>>( 2 );
    ret[ 0 ] = min;
    ret[ 1 ] = max;
    return ret;
  }
  std::vector<std::shared_ptr<ETP::Point<P,D>>> getPoints(){
    std::vector<std::shared_ptr<ETP::Point<P,D>>> ret;
    ret.push_back( edge1->getPointVal( 0 ) );
    ret.push_back( edge1->getPointVal( 1 ) );
    if( edge2->getPointVal( 0 ) == edge1->getPointVal( 0 ) || edge2->getPointVal( 0 ) == edge1->getPointVal( 1 ) ){
      ret.push_back( edge2->getPointVal( 1 ) );
    }else{
      ret.push_back( edge2->getPointVal( 0 ) );
    }
    return ret;
  }
  ETP::Point<D,D> calculateCentroid(){
    std::vector<std::shared_ptr<ETP::Point<P,D>>> pts = getPoints();
    return ETP::Point<D,D>( (D) pts[ 0 ]->x / 3.0 + (D) pts[ 1 ]->x / 3.0 + (D) pts[ 2 ]->x / 3.0,
                            (D) pts[ 0 ]->y / 3.0 + (D) pts[ 1 ]->y / 3.0 + (D) pts[ 2 ]->y / 3.0
   );
  }
  int getConnectedNodeIndex( std::shared_ptr<ETP::Triangle<P,D>> _other ){
    if( _other == nullptr ){ return -1; }
    for ( int i = 0; i < 3; i++ ){
      if( connectedNodes[ i ] == _other ){ return i; }
    }
    return -1;
  }
  int getReversedConnectedNodeIndex( std::shared_ptr<ETP::Triangle<P,D>> _connectedNode ){
    std::shared_ptr<ETP::Triangle<P,D>> otherConnectedNode;
    D lowerBound = getLowerBound( getConnectedNodeIndex(_connectedNode ) );
    std::shared_ptr<ETP::Edge<P,D>> e1 = getEdgeVal( getConnectedNodeIndex( _connectedNode ) );
    for( int k = 0; k < 3; k++ ){
      std::shared_ptr<ETP::Triangle<P,D>> kNode = getConnectedNode( k );
      if( kNode != nullptr && kNode != _connectedNode ){
        otherConnectedNode = kNode;
        lowerBound += getLowerBound( k );
        break;
      }
    }
    for( int l = 0; l < 3; l++ ){
      std::shared_ptr<ETP::Triangle<P,D>> lNode = _connectedNode->getConnectedNode( l );
      if( lNode != nullptr && lNode == otherConnectedNode /*&& _connectedNode->getLowerBound( l ) == lowerBound*/ ){
        std::shared_ptr<ETP::Edge<P,D>> e2 = getEdgeVal( getConnectedNodeIndex( otherConnectedNode ) );
        D innerAngle = getAngle( e1, e2 );
        if( _connectedNode->getLowerBound( l ) == lowerBound + innerAngle ){
          return l;
        }
      }
    }
  }
  int getLvl2ToLvl2Index( std::shared_ptr<ETP::Triangle<P,D>> _otherLvl2 ){
    for( int i = 0; i < 3; i++ ){
      std::shared_ptr<ETP::Triangle<P,D>> iNode = getConnectedNode( i );
      if( iNode != nullptr ){
        D iLowerBound = getLowerBound( i );
        for( int j = 0; j < 3; j++ ){
          std::shared_ptr<ETP::Triangle<P,D>> jNode = _otherLvl2->getConnectedNode( j );
          if( jNode != nullptr && jNode == iNode && _otherLvl2->getLowerBound( j ) < iLowerBound ){
            return i;
          }
        }
      }
    }
    return -1;
  }
  /*
  std::vector<std::shared_ptr<ETP::Triangle<P,D>>> getTrianglesToConnectedNode( int connectedNodeIdx ){
    std::vector<std::shared_ptr<ETP::Triangle<P,D>>> ret;
    std::shared_ptr<ETP::Triangle<P,D>> goal = connectedNodes[ connectedNodeIdx ];
    if( goal == nullptr ){ throw std::invalid_argument( "Triangle::getTrianglesToConnectedNode: search to null goal" ); }
    std::shared_ptr<ETP::Triangle<P,D>> tCurrent = this;
    bool endReached = false;
    while( endReached == false ){
      ret.push_back( tCurrent );
      if( tCurrent == goal ){
        endReached = true;
        break;
      }
      int nextIndex = tCurrent->getConnectedNodeIndex( goal );
      tCurrent = tCurrent->getAdjacentVal( nextIndex );
    }
    return ret;
  }
  */

  std::shared_ptr<ETP::Triangle<P,D>> getLvl1RootNode(){
    if( level != 1 ){ return nullptr; }
    for( size_t i = 0; i < 3; i++ ){
      if( connectedNodes[ i ] != nullptr ){
        return connectedNodes[ i ];
      }
    }
    return nullptr;
  }
  bool isPartOfLvl2Loop(){
    if( level != 2 ){ return false; }
    int ct = 0;
    std::shared_ptr<ETP::Triangle<P,D>> endPoint1 = nullptr;
    for ( int i = 0; i < 3; i++ ){
      if( endPoint1 == nullptr && connectedNodes[ i ] != nullptr ){
        endPoint1 = connectedNodes[ i ];
      }else if( endPoint1 != nullptr && endPoint1 == connectedNodes[ i ] ){
        return true;
      }
    }
    return false;
  }
  bool isPartoflvl2Ring(){
    if( level != 2 || connectedNodes[ 0 ] != nullptr || connectedNodes[ 1 ] != nullptr || connectedNodes[ 2 ] != nullptr ){ return false; }
    return true;
  }
  bool isOnSameLvl2Loop( std::shared_ptr<ETP::Triangle<P,D>> _other ){
    if( level != 2 || _other->getLevel() != 2 ){ return false; }
    std::shared_ptr<ETP::Triangle<P,D>> endPoint1 = nullptr;
    std::shared_ptr<ETP::Triangle<P,D>> endPoint2 = nullptr;
    for( int i = 0; i < 3; i++ ){
      std::shared_ptr<ETP::Triangle<P,D>> e1 = connectedNodes[ i ];
      if( endPoint1 == nullptr && e1 != nullptr ){
        endPoint1 = e1;
      }else if( endPoint1 != nullptr && ( e1 != nullptr && e1 != endPoint1 ) ){
        return false;
      }
      std::shared_ptr<ETP::Triangle<P,D>> e2 = _other->getConnectedNode( i );
      if( endPoint2 == nullptr && e2 != nullptr ){
        endPoint2 = e2;
      }else if( endPoint2 != nullptr && ( e2 != nullptr && e2 != endPoint2 ) ){
        return false;
      }
    }
    if( endPoint1 == endPoint2 ){
      return true;
    }
    return false;
  }
  bool haveSameLvl2CorridorEndpoints( std::shared_ptr<ETP::Triangle<P,D>> _other ){
    if( level != 2 || _other->getLevel() != 2 ){ return false; }
    std::shared_ptr<ETP::Triangle<P,D>> endPoint1 = nullptr;
    std::shared_ptr<ETP::Triangle<P,D>> endPoint2 = nullptr;
    std::shared_ptr<ETP::Triangle<P,D>> o_endPoint1 = nullptr;
    std::shared_ptr<ETP::Triangle<P,D>> o_endPoint2 = nullptr;
    for( int i = 0; i < 3; i++ ){
      std::shared_ptr<ETP::Triangle<P,D>> cn = connectedNodes[ i ];
      if( cn != nullptr ){
        if( endPoint1 == nullptr ){ endPoint1 = cn; }
        else if( cn == endPoint1 ){ return false; }
        else{ endPoint2 = cn; }
      }
      std::shared_ptr<ETP::Triangle<P,D>> cn2 = _other->getConnectedNode( i );
      if( cn2 != nullptr ){
        if( o_endPoint1 == nullptr ){ o_endPoint1 = cn2; }
        else if( cn2 == o_endPoint1 ){ return false; }
        else{ o_endPoint2 = cn2; }
      }
    }
    if( ( endPoint1 == o_endPoint1 && endPoint2 == o_endPoint2 ) || ( endPoint1 == o_endPoint2 && endPoint2 == o_endPoint1 ) ){
      return true;
    }
    return false;
  }

};// end of Class





}// end of namespace

#endif
