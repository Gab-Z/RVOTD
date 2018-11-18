#ifndef DEF_ETP_EDGE
#define DEF_ETP_EDGE

#include <iostream>
#include <vector>
#include <memory>
#include <cmath>

#include "./Point.h"

namespace ETP {

template <typename P,typename D>

class Edge {

  std::shared_ptr<ETP::Point<P,D>> point1;
  std::shared_ptr<ETP::Point<P,D>> point2;
  bool constrained;

public:
  Edge(){};
  Edge( std::shared_ptr<ETP::Point<P,D>> _p0, std::shared_ptr<ETP::Point<P,D>> _p1 ):point1( _p0 ), point2( _p1 ){};
  Edge( std::shared_ptr<ETP::Point<P,D>> _p0, std::shared_ptr<ETP::Point<P,D>> _p1, bool _constrained ){
    point1 = _p0,
    point2 = _p1,
    constrained = _constrained;
  };
  Edge( ETP::Point<P,D> _p0, ETP::Point<P,D> _p1, bool _constrained ){
    point1 = std::make_shared<ETP::Point<P,D>>( _p0 ),
    point2 = std::make_shared<ETP::Point<P,D>>( _p1 ),
    constrained = _constrained;
  };
  ~Edge(){};

  bool isConstrained(){
    return constrained;
  }

  void setContrained( bool _constrained ){
    constrained = _constrained;
  }

  std::shared_ptr<ETP::Point<P,D>>& getPointRef( const int& _i ){
    if( _i == 0 ){
      std::shared_ptr<Point<P,D>>& refPoint1 = point1;
      return refPoint1;
    }else{
      std::shared_ptr<Point<P,D>>& refPoint2 = point2;
      return refPoint2;
    }
  }

  std::shared_ptr<ETP::Point<P,D>> getPointVal( const int& _i ){
    if( _i == 0 ){
      return point1;
    }else{
      return point2;
    }
  }

  void setPoint( const size_t& _i, std::shared_ptr<ETP::Point<P,D>> _pt ){
    if( _i == 0 ){
      point1 = _pt;
    }else if( _i == 1 ){
      point2 = _pt;
    }
  }

  bool hasSamePoints( const std::shared_ptr<ETP::Point<P,D>>& _pt1, const std::shared_ptr<ETP::Point<P,D>>& _pt2 ){
    if( ( point1->equals( _pt1 ) && point2->equals( _pt2 ) ) || ( point1->equals( _pt2 ) && point2->equals( _pt1 ) ) ){
      return true;
    }
    return false;
  }

  bool hasSamePoints( const std::shared_ptr<ETP::Edge<P,D>>& _e ){
    std::shared_ptr<Point<P,D>>& e_pt1 = _e->getPointRef( 0 );
    std::shared_ptr<Point<P,D>>& e_pt2 = _e->getPointRef( 1 );
    if( ( point1->equals( e_pt1 ) && point2->equals( e_pt2 ) ) || ( point1->equals( e_pt2 ) && point2->equals( e_pt1 ) ) ){
      return true;
    }
    return false;
  }
  /*
  bool hasSamePoints( std::shared_ptr<ETP::Edge<P,D>> _e ){
    std::shared_ptr<Point<P,D>>& e_pt1 = _e->getPointRef( 0 );
    std::shared_ptr<Point<P,D>>& e_pt2 = _e->getPointRef( 1 );
    if( ( point1->equals( e_pt1 ) && point2->equals( e_pt2 ) ) || ( point1->equals( e_pt2 ) && point2->equals( e_pt1 ) ) ){
      return true;
    }
    return false;
  }
  */

  bool hasSamePoints( const P& _p1x, const P& _p1y, const P& _p2x, const P& _p2y ){
    ETP::Point<P,D>& e_pt1 = getPointRef( 0 );
    ETP::Point<P,D>& e_pt2 = getPointRef( 1 );
    if( ( point1->equals( _p1x, _p1y ) && point2->equals( _p2x, _p2y ) ) || ( point1->equals( _p2x, _p2y ) && point2->equals( _p1x, _p1y ) ) ){
      return true;
    }
    return false;
  }

  std::vector<P> getPositions(){
      std::vector<P> ret;
      ret.push_back( point1->x );
      ret.push_back( point1->y );
      ret.push_back( point2->x );
      ret.push_back( point2->y );
      return ret;
  }

  bool hasPoint( const std::shared_ptr<ETP::Point<P,D>>& _pt ){
    if( point1->equals( _pt ) || point2->equals( _pt ) ){
      return true;
    }
    return false;
  }

  D getLength(){
    return std::sqrt( (D) ( std::pow( point2->x - point1->x, 2 ) + std::pow( point2->y - point1->y, 2 ) ) );
  }





};

}

#endif
