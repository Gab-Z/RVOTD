#ifndef DEF_ETP_EDGE
#define DEF_ETP_EDGE

#include <iostream>
#include <vector>
#include <memory>
#include <cmath>

#include "./Point.h"

namespace ETP {

template <typename T>

class Edge {

  std::shared_ptr<ETP::Point<T>> point1;
  std::shared_ptr<ETP::Point<T>> point2;
  bool constrained;

public:
  Edge(){};
  Edge( std::shared_ptr<ETP::Point<T>> _p0, std::shared_ptr<ETP::Point<T>> _p1 ):point1( _p0 ), point2( _p1 ){};
  Edge( std::shared_ptr<ETP::Point<T>> _p0, std::shared_ptr<ETP::Point<T>> _p1, bool _constrained ){
    point1 = _p0,
    point2 = _p1,
    constrained = _constrained;
  };
  Edge( ETP::Point<T> _p0, ETP::Point<T> _p1, bool _constrained ){
    point1 = std::make_shared<ETP::Point<T>>( _p0 ),
    point2 = std::make_shared<ETP::Point<T>>( _p1 ),
    constrained = _constrained;
  };
  ~Edge(){};

  bool isConstrained(){
    return constrained;
  }

  void setContrained( bool _constrained ){
    constrained = _constrained;
  }

  std::shared_ptr<ETP::Point<T>>& getPointRef( const int& _i ){
    if( _i == 0 ){
      std::shared_ptr<Point<T>>& refPoint1 = point1;
      return refPoint1;
    }else{
      std::shared_ptr<Point<T>>& refPoint2 = point2;
      return refPoint2;
    }
  }

  std::shared_ptr<ETP::Point<T>> getPointVal( const int& _i ){
    if( _i == 0 ){
      return point1;
    }else{
      return point2;
    }
  }

  void setPoint( const size_t& _i, std::shared_ptr<ETP::Point<T>> _pt ){
    if( _i == 0 ){
      point1 = _pt;
    }else if( _i == 1 ){
      point2 = _pt;
    }
  }

  bool hasSamePoints( const std::shared_ptr<ETP::Point<T>>& _pt1, const std::shared_ptr<ETP::Point<T>>& _pt2 ){
    if( ( point1->equals( _pt1 ) && point2->equals( _pt2 ) ) || ( point1->equals( _pt2 ) && point2->equals( _pt1 ) ) ){
      return true;
    }
    return false;
  }

  bool hasSamePoints( const std::shared_ptr<ETP::Edge<T>>& _e ){
    std::shared_ptr<Point<T>>& e_pt1 = _e->getPointRef( 0 );
    std::shared_ptr<Point<T>>& e_pt2 = _e->getPointRef( 1 );
    if( ( point1->equals( e_pt1 ) && point2->equals( e_pt2 ) ) || ( point1->equals( e_pt2 ) && point2->equals( e_pt1 ) ) ){
      return true;
    }
    return false;
  }

  bool hasSamePoints( const T& _p1x, const T& _p1y, const T& _p2x, const T& _p2y ){
    ETP::Point<T>& e_pt1 = getPointRef( 0 );
    ETP::Point<T>& e_pt2 = getPointRef( 1 );
    if( ( point1->equals( _p1x, _p1y ) && point2->equals( _p2x, _p2y ) ) || ( point1->equals( _p2x, _p2y ) && point2->equals( _p1x, _p1y ) ) ){
      return true;
    }
    return false;
  }

  std::vector<T> getPositions(){
      std::vector<T> ret;
      ret.push_back( point1->x );
      ret.push_back( point1->y );
      ret.push_back( point2->x );
      ret.push_back( point2->y );
      return ret;
  }

  bool hasPoint( const std::shared_ptr<ETP::Point<T>>& _pt ){
    if( point1->equals( _pt ) || point2->equals( _pt ) ){
      return true;
    }
    return false;
  }

  T getLength(){
    return std::sqrt( std::pow( point2->x - point1->x, 2 ) + std::pow( point2->y - point1->y, 2 ) );
  }





};

}

#endif
