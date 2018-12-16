#ifndef ETP_POINT
#define ETP_POINT

#include <iostream>
#include <memory>
#include <nan.h>
#include <cmath>

namespace ETP {

template <typename P, typename D>

struct Point {

  P x;
  P y;

  Point(){}

  Point( P _x, P _y ): x( _x ), y( _y ){

  }

  ~Point(){}

  void set( P _x, P _y ){
    x = _x;
    y = _y;
  }

  bool equals( const std::shared_ptr<Point<P,D>>& _Point ){
    if( x == _Point->x && y == _Point->y ){
      return true;
    }
    return false;
  }
  bool equals( const Point<P,D>& _Point ){
    if( x == _Point.x && y == _Point.y ){
      return true;
    }
    return false;
  }
  bool equals( const P& _x, const P& _y ){
    if( x == _x && y == _y ){
      return true;
    }
    return false;
  }

  v8::Local<v8::Object> toV8Obj(){
    v8::Local<v8::Object> ret = Nan::New<v8::Object>();
    v8::Local<v8::String> xProp = Nan::New( "x" ).ToLocalChecked();
    v8::Local<v8::Value> xValue = Nan::New( x );
    ret->Set( xProp, xValue );
    v8::Local<v8::String> yProp = Nan::New( "y" ).ToLocalChecked();
    v8::Local<v8::Value> yValue = Nan::New( y );
    ret->Set( yProp, yValue );
    return ret;
  }

  Point perp(){
    P magnitude = std::sqrt( x * x + y * y );
    return Point( -y / magnitude, x / magnitude );
  }

  D dot( Point& other ){
    return (D)( x * other.x + y * other.y );
  }

  D dist( Point& other ){
    D dx = (D) other.x - (D) x;
    D dy = (D) other.y - (D) y;
    return std::pow( dx * dx + dy * dy, 0.5 );
  }
//((p1 - q1)^2 + (p2 - q2)^2)^(1/2)
/*
  Point operator* ( D _x ){
      return Point( _x * x, _x * y  );
  }
*/
  Point operator* ( const D _x ){
      return Point( _x * x, _x * y  );
  }

};

}

#endif
