#ifndef ETP_POINT
#define ETP_POINT

#include <iostream>
#include <memory>
#include <nan.h>

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


};

}

#endif
