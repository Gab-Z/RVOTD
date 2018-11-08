#ifndef DEF_ETP_UTILS
#define DEF_ETP_UTILS

#include <iostream>
#include <vector>
#include <nan.h>
#include <memory>
#include <limits>
#include <cmath>

#include "./core.h"

namespace ETP {

template<class C, typename T>
C normalisedVector( C _pt1, C _pt2 ){
  T vecx = _pt2.x - _pt1.x;
  T vecy = _pt2.y - _pt1.y;
  T vecL = std::abs( std::sqrt( std::pow( vecx, 2 ) + std::pow( vecy, 2 ) ) );
  return C( vecx / vecL, vecy / vecL  );
}

template<class C, typename T>
v8::Local<v8::Array> ptVecToJsArr( std::vector<C> _vec ){
  v8::Local<v8::Array> ret = Nan::New<v8::Array>();
  size_t vl = _vec.size();
  for( size_t i = 0; i < vl; i++ ){
    ret->Set( i * 2, Nan::New( _vec[ i ].x ) );
    ret->Set( i * 2 + 1, Nan::New( _vec[ i ].y ) );
  }
  return ret;
}

}

#endif
