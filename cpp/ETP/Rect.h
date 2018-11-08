#ifndef DEF_ETP_RECT
#define DEF_ETP_RECT

#include <algorithm>

namespace ETP{

template<typename T>
class Rect {

  T Left;
  T Top;
  T Right;
  T Bottom;

public:
  Rect(){}
  Rect( T _left, T _top, T _right, T _bottom ): Left( _left ), Top( _top ), Right( _right ), Bottom( _bottom ){}

/*
  template<class C>
  Rect( C _pt1, C _pt2 ){
    Left = std::min( _pt1.x, _pt2.x );
    Top = std::min( _pt1.y, _pt2.y );
    Right = std::max( _pt1.x, _pt2.x );
    Bottom = std::max( _pt1.y, _pt2.y );
  }
  */
  ~Rect(){}

  T top(){ return Top; }
  T left(){ return Left; }
  T right(){ return Right; }
  T bottom(){ return Bottom; }

  void sort(){
    if( Left > Right ){
      T tmpLeft = Right;
      T tmpRight = Left;
      Right = tmpLeft;
      Left = tmpRight;
    }
    if( Top > Bottom ){
      T tmpTop = Bottom;
      T tmpBottom = Top;
      Top = tmpBottom;
      Bottom = tmpTop;
    }
  }

  bool intersect( ETP::Rect<T> _rect ){
    //return a[0].x <= b[1].x && a[1].x >= b[0].x
      //      && a[0].y <= b[1].y && a[1].y >= b[0].y;
    return ( Left <= _rect.right() && Right >= _rect.left() && Top <= _rect.bottom() && Bottom >= _rect.top() );
  }


};

}



#endif
