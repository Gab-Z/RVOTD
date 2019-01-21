#ifndef DEF_ETP_SPACE
#define DEF_ETP_SPACE

#include <iostream>
#include <vector>
#include <iterator>
#include <cmath>
#include <memory>
#include <utility>
#include <limits>
#include <stdexcept>
#include <string>
#include <algorithm>

#include <nan.h>

#include "../clipper/clipper.hpp"
#include "../poly2tri/poly2tri.h"

#include "./Point.h"
#include "./Edge.h"
#include "./SearchNode.h"
#include "./Triangle.h"


#include "../utilz.h"

typedef int NanInt;
typedef float NanFloat;
typedef float NanDouble;
typedef int component;

namespace ETP {

using N = uint32_t;

template <typename P, typename D>
class TriangulationSpace {

  std::vector<std::shared_ptr<ETP::Point<P,D>>> points;
  std::vector<std::shared_ptr<ETP::Edge<P,D>>> edges;
  std::vector<std::shared_ptr<ETP::Triangle<P,D>>> triangles;
  P minX;
  P minY;
  P maxX;
  P maxY;
  P sectorWidth = 15;
  P sectorHeight = 15;
  std::vector<std::vector<std::vector<std::shared_ptr<ETP::Triangle<P,D>>>>> sectors;
  //unsigned int searchIndex = 0;

public:

  TriangulationSpace(){};

  v8::Local<v8::Array> getTriangles(){
    size_t tl = triangles.size();
    v8::Local<v8::Array> ret = Nan::New<v8::Array>( tl );
    for( size_t ti = 0; ti < tl; ti++ ){
      std::shared_ptr<ETP::Triangle<P,D>> tri = triangles[ ti ];
      ret->Set( ti, tri->toNanObject() );
    }
    return ret;
  }
  v8::Local<v8::Array> convertTriangles( std::vector<std::shared_ptr<ETP::Triangle<P,D>>> _triangles ){
    size_t tl = _triangles.size();
    v8::Local<v8::Array> ret = Nan::New<v8::Array>( tl );
    for( size_t ti = 0; ti < tl; ti++ ){
      std::shared_ptr<ETP::Triangle<P,D>>& tri = _triangles[ ti ];
      ret->Set( ti, tri->toNanObject() );
    }
    return ret;
  }
  v8::Local<v8::Array> getEdges(){
    size_t el = edges.size();
    v8::Local<v8::Array> ret = Nan::New<v8::Array>( el );
    for( size_t i = 0; i < el; i++ ){
      v8::Local<v8::Object> edge = Nan::New<v8::Object>();

      v8::Local<v8::String> constrProp = Nan::New( "constrained" ).ToLocalChecked();
      v8::Local<v8::Value> constrValue = Nan::New( edges[ i ]->isConstrained() );
      edge->Set( constrProp, constrValue );

      std::vector<P> _positions = edges[ i ]->getPositions();

      v8::Local<v8::String> p1xProp = Nan::New( "p1x" ).ToLocalChecked();
      v8::Local<v8::Value> p1xValue = Nan::New( _positions[ 0 ] );
      edge->Set( p1xProp, p1xValue );

      v8::Local<v8::String> p1yProp = Nan::New( "p1y" ).ToLocalChecked();
      v8::Local<v8::Value> p1yValue = Nan::New( _positions[ 1 ] );
      edge->Set( p1yProp, p1yValue );

      v8::Local<v8::String> p2xProp = Nan::New( "p2x" ).ToLocalChecked();
      v8::Local<v8::Value> p2xValue = Nan::New( _positions[ 2 ] );
      edge->Set( p2xProp, p2xValue );

      v8::Local<v8::String> p2yProp = Nan::New( "p2y" ).ToLocalChecked();
      v8::Local<v8::Value> p2yValue = Nan::New( _positions[ 3 ] );
      edge->Set( p2yProp, p2yValue );

      ret->Set( i, edge );
    }
    return ret;
  }
  v8::Local<v8::Object> getSectors(){
    v8::Local<v8::Object> ret = Nan::New<v8::Object>();
    v8::Local<v8::String> wProp = Nan::New( "width" ).ToLocalChecked();
    v8::Local<v8::Value> wValue = Nan::New( (int) sectorWidth );
    ret->Set( wProp, wValue );
    v8::Local<v8::String> hProp = Nan::New( "height" ).ToLocalChecked();
    v8::Local<v8::Value> hValue = Nan::New( (int) sectorHeight );
    ret->Set( hProp, hValue );

    v8::Local<v8::String> mxProp = Nan::New( "minX" ).ToLocalChecked();
    v8::Local<v8::Value> mxValue = Nan::New( (int) minX );
    ret->Set( mxProp, mxValue );
    v8::Local<v8::String> myProp = Nan::New( "minY" ).ToLocalChecked();
    v8::Local<v8::Value> myValue = Nan::New( (int) minY );
    ret->Set( myProp, myValue );

    v8::Local<v8::Array> sSectors = Nan::New<v8::Array>();
    size_t sl = sectors.size();
    for( size_t y = 0; y < sl; y ++ ){
      std::vector<std::vector<std::shared_ptr<ETP::Triangle<P,D>>>> row = sectors[ y ];
      v8::Local<v8::Array> sRow = Nan::New<v8::Array>();
      size_t rl = row.size();
      for( size_t x = 0; x < rl; x++ ){
        std::vector<std::shared_ptr<ETP::Triangle<P,D>>> cell = row[ x ];
        v8::Local<v8::Array> sCell = Nan::New<v8::Array>();
        size_t cl = cell.size();
        for( size_t j = 0; j < cl; j++ ){
          sCell->Set( j, Nan::New( cell[ j ]->getId() ) );
        }
        sRow->Set( x, sCell );
      }
      sSectors->Set( y, sRow );
    }
    v8::Local<v8::String> sProp = Nan::New( "sectors" ).ToLocalChecked();
    ret->Set( sProp, sSectors );
    return ret;
  }
  v8::Local<v8::Object> getBounds(){
    v8::Local<v8::Object> ret = Nan::New<v8::Object>();

    v8::Local<v8::String> minProp = Nan::New( "min" ).ToLocalChecked();
    v8::Local<v8::String> xProp = Nan::New( "x" ).ToLocalChecked();
    v8::Local<v8::String> yProp = Nan::New( "y" ).ToLocalChecked();
    v8::Local<v8::String> maxProp = Nan::New( "max" ).ToLocalChecked();

    v8::Local<v8::Object> min = Nan::New<v8::Object>();
    v8::Local<v8::Value> minxValue = Nan::New( (int) minX );
    v8::Local<v8::Value> minyValue = Nan::New( (int) minY );
    min->Set( xProp, minxValue );
    min->Set( yProp, minyValue );
    ret->Set( minProp, min );

    v8::Local<v8::Object> max = Nan::New<v8::Object>();
    v8::Local<v8::Value> maxxValue = Nan::New( (int) maxX );
    v8::Local<v8::Value> maxyValue = Nan::New( (int) maxY );
    max->Set( xProp, maxxValue );
    max->Set( yProp, maxyValue );
    ret->Set( maxProp, max );
    return ret;
  }
  v8::Local<v8::Array> listTrianglesIds( std::vector<std::shared_ptr<ETP::Triangle<P,D>>> _triangles ){
    size_t l = _triangles.size();
    v8::Local<v8::Array> ret =  Nan::New<v8::Array>( l );
    for( int i = 0; i < l; i++ ){
      ret->Set( i, Nan::New( _triangles[ i ]->getId() ) );
    }
    return ret;
  }
  template<class C>
  void buildFromPolyLines( std::vector<std::vector<C>> polygons ){
    //std::vector<C> boundingPolygon = std::vector<C>();
    //std::vector<std::vector<C>> inputPolygons = std::vector<std::vector<C>>();
    //std::vector<C> outputTriangles = std::vector<C>();
    size_t mal = polygons.size();

    //check if any edges intersect, if so return an error
    for( size_t pli = 0; pli < mal; pli ++ ){
      std::vector<C>& usedPoly = polygons[ pli ];
      size_t ul = usedPoly.size();
      for( size_t ui = 0; ui < ul; ui++ ){
        size_t ni = ui + 1 > ul - 1 ? 0 : ui + 1;
        C pt1 = usedPoly[ ui ];
        C pt2 = usedPoly [ ni ];
        for( size_t pli2 = 0; pli2 < mal; pli2 ++ ){
          std::vector<C>& usedPoly2 = polygons[ pli2 ];
          size_t ul2 = usedPoly2.size();
          for( size_t ui2 = 0; ui2 < ul2; ui2++ ){
            size_t ni2 = ui2 + 1 > ul2 - 1 ? 0 : ui2 + 1;
            C pt3 = usedPoly2[ ui2 ];
            C pt4 = usedPoly2[ ni2 ];
            if( doLinesIntersect<C>( pt1, pt2, pt3, pt4 ) ){
              throw std::invalid_argument( "TowerDefense::triangulate - edges cannot intersect !" );
            }
          }
        }
      }
    }
    std::vector<std::vector<C>> p2tPolys;
    //loop through c2t::Point vectors to build constrained edges and remove points on flat lines
    for( size_t pi = 0; pi < mal; pi ++ ){
      std::vector<C>& usedPoly = polygons[ pi ];
      if( usedPoly[ 0 ].x == usedPoly[ usedPoly.size() - 1 ].x && usedPoly[ 0 ].y == usedPoly[ usedPoly.size() - 1 ].y ){
        usedPoly.erase( usedPoly.end() - 1 );
      }
      p2tPolys.push_back( std::vector<C>() );
      std::vector<C>& renderPoly = p2tPolys[ p2tPolys.size() - 1 ];
      size_t ul = usedPoly.size();
      C vec0 = normalisedVector<C>( usedPoly[ 0 ], usedPoly[ 1 ] );
      size_t startI = 0;
      for( size_t ri = 0; ri < ul; ri++ ){
        size_t pri = ul - 1 - ri;
        size_t nextI = pri + 1 > ul - 1 ? 0 : pri + 1;
        C vec = normalisedVector<C>( usedPoly[ pri ], usedPoly[ nextI ] );
        if( vec.x == -vec0.x && vec.y == -vec0.y ){
          throw std::invalid_argument( "TowerDefense::triangulate - polygons edges can't reverse at 180° !" );
        }
        if( ( vec.x != vec0.x || vec.y != vec0.y ) && ( vec.x != -vec0.x || vec.y != -vec0.y ) ){
          startI = nextI;
          break;
        }
        if( pri == 1){
          //return Nan::ThrowError( Nan::New( "TowerDefense::triangulate - polygons points can't be all on the same line !" ).ToLocalChecked() );
          throw std::invalid_argument( "TowerDefense::triangulate - polygons points can't be all on the same line !" );
        }
      }

      for( size_t fi = 0; fi < ul; fi+=0 ){
        size_t ui = startI + fi > ul - 1 ? startI + fi - ul : startI + fi;
        size_t nextUi = ui + 1 > ul - 1 ? 0 : ui + 1;
        C sPt = usedPoly[ ui ];
        C fPt = usedPoly[ nextUi ];
        C fVec = normalisedVector<C>( sPt, fPt );
        int incr = 1;
        for( size_t ti = 1; ti < ul; ti++ ){
          size_t si = ui + ti > ul - 1 ? ui + ti - ul : ui + ti;
          size_t nextSi = si + 1 > ul - 1 ? 0 : si + 1;
          C sVec = normalisedVector<C>( usedPoly[ si ], usedPoly[ nextSi ] );
          if( fVec.x == -sVec.x && fVec.y == -sVec.y ){
            throw std::invalid_argument( "TowerDefense::triangulate - polygons edges can't reverse at 180° !" );
          }else if( fVec.x != sVec.x || fVec.y != sVec.y ){
            fPt = usedPoly[ si ];
            renderPoly.push_back( sPt );
            addEdge( sPt.x, sPt.y, fPt.x, fPt.y, true );
            incr = ti;
            break;
          }
        }
        fi = fi + incr;
      }
    }

    std::vector<p2t::Point*> boundingLine_trgt;
    std::vector<C>& boundingLine_src = p2tPolys[ 0 ];
    size_t nbPts = boundingLine_src.size();
    for( size_t i = 0; i < nbPts; i++ ){
      C uPt = boundingLine_src[ i ];
      boundingLine_trgt.push_back( new p2t::Point( uPt.x, uPt.y ) );
    }
    p2t::CDT _cdt( std::move( boundingLine_trgt ) );

    for( size_t i = 1; i < mal; i++ ){
      std::vector<C>& usedLine = p2tPolys[ i ];
      size_t nbPts = usedLine.size();
      std::vector<p2t::Point*> polyline;
      for( size_t ii = 0; ii < nbPts; ii++ ){
        C uPt = usedLine[ ii ];
        polyline.push_back( new p2t::Point( uPt.x, uPt.y ) );
      }
      _cdt.AddHole( std::move( polyline ) );
    }


    _cdt.Triangulate();
    std::vector<p2t::Triangle*> outputTriangles = _cdt.GetTriangles();

    //triangulate the polygon using bounding and holes coordinates

    /*
    c2t::clip2tri clip2tri;
    try{
      clip2tri.triangulate( inputPolygons, outputTriangles, boundingPolygon );
    }catch(  const std::invalid_argument& e  ){
      throw;
    }
    */
    /*
    //concat all polygon points to a single vector
    for( size_t ipl = 0; ipl < inputPolygons.size(); ipl ++ ){
      boundingPolygon.insert(
        boundingPolygon.end(),
        std::make_move_iterator(inputPolygons[ ipl ].begin()),
        std::make_move_iterator(inputPolygons[ ipl ].end())
      );
    }

    //search the output points of triangulation to find closest points in the inputed ponts. This is because process of c2t triangulation breaks number precision
    size_t pl = outputTriangles.size();
    size_t ipl = boundingPolygon.size();
    for( size_t tl = 0; tl < pl; tl ++ ){
      P maxDist = std::numeric_limits<D>::max();
      size_t idx = 0;
      C searchedPt = outputTriangles[ tl ];
      for( size_t b = 0; b < ipl; b++ ){
        P dist = std::abs( searchedPt.x - boundingPolygon[ b ].x ) + std::abs( searchedPt.y - boundingPolygon[ b ].y );
        if( dist < maxDist ){
          maxDist = dist;
          idx = b;
        }
      }
      outputTriangles[ tl ].x = boundingPolygon[ idx ].x;
      outputTriangles[ tl ].y = boundingPolygon[ idx ].y;
    }
    */

    //create unconstrained edges. Which are all edges created by the triangulation process
    size_t idl = outputTriangles.size();
    for( size_t id = 0; id < idl; id++ ){
      p2t::Triangle* t = outputTriangles[ id ];
      p2t::Point* pt0 = t->GetPoint( 0 );
      p2t::Point* pt1 = t->GetPoint( 1 );
      p2t::Point* pt2 = t->GetPoint( 2 );
      //C& p1 = outputTriangles[ id * 3 ];
      //C& p2 = outputTriangles[ id * 3 + 1 ];
      //C& p3 = outputTriangles[ id * 3 + 2 ];
      addTriangle(  addGetEdge( pt0->x, pt0->y, pt1->x, pt1->y, false ),
                    addGetEdge( pt1->x, pt1->y, pt2->x, pt2->y, false ),
                    addGetEdge( pt2->x, pt2->y, pt0->x, pt0->y, false )
      );
    }


    setAdjacentTriangles();
    setWidthsAndAngles();

    try{
      abstractLevel2();
    }catch( const std::invalid_argument& e ){
      throw;
    }

    assignSectors();

    //info.GetReturnValue().Set( self->space->getEdges() );
    //info.GetReturnValue().Set( self->space->getTriangles() );

  }
  void setSectorsDimensions( P _width, P _height ){
    sectorWidth = _width;
    sectorHeight = _height;
  }
  bool isNewEdge( const ETP::Point<P,D>& _pt1, const ETP::Point<P,D>& _pt2 ){
    for( size_t i = 0; i < edges.size(); i++ ){
      if( edges[ i ]->hasSamePoints( _pt1, _pt2 ) ){
        return false;
      }
    }
    return true;
  }
  bool isNewEdge( const P& _p1x, const P& _p1y, const P& _p2x, const P& _p2y ){
    for( size_t i = 0; i < edges.size(); i++ ){
      if( edges[ i ]->hasSamePoints( _p1x, _p1y, _p2x, _p2y ) ){
        return false;
      }
    }
    return true;
  }
  bool isNewPoint(  const std::shared_ptr<ETP::Point<P,D>>& _pt ){
    for( size_t i = 0; i < points.size(); i++ ){
      if( points[ i ]->equals( _pt ) ){
        return false;
      }
    }
    return true;
  }
  bool isNewPoint(  const P& _x, const P& _y ){
    for( size_t i = 0; i < points.size(); i++ ){
      if( points[ i ]->equals( _x, _y ) ){
        return false;
      }
    }
    return true;
  }
  std::shared_ptr<ETP::Point<P,D>> getPointByValue( std::shared_ptr<ETP::Point<P,D>> _pt ){
    for( size_t i = 0; i < points.size(); i++ ){
      if( points[ i ]->equals( _pt ) ){
        return points[ i ];
      }
    }
  }
  void addPoint( std::shared_ptr<ETP::Point<P,D>> _pt ){
    if( points.size() == 0 ){
      minX = _pt->x;
      minY = _pt->y;
      maxX = _pt->x;
      maxY = _pt->y;
    }else{
      minX = std::min( minX, _pt->x );
      minY = std::min( minY, _pt->y );
      maxX = std::max( maxX, _pt->x );
      maxY = std::max( maxY, _pt->y );
    }
    points.push_back( std::move( _pt ) );
  }
  void addEdge( std::shared_ptr<ETP::Edge<P,D>> _e ){
    edges.push_back( std::move( _e ) );
  }
  std::shared_ptr<ETP::Point<P,D>> addGetPoint( std::shared_ptr<ETP::Point<P,D>> _pt ){
    for( size_t i = 0; i < points.size(); i++ ){
      if( points[ i ]->equals( _pt ) ){
        return points[ i ];
      }
    }
    addPoint( std::move( _pt ) );
    return points[ points.size() - 1 ];
  }
  std::shared_ptr<ETP::Point<P,D>> addGetPoint( const P& _x, const P& _y ){
    for( size_t i = 0; i < points.size(); i++ ){
      if( points[ i ]->equals( _x, _y ) ){
        return points[ i ];
      }
    }
    addPoint( std::move( std::make_shared<ETP::Point<P,D>>( _x, _y ) ) );
    return points[ points.size() - 1 ];
  }
  std::shared_ptr<ETP::Edge<P,D>> addGetEdge( const P& _p1x, const P& _p1y, const P& _p2x, const P& _p2y, const bool& _constrained ){
    std::shared_ptr<ETP::Point<P,D>> p1 = addGetPoint( _p1x, _p1y );
    std::shared_ptr<ETP::Point<P,D>> p2 = addGetPoint( _p2x, _p2y );
    for( size_t i = 0; i < edges.size(); i++ ){
      if( edges[ i ]->hasSamePoints( p1, p2 ) ){
        return edges[ i ];
      }
    }
    edges.push_back( std::move( std::make_shared<ETP::Edge<P,D>>( p1, p2, _constrained ) ) );
    return edges[ edges.size() - 1 ];
  }
  void addEdge( const P& _p1x, const P& _p1y, const P& _p2x, const P& _p2y, const bool& _constrained  ){
    std::shared_ptr<ETP::Point<P,D>> p1 = addGetPoint( _p1x, _p1y );
    std::shared_ptr<ETP::Point<P,D>> p2 = addGetPoint( _p2x, _p2y );
    for( size_t i = 0; i < edges.size(); i++ ){
      if( edges[ i ]->hasSamePoints( p1, p2 ) ){
        return void();
      }
    }
    edges.push_back( std::move( std::make_shared<ETP::Edge<P,D>>( p1, p2, _constrained ) ) );
  }
  void addTriangle( std::shared_ptr<ETP::Edge<P,D>> _edge1, std::shared_ptr<ETP::Edge<P,D>> _edge2, std::shared_ptr<ETP::Edge<P,D>> _edge3 ){
    triangles.push_back( std::move( std::make_shared<ETP::Triangle<P,D>>( _edge1, _edge2, _edge3, triangles.size() ) ) );
  }
  void setAdjacentTriangles(){
    size_t tl = triangles.size();
    for( size_t i = 0; i < tl; i++ ){
      std::shared_ptr<ETP::Triangle<P,D>>& tri = triangles[ i ];
      for( int ei = 0; ei < 3; ei ++ ){
        std::shared_ptr<ETP::Edge<P,D>>& triEdge = tri->getEdgeRef( ei );
        if( triEdge->isConstrained() || tri->getAdjacentRef( ei ) != nullptr ){ continue; }
        for( size_t sti = 0; sti < tl; sti ++ ){
          std::shared_ptr<ETP::Triangle<P,D>>& searchedTri = triangles[ sti ];
          bool adjacentFound = false;
          if( searchedTri == tri ){ continue; }
          for( int sei = 0; sei  < 3; sei++ ){
            if( searchedTri->getEdgeRef( sei )->hasSamePoints( triEdge ) ){
              tri->setAdjacent( ( size_t ) ei, searchedTri );
              searchedTri->setAdjacent( (size_t ) sei, tri );
              adjacentFound = true;
              break;
            }
          }
          if( adjacentFound ){ break; }
        }
      }
    }
  }
  void setWidthsAndAngles(){
    size_t tl = triangles.size();
    for( size_t i = 0; i < tl; i++ ){
      std::shared_ptr<ETP::Triangle<P,D>> tri = triangles[ i ];
      tri->setWidth( 0, calculateWidth( tri, tri->getEdgeVal( 0 ), tri->getEdgeVal( 1 ) ) );
      tri->setWidth( 1, calculateWidth( tri, tri->getEdgeVal( 1 ), tri->getEdgeVal( 2 ) ) );
      tri->setWidth( 2, calculateWidth( tri, tri->getEdgeVal( 2 ), tri->getEdgeVal( 0 ) ) );

      tri->setAngle( 0, getAngle( tri->getEdgeVal( 0 ), tri->getEdgeVal( 1 ) ) );
      tri->setAngle( 1, getAngle( tri->getEdgeVal( 1 ), tri->getEdgeVal( 2 ) ) );
      tri->setAngle( 2, getAngle( tri->getEdgeVal( 2 ), tri->getEdgeVal( 0 ) ) );
    }
  }
  static D distanceBetween( const std::shared_ptr<ETP::Point<P,D>>& _pt, std::shared_ptr<ETP::Edge<P,D>> _edge ){
    std::shared_ptr<ETP::Point<P,D>> point1 = _edge->getPointVal( 0 );
    std::shared_ptr<ETP::Point<P,D>> point2 = _edge->getPointVal( 1 );
    if( point1->x == point2->x ){
      return (D) std::abs(point1->x - _pt->x);
    }
    D rise = (D) point2->y - (D) point1->y;
    D run = (D) point2->x - (D) point1->x;
    D intercept = (D) point1->y - ( rise / run ) * (D) point1->x;
    D a = rise;
    D b = -run;
    D c = run * intercept;
    return std::abs( a * (D) _pt->x + b * (D) _pt->y + c ) / ( std::sqrt( a * a + b * b ) );
  }
  static D distanceBetween( P _ptx, P _pty, P _edge1x, P _edge1y, P _edge2x, P _edge2y ){
    if( _edge1x == _edge2x ){
      return (D) std::abs( _edge1x - _ptx );
    }
    D rise = (D) _edge2y - (D) _edge1y;
    D run = (D) _edge2x - (D) _edge1x;
    D intercept = (D) _edge1y - ( rise / run ) * (D) _edge1x;
    D a = rise;
    D b = -run;
    D c = run * intercept;
    return std::abs( a * (D) _ptx + b * (D) _pty + c ) / ( std::sqrt( a * a + b * b ) );
  }
  static D searchWidth( const std::shared_ptr<ETP::Point<P,D>>& _pt, std::shared_ptr<ETP::Triangle<P,D>> _triangle, std::shared_ptr<ETP::Edge<P,D>> _edge, D _distance ){
    std::shared_ptr<ETP::Point<P,D>> point1 = _edge->getPointVal( 0 );
    std::shared_ptr<ETP::Point<P,D>> point2 = _edge->getPointVal( 1 );
    if( isObtuse( _pt, point1, point2 ) || isObtuse( _pt, point2, point1 ) ){
      return _distance;
    }
    D dist2 = distanceBetween( _pt, _edge );
    if( dist2 > _distance ){
      return _distance;
    }else if( _edge->isConstrained() ){
      return dist2;
    }else{
      std::shared_ptr<ETP::Triangle<P,D>> tri2 = _triangle->triangleOpposite( _edge );
      std::vector<std::shared_ptr<ETP::Edge<P,D>>> otherEdges = tri2->otherEdges( _edge );
      D dist3 = searchWidth( _pt, tri2, otherEdges[ 0 ], _distance );
      return searchWidth( _pt, tri2, otherEdges[ 1 ], dist3 );
    }
  }
  static D calculateWidth( std::shared_ptr<ETP::Triangle<P,D>> _triangle, std::shared_ptr<ETP::Edge<P,D>> _edge1, std::shared_ptr<ETP::Edge<P,D>> _edge2 ){
    if( _edge1->isConstrained() || _edge2->isConstrained() ){ return (D) 0.0f; }
    std::shared_ptr<ETP::Point<P,D>> ptC = vertexBetweenEdges( _edge1, _edge2 );
    std::shared_ptr<ETP::Edge<P,D>> edgc = _triangle->edgeOpposite( ptC );
    std::shared_ptr<ETP::Point<P,D>> ptA = _triangle->vertexOpposite( _edge1 );
    std::shared_ptr<ETP::Point<P,D>> ptB = _triangle->vertexOpposite( _edge2 );
    D d = std::min( _edge1->getLength(), _edge2->getLength() );
    if( isObtuse( ptC, ptA, ptB ) || isObtuse( ptC, ptB, ptA ) ){
      return d;
    }else if( edgc->isConstrained() ){
      return distanceBetween( ptC, edgc );
    }else{
      return searchWidth( ptC, _triangle, edgc, d );
    }
  }
  static bool isObtuse( std::shared_ptr<ETP::Point<P,D>> _pt1, std::shared_ptr<ETP::Point<P,D>> _pt2, std::shared_ptr<ETP::Point<P,D>> _pt3 ){
    if( getAngle( _pt1, _pt2, _pt3) >= (D) M_PI / 2.0f ){
      return true;
    }
    return false;
  }
  static D getAngle( std::shared_ptr<ETP::Point<P,D>> _pt1, std::shared_ptr<ETP::Point<P,D>> _vertex, std::shared_ptr<ETP::Point<P,D>> _pt2 ){
    P x1 = _pt1->x - _vertex->x;
    P y1 = _pt1->y - _vertex->y;
    P x2 = _pt2->x - _vertex->x;
    P y2 = _pt2->y - _vertex->y;
    P dot = x1 * x2 + y1 * y2;
    P det = x1 * y2 - y1 * x2;
    D angle = std::abs( std::atan2( (D) det, (D) dot ) );
    if( angle > (D) M_PI ){
      return ( (D) M_PI * 2.0 - angle );
    }
    return angle;
  }
  static D getAngle( std::shared_ptr<ETP::Edge<P,D>> _edge1, std::shared_ptr<ETP::Edge<P,D>> _edge2 ){
    if( _edge1 == nullptr || _edge2 == nullptr || _edge1->hasSamePoints( _edge2 ) ){ return (D) 0.0f; }
    std::shared_ptr<ETP::Point<P,D>> pt1;
    std::shared_ptr<ETP::Point<P,D>> vertex;
    std::shared_ptr<ETP::Point<P,D>> pt2;
    if( _edge1->getPointVal( 0 )->equals( _edge2->getPointVal( 0 ) ) ){
      vertex = _edge1->getPointVal( 0 );
      pt1 = _edge1->getPointVal( 1 );
      pt2 = _edge2->getPointVal( 1 );
    }else if( _edge1->getPointVal( 0 )->equals( _edge2->getPointVal( 1 ) ) ){
      vertex = _edge1->getPointVal( 0 );
      pt1 = _edge1->getPointVal( 1 );
      pt2 = _edge2->getPointVal( 0 );
    }else if( _edge1->getPointVal( 1 )->equals( _edge2->getPointVal( 0 ) ) ){
      vertex = _edge1->getPointVal( 1 );
      pt1 = _edge1->getPointVal( 0 );
      pt2 = _edge2->getPointVal( 1 );
    }else if( _edge1->getPointVal( 1 )->equals( _edge2->getPointVal( 1 ) ) ){
      vertex = _edge1->getPointVal( 1 );
      pt1 = _edge1->getPointVal( 0 );
      pt2 = _edge2->getPointVal( 0 );
    }else{
      return (D) 0.0f;
    }
    P x1 = pt1->x - vertex->x;
    P y1 = pt1->y - vertex->y;
    P x2 = pt2->x - vertex->x;
    P y2 = pt2->y - vertex->y;
    P dot = x1 * x2 + y1 * y2;
    P det = x1 * y2 - y1 * x2;
    D angle = std::abs( std::atan2( (D) det, (D) dot ) );
    if( angle > M_PI ){
      return ( M_PI * 2.0 - angle );
    }
    return angle;

  }
  /*
  unsigned int getSearchIndex(){
    searchIndex++;
    if( searchIndex == 100 ){ searchIndex = 0; }
    return searchIndex;
  }
  */
  static std::shared_ptr<ETP::Point<P,D>> vertexBetweenEdges( std::shared_ptr<ETP::Edge<P,D>> _edge1, std::shared_ptr<ETP::Edge<P,D>> _edge2 ){
    std::shared_ptr<ETP::Point<P,D>> e1pt1 = _edge1->getPointVal( 0 );
    if( e1pt1 == _edge2->getPointVal( 0 ) || e1pt1 == _edge2->getPointVal( 1 ) ){
      return e1pt1;
    }else{
      return _edge1->getPointVal( 1 );
    }
  }
  template<class C>
  C normalisedVector( C _pt1, C _pt2 ){
    D vecx = (D) ( _pt2.x - _pt1.x );
    D vecy = (D) ( _pt2.y - _pt1.y );
    D vecL = std::abs( std::sqrt( std::pow( vecx, 2 ) + std::pow( vecy, 2 ) ) );
    return C( vecx / vecL, vecy / vecL  );
  }
  template<class C, typename N>
  static N crossProduct( C a, C b ) {
    return (N) a.x * (N) b.y - (N) b.x * (N) a.y;
  }
  template<class C>
  static bool doBoundingBoxesIntersect( C _a0, C _a1, C _b0, C _b1 ) {
    C a0( std::min( _a0.x, _a1.x ), std::min( _a0.y, _a1.y ) );
    C a1( std::max( _a0.x, _a1.x ), std::max( _a0.y, _a1.y ) );
    C b0( std::min( _b0.x, _b1.x ), std::min( _b0.y, _b1.y ) );
    C b1( std::max( _b0.x, _b1.x ), std::max( _b0.y, _b1.y ) );
    return ( a0.x <= b1.x && a1.x >= b0.x && a0.y <= b1.y && a1.y >= b0.y );
  }
  template<class C>
  static bool isPointOnLine( C segPt1, C segPt2, C _pt ) {
    // Move the image, so that a.first is on (0|0)
    C bTmp( _pt.x - segPt1.x, _pt.y - segPt1.y );
    P r = crossProduct<C,P>( C( segPt2.x - segPt1.x, segPt2.y - segPt1.y ), bTmp );
    return ( std::abs( r ) < 0.000001 );
  }
  template<class C>
  static bool isPointRightOfLine( C segPt1, C segPt2, C _pt ) {
    // Move the image, so that a.first is on (0|0)
    //LineSegment aTmp = new LineSegment(new Point(0, 0), new Point( a.second.x - a.first.x, a.second.y - a.first.y ) );
    //Point bTmp = new Point(b.x - a.first.x, b.y - a.first.y);
    C bTmp( _pt.x - segPt1.x, _pt.y - segPt1.y );
    return crossProduct<C,P>( C( segPt2.x - segPt1.x, segPt2.y - segPt1.y ), bTmp ) < 0;
  }
  template<class C>
  static bool lineSegmentTouchesOrCrossesLine( C seg1Pt1, C seg1Pt2, C seg2Pt1, C seg2Pt2 ) {
    return  (   isPointOnLine<C>( seg1Pt1, seg1Pt2, seg2Pt1 )
            ||  isPointOnLine<C>( seg1Pt1, seg1Pt2, seg2Pt2 )
            || ( isPointRightOfLine<C>( seg1Pt1, seg1Pt2, seg2Pt1 ) ^ isPointRightOfLine<C>( seg1Pt1, seg1Pt2, seg2Pt2 ) )
    );
  }
  template<class C>
  static bool doLinesIntersect(  C seg1Pt1, C seg1Pt2, C seg2Pt1, C seg2Pt2  ) {
    if( seg1Pt1.equals( seg2Pt1 ) || seg1Pt1.equals( seg2Pt2 ) || seg1Pt2.equals( seg2Pt1 ) || seg1Pt2.equals( seg2Pt2 ) ){
      return false;
    }
    return (  doBoundingBoxesIntersect<C>( seg1Pt1, seg1Pt2, seg2Pt1, seg2Pt2 )
              && lineSegmentTouchesOrCrossesLine<C>(seg1Pt1, seg1Pt2, seg2Pt1, seg2Pt2)
              && lineSegmentTouchesOrCrossesLine<C>( seg2Pt1, seg2Pt2, seg1Pt1, seg1Pt2 )
    );
  }
  template<class C>
  static bool doLinesIntersectOrTouch(  C seg1Pt1, C seg1Pt2, C seg2Pt1, C seg2Pt2  ) {
    return (  doBoundingBoxesIntersect<C>( seg1Pt1, seg1Pt2, seg2Pt1, seg2Pt2 )
              && lineSegmentTouchesOrCrossesLine<C>(seg1Pt1, seg1Pt2, seg2Pt1, seg2Pt2)
              && lineSegmentTouchesOrCrossesLine<C>( seg2Pt1, seg2Pt2, seg1Pt1, seg1Pt2 )
    );
  }
  void assignSectors(){
    P width = maxX - minX;
    P height = maxY - minY;
    int nbX = (int) std::ceil( width / (P) sectorWidth );
    int nbY = (int) std::ceil( height / (P) sectorHeight );
    std::vector<std::vector<std::shared_ptr<ETP::Triangle<P,D>>>> row;
    for( size_t c = 0; c <= nbX; c++ ){
      std::vector<std::shared_ptr<ETP::Triangle<P,D>>> cell;
      row.push_back( cell );
    }
    for( size_t r = 0; r <= nbY; r++ ){
      sectors.push_back( row );
    }
    size_t nbTris = triangles.size();
    for( size_t i = 0; i < nbTris; i++ ){
      std::shared_ptr<ETP::Triangle<P,D>> t = triangles[ i ];
      std::vector<ETP::Point<P,D>> bound = t->getBoundingBox();
      P MinX = std::floor( ( bound[ 0 ].x - minX ) / (P) sectorWidth );
      P MinY = std::floor( ( bound[ 0 ].y - minY ) / (P) sectorHeight );
      P MaxX = std::floor( ( bound[ 1 ].x - minX ) / (P) sectorWidth );
      P MaxY = std::floor( ( bound[ 1 ].y - minY ) / (P) sectorHeight );
      std::vector<std::shared_ptr<ETP::Point<P,D>>> tPoints = t->getPoints();
      ETP::Point<P,D> pt1 = *tPoints[ 0 ].get();
      ETP::Point<P,D> pt2 = *tPoints[ 1 ].get();
      ETP::Point<P,D> pt3 = *tPoints[ 2 ].get();

      for( size_t y = MinY; y <= MaxY; y++ ){
        for( size_t x = MinX; x <= MaxX; x++ ){
          ETP::Point<P,D> gridPtMin = ETP::Point<P,D>( x * sectorWidth + minX,                  y * sectorHeight + minY );
          ETP::Point<P,D> gridPtMaxX = ETP::Point<P,D>( x * sectorWidth + sectorWidth + minX,   y * sectorHeight + minY );
          ETP::Point<P,D> gridPtMaxY = ETP::Point<P,D>( x * sectorWidth + minX,                 y * sectorHeight + sectorHeight + minY );
          ETP::Point<P,D> gridPtMaxXY = ETP::Point<P,D>( x * sectorWidth + sectorWidth + minX,  y * sectorHeight + sectorHeight + minY );
          if(     ( pt1.x >= gridPtMin.x && pt1.y >= gridPtMin.y && pt1.x <= gridPtMaxX.x && pt1.y <= gridPtMaxY.y )
              ||  ( pt2.x >= gridPtMin.x && pt2.y >= gridPtMin.y && pt2.x <= gridPtMaxX.x && pt2.y <= gridPtMaxY.y )
              ||  ( pt3.x >= gridPtMin.x && pt3.y >= gridPtMin.y && pt3.x <= gridPtMaxX.x && pt3.y <= gridPtMaxY.y )
              ||  ( t->isPointInside( gridPtMin ) == true )
              ||  ( t->isPointInside( gridPtMaxX ) == true )
              ||  ( t->isPointInside( gridPtMaxY ) == true )
              ||  ( t->isPointInside( gridPtMaxXY ) == true )

              ||  doLinesIntersectOrTouch<ETP::Point<P,D>>(  pt1, pt2, gridPtMin, gridPtMaxX )
              ||  doLinesIntersectOrTouch<ETP::Point<P,D>>(  pt1, pt2, gridPtMin, gridPtMaxY )
              ||  doLinesIntersectOrTouch<ETP::Point<P,D>>(  pt1, pt2, gridPtMaxX, gridPtMaxXY )
              ||  doLinesIntersectOrTouch<ETP::Point<P,D>>(  pt1, pt2, gridPtMaxY, gridPtMaxXY )

              ||  doLinesIntersectOrTouch<ETP::Point<P,D>>(  pt2, pt3, gridPtMin, gridPtMaxX )
              ||  doLinesIntersectOrTouch<ETP::Point<P,D>>(  pt2, pt3, gridPtMin, gridPtMaxY )
              ||  doLinesIntersectOrTouch<ETP::Point<P,D>>(  pt2, pt3, gridPtMaxX, gridPtMaxXY )
              ||  doLinesIntersectOrTouch<ETP::Point<P,D>>(  pt2, pt3, gridPtMaxY, gridPtMaxXY )

              ||  doLinesIntersectOrTouch<ETP::Point<P,D>>(  pt3, pt1, gridPtMin, gridPtMaxX )
              ||  doLinesIntersectOrTouch<ETP::Point<P,D>>(  pt3, pt1, gridPtMin, gridPtMaxY )
              ||  doLinesIntersectOrTouch<ETP::Point<P,D>>(  pt3, pt1, gridPtMaxX, gridPtMaxXY )
              ||  doLinesIntersectOrTouch<ETP::Point<P,D>>(  pt3, pt1, gridPtMaxY, gridPtMaxXY )
          ){
            sectors[ y ][ x ].push_back( t );
          }
        }
      }
    }
  }

  std::shared_ptr<ETP::Triangle<P,D>> getTriangleWithPoint( ETP::Point<P,D> _pt, D _scale ){
    D sectorY = std::floor( ( _pt.y - (D) minY * _scale) / ( (D) sectorHeight * _scale ) );
    if( sectorY >= (D) sectors.size() ){ return nullptr; }
    D sectorX = std::floor( ( _pt.x - (D) minX * _scale ) / ( (D) sectorWidth * _scale ) );
    if( sectorX >= (D) sectors[ sectorY ].size() ){ return nullptr; }
    std::vector<std::shared_ptr<ETP::Triangle<P,D>>> sector = sectors[ (int) sectorY ][ (int) sectorX ];
    size_t sl = sector.size();
    for( auto tri : sector ){
      if( tri->isPointInside( _pt, _scale ) ){ return tri; }
    }
    return nullptr;
  }
  std::shared_ptr<ETP::Triangle<P,D>> searchPointInAllTriangles( ETP::Point<P,D> _pt ){
    size_t tl = triangles.size();
    for( size_t i = 0; i < tl; i++ ){
      if( triangles[ i ]->isPointInside( _pt ) ){ return triangles[ i ]; }
    }
    return nullptr;
  }
  void collapseUnrootedTree( std::shared_ptr<ETP::Triangle<P,D>> t, component c ){
    std::vector<std::shared_ptr<ETP::Triangle<P,D>>> s;
    s.push_back( t );
    while( s.empty() == false ){
      std::shared_ptr<ETP::Triangle<P,D>> tCurrent = s.front();
      s.erase( s.begin() );
      tCurrent->setComponent( c );
      for( int i = 0; i < 3; i++ ){
        std::shared_ptr<ETP::Edge<P,D>> e = tCurrent->getEdgeVal( i );
        if( e->isConstrained() == true ){
          tCurrent->setAngle( i, 0.0f );
          tCurrent->setChoke( i, 0.0f );
        }else{
          tCurrent->setAngle( i, std::numeric_limits<D>::infinity() );
          tCurrent->setChoke( i, std::numeric_limits<D>::infinity() );
          std::shared_ptr<ETP::Triangle<P,D>> tNext = tCurrent->triangleOpposite( e );
          if( tNext->getComponent() == -1 ){
            s.insert( s.begin(), tNext );
          }
        }
      }
    }
  }
  void collapseRootedTree( std::shared_ptr<ETP::Triangle<P,D>> r, std::shared_ptr<ETP::Triangle<P,D>> t ){
    int c = r->getComponent();
    std::vector<std::shared_ptr<ETP::Triangle<P,D>>> s;
    s.push_back( t );
    std::vector<D> a;
    a.push_back( (D) 0.0f );
    std::vector<D> d;
    d.push_back( (D) 0.0f );
    std::vector<std::shared_ptr<ETP::Edge<P,D>>> edgesIn;
    edgesIn.push_back( r->getSharedEdgeWith( t ) );
    int lvl2RootToLvl1TreeIndex = r->getAdjacentIndex( t );
    while( s.empty() == false ){
      std::shared_ptr<ETP::Triangle<P,D>> tCurrent = s.back();
      s.pop_back();
      tCurrent->setComponent( c );
      D lastAngle = a.back();
      a.pop_back();
      D lastChoke = d.back();
      d.pop_back();
      std::shared_ptr<ETP::Edge<P,D>> edgeIn = edgesIn.back();
      edgesIn.pop_back();
      for( int i = 0; i < 3; i++ ){
          std::shared_ptr<ETP::Edge<P,D>> e =  tCurrent->getEdgeVal( i );
          std::shared_ptr<ETP::Triangle<P,D>> tCurrentNeighbour = tCurrent->getAdjacentVal( i );
          if( e->isConstrained() == true ){ continue; }
          if( e == edgeIn ){
            tCurrent->setConnectedNode( i, r );
            tCurrent->setIndexfromConnectedNode( i, lvl2RootToLvl1TreeIndex );
            continue;
          }
          D currAngle = lastAngle + getAngle( edgeIn, e );
          tCurrentNeighbour->setLowerBound( tCurrentNeighbour->getAdjacentIndex( tCurrent ), currAngle );
          D currChoke = std::min( lastChoke, calculateWidth( tCurrent, edgeIn, e ) );
          tCurrent->setChoke( i, currChoke );
          s.push_back( tCurrentNeighbour );
          a.push_back( currAngle );
          d.push_back( currChoke );
          edgesIn.push_back( e );
      }
    }
  }
  std::vector<std::shared_ptr<ETP::Triangle<P,D>>> abstractLevel0and1( component& c ){
    std::vector<std::shared_ptr<ETP::Triangle<P,D>>> q;
    std::vector<std::shared_ptr<ETP::Triangle<P,D>>> r;
    size_t tl = triangles.size();
    for( size_t ti = 0; ti < tl; ti++ ){
      std::shared_ptr<ETP::Triangle<P,D>> t = triangles[ ti ];
      t->setLevel( -1 );
      t->setComponent( -1 );
      t->setWidth( 0, calculateWidth( t, t->getEdgeVal( 0 ), t->getEdgeVal( 1 ) ) );
      t->setWidth( 1, calculateWidth( t, t->getEdgeVal( 1 ), t->getEdgeVal( 2 ) ) );
      t->setWidth( 2, calculateWidth( t, t->getEdgeVal( 2 ), t->getEdgeVal( 0 ) ) );
      int n = t->getNumConstrainedEdges();
      if( n == 3 ){
        t->setLevel( 0 );
        t->setComponent( c );
        c++;
        for( int i = 0; i < 3; i++ ){
          t->setAngle( i, 0.0f );
          t->setChoke( i, 0.0f );
          t->setAdjacent( i, nullptr );
        }
      }else if( n == 2 ){
        t->setLevel( 1 );
        for( int i = 0; i < 3; i++ ){
          std::shared_ptr<ETP::Edge<P,D>> e = t->getEdgeVal( i );
          if( e->isConstrained() == false ){
            q.push_back( t->triangleOpposite( e ) );
            break;
          }
        }
      }else{
        r.push_back( t );
      }
    }
    while( q.empty() == false ){
      std::shared_ptr<ETP::Triangle<P,D>> t = q.front();
      q.erase( q.begin() );
      if( t->getLevel() == -1 ){
        int n = t->getNumConstrainedEdges();
        int m = t->getNumAdjacentLevel( 1 );
        if( n + m  >= 2 ){
          t->setLevel( 1 );
          for( int i = 0; i < 3; i++ ){
            std::shared_ptr<ETP::Edge<P,D>> e = t->getEdgeVal( i );
            std::shared_ptr<ETP::Triangle<P,D>> tNext = t->triangleOpposite( e );
            if( e->isConstrained() == false && tNext->getLevel() == -1 ){
              q.push_back( tNext );
            }
          }
        }
        if( n + m == 3 ){
          collapseUnrootedTree( t, c );
          c++;
        }
      }
    }
    return r;
  }
  void abstractLevel2(){
    int c = 1;
    std::vector<std::shared_ptr<ETP::Triangle<P,D>>> q;
    q = abstractLevel0and1( c );
    for( size_t i = 0; i < q.size(); i++ ){
      std::shared_ptr<ETP::Triangle<P,D>> t = q[ i ];
      int n = t->getNumConstrainedEdges();
      int m = t->getNumAdjacentLevel( 1 );
      if( n + m == 0 && t->getLevel() == -1 ){
        abstractLevel3( t, c );
        c++;
      }
    }
    while( q.empty() == false ){
      std::shared_ptr<ETP::Triangle<P,D>> t = q.front();
      q.erase( q.begin() );
      if( t->getLevel() == -1 ){
        std::shared_ptr<ETP::Triangle<P,D>> tCurrent = t;
        while( tCurrent != nullptr ){
          tCurrent->setLevel( 2 );
          std::shared_ptr<ETP::Triangle<P,D>> tNext = nullptr;
          for( int i = 0; i < 3; i++ ){
            std::shared_ptr<ETP::Edge<P,D>> e = tCurrent->getEdgeVal( i );
            if( e->isConstrained() == true ){ continue; }
            //std::shared_ptr<ETP::Triangle<P,D>> tTemp = tCurrent->triangleOpposite( e );
            std::shared_ptr<ETP::Triangle<P,D>> tTemp = tCurrent->getAdjacentVal( i );
            if( tTemp->getLevel() == 1 ){
              collapseRootedTree( tCurrent, tTemp );
              //tCurrent->setAngle( i, 0.0f );
              //tCurrent->setChoke( i, 0.0f );
              //tCurrent->setAdjacent( i, nullptr );
            }else{
              if( tTemp->getLevel() == -1 ){
                tNext = tTemp;
              }
              //tCurrent->setAngle( i, std::numeric_limits<D>::infinity() );
              //tCurrent->setChoke( i, std::numeric_limits<D>::infinity() );
              //tCurrent->setLowerBound( i, 55.0 );
              //tCurrent->setAdjacent( i, nullptr );
            }
          }
          tCurrent = tNext;
        }
      }
    }
  }
  void abstractLevel3( std::shared_ptr<ETP::Triangle<P,D>> _t, int c ){
    std::vector<std::shared_ptr<ETP::Triangle<P,D>>> q = std::vector<std::shared_ptr<ETP::Triangle<P,D>>>();
    q.push_back( _t );
    while( q.empty() == false ){
      std::shared_ptr<ETP::Triangle<P,D>> t = q.front();
      q.erase( q.begin() );
      t->setLevel( 3 );
      t->setComponent( c );
      for( size_t i = 0; i < 3; i++ ){
        std::shared_ptr<ETP::Triangle<P,D>> tParent = t;
        std::shared_ptr<ETP::Triangle<P,D>> tCurrent = t->getAdjacentVal( i );
        if( tCurrent->getConnectedNode( tCurrent->getAdjacentIndex( tParent ) ) != nullptr ){ break; }
        std::shared_ptr<ETP::Edge<P,D>> edgeIn = t->getEdgeVal( i );
        std::shared_ptr<ETP::Edge<P,D>> edgeOut = nullptr;
        D lowerBoundCount = 0.0f;
        D chokeCount = edgeIn->getLength();
        for(;;){
          int constrainedAndLvl1Sum = tCurrent->getNumConstrainedEdges() + tCurrent->getNumAdjacentLevel( 1 );
          int currentToParentIndex = tCurrent->getAdjacentIndex( tParent );
          if( tCurrent->getNumConstrainedEdges() + tCurrent->getNumAdjacentLevel( 1 ) == 0 ){
            if( tCurrent->getLevel() == -1 ){
              q.push_back( tCurrent );
            }
            tCurrent->setConnectedNode( currentToParentIndex, t );
            tCurrent->setIndexfromConnectedNode( currentToParentIndex, i );
            break;
          }else{
            tCurrent->setLevel( 2 );
            tCurrent->setComponent( c );
            tCurrent->setConnectedNode( currentToParentIndex, t );
            tCurrent->setIndexfromConnectedNode( currentToParentIndex, i );
            std::shared_ptr<ETP::Triangle<P,D>> tNext;
            std::shared_ptr<ETP::Edge<P,D>> eNext;
            for( int j = 0; j < 3; j++ ){
              edgeOut = tCurrent->getEdgeVal( j );
              if( edgeOut->isConstrained() == true || edgeOut == edgeIn ){ continue; }
              std::shared_ptr<ETP::Triangle<P,D>> tNeighbour = tCurrent->getAdjacentVal( j );
              if( tNeighbour->getLevel() == 1 ){
                collapseRootedTree( tCurrent, tNeighbour );
              }else{
                tNext = tNeighbour;
                eNext = edgeOut;
                chokeCount = std::min( chokeCount, tCurrent->getWidthbetweenEdges( edgeIn, edgeOut ) );
                int neighbourToCurrentIdx = tNeighbour->getAdjacentIndex( tCurrent );
                tNeighbour->setChoke( neighbourToCurrentIdx, chokeCount );
                lowerBoundCount += getAngle( edgeIn, edgeOut );
                tNeighbour->setLowerBound( neighbourToCurrentIdx, lowerBoundCount );
              }
            }
            edgeIn = eNext;
            tParent = tCurrent;
            tCurrent = tNext;
          }
        }
      }
    }
  }



  std::vector<std::shared_ptr<ETP::Triangle<P,D>>> abstractTriangulationSearch( ETP::Point<P,D> _startPoint, ETP::Point<P,D> _goalPoint, D _radius, D _scale ){
    std::shared_ptr<ETP::Triangle<P,D>> _startTri = getTriangleWithPoint( _startPoint, _scale );
    std::shared_ptr<ETP::Triangle<P,D>> _goalTri = getTriangleWithPoint( _goalPoint, _scale );
    std::vector<std::shared_ptr<ETP::Triangle<P,D>>> ret;
    if( _startTri == nullptr || _goalTri == nullptr ){ return ret; }
    int sComp = _startTri->getComponent();
    int gComp = _goalTri->getComponent();
    if( sComp != gComp ){ return ret; }
    // on same triangle or within lvl 0 node
    if( _startTri == _goalTri ){
      ret.push_back( _startTri );
      return ret;
    }
    int sLvl = _startTri->getLevel();
    int gLvl = _goalTri->getLevel();

    //start and goal are connected lvl 1 and lvl 2 nodes
    if(  sLvl == 1 && gLvl == 2 ){ // start is lvl 1 and goal is lvl 2
      int connectedNodeIndex = _startTri->getConnectedNodeIndex( _goalTri );
      if( connectedNodeIndex > -1 ){
        return getTrianglesFromLvl1ToConnectedLvl2( _startTri, connectedNodeIndex );
      }
    }else if( sLvl == 2 && gLvl == 1 ){ // start is lvl 2 and goal is lvl 1
      int connectedNodeIndex = _goalTri->getConnectedNodeIndex( _startTri );
      if( connectedNodeIndex > -1 ){
        ret = getTrianglesFromLvl1ToConnectedLvl2( _goalTri, connectedNodeIndex );
        std::reverse( ret.begin(), ret.end() );
        return ret;
      }
    }
    //start and goal are inside the same lvl 1 tree
    if( sLvl == 1 && gLvl == 1 ){

      std::shared_ptr<ETP::Triangle<P,D>> startRoot = _startTri->getLvl1RootNode();
      std::shared_ptr<ETP::Triangle<P,D>> goalRoot = _goalTri->getLvl1RootNode();
      if(  startRoot == goalRoot && sComp == gComp  //start and goal belong to the same tree tree
      ){
        return searchLvl1Tree( _startTri, _goalTri );
      }
    }
    // start and goal belong either to the same level 2 loop or ring or to connected lvl 1 trees connected to it
    // includes also the case where start and goal belong to a level 2 corridor or to attached lvl 1 trees
    if( sLvl == 2 || ( sLvl == 1 && _startTri->getLvl1RootNode() != nullptr ) ){
      std::shared_ptr<ETP::Triangle<P,D>> startRoot = sLvl == 2 ? _startTri : _startTri->getLvl1RootNode();
      if( gLvl == 2 || ( gLvl == 1 && _goalTri->getLvl1RootNode() != nullptr ) ){
        std::shared_ptr<ETP::Triangle<P,D>> goalRoot = gLvl == 2 ? _goalTri : _goalTri->getLvl1RootNode();
        //return std::vector<std::shared_ptr<ETP::Triangle<P,D>>>{ startRoot, goalRoot };
        bool proceed = false;
        std::vector<std::shared_ptr<ETP::Triangle<P,D>>> bodySearch;
        if( startRoot->isPartoflvl2Ring() || startRoot->isOnSameLvl2Loop( goalRoot ) ){// inside a loop or ring
          proceed = true;
          return searchInsideLvl2RingOrLoop( _startPoint, _goalPoint, _startTri, _goalTri, _radius, _scale );
        }else if( startRoot->haveSameLvl2CorridorEndpoints( goalRoot )){// inside a corridor
          return searchInsideLvl2Corridor( _startPoint, _goalPoint, _startTri, _goalTri, _radius, _scale );
        }
        if( proceed == true ){
          if( _startTri != startRoot ){ // we start from a lvl 1 tree and addtriangles up to the lvl2 root
            std::vector<std::shared_ptr<ETP::Triangle<P,D>>> startPart = getTrianglesFromLvl1ToConnectedLvl2( _startTri, _startTri->getConnectedNodeIndex( startRoot ) );
            ret.insert(
              ret.end(),
              std::make_move_iterator( startPart.begin() ),
              std::make_move_iterator( startPart.end() - 1 )
           );
          }
          ret.insert(
            ret.end(),
            std::make_move_iterator( bodySearch.begin() ),
            std::make_move_iterator( bodySearch.end() )
          );
          if( _goalTri != goalRoot ){// goal is on a lvl 1 tree, we add triangles from the lvl 2 root to the goal
            std::vector<std::shared_ptr<ETP::Triangle<P,D>>> goalPart = getTrianglesFromLvl1ToConnectedLvl2( _goalTri, _goalTri->getConnectedNodeIndex( goalRoot ) );
            std::reverse( goalPart.begin(), goalPart.end() );
            ret.insert(
              ret.end(),
              std::make_move_iterator( goalPart.begin() + 1 ),
              std::make_move_iterator( goalPart.end() )
           );
          }
          return ret;
        }

      }
    }
    // perform classic Triangulation Reduction A*
    std::shared_ptr<ETP::UsedList<P,D>> usedList = std::make_shared<ETP::UsedList<P,D>>();

    std::vector<std::shared_ptr<ETP::searchEndPoint<P,D>>> goalEndPoints;
    if( gLvl == 3 ){
      goalEndPoints.push_back( std::make_shared<ETP::searchEndPoint<P,D>>( 0.0, _goalTri ) );
    }else if( gLvl == 2 ){
      std::shared_ptr<ETP::searchEndPoint<P,D>> lvl2EndPt = std::make_shared<ETP::searchEndPoint<P,D>>( 0.0, _goalTri );
      goalEndPoints.push_back( lvl2EndPt );
      std::shared_ptr<ETP::Triangle<P,D>> connectedNode_1 = nullptr;
      std::shared_ptr<ETP::Triangle<P,D>> connectedNode_2 = nullptr;
      int connectedNode_1_idx, connectedNode_2_idx;
      for( int i = 0; i < 3; i++ ){
        std::shared_ptr<ETP::Triangle<P,D>> connectedNode = _goalTri->getConnectedNode( i );
        if( connectedNode == nullptr ) continue;
        if( connectedNode_1 == nullptr ){
          connectedNode_1 = connectedNode;
          connectedNode_1_idx = i;
        }else{
          if( connectedNode != connectedNode_1 ){
            connectedNode_2 = connectedNode;
            connectedNode_2_idx = i;
          }else if( _goalTri->getLowerBound( connectedNode_1_idx ) > _goalTri->getLowerBound( i ) ){
            connectedNode_1 = connectedNode;
            connectedNode_1_idx = i;
          }
        }
      }
      int lvl3EndPt_1_To2_idx = _goalTri->getIndexfromConnectedNode( connectedNode_1_idx );
      goalEndPoints.push_back( std::make_shared<ETP::searchEndPoint<P,D>>( 0.0, connectedNode_1, lvl2EndPt, lvl3EndPt_1_To2_idx, getPointToEdgeLowerBound( _goalPoint, _goalTri->getEdgeVal( connectedNode_1_idx ), _radius, _scale ) + _goalTri->getLowerBound( connectedNode_1_idx ) * _radius ) );
      if( connectedNode_2 != nullptr ){
        int lvl3EndPt_2_To2_idx = _goalTri->getIndexfromConnectedNode( connectedNode_2_idx );
        goalEndPoints.push_back( std::make_shared<ETP::searchEndPoint<P,D>>( 0.0, connectedNode_2, lvl2EndPt, lvl3EndPt_2_To2_idx, getPointToEdgeLowerBound( _goalPoint, _goalTri->getEdgeVal( connectedNode_2_idx ), _radius, _scale ) + _goalTri->getLowerBound( connectedNode_2_idx ) * _radius ) );
      }
    }else if( gLvl == 1 ){
      std::shared_ptr<ETP::searchEndPoint<P,D>> lvl1EndPt = std::make_shared<ETP::searchEndPoint<P,D>>( 0.0, _goalTri );
      std::shared_ptr<ETP::Triangle<P,D>> lvl2 = _goalTri->getLvl1RootNode();
      int lvl1Tolvl2Idx = _goalTri->getConnectedNodeIndex( lvl2 );
      D lvl1ToLvl2LowerBound = getPointToEdgeLowerBound( _goalPoint, _goalTri->getEdgeVal( lvl1Tolvl2Idx ), _radius, _scale ) + _goalTri->getLowerBound( lvl1Tolvl2Idx ) * _radius;
      int lvl2ToLvl1Idx = _goalTri->getIndexfromConnectedNode( lvl1Tolvl2Idx );
      std::shared_ptr<ETP::searchEndPoint<P,D>> lvl2EndPt = std::make_shared<ETP::searchEndPoint<P,D>>( 0.0, lvl2, lvl1EndPt, lvl2ToLvl1Idx, lvl1ToLvl2LowerBound );
      goalEndPoints.push_back( lvl2EndPt );
      std::shared_ptr<ETP::Triangle<P,D>> connectedNode_1 = nullptr;
      std::shared_ptr<ETP::Triangle<P,D>> connectedNode_2 = nullptr;
      int connectedNode_1_idx, connectedNode_2_idx;
      for( int i = 0; i < 3; i++ ){
        std::shared_ptr<ETP::Triangle<P,D>> connectedNode = lvl2->getConnectedNode( i );
        if( connectedNode == nullptr ) continue;
        if( connectedNode_1 == nullptr ){
          connectedNode_1 = connectedNode;
          connectedNode_1_idx = i;
        }else{
          if( connectedNode != connectedNode_1 ){
            connectedNode_2 = connectedNode;
            connectedNode_2_idx = i;
          }else if( lvl2->getLowerBound( connectedNode_1_idx ) > lvl2->getLowerBound( i ) ){
            connectedNode_1 = connectedNode;
            connectedNode_1_idx = i;
          }
        }
        int lvl3EndPt_1_To2_idx = lvl2->getIndexfromConnectedNode( connectedNode_1_idx );
        std::shared_ptr<ETP::Edge<P,D>> lvl2entryEdge = lvl2->getEdgeVal( lvl2ToLvl1Idx );
        goalEndPoints.push_back( std::make_shared<ETP::searchEndPoint<P,D>>( 0.0, connectedNode_1, lvl2EndPt, lvl3EndPt_1_To2_idx, lvl1ToLvl2LowerBound + ( lvl2->getAngle( lvl2entryEdge, lvl2->getEdgeVal( connectedNode_1_idx ) ) + lvl2->getLowerBound( connectedNode_1_idx ) ) * _radius ) );
        if( connectedNode_2 != nullptr ){
          int lvl3EndPt_2_To2_idx = lvl2->getIndexfromConnectedNode( connectedNode_2_idx );
          goalEndPoints.push_back( std::make_shared<ETP::searchEndPoint<P,D>>( 0.0, connectedNode_2, lvl2EndPt, lvl3EndPt_2_To2_idx, lvl1ToLvl2LowerBound + ( lvl2->getAngle( lvl2entryEdge, lvl2->getEdgeVal( connectedNode_2_idx ) ) + lvl2->getLowerBound( connectedNode_1_idx ) ) * _radius ) );
        }
      }
    }
    std::vector<std::shared_ptr<ETP::Triangle<P,D>>> startPoints;
    if( sLvl == 3 ){
      setTriangleSearchNode( _startTri, ETP::SearchNode<P,D>( 0.0, 0.0, nullptr ), usedList );
      startPoints.push_back( _startTri );
    }else if( sLvl == 2 ){
      std::shared_ptr<ETP::Triangle<P,D>> connectedNode_1 = nullptr;
      std::shared_ptr<ETP::Triangle<P,D>> connectedNode_2 = nullptr;
      int connectedNode_1_idx, connectedNode_2_idx;
      for( int i = 0; i < 3; i++ ){
        std::shared_ptr<ETP::Triangle<P,D>> connectedNode = _startTri->getConnectedNode( i );
        if( connectedNode == nullptr ) continue;
        if( connectedNode_1 == nullptr ){
          connectedNode_1 = connectedNode;
          connectedNode_1_idx = i;
        }else{
          if( connectedNode != connectedNode_1 ){
            connectedNode_2 = connectedNode;
            connectedNode_2_idx = i;
          }else if( _startTri->getLowerBound( connectedNode_1_idx ) > _startTri->getLowerBound( i ) ){
            connectedNode_1 = connectedNode;
            connectedNode_1_idx = i;
          }
        }
      }
      int lvl3EndPt_1_To2_idx = _startTri->getIndexfromConnectedNode( connectedNode_1_idx );
      std::shared_ptr<ETP::Edge<P,D>> lvl3_1_entryEdge = connectedNode_1->getEdgeVal( lvl3EndPt_1_To2_idx );
      D lvl3_1_hVal = lvl3_1_entryEdge->minimumDistance( _goalPoint, _scale ) / _scale;
      D lvl3_1_gValue = std::max(
        lvl3_1_entryEdge->minimumDistance( _startPoint, _scale ) / _scale,
        getPointToEdgeLowerBound( _startPoint, _startTri->getEdgeVal( connectedNode_1_idx ), _radius, _scale ) + _startTri->getLowerBound( connectedNode_1_idx )  * _radius
      );
      setTriangleSearchNode( connectedNode_1, ETP::SearchNode<P,D>( lvl3_1_hVal, lvl3_1_gValue, _startTri, lvl3EndPt_1_To2_idx ), usedList );
      startPoints.push_back( connectedNode_1 );
      if( connectedNode_2 != nullptr ){
        int lvl3EndPt_2_To2_idx = _startTri->getIndexfromConnectedNode( connectedNode_2_idx );
        std::shared_ptr<ETP::Edge<P,D>> lvl3_2_entryEdge = connectedNode_2->getEdgeVal( lvl3EndPt_2_To2_idx );
        D lvl3_2_hVal = lvl3_2_entryEdge->minimumDistance( _goalPoint, _scale ) / _scale;
        D lvl3_2_gValue = std::max(
          lvl3_2_entryEdge->minimumDistance( _startPoint, _scale ) / _scale,
          getPointToEdgeLowerBound( _startPoint, _startTri->getEdgeVal( connectedNode_2_idx ), _radius, _scale ) + _startTri->getLowerBound( connectedNode_2_idx )  * _radius
        );
        setTriangleSearchNode( connectedNode_2, ETP::SearchNode<P,D>( lvl3_2_hVal, lvl3_2_gValue, _startTri, lvl3EndPt_2_To2_idx ), usedList );
        startPoints.push_back( connectedNode_2 );
      }
    }else if( sLvl == 1 ){
      std::shared_ptr<ETP::Triangle<P,D>> startLvl2 = _startTri->getLvl1RootNode();
      int lvl1To2Idx = _startTri->getConnectedNodeIndex( startLvl2 );
      int lvl2To1Idx = _startTri->getIndexfromConnectedNode( lvl1To2Idx );
      std::shared_ptr<ETP::Edge<P,D>> lvl2entryEdge = startLvl2->getEdgeVal( lvl2To1Idx );
      D lvl2HVal = lvl2entryEdge->minimumDistance( _goalPoint, _scale ) / _scale;
      D lvl1To2GVal = std::max( getPointToEdgeLowerBound( _startPoint,
                                                          _startTri->getEdgeVal( lvl1To2Idx ), _radius, _scale ) + _startTri->getLowerBound( lvl1To2Idx ) * _radius,
                                                          std::max(
                                                            lvl2entryEdge->minimumDistance( _startPoint, _scale ) / _scale,
                                                            lvl2HVal - _startPoint.dist( _goalPoint )
                                                          ));
      setTriangleSearchNode( startLvl2, ETP::SearchNode<P,D>( lvl2HVal, lvl1To2GVal, _startTri, lvl2To1Idx ), usedList );
      std::shared_ptr<ETP::Triangle<P,D>> connectedNode_1 = nullptr;
      std::shared_ptr<ETP::Triangle<P,D>> connectedNode_2 = nullptr;
      int connectedNode_1_idx, connectedNode_2_idx;
      for( int i = 0; i < 3; i++ ){
        std::shared_ptr<ETP::Triangle<P,D>> connectedNode = startLvl2->getConnectedNode( i );
        if( connectedNode == nullptr ) continue;
        if( connectedNode_1 == nullptr ){
          connectedNode_1 = connectedNode;
          connectedNode_1_idx = i;
        }else{
          if( connectedNode != connectedNode_1 ){
            connectedNode_2 = connectedNode;
            connectedNode_2_idx = i;
          }else if( startLvl2->getLowerBound( connectedNode_1_idx ) > startLvl2->getLowerBound( i ) ){
            connectedNode_1 = connectedNode;
            connectedNode_1_idx = i;
          }
        }
      }
      int lvl3EndPt_1_To2_idx = startLvl2->getIndexfromConnectedNode( connectedNode_1_idx );
      std::shared_ptr<ETP::Edge<P,D>> lvl3_1_entryEdge = connectedNode_1->getEdgeVal( lvl3EndPt_1_To2_idx );
      D lvl3_1_hVal = lvl3_1_entryEdge->minimumDistance( _goalPoint, _scale ) / _scale;
      D lvl3_1_gValue = std::max(
        lvl3_1_entryEdge->minimumDistance( _startPoint, _scale ) / _scale,
        std::max(
          lvl1To2GVal + ( startLvl2->getAngle( lvl2entryEdge, startLvl2->getEdgeVal( connectedNode_1_idx ) ) + startLvl2->getLowerBound( connectedNode_1_idx ) ) * _radius,
          lvl1To2GVal + ( lvl3_1_hVal - lvl2HVal )
        )
      );
      setTriangleSearchNode( connectedNode_1, ETP::SearchNode<P,D>( lvl3_1_hVal, lvl3_1_gValue, startLvl2, lvl3EndPt_1_To2_idx ), usedList );
      startPoints.push_back( connectedNode_1 );
      if( connectedNode_2 != nullptr ){
        int lvl3EndPt_2_To2_idx = startLvl2->getIndexfromConnectedNode( connectedNode_2_idx );
        std::shared_ptr<ETP::Edge<P,D>> lvl3_2_entryEdge = connectedNode_2->getEdgeVal( lvl3EndPt_2_To2_idx );
        D lvl3_2_hVal = lvl3_2_entryEdge->minimumDistance( _goalPoint, _scale ) / _scale;
        D lvl3_2_gValue = std::max(
          lvl3_2_entryEdge->minimumDistance( _startPoint, _scale ) / _scale,
          std::max(
            lvl1To2GVal + ( startLvl2->getAngle( lvl2entryEdge, startLvl2->getEdgeVal( connectedNode_2_idx ) ) + startLvl2->getLowerBound( connectedNode_2_idx ) ) * _radius,
            lvl1To2GVal + ( lvl3_2_hVal - lvl2HVal )
          )
        );
        setTriangleSearchNode( connectedNode_2, ETP::SearchNode<P,D>( lvl3_2_hVal, lvl3_2_gValue, startLvl2, lvl3EndPt_2_To2_idx ), usedList );
        startPoints.push_back( connectedNode_2 );
      }
    }
    ETP::searchResult<P,D> classicSearch = searchBetweenLvl3EndPoints( _startPoint, _goalPoint, startPoints, goalEndPoints, std::numeric_limits<D>::infinity(), std::numeric_limits<D>::infinity(), _radius, _scale, usedList );

    if( classicSearch.goalEndPoint != nullptr ){
      setSearchNodesFromEndPoint( classicSearch.goalEndPoint, usedList, true );
      std::vector<std::shared_ptr<ETP::Triangle<P,D>>> ret = getFunnel( _goalTri );
      ret[ ret.size() - 1 ]->setWidth( 0, classicSearch.goalEndPoint->gValue );
      return getFunnel( _goalTri );
    }else{
      return std::vector<std::shared_ptr<ETP::Triangle<P,D>>>();
    }

  }

  ETP::searchResult<P,D> searchBetweenLvl3EndPoints( ETP::Point<P,D> _startPoint, ETP::Point<P,D> _goalPoint, std::vector<std::shared_ptr<ETP::Triangle<P,D>>> _startPoints, std::vector<std::shared_ptr<ETP::searchEndPoint<P,D>>> _goalEndPoints, D _maxGValueToEndPoint, D _maxGValueToGoal, D _radius, D _scale, std::shared_ptr<ETP::UsedList<P,D>> usedList ){
    std::vector<std::shared_ptr<ETP::Triangle<P,D>>> openList;
    for( std::shared_ptr<ETP::Triangle<P,D>> startTri : _startPoints ){
      insertInOpenList( startTri->getSearchNode()->hValue, startTri, openList );
    }
    std::shared_ptr<ETP::searchEndPoint<P,D>> closestEndPoint = nullptr;
    D memLowerGValue;
    bool endPointFound = false;
    while( openList.empty() == false ){
      std::shared_ptr<ETP::Triangle<P,D>> tCurrent = openList.front();  openList.erase( openList.begin() );
      std::shared_ptr<ETP::SearchNode<P,D>> searchNodeCurrent = tCurrent->getSearchNode();
      D tCGVal = searchNodeCurrent->gValue;
      int comeFromIdx = searchNodeCurrent->comeFromIdx;
      for( int i = 0; i < 3; i++ ){
        if( i == comeFromIdx ) continue;
        std::shared_ptr<ETP::Triangle<P,D>> tNeighbour = tCurrent->getConnectedNode( i );
        if( tNeighbour == nullptr ) continue;
        std::shared_ptr<ETP::Edge<P,D>> tCurrentOuterEdge = tCurrent->getEdgeVal( i );
        int neighbourToCurrentIdx = tCurrent->getIndexfromConnectedNode( i );
        std::shared_ptr<ETP::Edge<P,D>> neighbEntryEdge = tNeighbour->getEdgeVal( neighbourToCurrentIdx );
        D tCurrentToNeighbLowerBound = tCGVal + tCurrent->getLowerBound( i ) * _radius;
        D startToClosestEdgeDist = neighbEntryEdge->minimumDistance( _startPoint, _radius, _scale ) / _scale;
        D hVal = neighbEntryEdge->minimumDistance( _goalPoint, _radius, _scale ) / _scale;
        D hValDiff = tCGVal + ( searchNodeCurrent->hValue - hVal );
        if( searchNodeCurrent->comeFrom == nullptr ){
          tCurrentToNeighbLowerBound += getPointToEdgeLowerBound(  _startPoint, tCurrentOuterEdge, _radius, _scale );
        }else /*if( tCurrent->getLevel() == 3 )*/{
          //tCurrentToNeighbLowerBound += tCurrent->getAngle( tCurrent->getEdgeVal( searchNodeCurrent->comeFromIdx ), tCurrentOuterEdge ) * _radius;
          tCurrentToNeighbLowerBound += getAngle( tCurrent->getEdgeVal( searchNodeCurrent->comeFromIdx ), tCurrentOuterEdge ) * _radius;
        }
        D endPointLowerBound = 0.0;
        D endPointGValToGoal = 0.0;
        std::shared_ptr<ETP::searchEndPoint<P,D>> goalEndPt = nullptr;
        for( std::shared_ptr<ETP::searchEndPoint<P,D>> gEndPt : _goalEndPoints ){
          if( gEndPt->triangle == tNeighbour ){
            goalEndPt = gEndPt;
            //endPointLowerBound += ( goalEndPt->entryEdge != nullptr ? tNeighbour->getAngle( neighbEntryEdge, goalEndPt->entryEdge ) * _radius : getPointToEdgeLowerBound( _goalPoint, neighbEntryEdge, _radius, _scale ) );
            endPointLowerBound += ( goalEndPt->entryEdge != nullptr ? getAngle( neighbEntryEdge, goalEndPt->entryEdge ) * _radius : getPointToEdgeLowerBound( _goalPoint, neighbEntryEdge, _radius, _scale ) );
            endPointGValToGoal += goalEndPt->gValueToSearchStart;

            /*
            openList.erase( std::remove_if( openList.begin(), openList.end(),
              [ gVal ]( std::shared_ptr<ETP::Triangle<P,D>> _t ){
                  if ( _t->getSearchNode()->gValue >= gVal ) {
                      return true;
                  }
                  return false;
              }
            ), openList.end() );
            */

            break;
          }
        }
        D gVal = tCurrentToNeighbLowerBound + endPointLowerBound;
        //D gVal = std::max( hValDiff, std::max( startToClosestEdgeDist, tCurrentToNeighbLowerBound + endPointLowerBound ) );
        D gValToGoal = gVal + endPointGValToGoal;
        if( ( endPointFound == true &&  gValToGoal >= memLowerGValue ) || gValToGoal >= _maxGValueToGoal ) continue;
        if( goalEndPt != nullptr ){
          if( endPointFound == false || goalEndPt->gValue > gVal ){
            goalEndPt->gValue = gValToGoal;//gVal;
            endPointFound = true;
            memLowerGValue = gValToGoal;//gVal;

          }
        }
        std::shared_ptr<ETP::SearchNode<P,D>> tNeighbSearchNode = tNeighbour->getSearchNode();
        if( tNeighbSearchNode == nullptr ){
          setTriangleSearchNode( tNeighbour, ETP::SearchNode<P,D>( hVal, gVal, tCurrent, neighbourToCurrentIdx ), usedList );
          if( goalEndPt == nullptr ){
            insertInOpenList( hVal, tNeighbour, openList );
          }else if( closestEndPoint == nullptr || gValToGoal < closestEndPoint->gValue ){
            closestEndPoint = goalEndPt;
            // return searchResult<P,D>( goalEndPt, usedList );
            //return searchResult<P,D>( closestEndPoint, usedList );
            //closestEndPoint->gValue = gVal;
          }
        }else if( tNeighbSearchNode->gValue > gVal ){
          //setTriangleSearchNode( tNeighbour, ETP::SearchNode<P,D>( hVal, gVal, tCurrent, neighbourToCurrentIdx ), usedList );
          tNeighbSearchNode->hValue = hVal;
          tNeighbSearchNode->gValue = gVal;
          tNeighbSearchNode->comeFrom = tCurrent;
          tNeighbSearchNode->comeFromIdx = neighbourToCurrentIdx;
          if( goalEndPt == nullptr ){
            removeFromOpenList( tNeighbour, openList );
            insertInOpenList( hVal, tNeighbour, openList );
          }else if( ( goalEndPt != nullptr && closestEndPoint == nullptr ) || gValToGoal < closestEndPoint->gValue ){
            closestEndPoint = goalEndPt;
            //closestEndPoint->gValue = gVal;
          }
        }
      }
    }
    return searchResult<P,D>( closestEndPoint, usedList );
  }


  std::vector<std::shared_ptr<ETP::Triangle<P,D>>> searchInsideLvl2RingOrLoop( ETP::Point<P,D> _startPoint, ETP::Point<P,D> _goalPoint, std::shared_ptr<ETP::Triangle<P,D>> _startTri, std::shared_ptr<ETP::Triangle<P,D>> _goalTri, D _radius, D _scale ){
    std::shared_ptr<ETP::UsedList<P,D>> usedList = std::make_shared<ETP::UsedList<P,D>>();
    D memGValue;
    std::shared_ptr<ETP::Triangle<P,D>> tCurrent = _startTri;
    if( _startTri->getLevel() == 1 ){
      std::shared_ptr<ETP::Triangle<P,D>> lvl2Root = _startTri->getLvl1RootNode();
      int startToLvl2Idx = _startTri->getConnectedNodeIndex( lvl2Root );
      //memGValue += getPointToEdgeLowerBound( _startPoint, _startTri->getEdgeVal( startToLvl2Idx ), _radius, _scale ) + _startTri->getLowerBound( startToLvl2Idx ) * radius;
      setTriangleSearchNode( lvl2Root, ETP::SearchNode<P,D>( 0.0, getPointToEdgeLowerBound( _startPoint, _startTri->getEdgeVal( startToLvl2Idx ), _radius, _scale ) + _startTri->getLowerBound( startToLvl2Idx ) * _radius, tCurrent, tCurrent->getIndexfromConnectedNode( startToLvl2Idx ) ), usedList );
      tCurrent = lvl2Root;
    }else{
      setTriangleSearchNode( _startTri, ETP::SearchNode<P,D>( 0.0, 0.0, nullptr ), usedList );
    }
    bool endIsLvl1;
    std::shared_ptr<ETP::Triangle<P,D>> lvl2EndPt;
    std::shared_ptr<ETP::Edge<P,D>> lv2EndEntry;
    D lvl1EndToRootLowerBound;
    if( _goalTri->getLevel() == 1 ){
      endIsLvl1 = true;
      lvl2EndPt = _goalTri->getLvl1RootNode();
      int idxToRoot = _goalTri->getConnectedNodeIndex( lvl2EndPt );
      int idxRootToEnd = _goalTri->getIndexfromConnectedNode( idxToRoot );
      std::shared_ptr<ETP::Edge<P,D>> lv2EndEntry = lvl2EndPt->getEdgeVal( idxRootToEnd );
      lvl1EndToRootLowerBound = _goalTri->getLowerBound( idxToRoot ) * _radius;
      setTriangleSearchNode( _goalTri, ETP::SearchNode<P,D>( 0.0, 0.0, lvl2EndPt, idxRootToEnd ), usedList );
    }else{
      endIsLvl1 = false;
      lvl2EndPt = _goalTri;
    }
    ETP::Point<D,D> goalPoint =  ETP::Point<D,D>( (D) _goalPoint.x, (D) _goalPoint.y );
    std::vector<std::shared_ptr<ETP::Triangle<P,D>>> openList;
    openList.push_back( tCurrent );
    bool goalReached = false;
    while( openList.empty() == false ){
      std::shared_ptr<ETP::Triangle<P,D>> tCurrent = openList.front();
      openList.erase( openList.begin() );
      ETP::SearchNode<P,D> snCurrent = *tCurrent->getSearchNode().get();
      int currentLvl = tCurrent->getLevel();
      for( int i = 0; i < 3; i++ ){
        if( currentLvl == 3 && tCurrent->getConnectedNode( i ) != tCurrent ){ continue; }
        std::shared_ptr<ETP::Triangle<P,D>> tNeighbour = tCurrent->getAdjacentVal( i );
        if( tNeighbour == nullptr || tNeighbour->getLevel() < 2 || tNeighbour->getSearchNode() != nullptr ){ continue; }
        D neighbGValue = snCurrent.gValue + tCurrent->getAngle( tCurrent->getEdgeVal( snCurrent.comeFromIdx ), tCurrent->getEdgeVal( i ) ) * _radius;
        if( goalReached == true && neighbGValue >= memGValue ){
          return getFunnel( _goalTri );
        }else if( tNeighbour == lvl2EndPt ){
          D gValToGoal = neighbGValue;
          if( endIsLvl1 ){
            gValToGoal += tNeighbour->getAngle( tNeighbour->getEdgeVal( tNeighbour->getAdjacentIndex( tCurrent ) ), lv2EndEntry ) * _radius + lvl1EndToRootLowerBound;
          }else{
            gValToGoal += getPointToEdgeLowerBound( _goalPoint, tNeighbour->getEdgeVal( tNeighbour->getAdjacentIndex( tCurrent ) ), _radius, _scale );
          }
          if( goalReached == false ){
            memGValue = gValToGoal;
            setTriangleSearchNode( tNeighbour, ETP::SearchNode<P,D>( 0.0, memGValue, tCurrent ), usedList );
            goalReached = true;
          }else if( gValToGoal < memGValue ){
            setTriangleSearchNode( tNeighbour, ETP::SearchNode<P,D>( 0.0, gValToGoal, tCurrent ), usedList );
            return getFunnel( _goalTri );
          }
        }else{
          D neighbHValue = tNeighbour->getCentroid().dist( goalPoint );
          setTriangleSearchNode( tNeighbour, ETP::SearchNode<P,D>( neighbHValue, neighbGValue, tCurrent ), usedList );
          insertInOpenList( neighbHValue, tNeighbour, openList );
        }
      }
    }
    return getFunnel( _goalTri );
  }

  std::vector<std::shared_ptr<ETP::Triangle<P,D>>> searchInsideLvl2Corridor( ETP::Point<P,D> _startPoint, ETP::Point<P,D> _goalPoint, std::shared_ptr<ETP::Triangle<P,D>> _startTri, std::shared_ptr<ETP::Triangle<P,D>> _goalTri, D _radius, D _scale ){
    std::shared_ptr<ETP::UsedList<P,D>> usedList = std::make_shared<ETP::UsedList<P,D>>();
    std::shared_ptr<ETP::Triangle<P,D>> startLvl2,
                                        goalLvl2,
                                        startLvl3,
                                        goalLvl3;
    std::shared_ptr<ETP::Edge<P,D>> startLvl1ToLvl2Edge, startLvl2ToLvl1Edge, startLvl2ToStartLvl3Edge, startLvl2ToGoalLvl3Edge,
                                    goalLvl1ToLvl2Edge, goalLvl2ToLvl1Edge, goalLvl2ToStartLvl3Edge, goalLvl2ToGoalLvl3Edge,
                                    startLvl3ToCorridorEdge,
                                    goalLvl3ToCorridorEdge;
    int lvl1StartToLvl2Idx, lvl2ToLvl1StartIdx,
        lvl2StartToLvl3StartIdx, lvl2StartToLvl3GoalIdx,
        lvl1GoalToLvl2Idx, lvl2ToLvl1GoalIdx,
        lvl2GoalToLvl3StartIdx, lvl2GoalToLvl3GoalIdx,
        startLvl3ToCorridorIdx, goalLvl3ToCorridorIdx;
    if( _startTri->getLevel() == 1 ){
      startLvl2 = _startTri->getLvl1RootNode();
      lvl1StartToLvl2Idx = _startTri->getConnectedNodeIndex( startLvl2 );
      lvl2ToLvl1StartIdx = _startTri->getIndexfromConnectedNode( lvl1StartToLvl2Idx );
      startLvl1ToLvl2Edge = _startTri->getEdgeVal( lvl1StartToLvl2Idx );
      startLvl2ToLvl1Edge = startLvl2->getEdgeVal( lvl2ToLvl1StartIdx );
      D hVal =  startLvl2ToLvl1Edge->minimumDistance( _goalPoint, _scale ) / _scale;
      D hValDiff = _startPoint.dist( _goalPoint ) / _scale - hVal;
      D dist = startLvl2ToLvl1Edge->minimumDistance( _startPoint, _scale ) / _scale;
      D lowerBound = getPointToEdgeLowerBound( _startPoint, startLvl1ToLvl2Edge, _radius, _scale ) + _startTri->getLowerBound( lvl1StartToLvl2Idx ) * _radius;
      D gVal = std::max( hValDiff, std::max( dist, lowerBound ) );
      setTriangleSearchNode( startLvl2, ETP::SearchNode<P,D>( hVal, gVal, _startTri, lvl2ToLvl1StartIdx ), usedList );
    }else{
      startLvl2 = _startTri;
      startLvl2ToLvl1Edge = nullptr;
    }
    if( _goalTri->getLevel() == 1 ){
      goalLvl2 = _goalTri->getLvl1RootNode();
      lvl1GoalToLvl2Idx = _goalTri->getConnectedNodeIndex( goalLvl2 );
      lvl2ToLvl1GoalIdx = _goalTri->getIndexfromConnectedNode( lvl1GoalToLvl2Idx );
      goalLvl1ToLvl2Edge = _goalTri->getEdgeVal( lvl1GoalToLvl2Idx );
      goalLvl2ToLvl1Edge = goalLvl2->getEdgeVal( lvl2ToLvl1GoalIdx );
    }else{
      goalLvl2 = _goalTri;
      goalLvl2ToLvl1Edge = nullptr;
    }
    std::shared_ptr<ETP::Triangle<P,D>> lvl3_1 = nullptr;
    std::shared_ptr<ETP::Triangle<P,D>> lvl3_2;
    for( int i = 0; i < 3; i++ ){
      std::shared_ptr<ETP::Triangle<P,D>> lv3ConnectedNode = startLvl2->getConnectedNode( i );
      if( lv3ConnectedNode == nullptr ) continue;
      if( lvl3_1 == nullptr ){
        lvl3_1 = lv3ConnectedNode;
      }else{
        lvl3_2 = lv3ConnectedNode;
        break;
      }
    }
    if( startLvl2->getLowerBound( startLvl2->getConnectedNodeIndex( lvl3_1 ) ) > goalLvl2->getLowerBound( goalLvl2->getConnectedNodeIndex( lvl3_1 ) ) ){
      goalLvl3 = lvl3_1;
      startLvl3 = lvl3_2;
    }else{
      goalLvl3 = lvl3_2;
      startLvl3 = lvl3_1;
    }

    lvl2StartToLvl3StartIdx = startLvl2->getConnectedNodeIndex( startLvl3 );
    lvl2StartToLvl3GoalIdx = startLvl2->getConnectedNodeIndex( goalLvl3 );
    lvl2GoalToLvl3StartIdx = goalLvl2->getConnectedNodeIndex( startLvl3 );
    lvl2GoalToLvl3GoalIdx = goalLvl2->getConnectedNodeIndex( goalLvl3 );
    startLvl3ToCorridorIdx = startLvl2->getIndexfromConnectedNode( lvl2StartToLvl3StartIdx );
    goalLvl3ToCorridorIdx = startLvl2->getIndexfromConnectedNode( lvl2StartToLvl3GoalIdx );

    startLvl2ToStartLvl3Edge = startLvl2->getEdgeVal( lvl2StartToLvl3StartIdx );
    startLvl2ToGoalLvl3Edge = startLvl2->getEdgeVal( lvl2StartToLvl3GoalIdx );
    goalLvl2ToStartLvl3Edge = goalLvl2->getEdgeVal( lvl2GoalToLvl3StartIdx );
    goalLvl2ToGoalLvl3Edge = goalLvl2->getEdgeVal( lvl2GoalToLvl3GoalIdx );

    startLvl3ToCorridorEdge = startLvl3->getEdgeVal( startLvl3ToCorridorIdx );
    goalLvl3ToCorridorEdge = goalLvl3->getEdgeVal( goalLvl3ToCorridorIdx );

    std::shared_ptr<ETP::Triangle<P,D>> tCurrent = startLvl2;
    std::shared_ptr<ETP::Edge<P,D>> edgeIn = startLvl2ToLvl1Edge;
    D memGValue = 0.0;
    for(;;){
      if( tCurrent == goalLvl2 ) break;
      int idxToNeighb = tCurrent->getConnectedNodeIndex( goalLvl3 );
      std::shared_ptr<ETP::Triangle<P,D>> tNeighbour = tCurrent->getAdjacentVal( idxToNeighb );
      std::shared_ptr<ETP::Edge<P,D>> edgeOut = tCurrent->getEdgeVal( idxToNeighb );
      D startToEdgeOutDist = edgeOut->minimumDistance( _startPoint, _scale ) / _scale;
      D hVal = edgeOut->minimumDistance( _goalPoint, _scale ) / _scale;
      D lowerBound = 0.0;
      D hValDiff;
      if( tNeighbour == goalLvl2 && goalLvl2ToLvl1Edge == nullptr ) lowerBound += getPointToEdgeLowerBound( _goalPoint, goalLvl2ToStartLvl3Edge, _radius, _scale );
      if( edgeIn != nullptr ){
        lowerBound += memGValue + tCurrent->getAngle( edgeIn, edgeOut ) * _radius;
        hValDiff = memGValue + ( tCurrent->getSearchNode()->hValue - hVal );
      }else{
        lowerBound += getPointToEdgeLowerBound( _startPoint, edgeOut, _radius, _scale );
        hValDiff = _startPoint.dist( _goalPoint ) - hVal;
      }
      D gValue = std::max( hValDiff, std::max( startToEdgeOutDist, lowerBound ) );
      setTriangleSearchNode( tNeighbour, ETP::SearchNode<P,D>( hVal, gValue, tCurrent, tNeighbour->getAdjacentIndex( tCurrent ) ), usedList );
      memGValue = gValue;
      edgeIn = edgeOut;
      tCurrent = tNeighbour;
    }
    std::shared_ptr<ETP::searchEndPoint<P,D>> goalEndPoint = nullptr;
    D innerCorridorGValue, startCorridorGValue, goalCorridorGValue;;
    if( goalLvl2ToLvl1Edge != nullptr ){
      D dist = _goalPoint.dist( _startPoint ) / _scale;
      D lvl2To1LowerBound = getPointToEdgeLowerBound( _goalPoint, goalLvl1ToLvl2Edge, _radius, _scale ) + _goalTri->getLowerBound( lvl1StartToLvl2Idx ) * _radius;
      std::shared_ptr<ETP::Triangle<P,D>> comeFrom = goalLvl2->getSearchNode()->comeFrom;
      D parentGVal = comeFrom->getSearchNode() != nullptr ? comeFrom->getSearchNode()->gValue : getPointToEdgeLowerBound( _startPoint, goalLvl2->getEdgeVal( goalLvl2->getAdjacentIndex( comeFrom ) ), _radius, _scale );
      innerCorridorGValue = std::max( dist, ( parentGVal + lvl2To1LowerBound + goalLvl2->getAngle( goalLvl2ToLvl1Edge, goalLvl2ToStartLvl3Edge ) * _radius ) );
      setTriangleSearchNode( _goalTri, ETP::SearchNode<P,D>( 0.0, innerCorridorGValue, goalLvl2, _goalTri->getConnectedNodeIndex( goalLvl2 ) ), usedList );

      //D lvl2Tolvl3LowerBound = ( goalLvl2->getAngle( goalLvl2ToLvl1Edge, goalLvl2ToGoalLvl3Edge ) + goalLvl2->getLowerBound( lvl2GoalToLvl3GoalIdx ) ) * _radius;
      goalCorridorGValue = lvl2To1LowerBound + ( goalLvl2->getAngle( goalLvl2ToLvl1Edge, goalLvl2ToGoalLvl3Edge ) + goalLvl2->getLowerBound( lvl2GoalToLvl3GoalIdx ) ) * _radius;
      //goalCorridorGValue = lvl2To1LowerBound + lvl2Tolvl3LowerBound;
      std::shared_ptr<ETP::searchEndPoint<P,D>> lvl1EndPoint = std::make_shared<ETP::searchEndPoint<P,D>>( 0.0, _goalTri, nullptr );
      std::shared_ptr<ETP::searchEndPoint<P,D>> lvl2EndPoint = std::make_shared<ETP::searchEndPoint<P,D>>( 0.0, goalLvl2, lvl1EndPoint, lvl2ToLvl1GoalIdx );
      goalEndPoint = std::make_shared<ETP::searchEndPoint<P,D>>( 0.0, goalLvl3, lvl2EndPoint, goalLvl3ToCorridorIdx, goalCorridorGValue );
    }else{
      //D lvl2Tolvl3LowerBound = getPointToEdgeLowerBound( _goalPoint, goalLvl2ToGoalLvl3Edge, _radius, _scale ) + goalLvl2->getLowerBound( lvl2GoalToLvl3GoalIdx ) * _radius;
      goalCorridorGValue = getPointToEdgeLowerBound( _goalPoint, goalLvl2ToGoalLvl3Edge, _radius, _scale ) + goalLvl2->getLowerBound( lvl2GoalToLvl3GoalIdx ) * _radius;
      //goalCorridorGValue = lvl2Tolvl3LowerBound;
      std::shared_ptr<ETP::searchEndPoint<P,D>> lvl2EndPoint = std::make_shared<ETP::searchEndPoint<P,D>>( 0.0, goalLvl2, nullptr, lvl2GoalToLvl3GoalIdx );
      goalEndPoint = std::make_shared<ETP::searchEndPoint<P,D>>( 0.0, goalLvl3, lvl2EndPoint, goalLvl3ToCorridorIdx, goalCorridorGValue );
      innerCorridorGValue = _goalTri->getSearchNode()->gValue;
    }
    if( startLvl2ToLvl1Edge != nullptr ){
      D dist = startLvl3ToCorridorEdge->minimumDistance( _startPoint, _scale ) / _scale;
      D lvl2Tolvl3LowerBound = ( startLvl2->getAngle( startLvl2ToLvl1Edge, startLvl2ToStartLvl3Edge ) + startLvl2->getLowerBound( lvl2StartToLvl3StartIdx ) ) * _radius;
      D lvl2GValue = startLvl2->getSearchNode()->gValue;
      D startToLvl3LowerBound = lvl2GValue + lvl2Tolvl3LowerBound;
      D hVal = startLvl3ToCorridorEdge->minimumDistance( _goalPoint, _scale ) / _scale;
      D hValDiff = lvl2GValue + ( startLvl2ToLvl1Edge->minimumDistance( _goalPoint, _scale ) / _scale ) - hVal;
      startCorridorGValue = std::max( startToLvl3LowerBound, std::max( dist, hValDiff ) );
      setTriangleSearchNode( startLvl3, ETP::SearchNode<P,D>( hVal, startCorridorGValue, startLvl2, startLvl2->getIndexfromConnectedNode( startLvl2->getConnectedNodeIndex( startLvl3 ) ) ), usedList );
    }else{
      D dist = startLvl3ToCorridorEdge->minimumDistance( _startPoint, _scale ) / _scale;
      D lvl2Tolvl3LowerBound = getPointToEdgeLowerBound( _startPoint, startLvl2ToStartLvl3Edge, _radius, _scale ) + startLvl2->getLowerBound( lvl2StartToLvl3StartIdx ) * _radius;
      startCorridorGValue = std::max( lvl2Tolvl3LowerBound, dist );
      setTriangleSearchNode( startLvl3, ETP::SearchNode<P,D>( dist, startCorridorGValue, startLvl2, startLvl3ToCorridorIdx ), usedList );
    }
    if( startCorridorGValue + goalCorridorGValue >= innerCorridorGValue ) return getFunnel( _goalTri );
    D maxSearchGValue = innerCorridorGValue - goalCorridorGValue;
    ETP::searchResult<P,D> classicSearch = searchBetweenLvl3EndPoints( _startPoint, _goalPoint, std::vector<std::shared_ptr<ETP::Triangle<P,D>>>{ startLvl3 }, std::vector<std::shared_ptr<ETP::searchEndPoint<P,D>>>{ goalEndPoint }, maxSearchGValue, innerCorridorGValue, _radius, _scale, usedList );

    if( classicSearch.goalEndPoint != nullptr ) setSearchNodesFromEndPoint( goalEndPoint, usedList, true );
  //return std::vector<std::shared_ptr<ETP::Triangle<P,D>>>{ triangles[ 0 ] };
    //if( classicSearch.goalEndPoint != nullptr ) return getFunnel( classicSearch.goalEndPoint->triangle );

    return getFunnel( _goalTri );
  }


  std::vector<std::shared_ptr<ETP::searchEndPoint<P,D>>> getTriangleEndPoints( std::shared_ptr<ETP::Triangle<P,D>> _tri, ETP::Point<P,D> _goalPoint, D _radius, D _scale ){
  std::vector<std::shared_ptr<ETP::searchEndPoint<P,D>>> ret;
  int lvl = _tri->getLevel();
  std::shared_ptr<ETP::searchEndPoint<P,D>> finalEndPt = std::make_shared<ETP::searchEndPoint<P,D>>( 0.0, _tri );
  if( lvl == 3 ){
    ret.push_back( std::make_shared<ETP::searchEndPoint<P,D>>( 0.0, _tri ) );
  }else if( lvl == 2 ){
    for( int i = 0; i < 3; i++ ){
      std::shared_ptr<ETP::Triangle<P,D>> connectedNode = _tri->getConnectedNode( i );
      if( connectedNode != nullptr ){
        bool isNew = true;
        D lowerBound = getPointToEdgeLowerBound( _goalPoint, _tri->getEdgeVal( i ), _radius, _scale ) + _tri->getLowerBound( i ) * _radius;
      //  D gValue = std::max(  )
        for ( auto ePt : ret ){
          if( connectedNode == ePt->triangle ){
            isNew = false;
            if( lowerBound < ePt->gValue ){
              ePt = std::make_shared<ETP::searchEndPoint<P,D>>( lowerBound, connectedNode, finalEndPt, _tri->getIndexfromConnectedNode( i ) );
            }
          }
        }
        if( isNew == true ){
          ret.push_back( std::make_shared<ETP::searchEndPoint<P,D>>( lowerBound, connectedNode, finalEndPt, _tri->getIndexfromConnectedNode( i ) ) );
        }
      }
    }
  }else if( lvl == 1 ){
    for( int i = 0; i < 3; i++ ){
      std::shared_ptr<ETP::Triangle<P,D>> connectedNode = _tri->getConnectedNode( i );
      D startLowerBound = getPointToEdgeLowerBound( _goalPoint, _tri->getEdgeVal( i ), _radius, _scale );
      if( connectedNode != nullptr ){
        std::shared_ptr<ETP::searchEndPoint<P,D>> lvl2EndPt = std::make_shared<ETP::searchEndPoint<P,D>>( startLowerBound + _tri->getLowerBound( i ) * _radius, connectedNode, finalEndPt );
        for( int j = 0; j < 3; j++ ){
          std::shared_ptr<ETP::Triangle<P,D>> lvl3Node = connectedNode->getConnectedNode( j );
          if( lvl3Node != nullptr ){
            bool isNew = true;
            D lowerBound = connectedNode->getLowerBound( j ) * _radius;
            size_t retl = ret.size();
            for( int k = 0; k < retl; k++ ){
              std::shared_ptr<ETP::searchEndPoint<P,D>> ePt = ret[ k ];
              if( lvl3Node == ePt->triangle ){
                isNew = false;
                if( lowerBound < ePt->gValue ){
                  ret[ k ] = std::make_shared<ETP::searchEndPoint<P,D>>( lvl2EndPt->gValue + lowerBound, lvl3Node, lvl2EndPt, connectedNode->getIndexfromConnectedNode( j ) );
                }
                break;
              }
            }
            if( isNew == true ){
              ret.push_back( std::make_shared<ETP::searchEndPoint<P,D>>( lvl2EndPt->gValue + lowerBound, lvl3Node, lvl2EndPt, connectedNode->getIndexfromConnectedNode( j ) ) );
            }
          }
        }
        break;
      }
    }
  }
  return ret;
  }
  void setSearchNodesFromEndPoint( std::shared_ptr<ETP::searchEndPoint<P,D>> endPoint, std::shared_ptr<ETP::UsedList<P,D>> usedList, bool isGoalPart ){
    std::shared_ptr<ETP::searchEndPoint<P,D>> currEndPt = endPoint;

    if( isGoalPart == true ){
      for(;;){
        //if( currEndPt == nullptr ) break;
        std::shared_ptr<ETP::searchEndPoint<P,D>> prevEndPt = currEndPt->parent;
        if( prevEndPt == nullptr ) break;
        std::shared_ptr<ETP::Triangle<P,D>> prevTri = prevEndPt->triangle;
        int prevLvl = prevTri->getLevel();
        if( prevLvl == 2 ){
          int lvl2Tolvl3Idx = prevTri->getIndexWithConnectedNodeIdx( currEndPt->parentIndex,  currEndPt->triangle );
          setTriangleSearchNode( prevTri, ETP::SearchNode<P,D>( 0.0, 0.0, currEndPt->triangle, lvl2Tolvl3Idx ), usedList );
          currEndPt = prevEndPt;
        }else if( prevLvl == 1 ){
          int lvl1To2Idx = prevTri->getIndexWithConnectedNodeIdx( currEndPt->parentIndex,  currEndPt->triangle );
          setTriangleSearchNode( prevTri, ETP::SearchNode<P,D>( 0.0, 0.0, currEndPt->triangle, lvl1To2Idx ), usedList );
          break;
        }
      }
    }else{
      for(;;){
        //if( currEndPt == nullptr ) break;
        std::shared_ptr<ETP::searchEndPoint<P,D>> prevEndPt = currEndPt->parent;
        if( prevEndPt == nullptr ) break;
        std::shared_ptr<ETP::Triangle<P,D>> prevTri = prevEndPt->triangle;
        int prevLvl = prevTri->getLevel();
        if( prevLvl == 2 ){
          std::shared_ptr<ETP::Triangle<P,D>> nextLvl2 = currEndPt->triangle->getAdjacentVal( prevTri->getIndexfromConnectedNode( prevTri->getConnectedNodeIndex( currEndPt->triangle ) ) );
          setTriangleSearchNode( currEndPt->triangle, ETP::SearchNode<P,D>( 0.0, 0.0, nextLvl2 ), usedList );
          for(;;){
            if( nextLvl2 == prevTri ) break;
            std::shared_ptr<ETP::Triangle<P,D>> secLvl2;
            for( int i = 0; i < 3; i++ ){
              std::shared_ptr<ETP::Triangle<P,D>> tmpNextLvl2 = nextLvl2->getAdjacentVal( i );
              if( tmpNextLvl2 != nullptr && tmpNextLvl2->getSearchNode() == nullptr && tmpNextLvl2->getLevel() == 2 ){
                secLvl2 = tmpNextLvl2;
                break;
              }
            }
            setTriangleSearchNode( nextLvl2, ETP::SearchNode<P,D>( 0.0, 0.0, secLvl2 ), usedList );
            if( secLvl2 == prevTri ) break;
            nextLvl2 = secLvl2;
          }
          currEndPt = prevEndPt;
          //prevEndPt = currEndPt->parent;
        }else if( prevLvl == 1 ){
          setTriangleSearchNode( currEndPt->triangle, ETP::SearchNode<P,D>( 0.0, 0.0, prevTri ), usedList );
          break;
        }
      }
    }
  }



  std::vector<std::shared_ptr<ETP::Triangle<P,D>>> getTrianglesFromLvl1ToConnectedLvl2( std::shared_ptr<ETP::Triangle<P,D>> _startTri, int connectedNodeIdx ){
      std::vector<std::shared_ptr<ETP::Triangle<P,D>>> ret;
      std::shared_ptr<ETP::Triangle<P,D>> goal = _startTri->getConnectedNode( connectedNodeIdx );
      if( goal == nullptr ){ throw std::invalid_argument( "Triangle::getTrianglesToConnectedNode: search to null goal" ); }
      std::shared_ptr<ETP::Triangle<P,D>> tCurrent = _startTri;
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

  std::vector<std::shared_ptr<ETP::Triangle<P,D>>> searchLvl1Tree( std::shared_ptr<ETP::Triangle<P,D>> _startTri, std::shared_ptr<ETP::Triangle<P,D>> _goalTri ){
    std::vector<std::shared_ptr<ETP::Triangle<P,D>>> openList;
    //ETP::UsedList<P,D> usedList = ETP::UsedList<P,D>();
    std::shared_ptr<ETP::UsedList<P,D>> usedList = std::make_shared<ETP::UsedList<P,D>>();
    ETP::Point<D,D> goalCentroid = _goalTri->getCentroid();
    setTriangleSearchNode( _startTri,ETP::SearchNode<P,D>( 0.0, 0.0 ), usedList );
    openList.push_back( _startTri );
    bool goalReached = false;
    while( openList.empty() == false && goalReached == false ){
      std::shared_ptr<ETP::Triangle<P,D>> tCurrent = openList.front();
      openList.erase( openList.begin() );
      ETP::SearchNode<P,D> snCurrent = *tCurrent->getSearchNode().get();
      ETP::Point<D,D> currentCentroid = tCurrent->getCentroid();
      for( int i = 0; i < 3; i++ ){
        std::shared_ptr<ETP::Triangle<P,D>> tNeighbour = tCurrent->getAdjacentVal( i );
        if( tNeighbour == nullptr || tNeighbour->getLevel() != 1 ){ continue; }
        ETP::Point<D,D> tNeighbCentroid = tNeighbour->getCentroid();
        if( tNeighbour == _goalTri ){
          setTriangleSearchNode( tNeighbour, ETP::SearchNode<P,D>( 0.0, snCurrent.gValue + currentCentroid.dist( tNeighbCentroid ), tCurrent ), usedList );
          goalReached = true;
          break;
        }
        std::shared_ptr<ETP::SearchNode<P,D>> snNeighbourPtr = tNeighbour->getSearchNode();
        bool neighbIsNew = false;
        if( snNeighbourPtr == nullptr ){
          D tnHVal = tNeighbCentroid.dist( goalCentroid );
          setTriangleSearchNode( tNeighbour, ETP::SearchNode<P,D>( tnHVal, snCurrent.gValue + currentCentroid.dist( tNeighbCentroid ), tCurrent ), usedList );
          insertInOpenList( tnHVal, tNeighbour, openList );
        }
      }
    }
    if( goalReached == false ){
      return std::vector<std::shared_ptr<ETP::Triangle<P,D>>>();
    }else{
      return getFunnel( _goalTri );
    }
  }


  std::vector<std::shared_ptr<ETP::Triangle<P,D>>> getFunnel( std::shared_ptr<ETP::Triangle<P,D>> _goalTri ){
    //return std::vector<std::shared_ptr<ETP::Triangle<P,D>>> { triangles[ 0 ] };
    std::vector<std::shared_ptr<ETP::Triangle<P,D>>> ret;
    if( _goalTri->getLevel() == 2 || _goalTri->getLevel() == 3 ) ret.push_back( _goalTri );
    std::shared_ptr<ETP::Triangle<P,D>> tCurrent = _goalTri;
    for(;;){
      std::shared_ptr<ETP::SearchNode<P,D>> tSN = tCurrent->getSearchNode();
      if( tSN == nullptr ) break;
      std::shared_ptr<ETP::Triangle<P,D>> tComeFrom = tSN->comeFrom;
      if( tComeFrom == nullptr ) break;
      int tCurLvl = tCurrent->getLevel();
      int tComeFromLvl = tComeFrom->getLevel();
      if( ( tCurLvl == 3 && tComeFromLvl == 3 ) /*|| ( tCurLvl == 2 && tComeFromLvl == 3 )*/ ){

/*
        if( tComeFrom->getId() == 72 ){
          tCurrent->setLevel( tSN->comeFromIdx );
          return std::vector<std::shared_ptr<ETP::Triangle<P,D>>>{ tCurrent, tCurrent->getAdjacentVal( tSN->comeFromIdx ) ,tComeFrom };
        }
*/

        std::shared_ptr<ETP::Triangle<P,D>> tmpCurrent = tCurrent;
        for(;;){
          tmpCurrent = tmpCurrent == tCurrent ? tmpCurrent->getAdjacentVal( tSN->comeFromIdx ) :  tmpCurrent->getAdjacentVal( tmpCurrent->getConnectedNodeIndex( tComeFrom ) );
          ret.insert( ret.begin(), tmpCurrent );
          if( tmpCurrent ==  tComeFrom ) break;
        }
      }

      else if( tCurLvl == 2 && tComeFromLvl == 3 ){
        int lvl2To3Idx, lvl3To2Idx;
        if( tSN->comeFromIdx != -1 ){
          lvl2To3Idx = tSN->comeFromIdx;
          lvl3To2Idx = tCurrent->getIndexfromConnectedNode( lvl2To3Idx );
          //return std::vector<std::shared_ptr<ETP::Triangle<P,D>>>{ tCurrent, tCurrent->getAdjacentVal( lvl2To3Idx ) };
        }else{
          D memLowerBound;
          bool first = false;
          int tmpLvl2To3Idx, tmpLvl3To2Idx;
          for( int i = 0; i < 3; i++ ){
            if( tCurrent->getConnectedNode( i ) == tComeFrom ){
              if( first == false || tCurrent->getLowerBound( i ) < memLowerBound ){
                tmpLvl2To3Idx = i;
                tmpLvl3To2Idx = tCurrent->getIndexfromConnectedNode( i );
                memLowerBound = tCurrent->getLowerBound( i );
                first = true;
              }
            }
          }
          lvl2To3Idx = tmpLvl2To3Idx;
          lvl3To2Idx = tmpLvl3To2Idx;
        }

        std::shared_ptr<ETP::Triangle<P,D>> tmpCurrent = tCurrent->getAdjacentVal( lvl2To3Idx );
        for( int i = 0; ; i++ ){
          ret.insert( ret.begin(), tmpCurrent );
          if( tmpCurrent == tComeFrom ) break;

          for( int j = 0; j < 3; j++ ){
            if( tmpCurrent->getConnectedNode( j ) == tComeFrom && tmpCurrent->getIndexfromConnectedNode( j ) == lvl3To2Idx ){
              tmpCurrent = tmpCurrent->getAdjacentVal( j );
              break;
            }
          }
        }

      }



      else if( tCurLvl == 1 && tComeFromLvl == 2 ){
        for(;;){
          ret.insert( ret.begin(), tCurrent );
          if( tCurrent == tComeFrom ) break;
          tCurrent = tCurrent->getAdjacentVal( tCurrent->getConnectedNodeIndex( tComeFrom ) );
        }
      }else if( tCurLvl == 2 && tComeFromLvl == 2 ){
        ret.insert( ret.begin(), tComeFrom );
      }else if( tCurLvl == 3 && tComeFromLvl == 2 ){
        int lvl2To3Idx, lvl3To2Idx;
        if( tSN->comeFromIdx != -1 ){
          //lvl2To3Idx = tSN->comeFromIdx;
          lvl3To2Idx = tSN->comeFromIdx;//tComeFrom->getIndexfromConnectedNode( lvl2To3Idx );
        }else{
          D memLowerBound;
          bool first = false;
          int tmpLvl2To3Idx, tmpLvl3To2Idx;
          for( int i = 0; i < 3; i++ ){
            if( tComeFrom->getConnectedNode( i ) == tCurrent ){
              if( first == false || tComeFrom->getLowerBound( i ) < memLowerBound ){
                tmpLvl2To3Idx = i;
                tmpLvl3To2Idx = tComeFrom->getIndexfromConnectedNode( i );
                memLowerBound = tComeFrom->getLowerBound( i );
                first = true;
              }
            }
          }
          lvl2To3Idx = tmpLvl2To3Idx;
          lvl3To2Idx = tmpLvl3To2Idx;
        }

        std::shared_ptr<ETP::Triangle<P,D>> tmpCurrent = tComeFrom;
        for( int i = 0; ; i++ ){
          if( tmpCurrent == tCurrent ) break;
          ret.insert( ret.begin() + i , tmpCurrent );
          for( int j = 0; j < 3; j++ ){
            if( tmpCurrent->getConnectedNode( j ) == tCurrent && tmpCurrent->getIndexfromConnectedNode( j ) == lvl3To2Idx ){
              tmpCurrent = tmpCurrent->getAdjacentVal( j );
              break;
            }
          }
        }

      }else if( tCurLvl == 2 && tComeFromLvl == 1 ){
        std::shared_ptr<ETP::Triangle<P,D>> tmpCurrent = tComeFrom;
        for( int i = 0; ; i++ ){
          if( tCurrent == tmpCurrent ) break;
          ret.insert( ret.begin() + i , tmpCurrent );
          tmpCurrent = tmpCurrent->getAdjacentVal( tmpCurrent->getConnectedNodeIndex( tCurrent ) );
        }
      }
      tCurrent = ret.front();
    }
    return ret;
  }

  void setTriangleSearchNode( std::shared_ptr<ETP::Triangle<P,D>> _tri, ETP::SearchNode<P,D> _searchNode, std::shared_ptr<ETP::UsedList<P,D>> _usedList ){
    _tri->setSearchNode( _searchNode );
    _usedList->push_back( _tri );
  }

  void insertInOpenList( D _hVal, std::shared_ptr<ETP::Triangle<P,D>> _tri, std::vector<std::shared_ptr<ETP::Triangle<P,D>>>& _openList ){
    size_t ol = _openList.size();
    for( size_t i = 0; i < ol; i++ ){
      if( _hVal < _openList[ i ]->getSearchNode()->hValue ){
        _openList.insert( _openList.begin() + i, _tri );
        return void();
      }
    }
    _openList.push_back( _tri );
  }
  void removeFromOpenList( std::shared_ptr<ETP::Triangle<P,D>> _tri, std::vector<std::shared_ptr<ETP::Triangle<P,D>>>& _openList ){
    size_t ol = _openList.size();
    for( size_t i = 0; i < ol; i++ ){
      if( _openList[ i ] == _tri ){
        _openList.erase( _openList.begin() + i );
        return void();
      }
    }
  }

  D getPointToEdgeLowerBound(  ETP::Point<P,D> _p, std::shared_ptr<ETP::Edge<P,D>> _edge, D _radius, D _scale ){
    return _edge->minimumDistance( _p, _radius, _scale ) / _scale;

    ETP::Point<P,D> bv = *_edge->getPointVal( 0 ).get();
    ETP::Point<D,D> v = ETP::Point<D,D>( (D) bv.x  * _scale, (D) bv.y  * _scale );
    ETP::Point<P,D> bw = *_edge->getPointVal( 1 ).get();
    ETP::Point<D,D> w = ETP::Point<D,D>( (D) bw.x * _scale, (D) bw.y * _scale );
    ETP::Point<D,D> p = ETP::Point<D,D>( (D) _p.x, (D) _p.y );
    D l2 = ( ( w.x - v.x ) * ( w.x - v.x ) ) + ( ( w.y - v.y ) * ( w.y - v.y ) );
    if ( l2 == 0.0 ) return p.dist( v );
    ETP::Point<D,D> vp = ETP::Point<D,D>( p.x - v.x, p.y - v.y );
    ETP::Point<D,D> vw = ETP::Point<D,D>( w.x - v.x, w.y - v.y );
    D t = vp.dot( vw ) / l2;
    D dist1, dist2;
    if( t < 0.0 ){ //find close to w
      ETP::Point<D,D> wv = ETP::Point<D,D>( v.x - w.x, v.y - w.y );
      D wvL = std::sqrt( wv.x * wv.x + wv.y * wv.y );
      ETP::Point<D,D> nwv = ETP::Point<D,D>( wv.x / wvL, wv.y / wvL );
      ETP::Point<D,D> pwv = ETP::Point<D,D>( w.x + nwv.x * _radius, w.y + nwv.y * _radius );
    //  return p.dist( pwv );
      dist1 = p.dist( pwv );
    }
    if( t > 1.0 ){ // find close to v
      ETP::Point<D,D> vw = ETP::Point<D,D>( w.x - v.x, w.y - v.y );
      D vwL = std::sqrt( vw.x * vw.x + vw.y * vw.y );
      ETP::Point<D,D> nvw = ETP::Point<D,D>( vw.x / vwL, vw.y / vwL );
      ETP::Point<D,D> pvw = ETP::Point<D,D>( v.x + nvw.x * _radius, v.y + nvw.y * _radius );
      //return p.dist( pvw );
      dist2 = p.dist( pvw );
    }
  //  D t = std::max( 0.0, std::min( 1.0, vp.dot( vw ) ) / l2 );
    ETP::Point<D,D> projection = ETP::Point<D,D>( v.x + t * ( w.x - v.x ),  v.y + t * ( w.y - v.y ) );
    D vproj = v.dist( projection );
    D wproj = w.dist( projection );
    ETP::Point<D,D> fartherP;
    D longestD;
    if( vproj >= wproj ){
      fartherP = v;
      longestD = vproj;
    }else{
      fartherP = w;
      longestD = wproj;
    }
    if( longestD >= _radius ){
      return p.dist( projection );
    }else{
      ETP::Point<D,D> fp = ETP::Point<D,D>( projection.x - fartherP.x, projection.y - fartherP.y );
      D fpL = std::sqrt( fp.x * fp.x + fp.y * fp.y );
      ETP::Point<D,D> nfp = ETP::Point<D,D>( fp.x / fpL, fp.y / fpL );
      ETP::Point<D,D> pfp = ETP::Point<D,D>( fartherP.x + nfp.x * _radius, fartherP.y + nfp.y * _radius );
      return p.dist( pfp );
    }
  }

  /*
  float minimum_distance(vec2 v, vec2 w, vec2 p) {
      // Return minimum distance between line segment vw and point p
      const float l2 = length_squared(v, w);  // i.e. |w-v|^2 -  avoid a sqrt
      if (l2 == 0.0) return distance(p, v);   // v == w case
      // Consider the line extending the segment, parameterized as v + t (w - v).
      // We find projection of point p onto the line.
      // It falls where t = [(p-v) . (w-v)] / |w-v|^2
      // We clamp t from [0,1] to handle points outside the segment vw.
      const float t = max(0, min(1, dot(p - v, w - v) / l2));
      const vec2 projection = v + t * (w - v);  // Projection falls on the segment
      return distance(p, projection);
  }
  */
};// end of class
}// end of namespace




#endif
