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
  std::shared_ptr<ETP::Triangle<P,D>> getTriangleWithPoint( ETP::Point<P,D> _pt ){
  //  if( triangles[ 0 ]->isPointInside( _pt ) ){ return triangles[ 0 ]; }
  //  return nullptr;

    P sectorY = std::floor( ( _pt.y - minY ) / sectorHeight );
    if( sectorY >= (P) sectors.size() ){ return nullptr; }
    P sectorX = std::floor( ( _pt.x - minX ) / sectorWidth );
    if( sectorX >= (P) sectors[ sectorY ].size() ){ return nullptr; }
    std::vector<std::shared_ptr<ETP::Triangle<P,D>>> sector = sectors[ sectorY ][ sectorX ];
    size_t sl = sector.size();
    for( size_t i = 0; i < sl; i++ ){
      if( sector[ i ]->isPointInside( _pt ) ){ return sector[ i ]; }
    }
    return nullptr;
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
  void abstractLevel3BACKUP( std::shared_ptr<ETP::Triangle<P,D>> _t, int c ){
    std::vector<std::shared_ptr<ETP::Triangle<P,D>>> q = std::vector<std::shared_ptr<ETP::Triangle<P,D>>>();
    q.push_back( _t );
    while( q.empty() == false ){
      std::shared_ptr<ETP::Triangle<P,D>> t = q.front();
      q.erase( q.begin() );
      t->setLevel( 3 );
      t->setComponent( c );
      for( size_t j = 0; j < 3; j++ ){
        //int previousTid = t->getId();
        std::shared_ptr<ETP::Triangle<P,D>> t1stNeighbour = t->getAdjacentVal( j );
        D chokeCount = t->getEdgeVal( j )->getLength();
        D lowerBoundCount = 0.0f;
        D memChoke = 0.0f;
        //D memLowerBound = 0.0f;
        if( t1stNeighbour->getNumConstrainedEdges() + t1stNeighbour->getNumAdjacentLevel( 1 ) == 0 ){
          if( t1stNeighbour->getLevel() == -1 ){
            q.push_back( t1stNeighbour );
          }
          t1stNeighbour->setConnectedNode( t1stNeighbour->getAdjacentIndex( t ), t );
          continue;
        }else{
          t1stNeighbour->setLevel( 2 );
          t1stNeighbour->setComponent( c );
          int indexToPrevious = t1stNeighbour->getAdjacentIndex( t );
          t1stNeighbour->setConnectedNode( indexToPrevious, t );
        }
        std::shared_ptr<ETP::Triangle<P,D>> tCurrent = t1stNeighbour;
        std::shared_ptr<ETP::Edge<P,D>> edgeIn = t->getEdgeVal( j );
        bool reachedLvl3 = false;
        while( reachedLvl3 == false ){
          for( size_t k = 0; k < 3; k++ ){
            std::shared_ptr<ETP::Edge<P,D>> edgeOut = tCurrent->getEdgeVal( k );
            if( edgeOut->isConstrained() == true || edgeOut == edgeIn ){ continue; }
            std::shared_ptr<ETP::Triangle<P,D>> tCurrentNeighbour = tCurrent->getAdjacentVal( k );
            //if( tCurrentNeighbour->getId() == previousTid ){ continue; }
            if( tCurrentNeighbour->getLevel() == 1 ){
              collapseRootedTree( tCurrent, tCurrentNeighbour );
            }else{
              int neighbN = tCurrentNeighbour->getNumConstrainedEdges();
              int neighbM = tCurrentNeighbour->getNumAdjacentLevel( 1 );
            //  memLowerBound = getAngle( edgeIn, edgeOut );

              tCurrentNeighbour->setComponent( c );
              tCurrentNeighbour->setConnectedNode( tCurrentNeighbour->getAdjacentIndex( tCurrent ), t );
              memChoke = tCurrent->getWidthbetweenEdges( edgeIn, edgeOut );
              lowerBoundCount += getAngle( edgeIn, edgeOut );
              tCurrentNeighbour->setLowerBound( tCurrentNeighbour->getEdgeIndex( edgeOut ), lowerBoundCount );

              if( neighbN + neighbM == 0 ){
                if( tCurrentNeighbour->getLevel() == -1 ){
                  q.push_back( tCurrentNeighbour );
                }
                tCurrentNeighbour->setConnectedNode( tCurrentNeighbour->getAdjacentIndex( tCurrent ), t );
                //tCurrentNeighbour->setLowerBound( tCurrentNeighbour->getEdgeIndex( edgeOut ), lowerBoundCount + memLowerBound );
                reachedLvl3 = true;
                break;
              }
              //previousTid = tCurrent->getId();
              tCurrentNeighbour->setLevel( 2 );
              tCurrent = tCurrentNeighbour;
              edgeIn = edgeOut;
              chokeCount = std::min( chokeCount, memChoke );

            }
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
          //return std::vector<std::shared_ptr<ETP::Triangle<P,D>>>{ triangles[ 0 ] };
          bodySearch =  searchInsideLvl2RingOrLoop( _startPoint, _goalPoint, _startTri, _goalTri, _radius, _scale );
        }else if( startRoot->haveSameLvl2CorridorEndpoints( goalRoot )){// inside a corridor
          //return std::vector<std::shared_ptr<ETP::Triangle<P,D>>>{ _startTri, _goalTri };
          //return std::vector<std::shared_ptr<ETP::Triangle<P,D>>>{ triangles[ 0 ] };
          bodySearch = searchInsideLvl2Corridor( _startPoint, _goalPoint, startRoot, goalRoot, _radius, _scale );
          return bodySearch;
          if( bodySearch.empty() == false ){
            proceed = true;
          }
        }

      //  if( startRoot->haveSameLvl2CorridorEndpoints( goalRoot )){ return std::vector<std::shared_ptr<ETP::Triangle<P,D>>>{ _startTri, _goalTri }; }

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

    //tmp
    //return searchClassic( _startPoint, _goalPoint, _startTri, _goalTri, _radius, _scale );
    ETP::searchResult<P,D> classicSearch = searchBetweenLvl3EndPoints( _startPoint, _goalPoint, std::make_shared<ETP::searchEndPoint<P,D>>( 0.0, _startTri, nullptr ),  std::vector<std::shared_ptr<ETP::searchEndPoint<P,D>>>{ std::make_shared<ETP::searchEndPoint<P,D>>( 0.0, _goalTri, nullptr ) }, std::numeric_limits<D>::infinity(), _radius, _scale );
    return classicSearch.getFunnel();


  }
  ETP::searchResult<P,D> searchBetweenLvl3EndPoints( ETP::Point<P,D> _startPoint, ETP::Point<P,D> _goalPoint, std::shared_ptr<ETP::searchEndPoint<P,D>> _startEndPoint, std::vector<std::shared_ptr<ETP::searchEndPoint<P,D>>> _goalEndPoints, D _maxGValue, D _radius, D _scale ){
    std::vector<std::shared_ptr<ETP::Triangle<P,D>>> openList;
    //ETP::UsedList<P,D> usedList = ETP::UsedList<P,D>();
    std::shared_ptr<ETP::UsedList<P,D>> usedList = std::make_shared<ETP::UsedList<P,D>>();
    std::shared_ptr<ETP::Triangle<P,D>> startTri = _startEndPoint->triangle;
    setTriangleSearchNode( startTri, ETP::SearchNode<P,D>( 0.0, 0.0, nullptr ), usedList );
    openList.push_back( startTri );
    D memLomerGValue;
    bool endPointFound = false;
    while( openList.empty() == false ){
      std::shared_ptr<ETP::Triangle<P,D>> tCurrent = openList.front();
      openList.erase( openList.begin() );
      std::shared_ptr<ETP::SearchNode<P,D>> searchNodeCurrent = tCurrent->getSearchNode();
      D tCGVal = searchNodeCurrent->gValue;
      for( int i = 0; i < 3; i++ ){
        std::shared_ptr<ETP::Triangle<P,D>> tNeighbour = tCurrent->getConnectedNode( i );
        std::shared_ptr<ETP::Edge<P,D>> tCurrentOuterEdge = tCurrent->getEdgeVal( i );
        if( tNeighbour == nullptr || /* tNeighbour == searchNodeCurrent->comeFrom || */ searchNodeCurrent->comeFromIdx == i ){ continue; }
        //int neighbourToCurrentIdx = tNeighbour->getConnectedNodeIndex( tCurrent, tCurrentOuterEdge );
        int neighbourToCurrentIdx = tCurrent->getIndexfromConnectedNode( i );
        std::shared_ptr<ETP::Edge<P,D>> neighbEntryEdge = tNeighbour->getEdgeVal( neighbourToCurrentIdx );

        D startToClosestEdgeDist = neighbEntryEdge->minimumDistance( _startPoint, _scale ) / _scale;
        D tCurrentToNeighbLowerBound = tCGVal + tCurrent->getLowerBound( i ) * _radius;

        D hVal = neighbEntryEdge->minimumDistance( _goalPoint, _scale ) / _scale;
        D hValDiff = tCGVal + searchNodeCurrent->hValue - hVal;

        if( tCurrent == startTri ){
          tCurrentToNeighbLowerBound += _startEndPoint->entryEdge != nullptr ? tCurrent->getAngle( _startEndPoint->entryEdge, tCurrentOuterEdge ) * _radius :  getPointToEdgeLowerBound(  _startPoint, tCurrentOuterEdge, _radius, _scale );
        }else if( tCurrent->getLevel() == 3 ){
          tCurrentToNeighbLowerBound += tCurrent->getAngle( tCurrent->getEdgeVal( searchNodeCurrent->comeFromIdx ), tCurrentOuterEdge ) * _radius;
        }

        D gVal = std::max( hValDiff, std::max( startToClosestEdgeDist, tCurrentToNeighbLowerBound ) );

        if( gVal >= _maxGValue || ( endPointFound == true && gVal >= memLomerGValue ) ){ continue; }
        bool neighbIsEndPoint = false;

        for( std::shared_ptr<ETP::searchEndPoint<P,D>> gEndPt : _goalEndPoints ){
          if( gEndPt->triangle == tNeighbour ){
            neighbIsEndPoint = true;
            //int neighbToCurrentIdx = tNeighbour->getConnectedNodeIndex( tCurrent,  );
            if( gEndPt->entryEdge != nullptr )
            gVal += gEndPt->entryEdge != nullptr ? tNeighbour->getAngle( neighbEntryEdge, gEndPt->entryEdge ) * _radius : getPointToEdgeLowerBound( _goalPoint, neighbEntryEdge, _radius, _scale );
            //gVal += gEndPt->gValue;

            if( endPointFound == false || gEndPt->gValue > gVal ){
              gEndPt->gValue = gVal;
              //remove triangles from openList which have a gValue higher than the endPoint's gValue
              openList.erase( std::remove_if( openList.begin(), openList.end(),
                [ gVal ]( std::shared_ptr<ETP::Triangle<P,D>> _t ){
                    if ( _t->getSearchNode()->gValue >= gVal ) {
                        return true;
                    }
                    return false;
                }
              ), openList.end() );
            }

            endPointFound = true;
            memLomerGValue = gVal;

            break;
          }// end of if( gEndPt->triangle == tNeighbour )
        }// end of for( std::shared_ptr<ETP::searchEndPoint<P,D>> gEndPt : _goalEndPoints )
        std::shared_ptr<ETP::SearchNode<P,D>> tNeighbSearchNode = tNeighbour->getSearchNode();
        if( tNeighbSearchNode == nullptr ){
          setTriangleSearchNode( tNeighbour, ETP::SearchNode<P,D>( hVal, gVal, tCurrent, neighbourToCurrentIdx ), usedList );
          if( neighbIsEndPoint == false ){ insertInOpenList( hVal, tNeighbour, openList ); }
        }else if( tNeighbSearchNode->gValue > gVal ){
          setTriangleSearchNode( tNeighbour, ETP::SearchNode<P,D>( hVal, gVal, tCurrent, neighbourToCurrentIdx ), usedList );
          if( neighbIsEndPoint == false ){
            removeFromOpenList( tNeighbour, openList );
            insertInOpenList( hVal, tNeighbour, openList );
          }
        }


      }//end of for( int i = 0; i < 3; i++ )


    }//end of while( openList.empty() == false )

    std::shared_ptr<ETP::searchEndPoint<P,D>> closestEndPoint = nullptr;
    for( std::shared_ptr<ETP::searchEndPoint<P,D>> testedClosestEndPoint : _goalEndPoints ){
      std::shared_ptr<ETP::SearchNode<P,D>> testedClosestSearchNode = testedClosestEndPoint->triangle->getSearchNode();
      //if( testedClosestSearchNode != nullptr && ( closestEndPoint == nullptr || closestEndPoint->gValue + closestEndPoint->triangle->getSearchNode()->gValue > testedClosestEndPoint->gValue + testedClosestEndPoint->gValue ) ){
      if( testedClosestSearchNode != nullptr && ( closestEndPoint == nullptr || _startEndPoint->gValueToSearchStart + testedClosestEndPoint->gValue + testedClosestEndPoint->gValueToSearchStart < _startEndPoint->gValueToSearchStart + closestEndPoint->gValue + closestEndPoint->gValueToSearchStart ) ){
        closestEndPoint = testedClosestEndPoint;
      }
    }
    return ETP::searchResult<P,D>( closestEndPoint, usedList );
  }// end of searchBetweenLvl3EndPoints

  std::vector<std::shared_ptr<ETP::Triangle<P,D>>> searchClassic( ETP::Point<P,D> _startPoint, ETP::Point<P,D> _goalPoint, std::shared_ptr<ETP::Triangle<P,D>> _startTri, std::shared_ptr<ETP::Triangle<P,D>> _goalTri, D _radius, D _scale ){
    //ETP::Point<D,D> goalCentroid = ETP::Point<D,D>( (D) _goalPoint.x / _scale, (D) _goalPoint.y / _scale );
    std::vector<std::shared_ptr<ETP::Triangle<P,D>>> openList;
    //ETP::UsedList<P,D> usedList = ETP::UsedList<P,D>();
    std::shared_ptr<ETP::UsedList<P,D>> usedList = std::make_shared<ETP::UsedList<P,D>>();
    std::vector<std::shared_ptr<ETP::searchEndPoint<P,D>>> goalEndPoints = getTriangleEndPoints( _goalTri, _goalPoint, _radius, _scale );
    setTriangleSearchNode( _startTri, ETP::SearchNode<P,D>( 0.0, 0.0, nullptr ), usedList );
    openList.push_back( _startTri );
    size_t gel = goalEndPoints.size();

    D memLomerGValue;
    bool endPointFound = false;

    while( openList.empty() == false ){
      std::shared_ptr<ETP::Triangle<P,D>> tCurrent = openList.front();
      openList.erase( openList.begin() );
      std::shared_ptr<ETP::SearchNode<P,D>> searchNodeCurrent = tCurrent->getSearchNode();
      D tCGVal = searchNodeCurrent->gValue;
      for( int i = 0; i < 3; i++ ){
        std::shared_ptr<ETP::Triangle<P,D>> tNeighbour = tCurrent->getConnectedNode( i );
        std::shared_ptr<ETP::Edge<P,D>> tCurrentOuterEdge = tCurrent->getEdgeVal( i );

        if( tNeighbour == nullptr || tNeighbour == searchNodeCurrent->comeFrom || searchNodeCurrent->comeFromIdx == i ){ continue; }
        int neighbourToCurrentIdx = tNeighbour->getConnectedNodeIndex( tCurrent, tCurrentOuterEdge );
        std::shared_ptr<ETP::Edge<P,D>> neighbEntryEdge = tNeighbour->getEdgeVal( neighbourToCurrentIdx );

        D startToClosestEdgeDist = neighbEntryEdge->minimumDistance( _startPoint, _scale ) / _scale;
        D tCurrentToNeighbLowerBound = tCGVal + tCurrent->getLowerBound( i ) * _radius;

        D hVal = neighbEntryEdge->minimumDistance( _goalPoint, _scale ) / _scale;
        D hValDiff = tCGVal + searchNodeCurrent->hValue - hVal;


        //D gVal = tCGVal + tCurrent->getLowerBound( i ) * _radius;

        if( tCurrent == _startTri ){
          tCurrentToNeighbLowerBound += getPointToEdgeLowerBound(  _startPoint, tCurrent->getEdgeVal( i ), _radius, _scale );
        }else if( tCurrent->getLevel() == 3 ){
          if( searchNodeCurrent->comeFrom != nullptr ){
            if( searchNodeCurrent->comeFrom->getLevel() == 3 ){
              //tCurrentToNeighbLowerBound += tCurrent->getAngle( tCurrent->getEdgeVal( tCurrent->getConnectedNodeIndex( searchNodeCurrent->comeFrom ) ), tCurrent->getEdgeVal( i ) ) * _radius;
              tCurrentToNeighbLowerBound += tCurrent->getAngle( tCurrent->getEdgeVal(  searchNodeCurrent->comeFromIdx ), tCurrent->getEdgeVal( i ) ) * _radius;
            }else if( searchNodeCurrent->comeFrom->getLevel() == 2 ){
              tCurrentToNeighbLowerBound += tCurrent->getAngle( tCurrent->getEdgeVal( searchNodeCurrent->comeFrom->getReversedConnectedNodeIndex( tCurrent ) ), tCurrent->getEdgeVal( i ) ) * _radius;
            }

          }
        }
        D gVal = std::max( hValDiff, std::max( startToClosestEdgeDist, tCurrentToNeighbLowerBound ) );

        if( endPointFound == true && gVal >= memLomerGValue ){ continue; }
        bool neighbIsEndPoint = false;

        for( std::shared_ptr<ETP::searchEndPoint<P,D>> gEndPt : goalEndPoints ){
          if( gEndPt->triangle == tNeighbour ){
            neighbIsEndPoint = true;
            //int neighbToCurrentIdx = tNeighbour->getConnectedNodeIndex( tCurrent,  );

            if( tNeighbour == _goalTri && tNeighbour->getLevel() == 3 ){
              gVal += getPointToEdgeLowerBound( _goalPoint, _goalTri->getEdgeVal( neighbourToCurrentIdx ), _radius, _scale );
            }else{
              gVal += tNeighbour->getAngle( tNeighbour->getEdgeVal( gEndPt->parentIndex ), tNeighbour->getEdgeVal( neighbourToCurrentIdx ) ) * _radius;
            }

            if( endPointFound == false ){
              //remove triangles from openList which have a gValue higher than the endPoint's gValue
              openList.erase( std::remove_if( openList.begin(), openList.end(),
                [ gVal ]( std::shared_ptr<ETP::Triangle<P,D>> _t ){
                    if ( _t->getSearchNode()->gValue >= gVal ) {
                        return true;
                    }
                    return false;
                }
              ), openList.end() );
            }

            endPointFound = true;
            memLomerGValue = gVal;
            break;
          }
        }
        //D hVal = tNeighbour->getCentroid().dist( goalCentroid );
        std::shared_ptr<ETP::SearchNode<P,D>> tNeighbSearchNode = tNeighbour->getSearchNode();
        if( tNeighbSearchNode == nullptr ){
          if( tCurrent == _startTri && tCurrent->getLevel() == 2 ){
            setTriangleSearchNode( tNeighbour, ETP::SearchNode<P,D>( hVal, gVal, tCurrent, tCurrent->getReversedConnectedNodeIndex( tNeighbour ) ), usedList );
          }else{
            setTriangleSearchNode( tNeighbour, ETP::SearchNode<P,D>( hVal, gVal, tCurrent ), usedList );
          }
          if( neighbIsEndPoint == false ){ insertInOpenList( hVal, tNeighbour, openList ); }
        }else if( tNeighbSearchNode->gValue > gVal ){
        //  tNeighbour->setSearchHValue( hVal );
        //  tNeighbour->setSearchGValue( gVal );
        //  tNeighbour->setSearchComeFrom( tCurrent );
        setTriangleSearchNode( tNeighbour, ETP::SearchNode<P,D>( hVal, gVal, tCurrent ), usedList );
          if( neighbIsEndPoint == false ){
            removeFromOpenList( tNeighbour, openList );
            insertInOpenList( hVal, tNeighbour, openList );
          }
        }
      }
    }

    std::shared_ptr<ETP::searchEndPoint<P,D>> closestEndPoint = nullptr;
    bool closestFound = false;
    //int closestSearchIndex;
    for( std::shared_ptr<ETP::searchEndPoint<P,D>> testedClosestEndPoint : goalEndPoints ){
      std::shared_ptr<ETP::SearchNode<P,D>> testedClosestSearchNode = testedClosestEndPoint->triangle->getSearchNode();
      if( testedClosestSearchNode != nullptr && ( closestEndPoint == nullptr || closestEndPoint->gValue + closestEndPoint->triangle->getSearchNode()->gValue > testedClosestEndPoint->gValue + testedClosestEndPoint->gValue ) ){
        closestEndPoint = testedClosestEndPoint;
        closestFound = true;
      }
    }
    if( closestFound == false ){
      //return std::vector<std::shared_ptr<ETP::Triangle<P,D>>>();
      return std::vector<std::shared_ptr<ETP::Triangle<P,D>>>{ triangles [1 ] };
    }
    bool endPointReached = false;
    std::vector<std::shared_ptr<ETP::Triangle<P,D>>> funnelEndPart = getFunnelToEndPoint( closestEndPoint );

//return std::vector<std::shared_ptr<ETP::Triangle<P,D>>>{ _startTri, closestEndPoint->triangle, closestEndPoint->triangle->getSearchNode()->comeFrom };

    std::shared_ptr<ETP::Triangle<P,D>> tCurrent = closestEndPoint->triangle;
    while( endPointReached == false ){
      std::shared_ptr<ETP::SearchNode<P,D>> tSN = tCurrent->getSearchNode();
      if( tSN == nullptr || tSN->comeFrom == nullptr  ){
        endPointReached = true;
        break;
      }
      std::shared_ptr<ETP::Triangle<P,D>> tComeFrom = tSN->comeFrom;
      std::shared_ptr<ETP::Triangle<P,D>> tNode = tComeFrom;
      bool partEnd = false;
      std::vector<std::shared_ptr<ETP::Triangle<P,D>>> tmpFunnel;
      while( partEnd == false ){
        if( tNode == tCurrent ){
          break;
        }
        tmpFunnel.push_back( tNode );
        if( tNode == tComeFrom ){
          int lowestDistIdx = -1;
          D memLowerBd;
          for( int i = 0; i < 3; i++ ){
            if( tNode->getConnectedNode( i ) == tCurrent ){
              if( lowestDistIdx == -1 || memLowerBd > tNode->getLowerBound( i ) ){
                lowestDistIdx = i;
                memLowerBd =  tNode->getLowerBound( i );
              }
            }
          }
          //tNode = tNode->getAdjacentVal( tNode->getConnectedNodeIndex( tCurrent, tCurrent->getEdgeVal( tSN->comeFromIdx ) ) );
          tNode = tNode->getAdjacentVal( lowestDistIdx );
        }else{
          tNode = tNode->getAdjacentVal( tNode->getConnectedNodeIndex( tCurrent ) );
        }
      }
      funnelEndPart.insert( funnelEndPart.begin(), tmpFunnel.begin(), tmpFunnel.end() );
      tCurrent = tComeFrom;
    }
    return funnelEndPart;
  }

  std::vector<std::shared_ptr<ETP::Triangle<P,D>>> searchClassicBACKUP( ETP::Point<P,D> _startPoint, ETP::Point<P,D> _goalPoint, std::shared_ptr<ETP::Triangle<P,D>> _startTri, std::shared_ptr<ETP::Triangle<P,D>> _goalTri, D _radius, D _scale ){
    ETP::Point<D,D> goalCentroid = ETP::Point<D,D>( (D) _goalPoint.x / _scale, (D) _goalPoint.y / _scale );
    std::vector<std::shared_ptr<ETP::Triangle<P,D>>> openList;
    //ETP::UsedList<P,D> usedList = ETP::UsedList<P,D>();
    std::shared_ptr<ETP::UsedList<P,D>> usedList = std::make_shared<ETP::UsedList<P,D>>();
    std::vector<std::shared_ptr<ETP::searchEndPoint<P,D>>> goalEndPoints = getTriangleEndPoints( _goalTri, _goalPoint, _radius, _scale );
    setTriangleSearchNode( _startTri, ETP::SearchNode<P,D>( 0.0, 0.0, nullptr ), usedList );
    openList.push_back( _startTri );
    size_t gel = goalEndPoints.size();

    D memLomerGValue;
    bool endPointFound = false;

    while( openList.empty() == false ){
      std::shared_ptr<ETP::Triangle<P,D>> tCurrent = openList.front();
      openList.erase( openList.begin() );
      std::shared_ptr<ETP::SearchNode<P,D>> searchNodeCurrent = tCurrent->getSearchNode();
      D tCGVal = searchNodeCurrent->gValue;
      for( int i = 0; i < 3; i++ ){
        std::shared_ptr<ETP::Triangle<P,D>> tNeighbour = tCurrent->getConnectedNode( i );

        if( tNeighbour == nullptr || tNeighbour == searchNodeCurrent->comeFrom || searchNodeCurrent->comeFromIdx == i ){ continue; }
        D gVal = tCGVal + tCurrent->getLowerBound( i ) * _radius;
        /*
        if( tCurrent->getLevel() == 3 && searchNodeCurrent->comeFrom != nullptr ){
          gVal += tCurrent->getAngle( tCurrent->getEdgeVal( tCurrent->getConnectedNodeIndex( searchNodeCurrent->comeFrom ) ), tCurrent->getEdgeVal( i ) );
        }
        */
        if( tCurrent == _startTri ){
          gVal += getPointToEdgeLowerBound(  _startPoint, tCurrent->getEdgeVal( i ), _radius, _scale );
        }else if( tCurrent->getLevel() == 3 ){
          if( searchNodeCurrent->comeFrom != nullptr ){
            if( searchNodeCurrent->comeFrom->getLevel() == 3 ){
              gVal += tCurrent->getAngle( tCurrent->getEdgeVal( tCurrent->getConnectedNodeIndex( searchNodeCurrent->comeFrom ) ), tCurrent->getEdgeVal( i ) ) * _radius;
            }else if( searchNodeCurrent->comeFrom->getLevel() == 2 ){
              gVal += tCurrent->getAngle( tCurrent->getEdgeVal( searchNodeCurrent->comeFrom->getReversedConnectedNodeIndex( tCurrent ) ), tCurrent->getEdgeVal( i ) ) * _radius;
            }

          }
        }
        if( endPointFound == true && gVal >= memLomerGValue ){ continue; }
        bool neighbIsEndPoint = false;
        for( auto gEndPt : goalEndPoints ){
          if( gEndPt->triangle == tNeighbour ){
            neighbIsEndPoint = true;
            int neighbToCurrentIdx = tNeighbour->getConnectedNodeIndex( tCurrent );
            if( tNeighbour == _goalTri && tNeighbour->getLevel() == 3 ){
              gVal += getPointToEdgeLowerBound( _goalPoint, _goalTri->getEdgeVal( neighbToCurrentIdx ), _radius, _scale );
            }else{
              gVal += tNeighbour->getAngle( tNeighbour->getEdgeVal( gEndPt->parentIndex ), tNeighbour->getEdgeVal( neighbToCurrentIdx ) ) * _radius;
            }
            if( endPointFound == false ){
              //remove triangles from openList which have a gValue higher than the endPoint's gValue
              openList.erase( std::remove_if( openList.begin(), openList.end(),
                [ gVal ]( std::shared_ptr<ETP::Triangle<P,D>> _t ){
                    if ( _t->getSearchNode()->gValue >= gVal ) {
                        return true;
                    }
                    return false;
                }
              ), openList.end() );
            }

            endPointFound = true;
            memLomerGValue = gVal;
            break;
          }
        }
        D hVal = tNeighbour->getCentroid().dist( goalCentroid );
        std::shared_ptr<ETP::SearchNode<P,D>> tNeighbSearchNode = tNeighbour->getSearchNode();
        if( tNeighbSearchNode == nullptr ){
          if( tCurrent == _startTri && tCurrent->getLevel() == 2 ){
            setTriangleSearchNode( tNeighbour, ETP::SearchNode<P,D>( hVal, gVal, tCurrent, tCurrent->getReversedConnectedNodeIndex( tNeighbour ) ), usedList );
          }else{
            setTriangleSearchNode( tNeighbour, ETP::SearchNode<P,D>( hVal, gVal, tCurrent ), usedList );

          }
          if( neighbIsEndPoint == false ){ insertInOpenList( hVal, tNeighbour, openList ); }
        }else if( tNeighbSearchNode->gValue > gVal ){
        //  tNeighbour->setSearchHValue( hVal );
        //  tNeighbour->setSearchGValue( gVal );
        //  tNeighbour->setSearchComeFrom( tCurrent );
        setTriangleSearchNode( tNeighbour, ETP::SearchNode<P,D>( hVal, gVal, tCurrent ), usedList );
          if( neighbIsEndPoint == false ){
            removeFromOpenList( tNeighbour, openList );
            insertInOpenList( hVal, tNeighbour, openList );
          }
        }
      }
    }

    std::shared_ptr<ETP::searchEndPoint<P,D>> closestEndPoint = nullptr;
    bool closestFound = false;
    int closestSearchIndex;
    for( int i = 0; i < gel; i++ ){
      std::shared_ptr<ETP::searchEndPoint<P,D>> testedClosestEndPoint = goalEndPoints[ i ];
      std::shared_ptr<ETP::SearchNode<P,D>> testedClosestSearchNode = testedClosestEndPoint->triangle->getSearchNode();
      if( testedClosestSearchNode != nullptr && ( closestEndPoint == nullptr || closestEndPoint->gValue + closestEndPoint->triangle->getSearchNode()->gValue > testedClosestEndPoint->gValue + testedClosestEndPoint->gValue ) ){
        closestEndPoint = testedClosestEndPoint;
        closestFound = true;
      }
    }
    if( closestFound == false ){
      return std::vector<std::shared_ptr<ETP::Triangle<P,D>>>();
    }
    bool endPointReached = false;
    std::vector<std::shared_ptr<ETP::Triangle<P,D>>> funnelEndPart = getFunnelToEndPoint( closestEndPoint );



    std::shared_ptr<ETP::Triangle<P,D>> tCurrent = closestEndPoint->triangle;
    while( endPointReached == false ){
      std::shared_ptr<ETP::SearchNode<P,D>> tSN = tCurrent->getSearchNode();
      if( tSN == nullptr || tSN->comeFrom == nullptr  ){
        endPointReached = true;
        break;
      }
      std::shared_ptr<ETP::Triangle<P,D>> tComeFrom = tSN->comeFrom;
      std::shared_ptr<ETP::Triangle<P,D>> tNode = tComeFrom;
      bool partEnd = false;
      std::vector<std::shared_ptr<ETP::Triangle<P,D>>> tmpFunnel;
      while( partEnd == false ){
        if( tNode == tCurrent ){
          break;
        }
        tmpFunnel.push_back( tNode );
        tNode = tNode->getAdjacentVal( tNode->getConnectedNodeIndex( tCurrent ) );

      }
      funnelEndPart.insert( funnelEndPart.begin(), tmpFunnel.begin(), tmpFunnel.end() );
      tCurrent = tComeFrom;
    }
    return funnelEndPart;
  }

  std::vector<std::shared_ptr<ETP::Triangle<P,D>>> searchInsideLvl2RingOrLoop( ETP::Point<P,D> _startPoint, ETP::Point<P,D> _goalPoint, std::shared_ptr<ETP::Triangle<P,D>> _startTri, std::shared_ptr<ETP::Triangle<P,D>> _goalTri, D _radius, D _scale ){
    std::vector<std::shared_ptr<ETP::Triangle<P,D>>> openList;
    //ETP::UsedList<P,D> usedList = ETP::UsedList<P,D>();
    std::shared_ptr<ETP::UsedList<P,D>> usedList = std::make_shared<ETP::UsedList<P,D>>();
    setTriangleSearchNode( _startTri, ETP::SearchNode<P,D>( 0.0, 0.0 ), usedList );
    ETP::Point<D,D> startPoint = _startTri->isPointInside( _startPoint ) ? ETP::Point<D,D>( (D) _startPoint.x, (D) _startPoint.y ) : _startTri->getCentroid();
    ETP::Point<D,D> goalPoint = _goalTri->isPointInside( _goalPoint ) ? ETP::Point<D,D>( (D) _goalPoint.x, (D) _goalPoint.y ) : _goalTri->getCentroid();
    openList.push_back( _startTri );
    bool goalReached = false;
    D memGValue;
    while( openList.empty() == false ){
      std::shared_ptr<ETP::Triangle<P,D>> tCurrent = openList.front();
      openList.erase( openList.begin() );
      ETP::SearchNode<P,D> snCurrent = *tCurrent->getSearchNode().get();
      ETP::Point<D,D> gValPoint = tCurrent == _startTri ? startPoint : tCurrent->getCentroid();
      int currentLvl = tCurrent->getLevel();
      for( int i = 0; i < 3; i++ ){
        if( currentLvl == 3 && tCurrent->getConnectedNode( i ) != tCurrent ){ continue; }
        std::shared_ptr<ETP::Triangle<P,D>> tNeighbour = tCurrent->getAdjacentVal( i );
        if( tNeighbour == nullptr || tNeighbour->getLevel() < 2 || tNeighbour->getSearchNode() != nullptr ){ continue; }
        ETP::Point<D,D> tNeighbCentroid = tNeighbour->getCentroid();
        D distToNeighbour = tNeighbour != _goalTri ? gValPoint.dist( tNeighbCentroid ) : gValPoint.dist( goalPoint );
        D neighbGValue = snCurrent.gValue + distToNeighbour;
        if( goalReached == true && neighbGValue >= memGValue ){
          return getFunnel( _goalTri );
        }else if( tNeighbour == _goalTri ){
          if( goalReached == false || neighbGValue < memGValue ){
            memGValue = neighbGValue;
            //tNeighbour->setSearchNode( ETP::SearchNode<P,D>( 0.0, memGValue, tCurrent ), usedList );
            setTriangleSearchNode( tNeighbour, ETP::SearchNode<P,D>( 0.0, memGValue, tCurrent ), usedList );
            if( goalReached == true ){
              return getFunnel( _goalTri );
            }
            goalReached = true;
          }
        }else{
          D neighbHValue = tNeighbCentroid.dist( goalPoint );
          setTriangleSearchNode( tNeighbour, ETP::SearchNode<P,D>( neighbHValue, neighbGValue, tCurrent ), usedList );
          insertInOpenList( neighbHValue, tNeighbour, openList );
        }
      }
    }
    return getFunnel( _goalTri );
  }

  std::vector<std::shared_ptr<ETP::Triangle<P,D>>> searchInsideLvl2RingOrLoopBACKUP( ETP::Point<P,D> _startPoint, ETP::Point<P,D> _goalPoint, std::shared_ptr<ETP::Triangle<P,D>> _startTri, std::shared_ptr<ETP::Triangle<P,D>> _goalTri, D _radius, D _scale ){

    std::vector<std::shared_ptr<ETP::Triangle<P,D>>> openList;
    //std::shared_ptr<ETP::UsedList<P,D>> usedList = std::make_shared<ETP::UsedList<P,D>>();
    //ETP::UsedList<P,D> usedList = ETP::UsedList<P,D>();
    std::shared_ptr<ETP::UsedList<P,D>> usedList = std::make_shared<ETP::UsedList<P,D>>();
    ETP::Point<D,D> goalCentroid = _goalTri->getCentroid();
    setTriangleSearchNode( _startTri, ETP::SearchNode<P,D>( 0.0, 0.0 ), usedList );
    openList.push_back( _startTri );
    bool goalReached = false;
    D memGValue;
    while( openList.empty() == false ){
      std::shared_ptr<ETP::Triangle<P,D>> tCurrent = openList.front();
      openList.erase( openList.begin() );
      ETP::SearchNode<P,D> snCurrent = *tCurrent->getSearchNode().get();
      ETP::Point<D,D> currentCentroid = tCurrent->getCentroid();
      int currentLvl = tCurrent->getLevel();
      for( int i = 0; i < 3; i++ ){
        if( currentLvl == 3 && tCurrent->getConnectedNode( i ) != tCurrent ){ continue; }
        std::shared_ptr<ETP::Triangle<P,D>> tNeighbour = tCurrent->getAdjacentVal( i );
        if( tNeighbour == nullptr || tNeighbour->getLevel() < 2 || tNeighbour->getSearchNode() != nullptr ){ continue; }
        //ETP::Point<D,D> tNeighbCentroid = tNeighbour->getCentroid();
        ETP::Point<D,D> tNeighbCentroid = tNeighbour->getCentroid();
        D neighbGValue = snCurrent.gValue;
        std::shared_ptr<ETP::Edge<P,D>> edgeOut = tCurrent->getEdgeVal( i );
        std::shared_ptr<ETP::Edge<P,D>> edgeIn = tCurrent->getEdgeVal( tCurrent->getAdjacentIndex( snCurrent.comeFrom ) );
        if( tCurrent == _startTri ){
          neighbGValue += getPointToEdgeLowerBound(  _startPoint, _startTri->getEdgeVal( i ), _radius, _scale );
        }else{

          neighbGValue += tCurrent->getAngle( edgeIn, edgeOut ) * _radius;
        }
        if( tNeighbour == _goalTri ){
          neighbGValue += getPointToEdgeLowerBound(  _goalPoint, edgeOut, _radius, _scale );
        }
        if( goalReached == true && neighbGValue >= memGValue ){
          return getFunnel( _goalTri );
        }else if( tNeighbour == _goalTri ){
          if( goalReached == false || neighbGValue < memGValue ){
            memGValue = neighbGValue;
            //tNeighbour->setSearchNode( ETP::SearchNode<P,D>( 0.0, memGValue, tCurrent ), usedList );
            setTriangleSearchNode( tNeighbour, ETP::SearchNode<P,D>( 0.0, memGValue, tCurrent ), usedList );
            if( goalReached == true ){
              return getFunnel( _goalTri );
            }
            goalReached = true;
          }
        }else{
          D neighbHValue = tNeighbCentroid.dist( goalCentroid );
          setTriangleSearchNode( tNeighbour, ETP::SearchNode<P,D>( neighbHValue, neighbGValue, tCurrent ), usedList );
          insertInOpenList( neighbHValue, tNeighbour, openList );
        }
      }
    }
  }

  std::vector<std::shared_ptr<ETP::Triangle<P,D>>> searchInsideLvl2CorridorTEST( ETP::Point<P,D> _startPoint, ETP::Point<P,D> _goalPoint, std::shared_ptr<ETP::Triangle<P,D>> _startTri, std::shared_ptr<ETP::Triangle<P,D>> _goalTri, D _radius, D _scale ){

  }
  std::vector<std::shared_ptr<ETP::Triangle<P,D>>> searchInsideLvl2Corridor( ETP::Point<P,D> _startPoint, ETP::Point<P,D> _goalPoint, std::shared_ptr<ETP::Triangle<P,D>> _startTri, std::shared_ptr<ETP::Triangle<P,D>> _goalTri, D _radius, D _scale ){
    std::shared_ptr<ETP::Triangle<P,D>> endPoint = nullptr;
    std::shared_ptr<ETP::Triangle<P,D>> otherEndPoint;
    D sDist;

    for( int i = 0; i < 3; i++ ){
      std::shared_ptr<ETP::Triangle<P,D>> cn = _startTri->getConnectedNode( i );
      if( cn != nullptr ){
        if( endPoint == nullptr ){
          endPoint = cn;
          sDist = _startTri->getLowerBound( i ) * _radius /*+ ( startPtInStartTri ?  getPointToEdgeLowerBound( _startPoint, _startTri->getEdgeVal( i ), _radius, _scale ) : 0.0 ) */;
        }else{
          otherEndPoint = cn;
          break;
        }
      }
    }
    int gEndPointIdx = _goalTri->getConnectedNodeIndex( endPoint );
    D gDist = _goalTri->getLowerBound( gEndPointIdx ) * _radius /*+ ( goalPtInGoalTri ? getPointToEdgeLowerBound( _goalPoint, _goalTri->getEdgeVal( gEndPointIdx ), _radius, _scale ) : 0.0 )*/;
    std::shared_ptr<ETP::Triangle<P,D>> searchStart;
    std::shared_ptr<ETP::Triangle<P,D>> searchGoal;
    bool reverse;
    ETP::Point<P,D> searchStartPoint;
    ETP::Point<P,D> searchGoalPoint;
    if( sDist >= gDist ){
      searchStart = _startTri;
      searchGoal = _goalTri;
      reverse = true;
      searchStartPoint = _startPoint;
      searchGoalPoint = _goalPoint;
    }else{
      searchStart = _goalTri;
      searchGoal = _startTri;
      reverse = false;
      D tmpDist = gDist;
      gDist = sDist;
      sDist = tmpDist;
      searchStartPoint = _goalPoint;
      searchGoalPoint = _startPoint;
    }
  //  return std::vector<std::shared_ptr<ETP::Triangle<P,D>>>{ searchStart, searchGoal, endPoint };

    bool goalReached = false;
    std::vector<std::shared_ptr<ETP::Triangle<P,D>>> openList;
    openList.push_back( searchStart );
    D startToGoalGValue;
    while( goalReached == false ){
      std::shared_ptr<ETP::Triangle<P,D>> tCurrent = openList.front();
      std::shared_ptr<ETP::Triangle<P,D>> tNeighbour = tCurrent->getAdjacentVal( tCurrent->getConnectedNodeIndex( endPoint ) );
      if( tNeighbour == searchGoal ){
        openList.insert( openList.begin(), tNeighbour );
        goalReached = true;
        startToGoalGValue = sDist - gDist;
        if( searchStart->isPointInside( searchStartPoint, _scale ) == true ){
          startToGoalGValue += getPointToEdgeLowerBound( searchStartPoint, searchStart->getEdgeVal( searchStart->getConnectedNodeIndex( endPoint ) ), _radius, _scale );
        }
        if( searchGoal->isPointInside( searchGoalPoint, _scale ) == true ){
          startToGoalGValue += getPointToEdgeLowerBound( searchGoalPoint, searchGoal->getEdgeVal( searchGoal->getConnectedNodeIndex( otherEndPoint ) ), _radius, _scale );
        }
        break;
      }else if( tNeighbour->getLevel() == 3 ){
        //return std::vector<std::shared_ptr<ETP::Triangle<P,D>>>{ triangles[ 6 ] };
        break;
      }else{
        openList.insert( openList.begin(), tNeighbour );
      }
    }
    if( goalReached == false ){
    //  return std::vector<std::shared_ptr<ETP::Triangle<P,D>>>{ triangles[ 5 ] };
      return std::vector<std::shared_ptr<ETP::Triangle<P,D>>>();
    }else{

      int otherEndPtToEndPtIdx = otherEndPoint->getConnectedNodeIndex( endPoint );
      //int otherEndPtToEndPtIdx = searchStart->getReversedConnectedNodeIndex( otherEndPoint );

      D baseGValue = //otherEndPoint->getLowerBound( otherEndPtToEndPtIdx ) * _radius - startToGoalGValue;
      getPointToEdgeLowerBound( searchStartPoint, searchStart->getEdgeVal( searchStart->getConnectedNodeIndex( otherEndPoint ) ), _radius, _scale ) + searchStart->getLowerBound( searchStart->getConnectedNodeIndex( otherEndPoint ) )
      + getPointToEdgeLowerBound( searchGoalPoint, searchGoal->getEdgeVal( searchGoal->getConnectedNodeIndex( endPoint ) ), _radius, _scale ) + searchGoal->getLowerBound( searchGoal->getConnectedNodeIndex( endPoint ) );

      if( baseGValue >= startToGoalGValue ){
        if( reverse == true ){
          std::reverse( openList.begin(), openList.end( ) );
        }
        //return std::vector<std::shared_ptr<ETP::Triangle<P,D>>>{ triangles[ 1 ] };
        return openList;
      }


      std::vector<std::shared_ptr<ETP::Triangle<P,D>>> secOpenList;
      //ETP::UsedList<P,D> usedList = ETP::UsedList<P,D>();
      std::shared_ptr<ETP::UsedList<P,D>> usedList = std::make_shared<ETP::UsedList<P,D>>();
      setTriangleSearchNode( otherEndPoint, ETP::SearchNode<P,D>( 0.0, baseGValue, endPoint ), usedList );
      secOpenList.push_back( otherEndPoint );
      ETP::Point<D,D> goalCentroid = endPoint->getCentroid();
      while( secOpenList.empty() == false ){
        std::shared_ptr<ETP::Triangle<P,D>> tCurrent = secOpenList.front();
        secOpenList.erase( secOpenList.begin() );
        ETP::SearchNode<P,D> tsn = *tCurrent->getSearchNode().get();
        std::shared_ptr<ETP::Edge<P,D>> edgeIn;
        if( tCurrent->getLevel() == 3 ){
          edgeIn = tCurrent->getEdgeVal( tCurrent->getConnectedNodeIndex( tsn.comeFrom ) );
        }
        for( int i = 0; i < 3; i++ ){
          std::shared_ptr<ETP::Triangle<P,D>> tConnectedNode = tCurrent->getConnectedNode( i );
          if( tConnectedNode == nullptr || ( tConnectedNode == tsn.comeFrom && tCurrent == otherEndPoint && i == otherEndPtToEndPtIdx ) || ( tConnectedNode == tsn.comeFrom && tCurrent != otherEndPoint ) ){ continue; }
          D tCNgValue = tCurrent->getLowerBound( i ) * _radius + tsn.gValue;
          if( tCurrent->getLevel() == 3 ){
            tCNgValue += tCurrent->getAngle( tCurrent->getEdgeVal( i ), edgeIn ) * _radius;
          }
          std::shared_ptr<ETP::SearchNode<P,D>> tCsN = tConnectedNode->getSearchNode();
          if( tCNgValue >= startToGoalGValue || ( tCsN != nullptr && tCsN->gValue <= tCNgValue ) ){ continue; }
          D tCNhValue = tConnectedNode->getCentroid().dist( goalCentroid );
          if( tConnectedNode == endPoint ){
            tCNgValue += endPoint->getAngle( endPoint->getEdgeVal( endPoint->getConnectedNodeIndex( tCurrent ) ), endPoint->getEdgeVal( endPoint->getConnectedNodeIndex( otherEndPoint ) ) ) * _radius;
            if( tCsN == nullptr ){
              setTriangleSearchNode( tConnectedNode, ETP::SearchNode<P,D>( tCNhValue, tCNgValue, tCurrent ), usedList );
            }else if( tCsN->gValue > tCNgValue ){
              tCsN->hValue = tCNhValue;
              tCsN->gValue = tCNgValue;
              tCsN->comeFrom = tCurrent;
            }
          }else{
            //tConnectedNode->setSearchNode( ETP::SearchNode<P,D>( tCNhValue, tCNgValue, tCurrent ), usedList );
            setTriangleSearchNode( tConnectedNode, ETP::SearchNode<P,D>( tCNhValue, tCNgValue, tCurrent ), usedList );
            insertInOpenList( tCNhValue, tConnectedNode, secOpenList );
          }
        }
      }

      if( endPoint->getSearchNode() != nullptr && endPoint->getSearchNode()->gValue < startToGoalGValue ){
        setTriangleSearchNode( searchGoal, ETP::SearchNode<P,D>( 0.0, 0.0, endPoint ), usedList );
        setTriangleSearchNode( otherEndPoint, ETP::SearchNode<P,D>( 0.0, 0.0, searchStart ), usedList );
        std::vector<std::shared_ptr<ETP::Triangle<P,D>>> path = getFunnelWithIntermediateNodes( searchGoal, searchStart );
        if( searchGoal == _startTri ){
          std::reverse( path.begin(), path.end( ) );
        }
        //return std::vector<std::shared_ptr<ETP::Triangle<P,D>>>{ triangles[2 ] };
        return path;
      }else{
        if( reverse == true ){
          std::reverse( openList.begin(), openList.end( ) );
        }
        //return std::vector<std::shared_ptr<ETP::Triangle<P,D>>>{ triangles[ 0 ] };
        return openList;
      }
    }
  }

  std::vector<std::shared_ptr<ETP::searchEndPoint<P,D>>> getTriangleEndPointsTEST( std::shared_ptr<ETP::Triangle<P,D>> _tri, ETP::Point<P,D> _goalPoint, D _radius, D _scale, std::shared_ptr<ETP::UsedList<P,D>> usedList ){
    std::vector<std::shared_ptr<ETP::searchEndPoint<P,D>>> ret;
    int lvl = _tri->getLevel();
    if( lvl == 3 ){
      ret.push_back( std::make_shared<ETP::searchEndPoint<P,D>>( 0.0, _tri ) );
    }else if( lvl == 2 ){
      setTriangleSearchNode( _tri, ETP::SearchNode<P,D>( 0.0, 0.0, nullptr ), usedList );
      for( int i = 0; i < 3; i++ ){
        std::shared_ptr<ETP::Triangle<P,D>> connectedNode = _tri->getConnectedNode( i );
        if( connectedNode == nullptr ) continue;
        bool isNew = true;
        D lowerBound = getPointToEdgeLowerBound( _goalPoint, _tri->getEdgeVal( i ), _radius, _scale );
        for ( std::shared_ptr<ETP::searchEndPoint<P,D>> ePt : ret ){
          if( connectedNode == ePt->triangle ){
            isNew = false;
            if( lowerBound < ePt->gValueToSearchStart ){
              int idxFromEndPtToStart = _tri->getIndexfromConnectedNode( i );
              ePt = std::make_shared<ETP::searchEndPoint<P,D>>( 0.0, connectedNode, nullptr, idxFromEndPtToStart, lowerBound );
              ePt->entryEdge = connectedNode->getEdgeVal( idxFromEndPtToStart );
            }
          }
        }// end of for ( std::shared_ptr<ETP::searchEndPoint<P,D>> ePt : ret )
        if( isNew == true ){
          int idxFromEndPtToStart = _tri->getIndexfromConnectedNode( i );
          std::shared_ptr<ETP::searchEndPoint<P,D>> endPt = std::make_shared<ETP::searchEndPoint<P,D>>( 0.0, connectedNode, nullptr, idxFromEndPtToStart, lowerBound );
          endPt->entryEdge = connectedNode->getEdgeVal( idxFromEndPtToStart );
          ret.push_back( endPt );
        }
      }// end of for( int i = 0; i < 3; i++ )
    }else if( lvl == 1 ){
      setTriangleSearchNode( _tri, ETP::SearchNode<P,D>( 0.0, 0.0, nullptr ), usedList );
      for( int i = 0; i < 3; i++ ){
        std::shared_ptr<ETP::Triangle<P,D>> connectedNode = _tri->getConnectedNode( i );
        if( connectedNode == nullptr ) continue;
        D startLowerBound = getPointToEdgeLowerBound( _goalPoint, _tri->getEdgeVal( i ), _radius, _scale ) + _tri->getLowerBound( i ) * _radius;
        int indexFromLvl2ToLvl1 = _tri->getIndexfromConnectedNode( i );
        setTriangleSearchNode( connectedNode, ETP::SearchNode<P,D>( startLowerBound, 0.0, _tri, indexFromLvl2ToLvl1 ), usedList );
        for( int j = 0; j < 3; j++ ){
          std::shared_ptr<ETP::Triangle<P,D>> lvl3Node = connectedNode->getConnectedNode( j );
          if( lvl3Node == nullptr ) continue;
          D finalLowerBound = startLowerBound + connectedNode->getAngle( connectedNode->getEdgeVal( indexFromLvl2ToLvl1 ), connectedNode->getEdgeVal( j ) ) * _radius + connectedNode->getLowerBound( j ) * _radius;
          int idxFromLvl3ToLvl2 = connectedNode->getIndexfromConnectedNode( j );
          bool isNew = true;
          for ( std::shared_ptr<ETP::searchEndPoint<P,D>> ePt : ret ){
            if( lvl3Node == ePt->triangle ){
              isNew = false;
              if( finalLowerBound < ePt->gValueToSearchStart ){
                ePt = std::make_shared<ETP::searchEndPoint<P,D>>( 0.0, lvl3Node, nullptr, idxFromLvl3ToLvl2, finalLowerBound );
                ePt->entryEdge = lvl3Node->getEdgeVal( idxFromLvl3ToLvl2 );
              }
            }
          }
          if( isNew == true ){
            std::shared_ptr<ETP::searchEndPoint<P,D>> endPt = std::make_shared<ETP::searchEndPoint<P,D>>( 0.0, lvl3Node, nullptr, idxFromLvl3ToLvl2, finalLowerBound );
            endPt->entryEdge = lvl3Node->getEdgeVal( idxFromLvl3ToLvl2 );
            ret.push_back( endPt );
          }
        }
      }// end of for( int i = 0; i < 3; i++ )
    }
    return ret;
  }

  std::vector<std::shared_ptr<ETP::searchEndPoint<P,D>>> getTriangleEndPoints( std::shared_ptr<ETP::Triangle<P,D>> _tri, ETP::Point<P,D> _goalPoint, D _radius, D _scale ){
  std::vector<std::shared_ptr<ETP::searchEndPoint<P,D>>> ret;
  //  ETP::Triangle<P,D> tri = *_tri.get();
  std::shared_ptr<ETP::searchEndPoint<P,D>> finalEndPt = std::make_shared<ETP::searchEndPoint<P,D>>( 0.0, _tri );
  int lvl = _tri->getLevel();
  if( lvl == 3 ){
    ret.push_back( finalEndPt );
  }else if( lvl == 2 ){
    for( int i = 0; i < 3; i++ ){
      std::shared_ptr<ETP::Triangle<P,D>> connectedNode = _tri->getConnectedNode( i );
      if( connectedNode != nullptr ){
        bool isNew = true;
        D lowerBound = getPointToEdgeLowerBound( _goalPoint, _tri->getEdgeVal( i ), _radius, _scale );
        size_t retl = ret.size();
        for ( auto ePt : ret ){
          if( connectedNode == ePt->triangle ){
            isNew = false;
            if( lowerBound < ePt->gValue ){
              ePt = std::make_shared<ETP::searchEndPoint<P,D>>( lowerBound, connectedNode, finalEndPt, _tri->getReversedConnectedNodeIndex( connectedNode ) );
            }
          }
        }
        if( isNew == true ){
          ret.push_back( std::make_shared<ETP::searchEndPoint<P,D>>( lowerBound, connectedNode, finalEndPt, _tri->getReversedConnectedNodeIndex( connectedNode ) ) );
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
                  ret[ k ] = std::make_shared<ETP::searchEndPoint<P,D>>( lvl2EndPt->gValue + lowerBound, lvl3Node, lvl2EndPt );
                }
                break;
              }
            }
            if( isNew == true ){
              ret.push_back( std::make_shared<ETP::searchEndPoint<P,D>>( lvl2EndPt->gValue + lowerBound, lvl3Node, lvl2EndPt, connectedNode->getReversedConnectedNodeIndex( lvl3Node ) ) );
            }
          }
        }
        break;
      }
    }
  }
  return ret;
  }







  std::vector<std::shared_ptr<ETP::Triangle<P,D>>> abstractTriangulationSearch( ETP::Point<P,D> _startPoint, ETP::Point<P,D> _goalPoint, D _radius ){
    std::shared_ptr<ETP::Triangle<P,D>> _startTri = getTriangleWithPoint( _startPoint );
    std::shared_ptr<ETP::Triangle<P,D>> _goalTri = getTriangleWithPoint( _goalPoint );
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
      if(    ( startRoot == nullptr ) //start and goal belong to the same unrooted tree
          || ( startRoot == _goalTri->getLvl1RootNode() ) // start and goal belong to the same rooted tree
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
        bool proceed = false;
        std::vector<std::shared_ptr<ETP::Triangle<P,D>>> bodySearch;
        if( startRoot->isPartoflvl2Ring() || startRoot->isOnSameLvl2Loop( goalRoot ) ){// inside a loop or ring
          proceed = true;
          bodySearch =  searchInsideLvl2RingOrLoop( _startPoint, _goalPoint, _startTri, _goalTri, _radius );
        }else if( startRoot->haveSameLvl2CorridorEndpoints( goalRoot )){// inside a corridor
          //return std::vector<std::shared_ptr<ETP::Triangle<P,D>>>{ _startTri, _goalTri };
          bodySearch = searchInsideLvl2Corridor( _startPoint, _goalPoint, startRoot, goalRoot, _radius );
          if( bodySearch.empty() == false ){
            proceed = true;
          }
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

    //tmp
    return searchClassic( _startPoint, _goalPoint, _startTri, _goalTri, _radius );

  }
  std::vector<std::shared_ptr<ETP::Triangle<P,D>>> searchClassic( ETP::Point<P,D> _startPoint, ETP::Point<P,D> _goalPoint, std::shared_ptr<ETP::Triangle<P,D>> _startTri, std::shared_ptr<ETP::Triangle<P,D>> _goalTri, D _radius ){
    ETP::Point<D,D> goalCentroid = ETP::Point<D,D>( (D) _goalPoint.x, (D) _goalPoint.y );
    std::vector<std::shared_ptr<ETP::Triangle<P,D>>> openList;
    //ETP::UsedList<P,D> usedList = ETP::UsedList<P,D>();
    std::shared_ptr<ETP::UsedList<P,D>> usedList = std::make_shared<ETP::UsedList<P,D>>();
    std::vector<std::shared_ptr<ETP::searchEndPoint<P,D>>> goalEndPoints = getTriangleEndPoints( _goalTri, _goalPoint, _radius );
    setTriangleSearchNode( _startTri, ETP::SearchNode<P,D>( 0.0, 0.0, nullptr ), usedList );
    openList.push_back( _startTri );
    size_t gel = goalEndPoints.size();

    D memLomerGValue;
    bool endPointFound = false;

    while( openList.empty() == false ){
      std::shared_ptr<ETP::Triangle<P,D>> tCurrent = openList.front();
      openList.erase( openList.begin() );
      std::shared_ptr<ETP::SearchNode<P,D>> searchNodeCurrent = tCurrent->getSearchNode();
      D tCGVal = searchNodeCurrent->gValue;
      for( int i = 0; i < 3; i++ ){
        std::shared_ptr<ETP::Triangle<P,D>> tNeighbour = tCurrent->getConnectedNode( i );
        if( tNeighbour == nullptr || tNeighbour == searchNodeCurrent->comeFrom ){ continue; }
        D gVal = tCGVal + tCurrent->getLowerBound( i ) * _radius;
        /*
        if( tCurrent->getLevel() == 3 && searchNodeCurrent->comeFrom != nullptr ){
          gVal += tCurrent->getAngle( tCurrent->getEdgeVal( tCurrent->getConnectedNodeIndex( searchNodeCurrent->comeFrom ) ), tCurrent->getEdgeVal( i ) );
        }
        */
        if( tCurrent == _startTri ){
          gVal += getPointToEdgeLowerBound(  _startPoint, tCurrent->getEdgeVal( i ), _radius );
        }else if( tCurrent->getLevel() == 3 ){
          if( searchNodeCurrent->comeFrom != nullptr ){
            gVal += tCurrent->getAngle( tCurrent->getEdgeVal( tCurrent->getConnectedNodeIndex( searchNodeCurrent->comeFrom ) ), tCurrent->getEdgeVal( i ) ) * _radius;
          }
        }
        if( endPointFound == true && gVal >= memLomerGValue ){ continue; }
        bool neighbIsEndPoint = false;
        for( auto gEndPt : goalEndPoints ){
          if( gEndPt->triangle == tNeighbour ){
            neighbIsEndPoint = true;
            int neighbToCurrentIdx = tNeighbour->getConnectedNodeIndex( tCurrent );
            if( tNeighbour == _goalTri && tNeighbour->getLevel() == 3 ){
              gVal += getPointToEdgeLowerBound( _goalPoint, _goalTri->getEdgeVal( neighbToCurrentIdx ), _radius );
            }else{
              gVal += tNeighbour->getAngle( tNeighbour->getEdgeVal( gEndPt->parentIndex ), tNeighbour->getEdgeVal( neighbToCurrentIdx ) ) * _radius;
            }
            if( endPointFound == false ){
              //remove triangles from openList which have a gValue higher than the endPoint's gValue
              openList.erase( std::remove_if( openList.begin(), openList.end(),
                [ gVal ]( std::shared_ptr<ETP::Triangle<P,D>> _t ){
                    if ( _t->getSearchNode()->gValue >= gVal ) {
                        return true;
                    }
                    return false;
                }
              ), openList.end() );
            }

            endPointFound = true;
            memLomerGValue = gVal;
            break;
          }
        }
        D hVal = tNeighbour->getCentroid().dist( goalCentroid );
        std::shared_ptr<ETP::SearchNode<P,D>> tNeighbSearchNode = tNeighbour->getSearchNode();
        if( tNeighbSearchNode == nullptr ){
          setTriangleSearchNode( tNeighbour, ETP::SearchNode<P,D>( hVal, gVal, tCurrent ), usedList );
          if( neighbIsEndPoint == false ){ insertInOpenList( hVal, tNeighbour, openList ); }
        }else if( tNeighbSearchNode->gValue > gVal ){
          tNeighbour->setSearchHValue( hVal );
          tNeighbour->setSearchGValue( gVal );
          tNeighbour->setSearchComeFrom( tCurrent );
          if( neighbIsEndPoint == false ){
            removeFromOpenList( tNeighbour, openList );
            insertInOpenList( hVal, tNeighbour, openList );
          }
        }
      }
    }

    std::shared_ptr<ETP::searchEndPoint<P,D>> closestEndPoint = nullptr;
    bool closestFound = false;
    int closestSearchIndex;
    for( int i = 0; i < gel; i++ ){
      std::shared_ptr<ETP::searchEndPoint<P,D>> testedClosestEndPoint = goalEndPoints[ i ];
      std::shared_ptr<ETP::SearchNode<P,D>> testedClosestSearchNode = testedClosestEndPoint->triangle->getSearchNode();
      if( testedClosestSearchNode != nullptr && ( closestEndPoint == nullptr || closestEndPoint->gValue + closestEndPoint->triangle->getSearchNode()->gValue > testedClosestEndPoint->gValue + testedClosestEndPoint->gValue ) ){
        closestEndPoint = testedClosestEndPoint;
        closestFound = true;
      }
    }
    if( closestFound == false ){
      return std::vector<std::shared_ptr<ETP::Triangle<P,D>>>();
    }
    bool endPointReached = false;
    std::vector<std::shared_ptr<ETP::Triangle<P,D>>> funnelEndPart = getFunnelToEndPoint( closestEndPoint );
    std::shared_ptr<ETP::Triangle<P,D>> tCurrent = closestEndPoint->triangle;
    while( endPointReached == false ){
      std::shared_ptr<ETP::SearchNode<P,D>> tSN = tCurrent->getSearchNode();
      if( tSN == nullptr || tSN->comeFrom == nullptr){
        endPointReached = true;
        break;
      }
      std::shared_ptr<ETP::Triangle<P,D>> tComeFrom = tSN->comeFrom;
      std::shared_ptr<ETP::Triangle<P,D>> tNode = tComeFrom;
      bool partEnd = false;
      std::vector<std::shared_ptr<ETP::Triangle<P,D>>> tmpFunnel;
      while( partEnd == false ){
        if( tNode == tCurrent ){
          break;
        }
        tmpFunnel.push_back( tNode );
        tNode = tNode->getAdjacentVal( tNode->getConnectedNodeIndex( tCurrent ) );

      }
      funnelEndPart.insert( funnelEndPart.begin(), tmpFunnel.begin(), tmpFunnel.end() );
      tCurrent = tComeFrom;
    }
    return funnelEndPart;
  }
  std::vector<std::shared_ptr<ETP::Triangle<P,D>>> searchInsideLvl2RingOrLoop( ETP::Point<P,D> _startPoint, ETP::Point<P,D> _goalPoint, std::shared_ptr<ETP::Triangle<P,D>> _startTri, std::shared_ptr<ETP::Triangle<P,D>> _goalTri, D _radius ){
    std::vector<std::shared_ptr<ETP::Triangle<P,D>>> openList;
    //std::shared_ptr<ETP::UsedList<P,D>> usedList = std::make_shared<ETP::UsedList<P,D>>();
    //ETP::UsedList<P,D> usedList = ETP::UsedList<P,D>();
    std::shared_ptr<ETP::UsedList<P,D>> usedList = std::make_shared<ETP::UsedList<P,D>>();
    ETP::Point<D,D> goalCentroid = _goalTri->getCentroid();
    setTriangleSearchNode( _startTri, ETP::SearchNode<P,D>( 0.0, 0.0 ), usedList );
    openList.push_back( _startTri );
    bool goalReached = false;
    D memGValue;
    while( openList.empty() == false ){
      std::shared_ptr<ETP::Triangle<P,D>> tCurrent = openList.front();
      openList.erase( openList.begin() );
      ETP::SearchNode<P,D> snCurrent = *tCurrent->getSearchNode().get();
      ETP::Point<D,D> currentCentroid = tCurrent->getCentroid();
      int currentLvl = tCurrent->getLevel();
      for( int i = 0; i < 3; i++ ){
        if( currentLvl == 3 && tCurrent->getConnectedNode( i ) != tCurrent ){ continue; }
        std::shared_ptr<ETP::Triangle<P,D>> tNeighbour = tCurrent->getAdjacentVal( i );
        int neighbLvl = tNeighbour->getLevel();
        if( tNeighbour == nullptr || neighbLvl < 2 || tNeighbour->getSearchNode() != nullptr ){ continue; }
        //ETP::Point<D,D> tNeighbCentroid = tNeighbour->getCentroid();
        ETP::Point<D,D> tNeighbCentroid = tNeighbour->getCentroid();
        D neighbGValue = snCurrent.gValue;
        if( tCurrent == _startTri ){
          neighbGValue += getPointToEdgeLowerBound(  _startPoint, _startTri->getEdgeVal( i ), _radius );
        }else{
          std::shared_ptr<ETP::Edge<P,D>> edgeOut = tCurrent->getEdgeVal( i );
          std::shared_ptr<ETP::Edge<P,D>> edgeIn;
          for( int e = 0; e < 3; e++ ){
            std::shared_ptr<ETP::Edge<P,D>> te = tCurrent->getEdgeVal( e );
            if( te->isConstrained() == false && te != edgeOut ){
              edgeIn = te;
              break;
            }
          }
          neighbGValue += tCurrent->getAngle( edgeIn, edgeOut ) * _radius;
        }
        if( tNeighbour == _goalTri ){
          neighbGValue += getPointToEdgeLowerBound(  _goalPoint, _goalTri->getEdgeVal( _goalTri->getLvl2ToLvl2Index( tCurrent ) ), _radius );
        }
        if( goalReached == true && neighbGValue >= memGValue ){
          return getFunnel( _goalTri );
        }else if( tNeighbour == _goalTri ){
          if( goalReached == false || neighbGValue < memGValue ){
            memGValue = neighbGValue;
            //tNeighbour->setSearchNode( ETP::SearchNode<P,D>( 0.0, memGValue, tCurrent ), usedList );
            setTriangleSearchNode( tNeighbour, ETP::SearchNode<P,D>( 0.0, memGValue, tCurrent ), usedList );
            if( goalReached == true ){
              return getFunnel( _goalTri );
            }
            goalReached = true;
          }
        }else{
          D neighbHValue = tNeighbCentroid.dist( goalCentroid );
          setTriangleSearchNode( tNeighbour, ETP::SearchNode<P,D>( neighbHValue, neighbGValue, tCurrent ), usedList );
          insertInOpenList( neighbHValue, tNeighbour, openList );
        }
      }
    }
  }
  std::vector<std::shared_ptr<ETP::Triangle<P,D>>> searchInsideLvl2Corridor( ETP::Point<P,D> _startPoint, ETP::Point<P,D> _goalPoint, std::shared_ptr<ETP::Triangle<P,D>> _startTri, std::shared_ptr<ETP::Triangle<P,D>> _goalTri, D _radius ){

    std::shared_ptr<ETP::Triangle<P,D>> endPoint = nullptr;
    D sDist;
    int sEndPointIdx;
    bool startPtInStartTri = _startTri->isPointInside( _startPoint );
    bool goalPtInGoalTri = _goalTri->isPointInside( _goalPoint );

    for( int i = 0; i < 3; i++ ){
      std::shared_ptr<ETP::Triangle<P,D>> cn = _startTri->getConnectedNode( i );
      if( cn != nullptr ){
        endPoint = cn;
        sDist = _startTri->getLowerBound( i ) * _radius + ( startPtInStartTri ?  getPointToEdgeLowerBound( _startPoint, _startTri->getEdgeVal( i ), _radius ) : 0.0 );
        sEndPointIdx = i;
        break;
      }
    }
    int gEndPointIdx = _goalTri->getConnectedNodeIndex( endPoint );
    D gDist = _goalTri->getLowerBound( gEndPointIdx ) * _radius + ( goalPtInGoalTri ? getPointToEdgeLowerBound( _goalPoint, _goalTri->getEdgeVal( gEndPointIdx ), _radius ) : 0.0 );
    std::shared_ptr<ETP::Triangle<P,D>> searchStart;
    std::shared_ptr<ETP::Triangle<P,D>> searchGoal;
    bool reverse;
    if( sDist >= gDist ){
      searchStart = _startTri;
      searchGoal = _goalTri;
      reverse = true;
    }else{
      searchStart = _goalTri;
      searchGoal = _startTri;
      reverse = false;
      D tmpDist = gDist;
      gDist = sDist;
      sDist = tmpDist;
    }

    bool goalReached = false;
    std::vector<std::shared_ptr<ETP::Triangle<P,D>>> openList;
    openList.push_back( searchStart );
    D startToGoalGValue;
    while( goalReached == false ){
      std::shared_ptr<ETP::Triangle<P,D>> tCurrent = openList.front();
      std::shared_ptr<ETP::Triangle<P,D>> tNeighbour = tCurrent->getAdjacentVal( tCurrent->getConnectedNodeIndex( endPoint ) );
      if( tNeighbour == searchGoal ){
        openList.insert( openList.begin(), tNeighbour );
        goalReached = true;
        startToGoalGValue = sDist - gDist;
        break;
      }else if( tNeighbour->getLevel() == 3 ){
        break;
      }else{
        openList.insert( openList.begin(), tNeighbour );
      }
    }
    if( goalReached == false ){
      return std::vector<std::shared_ptr<ETP::Triangle<P,D>>>();
    }else{
      std::shared_ptr<ETP::Triangle<P,D>> otherEndPoint;
      for( int i = 0; i < 3; i++ ){
        std::shared_ptr<ETP::Triangle<P,D>> cNode = searchStart->getConnectedNode( i );
        if( cNode != nullptr && cNode != endPoint ){
          otherEndPoint = cNode;
          break;
        }
      }
      //int otherEndPtToEndPtIdx = otherEndPoint->getConnectedNodeIndex( endPoint );
      int otherEndPtToEndPtIdx = searchStart->getReversedConnectedNodeIndex( otherEndPoint );
      D baseGValue = otherEndPoint->getLowerBound( otherEndPtToEndPtIdx ) * _radius - startToGoalGValue;
      /*
      ///---///
      return std::vector<std::shared_ptr<ETP::Triangle<P,D>>>{ endPoint, otherEndPoint };
      ///---///
      triangles[ 0 ]->setAngle( 0, sDist );
      triangles[ 0 ]->setAngle( 1,  gDist );
      triangles[ 0 ]->setAngle( 1,  otherEndPoint->getLowerBound( otherEndPtToEndPtIdx ) );
      return std::vector<std::shared_ptr<ETP::Triangle<P,D>>>{triangles[0]};
      ///---///
      */
      if( baseGValue >= startToGoalGValue ){
        if( reverse == true ){
          std::reverse( openList.begin(), openList.end( ) );
        }
        return openList;
      }
      std::vector<std::shared_ptr<ETP::Triangle<P,D>>> secOpenList;
      //ETP::UsedList<P,D> usedList = ETP::UsedList<P,D>();
std::shared_ptr<ETP::UsedList<P,D>> usedList = std::make_shared<ETP::UsedList<P,D>>();
      setTriangleSearchNode( otherEndPoint, ETP::SearchNode<P,D>( 0.0, baseGValue, endPoint ), usedList );
      secOpenList.push_back( otherEndPoint );
      ETP::Point<D,D> goalCentroid = endPoint->getCentroid();
      while( secOpenList.empty() == false ){
        std::shared_ptr<ETP::Triangle<P,D>> tCurrent = secOpenList.front();
        secOpenList.erase( secOpenList.begin() );
        ETP::SearchNode<P,D> tsn = *tCurrent->getSearchNode().get();
        std::shared_ptr<ETP::Edge<P,D>> edgeIn;
        if( tCurrent->getLevel() == 3 ){
          edgeIn = tCurrent->getEdgeVal( tCurrent->getConnectedNodeIndex( tsn.comeFrom ) );
        }
        for( int i = 0; i < 3; i++ ){
          std::shared_ptr<ETP::Triangle<P,D>> tConnectedNode = tCurrent->getConnectedNode( i );
          if( tConnectedNode == nullptr || ( tConnectedNode == tsn.comeFrom && tCurrent == otherEndPoint && i == otherEndPtToEndPtIdx ) || ( tConnectedNode == tsn.comeFrom && tCurrent != otherEndPoint ) ){ continue; }
          D tCNgValue = tCurrent->getLowerBound( i ) * _radius + tsn.gValue;
          if( tCurrent->getLevel() == 3 ){
            tCNgValue += tCurrent->getAngle( tCurrent->getEdgeVal( i ), edgeIn ) * _radius;
          }
          std::shared_ptr<ETP::SearchNode<P,D>> tCsN = tConnectedNode->getSearchNode();
          if( tCNgValue >= startToGoalGValue || ( tCsN != nullptr && tCsN->gValue <= tCNgValue ) ){ continue; }
          D tCNhValue = tConnectedNode->getCentroid().dist( goalCentroid );
          if( tConnectedNode == endPoint ){
            tCNgValue += endPoint->getAngle( endPoint->getEdgeVal( endPoint->getConnectedNodeIndex( tCurrent ) ), endPoint->getEdgeVal( endPoint->getConnectedNodeIndex( otherEndPoint ) ) ) * _radius;
            if( tCsN == nullptr ){
              setTriangleSearchNode( tConnectedNode, ETP::SearchNode<P,D>( tCNhValue, tCNgValue, tCurrent ), usedList );
            }else if( tCsN->gValue > tCNgValue ){
              tCsN->hValue = tCNhValue;
              tCsN->gValue = tCNgValue;
              tCsN->comeFrom = tCurrent;
            }
          }else{
            //tConnectedNode->setSearchNode( ETP::SearchNode<P,D>( tCNhValue, tCNgValue, tCurrent ), usedList );
            setTriangleSearchNode( tConnectedNode, ETP::SearchNode<P,D>( tCNhValue, tCNgValue, tCurrent ), usedList );
            insertInOpenList( tCNhValue, tConnectedNode, secOpenList );
          }
        }
      }
      /*
      ///---///
      triangles[ 0 ]->setAngle( 0, startToGoalGValue );
      triangles[ 0 ]->setAngle( 1,  endPoint->getSearchNode()->gValue );
      return std::vector<std::shared_ptr<ETP::Triangle<P,D>>>{triangles[0]};
      ///---///
      */
      if( endPoint->getSearchNode() != nullptr && endPoint->getSearchNode()->gValue < startToGoalGValue ){
        setTriangleSearchNode( searchGoal, ETP::SearchNode<P,D>( 0.0, 0.0, endPoint ), usedList );
        setTriangleSearchNode( otherEndPoint, ETP::SearchNode<P,D>( 0.0, 0.0, searchStart ), usedList );
        std::vector<std::shared_ptr<ETP::Triangle<P,D>>> path = getFunnelWithIntermediateNodes( searchGoal, searchStart );
        if( searchGoal == _startTri ){
          std::reverse( path.begin(), path.end( ) );
        }
        return path;
      }else{
        if( reverse == true ){
          std::reverse( openList.begin(), openList.end( ) );
        }
        return openList;
      }
    }
  }
  std::vector<std::shared_ptr<ETP::searchEndPoint<P,D>>> getTriangleEndPoints( std::shared_ptr<ETP::Triangle<P,D>> _tri, ETP::Point<P,D> _goalPoint, D _radius ){
  std::vector<std::shared_ptr<ETP::searchEndPoint<P,D>>> ret;
  //  ETP::Triangle<P,D> tri = *_tri.get();
  std::shared_ptr<ETP::searchEndPoint<P,D>> finalEndPt = std::make_shared<ETP::searchEndPoint<P,D>>( 0.0, _tri );
  int lvl = _tri->getLevel();
  if( lvl == 3 ){
    ret.push_back( finalEndPt );
  }else if( lvl == 2 ){
    for( int i = 0; i < 3; i++ ){
      std::shared_ptr<ETP::Triangle<P,D>> connectedNode = _tri->getConnectedNode( i );
      if( connectedNode != nullptr ){
        bool isNew = true;
        D lowerBound = getPointToEdgeLowerBound( _goalPoint, _tri->getEdgeVal( i ), _radius );
        size_t retl = ret.size();
        for ( auto ePt : ret ){
          if( connectedNode == ePt->triangle ){
            isNew = false;
            if( lowerBound < ePt->gValue ){
              ePt = std::make_shared<ETP::searchEndPoint<P,D>>( lowerBound, connectedNode, finalEndPt, _tri->getReversedConnectedNodeIndex( connectedNode ) );
            }
          }
        }
        if( isNew == true ){
          ret.push_back( std::make_shared<ETP::searchEndPoint<P,D>>( lowerBound, connectedNode, finalEndPt, _tri->getReversedConnectedNodeIndex( connectedNode ) ) );
        }
      }
    }
  }else if( lvl == 1 ){
    for( int i = 0; i < 3; i++ ){
      std::shared_ptr<ETP::Triangle<P,D>> connectedNode = _tri->getConnectedNode( i );
      D startLowerBound = getPointToEdgeLowerBound( _goalPoint, _tri->getEdgeVal( i ), _radius );
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
                  ret[ k ] = std::make_shared<ETP::searchEndPoint<P,D>>( lvl2EndPt->gValue + lowerBound, lvl3Node, lvl2EndPt );
                }
                break;
              }
            }
            if( isNew == true ){
              ret.push_back( std::make_shared<ETP::searchEndPoint<P,D>>( lvl2EndPt->gValue + lowerBound, lvl3Node, lvl2EndPt, connectedNode->getReversedConnectedNodeIndex( lvl3Node ) ) );
            }
          }
        }
        break;
      }
    }
  }
  return ret;
  };




  std::vector<std::shared_ptr<ETP::Triangle<P,D>>> abstractTriangulationSearch( ETP::Point<P,D> _startPoint, ETP::Point<P,D> _goalPoint ){
    std::shared_ptr<ETP::Triangle<P,D>> _startTri = getTriangleWithPoint( _startPoint );
    std::shared_ptr<ETP::Triangle<P,D>> _goalTri = getTriangleWithPoint( _goalPoint );
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
      if(    ( startRoot == nullptr ) //start and goal belong to the same unrooted tree
          || ( startRoot == _goalTri->getLvl1RootNode() ) // start and goal belong to the same rooted tree
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
        bool proceed = false;
        std::vector<std::shared_ptr<ETP::Triangle<P,D>>> bodySearch;
        if( startRoot->isPartoflvl2Ring() || startRoot->isOnSameLvl2Loop( goalRoot ) ){// inside a loop or ring
          proceed = true;
          bodySearch =  searchInsideLvl2RingOrLoop( _startPoint, _goalPoint, _startTri, _goalTri );
        }else if( startRoot->haveSameLvl2CorridorEndpoints( goalRoot )){// inside a corridor
          //return std::vector<std::shared_ptr<ETP::Triangle<P,D>>>{ _startTri, _goalTri };
          bodySearch = searchInsideLvl2Corridor( _startPoint, _goalPoint, startRoot, goalRoot );
          if( bodySearch.empty() == false ){
            proceed = true;
          }
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

    //tmp
    return searchClassic( _startPoint, _goalPoint, _startTri, _goalTri );

  }
  std::vector<std::shared_ptr<ETP::Triangle<P,D>>> abstractTriangulationSearch( std::shared_ptr<ETP::Triangle<P,D>> _startTri, std::shared_ptr<ETP::Triangle<P,D>> _goalTri ){
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
      if(    ( startRoot == nullptr ) //start and goal belong to the same unrooted tree
          || ( startRoot == _goalTri->getLvl1RootNode() ) // start and goal belong to the same rooted tree
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
        bool proceed = false;
        std::vector<std::shared_ptr<ETP::Triangle<P,D>>> bodySearch;
        if( startRoot->isPartoflvl2Ring() || startRoot->isOnSameLvl2Loop( goalRoot ) ){// inside a loop or ring
          proceed = true;
          bodySearch =  searchInsideLvl2RingOrLoop( _startTri, _goalTri );
        }else if( startRoot->haveSameLvl2CorridorEndpoints( goalRoot )){// inside a corridor
          //return std::vector<std::shared_ptr<ETP::Triangle<P,D>>>{ _startTri, _goalTri };
          bodySearch = searchInsideLvl2Corridor( startRoot, goalRoot );
          if( bodySearch.empty() == false ){
            proceed = true;
          }
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

    //tmp
    return searchClassic( _startTri, _goalTri );

  }

  std::vector<std::shared_ptr<ETP::Triangle<P,D>>> searchClassic( ETP::Point<P,D> _startPoint, ETP::Point<P,D> _goalPoint, std::shared_ptr<ETP::Triangle<P,D>> _startTri, std::shared_ptr<ETP::Triangle<P,D>> _goalTri ){
    //ETP::Point<D,D> goalCentroid = _goalTri->getCentroid();
    ETP::Point<D,D> goalCentroid = ETP::Point<D,D>( (D) _goalPoint.x, (D) _goalPoint.y );
    std::vector<std::shared_ptr<ETP::Triangle<P,D>>> openList;
    //ETP::UsedList<P,D> usedList = ETP::UsedList<P,D>();
    std::shared_ptr<ETP::UsedList<P,D>> usedList = std::make_shared<ETP::UsedList<P,D>>();
    std::vector<std::shared_ptr<ETP::searchEndPoint<P,D>>> goalEndPoints = getTriangleEndPoints( _goalTri, _goalPoint );
    setTriangleSearchNode( _startTri, ETP::SearchNode<P,D>( 0.0, 0.0, nullptr ), usedList );
    openList.push_back( _startTri );
    size_t gel = goalEndPoints.size();

    D memLomerGValue;
    bool endPointFound = false;

    while( openList.empty() == false ){
      std::shared_ptr<ETP::Triangle<P,D>> tCurrent = openList.front();
      openList.erase( openList.begin() );
      std::shared_ptr<ETP::SearchNode<P,D>> searchNodeCurrent = tCurrent->getSearchNode();
      D tCGVal = searchNodeCurrent->gValue;
      for( int i = 0; i < 3; i++ ){
        std::shared_ptr<ETP::Triangle<P,D>> tNeighbour = tCurrent->getConnectedNode( i );
        if( tNeighbour == nullptr || tNeighbour == searchNodeCurrent->comeFrom ){ continue; }
        D gVal = tCGVal + tCurrent->getLowerBound( i );
        /*
        if( tCurrent->getLevel() == 3 && searchNodeCurrent->comeFrom != nullptr ){
          gVal += tCurrent->getAngle( tCurrent->getEdgeVal( tCurrent->getConnectedNodeIndex( searchNodeCurrent->comeFrom ) ), tCurrent->getEdgeVal( i ) );
        }
        */
        if( tCurrent == _startTri ){
          gVal += getPointToEdgeLowerBound(  _startPoint, tCurrent->getEdgeVal( i ) );
        }else if( tCurrent->getLevel() == 3 ){
          if( searchNodeCurrent->comeFrom != nullptr ){
            gVal += tCurrent->getAngle( tCurrent->getEdgeVal( tCurrent->getConnectedNodeIndex( searchNodeCurrent->comeFrom ) ), tCurrent->getEdgeVal( i ) );
          }
        }
        if( endPointFound == true && gVal >= memLomerGValue ){ continue; }
        bool neighbIsEndPoint = false;
        for( auto gEndPt : goalEndPoints ){
          if( gEndPt->triangle == tNeighbour ){
            neighbIsEndPoint = true;
            int neighbToCurrentIdx = tNeighbour->getConnectedNodeIndex( tCurrent );
            if( tNeighbour == _goalTri && tNeighbour->getLevel() == 3 ){
              gVal += getPointToEdgeLowerBound( _goalPoint, _goalTri->getEdgeVal( neighbToCurrentIdx ) );
            }else{
              gVal += tNeighbour->getAngle( tNeighbour->getEdgeVal( gEndPt->parentIndex ), tNeighbour->getEdgeVal( neighbToCurrentIdx ) );
            }
            if( endPointFound == false ){
              //remove triangles from openList which have a gValue higher than the endPoint's gValue
              openList.erase( std::remove_if( openList.begin(), openList.end(),
                [ gVal ]( std::shared_ptr<ETP::Triangle<P,D>> _t ){
                    if ( _t->getSearchNode()->gValue >= gVal ) {
                        return true;
                    }
                    return false;
                }
              ), openList.end() );
            }

            endPointFound = true;
            memLomerGValue = gVal;
            break;
          }
        }
        D hVal = tNeighbour->getCentroid().dist( goalCentroid );
        std::shared_ptr<ETP::SearchNode<P,D>> tNeighbSearchNode = tNeighbour->getSearchNode();
        if( tNeighbSearchNode == nullptr ){
          setTriangleSearchNode( tNeighbour, ETP::SearchNode<P,D>( hVal, gVal, tCurrent ), usedList );
          if( neighbIsEndPoint == false ){ insertInOpenList( hVal, tNeighbour, openList ); }
        }else if( tNeighbSearchNode->gValue > gVal ){
          tNeighbour->setSearchHValue( hVal );
          tNeighbour->setSearchGValue( gVal );
          tNeighbour->setSearchComeFrom( tCurrent );
          if( neighbIsEndPoint == false ){
            removeFromOpenList( tNeighbour, openList );
            insertInOpenList( hVal, tNeighbour, openList );
          }
        }
      }
    }

    std::shared_ptr<ETP::searchEndPoint<P,D>> closestEndPoint = nullptr;
    bool closestFound = false;
    int closestSearchIndex;
    for( int i = 0; i < gel; i++ ){
      std::shared_ptr<ETP::searchEndPoint<P,D>> testedClosestEndPoint = goalEndPoints[ i ];
      std::shared_ptr<ETP::SearchNode<P,D>> testedClosestSearchNode = testedClosestEndPoint->triangle->getSearchNode();
      if( testedClosestSearchNode != nullptr && ( closestEndPoint == nullptr || closestEndPoint->gValue + closestEndPoint->triangle->getSearchNode()->gValue > testedClosestEndPoint->gValue + testedClosestEndPoint->gValue ) ){
        closestEndPoint = testedClosestEndPoint;
        closestFound = true;
      }
    }
    if( closestFound == false ){
      return std::vector<std::shared_ptr<ETP::Triangle<P,D>>>();
    }
    bool endPointReached = false;
    std::vector<std::shared_ptr<ETP::Triangle<P,D>>> funnelEndPart = getFunnelToEndPoint( closestEndPoint );
    std::shared_ptr<ETP::Triangle<P,D>> tCurrent = closestEndPoint->triangle;
    while( endPointReached == false ){
      std::shared_ptr<ETP::SearchNode<P,D>> tSN = tCurrent->getSearchNode();
      if( tSN == nullptr || tSN->comeFrom == nullptr){
        endPointReached = true;
        break;
      }
      std::shared_ptr<ETP::Triangle<P,D>> tComeFrom = tSN->comeFrom;
      std::shared_ptr<ETP::Triangle<P,D>> tNode = tComeFrom;
      bool partEnd = false;
      std::vector<std::shared_ptr<ETP::Triangle<P,D>>> tmpFunnel;
      while( partEnd == false ){
        if( tNode == tCurrent ){
          break;
        }
        tmpFunnel.push_back( tNode );
        tNode = tNode->getAdjacentVal( tNode->getConnectedNodeIndex( tCurrent ) );

      }
      funnelEndPart.insert( funnelEndPart.begin(), tmpFunnel.begin(), tmpFunnel.end() );
      tCurrent = tComeFrom;
    }
    return funnelEndPart;
  }

  std::vector<std::shared_ptr<ETP::Triangle<P,D>>> searchClassic( std::shared_ptr<ETP::Triangle<P,D>> _startTri, std::shared_ptr<ETP::Triangle<P,D>> _goalTri ){
    ETP::Point<D,D> goalCentroid = _goalTri->getCentroid();
    std::vector<std::shared_ptr<ETP::Triangle<P,D>>> openList;
    //ETP::UsedList<P,D> usedList = ETP::UsedList<P,D>();
    std::shared_ptr<ETP::UsedList<P,D>> usedList = std::make_shared<ETP::UsedList<P,D>>();
    std::vector<std::shared_ptr<ETP::searchEndPoint<P,D>>> goalEndPoints = getTriangleEndPoints( _goalTri );
    setTriangleSearchNode( _startTri, ETP::SearchNode<P,D>( 0.0, 0.0, nullptr ), usedList );
    openList.push_back( _startTri );
    size_t gel = goalEndPoints.size();

    D memLomerGValue;
    bool endPointFound = false;

    while( openList.empty() == false ){
      std::shared_ptr<ETP::Triangle<P,D>> tCurrent = openList.front();
      openList.erase( openList.begin() );
      std::shared_ptr<ETP::SearchNode<P,D>> searchNodeCurrent = tCurrent->getSearchNode();
      D tCGVal = searchNodeCurrent->gValue;
      for( int i = 0; i < 3; i++ ){
        std::shared_ptr<ETP::Triangle<P,D>> tNeighbour = tCurrent->getConnectedNode( i );
        if( tNeighbour == nullptr || tNeighbour == searchNodeCurrent->comeFrom ){ continue; }
        D gVal = tCGVal + tCurrent->getLowerBound( i );
        if( tCurrent->getLevel() == 3 && searchNodeCurrent->comeFrom != nullptr ){
          gVal += tCurrent->getAngle( tCurrent->getEdgeVal( tCurrent->getConnectedNodeIndex( searchNodeCurrent->comeFrom ) ), tCurrent->getEdgeVal( i ) );
        }
        if( endPointFound == true && gVal >= memLomerGValue ){ continue; }
        bool neighbIsEndPoint = false;
        for( auto gEndPt : goalEndPoints ){
          if( gEndPt->triangle == tNeighbour ){
            neighbIsEndPoint = true;
            if( endPointFound == false ){
              //remove triangles from openList which have a gValue higher than the endPoint's gValue
              openList.erase( std::remove_if( openList.begin(), openList.end(),
                [ gVal ]( std::shared_ptr<ETP::Triangle<P,D>> _t ){
                    if ( _t->getSearchNode()->gValue >= gVal ) {
                        return true;
                    }
                    return false;
                }
              ), openList.end() );
            }

            endPointFound = true;
            memLomerGValue = gVal;
            break;
          }
        }
        D hVal = tNeighbour->getCentroid().dist( goalCentroid );
        std::shared_ptr<ETP::SearchNode<P,D>> tNeighbSearchNode = tNeighbour->getSearchNode();
        if( tNeighbSearchNode == nullptr ){
          setTriangleSearchNode( tNeighbour, ETP::SearchNode<P,D>( hVal, gVal, tCurrent ), usedList );
          if( neighbIsEndPoint == false ){ insertInOpenList( hVal, tNeighbour, openList ); }
        }else if( tNeighbSearchNode->gValue > gVal ){
          tNeighbour->setSearchHValue( hVal );
          tNeighbour->setSearchGValue( gVal );
          tNeighbour->setSearchComeFrom( tCurrent );
          if( neighbIsEndPoint == false ){
            removeFromOpenList( tNeighbour, openList );
            insertInOpenList( hVal, tNeighbour, openList );
          }
        }
      }
    }

    std::shared_ptr<ETP::searchEndPoint<P,D>> closestEndPoint = nullptr;
    bool closestFound = false;
    int closestSearchIndex;
    for( int i = 0; i < gel; i++ ){
      std::shared_ptr<ETP::searchEndPoint<P,D>> testedClosestEndPoint = goalEndPoints[ i ];
      std::shared_ptr<ETP::SearchNode<P,D>> testedClosestSearchNode = testedClosestEndPoint->triangle->getSearchNode();
      if( testedClosestSearchNode != nullptr && ( closestEndPoint == nullptr || closestEndPoint->gValue + closestEndPoint->triangle->getSearchNode()->gValue > testedClosestEndPoint->gValue + testedClosestEndPoint->gValue ) ){
        closestEndPoint = testedClosestEndPoint;
        closestFound = true;
      }
    }
    if( closestFound == false ){
      return std::vector<std::shared_ptr<ETP::Triangle<P,D>>>();
    }
    bool endPointReached = false;
    std::vector<std::shared_ptr<ETP::Triangle<P,D>>> funnelEndPart = getFunnelToEndPoint( closestEndPoint );
    std::shared_ptr<ETP::Triangle<P,D>> tCurrent = closestEndPoint->triangle;
    while( endPointReached == false ){
      std::shared_ptr<ETP::SearchNode<P,D>> tSN = tCurrent->getSearchNode();
      if( tSN == nullptr || tSN->comeFrom == nullptr){
        endPointReached = true;
        break;
      }
      std::shared_ptr<ETP::Triangle<P,D>> tComeFrom = tSN->comeFrom;
      std::shared_ptr<ETP::Triangle<P,D>> tNode = tComeFrom;
      bool partEnd = false;
      std::vector<std::shared_ptr<ETP::Triangle<P,D>>> tmpFunnel;
      while( partEnd == false ){
        if( tNode == tCurrent ){
          break;
        }
        tmpFunnel.push_back( tNode );
        tNode = tNode->getAdjacentVal( tNode->getConnectedNodeIndex( tCurrent ) );

      }
      funnelEndPart.insert( funnelEndPart.begin(), tmpFunnel.begin(), tmpFunnel.end() );
      tCurrent = tComeFrom;
    }
    return funnelEndPart;
  }

  std::string  abstractTriangulationSearchString( std::shared_ptr<ETP::Triangle<P,D>> _startTri, std::shared_ptr<ETP::Triangle<P,D>> _goalTri ){
      std::vector<std::shared_ptr<ETP::Triangle<P,D>>> ret;
      if( _startTri == nullptr || _goalTri == nullptr ){
        //return ret;
        return std::string( "start or goal triangle is null" );
      }
      int sComp = _startTri->getComponent();
      int gComp = _goalTri->getComponent();
      if( sComp != gComp ){
        //return ret;
        return std::string( "on different Components" );
      }
      // on same triangle or within lvl 0 node
      if( _startTri == _goalTri ){
        //ret.push_back( _startTri );
        return std::string( "on same triangle or within lvl 0 node" );
      }
      int sLvl = _startTri->getLevel();
      int gLvl = _goalTri->getLevel();

      //start and goal are connected lvl 1 and lvl 2 nodes
      if( sLvl == 1 && gLvl == 2 ){ // start is lvl 1 and goal is lvl 2
        int connectedNodeIndex = _startTri->getConnectedNodeIndex( _goalTri );
        if( connectedNodeIndex > -1 ){
          //return _startTri->getTrianglesToConnectedNode( connectedNodeIdx );
          return std::string( "start and goal are connected lvl 1 and lvl 2 nodes - start is lvl 1 and goal is lvl 2" );
        }
      }else if( sLvl == 2 && gLvl == 1 ){ // start is lvl 2 and goal is lvl 1
        int connectedNodeIndex = _goalTri->getConnectedNodeIndex( _startTri );
        if( connectedNodeIndex > -1 ){
          //ret = _goalTri->getTrianglesToConnectedNode( connectedNodeIdx );
          //std::reverse(std::begin( ret ), std::end( ret ) );
          //return ret;
          return std::string( "start and goal are connected lvl 1 and lvl 2 nodes - start is lvl 2 and goal is lvl 1" );
        }
      }
      //start and goal are inside the same lvl 1 tree
      if( sLvl == 1 && gLvl == 1 ){
        std::shared_ptr<ETP::Triangle<P,D>> startRoot = _startTri->getLvl1RootNode();
        if(    ( startRoot == nullptr ) //start and goal belong to the same unrooted tree
            || ( startRoot == _goalTri->getLvl1RootNode() ) // start and goal belong to the same rooted tree
        ){
          return std::string( "start and goal are inside the same lvl 1 tree" );
        }
      }
      // start and goal belong either to the same level 2 loop or ring or to connected lvl 1 trees connected to it
      // includes also the case where start and goal belong to a level 2 corridor or to attached lvl 1 trees
      if( sLvl == 2 || ( sLvl == 1 && _startTri->getLvl1RootNode() != nullptr ) ){
        std::shared_ptr<ETP::Triangle<P,D>> startRoot = sLvl == 2 ? _startTri : _startTri->getLvl1RootNode();
        if( gLvl == 2 || ( gLvl == 1 && _goalTri->getLvl1RootNode() != nullptr ) ){
          std::shared_ptr<ETP::Triangle<P,D>> goalRoot = gLvl == 2 ? _goalTri : _goalTri->getLvl1RootNode();
          if( startRoot->isPartoflvl2Ring() || startRoot->isOnSameLvl2Loop( goalRoot ) ){// inside a loop or ring
            if( sLvl == 1 ){ // add triangles between start Lvl1 triangle and its lvl2 root

            }
            if( gLvl == 1 ){// add triangles between goal Lvl1 triangle and its lvl2 root

            }
            // add triangles between the two connected lvl 2 triangles
            return std::string( "start and goal belong either to the same level 2 loop or ring or to connected lvl 1 trees connected to it" );
          }else if( startRoot->isOnSameLvl2Corridor( goalRoot )){// inside a corridor
            return std::string( "start and goal belong to the same level 2 inside a corridor or to lvl 1 trees connected to it" );
          }
        }
      }
      // perform classic Triangulation Reduction A*
      return std::string( "start and goal match to no special case and a class TRA* is needed" );
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

  std::vector<std::shared_ptr<ETP::Triangle<P,D>>> searchInsideLvl2RingOrLoop( ETP::Point<P,D> _startPoint, ETP::Point<P,D> _goalPoint, std::shared_ptr<ETP::Triangle<P,D>> _startTri, std::shared_ptr<ETP::Triangle<P,D>> _goalTri ){

    std::vector<std::shared_ptr<ETP::Triangle<P,D>>> openList;
    //std::shared_ptr<ETP::UsedList<P,D>> usedList = std::make_shared<ETP::UsedList<P,D>>();
    //ETP::UsedList<P,D> usedList = ETP::UsedList<P,D>();
    std::shared_ptr<ETP::UsedList<P,D>> usedList = std::make_shared<ETP::UsedList<P,D>>();
    ETP::Point<D,D> goalCentroid = _goalTri->getCentroid();
    setTriangleSearchNode( _startTri, ETP::SearchNode<P,D>( 0.0, 0.0 ), usedList );
    openList.push_back( _startTri );
    bool goalReached = false;
    D memGValue;
    while( openList.empty() == false ){
      std::shared_ptr<ETP::Triangle<P,D>> tCurrent = openList.front();
      openList.erase( openList.begin() );
      ETP::SearchNode<P,D> snCurrent = *tCurrent->getSearchNode().get();
      ETP::Point<D,D> currentCentroid = tCurrent->getCentroid();
      int currentLvl = tCurrent->getLevel();
      for( int i = 0; i < 3; i++ ){
        if( currentLvl == 3 && tCurrent->getConnectedNode( i ) != tCurrent ){ continue; }
        std::shared_ptr<ETP::Triangle<P,D>> tNeighbour = tCurrent->getAdjacentVal( i );
        int neighbLvl = tNeighbour->getLevel();
        if( tNeighbour == nullptr || neighbLvl < 2 || tNeighbour->getSearchNode() != nullptr ){ continue; }
        //ETP::Point<D,D> tNeighbCentroid = tNeighbour->getCentroid();
        ETP::Point<D,D> tNeighbCentroid = tNeighbour->getCentroid();
        D neighbGValue = snCurrent.gValue;
        if( tCurrent == _startTri ){
          neighbGValue += getPointToEdgeLowerBound(  _startPoint, _startTri->getEdgeVal( i ) );
        }else{
          std::shared_ptr<ETP::Edge<P,D>> edgeOut = tCurrent->getEdgeVal( i );
          std::shared_ptr<ETP::Edge<P,D>> edgeIn;
          for( int e = 0; e < 3; e++ ){
            std::shared_ptr<ETP::Edge<P,D>> te = tCurrent->getEdgeVal( e );
            if( te->isConstrained() == false && te != edgeOut ){
              edgeIn = te;
              break;
            }
          }
          neighbGValue += tCurrent->getAngle( edgeIn, edgeOut );
        }
        if( tNeighbour == _goalTri ){
          neighbGValue += getPointToEdgeLowerBound(  _goalPoint, _goalTri->getEdgeVal( _goalTri->getLvl2ToLvl2Index( tCurrent ) ) );
        }
        if( goalReached == true && neighbGValue >= memGValue ){
          return getFunnel( _goalTri );
        }else if( tNeighbour == _goalTri ){
          if( goalReached == false || neighbGValue < memGValue ){
            memGValue = neighbGValue;
            //tNeighbour->setSearchNode( ETP::SearchNode<P,D>( 0.0, memGValue, tCurrent ), usedList );
            setTriangleSearchNode( tNeighbour, ETP::SearchNode<P,D>( 0.0, memGValue, tCurrent ), usedList );
            if( goalReached == true ){
              return getFunnel( _goalTri );
            }
            goalReached = true;
          }
        }else{
          D neighbHValue = tNeighbCentroid.dist( goalCentroid );
          setTriangleSearchNode( tNeighbour, ETP::SearchNode<P,D>( neighbHValue, neighbGValue, tCurrent ), usedList );
          insertInOpenList( neighbHValue, tNeighbour, openList );
        }
      }
    }
  }

  std::vector<std::shared_ptr<ETP::Triangle<P,D>>> searchInsideLvl2RingOrLoop( std::shared_ptr<ETP::Triangle<P,D>> _startTri, std::shared_ptr<ETP::Triangle<P,D>> _goalTri ){
    std::vector<std::shared_ptr<ETP::Triangle<P,D>>> openList;
    //std::shared_ptr<ETP::UsedList<P,D>> usedList = std::make_shared<ETP::UsedList<P,D>>();
    //ETP::UsedList<P,D> usedList = ETP::UsedList<P,D>();
    std::shared_ptr<ETP::UsedList<P,D>> usedList = std::make_shared<ETP::UsedList<P,D>>();
    ETP::Point<D,D> goalCentroid = _goalTri->getCentroid();
    setTriangleSearchNode( _startTri, ETP::SearchNode<P,D>( 0.0, 0.0 ), usedList );
    openList.push_back( _startTri );
    bool goalReached = false;
    D memGValue;
    while( openList.empty() == false ){
      std::shared_ptr<ETP::Triangle<P,D>> tCurrent = openList.front();
      openList.erase( openList.begin() );
      ETP::SearchNode<P,D> snCurrent = *tCurrent->getSearchNode().get();
      ETP::Point<D,D> currentCentroid = tCurrent->getCentroid();
      int currentLvl = tCurrent->getLevel();
      for( int i = 0; i < 3; i++ ){
        if( currentLvl == 3 && tCurrent->getConnectedNode( i ) != tCurrent ){ continue; }
        std::shared_ptr<ETP::Triangle<P,D>> tNeighbour = tCurrent->getAdjacentVal( i );
        int neighbLvl = tNeighbour->getLevel();
        if( tNeighbour == nullptr || neighbLvl < 2 || tNeighbour->getSearchNode() != nullptr ){ continue; }
        //ETP::Point<D,D> tNeighbCentroid = tNeighbour->getCentroid();
        ETP::Point<D,D> tNeighbCentroid = tNeighbour->getCentroid();
        D neighbGValue = snCurrent.gValue + currentCentroid.dist( tNeighbCentroid );
        if( goalReached == true && neighbGValue >= memGValue ){
          return getFunnel( _goalTri );
        }else if( tNeighbour == _goalTri ){
          if( goalReached == false || neighbGValue < memGValue ){
            memGValue = neighbGValue;
            //tNeighbour->setSearchNode( ETP::SearchNode<P,D>( 0.0, memGValue, tCurrent ), usedList );
            setTriangleSearchNode( tNeighbour, ETP::SearchNode<P,D>( 0.0, memGValue, tCurrent ), usedList );
            if( goalReached == true ){
              return getFunnel( _goalTri );
            }
            goalReached = true;
          }
        }else{
          D neighbHValue = tNeighbCentroid.dist( goalCentroid );
          setTriangleSearchNode( tNeighbour, ETP::SearchNode<P,D>( neighbHValue, neighbGValue, tCurrent ), usedList );
          insertInOpenList( neighbHValue, tNeighbour, openList );
        }
      }
    }
  }

  std::vector<std::shared_ptr<ETP::Triangle<P,D>>> searchInsideLvl2Corridor( ETP::Point<P,D> _startPoint, ETP::Point<P,D> _goalPoint, std::shared_ptr<ETP::Triangle<P,D>> _startTri, std::shared_ptr<ETP::Triangle<P,D>> _goalTri ){

    std::shared_ptr<ETP::Triangle<P,D>> endPoint = nullptr;
    D sDist;
    int sEndPointIdx;
    bool startPtInStartTri = _startTri->isPointInside( _startPoint );
    bool goalPtInGoalTri = _goalTri->isPointInside( _goalPoint );

    for( int i = 0; i < 3; i++ ){
      std::shared_ptr<ETP::Triangle<P,D>> cn = _startTri->getConnectedNode( i );
      if( cn != nullptr ){
        endPoint = cn;
        sDist = _startTri->getLowerBound( i ) + ( startPtInStartTri ?  getPointToEdgeLowerBound( _startPoint, _startTri->getEdgeVal( i ) ) : 0.0 );
        sEndPointIdx = i;
        break;
      }
    }
    int gEndPointIdx = _goalTri->getConnectedNodeIndex( endPoint );
    D gDist = _goalTri->getLowerBound( gEndPointIdx ) + ( goalPtInGoalTri ? getPointToEdgeLowerBound( _goalPoint, _goalTri->getEdgeVal( gEndPointIdx ) ) : 0.0 );
    std::shared_ptr<ETP::Triangle<P,D>> searchStart;
    std::shared_ptr<ETP::Triangle<P,D>> searchGoal;
    bool reverse;
    if( sDist >= gDist ){
      searchStart = _startTri;
      searchGoal = _goalTri;
      reverse = true;
    }else{
      searchStart = _goalTri;
      searchGoal = _startTri;
      reverse = false;
      D tmpDist = gDist;
      gDist = sDist;
      sDist = tmpDist;
    }

    bool goalReached = false;
    std::vector<std::shared_ptr<ETP::Triangle<P,D>>> openList;
    openList.push_back( searchStart );
    D startToGoalGValue;
    while( goalReached == false ){
      std::shared_ptr<ETP::Triangle<P,D>> tCurrent = openList.front();
      std::shared_ptr<ETP::Triangle<P,D>> tNeighbour = tCurrent->getAdjacentVal( tCurrent->getConnectedNodeIndex( endPoint ) );
      if( tNeighbour == searchGoal ){
        openList.insert( openList.begin(), tNeighbour );
        goalReached = true;
        startToGoalGValue = sDist - gDist;
        break;
      }else if( tNeighbour->getLevel() == 3 ){
        break;
      }else{
        openList.insert( openList.begin(), tNeighbour );
      }
    }
    if( goalReached == false ){
      return std::vector<std::shared_ptr<ETP::Triangle<P,D>>>();
    }else{
      std::shared_ptr<ETP::Triangle<P,D>> otherEndPoint;
      for( int i = 0; i < 3; i++ ){
        std::shared_ptr<ETP::Triangle<P,D>> cNode = searchStart->getConnectedNode( i );
        if( cNode != nullptr && cNode != endPoint ){
          otherEndPoint = cNode;
          break;
        }
      }
      //int otherEndPtToEndPtIdx = otherEndPoint->getConnectedNodeIndex( endPoint );
      int otherEndPtToEndPtIdx = searchStart->getReversedConnectedNodeIndex( otherEndPoint );
      D baseGValue = otherEndPoint->getLowerBound( otherEndPtToEndPtIdx ) - startToGoalGValue;
      /*
      ///---///
      return std::vector<std::shared_ptr<ETP::Triangle<P,D>>>{ endPoint, otherEndPoint };
      ///---///
      triangles[ 0 ]->setAngle( 0, sDist );
      triangles[ 0 ]->setAngle( 1,  gDist );
      triangles[ 0 ]->setAngle( 1,  otherEndPoint->getLowerBound( otherEndPtToEndPtIdx ) );
      return std::vector<std::shared_ptr<ETP::Triangle<P,D>>>{triangles[0]};
      ///---///
      */
      if( baseGValue >= startToGoalGValue ){
        if( reverse == true ){
          std::reverse( openList.begin(), openList.end( ) );
        }
        return openList;
      }
      std::vector<std::shared_ptr<ETP::Triangle<P,D>>> secOpenList;
      //ETP::UsedList<P,D> usedList = ETP::UsedList<P,D>();
std::shared_ptr<ETP::UsedList<P,D>> usedList = std::make_shared<ETP::UsedList<P,D>>();
      setTriangleSearchNode( otherEndPoint, ETP::SearchNode<P,D>( 0.0, baseGValue, endPoint ), usedList );
      secOpenList.push_back( otherEndPoint );
      ETP::Point<D,D> goalCentroid = endPoint->getCentroid();
      while( secOpenList.empty() == false ){
        std::shared_ptr<ETP::Triangle<P,D>> tCurrent = secOpenList.front();
        secOpenList.erase( secOpenList.begin() );
        ETP::SearchNode<P,D> tsn = *tCurrent->getSearchNode().get();
        std::shared_ptr<ETP::Edge<P,D>> edgeIn;
        if( tCurrent->getLevel() == 3 ){
          edgeIn = tCurrent->getEdgeVal( tCurrent->getConnectedNodeIndex( tsn.comeFrom ) );
        }
        for( int i = 0; i < 3; i++ ){
          std::shared_ptr<ETP::Triangle<P,D>> tConnectedNode = tCurrent->getConnectedNode( i );
          if( tConnectedNode == nullptr || ( tConnectedNode == tsn.comeFrom && tCurrent == otherEndPoint && i == otherEndPtToEndPtIdx ) || ( tConnectedNode == tsn.comeFrom && tCurrent != otherEndPoint ) ){ continue; }
          D tCNgValue = tCurrent->getLowerBound( i ) + tsn.gValue;
          if( tCurrent->getLevel() == 3 ){
            tCNgValue += tCurrent->getAngle( tCurrent->getEdgeVal( i ), edgeIn );
          }
          std::shared_ptr<ETP::SearchNode<P,D>> tCsN = tConnectedNode->getSearchNode();
          if( tCNgValue >= startToGoalGValue || ( tCsN != nullptr && tCsN->gValue <= tCNgValue ) ){ continue; }
          D tCNhValue = tConnectedNode->getCentroid().dist( goalCentroid );
          if( tConnectedNode == endPoint ){
            tCNgValue += endPoint->getAngle( endPoint->getEdgeVal( endPoint->getConnectedNodeIndex( tCurrent ) ), endPoint->getEdgeVal( endPoint->getConnectedNodeIndex( otherEndPoint ) ) );
            if( tCsN == nullptr ){
              setTriangleSearchNode( tConnectedNode, ETP::SearchNode<P,D>( tCNhValue, tCNgValue, tCurrent ), usedList );
            }else if( tCsN->gValue > tCNgValue ){
              tCsN->hValue = tCNhValue;
              tCsN->gValue = tCNgValue;
              tCsN->comeFrom = tCurrent;
            }
          }else{
            //tConnectedNode->setSearchNode( ETP::SearchNode<P,D>( tCNhValue, tCNgValue, tCurrent ), usedList );
            setTriangleSearchNode( tConnectedNode, ETP::SearchNode<P,D>( tCNhValue, tCNgValue, tCurrent ), usedList );
            insertInOpenList( tCNhValue, tConnectedNode, secOpenList );
          }
        }
      }
      /*
      ///---///
      triangles[ 0 ]->setAngle( 0, startToGoalGValue );
      triangles[ 0 ]->setAngle( 1,  endPoint->getSearchNode()->gValue );
      return std::vector<std::shared_ptr<ETP::Triangle<P,D>>>{triangles[0]};
      ///---///
      */
      if( endPoint->getSearchNode() != nullptr && endPoint->getSearchNode()->gValue < startToGoalGValue ){
        setTriangleSearchNode( searchGoal, ETP::SearchNode<P,D>( 0.0, 0.0, endPoint ), usedList );
        setTriangleSearchNode( otherEndPoint, ETP::SearchNode<P,D>( 0.0, 0.0, searchStart ), usedList );
        std::vector<std::shared_ptr<ETP::Triangle<P,D>>> path = getFunnelWithIntermediateNodes( searchGoal, searchStart );
        if( searchGoal == _startTri ){
          std::reverse( path.begin(), path.end( ) );
        }
        return path;
      }else{
        if( reverse == true ){
          std::reverse( openList.begin(), openList.end( ) );
        }
        return openList;
      }
    }
  }


  std::vector<std::shared_ptr<ETP::Triangle<P,D>>> searchInsideLvl2Corridor( std::shared_ptr<ETP::Triangle<P,D>> _startTri, std::shared_ptr<ETP::Triangle<P,D>> _goalTri ){
    std::shared_ptr<ETP::Triangle<P,D>> endPoint = nullptr;
    D sDist;
    int sEndPointIdx;
    for( int i = 0; i < 3; i++ ){
      std::shared_ptr<ETP::Triangle<P,D>> cn = _startTri->getConnectedNode( i );
      if( cn != nullptr ){
        endPoint = cn;
        sDist = _startTri->getLowerBound( i );
        sEndPointIdx = i;
        break;
      }
    }
    int gEndPointIdx = _goalTri->getConnectedNodeIndex( endPoint );
    D gDist = _goalTri->getLowerBound( gEndPointIdx );
    std::shared_ptr<ETP::Triangle<P,D>> searchStart;
    std::shared_ptr<ETP::Triangle<P,D>> searchGoal;
    bool reverse;
    if( sDist >= gDist ){
      searchStart = _startTri;
      searchGoal = _goalTri;
      reverse = true;
    }else{
      searchStart = _goalTri;
      searchGoal = _startTri;
      reverse = false;
      D tmpDist = gDist;
      gDist = sDist;
      sDist = tmpDist;
    }
    bool goalReached = false;
    std::vector<std::shared_ptr<ETP::Triangle<P,D>>> openList;
    openList.push_back( searchStart );
    D startToGoalGValue;
    while( goalReached == false ){
      std::shared_ptr<ETP::Triangle<P,D>> tCurrent = openList.front();
      std::shared_ptr<ETP::Triangle<P,D>> tNeighbour = tCurrent->getAdjacentVal( tCurrent->getConnectedNodeIndex( endPoint ) );
      if( tNeighbour == searchGoal ){
        openList.insert( openList.begin(), tNeighbour );
        goalReached = true;
        startToGoalGValue = sDist - gDist;
        break;
      }else if( tNeighbour->getLevel() == 3 ){
        break;
      }else{
        openList.insert( openList.begin(), tNeighbour );
      }
    }
    if( goalReached == false ){
      return std::vector<std::shared_ptr<ETP::Triangle<P,D>>>();
    }else{
      std::shared_ptr<ETP::Triangle<P,D>> otherEndPoint;
      for( int i = 0; i < 3; i++ ){
        std::shared_ptr<ETP::Triangle<P,D>> cNode = searchStart->getConnectedNode( i );
        if( cNode != nullptr && cNode != endPoint ){
          otherEndPoint = cNode;
          break;
        }
      }
      int otherEndPtToEndPtIdx = otherEndPoint->getConnectedNodeIndex( endPoint );
      D baseGValue = otherEndPoint->getLowerBound( otherEndPtToEndPtIdx ) - startToGoalGValue;
      if( baseGValue >= startToGoalGValue ){
        if( reverse == true ){
          std::reverse( openList.begin(), openList.end( ) );
        }
        return openList;
      }
      std::vector<std::shared_ptr<ETP::Triangle<P,D>>> secOpenList;
      //std::shared_ptr<ETP::UsedList<P,D>> usedList = std::make_shared<ETP::UsedList<P,D>>();
      //ETP::UsedList<P,D> usedList = ETP::UsedList<P,D>();
std::shared_ptr<ETP::UsedList<P,D>> usedList = std::make_shared<ETP::UsedList<P,D>>();
      setTriangleSearchNode( otherEndPoint, ETP::SearchNode<P,D>( 0.0, baseGValue, endPoint ), usedList );
      secOpenList.push_back( otherEndPoint );
      ETP::Point<D,D> goalCentroid = endPoint->getCentroid();
      while( secOpenList.empty() == false ){
        std::shared_ptr<ETP::Triangle<P,D>> tCurrent = secOpenList.front();
        secOpenList.erase( secOpenList.begin() );
        ETP::SearchNode<P,D> tsn = *tCurrent->getSearchNode().get();
        std::shared_ptr<ETP::Edge<P,D>> edgeIn;
        if( tCurrent->getLevel() == 3 ){
          edgeIn = tCurrent->getEdgeVal( tCurrent->getConnectedNodeIndex( tsn.comeFrom ) );
        }
        for( int i = 0; i < 3; i++ ){
          std::shared_ptr<ETP::Triangle<P,D>> tConnectedNode = tCurrent->getConnectedNode( i );
          if( tConnectedNode == nullptr || ( tConnectedNode == tsn.comeFrom && ( tCurrent != otherEndPoint || i != otherEndPtToEndPtIdx ) ) ){ continue; }
          D tCNgValue = tCurrent->getLowerBound( i ) + tsn.gValue;
          if( tCurrent->getLevel() == 3 ){
            tCNgValue += tCurrent->getAngle( tCurrent->getEdgeVal( i ), edgeIn );
          }
          std::shared_ptr<ETP::SearchNode<P,D>> tCsN = tConnectedNode->getSearchNode();
          if( tCNgValue >= startToGoalGValue || ( tCsN != nullptr && tCsN->gValue <= tCNgValue ) ){ continue; }
          D tCNhValue = tConnectedNode->getCentroid().dist( goalCentroid );
          if( tConnectedNode == endPoint ){
            if( tCsN == nullptr ){
              setTriangleSearchNode( tConnectedNode, ETP::SearchNode<P,D>( tCNhValue, tCNgValue, tCurrent ), usedList );
            }else if( tCsN->gValue > tCNgValue ){
              tCsN->hValue = tCNhValue;
              tCsN->gValue = tCNgValue;
              tCsN->comeFrom = tCurrent;
            }
          }else{
            //tConnectedNode->setSearchNode( ETP::SearchNode<P,D>( tCNhValue, tCNgValue, tCurrent ), usedList );
            setTriangleSearchNode( tConnectedNode, ETP::SearchNode<P,D>( tCNhValue, tCNgValue, tCurrent ), usedList );
            insertInOpenList( tCNhValue, tConnectedNode, secOpenList );
          }
        }
      }
      if( endPoint->getSearchNode() != nullptr ){
        setTriangleSearchNode( searchGoal, ETP::SearchNode<P,D>( 0.0, 0.0, endPoint ), usedList );
        setTriangleSearchNode( otherEndPoint, ETP::SearchNode<P,D>( 0.0, 0.0, searchStart ), usedList );
        std::vector<std::shared_ptr<ETP::Triangle<P,D>>> path = getFunnelWithIntermediateNodes( searchGoal, searchStart );
        if( searchGoal == _startTri ){
          std::reverse( path.begin(), path.end( ) );
        }
        return path;
      }else{
        if( reverse == true ){
          std::reverse( openList.begin(), openList.end( ) );
        }
        return openList;
      }
    }
  }



  std::vector<std::shared_ptr<ETP::Triangle<P,D>>> getFunnelWithIntermediateNodes( std::shared_ptr<ETP::Triangle<P,D>> _goalTri, std::shared_ptr<ETP::Triangle<P,D>> _startTri ){
    std::vector<std::shared_ptr<ETP::Triangle<P,D>>> ret;
    ret.push_back( _goalTri );
    bool startReached = false;
    while( startReached == false ){
      std::shared_ptr<ETP::Triangle<P,D>> tCurrent = ret.front();
      if( tCurrent == _startTri ){ break; }
      std::shared_ptr<ETP::SearchNode<P,D>> tSN = tCurrent->getSearchNode();
      if( tSN == nullptr ){
        startReached = true;
        break;
      }
      std::shared_ptr<ETP::Triangle<P,D>> tComeFrom = tSN->comeFrom;
      if( tComeFrom == nullptr ){
        startReached = true;
        break;
      }
      int tCurLvl = tCurrent->getLevel();
      int tComFrLvl = tComeFrom->getLevel();
      if( ( ( tCurLvl == 3 || tCurLvl == 2 ) && tComFrLvl == 3 ) || ( tCurLvl == 1 && tComFrLvl == 2 ) ){
        bool comeFromReached = false;
        std::shared_ptr<ETP::Triangle<P,D>> tNode = tCurrent;
        while( comeFromReached == false ){
          std::shared_ptr<ETP::Triangle<P,D>> tNext = tNode->getAdjacentVal( tNode->getConnectedNodeIndex( tComeFrom ) );
          ret.insert( ret.begin(), tNext );
          if( tNext == tComeFrom ){ break; }
          tNode = tNext;
        }
      }else if( tCurLvl == 3 && tComFrLvl == 2 ){
        std::vector<std::shared_ptr<ETP::Triangle<P,D>>> lvl2Part;
        std::shared_ptr<ETP::Triangle<P,D>> tNode = tComeFrom;
        lvl2Part.push_back( tNode );
        bool currentReached = false;
        while( currentReached == false ){
          std::shared_ptr<ETP::Triangle<P,D>> tn = lvl2Part.front();
          std::shared_ptr<ETP::Triangle<P,D>> tNxt = tn->getAdjacentVal( tn->getConnectedNodeIndex( tCurrent ) );
          if( tNxt == tCurrent ){ break; }
          lvl2Part.insert( lvl2Part.begin(), tNxt );
        }
        ret.insert( ret.begin(), lvl2Part.begin(), lvl2Part.end() );
      }else if( tCurLvl == 2 && tComFrLvl == 1 ){
        std::vector<std::shared_ptr<ETP::Triangle<P,D>>> lvl1Part;
        lvl1Part.push_back( tComeFrom );
        bool lvl2Reached = false;
        while( lvl2Reached == false ){
          std::shared_ptr<ETP::Triangle<P,D>> tn = lvl1Part.front();
          std::shared_ptr<ETP::Triangle<P,D>> tNxt = tn->getAdjacentVal( tn->getConnectedNodeIndex( tCurrent ) );
          if( tNxt == tCurrent ){ break; }
          lvl1Part.insert( lvl1Part.begin(), tNxt );
        }
        std::reverse( lvl1Part.begin(), lvl1Part.end() );
        ret.insert( ret.begin(), lvl1Part.begin(), lvl1Part.end() );
      }
    }
    return ret;
  }
  std::vector<std::shared_ptr<ETP::Triangle<P,D>>> getFunnel( std::shared_ptr<ETP::Triangle<P,D>> _goalTri ){
    std::vector<std::shared_ptr<ETP::Triangle<P,D>>> ret;
    ret.push_back( _goalTri );
    bool startReached = false;
    while( startReached == false ){
      std::shared_ptr<ETP::Triangle<P,D>> tCurrent = ret.front();
      std::shared_ptr<ETP::SearchNode<P,D>> tSN = tCurrent->getSearchNode();
      if( tSN == nullptr ){
        startReached = true;
        break;
      }
      std::shared_ptr<ETP::Triangle<P,D>> tComeFrom = tSN->comeFrom;
      if( tComeFrom == nullptr ){
        startReached = true;
        break;
      }
      ret.insert( ret.begin(), tComeFrom );
    }
    return ret;
  }

  std::vector<std::shared_ptr<ETP::Triangle<P,D>>> getFunnelToEndPoint( std::shared_ptr<ETP::searchEndPoint<P,D>> _endPoint ){
    std::vector<std::shared_ptr<ETP::Triangle<P,D>>> ret;
    ret.push_back( _endPoint->triangle );
    std::shared_ptr<ETP::searchEndPoint<P,D>> epCurrent = _endPoint;
    bool ended = false;
    while( ended == false ){
      std::shared_ptr<ETP::searchEndPoint<P,D>> parent = epCurrent->parent;
      if( parent == nullptr ){ break; }
      std::shared_ptr<ETP::Triangle<P,D>> tCurrent = parent->triangle;
      std::shared_ptr<ETP::Triangle<P,D>> endPt = epCurrent->triangle;
      std::vector<std::shared_ptr<ETP::Triangle<P,D>>> path;
      path.push_back( tCurrent );
      bool endPath = false;
      while( endPath == false ){
        std::shared_ptr<ETP::Triangle<P,D>> tNext = tCurrent->getAdjacentVal( tCurrent->getConnectedNodeIndex( endPt ) );
        if( tNext == endPt ){ break; }
        path.push_back( tNext );
        tCurrent = tNext;
      }
      epCurrent = parent;
      ret.insert( ret.begin(), path.begin(), path.end() );
    }
    std::reverse( ret.begin(), ret.end() );
    return ret;
  }

  void setTriangleSearchNode( std::shared_ptr<ETP::Triangle<P,D>> _tri, ETP::SearchNode<P,D> _searchNode, std::shared_ptr<ETP::UsedList<P,D>> _usedList ){
    _tri->setSearchNode( _searchNode );
    _usedList->push_back( _tri );
  }

  std::vector<std::shared_ptr<ETP::searchEndPoint<P,D>>> getTriangleEndPoints( std::shared_ptr<ETP::Triangle<P,D>> _tri, ETP::Point<P,D> _goalPoint ){
  std::vector<std::shared_ptr<ETP::searchEndPoint<P,D>>> ret;
//  ETP::Triangle<P,D> tri = *_tri.get();
  std::shared_ptr<ETP::searchEndPoint<P,D>> finalEndPt = std::make_shared<ETP::searchEndPoint<P,D>>( 0.0, _tri );
  int lvl = _tri->getLevel();
  if( lvl == 3 ){
    ret.push_back( finalEndPt );
  }else if( lvl == 2 ){
    for( int i = 0; i < 3; i++ ){
      std::shared_ptr<ETP::Triangle<P,D>> connectedNode = _tri->getConnectedNode( i );
      if( connectedNode != nullptr ){
        bool isNew = true;
        D lowerBound = getPointToEdgeLowerBound( _goalPoint, _tri->getEdgeVal( i ) );
        size_t retl = ret.size();
        for ( auto ePt : ret ){
          if( connectedNode == ePt->triangle ){
            isNew = false;
            if( lowerBound < ePt->gValue ){
              ePt = std::make_shared<ETP::searchEndPoint<P,D>>( lowerBound, connectedNode, finalEndPt, _tri->getReversedConnectedNodeIndex( connectedNode ) );
            }
          }
        }
        if( isNew == true ){
          ret.push_back( std::make_shared<ETP::searchEndPoint<P,D>>( lowerBound, connectedNode, finalEndPt, _tri->getReversedConnectedNodeIndex( connectedNode ) ) );
        }
      }
    }
  }else if( lvl == 1 ){
    for( int i = 0; i < 3; i++ ){
      std::shared_ptr<ETP::Triangle<P,D>> connectedNode = _tri->getConnectedNode( i );
      D startLowerBound = getPointToEdgeLowerBound( _goalPoint, _tri->getEdgeVal( i ) );
      if( connectedNode != nullptr ){
        std::shared_ptr<ETP::searchEndPoint<P,D>> lvl2EndPt = std::make_shared<ETP::searchEndPoint<P,D>>( startLowerBound + _tri->getLowerBound( i ), connectedNode, finalEndPt );
        for( int j = 0; j < 3; j++ ){
          std::shared_ptr<ETP::Triangle<P,D>> lvl3Node = connectedNode->getConnectedNode( j );
          if( lvl3Node != nullptr ){
            bool isNew = true;
            D lowerBound = connectedNode->getLowerBound( j );
            size_t retl = ret.size();
            for( int k = 0; k < retl; k++ ){
              std::shared_ptr<ETP::searchEndPoint<P,D>> ePt = ret[ k ];
              if( lvl3Node == ePt->triangle ){
                isNew = false;
                if( lowerBound < ePt->gValue ){
                  ret[ k ] = std::make_shared<ETP::searchEndPoint<P,D>>( lvl2EndPt->gValue + lowerBound, lvl3Node, lvl2EndPt );
                }
                break;
              }
            }
            if( isNew == true ){
              ret.push_back( std::make_shared<ETP::searchEndPoint<P,D>>( lvl2EndPt->gValue + lowerBound, lvl3Node, lvl2EndPt, connectedNode->getReversedConnectedNodeIndex( lvl3Node ) ) );
            }
          }
        }
        break;
      }
    }
  }
  return ret;
};

  std::vector<std::shared_ptr<ETP::searchEndPoint<P,D>>> getTriangleEndPoints( std::shared_ptr<ETP::Triangle<P,D>> _tri ){
    std::vector<std::shared_ptr<ETP::searchEndPoint<P,D>>> ret;
  //  ETP::Triangle<P,D> tri = *_tri.get();
    std::shared_ptr<ETP::searchEndPoint<P,D>> finalEndPt = std::make_shared<ETP::searchEndPoint<P,D>>( 0.0, _tri );
    int lvl = _tri->getLevel();
    if( lvl == 3 ){
      ret.push_back( finalEndPt );
    }else if( lvl == 2 ){
      for( int i = 0; i < 3; i++ ){
        std::shared_ptr<ETP::Triangle<P,D>> connectedNode = _tri->getConnectedNode( i );
        if( connectedNode != nullptr ){
          bool isNew = true;
          D lowerBound = _tri->getLowerBound( i );
          size_t retl = ret.size();
          for (int j = 0; j < retl; j++ ){
            std::shared_ptr<ETP::searchEndPoint<P,D>> ePt = ret[ j ];
            if( connectedNode == ePt->triangle ){
              isNew = false;
              if( lowerBound < ePt->gValue ){
                ret[ j ] = std::make_shared<ETP::searchEndPoint<P,D>>( lowerBound, connectedNode, finalEndPt );
              }
            }
          }
          if( isNew == true ){
            ret.push_back( std::make_shared<ETP::searchEndPoint<P,D>>( lowerBound, connectedNode, finalEndPt ) );
          }
        }
      }
    }else if( lvl == 1 ){
      for( int i = 0; i < 3; i++ ){
        std::shared_ptr<ETP::Triangle<P,D>> connectedNode = _tri->getConnectedNode( i );
        if( connectedNode != nullptr ){
          std::shared_ptr<ETP::searchEndPoint<P,D>> lvl2EndPt = std::make_shared<ETP::searchEndPoint<P,D>>( _tri->getLowerBound( i ), connectedNode, finalEndPt );
          for( int j = 0; j < 3; j++ ){
            std::shared_ptr<ETP::Triangle<P,D>> lvl3Node = connectedNode->getConnectedNode( j );
            if( lvl3Node != nullptr ){
              bool isNew = true;
              D lowerBound = connectedNode->getLowerBound( j );
              size_t retl = ret.size();
              for( int k = 0; k < retl; k++ ){
                std::shared_ptr<ETP::searchEndPoint<P,D>> ePt = ret[ k ];
                if( lvl3Node == ePt->triangle ){
                  isNew = false;
                  if( lowerBound < ePt->gValue ){
                    ret[ k ] = std::make_shared<ETP::searchEndPoint<P,D>>( lvl2EndPt->gValue + lowerBound, lvl3Node, lvl2EndPt );
                  }
                  break;
                }
              }
              if( isNew == true ){
                ret.push_back( std::make_shared<ETP::searchEndPoint<P,D>>( lvl2EndPt->gValue + lowerBound, lvl3Node, lvl2EndPt ) );
              }
            }
          }
          break;
        }
      }
    }
    return ret;
  };
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
  D getPointToEdgeLowerBound(  ETP::Point<P,D> vertex, std::shared_ptr<ETP::Edge<P,D>> _edge ){
    ETP::Point<P,D> pt1 = *_edge->getPointVal( 0 ).get();
    ETP::Point<P,D> pt2 = *_edge->getPointVal( 1 ).get();

    /*
    ETP::Point<P,D> ept1;
    ETP::Point<P,D> ept2;
    if( vertex.dist( pt1 ) < vertex.dist( pt2 ) ){
      pt1 = pt1;
      pt2 = pt2;
    }else{
      pt1 = ept2;
      pt2 = ept1;
    }
    */

    P x1 = vertex.x - pt1.x;
    P y1 = vertex.y - pt1.y;
    P x2 = pt2.x - pt1.x;
    P y2 = pt2.y - pt1.y;
    P dot = x1 * x2 + y1 * y2;
    P det = x1 * y2 - y1 * x2;
    D a1Angle = std::abs( std::atan2( (D) det, (D) dot ) );
    if( a1Angle > M_PI ){
      a1Angle =  ( M_PI * 2.0 - a1Angle );
    }
    //return a1Angle;


    x1 = vertex.x - pt2.x;
    y1 = vertex.y - pt2.y;
    x2 = pt1.x - pt2.x;
    y2 = pt1.y - pt2.y;
    dot = x1 * x2 + y1 * y2;
    det = x1 * y2 - y1 * x2;
    D a2Angle = std::abs( std::atan2( (D) det, (D) dot ) );
    if( a2Angle > M_PI ){
      a2Angle =  ( M_PI * 2.0 - a2Angle );
    }
    if( a1Angle < a2Angle ){
      return a1Angle;
    }
    return a2Angle;
  }
  D getPointToEdgeLowerBound(  ETP::Point<P,D> _p, std::shared_ptr<ETP::Edge<P,D>> _edge, D _radius ){
    ETP::Point<P,D> bv = *_edge->getPointVal( 0 ).get();
    ETP::Point<D,D> v = ETP::Point<D,D>( (D) bv.x, (D) bv.y );
    ETP::Point<P,D> bw = *_edge->getPointVal( 1 ).get();
    ETP::Point<D,D> w = ETP::Point<D,D>( (D) bw.x, (D) bw.y );
    ETP::Point<D,D> p = ETP::Point<D,D>( (D) _p.x, (D) _p.y );
    D l2 = ( ( w.x - v.x ) * ( w.x - v.x ) ) + ( ( w.y - v.y ) * ( w.y - v.y ) );
    if ( l2 == 0.0 ) return p.dist( v );
    ETP::Point<P,D> vp = ETP::Point<P,D>( p.x - v.x, p.y - v.y );
    ETP::Point<P,D> vw = ETP::Point<P,D>( w.x - v.x, w.y - v.y );
    D t = vp.dot( vw ) / l2;
    if( t < 0.0 ){ //find close to w
      ETP::Point<D,D> wv = ETP::Point<D,D>( v.x - w.x, v.y - w.y );
      D wvL = std::sqrt( wv.x * wv.x + wv.y * wv.y );
      ETP::Point<D,D> nwv = ETP::Point<D,D>( wv.x / wvL, wv.y / wvL );
      ETP::Point<D,D> pwv = ETP::Point<D,D>( w.x + nwv.x * _radius, w.y + nwv.y * _radius );
      return p.dist( pwv );
    }
    if( t > 1.0 ){ // find close to v
      ETP::Point<D,D> vw = ETP::Point<D,D>( w.x - v.x, w.y - v.y );
      D vwL = std::sqrt( vw.x * vw.x + vw.y * vw.y );
      ETP::Point<D,D> nvw = ETP::Point<D,D>( vw.x / vwL, vw.y / vwL );
      ETP::Point<D,D> pvw = ETP::Point<D,D>( v.x + nvw.x * _radius, v.y + nvw.y * _radius );
      return p.dist( pvw );
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
  D getPointToEdgeLowerBound(  ETP::Point<P,D> _p, std::shared_ptr<ETP::Edge<P,D>> _edge, D _radius, D _scale ){
    return _edge->minimumDistance( _p, _scale ) / _scale;


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
