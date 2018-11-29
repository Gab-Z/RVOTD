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
      std::shared_ptr<ETP::Triangle<P,D>>& tri = triangles[ ti ];
      v8::Local<v8::Object> retTri = Nan::New<v8::Object>();

      v8::Local<v8::String> idProp = Nan::New( "id" ).ToLocalChecked();
      v8::Local<v8::Value> idValue = Nan::New( tri->getId() );
      retTri->Set( idProp, idValue );

      std::shared_ptr<ETP::Edge<P,D>>& e1 = tri->getEdgeRef( 0 );
      std::shared_ptr<ETP::Edge<P,D>>& e2 = tri->getEdgeRef( 1 );
      std::shared_ptr<ETP::Edge<P,D>>& e3 = tri->getEdgeRef( 2 );

      std::vector<P> positions = tri->getPositions();

      v8::Local<v8::Array> posArr = Nan::New<v8::Array>( 6 );
        posArr->Set( 0, Nan::New( (NanInt) positions[ 0 ] ) );
        posArr->Set( 1, Nan::New( (NanInt) positions[ 1 ] ) );
        posArr->Set( 2, Nan::New( (NanInt) positions[ 2 ] ) );
        posArr->Set( 3, Nan::New( (NanInt) positions[ 3 ] ) );
        posArr->Set( 4, Nan::New( (NanInt) positions[ 4 ] ) );
        posArr->Set( 5, Nan::New( (NanInt) positions[ 5 ] ) );
      v8::Local<v8::String> posProp = Nan::New( "positions" ).ToLocalChecked();
      retTri->Set( posProp, posArr );

      v8::Local<v8::Array> constrArr = Nan::New<v8::Array>( 3 );
        constrArr->Set( 0, Nan::New( e1->isConstrained() ) );
        constrArr->Set( 1, Nan::New( e2->isConstrained() ) );
        constrArr->Set( 2, Nan::New( e3->isConstrained() ) );
      v8::Local<v8::String> constrProp = Nan::New( "constrained" ).ToLocalChecked();
      retTri->Set( constrProp, constrArr );

      v8::Local<v8::Array> adjacentArr = Nan::New<v8::Array>( 3 );
      //std::shared_ptr<ETP::Triangle<P,D>>& a1 = tri->getAdjacentRef( 0 );
      //std::shared_ptr<ETP::Triangle<P,D>>& a2 = tri->getAdjacentRef( 1 );
      //std::shared_ptr<ETP::Triangle<P,D>>& a3 = tri->getAdjacentRef( 2 );
      std::shared_ptr<ETP::Triangle<P,D>> a1 = tri->triangleOpposite( tri->getEdgeVal( 0 ) );
      std::shared_ptr<ETP::Triangle<P,D>> a2 = tri->triangleOpposite( tri->getEdgeVal( 1 ) );
      std::shared_ptr<ETP::Triangle<P,D>> a3 = tri->triangleOpposite( tri->getEdgeVal( 2 ) );
        adjacentArr->Set( 0, Nan::New( a1 != nullptr ? (int) a1->getId() : -1 ) );
        adjacentArr->Set( 1, Nan::New( a2 != nullptr ? (int) a2->getId() : -1 ) );
        adjacentArr->Set( 2, Nan::New( a3 != nullptr ? (int) a3->getId() : -1 ) );
      v8::Local<v8::String> adjacentProp = Nan::New( "adjacents" ).ToLocalChecked();
      retTri->Set( adjacentProp, adjacentArr );

      v8::Local<v8::Array> widthsArr = Nan::New<v8::Array>( 3 );
      widthsArr->Set( 0, Nan::New( (NanFloat) tri->getWidth( 0 ) ) );
      widthsArr->Set( 1, Nan::New( (NanFloat) tri->getWidth( 1 ) ) );
      widthsArr->Set( 2, Nan::New( (NanFloat) tri->getWidth( 2 ) ) );
      v8::Local<v8::String> widthsProp = Nan::New( "widths" ).ToLocalChecked();
      retTri->Set( widthsProp, widthsArr );

      v8::Local<v8::Array> anglesArr = Nan::New<v8::Array>( 3 );
      anglesArr->Set( 0, Nan::New( (NanFloat) tri->getAngle( 0 ) ) );
      anglesArr->Set( 1, Nan::New( (NanFloat) tri->getAngle( 1 ) ) );
      anglesArr->Set( 2, Nan::New( (NanFloat) tri->getAngle( 2 ) ) );
      v8::Local<v8::String> anglesProp = Nan::New( "angles" ).ToLocalChecked();
      retTri->Set( anglesProp, anglesArr );

      v8::Local<v8::Array> lowerBoundsArr = Nan::New<v8::Array>( 3 );
      lowerBoundsArr->Set( 0, Nan::New( (NanFloat) tri->getLowerBound( 0 ) ) );
      lowerBoundsArr->Set( 1, Nan::New( (NanFloat) tri->getLowerBound( 1 ) ) );
      lowerBoundsArr->Set( 2, Nan::New( (NanFloat) tri->getLowerBound( 2 ) ) );
      v8::Local<v8::String> lowerBoundsProp = Nan::New( "lowerBounds" ).ToLocalChecked();
      retTri->Set( lowerBoundsProp, lowerBoundsArr );

      v8::Local<v8::String> lvlProp = Nan::New( "level" ).ToLocalChecked();
      v8::Local<v8::Value> lvlVal = Nan::New( tri->getLevel() );
      retTri->Set( lvlProp, lvlVal );

      v8::Local<v8::String> componentProp = Nan::New( "component" ).ToLocalChecked();
      v8::Local<v8::Value> componentVal = Nan::New( tri->getComponent() );
      retTri->Set( componentProp, componentVal );

      v8::Local<v8::Array> nodesArr = Nan::New<v8::Array>( 3 );
      std::shared_ptr<ETP::Triangle<P,D>> n1 = tri->getConnectedNode( 0 );
      std::shared_ptr<ETP::Triangle<P,D>> n2 = tri->getConnectedNode( 1 );
      std::shared_ptr<ETP::Triangle<P,D>> n3 = tri->getConnectedNode( 2 );
        nodesArr->Set( 0, Nan::New( n1 != nullptr ? (int) n1->getId() : -1 ) );
        nodesArr->Set( 1, Nan::New( n2 != nullptr ? (int) n2->getId() : -1 ) );
        nodesArr->Set( 2, Nan::New( n3 != nullptr ? (int) n3->getId() : -1 ) );
      v8::Local<v8::String> nodesProp = Nan::New( "nodes" ).ToLocalChecked();
      retTri->Set( nodesProp, nodesArr );

      v8::Local<v8::String> centroidProp = Nan::New( "centroid" ).ToLocalChecked();
      v8::Local<v8::Object> centroidOb = Nan::New<v8::Object>();
      v8::Local<v8::String> xProp = Nan::New( "x" ).ToLocalChecked();
      v8::Local<v8::String> yProp = Nan::New( "y" ).ToLocalChecked();
      ETP::Point<D,D> centroid = tri->getCentroid();
      centroidOb->Set( xProp, Nan::New( centroid.x ) );
      centroidOb->Set( yProp, Nan::New( centroid.y ) );
      retTri->Set( centroidProp, centroidOb );

      ret->Set( ti, retTri );
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
            std::shared_ptr<ETP::Triangle<P,D>> tTemp = tCurrent->triangleOpposite( e );
            if( tTemp->getLevel() == 1 ){
              collapseRootedTree( tCurrent, tTemp );
              tCurrent->setAngle( i, 0.0f );
              tCurrent->setChoke( i, 0.0f );
              tCurrent->setAdjacent( i, nullptr );
            }else{
              if( tTemp->getLevel() == -1 ){
                tNext = tTemp;
              }
              tCurrent->setAngle( i, std::numeric_limits<D>::infinity() );
              tCurrent->setChoke( i, std::numeric_limits<D>::infinity() );
              //tCurrent->setLowerBound( i, 55.0 );
              tCurrent->setAdjacent( i, nullptr );
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
              try{
                collapseRootedTree( tCurrent, tCurrentNeighbour );
              }catch( const std::invalid_argument& err ){
                throw;
              }
              continue;
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
          bodySearch = searchInsideLvl2Corridor( _startTri, _goalTri );
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
              std::make_move_iterator( startPart.end() )
           );
          }
          ret.insert(
            ret.end(),
            std::make_move_iterator( bodySearch.begin() ),
            std::make_move_iterator( bodySearch.end() )
          );
          if( _goalTri != goalRoot ){// goal is on a lvl 1 tree, we add triangles from the lvl 2 root to the goal
            std::vector<std::shared_ptr<ETP::Triangle<P,D>>> goalPart = getTrianglesFromLvl1ToConnectedLvl2( _goalTri, _goalTri->getConnectedNodeIndex( goalRoot ) );
            ret.insert(
              ret.end(),
              std::make_move_iterator( goalPart.crend() ),
              std::make_move_iterator( goalPart.crbegin() )
           );
          }
          return ret;
        }

      }
    }
    // perform classic Triangulation Reduction A*

    //tmp
    return std::vector<std::shared_ptr<ETP::Triangle<P,D>>>();

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
    std::vector<std::shared_ptr<ETP::SearchNode<P,D>>> usedList;
    ETP::Point<D,D> goalCentroid = _goalTri->getCentroid();
    _startTri->setSearchNode( ETP::SearchNode<P,D>( 0.0, 0.0 ), usedList );
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
          tNeighbour->setSearchNode( ETP::SearchNode<P,D>( 0.0, snCurrent.gValue + currentCentroid.dist( tNeighbCentroid ), tCurrent ), usedList );
          goalReached = true;
          break;
        }
        std::shared_ptr<ETP::SearchNode<P,D>> snNeighbourPtr = tNeighbour->getSearchNode();
        bool neighbIsNew = false;
        if( snNeighbourPtr == nullptr ){
          D tnHVal = tNeighbCentroid.dist( goalCentroid );
          tNeighbour->setSearchNode( ETP::SearchNode<P,D>( tnHVal, snCurrent.gValue + currentCentroid.dist( tNeighbCentroid ), tCurrent ), usedList );
          insertInOpenList( tnHVal, tNeighbour, openList );
        }
      }
    }
    if( goalReached == false ){
      resetUsedList( std::move( usedList ) );
      return std::vector<std::shared_ptr<ETP::Triangle<P,D>>>();
    }else{
      return getFunnel( _goalTri, std::move( usedList ) );
    }
  }

  std::vector<std::shared_ptr<ETP::Triangle<P,D>>> searchInsideLvl2RingOrLoop( std::shared_ptr<ETP::Triangle<P,D>> _startTri, std::shared_ptr<ETP::Triangle<P,D>> _goalTri ){
    std::vector<std::shared_ptr<ETP::Triangle<P,D>>> openList;
    std::vector<std::shared_ptr<ETP::SearchNode<P,D>>> usedList;
    ETP::Point<D,D> goalCentroid = _goalTri->getCentroid();
    _startTri->setSearchNode( ETP::SearchNode<P,D>( 0.0, 0.0 ), usedList );
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
          return getFunnel( _goalTri, std::move( usedList ) );
        }else if( tNeighbour == _goalTri ){
          if( goalReached == false || neighbGValue < memGValue ){
            memGValue = neighbGValue;
            tNeighbour->setSearchNode( ETP::SearchNode<P,D>( 0.0, memGValue, tCurrent ), usedList );
            if( goalReached == true ){
              return getFunnel( _goalTri, std::move( usedList ) );
            }
            goalReached = true;
          }
        }else{
          D neighbHValue = tNeighbCentroid.dist( goalCentroid );
          tNeighbour->setSearchNode( ETP::SearchNode<P,D>( neighbHValue, neighbGValue, tCurrent ), usedList );
          insertInOpenList( neighbHValue, tNeighbour, openList );
        }
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
      reverse = false;
    }else{
      searchStart = _goalTri;
      searchGoal = _startTri;
      reverse = true;
    }
    bool endReached = false;
    std::vector<std::shared_ptr<ETP::Triangle<P,D>>> openList;
    openList.push_back( searchStart );
    while( endReached == false ){
      std::shared_ptr<ETP::Triangle<P,D>> tCurrent = openList.front();
      std::shared_ptr<ETP::Triangle<P,D>> tNeighbour = tCurrent->getAdjacentVal( tCurrent->getConnectedNodeIndex( endPoint ) );
      if( tNeighbour == searchGoal ){
        openList.insert( openList.begin(), tNeighbour );
        endReached = true;
        break;
      }else if( tNeighbour->getLevel() == 3 ){
        break;
      }else{
        openList.insert( openList.begin(), tNeighbour );
      }
    }
    if( endReached == false ){
      return std::vector<std::shared_ptr<ETP::Triangle<P,D>>>();
    }else{
      if( reverse == true ){
        std::reverse( openList.begin(), openList.end( ) );
      }
      return openList;
    }
  }

  std::vector<std::shared_ptr<ETP::Triangle<P,D>>> getFunnel( std::shared_ptr<ETP::Triangle<P,D>> _goalTri, std::vector<std::shared_ptr<ETP::SearchNode<P,D>>> _usedList ){
    std::vector<std::shared_ptr<ETP::Triangle<P,D>>> ret;
    ret.push_back( _goalTri );
    bool startReached = false;
    while( startReached == false ){
      std::shared_ptr<ETP::Triangle<P,D>> tCurrent = ret.front();
      std::shared_ptr<ETP::Triangle<P,D>> tComeFrom = tCurrent->getSearchNode()->comeFrom;
      if( tComeFrom == nullptr ){
        startReached = true;
        break;
      }
      ret.insert( ret.begin(), tComeFrom );
    }
    resetUsedList( std::move( _usedList ) );
    return ret;
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
  void resetUsedList( std::vector<std::shared_ptr<ETP::SearchNode<P,D>>> _usedList ){
    size_t l = _usedList.size();
    for( size_t i = 0; i < l; i++ ){
      _usedList[ i ] = nullptr;
    }
  }
};// end of class
}// end of namespace




#endif
