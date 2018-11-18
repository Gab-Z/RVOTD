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


#include <nan.h>

#include "../clipper/clipper.hpp"
#include "../poly2tri/poly2tri.h"

#include "./Point.h"
#include "./Edge.h"
#include "./Triangle.h"
#include "./Component.h"

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


public:

  TriangulationSpace(){};
/*
  void buildSpace( std::vector<std::vector<ETP::Point<P,D>>> _polygon ){
    std::vector<N> indices = mapbox::earcut<N>( _polygon );

    for( size_t i = 0; i < _polygon.size(); i++ ){
      std::vector<ETP::Point<P,D>> polyObj = _polygon[ i ];
      size_t polyObjEnd = polyObj[ 0 ].equals( polyObj[ polyObj.size() - 1 ] ) ? polyObj.size() - 1 : polyObj.size();
      //size_t polyObjEnd = polyObj.size();

      for( size_t ii = 0; ii < polyObjEnd; ii++ ){
        size_t iNext = (ii < (polyObjEnd - 1)) ?( ii + 1) : 0;
          edges.push_back( std::move( std::make_shared<ETP::Edge<P,D>>( polyObj[ ii ], polyObj[ iNext ], true ) ) );
      }

    }

//    for( size_t j = 0; j < _polygon.size(); j++ ){
//      std::vector<ETP::Point<P,D>>& polyObj = _polygon[ j ];
//      size_t polyObjEnd = polyObj[ 0 ].equals( polyObj[ (int) polyObj.size() - 1 ] ) ? (int) polyObj.size() - 1 : polyObj.size();
//      for( size_t jj = 0; jj < polyObjEnd; jj++ ){
//        points.push_back( std::make_shared<ETP::Point<P,D>>( polyObj[ jj ] ) );
//      }
//    }


    for( size_t el = 0; el < edges.size(); el++ ){
      bool point1IsNew = true;
      bool point2IsNew = true;
      for( size_t pi = 0; pi < points.size(); pi++ ){
        std::shared_ptr<ETP::Point<P,D>>& testPoint = points[ pi ];
        if( testPoint->equals( edges[ el ]->getPointRef( 0 ) ) ){
          point1IsNew = false;
        }
        if( testPoint->equals( edges[ el ]->getPointRef( 1 ) ) ){
          point2IsNew = false;
        }
      }
      if( point1IsNew ){
        points.push_back( edges[ el ]->getPointVal( 0 ) );
      }
      if( point2IsNew ){
        points.push_back( edges[ el ]->getPointVal( 1 ) );
      }
    }
    for( size_t n = 0; n < indices.size() / 3; n++ ){
      std::shared_ptr<ETP::Point<P,D>>& pt1 = points[ indices[ n * 3 ] ];
      std::shared_ptr<ETP::Point<P,D>>& pt2 = points[ indices[ n * 3 + 1 ] ];
      std::shared_ptr<ETP::Point<P,D>>& pt3 = points[ indices[ n * 3 + 2 ] ];
      bool newEdge1 = true;
      bool newEdge2 = true;
      bool newEdge3 = true;
      for( size_t e = 0; e < edges.size(); e++ ){
        if( edges[ e ]->hasSamePoints( points[ indices[ n * 3 ] ], points[ indices[ n * 3 + 1 ] ] ) ){
          newEdge1 = false;
        }
        if( edges[ e ]->hasSamePoints( points[ indices[ n * 3 + 1 ] ], points[ indices[ n * 3 + 2 ] ] ) ){
          newEdge2 = false;
        }
        if( edges[ e ]->hasSamePoints( points[ indices[ n * 3 + 2 ] ], points[ indices[ n * 3 ] ] ) ){
          newEdge3 = false;
        }
      }
      if( newEdge1 ){
        edges.push_back( std::make_shared<Edge<P,D>>( points[ indices[ n * 3 ] ], points[ indices[ n * 3 + 1 ] ], false  ) );
      }
      if( newEdge2 ){
        edges.push_back( std::make_shared<Edge<P,D>>( points[ indices[ n * 3 + 1 ] ], points[ indices[ n * 3 + 2 ] ], false  ) );
      }
      if( newEdge3 ){
        edges.push_back( std::make_shared<Edge<P,D>>( points[ indices[ n * 3 + 2 ] ], points[ indices[ n * 3 ] ], false  ) );
      }

    }

  }
*/

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


    //info.GetReturnValue().Set( self->space->getEdges() );
    //info.GetReturnValue().Set( self->space->getTriangles() );

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
    points.push_back( std::move( _pt ) );
    return points[ points.size() - 1 ];
  }
  std::shared_ptr<ETP::Point<P,D>> addGetPoint( const P& _x, const P& _y ){
    for( size_t i = 0; i < points.size(); i++ ){
      if( points[ i ]->equals( _x, _y ) ){
        return points[ i ];
      }
    }
    points.push_back( std::move( std::make_shared<ETP::Point<P,D>>( _x, _y ) ) );
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
/*
  static D getAngle( std::shared_ptr<ETP::Edge<P,D>> _edge1, std::shared_ptr<ETP::Edge<P,D>> _edge2 ){
    if( _edge1 == nullptr || _edge2 == nullptr ){ return 0.0f; }
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
    }
    D Pv1 = std::sqrt( (D) ( std::pow( vertex->x - pt1->x , 2 ) + std::pow( vertex->y - pt1->y, 2 ) ) );
    D Pv2 = std::sqrt( (D) ( std::pow( vertex->x - pt2->x , 2 ) + std::pow( vertex->y - pt2->y, 2 ) ) );
    D P12 = std::sqrt( (D) ( std::pow( pt1->x - pt2->x , 2 ) + std::pow( pt1->y - pt2->y, 2 ) ) );
    return std::acos( ( std::pow( Pv1, 2 ) + std::pow( Pv2, 2 ) - std::pow( P12, 2 ) ) / ( 2.0 * Pv1 * Pv2 ) );
  }
  */
  static D getAngle( std::shared_ptr<ETP::Edge<P,D>> _edge1, std::shared_ptr<ETP::Edge<P,D>> _edge2 ){
    if( _edge1 == nullptr || _edge2 == nullptr || _edge1->hasSamePoints( _edge2 )){ return (D) 0.0f; }
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
          tCurrent->setLowerBound( i, currAngle );
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
    int ct = 0;
    while( q.empty() == false ){
      ct++;
      std::shared_ptr<ETP::Triangle<P,D>> t = q.front();
      q.erase( q.begin(), q.begin() + 1 );
      t->setLevel( 3 );
      t->setComponent( c );
      for( size_t j = 0; j < 3; j++ ){
        int previousTid = t->getId();
        std::shared_ptr<ETP::Triangle<P,D>> t1stNeighbour = t->getAdjacentVal( j );
        D chokeCount = t->getEdgeVal( j )->getLength();
        D lowerBoundCount = 0.0f;
        D memChoke = 0.0f;
        D memLowerBound = 0.0f;
        if( t1stNeighbour->getNumConstrainedEdges() + t1stNeighbour->getNumAdjacentLevel( 1 ) == 0 ){
          if( t1stNeighbour->getLevel() == -1 ){
            q.push_back( t1stNeighbour );
          }
          t1stNeighbour->setConnectedNode( t1stNeighbour->getAdjacentIndex( t ), t );
          continue;
        }else{
          t1stNeighbour->setLevel( 2 );
          t1stNeighbour->setComponent( c );
          t1stNeighbour->setConnectedNode( t1stNeighbour->getAdjacentIndex( t ), t );
        }
        std::shared_ptr<ETP::Triangle<P,D>> tCurrent = t1stNeighbour;
        std::shared_ptr<ETP::Edge<P,D>> edgeIn = tCurrent->getEdgeVal( j );
        bool reachedLvl3 = false;
        while( reachedLvl3 == false ){
          for( size_t k = 0; k < 3; k++ ){
            std::shared_ptr<ETP::Edge<P,D>> edgeOut = tCurrent->getEdgeVal( k );
            if( edgeOut->isConstrained() == true ){ continue; }
            std::shared_ptr<ETP::Triangle<P,D>> tCurrentNeighbour = tCurrent->getAdjacentVal( k );
            if( tCurrentNeighbour->getId() == previousTid ){ continue; }
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
              if( neighbN + neighbM == 0 ){
                if( tCurrentNeighbour->getLevel() == -1 ){
                  q.push_back( tCurrentNeighbour );
                }
                tCurrentNeighbour->setConnectedNode( tCurrentNeighbour->getAdjacentIndex( tCurrent ), t );
                reachedLvl3 = true;
                break;
              }
              tCurrentNeighbour->setLevel( 2 );
              tCurrentNeighbour->setComponent( c );
              tCurrentNeighbour->setConnectedNode( tCurrentNeighbour->getAdjacentIndex( tCurrent ), t );
              memChoke = tCurrent->getWidthbetweenEdges( edgeIn, edgeOut );
              memLowerBound = getAngle( edgeIn, edgeOut );
              previousTid = tCurrent->getId();
              tCurrent = tCurrentNeighbour;
              edgeIn = edgeOut;
              chokeCount = std::min( chokeCount, memChoke );
              lowerBoundCount += memLowerBound;
            }
          }
        }
      }
    }
  }


};// end of class
}// end of namespace




#endif
