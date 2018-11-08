#include "./triangulationAbstraction.h"

void ETP::collapseUnrootedTree( std::shared_ptr<ETP::Triangle> t, std::shared_ptr<ETP::Component> c ){
  std::stack<std::shared_ptr<ETP::Triangle>> s;//Stack
  s.push( t );
  while( !s.empty() ){
    std::shared_ptr<ETP::Triangle> tCurrent = s.top();
    s.pop();
    tCurrent->setComponent( c );
    for( size_t i = 0; i < 3; i++ ){
      ETP::Edge e = ETP::getEdge( tCurrent, i );
      tCurrent->setAdjacent( i, nullptr );
      if( e->constrained() ){
        tCurrent->setAngle( i, 0 );
        tCurrent->setChoke( i, 0 );
      }else{
        tCurrent->setAngle( i, std::numeric_limits<float>::infinity() );
        tCurrent->setChoke( i, std::numeric_limits<float>::infinity() );
        std::shared_ptr<ETP::Triangle> tNext = getTriangleAccross( tCurrent, e );
        if( tNext->getComponent() == nullptr ){
          s.push( tNext );
        }
      }
    }
  }
}
