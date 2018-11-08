#ifndef DEF_TowerDefense
#define DEF_TowerDefense

#include <iostream>
#include <vector>
#include <nan.h>
#include <memory>
#include <limits>
#include <cmath>
#include <stdexcept>

#include "./clip2tri-master/clipper/clipper.hpp"
#include "./clip2tri-master/poly2tri/poly2tri.h"
#include "./clip2tri-master/clip2tri/clip2tri.h"

#include "./RVO2/RVOSimulator.h"
#include "./RVO2/Vector2.h"
#include "./ETP/Point.h"
#include "./ETP/Edge.h"
#include "./ETP/Triangle.h"
#include "./ETP/TriangulationSpace.h"

#include "./utils.h"


class TowerDefense  : public Nan::ObjectWrap {

  std::unique_ptr<RVO::RVOSimulator> sim;
  std::unique_ptr<ETP::TriangulationSpace<float>> space;


  public:

    static NAN_MODULE_INIT( Init );
    static NAN_METHOD( New );
    static Nan::Persistent<v8::FunctionTemplate> constructor;
    static NAN_METHOD( get );

    static NAN_METHOD( setTimeStep );
    static NAN_METHOD( setAgentDefaults );
    static NAN_METHOD( addAgent );
    static NAN_METHOD( addObstacle );
    static NAN_METHOD( processObstacles );
    static NAN_METHOD( getGlobalTime );
    static NAN_METHOD( getNumAgents );
    static NAN_METHOD( getAgentPosition );
    static NAN_METHOD( getAgentRadius );
    static NAN_METHOD( setAgentPrefVelocity );
    static NAN_METHOD( doStep );

    //static NAN_METHOD( triangulate );
    //static NAN_METHOD( buildTriangulation );

    static NAN_METHOD( clip2triangulate );

    static NAN_METHOD( tryClipper );
    static NAN_METHOD( tryPoly2Tri );


    void init();



};

#endif
