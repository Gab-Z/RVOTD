#include "./core.h"

Nan::Persistent<v8::FunctionTemplate> TowerDefense::constructor;

NAN_MODULE_INIT( TowerDefense::Init ) {
  v8::Local<v8::FunctionTemplate> ctor = Nan::New<v8::FunctionTemplate>( TowerDefense::New );
  constructor.Reset( ctor );
  ctor->InstanceTemplate()->SetInternalFieldCount( 1 );
  ctor->SetClassName( Nan::New( "TowerDefense" ).ToLocalChecked() );

  Nan::SetPrototypeMethod( ctor, "get", get );

  Nan::SetPrototypeMethod( ctor, "setTimeStep", setTimeStep );
  Nan::SetPrototypeMethod( ctor, "setAgentDefaults", setAgentDefaults );
  Nan::SetPrototypeMethod( ctor, "addAgent", addAgent );
  Nan::SetPrototypeMethod( ctor, "addObstacle", addObstacle );
  Nan::SetPrototypeMethod( ctor, "processObstacles", processObstacles );
  Nan::SetPrototypeMethod( ctor, "getGlobalTime", getGlobalTime );
  Nan::SetPrototypeMethod( ctor, "getNumAgents", getNumAgents );
  Nan::SetPrototypeMethod( ctor, "getAgentPosition", getAgentPosition );
  Nan::SetPrototypeMethod( ctor, "getAgentRadius", getAgentRadius );
  Nan::SetPrototypeMethod( ctor, "setAgentPrefVelocity", setAgentPrefVelocity );
  Nan::SetPrototypeMethod( ctor, "doStep", doStep );
  Nan::SetPrototypeMethod( ctor, "clip2triangulate", clip2triangulate );

  Nan::SetPrototypeMethod( ctor, "tryClipper", tryClipper );
  Nan::SetPrototypeMethod( ctor, "tryPoly2Tri", tryPoly2Tri );

  Nan::SetPrototypeMethod( ctor, "testAngle", testAngle );
  Nan::SetPrototypeMethod( ctor, "getTriangleId", getTriangleId );

  Nan::SetPrototypeMethod( ctor, "getSectors", getSectors );
  Nan::SetPrototypeMethod( ctor, "testTRAStar", testTRAStar );

  Nan::SetPrototypeMethod( ctor, "testTRAStarScale", testTRAStarScale );


  //Nan::SetPrototypeMethod( ctor, "triangulate", triangulate );
  //Nan::SetPrototypeMethod( ctor, "buildTriangulation", buildTriangulation );

  target->Set( Nan::New( "TowerDefense" ).ToLocalChecked(), ctor->GetFunction());

}

NAN_METHOD( TowerDefense::New ) {
  // throw an error if constructor is called without new keyword
  if( ! info.IsConstructCall() ) {
    return Nan::ThrowError( Nan::New( "TowerDefense::New - called without new keyword" ).ToLocalChecked() );
  }
  TowerDefense* towerDef = new TowerDefense();
  towerDef->init();
  towerDef->Wrap( info.Holder() );
  info.GetReturnValue().Set( info.Holder() );
}

NAN_METHOD( TowerDefense::get ){
  TowerDefense* self = Nan::ObjectWrap::Unwrap<TowerDefense>(info.This());
  v8::Local<v8::Object> ret = Nan::New<v8::Object>();
  //v8::Local<v8::String> dataProp = Nan::New( "prop" ).ToLocalChecked();
//  v8::Local<v8::Value> dataValue = Nan::New( "test is working !" ).ToLocalChecked();
  //ret->Set( dataProp, dataValue );

  std::vector<int> first;
  for (int i = 0; i < 10; i++ ){
    first.push_back( i );
  }
  std::vector<int> second;
  while( first.empty() == false ){
    int in = first.front();
    first.erase( first.begin() );
    second.push_back( in );
  }
  size_t l = second.size();
  v8::Local<v8::Array> arr = Nan::New<v8::Array>();
  for( size_t i = 0; i < l; i++ ){
    arr->Set( i, Nan::New( second[ i ] ) );
  }
  v8::Local<v8::String> arrProp = Nan::New( "array" ).ToLocalChecked();
  ret->Set( arrProp, arr );
  info.GetReturnValue().Set( ret );
}

void TowerDefense::init(){
  sim = std::make_unique<RVO::RVOSimulator>();
  space = std::make_unique<ETP::TriangulationSpace<tPos,tDecim>>();
}

NAN_METHOD( TowerDefense::setTimeStep ){
  if( info.Length() != 1 ) {
    return Nan::ThrowError( Nan::New( "TowerDefense::setTimeStep - expected 1 argument: number timeStep" ).ToLocalChecked() );
  }
  if( !info[ 0 ]->IsNumber() ){
    return Nan::ThrowError( Nan::New( "TowerDefense::setTimeStep - expected argument [ 0 ] to be a number" ).ToLocalChecked() );
  }
  TowerDefense* self = Nan::ObjectWrap::Unwrap<TowerDefense>( info.This() );
  self->sim->setTimeStep( ( float ) info[ 0 ]->NumberValue() );
}

NAN_METHOD( TowerDefense::setAgentDefaults ){
  TowerDefense* self = Nan::ObjectWrap::Unwrap<TowerDefense>( info.This() );
  if( info.Length() < 6 || info.Length() > 7 ) {
    return Nan::ThrowError( Nan::New( "TowerDefense::setAgentDefaults - expected 6 arguments + 1 optional" ).ToLocalChecked() );
  }
  if( !info[ 0 ]->IsNumber() ){ return Nan::ThrowError( Nan::New( "TowerDefense::setAgentDefaults - expected argument [ 0 ] neighborDist to be a number" ).ToLocalChecked() ); }
  if( !info[ 1 ]->IsNumber() ){ return Nan::ThrowError( Nan::New( "TowerDefense::setAgentDefaults - expected argument [ 1 ] maxNeighbors to be an integer" ).ToLocalChecked() ); }
  if( !info[ 2 ]->IsNumber() ){ return Nan::ThrowError( Nan::New( "TowerDefense::setAgentDefaults - expected argument [ 2 ] timeHorizon to be a number" ).ToLocalChecked() ); }
  if( !info[ 3 ]->IsNumber() ){ return Nan::ThrowError( Nan::New( "TowerDefense::setAgentDefaults - expected argument [ 3 ] timeHorizonObst to be a number" ).ToLocalChecked() ); }
  if( !info[ 4 ]->IsNumber() ){ return Nan::ThrowError( Nan::New( "TowerDefense::setAgentDefaults - expected argument [ 4 ] radius to be a number" ).ToLocalChecked() ); }
  if( !info[ 5 ]->IsNumber() ){ return Nan::ThrowError( Nan::New( "TowerDefense::setAgentDefaults - expected argument [ 5 ] maxSpeed to be a number" ).ToLocalChecked() ); }
  if( info.Length() == 7 && !info[ 6 ]->IsArray() ){ return Nan::ThrowError( Nan::New( "TowerDefense::setAgentDefaults - expected argument [ 6 ] velocity to be an array" ).ToLocalChecked() ); }

  if( info.Length() == 6 ){
    self->sim->setAgentDefaults( ( float ) info[ 0 ]->NumberValue(),//neighborDist
                                ( size_t ) info[ 1 ]->IntegerValue(),//maxNeighbors
                                ( float ) info[ 2 ]->NumberValue(),//timeHorizon
                                ( float ) info[ 3 ]->NumberValue(),//timeHorizonObst
                                ( float ) info[ 4 ]->NumberValue(),//radius
                                ( float ) info[ 5 ]->NumberValue()//maxSpeed
   );
 }else if( info.Length() == 7 ){
    v8::Local<v8::Array> arrVelocity = v8::Local<v8::Array>::Cast( info[ 6 ] );
    if( arrVelocity->Length() != 2 ){
      return Nan::ThrowError( Nan::New( "TowerDefense::setAgentDefaults - expected argument [ 6 ] to be an array of length 2" ).ToLocalChecked() );
    }
    if( !arrVelocity->Get( 0 )->IsNumber() || !arrVelocity->Get( 1 )->IsNumber() ){
      return Nan::ThrowError( Nan::New( "TowerDefense::setAgentDefaults - expected argument [ 6 ] to be an array of 2 numbers x & y" ).ToLocalChecked() );
    }
    self->sim->setAgentDefaults( ( float ) info[ 0 ]->NumberValue(),//neighborDist
                                ( size_t ) info[ 1 ]->IntegerValue(),//maxNeighbors
                                ( float ) info[ 2 ]->NumberValue(),//timeHorizon
                                ( float ) info[ 3 ]->NumberValue(),//timeHorizonObst
                                ( float ) info[ 4 ]->NumberValue(),//radius
                                ( float ) info[ 5 ]->NumberValue(),//maxSpeed
                                RVO::Vector2(   ( float ) arrVelocity->Get( 0 )->NumberValue(),//velocity x
                                                ( float ) arrVelocity->Get( 1 )->NumberValue()//velocity y
                                )
   );
 }
}

NAN_METHOD( TowerDefense::addAgent ){
  if( info.Length() < 2 ){
    return Nan::ThrowError( Nan::New( "TowerDefense::addAgent - expected at least 2 arguments numbers x & y" ).ToLocalChecked() );
  }
  if( info.Length() > 2 && ( info.Length() != 8 && info.Length() != 9 ) ){
    return Nan::ThrowError( Nan::New( "TowerDefense::addAgent - expected 9 arguments if setting defaults" ).ToLocalChecked() );
  }
  if( !info[ 0 ]->IsNumber() ){ return Nan::ThrowError( Nan::New( "TowerDefense::addAgent - expected argument [ 0 ] x to be a number" ).ToLocalChecked() ); }
  if( !info[ 1 ]->IsNumber() ){ return Nan::ThrowError( Nan::New( "TowerDefense::addAgent - expected argument [ 1 ] y to be a number" ).ToLocalChecked() ); }
  if( info.Length() > 2 ){
    if( !info[ 2 ]->IsNumber() ){ return Nan::ThrowError( Nan::New( "TowerDefense::addAgent - expected argument [ 2 ] neighborDist to be a number" ).ToLocalChecked() ); }
    if( !info[ 3 ]->IsNumber() ){ return Nan::ThrowError( Nan::New( "TowerDefense::addAgent - expected argument [ 3 ] maxNeighbors to be an integer" ).ToLocalChecked() ); }
    if( !info[ 4 ]->IsNumber() ){ return Nan::ThrowError( Nan::New( "TowerDefense::addAgent - expected argument [ 4 ] timeHorizon to be a number" ).ToLocalChecked() ); }
    if( !info[ 5 ]->IsNumber() ){ return Nan::ThrowError( Nan::New( "TowerDefense::addAgent - expected argument [ 5 ] timeHorizonObst to be a number" ).ToLocalChecked() ); }
    if( !info[ 6 ]->IsNumber() ){ return Nan::ThrowError( Nan::New( "TowerDefense::addAgent - expected argument [ 6 ] radius to be a number" ).ToLocalChecked() ); }
    if( !info[ 7 ]->IsNumber() ){ return Nan::ThrowError( Nan::New( "TowerDefense::addAgent - expected argument [ 7 ] maxSpeed to be a number" ).ToLocalChecked() ); }
    if( info.Length() == 9 && ( !info[ 8 ]->IsArray() ) ){ return Nan::ThrowError( Nan::New( "TowerDefense::addAgent - expected argument [ 8 ] velocity to be an array" ).ToLocalChecked() ); }
  }
  TowerDefense* self = Nan::ObjectWrap::Unwrap<TowerDefense>( info.This() );
  if( info.Length() == 2 ){
    self->sim->addAgent( RVO::Vector2( ( float ) info[ 0 ]->NumberValue(), ( float ) info[ 1 ]->NumberValue() ) );
  }else if( info.Length() == 8 ){
    self->sim->addAgent(
      RVO::Vector2( ( float ) info[ 0 ]->NumberValue(), ( float ) info[ 1 ]->NumberValue() ),
      ( float ) info[ 2 ]->NumberValue(),//neighborDist
      ( size_t ) info[ 3 ]->IntegerValue(),//maxNeighbors
      ( float ) info[ 4 ]->NumberValue(),//timeHorizon
      ( float ) info[ 5 ]->NumberValue(),//timeHorizonObst
      ( float ) info[ 6 ]->NumberValue(),//radius
      ( float ) info[ 7 ]->NumberValue()//maxSpeed
    );
  }else if( info.Length() == 9 ){
    v8::Local<v8::Array> arrVelocity = v8::Local<v8::Array>::Cast( info[ 8 ] );
    if( arrVelocity->Length() != 2 ){
      return Nan::ThrowError( Nan::New( "TowerDefense::addAgent - expected argument [ 8 ] to be an array of length 2" ).ToLocalChecked() );
    }
    if( !arrVelocity->Get( 0 )->IsNumber() || !arrVelocity->Get( 1 )->IsNumber() ){
      return Nan::ThrowError( Nan::New( "TowerDefense::addAgent - expected argument [ 8 ] to be an array of 2 numbers x & y" ).ToLocalChecked() );
    }
    self->sim->addAgent(
      RVO::Vector2( ( float ) info[ 0 ]->NumberValue(), ( float ) info[ 1 ]->NumberValue() ),
      ( float ) info[ 2 ]->NumberValue(),//neighborDist
      ( size_t ) info[ 3 ]->IntegerValue(),//maxNeighbors
      ( float ) info[ 4 ]->NumberValue(),//timeHorizon
      ( float ) info[ 5 ]->NumberValue(),//timeHorizonObst
      ( float ) info[ 6 ]->NumberValue(),//radius
      ( float ) info[ 7 ]->NumberValue(),//maxSpeed
      RVO::Vector2( ( float ) arrVelocity->Get( 0 )->NumberValue(),//velocity x
                    ( float ) arrVelocity->Get( 1 )->NumberValue()//velocity y
      )
    );
  }
}

NAN_METHOD( TowerDefense::addObstacle ){
  if( !info[ 0 ]->IsArray() ){
    return Nan::ThrowError( Nan::New( "TowerDefense::addObstacle - expected argument 0 vertices to be an array" ).ToLocalChecked() );
  }
  TowerDefense* self = Nan::ObjectWrap::Unwrap<TowerDefense>( info.This() );
  v8::Local<v8::Array> arrVertices = v8::Local<v8::Array>::Cast( info[ 0 ] );
  int l = arrVertices->Length();
  std::vector< RVO::Vector2 > vertices;
  for( int i = 0; i < l; i++ ){
    v8::Local<v8::Array> arrVertice = v8::Local<v8::Array>::Cast( arrVertices->Get( i ) );
    if( arrVertice->Length() != 2 ){
      return Nan::ThrowError( Nan::New( "TowerDefense::addObstacle - expected array elements to ben arrays of length 2" ).ToLocalChecked() );
    }
    if( !arrVertice->Get( 0 )->IsNumber() || !arrVertice->Get( 1 )->IsNumber() ){
      return Nan::ThrowError( Nan::New( "TowerDefense::addObstacle - expected array elements to ben arrays of two numbers" ).ToLocalChecked() );
    }
    vertices.push_back( RVO::Vector2( ( float ) arrVertice->Get( 0 )->NumberValue(), ( float ) arrVertice->Get( 1 )->NumberValue() ) );
  }
  self->sim->addObstacle( vertices );
}

NAN_METHOD( TowerDefense::processObstacles ){
  TowerDefense* self = Nan::ObjectWrap::Unwrap<TowerDefense>( info.This() );
  self->sim->processObstacles();
  info.GetReturnValue().Set( Nan::Undefined() );
}

NAN_METHOD( TowerDefense::getGlobalTime ){
  TowerDefense* self = Nan::ObjectWrap::Unwrap<TowerDefense>( info.This() );
  info.GetReturnValue().Set( self->sim->getGlobalTime() );
}

NAN_METHOD( TowerDefense::getNumAgents ){
  TowerDefense* self = Nan::ObjectWrap::Unwrap<TowerDefense>( info.This() );
  info.GetReturnValue().Set( ( int ) self->sim->getNumAgents() );
}

NAN_METHOD( TowerDefense::getAgentPosition ){
  if( info.Length() != 1 ){
    return Nan::ThrowError( Nan::New( "TowerDefense::getAgentPosition - expected 1 argument agentNo" ).ToLocalChecked() );
  }
  if( !info[ 0 ]->IsNumber() ){
    return Nan::ThrowError( Nan::New( "TowerDefense::getAgentPosition - expected argument [ 0 ] agentNo to be an integer" ).ToLocalChecked() );
  }
  TowerDefense* self = Nan::ObjectWrap::Unwrap<TowerDefense>( info.This() );
  const RVO::Vector2& agentPos = self->sim->getAgentPosition( (size_t ) info[ 0 ]->IntegerValue() );
  v8::Local<v8::Object> AgentPosOb = Nan::New<v8::Object>();
  v8::Local<v8::String> xProp = Nan::New( "x" ).ToLocalChecked();
  v8::Local<v8::Value> xValue = Nan::New( agentPos.x() );
  AgentPosOb->Set( xProp, xValue );
  v8::Local<v8::String> yProp = Nan::New( "y" ).ToLocalChecked();
  v8::Local<v8::Value> yValue = Nan::New( agentPos.y() );
  AgentPosOb->Set( yProp, yValue );
  info.GetReturnValue().Set( AgentPosOb );
}

NAN_METHOD( TowerDefense::getAgentRadius ){
  TowerDefense* self = Nan::ObjectWrap::Unwrap<TowerDefense>( info.This() );
  info.GetReturnValue().Set( Nan::Undefined() );
}
//ICI
NAN_METHOD( TowerDefense::setAgentPrefVelocity ){
  if( info.Length() != 2 ){
    return Nan::ThrowError( Nan::New( "TowerDefense::setAgentPrefVelocity - expected 2 arguments" ).ToLocalChecked() );
  }
  if( !info[ 0 ]->IsNumber() ){
    return Nan::ThrowError( Nan::New( "TowerDefense::setAgentPrefVelocity - expected argument [ 0 ] agentNo to be a number" ).ToLocalChecked() );
  }
  if( !info[ 1 ]->IsArray() ){
    return Nan::ThrowError( Nan::New( "TowerDefense::setAgentPrefVelocity - expected argument [ 1 ] prefVelocity to be an array" ).ToLocalChecked() );
  }
  TowerDefense* self = Nan::ObjectWrap::Unwrap<TowerDefense>( info.This() );
  v8::Local<v8::Array> arrVelocity = v8::Local<v8::Array>::Cast( info[ 1 ] );
  if( arrVelocity->Length() != 2 ){
    return Nan::ThrowError( Nan::New( "TowerDefense::setAgentPrefVelocity - expected argument [ 1 ] prefVelocity to be an array of length 2" ).ToLocalChecked() );
  }
  if( !arrVelocity->Get( 0 )->IsNumber() || !arrVelocity->Get( 1 )->IsNumber() ){
    return Nan::ThrowError( Nan::New( "TowerDefense::setAgentPrefVelocity - expected argument [ 1 ] prefVelocity to be an array of two numbers x & y" ).ToLocalChecked() );
  }
  self->sim->setAgentPrefVelocity( ( size_t ) info[ 0 ]->IntegerValue(), RVO::normalize( RVO::Vector2( ( float ) arrVelocity->Get( 0 )->NumberValue(), ( float ) arrVelocity->Get( 1 )->NumberValue() ) ) );
}

NAN_METHOD( TowerDefense::doStep ){
  TowerDefense* self = Nan::ObjectWrap::Unwrap<TowerDefense>( info.This() );
  self->sim->doStep();
}

//using Coord = double;
//using N = uint32_t;
//using Point = ETP::Point<Coord>;

NAN_METHOD( TowerDefense::clip2triangulate ){
  // expected array [ [[pt1x,pt1y],[pt2x,pt2y]...],[[pt1x,pt1y],[pt2x,pt2y]...],...] where first nested array represents bounding polygon and following arrays the holes
  if( info.Length() < 1 ){
    return Nan::ThrowError( Nan::New( "TowerDefense::triangulate - expected 1 argument" ).ToLocalChecked() );
  }
  if( !info[ 0 ]->IsArray() ){
    return Nan::ThrowError( Nan::New( "TowerDefense::triangulate - expected argument [ 0 ] to be an array" ).ToLocalChecked() );
  }
  TowerDefense* self = Nan::ObjectWrap::Unwrap<TowerDefense>( info.This() );
  typedef p2t::Point usedPoint;
  typedef double Coord;
  v8::Local<v8::Array> mainArray = v8::Local<v8::Array>::Cast( info[ 0 ] );
  size_t mal = mainArray->Length();
  std::vector<std::vector<usedPoint>> polygons;

  // copy positions stored in info to std::vector<std::vector<c2t::Point>> polygons
  for( size_t ai = 0; ai < mal; ai++ ){
    v8::Local<v8::Array> polyArray = v8::Local<v8::Array>::Cast( mainArray->Get( ai ) );
    size_t pal = polyArray->Length();
    polygons.push_back( std::vector<usedPoint>() );
    std::vector<usedPoint>& usedPoly = polygons[ polygons.size() - 1 ];
    for( size_t pi = 0; pi < pal; pi ++ ){
      v8::Local<v8::Array> ptArr = v8::Local<v8::Array>::Cast( polyArray->Get( pi ) );
      usedPoly.push_back( usedPoint( (Coord) ptArr->Get( 0 )->NumberValue(), (Coord) ptArr->Get( 1 )->NumberValue() ) );
    }
  }
  try{
    self->space->buildFromPolyLines<usedPoint>( polygons );
  }catch( const std::invalid_argument& e ){
    return Nan::ThrowError( Nan::New( e.what() ).ToLocalChecked() );
  }

  info.GetReturnValue().Set( self->space->getTriangles() );

}

NAN_METHOD( TowerDefense::tryClipper ){
  TowerDefense* self = Nan::ObjectWrap::Unwrap<TowerDefense>( info.This() );

  ClipperLib::Paths subj( 1 ), clip( 1 ), solution;

  subj[ 0 ].push_back( ClipperLib::IntPoint( 0.0, 0.0 ) );
  subj[ 0 ].push_back( ClipperLib::IntPoint( 100.0, 0.0 ) );
  subj[ 0 ].push_back( ClipperLib::IntPoint( 100.0, 100.0 ) );
  subj[ 0 ].push_back( ClipperLib::IntPoint( 0.0, 100.0 ) );

  clip[ 0 ].push_back( ClipperLib::IntPoint( 50.0, 50.0 ) );
  clip[ 0 ].push_back( ClipperLib::IntPoint( 150.0, 50.0 ) );
  clip[ 0 ].push_back( ClipperLib::IntPoint( 150.0, 80.0 ) );
  clip[ 0 ].push_back( ClipperLib::IntPoint( 50.0, 80.0 ) );

  ClipperLib::ClipType clipType;
  int inputType = info[ 0 ]->NumberValue();
  if( inputType == 0 ){
    clipType = ClipperLib::ctIntersection;
  }else if( inputType == 1 ){
    clipType = ClipperLib::ctUnion;
  }else if( inputType == 2 ){
    clipType = ClipperLib::ctDifference;
  }else if( inputType == 3 ){
    clipType = ClipperLib::ctXor;
  }

  ClipperLib::Clipper c;
  c.AddPaths(subj, ClipperLib::ptSubject, true);
  c.AddPaths(clip, ClipperLib::ptClip, true);
  c.Execute(clipType, solution, ClipperLib::pftNonZero, ClipperLib::pftNonZero);

  size_t sl = solution.size();
  v8::Local<v8::Array> ret = Nan::New<v8::Array>();
  for( size_t i = 0; i < sl; i++ ){
    ClipperLib::Path& sLine  = solution[ i ];
    size_t ls = sLine.size();
    v8::Local<v8::Array> rLine = Nan::New<v8::Array>();
    for( size_t y = 0; y < ls; y++ ){
      int _x = sLine[ y ].X;
      int _y = sLine[ y ].Y;
      rLine->Set( y * 2, Nan::New(  _x) );
      rLine->Set( y * 2 + 1, Nan::New( _y ) );
    }
    ret->Set( i, rLine );
  }
  info.GetReturnValue().Set( ret );
}

NAN_METHOD( TowerDefense::tryPoly2Tri ){
  TowerDefense* self = Nan::ObjectWrap::Unwrap<TowerDefense>( info.This() );
  std::vector<p2t::Point*> polyline;
  polyline.push_back( new p2t::Point( 10.0, 10.0 ) );
  polyline.push_back( new p2t::Point( 25.0, 25.0 ) );
  polyline.push_back( new p2t::Point( 50.0, 10.0 ) );
  polyline.push_back( new p2t::Point( 100.0, 10.0 ) );
  polyline.push_back( new p2t::Point( 100.0, 50.0 ) );
  polyline.push_back( new p2t::Point( 80.0, 50.0 ) );

  polyline.push_back( new p2t::Point( 60.0, 30.0 ) );
  polyline.push_back( new p2t::Point( 40.0, 50.0 ) );

  polyline.push_back( new p2t::Point( 10.0, 50.0 ) );

  p2t::CDT cdt( polyline );
  //p2t::Sweep sweep( sweepctx );
  cdt.Triangulate();
  std::vector<p2t::Triangle*> tris = cdt.GetTriangles();
  v8::Local<v8::Array> ret = Nan::New<v8::Array>();
  size_t l = tris.size();
  for( size_t i = 0; i < l; i++ ){
    v8::Local<v8::Array> triRet = Nan::New<v8::Array>();
    p2t::Triangle* t = tris[ i ];
    p2t::Point* pt0 = t->GetPoint( 0 );
    triRet->Set( 0,     Nan::New(pt0->x) );
    triRet->Set( 1, Nan::New(pt0->y) );
    p2t::Point* pt1 = t->GetPoint( 1 );
    triRet->Set( 2, Nan::New(pt1->x) );
    triRet->Set( 3, Nan::New(pt1->y) );
    p2t::Point* pt2 = t->GetPoint( 2 );
    triRet->Set( 4, Nan::New(pt2->x) );
    triRet->Set( 5, Nan::New(pt2->y) );
    ret->Set( i, triRet );
  }
  info.GetReturnValue().Set( ret );
}

NAN_METHOD( TowerDefense::testAngle ){
  TowerDefense* self = Nan::ObjectWrap::Unwrap<TowerDefense>( info.This() );
  int e1pt1x = info[ 0 ]->IntegerValue();
  int e1pt1y = info[ 1 ]->IntegerValue();
  int e1pt2x = info[ 2 ]->IntegerValue();
  int e1pt2y = info[ 3 ]->IntegerValue();

  int e2pt1x = info[ 4 ]->IntegerValue();
  int e2pt1y = info[ 5 ]->IntegerValue();
  int e2pt2x = info[ 6 ]->IntegerValue();
  int e2pt2y = info[ 7 ]->IntegerValue();

  int x1, y1, x2, y2;
  /*
  if( e1pt1x == e2pt1x &&  e1pt1y == e2pt1y ){

    x1 = e1pt1x - e1pt2x;
    y1 = e1pt1y - e1pt2y;

    x2 = e2pt2x - e2pt1x;
    y2 = e2pt2y - e2pt1y;
  }else if( e1pt1x == e2pt2x &&  e1pt1y == e2pt2y){
    x1 = e1pt1x - e1pt2x;
    y1 = e1pt1y - e1pt2y;

    x2 = e2pt1x - e2pt2x;
    y2 = e2pt1y - e2pt2y;
  }else if( e1pt2x == e2pt1x &&  e1pt2y == e2pt1y ){
    x1 = e1pt2x - e1pt1x;
    y1 = e1pt2y - e1pt1y;

    x2 = e2pt2x - e2pt1x;
    y2 = e2pt2y - e2pt1y;
  }else if( e1pt2x == e2pt2x &&  e1pt2y == e2pt2y ){
    x1 = e1pt2x - e1pt1x;
    y1 = e1pt2y - e1pt1y;

    x2 = e2pt1x - e2pt2x;
    y2 = e2pt1y - e2pt2y;
  }
  x1 = e1pt1x - e1pt2x;
  y1 = e1pt1y - e1pt2y;

  x2 = e2pt2x - e2pt1x;
  y2 = e2pt2y - e2pt1y;
  */
  /*
  D Pv1 = std::sqrt( (D) ( std::pow( vertex->x - pt1->x , 2 ) + std::pow( vertex->y - pt1->y, 2 ) ) );
  D Pv2 = std::sqrt( (D) ( std::pow( vertex->x - pt2->x , 2 ) + std::pow( vertex->y - pt2->y, 2 ) ) );
  D P12 = std::sqrt( (D) ( std::pow( pt1->x - pt2->x , 2 ) + std::pow( pt1->y - pt2->y, 2 ) ) );
  return std::acos( ( std::pow( Pv1, 2 ) + std::pow( Pv2, 2 ) - std::pow( P12, 2 ) ) / ( 2.0 * Pv1 * Pv2 ) );
  */
  x1 = e1pt2x - e1pt1x;
  y1 = e1pt2y - e1pt1y;

  x2 = e2pt2x - e2pt1x;
  y2 = e2pt2y - e2pt1y;
  //P x2 = pt2->x - vertex->x;
  //P y2 = pt2->y - vertex->y;
  int dot = x1 * x2 + y1 * y2;
  int det = x1 * y2 - y1 * x2;
  float angle = std::abs( (float) std::atan2( (float) det, (float) dot ) );
  info.GetReturnValue().Set( Nan::New( angle ) );

}

NAN_METHOD( TowerDefense::getTriangleId ){
  //TowerDefense* self = Nan::ObjectWrap::Unwrap<TowerDefense>( info.This() );
  //info.GetReturnValue().Set( self->space->getSectors() );
  //info.GetReturnValue().Set( self->space->getBounds() );

  TowerDefense* self = Nan::ObjectWrap::Unwrap<TowerDefense>( info.This() );
  tPos x = info[ 0 ]->IntegerValue();
  tPos y = info[ 1 ]->IntegerValue();
  tDecim scale = info[ 2 ]->NumberValue();
  std::shared_ptr<ETP::Triangle<tPos,tDecim>> t = self->space->getTriangleWithPoint( ETP::Point<tPos,tDecim>( x, y ), scale );
  if( t == nullptr ){
    info.GetReturnValue().Set( Nan::New( -1 ) );
  }else{
    info.GetReturnValue().Set( Nan::New( t->getId() ) );
  }

}

NAN_METHOD( TowerDefense::getSectors ){
  TowerDefense* self = Nan::ObjectWrap::Unwrap<TowerDefense>( info.This() );
  info.GetReturnValue().Set( self->space->getSectors() );
}

NAN_METHOD( TowerDefense::testTRAStar ){
  TowerDefense* self = Nan::ObjectWrap::Unwrap<TowerDefense>( info.This() );
  tPos sx = info[ 0 ]->IntegerValue();
  tPos sy = info[ 1 ]->IntegerValue();
  std::shared_ptr<ETP::Triangle<tPos,tDecim>> st = self->space->getTriangleWithPoint( ETP::Point<tPos,tDecim>( sx, sy ) );
  tPos gx = info[ 2 ]->IntegerValue();
  tPos gy = info[ 3 ]->IntegerValue();
  std::shared_ptr<ETP::Triangle<tPos,tDecim>> gt = self->space->getTriangleWithPoint( ETP::Point<tPos,tDecim>( gx, gy ) );
  //std::string ret = self->space->abstractTriangulationSearch( st, gt );
  //std::vector<std::shared_ptr<ETP::Triangle<tPos,tDecim>>> ret = self->space->abstractTriangulationSearch( st, gt );
  std::vector<std::shared_ptr<ETP::Triangle<tPos,tDecim>>> ret = self->space->abstractTriangulationSearch( ETP::Point<tPos,tDecim>( sx, sy ),  ETP::Point<tPos,tDecim>( gx, gy ), 2.0 );

  v8::Local<v8::Array> ids =  self->space->listTrianglesIds( ret );
  v8::Local<v8::Array> tris = self->space->convertTriangles( ret );

  v8::Local<v8::Object> ob = Nan::New<v8::Object>();
  v8::Local<v8::String> triProp = Nan::New( "triangles" ).ToLocalChecked();
  v8::Local<v8::String> idProp = Nan::New( "ids" ).ToLocalChecked();
  ob->Set( triProp, tris );
  ob->Set( idProp, ids );

  info.GetReturnValue().Set( ob );
}

NAN_METHOD( TowerDefense::testTRAStarScale ){
  TowerDefense* self = Nan::ObjectWrap::Unwrap<TowerDefense>( info.This() );
  tPos sx = info[ 0 ]->IntegerValue();
  tPos sy = info[ 1 ]->IntegerValue();

  tPos gx = info[ 2 ]->IntegerValue();
  tPos gy = info[ 3 ]->IntegerValue();
  tDecim scale = info[ 4 ]->NumberValue();
  std::shared_ptr<ETP::Triangle<tPos,tDecim>> st = self->space->getTriangleWithPoint( ETP::Point<tPos,tDecim>( sx, sy ), scale );
  std::shared_ptr<ETP::Triangle<tPos,tDecim>> gt = self->space->getTriangleWithPoint( ETP::Point<tPos,tDecim>( gx, gy ), scale );
  //std::string ret = self->space->abstractTriangulationSearch( st, gt );
  //std::vector<std::shared_ptr<ETP::Triangle<tPos,tDecim>>> ret = self->space->abstractTriangulationSearch( st, gt );
  std::vector<std::shared_ptr<ETP::Triangle<tPos,tDecim>>> ret = self->space->abstractTriangulationSearch( ETP::Point<tPos,tDecim>( sx, sy ),  ETP::Point<tPos,tDecim>( gx, gy ), 2.0, scale );

  v8::Local<v8::Array> ids =  self->space->listTrianglesIds( ret );
  v8::Local<v8::Array> tris = self->space->convertTriangles( ret );

  v8::Local<v8::Object> ob = Nan::New<v8::Object>();
  v8::Local<v8::String> triProp = Nan::New( "triangles" ).ToLocalChecked();
  v8::Local<v8::String> idProp = Nan::New( "ids" ).ToLocalChecked();
  ob->Set( triProp, tris );
  ob->Set( idProp, ids );

  info.GetReturnValue().Set( ob );
}
