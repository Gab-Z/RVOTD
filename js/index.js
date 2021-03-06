/*
window.eval = function () {
  throw new Error(`Sorry, this app does not support window.eval().`)
}
*/
const ipc = require( 'electron' ).ipcRenderer;
const bindings = require( "bindings" );
const towerdef = bindings( "towerdef" );
const td = new towerdef.TowerDefense();
console.log( td.get() );
/*
td.setTimeStep( 1 );
td.setAgentDefaults( 60, 10, 200,300, 30, 25 );

var agentCounter = 0;
var goals = [
  { x:400, y: 200 },
  { x:0, y:200 },
  {x:200, y: 400 },
  {x:200,y: 0},

  { x:400, y: 400 },
  { x:0, y:400 },
  {x:0, y: 0 },
  {x:400,y: 0}
];


addAgent( 0, 200, [ 1, 0 ] );
addAgent( 400, 200, [ -1, 0 ] );
addAgent( 200, 0, [ 0, 1 ] );
addAgent( 200, 400, [ 0, -1 ] );

addAgent( 0, 0, [ 1, 1 ] );
addAgent( 400, 0, [ -1, 1 ] );
addAgent( 400, 400, [ -1, -1 ] );
addAgent( 0, 400, [ 1, -1 ] );
let wall = addBlock( 50, 50, 200, 200 ) ;
//wall.style.transformOrigin = "50% 50%"
//wall.style.transform = "translate( 175px, 175px )"

var ct = 0;
var maxct = 10000000;

function move( ){
  td.doStep();
  for( let i = 0; i < agentCounter; i++ ){
    let agentPos = td.getAgentPosition( i );
    document.getElementById( "agent"+ i  ).style.transform = "translate( "+agentPos.x+"px, "+agentPos.y+"px )";
    let goal = goals[ i ];
    let angle = Math.atan2( goal.y - agentPos.y,  goal.x - agentPos.x );
    angle = (angle - Math.PI/72) + Math.random() * (Math.PI/72);
    td.setAgentPrefVelocity( i, [ Math.cos( angle ), Math.sin(angle) ] );
  }
  window.requestAnimationFrame( move );
}

function createSprite( id ){
  let sp = document.createElement( "div" );
  sp.style.width = "30px";
  sp.style.height = "30px";
  sp.style.backgroundColor = "black";
  sp.style.borderRadius = "50%";
  sp.id = id;
  sp.style.position = "absolute";
  return sp;
}
function addAgent( x, y, prefVelocity){
  td.addAgent( x, y );
  let agentId="agent"+agentCounter;
  document.body.appendChild( createSprite( agentId ) );
  document.getElementById( agentId ).style.transform = "translate( "+x+"px, "+y+"px )";
  td.setAgentPrefVelocity( agentCounter, prefVelocity );
  agentCounter++;
}
document.body.addEventListener( "click", e => { move(); } );

function addBlock( w, h, x, y ){
  let arr = [];
  arr.push( [ x - w/2, y - h/2 ],
            [ x + w/2, y - h/2 ],
            [ x + w/2, y + h/2 ],
            [ x - w/2, y + h/2 ]
  );
  td.addObstacle( arr );
  let sp =  document.body.appendChild( document.createElement( "div" ) );
  sp.style.margin = 0;
  sp.style.padding = 0;
  sp.style.width = w+"px";
  sp.style.height = h+"px";
  sp.style.backgroundColor = "red";
  sp.style.position = "absolute";
  sp.style.transform = "translate( "+ x +"px, "+ y +"px )"
  td.processObstacles();
  return sp;

}
*/
var displayScale = 0.5,
    agentRadius = 15;

var setCanvas = ( id, background )=>{
  let cv  = document.body.appendChild( document.createElement( "canvas" ) );
  cv.id = id;
  cv.style.position = "absolute";
  cv.style.width = "1000px";
  cv.style.height = "1000px";
  cv.style.margin = 0;
  cv.style.padding = 0;
  cv.width = 1000;
  cv.height = 1000;
  if( background ) cv.style.backgroundColor = background;
}
var agents = [];
var setAgent = ( x, y, goal, options )=>{
  td.addAgent( x, y );
  let agentOb = { x: x, y: y, goal: goal };
  if( options ){agentOb.options = options;}
  let angle = Math.atan2( agentOb.goal.y - agentOb.y,  agentOb.goal.x - agentOb.x );
  angle = ( angle - Math.PI / 72 ) + Math.random() * ( Math.PI / 72 );
  td.setAgentPrefVelocity( agents.length, [ Math.cos( angle ), Math.sin(angle) ] );
  agents.push( agentOb );
}
var blocks = [];
var setBlock = ( x, y, w, h )=>{
  let arrVertices = [
    [ x - w/2, y - h/2 ],
    [ x + w/2, y - h/2 ],
    [ x + w/2, y + h/2 ],
    [ x - w/2, y + h/2 ]
  ];
  let blockOb = { x: x, y: y, w: w, h: h };
  blocks.push(blockOb );
  td.addObstacle( arrVertices );
}
var drawFrame = ()=>{
  let cv = document.getElementById( "canvas" ),
      ctx = cv.getContext( "2d" );
  ctx.clearRect( 0, 0, cv.width, cv.height );
  drawBlocks();
  drawAgents();
}
var drawAgents = ()=>{
  let cv = document.getElementById( "canvas" ),
      ctx = cv.getContext( "2d" ),
      ds = displayScale;
      ctx.lineWidth = 0.5;
  agents.forEach( ( agentOb, index )=>{
    ctx.beginPath();
    ctx.moveTo(( agentOb.x+agentOb.w/2)*ds, agentOb.y * ds);
    ctx.arc( agentOb.x * ds, agentOb.y * ds, agentRadius * ds, 0, 2*Math.PI );
    ctx.stroke();
    if( agentOb.options && agentOb.options.color ){
      ctx.fillStyle = agentOb.options.color;
      ctx.fill();
    }
  } );
}
var drawBlocks = ()=>{
  let cv = document.getElementById( "canvas" ),
      ctx = cv.getContext( "2d" ),
      ds = displayScale;
  blocks.forEach( ( blockOb, index )=>{
    ctx.beginPath();
    ctx.moveTo( (blockOb.x - ( blockOb.w / 2 )) * ds, (blockOb.y - ( blockOb.h / 2 )) * ds );
    ctx.rect( (blockOb.x - ( blockOb.w / 2 )) * ds, (blockOb.y - ( blockOb.h / 2 )) * ds, blockOb.w * ds, blockOb.h * ds );
    ctx.stroke();
  } );
}
var updateAgents = ()=>{
  td.doStep();
  agents.forEach( ( agentOb, index )=>{
    let agentPos = td.getAgentPosition( index );
    agentOb.x = agentPos.x;
    agentOb.y = agentPos.y;
    let angle = Math.atan2( agentOb.goal.y - agentOb.y,  agentOb.goal.x - agentOb.x );
    angle = (angle - Math.PI/72) + Math.random() * (Math.PI/72);
    td.setAgentPrefVelocity( index, [ Math.cos( angle ), Math.sin(angle) ] );
  });
}
var loop = ()=>{
  updateAgents();
  drawFrame();
  window.requestAnimationFrame( loop );
}

setCanvas( "canvas", "#f0d5ba" );
setCanvas( "tmpCanvas" );
document.getElementById( "tmpCanvas" ).style.pointerEvents = "none";

polyMaker.createUI( document.body, document.getElementById( "canvas" ) );



td.setTimeStep( 2 );
td.setAgentDefaults( 80, 40, 100,100, agentRadius, 25 );
/*
let agentsData = [
  { x: 0,y: 0, goal: { x: 1000, y: 0 } },
  { x: 0,y: 40, goal: { x: 1000, y: 40 } },
  { x: 0,y: 80, goal: { x: 1000, y: 80 } },
  { x: 0,y: 120, goal: { x: 1000, y: 120 } },
  { x: 0,y: 160, goal: { x: 1000, y: 160 } },
  { x: 0,y: 200, goal: { x: 1000, y: 200 } },
  { x: 0,y: 240, goal: { x: 1000, y: 240 } },
  { x: 0,y: 280, goal: { x: 1000, y: 280 } },
  { x: 0,y: 320, goal: { x: 1000, y: 320 } },
  { x: 0,y: 360, goal: { x: 1000, y: 360 } },
  { x: 1000,y: 0, goal: { x: 0, y: 0 } },
  { x: 1000,y: 30, goal: { x: 0, y: 30 } },
  { x: 1000,y: 60, goal: { x: 0, y: 60 } },
  { x: 1000,y: 90, goal: { x: 0, y: 90 } },
  { x: 1000,y: 120, goal: { x: 0, y: 120 } },
  { x: 1000,y: 150, goal: { x: 0, y: 150 } },
  { x: 1000,y: 180, goal: { x: 0, y: 180 } },
  { x: 1000,y: 210, goal: { x: 0, y: 210 } },
  { x: 1000,y: 240, goal: { x: 0, y: 240 } },
  { x: 1000,y: 270, goal: { x: 0, y: 270 } }
];
*/


let agentsData = [];

for(let r = 0; r < 20; r++ ){
  for(let c = 0; c < 10;c++ ){
    agentsData.push( { x: c * agentRadius * 2, y: r * agentRadius * 2, goal: { x: 1000, y: r * agentRadius * 2 }, options:{ color: "#a071e3" } } );
  }
}
for(let r = 0; r < 20; r++ ){
  for(let c = 0; c < 10;c++ ){
    agentsData.push( { x: 1000 - (c * agentRadius * 2), y: r * agentRadius * 2, goal: { x: 0, y: r * agentRadius * 2 }, options:{ color: "#eea321" } } );
  }
}

agentsData.forEach( ( agentData, id ) => {

  setAgent( agentData.x, agentData.y, agentData.goal,agentData.options );
} );


for(let i = 0; i < 10; i++ ){
  setBlock( 500, i * 100 + 50, 50, 50 )
}
for(let i = 0; i < 10; i++ ){
  setBlock( 350, i * 100, 25, 25 )
}
for(let i = 0; i < 10; i++ ){
  setBlock( 650, i * 100, 25, 25 )
}
td.processObstacles();
/*
drawFrame();
document.body.addEventListener( "click", e => {
  loop();
} );
*/

const inverseSearchBut = document.body.appendChild( document.createElement("button"));
inverseSearchBut.style.position = "absolute";
inverseSearchBut.textContent = "inverse search";
inverseSearchBut.style.top = 0;
inverseSearchBut.style.left = 0;

var revertClickInfo = { start:{}, end:{} };
inverseSearchBut.addEventListener( "click", e => {
  if( revertClickInfo.start != {} && revertClickInfo.end != {} ){
    let tmpEnd = revertClickInfo.end;
    revertClickInfo.end = revertClickInfo.start;
    revertClickInfo.start = tmpEnd;
    let SCALE =1;
    let funnel = td.testTRAStarScale( revertClickInfo.start.x, revertClickInfo.start.y, revertClickInfo.end.x, revertClickInfo.end.y, SCALE );
    console.log(  stringIndentifyObject( funnel ));
    console.log( "end set at x:"+revertClickInfo.end.x+"/y:"+revertClickInfo.end.y+"/startId:"+td.getTriangleId(revertClickInfo.start.x, revertClickInfo.start.y, SCALE)+"/endId:"+td.getTriangleId(revertClickInfo.end.x, revertClickInfo.end.y, SCALE));
    if( funnel.triangles[ 0 ].positions )fillTriangles( revertClickInfo, funnel.triangles, document.getElementById( "tmpCanvas" ), SCALE, "rgba(235, 221, 56, 0.5)" );
  }
})


const triangulationBut = document.body.appendChild( document.createElement("button"));
triangulationBut.style.position = "absolute";
triangulationBut.textContent = "triangulate";
triangulationBut.style.top = 0;
triangulationBut.style.left = "920px";
var clickInfo = { start:{}, end:{}, first: true };

triangulationBut.addEventListener( "click", e => {
  if( !polyMaker.array ){ return false;}
  let triangulation = td.clip2triangulate( polyMaker.array );
  console.log( JSON.stringify( triangulation ) );
  let SCALE =1;
  drawTriangles( triangulation, document.getElementById( "canvas" ), SCALE );
  clickInfo = { start:{}, end:{}, first: true };
  document.getElementById( "canvas" ).addEventListener( "mouseup", e =>{
    let x =  e.clientX,
        y = e.clientY,
        triId = td.getTriangleId(x,y, SCALE);
    if( triId < 0 ){ return 0; }
    if( clickInfo.first == true ){
      clickInfo.start.x = x;
      clickInfo.start.y = y;
      clickInfo.first = false,
      cv = document.getElementById( "tmpCanvas" ),
      cv.getContext( "2d" ).clearRect( 0, 0, cv.width, cv.height );
      revertClickInfo = { start:{}, end:{} };
    console.log( "start set at x:"+x+"/y:"+y+"/id:"+td.getTriangleId(x,y, SCALE));

    }else{
        clickInfo.end.x = x;
        clickInfo.end.y = y;

        let funnel = td.testTRAStarScale( clickInfo.start.x, clickInfo.start.y, clickInfo.end.x, clickInfo.end.y, SCALE );

        //console.log( JSON.stringify( funnel ));
        console.log(  stringIndentifyObject( funnel ));

        console.log( "end set at x:"+x+"/y:"+y+"/startId:"+td.getTriangleId(clickInfo.start.x, clickInfo.start.y, SCALE)+"/endId:"+td.getTriangleId(clickInfo.end.x, clickInfo.end.y, SCALE));
        if( funnel.triangles[ 0 ].positions )fillTriangles( clickInfo, funnel.triangles, document.getElementById( "tmpCanvas" ), SCALE, "rgba(235, 221, 56, 0.5)" );

        revertClickInfo = clickInfo;
        clickInfo = { start:{}, end:{}, first: true };
        //
    }
  })
})

let polypts = [
  [ [ 10, 10 ],  [ 90, 10 ], //top
    [ 90, 40 ], [ 110, 20 ], [ 115, 45 ], [ 110, 50 ], [ 90, 50 ],

  [ 90, 90 ], [ 70, 90 ],[ 70, 95 ], [ 90, 95 ],[ 90, 110], [ 70, 110 ], [ 60, 115 ], [ 50, 125 ], [ 50, 90 ],  [ 10, 90 ]   ],
  [ [ 20, 20 ], [ 30, 20 ], [ 30, 30 ], [ 20, 30 ] ],
  [ [ 50, 20 ], [ 60, 20 ], [ 60, 30 ], [ 50, 30 ] ],
  [ [ 20, 60 ], [ 30, 60 ], [ 30, 70 ], [ 20, 70 ] ],
  [ [ 70, 70 ], [ 85, 70], [ 85, 80 ], [ 70, 80 ] ],
  [ [ 88, 43 ], [ 105, 43 ], [ 105, 45 ],[ 88, 46 ] ],
  [[ 50, 50], [ 60, 50], [ 60,60 ], [50, 60 ]]
];







/*
let indices =  td.triangulate( polypts );
console.log( JSON.stringify( indices ) );
drawPoly( polypts, indices, document.getElementById( "canvas" ), 3 );

///////

console.log( "##### starting triangulation #####" );
let triangulation = td.buildTriangulation( polypts );
console.log( "##### ending triangulation #####" );
console.log( "##### LENGTH : " +  triangulation.length );
console.log( JSON.stringify( triangulation ) );
drawTriangulation( triangulation, document.getElementById( "canvas" ), 3 );
*/

/* HERE GOOD PART FOR TRIANGULATION SEARCH
let triangulation = td.clip2triangulate( polypts );
console.log( JSON.stringify( triangulation ) );
let SCALE =8;
drawTriangles( triangulation, document.getElementById( "canvas" ), SCALE );
drawWidths( triangulation,[ { triIndex: ( Math.floor( Math.random() * triangulation.length ) ), edgeIndex: ( Math.floor( Math.random() * 3 ) ) } ], document.getElementById( "canvas" ), SCALE );
*/

//drawSectors( td.getSectors(), document.getElementById( "canvas" ), SCALE )

/*
 let  pt1 = [ 20, 20 ],
      pt2 = [ 20, 40],
      pt3 = [ 20, 20],
      pt4 = [ 40, 30];

let a = td.testAngle( pt1[ 0 ], pt1[ 1 ], pt2[ 0 ], pt2[ 1 ], pt3[ 0 ], pt3[ 1 ], pt4[ 0 ], pt4[ 1 ] )
console.log( a + " / " + a / Math.PI*180 );
drawAngle( pt1[ 0 ], pt1[ 1 ], pt2[ 0 ], pt2[ 1 ], pt3[ 0 ], pt3[ 1 ], pt4[ 0 ], pt4[ 1 ], document.getElementById( "canvas" ), 4.5 );
*/


//console.log( JSON.stringify( td.getLevel3Paths() ) );
//let clip = td.tryClipper( 3 );
//let clip = td.tryPoly2Tri();
//console.log( JSON.stringify( clip ) );
 //drawClip( clip[ 0 ], document.getElementById( "canvas" ), SCALE )
// drawP2T( clip, document.getElementById( "canvas" ), SCALE )
/*
document.getElementById( "canvas" ).addEventListener( "click", e =>{
//  console.log( e.clientX / SCALE + " / " + e.clientY / SCALE )
  let id = td.getTriangleId( Math.round(e.clientX / SCALE), Math.round(e.clientY / SCALE ) );
//  console.log( "id : " + id );
console.log( JSON.stringify( id ) );
})
*/
/* HERE MOUSE EVENTS FOR SEARCH
clickInfo = { start:{}, end:{}, first: true };
document.getElementById( "canvas" ).addEventListener( "mouseup", e =>{

  let x =  e.clientX,
      y = e.clientY,
      triId = td.getTriangleId(x,y, SCALE);
  if( triId < 0 ){ return 0; }
  if( clickInfo.first == true ){
    clickInfo.start.x = x;
    clickInfo.start.y = y;


    clickInfo.first = false,
    cv = document.getElementById( "tmpCanvas" ),
    cv.getContext( "2d" ).clearRect( 0, 0, cv.width, cv.height );
  console.log( "start set at x:"+x+"/y:"+y+"/id:"+td.getTriangleId(x,y, SCALE));

  }else{
      clickInfo.end.x = x;
      clickInfo.end.y = y;

      let funnel = td.testTRAStarScale( clickInfo.start.x, clickInfo.start.y, clickInfo.end.x, clickInfo.end.y, SCALE );

      console.log( JSON.stringify( funnel ));
      console.log( "end set at x:"+x+"/y:"+y+"/startId:"+td.getTriangleId(clickInfo.start.x, clickInfo.start.y, SCALE)+"/endId:"+td.getTriangleId(clickInfo.end.x, clickInfo.end.y, SCALE));
      if( funnel.triangles[ 0 ].positions )fillTriangles( funnel.triangles, document.getElementById( "tmpCanvas" ), SCALE, "rgba(235, 221, 56, 0.5)" );
      clickInfo = { start:{}, end:{}, first: true };
      //
  }
})

*/


/*
let satPoly1 = [ 0.0, 0.0,    10.0, 0.0,  5.0, 10.0 ],
    satPoly2 = [ 8.0, 0.0,    15.0, 0.0,  15.0, 10.0, 0.0, 10.0 ],
    collide = td.testSAT( satPoly1, satPoly2 );
    console.log( "SAT : " + collide );
*/
function drawSectors( sectors, _cv, _scale ){
  let nbRow = sectors.sectors.length,
      nbCol = sectors.sectors[ 0 ].length,
      ctx = _cv.getContext( "2d" ),
      scale = _scale || 1;
  ctx.beginPath();
  for( let y = 0; y < nbRow; y++ ){
    ctx.moveTo( sectors.minX * scale, ( sectors.minY + y * sectors.height ) * scale );
    ctx.lineTo( (sectors.minX + nbCol * sectors.width) * scale, ( sectors.minY + y * sectors.height ) * scale );
  }
  for( let x = 0; x < nbCol; x++ ){
    ctx.moveTo( ( sectors.minX + x * sectors.width ) * scale,  sectors.minY * scale );
    ctx.lineTo( ( sectors.minX + x * sectors.width ) * scale, ( sectors.minY + nbRow * sectors.height ) * scale );
  }
  ctx.stroke();
  ctx.font = '11px serif';
  ctx.textBaseline = "top";
  ctx.textAlign = "left"
  sectors.sectors.forEach( ( row, y )=>{
    row.forEach( ( cell, x )=>{
      let str = "";
      cell.forEach( id =>{
        str += id+",";
      })
      ctx.strokeText( str, ( x * sectors.width + sectors.minX ) * scale, ( y * sectors.height + sectors.minY ) * scale );
    })
  })
}

function getTriCentroid(_coord){
  let coord = [ [_coord[ 0 ], _coord[ 1 ] ], [ _coord[ 2 ], _coord[ 3 ] ], [ _coord[ 4 ], _coord[ 5 ] ] ];
	var center = coord.reduce(function (x,y) {
		return [x[0] + y[0]/coord.length, x[1] + y[1]/coord.length]
	}, [0,0])
	return center;
}

function drawTriangles( _triangulation, _cv, _scale ){
  let ctx = _cv.getContext( "2d" ),
      scale = _scale || 1,
      nbTri = _triangulation.length;
  ctx.font = '9px serif';
  ctx.textBaseline = "middle";
  ctx.textAlign = "center"
  _triangulation.forEach( tri =>{
    let p = tri.positions,
        c = tri.constrained;
    ctx.strokeStyle = ( c[ 0 ] ? "red" : "black" );
    ctx.beginPath();
    ctx.moveTo( p[ 0 ] * scale,  p[ 1 ] * scale );
    ctx.lineTo( p[ 2 ] * scale,  p[ 3 ] * scale );
    ctx.stroke();
    ctx.strokeStyle = ( c[ 1 ] ? "red" : "black" );
    ctx.beginPath();
    ctx.moveTo( p[ 2 ] * scale,  p[ 3 ] * scale );
    ctx.lineTo( p[ 4 ] * scale,  p[ 5 ] * scale );
    ctx.stroke();
    ctx.strokeStyle = ( c[ 2 ] ? "red" : "black" );
    ctx.beginPath();
    ctx.moveTo( p[ 4 ] * scale,  p[ 5 ] * scale );
    ctx.lineTo( p[ 0 ] * scale,  p[ 1 ] * scale );
    ctx.stroke();
    let center = getTriCentroid( p );
    ctx.strokeStyle = "black"
    ctx.strokeText( tri.id +"-("+tri.level+")", center[ 0 ]*scale, center[ 1 ]*scale);
    let edgeC1x = p[ 0 ] * scale + ((p[ 2 ] * scale - p[ 0 ] * scale) / 2),
        edgeC1y = p[ 1 ] * scale + ((p[ 3 ] * scale - p[ 1 ] * scale) / 2),
        edgeC2x = p[ 2 ] * scale + ((p[ 4 ] * scale - p[ 2 ] * scale) / 2),
        edgeC2y = p[ 3 ] * scale + ((p[ 5 ] * scale - p[ 3 ] * scale) / 2),
        edgeC3x = p[ 4 ] * scale + ((p[ 0 ] * scale - p[ 4 ] * scale) / 2),
        edgeC3y = p[ 5 ] * scale + ((p[ 1 ] * scale - p[ 5 ] * scale) / 2),
        cToE = 0.75,
        to1x = (center[ 0 ]*scale) + ( edgeC1x - center[ 0 ]*scale ) * cToE,
        to1y = (center[ 1 ]*scale) + ( edgeC1y - center[ 1 ]*scale ) * cToE,
        to2x = (center[ 0 ]*scale) + ( edgeC2x - center[ 0 ]*scale ) * cToE,
        to2y = (center[ 1 ]*scale) + ( edgeC2y - center[ 1 ]*scale ) * cToE,
        to3x = (center[ 0 ]*scale) + ( edgeC3x - center[ 0 ]*scale ) * cToE,
        to3y = (center[ 1 ]*scale) + ( edgeC3y - center[ 1 ]*scale ) * cToE,
        a = tri.lowerBounds,
        n = tri.nodes;
    ctx.strokeStyle = "red";

    if( a[ 0 ] > -1 ){ ctx.strokeText( a[ 0 ].toFixed(2)/*+"-"+n[ 0 ]*/, to1x, to1y );  }
    if( a[ 1 ] >  -1 ){ ctx.strokeText( a[ 1 ].toFixed(2)/*+"-"+n[ 1 ]*/, to2x, to2y );  }
    if( a[ 2 ] >  -1 ){ ctx.strokeText( a[ 2 ].toFixed(2)/*+"-"+n[ 2 ]*/, to3x, to3y );  }

    ctx.strokeStyle = "black";

    let centroid = tri.centroid;
    ctx.beginPath();
    ctx.moveTo( centroid.x * scale, centroid.y * scale );
    ctx.arc( centroid.x * scale, centroid.y * scale, 5, 0, 2*Math.PI)
  //  ctx.fill();


  })
}

function fillTriangles( _points, _triangulation, _cv, _scale, _color ){
  let ctx = _cv.getContext( "2d" ),
      scale = _scale || 1,
      nbTri = _triangulation.length;
  ctx.clearRect( 0, 0, _cv.width, _cv.height );

  _triangulation.forEach( tri =>{
    let p = tri.positions,
        c = tri.constrained;
    ctx.fillStyle = _color;
    ctx.beginPath();
    ctx.moveTo( p[ 0 ] * scale,  p[ 1 ] * scale );
    ctx.lineTo( p[ 2 ] * scale,  p[ 3 ] * scale );

    ctx.lineTo( p[ 4 ] * scale,  p[ 5 ] * scale );
    ctx.closePath();
    ctx.fill();
  })

  ctx.beginPath();
  let crossSize = 3;
  ctx.strokeStyle = "rgb(44, 111, 0)";
  ctx.moveTo( _points.start.x - crossSize, _points.start.y - crossSize );
  ctx.lineTo( _points.start.x + crossSize, _points.start.y + crossSize );
  ctx.moveTo( _points.start.x + crossSize, _points.start.y - crossSize );
  ctx.lineTo( _points.start.x - crossSize, _points.start.y + crossSize );
  ctx.stroke();
  ctx.beginPath();
  ctx.strokeStyle = "#e61818";
  ctx.moveTo( _points.end.x - crossSize, _points.end.y - crossSize );
  ctx.lineTo( _points.end.x + crossSize, _points.end.y + crossSize );
  ctx.moveTo( _points.end.x + crossSize, _points.end.y - crossSize );
  ctx.lineTo( _points.end.x - crossSize, _points.end.y + crossSize );
  ctx.stroke();
  ctx.strokeStyle = "black";

}

function drawWidths( _triangulation, _indices, _cv, _scale ){
  let ctx = _cv.getContext( "2d" ),
      scale = _scale || 1;
  _indices.forEach( indice =>{
    let tri = _triangulation[ indice.triIndex ],
        p = tri.positions,
        pId = indice.edgeIndex < 2 ?  indice.edgeIndex * 2 + 2 : 0,
        px = p[ pId ] * scale,
        py = p[ pId + 1 ] * scale,
        width = tri.widths[ indice.edgeIndex ] * scale;
    ctx.strokeStyle = "green";
    ctx.beginPath();
    ctx.moveTo( px + width, py );
    ctx.arc( px, py, width, 0, 2*Math.PI)
    ctx.stroke();
    ctx.beginPath();
    ctx.moveTo( px + width, py );
    ctx.arc( px, py, 10, 0, 2*Math.PI);
    ctx.fill();

  } )
}

function drawAngle( x1, y1, x2, y2, x3, y3, x4, y4, _cv, _scale ){
  let ctx = _cv.getContext( "2d" ),
      scale = _scale || 1;
  ctx.beginPath();
  ctx.moveTo( x1 * scale, y1 * scale );
  ctx.lineTo( x2 * scale, y2 * scale );
  ctx.moveTo( x3 * scale, y3 * scale );
  ctx.lineTo( x4 * scale, y4 * scale );
  ctx.stroke();
  let angle = td.testAngle( x1, y1, x2, y2, x3, y3, x4, y4 );
  ctx.beginPath();
  ctx.moveTo( x2 * scale + 20, y2 * scale );
  let angle0 = td.testAngle( x1, y1, x2, y2, x1, y1, x1 +10, y1  );
  //ctx.moveTo( x2 * scale + 20, y2 * scale );
  ctx.arc( x2 * scale, y2 * scale, 20, angle0 + Math.PI, angle + Math.PI);
  ctx.stroke();
}

function drawTriPoints( _triangulation, _cv, _scale ){
  console.log( "draw" );
  let ctx = _cv.getContext( "2d" ),
      scale = _scale || 1,
      nbPoints = _triangulation.length / 2,
      nbTris = nbPoints / 3;
  for( let i = 0; i < nbTris; i ++ ){
    let triIdx = i * 6,
        p0x = _triangulation[ triIdx ],
        p0y = _triangulation[ triIdx + 1 ],
        p1x = _triangulation[ triIdx + 2 ],
        p1y = _triangulation[ triIdx + 3 ],
        p2x =_triangulation[ triIdx + 4 ],
        p2y =_triangulation[ triIdx + 5 ];
    if( !p0x || ! p0y || !p1x || !p1y || !p2x || !p2y ) continue;
    ctx.beginPath();
    ctx.moveTo( p0x * scale, p0y * scale );
    ctx.lineTo( p1x * scale, p1y * scale );
    ctx.lineTo( p2x * scale, p2y * scale );
    ctx.lineTo( p0x * scale, p0y * scale );
    ctx.strokeStyle = "black";
    ctx.stroke();
  }
}

function drawClip( _arr, _cv, _scale ){
  let ctx = _cv.getContext( "2d" ),
      scale = _scale || 1,
      nbPts = _arr.length / 2;
  ctx.beginPath();
  ctx.moveTo( _arr[ 0 ] * scale,  _arr[ 1 ] * scale );
  for( let i = 1; i < nbPts; i++ ){
    ctx.lineTo( _arr[ i * 2 ] * scale,  _arr[ i * 2 + 1 ] * scale );
  }
  ctx.closePath();
  ctx.stroke();
}
function drawP2T( _arr, _cv, _scale ){
  let ctx = _cv.getContext( "2d" ),
      scale = _scale || 1,
      nbPts = _arr.length / 2;
  _arr.forEach( tri =>{
    ctx.beginPath();
    ctx.moveTo( tri[ 0 ] * scale,  tri[ 1 ] * scale );
    ctx.lineTo( tri[ 2 ] * scale,  tri[ 3 ] * scale );
    ctx.lineTo( tri[ 4 ] * scale,  tri[ 5 ] * scale );
    ctx.closePath();
    ctx.stroke();
  } )
}

function stringIndentifyObject( ob, _tabLvl ){
  let str = "";
  let tabLvl = _tabLvl || 0;
  let tabStr = "";
  for( let t = 0; t < tabLvl; t++ ){
    tabStr += "\t";
  }
  tabLvl++;
  if( Array.isArray( ob ) && !Array.isArray( ob[ 0 ] ) && !isObject( ob[ 0 ] ) ){
    str += JSON.stringify( ob );
  }else if( Array.isArray( ob ) && !Array.isArray( ob[ 0 ] ) && isObject( ob[ 0 ] ) ){
    str += tabStr + "[\n";
    ob.forEach( arrOb => {
      str += tabStr + stringIndentifyObject( arrOb, tabLvl );
    })
    str += "\n" + tabStr + "]";
  }else if( Array.isArray( ob ) && Array.isArray( ob[ 0 ] ) ){
    str += tabStr + "[\n";
    ob.forEach( arrEl => {
      str += tabStr + stringIndentifyObject( arrEl, tabLvl );
    })
    str += "\n" + tabStr + "]";
  }else if( !Array.isArray( ob ) && isObject( ob ) ){
    str += tabStr + "{\n";
    let oEntries = Object.entries( ob );
    let l = oEntries.length;
    oEntries.forEach( ( entry, i ) => {
      str += tabStr + entry[ 0 ] + " : " + stringIndentifyObject( entry[ 1 ], tabLvl ) + ( i < l - 1 ? "\n" : "");
    })
    str += "\n" + tabStr + "}";
  }else{
    str += JSON.stringify( ob );
  }
  return str;
}
function isObject( ob ){
  for( let k in ob ){
    if( ob.hasOwnProperty( k ) ) return true;
  }
  return false;
}
