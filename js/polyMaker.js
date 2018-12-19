const polyMaker = {};
polyMaker.createUI = ( parent, canvas ) => {
  let div = parent.appendChild( document.createElement( "div" ) );
  div.id = "polyMakerCont";
  div.style.position = "absolute";
  div.style.left = "1000px";
  div.style.right = 0;
  div.style.top = 0;
  div.style.bottom = 0;
  div.style.backgroundColor = "#e0cdb8";

  polyMaker.canvas = canvas;

  let shapesCont = polyMaker.addUiContainer( "polyMakerShapesCont" ),
      butCont = polyMaker.addUiContainer( "polyMakerButCont" ),
      stackCont = polyMaker.addUiContainer( "polyMakerStackCont" ),
      JSONCont = polyMaker.addUiContainer( "polyMakerJSONCont" ),
      stackLineCont = polyMaker.addUiContainer( "polyMakerStackLineCont", document.getElementById( "polyMakerStackCont" ) ),
      stackPointCont = polyMaker.addUiContainer( "polyMakerStackPointCont", document.getElementById( "polyMakerStackCont" ) ),
      startOutline = polyMaker.addButton( "Start outline", polyMaker.startOutline ),
      startObstacle = polyMaker.addButton( "Start obstacle", polyMaker.startObstacle ),
      closeLine = polyMaker.addButton( "close line", polyMaker.closeLine ),
      toJSON = polyMaker.addButton( "to JSON", polyMaker.toJSON );
      setArrayBut = polyMaker.addButton( "setArray", polyMaker.setArray );
  stackCont.classList.add( "frwsss" );
  stackLineCont.style.width = "45%";
  stackPointCont.style.width = "45%";
  polys.forEach( ( polyArray, idx ) => {
    polyMaker.createStackLine( shapesCont, "ployMakerShape_" + idx, "Poly " + idx, polyMaker.loadFromArray );
  })
}
polyMaker.addUiContainer = ( id, parent ) => {
  let _parent = parent || document.getElementById( "polyMakerCont" ),
      butCont = _parent.appendChild( document.createElement( "div" ) );
  butCont.id = id;
  butCont.style.position = "relative";
  butCont.style.margin = "5px";
  butCont.style.padding = "5px";
  butCont.style.border = "1px solid #8b9897";
  return butCont;
}
polyMaker.addButton = ( txt, callback ) => {
  let but = document.createElement( "button" );
  but.textContent = txt;
  document.getElementById( "polyMakerButCont" ).appendChild( but );
  but.addEventListener( "click", callback, false );
  return but;
}
polyMaker.startOutline = () => {
  if( !document.getElementById( "stackPolyLine" ) ){
    let polyLineDiv = polyMaker.createStackLine( document.getElementById( "polyMakerStackLineCont" ), "stackPolyLine", "Outline", polyMaker.clickPolyLine, true ),
        polyOb = {
          id: "stackPolyLine",
          points: []
        }
    polyMaker.polyLines.push( polyOb );
  }
  let polyOb = polyMaker.polyLines.find( ob => ob.id == "stackPolyLine" );
  polyMaker.activePoly = polyOb;
  polyMaker.startLine();
};
polyMaker.startObstacle = () => {
  let obstacle = {  id: "stackObstacleLine" + polyMaker.obstaclesCount,
                    points: []
                };
  polyMaker.obstaclesCount++
  polyMaker.polyLines.push( obstacle );
  polyMaker.activePoly = obstacle;
  let obstacleLineDiv = polyMaker.createStackLine( document.getElementById( "polyMakerStackLineCont" ), obstacle.id, "Obstacle : " + obstacle.id, polyMaker.clickPolyLine );
  polyMaker.startLine();
};
polyMaker.startLine = () => {
  polyMaker.updatePointsStack();
  polyMaker.canvas.addEventListener( "mousedown", polyMaker.startDragPoint, false );
  //polyMaker.startDragPoint()
};
polyMaker.closeLine = () => {
  polyMaker.activePoly.closed = true;
  polyMaker.activePoly = false;

  let cv = polyMaker.canvas;
  cv.removeEventListener( "mousemove", polyMaker.dragPoint );
  //cv.removeEventListener( "mouseup", polyMaker.endDragPoint );
  cv.removeEventListener( "mousedown", polyMaker.startDragPoint );
  cv.removeEventListener( "click", polyMaker.addPoint );
  polyMaker.drawAll();
};
polyMaker.createStackLine = ( parent, id, name, clickFunc, first ) => {
  let div = document.createElement( "div" );
  if( first ){
    parent.insertBefore( div, parent.childNodes[ 0 ] );
  }else{
    parent.appendChild( div );
  }
  div.id = id;
  div.style.width = "100%";
  div.style.height = "15px";
  div.style.border = "1px solid #eaece8";
  let nameCont = div.appendChild( document.createElement( "span" ) );
  nameCont.textContent = name;
  if( clickFunc ){
    div.addEventListener( "click", clickFunc )
  }

}
polyMaker.clickPolyLine = e => {
  let id = e.currentTarget.id,
      polyOb = polyMaker.polyLines.find( ob => ob.id == id );
  polyMaker.activePoly = polyOb;
  polyMaker.startLine();
}
polyMaker.updatePointsStack = () => {
  let pointsStack = document.getElementById( "polyMakerStackPointCont" );
  pointsStack.innerHTML = "";
  polyMaker.activePoly.points.forEach( ( pt, idx ) => {
    polyMaker.createStackLine( pointsStack, polyMaker.activePoly.id + "_pt_" + idx, "x : " + pt[ 0 ] + " / y : " + pt[ 1 ] );
  })
}
polyMaker.startDragPoint = e => {
  let cv = polyMaker.canvas;
  /*
  cv.removeEventListener( "mousedown", polyMaker.startDragPoint );
  cv.addEventListener( "mousemove", polyMaker.dragPoint );
  cv.addEventListener( "mouseup", polyMaker.endDragPoint );
  */
  cv.removeEventListener( "mousedown", polyMaker.startDragPoint );
  cv.addEventListener( "mousemove", polyMaker.dragPoint );
  cv.addEventListener( "click", polyMaker.addPoint );
  //cv.addEventListener( "mousedown", polyMaker.addPoint );
}
polyMaker.dragPoint = e => {
  polyMaker.memPt = { x: e.clientX, y:  e.clientY };
  window.requestAnimationFrame( ()=>{
    let cv = polyMaker.canvas;
    if( polyMaker.activePoly.points.length > 0 ){
      polyMaker.drawAll();
      let lastPoint = polyMaker.activePoly.points[ polyMaker.activePoly.points.length - 1 ],
          ctx = cv.getContext( "2d" ),
          _scale = polyMaker.scale || 1;
      ctx.beginPath();
      ctx.moveTo( lastPoint[ 0 ] * _scale, lastPoint[ 1 ] * _scale );
      ctx.lineTo( polyMaker.memPt.x, polyMaker.memPt.y );
      ctx.stroke();
      console.log("drag")
    }
  })
}
polyMaker.addPoint = e => {
  let cv = e.currentTarget;
  polyMaker.activePoly.points.push( [ e.clientX, e.clientY ] );
  polyMaker.updatePointsStack();
  polyMaker.drawAll();
}
polyMaker.endDragPoint = e => {
  let cv = e.currentTarget;
  cv.removeEventListener( "mousemove", polyMaker.dragPoint );
  cv.removeEventListener( "mouseup", polyMaker.endDragPoint );
  cv.addEventListener( "mousedown", polyMaker.startDragPoint );
  polyMaker.activePoly.points.push( [ e.clientX, e.clientY ] );
  polyMaker.updatePointsStack();
  polyMaker.drawAll();
}
polyMaker.drawAll = () => {
  let ctx = polyMaker.canvas.getContext( "2d" );
  ctx.clearRect( 0, 0, polyMaker.canvas.width, polyMaker.canvas.height );
  polyMaker.polyLines.forEach( polyOb => { polyMaker.drawPoly( polyOb ) } );
}
polyMaker.drawPoly = ob => {
  let pts = ob.points,
      ctx = polyMaker.canvas.getContext( "2d" ),
      _scale = polyMaker.scale || 1;
  ctx.beginPath();
  ctx.moveTo( pts[ 0 ][ 0 ] * _scale, pts[ 0 ][ 1 ] * _scale );
  for( let i = 1; i < pts.length; i++ ){
    ctx.lineTo( pts[ i ][ 0 ] * _scale, pts[ i ][ 1 ] * _scale );
  }
  if( ob.closed ) ctx.closePath();
  ctx.stroke();
}
polyMaker.toArray = () => {
  let ret = [];
  polyMaker.polyLines.forEach( polyOb => {
    ret.push( polyOb.points );
  })
  return ret;
}
polyMaker.toJSON = () => {
  document.getElementById( "polyMakerJSONCont" ).innerHTML = JSON.stringify( polyMaker.toArray() );
}
polyMaker.setArray = () => {
  polyMaker.array = polyMaker.toArray();
}
polyMaker.loadFromArray  = e => {
  let idx = parseFloat( e.currentTarget.id.replace("ployMakerShape_","") ),
      array = polys[ idx ];
  polyMaker.memPt = { x: null, y: null };
  polyMaker.polyLines = [];
  polyMaker.obstaclesCount = 0;
  polyMaker.activePoly = false;
  polyMaker.array = false;
  polyMaker.polyLines.push({
    id: "stackPolyLine",
    points: array[ 0 ],
    closed: true
  })
  for(let i = 1; i < array.length; i++ ){
    polyMaker.polyLines.push({
      id: "stackObstacleLine" + i - 1,
      points: array[ i ],
      closed: true
    })
    polyMaker.obstaclesCount++;
  }
  polyMaker.setArray();
  polyMaker.drawAll();
}

polyMaker.memPt = { x: null, y: null };
polyMaker.polyLines = [];
polyMaker.obstaclesCount = 0;
polyMaker.activePoly = false;
polyMaker.array = false;
