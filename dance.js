var svgDocument;
var xmlns="http://www.w3.org/2000/svg"
var xlinkns = "http://www.w3.org/1999/xlink"
function startup(evt) {
  O=evt.target
  svgDocument=O.ownerDocument
  O.setAttribute("onmousedown","run=!run;dt()")
  greedy_cones = svgDocument.getElementById("greedy_cones")
  cones = svgDocument.getElementById("cones")
  itext    = svgDocument.getElementById("itext")
}

colors = ['red' , 'green' , 'blue'] ;

%s

run    = true ;
t      = 0 ;

function dt(){
  if (!run) return ;
  itext.lastChild.nodeValue = 'Iteration '+t;
  iter() ;
  t = t+1 ;
  if (t>=data.length) {
    t = 0 ;
    for (i=0; i<4; i++){
      these_cones = svgDocument.getElementById('cones'+i) ;
      while(these_cones.hasChildNodes()) these_cones.removeChild(these_cones.firstChild);
    }
  }
  window.setTimeout("dt()",2)
}
function iter(){
  for (i=0; i<data[t].length; i++){
    dcone = data[t][i] ;
    cone  = svgDocument.getElementById(dcone[0]+'_'+dcone[1]+'_'+dcone[3]) ;
    if (!dcone[2] && cone) {
      cone.parentNode.removeChild(cone) ;
    }
    if ( dcone[2]) {
      if (!cone) {
	cone = svgDocument.createElementNS(xmlns,'use') ;
	cone.setAttributeNS(xlinkns,'xlink:href','#m'+dcone[3]) ;
	cone.setAttributeNS(null,'transform','translate('+dcone[0]+' '+dcone[1]+')') ;
	cone.setAttributeNS(null,'class','m'+dcone[3]) ;
	cone.setAttributeNS(null,'id',dcone[0]+'_'+dcone[1]+'_'+dcone[3]) ;
      }
      cone.setAttribute('stroke',colors[dcone[2]-1]) ;
      cone.setAttribute('fill','none') ;
      these_cones = svgDocument.getElementById('cones'+dcone[3]) ;
      these_cones.appendChild(cone) ;
    }
  }
}
function toggle(id){
  run = !run ;
  instance = svgDocument.getElementById('cones'+id[1]) ;
  if (instance.getAttribute("visibility")=="visible"){
    instance.setAttributeNS(null,'visibility','hidden') ;
  }
  else {
    instance.setAttributeNS(null,'visibility','visible') ;
  }
}
