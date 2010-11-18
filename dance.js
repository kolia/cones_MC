var svgDocument;
var xmlns="http://www.w3.org/2000/svg"
var xlinkns = "http://www.w3.org/1999/xlink"
function startup(evt) {
  O=evt.target
  svgDocument=O.ownerDocument
  O.setAttribute("onmousedown","run=!run;dt()")
  cones = svgDocument.getElementById("cones")
  itext = svgDocument.getElementById("itext")
}

colors = ['red' , 'green' , '#66FF66'] ;



run    = true ;
t      = 0 ;

function dt(){
  if (!run) return ;
  itext.lastChild.nodeValue = 'Iteration '+t;
  iter() ;
  t = t+1 ;
  if (t>=data.length) {
    t = 0 ;
    while(cones.hasChildNodes()) cones.removeChild(cones.firstChild);
  }
  window.setTimeout("dt()",20)
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
	cone.setAttributeNS(null,'id',dcone[0]+'_'+dcone[1]+'_'+dcone[3]) ;
      }
      cone.setAttribute('stroke',colors[dcone[2]-1]) ;
      cones.appendChild(cone) ;
    }
  }
}
