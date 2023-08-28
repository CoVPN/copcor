# extract assay name from marker names, which include Day, e.g.
marker.name.to.assay=function(a) {
  
  if (startsWith(a,"Day")) {
    # Day22pseudoneutid50 => pseudoneutid50 
    sub("Day[[0123456789]+", "", a)
    
  } else if (startsWith(a,"BD")) {
    # BD29pseudoneutid50 => pseudoneutid50
    sub("BD[[0123456789]+", "", a)
    
  } else if (contain(a,"overBD1")) {
    # DeltaBD29overBD1pseudoneutid50 => pseudoneutid50
    sub("DeltaBD[[0123456789]+overBD1", "", a)    
    
  } else if (contain(a,"overB")) {
    # Delta22overBpseudoneutid50 => pseudoneutid50
    sub("Delta[[0123456789]+overB", "", a)
    
  } else if (contain(a,"over")) {
    sub("Delta[[0123456789]+over[[0123456789]+", "", a)
    
  } else stop("marker.name.to.assay: not sure what to do")
}


