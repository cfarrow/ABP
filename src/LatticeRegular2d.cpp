#include "LatticeRegular2d.h"

/* Check the x=0 and x=length-1 rows for sites on the largest cluster */
bool LatticeRegular2d::isSpanning(size_t a) {
	if( a > 1) { a = 1; }
	size_t site;
	bool top, bottom;
	
	top = bottom = 0;
	
     for( size_t i=0; i < length; i++ ) {
         if(top&&bottom) break;
         site = i*length;
         if( isActive(site) == a && Clusters[ cluster_label[site] ] == getMaxMass() ) 
             top = 1;
         site = (i+1)*length - 1;
         if( isActive(site) == a && Clusters[ cluster_label[site] ] == getMaxMass() ) 
             bottom=1;
     }

     return (top&&bottom);
}

