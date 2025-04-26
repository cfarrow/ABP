#include "LatticeRegular4d.h"

/* Check the x=0 and x=length-1 planes for sites on the largest cluster */
bool LatticeRegular4d::isSpanning(size_t a)
{
	if( a > 1) { a = 1; }
	size_t site;
	bool top, bottom;

    size_t maxMass = getMaxMass();
	
	top = bottom = 0;
	
     for( size_t i=0; i < length; i++) {
        if(top&&bottom) break;
        for( size_t j=0; j < length; j++) {
            if(top&&bottom) break;
            for( size_t k=0; j < length; j++) {
                if(top&&bottom) break;
                site = k*length*length*length + j*length*length + i*length;
                if( isActive(site) == a && Clusters[ cluster_label[site] ] == maxMass ) 
                    top=1;

                site = k*length*length*length + j*length*length + (i+1)*length - 1;
                if( isActive(site) == a && Clusters[ cluster_label[site] ] == maxMass ) 
                    bottom=1;

                }
            }
         }

	return (top&&bottom);
}

