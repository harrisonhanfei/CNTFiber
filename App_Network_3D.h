//====================================================================================
//SOFTWARE:	Network of Eelectrically Conductive Nanocomposites (NECN)
//CODE FILE:	App_Network_3D.h
//OBJECTIVE:	Create a 3D nanotube netwok
//AUTHOR:		Fei Han; Angel Mora
//E-MAIL:			fei.han@kaust.edu.sa	;	angel.mora@kaust.edu.sa
//====================================================================================

#ifndef APPNETWORK3D_H
#define APPNETWORK3D_H

#include "time.h"
#include "Hns.h"
using namespace hns;

#include "Input_Reader.h"
#include "GenNetwork.h"
#include "Tecplot_Export.h"

//---------------------------------------------------------------------------
class App_Network_3D
{
	public:
		//Data Member
		vector<Point_3D> cnps;			//Define 3D point verctor of nanotuber points
		vector<double> cnts_radius;		//Define the radius of every nanotube in the network

		//Constructor
		App_Network_3D(){};

		//Member Functions
		int Create_conductive_network_3D(Input *Init)const;
};
//---------------------------------------------------------------------------
#endif
//===========================================================================

