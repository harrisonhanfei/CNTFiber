//====================================================================================
//SOFTWARE:	Network of Eelectrically Conductive Nanocomposites (NECN)
//CODE FILE:	App_Network_3D.cpp
//OBJECTIVE:	Create conductive nanotube network in 3D separated by backbone paths, dead branches and isolated clusters
//AUTHOR:		Fei Han; Angel Mora
//E-MAIL:			fei.han@kaust.edu.sa	;	angel.mora@kaust.edu.sa
//====================================================================================

#include "App_Network_3D.h"

//Generate 3D conductive nanotube network separated by backbone paths, dead branches and isolated clusters
int App_Network_3D::Create_conductive_network_3D(Input *Init)const
{
	//Time markers for total simulation
	time_t ct0, ct1;

	//-----------------------------------------------------------------------------------------------------------------------------------------
	//Network Generation with overlapping
	ct0 = time(NULL);
	
	vector<vector<Point_3D> > cnts_points;	//define two-dimensional vector of three-dimensional points for storing the CNT network
    vector<double> cnts_radius;						//define the radius of each nanotube in the network

	hout << "-_- To generate networks with overlapping......"<<endl;
	GenNetwork *Genet = new GenNetwork;
	if(Genet->Generate_geometric_networks(Init->geom_rve, Init->cluster_geo, Init->nanotube_geo, cnts_points, cnts_radius)==0) return 0;

	//-----------------------------------------------------------------------------------------------------------------------------------------
	//Checking the angle between two segments in one nanotube (if less than PI/2, provide an alarm)
//	if(Genet->CNTs_quality_testing(cnts_points)==0) return 0;

	//-----------------------------------------------------------------------------------------------------------------------------------------
	//A new class of Tecplot_Export
	Tecplot_Export *Tecexpt = new Tecplot_Export;
	
	//-----------------------------------------------------------------------------------------------------------------------------------------
	//The geometric structure of CNT network (by threads in Tecplot)
	if(Init->geom_rve.shape=="Cuboid")
	{
		struct cuboid cub;														//Generate a cuboid for RVE
		cub.poi_min = Init->geom_rve.origin;
		cub.len_x = Init->geom_rve.len_x;
		cub.wid_y = Init->geom_rve.wid_y;
		cub.hei_z = Init->geom_rve.hei_z;
		cub.volume = cub.len_x*cub.wid_y*cub.hei_z;

		//The geometric structure of CNT network (by threads in Tecplot)
		if(Tecexpt->Export_network_threads(cub, cnts_points)==0) return 0;
	}
	else if(Init->geom_rve.shape=="Cylinder")
	{
		struct cylinder cyl;														//Generate a cuboid for RVE
		cyl.bottom_center = Init->geom_rve.origin;
		cyl.length = Init->geom_rve.cyl_length;
		cyl.radius = Init->geom_rve.cyl_radius;

		//-----------------------------------------------------------------------------------------------------------------------------------------
		//The geometric structure of CNT network (by threads in Tecplot)
//		if(Tecexpt->Export_network_threads(cyl, cnts_points)==0) return 0;
		
		//The geometric structure of CNT network (by tetrahedron meshes in Tecplot) //Attention: little parts of nanotube volumes out of the cylinder
//		if(Tecexpt->Export_cnt_network_meshes(cyl, cnts_points, cnts_radius)==0) return 0;

		//-----------------------------------------------------------------------------------------------------------------------------------------
//		struct cylinder cyl1, cyl2;														//Generate cylinders for the shell
//		cyl1.bottom_center = Init->geom_rve.origin;
//		cyl1.length = Init->geom_rve.cyl_length;
//		cyl1.radius = Init->geom_rve.cyl_hollow_rad;
//		cyl2.bottom_center = cyl1.bottom_center;
//		cyl2.length = cyl1.length;
//		cyl2.radius = Init->geom_rve.cyl_hollow_rad + Init->geom_rve.cyl_shell_thick;

		//The geometric structure of CNT network (by tetrahedron meshes in Tecplot) and the outside thick-shell (by brick meshes in Tecplot)
//		if(Tecexpt->Export_cnt_network_meshes(cyl, cyl1, cyl2, cnts_points, cnts_radius)==0) return 0;

		//-----------------------------------------------------------------------------------------------------------------------------------------
		//Seperating several parts of cnts
		vector<vector<Point_3D> > temp_assign;
		vector<vector<vector<Point_3D> > > sevcnts_pois((int)Init->geom_rve.stretch_cross_z.size()+1, temp_assign);	//define the cnt group on several sections (z direction)
		if(Genet->Seperate_several_cnts(cnts_points, Init->geom_rve.stretch_cross_z, sevcnts_pois)==0) return 0;

		//-----------------------------------------------------------------------------------------------------------------------------------------
		//Export the several zones of CNT fiber
		if(Tecexpt->Export_several_cnt_meshes(cyl, sevcnts_pois, cnts_radius)==0) return 0;

		//-----------------------------------------------------------------------------------------------------------------------------------------
		//Translate cnts and then merge them to one verctor storage
		vector<double> delta_dist((int)Init->geom_rve.stretch_dist.size()+1, 0.0);
		double sum_dist = 0.0;
		for(int i=0; i<(int)Init->geom_rve.stretch_dist.size(); i++) sum_dist += Init->geom_rve.stretch_dist[i];
		double temp_dist = 0.0;
		for(int i=0; i<=(int)Init->geom_rve.stretch_dist.size(); i++)
		{
			delta_dist[i] = (temp_dist-0.5*sum_dist)/Init->geom_rve.stretch_steps;
			if(i<(int)Init->geom_rve.stretch_dist.size()) temp_dist += Init->geom_rve.stretch_dist[i];
		}

		for(int i=0; i<=Init->geom_rve.stretch_steps; i++)
		{
			//Translate several parts of cnts
			if(i>0) if(Genet->Translation_several_cnts(delta_dist, sevcnts_pois)==0) return 0;
			
			vector<vector<Point_3D> > temp_cnts;
			for(int j=0; j<(int)sevcnts_pois.size(); j++)
				for(int k=0; k<(int)sevcnts_pois[j].size(); k++) 
					temp_cnts.push_back(sevcnts_pois[j][k]);

			//-----------------------------------------------------------------------------------
			//The several CNT networks (by threads in Tecplot)
//			stringstream str_thread;
//			if(i<10)	str_thread << "CNT_Fiber_Multiple_Wires_000" << i << ".dat";
//			else if (i<100) str_thread << "CNT_Fiber_Multiple_Wires_00" << i << ".dat";
//			else if (i<1000)	str_thread << "CNT_Fiber_Multiple_Wires_0" << i << ".dat";
//			else str_thread << "CNT_Fiber_Multiple_Wires_" << i << ".dat";
//			if(Tecexpt->Export_several_threads(str_thread.str(), cyl, temp_cnts)==0) return 0;
		
			//The several CNT networks (by tetrahedron meshes in Tecplot) //Attention: little parts of nanotube volumes out of the cylinder
			stringstream str_mesh;
			if(i<10)	str_mesh << "CNT_Fiber_Multiple_Meshes_000" << i << ".dat";
			else if (i<100) str_mesh << "CNT_Fiber_Multiple_Meshes_00" << i << ".dat";
			else if (i<1000)	str_mesh << "CNT_Fiber_Multiple_Meshes_0" << i << ".dat";
			else str_mesh << "CNT_Fiber_Multiple_Meshes_" << i << ".dat";
			if(Tecexpt->Export_several_cnt_meshes(str_mesh.str(), cyl, temp_cnts, cnts_radius)==0) return 0;
		}
	}

	//-----------------------------------------------------------------------------------------------------------------------------------------
	//输出纳米管网络几何(Template for Sofiane)
//	if(Export_cnt_fibers(cnts_points)==0) return 0;

	//Generate the nodes and quadrilateral elements for CNT constructions
//	vector<vector<Node> > cnts_nodes;
//	vector<vector<Element> > cnts_eles;

//	if(Generate_cnts_nodes_elements(cnts_nodes, cnts_eles, cnts_points)==0) return 0;

	//-----------------------------------------------------------------------------------------------------------------------------------------
	//输出纳米管网络几何构造(四面体网格Tecplot) //注释在RVE表面处还没有考虑表面要切掉纳米管的一部分
//	if(Export_cnt_networks_meshes(cell_geo, cnts_points, ellips)==0) return 0;

	//-----------------------------------------------------------------------------------------------------------------------------------------
	//记录纳米管点信息(Point的flag==0表示纳米管的起始点; flag>0表示纳米管内的点, 按序依次编号)
//	if(Record_cnt_points_information(cnts_points)==0) return 0;

	//-----------------------------------------------------------------------------------------------------------------------------------------
	//Delete class points
	delete Genet;
	delete Tecexpt;
	//-----------------------------------------------------------------------------------------------------------------------------------------
	ct1 = time(NULL);
	hout << "Network generation time: "<<(int)(ct1-ct0)<<" secs."<<endl;
	hout << "^_^ End of network generation with overlapping."<<endl<<endl;

	return 1;
}
//===========================================================================
