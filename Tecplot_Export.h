//====================================================================================
//SOFTWARE:	Network of Eelectrically Conductive Nanocomposites (NECN)
//CODE FILE:	Tecplot_Export.h
//OBJECTIVE:	To export the 3D geometric images through Tecplot data files 
//AUTHOR:		Fei Han; Angel Mora
//E-MAIL:			fei.han@kaust.edu.sa	;	angel.mora@kaust.edu.sa
//====================================================================================

#ifndef  TECPLOTEXPORT_H
#define TECPLOTEXPORT_H

#include "Geometry_3D.h"
#include "Fem_3D.h"
#include "Input_Reader.h"
#include "GenNetwork.h"

//-------------------------------------------------------
class Tecplot_Export
{
	public:
		//Data Member
		
		//Constructor
		Tecplot_Export(){};

		//Member Functions
		//The geometric structure of CNT network (by quadrilaterial elements in Tecplot)
		int Export_cnt_network_meshes(const struct cylinder &cyl, const vector<vector<Point_3D> > &cnts_points, const vector<double> &cnts_radius)const;
		//The geometric structure of CNT network (by tetrahedron meshes in Tecplot) and the outside thick-shell (by brick meshes in Tecplot)
		int Export_cnt_network_meshes(const struct cylinder &cyl, const struct cylinder &cyl1, const struct cylinder &cyl2, const vector<vector<Point_3D> > &cnts_points, 
														  const vector<double> &cnts_radius)const;
		//The top and bottom of CNT networks (by tetrahedron meshes in Tecplot) //Attention: little parts of nanotube volumes out of the cylinder
		int Export_top_bottom_cnt_meshes(const string &output_file_name, const struct cylinder &cyl, const vector<vector<Point_3D> > &cnts_toppois, const vector<vector<Point_3D> > &cnts_botpois, const vector<double> &cnts_radius)const;
		int Export_top_bottom_cnt_meshes(const string &output_file_name, const struct cylinder &cyl, const vector<vector<Point_3D> > &cnts_points, const vector<double> &cnts_radius)const;
		//The geometric structure of CNT network (by threads in Tecplot) in a cuboid
		int Export_network_threads(const struct cuboid &cub, const vector<vector<Point_3D> > &cnts_points)const;
		//The geometric structure of CNT network (by threads in Tecplot) in a cylinder
		int Export_network_threads(const struct cylinder &cyl, const vector<vector<Point_3D> > &cnts_points)const;
		//The top and bottom of CNT networks (by threads in Tecplot)
		int Export_top_bottom_threads(const string &output_file_name, const struct cylinder &cyl, const vector<vector<Point_3D> > &cnts_toppois, const vector<vector<Point_3D> > &cnts_botpois)const;
		int Export_top_bottom_threads(const string &output_file_name, const struct cylinder &cyl, const vector<vector<Point_3D> > &cnts_points)const;

	private:
		//Export a 3D cuboid
		int Export_cuboid(ofstream &otec, const struct cuboid &cub)const;
		//Export a 3D cylinder
		int Export_cylinder(ofstream &otec, const struct cylinder &cell)const;
		//Export 3D nanotube threads
		int Export_nano_threads(ofstream &otec, const vector<vector<Point_3D> > &cnts_points)const;
		//输出纳米管线网格多Zones in Tecplot
		int Export_cnts_meshes_multizones(ofstream &otec, const struct cylinder &cyl, const vector<vector<Node> > &nodes, const vector<vector<Element> > &eles)const;
		//输出纳米管线网格单Zone in Tecplot
		int Export_cnts_meshes_singlezone(ofstream &otec, const struct cylinder &cyl, const vector<vector<Node> > &nodes, const vector<vector<Element> > &eles)const;
		//Export the thick shell of the cylinder
		int Export_cylinder_shell(ofstream &otec, const struct cylinder &cyl1, const struct cylinder &cyl2)const;
};
//-------------------------------------------------------
#endif
//===========================================================================