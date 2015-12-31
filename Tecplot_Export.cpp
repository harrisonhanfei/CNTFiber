//====================================================================================
//SOFTWARE:	Network of Eelectrically Conductive Nanocomposites (NECN)
//CODE FILE:	Tecplot_Export.h
//OBJECTIVE:	To export the 3D geometric images through Tecplot data files
//AUTHOR:		Fei Han; Angel Mora
//E-MAIL:			fei.han@kaust.edu.sa	;	angel.mora@kaust.edu.sa
//====================================================================================

#include "Tecplot_Export.h"

//---------------------------------------------------------------------------
//The geometric structure of CNT network (by threads in Tecplot) in a cuboid
int Tecplot_Export::Export_network_threads(const struct cuboid &cub, const vector<vector<Point_3D> > &cnts_points)const
{
	ofstream otec("CNT_Wires.dat");
	otec << "TITLE = CNT_Wires" << endl;
	otec << "VARIABLES = X, Y, Z" << endl;

	//---------------------------------------------------------------------------
	//Export a 3D cuboid
	if(Export_cuboid(otec, cub)==0) return 0;
	
	//---------------------------------------------------------------------------
	//Export 3D nanotube threads
	if(Export_nano_threads(otec, cnts_points)==0) return 0;

	//---------------------------------------------------------------------------
	otec.close();

	return 1;
}
//---------------------------------------------------------------------------
//The geometric structure of CNT network (by threads in Tecplot) in a cylinder
int Tecplot_Export::Export_network_threads(const struct cylinder &cyl, const vector<vector<Point_3D> > &cnts_points)const
{
	ofstream otec("CNT_Wires.dat");
	otec << "TITLE = CNT_Wires" << endl;
	otec << "VARIABLES = X, Y, Z" << endl;

	//---------------------------------------------------------------------------
	//Export a 3D cylinder
	if(Export_cylinder(otec, cyl)==0) return 0;
	
	//---------------------------------------------------------------------------
	//Export 3D nanotube threads
	if(Export_nano_threads(otec, cnts_points)==0) return 0;

	//---------------------------------------------------------------------------
	otec.close();

	return 1;
}
//---------------------------------------------------------------------------
//Export a 3D cuboid
int Tecplot_Export::Export_cuboid(ofstream &otec, const struct cuboid &cell)const
{
	otec << "ZONE N=" << 8 << ", E=" << 1 << ", F=FEPOINT, ET=BRICK" << endl;
	double cell_x[2] = {cell.poi_min.x, cell.poi_min.x+cell.len_x};
	double cell_y[2] = {cell.poi_min.y, cell.poi_min.y+cell.wid_y};
	double cell_z[2] = {cell.poi_min.z, cell.poi_min.z+cell.hei_z};
	for(int i=0; i<2; i++)
		for(int j=0; j<2; j++)
			for(int k=0; k<2; k++)
			{
				otec << cell_x[i] << "  " << cell_y[j] << "  " << cell_z[k] << endl;
			}
	
	otec << "1 2 4 3 5 6 8 7" << endl;
	otec << endl << endl;

	return 1;
}
//---------------------------------------------------------------------------
//Export a 3D cylinder
int Tecplot_Export::Export_cylinder(ofstream &otec, const struct cylinder &cyl)const
{
	//Define the number of nodes in top or bottom circles of the cylinder
	const int num_sec = 360;	

	vector<Node> nod_temp;
	vector<Element> ele_temp;

	//Define a class of GenNetwork for calling for the function below
	GenNetwork *Gentemp = new GenNetwork;
	Gentemp->Generate_cylinder_nodes_elements(cyl, num_sec, nod_temp, ele_temp);
	delete Gentemp;

	//---------------------------------------------------------------------------
	//Export the mesh of the cylinder
	otec << "ZONE N=" << (int)nod_temp.size() << ", E=" << (int)ele_temp.size() << ", F=FEPOINT, ET=TETRAHEDRON" << endl;
	for(int i=0; i<(int)nod_temp.size(); i++)
	{
		otec << nod_temp[i].x << "  " << nod_temp[i].y << "  " << nod_temp[i].z << endl;
	}
	otec << endl;
	for(int i=0; i<(int)ele_temp.size(); i++)
	{
		otec	<< ele_temp[i].nodes_id[0]+1 << "  " << ele_temp[i].nodes_id[1]+1 << "  " 
				<< ele_temp[i].nodes_id[2]+1 << "  " << ele_temp[i].nodes_id[3]+1 << endl;
	}
	otec	<< endl;

	otec << "ZONE N=" << 8 << ", E=" << 2 << ", F=FEPOINT, ET=BRICK" << endl;
	otec << nod_temp[num_sec/3+1].x << "  " << nod_temp[num_sec/3+1].y << "  " << nod_temp[num_sec/3+1].z << endl;
	otec << nod_temp[num_sec+1+num_sec/3+1].x << "  " << nod_temp[num_sec+1+num_sec/3+1].y << "  " << nod_temp[num_sec+1+num_sec/3+1].z << endl;
	otec << nod_temp[num_sec/3+1].x - 1E-6 << "  " << nod_temp[num_sec/3+1].y - 1E-6 << "  " << nod_temp[num_sec/3+1].z << endl;
	otec << nod_temp[num_sec+1+num_sec/3+1].x - 1E-6 << "  " << nod_temp[num_sec+1+num_sec/3+1].y - 1E-6 << "  " << nod_temp[num_sec+1+num_sec/3+1].z << endl;
	otec << nod_temp[5*num_sec/6+1].x << "  " << nod_temp[5*num_sec/6+1].y << "  " << nod_temp[5*num_sec/6+1].z << endl;
	otec << nod_temp[num_sec+1+5*num_sec/6+1].x << "  " << nod_temp[num_sec+1+5*num_sec/6+1].y << "  " << nod_temp[num_sec+1+5*num_sec/6+1].z << endl;
	otec << nod_temp[5*num_sec/6+1].x  + 1E-6 << "  " << nod_temp[5*num_sec/6+1].y  + 1E-6 << "  " << nod_temp[5*num_sec/6+1].z << endl;
	otec << nod_temp[num_sec+1+5*num_sec/6+1].x  + 1E-6 << "  " << nod_temp[num_sec+1+5*num_sec/6+1].y  + 1E-6 << "  " << nod_temp[num_sec+1+5*num_sec/6+1].z << endl;
	otec << "1 2 4 3 1 2 4 3" << endl;
	otec << "5 6 8 7 5 6 8 7" << endl;
	otec << endl << endl;

	return 1;
}
//---------------------------------------------------------------------------
//Export 3D nanotube threads
int Tecplot_Export::Export_nano_threads(ofstream &otec, const vector<vector<Point_3D> > &cnts_points)const
{
	for(int i=0; i<(int)cnts_points.size(); i++)
	{
		otec << "ZONE T=\"Line\"" << endl;
		otec << "i=1," << "j=" << (int)cnts_points[i].size() << ", f=point" << endl;
		for (int j=0; j<(int)cnts_points[i].size(); j++)
		{
			otec << cnts_points[i][j].x << "  " << cnts_points[i][j].y << "  " << cnts_points[i][j].z << endl;
		}
		otec << endl << endl;
	}

	return 1;
}
//---------------------------------------------------------------------------
//The geometric structure of CNT network (by tetrahedron meshes in Tecplot)
int Tecplot_Export::Export_cnt_network_meshes(const struct cylinder &cyl, const vector<vector<Point_3D> > &cnts_points, const vector<double> &cnts_radius)const
{
	//生成用于表示纳米管的节点及四面体网格
	vector<vector<Node> > cnts_nodes;
	vector<vector<Element> > cnts_eles;

	//Define a class of GenNetwork for calling for the function below
	GenNetwork *Gentemp = new GenNetwork;
	if(Gentemp->Generate_cnts_nodes_elements(cnts_nodes, cnts_eles, cnts_points, cnts_radius)==0) return 0;
	delete Gentemp;

	//输出纳米管线网格多Zones in Tecplot
//	ofstream otec_multip("CNT_Meshes_Multizones.dat");
//	if(Export_cnts_meshes_multizones(otec_multip, cyl, cnts_nodes, cnts_eles)==0) return 0;
//	otec_multip.close();

	//输出纳米管线网格单Zone in Tecplot
	ofstream otec_single("CNT_Meshes_Singlezone.dat");
	if(Export_cnts_meshes_singlezone(otec_single, cyl, cnts_nodes, cnts_eles)==0) return 0;
	otec_single.close();

	return 1;
}
//---------------------------------------------------------------------------
//The geometric structure of CNT network (by tetrahedron meshes in Tecplot)
int Tecplot_Export::Export_cnt_network_meshes(const struct cylinder &cyl, const struct cylinder &cyl1, const struct cylinder &cyl2, const vector<vector<Point_3D> > &cnts_points, 
																		   const vector<double> &cnts_radius)const
{
	//生成用于表示纳米管的节点及四面体网格
	vector<vector<Node> > cnts_nodes;
	vector<vector<Element> > cnts_eles;

	//Define a class of GenNetwork for calling for the function below
	GenNetwork *Gentemp = new GenNetwork;
	if(Gentemp->Generate_cnts_nodes_elements(cnts_nodes, cnts_eles, cnts_points, cnts_radius)==0) return 0;
	delete Gentemp;

	//输出纳米管线网格多Zones in Tecplot
	ofstream otec_multip("CNT_Meshes_Multizones.dat");
	if(Export_cnts_meshes_multizones(otec_multip, cyl, cnts_nodes, cnts_eles)==0) return 0;
	otec_multip.close();

	//输出纳米管线网格单Zone in Tecplot
	ofstream otec_single("CNT_Meshes_Singlezone.dat");
	if(Export_cnts_meshes_singlezone(otec_single, cyl, cnts_nodes, cnts_eles)==0) return 0;
	otec_single.close();

	//Export the thick shell of the cylinder
	ofstream otec_shell("CNT_Meshes_Multizones.dat", ios::app);
	if(Export_cylinder_shell(otec_shell, cyl1, cyl2)==0) return 0;
	otec_shell.close();

	return 1;
}
//---------------------------------------------------------------------------
//输出纳米管线网格多Zones in Tecplot
int Tecplot_Export::Export_cnts_meshes_multizones(ofstream &otec, const struct cylinder &cyl, const vector<vector<Node> > &nodes, const vector<vector<Element> > &eles)const
{
	//纳米管线的数量
	int cnts_account = (int)nodes.size();
	if(cnts_account!=(int)eles.size()) { hout << "节点和单元显示的纳米管线数量不一致！ 请检查！" << endl; return 0; }
	
	otec << "TITLE = CNT_Meshes_Multizones" << endl;
	otec << "VARIABLES = X, Y, Z" << endl;	
	
	//---------------------------------------------------------------------------
	//Export a 3D cylinder
	if(Export_cylinder(otec, cyl)==0) return 0;

	//---------------------------------------------------------------------------
	//输出纳米管网格
	for(int i=0; i<cnts_account; i++)
	{
		otec << "ZONE N=" << (int)nodes[i].size() << ", E=" << (int)eles[i].size() << ", F=FEPOINT, ET=TETRAHEDRON" << endl;
		for(int j=0; j<(int)nodes[i].size(); j++)
		{
			otec << nodes[i][j].x << "  " << nodes[i][j].y << "  " << nodes[i][j].z << endl;
		}
		otec << endl;
		for(int j=0; j<(int)eles[i].size(); j++)
		{
			otec	<< eles[i][j].nodes_id[0]+1 << "  " << eles[i][j].nodes_id[1]+1 << "  " 
					<< eles[i][j].nodes_id[2]+1 << "  " << eles[i][j].nodes_id[3]+1 << endl;
		}
		otec << endl << endl;
	}

	return 1;
}
//---------------------------------------------------------------------------
//输出纳米管线网格单Zone in Tecplot
int Tecplot_Export::Export_cnts_meshes_singlezone(ofstream &otec, const struct cylinder &cyl, const vector<vector<Node> > &nodes, const vector<vector<Element> > &eles)const
{
	//纳米管线的数量
	int cnts_account = (int)nodes.size();
	if(cnts_account!=(int)eles.size()) { hout << "节点和单元显示的纳米管线数量不一致！ 请检查！" << endl; return 0; }
	
	otec << "TITLE = CNT_Meshes_Singlezone" << endl;
	otec << "VARIABLES = X, Y, Z" << endl;	
	
	//---------------------------------------------------------------------------
	//Export a 3D cylinder
//	if(Export_cylinder(otec, cyl)==0) return 0;

	//---------------------------------------------------------------------------
	//输出纳米管网格
	int nodes_num = 0;
	int eles_num = 0;

	for(int i=0; i<cnts_account; i++)
	{
		nodes_num +=  (int)nodes[i].size();
		eles_num += (int)eles[i].size();
	}
		
	otec << "ZONE N=" << nodes_num << ", E=" << eles_num << ", F=FEPOINT, ET=TETRAHEDRON" << endl;
	for(int i=0; i<cnts_account; i++)
	{		
		for(int j=0; j<(int)nodes[i].size(); j++)
		{
			otec << nodes[i][j].x << "  " << nodes[i][j].y << "  " << nodes[i][j].z << endl;
		}
	}
	otec << endl;

	nodes_num = 0;
	for(int i=0; i<cnts_account; i++)
	{
		if(i!=0)  nodes_num +=  (int)nodes[i-1].size();
		for(int j=0; j<(int)eles[i].size(); j++)
		{
			otec	<< eles[i][j].nodes_id[0]+1+nodes_num << "  " << eles[i][j].nodes_id[1]+1+nodes_num << "  " 
					<< eles[i][j].nodes_id[2]+1+nodes_num << "  " << eles[i][j].nodes_id[3]+1+nodes_num << endl;
		}
	}
	otec << endl << endl;

	return 1;
}
//---------------------------------------------------------------------------
//Export the thick shell of the cylinder
int Tecplot_Export::Export_cylinder_shell(ofstream &otec, const struct cylinder &cyl1, const struct cylinder &cyl2)const
{
	//Define the number of nodes in top or bottom circles of the cylinder
	const int num_sec = 360;	

	vector<Node> nod_temp;
	vector<Element> ele_temp;

	//Define a class of GenNetwork for calling for the function below
	GenNetwork *Gentemp = new GenNetwork;
	Gentemp->Generate_cylshell_nodes_elements(cyl1, cyl2, num_sec, nod_temp, ele_temp);
	delete Gentemp;

	//---------------------------------------------------------------------------
	//Export the mesh of the cylinder
	otec << "ZONE N=" << (int)nod_temp.size() << ", E=" << (int)ele_temp.size() << ", F=FEPOINT, ET=BRICK" << endl;
	for(int i=0; i<(int)nod_temp.size(); i++)
	{
		otec << nod_temp[i].x << "  " << nod_temp[i].y << "  " << nod_temp[i].z << endl;
	}
	otec << endl;
	for(int i=0; i<(int)ele_temp.size(); i++)
	{
		otec	<< ele_temp[i].nodes_id[0]+1 << "  " << ele_temp[i].nodes_id[1]+1 << "  " 
				<< ele_temp[i].nodes_id[2]+1 << "  " << ele_temp[i].nodes_id[3]+1 << endl
				<< ele_temp[i].nodes_id[4]+1 << "  " << ele_temp[i].nodes_id[5]+1 << "  " 
				<< ele_temp[i].nodes_id[6]+1 << "  " << ele_temp[i].nodes_id[7]+1 << endl;
	}

	otec << "ZONE N=" << 8 << ", E=" << 2 << ", F=FEPOINT, ET=BRICK" << endl;
	otec << nod_temp[2*num_sec/3].x << "  " << nod_temp[2*num_sec/3].y << "  " << nod_temp[2*num_sec/3].z << endl;
	otec << nod_temp[2*num_sec+2*num_sec/3].x << "  " << nod_temp[2*num_sec+2*num_sec/3].y << "  " << nod_temp[2*num_sec+2*num_sec/3].z << endl;
	otec << nod_temp[2*num_sec/3].x - 1E-6 << "  " << nod_temp[2*num_sec/3].y - 1E-6 << "  " << nod_temp[2*num_sec/3].z << endl;
	otec << nod_temp[2*num_sec+2*num_sec/3].x - 1E-6 << "  " << nod_temp[2*num_sec+2*num_sec/3].y - 1E-6 << "  " << nod_temp[2*num_sec+2*num_sec/3].z << endl;
	otec << nod_temp[5*num_sec/3].x << "  " << nod_temp[5*num_sec/3].y << "  " << nod_temp[5*num_sec/3].z << endl;
	otec << nod_temp[2*num_sec+5*num_sec/3].x << "  " << nod_temp[2*num_sec+5*num_sec/3].y << "  " << nod_temp[2*num_sec+5*num_sec/3].z << endl;
	otec << nod_temp[5*num_sec/3].x  + 1E-6 << "  " << nod_temp[5*num_sec/3].y  + 1E-6 << "  " << nod_temp[5*num_sec/3].z << endl;
	otec << nod_temp[2*num_sec+5*num_sec/3].x  + 1E-6 << "  " << nod_temp[2*num_sec+5*num_sec/3].y  + 1E-6 << "  " << nod_temp[2*num_sec+5*num_sec/3].z << endl;
	otec << "1 2 4 3 1 2 4 3" << endl;
	otec << "5 6 8 7 5 6 8 7" << endl;
	otec << endl << endl;

	otec << "ZONE N=" << 8 << ", E=" << 2 << ", F=FEPOINT, ET=BRICK" << endl;
	otec << nod_temp[2*num_sec/3+1].x << "  " << nod_temp[2*num_sec/3+1].y << "  " << nod_temp[2*num_sec/3+1].z << endl;
	otec << nod_temp[2*num_sec+2*num_sec/3+1].x << "  " << nod_temp[2*num_sec+2*num_sec/3+1].y << "  " << nod_temp[2*num_sec+2*num_sec/3+1].z << endl;
	otec << nod_temp[2*num_sec/3+1].x - 1E-6 << "  " << nod_temp[2*num_sec/3+1].y - 1E-6 << "  " << nod_temp[2*num_sec/3+1].z << endl;
	otec << nod_temp[2*num_sec+2*num_sec/3+1].x - 1E-6 << "  " << nod_temp[2*num_sec+2*num_sec/3+1].y - 1E-6 << "  " << nod_temp[2*num_sec+2*num_sec/3+1].z << endl;
	otec << nod_temp[5*num_sec/3+1].x << "  " << nod_temp[5*num_sec/3+1].y << "  " << nod_temp[5*num_sec/3+1].z << endl;
	otec << nod_temp[2*num_sec+5*num_sec/3+1].x << "  " << nod_temp[2*num_sec+5*num_sec/3+1].y << "  " << nod_temp[2*num_sec+5*num_sec/3+1].z << endl;
	otec << nod_temp[5*num_sec/3+1].x  + 1E-6 << "  " << nod_temp[5*num_sec/3+1].y  + 1E-6 << "  " << nod_temp[5*num_sec/3+1].z << endl;
	otec << nod_temp[2*num_sec+5*num_sec/3+1].x  + 1E-6 << "  " << nod_temp[2*num_sec+5*num_sec/3+1].y  + 1E-6 << "  " << nod_temp[2*num_sec+5*num_sec/3+1].z << endl;
	otec << "1 2 4 3 1 2 4 3" << endl;
	otec << "5 6 8 7 5 6 8 7" << endl;
	otec << endl << endl;

	return 1;
}
//---------------------------------------------------------------------------
//The top and bottom of CNT networks (by threads in Tecplot)
int Tecplot_Export::Export_top_bottom_threads(const string &output_file_name, const struct cylinder &cyl, const vector<vector<Point_3D> > &cnts_toppois, const vector<vector<Point_3D> > &cnts_botpois)const
{
	ofstream otec(output_file_name.c_str());
	otec << "TITLE = CNT_Wires" << endl;
	otec << "VARIABLES = X, Y, Z" << endl;

	//---------------------------------------------------------------------------
	//Export a 3D cylinder
//	if(Export_cylinder(otec, cyl)==0) return 0;
	
	//---------------------------------------------------------------------------
	//Export 3D nanotube threads
	if(Export_nano_threads(otec, cnts_toppois)==0) return 0;
	if(Export_nano_threads(otec, cnts_botpois)==0) return 0;

	//---------------------------------------------------------------------------
	otec.close();

	return 1;
}
//---------------------------------------------------------------------------
//The top and bottom of CNT networks (by threads in Tecplot)
int Tecplot_Export::Export_top_bottom_threads(const string &output_file_name, const struct cylinder &cyl, const vector<vector<Point_3D> > &cnts_points)const
{
	ofstream otec(output_file_name.c_str());
	otec << "TITLE = CNT_Wires" << endl;
	otec << "VARIABLES = X, Y, Z" << endl;

	//---------------------------------------------------------------------------
	//Export a 3D cylinder
//	if(Export_cylinder(otec, cyl)==0) return 0;
	
	//---------------------------------------------------------------------------
	//Export 3D nanotube threads
	if(Export_nano_threads(otec, cnts_points)==0) return 0;

	//---------------------------------------------------------------------------
	otec.close();

	return 1;
}
//---------------------------------------------------------------------------
//The top and bottom of CNT networks (by tetrahedron meshes in Tecplot) //Attention: little parts of nanotube volumes out of the cylinder
int Tecplot_Export::Export_top_bottom_cnt_meshes(const string &output_file_name, const struct cylinder &cyl, const vector<vector<Point_3D> > &cnts_toppois, const vector<vector<Point_3D> > &cnts_botpois, const vector<double> &cnts_radius)const
{
	//输出纳米管线网格单Zone in Tecplot
	ofstream otec_tbm(output_file_name.c_str());
	//---------------------------------------------------------------------------
	//生成用于表示纳米管的节点及四面体网格
	vector<vector<Node> > cnts_nodes;
	vector<vector<Element> > cnts_eles;

	//Define a class of GenNetwork for calling for the function below
	GenNetwork *Gentemp = new GenNetwork;
	if(Gentemp->Generate_cnts_nodes_elements(cnts_nodes, cnts_eles, cnts_toppois, cnts_radius)==0) return 0;
	if(Export_cnts_meshes_singlezone(otec_tbm, cyl, cnts_nodes, cnts_eles)==0) return 0;

	//---------------------------------------------------------------------------
	cnts_nodes.clear();
	cnts_eles.clear();

	//Define a class of GenNetwork for calling for the function below
	if(Gentemp->Generate_cnts_nodes_elements(cnts_nodes, cnts_eles, cnts_botpois, cnts_radius)==0) return 0;
	if(Export_cnts_meshes_singlezone(otec_tbm, cyl, cnts_nodes, cnts_eles)==0) return 0;

	//---------------------------------------------------------------------------
	otec_tbm.close();

	delete Gentemp;

	return 1;
}
//---------------------------------------------------------------------------
//The top and bottom of CNT networks (by tetrahedron meshes in Tecplot) //Attention: little parts of nanotube volumes out of the cylinder
int Tecplot_Export::Export_top_bottom_cnt_meshes(const string &output_file_name, const struct cylinder &cyl, const vector<vector<Point_3D> > &cnts_points, const vector<double> &cnts_radius)const
{
	//输出纳米管线网格单Zone in Tecplot
	ofstream otec_tbm(output_file_name.c_str());
	//---------------------------------------------------------------------------
	//生成用于表示纳米管的节点及四面体网格
	vector<vector<Node> > cnts_nodes;
	vector<vector<Element> > cnts_eles;

	//Define a class of GenNetwork for calling for the function below
	GenNetwork *Gentemp = new GenNetwork;
	if(Gentemp->Generate_cnts_nodes_elements(cnts_nodes, cnts_eles, cnts_points, cnts_radius)==0) return 0;
	delete Gentemp;

	if(Export_cnts_meshes_singlezone(otec_tbm, cyl, cnts_nodes, cnts_eles)==0) return 0;
	//---------------------------------------------------------------------------
	otec_tbm.close();

	return 1;
}
//===========================================================================
