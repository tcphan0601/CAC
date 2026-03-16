#include <string.h>
#include "dump_tecplot.h"
#include "compute.h"
#include "domain.h"
#include "atom.h"
#include "element.h"
#include "modify.h"
#include "update.h"
#include "group.h"
#include "memory.h"
#include "error.h"
#include "comm.h"
#include "universe.h"

using namespace CAC_NS;

#define ONELINE 256
#define DELTA 1048576

// same as in element.cpp

enum{QUADRILATERAL,TRIANGLE,HEXAHEDRON,TETRAHEDRON,OCTAHEDRON,PYRAMID,WEDGE};

/*--------------------------------------------------------------------------------*/

DumpTecplot::DumpTecplot(CAC *cac, int narg, char **arg) : Dump(cac,narg,arg)
{
<<<<<<< HEAD:V2.4.0-a1/src/dump_tecplot.cpp
  if (narg < 5) error->all(FLERR,"Illegal dump tecplot command");
  nevery = universe->inumeric(FLERR,arg[3]);

  if (domain->dimension == 3) 
    npe_connect = 8;
  else 
    npe_connect = 4;

  scale_flag = 0;
  image_flag = 0;
  buffer_allow = 1;
  buffer_flag = 1;
  clearstep = 1;


  int iarg = 5;
  stress_flag = 0;
  force_flag = 0;
  displace_flag = 0;
  average_flag = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"stress") == 0) stress_flag = 1;
    else if (strcmp(arg[iarg],"force") == 0) force_flag = 1;
    else if (strcmp(arg[iarg],"displace") == 0) displace_flag = 1;
    else if (strcmp(arg[iarg],"ave") == 0) {
      if (element->element_cluster_flag)
        average_flag = 1;
      else 
        error->warning(FLERR,"No elements with multiple atoms per unit cell, ave keyword is ignored");
    } else error->all(FLERR,"Illegal dump tecplot command");
    iarg++;
  }

  int n;

  if (stress_flag) {
    n = strlen(id) + strlen("_DUMP_STRESS") + 1;
    id_stress = new char[n];
    strcpy(id_stress,id);
    strcat(id_stress,"_DUMP_STRESS");
    char **newarg = new char*[4];
    newarg[0] = id_stress;
    newarg[1] = arg[1];
    newarg[2] = (char *) "stress/atom";
    newarg[3] = (char *) "NULL";
    modify->add_compute(4,newarg);
    delete [] newarg;
  }

  if (displace_flag) {
    n = strlen(id) + strlen("_DUMP_DISPLACE") + 1;
    id_displace = new char[n];
    strcpy(id_displace,id);
    strcat(id_displace,"_DUMP_DISPLACE");
    char **newarg = new char*[3];
    newarg[0] = id_displace;
    newarg[1] = arg[1];
    newarg[2] = (char *) "displace/atom";
    modify->add_compute(3,newarg);
    delete [] newarg;
  }

  if (average_flag) {
    if (stress_flag) comm_elem_forward += 6*element->max_npe;
    if (force_flag) comm_elem_forward += 6*element->max_npe;
    if (displace_flag) comm_elem_forward += 4*element->max_npe;
  }
=======
	if (narg < 5) error->all(FLERR,"Illegal dump cac command");
	nevery = universe->inumeric(FLERR,arg[3]);
	scale_flag = 0;
	image_flag = 0;
	buffer_allow = 1;
	buffer_flag = 1;
	clearstep = 1;


	int iarg = 5;
	stress_flag = 0;
	force_flag = 0;
	displace_flag = 0;
	average_flag = 0;
	while (iarg < narg) {
		if (strcmp(arg[iarg],"stress") == 0) stress_flag = 1;
		else if (strcmp(arg[iarg],"force") == 0) force_flag = 1;
		else if (strcmp(arg[iarg],"displace") == 0) displace_flag = 1;
		else if (strcmp(arg[iarg],"ave") == 0) {
			if (element->element_cluster_flag)
				average_flag = 1;
			else
				error->warning(FLERR,"No elements with multiple atoms per unit cell, ave keyword is ignored");
		} else error->all(FLERR,"Illegal dump tecplot command");
		iarg++;
	}

	int n;

	if (stress_flag) {
		n = strlen(id) + strlen("_DUMP_STRESS") + 1;
		id_stress = new char[n];
		strcpy(id_stress,id);
		strcat(id_stress,"_DUMP_STRESS");
		char **newarg = new char*[4];
		newarg[0] = id_stress;
		newarg[1] = arg[1];
		newarg[2] = (char *) "stress/atom";
		newarg[3] = (char *) "NULL";
		modify->add_compute(4,newarg);
		delete [] newarg;
	}

	if (displace_flag) {
		n = strlen(id) + strlen("_DUMP_DISPLACE") + 1;
		id_displace = new char[n];
		strcpy(id_displace,id);
		strcat(id_displace,"_DUMP_DISPLACE");
		char **newarg = new char*[3];
		newarg[0] = id_displace;
		newarg[1] = arg[1];
		newarg[2] = (char *) "displace/atom";
		modify->add_compute(3,newarg);
		delete [] newarg;
	}

	if (average_flag) {
		if (stress_flag) comm_elem_forward += 6*element->max_npe;
		if (force_flag) comm_elem_forward += 6*element->max_npe;
		if (displace_flag) comm_elem_forward += 4*element->max_npe;
	}
>>>>>>> 35213ed516bb7dded2274611730a1df3fc3c02a5:V2.4.0-a1/src/dump_tec_plot.cpp
}

/*--------------------------------------------------------------------------------*/

DumpTecplot::~DumpTecplot()
{
	if (modify->ncompute && stress_flag) modify->delete_compute(id_stress);
	if (modify->ncompute && displace_flag) modify->delete_compute(id_displace);

}

/*--------------------------------------------------------------------------------*/

<<<<<<< HEAD:V2.4.0-a1/src/dump_tecplot.cpp
void DumpTecplot::init_style() 
{

  atom_size_one = node_size_one = 
    3 + 6*(stress_flag + force_flag) + 4*displace_flag;
  elem_size_one = element->max_npe+1;

  int n;

  // format for writing atom and node section 

  delete [] format;
  char *str;
  str = (char *) "%g %g %g";
  n = strlen(str);
  format = new char[n];
  strcpy(format,str);

  if (stress_flag) {
    delete [] format_stress;
    str = (char *) " %g %g %g %g %g %g";
    n = strlen(str);
    format_stress = new char[n];
    strcpy(format_stress,str);
  }

  if (force_flag) {
    delete [] format_force;
    str = (char *) " %g %g %g %g %g %g";
    n = strlen(str);
    format_force = new char[n];
    strcpy(format_force,str);
  }

  if (displace_flag) {
    delete [] format_displace;
    str = (char *) " %g %g %g %g";
    n = strlen(str);
    format_displace = new char[n];
    strcpy(format_displace,str);
  }

  // format for writing element section (node connectivity)

  delete [] format_elem;
  str = (char *) TAGINT_FORMAT " ";
  n = strlen(str) + 2;
  format_elem = new char[n];
  strcpy(format_elem, str);

  // setup column string 

  delete [] columns; 
  columns = new char[500];
  columns = strcpy(columns,"\"x\", \"y\", \"z\"");
  if (stress_flag) 
    strcat(columns,", \"sxx\", \"syy\", \"szz\", \"sxy\", \"sxz\", \"syz\"");
  if (force_flag) 
    strcat(columns,", \"vx\", \"vy\", \"vz\", \"fx\", \"fy\", \"fz\"");
  if (displace_flag) 
    strcat(columns,", \"dx\", \"dy\", \"dz\", \"d\"");

  header_choice = &DumpTecplot::header_item;
  atom_header_choice = &DumpTecplot::atom_header_item;
  elem_header_choice = &DumpTecplot::elem_header_item;
  node_header_choice = &DumpTecplot::node_header_item;

  pack_atom_choice = &DumpTecplot::pack_atom_noscale_noimage;
  pack_node_choice = &DumpTecplot::pack_node_noscale_noimage;
  pack_elem_choice = &DumpTecplot::pack_elem_noscale_noimage;

  convert_atom_choice = &DumpTecplot::convert_atom_noimage;
  convert_node_choice = &DumpTecplot::convert_node_noimage;
  convert_elem_choice = &DumpTecplot::convert_elem_noimage;

  if (buffer_flag == 1) {
    write_atom_choice = &DumpTecplot::write_string;
    write_node_choice = &DumpTecplot::write_string;
    write_elem_choice = &DumpTecplot::write_string;
  } else {
    write_atom_choice = &DumpTecplot::write_lines_noimage;
    write_elem_choice = &DumpTecplot::write_elem_lines_noimage;
    write_node_choice = &DumpTecplot::write_lines_noimage;
  }

  // set stress/displace compute pointer

  if (stress_flag) {
    int icompute = modify->find_compute(id_stress);
    if (icompute < 0) error->all(FLERR,"stress/atom ID for dump tecplot does not exist");
    stress = modify->compute[icompute];
  } else stress = NULL;

  if (displace_flag) {
    int icompute = modify->find_compute(id_displace);
    if (icompute < 0) error->all(FLERR,"displace/atom ID for dump tecplot does not exist");
    displace = modify->compute[icompute];
  } else displace = NULL;
=======
void DumpTecPlot::init_style()
{

	atom_size_one = node_size_one =
		3 + 6*(stress_flag + force_flag) + 4*displace_flag;
	elem_size_one = element->max_npe+1;

	int n;

	// format for writing atom and node section

	delete [] format;
	char *str;
	str = (char *) "%g %g %g";
	n = strlen(str)+2;
	format = new char[n];
	strcpy(format,str);

	if (stress_flag) {
		delete [] format_stress;
		str = (char *) " %g %g %g %g %g %g";
		n = strlen(str)+2;
		format_stress = new char[n];
		strcpy(format_stress,str);
	}

	if (force_flag) {
		delete [] format_force;
		str = (char *) " %g %g %g %g %g %g";
		n = strlen(str)+2;
		format_force = new char[n];
		strcpy(format_force,str);
	}

	if (displace_flag) {
		delete [] format_displace;
		str = (char *) " %g %g %g %g";
		n = strlen(str)+2;
		format_displace = new char[n];
		strcpy(format_displace,str);
	}

	// format for writing element section (node connectivity)

	delete [] format_elem;
	str = (char *) TAGINT_FORMAT " ";
	n = strlen(str) + 2;
	format_elem = new char[n];
	strcpy(format_elem, str);

	// setup column string

	delete [] columns;
	columns = new char[500];
	columns = strcpy(columns,"\"x\", \"y\", \"z\"");
	if (stress_flag)
		strcat(columns,", \"sxx\", \"syy\", \"szz\", \"sxy\", \"sxz\", \"syz\"");
	if (force_flag)
		strcat(columns,", \"vx\", \"vy\", \"vz\", \"fx\", \"fy\", \"fz\"");
	if (displace_flag)
		strcat(columns,", \"dx\", \"dy\", \"dz\", \"d\"");

	header_choice = &DumpTecPlot::header_item;
	atom_header_choice = &DumpTecPlot::atom_header_item;
	elem_header_choice = &DumpTecPlot::elem_header_item;
	node_header_choice = &DumpTecPlot::node_header_item;

	pack_atom_choice = &DumpTecPlot::pack_atom_noscale_noimage;
	pack_node_choice = &DumpTecPlot::pack_node_noscale_noimage;
	pack_elem_choice = &DumpTecPlot::pack_elem_noscale_noimage;

	convert_atom_choice = &DumpTecPlot::convert_atom_noimage;
	convert_node_choice = &DumpTecPlot::convert_node_noimage;
	convert_elem_choice = &DumpTecPlot::convert_elem_noimage;

	if (buffer_flag == 1) {
		write_atom_choice = &DumpTecPlot::write_string;
		write_node_choice = &DumpTecPlot::write_string;
		write_elem_choice = &DumpTecPlot::write_string;
	} else {
		write_atom_choice = &DumpTecPlot::write_lines_noimage;
		write_elem_choice = &DumpTecPlot::write_elem_lines_noimage;
		write_node_choice = &DumpTecPlot::write_lines_noimage;
	}

	// set stress/displace compute pointer

	if (stress_flag) {
		int icompute = modify->find_compute(id_stress);
		if (icompute < 0) error->all(FLERR,"stress/atom ID for dump tecplot does not exist");
		stress = modify->compute[icompute];
	} else stress = NULL;

	if (displace_flag) {
		int icompute = modify->find_compute(id_displace);
		if (icompute < 0) error->all(FLERR,"displace/atom ID for dump tecplot does not exist");
		displace = modify->compute[icompute];
	} else displace = NULL;
>>>>>>> 35213ed516bb7dded2274611730a1df3fc3c02a5:V2.4.0-a1/src/dump_tec_plot.cpp
}

/*------------------------------------------------------------------------------------------*/

void DumpTecplot::header_item(bigint ndump)
{
	fprintf(fp, "title = \"TIMESTEP " BIGINT_FORMAT "\"\n",update->ntimestep);
	fprintf(fp, "variables = %s\n",columns);
}

/*------------------------------------------------------------------------------------------*/

void DumpTecplot::node_header_item(bigint ndump)
{
<<<<<<< HEAD:V2.4.0-a1/src/dump_tecplot.cpp
  fprintf(fp,"zone t = \"Coarse Element\",");
  fprintf(fp," n = " BIGINT_FORMAT " e = " BIGINT_FORMAT,ndump,nelem_total);
  if (domain->dimension == 3) 
    fprintf(fp," datapacking = point, zonetype = febrick\n");
  else 
    fprintf(fp," datapacking = point, zonetype = fequadrilateral\n");
=======
	fprintf(fp,"zone t = \"Coarse Element\",");
	fprintf(fp," n = " BIGINT_FORMAT " e = " BIGINT_FORMAT,ndump,nelem_total);
	fprintf(fp," datapacking = point, zonetype = febrick\n");
>>>>>>> 35213ed516bb7dded2274611730a1df3fc3c02a5:V2.4.0-a1/src/dump_tec_plot.cpp
}

/*------------------------------------------------------------------------------------------*/

void DumpTecplot::atom_header_item(bigint ndump)
{
<<<<<<< HEAD:V2.4.0-a1/src/dump_tecplot.cpp
  fprintf(fp,"zone t = \"Discrete Atoms\", datapacking = point\n");
=======
	fprintf(fp,"zone t = \"Discrete Atoms\", f = point\n");
>>>>>>> 35213ed516bb7dded2274611730a1df3fc3c02a5:V2.4.0-a1/src/dump_tec_plot.cpp
}

/*------------------------------------------------------------------------------------------*/

void DumpTecplot::elem_header_item(bigint ndump)
{
	// no header for element section (node connectivity)
}

/*-------------------------------------------------------------------------------------------*/

void DumpTecplot::write_header(bigint ndump)
{
	if (multiproc) (this->*header_choice)(ndump);
	else if (filewriter) (this->*header_choice)(ndump);
}

/*-------------------------------------------------------------------------------------------*/

void DumpTecplot::write_elem_header(bigint ndump)
{
	if (multiproc) (this->*elem_header_choice)(ndump);
	else if (filewriter) (this->*elem_header_choice)(ndump);
}
/*-------------------------------------------------------------------------------------------*/

void DumpTecplot::write_atom_header(bigint ndump)
{
	if (multiproc) (this->*atom_header_choice)(ndump);
	else if (filewriter) (this->*atom_header_choice)(ndump);
}

/*-------------------------------------------------------------------------------------------*/

void DumpTecplot::write_node_header(bigint ndump)
{
	if (multiproc) (this->*node_header_choice)(ndump);
	else if (filewriter) (this->*node_header_choice)(ndump);
}

/*-----------------------------------------------------------------------------------------*/

void DumpTecplot::pack_elem(tagint *ids)
{
	(this->*pack_elem_choice)(ids);
}

/*-----------------------------------------------------------------------------------------*/

void DumpTecplot::pack_node(tagint *ids)
{
	(this->*pack_node_choice)(ids);
}

/*-----------------------------------------------------------------------------------------*/

void DumpTecplot::pack_atom(tagint *ids)
{
	(this->*pack_atom_choice)(ids);
}

/*------------------------------------------------------------------------------------------*/

void DumpTecplot::pack_elem_noscale_noimage(tagint *ids)
{
<<<<<<< HEAD:V2.4.0-a1/src/dump_tecplot.cpp
  int *mask = element->mask;
  int nlocal = element->nlocal;
  tagint *tag = element->tag;
  tagint **nodetag = element->nodetag;

  int *npe = element->npe;
  int *etype = element->etype;
  int **element_clusters = element->element_clusters;
  int *element_shape_ids = element->element_shape_ids;
  int m = 0;
  int inpe,itype;

  element->reset_node_tags(groupbit,average_flag);

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      if (average_flag && tag[i] != element_clusters[i][0]) 
        continue;

      itype = etype[i];
      if (element_shape_ids[itype] == QUADRILATERAL) {
        for (int j = 0; j < 4; j++) 
          ebuf[m++] = nodetag[i][j];

      } else if (element_shape_ids[itype] == TRIANGLE) {

        // node connectivity scheme: 1 2 3 3
        
        for (int j = 0; j < 3; j++) 
          ebuf[m++] = nodetag[i][j];
        ebuf[m++] = nodetag[i][2];

      } else if (element_shape_ids[itype] == HEXAHEDRON) {

        for (int j = 0; j < 8; j++) 
          ebuf[m++] = nodetag[i][j];

      } else if (element_shape_ids[itype] == PYRAMID) {

        // node connectivity scheme: 1 2 3 4 5 5 5 5

        for (int j = 0; j < 4; j++) 
          ebuf[m++] = nodetag[i][j];
        for (int j = 4; j < 8; j++)
          ebuf[m++] = nodetag[i][4];

      } else if (element_shape_ids[itype] == TETRAHEDRON) {

        // node connectivity scheme: 1 2 3 3 4 4 4 4

        for (int j = 0; j < 3; j++) 
          ebuf[m++] = nodetag[i][j];
        ebuf[m++] = nodetag[i][2];
        for (int j = 4; j < 8; j++) 
          ebuf[m++] = nodetag[i][3];
      }
    }

} 
=======
	int *mask = element->mask;
	int nlocal = element->nlocal;
	tagint *tag = element->tag;
	tagint **nodetag = element->nodetag;

	int *npe = element->npe;
	int *etype = element->etype;
	int **element_clusters = element->element_clusters;
	int m = 0;
	int inpe;

	// determine maximum number of node tags for this proc

	tagint maxtag = 0;
	if (average_flag) {
		for (int i = 0; i < nlocal; i++)
			if (mask[i] & groupbit &&
			    tag[i] == element_clusters[i][0])
				maxtag += npe[etype[i]];
	} else {
		for (int i = 0; i < nlocal; i++)
			if (mask[i] & groupbit) maxtag += npe[etype[i]];
	}

	tagint maxtag_sum;
	MPI_Scan(&maxtag,&maxtag_sum,1,MPI_CAC_TAGINT,MPI_SUM,world);

	// itag = 1st tag that my nodes in group should use

	tagint itag = maxtag_sum - maxtag + 1;

	if (average_flag) {
		for (int i = 0; i < nlocal; i++)
			if (mask[i] & groupbit &&
			    tag[i] == element_clusters[i][0]) {
				inpe = npe[etype[i]];
				ebuf[m++] = inpe;
				for (int j = 0; j < inpe; j++)
					ebuf[m++] = itag++;
			}
	} else {
		for (int i = 0; i < nlocal; i++)
			if (mask[i] & groupbit) {
				inpe = npe[etype[i]];
				ebuf[m++] = inpe;
				for (int j = 0; j < inpe; j++)
					ebuf[m++] = itag++;
			}
	}
}
>>>>>>> 35213ed516bb7dded2274611730a1df3fc3c02a5:V2.4.0-a1/src/dump_tec_plot.cpp

/*------------------------------------------------------------------------------------------*/

void DumpTecplot::pack_atom_noscale_noimage(tagint *ids)
{
	double **s,**d;
	if (stress_flag) {
		if (stress->invoked_peratom != update->ntimestep) {
			stress->compute_peratom();
			stress->addstep(update->ntimestep+nevery);
		}
		s = stress->array_atom;
	}
	if (displace_flag) {
		if (displace->invoked_peratom != update->ntimestep) {
			displace->compute_peratom();
			displace->addstep(update->ntimestep+nevery);
		}
		d = displace->array_atom;
	}
	int m,n;
	int *mask = atom->mask;
	double **x = atom->x;
	double **v = atom->v;
	double **f = atom->f;

	int nlocal = atom->nlocal;

	m = n = 0;
	for (int i = 0; i < nlocal; i++) {
		if (mask[i] & groupbit) {
			abuf[m++] = x[i][0];
			abuf[m++] = x[i][1];
			abuf[m++] = x[i][2];
			if (stress_flag) {
				abuf[m++] = s[i][0];
				abuf[m++] = s[i][1];
				abuf[m++] = s[i][2];
				abuf[m++] = s[i][3];
				abuf[m++] = s[i][4];
				abuf[m++] = s[i][5];
			}
			if (force_flag) {
				abuf[m++] = v[i][0];
				abuf[m++] = v[i][1];
				abuf[m++] = v[i][2];
				abuf[m++] = f[i][0];
				abuf[m++] = f[i][1];
				abuf[m++] = f[i][2];
			}
			if (displace_flag) {
				abuf[m++] = d[i][0];
				abuf[m++] = d[i][1];
				abuf[m++] = d[i][2];
				abuf[m++] = d[i][3];
			}
			if (ids) ids[n++] = i;
		}
	}
}


/*------------------------------------------------------------------------------------------*/

void DumpTecplot::pack_node_noscale_noimage(tagint *ids)
{
	double ***nodes,***noded;
	if (stress_flag) {
		if (stress->invoked_peratom != update->ntimestep) {
			stress->compute_peratom();
			stress->addstep(update->ntimestep+nevery);
		}
		nodes = stress->array_node;
	}
	if (displace_flag) {
		if (displace->invoked_peratom != update->ntimestep) {
			displace->compute_peratom();
			displace->addstep(update->ntimestep+nevery);
		}
		noded = displace->array_node;
	}

	int m,n,ii,i,j,k;
	int *mask = element->mask;
	double ***nodex = element->nodex;
	double ***nodev = element->nodev;
	double ***nodef = element->nodef;
	int nlocal = element->nlocal;
	int *npe = element->npe;
	int *etype = element->etype;

	m = n = 0;

	if (average_flag) {
		int **element_clusters = element->element_clusters;
		tagint *tag = element->tag;
		int *apc = element->apc;

		// forward comm of ghost stress or force/velocity

		comm->forward_comm_dump(this);

		// build element map to find elements in clusters

		int mapflag = 0;
		if (element->map_style == 0) {
			mapflag = 1;
			element->map_init();
			element->map_set();
		}

		double ix,iy,iz,s0,s1,s2,s3,s4,s5;
		double vx,vy,vz,fx,fy,fz,d0,d1,d2,d3;
		for (i = 0; i < nlocal; i++) {
			if (tag[i] != element_clusters[i][0]) continue;
			int iapc = apc[etype[i]];
			int inpe = npe[etype[i]];
			for (j = 0; j < inpe; j++) {
				ix = nodex[i][j][0];
				iy = nodex[i][j][1];
				iz = nodex[i][j][2];
				if (stress_flag) {
					s0 = nodes[i][j][0];
					s1 = nodes[i][j][1];
					s2 = nodes[i][j][2];
					s3 = nodes[i][j][3];
					s4 = nodes[i][j][4];
					s5 = nodes[i][j][5];
				}
				if (force_flag) {
					vx = nodev[i][j][0];
					vy = nodev[i][j][1];
					vz = nodev[i][j][2];
					fx = nodef[i][j][0];
					fy = nodef[i][j][1];
					fz = nodef[i][j][2];
				}
				if (displace_flag) {
					d0 = noded[i][j][0];
					d1 = noded[i][j][1];
					d2 = noded[i][j][2];
					d3 = noded[i][j][3];
				}
				for (k = 1; k < iapc; k++) {
					ii = element->map(element_clusters[i][k]);
					ix += nodex[ii][j][0];
					iy += nodex[ii][j][1];
					iz += nodex[ii][j][2];
					if (stress_flag) {
						s0 += nodes[ii][j][0];
						s1 += nodes[ii][j][1];
						s2 += nodes[ii][j][2];
						s3 += nodes[ii][j][3];
						s4 += nodes[ii][j][4];
						s5 += nodes[ii][j][5];
					}
					if (force_flag) {
						vx += nodev[ii][j][0];
						vy += nodev[ii][j][1];
						vz += nodev[ii][j][2];
						fx += nodef[ii][j][0];
						fy += nodef[ii][j][1];
						fz += nodef[ii][j][2];
					}
					if (displace_flag) {
						d0 += noded[ii][j][0];
						d1 += noded[ii][j][1];
						d2 += noded[ii][j][2];
						d3 += noded[ii][j][3];
					}
				}
				nbuf[m++] = ix/iapc;
				nbuf[m++] = iy/iapc;
				nbuf[m++] = iz/iapc;
				if (stress_flag) {
					nbuf[m++] = s0/iapc;
					nbuf[m++] = s1/iapc;
					nbuf[m++] = s2/iapc;
					nbuf[m++] = s3/iapc;
					nbuf[m++] = s4/iapc;
					nbuf[m++] = s5/iapc;
				}
				if (force_flag) {
					nbuf[m++] = vx/iapc;
					nbuf[m++] = vy/iapc;
					nbuf[m++] = vz/iapc;
					nbuf[m++] = fx/iapc;
					nbuf[m++] = fy/iapc;
					nbuf[m++] = fz/iapc;
				}
				if (displace_flag) {
					nbuf[m++] = d0/iapc;
					nbuf[m++] = d1/iapc;
					nbuf[m++] = d2/iapc;
					nbuf[m++] = d3/iapc;
				}
			}
		}

		// clear map if needed

		if (mapflag) {
			element->map_delete();
			element->map_style = 0;
		}

	} else
		for (i = 0; i < nlocal; i++)
			if (mask[i] & groupbit) {
				for (j = 0; j < npe[etype[i]]; j++) {
					nbuf[m++] = nodex[i][j][0];
					nbuf[m++] = nodex[i][j][1];
					nbuf[m++] = nodex[i][j][2];
					if (stress_flag) {
						nbuf[m++] = nodes[i][j][0];
						nbuf[m++] = nodes[i][j][1];
						nbuf[m++] = nodes[i][j][2];
						nbuf[m++] = nodes[i][j][3];
						nbuf[m++] = nodes[i][j][4];
						nbuf[m++] = nodes[i][j][5];
					}
					if (force_flag) {
						nbuf[m++] = nodev[i][j][0];
						nbuf[m++] = nodev[i][j][1];
						nbuf[m++] = nodev[i][j][2];
						nbuf[m++] = nodef[i][j][0];
						nbuf[m++] = nodef[i][j][1];
						nbuf[m++] = nodef[i][j][2];
					}
					if (displace_flag) {
						nbuf[m++] = noded[i][j][0];
						nbuf[m++] = noded[i][j][1];
						nbuf[m++] = noded[i][j][2];
						nbuf[m++] = noded[i][j][3];
					}
				}
				if (ids) ids[n++] = i;
			}
}

/*--------------------------------------------------------------------------------------------*/

int DumpTecplot::convert_atom_string(int n, int size_one, double *mybuf)
{
	return (this->*convert_atom_choice)(n,size_one,mybuf);
}

/*--------------------------------------------------------------------------------------------*/

int DumpTecplot::convert_node_string(int n, int size_one, double *mybuf)
{
	return (this->*convert_node_choice)(n,size_one,mybuf);
}

/*--------------------------------------------------------------------------------------------*/

int DumpTecplot::convert_elem_string(int n, int size_one, double *mybuf)
{
	return (this->*convert_elem_choice)(n,size_one,mybuf);
}

/*--------------------------------------------------------------------------------------------*/

int DumpTecplot::convert_node_noimage(int n, int size_one, double *mybuf)
{
	int offset = 0;
	int m = 0;
	for (int i = 0; i < n; i++) {
		if (offset + ONELINE > maxsnbuf) {
			if ((bigint) maxsnbuf + DELTA > MAXSMALLINT) return -1;
			maxsnbuf += DELTA;
			memory->grow(snbuf,maxsnbuf,"dump:snbuf");
		}
		offset += sprintf(&snbuf[offset],format,
		                  mybuf[m],mybuf[m+1],mybuf[m+2]);
		m += 3;
		if (stress_flag) {
			offset += sprintf(&snbuf[offset],format_stress,
			                  mybuf[m],mybuf[m+1],mybuf[m+2],
			                  mybuf[m+3],mybuf[m+4],mybuf[m+5]);
			m += 6;
		}
		if (force_flag) {
			offset += sprintf(&snbuf[offset],format_force,
			                  mybuf[m],mybuf[m+1],mybuf[m+2],
			                  mybuf[m+3],mybuf[m+4],mybuf[m+5]);
			m += 6;
		}
		if (displace_flag) {
			offset += sprintf(&snbuf[offset],format_displace,
			                  mybuf[m],mybuf[m+1],mybuf[m+2],mybuf[m+3]);
			m += 4;
		}
		offset += sprintf(&snbuf[offset],"\n");
	}
	return offset;
}
/*--------------------------------------------------------------------------------------------*/

int DumpTecplot::convert_atom_noimage(int n, int size_one, double *mybuf)
{
	int offset = 0;
	int m = 0;
	for (int i = 0; i < n; i++) {
		if (offset + ONELINE > maxsabuf) {
			if ((bigint) maxsabuf + DELTA > MAXSMALLINT) return -1;
			maxsabuf += DELTA;
			memory->grow(sabuf,maxsabuf,"dump:sabuf");
		}
		offset += sprintf(&sabuf[offset],format,
		                  mybuf[m],mybuf[m+1],mybuf[m+2]);
		m += 3;
		if (stress_flag) {
			offset += sprintf(&sabuf[offset],format_stress,
			                  mybuf[m],mybuf[m+1],mybuf[m+2],
			                  mybuf[m+3],mybuf[m+4],mybuf[m+5]);
			m += 6;
		}
		if (force_flag) {
			offset += sprintf(&sabuf[offset],format_force,
			                  mybuf[m],mybuf[m+1],mybuf[m+2],
			                  mybuf[m+3],mybuf[m+4],mybuf[m+5]);
			m += 6;
		}
		if (displace_flag) {
			offset += sprintf(&sabuf[offset],format_displace,
			                  mybuf[m],mybuf[m+1],mybuf[m+2],mybuf[m+3]);
			m += 4;
		}
		offset += sprintf(&sabuf[offset],"\n");
	}
	return offset;
}

/*----------------------------------------------------------------------------------------------------*/

int DumpTecplot::convert_elem_noimage(int n, int size_one, double *mybuf)
{
<<<<<<< HEAD:V2.4.0-a1/src/dump_tecplot.cpp
  int offset = 0;
  int m = 0;
  for (int i = 0; i < n; i++) {
    if (offset + ONELINE > maxsebuf) {
      if ((bigint) maxsebuf + DELTA > MAXSMALLINT) return -1;
      maxsebuf += DELTA;
      memory->grow(sebuf,maxsebuf,"dump:sebuf");
    }
    for (int j = 0; j < npe_connect; j++)
      offset += sprintf(&sebuf[offset],format_elem,
          static_cast<tagint>(mybuf[m++]));

    offset += sprintf(&sebuf[offset],"\n");
  }
  return offset;
}	
=======
	int offset = 0;
	int m = 0;
	int k;
	for (int i = 0; i < n; i++) {
		if (offset + ONELINE > maxsebuf) {
			if ((bigint) maxsebuf + DELTA > MAXSMALLINT) return -1;
			maxsebuf += DELTA;
			memory->grow(sebuf,maxsebuf,"dump:sebuf");
		}
		k = static_cast<int>(mybuf[m++]);
		for (int j = 0; j < k; j++)
			offset += sprintf(&sebuf[offset],format_elem,
			                  static_cast<tagint>(mybuf[m++]));

		offset += sprintf(&sebuf[offset],"\n");
	}
	return offset;
}
>>>>>>> 35213ed516bb7dded2274611730a1df3fc3c02a5:V2.4.0-a1/src/dump_tec_plot.cpp

/*----------------------------------------------------------------------------------------------------*/

void DumpTecplot::write_node_data(int n, double *mybuf)
{
	(this->*write_node_choice)(n,mybuf);
}

/*----------------------------------------------------------------------------------------------------*/

void DumpTecplot::write_atom_data(int n, double *mybuf)
{
	(this->*write_atom_choice)(n,mybuf);
}

/*----------------------------------------------------------------------------------------------------*/

void DumpTecplot::write_elem_data(int n, double *mybuf)
{
	(this->*write_elem_choice)(n,mybuf);
}

/*---------------------------------------------------------------------------------------------------*/

void DumpTecplot::write_string(int n, double *mybuf)
{
	int m;
	fwrite(mybuf, sizeof(char),n,fp);
}
/*-------------------------------------------------------------------------------------------------*/

void DumpTecplot::write_lines_noimage(int n, double *mybuf)
{
	int m = 0;
	for (int i = 0; i < n; i++) {
		fprintf(fp,format,mybuf[m],mybuf[m+1],mybuf[m+2]);
		m += 3;
		if (stress_flag) {
			fprintf(fp,format_stress,mybuf[m],mybuf[m+1],
			        mybuf[m+2],mybuf[m+3],mybuf[m+4],mybuf[m+5]);
			m += 6;
		}
		if (force_flag) {
			fprintf(fp,format_force,mybuf[m],mybuf[m+1],
			        mybuf[m+2],mybuf[m+3],mybuf[m+4],mybuf[m+5]);
			m += 6;
		}
		if (displace_flag) {
			fprintf(fp,format_displace,mybuf[m],mybuf[m+1],
			        mybuf[m+2],mybuf[m+3]);
			m += 4;
		}
		fprintf(fp,"\n");
	}
}

/*----------------------------------------------------------------------------------------------------*/

void DumpTecplot::write_elem_lines_noimage(int n, double *mybuf)
{
<<<<<<< HEAD:V2.4.0-a1/src/dump_tecplot.cpp
  int m = 0;
  int k;
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < npe_connect; j++)
      fprintf(fp,format_elem, static_cast<tagint>(mybuf[m++]));
    fprintf(fp,"\n");
  }
=======
	int m = 0;
	int k;
	for (int i = 0; i < n; i++) {
		k = static_cast<int> (mybuf[m++]);
		for (int j = 1; j <= k; j++)
			fprintf(fp,format_elem, static_cast<tagint>(mybuf[m++]));
		fprintf(fp,"\n");
	}
>>>>>>> 35213ed516bb7dded2274611730a1df3fc3c02a5:V2.4.0-a1/src/dump_tec_plot.cpp
}

/* ---------------------------------------------------------------------- */

int DumpTecplot::pack_elem_forward_comm(int n, int *list, double *buf, int pbc_flag, int *pbc)
{
	int ii,i,j,m;
	int npe = element->max_npe;
	m = 0;

	if (stress_flag) {
		double ***nodes = stress->array_node;
		for (ii = 0; ii < n; ii++) {
			i = list[ii];
			for (j = 0; j < npe; j++) {
				buf[m++] = nodes[i][j][0];
				buf[m++] = nodes[i][j][1];
				buf[m++] = nodes[i][j][2];
				buf[m++] = nodes[i][j][3];
				buf[m++] = nodes[i][j][4];
				buf[m++] = nodes[i][j][5];
			}
		}
	}
	if (force_flag) {
		double ***nodef = element->nodef;
		double ***nodev = element->nodev;
		for (ii = 0; ii < n; ii++) {
			i = list[ii];
			for (j = 0; j < npe; j++) {
				buf[m++] = nodev[i][j][0];
				buf[m++] = nodev[i][j][1];
				buf[m++] = nodev[i][j][2];
				buf[m++] = nodef[i][j][0];
				buf[m++] = nodef[i][j][1];
				buf[m++] = nodef[i][j][2];
			}
		}
	}
	if (displace_flag) {
		double ***noded = displace->array_node;
		for (ii = 0; ii < n; ii++) {
			i = list[ii];
			for (j = 0; j < npe; j++) {
				buf[m++] = noded[i][j][0];
				buf[m++] = noded[i][j][1];
				buf[m++] = noded[i][j][2];
				buf[m++] = noded[i][j][3];
			}
		}
	}

	return m;
}

/* ---------------------------------------------------------------------- */

void DumpTecplot::unpack_elem_forward_comm(int n, int first, double *buf)
{
	int i,j,m,last;
	last = first + n;
	int npe = element->max_npe;
	m = 0;

	if (stress_flag) {
		double ***nodes = stress->array_node;
		for (i = first; i < last; i++)
			for (j = 0; j < npe; j++) {
				nodes[i][j][0] = buf[m++];
				nodes[i][j][1] = buf[m++];
				nodes[i][j][2] = buf[m++];
				nodes[i][j][3] = buf[m++];
				nodes[i][j][4] = buf[m++];
				nodes[i][j][5] = buf[m++];
			}
	}
	if (force_flag) {
		double ***nodef = element->nodef;
		double ***nodev = element->nodev;
		for (i = first; i < last; i++)
			for (j = 0; j < npe; j++) {
				nodev[i][j][0] = buf[m++];
				nodev[i][j][1] = buf[m++];
				nodev[i][j][2] = buf[m++];
				nodef[i][j][0] = buf[m++];
				nodef[i][j][1] = buf[m++];
				nodef[i][j][2] = buf[m++];
			}
	}
	if (displace_flag) {
		double ***noded = displace->array_node;
		for (i = first; i < last; i++)
			for (j = 0; j < npe; j++) {
				noded[i][j][0] = buf[m++];
				noded[i][j][1] = buf[m++];
				noded[i][j][2] = buf[m++];
				noded[i][j][3] = buf[m++];
			}
	}
}

/*----------------------------------------------------------------------------------------------------*/

int DumpTecplot::count_elements()
{
	if (average_flag) return element->count_element_clusters();
	else if (igroup == 0 || igroup == 1) return element->nlocal;
	int m = 0;
	int *mask = element->mask;
	for (int i = 0; i < element->nlocal; i++)
		if (mask[i] & groupbit) m++;
	return m;
}

/*----------------------------------------------------------------------------------------------------*/

int DumpTecplot::count_nodes()
{
	int n = 0;
	int *etype = element->etype;
	int *npe = element->npe;
	for (int i = 0; i < element->nlocal; i++)
		n += npe[etype[i]];
	return n;
}
