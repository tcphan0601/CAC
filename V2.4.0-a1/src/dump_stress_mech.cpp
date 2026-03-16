#include <string.h>
#include "dump_stress_mech.h"
#include "domain.h"
#include "atom.h"
#include "element.h"
#include "update.h"
#include "group.h"
#include "memory.h"
#include "error.h"
#include "comm.h"
#include "modify.h"
#include "compute.h"

using namespace CAC_NS;

#define ONELINE 256
#define DELTA 1048576

enum {UNKNOWN,FCC,HCP,BCC,ICOS,OTHER};

/*---------------------------------------------
   temporary dump command to output stress from stress/mech command
   will remove when find a more general way to output
   -----------------------------------*/

DumpStressMech::DumpStressMech(CAC *cac, int narg, char **arg) : Dump(cac,narg,arg)
{

	if (narg != 6) error->all(FLERR,"Illegal dump stress/mech command");
	int compute_not_exist = 1;

	for (int i = 0; i < modify->ncompute; i++) {
		if (strcmp(modify->compute[i]->id,arg[5]) == 0) {
			compute = modify->compute[i];
			if (strcmp(compute->style,"stress/mech") == 0)
				compute_not_exist = 0;
			break;
		}
	}

	if (compute_not_exist) error->all(FLERR,"Invalid compute stress/mech id for dump stress/mech command");

	scale_flag = 0;
	image_flag = 0;
	buffer_allow = 1;
	buffer_flag = 1;
	format_default = NULL;
	format = NULL;
	clearstep = 1;
}

/*--------------------------------------------------------------------------------*/

void DumpStressMech::init_style()
{

	atom_size_one = 6;
	elem_size_one = 0;
	node_size_one = 0;

	int n;

	// setup line format

	delete [] format;
	char *str;
	str = (char *) "%g %g %g %g %g %g";
	n = strlen(str) + 2;
	format = new char[n];
	strcpy(format,str);
	strcat(format,"\n");

	// setup column string

	int d = (int) compute->compute_scalar();
	if (d == 0)
		columns = (char *) "\"x\", \"y\", \"z\", \"sxx\", \"sxy\", \"sxz\"";
	else if (d == 1)
		columns = (char *) "\"x\", \"y\", \"z\", \"syx\", \"syy\", \"syz\"";
	else
		columns = (char *) "\"x\", \"y\", \"z\", \"szx\", \"szy\", \"szz\"";

	header_choice = &DumpStressMech::header_item;

	pack_atom_choice = &DumpStressMech::pack_surfaces;

	convert_atom_choice = &DumpStressMech::convert_surfaces;

	if (buffer_flag == 1)
		write_atom_choice = &DumpStressMech::write_string;
	else
		write_atom_choice = &DumpStressMech::write_surfaces;


}
/*------------------------------------------------------------------------------------------*/

void DumpStressMech::header_item(bigint ndump)
{
	fprintf(fp,"title=\"TIMESTEP " BIGINT_FORMAT "\"\n",update->ntimestep);
	fprintf(fp,"variables=%s\n",columns);
}

/*-------------------------------------------------------------------------------------------*/

void DumpStressMech::write_header(bigint ndump)
{
	if (multiproc) (this->*header_choice)(ndump);
	else if (filewriter) (this->*header_choice)(ndump);
}

/*-------------------------------------------------------------------------------------------*/

void DumpStressMech::write_elem_header(bigint ndump)
{
}

/*-------------------------------------------------------------------------------------------*/

void DumpStressMech::write_atom_header(bigint ndump)
{
	fprintf(fp,"zone t = \"Surface Center\", f = point\n");
}

/*-------------------------------------------------------------------------------------------*/

void DumpStressMech::write_node_header(bigint ndump)
{
}

/*-----------------------------------------------------------------------------------------*/

void DumpStressMech::pack_elem(tagint *ids)
{
}

/*-----------------------------------------------------------------------------------------*/

void DumpStressMech::pack_node(tagint *ids)
{
}

/*-----------------------------------------------------------------------------------------*/

void DumpStressMech::pack_atom(tagint *ids)
{
	(this->*pack_atom_choice)(ids);
}

/*-------------------------------------
   pack surface coordinate and stress value
   only node 0 output
   -----------------------------------------------------*/


void DumpStressMech::pack_surfaces(tagint *ids)
{

	// invoke compute stress/mech

	if (compute->invoked_array != update->ntimestep)
		compute->compute_array();

	int f = 0;
	if (me) return;

	int nsurface = compute->size_array_rows;
	double **array = compute->array;

	int m = 0;

	for (int i = 0; i < nsurface; i++) {
		abuf[m++] = array[i][0];
		abuf[m++] = array[i][1];
		abuf[m++] = array[i][2];
		abuf[m++] = array[i][3];
		abuf[m++] = array[i][4];
		abuf[m++] = array[i][5];
	}

}

/*--------------------------------------------------------------------------------------------*/

int DumpStressMech::convert_node_string(int n,int size_one,double *mybuf)
{
	return 0;
}

/*--------------------------------------------------------------------------------------------*/
int DumpStressMech::convert_atom_string(int n,int size_one,double *mybuf)
{
	return (this->*convert_atom_choice)(n,size_one,mybuf);
}


/*--------------------------------------------------------------------------------------------*/

int DumpStressMech::convert_elem_string(int n,int size_one,double *mybuf)
{
	return 0;
}


/*--------------------------------------------------------------------------------------------*/

int DumpStressMech::convert_surfaces(int n,int size_one,double *mybuf)
{
	int offset = 0;
	int m = 0;
	for (int i = 0; i < n; i++) {
		if (offset + ONELINE > maxsabuf) {
			if ((bigint) maxsabuf + DELTA > MAXSMALLINT) return -1;
			maxsabuf += DELTA;
			memory->grow(sabuf, maxsabuf,"dump:sabuf");
		}
		offset += sprintf(&sabuf[offset],format,
		                  mybuf[m],mybuf[m+1],mybuf[m+2],mybuf[m+3],mybuf[m+4],mybuf[m+5]);
		m += 6;
	}
	return offset;
}

/*----------------------------------------------------------------------------------------------------*/

void DumpStressMech::write_node_data(int n, double *mybuf)
{

}

/*----------------------------------------------------------------------------------------------------*/

void DumpStressMech::write_atom_data(int n, double *mybuf)
{
	(this->*write_atom_choice)(n,mybuf);
}

/*----------------------------------------------------------------------------------------------------*/

void DumpStressMech::write_elem_data(int n, double *mybuf)
{
}

/*---------------------------------------------------------------------------------------------------*/

void DumpStressMech::write_string(int n, double *mybuf)
{
	int m;
	fwrite(mybuf, sizeof(char),n,fp);
}
/*-------------------------------------------------------------------------------------------------*/

void DumpStressMech::write_surfaces(int n, double *mybuf)
{
	int m = 0;
	for (int i = 0; i < n; i++) {
		fprintf(fp,format,mybuf[m+1],mybuf[m+2],mybuf[m+3],mybuf[m+4],mybuf[m+5]);
		m += 6;
	}
}

/*-----------------------------------------------
   Counting number of surface this proc contribute to dump
   --------------------------------------------------*/

int DumpStressMech::count_atoms()
{
	if (compute->invoked_array != update->ntimestep)
		compute->compute_array();

	if (me == 0)
		return compute->size_array_rows;
	else return 0;
}

/*-------------------------------------------------------------------------------------------------*/

int DumpStressMech::count_elements()
{
	return 0;
}


/*-------------------------------------------------------------------------------------------------*/

int DumpStressMech::count_nodes()
{
	return 0;
}

/*-------------------------------------------------------------------------------------------------*/

int DumpStressMech::count_intpl_atoms()
{
	return 0;
}
