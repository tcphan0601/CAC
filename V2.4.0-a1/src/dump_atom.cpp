#include <string.h>
#include "dump_atom.h"
#include "domain.h"
#include "atom.h"
#include "element.h"
#include "modify.h"
#include "update.h"
#include "group.h"
#include "memory.h"
#include "error.h"
#include "comm.h"
#include "compute.h"
#include "universe.h"

using namespace CAC_NS;

#define ONELINE 256
#define DELTA 1048576

/*--------------------------------------------------------------------------------*/

DumpAtom::DumpAtom(CAC *cac, int narg, char **arg) : Dump(cac,narg,arg)
{
	if (narg < 5) error->all(FLERR,"Illegal dump atom command");
	nevery = universe->inumeric(FLERR,arg[3]);

	// will remove optional args after added dump custom style

	int iarg = 5;
	stress_flag = 0;
	force_flag = 0;
	displace_flag = 0;
	while (iarg < narg) {
		if (strcmp(arg[iarg],"stress") == 0) stress_flag = 1;
		else if (strcmp(arg[iarg],"force") == 0) force_flag = 1;
		else if (strcmp(arg[iarg],"displace") == 0) displace_flag = 1;
		else error->all(FLERR,"Illegal dump atom command");
		iarg++;
	}

	scale_flag = 0;
	image_flag = 0;
	buffer_allow = 1;
	buffer_flag = 1;

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
}

/*--------------------------------------------------------------------------------*/

DumpAtom::~DumpAtom()
{
	// check ncompute in case all fixes have already been deleted

	if (modify->ncompute && stress_flag) modify->delete_compute(id_stress);
	if (modify->ncompute && displace_flag) modify->delete_compute(id_displace);
}
/*--------------------------------------------------------------------------------*/

void DumpAtom::init_style()
{
	atom_size_one = 5 + 6*(stress_flag + force_flag) + 4*displace_flag;
	if (image_flag) atom_size_one += 3;
	elem_size_one = atom_size_one*element->max_nintpl;
	node_size_one = 0;

	int n;

	delete [] format;
	char *str;
	str = (char *) TAGINT_FORMAT " %d %g %g %g";
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

	// setup boundary string

	domain->boundary_string(boundstr);

	// setup column string

	delete [] columns;
	columns = new char[500];
	columns = strcpy(columns,"id type x y z");
	if (stress_flag)
		strcat(columns," sxx syy szz sxy sxz syz");
	if (force_flag)
		strcat(columns," vx vy vz fx fy fz");
	if (displace_flag)
		strcat(columns," dx dy dz d");

	header_choice = &DumpAtom::header_item;

	pack_atom_choice = &DumpAtom::pack_atom_noscale_noimage;
	pack_elem_choice = &DumpAtom::pack_elem_noscale_noimage;

	convert_atom_choice = &DumpAtom::convert_atom_noimage;
	convert_elem_choice = &DumpAtom::convert_elem_noimage;

	if (buffer_flag == 1) {
		write_atom_choice = &DumpAtom::write_string;
		write_elem_choice = &DumpAtom::write_string;
	}
	else {
		write_atom_choice = &DumpAtom::write_atom_lines_noimage;
		write_elem_choice = &DumpAtom::write_elem_lines_noimage;
	}

	// set stress compute pointer

	if (stress_flag) {
		int icompute = modify->find_compute(id_stress);
		if (icompute < 0) error->all(FLERR,"stress/atom ID for dump atom does not exist");
		stress = modify->compute[icompute];
	} else stress = NULL;

	if (displace_flag) {
		int icompute = modify->find_compute(id_displace);
		if (icompute < 0) error->all(FLERR,"displace/atom ID for dump tecplot does not exist");
		displace = modify->compute[icompute];
	} else displace = NULL;
}

/*------------------------------------------------------------------------------------------*/

void DumpAtom::header_item(bigint ndump)
{
	fprintf(fp,"ITEM: TIMESTEP\n");

	fprintf(fp,BIGINT_FORMAT "\n",update->ntimestep);

	fprintf(fp,"ITEM: NUMBER OF ATOMS\n");
	fprintf(fp,BIGINT_FORMAT "\n",ndump);
	fprintf(fp,"ITEM: BOX BOUNDS %s\n",boundstr);

	fprintf(fp,"%-1.16e %-1.16e\n",boxxlo,boxxhi);
	fprintf(fp,"%-1.16e %-1.16e\n",boxylo,boxyhi);
	fprintf(fp,"%-1.16e %-1.16e\n",boxzlo,boxzhi);
	fprintf(fp,"ITEM: ATOMS %s\n",columns);

}

/*-------------------------------------------------------------------------------------------*/

void DumpAtom::write_header(bigint ndump)
{
	if (multiproc) (this->*header_choice)(ndump);
	else if (filewriter) (this->*header_choice)(ndump);
}

/*-------------------------------------------------------------------------------------------*/

void DumpAtom::write_elem_header(bigint ndump)
{
}
/*-------------------------------------------------------------------------------------------*/

void DumpAtom::write_atom_header(bigint ndump)
{
}

/*-------------------------------------------------------------------------------------------*/

void DumpAtom::write_node_header(bigint ndump)
{
}

/*-----------------------------------------------------------------------------------------*/

void DumpAtom::pack_elem(tagint *ids)
{
	(this->*pack_elem_choice)(ids);
}

/*-----------------------------------------------------------------------------------------*/

void DumpAtom::pack_node(tagint *ids)
{
}

/*-----------------------------------------------------------------------------------------*/

void DumpAtom::pack_atom(tagint *ids)
{
	(this->*pack_atom_choice)(ids);
}

/*------------------------------------------------------------------------------------------*/

void DumpAtom::pack_elem_noscale_noimage(tagint *ids)
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

	tagint maxtag,maxtag_all,*tag;
	if (atom->natoms) {
		tag = atom->tag;
		maxtag = 0;
		for (int i = 0; i < atom->nlocal; i++)
			maxtag = MAX(maxtag,tag[i]);
		MPI_Allreduce(&maxtag,&maxtag_all,1,MPI_CAC_TAGINT,MPI_MAX,world);
	} else maxtag_all = 0;

	int nlocal = element->nlocal;
	double ***nodex = element->nodex;
	double ***nodev = element->nodev;
	double ***nodef = element->nodef;
	int *nintpl = element->nintpl;
	int *mask = element->mask;
	int *etype = element->etype;
	int *ctype = element->ctype;
	tag = element->tag;
	int *npe = element->npe;
	double ***shape_array = element->shape_array;
	double x,y,z,tmp;
	double vx,vy,vz,fx,fy,fz;
	int m = 0;
	int ietype;

	//bigint needtag,needtag_sum;
	//needtag = 0;
	//for (int i = 0; i < nlocal; i++)
	//  if (mask[i] & groupbit ) needtag += nintpl[etype[i]];

	//MPI_Scan(&needtag,&needtag_sum,1,MPI_CAC_BIGINT,MPI_SUM,world);

	// itag = 1st new tag that my interpolated atoms should use

	//bigint itag = maxtag_all + needtag_sum - needtag + 1;
	tagint itag = maxtag_all + 1;
	int maxintpl = element->maxintpl;
	for (int i = 0; i < nlocal; i++) {
		ietype = etype[i];
		if (mask[i] & groupbit)
			for (int iintpl = 0; iintpl < nintpl[ietype]; iintpl++) {
				ebuf[m++] = itag + (tag[i]-1)*maxintpl + iintpl;
				ebuf[m++] = ctype[i];
				x = y = z = 0.0;
				for (int inode = 0; inode < npe[ietype]; inode++) {
					x += shape_array[ietype][iintpl][inode]*nodex[i][inode][0];
					y += shape_array[ietype][iintpl][inode]*nodex[i][inode][1];
					z += shape_array[ietype][iintpl][inode]*nodex[i][inode][2];
				}
				ebuf[m++] = x;
				ebuf[m++] = y;
				ebuf[m++] = z;

				if (stress_flag) {
					for (int j = 0; j < 6; j++) {
						tmp = 0.0;
						for (int inode = 0; inode < npe[ietype]; inode++)
							tmp += shape_array[ietype][iintpl][inode]*nodes[i][inode][j];
						ebuf[m++] = tmp;
					}
				}
				if (force_flag) {
					vx = vy = vz = fx = fy = fz = 0;
					for (int inode = 0; inode < npe[ietype]; inode++) {
						vx += shape_array[ietype][iintpl][inode]*nodev[i][inode][0];
						vy += shape_array[ietype][iintpl][inode]*nodev[i][inode][1];
						vz += shape_array[ietype][iintpl][inode]*nodev[i][inode][2];
						fx += shape_array[ietype][iintpl][inode]*nodef[i][inode][0];
						fy += shape_array[ietype][iintpl][inode]*nodef[i][inode][1];
						fz += shape_array[ietype][iintpl][inode]*nodef[i][inode][2];
					}
					ebuf[m++] = vx;
					ebuf[m++] = vy;
					ebuf[m++] = vz;
					ebuf[m++] = fx;
					ebuf[m++] = fy;
					ebuf[m++] = fz;
				}
				if (displace_flag) {
					for (int j = 0; j < 4; j++) {
						tmp = 0.0;
						for (int inode = 0; inode < npe[ietype]; inode++)
							tmp += shape_array[ietype][iintpl][inode]*noded[i][inode][j];
						ebuf[m++] = tmp;
					}
				}
			}
	}
}

/*------------------------------------------------------------------------------------------*/

void DumpAtom::pack_atom_noscale_noimage(tagint *ids)
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

	int *type = atom->type;
	tagint *tag = atom->tag;
	int nlocal = atom->nlocal;

	m = n = 0;
	for (int i = 0; i < nlocal; i++) {
		if (mask[i] & groupbit) {
			abuf[m++] = tag[i];
			abuf[m++] = type[i];
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

/*--------------------------------------------------------------------------------------------*/
int DumpAtom::convert_node_string(int n,int size_one,double *mybuf)
{
	return 0;
}

/*--------------------------------------------------------------------------------------------*/
int DumpAtom::convert_atom_string(int n,int size_one,double *mybuf)
{
	return (this->*convert_atom_choice)(n,size_one,mybuf);
}

/*--------------------------------------------------------------------------------------------*/

int DumpAtom::convert_elem_string(int n,int size_one,double *mybuf)
{
	return (this->*convert_elem_choice)(n,size_one,mybuf);
}


/*--------------------------------------------------------------------------------------------*/

int DumpAtom::convert_atom_noimage(int n,int size_one,double *mybuf)
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
		                  static_cast<tagint> (mybuf[m]),
		                  static_cast<int> (mybuf[m+1]),
		                  mybuf[m+2],mybuf[m+3],mybuf[m+4]);
		m += 5;
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

/*--------------------------------------------------------------------------------------------*/

int DumpAtom::convert_elem_noimage(int n, int size_one, double *mybuf)
{
	int offset = 0;
	int m = 0;
	int *etype = element->etype;
	for (int i = 0; i < n; i++)
		for (int j = 0; j < element->nintpl[etype[i]]; j++) {
			if (offset + ONELINE > maxsebuf) {
				if ((bigint) maxsebuf + DELTA > MAXSMALLINT) return -1;
				maxsebuf += DELTA;
				memory->grow(sebuf,maxsebuf,"dump:sebuf");
			}
			offset += sprintf(&sebuf[offset],format,
			                  static_cast<tagint> (mybuf[m]),
			                  static_cast<int> (mybuf[m+1]),
			                  mybuf[m+2],mybuf[m+3],mybuf[m+4]);
			m += 5;
			if (stress_flag) {
				offset += sprintf(&sebuf[offset],format_stress,
				                  mybuf[m],mybuf[m+1],mybuf[m+2],
				                  mybuf[m+3],mybuf[m+4],mybuf[m+5]);
				m += 6;
			}
			if (force_flag) {
				offset += sprintf(&sebuf[offset],format_force,
				                  mybuf[m],mybuf[m+1],mybuf[m+2],
				                  mybuf[m+3],mybuf[m+4],mybuf[m+5]);
				m += 6;
			}
			if (displace_flag) {
				offset += sprintf(&sebuf[offset],format_displace,
				                  mybuf[m],mybuf[m+1],mybuf[m+2],mybuf[m+3]);
				m += 4;
			}
			offset += sprintf(&sebuf[offset],"\n");
		}
	return offset;
}


/*----------------------------------------------------------------------------------------------------*/

void DumpAtom::write_node_data(int n, double *mybuf)
{

}

/*----------------------------------------------------------------------------------------------------*/

void DumpAtom::write_atom_data(int n, double *mybuf)
{
	(this->*write_atom_choice)(n,mybuf);
}

/*----------------------------------------------------------------------------------------------------*/

void DumpAtom::write_elem_data(int n, double *mybuf)
{
	(this->*write_elem_choice)(n,mybuf);
}

/*---------------------------------------------------------------------------------------------------*/

void DumpAtom::write_string(int n, double *mybuf)
{
	int m;
	fwrite(mybuf, sizeof(char),n,fp);
}
/*-------------------------------------------------------------------------------------------------*/

void DumpAtom::write_atom_lines_noimage(int n, double *mybuf)
{
	int m = 0;
	for (int i = 0; i < n; i++) {
		fprintf(fp,format,
		        static_cast<tagint> (mybuf[m]),
		        static_cast<int> (mybuf[m+1]),
		        mybuf[m+2],mybuf[m+3],mybuf[m+4]);
		m += 5;
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

/*-------------------------------------------------------------------------------------------------*/

void DumpAtom::write_elem_lines_noimage(int n, double *mybuf)
{
	int m = 0;
	int *etype = element->etype;
	for (int i = 0; i < n; i++)
		for (int j = 0; j < element->nintpl[etype[i]]; j++) {
			fprintf(fp,format,
			        static_cast<tagint> (mybuf[m]),
			        static_cast<int> (mybuf[m+1]),
			        mybuf[m+2],mybuf[m+3],mybuf[m+4]);
			m += 5;
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
