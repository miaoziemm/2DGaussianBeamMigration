/****************************************************
MPIFrame: MPI Frame for Gaussian Beam Migration (include Segy File Read)
Author: Chang Zhimiao
Email: changzm20@mails.jlu.edu.cn
Date: 2023.06.15
*****************************************************/

#include <stdio.h>
#include "SegyRead.h"
#include "split.h"
#include "par.h"
#include "mpi.h"
#include "common.h"
#include "alloc.h"
#include "main.h"
#include "cbm.h"

/* MPI Handle */
#define MASTER 0
#define REQUEST 11
#define PARAM 22
#define COORDX 33
#define COORDY 44
#define COORDZ 55
#define DATAX 66
#define DATAY 77
#define DATAZ 88

void parBcastSet(int nx, int nz,
				 float dx, float dz, float fx, float fz,
				 float vpavg, float vpmin,
				 int nt, float dt, int nstotal, int ntrps_max,
				 parBcast *pbc);

void parBcastGet(int *nx, int *nz,
				 float *dx, float *dz,
				 float *fx, float *fz,
				 float *vpavg, float *vpmin,
				 int *nt, float *dt, int *nstotal, int *ntrps_max,
				 parBcast *pbc);

int master(int argc, char **argv);
int slave(int argc, char **argv);
int xargc;
char **xargv;

int main(int argc, char **argv)
{

	/*   MPI init   */
	int rank, np;
	double starttime, endtime;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &np);
	starttime = MPI_Wtime();
	/*****************************************/
	if (np < 2)
		err("error: np at least 2!\n");
	if (rank == 0)
	{
		master(argc, argv); /* master process */
	}
	else
	{
		slave(argc, argv); /* slave process */
	}
	/* MPI finished  */
	endtime = MPI_Wtime();
	if (rank == 0)
		fprintf(stderr, "program takes %f seconds\n", endtime - starttime);
	MPI_Finalize();
	/****************************************/

	return (0);
}

int master(int argc, char **argv)
{
	int i, j, k, l, m, n;
	int ii, jj, kk, ll, mm, nn;
	int segyCount;
	int segyFileNum;
	char **index;	// index file name--segy files
	int nx;			/* number of horizontal samples in input model 	  */
	int nz;			/* number of vertical samples in input model      */
	int nt;			/* number of time samples in input date		  */
	int ix, iy, iz; /* index in x y z				  */
	int is, ir;		/* index in source and receiver		 	  */

	float dx;	  /* model spacing in input models		  */
	float dz;	  /* model spacing in input models		  */
	float dt;	  /* time sampling interval in input data	          */
	float fx, fz; /* first sample point in models			  */
	float zs;	  /* z coordinate of source			  */
	float zr;	  /* z coordinate of receiver			  */
	double vpavg; /* average velocity in input vpfile		  */
	float vpmin;  /* minimum velocity at surface in input vpfile	  */

	char *vpfile = ""; /* P wave velocity filename            */
	char *indexFile = "";
	FILE *vpfp; /* P wave velocity file pointer        */
				// FILE *fpdat;		/* file pointer to seismic data */
	FILE *indexfp;

	long ntrtotal; /* total trace number in input seismic data */
	int nstotal;   /* total shot number in input seismic data */
	int ntrps_max; /* maximum trace number in single shot gathers */
	int *ntrps;	   /* array[nstotal] containing ntr of every shot */
	coord *scd;
	int sx, gx, sy, gy;	  /* coordinate value in trace headers */
	float scalel, scalco; /* scales of elevation and coordinate */
	int ntr;			  /* trace number in a shot gather */
	int irug;			  /* index of rugged topography */
	int **dataTemp;
	int **dataTempTo;

	char *elevfile = "";
	FILE *elevfp;

	int direction; /* direction of migration */

	int endian = 0; /* endian of input data */

	/*   MPI init   */
	int rank, np, rank_req;
	MPI_Status status;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &np);

	initargs(argc, argv);

	/* Model par */
	if (!getparint("nx", &nx))
		nx = 100;
	if (!getparint("nz", &nz))
		nz = 100;
	if (!getparfloat("dx", &dx))
		dx = 10.0;
	if (!getparfloat("dz", &dz))
		dz = 10.0;
	if (!getparfloat("fx", &fx))
		fx = 0.0;
	if (!getparfloat("fz", &fz))
		fz = 0.0;
	if (!getparstring("vpfile", &vpfile))
		vpfile = "vpfile";
	if (!getparint("irug", &irug))
		irug = 0;
	if (!getparfloat("zs", &zs))
		zs = 0.0;
	if (!getparfloat("zr", &zr))
		zr = 0.0;
	if (!getparint("SeisFileNum", &segyFileNum))
		segyFileNum = 1;
	if (!getparstring("index", &indexFile))
		indexFile = "Index";
	if (!getparstring("elev", &elevfile))
		elevfile = "elev";
	if (!getparint("direction", &direction))
		direction = 0;
	if (!getparint("endian", &endian))
		endian = 0;
	/* Data par */
	indexfp = fopen(indexFile, "r");
	index = IndexProcess(indexfp, segyFileNum);
	fclose(indexfp);

	/* Data par end */
	Model vm; /* velocity model */
	vm.nx = nx;
	vm.nz = nz;
	vm.dx = dx;
	vm.dz = dz;
	vm.fx = fx;
	vm.fz = fz;
	vm.vp = alloc2float(nz, nx);
	/* Model par end */

	/* Model par init */
	vpfp = fopen(vpfile, "rb");
	if (vpfp == NULL)
		err("error: can't open vpfile!\n");
	vpavg = 0.0;
	vpmin = 1000000.0;

	for (ix = 0; ix < nx; ix++)
	{
		for (iz = 0; iz < nz; iz++)
		{
			fread(&vm.vp[ix][iz], sizeof(float), 1, vpfp);
			vpavg = vpavg + vm.vp[ix][iz] / (nx * nz);
			if (vm.vp[ix][iz] < vpmin)
				vpmin = vm.vp[ix][iz];
		}
	}
	vm.vpavg = vpavg;
	vm.vpmin = vpmin;
	printf("rank::%d vpavg=%f vpmin=%f\n", rank, vm.vpavg, vm.vpmin);
	fclose(vpfp);
	/* Model init end */

	/* Model Bcast end */
	FILE *dataIn;
	dataIn = fopen(index[0], "rb");
	getShotInfo_su(dataIn, &nt, &dt, &ntrtotal, &nstotal, &ntrps_max, &scalel, &scalco, endian);
	fclose(dataIn);

	/* Model Bcast */
	parBcast *pbc;
	pbc = (parBcast *)malloc(sizeof(parBcast));
	parBcastSet(nx, nz, dx, dz, fx, fz, vm.vpavg, vm.vpmin, nt, dt, nstotal, ntrps_max, pbc);
	MPI_Bcast(pbc, sizeof(parBcast), MPI_BYTE, MASTER, MPI_COMM_WORLD);
	free(pbc);
	smooth2d(nz,nx,5,5,vm.vp);
	for (i = 0; i < nx; i++)
	{
		MPI_Bcast(&(vm.vp[i][0]), nz,
				  MPI_FLOAT, MASTER, MPI_COMM_WORLD);
	}

	parSendRecv psr;
	FILE *datain;
	segy *tr;

	/* Segy data loop */
	for (segyCount = 0; segyCount < segyFileNum; segyCount++)
	{
		if (segyCount != 0)
		{
			/* Model Bcast end */
			FILE *dataIn;
			dataIn = fopen(index[segyCount], "rb");
			getShotInfo_su(dataIn, &nt, &dt, &ntrtotal, &nstotal, &ntrps_max, &scalel, &scalco, endian);
			fclose(dataIn);

			/* Model Bcast */
			parBcast *pbc;
			pbc = (parBcast *)malloc(sizeof(parBcast));
			parBcastSet(nx, nz, dx, dz, fx, fz, vpavg, vpmin, nt, dt, nstotal, ntrps_max, pbc);
			MPI_Bcast(pbc, sizeof(parBcast), MPI_BYTE, MASTER, MPI_COMM_WORLD);
			free(pbc);
		}

		dataTemp = alloc2int(nt, ntrps_max);
		dataTempTo = alloc2int(nt, ntrps_max);
		tr = (segy *)malloc(sizeof(segy));

		// datain = fopen(index[segyCount], "rb");
		// if (datain == NULL)
		// 	err("error: can't open SeismicDataFile!\n");
		// getShotInfo(datain, &nt, &dt, &ntrtotal, &nstotal, &ntrps_max, &scalel, &scalco);
		// fclose(datain);

		ntrps = alloc1int(nstotal);
		scd = (coord *)malloc(sizeof(coord) * nstotal);
		datain = fopen(index[segyCount], "rb");

		if (fread(tr, 240, 1, datain) != 1)
			err("error: can't read trace header!\n");
		ir = 1;
		is = 0;
		sx = tr->sx;
		sy = tr->sy;
		gx = tr->gx;
		gy = tr->gy;

		/* apply coordinate scale */
		scd[0].x = (float)(tr->sx) * scalco;
		scd[0].y = (float)(tr->sy) * scalco;
		// printf("scalco=%f\n", scalco);
		if (irug)
			scd[0].z = (float)(-tr->selev + tr->sdepth) * scalel;
		else
			scd[0].z = zs;
		fseek(datain, nt * sizeof(float), SEEK_CUR);
		while (fread(tr, 240, 1, datain) == 1)
		{
			if (tr->sx != sx || tr->sy != sy)
			{
				ntrps[is] = ir;
				is++;
				sx = tr->sx;
				sy = tr->sy;
				scd[is].x = (float)(tr->sx) * scalco;
				scd[is].y = (float)(tr->sy) * scalco;
				if (irug)
					scd[is].z = (float)(-tr->selev + tr->sdepth) * scalel;
				else
					scd[is].z = zs;
				ir = 0;
			}
			ir++;

			fseek(datain, nt * sizeof(float), SEEK_CUR);
		}
		ntrps[is] = ir;
		fclose(datain);

		/* Shot par */
		Shot ss;
		ss.nt = nt;
		ss.dt = dt;
		ss.irug = irug;
		ss.zs = zs;
		ss.zr = zr; // for rugged
		ss.coordx = alloc1float(ntrps_max);
		ss.coordy = alloc1float(ntrps_max);
		if (irug)
			ss.coordz = alloc1float(ntrps_max);
		else
			ss.coordz = NULL;				  /* alloc in rugged */
		ss.data = alloc2float(nt, ntrps_max); /* alloc data space based on the maximum trace number */
		datain = fopen(index[segyCount], "rb");

		// for (is = 0; is < nstotal; is++)
		for (is = 0; is < nstotal; is++)
		{
			ntr = ntrps[is]; /* ntrps is the trace number for each shot */
			psr.is = is;
			psr.ntr = ntr;
			psr.xs = scd[is].x;
			psr.ys = scd[is].y;
			psr.zs = scd[is].z;

			/* Read in shot */
			for (ir = 0; ir < ntr; ir++)
			{
				// printf("scalco=%f\n", scalco);
				fread(tr, 240, 1, datain);

				fread(&ss.data[ir][0], sizeof(int), nt, datain);

				ss.coordx[ir] = (float)tr->gx * scalco;
				ss.coordy[ir] = (float)tr->gy * scalco;

				if (irug)
					ss.coordz[ir] = -(float)tr->gelev * scalel;
			}

			/* Shot Bcast */
			MPI_Recv(&rank_req, 1, MPI_INT, MPI_ANY_SOURCE,
					 REQUEST, MPI_COMM_WORLD, &status);
			/* send parameters */
			MPI_Send(&psr, sizeof(parSendRecv), MPI_BYTE, rank_req,
					 PARAM, MPI_COMM_WORLD);
			/*  recerve shot data */
			MPI_Send(&(ss.coordx[0]), ntr, MPI_FLOAT, rank_req,
					 COORDX, MPI_COMM_WORLD);
			MPI_Send(&(ss.coordy[0]), ntr, MPI_FLOAT, rank_req,
					 COORDY, MPI_COMM_WORLD);
			if (irug)
				MPI_Send(&(ss.coordz[0]), ntr, MPI_FLOAT, rank_req,
						 COORDZ, MPI_COMM_WORLD);
			MPI_Send(&(ss.data[0][0]), nt * ntr, MPI_FLOAT, rank_req,
					 DATAZ, MPI_COMM_WORLD);

		} // end shot loop
		free(tr);
		free(scd);
		free(ntrps);
		free2int(dataTemp);
		free2int(dataTempTo);
		if (segyCount != segyFileNum - 1)
		{
			for (ir = 1; ir < np; ir++)
			{
				MPI_Recv(&rank_req, 1, MPI_INT, MPI_ANY_SOURCE,
						 REQUEST, MPI_COMM_WORLD, &status);
				psr.is = -1;
				psr.ntr = -1;
				/* send parameters */
				MPI_Send(&psr, sizeof(parSendRecv), MPI_BYTE, rank_req,
						 PARAM, MPI_COMM_WORLD);
			}
		}
	}
	// save rank index
	FILE *rankIndex;
	rankIndex = fopen("rankIndex", "w");
	for (ir = 1; ir < np; ir++)
		fprintf(rankIndex, "./image/imgFile_%d\n", ir);
	fprintf(rankIndex, "END\n");
	fclose(rankIndex);

	printf("rank %d finished\n", rank);

	fclose(datain);
	for (ir = 1; ir < np; ir++)
	{
		/* receive request */
		MPI_Recv(&rank_req, 1, MPI_INT, MPI_ANY_SOURCE,
				 REQUEST, MPI_COMM_WORLD, &status);
		psr.is = -1;
		psr.ntr = -1;
		/* send parameters */
		MPI_Send(&psr, sizeof(parSendRecv), MPI_BYTE, rank_req,
				 PARAM, MPI_COMM_WORLD);
		fprintf(stderr, "rank %d finished\n", rank_req);
	}
	fprintf(stderr, "migration finished\n");

	// free1float(scd);
	// free1int(ntrps);
	// free2int(dataTemp);
	// free2int(dataTempTo);
	// free2(index);
	// free3float(vm.vp);
	// TODO:check free

	return 0;
} // end master

int slave(int argc, char **argv)
{
	int SegyCount = 0;
	int ii, i, j, k, it, ixx;
	int Count;
	int beamDirection = 0;
	int nx, nz, nt, ix, iy, iz, is;
	int nstotal;
	int ntrps_max;
	int ntr;
	float xs;
	float fx, fz;
	float dxline, dyline;
	int ixline, iyline;
	float fx_min, fx_max;	// the minimum and maximum coordinates of velocity field//
	float sx, sy, gx, gy;	// gx,gy are the minimum coordinates of one shot gater//
	float apx_min, apx_max; // the min and max coordinates of aperture of one shot//
	float aperx;			// aperture of crossline and inline//
	float bwh;
	int iaperx_min, iaperx_max; // min and max location of aperture in the velocity field//
	int naperx;
	int cx_min, cx_max;		// the minimum and maximum coordinates of receivers of one shot//
	float fcx_min, fcx_max; // minimum coordinates of receivers in the temporary velocity field //
	int nxline;				/* inline and crossline samples */
	int live, dead;
	float dx, dz, dt, fmin, fmax;
	float vpavg, vpmin;
	float amax, amin;
	float zr, zs;
	int nangle, iangle;
	float **tempv, **tempg;
	float **f;
	float *elevtemp;
	float *coord;

	int irug; /* index of rugged topography */
	int direction;

	char imgFileHandle[1024], angFileHandle[1024];
	FILE *imgFile;
	FILE *angfp;

	char *elevfile = "";
	FILE *elevfp;
	float *elev;

	float **g, **finalg;
	float **tang, **fang;

	/*   MPI init   */
	int rank, np;
	MPI_Status status;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &np);

	initargs(argc, argv);
	// TODO:par read
	if (!getparint("SeisFileNum", &SegyCount))
		SegyCount = 1;
	if (!getparint("nx", &nx))
		nx = 100;
	if (!getparint("nz", &nz))
		nz = 100;
	if (!getparfloat("dx", &dx))
		dx = 10.0;
	if (!getparfloat("dz", &dz))
		dz = 10.0;
	if (!getparfloat("fx_min", &fx_min))
		fx_min = 0.0;
	if (!getparfloat("fx_max", &fx_max))
		fx_max = fx_min + nx * dx;
	if (!getparfloat("dxline", &dxline))
		dxline = 20.0;
	if (!getparfloat("aperx", &aperx))
		aperx = 5000.0;
	if (!getparfloat("fmin", &fmin))
		fmin = 0.080 / dt;
	if (!getparfloat("fmax", &fmax))
		fmax = 5.0 * fmin;
	if (!getparfloat("amax", &amax))
		amax = 46.0;
	amin = -amax;
	if (!getparint("nangle", &nangle))
		nangle = 46;
	if (!getparfloat("zs", &zs))
		zs = 0.;
	if (!getparfloat("zr", &zr))
		zr = 0.;
	if (!getparint("irug", &irug))
		irug = 0;
	if (!getparint("beamDirection", &beamDirection))
		beamDirection = 0;
	if (!getparstring("elev", &elevfile))
		elevfile = "elev";

	if (!getparint("direction", &direction))
		direction = 0;
	parBcast *pbc;
	pbc = (parBcast *)malloc(sizeof(parBcast));
	MPI_Bcast(pbc, sizeof(parBcast), MPI_BYTE, MASTER, MPI_COMM_WORLD);
	parBcastGet(&nx, &nz, &dx, &dz, &fx, &fz,
				&vpavg, &vpmin, &nt, &dt, &nstotal, &ntrps_max, pbc);
	free(pbc);
	printf("rank::%d nx=%d nz=%d\n", rank, nx, nz);
	printf("rank::%d dx=%f dz=%f\n", rank, dx, dz);
	printf("rank::%d fx=%f fz=%f\n", rank, fx, fz);
	printf("rank::%d vpavg=%f vpmin=%f\n", rank, vpavg, vpmin);
	printf("rank::%d nt=%d dt=%f nstotal=%d ntrps_max=%d\n", rank, nt, dt, nstotal, ntrps_max);

	sprintf(imgFileHandle, "./image/imgFile_%d", rank);
	sprintf(angFileHandle, "./ang/angFile_%d", rank);
	// TODO:save file
	imgFile = fopen(imgFileHandle, "wb");
	angfp = fopen(angFileHandle, "wb");

	g = alloc2float(nz, nx);
	finalg = alloc2float(nz, nx);
	tang = alloc2float(nz, nangle);
	fang = alloc2float(nz, nangle);

	if (irug == 1)
	{
		elev = alloc1float(nx);
		elevfp = fopen(elevfile, "rb");
		fread(elev, sizeof(float), nx, elevfp);
		fclose(elevfp);
	}

	for (ix = 0; ix < nx; ix++) //  zero the migrated section
		for (iz = 0; iz < nz; iz++)
			g[ix][iz] = 0.0, finalg[ix][iz] = 0.0;

	//  zero the migrated section
	for (iangle = 0; iangle < nangle; iangle++)
		for (iz = 0; iz < nz; iz++)
			tang[iangle][iz] = 0.0, fang[iangle][iz] = 0.0;

	Model vm;
	vm.nx = nx;
	vm.nz = nz;
	vm.dx = dx;
	vm.dz = dz;
	vm.fx = fx;
	vm.fz = fz;
	vm.vpavg = vpavg;
	vm.vpmin = vpmin;
	vm.vp = alloc2float(nz, nx);
	for (i = 0; i < nx; i++)
	{
		MPI_Bcast(&(vm.vp[i][0]), nz,
				  MPI_FLOAT, MASTER, MPI_COMM_WORLD);
	}

	Shot ss; // receiver coordinate
	ss.nt = nt;
	ss.dt = dt;
	ss.irug = irug;
	ss.zs = zs;
	ss.zr = zr;
	ss.coordx = alloc1float(ntrps_max);
	ss.coordy = alloc1float(ntrps_max);

	float *originCoordx;
	float *originCoordy;
	int *reCoordx, reXs;
	int *reCoordy, reYs;

	originCoordx = alloc1float(ntrps_max);
	originCoordy = alloc1float(ntrps_max);
	reCoordx = alloc1int(ntrps_max);
	reCoordy = alloc1int(ntrps_max);

	if (irug)
		ss.coordz = alloc1float(ntrps_max);
	else
		ss.coordz = NULL;
	ss.data = alloc2float(nt, ntrps_max);

	bwh = 1.0 * vpavg / fmin;

	/* random numbers used to denote live and dead cells */
	live = 98765;
	dead = 1 + (int)(1.0e7 * franuni());
	parSendRecv psr; // shot coordinate

	for (Count = 0; Count < SegyCount; Count++)
	{

		if (Count != 0)
		{
			parBcast *pbc;
			pbc = (parBcast *)malloc(sizeof(parBcast));
			MPI_Bcast(pbc, sizeof(parBcast), MPI_BYTE, MASTER, MPI_COMM_WORLD);
			parBcastGet(&nx, &nz, &dx, &dz, &fx, &fz,
						&vpavg, &vpmin, &nt, &dt, &nstotal, &ntrps_max, pbc);
			free(pbc);
			printf("rank::%d nx=%d nz=%d\n", rank, nx, nz);
			printf("rank::%d dx=%f dz=%f\n", rank, dx, dz);
			printf("rank::%d fx=%f fz=%f\n", rank, fx, fz);
			printf("rank::%d vpavg=%f vpmin=%f\n", rank, vpavg, vpmin);
			printf("rank::%d nt=%d dt=%f nstotal=%d ntrps_max=%d\n", rank, nt, dt, nstotal, ntrps_max);
		}

		for (;;)
		{
			
			MPI_Send(&rank, 1, MPI_INT, MASTER, REQUEST, MPI_COMM_WORLD);
			MPI_Recv(&psr, sizeof(parSendRecv), MPI_BYTE, MASTER,
					 PARAM, MPI_COMM_WORLD, &status);

			if (psr.ntr < 1)
				break; /* job alloc finished */

			ss.is = is = psr.is;
			ss.ntr = ntr = psr.ntr;
			ss.xs = xs = psr.xs;
			ss.zs = zs = psr.zs;
			sx=xs;
			MPI_Recv(&(ss.coordx[0]), ntr, MPI_FLOAT, MASTER,
					 COORDX, MPI_COMM_WORLD, &status);

			MPI_Recv(&(ss.coordy[0]), ntr, MPI_FLOAT, MASTER,
					 COORDY, MPI_COMM_WORLD, &status);
			// 找到ss.coordx的最大值cx_maxf和最小值cx_minf

if(direction==1)
{
	for(i=0;i<ntr;i++)
	{
		ss.coordx[i]=ss.coordy[i];
	}
	sx=psr.ys;
}


			float cx_maxf = ss.coordx[0];
			float cx_minf = ss.coordx[0];
			for (i = 0; i < ntr; i++)
			{
				if (ss.coordx[i] > cx_maxf)
					cx_maxf = ss.coordx[i];
				if (ss.coordx[i] < cx_minf)
					cx_minf = ss.coordx[i];
			}

			cx_max = (int)cx_maxf;
			cx_min = (int)cx_minf;
			printf("rank::%d cx_max=%d cx_min=%d ntr=%d\n", rank, cx_max, cx_min,ntr);

			if (irug)
				MPI_Recv(&(ss.coordz[0]), ntr, MPI_FLOAT, MASTER,
						 COORDZ, MPI_COMM_WORLD, &status);

			MPI_Recv(&(ss.data[0][0]), nt * ntr, MPI_FLOAT, MASTER,
					 DATAZ, MPI_COMM_WORLD, &status);

			apx_min = MIN(sx, cx_min);
			apx_max = MAX(sx, cx_max);

			// aperture of one shot, the imaging area//
			apx_min = MAX(fx_min, apx_min - aperx);
			apx_max = MIN(fx_max, apx_max + aperx);

			nxline = NINT(cx_max - cx_min) / dxline + 1; // the record
			nxline=ntr;
			// nxline = ntr;
			f = alloc2float(nt, nxline);
			for (ix = 0; ix < nxline; ix++)
				for (it = 0; it < nt; it++)
					f[ix][it] = 0.0;
	// printf("cx_min=%d cx_max=%d\n",cx_min,cx_max);
			// 将ss.data[nt][ntr]转换为f[nt][nxline]
			for (ix = 0; ix < ntr; ix++)
			{
				// ixx = NINT((ss.coordx[ix] - cx_min) / dxline);
				for (it = 0; it < nt; it++)
					f[ix][it] = ss.data[ix][it];
			}
			
			
			// printf("rank::%d fx_min=%f fx_max=%f\n",rank,fx_min,fx_max);
			iaperx_min = NINT((apx_min - fx_min) / dx);
			iaperx_max = NINT((apx_max - fx_min) / dx);
			naperx = (iaperx_max - iaperx_min);
			tempv = alloc2float(nz, naperx);
			tempg = alloc2float(nz, naperx);
			elevtemp = alloc1float(naperx);
			coord = alloc1float(ntr);
			fcx_min = cx_min - iaperx_min * dx - fx_min;
			fcx_max = cx_max - iaperx_min * dx - fx_min;
			// printf("iperx_min=%d iaperx_max=%d\n",iaperx_min,iaperx_max);
		// printf("fcx_min=%f fcx_max=%f\n",fcx_min,fcx_max);
			sx = sx - iaperx_min * dx - fx_min;
			// printf("fcx_min=%f fcx_max=%f\n",fcx_min,fcx_max);
			// printf("sx=%f\n",sx);
			// printf("naperx=%d\n",naperx);
			
			for (ix = 0; ix < naperx; ix++)
			{
				// printf("ix+ia=%d,ix=%d\n",ix+iaperx_min,ix);
				if(irug==1)
				{
					elevtemp[ix] = elev[ix + iaperx_min];
				}
				
				for (iz = 0; iz < nz; iz++)
				{
					tempv[ix][iz] = vm.vp[iaperx_min + ix][iz];
					tempg[ix][iz] = 0.0;
				}
			}

			for(it=0;it<ntr;it++)
			{
				coord[it]=ss.coordx[it]-fx_min-iaperx_min * dx;
			}
				

			bwh = 1.0 * vpavg / fmin;

			if (irug == 0)
			{
				csmiggb(-1, bwh, fmin, fmax, amin, amax, live, dead, nt, dt, sx, fcx_min, fcx_max,
						naperx, dx, ntr, dxline, nz, dz, f, tempv, tempg, coord);
			}
			else
			{
				csmiggb_irr(-1, bwh, fmin, fmax, amin, amax, live, dead, nt, dt, sx, fcx_min, fcx_max,
							naperx, dx, ntr, dxline, nz, dz, f, tempv, tempg, coord, elevtemp,elev ,nxline, dxline);
			}

			// TODO:有关每炮地震记录的叠加
			for (ix = 0; ix < naperx; ix++)
				for (iz = 0; iz < nz; iz++)
					g[iaperx_min + ix][iz] += tempg[ix][iz];
printf("rank::%d,shot %d finished\n", rank, is);

			free2float(tempv);
			free2float(tempg);
			free2float(f);
			free1float(elevtemp);
			free1float(coord);
		}

		for (iangle = 0; iangle < nangle; iangle++)
			for (iz = 0; iz < nz; iz++)
				fang[iangle][iz] += tang[iangle][iz];
		for (ix = 0; ix < nx; ix++)
			for (iz = 0; iz < nz; iz++)
				finalg[ix][iz] += g[ix][iz];
		printf("np=%d,job finished\n", rank);
	}

	// 将三维的g写入文件
	for (ix = 0; ix < nx; ix++)
		for (iz = 0; iz < nz; iz++)
			fwrite(&finalg[ix][iz], sizeof(float), 1, imgFile);

	// 将三维的fang写入文件
	for (iangle = 0; iangle < nangle; iangle++)
		for (iz = 0; iz < nz; iz++)
			fwrite(&fang[iangle][iz], sizeof(float), 1, angfp);

	printf("np=%d,file write finished\n", rank);

	free2float(g);
	free2float(finalg);
	free2float(tang);
	free2float(fang);
	fclose(angfp);
	fclose(imgFile);

	// TODO:check free

	return 0;
} // end slave

void parBcastSet(int nx, int nz,
				 float dx, float dz, float fx, float fz,
				 float vpavg, float vpmin,
				 int nt, float dt, int nstotal, int ntrps_max,
				 parBcast *pbc)
{
	pbc->nx = nx;

	pbc->nz = nz;
	pbc->dx = dx;

	pbc->dz = dz;
	pbc->fx = fx;

	pbc->fz = fz;

	pbc->vpavg = vpavg;
	pbc->vpmin = vpmin;
	pbc->nt = nt;
	pbc->dt = dt;
	pbc->nstotal = nstotal;
	pbc->ntrps_max = ntrps_max;
} // end parBcastSet
void parBcastGet(int *nx, int *nz,
				 float *dx, float *dz,
				 float *fx, float *fz,
				 float *vpavg, float *vpmin,
				 int *nt, float *dt, int *nstotal, int *ntrps_max,
				 parBcast *pbc)
{
	*nx = pbc->nx;

	*nz = pbc->nz;
	*dx = pbc->dx;

	*dz = pbc->dz;
	*fx = pbc->fx;

	*fz = pbc->fz;

	*vpavg = pbc->vpavg;
	*vpmin = pbc->vpmin;
	*nt = pbc->nt;
	*dt = pbc->dt;
	*nstotal = pbc->nstotal;
	*ntrps_max = pbc->ntrps_max;
} // end parBcastGet
