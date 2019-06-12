#include <stdio.h>
#include <string.h>
#include <iostream>
#include <sys/stat.h>
#include <unistd.h>
#include <string>
#include <fstream>
#include "main.h"
#include "stdlib.h"
#include "math.h"
#include <cassert>
#include <cmath>
#include <vector>
#include <algorithm>
#include <list>
#include <time.h>
#include "func.h"
//#include <mpi.h>
#include <omp.h>
#include <ctime>
//#include "input.h"

void help(char *infile) {
	fprintf(screen,
		"\nList of command line options which MUST be specified\n\n"
		"start(>=0)	: start at this timestep (-start) (default:%d)\n"
		"stop(>=0)	: stop at this timestep (-stop) (default:%d)\n"
		"output(>=0)	: when to print timestep to screen (-output) (default:%d)\n"
		"stepsize(>=0)	: timestep between data saved in file (-stepsize) (default:%d)\n"
		"LB_start(Ang)	: start of liquid boundary (-LB_start) (default:%f)\n"
		"LB_stop(Ang)	: stop of liquid boundary (-LB_stop) (default:%f)\n"
		"iLB(Ang)	: stepsize of liquid boundary (-iLB) (default:%f)\n"
		"VB_start(Ang)	: start of vapor boundary (-VB_start) (default:%f)\n"
		"VB_start(ang)	: stop of vapor boundary (-VB_stop) (default:%f)\n"
		"iVB(Ang)	: stepsize of vapor boundary (-iVB) (default:%f)\n"
		"vmin(Ang/fs)	: minimum velocity for distribution (-vmin) (default:%f)\n"
		"vmax(Ang/fs)	: maximum velocity for distribution (-vmax) (default:%f)\n"
		"iV		: stepsize for velocity distribution (-iV) (default:%d)\n"
		"string		: append name to generated output files  (-string)\n"
		"vflag		: calculate velocity distribution of the fluxes (-vflag) (default:%d)\n\n",
		start,stop,output,step_size,LB_start,LB_stop,iLB,VB_start,VB_stop,iVB,vmin,vmax,iV,vflag
	);

};

int main(int narg, char **arg)
{

// set initial values
	start = 0;
	stop = 0;
	output = 2000;
	step_size = 500;
	LB_start = 0.0;
	LB_stop = 0.0;
	iLB = 1.0;
	VB_start = 0.0;
	VB_stop = 0.0;
	iVB = 1.0;
	threads = 1;
	vmin = -0.01;
	vmax = -vmin;
	iV = 500;
	data = 1;
	STRflag = 0;
	vflag = 0;
	nflag = 0;
	NATOMS = 0;
	SBYTES = 0;

	lineD = new char[91];
	strcpy(lineD,"==========================================================================================");
	lineS = new char[91];
	strcpy(lineS,"------------------------------------------------------------------------------------------");

	int flag_warning = 0;	
	int iarg = 2;
	if ( strcmp(arg[1],"-help") == 0 || strcmp(arg[1],"-h") == 0) {
		help(arg[1]);
		return 1;
	};


	while (iarg < narg) {
		if ( strcmp(arg[iarg],"-start") == 0 ) {
			if (iarg+2 > narg) {fprintf(screen,"Missing argument value of '%s'\n",arg[iarg]); return 1;}
			double starti = atof(arg[iarg+1]);
			if (starti!=floor(starti)) {
				fprintf(screen,"WARNING: Start value changed to %d!\n",(int)floor(starti));
 				start = (int)floor(starti);
			} else {
				start = (int)starti;
			};
			iarg += 2;
		} else if ( strcmp(arg[iarg],"-stop") == 0 ) {
			if (iarg+2 > narg) {fprintf(screen,"Missing argument value of '%s'\n",arg[iarg]); return 1;}
			double stopi = atof(arg[iarg+1]);
			if (stopi!=floor(stopi)) {
				fprintf(screen,"WARNING: Stop value changed to %d!\n",(int)floor(stopi));
 				stop = (int)floor(stopi);
			} else {
				stop = (int)stopi;
			};
			iarg += 2;
		} else if ( strcmp(arg[iarg],"-output") == 0 ) {
			if (iarg+2 > narg) {fprintf(screen,"Missing argument value of '%s'\n",arg[iarg]); return 1;}
			double outputi = atof(arg[iarg+1]);
			if (outputi!=floor(outputi)) {
				fprintf(screen,"WARNING: Output value changed to %d!\n",(int)floor(outputi));
 				output = (int)floor(outputi);
			} else {
				output = (int)outputi;
			};
			iarg += 2;
		} else if ( strcmp(arg[iarg],"-stepsize") == 0 ) {
			if (iarg+2 > narg) {fprintf(screen,"Missing argument value of '%s'\n",arg[iarg]); return 1;}
			double step_sizei = atof(arg[iarg+1]);
			if (step_sizei!=floor(step_sizei)) {
				fprintf(screen,"WARNING: step_size value changed to %d!\n",(int)floor(step_sizei));
 				step_size = (int)floor(step_sizei);
			} else {
				step_size = (int)step_sizei;
			};
			iarg += 2;
		} else if ( strcmp(arg[iarg],"-output_every") == 0 ) {
			if (iarg+2 > narg) {fprintf(screen,"Missing argument value of '%s'\n",arg[iarg]); return 1;}
			double output_everyi = atof(arg[iarg+1]);
			output_every = (int)output_everyi;
			output_every_flag = 1;
			iarg += 2;
		} else if ( strcmp(arg[iarg],"-LB_start") == 0 ) {
			if (iarg+2 > narg) {fprintf(screen,"Missing argument value of '%s'\n",arg[iarg]); return 1;}
			LB_start = atof(arg[iarg+1]);
			iarg += 2;
		} else if ( strcmp(arg[iarg],"-LB_stop") == 0 ) {
			if (iarg+2 > narg) {fprintf(screen,"Missing argument value of '%s'\n",arg[iarg]); return 1;}
			LB_stop = atof(arg[iarg+1]);
			iarg += 2;
		} else if ( strcmp(arg[iarg],"-iLB") == 0 ) {
			if (iarg+2 > narg) {fprintf(screen,"Missing argument value of '%s'\n",arg[iarg]); return 1;}
			iLB = atof(arg[iarg+1]);
			iarg += 2;
		} else if ( strcmp(arg[iarg],"-VB_start") == 0 ) {
			if (iarg+2 > narg) {fprintf(screen,"Missing argument value of '%s'\n",arg[iarg]); return 1;}
			VB_start = atof(arg[iarg+1]);
			iarg += 2;
		} else if ( strcmp(arg[iarg],"-VB_stop") == 0 ) {
			if (iarg+2 > narg) {fprintf(screen,"Missing argument value of '%s'\n",arg[iarg]); return 1;}
			VB_stop = atof(arg[iarg+1]);
			iarg += 2;
		} else if ( strcmp(arg[iarg],"-iVB") == 0 ) {
			if (iarg+2 > narg) {fprintf(screen,"Missing argument value of '%s'\n",arg[iarg]); return 1;}
			iVB = atof(arg[iarg+1]);
			iarg += 2;
		} else if ( strcmp(arg[iarg],"-vmin") == 0 ) {
			if (iarg+2 > narg) {fprintf(screen,"Missing argument value of '%s'\n",arg[iarg]); return 1;}
			vmin = atof(arg[iarg+1]);
			iarg += 2;
		} else if ( strcmp(arg[iarg],"-vmax") == 0 ) {
			if (iarg+2 > narg) {fprintf(screen,"Missing argument value of '%s'\n",arg[iarg]); return 1;}
			vmax = atof(arg[iarg+1]);
			iarg += 2;
		} else if ( strcmp(arg[iarg],"-iV") == 0 ) {
			if (iarg+2 > narg) {fprintf(screen,"Missing argument value of '%s'\n",arg[iarg]); return 1;}
			iV = atoi(arg[iarg+1]);
			iarg += 2;
		} else if ( strcmp(arg[iarg],"-string") == 0 ) {
			if (iarg+2 > narg) {fprintf(screen,"Missing argument value of '%s'\n",arg[iarg]); return 1;}
			STR = new char[strlen(arg[iarg+1])+1];
			strcpy(STR,arg[iarg+1]);
			STRflag = 1;
			iarg += 2;
		} else if ( strcmp(arg[iarg],"-data") == 0 ) {
			if (iarg+2 > narg) {fprintf(screen,"Missing argument value of '%s'\n",arg[iarg]); return 1;}
			data = atoi(arg[iarg+1]);
			iarg += 2;
		} else if ( strcmp(arg[iarg],"-vflag") == 0 ) {
			if (iarg+2 > narg) {fprintf(screen,"Missing argument value of '%s'\n",arg[iarg]); return 1;}
			vflag = atoi(arg[iarg+1]);
			iarg += 2;
		} else if ( strcmp(arg[iarg],"-natoms") == 0 ) {
			if (iarg+2 > narg) {fprintf(screen,"Missing argument value of '%s'\n",arg[iarg]); return 1;}
			nflag = 1;
			NATOMS = atoi(arg[iarg+1]);
			SBYTES = NATOMS;
			iarg += 2;
//		} else if ( strcmp(arg[iarg],"-threads") == 0 ) {
//			if (iarg+2 > narg) {fprintf(screen,"Missing argument value of '%s'\n",arg[iarg]); return 1;}
//			threads = atoi(arg[iarg+1]);
//			iarg += 2;
		} else if ( strcmp(arg[iarg],"-help") == 0 || strcmp(arg[iarg],"-h") == 0 ) {
			help(arg[1]);
			return 1;
		} else { fprintf(screen,"Illegal name: '%s'\n",arg[iarg]); return 1;} 
	}


	if (!output_every_flag)	output_every = step_size * output;



	time_t now = time(0);
	tm *ltm = localtime(&now);
	int year = 1900 + ltm->tm_year;
	int month = 1 + ltm->tm_mon;
	int day = ltm->tm_mday;
	int hour = ltm->tm_hour;
	int minute = ltm->tm_min;
	int second = ltm->tm_sec;


	fprintf(screen,"\n%s\n\n",lineD);
	fprintf(screen,"%-37s%s","","POST-PROCESSING");
	fprintf(screen,"\n\n%s\n",lineD);
	fprintf(screen,"Date: %02d/%02d/%04d (d/m/y)\n",day,month,year);
	fprintf(screen,"Time: %02d:%02d:%02d (hr:min:sec)",hour,minute,second);
	fprintf(screen,"\n%s\n",lineD);

	clock_t starttime, stoptime;
	double totaltime;
	starttime = clock();

	int RUN = 0;


	if (data==1) {
		imass=0; imol=1; iid=2; itype=3; ix=4; iy=5; iz=6; ivx=7; ivy=8; ivz=9; ike=10; ipe=11; is1=12; is2=13; is3=14; is4=15; is5=16; is6=17;
	} else if (data==2) {
		imass=0; iid=1; ix=2; iy=3; iz=4;
	} else {

/*
		char *str = new char[4];
		fprintf(screen,"Number of atom data to define:\n");
		fgets(str,4,stdin);
		fprintf(screen,"Define atom data e.g. 'x', 'type', ...\n");
		char **strdata = new char*[6];
	//	fgets(*strdata,100,stdin);

		for (int i = 0; i < atoi(str); i++) {
			fprintf(screen,"%d: ",i);
			fgets(strdata[i],10,stdin);
		};
	
		for (int i = 0; i < atoi(str); i++) {
			fprintf(screen,"%s\n",strdata[i]);
		};


		delete [] str;
		delete [] strdata;
*/
		return 1;		
	
	};

	FILE *fp = fopen(arg[1],"rb");
	if (fp == NULL) {
		perror("Error opening file");
		return 1;
	};

	fprintf(screen,"\n%s\n",lineD);
	fprintf(screen,"ATOM DATA");
	fprintf(screen,"\n%s\n",lineS);
	fprintf(screen,"mass:%d - mol:%d - id:%d - type:%d\n",imass,imol,iid,itype);
	fprintf(screen,"x:%d - y:%d - z:%d\n",ix,iy,iz);
	fprintf(screen,"vx:%d - vy:%d - vz:%d - ke:%d - pe:%d\n",ivx,ivy,ivz,ike,ipe);
	fprintf(screen,"s1:%d - s2:%d - s3:%d - s4:%d - s5:%d - s6:%d\n",is1,is2,is3,is4,is5,is6);
	fprintf(screen,"\n%s\n",lineD);

	fprintf(screen,"\n%s\n",lineD);
	fprintf(screen,"INPUT");
	fprintf(screen,"\n%s\n",lineS);

	fprintf(screen,"%-25s%-25s%s\n","NAME","COMMAND","VALUE");
	fprintf(screen,"%-25s%-25s%d\n","Start","-start",start);
	fprintf(screen,"%-25s%-25s%d\n","Stop","-stop",stop);
	fprintf(screen,"%-25s%-25s%d\n","Output","-output",output);
	fprintf(screen,"%-25s%-25s%d\n","Step size","-stepsize",step_size);
	fprintf(screen,"%-25s%-25s%d\n","Output every","",output_every);
	fprintf(screen,"%-25s%-25s%.4f\n","Liquid bnd start","-LB_start",LB_start);
	fprintf(screen,"%-25s%-25s%.4f\n","Liquid bnd stop","-LB_stop",LB_stop);
	fprintf(screen,"%-25s%-25s%.4f\n","Liquid bnd step","-iLB",iLB);
	fprintf(screen,"%-25s%-25s%.4f\n","Vapor bnd start","-LB_start",VB_start);
	fprintf(screen,"%-25s%-25s%.4f\n","Vapor bnd stop","-LB_stop",VB_stop);
	fprintf(screen,"%-25s%-25s%.4f\n","Vapor bnd step","-iVB",iVB);
	fprintf(screen,"%-25s%-25s%s\n","String","-string",STR);
	fprintf(screen,"%-25s%-25s%d\n","Data","-data",data);
	fprintf(screen,"%-25s%-25s%d\n","vflag","-vflag",vflag);
	if (vflag) fprintf(screen,"%-25s%-25s%f\n","vmin","-vmin",vmin);
	if (vflag) fprintf(screen,"%-25s%-25s%f\n","vmax","-vmax",vmax);
	if (vflag) fprintf(screen,"%-25s%-25s%d\n","Velocity step","-iV",iV);


	fprintf(screen,"%s\n",lineD);

	int step = 0;

	if (!nflag) {

	fprintf(screen,"\nFind largest molecule ID and number of atoms\nTime steps: ");

	while(1) {

		fread(&ntimestep,sizeof(bigint),1,fp);	
	
		if (feof(fp)) break;

		if ((step%output_every)==0) fprintf(screen," " BIGINT_FORMAT, ntimestep);
		fflush(stdout);		
		step += step_size;
			
		fread(&natoms,sizeof(bigint),1,fp);		
		fread(&triclinic,sizeof(int),1,fp);		
		fread(&boundary[0][0],6*sizeof(int),1,fp);	
    fread(&xlo,sizeof(double),1,fp);
    fread(&xhi,sizeof(double),1,fp);
 	  fread(&ylo,sizeof(double),1,fp);
    fread(&yhi,sizeof(double),1,fp);
    fread(&zlo,sizeof(double),1,fp);
    fread(&zhi,sizeof(double),1,fp);
    if (triclinic) {
			fread(&xy,sizeof(double),1,fp);
			fread(&xz,sizeof(double),1,fp);
			fread(&yz,sizeof(double),1,fp);
    }
    fread(&size_one,sizeof(int),1,fp);
    fread(&nchunk,sizeof(int),1,fp);
	
		if (natoms > NATOMS) NATOMS = natoms;

		if (ntimestep >= start && ntimestep <= stop) {
		
		for (int i = 0; i < nchunk; i++) {
			fread(&n,sizeof(int),1,fp);				

			if (n > maxbuf) {
				if (buf) delete [] buf;
				buf = new double[n];
				maxbuf = n;
			}

			fread(buf,sizeof(double),n,fp);
				
			for (int S = 0; S < n/size_one; S++) {
				if (buf[iid+size_one*S]>SBYTES) SBYTES = buf[iid+size_one*S];
			};

		};
	

		} else if (ntimestep > stop) {

			break;

		} else {

			for (int i = 0; i < nchunk; i++) {
				fread(&n,sizeof(int),1,fp);				

				if (n > maxbuf) {
					if (buf) delete [] buf;
					buf = new double[n];
					maxbuf = n;
				}

				fread(buf,sizeof(double),n,fp);

			}; //end of chunk
		};

	};//end of while
	};

	fprintf(screen,"\nLargest molecule ID: %d\n",SBYTES);
	fprintf(screen,"Number of atoms: %d\n\n",NATOMS);

//	int NATOMS = SBYTES;	
//	SBYTES = 10000000;

	SBYTES += 1;
	NATOMS += 1;
	
	indexID = new double[NATOMS];
	indices = new int[SBYTES];
	JevapL_oldlist = new int[SBYTES];	JevapL_newlist = new int[SBYTES];
	JcondL_oldlist = new int[SBYTES];	JcondL_newlist = new int[SBYTES];
	JcollL_oldlist = new int[SBYTES];	JcollL_newlist = new int[SBYTES];
	JoutL_oldlist  = new int[SBYTES];	JoutL_newlist  = new int[SBYTES];

	int MAX = 0;
	double max_evap = 0.0;
	int count_evap = 0;


	double Liquid[4] = {-1000.0, 20.0, 730.0, 1000.0}; //set particles in this region to liquidflag

	if (!STRflag) {
			STR = new char[1];
			strcpy(STR,"");
	};

	char *filename = new char[40+strlen(STR)+1];
	sprintf(filename,"Output_evap_cond_%02d_%02d_%04d%s.txt",day,month,year,STR);

	if (fexists(filename)) sprintf(filename,"Output_evap_cond_%02d_%02d_%04d_%02d%02d%02d%s.txt",day,month,year,hour,minute,second,STR);

	FILE *fpdatcoeff = fopen(filename,"w");

	char *filename1 = new char[40+strlen(STR)+1];
	sprintf(filename1,"Output_flux_time_%02d_%02d_%04d%s.txt",day,month,year,STR);

	if (fexists(filename1)) sprintf(filename1,"Output_flux_time_%02d_%02d_%04d_%02d%02d%02d%s.txt",day,month,year,hour,minute,second,STR);

	FILE *fpfluxtime = fopen(filename1,"w");

	fprintf(screen,"Created output file(s):\n");
	fprintf(screen,"> %s\n",filename);
	fprintf(screen,"> %s\n",filename1);

	delete [] filename;
	delete [] filename1;

	fprintf(fpdatcoeff,"Liquid_bnd Vapor_bnd L_Jevap L_Jout L_Jevap/Jout L_Jcond L_Jcoll L_Jcond/Jcoll\n");

	fseek(fp,0,SEEK_SET);

	JL = new int[4];

	for (double iLB_start = LB_start; iLB_start <= LB_stop; iLB_start += iLB) {

		int Vel_size = 3 * iV;
//		int Vel_size = 3 * iV * (ceil((VB_stop-iLB_start)/iVB)+1);
//		fprintf(screen,"size: %d\n",Vel_size);

		if (vel_Jout) delete [] vel_Jout;
		if (vel_Jcoll) delete [] vel_Jcoll;
		if (vel_Jevap) delete [] vel_Jevap;
		if (vel_Jcond) delete [] vel_Jcond;

		vel_Jout = new int[Vel_size];
		vel_Jcoll = new int[Vel_size];
		vel_Jevap = new int[Vel_size];
		vel_Jcond = new int[Vel_size];

//		for (int ivel = 0; ivel < Vel_size; ivel++) {
//			vel_Jout[ivel] = 0;
//			vel_Jcoll[ivel] = 0;
//			vel_Jevap[ivel] = 0;
//			vel_Jcond[ivel] = 0;
//		}

	//	VB_tic = 0;

// Create file for velocity distribution. New file will be created for each liquid boundary position (LB).
		if (vflag) {	

			char *filename = new char[40+strlen(STR)+1];
			sprintf(filename,"vel_evap_cond_%.2f_%02d_%02d_%04d%s.txt",iLB_start,day,month,year,STR);

			if (fexists(filename)) sprintf(filename,"vel_evap_cond_%.2f_%02d_%02d_%04d_%02d%02d%02d%s.txt",iLB_start,day,month,year,hour,minute,second,STR);
	
//			fprintf(screen,"Created output file(s):\n");
			fprintf(screen,"> %s\n\n",filename);
	
			fpdatvel = fopen(filename,"w");
			delete [] filename;

			fprintf(fpdatvel,"vmin vmax iV\n");
			fprintf(fpdatvel,"%E %E %d\n",vmin,vmax,iV);
			fprintf(fpdatvel,"LB VB Jout Jcoll Jevap Jcond\n");

		};


		for (double iVB_start = VB_start; iVB_start <= VB_stop; iVB_start += iVB) {

			fseek(fp, 0,SEEK_SET);

			for (int ivel = 0; ivel < Vel_size; ivel++) {
				vel_Jout[ivel] = 0;
				vel_Jcoll[ivel] = 0;
				vel_Jevap[ivel] = 0;
				vel_Jcond[ivel] = 0;
			};

			for (int t = 0; t < SBYTES; t++) {
				indices[t] = 0;
				JevapL_oldlist[t] = 0;	JevapL_newlist[t] = 0;
				JcondL_oldlist[t] = 0;	JcondL_newlist[t] = 0;
				JcollL_oldlist[t] = 0;	JcollL_newlist[t] = 0;
				JoutL_oldlist[t] = 0;  	JoutL_newlist[t] = 0;
			};

			JL[0] = 0; JL[1] = 0; JL[2] = 0;JL[3] = 0;

//		double JL1 = 0.0;
//		double JL2 = 0.0;
//		double JL3 = 0.0;
//		double JL4 = 0.0;

		double Left_boundary[4] = {iLB_start, iVB_start, iVB_start, iVB_start}; //Jcond,Jevap,Jcoll,Jout
		step = 0;

		RUN++;
		fprintf(screen,"\n%s\n",lineD);
		fprintf(screen,"%-42sRUN %d","",RUN);
//		fprintf(screen,"\n%s\n",lineS);
//		fprintf(screen,"Liquid-boundary: %E Ang -- Vapor-boundary: %E Ang\n",Left_boundary[0],Left_boundary[1]);
		fprintf(screen,"\n%s\n",lineS);
		fprintf(screen,"Time steps:");

		while(1) {

			fread(&ntimestep,sizeof(bigint),1,fp);	
	
			if (feof(fp)) {
					fprintf(screen,"\n%s\n",lineS);
					break;
			}
	
			if ((step%output_every)==0) fprintf(screen," " BIGINT_FORMAT, ntimestep);
			step += step_size;

			fflush(stdout);
			
			fread(&natoms,sizeof(bigint),1,fp);		
			fread(&triclinic,sizeof(int),1,fp);		
			fread(&boundary[0][0],6*sizeof(int),1,fp);	
      fread(&xlo,sizeof(double),1,fp);
      fread(&xhi,sizeof(double),1,fp);
  	  fread(&ylo,sizeof(double),1,fp);
	    fread(&yhi,sizeof(double),1,fp);
      fread(&zlo,sizeof(double),1,fp);
      fread(&zhi,sizeof(double),1,fp);
      if (triclinic) {
				fread(&xy,sizeof(double),1,fp);
				fread(&xz,sizeof(double),1,fp);
				fread(&yz,sizeof(double),1,fp);
      }
      fread(&size_one,sizeof(int),1,fp);
      fread(&nchunk,sizeof(int),1,fp);

	
			if (ntimestep <= start) {
				int rowID = 0;

				for (int i = 0; i < nchunk; i++) {
					fread(&n,sizeof(int),1,fp);				

					if (n > maxbuf) {
						if (buf) delete [] buf;
						buf = new double[n];
						maxbuf = n;
					}

					fread(buf,sizeof(double),n,fp);
				
						for (int S = 0; S < n/size_one; S++) {

							int row = buf[iid+size_one*S];
							if (buf[iz+size_one*S] < Left_boundary[0]) indices[row] = -2; //liquid region
							if (buf[iz+size_one*S] > Left_boundary[1]) indices[row] = -1; //vapor region

						};

				}; // end of chunk
				
			}	else if (ntimestep > start && ntimestep <= stop) { 
			
				int rowID = 0;

				for (int i = 0; i < nchunk; i++) {
					fread(&n,sizeof(int),1,fp);

					if (n > maxbuf) {
						if (buf) delete [] buf;
						buf = new double[n];
						maxbuf = n;
					}

					fread(buf,sizeof(double),n,fp);	

						int id;
						double z;

//						int JL_prv[4] = {0,0,0,0};
//						int JR_prv[4] = {0,0,0,0};

						for (int S = 0; S < n/size_one; S++) {

							indexID[rowID] = buf[iid+size_one*S];

							id = buf[iid+size_one*S];
							z = buf[iz+size_one*S];
							
							flux(S,Left_boundary);

							if (z < Left_boundary[0]) indices[id] = -2;//liquid 
							if (z > Left_boundary[1]) indices[id] = -1;//vapor
			
							rowID += 1;

						}; // end of for loop

//						JL[0] += JL_prv[0];
//						JL[1] += JL_prv[1];
//						JL[2] += JL_prv[2];
//						JL[3] += JL_prv[3];

				}; //end of chunk
				//																											time step,LB start, VB start, Jcond, Jevap, Jcoll, Jout 
//				fprintf(fpfluxtime,BIGINT_FORMAT " %E %E %E %E %E %E\n",ntimestep,iLB_start,iVB_start,JL[0]-JL1,JL[1]-JL2,JL[2]-JL3,JL[3]-JL4);
//				JL1=JL[0];JL2=JL[1];JL3=JL[2];JL4=JL[3];


				for (int t1 = 0; t1 < natoms; t1++) {
					int t = indexID[t1];
					JevapL_oldlist[t] = JevapL_newlist[t];
					JevapL_newlist[t] = 0;
					JcondL_oldlist[t] = JcondL_newlist[t];
					JcondL_newlist[t] = 0;
					JcollL_oldlist[t] = JcollL_newlist[t];
					JcollL_newlist[t] = 0;
					JoutL_oldlist[t] = JoutL_newlist[t];
					JoutL_newlist[t] = 0;
				};

//				};

			} else {

				fprintf(screen,"\n%s\n",lineS);
				break; 

			};
	
		}; // end of while
		fprintf(screen,"%-13s%-13s%-10s%-10s%-13s%-10s%-10s%s\n","Liquid bnd","Vapor bnd","Jevap","Jout","Jevap/Jout","Jcond","Jcoll","Jcond/Jcoll");
		fprintf(screen,"%-13.4E%-13.4E%-10d%-10d%-13.4E%-10d%-10d%-10.4E\n",Left_boundary[0],Left_boundary[1],JL[1],JL[3],JL[1]/(double)JL[3],JL[0],JL[2],JL[0]/(double)JL[2]);

		// liquid, vapor, Jevap, Jout, Jevap/Jout, Jcond, Jcoll, Jcond/Jcoll
		fprintf(fpdatcoeff,"%E %E %d %d %E %d %d %E\n",Left_boundary[0],Left_boundary[1],JL[1],JL[3],JL[1]/(double)JL[3],JL[0],JL[2],JL[0]/(double)JL[2]);

		fflush(fpdatcoeff);

		totaltime = (clock()-starttime)/(double)CLOCKS_PER_SEC;
		int hours = floor(totaltime/3600.0);
		int minutes = floor((totaltime - (hours*3600))/60.0);
		double seconds = totaltime - (hours*3600.0) - (minutes*60.0);

		fprintf(screen,"\n%s",lineS);
		fprintf(screen,"\nTotal time: %f seconds - [%02d:%02d:%f (hr:min:sec)]",totaltime,hours,minutes,seconds);
		fprintf(screen,"\n%s\n\n",lineD);

//		VB_tic += 1;

		if (vflag) {	

//			for (int iv1 = 0; iv1 < VB_tic; iv1++ ) {
				for (int iv2 = 0; iv2 < iV; iv2++) {
					fprintf(fpdatvel,"%E %E %d %d %d %d %d %d %d %d %d %d %d %d\n",iLB_start, iVB_start,
					vel_Jout[iv2],vel_Jout[iv2+iV],vel_Jout[iv2+2*iV],
					vel_Jcoll[iv2],vel_Jcoll[iv2+iV],vel_Jcoll[iv2+2*iV],
					vel_Jevap[iv2],vel_Jevap[iv2+iV],vel_Jevap[iv2+2*iV],
					vel_Jcond[iv2],vel_Jcond[iv2+iV],vel_Jcond[iv2+2*iV]
					);
				};
//			};

		};

		fflush(fpdatvel);

  	}; //end for loop VB_start

		if (vflag) fclose(fpdatvel); // Close velocity distribution file

	}; // end for loop LB_start


	if (buf) delete [] buf;
	if (vel_Jout) delete [] vel_Jout;
	if (vel_Jcoll) delete [] vel_Jcoll;
	if (vel_Jevap) delete [] vel_Jevap;
	if (vel_Jcond) delete [] vel_Jcond;

	delete [] indices;
	delete [] indexID;
	delete []	JevapL_oldlist;
	delete [] JevapL_newlist;
	delete [] JcondL_oldlist;	
	delete [] JcondL_newlist;
	delete []	JcollL_oldlist;
	delete []	JcollL_newlist;
	delete []	JoutL_oldlist;
	delete []	JoutL_newlist;
	delete [] STR;

	delete [] JL;

	fclose(fpdatcoeff);
	fclose(fpfluxtime);

	return 0;

} // end of main


