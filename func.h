//===============flux=========================
void flux(int S, double *Left_bnd)
{
	int iVX,iVY,iVZ;

	int id = buf[iid+size_one*S];
	double z = buf[iz+size_one*S];
	double dz = 20.0;
	
	double z_left[4] = {Left_bnd[0] + 10.0, Left_bnd[1] - 10.0, 
											Left_bnd[2] + 10.0, Left_bnd[3] - 10.0};

//	int pos;
	int sign_condL, sign_collL, sign_evapL, sign_outL;
//	int sign_condR, sign_collR, sign_evapR, sign_outR;

	sign_outL  = (z > Left_bnd[3]) - 2; // sign -1 or -2
	sign_collL = (z > Left_bnd[2]) - 2; // sign
//	sign_evapL = (z > Left_bnd[1]) - 2; // sign
	sign_condL = (z > Left_bnd[0]) - 2; // sign

	if (z>Left_bnd[0]-20 && z<Left_bnd[1]+20) fprintf(fpposition,"%d %E\n",id,z);


	if (vflag) {
		double VX = buf[ivx+size_one*S];
		double VY = buf[ivy+size_one*S];
		double VZ = buf[ivz+size_one*S];

		double iivx = floor((VX-vmin)*iV/(vmax-vmin));
		double iivy = floor((VY-vmin)*iV/(vmax-vmin));
		double iivz = floor((VZ-vmin)*iV/(vmax-vmin));

		iVX = (int)sqrt(iivx*iivx);
		iVY = (int)sqrt(iivy*iivy);
		iVZ = (int)sqrt(iivz*iivz);
	};

	//====Jout====//
	if (JoutL_oldlist[id] != 0) {
		if (JoutL_oldlist[id] == sign_outL) { //if old sign == new sign => atom did not cross fzhi
//			if (z <= Left_bnd[3]) {
				JoutL_newlist[id] = sign_outL;		
//			}
		} else {
//		if (JoutL_oldlist[id] != sign_outL) {
				if (z < (Left_bnd[3]+dz)) {
				JL[3]++;
				fprintf(fpout,"%d\n",id);
				if (vflag) {
					vel_Jout[iVX] += 1.0;
					vel_Jout[iVY+iV] += 1.0;
					vel_Jout[iVZ+2*iV] += 1.0;
				};
//				JoutL_newlist[id] = 0;		
				if (indices[id] == -2) { //flag[pos_flag+1] == -2) {
					JL[1]++; // old sign /= new sign
					fprintf(fpevap,"%d\n",id);
					if (vflag) {
						vel_Jevap[iVX] += 1.0;
						vel_Jevap[iVY+iV] += 1.0;
						vel_Jevap[iVZ+2*iV] += 1.0;
					};
				};
				};
		}
	} else if (z < Left_bnd[3]) {
		JoutL_newlist[id] = sign_outL;		
	}

/*
	//====Jevap====//
	if (JevapL_oldlist[id] != 0) {
		if (JevapL_oldlist[id] == sign_evapL) { //if old sign == new sign => atom did not cross fzhi
//			if (z <= Left_bnd[1]) {
					JevapL_newlist[id] = sign_evapL;
//			}
		} else {
//			if (JevapL_oldlist[id] != sign_evapL) {
			if (indices[id] == -2) { //flag[pos_flag+1] == -2) {
				JL[1] += 1; // old sign /= new sign
				vel_Jevap[iVX+VB_tic*3*iV] += 1.0;
				vel_Jevap[iVY+iV+VB_tic*3*iV] += 1.0;
				vel_Jevap[iVZ+2*iV+VB_tic*3*iV] += 1.0;
				JevapL_newlist[id] = 0;
			};
		}
	} else if (z <= Left_bnd[1]) {
			JevapL_newlist[id] = sign_evapL;
	};
*/

	//====Jcoll====//
	if (JcollL_oldlist[id] != 0) {
		if (JcollL_oldlist[id] == sign_collL) { //if old sign == new sign => atom did not cross fzhi
//			if (z >= Left_bnd[2] && z <= z_left[2]) {
				JcollL_newlist[id] = sign_collL;		
//			}
	} else {
			if (z>(Left_bnd[2]-dz)) {
			JL[2]++;
			fprintf(fpcoll,"%d\n",id);
			//fprintf(fpposition,"%d %E\n",id,z);
			if (vflag) {
				vel_Jcoll[iVX] += 1.0;
				vel_Jcoll[iVY+iV] += 1.0;
				vel_Jcoll[iVZ+2*iV] += 1.0;
			};
			};
//			JcollL_newlist[id] = 0;		
		}
	} else if (z >= Left_bnd[2]) {
			JcollL_newlist[id] = sign_collL;		
	}



	//====Jcond====//
	if (JcondL_oldlist[id] != 0) {
			if (JcondL_oldlist[id] == sign_condL) { //if old sign == new sign => atom did not cross fzhi
//				if (z >= Left_bnd[0] && z <= z_left[0]) {
					JcondL_newlist[id] = sign_condL;
//			}
		} else {
			if (z>(Left_bnd[0]-dz)) {
			if (indices[id] == -1) { //flag[pos_flag+1] == -1) {
				JL[0]++; // old sign /= new sign
				fprintf(fpcond,"%d\n",id);
				if (vflag) {
					vel_Jcond[iVX] += 1.0;
					vel_Jcond[iVY+iV] += 1.0;
					vel_Jcond[iVZ+2*iV] += 1.0;
				};
//				JcondL_newlist[id] = 0;
			};
			};
		}
	} else if (z >= Left_bnd[0]) {
			JcondL_newlist[id] = sign_condL;
	}

}
