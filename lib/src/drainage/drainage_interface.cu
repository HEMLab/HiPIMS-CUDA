// ======================================================================================
// Name                :    Drainage model (first version)
// Description         :    This drainage model is developed by coupling TPA and 2D SWEs
// ======================================================================================
// Version             :    1.00
// Authors             :    Qian Li                                                   
//                           PhD candidate in Newcastle University
// Create Time         :    2018/4/30
// Update Time         :    
// ======================================================================================
// Copyright @	Qian Li . All rights reserved.
// ======================================================================================


#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "drainage_class.h"
#include "drainage_functions.h"
#include "drainage_interface.h"
#include <cuda.h>
#include <iostream>
#include <cstdio>
#include <fstream>
#include <cmath>
#include <string>
#include <sstream>
#include <map>
#include <thrust/device_ptr.h>
#include <thrust/reduce.h>
#include "mapped_field.h"
#include "cuda_mapped_field.h"




namespace GC {

	//--------global functions list------------------------------------------------------------------------------------------------------////////////////////////////
	__global__ void surfDrainQKernal_calculation( int num_Junc,
		 double g,
		 double pi,
		 double dT,
		 double dx,
		 double *Qe,
		double *Qdrainage,
		double *juncMaxDept,
		double *juncR,
		double *juncArea,
		double *hj,
		double *zbJunc,
		int *index,
		double *zbSurf,
		double *hs,
		double *QSurf);
		__global__ void surfH_limilatorKernal(unsigned int num_Girds, double *h_surf_device);

	  __device__ void QeLimitor(double &QeD, double dx, double hs, double juncArea, double hj, double dT, double maxDept);






	//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	void surfDrainQ_calculation(double dT, Input_information inputKeyPtr, Junc_variables &juncVariables,
		Const_variables constVariables, cuFvMappedField<Scalar, on_cell>& z, cuFvMappedField<Scalar, on_cell>& h, cuFvMappedField<Scalar, on_cell>& Q_surface,
		int blocksPerGrid, int threadsPerBlock) {
		auto mesh = h.mesh;
		auto cell_volumes = mesh->cell_volumes.host_ptr();
		double dx_surface = sqrt(cell_volumes[0]);

		surfDrainQKernal_calculation << <blocksPerGrid, threadsPerBlock >> > (inputKeyPtr.num_Junc,
			constVariables.g,
			constVariables.pi,
			dT,
			dx_surface,
			juncVariables.Qe_device,
			juncVariables.Qexchange_device,
			juncVariables.maxDepth_device,
			juncVariables.j_radius_device,
			juncVariables.junc_Area_device,
			juncVariables.h_j_device,
			juncVariables.zb_j_device,
			juncVariables.surfGrid_device,
			z.data.dev_ptr(),
			h.data.dev_ptr(),
			Q_surface.data.dev_ptr());
	}



	

	//----------------------------------------------------------------------------------------------------------------------
	__global__ void surfDrainQKernal_calculation( int num_Junc,
		 double g,
		 double pi,
		 double dT,
		 double dx,
		 double *Qe,
		double *Qdrainage,
		double *juncMaxDept,
		double *juncR,
		double *juncArea,
		double *hj,
		double *zbJunc,
		int *index,
		double *zbSurf,
		double *hs,
		double *QSurf) {

		// hs : water depth on the surface ground
		// hj: water depth of the junction in the drainage system
		// ind: the corresponding surface grid index
		// Q_s: exchange flow rate for the surface ground
		// Q_d: exchange flow rate for the drainage system


		int tid = blockDim.x * blockIdx.x + threadIdx.x;
		double kesi=1e-10;
		double Ci = 1.0;
	



		double tempSingleQ;
	    double tempAllQ[3];

          int J1_allIndex[25]={3255, 3256, 3257, 3258, 3259,
	  	                  3154, 3155, 3156, 3157, 3158,
	  	                  3053, 3054, 3055, 3056, 3057,
	  	                  2952, 2953, 2954, 2955, 2956,
	  	                  2851, 2852, 2853, 2854, 2855};

	  int J2_allIndex[25]={9214, 9215, 9216, 9217, 9218,
	  	                  9113, 9114, 9115, 9116, 9117,
	  	                  9012, 9013, 9014, 9015, 9016,
	  	                  8911, 8912, 8913, 8914, 8915,
	  	                  8810, 8811, 8812, 8813, 8814};

	  int J3_allIndex[25]={10650, 10651, 10652, 10653, 10654, 10655, 10656, 10657, 10658, 
	  	                   10532, 10533, 10534, 10535, 10536, 10537, 10538, 10539, 
	  	                   10413, 10414, 10415, 10416, 10417, 10418, 10419, 10420};



	 /* int J1_allIndex[25]={2947, 2948, 2949, 2050, 2951,
	  	                  2846, 2847, 2848, 2849, 2850,
	  	                  2745, 2746, 2747, 2748, 2749,
	  	                  2644, 2645, 2646, 2647, 2648,
	  	                  2543, 2544, 2545, 2546, 2547};

	  int J2_allIndex[25]={8906, 8907, 8908, 8909, 8910,
	  	                  8805, 8806, 8807, 8808, 8809,
	  	                  8704, 8705, 8706, 8707, 8708,
	  	                  8603, 8604, 8605, 8606, 8607,
	  	                  8502, 8503, 8504, 8505, 8506};

	  int J3_allIndex[25]={10343, 10344, 10345, 10346, 10347, 10348, 10349, 10350, 
	  	                   10224, 10225, 10226, 10227, 10228, 10229, 10230, 10231, 
	  	                   10105, 10106, 10107, 10108, 10109, 10110, 10111, 10112};*/

	  int junctionCellSize=25;
     // printf("J1_allIndex[0]=%d\n", J1_allIndex[0]);

		while (tid< num_Junc)
		{
			Qe[tid]=0.00;
			//printf("tid= %d, hj= %f, hs=%f, zb_j= %f, zbSurf= %f, juncMaxDept= %f, junc2surfGrid= %d\n", tid, hj[tid], hs[index[tid]], zbJunc[tid], zbSurf[index[tid]], juncMaxDept[tid], index[tid]);
			// if (zbJunc[tid]+juncMaxDept[tid]!=zbSurf[index[tid]])
			// {
			// 	printf("Wrong elevation relationship: the %d junction invert elevation and its corresponding surface elevation.\n", tid);
			// 	printf("zbJunc[%d]=%f, juncMaxDep[%d]=%f, index[%d]=%d, zbSurf[%d]=%f\n", tid, zbJunc[tid], tid, juncMaxDept[tid], tid, index[tid], index[tid], zbSurf[index[tid]]);
			// 	break;
			// }


            tempAllQ[tid]=0.0;
			if (tid==0)
			{
				//printf("junction index: %d\n", tid);
				for (int k=0; k<junctionCellSize; k++)
				{
				    tempSingleQ =  2.0 / 3.0*Ci*pi*0.06*pow(2 * g, 0.5)*pow(hs[J1_allIndex[k]], 3.0 / 2.0); // positive: from surface to drainage
                    //QeLimitor(tempSingleQ,  dx,  hs[J1_allIndex[k]],  juncArea[tid],  hj[tid], dT, juncMaxDept[tid]);
                    if (tempSingleQ>0.0 && hj[tid] < juncMaxDept[tid]) 
					{
						if (tempSingleQ>hs[J1_allIndex[k]]*dx*dx/dT)
						{
							tempSingleQ=hs[J1_allIndex[k]]*dx*dx/dT;
						}
						//printf("bbbbb\n");
					}
                    //printf("tid=%d, index=%d, hs=%f, tempSingleQ=%f\n", tid, J1_allIndex[k], hs[J1_allIndex[k]], tempSingleQ );
                    tempAllQ[tid]=tempAllQ[tid]+tempSingleQ;
                    QSurf[J1_allIndex[k]]=-tempSingleQ/pow(dx,2);
				}
				Qe[tid]=tempAllQ[tid];
				//printf("tid=%d, tempAllQ[0]=%f\n", tid, tempAllQ[tid]);
			}


			else if (tid==1)
			{
				//printf("junction index: %d\n", tid);
				for (int k=0; k<junctionCellSize; k++)
				{
				    tempSingleQ = 2.0 / 3.0*Ci*pi*0.06*pow(2 * g, 0.5)*pow(hs[J2_allIndex[k]], 3.0 / 2.0); // positive: from surface to drainage
				    if (tempSingleQ>0.0 && hj[tid] < juncMaxDept[tid]) 
					{
						if (tempSingleQ>hs[J2_allIndex[k]]*dx*dx/dT)
						{
							tempSingleQ=hs[J2_allIndex[k]]*dx*dx/dT;
						}
						//printf("ccccc\n");
					}
                    //printf("tid=%d, index=%d, hs=%f, tempSingleQ=%f\n", tid, J2_allIndex[k], hs[J2_allIndex[k]], tempSingleQ );
                    tempAllQ[tid]=tempAllQ[tid]+tempSingleQ;
                    QSurf[J2_allIndex[k]]=-tempSingleQ/pow(dx,2);
				}
				Qe[tid]=tempAllQ[tid];
				//printf("tid=%d, tempAllQ[tid]=%f\n", tid, tempAllQ[tid]);

			}


			else
			{
				//printf("junction index: %d\n", tid);
				for (int k=0; k<junctionCellSize; k++)
				{
				    tempSingleQ = 2.0 / 3.0*Ci*pi*0.06*pow(2 * g, 0.5)*pow(hs[J3_allIndex[k]], 3.0 / 2.0); // positive: from surface to drainage
				    if (tempSingleQ>0.0 && hj[tid] < juncMaxDept[tid]) 
					{
						if (tempSingleQ>hs[J3_allIndex[k]]*dx*dx/dT)
						{
							tempSingleQ=hs[J3_allIndex[k]]*dx*dx/dT;
						}
						//printf("ddddd\n");
					}
                    //printf("tid=%d, index=%d, hs=%f, tempSingleQ=%f\n", tid, J3_allIndex[k], hs[J3_allIndex[k]], tempSingleQ );
                    tempAllQ[tid]=tempAllQ[tid]+tempSingleQ;
                    QSurf[J3_allIndex[k]]=-tempSingleQ/pow(dx,2);
				}
				Qe[tid]=tempAllQ[tid];
				//printf("tid=%d, tempAllQ[2]=%f\n", tid, tempAllQ[tid]);

			}





		/*	if (hs[index[tid]]<kesi)
			{
				if (hj[tid]<=juncMaxDept[tid])
				{
					Qe[tid]=-0;
					
				}
				if (hj[tid]>juncMaxDept[tid])// 2.2
				{
					Qe[tid]=-Ci*juncArea[tid] * pow(2 * g, 0.5)*pow((hj[tid] - (hs[index[tid]] + juncMaxDept[tid])), 0.5);
					
				}
			
			}
			else // surface water depth >0
			{
				
				if (hj[tid]<=juncMaxDept[tid])  // 2.1.1
				{
					
					Qe[tid] = 2.0 / 3.0*Ci*pi*juncR[tid]*2.0*pow(2 * g, 0.5)*pow(hs[index[tid]], 3.0 / 2.0); // positive: from surface to drainage
					
					
				}
				else if (hj[tid]>juncMaxDept[tid] && hj[tid]<hs[index[tid]]+juncMaxDept[tid])
				{
					
					if (hs[index[tid]]<=juncArea[tid]/(pi*2*juncR[tid]))// 2.1.2-(1)
					{
						
						Qe[tid] = Ci*pi*juncR[tid] * 2 * pow(2 * g, 0.5)*hs[index[tid]] * pow((hs[index[tid]] + juncMaxDept[tid] - hj[tid]), 0.5); // positive: from surface to drainage
					    
					}
					else //2.1.2-(2)
					{
						
						Qe[tid] = Ci*juncArea[tid] * pow(2 * g, 0.5)*pow((hs[index[tid]] + juncMaxDept[tid] - hj[tid]), 0.5); // positive: from surface to drainage
					    
					}

				}
				else if (hj[tid]==hs[index[tid]]+juncMaxDept[tid])
				{
					Qe[tid]=0;
				    
				}
				else if (hj[tid] > hs[index[tid]]+juncMaxDept[tid]) // hj>hs[index[tid]] //2.2
				{
					Qe[tid]=-Ci*juncArea[tid] * pow(2 * g, 0.5)*pow((hj[tid] - (hs[index[tid]] + juncMaxDept[tid])), 0.5);
				    
				}
				else
				{
				}
			}
			*/
				
		
			// // Qdrainage[tid] = Qe[tid];
			// // QSurf[index[tid]] = -Qe[tid];
			// // printf("tid=%d, Qdrainage= %f\n", tid, Qe[tid]);
			// // printf("tid=%d, index= %d, Qsurf= %f\n", tid, index[tid], QSurf[index[tid]]);
			
			//QeLimitor(Qe[tid],  dx,  hs[index[tid]],  juncArea[tid],  hj[tid], dT, juncMaxDept[tid]);
			Qdrainage[tid] = Qe[tid]/juncArea[tid];
			//QSurf[index[tid]] = -Qe[tid]/pow(dx,2);
			// if (tid<10)
			// {
			// 	printf("Qe[%d]=%f\n", tid, Qe[tid]);
			// 	printf("juncArea[%d]=%f, Qdrainage[%d]= %f (m/s)\n", tid, juncArea[tid], tid, Qdrainage[tid]);
		 //        //printf("tid=%d, index= %d, Qsurf[index[%d]]= %f (m/s)\n", tid, index[tid], index[tid], QSurf[index[tid]]);

			// }
          
			tid += blockDim.x * gridDim.x;
		}

		//delete[] Qe;
	}


	//----------------------------------------------------------------------------------------------------------------------
    __device__ void QeLimitor(double &QeD, double dx, double hs, double juncArea, double hj, double dT, double maxDept){

		double Qmax1, Qmax2, Qmax3, Qmax4;

		// water from surface to drainage system condition 1:
		if (QeD>0 && hj < maxDept) 
		{
			if (QeD>hs*dx*dx/dT)
			{
				QeD=hs*dx*dx/dT;
			}
			//printf("aaaa\n");
		}

		// water from surface to drainage system condition 2:
		if (QeD>0 && hj >= maxDept && hj < hs+maxDept)
		{
			Qmax1 = (hs+maxDept-hj)/(dT*(1/juncArea+1/pow(dx,2)));
			if (QeD> Qmax1)
			{
				QeD = Qmax1;
			}
		}

		// water from drainage system to water surface
		if(QeD<0 && hj>hs+maxDept)
		{
			Qmax2 = abs(hs+maxDept-hj)/(dT*(1/juncArea+1/pow(dx,2)));
			Qmax3 = hj*juncArea/dT;
			if (Qmax2 < Qmax3)
				Qmax4=Qmax2;
			else
				Qmax4=Qmax3;
			if (abs(QeD)>Qmax4)
				QeD=-Qmax4;
			
		}
		

		//printf("QeLimitor: hs= %f, dx=%f, juncArea=%f, hj=%f, Qe=%f\n",hs, dx, juncArea, hj, QeD);
	}

    //----------------------------------------------------------------------------------------------------------------------

	void surfH_limitator(cuFvMappedField<Scalar, on_cell>& h, int blocksPerGrid, int threadsPerBlock) {

		surfH_limilatorKernal << <blocksPerGrid, threadsPerBlock >> > (h.data.size(), h.data.dev_ptr());


	}

	//----------------------------------------------------------------------------------------------------------------------

	__global__ void surfH_limilatorKernal(unsigned int num_Girds, double *h_surf_device) {

		unsigned int index = blockDim.x * blockIdx.x + threadIdx.x;
		while (index  < num_Girds) {
			if (h_surf_device[index] < 0)
			{
				printf("Water depth in surfacel %d cell becomes negative\n ", index);
			}
			// if (index==154 || index==156 || index==196 || index==198)
        	//printf("surface waterdepth: cell %d, water: %f\n ", index, h_surf_device[index]);
			index += blockDim.x * gridDim.x;
		}

	}




}







