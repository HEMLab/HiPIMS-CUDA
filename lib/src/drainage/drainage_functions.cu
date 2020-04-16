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
#include "Vector.h"
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
//
namespace GC{

	
//--------global functions------------------------------------------------------------------------------------------------------
__global__ void pipe_boundaryKernal (int N, double *h_p_device, double *h_j_device, double *qx_j_device, double *qy_j_device,
	                          int *P_boundType_device,
	                          int *upsPJ_device, int *downsPJ_device, int *downsPO_device, double *upNorm1_device,
							  double *upNorm2_device, double *downNorm1_device, double *downNorm2_device,
	                          double *p_ups_h_device, double *p_ups_Q_device, double *p_downs_h_device,
							  double *p_downs_Q_device);

__global__ void pipe_calculationKernal (int CN, int num_Pipe, double T, double dt, double a_speed, 
	                              int *upsBCell_device, int*downsBCell_device, int *upsPJ_device, int *downsPJ_device,
								  int *downsPO_device, int *P_outfBound, 
								  double *P_diameter_device, double *dx_cell_device, double *zb_j_device, double *pipe_MN_device, 
								  double *zb_p_device, double *h_p_device, double *Q_p_device, 
								  double *sita_p_device, double *A_p_device, double *An_p_device, double *Qn_p_device,
								  double *pipe_flux_ups1, double *pipe_flux_ups2, double *pipe_flux_downs1,
								  double *pipe_flux_downs2, double *p_ups_h_device, double *p_ups_Q_device, 
								  double *p_downs_h_device, double *p_downs_Q_device,
								  double *downNorm1_device, double *downNorm2_device,
								  int *interval_ptr_device, double *T_series, double *h_series, double *Q_series,
								  int *T_start, int *h_start, int *Q_start, int *intervalLen, int* outfSurfIndex,
								  double *surfWaterDepth, Vector *surfQ,
								  Const_variables constVariables);

__global__ void pipe_updateKernal(int CN, int num_Pipe, double a_speed, double dT, double *P_diameter_device, double *A_p_device, 
	                        double *An_p_device, double *Q_p_device,
	                        double *Qn_p_device, double *h_p_device, double *hn_p_device, double *sita_p_device, 
							double *sitan_p_device, int *downsBCell_device, double *pipe_Ap_device, 
							double *dx_cell_device, double *t_pipe_device,
							double *pipe_flux_downs1, Const_variables constVariables);

__global__ void junction_calculationKernal(int N, int width_upsJP, int width_downsJP, double dt, double *junc_MN_device, 
	                                 double *pipe_Ap_device, 
	                                 double *P_diameter_device, double *junc_Area_device,
	                                 double *upNorm1_device, double *upNorm2_device, double *downNorm1_device, 
									 double *downNorm2_device, int *upsJP_device, int *downsJP_device, 
									 double *pipe_flux_ups1, double *pipe_flux_ups2, double *pipe_flux_downs1,
									 double *pipe_flux_downs2, double *h_j_device,double *hn_j_device, 
									 double *qx_j_device, double *qxn_j_device, double *qy_j_device, 
									 double *qyn_j_device, int *j_surfGird, double* Qexchange_device,
									 double *Ax_junc, double *Ay_junc, 
									 double *hydrostaticX, double *hydrostaticY, 
									 double *pipe_F1, double *pipe_F2,
									 double *pipe_F3,
									 Const_variables constVariables);

__global__ void junction_updateKernal(int num_Junc, double dT,  double *h_j_device, double *hn_j_device, double *qx_j_device, 
	                            double *qxn_j_device, double *qy_j_device, double *qyn_j_device, 
								double *t_junc_device, double *junc_Area_device, Const_variables constVariables);



__global__ void surfDrainQKernal_calculation(Input_information inputKeyPtr, 
											 Const_variables constVariables,
											 double *Qdrainage,
	                                         double *juncMaxDept, 
											 double *junc_R, 
											 double *junc_area, 
											 double *hj,
	                                         double *zb_junction, 
											 int *index,
											 double *zb_surface, 
											 double *hs, 
											 double *Qsuface);


//------list of device functions-------------------------------------------------------------------------------------------

__device__ int findPipeOrder(int pipe_num, int *downstreamCell, int temp1);

__device__ void TJunctionBoundary(double a_speed, int pipe_index, int bound_pipe_index, double d1, double d2, double *h_p_device, double *A_p_device, 
	                              double *sita_p_device, double *P_diameter_device, double &zb_R, double &h_R, double &A_R, double &Q_R, 
								  double &V_R, double &c_R, double &I_R, double &Phi_R, Const_variables constVariables);

__device__ void flowVariables_calculation(double a_speed, double P_diameter, double sita, double A, double Q, double h, double &V,
	                                      double &c, double &I, double &Phi, Const_variables constVariables);

__device__ int boundaryCellJugement(int pipe_num, int *BoundaryCell, int temp1);

__device__ void pipeToJuncBoundary(int j, double a_speed, double *h_junc, double *Q_junc, double h_elev, double P_diameter, double &zb_b, double &h_b,
	                          double &Q_b, double &A_b, double &Phi_b, double &V_b, double &c_b, double &I_b, Const_variables constVariables);

__device__ void pipeToOutfBoundary(double T, int *outfSurfIndex,  double *surfWaterDepth,  Vector *surfQ,  
	                               double *downNorm1_device, double *downNorm2_device,  int pipe_index, int downsPO_device, int boundType, 
								   int *interval_ptr_device, double *T_series_device, double *h_series_device, double *Q_series_device, 
								   int *T_start_device, int *h_start_device, int *Q_start_device, int *intervalLen_device, double a_speed,
								   double P_diameter, double h_p, double Q_p, double A_p, double Phi_p, double V_p, double c_p, double I_p, 
								   double zb_p, double &h_b, double &Q_b, double &A_b, double &Phi_b, double &V_b, double &c_b, double &I_b, 
								   double &zb_b, Const_variables constVariables);

__device__ void SimpleDifference(double *timeSeries, double *goalSeries, double T, int timeStart, int goalStart, 
	                             int seriesLength, int &interval_ptr_device,  double &goal_value);

__device__ void HLLFluxCalculation(int tid, double P_diameter,double a_speed, double Phi_L,double Phi_R,double V_L, double V_R,
	                                double h_L, double h_R, double c_L, double c_R, double A_L, double A_R, 
									double Q_L, double Q_R, double I_L, double I_R, double *F, Const_variables constVariables);
    
__device__ void pipeFrictionUpdate(int index, double dt, double P_manning, double P_diameter, double F_n, double s2, 
	                                double A_p_device, double Q_p_device, double h_p_device, double &Qn_p_device,
									Const_variables constVariables);

__device__ void juncFrictionUpdate(int index, int junc_num, double dt, double J_manning, double Ax_junc, double Ay_junc, double qx_j_device, 
	                                double &qxn_j_device, double qy_j_device, double &qyn_j_device, double h_j_device,
									Const_variables constVariables);

//----------------------------------------------------------------------------------------------------------------------
void gaugeOutputData(double T, Drain_gauges &drainGauges, Pipe_variables pipeVariables, Junc_variables juncVariables){
	if (drainGauges.pGaugeNum>0)
	{
		pipeGaugeOutput(T, drainGauges, pipeVariables);
	}
	if (drainGauges.jGaugeNum>0)
	{
		juncGaugeOutput(T, drainGauges, juncVariables);
	}
}
//----------------------------------------------------------------------------------------------------------------------
void pipeGaugeOutput(double T, Drain_gauges &drainGauges, Pipe_variables pipeVariables){

	for (int i=0;i<drainGauges.pGaugeNum;i++)
	{
		drainGauges.pGaugeh[i]=pipeVariables.h_p_out[drainGauges.pGaugePosi[i]];
		drainGauges.pGaugeQ[i]=pipeVariables.Q_p_out[drainGauges.pGaugePosi[i]];
		//std::cout<<"pGaugeh["<<i<<"]= "<< drainGauges.pGaugeh[i]<<std::endl;
		//std::cout<<"pGaugeQ["<<i<<"]= "<< drainGauges.pGaugeQ[i]<<std::endl;
	}

	FILE * output;
	const char* directory = "output/";
	std::string name = std::string(directory) + "Pipe_gauges"+ ".txt";
	output = fopen(name.c_str(), "a+");
	fprintf(output, "%f ", T);
	for (int i=0;i<drainGauges.pGaugeNum;i++)
	{
		fprintf(output, "%f %f  ", drainGauges.pGaugeh[i],drainGauges.pGaugeQ[i]);
	}
	fprintf(output, "\n");	
	fclose(output);
}
//----------------------------------------------------------------------------------------------------------------------
void juncGaugeOutput(double T, Drain_gauges &drainGauges, Junc_variables juncVariables){

	for (int i=0;i<drainGauges.jGaugeNum;i++)
	{
		drainGauges.jGaugeh[i]=juncVariables.h_j_out[drainGauges.jGaugePosi[i]];	
		std::cout<<"jGaugeh["<<i<<"]= "<< drainGauges.jGaugeh[i]<<std::endl;
	}
	FILE * output;
	const char* directory = "output/";
	std::string name = std::string(directory) + "Junc_gauges"+ ".txt";
	output = fopen(name.c_str(), "a+");
	fprintf(output, "%f ", T);
	for (int i=0;i<drainGauges.jGaugeNum;i++)
	{
		fprintf(output, "%f  ", drainGauges.jGaugeh[i]);
	}
	fprintf(output, "\n");	
	fclose(output);
}
//----------------------------------------------------------------------------------------------------------------------
void gaugeSettings(Drain_gauges &drainGauges){
	obtainGaugeNumber(drainGauges);
	drainGauges.setHostMemory(drainGauges.pGaugeNum, drainGauges.jGaugeNum);
	obtainGaugePositions(drainGauges);
}
//----------------------------------------------------------------------------------------------------------------------
void obtainGaugePositions(Drain_gauges &drainGauges){
	int index_p, index_j;
	int temp; // interval between pipe_gauge_position to junc_gauge_position
	if (drainGauges.pGaugeNum==0 || drainGauges.pGaugeNum==1)
	{
		temp=1+1;
	}
	else
	{
		temp=1+drainGauges.pGaugeNum;
	}
	std::string line;
	std::ifstream InFile;
	InFile.open("input/drainage/gauge.dat");
	if (!InFile) 
		std::cout<<"No gauge file exists"<<std::endl;
	else
	{
		int n=1;  // ptr of the row when reading the file
		int gauge_int;
		
		while (getline(InFile,line))
		{
		
			std::istringstream strm (line); 
			if (drainGauges.pGaugeNum>0)
			{
				if (n>5 && n<=5+drainGauges.pGaugeNum)
				{
					strm >> gauge_int;
					index_p=gauge_int;
					strm >> gauge_int;
					drainGauges.pGaugePosi[index_p-1]=gauge_int;
				}	
			}
			if(drainGauges.jGaugeNum>0)
			{
				if (n>4+temp+1 && n<=4+temp+1+drainGauges.jGaugeNum)
				{
					strm >> gauge_int;
					index_j=gauge_int;
					strm >> gauge_int;
					drainGauges.jGaugePosi[index_j-1]=gauge_int;
				}
			}
			
			n++;
		}
	}

	// if (drainGauges.pGaugeNum>0)
	// {
	// 	for (int i=0;i<drainGauges.pGaugeNum;i++)
	// 	{
	// 		std::cout<<"pGaugePosi["<<i<<"]= "<<drainGauges.pGaugePosi[i]<<std::endl;
	// 	}
	// }
	// if (drainGauges.jGaugeNum>0)
	// {
	// 	for (int i=0;i<drainGauges.jGaugeNum;i++)
	// 	{
	// 		std::cout<<"jGaugePosi["<<i<<"]= "<<drainGauges.jGaugePosi[i]<<std::endl;
	// 	}

	// }
	InFile.close();
}



//----------------------------------------------------------------------------------------------------------------------
void obtainGaugeNumber(Drain_gauges &drainGauges){

	std::string line;
	std::ifstream InFile;
	InFile.open("input/drainage/gauge.dat");
	if (!InFile) 
		std::cout<<"No gauge file exists"<<std::endl;
	else
	{
		int n=1;  // ptr of the row when reading the file
		int gauge_int;
		
		while (getline(InFile,line))
		{
		
			std::istringstream strm (line); 
			if (n==2)
			{
				strm >> gauge_int;
				drainGauges.pGaugeNum=gauge_int;
			}
			if (n==4)
			{
				strm >> gauge_int;
				drainGauges.jGaugeNum=gauge_int;
			}
			n++;
		}
	}
	
	std::cout<<"pipe gauge number = "<<drainGauges.pGaugeNum<<std::endl;
	std::cout<<"junction gauge number = "<<drainGauges.jGaugeNum<<std::endl;
	InFile.close();
}




//----------------------------------------------------------------------------------------------------------------------
void outfallBoundaryData(Input_information inputKeyPtr,  
	                    Outf_attributes &outfAttributes,
						Drain_outfBounds &drainOutfBounds){

	// obtain the number for each different type of outfall boundaries
	seriesLengthData(drainOutfBounds);
	
	// set host and Device memory
	drainOutfBounds.setHostMemory(drainOutfBounds.length_T, drainOutfBounds.length_h, drainOutfBounds.length_Q, inputKeyPtr.num_Outfall);
	drainOutfBounds.setDeviceMemory(drainOutfBounds.length_T, drainOutfBounds.length_h, drainOutfBounds.length_Q, inputKeyPtr.num_Outfall);
	drainOutfBounds.initIntervalPtr(inputKeyPtr.num_Outfall);

	// read the data from outfall boundary files
	outfallBoundFile(inputKeyPtr, drainOutfBounds); // read intervalLen, TStart, hStart, QStart vectors
	if (drainOutfBounds.length_T>0)
	    TSeriesBoundFile(drainOutfBounds);
	if (drainOutfBounds.length_h>0)
	    hSeriesBoundFile(drainOutfBounds);
	if (drainOutfBounds.length_Q>0)
	    QSeriesBoundFile(drainOutfBounds);
	drainOutfBounds.copyCudaMemory(drainOutfBounds.length_T, drainOutfBounds.length_h, drainOutfBounds.length_Q, inputKeyPtr.num_Outfall);

}
//----------------------------------------------------------------------------------------------------------------------
void seriesLengthData(Drain_outfBounds &drainOutfBounds){
	// obtain the T, h, Q series length that can set meomry for their host and device vector
	std::string line;
	std::ifstream InFile;
	InFile.open("input/drainage/Drainage_boundary.dat");
	if (!InFile) 
		std::cout<<"No Drainage_boundary file exists"<<std::endl;
	else
	{
		int n=1;  // ptr of the row when reading the file
		int outB_int;
		
		while (getline(InFile,line))
		{
		
			std::istringstream strm (line); 
			if (n==2)
			{
				strm >> outB_int;
				drainOutfBounds.length_T=outB_int;
			}
			if (n==4)
			{
				strm >> outB_int;
				drainOutfBounds.length_h=outB_int;
			}
			if (n==6)
			{
				strm >> outB_int;
				drainOutfBounds.length_Q=outB_int;
			}
			n++;
		}
	}
	
	/*std::cout<<"length_T = "<<drainOutfBounds.length_T<<std::endl;
	std::cout<<"length_h = "<<drainOutfBounds.length_h<<std::endl;
	std::cout<<"length_Q = "<<drainOutfBounds.length_Q<<std::endl;*/
	InFile.close();

}

//----------------------------------------------------------------------------------------------------------------------
void outfallBoundFile(Input_information inputKeyPtr, Drain_outfBounds &drainOutfBounds){

	std::string line;
	std::ifstream InFile;
	InFile.open("input/drainage/Drainage_boundary.dat");
	if (!InFile) 
		std::cout<<"No Drainage_boundary file exists"<<std::endl;
	else
	{
		int n=1;  // ptr of the row when reading the file
		int num_ptr=0;  // ptr if the row of the vector data
		int outB_int;
		int intervalTitleLine=9;  // 9 is the title line "Iinterval length of each outfall boundary*******************" 
		int TStartTitleLine=9+inputKeyPtr.num_Outfall+3;
		int hStartTitleLine=TStartTitleLine+inputKeyPtr.num_Outfall+3;
		int QStartTitleLine=hStartTitleLine+inputKeyPtr.num_Outfall+3;
		/*std::cout<<"TStartTitleLine= "<< TStartTitleLine<<std::endl;
		std::cout<<"hStartTitleLine= "<< hStartTitleLine<<std::endl;
		std::cout<<"QStartTitleLine= "<< QStartTitleLine<<std::endl;*/
		while (getline(InFile,line))
		{
		
			std::istringstream strm (line); 

			if (n>intervalTitleLine && n<=intervalTitleLine+inputKeyPtr.num_Outfall){
				strm >> outB_int;
				drainOutfBounds.intervalLen[num_ptr]=outB_int;
				num_ptr++;

				if (num_ptr==inputKeyPtr.num_Outfall)
					num_ptr=0;  //set the ptr to be 0 for reading the next start points vector
			}

			if (n>TStartTitleLine && n<=TStartTitleLine+inputKeyPtr.num_Outfall){
				strm >> outB_int;
				drainOutfBounds.T_start[num_ptr]=outB_int;
				num_ptr++;

				if (num_ptr==inputKeyPtr.num_Outfall)
					num_ptr=0;  //set the ptr to be 0 for reading the next start points vector
			}

			if (n>hStartTitleLine && n<=hStartTitleLine+inputKeyPtr.num_Outfall){
				strm >> outB_int;
				drainOutfBounds.h_start[num_ptr]=outB_int;
				num_ptr++;

				if (num_ptr==inputKeyPtr.num_Outfall)
					num_ptr=0;  //set the ptr to be 0 for reading the next start points vector
			}
			if (n>QStartTitleLine && n<=QStartTitleLine+inputKeyPtr.num_Outfall){
				strm >> outB_int;
				drainOutfBounds.Q_start[num_ptr]=outB_int;
				num_ptr++;

				if (num_ptr==inputKeyPtr.num_Outfall)
					num_ptr=0;  //set the ptr to be 0 for reading the next start points vector
			}
		
			n++;
		}
	}

	//for(int i=0;i<inputKeyPtr.num_Outfall;i++)
	//{
	//	std::cout<<"interval["<<i<<"]= "<<drainOutfBounds.intervalLen[i]<<std::endl;
	//	std::cout<<"TStart["<<i<<"]= "<<drainOutfBounds.T_start[i]<<std::endl;
	//	std::cout<<"hStart["<<i<<"]= "<<drainOutfBounds.h_start[i]<<std::endl;
	//	std::cout<<"QStart["<<i<<"]= "<<drainOutfBounds.Q_start[i]<<std::endl;
	//}
	
	InFile.close();


}

//----------------------------------------------------------------------------------------------------------------------
void TSeriesBoundFile(Drain_outfBounds &drainOutfBounds){
	int n=0;
	double outfB_double;
	std::string line;
	std::ifstream InFile;
	InFile.open("input/drainage/Drainage_T_series.dat");
	if (!InFile) 
		std::cout<<"No Drainage_T_B file exists"<<std::endl;
	else
	{
		while (getline(InFile,line))
		{
		
			std::istringstream strm (line); 
			strm >> outfB_double;
			drainOutfBounds.T_series[n]=outfB_double;
			//std::cout<<"T_series["<<n<<"]= "<<drainOutfBounds.T_series[n]<<std::endl;
			n++;
			if (n==drainOutfBounds.length_T)
				break;
		}
	}
	//std::cout<<"in T_series file, n=: "<<n<<std::endl;
	InFile.close();
}



//----------------------------------------------------------------------------------------------------------------------
void hSeriesBoundFile(Drain_outfBounds &drainOutfBounds){

	// read h boundary data for outfalls
	int n=0;
	double outfB_double;
	std::string line;
	std::ifstream InFile;
	InFile.open("input/drainage/Drainage_h_series.dat");
	if (!InFile) 
		std::cout<<"No Drainage_h_B file exists"<<std::endl;
	else
	{
		while (getline(InFile,line))
		{
		
			std::istringstream strm (line); 
			//strm >> outfB_double;
			strm >> outfB_double;
			drainOutfBounds.h_series[n]=outfB_double;
			//std::cout<<"h_series["<<n<<"]= "<<drainOutfBounds.h_series[n]<<std::endl;
			n++;
			if (n==drainOutfBounds.length_h)
				break;
		}
	}
	InFile.close();
}

//----------------------------------------------------------------------------------------------------------------------
void QSeriesBoundFile(Drain_outfBounds &drainOutfBounds){
	// read h boundary data for outfalls
	int n=0;
	double outfB_double;
	std::string line;
	std::ifstream InFile;
	InFile.open("input/drainage/Drainage_Q_series.dat");
	if (!InFile) 
		std::cout<<"No Drainage_Q_B file exists"<<std::endl;
	else
	{
		while (getline(InFile,line))
		{
		
			std::istringstream strm (line); 
			//strm >> outfB_double;
			strm >> outfB_double;
			drainOutfBounds.Q_series[n]=outfB_double;
			//std::cout<<"Q_series["<<n<<"]= "<<drainOutfBounds.Q_series[n]<<std::endl;
			n++;
			if (n==drainOutfBounds.length_Q)
				break;
		}
	}
	InFile.close();

}


//----------------------------------------------------------------------------------------------------------------------
void readInputFile(char* Keywords[], Input_information &inputKeyPtr, 
	               Junc_attributes &juncAttributes, 
				   Outf_attributes &outfAttributes,
				   Pipe_attributes &pipeAttributes,
				   Pipe_variables &pipeVariables){
	// read input file

	// read the introduction file to obtain the basic information of this project including number of junctions and pipes
	introductionFile(Keywords, inputKeyPtr);

	// creat memory for all kinds of objects
	hostMemoryCreationForObjects(inputKeyPtr, juncAttributes, outfAttributes, pipeAttributes);

    // read junction file
	junctionFile(inputKeyPtr,juncAttributes);

	// read outfall file
	outfallFile(inputKeyPtr, outfAttributes);

	// read pipe file
	pipeFile(inputKeyPtr,pipeAttributes);

	// read pipe vertice file
	verticeFile(inputKeyPtr,pipeAttributes);

	// read pipe cells file
	pipeCellsDataFile(inputKeyPtr,pipeAttributes, pipeVariables);

	// read the surface grid index of junctions and outfalls
	NodeGridFile(inputKeyPtr, juncAttributes, outfAttributes);

}



//----------------------------------------------------------------------------------------------------------------------
void pipeCellsDataFile(Input_information inputKeyPtr, 
					   Pipe_attributes pipeAttributes,
					   Pipe_variables &pipeVariables ){
	
	int n=0;
	int totalCellNum_ptr, cellNum_ptr, cellWidth_ptr;
	

	std::string line;
	std::ifstream InFile;
	InFile.open("input/drainage/PipeCells.dat");
	if (!InFile) 
		std::cout<<"No PipeCells file exists"<<std::endl;
	while (getline(InFile,line))
	{
		n++;
		std::istringstream strm (line); 
		char token[1024];
		strm >> token;
		if (token[0]=='&')
		{
			totalCellNum_ptr=n+2;
		}

		if (token[0]=='-')
		{
			cellNum_ptr=n+2;
		}

		if (token[0]=='*')
		{
			cellWidth_ptr=n+2;
		}

	}
	/*std::cout<<"totalCellNum_ptr= "<<totalCellNum_ptr<<std::endl;
	std::cout<<"cellNum_ptr= "<<cellNum_ptr<<std::endl;
	std::cout<<"cellWidth_ptr= "<<cellWidth_ptr<<std::endl;*/
	InFile.close();

	getPipeCellsData(inputKeyPtr, pipeAttributes, pipeVariables, totalCellNum_ptr, cellNum_ptr, cellWidth_ptr);
}


//----------------------------------------------------------------------------------------------------------------------
void getPipeCellsData(Input_information inputKeyPtr, 
					   Pipe_attributes pipeAttributes,
					   Pipe_variables &pipeVariables,
					   int totalCellNum_ptr, int cellNum_ptr, int cellWidth_ptr){
	int P_int1, P_int2;
	double P_double;
	int n=0;
	pipeVariables.p_cellNum = new int[inputKeyPtr.num_Pipe];
	

	std::string line;
	std::ifstream InFile;
	InFile.open("input/drainage/PipeCells.dat");
	if (!InFile) 
		std::cout<<"No PipeCells file exists"<<std::endl;

	while (getline(InFile,line))
	{
		n++;
		std::istringstream strm (line); 
		
		if (n==totalCellNum_ptr)
		{
			strm >> P_int1;
			pipeVariables.totalCellNum = P_int1;
			pipeVariables.dx_cell = new double[pipeVariables.totalCellNum];
		}

		if (n>=cellNum_ptr && n< cellNum_ptr+inputKeyPtr.num_Pipe)
		{
			strm >> P_int1;
			strm >> P_int2;
			pipeVariables.p_cellNum[P_int1-1]=P_int2;
		}

		if (n>=cellWidth_ptr && n< cellWidth_ptr+pipeVariables.totalCellNum)
		{
			strm >> P_int1;
			strm >> P_double;
			pipeVariables.dx_cell[P_int1-1]=P_double;
		}
		
	}
	
	InFile.close();
	std::cout<<"total cell number: "<<pipeVariables.totalCellNum<<std::endl;
	// for (int i=0;i<10;i++)
	// {
	// 	std::cout<<"pipeVariables.p_cellNum["<<i<<"]= "<<pipeVariables.p_cellNum[i]<<std::endl;
	// }
	// for (int i=0;i<10;i++)
	// {
	// 	std::cout<<"pipeVariables.dx_cell["<<i<<"]= "<<pipeVariables.dx_cell[i]<<std::endl;
	// }

}

//----------------------------------------------------------------------------------------------------------------------
void NodeGridFile(Input_information inputKeyPtr,
				  Junc_attributes &juncAttributes,
				  Outf_attributes &outfAttributes){

	std::string line;
	std::ifstream InFile;
	InFile.open("input/drainage/nodeGridIndex.dat");
	if (!InFile) 
		std::cout<<"No nodeGridIndex file exists"<<std::endl;
	else
	{
		int n=1;  // ptr of the row when reading the file
		int num_ptr=0;  // ptr if the row of the vector data
		int node_int;
		
		while (getline(InFile,line))
		{
		
			std::istringstream strm (line); 
			if (n>1 && n<=1+inputKeyPtr.num_Junc)
			{
				strm >> node_int;
				strm >> node_int;
				strm >> node_int;
				strm >> node_int;
				juncAttributes.J_surfGrid[num_ptr]=node_int;
				num_ptr++;
				if (num_ptr==inputKeyPtr.num_Junc)
					num_ptr=0;  //set the ptr to be 0 for reading the next start points vector
			}
			if(inputKeyPtr.num_Outfall>0)
			{
				if (n>1+inputKeyPtr.num_Junc+1 && n<=1+inputKeyPtr.num_Junc+1+inputKeyPtr.num_Outfall)
				{
					strm >> node_int;
					strm >> node_int;
					strm >> node_int;
					strm >> node_int;
					outfAttributes.O_surfGrid[num_ptr]=node_int;
					num_ptr++;
				}

			}
		
			n++;
		}
	}
	

	InFile.close();
	/*for (int i=0;i<inputKeyPtr.num_Junc;i++)
		printf("junc_surfaceGrid[%d]=%d\n", i, juncAttributes.J_surfGrid[i]);
	for (int k=0; k<inputKeyPtr.num_Outfall;k++)
		printf("outf_surfaceGrid[%d]=%d\n", k, outfAttributes.O_surfGrid[k]);*/
}


//----------------------------------------------------------------------------------------------------------------------
void introductionFile(char* Keywords[], Input_information &inputKeyPtr){
    // the first time to open the file is to find the pointer position of each necessary data	
	firstOpenIntroductionFile(Keywords, inputKeyPtr);
	secondOpenIntroductionFile(inputKeyPtr);
}


//----------------------------------------------------------------------------------------------------------------------
void firstOpenIntroductionFile(char* Keywords[], Input_information &inputKeyPtr){
    int n=0;  // n means the row number of the key words
	int newsect=-1;

	std::string line;
	std::ifstream InFile;

	InFile.open("input/drainage/Introduction.dat");
	if (!InFile) 
		std::cout<<"no input file exist"<<std::endl;


	while (getline(InFile,line))
	{
		n++;
		std::istringstream strm (line); 
		char token[1024];
		strm >> token;
		if (token[0]=='[')
		{
			newsect = findMatch(token, Keywords);

			switch (newsect)
			{
				case 0:
					inputKeyPtr.p_aspeed=n;
					//std::cout<<"row of a_speed= "<<inputKeyPtr.p_aspeed<<std::endl;
					break;
			

				case 1:
					inputKeyPtr.p_juncNum=n;
					//std::cout<<"row of junction number= "<<inputKeyPtr.p_juncNum<<std::endl;
					break;

				case 2:
					inputKeyPtr.p_outfallNum=n;
					//std::cout<<"row of outfall number= "<<inputKeyPtr.p_outfallNum<<std::endl;
					break;

				case 3:
					inputKeyPtr.p_pipeNum=n;
					//std::cout<<"row of pipe data= "<<inputKeyPtr.p_pipeNum<<std::endl;
					break;

				//case 4:
				//	inputKeyPtr.p_outfalls=n;
				//	//std::cout<<"row of outfall data= "<<p_outfalls<<std::endl;
				//	break;

				//case 5:
				//	inputKeyPtr.p_pipes=n;
				//	//std::cout<<"row of pipe data= "<<p_pipes<<std::endl;
				//	break;

				//case 6:
				//	inputKeyPtr.p_xsections=n;
				//	//std::cout<<"row of xsection data= "<<p_xsections<<std::endl;
				//	break;

				//case 7:
				//	inputKeyPtr.p_vertices=n;
				//	//std::cout<<"row of pipe vertice data= "<<p_vertices<<std::endl;
				//	break;

				//case 8:
				//	inputKeyPtr.p_symbols=n;
				//	//std::cout<<"row of symbol data= "<<p_symbols<<std::endl;
				//	break;
			}
		}
	}

	InFile.close();

}



//----------------------------------------------------------------------------------------------------------------------
int  findMatch(char *s, char *keyword[])
//
//  Input:   s = character string
//           keyword = array of keyword strings
//  Output:  returns index of matching keyword or -1 if no match found  
//  Purpose: finds match between string and array of keyword strings.
//
{
   int i = 0;
   while (keyword[i] != NULL)
   {
      if (match(s, keyword[i])) return(i);
      i++;
   }
   return(-1);
}



//----------------------------------------------------------------------------------------------------------------------
int  match(char *str, char *substr)
//
//  Input:   str = character string being searched
//           substr = sub-string being searched for
//  Output:  returns 1 if sub-string found, 0 if not
//  Purpose: sees if a sub-string of characters appears in a string
//           (not case sensitive).
//
{
    int i,j;

    // --- fail if substring is empty
    if (!substr[0]) return(0);

    // --- skip leading blanks of str
    for (i = 0; str[i]; i++)
    {
        if (str[i] != ' ') break;
    }

    // --- check if substr matches remainder of str
    for (i = i,j = 0; substr[j]; i++,j++)
    {
        if (!str[i] || UCHAR(str[i]) != UCHAR(substr[j])) return(0);
    }
    return(1);
}

//----------------------------------------------------------------------------------------------------------------------
char UCHAR(char x){
	
//#define UCHAR(x) (((x) >= 'a' && (x) <= 'z') ? ((x)&~32) : (x))
	if ((x) >= 'a' && (x) <= 'z')
		return (x)&~32;
	else
		return (x);
}


//----------------------------------------------------------------------------------------------------------------------
void secondOpenIntroductionFile(Input_information &inputKeyPtr){
	// the second time to open the file is to read the necessary data information including juntion and pipe number
	// and creat memory for them 

	int n=0;
	std::string line;
	std::ifstream InFile;

	InFile.open("input/drainage/Introduction.dat");
	while (getline(InFile,line))
	{
		n++;
		std::istringstream strm (line); 
		if (n==(inputKeyPtr.p_aspeed+1))
		{
			strm >> inputKeyPtr.a_speed;
			std::cout<<"a_speed= "<<inputKeyPtr.a_speed<<std::endl;
		}
		else if (n==(inputKeyPtr.p_juncNum)+1)
		{
			strm >> inputKeyPtr.num_Junc;
			std::cout<<"junction_number= "<<inputKeyPtr.num_Junc<<std::endl;
		
		}
		else if (n==(inputKeyPtr.p_outfallNum+1))
		{
			strm >> inputKeyPtr.num_Outfall;
			std::cout<<"outfall_number= "<<inputKeyPtr.num_Outfall<<std::endl;
		}

		else if (n==(inputKeyPtr.p_pipeNum+1))
		{
			strm >> inputKeyPtr.num_Pipe;
			std::cout<<"pipe_number= "<<inputKeyPtr.num_Pipe<<std::endl;
		}
	
		else
		{
			//std::cout<<"nothing"<<std::endl;
		}


	}
	InFile.close();

}

//----------------------------------------------------------------------------------------------------------------------
void hostMemoryCreationForObjects(Input_information &inputKeyPtr, 
		                         Junc_attributes &juncAttributes,
								 Outf_attributes &outfAttributes,
				                 Pipe_attributes &pipeAttributes){
	// creat memory for all objects;
	juncAttributes.setHostMemory(inputKeyPtr.num_Junc);
	pipeAttributes.setHostMemory(inputKeyPtr.num_Pipe);
	outfAttributes.setHostMemory(inputKeyPtr.num_Outfall);

}

//----------------------------------------------------------------------------------------------------------------------
void junctionFile(Input_information &inputKeyPtr, 
		          Junc_attributes &juncAttributes){
	// read junction attributes from junction file
	int n=0;
	int i_j=0;
	int j_int;
	double j_double;
	std::string line;
	std::ifstream InFile;
	InFile.open("input/drainage/Junction.dat");
	if (!InFile) 
		std::cout<<"No junction data file exists"<<std::endl;
	while (getline(InFile,line))
	{
		n++;
		std::istringstream strm (line); 
		if (n>4 && n<=4+inputKeyPtr.num_Junc)
		{
			strm >> j_int;
			juncAttributes.J_index[i_j]=j_int;
			//std::cout<< "junction index["<<i_j<<"]= "<<juncAttributes.J_index[i_j]<< std::endl;

			strm >> j_double;
			juncAttributes.J_elev[i_j]=j_double;
			//std::cout<< "junction elevation["<<i_j<<"]= "<<juncAttributes.J_elev[i_j]<<std::endl;

			strm >> j_double;
			juncAttributes.J_maxDep[i_j]=j_double;
			//std::cout<< "junction maxim depth["<<i_j<<"]= "<<juncAttributes.J_maxDep[i_j]<<std::endl;

			strm >> j_double;
			juncAttributes.J_radius[i_j]=j_double;
			//std::cout<< "junction radius["<<i_j<<"]= "<<juncAttributes.J_radius[i_j]<<std::endl;

			strm >> j_double;
			juncAttributes.J_xcoor[i_j]=j_double;
			//std::cout<< "junction xcoor["<<i_j<<"]= "<<juncAttributes.J_xcoor[i_j]<<std::endl;

			strm >> j_double;
			juncAttributes.J_ycoor[i_j]=j_double;
			//std::cout<< "junction ycoor["<<i_j<<"]= "<<juncAttributes.J_ycoor[i_j]<<std::endl;

			strm >> j_double;
			juncAttributes.J_manning[i_j]=j_double;
			//std::cout<< "junction manning["<<i_j<<"]= "<<juncAttributes.J_manning[i_j]<<std::endl;

			strm >> j_double;
			juncAttributes.J_initWaterDepth[i_j]=j_double;
			//std::cout<< "junction initWaterDepth["<<i_j<<"]= "<<juncAttributes.J_initWaterDepth[i_j]<<std::endl;
			
			strm >> j_double;
			juncAttributes.J_initXFlowRate[i_j]=j_double;
			//std::cout<< "junction initXFlowRate["<<i_j<<"]= "<<juncAttributes.J_initXFlowRate[i_j]<<std::endl;

			strm >> j_double;
			juncAttributes.J_initYFlowRate[i_j]=j_double;
			//td::cout<< "junction initYFlowRate["<<i_j<<"]= "<<juncAttributes.J_initYFlowRate[i_j]<<std::endl;

			i_j++;
		}
	}
	InFile.close();
}


//----------------------------------------------------------------------------------------------------------------------
void outfallFile(Input_information &inputKeyPtr,
	             Outf_attributes &outfAttributes)
{
	if(inputKeyPtr.num_Outfall>0)
	{

		// read outfall attributes from outfall file
		int n=0;
		int i_O=0;
		int O_int;
		double O_double;
		std::string line;
		std::ifstream InFile;
		InFile.open("input/drainage/Outfall.dat");
		if (!InFile) 
			std::cout<<"No pipe data file exists"<<std::endl;
		while (getline(InFile,line))
		{
			n++;
			std::istringstream strm (line); 
			if (n>4 && n<=4+inputKeyPtr.num_Outfall)
			{
				strm >> O_int;
				outfAttributes.O_index[i_O]=O_int;
				//std::cout<< "outfall index["<<i_O<<"]= "<<outfAttributes.O_index[i_O]<< std::endl;

				strm >> O_double;
				outfAttributes.O_elev[i_O]=O_double;
				//std::cout<< "outfall elevation["<<i_O<<"]= "<<outfAttributes.O_elev[i_O]<< std::endl;

				strm >> O_double;
				outfAttributes.O_maxDep[i_O]=O_double;
				//std::cout<< "outfall O_maxDep["<<i_O<<"]= "<<outfAttributes.O_maxDep[i_O]<< std::endl;

				strm >> O_double;
				outfAttributes.O_xcoor[i_O]=O_double;
				//std::cout<< "outfall O_xcoor["<<i_O<<"]= "<<outfAttributes.O_xcoor[i_O]<< std::endl;

				strm >> O_double;
				outfAttributes.O_ycoor[i_O]=O_double;
				//std::cout<< "outfall O_ycoor["<<i_O<<"]= "<<outfAttributes.O_ycoor[i_O]<< std::endl;

				strm >> O_double;
				outfAttributes.O_manning[i_O]=O_double;
				//std::cout<< "outfall O_manning["<<i_O<<"]= "<<outfAttributes.O_manning[i_O]<< std::endl;

				strm >> O_int;
				outfAttributes.O_boundType[i_O]=O_int;
				//std::cout<< "outfall boundType["<<i_O<<"]= "<<outfAttributes.O_boundType[i_O]<< std::endl;
			
				i_O++;

			}
		
		}
		InFile.close();
	}
	else
	{
		std::cout<<"No outfall test\n";
	}

}
//----------------------------------------------------------------------------------------------------------------------
void pipeFile(Input_information &inputKeyPtr, 
		      Pipe_attributes &pipeAttributes){
	// read pipe attributes from pipe file
	int n=0;
	int i_p=0;
	int p_int;
	double p_double;
	std::string line;
	std::ifstream InFile;
	InFile.open("input/drainage/Pipe.dat");
	if (!InFile) 
		std::cout<<"No pipe data file exists"<<std::endl;
	while (getline(InFile,line))
	{
		n++;
		std::istringstream strm (line); 
		if (n>4 && n<=4+inputKeyPtr.num_Pipe)
		{
			strm >> p_int;
			pipeAttributes.P_index[i_p]=p_int;
			//std::cout<< "pipe index["<<i_p<<"] = "<<pipeAttributes.P_index[i_p]<< std::endl;

			strm >> p_int;
			pipeAttributes.P_inNode[i_p]=p_int;
			//std::cout<< "pipe upstream node["<<i_p<<"] = "<< pipeAttributes.P_inNode[i_p] << std::endl;

			strm >> p_int;
			pipeAttributes.P_outNode[i_p]=p_int;
			//std::cout<< "pipe downstream node["<<i_p<<"] = "<< pipeAttributes.P_outNode[i_p] << std::endl;

			strm>> p_int;
			pipeAttributes.P_outfall[i_p]=p_int;
			//std::cout<<"pipe outfall node["<<i_p<<"] = "<< pipeAttributes.P_outfall[i_p] << std::endl;

			strm>> p_int;
			pipeAttributes.P_outfBound[i_p]=p_int;
			//std::cout<<"pipe outfall type["<<i_p<<"] = "<< pipeAttributes.P_outfBound[i_p] << std::endl;

			strm >> p_double;
			pipeAttributes.P_length[i_p]=p_double;
			//std::cout<< "pipe length["<<i_p<<"] = "<< pipeAttributes.P_length[i_p] << std::endl;

			strm >> p_double;
			pipeAttributes.P_dx[i_p]=p_double;
			//std::cout<< "pipe dx["<<i_p<<"] = "<< pipeAttributes.P_dx[i_p] << std::endl;

			strm >> p_double;
			pipeAttributes.P_manning[i_p]=p_double;
			//std::cout<< "pipe manning coefficient["<<i_p<<"] = "<< pipeAttributes.P_manning[i_p] << std::endl;

			strm >> p_double;
			pipeAttributes.P_inletHeight[i_p]=p_double;
			//std::cout<< "pipe inletHeight["<<i_p<<"] = "<< pipeAttributes.P_inletHeight[i_p] << std::endl;

			strm >> p_double;
			pipeAttributes.P_outletHeight[i_p]=p_double;
			//std::cout<< "pipe outletHeight["<<i_p<<"] = "<< pipeAttributes.P_outletHeight[i_p] << std::endl;

			strm >> p_double;
			pipeAttributes.P_diameter[i_p]=p_double;
			//std::cout<< "pipe diameter["<<i_p<<"] = "<< pipeAttributes.P_diameter[i_p] << std::endl;

			strm >> p_double;
			pipeAttributes.P_initWaterDepth[i_p] = p_double;
			//std::cout<< "pipe initWaterDepth["<<i_p<<"] = "<< pipeAttributes.P_initWaterDepth[i_p] << std::endl;

			strm >> p_double;
			pipeAttributes.P_initFlowRate[i_p] = p_double;
			//std::cout<< "pipe initFlowRate["<<i_p<<"] = "<< pipeAttributes.P_initFlowRate[i_p] << std::endl;
			
			i_p++;
			
		}
	}
	InFile.close();
}



//----------------------------------------------------------------------------------------------------------------------
void verticeFile(Input_information &inputKeyPtr, 
		         Pipe_attributes &pipeAttributes){
	// read vertice data from vertice file
	int n=0;
	int i_p=0;
	//int p_int;
	double p_double;
	std::string line;
	std::ifstream InFile;
	InFile.open("input/drainage/Vertice.dat");
	if (!InFile) 
		std::cout<<"No Vertice data file exists"<<std::endl;
	while (getline(InFile,line))
	{
		n++;
		std::istringstream strm (line); 
		if (n>4 && n<=4+inputKeyPtr.num_Pipe)
		{
			strm >> p_double;
			pipeAttributes.P_coorUpx[i_p]=p_double;
			//if (i_p==0 || i_p==inputKeyPtr.num_Pipe-1)
			//  std::cout<< "P_coorUpx["<<i_p<<"]= "<< pipeAttributes.P_coorUpx[i_p] << std::endl;
				//printf("P_coorUpx[%d]= %.4lf\n",i_p, pipeAttributes.P_coorUpx[i_p] );
			
			strm >> p_double;
			pipeAttributes.P_coorUpy[i_p]=p_double;
			//if (i_p==0 || i_p==inputKeyPtr.num_Pipe-1)
			//std::cout<< "P_coorUpy["<<i_p<<"]= "<< pipeAttributes.P_coorUpy[i_p] << std::endl;
				//printf("P_coorUpx[%d]= %.4lf\n",i_p, pipeAttributes.P_coorUpy[i_p] );

			strm >> p_double;
			pipeAttributes.P_coorDownx[i_p]=p_double;
			//if (i_p==0 || i_p==inputKeyPtr.num_Pipe-1)
			//std::cout<< "P_coorCenx["<<i_p<<"]= "<< pipeAttributes.P_coorDownx[i_p] << std::endl;
				//printf("P_coorDownx[%d]= %.4lf\n",i_p, pipeAttributes.P_coorDownx[i_p] );

			strm >> p_double;
			pipeAttributes.P_coorDowny[i_p]=p_double;
			//if (i_p==0 || i_p==inputKeyPtr.num_Pipe-1)
			//	printf("P_coorDowny[%d]= %.4lf\n",i_p, pipeAttributes.P_coorDowny[i_p] );
			
			strm >> p_double;
			pipeAttributes.P_coorCenx[i_p]=p_double;
			//if (i_p==0 || i_p==inputKeyPtr.num_Pipe-1)
			//	printf("P_coorCenx[%d]= %.4lf\n",i_p, pipeAttributes.P_coorCenx[i_p] );

			strm >> p_double;
			pipeAttributes.P_coorCeny[i_p]=p_double;
			//if (i_p==0 || i_p==inputKeyPtr.num_Pipe-1)
			//	printf("P_coorCeny[%d]= %.4lf\n",i_p, pipeAttributes.P_coorCeny[i_p] );

			i_p++;
			
		}
	}
	InFile.close();
}

//----------------------------------------------------------------------------------------------------------------------
void getPipeCellNumber(Input_information inputKeyPtr, 
					   Pipe_attributes pipeAttributes,
					   Pipe_variables &pipeVariables){
	// calculate the total cell number of pipe
	
	int i;
	double temp;
	int N_integer;
	
	pipeVariables.p_cellNum = new int[inputKeyPtr.num_Pipe];
	pipeVariables.totalCellNum = 0;

	for (i=0;i<inputKeyPtr.num_Pipe;i++)
	{
		// calculate the cell number in each pipe
		temp = pipeAttributes.P_length[i] / pipeAttributes.P_dx[i];
		N_integer = int (temp);
		if (temp - N_integer > 0.5)
		{
			N_integer = N_integer + 1;
		}
		
		pipeVariables.p_cellNum[i] = N_integer;
		/*if (i<10)
		    std::cout<< "cell num of each pipe: "<<pipeVariables.p_cellNum[i]<<std::endl;*/
		pipeVariables.totalCellNum = pipeVariables.totalCellNum + pipeVariables.p_cellNum[i];
		
	}
	
	/*for (i=0;i<inputKeyPtr.num_Pipe;i++)
	{
		std::cout<<"each pipe cell num: "<<pipeVariables.p_cellNum[i]<<std::endl;

	}
	std::cout<<"pipe total cell num: "<< pipeVariables.totalCellNum <<std::endl;*/

}

//----------------------------------------------------------------------------------------------------------------------
void allocMemory(Input_information inputKeyPtr,
	              Pipe_variables &pipeVariables,
				  Junc_variables &juncVariables,
				  Outf_variables &outfVariables){



 
	// allocate pipe variables memory space on CPU
	pipeVariables.setHostMemory(inputKeyPtr.num_Pipe, pipeVariables.totalCellNum);
	// allocate pipe variables memory space on GPU
	pipeVariables.setDeviceMemory(inputKeyPtr.num_Pipe, pipeVariables.totalCellNum);

    // allocate junction variables memory space on CPU
	juncVariables.setHostMemory(inputKeyPtr.num_Junc);
	// allocate junction variables memory space on GPU
	juncVariables.setDeviceMemory(inputKeyPtr.num_Junc);

	if (inputKeyPtr.num_Outfall>0)
	{    // allocate junction variables memory space on CPU
		outfVariables.setHostMemory(inputKeyPtr.num_Outfall);
		// allocate junction variables memory space on GPU
		outfVariables.setDeviceMemory(inputKeyPtr.num_Outfall);
	}
	

	
}


//----------------------------------------------------------------------------------------------------------------------
void preProcessing(Input_information inputKeyPtr,
	               Pipe_attributes pipeAttributes,
			   	   Pipe_variables &pipeVariables,
				   Junc_attributes juncAttributes,
				   Junc_variables &juncVariables,
				   Outf_attributes outfAttributes,
				   Outf_variables &outfVariables,
				   Const_variables constVariables){
	// preprocess the data before calculation

	// determin the each cell width in each pipe
    //pipeCellWidth(inputKeyPtr, pipeAttributes, pipeVariables);

	// to obtain the logic matrix: the upstream and downstream node connected by a pipe
	pipeToJuncMatrix(inputKeyPtr, pipeAttributes, pipeVariables);
	
	// to obtain the logic matrix: the upstream and downstream pipe connected by a junction
	juncToPipeMatrix(inputKeyPtr, pipeVariables, juncVariables);

	if (inputKeyPtr.num_Outfall>0)
	// to obtain the logic matrix: the  downstream pipe connected by an outfall
		outfToPipeMatrix(inputKeyPtr, pipeVariables, outfVariables);

	// to obtain the up_nor1, up_nor2, down_nor1, down_nor1, outside norm vector among the interface between pipe and junction
	normVec(inputKeyPtr, pipeAttributes, pipeVariables, juncAttributes, juncVariables);
	
	// calculate the bed elevation of all pipe cells
	pipeBedElevation(inputKeyPtr, pipeAttributes, pipeVariables);
	
	// find the identity of all upstream and downstream boundary cells and save them in two matixes respectively
	boundaryPipeCellIndex(inputKeyPtr, pipeVariables);
	
	// Area calculation including cross sectional area of pipe cell and junction area
	areaCalculation(inputKeyPtr, constVariables, pipeAttributes, pipeVariables, juncAttributes, juncVariables);

	
	
}


//----------------------------------------------------------------------------------------------------------------------
void pipeCellWidth(Input_information inputKeyPtr,
	               Pipe_attributes pipeAttributes,
			   	   Pipe_variables &pipeVariables)
{
	double temp;
	int sub_p=0;  // represents the total cell number in pipe matrix before the pointer position
	int N_integer;  // the integer part of the result: length/dx
	double N_decimal; // the decimal part of the result: length/dx

	
	
	// in this function,the dx of each cell is determined, expecially for the cells close to the pipe downstream part 
	for (int i=0;i<inputKeyPtr.num_Pipe;i++)
	{
		
		if (i==0)
		{
			for (int j=0;j<pipeVariables.p_cellNum[i];j++)
			{


				if (pipeAttributes.P_length[i] > pipeAttributes.P_dx[i])
				{
					pipeVariables.dx_cell[j]=pipeAttributes.P_dx[i];
					// calculate the width of the last cell in the first pipe
					if (j==pipeVariables.p_cellNum[i]-1)
					{
						temp=pipeAttributes.P_length[i] / pipeAttributes.P_dx[i];
						N_integer = int (temp);
						N_decimal=temp-N_integer;
						if (N_integer==pipeVariables.p_cellNum[i])
						{
							pipeVariables.dx_cell[j]=pipeAttributes.P_dx[i]+(pipeAttributes.P_length[i]-N_integer*pipeAttributes.P_dx[i]);
						}
						else
						{
							pipeVariables.dx_cell[j]=pipeAttributes.P_length[i]-N_integer*pipeAttributes.P_dx[i];

						}

				    }

				}
				else if (pipeAttributes.P_length[i] = pipeAttributes.P_dx[i])
				{
					pipeVariables.dx_cell[j]= pipeAttributes.P_dx[i];
				}
				else
				{
					printf("Warnning: the input dx of %d pipe is longer than its length.\n", i);
					pipeVariables.dx_cell[j]=pipeAttributes.P_length[i];
				}

			}
		}
		else
		{
			sub_p=sub_p+pipeVariables.p_cellNum[i-1];
			for (int j=0;j<pipeVariables.p_cellNum[i];j++)
			{
				if(pipeAttributes.P_length[i] > pipeAttributes.P_dx[i])
				{
					pipeVariables.dx_cell[sub_p+j]=pipeAttributes.P_dx[i];
					if (j==pipeVariables.p_cellNum[i]-1)
					{
						temp=pipeAttributes.P_length[i] / pipeAttributes.P_dx[i];
						N_integer = int (temp);
						N_decimal=temp-N_integer;
						if (N_integer==pipeVariables.p_cellNum[i])
						{
							pipeVariables.dx_cell[sub_p+j]=pipeAttributes.P_dx[i]+(pipeAttributes.P_length[i]-N_integer*pipeAttributes.P_dx[i]);
						}
						else
						{
							pipeVariables.dx_cell[sub_p+j]=pipeAttributes.P_length[i]-N_integer*pipeAttributes.P_dx[i];

						}

					}

				}
				else if (pipeAttributes.P_length[i] = pipeAttributes.P_dx[i])
				{
					pipeVariables.dx_cell[sub_p+j] = pipeAttributes.P_dx[i];
				}
				else
				{
					printf("Warnning: the input dx of % pipe is longer than its length.\n");
					pipeVariables.dx_cell[sub_p+j] = pipeAttributes.P_length[i];
				}
			}
		}

		
	}
	for (int i=0; i<pipeVariables.totalCellNum; i++)
	{
		std::cout<<"pipeVariables.dx_cell["<<i<<"]= "<<pipeVariables.dx_cell[i]<<std::endl;
	}
}


//----------------------------------------------------------------------------------------------------------------------
void pipeToJuncMatrix(Input_information inputKeyPtr,
	                              Pipe_attributes pipeAttributes,
								  Pipe_variables &pipeVariables){

// to obtain the logic matrix: the upstream and downstream node connected by a pipe
	for (int k=0; k<inputKeyPtr.num_Pipe; k++)
	{
		if (pipeAttributes.P_inNode[k]<0)
			std::cout<<"the upstream part of pipe "<<k<< " connects to an outfall "<< std::endl;
		else
			pipeVariables.ups_PJ[k]=pipeAttributes.P_inNode[k]-1;
		/*if (k<12)
		    std::cout<< "pipeVariables.ups_PJ["<<k<<"]= "<< pipeVariables.ups_PJ[k] << std::endl;*/
		if (pipeAttributes.P_outNode[k]>0)
		{
			pipeVariables.downs_PJ[k]=pipeAttributes.P_outNode[k]-1;
		    pipeVariables.downs_PO[k]=pipeAttributes.P_outfall[k];  //-9999
		}	
		else
		{
			if (pipeAttributes.P_outfall[k]>0)
			{
				pipeVariables.downs_PJ[k]=pipeAttributes.P_outNode[k];  // -9999
			    pipeVariables.downs_PO[k]=pipeAttributes.P_outfall[k]-1;
			}
			else
				std::cout<<"the downstream of "<<k<<" pipe is neither node nor outfall"<<std::endl;
		}
			
	/*	if (k<12)
		{
			 std::cout<< "pipeVariables.downs_PJ["<<k<<"]= "<< pipeVariables.downs_PJ[k] << std::endl;
			 std::cout<< "pipeVariables.downs_PO["<<k<<"]= "<< pipeVariables.downs_PO[k] << std::endl;
		}
		  */ 
		
	}


	/* const char* directory = "output/";
	int i;
	FILE * output;


    std::string name = std::string(directory) + "PJmatrix" + ".txt";
    output = fopen(name.c_str(), "w");

	
	for (i=0;i<inputKeyPtr.num_Pipe;i++)
	{
		fprintf(output, "  %d   %d  %d\n", pipeVariables.ups_PJ[i], pipeVariables.downs_PJ[i],  pipeVariables.downs_PO[i]);
	}

    fclose(output);*/



}



//----------------------------------------------------------------------------------------------------------------------
void juncToPipeMatrix(Input_information inputKeyPtr,
			   					  Pipe_variables &pipeVariables,
								  Junc_variables &juncVariables){
	// to obtain the logic matrix: the upstream and downstream pipes connected by a junction
	// to obtain the ups_JP and downs_JP matrix
	// the width of the JP matrix is determined by the maximum number of connected pipes to a single junction
	
	//int width_upsJP, width_downsJP;
	int pm,pn;


	std::map<int, int> m_up, m_down;  
    int val_up = 0, val_down = 0;  
    for (int k = 0; k < inputKeyPtr.num_Pipe; k++)  
    {  
		if (pipeVariables.ups_PJ[k] < 0)
		{
			std::cout<<"the upstream part of pipe "<<k<< " connects to an outfall "<< std::endl;
		}
		else
		{    m_up[pipeVariables.ups_PJ[k]]++;  
			 if (m_up[pipeVariables.ups_PJ[k]] > m_up[val_up])// at the beginning m[value]=0  
				val_up = pipeVariables.ups_PJ[k]; 
		}
		
		if (pipeVariables.downs_PJ[k] >= 0) // the valid value of pipeVariables.downs_PJ[k] starts from 0 to num_pipe-1
		{
			
			m_down[pipeVariables.downs_PJ[k]]++;
			if (m_down[pipeVariables.downs_PJ[k]] > m_down[val_down])
				val_down = pipeVariables.downs_PJ[k];
		}
    }  
	/*std::cout << "the element is:" << val_up<< ",maximum repeating time is :" << m_up[val_up] << std::endl; 
	std::cout << "the element is:" << val_down<< ",maximum repeating time is :" << m_down[val_down] << std::endl; */

	juncVariables.width_upsJP=m_up[val_up];
	juncVariables.width_downsJP=m_down[val_down];

	// the CPU and GPU memory can only be released here due to the width_upsJP and width_downsJP matrix
	juncVariables.ups_JP= new int [inputKeyPtr.num_Junc * juncVariables.width_upsJP];
	juncVariables.downs_JP= new int [inputKeyPtr.num_Junc * juncVariables.width_downsJP];
	cudaMalloc((void **) &juncVariables.upsJP_device, juncVariables.width_upsJP * inputKeyPtr.num_Junc*sizeof(int));
	cudaMalloc((void **) &juncVariables.downsJP_device, juncVariables.width_downsJP * inputKeyPtr.num_Junc*sizeof(int));

	for(int i=0;i<inputKeyPtr.num_Junc * juncVariables.width_upsJP;i++)
	{
		juncVariables.ups_JP[i]=-1;
		
	}
	for(int i=0;i<inputKeyPtr.num_Junc * juncVariables.width_downsJP;i++)
	{
		juncVariables.downs_JP[i]=-1;
		
	}

	for(int i=0;i<inputKeyPtr.num_Junc;i++)
	{
		pm=0;
	
		for(int k=0;k<inputKeyPtr.num_Pipe;k++)
		{
			if(pipeVariables.ups_PJ[k]==i)
			{
				juncVariables.ups_JP[pm+i * juncVariables.width_upsJP]=k;
				pm++;
			}
			
		}
	}

	for(int i=0;i<inputKeyPtr.num_Junc;i++)
	{
		pn=0;
	
		for(int k=0;k<inputKeyPtr.num_Pipe;k++)
		{
			if(pipeVariables.downs_PJ[k]==i)
			{
				juncVariables.downs_JP[pn+i * juncVariables.width_downsJP]=k;
				pn++;
			}
			
		}
	}
	
	//for(int i=0;i<inputKeyPtr.num_Junc*juncVariables.width_upsJP;i++)
	//{
	//	std::cout<<"ups_JP["<<i<<"]= "<<juncVariables.ups_JP[i]<<std::endl;
	//}
	////std::cout<<std::endl;
	//for(int i=0;i<inputKeyPtr.num_Junc*juncVariables.width_downsJP;i++)
	//{
	//	//if (i>=89*5 && i<=89*5+4)
	//	    std::cout<<"downs_JP["<<i<<"]= "<<juncVariables.downs_JP[i]<<std::endl;
	//}
	//

	/*const char* directory = "output/";

	FILE * output1;
    std::string name = std::string(directory) + "Ups_JPmatrix" + ".txt";
    output1 = fopen(name.c_str(), "w");
	for (int i=0;i<inputKeyPtr.num_Junc*juncVariables.width_upsJP;i++)
	{
		fprintf(output1, "  %d  \n", juncVariables.ups_JP[i]);
	}

    fclose(output1);

		

	FILE * output2;
    std::string name1 = std::string(directory) + "Downs_JPmatrix" + ".txt";
    output2 = fopen(name1.c_str(), "w");
    for (int i=0;i<inputKeyPtr.num_Junc*juncVariables.width_downsJP;i++)
	{
		fprintf(output2, "  %d  \n", juncVariables.downs_JP[i]);
	}

    fclose(output2);*/




}

//----------------------------------------------------------------------------------------------------------------------
void outfToPipeMatrix(Input_information inputKeyPtr,
			   			Pipe_variables &pipeVariables,
						Outf_variables &outfVariables){
	// to obtain the logic matrix: the downstream pipes connected by a outfall
	// to obtain the  and downs_OP matrix
	// the width of the OP matrix is determined by the maximum number of connected pipes to a single outfall
	
	//int width_upsJP, width_downsJP;
	int pn;


	std::map<int, int>  m_down;  
    int val_down = 0;  
    for (int k = 0; k < inputKeyPtr.num_Pipe; k++)  
    {  

		if (pipeVariables.downs_PO[k] >= 0)  // the valid value of pipeVariables.downs_PO starts from 0 to num_outfall-1
		{
			
			m_down[pipeVariables.downs_PO[k]]++;
			if (m_down[pipeVariables.downs_PO[k]] > m_down[val_down])
				val_down = pipeVariables.downs_PO[k];
		}
    }  
	//std::cout << "the element is:" << val_down<< ",maximum repeating time is :" << m_down[val_down] << std::endl; 

	
	outfVariables.width_downsOP=m_down[val_down];

	// the CPU and GPU memory can only be released here due to the width_downsOP matrix
	
	outfVariables.downs_OP= new int [inputKeyPtr.num_Outfall * outfVariables.width_downsOP];
	cudaMalloc((void **) &outfVariables.downsOP_device, outfVariables.width_downsOP * inputKeyPtr.num_Outfall*sizeof(int));


	for(int i=0;i<inputKeyPtr.num_Outfall * outfVariables.width_downsOP;i++)
	{
		outfVariables.downs_OP[i]=-1;
		
	}

	for(int i=0;i<inputKeyPtr.num_Outfall;i++)
	{
		pn=0;
	
		for(int k=0;k<inputKeyPtr.num_Pipe;k++)
		{
			if(pipeVariables.downs_PO[k]==i)
			{
				outfVariables.downs_OP[pn+i * outfVariables.width_downsOP]=k;
				pn++;
			}
			
		}
	}
	
	/*for(int i=0;i<inputKeyPtr.num_Outfall*outfVariables.width_downsOP;i++)
	{
		
		  std::cout<<"downs_OP["<<i<<"]= "<<outfVariables.downs_OP[i]<<std::endl;
	}
*/
	/*const char* directory = "output/";
	FILE * output2;
    std::string name = std::string(directory) + "Downs_OPmatrix" + ".txt";
    output2 = fopen(name.c_str(), "w");
    for (int i=0;i<inputKeyPtr.num_Outfall* outfVariables.width_downsOP;i++)
	{
		fprintf(output2, "  %d  \n", outfVariables.downs_OP[i]);
	}

    fclose(output2);*/

	
}
//----------------------------------------------------------------------------------------------------------------------
void normVec(Input_information inputKeyPtr,
						Pipe_attributes pipeAttributes,
			   			Pipe_variables &pipeVariables,
						Junc_attributes juncAttributes,
						Junc_variables &juncVariables){
// to obtain the up_nor1, up_nor2, down_nor1, down_nor1, outside norm vector among the interface between pipe and junction

	//double x1,x2;
	double normal1, normal2;

	// calculate the upstream norm vector

	for(int i=0; i<inputKeyPtr.num_Pipe; i++)
	{
		/*if (abs(pipeAttributes.P_coorUpAx[i]-pipeAttributes.P_coorUpBx[i])>0)
		{
			x1=-1*(pipeAttributes.P_coorUpAy[i]-pipeAttributes.P_coorUpBy[i])
				/(pipeAttributes.P_coorUpAx[i]-pipeAttributes.P_coorUpBx[i]);
			x2=-x1;
			if(x1*(pipeAttributes.P_coorUpBx[i]-juncAttributes.J_xcoor[pipeVariables.ups_PJ[i]])
			   +1*(pipeAttributes.P_coorUpBy[i]-juncAttributes.J_ycoor[pipeVariables.ups_PJ[i]])>0)
			{
				normal1=x1;
				normal2=1;
			}
			else
			{
				normal1=x2;
				normal2=-1;
			}
		}
		else
		{
			x1=1;
			x2=-1;

			if(x1*(pipeAttributes.P_coorUpBx[i]-juncAttributes.J_xcoor[pipeVariables.ups_PJ[i]])
				+0*(pipeAttributes.P_coorUpBy[i]-juncAttributes.J_ycoor[pipeVariables.ups_PJ[i]])>0)
			{
				normal1=x1;
				normal2=0;
			}
			else
			{
				normal1=x2;
				normal2=0;
			}
		}
*/

		normal1=pipeAttributes.P_coorCenx[i]-pipeAttributes.P_coorUpx[i];
		normal2=pipeAttributes.P_coorCeny[i]-pipeAttributes.P_coorUpy[i];
	
		if (normal1==0 && normal2==0)
			std::cout<<"error: no upstream noral vector in "<<i<<" pipe.\n";
		else
		{
			pipeVariables.up_nor1[i]=normal1/(sqrt(pow(normal1,2)+pow(normal2,2)));
		    pipeVariables.up_nor2[i]=normal2/(sqrt(pow(normal1,2)+pow(normal2,2)));
		}

		/*std::cout<<"up_nor1["<<i<<"]= "<<pipeVariables.up_nor1[i]<<std::endl;
		std::cout<<"up_nor2["<<i<<"]= "<<pipeVariables.up_nor2[i]<<std::endl;*/

		
		
		//std::cout<<std::endl;
	}


	// calculate the downstream norm vector

	for(int i=0; i<inputKeyPtr.num_Pipe; i++)

	{
		//std::cout<<"J_ycoor[downs_PJ["<<i<<"]= "<<J_ycoor[downs_PJ[i]]<<std::endl;
		//std::cout<<"downs_PJ["<<i<<"]= "<<downs_PJ[i]<<std::endl;
		/*if (abs(pipeAttributes.P_coorDownCx[i]-pipeAttributes.P_coorDownDx[i])>0)
		{
			x1=-1*(pipeAttributes.P_coorDownCy[i]-pipeAttributes.P_coorDownDy[i])
				/(pipeAttributes.P_coorDownCx[i]-pipeAttributes.P_coorDownDx[i]);
			x2=-x1;
			if(x1*(pipeAttributes.P_coorDownDx[i]-juncAttributes.J_xcoor[pipeVariables.downs_PJ[i]])
			   +1*(pipeAttributes.P_coorDownDy[i]-juncAttributes.J_ycoor[pipeVariables.downs_PJ[i]])>0)
			{
				normal1=x1;
				normal2=1;
			}
			else
			{
				normal1=x2;
				normal2=-1;
			}
		}
		else
		{
			x1=1;
			x2=-1;

			if(x1*(pipeAttributes.P_coorDownDx[i]-juncAttributes.J_xcoor[pipeVariables.downs_PJ[i]])
			  +0*(pipeAttributes.P_coorDownDy[i]-juncAttributes.J_ycoor[pipeVariables.downs_PJ[i]])>0)
			{
				normal1=x1;
				normal2=0;
			}
			else
			{
				normal1=x2;
				normal2=0;
			}
		}*/



		normal1=pipeAttributes.P_coorCenx[i]-pipeAttributes.P_coorDownx[i];
		normal2=pipeAttributes.P_coorCeny[i]-pipeAttributes.P_coorDowny[i];
		if (normal1==0 && normal2==0)
			std::cout<<"error: no downstream noral vector in "<<i<<" pipe.\n";
		else
		{
			pipeVariables.down_nor1[i]=normal1/(sqrt(pow(normal1,2)+pow(normal2,2)));
		    pipeVariables.down_nor2[i]=normal2/(sqrt(pow(normal1,2)+pow(normal2,2)));
		}
		
		
		
		// std::cout<<"down_nor1["<<i<<"]= "<<pipeVariables.down_nor1[i]<<std::endl;
		// std::cout<<"down_nor2["<<i<<"]= "<<pipeVariables.down_nor2[i]<<std::endl;
		
		
		//std::cout<<std::endl;
	}


}

//----------------------------------------------------------------------------------------------------------------------
void pipeBedElevation(Input_information inputKeyPtr,
								  Pipe_attributes pipeAttributes,
			   					  Pipe_variables &pipeVariables){
	// calculate the bed elevation of all pipe cells
	int sub_p=0;  // represents the total cell number in pipe matrix before the pointer position
	int i;
	double length;
	//double pipe1_Bzb=0.235;
	//double pipe2_Bzb=0.235;
	for(i=0;i<inputKeyPtr.num_Pipe;i++)
	{
		//P_inletHeight[i]=J_elev[ups_PJ[i]]+P_inletHeight[i];
		//P_outletHeight[i]=J_elev[downs_PJ[i]]+P_outletHeight[i];
		if (i>0)
		{
			sub_p=sub_p+pipeVariables.p_cellNum[i-1];
		}
		for (int j=0;j<pipeVariables.p_cellNum[i];j++)
		{
			//pipeVariables.dx_cell[sub_p+j]=pipeAttributes.P_dx[i];
			if (j==0)
				length=0.5*pipeVariables.dx_cell[sub_p+j];
			else
				length=length+0.5*pipeVariables.dx_cell[sub_p+j-1]+0.5*pipeVariables.dx_cell[sub_p+j];

			pipeVariables.zb_p[sub_p+j]=pipeAttributes.P_inletHeight[i]
			-((pipeAttributes.P_inletHeight[i]-pipeAttributes.P_outletHeight[i])/pipeAttributes.P_length[i])*length;
			
			if (sub_p+j<50)
			{
				std::cout<<"dx_cell["<<sub_p+j<<"]= "<<pipeVariables.dx_cell[sub_p+j]<<std::endl;
				std::cout<<"zb["<<sub_p+j<<"]= "<<pipeVariables.zb_p[sub_p+j]<<std::endl;
			}
				
		}
	}

	int T_node1=110;
	int T_node2=147;
	int L_node1=174;
	int L_node2=191;
	double slopp3=1.03;
	double slopp4=0.91;
	double slopp5=5.07;
	double slopp6=1.26;

	for (int k=104; k<218; k++)
	{
		if (k<=T_node1)
			pipeVariables.zb_p[k]=0.134;
		if (k>T_node1 && k<T_node2)
		{
			pipeVariables.zb_p[k]=0.134-((k-T_node1-1)*0.1+0.05)*slopp3/100.00;
		}
		if (k==T_node2)
		{
			pipeVariables.zb_p[k]=0.096405;
		}
		if (k>T_node2 && k<L_node1)
		{
			pipeVariables.zb_p[k]=0.095787-((k-T_node2-1)*0.1+0.05)*slopp4/100.00;
		}
		if (k==L_node1)
		{
			pipeVariables.zb_p[k]=0.071672;
		}
		if (k>L_node1 && k<L_node2)
		{
			pipeVariables.zb_p[k]=0.070762-((k-L_node1-1)*0.1+0.05)*slopp5/100.00;
		}
		if (k==L_node2)
		{
			pipeVariables.zb_p[k]=-0.012893;
		}
		if (k>L_node2)
		{
			pipeVariables.zb_p[k]=-0.015428-((k-L_node2-1)*0.1+0.05)*slopp6/100.00;
		}
	}

//	std::cout<<"totalCellNum=                "<<pipeVariables.totalCellNum<<std::endl;

}

//----------------------------------------------------------------------------------------------------------------------
void boundaryPipeCellIndex(Input_information inputKeyPtr,
			   			   Pipe_variables &pipeVariables){
// find the identity of all upstream and downstream boundary cells and save them in two matixes respectively
    int i;
	//CNum[]: each element represents the cell number of each pipe
	for (i=0;i<inputKeyPtr.num_Pipe;i++)
	{
		if (i==0)
		{
			pipeVariables.upstreamCell[i]=0;
			pipeVariables.downstreamCell[i]=pipeVariables.p_cellNum[i]-1;
			
		}
		else
		{
			pipeVariables.upstreamCell[i]=pipeVariables.upstreamCell[i-1]+pipeVariables.p_cellNum[i-1];
			pipeVariables.downstreamCell[i]=pipeVariables.downstreamCell[i-1]+pipeVariables.p_cellNum[i];

		}
	}
	
	// for(i=0;i<inputKeyPtr.num_Pipe;i++)
	// {
	// 	std::cout<<"pipeVariables.upstreamCell["<<i<<"]= "<<pipeVariables.upstreamCell[i]<<std::endl;
	// 	std::cout<<"pipeVariables.downstreamCell["<<i<<"]= "<<pipeVariables.downstreamCell[i]<<std::endl;
	// }
}

//----------------------------------------------------------------------------------------------------------------------
void areaCalculation(Input_information inputKeyPtr,
	                  Const_variables constVariables,
	                  Pipe_attributes pipeAttributes,
			   		  Pipe_variables &pipeVariables,
					  Junc_attributes juncAttributes,
					  Junc_variables &juncVariables){
// Area calculation including cross sectional area of pipe cell and junction area
	
	int i;
	
	// pipe cross sectional area

	for (i=0;i<inputKeyPtr.num_Pipe;i++)
	{
		pipeVariables.pipe_Ap[i]=constVariables.pi*pow(pipeAttributes.P_diameter[i],2)/4.0;
	}
	
	// junction area
	for(i=0;i<inputKeyPtr.num_Junc;i++)
	{
		//juncVariables.junc_Area[i]=constVariables.pi*pow(juncAttributes.J_radius[i],2);
		juncVariables.junc_Area[i]=pow(juncAttributes.J_radius[i],2);
	}
}


//----------------------------------------------------------------------------------------------------------------------

void initFlow(Input_information inputKeyPtr,
	           Pipe_attributes pipeAttributes,
			   Pipe_variables &pipeVariables,
			   Junc_attributes juncAttributes,
			   Junc_variables &juncVariables){
	
	int i;
	int k,sub=0;
	// junction initialization

	for (i=0;i<inputKeyPtr.num_Junc;i++)
	{
	 
		juncVariables.h_j[i]=juncAttributes.J_initWaterDepth[i];
		juncVariables.qx_j[i]=juncAttributes.J_initXFlowRate[i];
		juncVariables.qy_j[i]=juncAttributes.J_initYFlowRate[i];
		
			//std::cout<<"juncVariables.h_j["<<i<<"]= "<<juncVariables.h_j[i]<<std::endl;
			//std::cout<<"juncVariables.qx_j["<<i<<"]= "<<juncVariables.qx_j[i]<<std::endl;
			//std::cout<<"juncVariables.qy_j["<<i<<"]= "<<juncVariables.qy_j[i]<<std::endl;
			//std::cout<<std::endl;
		
	}

	

	// pipe initialization

	for (i=0;i<inputKeyPtr.num_Pipe;i++)
	{
		if (i==0)
		{
			for(k=0;k<pipeVariables.p_cellNum[i];k++)	
			{
				pipeVariables.h_p[k]=pipeAttributes.P_initWaterDepth[i];
				pipeVariables.Q_p[k]=pipeAttributes.P_initFlowRate[i];
				pipeVariables.sita_p[k]=2*acos(1-2.0*pipeVariables.h_p[k]/pipeAttributes.P_diameter[i]);
				pipeVariables.A_p[k]=1/8.0*(pipeVariables.sita_p[k]-sin(pipeVariables.sita_p[k]))
									   *pipeAttributes.P_diameter[i]*pipeAttributes.P_diameter[i];
			}

		}
		else
		{
			sub=sub+pipeVariables.p_cellNum[i-1];
			for(k=0;k<pipeVariables.p_cellNum[i];k++)
			{
				pipeVariables.h_p[sub+k]=pipeAttributes.P_initWaterDepth[i];
				pipeVariables.Q_p[sub+k]=pipeAttributes.P_initFlowRate[i];
				pipeVariables.sita_p[sub+k]=2*acos(1-2.0*pipeVariables.h_p[sub+k]/pipeAttributes.P_diameter[i]);
				pipeVariables.A_p[sub+k]=1/8.0*(pipeVariables.sita_p[sub+k]-sin(pipeVariables.sita_p[sub+k]))
									   *pipeAttributes.P_diameter[i]*pipeAttributes.P_diameter[i];
                //printf("sub= %d, hp= %f, Q_p=%f, sita= %f, A_p=%f\n",sub,pipeVariables.h_p[sub],pipeVariables.Q_p[sub],pipeVariables.sita_p[sub],pipeVariables.A_p[sub]);

			}
		}
	}

	/*for (i=0;i<pipeVariables.totalCellNum;i++)
	{
		  printf("i= %d, hp= %f,  Q_p=%f, sita= %f, A_p=%f\n",i,pipeVariables.h_p[i], pipeVariables.Q_p[i],pipeVariables.sita_p[i],pipeVariables.A_p[i]);

	}*/

}

//----------------------------------------------------------------------------------------------------------------------
void copyMemory(Input_information inputKeyPtr,
	             Pipe_attributes pipeAttributes,
	             Pipe_variables &pipeVariables,
				 Junc_attributes juncAttributes,
	             Junc_variables &juncVariables,
				 Outf_attributes outfAttributes,
				 Outf_variables &outfVariables){

				//	 for (int i=0;i<4;i++)
				//	 std::cout<<" pipeAttributes.P_outfBound[= "<<i<<"]= "<< pipeAttributes.P_outfBound[i]<<std::endl;
    // copy pipe variables memory from host to device
	//void copyCudaMemory(int cellTotalNum, int num_Pipe, double *pipeD, double *pipeM)
	pipeVariables.copyCudaMemory(pipeVariables.totalCellNum, 
		inputKeyPtr.num_Pipe, pipeAttributes.P_diameter, pipeAttributes.P_manning, pipeAttributes.P_outfBound);

	// copy junction variables memory from host to device
	//void copyCudaMemory(int num_Junc, double *J_elev, double *J_manning)
	juncVariables.copyCudaMemory(inputKeyPtr.num_Junc, juncAttributes.J_elev, 
		juncAttributes.J_manning, juncAttributes.J_surfGrid, juncAttributes.J_radius,
		juncAttributes.J_maxDep);

	// 	void copyCudaMemory(int num_Outf, double *O_elev, double *O_manning, int *O_surfGrid, double *O_radius, double *maxDepth)
	// copy outfall variables memory from host to device
	if (inputKeyPtr.num_Outfall>0)
	{
		outfVariables.copyCudaMemory(inputKeyPtr.num_Outfall, outfAttributes.O_elev, outfAttributes.O_manning, 
		outfAttributes.O_surfGrid, outfAttributes.O_maxDep, outfAttributes.O_boundType);
	}
	
	
	
}

//---------------------------------------------------------------------------------------------------
void releaseMemory(Pipe_attributes &pipeAttributes,
					Pipe_variables &pipeVariables,
					Junc_attributes &juncAttributes,
					Junc_variables &juncVariables,
					Outf_attributes &outfAttributes,
					Outf_variables &ourfVariables,
					Drain_outfBounds &drainOutfBounds, 
					Drain_gauges &drainGauges){

	pipeAttributes.deleteHostMemory();
	pipeVariables.deleteHostMemory();
	pipeVariables.deleteDeviceMemory();
	juncAttributes.deleteHostMemory();
	juncVariables.deleteHostMemory();
	juncVariables.deleteDeviceMemory();
	outfAttributes.deleteHostMemory();
	ourfVariables.deleteDeviceMemory();
	drainOutfBounds.deleteHostMemory();
	drainOutfBounds.deleteDeviceMemory();
	drainGauges.deleteHostMemory();						
}

//---------------------------------------------------------------------------------------------------
__global__ void pipe_boundaryKernal (int N, double *h_p_device, double *h_j_device, double *qx_j_device, double *qy_j_device,
	                          int *P_boundType_device,
	                          int *upsPJ_device, int *downsPJ_device, int *downsPO_device, double *upNorm1_device,
							  double *upNorm2_device, double *downNorm1_device, double *downNorm2_device,
	                          double *p_ups_h_device, double *p_ups_Q_device, double *p_downs_h_device,
							  double *p_downs_Q_device){

						  
							 /* pipe_boundaryKernal(inputKeyPtr.num_Pipe,
								  juncVariables.h_j_device,
								  juncVariables.qx_j_device,
								  juncVariables.qy_j_device,
								  pipeVariables.upsPJ_device,
								  pipeVariables.downsPJ_device,
								  pipeVariables.upNorm1_device,
								  pipeVariables.upNorm2_device,
								  pipeVariables.downNorm1_device,
								  pipeVariables.downNorm2_device,
								  pipeVariables.p_ups_h_device,
								  pipeVariables.p_ups_Q_device,
								  pipeVariables.p_downs_h_device,
								  pipeVariables.p_downs_Q_device);*/
				
	int tid= blockDim.x * blockIdx.x + threadIdx.x; 
	
	//printf("upNorm1_device[tid]= %d  \n", upNorm1_device[0]);
	//printf("upNorm1_device[tid]= %d  \n",upNorm2_device[0]);
	/*while (tid<N)
	{
		printf("N= %d, h_j_device[%d]= %f\n",N, tid,  h_j_device[tid]);
		tid += blockDim.x * gridDim.x;
	}*/
	
	while (tid< N)
	{
		
		if (tid==0)
		{
			p_ups_h_device[0]=h_j_device[1];
			p_ups_Q_device[0]=qx_j_device[1];

			p_downs_h_device[0]=h_p_device[2];
			p_downs_Q_device[0]=0;
		}
		else if (tid==1)
		{
			p_ups_h_device[1]=h_j_device[0];
			p_ups_Q_device[1]=qx_j_device[0];

			p_downs_h_device[1]=h_p_device[3];
			p_downs_Q_device[1]=0;
		}
		else if (tid==2)
		{
			p_ups_h_device[2]=h_j_device[2];
			p_ups_Q_device[2]=-qy_j_device[2];

			p_downs_h_device[2]=0;
			p_downs_Q_device[2]=0;
		}
		else
		{
			p_ups_h_device[3]=h_j_device[2];
			p_ups_Q_device[3]=qx_j_device[2];

			p_downs_h_device[3]=0;
			p_downs_Q_device[3]=0;
		}

		//if (upsPJ_device[tid]<0)
		//{
		//	printf("error: the upstream node becomes an outfall");
		//}
		//else
		//{
		//	p_ups_h_device[tid]=h_j_device[upsPJ_device[tid]];
		//	p_ups_Q_device[tid]=qx_j_device[upsPJ_device[tid]]*(upNorm1_device[tid]*1+upNorm2_device[tid]*0)+
		//	                 qy_j_device[upsPJ_device[tid]]*(upNorm1_device[tid]*0+upNorm2_device[tid]*1) ;
		//}
  //       
		//if (downsPJ_device[tid]<0)
		//{
		//	if (downsPO_device[tid]<0)
		//	{
		//		printf("the downstram part of %d pipe is connected by neither a junction nor an outfall\n", tid);
		//	}
		//	else
		//	{
		//		
		//		//// BOUNDARY CONDITIONS FROM THE OUTFALLS
		//		//printf("P_boundType_device[%d]= %d\n", tid, P_boundType_device[tid]);
		//		//if (P_boundType_device[tid]==100 ||P_boundType_device[tid]==200) // 100 represents fixed boundary and 200 represents the continuous boundary
		//		//{
		//			p_downs_h_device[tid]=0.00;
		//			p_downs_Q_device[tid]=0.00;
		//			
		//		//}
		//		//else
		//		//{
		//		//	
		//		//	/*p_downs_h_device[tid]=h_O_device[downsPO_device[tid]];
		//		//    p_downs_Q_device[tid]=qx_O_device[downsPO_device[tid]]*(-downNorm1_device[tid]*1+(-downNorm2_device[tid])*0)+
		//	 //                  qy_O_device[downsPO_device[tid]]*(-downNorm1_device[tid]*0+(-downNorm2_device[tid])*1) ;*/
		//		//}

		//	}
		//}
		//else
		//{
		//	 p_downs_h_device[tid]=h_j_device[downsPJ_device[tid]];
		//	 p_downs_Q_device[tid]=qx_j_device[downsPJ_device[tid]]*(-downNorm1_device[tid]*1+(-downNorm2_device[tid])*0)+
		//	                   qy_j_device[downsPJ_device[tid]]*(-downNorm1_device[tid]*0+(-downNorm2_device[tid])*1) ;
		//}
		
	
		
		

		 tid += blockDim.x * gridDim.x;
		 
	}

}


//----------------------------------------------------------------------------------------------------------------------
__global__ void pipe_calculationKernal (int CN, int num_Pipe, double T, double dt, double a_speed, 
	                              int *upsBCell_device, int*downsBCell_device, int *upsPJ_device, int *downsPJ_device,
								  int *downsPO_device, int *P_outfBound, 
								  double *P_diameter_device, double *dx_cell_device, double *zb_j_device, double *pipe_MN_device, 
								  double *zb_p_device, double *h_p_device, double *Q_p_device, 
								  double *sita_p_device, double *A_p_device, double *An_p_device, double *Qn_p_device,
								  double *pipe_flux_ups1, double *pipe_flux_ups2, double *pipe_flux_downs1,
								  double *pipe_flux_downs2, double *p_ups_h_device, double *p_ups_Q_device, 
								  double *p_downs_h_device, double *p_downs_Q_device,
								  double *downNorm1_device, double *downNorm2_device,
								  int *interval_ptr_device, double *T_series, double *h_series, double *Q_series,
								  int *T_start, int *h_start, int *Q_start, int *intervalLen, int* outfSurfIndex,
								  double *surfWaterDepth, Vector *surfQ,
								  Const_variables constVariables){

//	double mass_pipe=0;  // to verify mass conservation
	int j;
	double zb_L, zb_R, sita_L, sita_R, A_L, A_R, Q_L, Q_R, h_L, h_R, V_L, V_R;
	double T_L, T_R, c_L, c_R, H_L, H_R, I_L, I_R, Phi_L, Phi_R;
	double zbfe, zbfw;
	double Fe[2], Fw[2];
	double s1, s2; // source term
	double dAdt, dqxdt;
	double g=9.81;
	
	int T_node_position1=110;
	int T_node_position2=147;
	int L_node_position1=174;
	int L_node_position2=191;
	double T_valveOn=200.0;
	int valve_position=106;

	/*int T_node_position1=2;
	int T_node_position2=3;
	double T_valveOn=500.0;
	int valve_position=2;*/

	// N : number of all pipes; M: cell number of a single pipe
	int tid= blockDim.x * blockIdx.x + threadIdx.x;  

	while (tid< CN)
	{

       	//double surfQx = surfQ[5].x;
	    //double surfQy = surfQ[5].y;
        //printf("surfQx[5]=%f, surfQy[5]=%f\n", surfQx, surfQy);

		//i=(tid%N);    // i denotes the colome of the 2D matrix (the order of cells)
        //j=(tid/N);    // j denotes the row of the 2D matrix (the order of pipes)

		
		j=findPipeOrder(num_Pipe,downsBCell_device,tid);

		// calculate the east face

        // left hand of east face
		zb_L=zb_p_device[tid];
		sita_L=sita_p_device[tid];
		A_L=A_p_device[tid];
		Q_L=Q_p_device[tid];
		h_L=h_p_device[tid];
		//printf("tid= %d, zb_L= %f\n",tid,zb_L);
		//printf("tid= %d, sita_L= %f\n",tid,sita_L);
		//printf("tid= %d, A_L= %f\n",tid,A_L);
		//printf("tid= %d, Q_L= %f\n",tid,Q_L);
		//printf("tid= %d, h_L= %f\n",tid,h_L);



		flowVariables_calculation(a_speed, P_diameter_device[j], sita_L, A_L, Q_L, h_L, V_L, c_L, I_L, Phi_L, constVariables);
		
		// right hand of east face

		if(boundaryCellJugement(num_Pipe,downsBCell_device,tid))
		{
			if (j==0)
			{
				
				TJunctionBoundary(a_speed, j, T_node_position1, P_diameter_device[2], P_diameter_device[j], h_p_device, A_p_device, sita_p_device, P_diameter_device, zb_R, h_R, A_R, Q_R, V_R,
								                    c_R, I_R, Phi_R, constVariables);
			}
			else if (j==1)
			{
				TJunctionBoundary(a_speed, j, T_node_position2, P_diameter_device[2], P_diameter_device[j], h_p_device, A_p_device, sita_p_device, P_diameter_device, zb_R, h_R, A_R, Q_R, V_R,
								                    c_R, I_R, Phi_R, constVariables);
			}
			else
			{
				// the downstream of the pipe is connected by a junction
				if(downsPJ_device[j]>=0)  // the valid value of downsPJ_device starts from 0 to num_Pipe-1
				{
					pipeToJuncBoundary(j, a_speed, p_downs_h_device, p_downs_Q_device, zb_j_device[downsPJ_device[j]], P_diameter_device[j], zb_R, h_R, Q_R, A_R, Phi_R,
					          V_R, c_R, I_R, constVariables);
				

				}
				else
				{
					//printf("BoundOutfall calculation: %d\n",j);
					if (P_outfBound[j]>=0 && downsPO_device[j]>=0)
					{
						pipeToOutfBoundary(T, outfSurfIndex, surfWaterDepth, surfQ, downNorm1_device, downNorm2_device, j, downsPO_device[j], 
							               P_outfBound[j], interval_ptr_device, T_series, h_series, Q_series, 
							               T_start, h_start, Q_start, intervalLen, a_speed, P_diameter_device[j], 
							               h_L, Q_L, A_L, Phi_L, V_L, c_L, I_L, zb_L, h_R, Q_R, A_R, Phi_R, V_R, c_R, I_R, zb_R,
										   constVariables);
			
					
					}
					else
					{
						printf("the downstream of pipe %d is connected by a outfall with a wrong outfall boundary type \n");
					}
					
					/*printf("tid= %d, phi_L=%f, V_L=%f, V_R=%f, h_L=%f, h_R=%f, c_L=%f, c_R=%f, A_L=%f, A_R=%f, Q_L=%f, Q_R=%f, I_L=%f,I_R=%f\n",tid ,Phi_L, V_L, V_R, h_L, h_R, c_L, c_R, A_L, A_R, Q_L, Q_R,
										I_L, I_R);*/
					

					
				}
				//printf("downstream: j= %d, I_R= %f, h= %f, h_b= %f, Q=%f Q_b=%f \n",j,I_R, p_downs_h_device[j], h_R, p_downs_Q_device[j], Q_R);
			
			}




		
		}
		else
		{
			zb_R=zb_p_device[tid+1];
			sita_R=sita_p_device[tid+1];
			A_R=A_p_device[tid+1];
			Q_R=Q_p_device[tid+1];
			h_R=h_p_device[tid+1];

			flowVariables_calculation(a_speed, P_diameter_device[j], sita_R, A_R, Q_R, h_R, V_R, c_R, I_R, Phi_R, constVariables);
		}

		if(tid==valve_position && T>T_valveOn)
		{
			//printf("aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa\n");
			zb_R=zb_p_device[valve_position];
			sita_R=sita_p_device[valve_position];
			A_R=A_p_device[valve_position];
			Q_R= - Q_p_device[valve_position];
			h_R=h_p_device[valve_position];
			flowVariables_calculation(a_speed, P_diameter_device[j], sita_R, A_R, Q_R, h_R, V_R, c_R, I_R, Phi_R, constVariables);
		}
        
		zbfe=(zb_L+zb_R)/2.0;
	
		HLLFluxCalculation(tid, P_diameter_device[j], a_speed, Phi_L, Phi_R, V_L, V_R, h_L, h_R, c_L, c_R, A_L, A_R, Q_L, Q_R,
			                I_L, I_R, Fe, constVariables);
	    
		// 	if (tid==17)
		// {
		// 	printf("***********************east face*******************:\n");
		// 	printf("diameter=%f\n",P_diameter_device[j]);
		// 	printf("zb_L=%f, zb_R=%f, zbfe=%f\n", zb_L, zb_R, zbfe);
		// 	printf("tid= %d, phi_L=%f, phi_R=%f, V_L=%f, V_R=%f, h_L=%f, h_R=%f, c_L=%f, c_R=%f, A_L=%f, A_R=%f, Q_L=%f, Q_R=%f, I_L=%f,I_R=%f, Fe0=%f, Fe1=%f\n",tid ,Phi_L, Phi_R, V_L, V_R, h_L, h_R, c_L, c_R, A_L, A_R, Q_L, Q_R,
		// 	                I_L, I_R, Fe[0], Fe[1]);
		// 	printf("\n");

		// }
		

		//printf("tid= %d, Fe1=%f\n",tid, Fe[1]);
		if (boundaryCellJugement(num_Pipe,downsBCell_device,tid))
		{
			pipe_flux_downs1[j]=Fe[0];
			pipe_flux_downs2[j]=Fe[1];
	/*		printf("pipe_flux_downs1[%d]= %f\n",tid, pipe_flux_downs1[j]);
			printf("pipe_flux_downs2[%d]= %f\n",tid, pipe_flux_downs2[j]);*/
		}

		// west face calculation

		// right hand of west face

        zb_R=zb_L;
		sita_R=sita_L;
		A_R=A_L;
		Q_R=Q_L;
		h_R=h_L; 
		V_R=V_L;
		T_R=T_L;
		c_R=c_L;
		H_R=H_L,
		I_R=I_L;
		Phi_R=Phi_L;

		if (boundaryCellJugement(num_Pipe,upsBCell_device,tid))
		{
			// boundary conditions

			pipeToJuncBoundary(j, a_speed, p_ups_h_device, p_ups_Q_device, zb_j_device[upsPJ_device[j]], P_diameter_device[j], zb_L, h_L, Q_L, A_L, Phi_L,
				          V_L, c_L, I_L, constVariables);
		    //printf("upstream: j= %d, I_L= %f, h=%f,hB=%f Q=%f, QB=%f\n", j, I_L, p_ups_h_device[j], h_L, p_ups_Q_device[j], Q_L);
			
		}
		else
		{
			zb_L=zb_p_device[tid-1];
			sita_L=sita_p_device[tid-1];
			A_L=A_p_device[tid-1];
			Q_L=Q_p_device[tid-1];
			h_L=h_p_device[tid-1];

			flowVariables_calculation(a_speed, P_diameter_device[j], sita_L, A_L, Q_L, h_L, V_L, c_L, I_L, Phi_L, constVariables);

			
		}
		if(tid==valve_position+1 && T>T_valveOn)
		{
			
			zb_L=zb_p_device[valve_position+1];
			sita_L=sita_p_device[valve_position+1];
			A_L=A_p_device[valve_position+1];
			Q_L= - Q_p_device[valve_position+1];
			h_L=h_p_device[valve_position+1];
			flowVariables_calculation(a_speed, P_diameter_device[j], sita_L, A_L, Q_L, h_L, V_L, c_L, I_L, Phi_L, constVariables);
		}
		
		zbfw=(zb_L+zb_R)/2.0;
        HLLFluxCalculation(tid, P_diameter_device[j], a_speed, Phi_L, Phi_R, V_L, V_R, h_L, h_R, c_L, c_R, A_L, A_R, Q_L, Q_R,
			                I_L, I_R, Fw, constVariables); 
		
		// 	if (tid==17)
		// {
		// 	printf("***********************west face*******************:\n");
		// 	printf("tid= %d, V_L=%f, V_R=%f, h_L=%f, h_R=%f, c_L=%f, c_R=%f, A_L=%f, A_R=%f, Q_L=%f, Q_R=%f, I_L=%f,I_R=%f, Fe0=%f, Fe1=%f\n",tid ,V_L, V_R, h_L, h_R, c_L, c_R, A_L, A_R, Q_L, Q_R,
		// 	                I_L, I_R, Fw[0], Fw[1]);
		// }
		
		if (boundaryCellJugement(num_Pipe,upsBCell_device,tid))
		{
			pipe_flux_ups1[j]=Fw[0];
			pipe_flux_ups2[j]=Fw[1];
		/*	printf("pipe_flux_ups1[%d]=%f\n", tid, pipe_flux_ups1[j]);
			printf("pipe_flux_ups2[%d]=%f\n", tid, pipe_flux_ups2[j]);*/
		}

		//printf("j= %d, pipe_flux_ups1= %f, pipe_flux_ups2= %f\n",j,pipe_flux_ups1[j],pipe_flux_ups2[j]);
		//printf("j= %d, pipe_flux_downs1= %f, pipe_flux_downs2= %f\n",j,pipe_flux_downs1[j],pipe_flux_downs2[j]);

		// calculate the bed elevation source term
		s1=0;
		s2=-g*A_p_device[tid]*(zbfe-zbfw)/dx_cell_device[tid];
		//printf("s2= %f\n",s2);

		// calculate new value of A and Q
		dAdt=(Fe[0]-Fw[0])/dx_cell_device[tid]-s1;
	    An_p_device[tid]=A_p_device[tid]-dt*dAdt;


   //     if (tid==51 || tid==103 || tid==110 || tid==147){
			// printf("tid= %d, Fe[0]= %f, Fw[0]= %f\n",tid, Fe[0],Fw[0]);
			// printf("tid= %d, Fe[1]= %f, Fw[1]= %f\n",tid, Fe[1],Fw[1]);
			// printf("tid= %d, dAdt= %f\n",tid, dAdt);
			// printf("tid= %d, An_p_device[%d]= %f, A_p_device[%d]= %f\n",tid, tid, An_p_device[tid], tid, A_p_device[tid]);
			// printf("\n");
	  //   }


	 //    if (tid==2)
		// {
			
		// 	An_p_device[tid]=An_p_device[tid]-dt/dx_cell_device[tid]*(-pipe_flux_downs1[0]);
		// 	printf("tid=%d, Tjunction: An=%f, A=%f, dx=%f, flux=%f\n ", tid, An_p_device[tid], A_p_device[tid], dx_cell_device[tid], pipe_flux_downs1[0]);
		// 	double massChange2=(An_p_device[tid]-A_p_device[tid])*dx_cell_device[tid];
		// 	printf("mass changed: tid=%d, dx=%f, (An-A)*dx=%f\n", tid, dx_cell_device[tid], massChange2);
		// }
		// if (tid==3)
		// {
			
		// 	An_p_device[tid]=An_p_device[tid]-dt/dx_cell_device[tid]*(-pipe_flux_downs1[1]);
		// 	printf("tid=%d Tjunction: An=%f, A=%f, dx=%f, flux=%f\n ", tid, An_p_device[tid], A_p_device[tid], dx_cell_device[tid], pipe_flux_downs1[0]);
		// 	double massChange3=(An_p_device[tid]-A_p_device[tid])*dx_cell_device[tid];
		// 	printf("mass changed: tid=%d, dx=%f, (An-A)*dx=%f\n", tid, dx_cell_device[tid], massChange3);
		// }

		//dqxdt=(Fe[1]-Fw[1])/dx_cell_device[tid]-s2;
		//Qn_p_device[tid]=Q_p_device[tid]-dt*dqxdt;

	 //   		if (tid==17)	
		// {
		// 	printf("tid= %d, zbfe= %f, zbfw=%f, Fe[0]= %f, Fe[1]=%f, Fw[0]=%f, Fw[1]=%f\n", tid, zbfe, zbfw, Fe[0], Fe[1], Fw[0], Fw[1]);
		// 	printf("***********************integration*******************:\n");
		// 	printf("tid= %d, dx_cell_device=%f, s1=%f, A_p_device=%f, An_p_device=%f, dqxdt=%f, s2=%f, Q_p_device=%f, Qn_p_device=%f\n",tid ,
		// 		dx_cell_device[tid], s1, A_p_device[tid], An_p_device[tid], dqxdt, s2, Q_p_device[tid], Qn_p_device[tid]);
		// }


		dqxdt=(Fe[1]-Fw[1])/dx_cell_device[tid];
		pipeFrictionUpdate(tid, dt, pipe_MN_device[j], P_diameter_device[j], dqxdt, s2,  A_p_device[tid], Q_p_device[tid],
			                 h_p_device[tid], Qn_p_device[tid], constVariables);
		//printf("tid = %d, Q_p= %f, Qn_p= %f\n", tid,Q_p_device[tid], Qn_p_device[tid]);


		// if (tid<50)
		// {
		//     printf("tid= %d, Fe[0]= %f, Fw[0]= %f\n",tid, Fe[0],Fw[0]);
		// 	printf("tid= %d, Fe[1]= %f, Fw[1]= %f\n",tid, Fe[1],Fw[1]);
		// 	printf("tid= %d, dAdt= %f, dqxdt= %f\n",tid, dAdt,dqxdt);
		// 	printf("tid= %d, An_p_device[%d]= %f, A_p_device[%d]= %f\n",tid, tid, An_p_device[tid], tid, A_p_device[tid]);
		// 	printf("tid= %d, Qn_p_device[%d]= %f, Q_p_device[%d]= %f\n",tid, tid, Qn_p_device[tid], tid, Q_p_device[tid]);
			
		// 	printf("\n");
	 //    }
		// if (tid==23)
		// {
		// 	 printf("A_p[%d]=%f, An_p[%d]=%f\n", tid, A_p_device[tid], tid, An_p_device[tid]);
  //            printf("Q_p[%d]=%f, Qn_p[%d]=%f\n", tid, Q_p_device[tid], tid, Qn_p_device[tid]);

		// }
       
		
	    tid += blockDim.x * gridDim.x;
		
	}
}

//---------------------------------------------------------------------------------------------------
__device__ void TJunctionBoundary(double a_speed, int pipe_index, int bound_pipe_index, double d1, double d2, double *h_p_device, double *A_p_device, 
	                              double *sita_p_device, double *P_diameter_device, double &zb_R, double &h_R, double &A_R, double &Q_R, 
								  double &V_R, double &c_R, double &I_R, double &Phi_R, Const_variables constVariables){
	//zb_R=0.235; // this should be noticed
	double sita_R;
	int j=pipe_index;  // current pipe index
	int i=bound_pipe_index; // boundary pipe cell index
	h_R=h_p_device[i];  
	//printf("i=%d, h_R=%f",i, h_R);
	//double d1=0.19;
	//double d2=0.19;  // !!!!!!!!!!!!!!!!!!!!!!!!
	double Ap1=constVariables.pi*d1*d1/4.0; // this should be noticed
	double Ap2=constVariables.pi*d2*d2/4.0;
	double offset=d1/2.0-d2/2.0;
	if (j==0)
	{
		zb_R=0.19416;
		//zb_R=0.235;
	}
		
	if (j==1)
	{
		zb_R=0.095787;
		//zb_R=0.235;
	}
		
	
	if (offset<0)
	{
		printf("Error: the offset of the T junction node contecting pipe %d becomes negative\n", j);
	}
	else if (offset ==0)
	{
		sita_R=sita_p_device[i];
		A_R=A_p_device[i];
		Q_R=0;
		flowVariables_calculation(a_speed, P_diameter_device[j], sita_R, A_R, Q_R, h_R, V_R, c_R, I_R, Phi_R, constVariables);
	}
	else
	{
		if (h_R<=offset)
		{
			//printf("*$$$$$$$$$$$$$$$$$$$$$$$$$$$*: \n");
			h_R=0;
			sita_R=0;
			A_R=0;
			Q_R=0;
			flowVariables_calculation(a_speed, P_diameter_device[j], sita_R, A_R, Q_R, h_R, V_R, c_R, I_R, Phi_R, constVariables);
		}
		else if (h_R<=P_diameter_device[j]+offset && h_R > offset)
		{
			//printf("******************************: \n");
			h_R=h_p_device[i]-offset;
			sita_R=2*acos(1-2.0*h_R/P_diameter_device[j]);
			A_R=1/8.0*(sita_R-sin(sita_R))*P_diameter_device[j]*P_diameter_device[j];
			Q_R=0;
			flowVariables_calculation(a_speed, P_diameter_device[j], sita_R, A_R, Q_R, h_R, V_R, c_R, I_R, Phi_R, constVariables);
			//printf("aaa\n");
		}
		else
		{
			Q_R=0.0;
			V_R=0.0;
			sita_R=2*constVariables.pi;
			//A_R=constVariables.pi*P_diameter_device[j]*P_diameter_device[j]/4.0;
			c_R=a_speed;
					
			if (A_p_device[i]>Ap1)  // this should be noticed, means p1 become pressurized
			{
				//printf("111111\n");
				//printf("h_R= %f, A_p_device[%d]=%f\n", h_R, i, A_p_device[i]);
				double H1 =(pow(a_speed,2)/constVariables.g)*((A_p_device[i]-Ap1)/Ap1);
				double H2 =H1+offset;
				I_R=constVariables.pi/4.0*constVariables.g*pow(d2,2)*(H2+d2/2.0);
				A_R=Ap2+Ap2*constVariables.g*H2/pow(a_speed,2);
				//printf("H1= %f, H2=%f, I_R=%f, A_R=%f\n", H1, H2, I_R, A_R);
			
			}
			else
			{
				//printf("222222\n");
				if (h_R>d1)
					printf("error: in the T node, the water in the side pipe cell % geos beyongd its diameter\n",i);
				else
				{
					double H2=h_R-d2-offset;
				    if (H2<0)
						printf("error: the H of pipe %d becomes negative: H: %f, h_R: %f, d2: %f, offset: %f\n", j, H2, h_R, d2, offset);
					else
					{
						I_R=constVariables.pi *constVariables.g*d2*d2/4.0*(H2+d2/2.0);
						A_R=Ap2+Ap2*constVariables.g*H2/pow(a_speed,2);
						//printf("h_R= %f, A_p_device[%d]=%f\n", h_R, i, A_p_device[i]);
						//printf("H2=%f, I_R=%f, A_R=%f\n", H2, I_R, A_R);
					}
					
				}
				//if (I_R<constVariables.pi*constVariables.g*(d2/2.0)*d2*d2/4.0)
				//	printf("error: in the T node the pipe %d 's boundary pressure is smaller than its hydrostatic pressure\n",j);
				//double H2=4.0*I_R/(constVariables.pi*constVariables.g*0.085*0.085)-0.5*0.085;
			
			}
			Phi_R=a_speed*(log(A_R)-log(Ap2))+6.41*sqrt(constVariables.g*P_diameter_device[j]/8.0)*sin(2*constVariables.pi/4.0);
			//printf("bbb: h_R= %f, A=%f, c=%f, I=%f, Phi_R=%f, sita=%f\n", h_R, A_R, c_R, I_R, Phi_R, sita_p_device[i]);

		}
		//printf("pipe_index:%d, pipe_boundary_index:%d, h_boundary=%f, pipe_d=%f,  sita=%f, A=%f\n", j, i, h_p_device[i], P_diameter_device[j], sita_R, A_R);

	}

}


//---------------------------------------------------------------------------------------------------
__device__ int findPipeOrder(int pipe_num, int *downstreamCell, int temp1){
	int i;
	for(i=0;i<pipe_num;i++)
	{
		if (temp1<=downstreamCell[i])
			break;
	}
	return i;

}

//---------------------------------------------------------------------------------------------------
__device__ void flowVariables_calculation(double a_speed, double P_diameter, double sita, double A, double Q, double h, double &V,
	                                      double &c, double &I, double &Phi, Const_variables constVariables){

	double T, H;
	double Ap=constVariables.pi*pow(P_diameter,2)/4.0;

	if (A<constVariables.kesi)
	{
		V=0.0;
	}
	else 
	{
		V=Q/A;
	}

	T=P_diameter*sin(sita/2.0);

	if (T<constVariables.kesi && h<constVariables.kesi)
	{
		c=0.0;
	}   
    else
	{
		c=sqrt(constVariables.g*A/T);
        c=constVariables.MIN(c,a_speed);
	}

	if (A>Ap)
	{
		H =(pow(a_speed,2)/constVariables.g)*((A-Ap)/Ap);
		I=constVariables.pi/4.0*constVariables.g*pow(P_diameter,2)*(H+P_diameter/2.0);
	}
	else
	{
		H=0.0;
		I=1.0/24.0*(3.0*sin(sita/2.0)-pow(sin(sita/2.0),3)-3.0*(sita/2.0) * cos(sita/2.0))*constVariables.g*pow(P_diameter,3);
	}
		
	if (A<=Ap)
	{
		Phi=6.41*sqrt(constVariables.g*P_diameter/8.0)*sin(sita/4.0);
	}   
    else
	{
		Phi=a_speed*(log(A)-log(Ap))+6.41*sqrt(constVariables.g*P_diameter/8.0)*sin(2*constVariables.pi/4.0);
	}


}


//---------------------------------------------------------------------------------------------------
__device__ int boundaryCellJugement(int pipe_num, int *BoundaryCell, int temp1){
	
	
	int i;
	for(i=0;i<pipe_num;i++)
	{
		if (temp1==BoundaryCell[i])
			return 1;
	}
	return 0;
}

//---------------------------------------------------------------------------------------------------
__device__ void pipeToJuncBoundary(int j, double a_speed, double *h_junc, double *Q_junc, double h_elev, double P_diameter, double &zb_b, double &h_b,
	                          double &Q_b, double &A_b, double &Phi_b, double &V_b, double &c_b, double &I_b, Const_variables constVariables){

    double Ap=constVariables.pi*pow(P_diameter,2)/4.0;
	double sita_b, T_b;
	double H_junction;


	zb_b=h_elev;
	h_b=h_junc[j];
	Q_b=Q_junc[j]*P_diameter;

	if (h_b>=P_diameter)
	{
		A_b=Ap;   // A_b will be revised in the next several steps
        sita_b=2.0*constVariables.pi;
	}  
	else
	{
		sita_b=2.0*acos(1-2.0*h_b/P_diameter);
		A_b=1/8.0*(sita_b-sin(sita_b))*P_diameter*P_diameter;
	}
    
	if (A_b<constVariables.kesi)
	{
		V_b=0.0;
	}
    else
	{
		V_b=Q_b/A_b;   // V_L will be revised in the next several steps if h_junc>P_diameter
	}
     
    T_b=P_diameter*sin(sita_b/2.0);

	if (T_b<constVariables.kesi && h_b<constVariables.kesi)
	{
		c_b=0.0;
	}
    else
	{
		c_b=sqrt(constVariables.g*A_b/T_b);
        c_b=constVariables.MIN(c_b,a_speed);
	}
   
	if (h_junc[j]>=P_diameter)
	{
		I_b=constVariables.g*(h_junc[j]-P_diameter/2.0)*Ap;
		H_junction=I_b*4.0/(constVariables.pi*constVariables.g*pow(P_diameter,2))-P_diameter/2.0;
		A_b=H_junction*constVariables.g*Ap/pow(a_speed,2)+Ap;
		V_b=Q_b/A_b;  // at this point, V_L is revised by the new H_Junction deduced from pressre in the junction 
	}
    else
	{
		I_b=1/24.0*(3.0*sin(sita_b/2.0)-pow(sin(sita_b/2.0),3)-3*(sita_b/2.0) 
				   *cos(sita_b/2.0))*constVariables.g*pow(P_diameter,3);
	}


	if (A_b<=Ap)
	{
		Phi_b=6.41*sqrt(constVariables.g*P_diameter/8.0)*sin(sita_b/4.0);
	}
	else
	{
		Phi_b=a_speed*(log(A_b)-log(Ap))+6.41*sqrt(constVariables.g*P_diameter/8.0)*sin(2*constVariables.pi/4.0);
	}
}

//---------------------------------------------------------------------------------------------------
__device__ void pipeToOutfBoundary(double T, int *outfSurfIndex,  double *surfWaterDepth,  Vector *surfQ,
	                               double *downNorm1_device, double *downNorm2_device,  int pipe_index, int downsPO_device, int boundType, 
								   int *interval_ptr_device, double *T_series_device, double *h_series_device, double *Q_series_device, 
								   int *T_start_device, int *h_start_device, int *Q_start_device, int *intervalLen_device, double a_speed,
								   double P_diameter, double h_p, double Q_p, double A_p, double Phi_p, double V_p, double c_p, double I_p, 
								   double zb_p, double &h_b, double &Q_b, double &A_b, double &Phi_b, double &V_b, double &c_b, double &I_b, 
								   double &zb_b, Const_variables constVariables){

    // pipe_index means the current pipe that connects an outfall
	// downsPO_device means the current outfall index

	//printf("boundType= %d\n",boundType);
	/*for (int i=0;i<5;i++)
	{
		printf("the %d outfdall and its surface water grid: %d and surface water %f\n", i, outfSurfIndex[i], surfWaterDepth[outfSurfIndex[i]]);
	    printf("Qx: %f  Qy:%f\n",  surfQx[outfSurfIndex[i]], surfQy[outfSurfIndex[i]]);
	}*/
	 if (T>200)
	 {
	 	if (pipe_index==2)
	 	{
	 		h_b=h_p;
			Q_b= - Q_p;
			A_b=A_p;
			Phi_b=Phi_p;
			V_b= - V_p;
			c_b=c_p;
			I_b=I_p;
			zb_b=zb_p;

	 	}
	 	if (pipe_index==3)
	 	{
	 		  zb_b=zb_p;
	 		  h_b=0.0;
	 		  Q_b=0.0;
			  double sita_b = 2*acos(1-2.0*h_b/P_diameter);
			  A_b = 1/8.0*(sita_b-sin(sita_b)) *P_diameter*P_diameter;
			  flowVariables_calculation(a_speed, P_diameter, sita_b, A_b, Q_b, h_b, V_b, c_b, I_b, Phi_b, constVariables);

	 	}

	 }
     else
     {
     		switch (boundType){
			  // 'rigid' boundary
			  case 100 :{
				h_b=h_p;
				Q_b= - Q_p;
				A_b=A_p;
				Phi_b=Phi_p;
				V_b= - V_p;
				c_b=c_p;
				I_b=I_p;
				zb_b=zb_p;

				break;
			   }

			  // 'open' boundary
			  case 200 :{
				h_b=h_p;
				Q_b=Q_p;
				A_b=A_p;
				Phi_b=Phi_p;
				V_b=V_p;
				c_b=c_p;
				I_b=I_p;
				zb_b=zb_p;
				break;
			   }
			 // 'h given' boundary
			  case 300 :{
				  /* printf("case 300:\n");
				   printf("T= %f, downsPO_device=%d, T_start_device=%d, h_start_device=%d, intervalLen_device=%d, interval_ptr_device=%d\n",
					     T, downsPO_device, T_start_device[downsPO_device],  h_start_device[downsPO_device], intervalLen_device[downsPO_device],
						  interval_ptr_device[downsPO_device]);*/
		          zb_b=zb_p;
				  if (T_start_device[downsPO_device]<0 || h_start_device[downsPO_device]<0 || intervalLen_device[downsPO_device]<0 || downsPO_device<0 )
					  printf("The %d outfall boundary is wrong (Boundary type: 300\n", downsPO_device);
				  else
				      SimpleDifference(T_series_device, h_series_device, T, T_start_device[downsPO_device], h_start_device[downsPO_device], 
					                   intervalLen_device[downsPO_device], interval_ptr_device[downsPO_device], h_b);
				 
				/*  printf("T= %f, T_start_device[%d]=%f, h_start_device[%d]=%f, intervalLen_device[%d]=%f, interval_ptr_device[%d]=%f, hb=%f\n",
					     T, downsPO_device, T_start_device[downsPO_device], downsPO_device, h_start_device[downsPO_device], downsPO_device, intervalLen_device[downsPO_device],
						 downsPO_device, interval_ptr_device[downsPO_device], h_b);*/
				
				  Q_b=Q_p;
				  double sita_b = 2*acos(1-2.0*h_b/P_diameter);
				  A_b = 1/8.0*(sita_b-sin(sita_b)) *P_diameter*P_diameter;
				  flowVariables_calculation(a_speed, P_diameter, sita_b, A_b, Q_b, h_b, V_b, c_b, I_b, Phi_b, constVariables);
				 // printf("P_diameter=%f, sita_b= %f, A_b=%f, Q_b=%f, V_b=%f, c_b=%f, I_b=%f, phi_b=%f\n", P_diameter, sita_b, A_b, Q_b, V_b, c_b, I_b, Phi_b);
				 // printf("P_diameter=%f, h_p= %f, A_p=%f, Q_p=%f, V_p=%f, c_p=%f, I_p=%f, phi_p=%f\n", P_diameter, h_p, A_p, Q_p, V_p, c_p, I_p, Phi_p);

				  break;

			   }
		      // 'Q given' boundary
			  case 400 :{
			  	  zb_b=zb_p;
				  
				  if (T_start_device[downsPO_device]<0 || Q_start_device[downsPO_device]<0 || intervalLen_device[downsPO_device]<0 || downsPO_device<0 )
					  printf("The %d outfall boundary is wrong (Boundary type: 400\n", downsPO_device);
				  else
				      SimpleDifference(T_series_device, Q_series_device, T, T_start_device[downsPO_device], Q_start_device[downsPO_device], 
					                   intervalLen_device[downsPO_device], interval_ptr_device[downsPO_device], Q_b);
			      
				  h_b=h_p;
				  double sita_b = 2*acos(1-2.0*h_b/P_diameter);
				  A_b = 1/8.0*(sita_b-sin(sita_b)) *P_diameter*P_diameter;
				  flowVariables_calculation(a_speed, P_diameter, sita_b, A_b, Q_b, h_b, V_b, c_b, I_b, Phi_b, constVariables);
				 // printf("P_diameter=%f, sita_b= %f, A_b=%f, Q_b=%f, V_b=%f, c_b=%f, I_b=%f, phi_b=%f\n", P_diameter, sita_b, A_b, Q_b, V_b, c_b, I_b, Phi_b);
				 // printf("P_diameter=%f, h_p= %f, A_p=%f, Q_p=%f, V_p=%f, c_p=%f, I_p=%f, phi_p=%f\n", P_diameter, h_p, A_p, Q_p, V_p, c_p, I_p, Phi_p);

				  break;

			   }
			  // 'hQ given boundary'
			  case 500 :{
				 // printf("Q_series_device[%d]= %f, Q_series_device[%d]=%f, Q_series_device[%d]=%f, Q_series_device[%d]=%f\n",0, Q_series_device[0], 1, Q_series_device[1], 2, Q_series_device[2], 3, Q_series_device[3]);
				  if (T_start_device[downsPO_device]<0 || h_start_device[downsPO_device]<0 || Q_start_device[downsPO_device]<0 ||intervalLen_device[downsPO_device]<0 || downsPO_device<0 )
					  printf("The %d outfall boundary is wrong (Boundary type: 500)\n", downsPO_device);
				  else
				  {
					  SimpleDifference(T_series_device, h_series_device, T, T_start_device[downsPO_device], h_start_device[downsPO_device], 
					                   intervalLen_device[downsPO_device], interval_ptr_device[downsPO_device], h_b);
					  SimpleDifference(T_series_device, Q_series_device, T, T_start_device[downsPO_device], Q_start_device[downsPO_device], 
					                   intervalLen_device[downsPO_device], interval_ptr_device[downsPO_device], Q_b);
				  }
				  zb_b=zb_p;
				  double sita_b = 2*acos(1-2.0*h_b/P_diameter);
				  A_b = 1/8.0*(sita_b-sin(sita_b)) *P_diameter*P_diameter;
				  flowVariables_calculation(a_speed, P_diameter, sita_b, A_b, Q_b, h_b, V_b, c_b, I_b, Phi_b, constVariables);

				   break;

			    }

			  case 600 :{
				  h_b=surfWaterDepth[outfSurfIndex[downsPO_device]];
				  //Q_b=0;
				  double surfQx = surfQ[outfSurfIndex[downsPO_device]].x;
				  double surfQy = surfQ[outfSurfIndex[downsPO_device]].y;
				  Q_b=surfQx*(-downNorm1_device[pipe_index]*1+(-downNorm2_device[pipe_index])*0)+
					  surfQy*(-downNorm1_device[pipe_index]*0+(-downNorm2_device[pipe_index])*1);
				/*  printf("outfall: %d, surfQx=%f, surfQy=%f, downNorm1_device=%f, downNorm2_device=%f\n  ", downsPO_device,
					    surfQx[outfSurfIndex[downsPO_device]], surfQy[outfSurfIndex[downsPO_device]], downNorm1_device[pipe_index], downNorm2_device[pipe_index]);*/

				  double sita_b = 2*acos(1-2.0*h_b/P_diameter);
				  A_b = 1/8.0*(sita_b-sin(sita_b)) *P_diameter*P_diameter;
				  flowVariables_calculation(a_speed, P_diameter, sita_b, A_b, Q_b, h_b, V_b, c_b, I_b, Phi_b, constVariables);
				  zb_b=zb_p;
				  //printf("Boundary type:600 , downsPO_device= %d, outfSurfIndex=%d\n", downsPO_device, outfSurfIndex[downsPO_device]);
				   break;
				}

			}

     }



}

//---------------------------------------------------------------------------------------------------
__device__ void SimpleDifference(double *timeSeries, double *goalSeries, double T, int timeStart, int goalStart, 
	                             int seriesLength, int &interval_ptr_device,  double &goal_value){


	 
    int timeEnd=seriesLength-1;
	int goalEnd=seriesLength-1;
	int n=interval_ptr_device;

		  if (T>=timeSeries[timeEnd])
		  {
			  goal_value=goalSeries[goalEnd];
			 
		  }
		  else
		  {

			  if (T>=timeSeries[timeStart+n])
			  {

				  for (int i=timeStart+n+1;i<seriesLength;i++)
				  {
					  if (T<timeSeries[i])
					  {
						  n=i-timeStart;
						  break;
					  }
				  }
			  }
			  double T_left=timeSeries[timeStart+n-1];
			  double T_right=timeSeries[timeStart+n];
			  double g_left=goalSeries[goalStart+n-1];
			  double g_right=goalSeries[goalStart+n];
		      goal_value=(g_right-g_left)*(T-T_left)/(T_right-T_left) + g_left;
			  //printf("n= %d, start= %d, T_L=%f, T_R=%f, h_L=%f, h_R=%f\n", n, timeStart, T_left, T_right, g_left, g_right);
		  }
		  interval_ptr_device=n;
		 // printf("T= %f, differetial_value=%f\n", T, goal_value);



}

//---------------------------------------------------------------------------------------------------
__device__ void HLLFluxCalculation(int tid, double P_diameter,double a_speed, double Phi_L,double Phi_R,double V_L, double V_R,
	                                double h_L, double h_R, double c_L, double c_R, double A_L, double A_R, 
									double Q_L, double Q_R, double I_L, double I_R, double *F, Const_variables constVariables){
    

	double Phi_star, sita_star, c_star, V_star;
	double S_L, S_R, F_L1, F_L2, F_R1, F_R2;

    if (h_L<constVariables.kesi && h_R<constVariables.kesi)
	{
		F[0]=0;
		F[1]=0;
		//printf("tid=%d F[0]= %f, F[1]=%f \n",tid, F[0], F[1]);
	}
	else
	{
		Phi_star=0.5*(Phi_L+Phi_R)+0.5*(V_L-V_R);
		sita_star=4*asin(Phi_star/(6.41*sqrt(constVariables.g*P_diameter/8.0)));
	
		if (Phi_star>=(6.41*sqrt(constVariables.g*P_diameter/8.0)))
		{
			c_star=a_speed;
		}
		else
		{
			c_star=sqrt((constVariables.g*P_diameter*(sita_star-sin(sita_star)))/(8.0*sin(sita_star/2.0)));
			c_star=constVariables.MIN(c_star,a_speed);
		}

		 V_star=0.5*(V_L+V_R)+0.5*(Phi_L-Phi_R);

		if (h_L>constVariables.kesi && h_R>constVariables.kesi)
		{
			S_L=constVariables.MIN(V_L-c_L,V_star-c_star);
			S_R=constVariables.MAX(V_star+c_star,V_R+c_R);
			F_L1=Q_L;
			F_L2=Q_L*V_L+I_L;
			F_R1=Q_R;
			F_R2=Q_R*V_R+I_R;

		}

	

		if (h_R<constVariables.kesi && h_L>constVariables.kesi)
		{
			S_L=V_L-c_L;
			S_R=V_L+Phi_L;
			F_L1=Q_L;
			F_L2=Q_L*V_L+I_L;
			F_R1=0;
			F_R2=0;
		}
		
		if (h_L<constVariables.kesi && h_R>constVariables.kesi)
		{
			S_L=V_R-Phi_R;
			S_R=V_R+c_R;
			F_L1=0;
			F_L2=0;
			F_R1=Q_R;
			F_R2=Q_R*V_R+I_R;
		}

		if (S_L>=0)
		{
			F[0]=F_L1;
			F[1]=F_L2;
		}
		else
		{
			if (S_R>=0)
			{
				F[0]=(S_R*F_L1-S_L*F_R1+S_L*S_R*(A_R-A_L))/(S_R-S_L);
				F[1]=(S_R*F_L2-S_L*F_R2+S_L*S_R*(Q_R-Q_L))/(S_R-S_L);
			}
			else
			{
				F[0]=F_R1;
				F[1]=F_R2;
			}
		}		
	//printf("tid= %d, phi_star=%f, c_star=%f,  sita_star=%f, V_star=%f, s_L=%f, s_R=%f,\n",tid, Phi_star, c_star, sita_star, V_star,S_L, S_R);
	}
}


//---------------------------------------------------------------------------------------------------
__device__ void pipeFrictionUpdate(int index, double dt, double P_manning, double P_diameter, double F_n, double s2, 
	                                double A_p_device, double Q_p_device, double h_p_device, double &Qn_p_device,
									Const_variables constVariables){

	double B;
	double SGN;
	double numerator_p;
	double denominator_p;
	double Ap=constVariables.pi*P_diameter*P_diameter/4.0;
	double sita, P_n;


	B=-Q_p_device + dt*F_n - dt*s2;

	if (index==174|| index==191)
		P_manning=P_manning*1.1;

	if (dt<constVariables.kesi)
	{
		Qn_p_device=Q_p_device;

	}
	else
	{
		if (B>0)
		{
			SGN=-1.0;
		}
		else if (B<=0)
		{
			SGN=1.0;
		}

		if (A_p_device<=Ap)
			 sita=2*acos(1-2*h_p_device/P_diameter);
		else
			sita=2*constVariables.pi;
		P_n=0.5*sita*P_diameter;



		if (P_n<constVariables.kesi)
		{
			Qn_p_device=Q_p_device-dt*F_n + dt*s2;
		}
		else
		{
			numerator_p=-pow(A_p_device,7.0/3.0)+sqrt(pow(A_p_device,14.0/3.0)+4*dt*abs(B)*constVariables.g*pow(P_manning,2)*pow(P_n,4.0/3.0)*pow(A_p_device,7.0/3.0));
			denominator_p=2.0*dt*SGN*constVariables.g*pow(P_manning,2)*pow(P_n,(4/3.0));
			Qn_p_device=numerator_p/denominator_p;
		}

		if (abs(B)<constVariables.kesi)
		{
			Qn_p_device=Q_p_device-dt*F_n+dt*s2;
		}

	}
	// if (index==51 || index==103)
 //    {
 //    	printf("index= %d, numerator=%f,  denominator= %f \n",index, numerator_p,denominator_p);
	//     printf("index=%d, B= %f, P_manning= %f, Pn= %f, SGN= %f\n", index, B, P_manning, P_n, SGN );
	//     printf("index=%d, Fn=%f, sb=%f, Qx_pre=%f, Q_new=%f\n", index, F_n, s2, Q_p_device, Qn_p_device );

 //    }
	

}

//----------------------------------------------------------------------------------------------------------------------
__global__ void pipe_updateKernal(int CN, int num_Pipe, double a_speed, double dT, double *P_diameter_device, double *A_p_device, 
	                        double *An_p_device, double *Q_p_device,
	                        double *Qn_p_device, double *h_p_device, double *hn_p_device, double *sita_p_device, 
							double *sitan_p_device, int *downsBCell_device, double *pipe_Ap_device, 
							double *dx_cell_device, double *t_pipe_device,
							double *pipe_flux_downs1, Const_variables constVariables){
	
	/*Input_information inputKeyPtr;
							  Pipe_variables pipeVariables;
							  Junc_variables juncVariables;
							  pipe_updateKernal(pipeVariables.totalCellNum, inputKeyPtr.num_Pipe,
								          inputKeyPtr.a_speed, pipeVariables.P_diameter_device,
										  pipeVariables.A_p_device, pipeVariables.An_p_device,
										  pipeVariables.Q_p_device, pipeVariables.Qn_p_device,
										  pipeVariables.h_p_device, pipeVariables.hn_p_device,
										  pipeVariables.sita_p_device, pipeVariables.sitan_p_device,
										  pipeVariables.downsBCell_device, pipeVariables.pipe_Ap_device,
										  pipeVariables.dx_cell_device, pipeVariables.t_pipe_device,
										  constVariables);
*/
	// update A, Q, h and sita  
    int i,j;
	//double Ap=pi*pow(P_diameter,2)/4.0;
	double sita_n, sita_k, sita_q;
	double c_current;
	int steps_sita;
	int tid= blockDim.x * blockIdx.x + threadIdx.x;  

	double kesi=constVariables.kesi;
	double g=constVariables.g;
	double pi=constVariables.pi;
	int Maxtrials=constVariables.Maxtrials;

	while (tid< CN)
	{
		
		j=findPipeOrder(num_Pipe,downsBCell_device,tid);


		if (tid==110)
		{
			
			An_p_device[tid]=An_p_device[tid]-dT/dx_cell_device[tid]*(-pipe_flux_downs1[0]);
			//printf("Tjunction: An=%f, A=%f, dx=%f, flux=%f\n ", An_p_device[tid], A_p_device[tid], dx_cell_device[tid], pipe_flux_downs1[0]*10000000);
			//double massChange2=(An_p_device[tid]-A_p_device[tid])*dx_cell_device[tid];
			//printf("mass changed: tid=%d, (An-A)*dx=%f\n", tid, massChange2);
		}
		if (tid==147)
		{
			
			An_p_device[tid]=An_p_device[tid]-dT/dx_cell_device[tid]*(-pipe_flux_downs1[1]);
			//printf("Tjunction: An=%f, A=%f, dx=%f, flux=%f\n ", An_p_device[tid], A_p_device[tid], dx_cell_device[tid], pipe_flux_downs1[0]*10000000);
			//double massChange3=(An_p_device[tid]-A_p_device[tid])*dx_cell_device[tid];
			//printf("mass changed: tid=%d, (An-A)*dx=%f\n", tid, massChange3);
		}




		
		A_p_device[tid]=An_p_device[tid];
        Q_p_device[tid]=Qn_p_device[tid];
		//printf("tid= %d, A_p_device= %f, Q= %f\n",tid,A_p_device[tid], Q_p_device[tid]);







		steps_sita=0;
		if (A_p_device[tid]<kesi)
		{
			sitan_p_device[tid]=0;
		}
		else if (A_p_device[tid]>= pipe_Ap_device[j])
		{
			sitan_p_device[tid]=2*pi;
		}
		else
		{
			if (sita_p_device[tid]<kesi)
			{
				sita_k=1.57;
			}
			else if (abs(sita_p_device[tid]-2*pi)<1e-6)
			{
				sita_k=2*pi-0.5;
			}
			else
			{
				sita_k=sita_p_device[tid];
			}

			while (steps_sita<Maxtrials)
			{
				steps_sita=steps_sita+1;
				sita_q=sita_k-(sita_k-sin(sita_k)-2*pi*A_p_device[tid]/pipe_Ap_device[j])/(1-cos(sita_k));

				if (abs(sita_q-sita_k)<=0.001*abs(sita_k))
				{
					sitan_p_device[tid]=sita_q;
					break;
				}
				sita_k=sita_q;
			}
		}

		hn_p_device[tid]=0.5*(1-cos(sitan_p_device[tid]/2.0))*P_diameter_device[j];
	
		sita_p_device[tid]=sitan_p_device[tid];
	    h_p_device[tid]=hn_p_device[tid];

		//if (tid==51 || tid==103 || tid==110 || tid==147)
		//    printf("tid= %d, h= %f, Q= %f\n",tid,h_p_device[tid], Q_p_device[tid]);


		//	printf("hn_p_device[%d]= %f\n",tid,hn_p_device[tid]);

		// calculate the dt for each pipe cell
		//if (An_p_device[tid]>=(pi*pow(P_diameter_device[j],2)/4.0))
		//{
		//	t_pipe_device[tid]=(dx_cell_device[tid]/2.0)/(abs(Qn_p_device[tid]/An_p_device[tid])+a_speed);
		//}
		//else
		//{
		//	// sqrt((g*P_diameter*(sita_star-sin(sita_star)))/(8.0*sin(sita_star/2.0)))
		//	c_current=sqrt((g*P_diameter_device[j]*(sita_p_device[tid]-sin(sita_p_device[tid])))/(8.0*sin(sita_p_device[tid]/2.0)));
		//	t_pipe_device[tid]=(dx_cell_device[tid]/2.0)/(abs(Qn_p_device[tid]/An_p_device[tid])+c_current);
		//}
		if (h_p_device[tid]>kesi)
		{
			if (An_p_device[tid]>=(pi*pow(P_diameter_device[j],2)/4.0))
			{
				t_pipe_device[tid]=(dx_cell_device[tid]/2.0)/(abs(Qn_p_device[tid]/An_p_device[tid])+a_speed);
			/*		if (tid==1124)
				{
					printf("^^pipe %d, time: %f, dx=%f, Q=%f, A=%f, c=%f\n", tid, t_pipe_device[tid], dx_cell_device[tid], Qn_p_device[tid],  An_p_device[tid], a_speed);

				}*/
			}
			else
			{
				c_current=sqrt((g*P_diameter_device[j]*(sita_p_device[tid]-sin(sita_p_device[tid])))/(8.0*sin(sita_p_device[tid]/2.0)));
				//printf("c_current[%d]=%f\n", tid, c_current);
				t_pipe_device[tid]=(dx_cell_device[tid]/2.0)/(abs(Qn_p_device[tid]/An_p_device[tid])+c_current);
				/*if (tid==1124)
				{
					printf("**pipe %d, time: %f, dx=%f, Q=%f, A=%f, c=%f\n", tid, t_pipe_device[tid], dx_cell_device[tid], Qn_p_device[tid],  An_p_device[tid], c_current);

				}*/
			/*	if (t_pipe_device[tid]<0.1)
				{
					printf("pipe %d, time: %f, dx=%f, Q=%f, A=%f, c=%f\n", tid, t_pipe_device[tid], dx_cell_device[tid], Qn_p_device[tid],  An_p_device[tid], c_current);
				}*/
				/*if (t_pipe_device[tid]<0.3)
				printf("t_pipe_device[%d]=%f\n", tid, t_pipe_device[tid]);*/
			}


		}
		else
		{
			t_pipe_device[tid]=0.005;
		}

		//printf("tid= %d, c_current[%d]= %f\n", tid,tid, c_current );
		//printf("tid= %d, t_pipe_device[%d]= %f\n", tid,tid, t_pipe_device[tid] );
		//if (100<tid<200)
		//printf("tid= %d, t_pipe_device= %f\n", tid,t_pipe_device[tid] );
		//printf("tid= %d, h_pipe_device[tid]=%f\n",tid, h_p_device[tid]);
		tid += blockDim.x * gridDim.x;
		
      
	}
	
	
	
	
}

//----------------------------------------------------------------------------------------------------------------------
__global__ void junction_calculationKernal(int N, int width_upsJP, int width_downsJP, double dt, double *junc_MN_device, 
	                                 double *pipe_Ap_device, 
	                                 double *P_diameter_device, double *junc_Area_device,
	                                 double *upNorm1_device, double *upNorm2_device, double *downNorm1_device, 
									 double *downNorm2_device, int *upsJP_device, int *downsJP_device, 
									 double *pipe_flux_ups1, double *pipe_flux_ups2, double *pipe_flux_downs1,
									 double *pipe_flux_downs2, double *h_j_device,double *hn_j_device, 
									 double *qx_j_device, double *qxn_j_device, double *qy_j_device, 
									 double *qyn_j_device, int *j_surfGird, double* Qexchange_device,
									 double *Ax_junc, double *Ay_junc, 
									 double *hydrostaticX, double *hydrostaticY, 
									 double *pipe_F1, double *pipe_F2,
									 double *pipe_F3,
									 Const_variables constVariables){

// __global__ void junction_calculationKernal(int N, int width_upsJP, int width_downsJP, double dt, double *junc_MN_device, 
// 	                                 double *pipe_Ap_device, 
// 	                                 double *P_diameter_device, double *junc_Area_device,
// 	                                 double *upNorm1_device, double *upNorm2_device, double *downNorm1_device, 
// 									 double *downNorm2_device, int *upsJP_device, int *downsJP_device, 
// 									 double *pipe_flux_ups1, double *pipe_flux_ups2, double *pipe_flux_downs1,
// 									 double *pipe_flux_downs2, double *h_j_device,double *hn_j_device, 
// 									 double *qx_j_device, double *qxn_j_device, double *qy_j_device, 
// 									 double *qyn_j_device, int *j_surfGird, double* Qexchange_device,
// 									 double *Ax_junc, double *Ay_junc, 
// 									 double *hydrostaticX, double *hydrostaticY, 
// 									 double *pipe_F1, double *pipe_F2,
// 									 double *pipe_F3,
// 									 Const_variables constVariables);



	int tid= blockDim.x * blockIdx.x + threadIdx.x;
//printf("N=%d, width_upsJP=%d, width_downsJp=%d, dt=%f\n",N, width_upsJP, width_downsJP, dt);
//for(int g=0;g<N;g++)
//{
//	printf("junc_MN_device[0]= %f, pipe_Ap_device[0]= %f\n", upNorm1_device[0], upNorm2_device[0]);
//}


	/*
	Input_information inputKeyPtr;
	Pipe_variables pipeVariables;
	Junc_variables juncVariables;
	junction_calculationKernal (inputKeyPtr.num_Junc, juncVariables.width_upsJP,
						  juncVariables.width_downsJP, dT,
						  juncVariables.junc_MN_device, pipeVariables.pipe_Ap_device,
						  pipeVariables.P_diameter_device, juncVariables.junc_Area_device,
						  pipeVariables.upNorm1_device, pipeVariables.upNorm2_device,
						  pipeVariables.downNorm1_device, pipeVariables.downNorm2_device, 
						  juncVariables.upsJP_device, juncVariables.downsJP_device,
						  pipeVariables.pipe_flux_ups1, pipeVariables.pipe_flux_ups2,
						  pipeVariables.pipe_flux_downs1, pipeVariables.pipe_flux_downs2,
						  juncVariables.h_j_device, juncVariables.hn_j_device,
						  juncVariables.qx_j_device, juncVariables.qxn_j_device,
						  juncVariables.qy_j_device, juncVariables.qyn_j_device,
						  juncVariables.surfGrid_device, constVariables);*/
	
										
	// double *pipe_F1= new double[N];
	// double *pipe_F2= new double[N];
	// double *pipe_F3= new double[N];
	// double *hydrostaticX= new double[N];
	// double *hydrostaticY= new double[N];
	// double *Ax_junc= new double[N];
	// double *Ay_junc= new double[N];
	
	//printf("downNorm1_device[0]= %f, downNorm2_device[0]= %f\n", downNorm1_device[0], downNorm2_device[0]);
	double sita_pipe, A_pipe, Pressure_pipe;
	int i, j;
	int sub;

	double g=constVariables.g;
	
	// for (i=0;i<N;i++)
	// {
	// 	pipe_F1[i]=0.00;
	// 	pipe_F2[i]=0.00;
	// 	pipe_F3[i]=0.00;
	// 	hydrostaticX[i]=0.00;
	// 	hydrostaticY[i]=0.00;
	// 	//printf("pipe_F1[%d]=%f\n",i,pipe_F1[i]);
	// }
	

	//printf("g=%f\n", constVariables.g);
	while (tid<N)   // in this function, N is the number of junction
	{

		i=tid;
		pipe_F1[i]=0.00;
		pipe_F2[i]=0.00;
		pipe_F3[i]=0.00;
		hydrostaticX[i]=0.00;
		hydrostaticY[i]=0.00;
		Ax_junc[i]=0.00;
		Ay_junc[i]=0.00;
		
		//printf("j_surfIndex[%d]=%d, juncVariables.Qexchange_device[%d]=%f\n,", tid, j_surfGird[tid], tid, Qexchange_device[tid]);
	
		// i means the order of the junction, j means the order of the connecting pipes
		for (j=0; j<width_upsJP; j++)
		{
			sub=i*width_upsJP+j;
			//printf("sub=%d , upsJP_device[sub]=%d \n",sub,upsJP_device[sub]);
			//printf("sub=%d , downsJP_device[sub]=%d \n",sub,downsJP_device[sub]);
			if (upsJP_device[sub]>=0)
			{
				
				pipe_F1[i]=pipe_F1[i]+pipe_flux_ups1[upsJP_device[sub]];
				pipe_F2[i]=pipe_F2[i]+pipe_flux_ups2[upsJP_device[sub]]*
					(upNorm1_device[upsJP_device[sub]]*1+upNorm2_device[upsJP_device[sub]]*0);
				pipe_F3[i]=pipe_F3[i]+pipe_flux_ups2[upsJP_device[sub]]*
					(upNorm1_device[upsJP_device[sub]]*0+upNorm2_device[upsJP_device[sub]]*1);	
				//printf("i=%d, sub=%d,  pipe_flux_ups1=%f, pipe_flux_ups2= %f\n", i, sub,pipe_flux_ups1[upsJP_device[sub]], pipe_flux_ups2[upsJP_device[sub]]);
				
				// pressure calculation
				if (h_j_device[i]<P_diameter_device[upsJP_device[sub]])
				{
					sita_pipe=2*acos(1-2*h_j_device[i]/P_diameter_device[upsJP_device[sub]]);
					A_pipe=1/8.0*(sita_pipe-sin(sita_pipe))*P_diameter_device[upsJP_device[sub]]*P_diameter_device[upsJP_device[sub]];
					Pressure_pipe=1/24.0*(3.0*sin(sita_pipe/2.0)-pow(sin(sita_pipe/2.0),3)-3.0*(sita_pipe/2.0)
                     *cos(sita_pipe/2.0))*g*pow(P_diameter_device[upsJP_device[sub]],3);
					
				}
				else
				{
					Pressure_pipe=g*(h_j_device[i]-P_diameter_device[upsJP_device[sub]]/2.0)*pipe_Ap_device[upsJP_device[sub]];
				}

				
				hydrostaticX[i]=hydrostaticX[i]+Pressure_pipe*
					(upNorm1_device[upsJP_device[sub]]*1+upNorm2_device[upsJP_device[sub]]*0);
				hydrostaticY[i]=hydrostaticY[i]+Pressure_pipe*
					(upNorm1_device[upsJP_device[sub]]*0+upNorm2_device[upsJP_device[sub]]*1);
				
			}
		}


		for (j=0; j<width_downsJP; j++)
		{
			sub=i*width_downsJP+j;
			if (downsJP_device[sub]>=0)
			{
				pipe_F1[i]=pipe_F1[i]+(-pipe_flux_downs1[downsJP_device[sub]]);  // noticed!!!!
				pipe_F2[i]=pipe_F2[i]+pipe_flux_downs2[downsJP_device[sub]]*
					(downNorm1_device[downsJP_device[sub]]*1+downNorm2_device[downsJP_device[sub]]*0);
				pipe_F3[i]=pipe_F3[i]+pipe_flux_downs2[downsJP_device[sub]]*
					(downNorm1_device[downsJP_device[sub]]*0+downNorm2_device[downsJP_device[sub]]*1);
			

				// pressure calculation

				if (h_j_device[i]<P_diameter_device[downsJP_device[sub]])
				{
					sita_pipe=2*acos(1-2*h_j_device[i]/P_diameter_device[downsJP_device[sub]]);
					A_pipe=1/8.0*(sita_pipe-sin(sita_pipe))*P_diameter_device[downsJP_device[sub]]*P_diameter_device[downsJP_device[sub]];
					Pressure_pipe=1/24.0*(3.0*sin(sita_pipe/2.0)-pow(sin(sita_pipe/2.0),3)-3.0*(sita_pipe/2.0)
                     *cos(sita_pipe/2.0))*g*pow(P_diameter_device[downsJP_device[sub]],3);
				}
				else
				{
					Pressure_pipe=g*(h_j_device[i]-P_diameter_device[downsJP_device[sub]]/2.0)*pipe_Ap_device[downsJP_device[sub]];
				}

				hydrostaticX[i]=hydrostaticX[i]+Pressure_pipe*
					(downNorm1_device[downsJP_device[sub]]*1+downNorm2_device[downsJP_device[sub]]*0);
				hydrostaticY[i]=hydrostaticY[i]+Pressure_pipe*
					(downNorm1_device[downsJP_device[sub]]*0+downNorm2_device[downsJP_device[sub]]*1);
			
				//printf("i= %d, j= %d, sub= %d, Pressure_pipe= %f\n", i,j,sub,Pressure_pipe);
			
			}
		}

		hydrostaticX[tid]=-hydrostaticX[tid];
		hydrostaticY[tid]=-hydrostaticY[tid];
	

		hn_j_device[tid]=h_j_device[tid]-dt/junc_Area_device[tid]*(pipe_F1[tid]);	

       // printf("before adding surface water: hn_j[%d]=%f \n", tid, hn_j_device[tid]);
       // add the water from surface
		hn_j_device[tid]=hn_j_device[tid]+(Qexchange_device[tid]*dt);
      //   printf("after adding surface water: hn_j[%d]=%f \n", tid, hn_j_device[tid]);

		//qxn_j_device[tid]=qx_j_device[tid]-dt/junc_Area_device[tid]*(pipe_F2[tid]+hydrostaticX[tid]); 	
		//qyn_j_device[tid]=qy_j_device[tid]-dt/junc_Area_device[tid]*(pipe_F3[tid]+hydrostaticY[tid]);
		
		Ax_junc[tid]=-(pipe_F2[tid]+hydrostaticX[tid])/junc_Area_device[tid];
		Ay_junc[tid]=-(pipe_F3[tid]+hydrostaticY[tid])/junc_Area_device[tid];

		juncFrictionUpdate(tid, N, dt, junc_MN_device[tid], Ax_junc[tid], Ay_junc[tid], qx_j_device[tid], 
			               qxn_j_device[tid], qy_j_device[tid], qyn_j_device[tid], h_j_device[tid], constVariables);

		
		/*printf("tid= %d, pipe_F1[tid]= %f\n",tid, pipe_F1[tid]);
		printf("tid= %d, pipe_F2[tid]= %f\n", tid,pipe_F2[tid]);
		printf("tid= %d, pipe_F3[tid]= %f\n", tid,pipe_F3[tid]);
		printf("tid= %d, hydrostaticX[tid]= %f\n", tid,hydrostaticX[tid]);
		printf("tid= %d, hydrostaticY[tid]= %f\n", tid,hydrostaticY[tid]);
		printf("tid= %d, junc_Area_device[tid]= %f\n", tid,junc_Area_device[tid]);*/
		if (tid<500)
		{
		printf("tid= %d, hn[%d]= %f\n", tid,tid, hn_j_device[tid]);
		printf("tid= %d, qxn[%d]= %f\n", tid,tid,qxn_j_device[tid]);
		printf("tid= %d, qyn[%d]= %f\n", tid,tid, qyn_j_device[tid]);
		}
		/*printf("tid= %d, j_surfGird[tid]= %d\n", tid,j_surfGird[tid]);
		printf("*****\n");*/
		
		tid += blockDim.x * gridDim.x;
	}

	// delete[] pipe_F1;
	// delete[] pipe_F2;
	// delete[] pipe_F3;
	// delete[] hydrostaticX;
	// delete[] hydrostaticY;
	// delete[] Ax_junc;
	// delete[] Ay_junc;
	
}

//----------------------------------------------------------------------------------------------------------------------
__device__ void juncFrictionUpdate(int index, int junc_num, double dt, double J_manning, double Ax_junc, double Ay_junc, double qx_j_device, 
	                                double &qxn_j_device, double qy_j_device, double &qyn_j_device, double h_j_device,
									Const_variables constVariables){
									
	double mx, my;
	//double numerator_x, denominator_x, numerator_y, denominator_y;
	double temp;
	double g=constVariables.g;

	//printf("friction: hj[%d]=%f\n", index, h_j_device);
	if (dt<constVariables.kesi)
	{
		qxn_j_device=qx_j_device;
		qyn_j_device=qy_j_device;

	}
	else
	{
		mx=qx_j_device+dt*Ax_junc;
		my=qy_j_device+dt*Ay_junc;

		if (h_j_device<constVariables.kesi)
		{
			qxn_j_device=qx_j_device+dt*Ax_junc;
			qyn_j_device=qy_j_device+dt*Ay_junc;
			//printf("111: junc %d, qx[%d]=%f, qy[%d]=%f\n", index, index, qxn_j_device, index, qyn_j_device);
		}
		else
		{
			temp=dt*g*J_manning*J_manning*pow(h_j_device, (-4.0/3.0))*sqrt(pow((mx/h_j_device),2)+pow((my/h_j_device),2));
			if (temp<constVariables.kesi)
			{
				qxn_j_device=mx;
				qyn_j_device=my;
				//printf("222: junc=%d, temp=%f, mx=%f, my=%f\n", index, temp, mx, my);
				//printf("222: junc %d, qx[%d]=%f, qy[%d]=%f\n", index, index, qxn_j_device, index, qyn_j_device);
			}
			else
			{
				qxn_j_device=(mx-mx*sqrt(1.0+4.0*temp))/(-2.0*temp);
				qyn_j_device=(my-my*sqrt(1.0+4.0*temp))/(-2.0*temp);
				//printf("333: junc=%d, temp=%f, mx=%f, my=%f\n", index, temp, mx, my);
				//printf("333: junc %d, qx[%d]=%f, qy[%d]=%f\n", index, index, qxn_j_device, index, qyn_j_device);
			}
		}


	}
	
	// mx=qx_j_device+dt*Ax_junc;
	// my=qy_j_device+dt*Ay_junc;

	// numerator_x=(mx-mx*sqrt(1+4.0*dt*g*pow(J_manning,2)*pow(hn_j_device,-7.0/3.0)*sqrt(mx*mx+my*my)))*pow(hn_j_device,7.0/3.0);
	// denominator_x=10e-10-2*dt*g*sqrt(mx*mx+my*my);

	// qxn_j_device=numerator_x/denominator_x;

	// numerator_y=(my-my*sqrt(1+4.0*dt*g*pow(J_manning,2)*pow(hn_j_device,-7.0/3.0)*sqrt(mx*mx+my*my)))*pow(hn_j_device,7.0/3.0);
	// denominator_y=denominator_x;

	// qyn_j_device=numerator_y/denominator_y;

	//printf('')	
}




//----------------------------------------------------------------------------------------------------------------------
__global__ void junction_updateKernal(int num_Junc, double dT,  double *h_j_device, double *hn_j_device, double *qx_j_device, 
	                            double *qxn_j_device, double *qy_j_device, double *qyn_j_device, 
								double *t_junc_device, double *junc_Area_device, Const_variables constVariables){

	/*Input_information inputKeyPtr;
	Pipe_variables pipeVariables;
	Junc_variables juncVariables;
	junction_updateKernal(inputKeyPtr.num_Junc, juncVariables.h_j_device,
					juncVariables.hn_j_device, juncVariables.qx_j_device,
					juncVariables.qxn_j_device, juncVariables.qy_j_device,
					juncVariables.qyn_j_device, juncVariables.t_junc_device,
					juncVariables.junc_Area_device, constVariables);*/

	int tid= blockDim.x * blockIdx.x + threadIdx.x;
	double ux, uy, c_sound;
	double dt1, dt2, dl;

	double kesi=constVariables.kesi;
	double g=constVariables.g;

	

	while (tid< num_Junc)
	{
		h_j_device[tid]=hn_j_device[tid];
	    qx_j_device[tid]=qxn_j_device[tid];
		qy_j_device[tid]=qyn_j_device[tid];

		// calculate dt for junction cells
		/*if (h_j_device[tid]<kesi)
		{
			u_j=0;
		}
		else
		{
			u_j=sqrt(pow(qx_j_device[tid]/h_j_device[tid],2)+pow(qy_j_device[tid]/h_j_device[tid],2));
		}
		t_junc_device[tid]=sqrt(junc_Area_device[tid])/(u_j+sqrt(g*h_j_device[tid]));*/
		//printf("tid= %d, t_junc_device= %f\n",tid,t_junc_device[tid] );
		if (h_j_device[tid]>kesi)
		{
			dl=sqrt(junc_Area_device[tid]);
			ux=qx_j_device[tid]/h_j_device[tid];
			uy=qy_j_device[tid]/h_j_device[tid];
			c_sound=sqrt(h_j_device[tid]*g);
			dt1=dl/(c_sound + fabs(ux));
			dt2=dl/(c_sound + fabs(uy));
			t_junc_device[tid]=constVariables.MIN(dt1,dt2);

		}
		else
		{
			t_junc_device[tid]=0.01;
		}
	//	printf("tid= %d, t_junc_device[tid]=%f \n",tid,t_junc_device[tid] );



		tid += blockDim.x * gridDim.x;

	}
}


//----------------------------------------------------------------------------------------------------------------------
void pipe_boundary(Input_information inputKeyPtr,
	               Pipe_variables &pipeVariables,
				   Junc_variables &juncVariables,
				   Outf_variables &outfVariables,
				   int blocksPerGrid, 
				   int threadsPerBlock){

	 pipe_boundaryKernal <<<blocksPerGrid,threadsPerBlock>>> (inputKeyPtr.num_Pipe,
		                                                        pipeVariables.h_p_device,
																juncVariables.h_j_device,
																juncVariables.qx_j_device,
																juncVariables.qy_j_device,
																pipeVariables.P_Obound_device,
																pipeVariables.upsPJ_device,
																pipeVariables.downsPJ_device,
																pipeVariables.downsPO_device,
																pipeVariables.upNorm1_device,
																pipeVariables.upNorm2_device,
																pipeVariables.downNorm1_device,
																pipeVariables.downNorm2_device,
																pipeVariables.p_ups_h_device,
																pipeVariables.p_ups_Q_device,
																pipeVariables.p_downs_h_device,
													  			pipeVariables.p_downs_Q_device);
}



//----------------------------------------------------------------------------------------------------------------------
void pipe_calculation(double T, 
	                  double dT,
					  double *surfWaterDepth,
					  Vector *surfQ_device,
					  Input_information inputKeyPtr,
					  Pipe_variables &pipeVariables,
					  Junc_variables &juncVariables,
					  Outf_variables &outfVariables,
					  Drain_outfBounds &drainOutfBounds,
					  Const_variables constVariables,
					  int blocksPerGrid, 
				      int threadsPerBlock){




			  pipe_calculationKernal <<< blocksPerGrid,threadsPerBlock>>>(pipeVariables.totalCellNum, inputKeyPtr.num_Pipe,
								                T, dT, inputKeyPtr.a_speed, 
												pipeVariables.upsBCell_device, pipeVariables.downsBCell_device,
												pipeVariables.upsPJ_device, pipeVariables.downsPJ_device,
												pipeVariables.downsPO_device, pipeVariables.P_Obound_device,
												pipeVariables.P_diameter_device, pipeVariables.dx_cell_device,
												juncVariables.zb_j_device, pipeVariables.P_manning_device, 
												pipeVariables.zb_p_device, pipeVariables.h_p_device, 
												pipeVariables.Q_p_device, pipeVariables.sita_p_device, 
												pipeVariables.A_p_device, pipeVariables.An_p_device, 
												pipeVariables.Qn_p_device, pipeVariables.pipe_flux_ups1,
												pipeVariables.pipe_flux_ups2, pipeVariables.pipe_flux_downs1,
												pipeVariables.pipe_flux_downs2, pipeVariables.p_ups_h_device,
												pipeVariables.p_ups_Q_device, pipeVariables.p_downs_h_device,
												pipeVariables.p_downs_Q_device, 
												pipeVariables.downNorm1_device,
                                                pipeVariables.downNorm2_device,
												drainOutfBounds.interval_ptr_device,  
												drainOutfBounds.T_series_device,
												drainOutfBounds.h_series_device,
												drainOutfBounds.Q_series_device,
												drainOutfBounds.T_start_device,
												drainOutfBounds.h_start_device,
												drainOutfBounds.Q_start_device,
												drainOutfBounds.intervalLen_device,
												outfVariables.surfGrid_device,
												surfWaterDepth, 
												surfQ_device,
												constVariables);


}			

//----------------------------------------------------------------------------------------------------------------------
void pipe_update(Input_information inputKeyPtr,
				Pipe_variables &pipeVariables,
				Const_variables constVariables,
				double dT,
				int blocksPerGrid, 
				int threadsPerBlock){


    pipe_updateKernal <<< blocksPerGrid,threadsPerBlock>>> (pipeVariables.totalCellNum, inputKeyPtr.num_Pipe,
														  inputKeyPtr.a_speed, dT, pipeVariables.P_diameter_device,
														  pipeVariables.A_p_device, pipeVariables.An_p_device,
														  pipeVariables.Q_p_device, pipeVariables.Qn_p_device,
														  pipeVariables.h_p_device, pipeVariables.hn_p_device,
														  pipeVariables.sita_p_device, pipeVariables.sitan_p_device,
														  pipeVariables.downsBCell_device, pipeVariables.pipe_Ap_device,
														  pipeVariables.dx_cell_device, pipeVariables.t_pipe_device,
														  pipeVariables.pipe_flux_downs1, constVariables);


}



//----------------------------------------------------------------------------------------------------------------------
void junction_calculation(double dT,
						  Input_information inputKeyPtr,
						  Pipe_variables &pipeVariables,
						  Junc_variables &juncVariables,
						  Const_variables constVariables,
						  int blocksPerGrid, 
						  int threadsPerBlock){


	// junction_calculationKernal  <<<blocksPerGrid,threadsPerBlock>>> (inputKeyPtr.num_Junc, juncVariables.width_upsJP,
	// 																	juncVariables.width_downsJP, dT,
	// 																	juncVariables.junc_MN_device, pipeVariables.pipe_Ap_device,
	// 																	pipeVariables.P_diameter_device, juncVariables.junc_Area_device,
	// 																	pipeVariables.upNorm1_device, pipeVariables.upNorm2_device,
	// 																	pipeVariables.downNorm1_device, pipeVariables.downNorm2_device, 
	// 																	juncVariables.upsJP_device, juncVariables.downsJP_device,
	// 																	pipeVariables.pipe_flux_ups1, pipeVariables.pipe_flux_ups2,
	// 																	pipeVariables.pipe_flux_downs1, pipeVariables.pipe_flux_downs2,
	// 																	juncVariables.h_j_device, juncVariables.hn_j_device,
	// 																	juncVariables.qx_j_device, juncVariables.qxn_j_device,
	// 																	juncVariables.qy_j_device, juncVariables.qyn_j_device,
	// 																	juncVariables.surfGrid_device, juncVariables.Qexchange_device,
	// 																	constVariables);

	junction_calculationKernal  <<<blocksPerGrid,threadsPerBlock>>> (inputKeyPtr.num_Junc, juncVariables.width_upsJP,
																		juncVariables.width_downsJP, dT,
																		juncVariables.junc_MN_device, pipeVariables.pipe_Ap_device,
																		pipeVariables.P_diameter_device, juncVariables.junc_Area_device,
																		pipeVariables.upNorm1_device, pipeVariables.upNorm2_device,
																		pipeVariables.downNorm1_device, pipeVariables.downNorm2_device, 
																		juncVariables.upsJP_device, juncVariables.downsJP_device,
																		pipeVariables.pipe_flux_ups1, pipeVariables.pipe_flux_ups2,
																		pipeVariables.pipe_flux_downs1, pipeVariables.pipe_flux_downs2,
																		juncVariables.h_j_device, juncVariables.hn_j_device,
																		juncVariables.qx_j_device, juncVariables.qxn_j_device,
																		juncVariables.qy_j_device, juncVariables.qyn_j_device,
																		juncVariables.surfGrid_device, juncVariables.Qexchange_device,
																		juncVariables.Ax_junc_device, juncVariables.Ay_junc_device, 
																		juncVariables.hydroPX_device, juncVariables.hydroPY_device,
																		juncVariables.pipe_F1_device, juncVariables.pipe_F2_device,
																		juncVariables.pipe_F3_device,
																		constVariables);

}


//----------------------------------------------------------------------------------------------------------------------
void junction_update(Input_information inputKeyPtr,
					 Junc_variables &juncVariables,
					 Const_variables constVariables,
					 double dT, 
					 int blocksPerGrid, 
					 int threadsPerBlock){

	junction_updateKernal <<<blocksPerGrid,threadsPerBlock>>> (inputKeyPtr.num_Junc, dT, juncVariables.h_j_device,
															 juncVariables.hn_j_device, juncVariables.qx_j_device,
															 juncVariables.qxn_j_device, juncVariables.qy_j_device,
															 juncVariables.qyn_j_device, juncVariables.t_junc_device,
															 juncVariables.junc_Area_device, constVariables);


}


//----------------------------------------------------------------------------------------------------------------------
void adaptiveDrainageTime(double T, double Tout,
	                 Input_information inputKeyPtr,
					 Junc_variables &juncVariables,
					 Pipe_variables &pipeVariables,
					 Const_variables constVariables,
					 double &dT_adjust){

	double dT_pipeAdjust, dT_juncAdjust;
	
		
	adaptivePipeTime(dT_pipeAdjust, pipeVariables, constVariables);	
			
	//adaptiveJuncTime(dT_juncAdjust, inputKeyPtr,
	//				 juncVariables, constVariables);
	

	//dT_adjust=constVariables.fmin(dT_pipeAdjust, dT_juncAdjust)*constVariables.CFL;
	dT_adjust=dT_pipeAdjust*constVariables.CFL;

	// if (T+dT_adjust>Tout)
	// {
	// 	dT_adjust=Tout-T;
	// }
	printf("next time step adjusted by CFL: \n");	
	//printf("dT_pipeAdjust= %f, dT_juncAdjust= %f, dT_adjust=%f\n",dT_pipeAdjust, dT_juncAdjust, dT_adjust);
	printf("dT_pipeAdjust= %f, dT_adjust=%f\n",dT_pipeAdjust, dT_adjust);
}


//----------------------------------------------------------------------------------------------------------------------
void adaptivePipeTime( double &dT_pipeAdjust,
					   Pipe_variables &pipeVariables,
	                   Const_variables constVariables){
	
	//cudaMemcpy(pipeVariables.t_pipe, pipeVariables.t_pipe_device, pipeVariables.totalCellNum*sizeof(double), cudaMemcpyDeviceToHost);
	int i=0;
	double temp=3e35;
	while (i< pipeVariables.totalCellNum)
	{
		
		temp=constVariables.fmin(temp, pipeVariables.t_pipe[i]);
		//printf("tid= %d, temp= %f, dt_pipe_device=%f\n", i, temp,pipeVariables.t_pipe[i]);
		i=i+1;
	}
	dT_pipeAdjust=temp;
	//printf("dT_pipeAdjust=%f\n", dT_pipeAdjust);
}

//----------------------------------------------------------------------------------------------------------------------
void adaptiveJuncTime(double &dT_juncAdjust,
	                       Input_information inputKeyPtr,
					       Junc_variables &juncVariables,
	                       Const_variables constVariables){

    //cudaMemcpy(juncVariables.t_junc, juncVariables.t_junc_device,  inputKeyPtr.num_Junc*sizeof(double), cudaMemcpyDeviceToHost);
	double temp=3e35;
	int i=0;
	while (i<inputKeyPtr.num_Junc)
	{
		
		temp=constVariables.fmin(temp, juncVariables.t_junc[i]);
	//	printf("tid= %d, temp= %f, dt_junc_device=%f\n", i, temp, juncVariables.t_junc[i]);
		i=i+1;
	}

	dT_juncAdjust=temp;
	//printf("dT_juncAdjust=%f\n", dT_juncAdjust);

}


// //----------------------------------------------------------------------------------------------------------------------
// void outputResults(double timeInterval, double &outputTime, 
// 	               double T, double &dT, double Tout,
// 				   Input_information inputKeyPtr,
// 	               Pipe_attributes pipeAttributes, 
// 				   Pipe_variables &pipeVariables,
// 				   Junc_attributes juncAttributes,
// 				   Junc_variables &juncVariables,
// 				   Const_variables constVariables){

// 	// in this function, dT is adjusted to reach each output time
// 	// here the dT represents the time step in the next time level
// 	double mass;
// 	mass=pipeVariables.mass+juncVariables.mass;
// 	if (T+dT>outputTime)
// 	{
// 		dT=outputTime-T;
// 	}
	
// 	if (T==outputTime)
// 	{
// 		//cudaMemcpy(pipeVariables.h_p_out, pipeVariables.h_p_device,  pipeVariables.totalCellNum*sizeof(double), cudaMemcpyDeviceToHost);
// 	    //cudaMemcpy(juncVariables.h_j_out, juncVariables.h_j_device,  inputKeyPtr.num_Junc*sizeof(double), cudaMemcpyDeviceToHost);

// 		writePipeResults(T, inputKeyPtr, pipeAttributes, pipeVariables);
// 		writeJuncResults(T, inputKeyPtr, juncAttributes, juncVariables);

// 		outputTime=outputTime+timeInterval;
// 		dT=constVariables.dT_Default;
// 	}

// 	writeTimeDataFile("Time",T, dT, Tout);
// 	writeMassDataFile("Mass",T, dT, Tout, mass);
// 	printf("next time step adjusted by report time: %f\n", dT);
	
// }
//----------------------------------------------------------------------------------------------------------------------
void writeMassDataFile(const char* filename, double T, double dT, double Tout, 
	                   double mass){
	
	FILE * output;
	const char* directory = "output/";
	std::string name = std::string(directory) + "Drainage_"+ std::string(filename) + ".txt";
	output = fopen(name.c_str(), "a");
	
	fprintf(output, "%f            %f \n", T, mass);
	
	if (T==Tout)
		fprintf(output, "*********************************************************\n", T, dT);
	
	fclose(output);

}

//----------------------------------------------------------------------------------------------------------------------
void writeTimeDataFile(const char* filename, double T, double dT, double Tout){
						  
	
	FILE * output;
	const char* directory = "output/";
	std::string name = std::string(directory) + "Drainage_"+ std::string(filename) + ".txt";
	output = fopen(name.c_str(), "a");
	
	fprintf(output, "%f            %f \n", T, dT);
	
	if (T==Tout)
		fprintf(output, "*********************************************************\n", T, dT);
	
	fclose(output);

}

//----------------------------------------------------------------------------------------------------------------------

void writePipeResults(double T, 
	                  Input_information inputKeyPtr,
	                  Pipe_attributes pipeAttributes, 
	                  Pipe_variables &pipeVariables){

	//const char* directory = "output/";
	//int i, j, preSumIndex=0;
	//FILE * output;

	///*Pipe_attributes pipeAttributes;
	// Pipe_variables pipeVariables;
	// Input_information inputKeyPtr;*/

 //   std::ostringstream out_time;
	//out_time << t;

	//// write the pipe water depth file
 //   std::string name = std::string(directory) + "DrainagePipe_"+ std::string(filename) + "_" + out_time.str() + ".txt";
 //   output = fopen(name.c_str(), "w");
	//fprintf(output, "Pipe index         Cell Number     Water depth in each cell \n");
	//
	//for (i=0;i<inputKeyPtr.num_Pipe;i++)
	//{
	//	fprintf(output, "  %d   %d  ", pipeAttributes.P_index[i], pipeVariables.p_cellNum[i]);
	//	for (j=0;j<pipeVariables.p_cellNum[i];j++)
	//	{
	//		fprintf(output, "    %f   ", pipeVariables.h_p_out[preSumIndex+j]);
	//	}
	//	fprintf(output, "\n");
	//	preSumIndex=preSumIndex+pipeVariables.p_cellNum[i];

	//}

	pipeWaterDepthFile("h", T, inputKeyPtr, pipeAttributes, pipeVariables);
	pipeFlowRateFile("Q", T, inputKeyPtr, pipeAttributes, pipeVariables);
	//pipeCrossSectionFile("A", T, inputKeyPtr, pipeAttributes, pipeVariables);



    //fclose(output);

}

//----------------------------------------------------------------------------------------------------------------------
void pipeWaterDepthFile(const char* filename, double t, 
	                  Input_information inputKeyPtr,
	                  Pipe_attributes pipeAttributes, 
	                  Pipe_variables &pipeVariables){

	const char* directory = "output/";
	int i, j, preSumIndex=0;
	FILE * output;

	/*Pipe_attributes pipeAttributes;
	 Pipe_variables pipeVariables;
	 Input_information inputKeyPtr;*/

    std::ostringstream out_time;
	out_time << t;

	// write the pipe water depth file
    std::string name = std::string(directory) + "DrainagePipe_"+ std::string(filename) + "_" + out_time.str() + ".txt";
    output = fopen(name.c_str(), "w");
	fprintf(output, "Pipe index         Cell Number     Water depth in each cell \n");
	
	for (i=0;i<inputKeyPtr.num_Pipe;i++)
	{
		fprintf(output, "  %d   %d  ", pipeAttributes.P_index[i], pipeVariables.p_cellNum[i]);
		for (j=0;j<pipeVariables.p_cellNum[i];j++)
		{
			fprintf(output, "    %f   ", pipeVariables.h_p_out[preSumIndex+j]);
		}
		fprintf(output, "\n");
		preSumIndex=preSumIndex+pipeVariables.p_cellNum[i];

	}
	fclose(output);
}


//----------------------------------------------------------------------------------------------------------------------
void pipeCrossSectionFile(const char* filename, double t, 
	                  Input_information inputKeyPtr,
	                  Pipe_attributes pipeAttributes, 
	                  Pipe_variables &pipeVariables){

	const char* directory = "output/";
	int i, j, preSumIndex=0;
	FILE * output;

	/*Pipe_attributes pipeAttributes;
	 Pipe_variables pipeVariables;
	 Input_information inputKeyPtr;*/

    std::ostringstream out_time;
	out_time << t;

	// write the pipe water depth file
    std::string name = std::string(directory) + "DrainagePipe_"+ std::string(filename) + "_" + out_time.str() + ".txt";
    output = fopen(name.c_str(), "w");
	fprintf(output, "Pipe index         Cell Number     Water depth in each cell \n");
	
	for (i=0;i<inputKeyPtr.num_Pipe;i++)
	{
		/*fprintf(output, " %d  %f   %f  %f  %f\n",i,  pipeVariables.pipe_flux_ups1_host[i], pipeVariables.pipe_flux_ups2_host[i], 
			                                    pipeVariables.pipe_flux_downs1_host[i], pipeVariables.pipe_flux_downs2_host[i]);*/
		/*for (j=0;j<pipeVariables.p_cellNum[i];j++)
		{
			fprintf(output, "    %f   ", pipeVariables.A_p_out[preSumIndex+j]);
		}
		fprintf(output, "\n");
		preSumIndex=preSumIndex+pipeVariables.p_cellNum[i];*/

	}
	fclose(output);
}

//----------------------------------------------------------------------------------------------------------------------
void pipeFlowRateFile(const char* filename, double t, 
	                  Input_information inputKeyPtr,
	                  Pipe_attributes pipeAttributes, 
	                  Pipe_variables &pipeVariables){

	//cudaMemcpy(pipeVariables.Q_p_out, pipeVariables.Q_p_device, pipeVariables.totalCellNum*sizeof(double), cudaMemcpyDeviceToHost);


	const char* directory = "output/";
	int i, j, preSumIndex=0;
	FILE * output;

	/*Pipe_attributes pipeAttributes;
	 Pipe_variables pipeVariables;
	 Input_information inputKeyPtr;*/

    std::ostringstream out_time;
	out_time << t;

	// write the pipe water depth file
    std::string name = std::string(directory) + "DrainagePipe_"+ std::string(filename) + "_" + out_time.str() + ".txt";
    output = fopen(name.c_str(), "w");
	fprintf(output, "Pipe index         Cell Number     Flow Rate in each cell \n");
	
	for (i=0;i<inputKeyPtr.num_Pipe;i++)
	{
		fprintf(output, "  %d   %d  ", pipeAttributes.P_index[i], pipeVariables.p_cellNum[i]);
		for (j=0;j<pipeVariables.p_cellNum[i];j++)
		{
			fprintf(output, "    %f   ", pipeVariables.Q_p_out[preSumIndex+j]);
		}
		fprintf(output, "\n");
		preSumIndex=preSumIndex+pipeVariables.p_cellNum[i];

	}
	fclose(output);

}


//----------------------------------------------------------------------------------------------------------------------
void writeJuncResults(double T, 
	                  Input_information inputKeyPtr,
	                  Junc_attributes juncAttributes, 
	                  Junc_variables &juncVariables){

	
	juncWaterDepthFile("h", T, inputKeyPtr, juncAttributes, juncVariables);
	juncFlowRateFile("Q", T, inputKeyPtr, juncAttributes, juncVariables);
}

//----------------------------------------------------------------------------------------------------------------------
void juncWaterDepthFile(const char* filename, double t, 
	                  Input_information inputKeyPtr,
	                  Junc_attributes juncAttributes, 
	                  Junc_variables &juncVariables){

	const char* directory = "output/";
	int i;
	FILE * output;

	/*Pipe_attributes pipeAttributes;
	 Pipe_variables pipeVariables;
	 Input_information inputKeyPtr;*/

    std::ostringstream out_time;
	out_time << t;
    std::string name = std::string(directory) + "DrainageJunc_"+ std::string(filename) + "_" + out_time.str() + ".txt";
    output = fopen(name.c_str(), "w");
	fprintf(output, "Junction index        Water depth in the junction \n");
	
	for (i=0;i<inputKeyPtr.num_Junc;i++)
	{
		fprintf(output, "  %d   %f  \n", juncAttributes.J_index[i], juncVariables.h_j_out[i]);
	}

    fclose(output);
}


//----------------------------------------------------------------------------------------------------------------------
void juncFlowRateFile(const char* filename, double t, 
	                  Input_information inputKeyPtr,
	                  Junc_attributes juncAttributes, 
	                  Junc_variables &juncVariables){

	cudaMemcpy(juncVariables.qx_j_out, juncVariables.qx_j_device, inputKeyPtr.num_Junc*sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(juncVariables.qy_j_out, juncVariables.qy_j_device, inputKeyPtr.num_Junc*sizeof(double), cudaMemcpyDeviceToHost);

    const char* directory = "output/";
	int i;
	FILE * output;

	/*Pipe_attributes pipeAttributes;
	 Pipe_variables pipeVariables;
	 Input_information inputKeyPtr;*/

    std::ostringstream out_time;
	out_time << t;
    std::string name = std::string(directory) + "DrainageJunc_"+ std::string(filename) + "_" + out_time.str() + ".txt";
    output = fopen(name.c_str(), "w");
	fprintf(output, "Junction index       X- Flow Rate    Y- Flow Rate \n");
	
	for (i=0;i<inputKeyPtr.num_Junc;i++)
	{
		fprintf(output, "  %d   %f  %f\n", juncAttributes.J_index[i], juncVariables.qx_j_out[i], juncVariables.qy_j_out[i]);
	}

    fclose(output);
}



}