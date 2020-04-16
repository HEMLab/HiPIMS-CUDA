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




#ifndef DRAINAGE_FUNCTIONS_H
#define DRAINAGE_FUNCTIONS_H
//#include "cuda_mapped_field.h"
//#include "Scalar.h"
#include "Vector.h"
//#include "Tensor.h"
//#include "Flag.h"

namespace GC{

	//--------host function list---------------------------------------------------------------------------------------------------------

char UCHAR(char x);
void readInputFile(char* Keywords[], Input_information &inputKeyPtr, 
	               Junc_attributes &juncAttributes, 
				   Outf_attributes &outfAttributes,
				   Pipe_attributes &pipeAttributes,
				   Pipe_variables &pipeVariables);
void introductionFile(char* Keywords[], Input_information &inputKeyPtr);
void firstOpenIntroductionFile(char* Keywords[], Input_information &inputKeyPtr);
int  match(char *str, char *substr);
int  findMatch(char *s, char *keyword[]);
void secondOpenIntroductionFile(Input_information &inputKeyPtr);
void hostMemoryCreationForObjects(Input_information &inputKeyPtr, 
		                         Junc_attributes &juncAttributes,
								 Outf_attributes &outfAttributes,
				                 Pipe_attributes &pipeAttributes);
void junctionFile(Input_information &inputKeyPtr, 
		          Junc_attributes &juncAttributes);
void outfallFile(Input_information &inputKeyPtr,
	             Outf_attributes &outfAttributes);
void pipeFile(Input_information &inputKeyPtr, 
		      Pipe_attributes &pipeAttributes);
void verticeFile(Input_information &inputKeyPtr, 
		         Pipe_attributes &pipeAttributes);
void pipeCellsDataFile(Input_information inputKeyPtr, 
					   Pipe_attributes pipeAttributes,
					   Pipe_variables &pipeVariables );
void NodeGridFile(Input_information inputKeyPtr,
				  Junc_attributes &juncAttributes,
				  Outf_attributes &outfAttributes);
void getPipeCellNumber(Input_information inputKeyPtr, 
				   Pipe_attributes pipeAttributes,
				   Pipe_variables &pipeVariables);
void getPipeCellsData(Input_information inputKeyPtr, 
					   Pipe_attributes pipeAttributes,
					   Pipe_variables &pipeVariables,
					   int totalCellNum_ptr, int cellNum_ptr, int cellWidth_ptr);
void allocMemory(Input_information inputKeyPtr,
	              Pipe_variables &pipeVariables,
				  Junc_variables &juncVariables,
				  Outf_variables &outfVariables);
void preProcessing(Input_information inputKeyPtr,
	               Pipe_attributes pipeAttributes,
			   	   Pipe_variables &pipeVariables,
				   Junc_attributes juncAttributes,
				   Junc_variables &juncVariables,
				   Outf_attributes outfAttributes,
				   Outf_variables &outfVariables,
				   Const_variables constVariables);
void pipeCellWidth(Input_information inputKeyPtr,
	               Pipe_attributes pipeAttributes,
			   	   Pipe_variables &pipeVariables);
void pipeToJuncMatrix(Input_information inputKeyPtr,
	                              Pipe_attributes pipeAttributes,
								  Pipe_variables &pipeVariables);
void juncToPipeMatrix(Input_information inputKeyPtr,
			   					  Pipe_variables &pipeVariables,
								  Junc_variables &juncVariables);
void outfToPipeMatrix(Input_information inputKeyPtr,
			   			Pipe_variables &pipeVariables,
						Outf_variables &outfVariables);
void normVec(Input_information inputKeyPtr,
						Pipe_attributes pipeAttributes,
			   			Pipe_variables &pipeVariables,
						Junc_attributes juncAttributes,
						Junc_variables &juncVariables);
void pipeBedElevation(Input_information inputKeyPtr,
					Pipe_attributes pipeAttributes,
			   		Pipe_variables &pipeVariables);
void boundaryPipeCellIndex(Input_information inputKeyPtr,
			   			   Pipe_variables &pipeVariables);
void areaCalculation(Input_information inputKeyPtr,
	                  Const_variables constVariables,
	                  Pipe_attributes pipeAttributes,
			   		  Pipe_variables &pipeVariables,
					  Junc_attributes juncAttributes,
					  Junc_variables &juncVariables);

void initFlow(Input_information inputKeyPtr,
	           Pipe_attributes pipeAttributes,
			   Pipe_variables &pipeVariables,
			   Junc_attributes juncAttributes,
			   Junc_variables &juncVariables);
void copyMemory(Input_information inputKeyPtr,
	             Pipe_attributes pipeAttributes,
	             Pipe_variables &pipeVariables,
				 Junc_attributes juncAttributes,
	             Junc_variables &juncVariables,
				 Outf_attributes outfAttributes,
				 Outf_variables &outfVariables);
void releaseMemory(Pipe_attributes &pipeAttributes,
					Pipe_variables &pipeVariables,
					Junc_attributes &juncAttributes,
					Junc_variables &juncVariables,
					Outf_attributes &outfAttributes,
					Outf_variables &ourfVariables,
					Drain_outfBounds &drainOutfBounds, 
					Drain_gauges &drainGauges);


void pipe_boundary(Input_information inputKeyPtr,
	               Pipe_variables &pipeVariables,
				   Junc_variables &juncVariables,
				   Outf_variables &outfVariables,
				   int blocksPerGrid, 
				   int threadsPerBlock);
void pipe_calculation(double T, 
	                  double dT,
					  double *surfWaterDepth,
					  Vector *surfQx_device,
					  Input_information inputKeyPtr,
					  Pipe_variables &pipeVariables,
					  Junc_variables &juncVariables,
					  Outf_variables &outfVariables,
					  Drain_outfBounds &drainOutfBounds,
					  Const_variables constVariables,
					  int blocksPerGrid, 
				      int threadsPerBlock);
void pipe_update(Input_information inputKeyPtr,
				Pipe_variables &pipeVariables,
				Const_variables constVariables,
				double dT,
				int blocksPerGrid, 
				int threadsPerBlock);
void junction_calculation(double dT,
						Input_information inputKeyPtr,
						Pipe_variables &pipeVariables,
						Junc_variables &juncVariables,
						Const_variables constVariables,
						int blocksPerGrid, 
						int threadsPerBlock);
void junction_update(Input_information inputKeyPtr,
					 Junc_variables &juncVariables,
					 Const_variables constVariables,
					 double dT, 
					 int blocksPerGrid, 
					 int threadsPerBlock);
void adaptiveDrainageTime(double T, double Tout,
	                 Input_information inputKeyPtr,
					 Junc_variables &juncVariables,
					 Pipe_variables &pipeVariables,
					 Const_variables constVariables,
					 double &dT_adjust);
void adaptivePipeTime( double &dT_pipeAdjust,
					   Pipe_variables &pipeVariables,
	                   Const_variables constVariables);
void adaptiveJuncTime(double &dT_juncAdjust,
	                       Input_information inputKeyPtr,
					       Junc_variables &juncVariables,
	                       Const_variables constVariables);
void outputResults(double timeInterval, double &outputTime, 
	               double T, double &dT, double Tout,
				   Input_information inputKeyPtr,
	               Pipe_attributes pipeAttributes, 
				   Pipe_variables &pipeVariables,
				   Junc_attributes juncAttributes,
				   Junc_variables &juncVariables,
				   Const_variables constVariables);
void writePipeResults(double T, 
	                  Input_information inputKeyPtr,
	                  Pipe_attributes pipeAttributes, 
	                  Pipe_variables &pipeVariables);
void pipeCrossSectionFile(const char* filename, double t, 
	                  Input_information inputKeyPtr,
	                  Pipe_attributes pipeAttributes, 
	                  Pipe_variables &pipeVariables);
void pipeFlowRateFile(const char* filename, double t, 
	                  Input_information inputKeyPtr,
	                  Pipe_attributes pipeAttributes, 
	                  Pipe_variables &pipeVariables);
void pipeWaterDepthFile(const char* filename, double t, 
	                  Input_information inputKeyPtr,
	                  Pipe_attributes pipeAttributes, 
	                  Pipe_variables &pipeVariables);
void writeJuncResults(double T, 
	                  Input_information inputKeyPtr,
	                  Junc_attributes juncAttributes, 
	                  Junc_variables &juncVariables);
void juncWaterDepthFile(const char* filename, double t, 
	                  Input_information inputKeyPtr,
	                  Junc_attributes juncAttributes, 
	                  Junc_variables &juncVariables);
void juncFlowRateFile(const char* filename, double t, 
	                  Input_information inputKeyPtr,
	                  Junc_attributes juncAttributes, 
	                  Junc_variables &juncVariables);
void writeTimeDataFile(const char* filename, double T, double dT, double Tout);
void writeMassDataFile(const char* filename, double T, double dT, double Tout, 
	                   double mass);

void outfallBoundaryData(Input_information inputKeyPtr,  
	                    Outf_attributes &outfAttributes,
						Drain_outfBounds &drainOutfBounds);
void seriesLengthData(Drain_outfBounds &drainOutfBounds);
void outfallBoundFile(Input_information inputKeyPtr, Drain_outfBounds &drainOutfBounds);
void TSeriesBoundFile(Drain_outfBounds &drainOutfBounds);
void hSeriesBoundFile(Drain_outfBounds &drainOutfBounds);
void QSeriesBoundFile(Drain_outfBounds &drainOutfBounds);
void gaugeSettings(Drain_gauges &drainGauges);
void obtainGaugeNumber(Drain_gauges &drainGauges);
void obtainGaugePositions(Drain_gauges &drainGauges);
void gaugeOutputData(double T, Drain_gauges &drainGauges, Pipe_variables pipeVariables, Junc_variables juncVariables);
void pipeGaugeOutput(double T, Drain_gauges &drainGauges, Pipe_variables pipeVariables);
void juncGaugeOutput(double T, Drain_gauges &drainGauges, Junc_variables juncVariables);

void timeAdjustForHipims(double dt_s, double dt_c, double &dt_pre, double &dt_new);
void surfDrainQ_calculation(Input_information inputKeyPtr,Junc_variables &juncVariables, 
	                        Const_variables constVariables, double *zb_surface, double *hs,
							double *Qsuface, int blocksPerGrid, int threadsPerBlock);
	


}

#endif