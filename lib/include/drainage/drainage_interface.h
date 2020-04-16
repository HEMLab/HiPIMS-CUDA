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




#ifndef DRAINAGE_INTERFACE_H
#define DRAINAGE_INTERFACE_H
#include "cuda_mapped_field.h"
#include "Scalar.h"
#include "Vector.h"
#include "Tensor.h"
#include "Flag.h"


namespace GC {



	void surfDrainQ_calculation(double dT, Input_information inputKeyPtr, Junc_variables &juncVariables,
		Const_variables constVariables, cuFvMappedField<Scalar, on_cell>& z, cuFvMappedField<Scalar, on_cell>& h, cuFvMappedField<Scalar, on_cell>& Q_surface,
		int blocksPerGrid, int threadsPerBlock);
	void surfH_limitator(cuFvMappedField<Scalar, on_cell>& h, int blocksPerGrid, int threadsPerBlock);


}

#endif