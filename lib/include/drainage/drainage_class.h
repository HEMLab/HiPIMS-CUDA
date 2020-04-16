
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



#ifndef DRAINAGE_CLASS_H
#define DRAINAGE_CLASS_H


namespace GC{
	//--------class definition---------------------------------------------------------------------------------------------------------

	class Input_information
	{
	public:
		int p_aspeed;
		int p_juncNum;
		int p_pipeNum;
		int p_outfallNum;
	 /*   int p_junctions;
		int p_outfalls;
		int p_pipes;
		int p_xsections;
		int p_vertices;    
		int p_symbols; */
	
		Input_information()
		{
			p_aspeed=-1;
			p_juncNum=-1;
			p_pipeNum=-1;
			p_outfallNum=-1;
			/*p_junctions=-1;
			p_outfalls=-1;
			p_pipes=-1;
			p_xsections=-1;
			p_vertices=-1;
			p_symbols=-1;*/
		
	

		}
	public:
		double a_speed;
		int num_Junc;
		int num_Pipe;
		int num_Outfall;

	};


	//******************************************************************************************************************* 
	class Pipe_attributes
	{
	public:
		// variables of pipe attributes on CPU
		int *P_index;
		int *P_inNode;
		int *P_outNode;
		int *P_outfall;
		int *P_outfBound;
		double *P_length;
		double *P_dx;
		double *P_manning;
		double *P_diameter;
		double *P_inletHeight;
		double *P_outletHeight;
		double *P_coorUpx;
		double *P_coorUpy;
		double *P_coorCenx;
		double *P_coorCeny;
		double *P_coorDownx;
		double *P_coorDowny;
		double *P_coorUpAx;
		double *P_coorUpAy;
		double *P_coorUpBx;
		double *P_coorUpBy;
		double *P_coorDownCx;
		double *P_coorDownCy;
		double *P_coorDownDx;
		double *P_coorDownDy;
		double *P_initWaterDepth;
		double *P_initFlowRate;
	
		void setHostMemory(int num_Pipe)
		{
			P_index = new int[num_Pipe];
			P_inNode = new int[num_Pipe];
			P_outNode = new int[num_Pipe];
			P_outfall= new int[num_Pipe];
			P_outfBound = new int[num_Pipe];
			P_length = new double[num_Pipe];
			P_dx = new double[num_Pipe];
			P_manning = new double[num_Pipe];
			P_diameter = new double[num_Pipe];
			P_inletHeight = new double[num_Pipe];
			P_outletHeight = new double[num_Pipe];
			P_coorUpx = new double[num_Pipe];
			P_coorUpy = new double[num_Pipe];
			P_coorCenx = new double[num_Pipe];
			P_coorCeny = new double[num_Pipe];
			P_coorDownx = new double[num_Pipe];
			P_coorDowny = new double[num_Pipe];
			P_coorUpAx = new double[num_Pipe];
			P_coorUpAy = new double[num_Pipe];
			P_coorUpBx = new double[num_Pipe];
			P_coorUpBy = new double[num_Pipe];
			P_coorDownCx = new double[num_Pipe];
			P_coorDownCy = new double[num_Pipe];
			P_coorDownDx = new double[num_Pipe];
			P_coorDownDy = new double[num_Pipe];
			P_initWaterDepth = new double[num_Pipe];
			P_initFlowRate = new double[num_Pipe];
		
		}


		void deleteHostMemory()
		{
			delete[] P_index;
			delete[] P_inNode;
			delete[] P_outNode;
			delete[] P_outfall;
			delete[] P_outfBound;
			delete[] P_length;
			delete[] P_dx;
			delete[] P_manning;
			delete[] P_diameter;
			delete[] P_inletHeight;
			delete[] P_outletHeight;
			delete[] P_coorUpx;
			delete[]P_coorUpy;
			delete[]P_coorCenx;
			delete[]P_coorCeny;
			delete[]P_coorDownx;
			delete[]P_coorDowny;
			delete[] P_coorUpAx;
			delete[] P_coorUpAy;
			delete[] P_coorUpBx;
			delete[] P_coorUpBy;
			delete[] P_coorDownCx;
			delete[] P_coorDownCy;
			delete[] P_coorDownDx;
			delete[] P_coorDownDy;
			delete[] P_initWaterDepth;
			delete[] P_initFlowRate;
		}

	
	};

	//******************************************************************************************************************* 
	class Pipe_variables
	{
	public:

		int totalCellNum;
		double mass;

		double *pipe_Ap;
		int *p_cellNum;  // the cell num of each pipe
		int *ups_PJ;
		int *downs_PJ;
		int *downs_PO;
		double *up_nor1;
		double *up_nor2;
		double *down_nor1;
		double* down_nor2; 
		int *upstreamCell;
		int *downstreamCell;
		double *dx_cell;
		double *h_p_out;
		double *A_p_out;
		double *Q_p_out;

		double *h_p;
		double *Q_p;
		double *sita_p;
		double *A_p;
		double *zb_p;

		double *zb_p_device;
		double *h_p_device;
		double *Q_p_device;
		double *sita_p_device;
		double *A_p_device;
		double *hn_p_device;
		double *Qn_p_device;
		double *An_p_device;
		double *sitan_p_device;

		double *pipe_flux_ups1;
		double *pipe_flux_ups2;
		double *pipe_flux_downs1;
		double *pipe_flux_downs2;
		double *pipe_flux_ups1_host;
		double *pipe_flux_ups2_host;
		double *pipe_flux_downs1_host;
		double *pipe_flux_downs2_host;

		double *p_ups_h_device;
		double *p_ups_Q_device;
		double *p_downs_h_device;
		double *p_downs_Q_device;

		int *upsPJ_device;
		int *downsPJ_device;
		int *downsPO_device;
		double *upNorm1_device;
		double *upNorm2_device;
		double *downNorm1_device;
		double *downNorm2_device;

		int *upsBCell_device;
		int *downsBCell_device;
	
		double *P_diameter_device;
		double *P_manning_device;
		double *pipe_Ap_device;
		double *dx_cell_device;
		double *t_pipe_device;   // time step of each pipe cell on device
		double *t_pipe; // time step of each pipe cell on host
	
		int *P_Obound_device; // bound type of the downstream outfall 

		void setHostMemory(int num_Pipe, int num_allPipeCell)
		{
			ups_PJ = new int[num_Pipe];
			downs_PJ = new int[num_Pipe];
			downs_PO = new int[num_Pipe];
			P_Obound_device = new int[num_Pipe];

			up_nor1 = new double[num_Pipe];
			up_nor2 = new double[num_Pipe];
			down_nor1 = new double[num_Pipe];
			down_nor2 = new double[num_Pipe];
			upstreamCell = new int[num_Pipe];
			downstreamCell = new int[num_Pipe];
			pipe_Ap = new double[num_Pipe];

			pipe_flux_ups1_host = new double[num_Pipe];
			pipe_flux_ups2_host = new double[num_Pipe];
			pipe_flux_downs1_host = new double[num_Pipe];
			pipe_flux_downs2_host = new double[num_Pipe];


			h_p_out=new double[num_allPipeCell];
			A_p_out=new double[num_allPipeCell];
			Q_p_out=new double[num_allPipeCell];
			//dx_cell = new double[num_allPipeCell];

			h_p = new double[num_allPipeCell];
			Q_p = new double[num_allPipeCell];
			sita_p = new double[num_allPipeCell];
			A_p = new double[num_allPipeCell];
			zb_p = new double[num_allPipeCell];
			t_pipe = new double[num_allPipeCell];
		

		}

		void setDeviceMemory(int num_Pipe, int num_allPipeCell)
		{
			cudaMalloc((void **) &zb_p_device, num_allPipeCell*sizeof(double));
			cudaMalloc((void **) &h_p_device, num_allPipeCell*sizeof(double));
			cudaMalloc((void **) &Q_p_device, num_allPipeCell*sizeof(double));
			cudaMalloc((void **) &sita_p_device, num_allPipeCell*sizeof(double));
			cudaMalloc((void **) &A_p_device, num_allPipeCell*sizeof(double));
			cudaMalloc((void **) &hn_p_device, num_allPipeCell*sizeof(double));
			cudaMalloc((void **) &Qn_p_device, num_allPipeCell*sizeof(double));
			cudaMalloc((void **) &An_p_device, num_allPipeCell*sizeof(double));
			cudaMalloc((void **) &sitan_p_device, num_allPipeCell*sizeof(double));

			cudaMalloc((void **) &pipe_flux_ups1, num_Pipe*sizeof(double));
			cudaMalloc((void **) &pipe_flux_ups2, num_Pipe*sizeof(double));
			cudaMalloc((void **) &pipe_flux_downs1, num_Pipe*sizeof(double));
			cudaMalloc((void **) &pipe_flux_downs2, num_Pipe*sizeof(double));

			cudaMalloc((void **) &p_ups_h_device, num_Pipe*sizeof(double));
			cudaMalloc((void **) &p_ups_Q_device, num_Pipe*sizeof(double));
			cudaMalloc((void **) &p_downs_h_device, num_Pipe*sizeof(double));
			cudaMalloc((void **) &p_downs_Q_device, num_Pipe*sizeof(double));

			cudaMalloc((void **) &upsPJ_device, num_Pipe*sizeof(int));
			cudaMalloc((void **) &downsPJ_device, num_Pipe*sizeof(int));
			cudaMalloc((void **) &downsPO_device, num_Pipe*sizeof(int));
			cudaMalloc((void **) &upNorm1_device, num_Pipe*sizeof(double));
			cudaMalloc((void **) &upNorm2_device, num_Pipe*sizeof(double));
			cudaMalloc((void **) &downNorm1_device, num_Pipe*sizeof(double));
			cudaMalloc((void **) &downNorm2_device, num_Pipe*sizeof(double));
			cudaMalloc((void **) &P_Obound_device, num_Pipe*sizeof(int));

			cudaMalloc((void **) &upsBCell_device,num_Pipe*sizeof(int));
			cudaMalloc((void **) &downsBCell_device,num_Pipe*sizeof(int));

			cudaMalloc((void **) &P_diameter_device,num_Pipe*sizeof(double));
			cudaMalloc((void **) &P_manning_device,num_Pipe*sizeof(double));
			cudaMalloc((void **) &pipe_Ap_device,num_Pipe*sizeof(double));

			cudaMalloc((void **) &dx_cell_device,num_allPipeCell*sizeof(double));
			cudaMalloc((void **) &t_pipe_device,num_allPipeCell*sizeof(double));
		
		}
	
		void deleteHostMemory()
		{
			delete[] ups_PJ;
			delete[] downs_PJ;
			delete[] downs_PO;

			delete[] up_nor1;
			delete[] up_nor2;
			delete[] down_nor1;
			delete[] down_nor2;
			delete[] upstreamCell;
			delete[] downstreamCell;
			delete[] pipe_Ap;
			delete[] h_p_out;
			delete[] A_p_out;
			delete[] Q_p_out;
			delete[] dx_cell;

			delete[] h_p;
			delete[] Q_p;
			delete[] sita_p;
			delete[] A_p;
			delete[] zb_p;
			delete[] t_pipe;
		}

		void deleteDeviceMemory()
		{
			cudaFree(zb_p_device);
			cudaFree(h_p_device);
			cudaFree(Q_p_device);
			cudaFree(sita_p_device);
			cudaFree(A_p_device);
			cudaFree(hn_p_device);
			cudaFree(Qn_p_device);
			cudaFree(An_p_device);
			cudaFree(sitan_p_device);

			cudaFree(pipe_flux_ups1);
			cudaFree(pipe_flux_ups2);
			cudaFree(pipe_flux_downs1);
			cudaFree(pipe_flux_downs2);

			cudaFree(p_ups_h_device);
			cudaFree(p_ups_Q_device);
			cudaFree(p_downs_h_device);
			cudaFree(p_downs_Q_device);

			cudaFree(upsPJ_device);
			cudaFree(downsPJ_device);
			cudaFree(downsPO_device);
			cudaFree(upNorm1_device);
			cudaFree(upNorm2_device);
			cudaFree(downNorm1_device);
			cudaFree(downNorm2_device);

			cudaFree(upsBCell_device);
			cudaFree(downsBCell_device);

			cudaFree(P_diameter_device);
			cudaFree(P_manning_device);
			cudaFree(pipe_Ap_device);

			cudaFree(dx_cell_device);
			cudaFree(t_pipe_device);
	
			cudaFree(P_Obound_device);
		
		}

		void copyCudaMemory(int cellTotalNum, int num_Pipe, double *pipeD, double *pipeM, int *P_outfBoundType)
		{
			cudaMemcpy(zb_p_device, zb_p, cellTotalNum*sizeof(double), cudaMemcpyHostToDevice);
			cudaMemcpy(h_p_device, h_p, cellTotalNum*sizeof(double), cudaMemcpyHostToDevice);
			cudaMemcpy(Q_p_device, Q_p, cellTotalNum*sizeof(double), cudaMemcpyHostToDevice);
			cudaMemcpy(sita_p_device, sita_p, cellTotalNum*sizeof(double), cudaMemcpyHostToDevice);
			cudaMemcpy(A_p_device, A_p, cellTotalNum*sizeof(double), cudaMemcpyHostToDevice);

			cudaMemcpy(upsPJ_device, ups_PJ, num_Pipe*sizeof(int), cudaMemcpyHostToDevice);
			cudaMemcpy(downsPJ_device, downs_PJ, num_Pipe*sizeof(int), cudaMemcpyHostToDevice);
			cudaMemcpy(downsPO_device, downs_PO, num_Pipe*sizeof(int), cudaMemcpyHostToDevice);
			cudaMemcpy(upNorm1_device, up_nor1, num_Pipe*sizeof(double), cudaMemcpyHostToDevice);
			cudaMemcpy(upNorm2_device, up_nor2, num_Pipe*sizeof(double), cudaMemcpyHostToDevice);
			cudaMemcpy(downNorm1_device, down_nor1, num_Pipe*sizeof(double), cudaMemcpyHostToDevice);
			cudaMemcpy(downNorm2_device, down_nor2, num_Pipe*sizeof(double), cudaMemcpyHostToDevice);

			cudaMemcpy(upsBCell_device,upstreamCell, num_Pipe*sizeof(int), cudaMemcpyHostToDevice);
			cudaMemcpy(downsBCell_device,downstreamCell, num_Pipe*sizeof(int), cudaMemcpyHostToDevice);

			cudaMemcpy(P_diameter_device,pipeD, num_Pipe*sizeof(double), cudaMemcpyHostToDevice);
			cudaMemcpy(P_manning_device,pipeM, num_Pipe*sizeof(double), cudaMemcpyHostToDevice);
			cudaMemcpy(pipe_Ap_device,pipe_Ap, num_Pipe*sizeof(double), cudaMemcpyHostToDevice);
			cudaMemcpy(dx_cell_device,dx_cell, cellTotalNum*sizeof(double), cudaMemcpyHostToDevice);
	
			cudaMemcpy(P_Obound_device,P_outfBoundType, num_Pipe*sizeof(int), cudaMemcpyHostToDevice);
		}

	};

	//******************************************************************************************************************* 
	class Junc_attributes
	{
	public:

		int *J_index;
		double *J_elev;
		double *J_maxDep;
		double *J_radius;
		double *J_xcoor;
		double *J_ycoor;
		double *J_manning;
		double *J_initWaterDepth;
		double *J_initXFlowRate;
		double *J_initYFlowRate;
		int *J_surfGrid;
	

		void setHostMemory(int num_Junc)
		{
			J_index = new int[num_Junc];
			J_elev = new double[num_Junc];
			J_maxDep = new double[num_Junc];
			J_radius = new double[num_Junc];
			J_xcoor = new double[num_Junc];
			J_ycoor = new double[num_Junc];
			J_manning = new double[num_Junc];
			J_initWaterDepth = new double[num_Junc];
			J_initXFlowRate = new double[num_Junc];
			J_initYFlowRate = new double[num_Junc];
			J_surfGrid = new int[num_Junc];
		}

		void deleteHostMemory()
		{
			delete[] J_index;
			delete[] J_elev;
			delete[] J_maxDep;
			delete[] J_radius;
			delete[] J_xcoor;
			delete[] J_ycoor;
			delete[] J_manning;
			delete[] J_initWaterDepth;
			delete[] J_initXFlowRate;
			delete[] J_initYFlowRate;
			delete[] J_surfGrid;
		}
	};

	//******************************************************************************************************************* 
	class Junc_variables 
	{
	public:

		int width_upsJP;
		int width_downsJP;
		double mass;

		double *h_j;
		double *qx_j;
		double *qy_j;
		double *junc_Area;
		double *h_j_out;
		double *qx_j_out;
		double *qy_j_out;

		double *h_j_device;
		double *hn_j_device;
		double *qx_j_device;
		double *qxn_j_device;
		double *qy_j_device;
		double *qyn_j_device;
		double *zb_j_device;
		double *junc_Area_device;
		double *t_junc_device;
		double *t_junc;
		double *junc_MN_device;
		double *j_radius_device;
		double *maxDepth_device;
		double *Qexchange_device;
		double *Qe_device; // Qe is for calculating the exchanage flow rate at junctions, Qexchange is the result

		double *Ax_junc_device;
		double *Ay_junc_device;
		double *hydroPX_device;
		double *hydroPY_device;
		double *pipe_F1_device;
		double *pipe_F2_device;
		double *pipe_F3_device;

		// the memory of the four following variables is created in the juncToPipeMatrix function
		int *ups_JP;
		int *downs_JP;
		int *upsJP_device;
		int *downsJP_device;

		int *surfGrid_device;
	

	
		void setHostMemory(int num_Junc)
		{
			h_j = new double[num_Junc];
			qx_j = new double[num_Junc];
			qy_j = new double[num_Junc];
			junc_Area = new double[num_Junc];
			h_j_out=new double[num_Junc];
			qx_j_out=new double[num_Junc];
			qy_j_out=new double[num_Junc];
			t_junc=new double[num_Junc];
		}

		void setDeviceMemory(int num_Junc)
		{
			cudaMalloc((void **) &h_j_device, num_Junc*sizeof(double));
			cudaMalloc((void **) &hn_j_device, num_Junc*sizeof(double));
			cudaMalloc((void **) &qx_j_device, num_Junc*sizeof(double));
			cudaMalloc((void **) &qxn_j_device, num_Junc*sizeof(double));
			cudaMalloc((void **) &qy_j_device, num_Junc*sizeof(double));
			cudaMalloc((void **) &qyn_j_device, num_Junc*sizeof(double));

			cudaMalloc((void **) &zb_j_device,num_Junc*sizeof(double));
			cudaMalloc((void **) &junc_Area_device,num_Junc*sizeof(double));
			cudaMalloc((void **) &t_junc_device,num_Junc*sizeof(double));
			cudaMalloc((void **) &junc_MN_device,num_Junc*sizeof(double));
			cudaMalloc((void **) &j_radius_device,num_Junc*sizeof(double));
			cudaMalloc((void **) &maxDepth_device,num_Junc*sizeof(double));
			cudaMalloc((void **) &Qexchange_device,num_Junc*sizeof(double));
			cudaMalloc((void **) &Qe_device,num_Junc*sizeof(double));

			cudaMalloc((void **) &surfGrid_device,num_Junc*sizeof(int));

			cudaMalloc((void **) &Ax_junc_device,num_Junc*sizeof(double));
			cudaMalloc((void **) &Ay_junc_device,num_Junc*sizeof(double));
			cudaMalloc((void **) &hydroPX_device,num_Junc*sizeof(double));
			cudaMalloc((void **) &hydroPY_device,num_Junc*sizeof(double));
			cudaMalloc((void **) &pipe_F1_device,num_Junc*sizeof(double));
			cudaMalloc((void **) &pipe_F2_device,num_Junc*sizeof(double));
			cudaMalloc((void **) &pipe_F3_device,num_Junc*sizeof(double));

		}

		void deleteHostMemory()
		{
			delete[] h_j;
			delete[] qx_j;
			delete[] qy_j;
			delete[] junc_Area;
			delete[] h_j_out;
			delete[] qx_j_out;
			delete[] qy_j_out;

			delete[] ups_JP;
			delete[] downs_JP;
			delete[] t_junc;
	
		
		}

		void deleteDeviceMemory()
		{
			cudaFree(h_j_device);
			cudaFree(hn_j_device);
			cudaFree(qx_j_device);
			cudaFree(qxn_j_device);
			cudaFree(qy_j_device);
			cudaFree(qyn_j_device);

			cudaFree(zb_j_device);
			cudaFree(junc_Area_device);
			cudaFree(t_junc_device);
			cudaFree(junc_MN_device);

			cudaFree(upsJP_device);
			cudaFree(downsJP_device);
			cudaFree(surfGrid_device);
			cudaFree(j_radius_device);
			cudaFree(maxDepth_device);
			cudaFree(Qexchange_device);
			cudaFree(Qe_device);

			cudaFree(Ax_junc_device);
			cudaFree(Ay_junc_device);
			cudaFree(hydroPX_device);
			cudaFree(hydroPY_device);
			cudaFree(pipe_F1_device);
			cudaFree(pipe_F2_device);
			cudaFree(pipe_F3_device);
		
		}

		void copyCudaMemory(int num_Junc, double *J_elev, double *J_manning, int *J_surfGrid, double *J_radius, double *maxDepth)
		{
			//printf("width_upsJP=%d\n",width_upsJP);
			cudaMemcpy(h_j_device, h_j, num_Junc*sizeof(double), cudaMemcpyHostToDevice);
			cudaMemcpy(qx_j_device, qx_j, num_Junc*sizeof(double), cudaMemcpyHostToDevice);
			cudaMemcpy(qy_j_device, qy_j, num_Junc*sizeof(double), cudaMemcpyHostToDevice);

			cudaMemcpy(upsJP_device, ups_JP, num_Junc*width_upsJP*sizeof(int), cudaMemcpyHostToDevice);
			cudaMemcpy(downsJP_device, downs_JP, num_Junc*width_downsJP*sizeof(int), cudaMemcpyHostToDevice);

			cudaMemcpy(zb_j_device,J_elev, num_Junc*sizeof(double), cudaMemcpyHostToDevice);
			cudaMemcpy(junc_Area_device,junc_Area, num_Junc*sizeof(double), cudaMemcpyHostToDevice);
			cudaMemcpy(junc_MN_device,J_manning, num_Junc*sizeof(double), cudaMemcpyHostToDevice);
			cudaMemcpy(j_radius_device,J_radius, num_Junc*sizeof(double), cudaMemcpyHostToDevice);
			cudaMemcpy(maxDepth_device,maxDepth, num_Junc*sizeof(double), cudaMemcpyHostToDevice);

			cudaMemcpy(surfGrid_device,J_surfGrid, num_Junc*sizeof(int), cudaMemcpyHostToDevice);
		}
	
	};
	//******************************************************************************************************************* 
	class Outf_attributes
	{
	public:

		int *O_index;
		double *O_elev;
		double *O_maxDep;
		double *O_xcoor;
		double *O_ycoor;
		double *O_manning;
		double *O_WaterDepth;
		double *O_XFlowRate;
		double *O_YFlowRate;
		int *O_surfGrid;
		int *O_boundType;

		void setHostMemory(int num_Outf)
		{
			O_index = new int[num_Outf];
			O_elev = new double[num_Outf];
			O_maxDep = new double[num_Outf];
			O_xcoor = new double[num_Outf];
			O_ycoor = new double[num_Outf];
			O_manning = new double[num_Outf];
			O_WaterDepth = new double[num_Outf];
			O_XFlowRate = new double[num_Outf];
			O_YFlowRate = new double[num_Outf];
			O_surfGrid = new int[num_Outf];
			O_boundType = new int[num_Outf];
		}

		void deleteHostMemory()
		{
			delete[] O_index;
			delete[] O_elev;
			delete[] O_maxDep;
			delete[] O_xcoor;
			delete[] O_ycoor;
			delete[] O_manning;
			delete[] O_WaterDepth;
			delete[] O_XFlowRate;
			delete[] O_YFlowRate;
			delete[] O_surfGrid;
			delete[] O_boundType;
		}
	};



	//******************************************************************************************************************* 
	class Outf_variables 
	{
	public:

	
		int width_downsOP;
		double mass;

		double *h_O;
		double *qx_O;
		double *qy_O;
		double *Outf_Area;
		double *h_O_out; // transffered to the CPU to calculate mass 

		double *h_O_device;
		double *hn_O_device;
		double *qx_O_device;
		double *qxn_O_device;
		double *qy_O_device;
		double *qyn_O_device;
		double *zb_O_device;
		double *Outf_Area_device;
		double *t_Outf_device;
		double *t_Outf;
		double *Outf_MN_device;
		double *Outf_radius_device;
		double *maxDepth_device;
		double *Qexchange_device;

		// the memory of the four following variables is created in the juncToPipeMatrix function
	
		int *downs_OP;
	
		int *downsOP_device;

		int *surfGrid_device;
	
		int *boundType_device;
	
		void setHostMemory(int num_Outf)
		{
			h_O = new double[num_Outf];
			qx_O = new double[num_Outf];
			qy_O = new double[num_Outf];
			Outf_Area = new double[num_Outf];
			h_O_out=new double[num_Outf];
			t_Outf=new double[num_Outf];
		}

		void setDeviceMemory(int num_Outf)
		{
			cudaMalloc((void **) &h_O_device, num_Outf*sizeof(double));
			cudaMalloc((void **) &hn_O_device, num_Outf*sizeof(double));
			cudaMalloc((void **) &qx_O_device, num_Outf*sizeof(double));
			cudaMalloc((void **) &qxn_O_device, num_Outf*sizeof(double));
			cudaMalloc((void **) &qy_O_device, num_Outf*sizeof(double));
			cudaMalloc((void **) &qyn_O_device, num_Outf*sizeof(double));

			cudaMalloc((void **) &zb_O_device,num_Outf*sizeof(double));
			cudaMalloc((void **) &Outf_Area_device,num_Outf*sizeof(double));
			cudaMalloc((void **) &t_Outf_device,num_Outf*sizeof(double));
			cudaMalloc((void **) &Outf_MN_device,num_Outf*sizeof(double));
			cudaMalloc((void **) &Outf_radius_device,num_Outf*sizeof(double));
			cudaMalloc((void **) &maxDepth_device,num_Outf*sizeof(double));
			cudaMalloc((void **) &Qexchange_device,num_Outf*sizeof(double));

			cudaMalloc((void **) &surfGrid_device,num_Outf*sizeof(int));
			cudaMalloc((void **) &boundType_device,num_Outf*sizeof(int));

		}

		void deleteHostMemory()
		{
			delete[] h_O;
			delete[] qx_O;
			delete[] qy_O;
			delete[] Outf_Area;
			delete[] h_O_out;

			//delete[] ups_JP;
			delete[] downs_OP;
			delete[] t_Outf;
	
		
		}

		void deleteDeviceMemory()
		{
			cudaFree(h_O_device);
			cudaFree(hn_O_device);
			cudaFree(qx_O_device);
			cudaFree(qxn_O_device);
			cudaFree(qy_O_device);
			cudaFree(qyn_O_device);

			cudaFree(zb_O_device);
			cudaFree(Outf_Area_device);
			cudaFree(t_Outf_device);
			cudaFree(Outf_MN_device);

			//cudaFree(upsJP_device);
			cudaFree(downsOP_device);
			cudaFree(surfGrid_device);
			cudaFree(Outf_radius_device);
			cudaFree(maxDepth_device);
			cudaFree(Qexchange_device);
			cudaFree(boundType_device);
		
		}

		void copyCudaMemory(int num_Outf, double *O_elev, double *O_manning, int *O_surfGrid, double *maxDepth, int *O_boundType)
		{
	
			cudaMemcpy(h_O_device, h_O, num_Outf*sizeof(double), cudaMemcpyHostToDevice);
			cudaMemcpy(qx_O_device, qx_O, num_Outf*sizeof(double), cudaMemcpyHostToDevice);
			cudaMemcpy(qy_O_device, qy_O, num_Outf*sizeof(double), cudaMemcpyHostToDevice);

		
			cudaMemcpy(downsOP_device, downs_OP, num_Outf*width_downsOP*sizeof(int), cudaMemcpyHostToDevice);

			cudaMemcpy(zb_O_device,O_elev, num_Outf*sizeof(double), cudaMemcpyHostToDevice);
			cudaMemcpy(Outf_Area_device,Outf_Area, num_Outf*sizeof(double), cudaMemcpyHostToDevice);
			cudaMemcpy(Outf_MN_device,O_manning, num_Outf*sizeof(double), cudaMemcpyHostToDevice);
			cudaMemcpy(maxDepth_device,maxDepth, num_Outf*sizeof(double), cudaMemcpyHostToDevice);

			cudaMemcpy(surfGrid_device,O_surfGrid, num_Outf*sizeof(int), cudaMemcpyHostToDevice);
			cudaMemcpy(boundType_device,O_boundType, num_Outf*sizeof(int), cudaMemcpyHostToDevice);

		}
	
	};

    //******************************************************************************************************************* 
class Drain_outfBounds
{
public:
	
	double *T_series;
	double *T_series_device;
	double *h_series;
	double *h_series_device;
	double *Q_series;
	double *Q_series_device;
	int length_T;
	int length_h;
	int length_Q;
	int *T_start;
	int *T_start_device;
	int *h_start;
	int *h_start_device;
	int *Q_start;
	int *Q_start_device;
	int *intervalLen;
	int *intervalLen_device;
	int *interval_ptr;
	int *interval_ptr_device;

	void setHostMemory(int length_T, int length_h, int length_Q, int outfall_num)
	{
		
		T_series = new double[length_T];
		h_series = new double[length_h];
		Q_series = new double[length_Q];
		T_start = new int[outfall_num];
		h_start = new int[outfall_num];
		Q_start = new int[outfall_num];
		intervalLen = new int[outfall_num];
		interval_ptr = new int[outfall_num];
		
	}

	void setDeviceMemory(int length_T, int length_h, int length_Q, int outfall_num)
	{
		
		cudaMalloc((void **) &T_series_device, length_T*sizeof(double));
		cudaMalloc((void **) &h_series_device, length_h*sizeof(double));
		cudaMalloc((void **) &Q_series_device, length_Q*sizeof(double));
		cudaMalloc((void **) &T_start_device, outfall_num*sizeof(int));
		cudaMalloc((void **) &h_start_device, outfall_num*sizeof(int));
		cudaMalloc((void **) &Q_start_device, outfall_num*sizeof(int));
		cudaMalloc((void **) &intervalLen_device, outfall_num*sizeof(int));
		cudaMalloc((void **) &interval_ptr_device, outfall_num*sizeof(int));

	}


	void deleteHostMemory()
	{
	
		delete[] T_series;
		delete[] h_series;
		delete[] Q_series;
		delete[] T_start;
		delete[] h_start;
		delete[] Q_start;
		delete[] intervalLen;
		delete[] interval_ptr;
	}

	void deleteDeviceMemory()
	{
		
		cudaFree(T_series_device);
		cudaFree(h_series_device);
		cudaFree(Q_series_device);
		cudaFree(T_start_device);
		cudaFree(h_start_device);
		cudaFree(Q_start_device);
		cudaFree(intervalLen_device);
		cudaFree(interval_ptr_device);
	}


	void copyCudaMemory(int length_T, int length_h, int length_Q, int outfall_num)
	{
		
		cudaMemcpy(T_series_device, T_series, length_T*sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(h_series_device, h_series, length_h*sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(Q_series_device, Q_series, length_Q*sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(T_start_device, T_start, outfall_num*sizeof(int), cudaMemcpyHostToDevice);
		cudaMemcpy(h_start_device, h_start, outfall_num*sizeof(int), cudaMemcpyHostToDevice);
		cudaMemcpy(Q_start_device, Q_start, outfall_num*sizeof(int), cudaMemcpyHostToDevice);
		cudaMemcpy(intervalLen_device, intervalLen, outfall_num*sizeof(int), cudaMemcpyHostToDevice);
		cudaMemcpy(interval_ptr_device, interval_ptr, outfall_num*sizeof(int), cudaMemcpyHostToDevice);
	}

public:
	
	void initIntervalPtr(int outfall_num)
	{
		//interval_ptr=0;
		for (int i=0;i<outfall_num;i++)
		{
			interval_ptr[i]=0;
		}
	}


};

//******************************************************************************************************************* 
class Drain_gauges
{
public: 
	int pGaugeNum;
	int jGaugeNum;
	double *pGaugeh;
	double *pGaugeQ;
	double *jGaugeh;
	int *pGaugePosi;
	int *jGaugePosi;
	void setHostMemory(int pGaugeNum, int jGaugeNum)
	{
		pGaugeh = new double[pGaugeNum];
		pGaugeQ = new double[pGaugeNum];
		jGaugeh = new double[jGaugeNum];
		pGaugePosi = new int[pGaugeNum];
		jGaugePosi = new int[jGaugeNum];
	}
	void deleteHostMemory()
	{
		delete[] pGaugeh;
		delete[] pGaugeQ;
		delete[] jGaugeh;
		delete[] pGaugePosi;
		delete[] jGaugePosi;
	}
};



	//******************************************************************************************************************* 
	class Const_variables
	{
	public:
		double pi;
		double g;
		double kesi;
		double CFL;
		double dT_Default;
		int Maxtrials;
	
		Const_variables()
		{
			pi=3.1415926;
			g=9.81;
			kesi=1.0e-10;
			Maxtrials=15;
			CFL=0.8;
			dT_Default=0.005;
		}
		__device__ double MIN(double x, double y)
		{
			// MIN(x, y) (((x) < (y)) ? (x) : (y))
			if (x<y)
				return x;
			else
				return y;

		}

		__device__ double MAX(double x, double y)
		{
			if (x>y)
				return x;
			else
				return y;
		}

		double fmin(double x, double y)
		{
			if (x<y)
				return x;
			else
				return y;

		}
	};


	


}

#endif











