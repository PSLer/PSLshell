%%PSLshell
%%Date: 	2025.09.30
%%Author:	Junpeng Wang (junpeng.wang@tum.de; junwa@dtu.dk)

%%This repository is the associated software of the paper: 
%% "Topology-aware Stress Analysis in Shell Structures"
%% by Junpeng Wang, Yingjian Liu, Jun Wu, and Rüdiger Westermann.
%% Preprint can be found at https://papers.ssrn.com/sol3/papers.cfm?abstract_id=5776659

clear all; clc;
GlobalVariables();

%%1. Import Data
tStart = tic;
stressfileName = './data/cylinder_1_Q4';
ImportStressFields(strcat(stressfileName, '.TSV'));
disp(['Import Stress Field Costs: ' sprintf('%10.3g',toc(tStart)) 's']);

%%2. Pre-process Data
tStart = tic;
PreProcessData();
disp(['Pre-process Data Costs: ' sprintf('%10.3g',toc(tStart)) 's']);
figure; ShowProblemDescription();

if 1 %% PSLs Generation
	%%3. Seeding PSLs
	tStart = tic;
	profile = struct(...
		'PSLsDensityCtrl',		10,				...	%% minDim/PSLsDensityCtrl ('\omega in the paper')
		'PSLtypeCtrl',			[0 0], 			... %% [major, minor], 0==generate, 1==skip
		'topologyAware',		1, 				... %% Enable Topology Analysis
		'distanceMetric',		'Euclidean', 	... %% 'Euclidean', 'Geodesic'
		'seedSparsityCtrl',		1,				... %% Every 'value' element is sampled
		'maxIts',				1000			... %% Maximum Integration Steps in a Single Run of PSL Tracing
	);
	%% =================================Settings for Paper replicability=================================
	%% #Planar Plate
	% profile = struct('PSLsDensityCtrl', 30, 'PSLtypeCtrl', [0 0], 'topologyAware', 1, 'distanceMetric', 'Euclidean', 'seedSparsityCtrl', 1, 'maxIts', 1000);
	%% #Cylinder; #Dome
	% profile = struct('PSLsDensityCtrl', 10, 'PSLtypeCtrl', [0 0], 'topologyAware', 1, 'distanceMetric', 'Euclidean', 'seedSparsityCtrl', 1, 'maxIts', 1000);
	%% #Polyhedra-like Structure
	% profile = struct('PSLsDensityCtrl', 15, 'PSLtypeCtrl', [0 0], 'topologyAware', 1, 'distanceMetric', 'Euclidean', 'seedSparsityCtrl', 1, 'maxIts', 1000);
	%% #Wing
	% profile = struct('PSLsDensityCtrl', 3, 'PSLtypeCtrl', [0 0], 'topologyAware', 1, 'distanceMetric', 'Geodesic', 'seedSparsityCtrl', 1, 'maxIts', 2000);
	%% #Wind turbine blade
	% profile = struct('PSLsDensityCtrl', 5, 'PSLtypeCtrl', [0 0], 'topologyAware', 1, 'distanceMetric', 'Geodesic', 'seedSparsityCtrl', 1, 'maxIts', 2000);
	%% #TPMS
	% profile = struct('PSLsDensityCtrl', 40, 'PSLtypeCtrl', [0 0], 'topologyAware', 1, 'distanceMetric', 'Euclidean', 'seedSparsityCtrl', 1, 'maxIts', 2000);    
	SeedingPSLs(profile); 
	disp(['Seeding PSL Costs: ' sprintf('%10.3g',toc(tStart)) 's']);
	
	%%3. Vis.
	ShowPSLs('None', 'None', 1.0, 20);
	% ShowPrincipalStressDirectionsByArrowsGlobally(); %%Show Principal Stress Directions via Arrows
else %% Check Mesh Normals
	%%3. Evaluate Element Normals
	faceNormals = EvaluateMeshNormals();
	figure; ShowMeshNormals(faceNormals);
	% optFlipNormal = 1;
	% if optFlipNormal %% If the normals are not outward
		% stressfileNameNew = strcat(stressfileName, '_new');
		% ExportStressFieldWithFlippedNormals(strcat(stressfileNameNew, '.TSV'));
	% end	
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Source
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function GlobalVariables()
	global boundingBox_; boundingBox_ = [];
	global meshType_; meshType_ = [];
	global frameType_; frameType_ = [];
	global meshOrder_; meshOrder_ = [];
	global numNodes_; numNodes_ = 0;
	global nodeCoords_; nodeCoords_ = [];
	global numEles_; numEles_ = 0;
	global eNodMat_; eNodMat_ = [];
	global eleCentroids_; eleCentroids_ = [];
	global meshTypeMap_; meshTypeMap_ = [];
	global localStressField_; localStressField_ = [];
	global globalStressField_; globalStressField_ = [];
	global seedPointsRef_; seedPointsRef_ = [];
	global silhouetteStruct_; silhouetteStruct_ = [];
	global eleCharacterSizeList_; eleCharacterSizeList_ = [];
	global shellEleNormals_; shellEleNormals_ = [];
	global shellCurvatureScaling_; shellCurvatureScaling_ = [];
		
	global nodStruct_; nodStruct_ = [];
	global eleStruct_; eleStruct_ = [];
	global boundaryElements_; boundaryElements_ = [];
	
	%% 2. Macros
	global refVec_; refVec_ = [1 0 0]; refVec_ = refVec_/norm(refVec_);
	global refVecFallback_; refVecFallback_ = [0 0 1]; refVecFallback_ = refVecFallback_/norm(refVecFallback_);
	global tolRefVecFallback_; tolRefVecFallback_ = 1.0e-6;
	
	%% 3. Algorithm Control
	%% True: Directly Conduct Stress Interpolation in Local Frame; 
	%% False: Convert Per-element Stress State to Global Frame, Conduct Element Interpolation, Convert back to Local Frame
	global tracingSche_; tracingSche_ = 'RK2'; %% 'Euler' 'RK2'
	global approxiInterp_; approxiInterp_ = 0; 
	global scalingStepSize4PSL_; scalingStepSize4PSL_ = 5;
	global limitSteps_; limitSteps_ = 2000; 
	% %% 2.2 %% Tracing PSL stops when the angle deviation between the neighboring tangents is larger than permittedMaxAdjacentTangentAngleDeviation_
	global permittedMaxAdjacentTangentAngleDeviation_; permittedMaxAdjacentTangentAngleDeviation_ = 45; %% in degree
	
	% %% 3. Result
	global majorPSLpool_; majorPSLpool_ = PrincipalStressLineStruct();
	global minorPSLpool_; minorPSLpool_ = PrincipalStressLineStruct();
	global candidateElementsIncDegePts_; candidateElementsIncDegePts_ = [];
	global degePts_; degePts_ = [];
end


function ImportStressFields(fileName)
	global frameType_ meshOrder_ numNodes_ nodeCoords_ numEles_ eNodMat_ meshTypeMap_;
	global localStressField_ globalStressField_ loadingCond_ fixingCond_ nodeWiseStressField_;
	global refVec_ refVecFallback_;
	
	%%1. Read mesh and cartesian stress field
	fid = fopen(fileName, 'r');
	%%1.1 Mesh
	% fgetl(fid);
	tmp = fscanf(fid, '%s', 1);
	dataVersionID = fscanf(fid, '%f', 1);
	tmp = fscanf(fid, '%s %s %s', 3);
	dataType = fscanf(fid, '%s', 1);
	switch dataType
		case 'NODE', nodeWiseStressField_ = 1;
		case 'ELEMENT', nodeWiseStressField_ = 0; %%Element-wise Stress Data
		otherwise, error('Un-supported Stress Data!');
	end
	domainType = fscanf(fid, '%s', 1);
	if ~strcmp(domainType, 'Shell'), error('Un-supported Data!'); end
	%meshType_ = fscanf(fid, '%s', 1);
	%if ~(strcmp(meshType_, 'Quad') || strcmp(meshType_, 'Tri') || strcmp(meshType_, 'Hybrid')), error('Un-supported Mesh!'); end
	meshOrder_ = fscanf(fid, '%d', 1);
	if ~(1==meshOrder_ || 2==meshOrder_), error('Un-supported Mesh!'); end
	tmp = fscanf(fid, '%s', 1); frameType_ = fscanf(fid, '%s', 1);
	if strcmp(frameType_, 'LOCAL')
		tmp = fscanf(fid, '%s %s', 2); refVec_ = fscanf(fid, '%e %e %e', [3, 1])';
		tmp = fscanf(fid, '%s %s', 2); refVecFallback_ = fscanf(fid, '%e %e %e', [3, 1])';
	end
	startReadingVertices = fscanf(fid, '%s', 1);
	if ~strcmp(startReadingVertices, 'Vertices:'), error('Un-supported Data!'); end
	numNodes_ = fscanf(fid, '%d', 1);
	nodeCoords_ = fscanf(fid, '%e %e %e', [3, numNodes_])'; 
	startReadingElements = fscanf(fid, '%s', 1);
	if ~strcmp(startReadingElements, 'Elements:'), error('Un-supported Data!'); end
	numEles_ = fscanf(fid, '%d', 1);
	switch meshOrder_
		case 1
			eNodMat_ = NaN(numEles_, 4);
			meshTypeMap_ = repmat("Q4", numEles_,1);
			for ii=1:numEles_
				iNumNodesPerEle = fscanf(fid, '%d', 1);
				switch iNumNodesPerEle
					case 3
						eNodMat_(ii,1:3) = fscanf(fid, '%d %d %d', [3, 1])';
						meshTypeMap_(ii) = "T3";
					case 4
						eNodMat_(ii,1:4) = fscanf(fid, '%d %d %d %d', [4, 1])';
				end
			end					
		case 2
			eNodMat_ = NaN(numEles_, 8);
			meshTypeMap_ = repmat("Q8", numEles_,1);
			for ii=1:numEles_
				iNumNodesPerEle = fscanf(fid, '%d', 1);
				switch iNumNodesPerEle
					case 6
						eNodMat_(ii,1:6) = fscanf(fid, '%d %d %d %d %d %d', [6, 1])';
						meshTypeMap_(ii) = "T6";
					case 8
						eNodMat_(ii,1:8) = fscanf(fid, '%d %d %d %d %d %d %d %d', [8, 1])';
				end
			end					
	end
	%%1.2 Stress Related
	startReadingLoads = fscanf(fid, '%s %s', 2); 
	if ~strcmp(startReadingLoads, 'NodeForces:'), error('Un-supported Data!'); end
	numLoadedNodes = fscanf(fid, '%d', 1);
	if numLoadedNodes>0, loadingCond_ = fscanf(fid, '%d %e %e %e %e %e %e', [7, numLoadedNodes])'; else, loadingCond_ = []; end
	
	startReadingFixations = fscanf(fid, '%s %s', 2);
	if ~strcmp(startReadingFixations, 'FixedNodes:'), error('Un-supported Data!'); end
	numFixedNodes = fscanf(fid, '%d', 1);
	if numFixedNodes>0, fixingCond_ = fscanf(fid, '%d %d %d %d %d %d %d', [7, numFixedNodes])'; else, fixingCond_ = []; end
	
	startReadingStress = fscanf(fid, '%s %s', 2); 
	if ~strcmp(startReadingStress, 'CartesianStress:'), error('Un-supported Data!'); end
	numValidNods = fscanf(fid, '%d', 1);
	switch frameType_
		case 'LOCAL'
			%% Per-row: sigma_xx, sigma_yy, sigma_xy
			localStressField_ = fscanf(fid, '%e %e %e', [3, numValidNods])';
		case 'GLOBAL'
			%% Per-row: sigma_xx, sigma_yy, sigma_zz, sigma_yz, sigma_zx, sigma_xy
			globalStressField_ = fscanf(fid, '%e %e %e %e %e %e', [6, numValidNods])';
	end
	fclose(fid);	
end


function PreProcessData()
	global boundingBox_ frameType_ meshOrder_;
	global numNodes_ nodeCoords_ numEles_ eNodMat_ eleCentroids_ meshTypeMap_ voidEleMap_;
	global localStressField_ globalStressField_ localStressFieldPerEle_ globalStressFieldPerEle_;
	global silhouetteStruct_ eleCharacterSizeList_ refScalingSize_ shellEleNormals_ shellCurvatureScaling_;
	global nodStruct_ eleStruct_ boundaryElements_; 
	global triangularizedMesh_ triangularizedMeshGraph_;
	
	%%1. Extract Silhouette for Vis.
	boundingBox_ = [min(nodeCoords_, [], 1); max(nodeCoords_, [], 1)];
	silhouetteStruct_.vertices = nodeCoords_;
	silhouetteStruct_.faces = eNodMat_;
	if 2==meshOrder_
		silhouetteStruct_.faces('Q8'==meshTypeMap_,:) = silhouetteStruct_.faces('Q8'==meshTypeMap_,[1 5 2 6 3 7 4 8]);
		silhouetteStruct_.faces('T6'==meshTypeMap_,:) = silhouetteStruct_.faces('T6'==meshTypeMap_,[1 4 2 5 3 6 7 8]);	
	end
	eleCentroids_ = zeros(numEles_, 3);
	for ii=1:numEles_
		iEleType = meshTypeMap_(ii);
		switch iEleType
			case 'T3'
				iNumNodesPerEle = 3;
				iParas = [1 1]/3;
			case 'Q4'
				iNumNodesPerEle = 4;
				iParas = [0 0];
			case 'T6'
				iNumNodesPerEle = 6;
				iParas = [1 1]/3;
			case 'Q8'
				iNumNodesPerEle = 8;
				iParas = [0 0];
		end
		iSF = ShapeFunction(iParas, iEleType);
		eleCentroids_(ii,:) = iSF * nodeCoords_(eNodMat_(ii,1:iNumNodesPerEle),:);
	end
	
	%%2. Relate the Shared Element for each Vertex Node in unstructured shell mesh
	iNodStruct = struct('adjacentEles', []);
	nodStruct_ = repmat(iNodStruct, numNodes_, 1);
	for ii=1:numEles_
		switch meshTypeMap_(ii)
			case 'T3'
				iNumVertexNodesPerEle = 3;
			case 'Q4'
				iNumVertexNodesPerEle = 4;
			case 'T6'
				iNumVertexNodesPerEle = 3;
			case 'Q8'
				iNumVertexNodesPerEle = 4;
		end
		for jj=1:iNumVertexNodesPerEle
			nodStruct_(eNodMat_(ii,jj)).adjacentEles(1,end+1) = ii;
		end	
	end
	
	%%3. Relate the Shared Element of each Edge
	iEleStruct = struct('type', 0, 'edge2nextEle', zeros(1,4), 'edgeIDsOfNextEle', zeros(1,4), 'vtxEdgesGlobal', zeros(4,2));
	eleStruct_ = repmat(iEleStruct, numEles_, 1);
	for ii=1:numEles_
		switch meshTypeMap_(ii)
			case 'T3'
				eleStruct_(ii).type = 3; %%T3
				iNumVertexEdgesPerEle = 3;
				iNodesPerEle = eNodMat_(ii,1:3);
				edgePerEle = iNodesPerEle([1 2; 2 3; 3 1]);
				eleStruct_(ii).edge2nextEle(end) = [];
				eleStruct_(ii).edgeIDsOfNextEle(end) = [];
				eleStruct_(ii).vtxEdgesGlobal = edgePerEle;
			case 'Q4'
				eleStruct_(ii).type = 4; %%Q4
				iNumVertexEdgesPerEle = 4;
				iNodesPerEle = eNodMat_(ii,1:4);
				edgePerEle = iNodesPerEle([1 2; 2 3; 3 4; 4 1]);
				eleStruct_(ii).vtxEdgesGlobal = edgePerEle;
			case 'T6'
				eleStruct_(ii).type = 6; %%T6
				iNumVertexEdgesPerEle = 3;
				iNodesPerEle = eNodMat_(ii,1:3);
				edgePerEle = iNodesPerEle([1 2; 2 3; 3 1]);
				eleStruct_(ii).edge2nextEle(end) = [];
				eleStruct_(ii).edgeIDsOfNextEle(end) = [];
				eleStruct_(ii).vtxEdgesGlobal(end,:) = [];
				eleStruct_(ii).vtxEdgesGlobal = edgePerEle;
			case 'Q8'
				eleStruct_(ii).type = 8; %%Q8
				iNumVertexEdgesPerEle = 4;
				iNodesPerEle = eNodMat_(ii,1:4);
				edgePerEle = iNodesPerEle([1 2; 2 3; 3 4; 4 1]);
				eleStruct_(ii).vtxEdgesGlobal = edgePerEle;
		end
		for jj=1:iNumVertexEdgesPerEle
			iEdge = edgePerEle(jj,:);
			eleGroup1 = nodStruct_(iEdge(1)).adjacentEles;
			eleGroup2 = nodStruct_(iEdge(2)).adjacentEles;
			sharedEles = intersect(eleGroup1, eleGroup2);
			numSharedEles = numel(sharedEles);				
			if 2==numSharedEles
				nextEle = setdiff(sharedEles, ii);
				if 1~=numel(nextEle)
					error('Un-expected Error in Mesh Structure 1!');
				end
				eleStruct_(ii).edge2nextEle(jj) = nextEle;
				switch meshTypeMap_(nextEle)
					case 'T3'
						iNodesNextEle = eNodMat_(nextEle,1:3);
						edgesNextEle = iNodesNextEle([1 2; 2 3; 3 1]);
					case 'Q4'
						iNodesNextEle = eNodMat_(nextEle,1:4);
						edgesNextEle = iNodesNextEle([1 2; 2 3; 3 4; 4 1]);
					case 'T6'
						iNodesNextEle = eNodMat_(nextEle,1:3);
						edgesNextEle = iNodesNextEle([1 2; 2 3; 3 1]);					
					case 'Q8'
						iNodesNextEle = eNodMat_(nextEle,1:4);
						edgesNextEle = iNodesNextEle([1 2; 2 3; 3 4; 4 1]);					
				end
				edgeInCurrentEle = sort(iEdge);
				edgesNextEle = sort(edgesNextEle,2);
				[~, localEdgeIDinNextEle] = intersect(edgesNextEle, edgeInCurrentEle, 'rows');
				eleStruct_(ii).edgeIDsOfNextEle(jj) = localEdgeIDinNextEle;
			elseif 1==numSharedEles
				if sharedEles~=ii
					error('Un-expected Error in Mesh Structure 2!');
				end
			else
				error('Un-expected Error in Mesh Structure 3!');
			end
		end		
	end
	
	%%3.1 Extract Boundary Elements for Testing
	eleState = zeros(numEles_,1);
	for ii=1:numEles_
		if ~isempty(find(0==eleStruct_(ii).edge2nextEle))
			eleState(ii) = 1;
		end
	end
	boundaryElements_ = find(1==eleState);	
	
	%%4. Evaluate the Character Size of each Element (used in PSL tracing process)
	%% Roughly Minimum Edge
	eleCharacterSizeList_ = zeros(numEles_,1);
	for ii=1:numEles_
		switch meshTypeMap_(ii)
			case 'T3'
				iNodesPerEle = eNodMat_(ii,1:3);
				edgePerEle = [1 2; 2 3; 3 1];
			case 'Q4'
				iNodesPerEle = eNodMat_(ii,1:4);
				edgePerEle = [1 2; 2 3; 3 4; 4 1];
			case 'T6'
				iNodesPerEle = eNodMat_(ii,1:3);
				edgePerEle = [1 2; 2 3; 3 1];
			case 'Q8'
				iNodesPerEle = eNodMat_(ii,1:4);
				edgePerEle = [1 2; 2 3; 3 4; 4 1];
		end
		iEleNodeCoords = nodeCoords_(iNodesPerEle,:);
		iEleEdgeLengths = vecnorm(iEleNodeCoords(edgePerEle(:,2),:)-iEleNodeCoords(edgePerEle(:,1),:),2,2);
		eleCharacterSizeList_(ii) = min(iEleEdgeLengths);
	end
	refScalingSize_ = mean(eleCharacterSizeList_);	
	
	%%5. Evaluate the Curvature of each higher-order Element (used in PSL tracing process, scale the tracing step size)
	%% Rougly Maximum Directional Deviation of the Sampled Shell Normals per Element
	%% Scaling range: 1-5
	maxScalingFactor = 5;
	shellEleNormals_ = struct('normals', [], 'paras', []);
	shellEleNormals_ = repmat(shellEleNormals_, numEles_, 1);
	shellCurvatureScaling_ = ones(numEles_, 1);
	if meshOrder_>1
		for ii=1:numEles_
			switch meshTypeMap_(ii)
				case 'T6'
					iEleNodes = nodeCoords_(eNodMat_(ii,1:6),:);
					sampledPts = [
						1/3		1/3
						1/6		1/6
						2/3		1/6
						1/6		2/3
					];
					[faceNormList, ~, ~] = ComputeNormalsAtGivenPosition(sampledPts, iEleNodes, 'T6');
				case 'Q8'
					iEleNodes = nodeCoords_(eNodMat_(ii,:),:);
					sampledPts = [ %%Natural Coordinates
						 0.0	 			 0.0
						-0.5773502691896258	-0.5773502691896258
						 0.5773502691896258	-0.5773502691896258
						 0.5773502691896258	 0.5773502691896258
						-0.5773502691896258	 0.5773502691896258
					];
					[faceNormList, ~, ~] = ComputeNormalsAtGivenPosition(sampledPts, iEleNodes, 'Q8');
			end
			shellEleNormals_(ii).normal = faceNormList;
			shellEleNormals_(ii).paras = sampledPts;
			
			refNorm = faceNormList(1,:);
			testNorm = faceNormList(2:end,:);
			dirEval = acos(sum(refNorm .* testNorm, 2));
			iMaxDirDev = max(dirEval);
			if iMaxDirDev > 45
				warning(['Element ', sprintf('%4i',ii), ' is highly curved!']);
			end
			shellCurvatureScaling_(ii) = max(iMaxDirDev, 1);
		end
	end
	maxDirectionalDev = max(shellCurvatureScaling_);
	if maxDirectionalDev > 1
		tmp = 1 + (maxScalingFactor-1)*(shellCurvatureScaling_-1)/(maxDirectionalDev-1);
		shellCurvatureScaling_ = tmp;
	end
	
	%%6. Pre-compute per-element Node Stress Tensors in Global
	localStressFieldPerEle_ = cell(numEles_,1);
	globalStressFieldPerEle_ = cell(numEles_,1);
	voidEleMap_ = false(numEles_,1);
	for ii=1:numEles_
		iEleType = meshTypeMap_(ii);
		switch iEleType
			case 'T3'
				iEleNodes = eNodMat_(ii,1:3);
				parasNodes = [0.0 0.0; 1.0 0.0; 0.0 1.0];
				iEleNodeCoords = nodeCoords_(iEleNodes,:);
				nodeLocalFrames = zeros(3,3,3);
				for jj=1:3
					[nodeLocalFrames(:,:,jj), ~, ~, ~] = ComputeLocalFrameAtGivenPosition(parasNodes(jj,:), iEleNodeCoords, iEleType);
				end				
			case 'Q4'
				iEleNodes = eNodMat_(ii,1:4);
				parasNodes = [-1.0 -1.0; 1.0 -1.0; 1.0 1.0; -1.0 1.0];
				iEleNodeCoords = nodeCoords_(iEleNodes,:);	
				nodeLocalFrames = zeros(3,3,4);
				for jj=1:4
					[nodeLocalFrames(:,:,jj), ~, ~, ~] = ComputeLocalFrameAtGivenPosition(parasNodes(jj,:), iEleNodeCoords, iEleType);
				end						
			case 'T6'
				iEleNodes = eNodMat_(ii,1:6);
				parasNodes = [0.0 0.0; 1.0 0.0; 0.0 1.0; 0.5 0.0; 0.5 0.5; 0.0 0.5];
				iEleNodeCoords = nodeCoords_(iEleNodes,:);	
				nodeLocalFrames = zeros(3,3,6);
				for jj=1:6
					[nodeLocalFrames(:,:,jj), ~, ~, ~] = ComputeLocalFrameAtGivenPosition(parasNodes(jj,:), iEleNodeCoords, iEleType);
				end				
			case 'Q8'
				iEleNodes = eNodMat_(ii,1:8);
				parasNodes = [-1.0 -1.0; 1.0 -1.0; 1.0 1.0; -1.0 1.0; 0.0 -1.0; 1.0 0.0; 0.0 1.0; -1.0 0.0];
				iEleNodeCoords = nodeCoords_(iEleNodes,:);				
				nodeLocalFrames = zeros(3,3,8);
				for jj=1:8
					[nodeLocalFrames(:,:,jj), ~, ~, ~] = ComputeLocalFrameAtGivenPosition(parasNodes(jj,:), iEleNodeCoords, iEleType);
				end							
		end
		switch frameType_
			case 'LOCAL'
				iEleStressesLocal = localStressField_(iEleNodes,:);
				iEleStressesGlobal = GlobalFrame2Local_StressTensor(iEleStressesLocal, nodeLocalFrames, 0);	
			case 'GLOBAL'
				iEleStressesGlobal = globalStressField_(iEleNodes,:);
				iEleStressesLocal = GlobalFrame2Local_StressTensor(iEleStressesGlobal, nodeLocalFrames, 1);					
		end
		localStressFieldPerEle_{ii} = iEleStressesLocal;
		globalStressFieldPerEle_{ii} = iEleStressesGlobal;
		if ~isempty(find(0==ComputeVonMisesStressLocalFrame(iEleStressesLocal)))
			voidEleMap_(ii) = true;
		end
	end	

	%%7. Triangulation for Geodesic Distance Computing
	elesQ8 = find('Q8'==meshTypeMap_);
	if ~isempty(elesQ8)
		numQ8eles = numel(elesQ8);
		newlyInsertedNodes = zeros(numQ8eles, 3);
		Q8SFctr = ShapeFunction([0 0], 'Q8');
		for ii=1:numQ8eles
			iEleNodes = eNodMat_(elesQ8(ii),1:8);
			iEleNodeCoords = nodeCoords_(iEleNodes,:);
			newlyInsertedNodes(ii,:) = Q8SFctr * iEleNodeCoords;
		end
		mapQ8 = zeros(numEles_,1);
		mapQ8(elesQ8,1) = (1:numel(elesQ8))';
	else
		newlyInsertedNodes = [];
		mapQ8 = (1:numEles_)';
	end
	triangularizedNodeCoords = [nodeCoords_; newlyInsertedNodes];
	triangularizedElementsTmp = struct('arr', []);
	triangularizedElementsTmp = repmat(triangularizedElementsTmp, numEles_, 1);
	
	for ii=1:numEles_
		iEleType = meshTypeMap_(ii);
		switch iEleType
			case 'T3'
				iEleNodes = eNodMat_(ii,1:3);
			case 'Q4'
				iEleNodes = eNodMat_(ii,1:4);
				iEleNodes = iEleNodes([1 2 3  3 4 1]);
			case 'T6'
				iEleNodes = eNodMat_(ii,1:6);
				iEleNodes = iEleNodes([1 4 6  4 2 5  6 5 3  4 5 6]);
			case 'Q8'
				iEleNodes = zeros(1,9);
				iEleNodes(1:8) = eNodMat_(ii,1:8);
				iEleNodes(9) = numNodes_+mapQ8(ii);
				iEleNodes = iEleNodes([1 5 8  9 8 5  5 2 9  6 9 2  9 6 7  3 7 6  8 9 4  7 4 9]);
		end
		triangularizedElementsTmp(ii).arr = iEleNodes;
	end
	triangularizedElements = [triangularizedElementsTmp.arr]';
	numTriangularizedEles = numel(triangularizedElements)/3;
	triangularizedElements = reshape(triangularizedElements, 3, numTriangularizedEles)';
	triangularizedMesh_ = triangulation(triangularizedElements, triangularizedNodeCoords);
	triangularizedMeshEdges = edges(triangularizedMesh_);
	weights = vecnorm(triangularizedNodeCoords(triangularizedMeshEdges(:,1),:) - triangularizedNodeCoords(triangularizedMeshEdges(:,2),:), 2, 2);
	triangularizedMeshGraph_ = graph(triangularizedMeshEdges(:,1), triangularizedMeshEdges(:,2), weights);
end


function SeedingPSLs(profile)
	global boundingBox_ eleCentroids_;
	global mergingThreshold_ distanceMetric_ limitSteps_
	global iniSeedPoints_ seedPointsValence_ seedPoints_ seedPointsRef_
	global majorPSLpool_ minorPSLpool_ majorCoordList_ minorCoordList_
	
	%%1. Initialization
	limitSteps_ = profile.maxIts;
	distanceMetric_ = profile.distanceMetric;
	iniSeedPoints_ = SetupSeeds(profile.seedSparsityCtrl);
	seedPoints_ = iniSeedPoints_;
	seedPointsRef_ = [seedPoints_(:,1) eleCentroids_(seedPoints_(:,1),:)];
	numSeedPoints = size(iniSeedPoints_,1);
	assert(2==numel(profile.PSLtypeCtrl));
	seedPointsValence_ = repmat(profile.PSLtypeCtrl(:)', numSeedPoints, 1);
	mergingThreshold_ = min(boundingBox_(2,:)-boundingBox_(1,:))/profile.PSLsDensityCtrl;
	if 0==mergingThreshold_ %%2D
		mergingThreshold_ = max(boundingBox_(2,:)-boundingBox_(1,:))/profile.PSLsDensityCtrl;
	end
	
	majorPSLpool_ = PrincipalStressLineStruct();
	minorPSLpool_ = PrincipalStressLineStruct();
	majorCoordList_ = [];	
	minorCoordList_ = [];
	% ShowSeedingProcess(0.75);	
	[~, startEle] = min(vecnorm(sum(boundingBox_, 1)/2-eleCentroids_, 2, 2));
	startCoord = eleCentroids_(startEle,:);
	
	%%2. Topology Analysis
	if profile.topologyAware
		StressTopologyAnalysis();
		EmbedTopologicalSkeleton();
	end
	% ShowSeedingProcess(0.75);	
	%%3. Seeding (ref. 3D-TSV)
	its = 0;
	looper = sum(seedPointsValence_(:));
	while looper<2*numSeedPoints
		its = its + 1;
		valenceMetric = sum(seedPointsValence_,2);
		unFinishedSpps = find(valenceMetric<2);
		% spp = unFinishedSpps(1);
		if 0==looper
			[~, spp] = min(vecnorm(startCoord-seedPointsRef_(:,end-2:end),2,2));
		else
			tmp	= seedPointsValence_(unFinishedSpps,:);
			tmp = find(1==sum(tmp,2));
			if ~isempty(tmp), unFinishedSpps = unFinishedSpps(tmp); end
			[~, tarPos] = min(vecnorm(startCoord-seedPointsRef_(unFinishedSpps,end-2:end),2,2));
			spp = unFinishedSpps(tarPos);
		end
		valences = seedPointsValence_(spp,:);
		iSeed = seedPoints_(spp,:);
		if 0==valences(1)
			seedPointsValence_(spp,1) = 1;
			majorPSL = CreatePrincipalStressLine(iSeed, 'MAJOR');
			if 0==majorPSL.numIntergPts
				looper = sum(seedPointsValence_(:)); 
				disp([' Iteration.: ' sprintf('%4i',its) ' Progress.: ' sprintf('%6i',looper) ...
					' Total.: ' sprintf('%6i',2*numSeedPoints)]); continue; 			
			end
			majorPSLpool_(end+1,1) = majorPSL;
			majorCoordList_(end+1:end+majorPSL.numIntergPts,:) = majorPSL.phyCoordList;
			sppsEmptyMajorValence = find(0==seedPointsValence_(:,1));
			if ~isempty(sppsEmptyMajorValence)
				[potentialDisListMajor, potentialPosListMajor, potentialPosListMajorPara] = GetDisListOfPointList2Curve(...	
					seedPointsRef_(sppsEmptyMajorValence,:), [majorPSL.eleIndexList majorPSL.phyCoordList], [majorPSL.eleIndexList majorPSL.paraCoordList]);
				potentialSolidSppsMajor = find(potentialDisListMajor<mergingThreshold_);
				if ~isempty(potentialSolidSppsMajor)
					spps2BeMerged = sppsEmptyMajorValence(potentialSolidSppsMajor);							
					seedPointsRef_(spps2BeMerged,:) = potentialPosListMajor(potentialSolidSppsMajor,:);
					seedPoints_(spps2BeMerged,:) = potentialPosListMajorPara(potentialSolidSppsMajor,:);
					seedPointsValence_(spps2BeMerged,1) = 1;											
					modifiedMinorValences = HighCurvatureModification(spps2BeMerged, 'MINOR');
					seedPointsValence_(modifiedMinorValences,2) = 1;							
				end				
			end
		end
		
		if 0==valences(2)
			seedPointsValence_(spp,2) = 1;
			minorPSL = CreatePrincipalStressLine(iSeed, 'MINOR');
			if 0==minorPSL.numIntergPts
				looper = sum(seedPointsValence_(:)); 
				disp([' Iteration.: ' sprintf('%4i',its) ' Progress.: ' sprintf('%6i',looper) ...
					' Total.: ' sprintf('%6i',2*numSeedPoints)]); continue; 			
			end
			minorPSLpool_(end+1,1) = minorPSL;
			minorCoordList_(end+1:end+minorPSL.numIntergPts,:) = minorPSL.phyCoordList;
			sppsEmptyMinorValence = find(0==seedPointsValence_(:,2));
			if ~isempty(sppsEmptyMinorValence)
				[potentialDisListMinor, potentialPosListMinor, potentialPosListMinorPara] = GetDisListOfPointList2Curve(...	
					seedPointsRef_(sppsEmptyMinorValence,:), [minorPSL.eleIndexList minorPSL.phyCoordList], [minorPSL.eleIndexList minorPSL.paraCoordList]);
				potentialSolidSppsMinor = find(potentialDisListMinor<mergingThreshold_);
				if ~isempty(potentialSolidSppsMinor)
					spps2BeMerged = sppsEmptyMinorValence(potentialSolidSppsMinor);
					seedPointsRef_(spps2BeMerged,:) = potentialPosListMinor(potentialSolidSppsMinor,:);
					seedPoints_(spps2BeMerged,:) = potentialPosListMinorPara(potentialSolidSppsMinor,:);
					seedPointsValence_(spps2BeMerged,2) = 1;
					modifiedMajorValences = HighCurvatureModification(spps2BeMerged, 'MAJOR');
					seedPointsValence_(modifiedMajorValences,1) = 1;	
				end
			end
		end
		looper = sum(seedPointsValence_(:));
		disp([' Iteration.: ' sprintf('%4i',its) ' Progress.: ' sprintf('%6i',looper) ...
			' Total.: ' sprintf('%6i',2*numSeedPoints)]);	
		% ShowSeedingProcess(0.75);	
	end
	majorPSLpool_ = CompactStreamlines(majorPSLpool_, 5);
	minorPSLpool_ = CompactStreamlines(minorPSLpool_, 5);		
end


function ShowPSLs(stressCompBG, stressCompPSLs, lw, miniIntegralPts)
	%% Syntax:
	%% stressCompBG: %% 'None', 'Sigma_1', 'Sigma_2', 'Sigma_xx', 'Sigma_yy', 'Sigma_xy', 'Sigma_vM' ('Sigma_3', 'Sigma_zz', 'Sigma_xz', 'Sigma_yz', GLOBAL Stress Input) 
	%% stressCompPSLs: %% 'None', 'Sigma', 'Sigma_xx', 'Sigma_yy', 'Sigma_xy', 'Sigma_vM'
	%% ShowPSLs(stressCompBG, stressCompPSLs, lw, miniIntegralPts);
	global frameType_ numNodes_ silhouetteStruct_ boundingBox_ localStressField_ globalStressField_;
	global degePts_ majorPSLpool_ minorPSLpool_;

	%% Prepare PSLs for Vis
	tarIndice = [];
	for ii=1:numel(majorPSLpool_)
		if majorPSLpool_(ii).numIntergPts > miniIntegralPts, tarIndice(end+1,1) = ii; end
	end
	tarMajorPSLs = majorPSLpool_(tarIndice);	
	numTarMajorPSLs = numel(tarMajorPSLs);
	
	tarIndice = [];
	for ii=1:numel(minorPSLpool_)
		if minorPSLpool_(ii).numIntergPts > miniIntegralPts, tarIndice(end+1,1) = ii; end
	end
	tarMinorPSLs = minorPSLpool_(tarIndice);	
	numTarMinorPSLs = numel(tarMinorPSLs);
	numDege = numel(degePts_);
	degePtPos = zeros(numDege,3);
	if numDege>0
		for ii=1:numDege
			degePtPos(ii,:) = degePts_(ii).phyCoord;
			if 1==ii
				majorPSLs_TS = degePts_(ii).majorSkeletons;
				minorPSLs_TS = degePts_(ii).minorSkeletons;
			else
				majorPSLs_TS(end+1:end+numel(degePts_(ii).majorSkeletons),1) = degePts_(ii).majorSkeletons;
				minorPSLs_TS(end+1:end+numel(degePts_(ii).minorSkeletons),1) = degePts_(ii).minorSkeletons;
			end
		end
		tarIndiceTS = [];
		for ii=1:numel(majorPSLs_TS)
			if majorPSLs_TS(ii).numIntergPts > 2, tarIndiceTS(end+1,1) = ii; end
		end
		tarMajorPSLs_TS = majorPSLs_TS(tarIndiceTS);	
		numTarMajorPSLs_TS = numel(tarMajorPSLs_TS);
		tarIndiceTS = [];
		for ii=1:numel(minorPSLs_TS)
			if minorPSLs_TS(ii).numIntergPts > 2, tarIndiceTS(end+1,1) = ii; end
		end
		tarMinorPSLs_TS = minorPSLs_TS(tarIndiceTS);
		numTarMinorPSLs_TS = numel(tarMinorPSLs_TS);
	else
		tarMajorPSLs_TS = []; tarMinorPSLs_TS = [];
		numTarMajorPSLs_TS = 0; numTarMinorPSLs_TS = 0;
	end
	
	if 0==numTarMajorPSLs && 0==numTarMinorPSLs && 0==numTarMajorPSLs_TS && 0==numTarMinorPSLs_TS, return; end

	%% Setup Color Scheme
	switch stressCompBG
		case 'Sigma_1' %%Major Principal Stress
			switch frameType_
				case 'LOCAL'
					tmp = ComputePrincipalStressLocalFrame(localStressField_);
					cVal = tmp(:,4);					
				case 'GLOBAL'
					tmp = ComputePrincipalStressLocalFrame(globalStressField_);
					cVal = tmp(:,9);						
			end		
		case 'Sigma_2' %%Minor (Medium for global) Principal Stress
			switch frameType_
				case 'LOCAL'
					tmp = ComputePrincipalStressLocalFrame(localStressField_);
					cVal = tmp(:,1);				
				case 'GLOBAL'
					tmp = ComputePrincipalStressLocalFrame(globalStressField_);
					cVal = tmp(:,5);
			end		
		case 'Sigma_3' %%Minor Principal Stress
			if strcmp(frameType_, 'GLOBAL')
				tmp = ComputePrincipalStressLocalFrame(globalStressField_);
				cVal = tmp(:,1);
			else
				warning('Un-supported of the local stres field input!'); return;
			end
		case 'Sigma_vM' %%von Mises Stress
			switch frameType_
				case 'LOCAL'
					cVal = ComputeVonMisesStressLocalFrame(localStressField_);
				case 'GLOBAL'
					cVal = ComputeVonMisesStressLocalFrame(globalStressField_);
			end		
		case 'Sigma_xx'
			switch frameType_
				case 'LOCAL'
					cVal = localStressField_(:,1);
				case 'GLOBAL'
					cVal = globalStressField_(:,1);
			end		
		case 'Sigma_yy'
			switch frameType_
				case 'LOCAL'
					cVal = localStressField_(:,2);
				case 'GLOBAL'
					cVal = globalStressField_(:,2);
			end		
		case {'Sigma_xy', 'Sigma_yx'}
			switch frameType_
				case 'LOCAL'
					cVal = localStressField_(:,3);
				case 'GLOBAL'
					cVal = globalStressField_(:,6);
			end
		case 'Sigma_zz'
			if strcmp(frameType_, 'GLOBAL')
				cVal = globalStressField_(:,3);
			else
				warning('Un-supported of the local stres field input!'); return;
			end			
		case {'Sigma_xz', 'Sigma_zx'}
			if strcmp(frameType_, 'GLOBAL')
				cVal = globalStressField_(:,5);
			else
				warning('Un-supported of the local stres field input!'); return;
			end	
		case {'Sigma_yz', 'Sigma_zy'}
			if strcmp(frameType_, 'GLOBAL')
				cVal = globalStressField_(:,4);
			else
				warning('Un-supported of the local stres field input!'); return;
			end			
		case 'None'
			cVal = ones(numNodes_,1);			
		otherwise
			warning('Wrong Input!');
	end
	minValCol = min(cVal);
	
	if ~strcmp(stressCompBG, 'None')
		stressCompPSLs = 'None';
	end
	color4MajorPSLs = struct('arr', []); color4MajorPSLs = repmat(color4MajorPSLs, numTarMajorPSLs, 1);
	color4MinorPSLs = struct('arr', []); color4MinorPSLs = repmat(color4MinorPSLs, numTarMinorPSLs, 1);
	color4MajorPSLs_TS = struct('arr', []); color4MajorPSLs_TS = repmat(color4MajorPSLs_TS, numTarMajorPSLs_TS, 1);
	color4MinorPSLs_TS = struct('arr', []); color4MinorPSLs_TS = repmat(color4MinorPSLs_TS, numTarMinorPSLs_TS, 1);	
	switch stressCompPSLs
		case 'None'
			for ii=1:numTarMajorPSLs
				color4MajorPSLs(ii).arr = minValCol*ones(1, tarMajorPSLs(ii).numIntergPts);
			end		
			for ii=1:numTarMinorPSLs
				color4MinorPSLs(ii).arr = minValCol*ones(1, tarMinorPSLs(ii).numIntergPts);
			end
			for ii=1:numTarMajorPSLs_TS
				color4MajorPSLs_TS(ii).arr = minValCol*ones(1, tarMajorPSLs_TS(ii).numIntergPts);
			end		
			for ii=1:numTarMinorPSLs_TS
				color4MinorPSLs_TS(ii).arr = minValCol*ones(1, tarMinorPSLs_TS(ii).numIntergPts);
			end			
		case 'Sigma'
			for ii=1:numTarMajorPSLs
				color4MajorPSLs(ii).arr = tarMajorPSLs(ii).localPrincipalStressList(:,4)';
			end	
			for ii=1:numTarMinorPSLs
				color4MinorPSLs(ii).arr = tarMinorPSLs(ii).localPrincipalStressList(:,1)';
			end
			for ii=1:numTarMajorPSLs_TS
				color4MajorPSLs_TS(ii).arr = tarMajorPSLs_TS(ii).localPrincipalStressList(:,4)';
			end	
			for ii=1:numTarMinorPSLs_TS
				color4MinorPSLs_TS(ii).arr = tarMinorPSLs_TS(ii).localPrincipalStressList(:,1)';
			end
			cVal = ones(size(cVal)) * min([[color4MajorPSLs.arr], [color4MinorPSLs.arr]]);
		case 'Sigma_xx'
			for ii=1:numTarMajorPSLs
				color4MajorPSLs(ii).arr = tarMajorPSLs(ii).localStressTensorList(:,1)';
			end
			for ii=1:numTarMinorPSLs
				color4MinorPSLs(ii).arr = tarMinorPSLs(ii).localStressTensorList(:,1)';
			end
			for ii=1:numTarMajorPSLs_TS
				color4MajorPSLs_TS(ii).arr = tarMajorPSLs_TS(ii).localStressTensorList(:,1)';
			end
			for ii=1:numTarMinorPSLs_TS
				color4MinorPSLs_TS(ii).arr = tarMinorPSLs_TS(ii).localStressTensorList(:,1)';
			end
			cVal = ones(size(cVal)) * min([[color4MajorPSLs.arr], [color4MinorPSLs.arr]]);
		case 'Sigma_yy'
			for ii=1:numTarMajorPSLs
				color4MajorPSLs(ii).arr = tarMajorPSLs(ii).localStressTensorList(:,2)';
			end	
			for ii=1:numTarMinorPSLs
				color4MinorPSLs(ii).arr = tarMinorPSLs(ii).localStressTensorList(:,2)';
			end
			for ii=1:numTarMajorPSLs_TS
				color4MajorPSLs_TS(ii).arr = tarMajorPSLs_TS(ii).localStressTensorList(:,2)';
			end	
			for ii=1:numTarMinorPSLs_TS
				color4MinorPSLs_TS(ii).arr = tarMinorPSLs_TS(ii).localStressTensorList(:,2)';
			end
			cVal = ones(size(cVal)) * min([[color4MajorPSLs.arr], [color4MinorPSLs.arr]]);
		case {'Sigma_xy', 'Sigma_yx'}
			for ii=1:numTarMajorPSLs
				color4MajorPSLs(ii).arr = tarMajorPSLs(ii).localStressTensorList(:,3)';
			end	
			for ii=1:numTarMinorPSLs
				color4MinorPSLs(ii).arr = tarMinorPSLs(ii).localStressTensorList(:,3)';
			end
			for ii=1:numTarMajorPSLs_TS
				color4MajorPSLs_TS(ii).arr = tarMajorPSLs_TS(ii).localStressTensorList(:,3)';
			end	
			for ii=1:numTarMinorPSLs_TS
				color4MinorPSLs_TS(ii).arr = tarMinorPSLs_TS(ii).localStressTensorList(:,3)';
			end	
			cVal = ones(size(cVal)) * min([[color4MajorPSLs.arr], [color4MinorPSLs.arr]]);
		case 'Sigma_vM'
			for ii=1:numTarMajorPSLs
				color4MajorPSLs(ii).arr = tarMajorPSLs(ii).localVonMisesStressList';
			end			
			for ii=1:numTarMinorPSLs
				color4MinorPSLs(ii).arr = tarMinorPSLs(ii).localVonMisesStressList';
			end
			for ii=1:numTarMajorPSLs_TS
				color4MajorPSLs_TS(ii).arr = tarMajorPSLs_TS(ii).localVonMisesStressList';
			end			
			for ii=1:numTarMinorPSLs_TS
				color4MinorPSLs_TS(ii).arr = tarMinorPSLs_TS(ii).localVonMisesStressList';
			end
			cVal = ones(size(cVal)) * min([[color4MajorPSLs.arr], [color4MinorPSLs.arr]]);
		otherwise
			error('Wrong Input!');			
	end	
	
	
	%%Draw
	figure; 
	handleSilhouette = patch(silhouetteStruct_, 'FaceVertexCData', cVal); hold('on');

	lineWidthTube = lw*min(boundingBox_(2,:)-boundingBox_(1,:))/100;
	if 0==lineWidthTube %%2D
		majorPSLsGraph.vertices = [];
		majorPSLsGraph.faces = [];
		majorPSLsGraphColor = [];
		for jj=1:numel(tarMajorPSLs)
			if 1==jj
				majorPSLsGraph.vertices = tarMajorPSLs(jj).phyCoordList;
				majorPSLsGraph.faces = [1:tarMajorPSLs(jj).numIntergPts-1; 2:tarMajorPSLs(jj).numIntergPts]';
				majorPSLsGraphColor = color4MajorPSLs(jj).arr(:);
			else
				existingVertices = size(majorPSLsGraph.vertices,1);
				majorPSLsGraph.vertices(end+1:end+tarMajorPSLs(jj).numIntergPts,:) = tarMajorPSLs(jj).phyCoordList;
				majorPSLsGraph.faces(end+1:end+tarMajorPSLs(jj).numIntergPts-1,:) = existingVertices + [1:tarMajorPSLs(jj).numIntergPts-1; 2:tarMajorPSLs(jj).numIntergPts]';
				majorPSLsGraphColor(end+1:end+tarMajorPSLs(jj).numIntergPts,:) = color4MajorPSLs(jj).arr(:);
			end
		end
		if numTarMajorPSLs_TS > 0
			majorPSLsGraph_TS.vertices = [];
			majorPSLsGraph_TS.faces = [];
			majorPSLsGraphColor_TS = [];
			for jj=1:numel(tarMajorPSLs_TS)
				if 1==jj
					majorPSLsGraph_TS.vertices = tarMajorPSLs_TS(jj).phyCoordList;
					majorPSLsGraph_TS.faces = [1:tarMajorPSLs_TS(jj).numIntergPts-1; 2:tarMajorPSLs_TS(jj).numIntergPts]';
					majorPSLsGraphColor_TS = color4MajorPSLs_TS(jj).arr(:);
				else
					existingVertices = size(majorPSLsGraph_TS.vertices,1);
					majorPSLsGraph_TS.vertices(end+1:end+tarMajorPSLs_TS(jj).numIntergPts,:) = tarMajorPSLs_TS(jj).phyCoordList;
					majorPSLsGraph_TS.faces(end+1:end+tarMajorPSLs_TS(jj).numIntergPts-1,:) = existingVertices + [1:tarMajorPSLs_TS(jj).numIntergPts-1; 2:tarMajorPSLs_TS(jj).numIntergPts]';
					majorPSLsGraphColor_TS(end+1:end+tarMajorPSLs_TS(jj).numIntergPts,:) = color4MajorPSLs_TS(jj).arr(:);
				end
			end			
		end
		minorPSLsGraph.vertices = [];
		minorPSLsGraph.faces = [];
		minorPSLsGraphColor = [];
		for jj=1:numel(tarMinorPSLs)
			if 1==jj
				minorPSLsGraph.vertices = tarMinorPSLs(jj).phyCoordList;
				minorPSLsGraph.faces = [1:tarMinorPSLs(jj).numIntergPts-1; 2:tarMinorPSLs(jj).numIntergPts]';
				minorPSLsGraphColor = color4MinorPSLs(jj).arr(:);
			else
				existingVertices = size(minorPSLsGraph.vertices,1);
				minorPSLsGraph.vertices(end+1:end+tarMinorPSLs(jj).numIntergPts,:) = tarMinorPSLs(jj).phyCoordList;
				minorPSLsGraph.faces(end+1:end+tarMinorPSLs(jj).numIntergPts-1,:) = existingVertices + [1:tarMinorPSLs(jj).numIntergPts-1; 2:tarMinorPSLs(jj).numIntergPts]';
				minorPSLsGraphColor(end+1:end+tarMinorPSLs(jj).numIntergPts,:) = color4MinorPSLs(jj).arr(:);
			end
		end
		if numTarMinorPSLs_TS > 0
			minorPSLsGraph_TS.vertices = [];
			minorPSLsGraph_TS.faces = [];
			minorPSLsGraphColor_TS = [];
			for jj=1:numel(tarMinorPSLs_TS)
				if 1==jj
					minorPSLsGraph_TS.vertices = tarMinorPSLs_TS(jj).phyCoordList;
					minorPSLsGraph_TS.faces = [1:tarMinorPSLs_TS(jj).numIntergPts-1; 2:tarMinorPSLs_TS(jj).numIntergPts]';
					minorPSLsGraphColor_TS = color4MinorPSLs_TS(jj).arr(:);
				else
					existingVertices = size(minorPSLsGraph_TS.vertices,1);
					minorPSLsGraph_TS.vertices(end+1:end+tarMinorPSLs_TS(jj).numIntergPts,:) = tarMinorPSLs_TS(jj).phyCoordList;
					minorPSLsGraph_TS.faces(end+1:end+tarMinorPSLs_TS(jj).numIntergPts-1,:) = existingVertices + [1:tarMinorPSLs_TS(jj).numIntergPts-1; 2:tarMinorPSLs_TS(jj).numIntergPts]';
					minorPSLsGraphColor_TS(end+1:end+tarMinorPSLs_TS(jj).numIntergPts,:) = color4MinorPSLs_TS(jj).arr(:);
				end
			end		
		end
		handleMajorPSL = []; 
		if ~isempty(majorPSLsGraphColor)
			hold('on');
			handleMajorPSL = patch(majorPSLsGraph, 'FaceVertexCData', majorPSLsGraphColor);
		end
		handleMinorPSL = []; 
		if ~isempty(minorPSLsGraphColor)
			hold('on');
			handleMinorPSL = patch(minorPSLsGraph, 'FaceVertexCData', minorPSLsGraphColor);
		end
		handleMajorPSL_TS = [];
		if ~isempty(majorPSLsGraphColor_TS)
			hold('on');
			handleMajorPSL_TS = patch(majorPSLsGraph_TS, 'FaceVertexCData', majorPSLsGraphColor_TS);
		end
		handleMinorPSL_TS = []; 
		if ~isempty(minorPSLsGraphColor_TS)
			hold('on');
			handleMinorPSL_TS = patch(minorPSLsGraph_TS, 'FaceVertexCData', minorPSLsGraphColor_TS);
		end
		handleDegPot = [];
		if numDege > 0
			handleDegPot = plot3(degePtPos(:,1), degePtPos(:,2), degePtPos(:,3), 'o', 'Color', [112 48 160]/255, 'lineWidth', 5, 'MarkerSize', 20);
		end

		if strcmp(stressCompBG, 'None')
			if strcmp(stressCompPSLs, 'None')
				set(handleMajorPSL, 'FaceColor', 'None', 'EdgeColor', [252 141 98]/255, 'lineWidth', 3);
				set(handleMinorPSL, 'FaceColor', 'None', 'EdgeColor', [102 194 165]/255, 'lineWidth', 3); %[102 194 165]/255
				set(handleMajorPSL_TS, 'FaceColor', 'None', 'EdgeColor', [252 141 98]/255, 'lineWidth', 6);
				set(handleMinorPSL_TS, 'FaceColor', 'None', 'EdgeColor', [102 194 165]/255, 'lineWidth', 6);					
			else
				colormap('parula'); 
				%%Colorbar
				cb = colorbar('Location', 'east', 'AxisLocation','in');
				t=get(cb,'Limits'); 
				set(cb,'Ticks',linspace(t(1),t(2),5));
				L=cellfun(@(x)sprintf('%.2e',x),num2cell(linspace(t(1),t(2),5)),'Un',0); 
				set(cb,'xticklabel',L);
				set(handleMajorPSL, 'EdgeColor', 'interp', 'lineWidth', 3);
				set(handleMinorPSL, 'EdgeColor', 'interp', 'lineWidth', 3);
				set(handleMajorPSL_TS, 'EdgeColor', 'interp', 'lineWidth', 6);
				set(handleMinorPSL_TS, 'EdgeColor', 'interp', 'lineWidth', 6);					
			end
			set(handleSilhouette, 'faceColor', [0.9 0.9 0.9], 'FaceAlpha', 1.0, 'EdgeColor', 'None');
		else
			colormap('parula'); 
			%%Colorbar
			cb = colorbar('Location', 'east', 'AxisLocation','in');
			t=get(cb,'Limits'); 
			set(cb,'Ticks',linspace(t(1),t(2),5));
			L=cellfun(@(x)sprintf('%.2e',x),num2cell(linspace(t(1),t(2),5)),'Un',0); 
			set(cb,'xticklabel',L);			
			set(handleSilhouette, 'FaceColor', 'interp', 'FaceAlpha', 1.0, 'EdgeColor', 'none');	
			set(handleMajorPSL, 'FaceColor', 'None', 'EdgeColor', [252 141 98]/255, 'lineWidth', 3);
			set(handleMinorPSL, 'FaceColor', 'None', 'EdgeColor', [102 194 165]/255, 'lineWidth', 3);
			set(handleMajorPSL_TS, 'FaceColor', 'None', 'EdgeColor', [252 141 98]/255, 'lineWidth', 6);
			set(handleMinorPSL_TS, 'FaceColor', 'None', 'EdgeColor', [102 194 165]/255, 'lineWidth', 6);	
		end		
	else
		handleMajorPSL = []; 
		[gridX, gridY, gridZ, gridC, ~] = ExpandPSLs2Tubes(tarMajorPSLs, color4MajorPSLs, lineWidthTube);
		if ~isempty(gridX)
			hold('on'); 
			handleMajorPSL = surf(gridX, gridY, gridZ, gridC);
			shading('interp');
		end
		handleMajorPSL_TS = []; 
		[gridX, gridY, gridZ, gridC, ~] = ExpandPSLs2Tubes(tarMajorPSLs_TS, color4MajorPSLs_TS, 1.5*lineWidthTube);
		if ~isempty(gridX)
			hold('on'); 
			handleMajorPSL_TS = surf(gridX, gridY, gridZ, gridC);
			shading('interp');
		end
		handleMinorPSL = []; 
		[gridX, gridY, gridZ, gridC, ~] = ExpandPSLs2Tubes(tarMinorPSLs, color4MinorPSLs, lineWidthTube);
		if ~isempty(gridX)
			hold('on');
			handleMinorPSL = surf(gridX, gridY, gridZ, gridC);
			shading('interp');
		end
		handleMinorPSL_TS = []; 
		[gridX, gridY, gridZ, gridC, ~] = ExpandPSLs2Tubes(tarMinorPSLs_TS, color4MinorPSLs_TS, 1.5*lineWidthTube);
		if ~isempty(gridX)
			hold('on');
			handleMinorPSL_TS = surf(gridX, gridY, gridZ, gridC);
			shading('interp');
		end
		handleDegPot = [];
		if numDege > 0
			[patchX, patchY, patchZ] = GetSpheresOfPoints(degePtPos, 2.5*lineWidthTube);
			handleDegPot = surf(patchX, patchY, patchZ, min([[color4MajorPSLs.arr], [color4MinorPSLs.arr]])*ones(size(patchZ)));
		end
		if strcmp(stressCompBG, 'None')
			if strcmp(stressCompPSLs, 'None')
				set(handleMajorPSL, 'FaceColor', [252 141 98]/255, 'FaceAlpha', 1, 'EdgeAlpha', 0);
				set(handleMinorPSL, 'FaceColor', [102 194 165]/255, 'FaceAlpha', 1, 'EdgeAlpha', 0);
				set(handleMajorPSL_TS, 'FaceColor', [252 141 98]/255, 'FaceAlpha', 1, 'EdgeAlpha', 0);
				set(handleMinorPSL_TS, 'FaceColor', [102 194 165]/255, 'FaceAlpha', 1, 'EdgeAlpha', 0);							
			else
				colormap('parula'); 
				%%Colorbar
				cb = colorbar('Location', 'east', 'AxisLocation','in');
				t=get(cb,'Limits'); 
				set(cb,'Ticks',linspace(t(1),t(2),5));
				L=cellfun(@(x)sprintf('%.2e',x),num2cell(linspace(t(1),t(2),5)),'Un',0); 
				set(cb,'xticklabel',L);							
			end
			set(handleSilhouette, 'faceColor', [0.9 0.9 0.9], 'FaceAlpha', 1.0, 'EdgeColor', 'None');
			set(handleDegPot, 'FaceColor', [112 48 160]/255, 'FaceAlpha', 1, 'EdgeAlpha', 0);
		else
			colormap('parula'); 
			%%Colorbar
			cb = colorbar('Location', 'east', 'AxisLocation','in');
			t=get(cb,'Limits'); 
			set(cb,'Ticks',linspace(t(1),t(2),5));
			L=cellfun(@(x)sprintf('%.2e',x),num2cell(linspace(t(1),t(2),5)),'Un',0); 
			set(cb,'xticklabel',L);			
			set(handleSilhouette, 'FaceColor', 'interp', 'FaceAlpha', 1.0, 'EdgeColor', 'none');	
			
			set(handleMajorPSL, 'FaceColor', [252 141 98]/255, 'FaceAlpha', 1, 'EdgeAlpha', 0);
			set(handleMinorPSL, 'FaceColor', [102 194 165]/255, 'FaceAlpha', 1, 'EdgeAlpha', 0);
			set(handleMajorPSL_TS, 'FaceColor', [252 141 98]/255, 'FaceAlpha', 1, 'EdgeAlpha', 0);
			set(handleMinorPSL_TS, 'FaceColor', [102 194 165]/255, 'FaceAlpha', 1, 'EdgeAlpha', 0);	
			set(handleDegPot, 'FaceColor', [112 48 160]/255, 'FaceAlpha', 1, 'EdgeAlpha', 0);	%%[112 48 160]/255	
		end
	end
	
	% %%Lighting, Reflection
	if 0==lineWidthTube
		view(2);
		xlabel('X'); ylabel('Y');
	else
		view(3);
		lighting('gouraud');
		material('dull');
		camlight('headlight', 'infinite');
		xlabel('X'); ylabel('Y'); zlabel('Z');
	end
	axis('equal', 'tight', 'off');
	set(gca, 'FontName', 'Times New Roman', 'FontSize', 30);
end


function oPSLs = CompactStreamlines(iPSLs, truncatedThreshold)
	tarIndice = [];
	for ii=1:numel(iPSLs)
		if iPSLs(ii).numIntergPts > truncatedThreshold
			tarIndice(end+1,1) = ii;
		end
	end
	oPSLs = iPSLs(tarIndice);
	if isempty(oPSLs), oPSLs = []; end
end

function iniSeedPoints = SetupSeeds(seedSparsityCtrl)
	global numEles_ voidEleMap_;
	global meshTypeMap_;
	iniSeedPoints = zeros(numEles_, 3);
	for ii=1:numEles_
		iEleType = meshTypeMap_(ii);
		switch iEleType
			case {'T3', 'T6'}
				iniSeedPoints(ii,:) = [ii 1/3 1/3];
			case {'Q4', 'Q8'}
				iniSeedPoints(ii,:) = [ii 0 0];
		end
	end
	iniSeedPoints = iniSeedPoints(1:seedSparsityCtrl:end,:);
	iniSeedPoints(find(voidEleMap_),:) = [];
end

function EmbedTopologicalSkeleton()
	global seedPoints_;
	global seedPointsValence_;
	global seedPointsRef_;
	global majorPSLpool_ minorPSLpool_
	global majorCoordList_ minorCoordList_
	global mergingThreshold_;
	global degePts_;
	
	%%Move Topological Skeleton to Major and Minor PSL Pools
	if ~isempty(degePts_)
		for ii=1:numel(degePts_)
			for jj=1:numel(degePts_(ii).majorSkeletons)
				majorPSLpool_(end+1,1) = degePts_(ii).majorSkeletons(jj);
				majorCoordList_(end+1:end+degePts_(ii).majorSkeletons(jj).numIntergPts,:) = degePts_(ii).majorSkeletons(jj).phyCoordList;
			end
			for jj=1:numel(degePts_(ii).minorSkeletons)
				minorPSLpool_(end+1,1) = degePts_(ii).minorSkeletons(jj);
				minorCoordList_(end+1:end+degePts_(ii).minorSkeletons(jj).numIntergPts,:) = degePts_(ii).minorSkeletons(jj).phyCoordList;
			end			
		end
	end
	numMajorPSLs = numel(majorPSLpool_);
	for ii=1:numMajorPSLs
		majorPSL = majorPSLpool_(ii);
		if majorPSL.numIntergPts>0					
			sppsEmptyMajorValence = find(0==seedPointsValence_(:,1));
            if ~isempty(sppsEmptyMajorValence)
				[potentialDisListMajor, potentialPosListMajor, potentialPosListMajorPara] = GetDisListOfPointList2Curve(...	
					seedPointsRef_(sppsEmptyMajorValence,:), [majorPSL.eleIndexList majorPSL.phyCoordList], [majorPSL.eleIndexList majorPSL.paraCoordList]);
				potentialSolidSppsMajor = find(potentialDisListMajor<=mergingThreshold_);
				if ~isempty(potentialSolidSppsMajor)
					spps2BeMerged = sppsEmptyMajorValence(potentialSolidSppsMajor);							
					seedPointsRef_(spps2BeMerged,:) = potentialPosListMajor(potentialSolidSppsMajor,:);
					seedPoints_(spps2BeMerged,:) = potentialPosListMajorPara(potentialSolidSppsMajor,:);
					seedPointsValence_(spps2BeMerged,1) = 1;											
					modifiedMinorValences = HighCurvatureModification(spps2BeMerged, 'MINOR');
					seedPointsValence_(modifiedMinorValences,2) = 1;							
				end
			end
		end
	end
	numMinorPSLs = numel(minorPSLpool_);
	for ii=1:numMinorPSLs
		minorPSL = minorPSLpool_(ii);
		if minorPSL.numIntergPts>0	
			sppsEmptyMinorValence = find(0==seedPointsValence_(:,2));
            if ~isempty(sppsEmptyMinorValence)
				[potentialDisListMinor, potentialPosListMinor, potentialPosListMinorPara] = GetDisListOfPointList2Curve(...	
					seedPointsRef_(sppsEmptyMinorValence,:), [minorPSL.eleIndexList minorPSL.phyCoordList], [minorPSL.eleIndexList minorPSL.paraCoordList]);
				potentialSolidSppsMinor = find(potentialDisListMinor<=mergingThreshold_);
				if ~isempty(potentialSolidSppsMinor)
					spps2BeMerged = sppsEmptyMinorValence(potentialSolidSppsMinor);
					seedPointsRef_(spps2BeMerged,:) = potentialPosListMinor(potentialSolidSppsMinor,:);
					seedPoints_(spps2BeMerged,:) = potentialPosListMinorPara(potentialSolidSppsMinor,:);
					seedPointsValence_(spps2BeMerged,2) = 1;
					modifiedMajorValences = HighCurvatureModification(spps2BeMerged, 'MAJOR');
					seedPointsValence_(modifiedMajorValences,1) = 1;						
				end
			end
		end
	end
end

function [potentialDisList, potentialPosList, potentialPosListPara] = GetDisListOfPointList2Curve(srcPoints, tarCurve, tarCurvePara)
	disT = DistanceMap2PtCloud(srcPoints(:,end-2:end), tarCurve(:,end-2:end));
	[minVal, minValPos] = min(disT,[],1);
	potentialPosList = tarCurve(minValPos,:);
	potentialPosListPara = tarCurvePara(minValPos,:);
	potentialDisList = minVal(:);
end

function modifiedValences = HighCurvatureModification(spps2BeMerged, psDir)
	global majorCoordList_; 
	global minorCoordList_;		
	global seedPointsRef_;
	global seedPointsValence_;
	global mergingThreshold_;

	tarCurve = [];
	switch psDir
		case 'MAJOR'
			if isempty(majorCoordList_), modifiedValences = []; return; end
			tarCurve = majorCoordList_;
            spps2BeMerged = spps2BeMerged(0==seedPointsValence_(spps2BeMerged,1));
		case 'MINOR'
			if isempty(minorCoordList_), modifiedValences = []; return; end
			tarCurve = minorCoordList_;
            spps2BeMerged = spps2BeMerged(0==seedPointsValence_(spps2BeMerged,2));
	end
	srcPoints = seedPointsRef_(spps2BeMerged,end-2:end);
	disT = DistanceMap2PtCloud(srcPoints, tarCurve);	
	minVal = min(disT, [], 1);
	modifiedValences = spps2BeMerged(find(minVal<mergingThreshold_));
end

function parasOut = AvoidExactLocationAtVertexOrEdge4Robustness(parasIn, eleType, safetyScaling)
	xi = parasIn(1);
	eta = parasIn(2);
	tolEff = safetyScaling;
	if strcmp(eleType, 'T3') || strcmp(eleType, 'T6')
		if abs(xi)<tolEff
			xi = safetyScaling;
		elseif abs(xi-1)<tolEff
			xi = 1 - safetyScaling;
		end
		if abs(eta)<tolEff
			eta = safetyScaling;
		elseif abs(eta-1)<tolEff
			eta = 1 - safetyScaling;
		end
		if abs(1-xi-eta)<tolEff
			xi = xi*(1-safetyScaling);
			eta = eta*(1-safetyScaling);
		end			
	elseif strcmp(eleType, 'Q4') || strcmp(eleType, 'Q8')
		if xi>-1-tolEff &&  xi<-1+tolEff
			xi = -1+safetyScaling;
		elseif xi>1-tolEff && xi<1+tolEff
			xi = 1-safetyScaling;
		end
		if eta>-1-tolEff &&  eta<-1+tolEff
			eta = -1+safetyScaling;
		elseif eta>1-tolEff && eta<1+tolEff
			eta = 1-safetyScaling;
		end			
	else
		error('Un-supported Element Type!');
	end
	parasOut = [xi eta];
end

function inside = CheckWhetherGivenNaturalCoordinatesWithinDomain(paras, eleType, varargin)
	inside  = 0;
	xi = paras(1); eta = paras(2);
	if ~(isfinite(xi) && isfinite(eta))
		inside = false;
		error('Infinite Natural Coordinates!');
	end
	
	if 3==nargin
		tolEff = varargin{1};
	else
		tolEff = 0;
	end
	%tolEff = 50 * eps(scale);
	
	if strcmp(eleType, 'T3') || strcmp(eleType, 'T6')
		l1 = 1 - xi - eta;
		l2 = xi;
		l3 = eta;
		inside = (l1 >= -tolEff) & (l2 >= -tolEff) & (l3 >= -tolEff);	
	elseif strcmp(eleType, 'Q4') || strcmp(eleType, 'Q8')
		inside = (xi >= -1 - tolEff) && (xi <= 1 + tolEff) && (eta >= -1 - tolEff) && (eta <= 1 + tolEff);
	else
		error('Un-supported Element Type!');
	end
end

function [R, origin, t1, t2] = ComputeLocalFrameAtGivenPosition(paras, iEleNodes, elementType)
	global refVec_ refVecFallback_ tolRefVecFallback_;
	N = ShapeFunction(paras, elementType);
	[dNdxi, dNdeta] = DeShapeFunction(paras, elementType);
	%%Tangent Planes at "paras"
	t1 = dNdxi * iEleNodes;
	t2 = dNdeta * iEleNodes;
	origin = N*iEleNodes;
	nVec = cross(t1, t2);
	e3 = nVec ./ vecnorm(nVec,2,2);
	rVec = refVec_;
	tVec = rVec - (dot(rVec, e3))*e3;
	if norm(tVec) < tolRefVecFallback_
		rVec = refVecFallback_;
		tVec = rVec - (dot(rVec, e3))*e3;
		if norm(tVec) < tolRefVecFallback_
			error('Reference directions are parallel to normal; provide a custom ref.');
		end				
	end
	e1 = tVec / norm(tVec);
	e2 = cross(e3, e1);
	R = [e1(:), e2(:), e3(:)];
end

function [faceNormList, tVec1, tVec2] = ComputeNormalsAtGivenPosition(paras, iEleNodes, elementType)
	[dNdxi, dNdeta] = DeShapeFunction(paras, elementType);
	%%Tangent Planes at "paras"
	tVec1 = dNdxi * iEleNodes;
	tVec2 = dNdeta * iEleNodes;
	nVec = cross(tVec1, tVec2);
	nVec = nVec ./ vecnorm(nVec,2,2);
	faceNormList = nVec;
end

function psList = ComputePrincipalStressLocalFrame(stressTensors)
	numStressTensors = size(stressTensors,1);
	switch size(stressTensors,2)
		case 3 %%Local Frame
			%% Per-row: sigma_xx, sigma_yy, sigma_xy
			psList = zeros(numStressTensors,6);
			for ii=1:numStressTensors
				iStressTensor = stressTensors(ii,:);
				iStressTensor = iStressTensor([1 3; 3 2]);
				[eigenVec, eigenVal] = eig(iStressTensor);
				iPS = [eigenVal(1,1); eigenVec(:,1); eigenVal(2,2); eigenVec(:,2)]';
				psList(ii,:) = iPS;		
			end		
		case 6
			%% Per-row: sigma_xx, sigma_yy, sigma_zz, sigma_yz, sigma_zx, sigma_xy
			psList = zeros(numStressTensors,12);
			for ii=1:numStressTensors
				iStressTensor = stressTensors(ii,:);
				iStressTensor = iStressTensor([1 6 5; 6 2 4; 5 4 3]);
				[eigenVec, eigenVal] = eig(iStressTensor);
				iPS = [eigenVal(1,1); eigenVec(:,1); eigenVal(2,2); eigenVec(:,2); eigenVal(3,3); eigenVec(:,3)]';
				psList(ii,:) = iPS;		
			end		
	end
end

function iVonMisesStress = ComputeVonMisesStressLocalFrame(iStressTensor)
	switch size(iStressTensor,2)
		case 3 %%Local Frame
			%% Per-row: sigma_xx, sigma_yy, sigma_xy
			iVonMisesStress = sqrt(iStressTensor(:,1).^2 + iStressTensor(:,2).^2 - iStressTensor(:,1).*iStressTensor(:,2) + 3*iStressTensor(:,3).^2);	
		case 6
			%% Per-row: sigma_xx, sigma_yy, sigma_zz, sigma_yz, sigma_zx, sigma_xy
			iVonMisesStress = sqrt(0.5*((iStressTensor(:,1)-iStressTensor(:,2)).^2 + (iStressTensor(:,2)-iStressTensor(:,3)).^2 + (iStressTensor(:,3)...
					-iStressTensor(:,1)).^2 ) + 3*( iStressTensor(:,6).^2 + iStressTensor(:,4).^2 + iStressTensor(:,5).^2 ));			
	end
end

function e1paras = ConvertNaturalCoordsAtSharedEdge(element0, e0paras, edgeIdLocal0, type0, vtxEdgeGloal0, element, edgeIdLocal1, type1, vtxEdgeGloal1)
    xi0  = e0paras(1);  eta0 = e0paras(2);
    s0 = EdgeProgress(type0, edgeIdLocal0, xi0, eta0);
    if  isequal(vtxEdgeGloal0, vtxEdgeGloal1)
        s1 = s0;            % same local direction
    elseif isequal(vtxEdgeGloal0, fliplr(vtxEdgeGloal1))
        s1 = 1 - s0;        % opposite direction (typical for interior edges)
    else
        error('Shared edge mismatch: different end-node pair.');
    end
    [xi1, eta1] = NatFromProgress(type1, edgeIdLocal1, s1);
    e1paras = [xi1, eta1];
end

function s = EdgeProgress(t, edgeId, xi, eta)
	if 'T3'==t || 'T6'==t
        switch edgeId
            case 1, s = xi;          % along base 1->2
            case 2, s = eta;         % along hypotenuse 2->3
            case 3, s = 1 - eta;     % along left 3->1
            otherwise, error('tri edgeId must be 1..3');
        end	
	elseif 'Q4'==t || 'Q8'==t
       switch edgeId
           case 1, s = (xi  + 1)/2; % bottom
           case 2, s = (eta + 1)/2; % right
           case 3, s = (1 - xi)/2;  % top
           case 4, s = (1 - eta)/2; % left
           otherwise, error('quad edgeId must be 1..4');
       end
	else
		error('Un-supported Element Type!');
	end
end

function [xi,eta] = NatFromProgress(t, edgeId, s)
    s = max(0, min(1, s));  % tiny clamp
	if 'T3'==t || 'T6'==t
        switch edgeId
            case 1, xi = s;     eta = 0;
            case 2, xi = 1 - s; eta = s;
            case 3, xi = 0;     eta = 1 - s;
        end	
	elseif 'Q4'==t || 'Q8'==t
        switch edgeId
            case 1, xi = -1 + 2*s; eta = -1;
            case 2, xi =  1;      eta = -1 + 2*s;
            case 3, xi =  1 - 2*s;eta =  1;
            case 4, xi = -1;      eta =  1 - 2*s;
        end	
	else
		error('Un-supported Element Type!');
	end
end

function iPSL = CreatePrincipalStressLine(seed, tracingType, varargin)
	iPSL = PrincipalStressLineStruct();
	switch tracingType
		case 'MAJOR', psDir = [5 6];
		case 'MINOR', psDir = [2 3];
	end
	%%1. prepare for tracing
	[seed, phyCoord0, ~, localFrame0, t1, t2, localStressTensor0, localPrincipalStress0, localVonMisesStress0, opt] = PreparingForTracing(seed);
	if ~opt, return; end
	if 3==nargin
		localPrincipalStress0(psDir) = varargin{1}; %%Prescribed Principal Stress Direction at Degenerate Point
	end
	
	%%2. tracing PSL
	PSLphyCoordList = phyCoord0;
	PSLeleIndexList = seed(1);
	PSLparaCoordList = seed(2:3);
	PSLocalFrameList = localFrame0;
	PSLocalStressTensorList = localStressTensor0;
	PSLocalPrincipalStressList = localPrincipalStress0;
	PSLocalVonMisesStressList = localVonMisesStress0;
	
	e0 = PSLeleIndexList; para0 = PSLparaCoordList;
	v0 = localPrincipalStress0(1,psDir);	
	%%2.1 along first direction (v0)
	[phyCoordList, eleIndexList, paraCoordList, localFrameList, localStressTensorList, localPrincipalStressList, localVonMisesStressList] = ...
		Tracing(e0, para0, localFrame0, t1, t2, v0, psDir);
	numNewlyGeneratedIntergPts = numel(eleIndexList);
	PSLphyCoordList = [PSLphyCoordList; phyCoordList];
	PSLeleIndexList = [PSLeleIndexList; eleIndexList];
	PSLparaCoordList = [PSLparaCoordList; paraCoordList];
	PSLocalFrameList(:,:,end+1:end+numNewlyGeneratedIntergPts) = localFrameList;
	PSLocalStressTensorList = [PSLocalStressTensorList; localStressTensorList];
	PSLocalPrincipalStressList = [PSLocalPrincipalStressList; localPrincipalStressList];	
	PSLocalVonMisesStressList = [PSLocalVonMisesStressList; localVonMisesStressList];
	%%2.2 along second direction (-v0)
	[phyCoordList, eleIndexList, paraCoordList, localFrameList, localStressTensorList, localPrincipalStressList, localVonMisesStressList] = ...
		Tracing(e0, para0, localFrame0, t1, t2, -v0, psDir);
	numNewlyGeneratedIntergPts = numel(eleIndexList);
	if numNewlyGeneratedIntergPts > 1
		phyCoordList = flip(phyCoordList);
		eleIndexList = flip(eleIndexList);
		paraCoordList = flip(paraCoordList);
		localFrameList = flip(localFrameList,3);
		localStressTensorList = flip(localStressTensorList);
		localPrincipalStressList = flip(localPrincipalStressList);
		localVonMisesStressList = flip(localVonMisesStressList);
	end
	PSLphyCoordList = [phyCoordList; PSLphyCoordList];
	PSLeleIndexList = [eleIndexList; PSLeleIndexList];
	PSLparaCoordList = [paraCoordList; PSLparaCoordList];
	localFrameList(:,:,end+1:end+size(PSLocalFrameList,3)) = PSLocalFrameList; PSLocalFrameList = localFrameList;
	PSLocalStressTensorList = [localStressTensorList; PSLocalStressTensorList];
	PSLocalPrincipalStressList = [localPrincipalStressList; PSLocalPrincipalStressList];	
	PSLocalVonMisesStressList = [localVonMisesStressList; PSLocalVonMisesStressList];
	
	%%2.3 finish Tracing the current major PSL	
	iPSL.midPointPosition = size(phyCoordList,1)+1;
	iPSL.numIntergPts = size(PSLphyCoordList,1);
	iPSL.phyCoordList = PSLphyCoordList;
	iPSL.eleIndexList = PSLeleIndexList;
	iPSL.paraCoordList = PSLparaCoordList;
	iPSL.localFrames = PSLocalFrameList;
	iPSL.localStressTensorList = PSLocalStressTensorList;	
	iPSL.localPrincipalStressList = PSLocalPrincipalStressList;	
	iPSL.localVonMisesStressList = PSLocalVonMisesStressList;	
end

function [phyCoordList, eleIndexList, paraCoordList, localFrameList, localStressTensorList, localPrincipalStressList, ...
	localVonMisesStressList] = Tracing(e0, paraCoord0, localFrame0, t1, t2, v0, psDir)
	global tracingSche_ scalingStepSize4PSL_ limitSteps_;
	global meshTypeMap_ voidEleMap_ eleStruct_;
	global eleCharacterSizeList_ shellCurvatureScaling_;

	phyCoordList = zeros(limitSteps_,3);
	tangentList = zeros(limitSteps_,3);
	eleIndexList = zeros(limitSteps_,1);
	paraCoordList = zeros(limitSteps_,2);
	localFrameList = zeros(3,3,limitSteps_);
	localStressTensorList = zeros(limitSteps_,3);
	localPrincipalStressList = zeros(limitSteps_,6);
	localVonMisesStressList = zeros(limitSteps_,1);
	
	[e0, paraCoord0, localFrame0, t1, t2, v0] = SeedSecurityCheck4Robustness(e0, paraCoord0, localFrame0, t1, t2, v0);
	if voidEleMap_(e0), e0 = 0; end
	index = 0;
	reachingBoundary = false;
    eStart = eleStruct_(e0).edge2nextEle;
	switch tracingSche_
		case 'Euler'
			while e0 && index<=limitSteps_ && ~reachingBoundary
				e0Type = meshTypeMap_(e0);
				tracingStepSize = eleCharacterSizeList_(e0)/scalingStepSize4PSL_/shellCurvatureScaling_(e0);
				localCoordIncrement = tracingStepSize*v0;
				paraCoordIncrement = ComputeIncrementOfNaturalCoordFromIncrementOfLocalCoord(localCoordIncrement, localFrame0, t1(:), t2(:));
				paraCoord1 = paraCoord0 + paraCoordIncrement;
				inside_e0 = CheckWhetherGivenNaturalCoordinatesWithinDomain(paraCoord1, e0Type);
				eRef = e0;
				if ~inside_e0
					indexLocal = 0;
					[e0, e0Type, paraCoord1, e0_edgeIDlocal, optAtBoundary, ~] = EdgeSituationProcessing(e0, paraCoord0, paraCoord1, e0Type);
					if optAtBoundary, break; end
				end
				
				%%Advance Point
				index = index + 1;
				[iPhyCoord, iLocalFrame, t1, t2, iStressTensor, iVonMisesStress, iPrincipalStress] = RetrievePhyInfoAtUpdatedPt(e0, paraCoord1, e0Type);
				phyCoordList(index,:) = iPhyCoord;
				eleIndexList(index,1) = e0;
				paraCoordList(index,:) = paraCoord1;
				localFrameList(:,:,index) = iLocalFrame;
				localStressTensorList(index,:) = iStressTensor;
				localPrincipalStressList(index,:) = iPrincipalStress;
				localVonMisesStressList(index,1) = iVonMisesStress;		
				
				%% Temination Condition
				if inside_e0
					[v1, terminationCond] = BidirectionalFeatureProcessing(v0, iPrincipalStress(psDir), localFrame0, iLocalFrame);            
                else
					[v1, terminationCond] = BidirectionalFeatureProcessing_Advanced(v0, iPrincipalStress(psDir), localFrame0, iLocalFrame, eRef, e0_edgeIDlocal);				
                end
                
				tangentList(index,:) = GlobalFrame2Local_Vecs(v1, iLocalFrame, 0);
				if terminationCond, index = index-1; break; end
                if index > 100 && isempty(setdiff(e0, eStart)), break; end
				v0 = v1;
				paraCoord0 = paraCoord1;
				localFrame0 = iLocalFrame;
			end		
		case 'RK2'
			while e0 && index<=limitSteps_ && ~reachingBoundary
				e0Type = meshTypeMap_(e0);
				tracingStepSize = eleCharacterSizeList_(e0)/scalingStepSize4PSL_/shellCurvatureScaling_(e0);
				localCoordIncrement = tracingStepSize*v0;
				paraCoordIncrement = ComputeIncrementOfNaturalCoordFromIncrementOfLocalCoord(localCoordIncrement, localFrame0, t1(:), t2(:)); %%k1
				paraCoord1 = paraCoord0 + paraCoordIncrement;
				inside_e0 = CheckWhetherGivenNaturalCoordinatesWithinDomain(paraCoord1, e0Type);
				eRef = e0;
				if ~inside_e0
					[~, ~, ~, e0_edgeIDlocal, optAtBoundary, paraCoord1e0] = EdgeSituationProcessing(e0, paraCoord0, paraCoord1, e0Type);		
					if optAtBoundary, break; end
					k1 = paraCoord1e0 - paraCoord0;
				else
					k1 = paraCoordIncrement;
				end
				
				%%Evaluate at Mid-point
				midPot = paraCoord0 + k1/2;
				[~, iLocalFrame_midPot, t1_midPot, t2_midPot, ~, ~, iPrincipalStress_midPot] = ...
					RetrievePhyInfoAtUpdatedPt(eRef, midPot, e0Type);
				[v1_midPot, ~] = BidirectionalFeatureProcessing(v0, iPrincipalStress_midPot(psDir), localFrame0, iLocalFrame_midPot);
				localCoordIncrement2 = tracingStepSize * v1_midPot;
				paraCoordIncrement2 = ComputeIncrementOfNaturalCoordFromIncrementOfLocalCoord(localCoordIncrement2, iLocalFrame_midPot, t1_midPot(:), t2_midPot(:)); 
				paraCoord2 = paraCoord0 + paraCoordIncrement2;
				inside_e0 = CheckWhetherGivenNaturalCoordinatesWithinDomain(paraCoord2, e0Type);

				if ~inside_e0
					[e0, e0Type, paraCoord2, e0_edgeIDlocal, optAtBoundary] = EdgeSituationProcessing(e0, paraCoord0, paraCoord2, e0Type);
					if optAtBoundary, break; end
				end
				
				%%Advance Point
				index = index + 1;
				[iPhyCoord, iLocalFrame, t1, t2, iStressTensor, iVonMisesStress, iPrincipalStress] = RetrievePhyInfoAtUpdatedPt(e0, paraCoord2, e0Type);
				phyCoordList(index,:) = iPhyCoord;
				eleIndexList(index,1) = e0;
				paraCoordList(index,:) = paraCoord2;
				localFrameList(:,:,index) = iLocalFrame;
				localStressTensorList(index,:) = iStressTensor;
				localPrincipalStressList(index,:) = iPrincipalStress;
				localVonMisesStressList(index,1) = iVonMisesStress;
	
				%% Temination Condition
				if inside_e0
					[v1, terminationCond] = BidirectionalFeatureProcessing(v0, iPrincipalStress(psDir), localFrame0, iLocalFrame);
				else
					[v1, terminationCond] = BidirectionalFeatureProcessing_Advanced(v0, iPrincipalStress(psDir), localFrame0, iLocalFrame, eRef, e0_edgeIDlocal);
                end			
				tangentList(index,:) = GlobalFrame2Local_Vecs(v1, iLocalFrame, 0);
				if terminationCond, index = index-1; break; end
                if index > 100 && isempty(setdiff(e0, eStart)), break; end
				v0 = v1;
				paraCoord0 = paraCoord2;
				localFrame0 = iLocalFrame;				
			end		
	end

	phyCoordList = phyCoordList(1:index,:);
	tangentList = tangentList(1:index,:);
	eleIndexList = eleIndexList(1:index,1);
	paraCoordList = paraCoordList(1:index,:);
	localFrameList = localFrameList(:,:,1:index);
	localStressTensorList = localStressTensorList(1:index,:);
	localPrincipalStressList = localPrincipalStressList(1:index,:);
	localVonMisesStressList = localVonMisesStressList(1:index,1);		
end

function [v1, terminationCond] = BidirectionalFeatureProcessing(v0, v1Potential, localFrame0, localFrame1)
	global permittedMaxAdjacentTangentAngleDeviation_;
	terminationCond = 0;
	
	v0Global = GlobalFrame2Local_Vecs(v0, localFrame0, 0);
	v1Potentialglobal = GlobalFrame2Local_Vecs(v1Potential, localFrame1, 0);	
	
	angle1 = acos(v0Global*v1Potentialglobal');
	angle2 = acos(-v0Global*v1Potentialglobal');
	if angle1 < angle2
		v1 = v1Potential;
		if angle1 > permittedMaxAdjacentTangentAngleDeviation_/180*pi, terminationCond = 1; end
	else
		v1 = -v1Potential;
		if angle2 > permittedMaxAdjacentTangentAngleDeviation_/180*pi, terminationCond = 1; end
	end
end

function [targetDirection, terminationCond] = BidirectionalFeatureProcessing_Advanced(d0Local, v1Local, localFrame0, iLocalFrame, eRef, refEdgeID)
	terminationCond = 0;
	
	n0 = localFrame0(:,3);
	n1 = iLocalFrame(:,3);
	[alpha0, tVec] = CheckSignedDihedralAngle(eRef, refEdgeID, localFrame0, iLocalFrame);
	D0Global = GlobalFrame2Local_Vecs(d0Local, localFrame0, 0);
	V1Global = GlobalFrame2Local_Vecs(v1Local, iLocalFrame, 0);
	if abs(alpha0) <= pi/6
		angle1 = acos(V1Global*D0Global(:));
		angle2 = acos(-V1Global*D0Global(:));		
	else
		% turningAng = alpha0/pi*180
		T = [0 -tVec(3) tVec(2); tVec(3) 0 -tVec(1); -tVec(2) tVec(1) 0];
		W = eye(3) + sin(alpha0)*T + (1-cos(alpha0))*(T*T);
		
		D0W = W*D0Global(:);
		%%Numerical instability
		D0W = D0W - (dot(D0W, n0))* n1; D0W = D0W/norm(D0W);
		angle1 = acos(V1Global*D0W);
		angle2 = acos(-V1Global*D0W);	
	end
	
	if angle1 < angle2
		targetDirection = v1Local;
		if angle1 > 70/180*pi, terminationCond = 1; end
	else
		targetDirection = -v1Local;
		if angle2 > 70/180*pi, terminationCond = 1; end
	end	
end

function paraCoordIncrement = ComputeIncrementOfNaturalCoordFromIncrementOfLocalCoord(localCoordIncrement, localFrame, t1, t2)
	e1 = localFrame(:,1);
	e2 = localFrame(:,2);
	JsurfaceMapping = [sum(t1.*e1) sum(t2.*e1); sum(t1.*e2) sum(t2.*e2)];
	paraCoordIncrement = JsurfaceMapping \ localCoordIncrement(:);
	paraCoordIncrement = paraCoordIncrement';
end

function [alpha0, tVec] = CheckSignedDihedralAngle(e1, e1_edgeID, localFrame1, localFrame2)
	global eleStruct_;
	global nodeCoords_;
	
	tarEdgeVec = nodeCoords_(eleStruct_(e1).vtxEdgesGlobal(e1_edgeID,:),:);
	tVec = tarEdgeVec(2,:) - tarEdgeVec(1,:); tVec = tVec / norm(tVec);
	n1 = localFrame1(:,3);
	n2 = localFrame2(:,3);
	
	alpha0 = atan2(dot(tVec, cross(n1,n2)), dot(n1, n2));
end

function [e0, e0Type, paraCoord1, e0_edgeIDlocal, optAtBoundary, paraCoord1Fallback] = EdgeSituationProcessing(e0, paraCoord0, paraCoord1, e0Type)
	global eleStruct_ meshTypeMap_ voidEleMap_;
	tolBoundaryCloseness = 1.0e-2;
	optAtBoundary = false;
	%%Get Intersection
	[e0_edgeIDlocal, paraCoord1Fallback, stateIntersection] = DetectIntersectedEdgeInNaturalFrame(paraCoord0, paraCoord1, e0Type);
	if ~stateIntersection
		warning('No Intersection Found!')
		paraCoord0
		paraCoord1
	end
	e1 = eleStruct_(e0).edge2nextEle(e0_edgeIDlocal);
	if 0==e1%%paraCoord1 is not in the domain 
		%%Endpoint treatment when reaching boundary 
		if norm(paraCoord0-paraCoord1Fallback) > tolBoundaryCloseness
			paraCoord1 = paraCoord1Fallback;
			reachingBoundary = true;
		else %%Already at the boundary
			optAtBoundary = true;
			%break;
		end
	elseif voidEleMap_(e1) %%enter into void element
		optAtBoundary = true;
	else %%paraCoord1 is not in e0 but still in the domain, specificly in e1
		%%Shift computing stencils (natural coordinates) to next element (e1)
		%%Convert Natural Coordinates (NC) at Shared Edge of 2 Adjacent Elements
		e0_sharedEdge_NatCoord = paraCoord1Fallback;
		e1_edgeIDlocal = eleStruct_(e0).edgeIDsOfNextEle(e0_edgeIDlocal);
		e1_sharedEdge_NatCoord = ConvertNaturalCoordsAtSharedEdge(e0, e0_sharedEdge_NatCoord, e0_edgeIDlocal, e0Type, ...
			eleStruct_(e0).vtxEdgesGlobal(e0_edgeIDlocal,:), e1, e1_edgeIDlocal, meshTypeMap_(e1), ...
				eleStruct_(e1).vtxEdgesGlobal(e1_edgeIDlocal,:));
		%%***Maybe Problematic for Sharp Edge***
		e0 = e1;
		e0Type = meshTypeMap_(e0);
		paraCoord1 = e1_sharedEdge_NatCoord;			
	end	
end	

function [iPhyCoord, iLocalFrame, t1, t2, iStressTensor, iVonMisesStress, iPrincipalStress] = RetrievePhyInfoAtUpdatedPt(e0, paraCoord1, e0Type)
	global approxiInterp_;
	global nodeCoords_ eNodMat_;
	global localStressFieldPerEle_ globalStressFieldPerEle_;

	switch e0Type
		case 'T3', iEleNodes = eNodMat_(e0,1:3);
		case 'Q4', iEleNodes = eNodMat_(e0,1:4);
		case 'T6', iEleNodes = eNodMat_(e0,1:6);
		case 'Q8', iEleNodes = eNodMat_(e0,1:8);				
	end
	SF = ShapeFunction(paraCoord1, e0Type);
	iEleNodeCoords = nodeCoords_(iEleNodes,:);
	iPhyCoord = SF * iEleNodeCoords;
	[iLocalFrame, ~, t1, t2] = ComputeLocalFrameAtGivenPosition(paraCoord1, iEleNodeCoords, e0Type);
	if approxiInterp_
		iStressTensor = SF * localStressFieldPerEle_{e0};
	else
		iStressTensorGlobal = SF * globalStressFieldPerEle_{e0};
		iStressTensor = GlobalFrame2Local_StressTensor(iStressTensorGlobal, iLocalFrame, 1);
	end
	iVonMisesStress = ComputeVonMisesStressLocalFrame(iStressTensor);
	iPrincipalStress = ComputePrincipalStressLocalFrame(iStressTensor);	
end

function [dNds, dNdt] = DeShapeFunction(paras, elementType)
	numSamps = size(paras,1);
	s = paras(:,1);
	t = paras(:,2);	
	
	switch elementType
		case 'T3'
			dNds = zeros(numSamps,3); 
			dNdt = zeros(numSamps,3);
			
			dNds(:,1) = -ones(numSamps,1); 
			dNdt(:,1) = -ones(numSamps,1);
			dNds(:,2) = ones(numSamps,1);
			dNdt(:,2) = zeros(numSamps,1);
			dNds(:,3) = zeros(numSamps,1);
			dNdt(:,3) = ones(numSamps,1);			
		case 'T6'
			dNds = zeros(numSamps,6);
			dNdt = zeros(numSamps,6);
			
			dNds(:,1) = 4*s + 4*t - 3;
			dNdt(:,1) = 4*s + 4*t - 3;
			dNds(:,2) = 4*s - 1;
			dNdt(:,2) = zeros(numSamps,1);					
			dNds(:,3) = zeros(numSamps,1);
			dNdt(:,3) = 4*t - 1;						
			dNds(:,4) = 4 - 4*t - 8*s;
			dNdt(:,4) = -4*s;						
			dNds(:,5) = 4*t;
			dNdt(:,5) = 4*s;	
			dNds(:,6) = -4*t;
			dNdt(:,6) = 4 - 8*t - 4*s;			
		case 'Q4'
			dNds = zeros(numSamps,4);
			dNdt = zeros(numSamps,4);	

			dNds(:,1) = -0.25*(1-t);
			dNdt(:,1) = -0.25*(1-s);
			dNds(:,2) = 0.25*(1-t);
			dNdt(:,2) = -0.25*(1+s);					
			dNds(:,3) = 0.25*(1+t);
			dNdt(:,3) = 0.25*(1+s);						
			dNds(:,4) = -0.25*(1+t);
			dNdt(:,4) = 0.25*(1-s);			
		case 'Q8'
			dNds = zeros(numSamps,8);
			dNdt = zeros(numSamps,8);					

			dNds(:,1) = -(s/4 - 1/4).*(t - 1) - ((t - 1).*(s + t + 1))/4;
			dNdt(:,1) = -(s/4 - 1/4).*(t - 1) - (s/4 - 1/4).*(s + t + 1);	
			dNds(:,2) = ((t - 1).*(t - s + 1))/4 - (s/4 + 1/4).*(t - 1);
			dNdt(:,2) = (s/4 + 1/4).*(t - s + 1) + (s/4 + 1/4).*(t - 1);	
			dNds(:,3) = (s/4 + 1/4).*(t + 1) + ((t + 1).*(s + t - 1))/4;
			dNdt(:,3) = (s/4 + 1/4).*(t + 1) + (s/4 + 1/4).*(s + t - 1);	
			dNds(:,4) = (s/4 - 1/4).*(t + 1) + ((t + 1).*(s - t + 1))/4;
			dNdt(:,4) = (s/4 - 1/4).*(s - t + 1) - (s/4 - 1/4).*(t + 1);	
			dNds(:,5) = s.*(t - 1);
			dNdt(:,5) = s.^2 /2 - 1/2;	
			dNds(:,6) = 1/2 - t.^2 /2;
			dNdt(:,6) = -2*t.*(s/2 + 1/2);	
			dNds(:,7) = -s.*(t + 1);
			dNdt(:,7) = 1/2 - s.^2 /2;	
			dNds(:,8) = t.^2 /2 - 1/2;
			dNdt(:,8) = 2*t.*(s/2 - 1/2);		
		otherwise
			error('Un-supported Element Type!');
	end
end

function [edgeID, paraCoord1Fallback, successfulIntersection] = DetectIntersectedEdgeInNaturalFrame(p1, p2, eleType)
	edgeID = 0;
	paraCoord1Fallback = [];
	successfulIntersection = true;
	p1 = AvoidExactLocationAtVertexOrEdge4Robustness(p1, eleType, 1.0e-4);
	tangent = p2 - p1;
	tolEff = 1.0e6 * eps(1.0); 
	if norm(tangent)<1.0e2*tolEff
		error('p2 coincides with p1!');
	end
	%%p1p2 -> y = a*x + b
	if abs(p2(1)-p1(1)) < tolEff 
		opt = 1; %%Parallel to y-axis
	elseif abs(p2(2)-p1(2)) < tolEff
		opt = 2; %% Parallel to x-axis
	else
		opt = 3; %% Arbitrary
		a = (p2(2)-p1(2))/(p2(1)-p1(1));
		b = p1(2) - a*p1(1);
	end
	lowerP12X = min([p1(1) p2(1)]); upperP12X = max([p1(1) p2(1)]);
	lowerP12Y = min([p1(2) p2(2)]); upperP12Y = max([p1(2) p2(2)]);
	
	if strcmp(eleType, 'T3') || strcmp(eleType, 'T6')
		switch opt
			case 1
				if p2(2)>1-p2(1)
					edgeID = 2;
					paraTmp = [p2(1) 1-p2(1)];
					if norm(paraTmp-p1)>tolEff
						paraCoord1Fallback = paraTmp;
					else
						edgeID = 0;
					end
				elseif p2(2) < 0
					edgeID = 1;
					paraTmp = [p2(1) 0];
					if norm(paraTmp-p1)>tolEff
						paraCoord1Fallback = paraTmp;
					else
						edgeID = 0;
					end					
				else
					disp(eleType);
					disp([p1; p2]);
					error('Unexpected Error in Detecting Intersected Edge In Natural Frame!');					
				end
			case 2
				if p2(1)>1-p2(2)
					edgeID = 2;
					paraTmp = [1-p2(2) p2(2)];
					if norm(paraTmp-p1)<tolEff
						error('p2 coincides with p1!');
					end
					paraCoord1Fallback = paraTmp;						
				elseif p2(1) < 0
					edgeID = 3;
					paraTmp = [0 p2(2)];
					if norm(paraTmp-p1)<tolEff
						error('p2 coincides with p1!');
					end
					paraCoord1Fallback = paraTmp;						
				else
					disp(eleType);
					disp([p1; p2]);
					error('Unexpected Error in Detecting Intersected Edge In Natural Frame!');							
				end
			case 3
				for ii=1:3
					switch ii
						case 1
							xi = (-b)/a;
							paraTmp = [xi 0];
							if xi<0 || xi>1, continue; end
							if xi<=lowerP12X || xi>=upperP12X, continue; end %% Out of segment p1p2
							if norm(paraTmp-p1)<tolEff, continue; end
							edgeID = ii;
							paraCoord1Fallback = paraTmp;									
							return	
						case 2
							if -1==a, continue; end %%Parallel to edge 2
							xi = (1-b)/(a+1);
							eta = 1-xi;
							paraTmp = [xi eta];
							if xi<0 || xi>1, continue; end
							if eta<0 || eta>1, continue; end
							if xi<=lowerP12X || xi>=upperP12X, continue; end %% Out of segment p1p2
							if norm(paraTmp-p1)<tolEff, continue; end
							edgeID = ii;
							paraCoord1Fallback = paraTmp;	                                
							return
						case 3
							eta = b;
							paraTmp = [0 eta];
							if eta<0 || eta>1, continue; end
							if eta<=lowerP12Y || eta>=upperP12Y, continue; end %% Out of segment p1p2
							if norm(paraTmp-p1)<tolEff, continue; end
							edgeID = ii;
							paraCoord1Fallback = paraTmp;										
							return										
					end
				end
		end	
	elseif strcmp(eleType, 'Q4') || strcmp(eleType, 'Q8')
		switch opt
			case 1
				if p2(2)>1
					edgeID = 3;
					paraTmp = [p2(1) 1];
					if norm(paraTmp-p1)>tolEff
						paraCoord1Fallback = paraTmp;
					else
						edgeID = 0; 
					end						
				elseif p2(2)<-1
					edgeID = 1;
					paraTmp = [p2(1) -1];
					if norm(paraTmp-p1)>tolEff
						paraCoord1Fallback = paraTmp;
					else
						edgeID = 0; 
					end							
				else
					disp(eleType);
					disp([p1; p2]);
					error('Unexpected Error in Detecting Intersected Edge In Natural Frame!');
				end		
			case 2
				if p2(1)>1
					edgeID = 2;
					paraTmp = [1 p2(2)];
					if norm(paraTmp-p1)>tolEff
						paraCoord1Fallback = paraTmp;
					else
						edgeID = 0; 
					end						
				elseif p2(1)<-1
					edgeID = 4;
					paraTmp = [-1 p2(2)];
					if norm(paraTmp-p1)>tolEff
						paraCoord1Fallback = paraTmp;
					else
						edgeID = 0; 
					end							
				else
					disp(eleType);
					disp([p1; p2]);
					error('Unexpected Error in Detecting Intersected Edge In Natural Frame!');
				end				
			case 3
				for ii=1:4
					switch ii
						case 1
							xi = (-1-b)/a;
							paraTmp = [xi -1];
							if xi<-1 || xi>1, continue; end
							if xi<=lowerP12X || xi>=upperP12X, continue; end %% Out of segment p1p2
							if norm(paraTmp-p1)<tolEff, continue; end
							edgeID = ii;
							paraCoord1Fallback = paraTmp;									
							return								
						case 2
							eta = a*1+b;
							paraTmp = [1 eta];
							if eta<-1 || eta>1, continue; end
							if eta<=lowerP12Y || eta>=upperP12Y, continue; end %% Out of segment p1p2
							if norm(paraTmp-p1)<tolEff, continue; end
							edgeID = ii;
							paraCoord1Fallback = paraTmp;									
							return	
						case 3
							xi = (1-b)/a;
							paraTmp = [xi 1];
							if xi<-1 || xi>1, continue; end
							if xi<=lowerP12X || xi>=upperP12X, continue; end %% Out of segment p1p2
							if norm(paraTmp-p1)<tolEff, continue; end
							edgeID = ii;
							paraCoord1Fallback = paraTmp;									
							return								
						case 4
							eta = -a*1+b;
							paraTmp = [-1 eta];
							if eta<-1 || eta>1, continue; end
							if eta<=lowerP12Y || eta>=upperP12Y, continue; end %% Out of segment p1p2
							if norm(paraTmp-p1)<tolEff, continue; end
							edgeID = ii;
							paraCoord1Fallback = paraTmp;									
							return		
					end
				end
		end	
	else
		error('Un-supported Element type!');
	end
	
	if 0==edgeID
		successfulIntersection = false;
	end
end

function disT = DistanceMap2PtCloud(srcPoints, tarCurve)
	global triangularizedMesh_ triangularizedMeshGraph_
	global distanceMetric_;
	switch distanceMetric_
		case 'Euclidean'
			disT = (tarCurve(:,1) - srcPoints(:,1)').^2;
			disT = disT + (tarCurve(:,2) - srcPoints(:,2)').^2;
			disT = disT + (tarCurve(:,3) - srcPoints(:,3)').^2;
			disT = sqrt(disT);
		case 'Geodesic'
			idtarCurve = nearestNeighbor(triangularizedMesh_, tarCurve);
			idsrcPoints = nearestNeighbor(triangularizedMesh_, srcPoints);	
			[tarCurve_campact, ~, tarCurve_recoverMap] = unique(idtarCurve);
			[srcPoints_campact, ~, srcPoints_recoverMap] = unique(idsrcPoints);
			disT = distances(triangularizedMeshGraph_, tarCurve_campact, srcPoints_campact);
			disT = disT(tarCurve_recoverMap,:);
			disT = disT(:,srcPoints_recoverMap);
	end
end

function faceNormals = EvaluateMeshNormals()
	global numEles_ meshTypeMap_ nodeCoords_ eNodMat_
	faceNormals = zeros(numEles_,3);
	for ii=1:numEles_
		iEleType = meshTypeMap_(ii);
		switch iEleType
			case {'T3', 'T6'}
				iEleCoords = nodeCoords_(eNodMat_(ii,1:3),:);
				[R, origin, t1, t2] = ComputeLocalFrameAtGivenPosition([1 1]/3, iEleCoords, 'T3');
			case {'Q4', 'Q8'}
				iEleCoords = nodeCoords_(eNodMat_(ii,1:4),:);
				[R, origin, t1, t2] = ComputeLocalFrameAtGivenPosition([0 0], iEleCoords, 'Q4');
		end
		faceNormals(ii,:) = R(:,3)';
	end	
end

function [gridX, gridY, gridZ, gridC, gridIndices] = ExpandPSLs2Tubes(PSLs, colorSrc, r)
	%%Syntax
	%% [gridX, gridY, gridZ, gridC] = ExpandPSLs2Tubes(PSLs, colorSrc, r)
	gridX = [];
	gridY = [];
	gridZ = [];
	gridC = [];
	gridIndices = [];
	if isempty(PSLs), return; end
	n = 8; 
	numLines = length(PSLs);
	gridXYZ = zeros(3,n+1,1);
	gridC = zeros(n+1,1);
	gridIndices = struct('arr', []);
	gridIndices = repmat(gridIndices, numLines, 1);	
	for ii=1:numLines		
		curve = PSLs(ii).phyCoordList';
		npoints = size(curve,2);
		%deltavecs: average for internal points. first strecth for endpoitns.		
		dv = curve(:,[2:end,end])-curve(:,[1,1:end-1]);		
		%make nvec not parallel to dv(:,1)
		nvec=zeros(3,1); 
		[~,idx]=min(abs(dv(:,1))); 
		nvec(idx)=1;
		%precalculate cos and sing factors:
		cfact=repmat(cos(linspace(0,2*pi,n+1)),[3,1]);
		sfact=repmat(sin(linspace(0,2*pi,n+1)),[3,1]);
		%Main loop: propagate the normal (nvec) along the tube
		xyz = zeros(3,n+1,npoints+2);
		for k=1:npoints
			convec=cross(nvec,dv(:,k));
			convec=convec./norm(convec);
			nvec=cross(dv(:,k),convec);
			nvec=nvec./norm(nvec);
			%update xyz:
			xyz(:,:,k+1)=repmat(curve(:,k),[1,n+1]) + cfact.*repmat(r*nvec,[1,n+1]) + sfact.*repmat(r*convec,[1,n+1]);
        end
		%finally, cap the ends:
		xyz(:,:,1)=repmat(curve(:,1),[1,n+1]);
		xyz(:,:,end)=repmat(curve(:,end),[1,n+1]);
		gridIndices(ii).arr = size(gridXYZ,3) : size(gridXYZ,3)+size(xyz,3)-1;
		gridXYZ(:,:,end+1:end+npoints+2) = xyz;	
		color = colorSrc(ii).arr;	
		c = [color(1) color color(end)];
		c = repmat(c, n+1, 1);
		gridC(:,end+1:end+npoints+2) = c;
	end		
	gridX = squeeze(gridXYZ(1,:,:)); 
	gridX(:,1) = [];
	gridY = squeeze(gridXYZ(2,:,:)); 
	gridY(:,1) = [];
	gridZ = squeeze(gridXYZ(3,:,:)); 
	gridZ(:,1) = [];
	gridC(:,1) = [];
end

function [patchX, patchY, patchZ] = GetSpheresOfPoints(points, sphereRadius)
	numSeedPoints = size(points,1);
	[sphereX,sphereY,sphereZ] = sphere(10);
	sphereX = sphereRadius*sphereX;
	sphereY = sphereRadius*sphereY;
	sphereZ = sphereRadius*sphereZ;
	nn = size(sphereX,1);

	patchX = sphereX; ctrX = points(:,1);
	patchX = repmat(patchX, numSeedPoints, 1); ctrX = repmat(ctrX, 1, nn); ctrX = reshape(ctrX', numel(ctrX), 1);
	patchX = ctrX + patchX;
	
	patchY = sphereY; ctrY = points(:,2);
	patchY = repmat(patchY, numSeedPoints, 1); ctrY = repmat(ctrY, 1, nn); ctrY = reshape(ctrY', numel(ctrY), 1);
	patchY = ctrY + patchY;	
	
	patchZ = sphereZ; ctrZ = points(:,3);
	patchZ = repmat(patchZ, numSeedPoints, 1); ctrZ = repmat(ctrZ, 1, nn); ctrZ = reshape(ctrZ', numel(ctrZ), 1);
	patchZ = ctrZ + patchZ;
end

function coordsOut = GlobalFrame2Local_Coords(coordsIn, R, origin, opt)
	switch opt
		case 1 %%Global -> Local
			coordsLocal3x3 = R' * (coordsIn - origin)';
			coordsLocal = coordsLocal3x3(1:2,:)';
			coordsOut = coordsLocal;
		case 0 %%Local -> Global
			coordsLocal3x3 = [coordsIn, zeros(size(coordsIn,1), 1)]';
			coordsGloabl = origin(:) + R * coordsLocal3x3;
			coordsOut = coordsGloabl';
	end
end

function stressOut = GlobalFrame2Local_StressTensor(stressIn, R, opt)
	%% sigma_xx, sigma_yy, sigma_zz, sigma_yz, sigma_zx, sigma_xy (global)
	numComps = size(R,3);
	if numComps~=size(stressIn,1)
		error('Un-matched Dimensions!');
	end
	switch opt
		case 1 %%Global -> Local
			stressOut = zeros(numComps,3);
			for ii=1:numComps
				iR = R(:,:,ii);
				iStressIn = stressIn(ii,:);
				iStressLocal3x3 = iR' * iStressIn([1 6 5; 6 2 4; 5 4 3]) * iR;
				stressOut(ii,:) = iStressLocal3x3([1 5 2]);
			end
		case 0 %%Local -> Global
			stressOut = zeros(numComps,6);
			for ii=1:numComps
				iR = R(:,:,ii);
				iStressIn = stressIn(ii,:);
				iStressLocal3x3 = zeros(3,3);
				iStressLocal3x3(1,1:2) = iStressIn([1 3]);
				iStressLocal3x3(2,1:2) = iStressIn([3 2]);
				iStressGlobal3x3 = iR * iStressLocal3x3 * iR';
				stressOut(ii,:) = iStressGlobal3x3([1 5 9 6 3 2]);
			end
	end
end

function vecsOut = GlobalFrame2Local_Vecs(vecsIn, R, opt)
	switch opt
		case 1 %%Global -> Local
			vecsLocal3x3 = R' * vecsIn';
			vecsLocal = vecsLocal3x3(1:2,:)';
			vecsOut = vecsLocal;
		case 0 %%Local -> Global
			vecsLocal3x3 = [vecsIn, zeros(size(vecsIn,1), 1)]';
			vecsGlobal = (R * vecsLocal3x3)';
			vecsOut = vecsGlobal;
	end
end

function [seed0, phyCoord, localCoord, localFrame, t1, t2, localStressTensor, localPrincipalStress, localVonMisesStress, opt] = PreparingForTracing(startPoint)
	global nodeCoords_ eNodMat_ meshTypeMap_ voidEleMap_;
	global refScalingSize_ approxiInterp_;
	global localStressFieldPerEle_ globalStressFieldPerEle_;
	
	ele = startPoint(1);
	paras = startPoint(2:3); %% Natural Coordinates
	
	phyCoord = NaN; 
	localCoord = NaN;
	localFrame = NaN;
	t1 = NaN;
	t2 = NaN;
	localStressTensor = NaN; 
	localPrincipalStress = NaN; 
	localVonMisesStress = NaN;
	if voidEleMap_(ele), opt = false; return; end
	
	eleType = meshTypeMap_(ele);
	switch eleType
		case 'T3'
			eleOrder = 1;
			nodeListPerEle = 1:3;			
		case 'Q4'
			eleOrder = 1;
			nodeListPerEle = 1:4;
		case 'T6'
			eleOrder = 2;
			nodeListPerEle = 1:6;
		case 'Q8'
			eleOrder = 2;
			nodeListPerEle = 1:8;
	end
	% paras = AvoidExactLocationAtVertexOrEdge4Robustness(paras, eleType, 1.0e-4);
	state_withinit = CheckWhetherGivenNaturalCoordinatesWithinDomain(paras, eleType);
	if state_withinit
		iEleNodes = eNodMat_(ele,nodeListPerEle);
		iEleNodeCoords = nodeCoords_(iEleNodes,:);
		SF = ShapeFunction(paras, eleType);
		phyCoord = SF * iEleNodeCoords;

		[localFrame, origin0, t1, t2] = ComputeLocalFrameAtGivenPosition(paras, iEleNodeCoords, eleType);
		localCoord = [0 0];
			
		if approxiInterp_
			localStressTensor = SF * localStressFieldPerEle_(ele,:);
		else
			globalStressTensor = SF * globalStressFieldPerEle_{ele};
			localStressTensor = GlobalFrame2Local_StressTensor(globalStressTensor, localFrame, 1);
		end
		localPrincipalStress = ComputePrincipalStressLocalFrame(localStressTensor);
		localVonMisesStress = ComputeVonMisesStressLocalFrame(localStressTensor);
	end
	seed0 = [ele paras];
	opt = state_withinit;
end

function val = PrincipalStressLineStruct()
	val = struct(...
		'numIntergPts',					0,	...
		'midPointPosition',				0,	...		
		'phyCoordList',					[], ...
		'eleIndexList',					[], ...
		'paraCoordList',				[], ...
		'localFrames',					[],	...		
		'localStressTensorList',		[],	...
		'localVonMisesStressList',		[], ...
		'localPrincipalStressList', 	[], ...
		'strainEnergyApprox',			0	...	
	);	
end

function [e0out, paraCoord0out, localFrame0out, t1out, t2out, v0out] = SeedSecurityCheck4Robustness(e0in, paraCoord0in, localFrame0in, t1in, t2in, v0in)
	%% This function is introduced to counteract the situation where the seed point is exactly located at the element 
	%% edge or vertex meanwhile the initial direction is outward
	global meshTypeMap_;
	global eleStruct_;
	global nodeCoords_;
	global eNodMat_;
	
	e0out = e0in;
	paraCoord0out = paraCoord0in;
	localFrame0out = localFrame0in; 
	t1out = t1in;
	t2out = t2in;
	v0out = v0in;
	iEleType = meshTypeMap_(e0in);
	tolEff = 1.0e6 * eps(1.0);
	s = paraCoord0in(1); t = paraCoord0in(2);
	if strcmp(iEleType, 'T3') || strcmp(iEleType, 'T6')
		%%Edge-Vertex Check
		if s<tolEff || t<tolEff || s+t > 1-tolEff
			if norm([s t]) < tolEff || norm([s t]) > 1 - tolEff %% at Vertex
				% warning('Approaching Vertex!');
				[~, nodeID] = min(vecnorm([s t] - [0 0; 1 0; 0 1], 2, 2));
				[e0out, paraCoord0out, localFrame0out, t1out, t2out, v0out] = ...
					TransferToNextElement_StandardSharedNodeCase(nodeID, e0in, [s t], localFrame0in, t1in, t2in, v0in);
			else %% at Edge
				% warning('Approaching Edge!');
				if s < tolEff
					s = round(s);
					edgeID = 3;
				end
				if t < tolEff
					t = round(t);
					edgeID = 1;
				end
				if 1-s-t < tolEff
					edgeID = 2;
				end
				eNormal = DetermineExternalNormalOfGivenEdgeInLocalFrame(e0in, edgeID);
				if acos(eNormal*v0in') > pi/2 %%v0in is inward
					return
				else
					e1 = eleStruct_(e0in).edge2nextEle(edgeID);
					if 0==e1, return; end %%approaching boundary
					[e0out, paraCoord0out, localFrame0out, t1out, t2out, v0out] = ...
						TransferToNextElement_StandardSharedEdgeCase(e0in, edgeID, [s t], iEleType, e1, v0in, localFrame0in);						
				end
			end
		else
			return
		end
	elseif strcmp(iEleType, 'Q4') || strcmp(iEleType, 'Q8')
		%%Edge-Vertex Check
		if s < -1+tolEff || s > 1-tolEff || t < -1+tolEff || t > 1-tolEff
			if norm([s t]) > sqrt(2) - tolEff %% at Vertex
				% warning('Approaching Vertex!');
				[~, nodeID] = min(vecnorm([s t] - [-1 -1; 1 -1; 1 1; -1 1], 2, 2));
				[e0out, paraCoord0out, localFrame0out, t1out, t2out, v0out] = ...
					TransferToNextElement_StandardSharedNodeCase(nodeID, e0in, [s t], localFrame0in, t1in, t2in, v0in);				
			else %% at Edge
				% warning('Approaching Edge!');
				if s < -1+tolEff
					s = round(s);
					edgeID = 4;
				elseif s > 1-tolEff
					s = round(s);
					edgeID = 2;
				end
				if t < -1+tolEff
					t = round(t);
					edgeID = 1;
				elseif t > 1-tolEff
					t = round(t);
					edgeID = 3;
				end
				eNormal = DetermineExternalNormalOfGivenEdgeInLocalFrame(e0in, edgeID);
				if acos(eNormal*v0in') > pi/2 %%v0in is inward
					return
				else
					e1 = eleStruct_(e0in).edge2nextEle(edgeID);
					if 0==e1, return; end %%approaching boundary
					[e0out, paraCoord0out, localFrame0out, t1out, t2out, v0out] = ...
						TransferToNextElement_StandardSharedEdgeCase(e0in, edgeID, [s t], iEleType, e1, v0in, localFrame0in);			
				end
			end		
		else
			return
		end
	else
		error('Un-supported Element type!');
	end
end

function [e0out, paraCoord0out, localFrame0out, t1out, t2out, v0out] = ...
			TransferToNextElement_StandardSharedEdgeCase(e0in, edgeID, e0Para, e0Type, e1, v0in, localFrame0in)
	global nodeCoords_;
	global eNodMat_;
	global eleStruct_;
	global meshTypeMap_;
	
	e1_edgeIDlocal = eleStruct_(e0in).edgeIDsOfNextEle(edgeID);
	e1_sharedEdge_NatCoord = ConvertNaturalCoordsAtSharedEdge(e0in, e0Para, edgeID, e0Type, ...
		eleStruct_(e0in).vtxEdgesGlobal(edgeID,:), e1, e1_edgeIDlocal, meshTypeMap_(e1), ...
			eleStruct_(e1).vtxEdgesGlobal(e1_edgeIDlocal,:));
	v0inGlobal = GlobalFrame2Local_Vecs(v0in, localFrame0in, 0);
	switch meshTypeMap_(e1)
		case 'T3'
			iEleCoords = nodeCoords_(eNodMat_(e1,1:3),:);
		case 'T6'
			iEleCoords = nodeCoords_(eNodMat_(e1,1:6),:);
		case 'Q4'
			iEleCoords = nodeCoords_(eNodMat_(e1,1:4),:);
		case 'Q8'
			iEleCoords = nodeCoords_(eNodMat_(e1,1:8),:);	
	end
	[e1R, ~, e1t1, e2t2] = ComputeLocalFrameAtGivenPosition(e1_sharedEdge_NatCoord, iEleCoords, meshTypeMap_(e1));
	v0outLocal = GlobalFrame2Local_Vecs(v0inGlobal, e1R, 1);
	
	
	e0out = e1;
	paraCoord0out = e1_sharedEdge_NatCoord;
	localFrame0out = e1R; 
	t1out = e1t1;
	t2out = e2t2;
	v0out = v0outLocal;		
end

function eNormal = DetermineExternalNormalOfGivenEdgeInLocalFrame(e0, edgeID)
	global eleStruct_;
	global meshTypeMap_;
	global nodeCoords_;
	global eNodMat_;
	iEleType = meshTypeMap_(e0);
	vtxCoordsGlobal = eleStruct_(e0).vtxEdgesGlobal(edgeID,:);
	if strcmp(iEleType, 'T3') || strcmp(iEleType, 'T6')
		iEleCoords = nodeCoords_(eNodMat_(e0,1:3),:);
		[R, origin, ~, ~] = ComputeLocalFrameAtGivenPosition([1 1]/3, iEleCoords, 'T3');
	elseif strcmp(iEleType, 'Q4') || strcmp(iEleType, 'Q8')
		iEleCoords = nodeCoords_(eNodMat_(e0,1:4),:);
		[R, origin, ~, ~] = ComputeLocalFrameAtGivenPosition([0 0], iEleCoords, 'Q4');
	end
	vtxCoordsLocal = GlobalFrame2Local_Coords(nodeCoords_(vtxCoordsGlobal,:), R, origin, 1);
	edgeVecLocal = vtxCoordsLocal(2,:) - vtxCoordsLocal(1,:); edgeVecLocal = edgeVecLocal / norm(edgeVecLocal);
	normVecRaw = [-edgeVecLocal(2) edgeVecLocal(1)];
	refVec = sum(edgeVecLocal,1)/2; refVec = refVec/norm(refVec);
	devAng = acos(refVec*normVecRaw');
	if devAng < pi/2
		eNormal = normVecRaw;
	else
		eNormal = -normVecRaw;
	end
end

function [e0out, paraCoord0out, localFrame0out, t1out, t2out, v0out] = ...
	TransferToNextElement_StandardSharedNodeCase(nodeID, e0in, e0Para, localFrame0in, t1in, t2in, v0in)
	global nodeCoords_;
	global eNodMat_;
	global nodStruct_;
	global eleCentroids_;
	global meshTypeMap_;
	
	e0out = e0in;
	paraCoord0out = e0Para;
	localFrame0out = localFrame0in;
	t1out = t1in; 
	t2out = t2in;
	v0out = v0in;
	
	v0inGlobal = GlobalFrame2Local_Vecs(v0in, localFrame0in, 0);
	tarNodeGlobal = eNodMat_(e0in, nodeID);
	tarNodeCoord = nodeCoords_(tarNodeGlobal,:);
	adjElesSharingThisNode = nodStruct_(tarNodeGlobal).adjacentEles;
	numAdjEles = numel(adjElesSharingThisNode);
	refSize = vecnorm(tarNodeCoord-eleCentroids_,2,2)/4;
	
	inside = false;
	for ii=1:numAdjEles
		iEle = adjElesSharingThisNode(ii);
		iEleType = meshTypeMap_(iEle);
		iPerturbationPos = tarNodeCoord + v0inGlobal*refSize(ii);
		switch iEleType
			case 'T3'
				iEleNodes = eNodMat_(iEle, 1:3);
			case 'Q4'
				iEleNodes = eNodMat_(iEle, 1:4);
			case 'T6'
				iEleNodes = eNodMat_(iEle, 1:6);
			case 'Q8'
				iEleNodes = eNodMat_(iEle, 1:8);
		end
		iEleNodeCoordsPhy = nodeCoords_(iEleNodes,:);
		[~, inside] = NewtonRhapson_LocatingPosition3D(iEleNodeCoordsPhy, iPerturbationPos, iEleType);
		if inside
			e0out = iEle;
			jNodeLocal = find(0==(tarNodeGlobal-iEleNodes));
			switch iEleType
				case {'T3','T6'}
					switch jNodeLocal
						case 1
							paraCoord0out = [0 0];
						case 2
							paraCoord0out = [1 0];
						case 3
							paraCoord0out = [0 1];
					end
				case {'Q4','Q8'}
					switch jNodeLocal
						case 1
							paraCoord0out = [-1 -1];
						case 2
							paraCoord0out = [1 -1];
						case 3
							paraCoord0out = [1 1];
						case 4
							paraCoord0out = [-1 1];
					end				
			end
			[localFrame0out, origin, t1out, t2out] = ComputeLocalFrameAtGivenPosition(paraCoord0out, iEleNodeCoordsPhy, iEleType);
			v0out = GlobalFrame2Local_Vecs(v0inGlobal, localFrame0out, 1);
			return
		end			
	end
end	

function [paras, inside] = NewtonRhapson_LocatingPosition3D(eleNodeCoords, target, iEleType)
	% Project a 3D point 'target' onto an isoparametric shell element
	% (T3, T6, Q4, Q8). Returns natural coords [s t] and 'inside' flag
	% (true if converged and (s,t) lies in the element domain).
	%
	% Dependencies: ShapeFunction([s t], type), DeShapeFunction([s t], type)
	
	target = target(:).';            % 1x3
	maxIts = 30;
	tol_r  = 1e-10;                  % residual tolerance (length units)
	tol_x  = 1e-12;                  % step tolerance in (s,t)
	tol_g  = 1e-12;                  % gradient tolerance
	lambda0 = 0;                     % LM damping (0 = pure GN/Newton)
	project_to_domain = true;
	
	% --- initial guess
	switch iEleType
	case 'T3',  s = 1/3; t = 1/3;
	case 'T6',  s = 1/3; t = 1/3;
	case 'Q4',  s = 0.0; t = 0.0;
	case 'Q8',  s = 0.0; t = 0.0;
	otherwise, error('Unknown iEleType');
	end
	
	for it = 1:maxIts
		% shape functions and tangents
		N = ShapeFunction([s t], iEleType);                 % 1xMe
		[dNds, dNdt] = DeShapeFunction([s t], iEleType);    % 1xMe
	
		x   = N    * eleNodeCoords;       % 1x3
		xs  = dNds * eleNodeCoords;       % 1x3  = ∂x/∂s
		xt  = dNdt * eleNodeCoords;       % 1x3  = ∂x/∂t
	
		r   = x - target;                 % 1x3
		rn  = norm(r);
	
		% Converged by residual length?
		if rn <= tol_r, break; end
	
		% Gauss–Newton / LM in parameter space: (H + λI) Δ = -g
		H = [ dot(xs,xs), dot(xs,xt);
			dot(xt,xs), dot(xt,xt) ];
		g = [ dot(xs,r);  dot(xt,r) ];
	
		% Optionally adapt λ when H is near-singular
		lambda = lambda0;
		if rcond(H) < 1e-12
			lambda = max(lambda, 1e-12 * trace(H));
		end
	
		dx = - (H + lambda*eye(2)) \ g;   % 2x1 step in (s,t)
	
		% Backtracking line search on φ = 0.5||r||^2
		phi0  = 0.5*rn*rn;
		alpha = 1.0;
		for bt = 1:8
			s_try = s + alpha*dx(1);
			t_try = t + alpha*dx(2);
	
			if project_to_domain
				switch iEleType
					case {'T3','T6'}
						[s_try, t_try] = ProjectToTriangle([s_try, t_try], [0 0], [1 0], [0 1]);
					case {'Q4','Q8'}
						s_try = max(-1, min( 1, s_try));
						t_try = max(-1, min( 1, t_try));
					otherwise
						error('Unknown iEleType');
				end
			end
	
			% trial residual
			N_try = ShapeFunction([s_try t_try], iEleType);
			x_try = N_try * eleNodeCoords;
			r_try = x_try - target;
			phi   = 0.5 * dot(r_try, r_try);
	
			if phi <= phi0 * (1 - 1e-4*alpha)    % Armijo
				s = s_try; t = t_try;
				break;
			else
				alpha = 0.5*alpha;
			end
		end
	
		% Converged by small parameter step or small gradient?
		if norm(alpha*dx) <= tol_x || norm(g) <= tol_g
			break;
		end
	end
	
	paras = [s t];

	% Inside test: reference domain with small slack to accept boundary
	slack = 1e-9;
	switch iEleType
	case {'T3','T6'}
		inside = (s >= -slack) && (t >= -slack) && (s+t <= 1+slack);
	case {'Q4','Q8'}
		inside = (s >= -1-slack) && (s <= 1+slack) && (t >= -1-slack) && (t <= 1+slack);
	otherwise
		inside = false;
	end
	
	% Also require that iterations did not exhaust
	inside = inside && (it < maxIts);

end

function [sp,tp] = ProjectToTriangle(st, A, B, C)
% Orthogonal projection of point st=[s t] onto triangle ABC in R^2
    P = st(:); A2=A(:); B2=B(:); C2=C(:);
    v0 = B2 - A2; v1 = C2 - A2; v2 = P - A2;
    M  = [v0 v1];                             % 2x2
    uv = M \ v2;  u = uv(1); v = uv(2);
    if u >= 0 && v >= 0 && u+v <= 1
        Q = A2 + M*[u; v];
        sp = Q(1); tp = Q(2); return;
    end
    [pAB,dAB] = ProjOnSeg(P, A2, B2);
    [pAC,dAC] = ProjOnSeg(P, A2, C2);
    [pBC,dBC] = ProjOnSeg(P, B2, C2);
    [~,idx] = min([dAB,dAC,dBC]);
    if idx==1, Q=pAB; elseif idx==2, Q=pAC; else, Q=pBC; end
    sp = Q(1); tp = Q(2);
end

function [p,d2] = ProjOnSeg(P,A,B)
    AB = B - A; t = dot(P-A,AB)/max(dot(AB,AB),eps);
    t  = max(0,min(1,t));
    p  = A + t*AB;
    d2 = sum((P - p).^2);
end

function N = ShapeFunction(paras, elementType)
	numSamps = size(paras,1);
	s = paras(:,1);
	t = paras(:,2);
	switch elementType
		case 'T3'
			N = zeros(numSamps,3);
			N(:,1) = 1-s-t;
			N(:,2) = s;
			N(:,3) = t;				
		case 'T6'
			N = zeros(numSamps,6);
			N(:,1) = (1-s-t).*(2*(1-s-t)-1);
			N(:,2) = s.*(2*s-1);
			N(:,3) = t.*(2*t-1);
			N(:,4) = 4*s.*(1-s-t);
			N(:,5) = 4*s.*t;
			N(:,6) = 4*t.*(1-s-t);			
		case 'Q4'
			N = zeros(numSamps,4);
			N(:,1) = 0.25*(1-s).*(1-t);
			N(:,2) = 0.25*(1+s).*(1-t);
			N(:,3) = 0.25*(1+s).*(1+t);
			N(:,4) = 0.25*(1-s).*(1+t);			
		case 'Q8'
			N = zeros(numSamps,8);
			N(:,1) = 0.25*(1-s).*(1-t).*(-s-t-1);
			N(:,2) = 0.25*(1+s).*(1-t).*(s-t-1);
			N(:,3) = 0.25*(1+s).*(1+t).*(s+t-1);
			N(:,4) = 0.25*(1-s).*(1+t).*(-s+t-1);
			N(:,5) = 0.5*(1-s.^2).*(1-t);
			N(:,6) = 0.5*(1+s).*(1-t.^2);
			N(:,7) = 0.5*(1-s.^2).*(1+t);
			N(:,8) = 0.5*(1-s).*(1-t.^2);		
		otherwise
			error('Un-supported Element Type!');
	end
end

function ShowMeshNormals(faceNormals)
	global silhouetteStruct_ eleCentroids_
	hdEles = patch(silhouetteStruct_); hold('on');
	hdNormals = quiver3(eleCentroids_(:,1), eleCentroids_(:,2), eleCentroids_(:,3), faceNormals(:,1), faceNormals(:,2), faceNormals(:,3));
	set(hdEles, 'FaceColor', [65 174 118]/255, 'FaceAlpha', 1.0, 'EdgeColor', 'k');
	set(hdNormals, 'LineWidth', 0.5, 'Color', 'r', 'MaxHeadSize', 0.5);
	view(3);
	axis('equal', 'tight', 'off'); 	
end

function ShowPrincipalStressDirectionsByArrowsGlobally(varargin)
	global nodeCoords_ numEles_ eNodMat_ silhouetteStruct_;
	global meshTypeMap_ voidEleMap_;	
	global localStressFieldPerEle_ globalStressFieldPerEle_;
	
	if 1==nargin
		eleList = varargin{1};
		eleList = unique(eleList);
	else
		eleList = (1:numEles_)';
	end
	volidEles = find(voidEleMap_);
	eleList = setdiff(eleList, volidEles);
	
	numElesConsidered = numel(eleList);
	originList = zeros(numElesConsidered,3);
	dirMajorList = zeros(numElesConsidered,3);
	dirMinorList = zeros(numElesConsidered,3);
	
	opt_local2Global2Local = 1;
	for ii=1:numElesConsidered
		iEle = eleList(ii);
		switch meshTypeMap_(iEle)
			case 'T3'
				paras = [1 1]/3;
				nodeListPerEle = 1:3;
				iEleNodes = eNodMat_(iEle,nodeListPerEle);
				SF = ShapeFunction(paras, 'T3');
				originList(ii,:) = SF * nodeCoords_(iEleNodes, :);
				
				[iLocalFrame, ~, ~, ~] = ComputeLocalFrameAtGivenPosition(paras, nodeCoords_(iEleNodes, :), 'T3');
				iEleStressesLocal = localStressFieldPerEle_{iEle};
				
				if opt_local2Global2Local
					iEleStressesGlobal = globalStressFieldPerEle_{iEle};
					iStressTensorGlobal = SF * iEleStressesGlobal;
					iStressTensor = GlobalFrame2Local_StressTensor(iStressTensorGlobal, iLocalFrame, 1);
				else
					iStressTensor = SF * iEleStressesLocal;
				end
				
				iPSlocal = ComputePrincipalStressLocalFrame(iStressTensor);
				iMajor = iPSlocal(1,5:6); iMinor = iPSlocal(1,2:3);
				
				
				dirMajorList(ii,:) = GlobalFrame2Local_Vecs(iMajor, iLocalFrame, 0)';
				dirMinorList(ii,:) = GlobalFrame2Local_Vecs(iMinor, iLocalFrame, 0)';
			case 'Q4'
				paras = [0 0];
				nodeListPerEle = 1:4;
				iEleNodes = eNodMat_(iEle,nodeListPerEle);
				SF = ShapeFunction(paras, 'Q4');
				originList(ii,:) = SF * nodeCoords_(iEleNodes, :);
				
				[iLocalFrame, ~, ~, ~] = ComputeLocalFrameAtGivenPosition(paras, nodeCoords_(iEleNodes, :), 'Q4');
				iEleStressesLocal = localStressFieldPerEle_{iEle};
				
				if opt_local2Global2Local
					iEleStressesGlobal = globalStressFieldPerEle_{iEle};
					iStressTensorGlobal = SF * iEleStressesGlobal;
					iStressTensor = GlobalFrame2Local_StressTensor(iStressTensorGlobal, iLocalFrame, 1);				
				else
					iStressTensor = SF * iEleStressesLocal;
				end
				
				iPSlocal = ComputePrincipalStressLocalFrame(iStressTensor);
				iMajor = iPSlocal(1,5:6); iMinor = iPSlocal(1,2:3);
	
				dirMajorList(ii,:) = GlobalFrame2Local_Vecs(iMajor, iLocalFrame, 0)';
				dirMinorList(ii,:) = GlobalFrame2Local_Vecs(iMinor, iLocalFrame, 0)';				
			case 'T6'
				paras = [1 1]/3;
				nodeListPerEle = 1:6;
				iEleNodes = eNodMat_(iEle,nodeListPerEle);
				SF = ShapeFunction(paras, 'T6');
				originList(ii,:) = SF * nodeCoords_(iEleNodes, :);
				iEleStressesLocal = localStressFieldPerEle_{iEle};
				[iLocalFrame, ~, ~, ~] = ComputeLocalFrameAtGivenPosition(paras, nodeCoords_(iEleNodes, :), 'T6');
				
				if opt_local2Global2Local
					iEleStressesGlobal = globalStressFieldPerEle_{iEle};
					iStressTensorGlobal = SF * iEleStressesGlobal;
					iStressTensor = GlobalFrame2Local_StressTensor(iStressTensorGlobal, iLocalFrame, 1);					
				else
					iStressTensor = SF * iEleStressesLocal;
				end

				iPSlocal = ComputePrincipalStressLocalFrame(iStressTensor);
				iMajor = iPSlocal(1,5:6); iMinor = iPSlocal(1,2:3);	
				
				dirMajorList(ii,:) = GlobalFrame2Local_Vecs(iMajor, iLocalFrame, 0)';
				dirMinorList(ii,:) = GlobalFrame2Local_Vecs(iMinor, iLocalFrame, 0)';					
			case 'Q8'
				paras = [0 0];
				nodeListPerEle = 1:8;
				iEleNodes = eNodMat_(iEle,nodeListPerEle);
				SF = ShapeFunction(paras, 'Q8');
				originList(ii,:) = SF * nodeCoords_(iEleNodes, :);
				iEleStressesLocal = localStressFieldPerEle_{iEle};
				[iLocalFrame, ~, ~, ~] = ComputeLocalFrameAtGivenPosition(paras, nodeCoords_(iEleNodes, :), 'Q8');
				
				if opt_local2Global2Local
					iEleStressesGlobal = globalStressFieldPerEle_{iEle};
					iStressTensorGlobal = SF * iEleStressesGlobal;
					iStressTensor = GlobalFrame2Local_StressTensor(iStressTensorGlobal, iLocalFrame, 1);						
				else
					iStressTensor = SF * iEleStressesLocal;
				end
				
				iPSlocal = ComputePrincipalStressLocalFrame(iStressTensor);
				iMajor = iPSlocal(1,5:6); iMinor = iPSlocal(1,2:3);
				
				dirMajorList(ii,:) = GlobalFrame2Local_Vecs(iMajor, iLocalFrame, 0)';
				dirMinorList(ii,:) = GlobalFrame2Local_Vecs(iMinor, iLocalFrame, 0)';					
		end
	end


	elesVis.vertices = nodeCoords_;
	elesVis.faces = silhouetteStruct_.faces(eleList,:);
	
	figure;
	hd_ele = patch(elesVis); hold on;
	hd_SampledPts = plot3(originList(:,1), originList(:,2), originList(:,3), '.r', 'LineWidth', 2, 'MarkerSize', 10); hold on;
	hd_major = quiver3(originList(:,1), originList(:,2), originList(:,3), dirMajorList(:,1), dirMajorList(:,2), dirMajorList(:,3), 'LineWidth', 1, 'Color', [252 141 98]/255); hold on;
	hd_major = quiver3(originList(:,1), originList(:,2), originList(:,3), dirMinorList(:,1), dirMinorList(:,2), dirMinorList(:,3), 'LineWidth', 1, 'Color', [0 176 80]/255);
	set(hd_ele, 'FaceColor', [204 204 204]/255, 'FaceAlpha', 1.0, 'EdgeColor', 'k', 'LineWidth', 2);
	if 0==sum(nodeCoords_(:,3))
		view(2); %light	
	else
		view(3); %light	 
	end
	axis equal
end

function ShowProblemDescription()
	global nodeCoords_;
	global loadingCond_;
	global fixingCond_;
	global boundingBox_;
	global silhouetteStruct_;

	
	hd = patch(silhouetteStruct_); hold('on');
	set(hd, 'FaceColor', [65 174 118]/255, 'FaceAlpha', 1.0, 'EdgeColor', 'k');
    if ~isempty(loadingCond_)
		lB = 0.2; uB = 1.0;
		amps = vecnorm(loadingCond_(:,2:4),2,2);
		maxAmp = max(amps); minAmp = min(amps);
		if abs(minAmp-maxAmp)/(minAmp+maxAmp)<0.1
			scalingFac = 1;
		else
			if minAmp/maxAmp>lB/uB, lB = minAmp/maxAmp; end
			scalingFac = lB + (uB-lB)*(amps-minAmp)/(maxAmp-minAmp);
		end
		loadingDirVec = loadingCond_(:,2:4)./amps.*scalingFac;
		coordLoadedNodes = nodeCoords_(loadingCond_(:,1),:);
		amplitudesF = mean(boundingBox_(2,:)-boundingBox_(1,:))/5 * loadingDirVec;
		hold('on'); quiver3(coordLoadedNodes(:,1), coordLoadedNodes(:,2), coordLoadedNodes(:,3), amplitudesF(:,1), ...
			amplitudesF(:,2), amplitudesF(:,3), 0, 'Color', [255 127 0.0]/255, 'LineWidth', 2, 'MaxHeadSize', 1, 'MaxHeadSize', 1);
	end
    if ~isempty(fixingCond_)
		tarNodeCoord = nodeCoords_(fixingCond_(:,1),:);
		hold('on'); hd1 = plot3(tarNodeCoord(:,1), tarNodeCoord(:,2), tarNodeCoord(:,3), 'x', ...
			'color', [153 153 153]/255, 'LineWidth', 3, 'MarkerSize', 15);		
    end
	if ~isempty(find(0==(boundingBox_(2,:)-boundingBox_(1,:)))) %%2D
		view(2);
	else
		view(3);
	end
	axis(gca, 'equal'); 
	axis(gca, 'tight');
	axis(gca, 'off');
end

function StressTopologyAnalysis()
	%%1. Identify Candidate Elements for Topology Analysis
	IdentifyCandidateElements();
	
	%%2. Locate Degenerate Points	
	LocateDegeneratePoints();	
	
	%%3. Compute Topological Skeleton
	ComputeTopologicalSkeletons();
end

function ComputeTopologicalSkeletons()
	global degePts_;
	global nodeCoords_ eNodMat_ meshTypeMap_;
	global localStressFieldPerEle_ globalStressFieldPerEle_;
	global approxiInterp_;
	global tracingSche_; 
	
	prescribedTracingSche = tracingSche_;
	tracingSche_ = 'Euler';
	numDegPts = numel(degePts_);
	if 0==numDegPts, return; end
	
	%%1. Derivatives of Cartesian Stresses at Degenerate Points w.r.t. Cartesian Coordinates
	for ii=1:numDegPts
		iEle = degePts_(ii).eleIndex;
		iEleType = meshTypeMap_(iEle);
		iPara = degePts_(ii).paraCoord;
		[dNds, dNdt] = DeShapeFunction(iPara, iEleType);
		switch iEleType
			case 'T3'
				iEleNodes = eNodMat_(iEle,1:3);
				iEleCoords = nodeCoords_(iEleNodes,:);
				[iLocalFrame, ~, t1, t2] = ComputeLocalFrameAtGivenPosition(iPara, iEleCoords, iEleType);
				if approxiInterp_
					iEleStressTensorLocal = localStressFieldPerEle_{iEle};
				else
					iEleStressTensorGlobal = globalStressFieldPerEle_{iEle};
				end
			case 'Q4'
				iEleNodes = eNodMat_(iEle,1:4);
				iEleCoords = nodeCoords_(iEleNodes,:);
				[iLocalFrame, ~, t1, t2] = ComputeLocalFrameAtGivenPosition(iPara, iEleCoords, iEleType);
				if approxiInterp_
					iEleStressTensorLocal = localStressFieldPerEle_{iEle};
				else
					iEleStressTensorGlobal = globalStressFieldPerEle_{iEle};
				end				
			case 'T6'
				iEleNodes = eNodMat_(iEle,1:6);
				iEleCoords = nodeCoords_(iEleNodes,:);
				[iLocalFrame, ~, t1, t2] = ComputeLocalFrameAtGivenPosition(iPara, iEleCoords, iEleType);
				if approxiInterp_
					iEleStressTensorLocal = localStressFieldPerEle_{iEle};
				else
					iEleStressTensorGlobal = globalStressFieldPerEle_{iEle};
				end				
			case 'Q8'
				iEleNodes = eNodMat_(iEle,1:8);
				iEleCoords = nodeCoords_(iEleNodes,:);
				[iLocalFrame, ~, t1, t2] = ComputeLocalFrameAtGivenPosition(iPara, iEleCoords, iEleType);
				if approxiInterp_
					iEleStressTensorLocal = localStressFieldPerEle_{iEle};
				else
					iEleStressTensorGlobal = globalStressFieldPerEle_{iEle};
				end				
		end
		e1 = iLocalFrame(:,1);
		e2 = iLocalFrame(:,2);
		JsurfaceMapping = [sum(t1(:).*e1) sum(t2(:).*e1); sum(t1(:).*e2) sum(t2(:).*e2)];
		dNdxy = inv(JsurfaceMapping)' * [dNds; dNdt];
		if approxiInterp_
			dStressdxLocal = dNdxy(1,:) * iEleStressTensorLocal;
			dStressdyLocal = dNdxy(2,:) * iEleStressTensorLocal;
			df1dx = (dStressdxLocal(1)-dStressdxLocal(2))/2; a = df1dx;
			df1dy = (dStressdyLocal(1)-dStressdyLocal(2))/2; b = df1dy;
			df2dx = dStressdxLocal(3); c = df2dx;
			df2dy = dStressdyLocal(3); d = df2dy;
		else
			dStressdxGlobal = dNdxy(1,:) * iEleStressTensorGlobal;
			dStressdxLocal = GlobalFrame2Local_StressTensor(dStressdxGlobal, iLocalFrame, 1); %% D_sigma_xx_D_xx, D_sigma_yy_D_xx, D_sigma_xy_D_xx
			dStressdyGlobal = dNdxy(2,:) * iEleStressTensorGlobal;
			dStressdyLocal = GlobalFrame2Local_StressTensor(dStressdyGlobal, iLocalFrame, 1); %% D_sigma_xx_D_yy, D_sigma_yy_D_yy, D_sigma_xy_D_yy
			df1dx = (dStressdxLocal(1)-dStressdxLocal(2))/2; a = df1dx;
			df1dy = (dStressdyLocal(1)-dStressdyLocal(2))/2; b = df1dy;
			df2dx = dStressdxLocal(3); c = df2dx;
			df2dy = dStressdyLocal(3); d = df2dy;
		end
		degePts_(ii).df12dxy = [a b; c d];
		degePts_(ii).delta = a*d - b*c; %% Negative -> Trisector; Positive -> Wedge
		rawRoots = roots([d (c+2*b) (2*a-d) -c]);
		degePts_(ii).tangentList = rawRoots(0==imag(rawRoots));
	end
	
	%%2. Get Topological Skeleton
	for ii=1:numDegPts
		degePts_(ii).majorSkeletons = repmat(degePts_(ii).majorSkeletons, numel(degePts_(ii).tangentList), 1);
		degePts_(ii).minorSkeletons = repmat(degePts_(ii).minorSkeletons, numel(degePts_(ii).tangentList), 1);
		for jj=1:numel(degePts_(ii).tangentList)
			ijSeed = [degePts_(ii).eleIndex degePts_(ii).paraCoord];
			iniDir = [1 degePts_(ii).tangentList(jj)]; iniDir = iniDir/norm(iniDir);
			ijMajorPSL = CreatePrincipalStressLine(ijSeed, 'MAJOR', iniDir);
			ijMinorPSL = CreatePrincipalStressLine(ijSeed, 'MINOR', iniDir);
			degePts_(ii).majorSkeletons(jj) = ijMajorPSL;
			degePts_(ii).minorSkeletons(jj) = ijMinorPSL;	
		end
	end
	tracingSche_ = prescribedTracingSche;
end

function val = DegeneratePointStruct()
	PSL = PrincipalStressLineStruct();
	val = struct(	...
		'eleIndex',							0,	...
		'paraCoord',						[],	...		
		'phyCoord',							[],	...
		'cartesianStress',					[], ...
		'principalStress',					[],	...
		'directDegenerancyExtentMetric', 	[], ...
		'tangentList',						[],	...
		'df12dxy',							[],	...
		'abcd',								[],	...
		'delta',							0,	...
		'majorSkeletons',					PSL,...
		'minorSkeletons',					PSL ...
	);
end

function LocateDegeneratePoints()
	global nodeCoords_ eNodMat_ meshTypeMap_
	global candidateElementsIncDegePts_ degePts_
	global refScalingSize_;
	
	numCandiEles = numel(candidateElementsIncDegePts_);
	distilledCandidateElements = [];
	for ii=1:numCandiEles
		iEle = candidateElementsIncDegePts_(ii);
		iEleType = meshTypeMap_(iEle);
		switch iEleType
			case 'T3'
				[paras, isDegePt] = NewtonRhapson_LocatingDegPt(iEle, [0 0]);
				if isDegePt
					distilledCandidateElements(end+1,1:3) = [iEle paras(:)'];
				end				
			case 'Q4'
				[paras, isDegePt] = NewtonRhapson_LocatingDegPt(iEle, [0 0]);
				if isDegePt
					distilledCandidateElements(end+1,1:3) = [iEle paras(:)'];
				end						
			case 'T6'
				tmpDistilledCandiEle = [];
				[paras1, isDegePt1] = NewtonRhapson_LocatingDegPt(iEle, [0 0], 'T6_1');
				if isDegePt1
					tmpDistilledCandiEle(end+1,1:3) = [iEle paras1(:)'];
				end
				[paras2, isDegePt2] = NewtonRhapson_LocatingDegPt(iEle, [0 0], 'T6_2');
				if isDegePt2
					if isempty(tmpDistilledCandiEle)
						tmpDistilledCandiEle(end+1,1:3) = [iEle paras2(:)'];
					else
						if min(vecnorm(paras2(:)' - tmpDistilledCandiEle(:,2:3),2,2)) > 0.1 %%Additional Degenerate Point
							tmpDistilledCandiEle(end+1,1:3) = [iEle paras2(:)'];
						end
					end				
				end
				[paras3, isDegePt3] = NewtonRhapson_LocatingDegPt(iEle, [0 0], 'T6_3');
				if isDegePt3
					if isempty(tmpDistilledCandiEle)
						tmpDistilledCandiEle(end+1,1:3) = [iEle paras3(:)'];
					else
						if min(vecnorm(paras3(:)' - tmpDistilledCandiEle(:,2:3),2,2)) > 0.1 %%Additional Degenerate Points
							tmpDistilledCandiEle(end+1,1:3) = [iEle paras3(:)'];
						end
					end				
				end
				[paras4, isDegePt4] = NewtonRhapson_LocatingDegPt(iEle, [0 0], 'T6_4');
				if isDegePt4
					if isempty(tmpDistilledCandiEle)
						tmpDistilledCandiEle(end+1,1:3) = [iEle paras4(:)'];
					else
						if min(vecnorm(paras4(:)' - tmpDistilledCandiEle(:,2:3),2,2)) > 0.1 %%Additional Degenerate Points
							tmpDistilledCandiEle(end+1,1:3) = [iEle paras4(:)'];
						end						
					end
				end
				if ~isempty(tmpDistilledCandiEle)
					distilledCandidateElements(end+1:end+size(tmpDistilledCandiEle,1),:) = tmpDistilledCandiEle;
				end			
			case 'Q8'
				tmpDistilledCandiEle = [];
				[paras1, isDegePt1] = NewtonRhapson_LocatingDegPt(iEle, [0 0], 'Q8_1');
				if isDegePt1
					tmpDistilledCandiEle(end+1,1:3) = [iEle paras1(:)'];
				end
				[paras2, isDegePt2] = NewtonRhapson_LocatingDegPt(iEle, [0 0], 'Q8_2');
				if isDegePt2
					if isempty(tmpDistilledCandiEle)
						tmpDistilledCandiEle(end+1,1:3) = [iEle paras2(:)'];
					else
						if min(vecnorm(paras2(:)' - tmpDistilledCandiEle(:,2:3),2,2)) > 0.1 %%Additional Degenerate Point
							tmpDistilledCandiEle(end+1,1:3) = [iEle paras2(:)'];
						end
					end				
				end
				[paras3, isDegePt3] = NewtonRhapson_LocatingDegPt(iEle, [0 0], 'Q8_3');
				if isDegePt3
					if isempty(tmpDistilledCandiEle)
						tmpDistilledCandiEle(end+1,1:3) = [iEle paras3(:)'];
					else
						if min(vecnorm(paras3(:)' - tmpDistilledCandiEle(:,2:3),2,2)) > 0.1 %%Additional Degenerate Points
							tmpDistilledCandiEle(end+1,1:3) = [iEle paras3(:)'];
						end
					end				
				end
				[paras4, isDegePt4] = NewtonRhapson_LocatingDegPt(iEle, [0 0], 'Q8_4');
				if isDegePt4
					if isempty(tmpDistilledCandiEle)
						tmpDistilledCandiEle(end+1,1:3) = [iEle paras4(:)'];
					else
						if min(vecnorm(paras4(:)' - tmpDistilledCandiEle(:,2:3),2,2)) > 0.1 %%Additional Degenerate Points
							tmpDistilledCandiEle(end+1,1:3) = [iEle paras4(:)'];
						end						
					end
				end
				if ~isempty(tmpDistilledCandiEle)
					distilledCandidateElements(end+1:end+size(tmpDistilledCandiEle,1),:) = tmpDistilledCandiEle;
				end				
		end
	end
	if isempty(distilledCandidateElements), degePts_ = []; return; end
	
	%% Merge Degenerate Point(s) Located at the element Edge
	numDistilledCandiEles = size(distilledCandidateElements, 1);
	validDegeneratePointsMap = false(numDistilledCandiEles,1);
	refCoords = [];
	for ii=1:numDistilledCandiEles
		iEle = distilledCandidateElements(ii,1);
		iPara = distilledCandidateElements(ii,2:3);
		iEleType = meshTypeMap_(iEle);
		switch iEleType
			case 'T3'
				iEleNodes = eNodMat_(iEle,1:3);
			case 'Q4'
				iEleNodes = eNodMat_(iEle,1:4);
			case 'T6'
				iEleNodes = eNodMat_(iEle,1:6);
			case 'Q8'
				iEleNodes = eNodMat_(iEle,1:8);	
		end
		SF = ShapeFunction(iPara, iEleType);
		iEleNodeCoords = nodeCoords_(iEleNodes,:);
		iPhyCoord = SF * iEleNodeCoords;
		if 1==ii
			refCoords = iPhyCoord;
			validDegeneratePointsMap(ii) = true;
		else
			if min(vecnorm(iPhyCoord-refCoords, 2, 2)) > 1.0e-2*refScalingSize_ %%Not the Same Degenerate Point
				refCoords(end+1,:) = iPhyCoord;
				validDegeneratePointsMap(ii) = true;
			end
		end
	end
	distilledCandidateElements = distilledCandidateElements(find(validDegeneratePointsMap),:);
	
	%%Setup Degenerate Point Structure
	numDegePts = size(distilledCandidateElements,1);
	degePts_ = DegeneratePointStruct();
	degePts_ = repmat(degePts_, numDegePts, 1);
	for ii=1:numDegePts
		degePts_(ii).eleIndex = distilledCandidateElements(ii,1);
		degePts_(ii).paraCoord = distilledCandidateElements(ii,2:3);
		degePts_(ii).phyCoord = refCoords(ii,:);
	end
end

function IdentifyCandidateElements()
	global numEles_ meshTypeMap_ boundaryElements_ voidEleMap_;
	global localStressFieldPerEle_;
	global candidateElementsIncDegePts_;
	
	excludedElementMap = true(numEles_,1);
	
	for ii=1:numEles_
		if voidEleMap_(ii), continue; end
		eleStress = localStressFieldPerEle_{ii};
		switch meshTypeMap_(ii)
			case {'T3', 'Q4', 'T6'}
				opt = DegenrationMeasure_SignConsistency(eleStress);
			case 'Q8'
				stressAtCtr = ShapeFunction([0 0], 'Q8') * eleStress;		
				opt = DegenrationMeasure_BernsteinBoundingQuad([eleStress; stressAtCtr]);
		end
		excludedElementMap(ii) = opt;
	end
	candidateElementsIncDegePts_ = find(false==excludedElementMap);
	candidateElementsIncDegePts_ = setdiff(candidateElementsIncDegePts_, boundaryElements_);
end

function opt = DegenrationMeasure_SignConsistency(eleStress)
	f1 = (eleStress(:,1)-eleStress(:,2))/2;
	f2 = eleStress(:,3);
	discriminantf1 = (min(f1) > 0) || (max(f1) < 0);
	discriminantf2 = (min(f2) > 0) || (max(f2) < 0);
	opt = (discriminantf1 || discriminantf2);
end

function opt = DegenrationMeasure_BernsteinBoundingQuad(eleStress)
	f1 = (eleStress(:,1)-eleStress(:,2))/2;
	f2 = eleStress(:,3);	
	reOrdering = [1 8 4 5 9 7 2 6 3]';
	f1 = f1(reOrdering);
	f2 = f2(reOrdering);
	S1 = reshape(f1, 3, 3)';
	S2 = reshape(f2, 3, 3)';

	C = [1 0 0; -0.5 2 -0.5; 0 0 1];
	B1 = C * S1 * C';
	B2 = C * S2 * C';

	f1min = min(B1(:)); f1max = max(B1(:));
	f2min = min(B2(:)); f2max = max(B2(:));

	discriminantf1 = (f1min > 0) || (f1max < 0);
	discriminantf2 = (f2min > 0) || (f2max < 0);

	opt = (discriminantf1 || discriminantf2);
end

function [paras, isDegePt] = NewtonRhapson_LocatingDegPt(iEle, target, varargin)
	global nodeCoords_ eNodMat_ meshTypeMap_;
	global localStressFieldPerEle_ globalStressFieldPerEle_;
	global approxiInterp_;
	
	target = target(:)';
	
	maxIts = 30;
	tol_r  = 1e-5;                      % residual tolerance (stress units)
	tol_x  = 1e-12;                      % step tolerance in param space
	lambda0 = 0;                         % LM damping (0 = pure Newton)
	project_to_domain = true;
	iEleType = meshTypeMap_(iEle);	
	switch iEleType
		case 'T3'
			natCoordScope = [0 0; 1 0; 0 1];
			if approxiInterp_
				eleStressTensorsLocal = localStressFieldPerEle_{iEle};
			else
				vtxCoords = nodeCoords_(eNodMat_(iEle,1:3),:);
                eleStressTensorsGlobal = globalStressFieldPerEle_{iEle};
			end
		case 'Q4'
			natCoordScope = [-1 -1; 1 -1; 1 1; -1 1];
			if approxiInterp_
				eleStressTensorsLocal = localStressFieldPerEle_{iEle};
			else
				vtxCoords = nodeCoords_(eNodMat_(iEle,1:4),:);
                eleStressTensorsGlobal = globalStressFieldPerEle_{iEle};
			end
		case 'T6'
			natCoordScope = [0 0; 1 0; 0 1; 0.5	0; 0.5 0.5; 0 0.5];
			if 3==nargin
				subT6 = varargin{1};
				switch subT6
					case 'T6_1'
						natCoordScope = [0 0; 0.5 0; 0 0.5];
					case 'T6_2'
						natCoordScope = [0.5 0; 1 0; 0.5 0.5];
					case 'T6_3'
						natCoordScope = [0 0.5; 0.5 0.5; 0 1];
					case 'T6_4'
						natCoordScope = [0.5 0; 0.5 0.5; 0 0.5];
					otherwise
						error('Un-supported Division!');
				end
			end
			if approxiInterp_
				eleStressTensorsLocal = localStressFieldPerEle_{iEle};
			else
				vtxCoords = nodeCoords_(eNodMat_(iEle,1:6),:);
				eleStressTensorsGlobal = globalStressFieldPerEle_{iEle};
			end
		case 'Q8'
			natCoordScope = [-1 -1; 1 -1; 1 1; -1 1; 0 -1; 1 0; 0 1; -1 0];
			if 3==nargin
				subQ8 = varargin{1};
				switch subQ8
					case 'Q8_1'
						natCoordScope = [-1 -1; 0 -1; 0 0; -1 0];
					case 'Q8_2'
						natCoordScope = [0 -1; 1 -1; 1 0; 0 0];
					case 'Q8_3'
						natCoordScope = [0 0; 1 0; 1 1; 0 1];
					case 'Q8_4'
						natCoordScope = [-1 0; 0 0; 0 1; -1 1];
					otherwise
						error('Un-supported Division!');
				end				
			end
			if approxiInterp_
				eleStressTensorsLocal = localStressFieldPerEle_{iEle};
			else
				vtxCoords = nodeCoords_(eNodMat_(iEle,1:8),:);
				eleStressTensorsGlobal = globalStressFieldPerEle_{iEle};
			end
		otherwise
			error('Un-supported Element Type!');
	end
	bBox = [min(natCoordScope,[],1); max(natCoordScope,[],1)];
	ctr = sum(natCoordScope,1)/size(natCoordScope,1);
	s = ctr(1); t = ctr(2);
	
	for it = 1:maxIts
		% Shape functions and first derivatives
		N = ShapeFunction([s t], iEleType);        % 1xMe
		[dNds, dNdt] = DeShapeFunction([s t], iEleType); % 1xMe
	
		% q(s,t) and residual r = q - target
		if approxiInterp_
			%% q   = N    * vtxVals;                % 1x2
			q = N * [(eleStressTensorsLocal(:,1)-eleStressTensorsLocal(:,2))/2 eleStressTensorsLocal(:,3)];
		else
			iStressGlobal = N * eleStressTensorsGlobal;
			[iLocalFrame, ~, ~, ~] = ComputeLocalFrameAtGivenPosition([s t], vtxCoords, iEleType);
			iStressLocal = GlobalFrame2Local_StressTensor(iStressGlobal, iLocalFrame, 1);
			q = [(iStressLocal(1)-iStressLocal(2))/2 iStressLocal(3)];		
		end
		r   = q - target;                    % 1x2
		rn  = norm(r);
	
		% Converged by residual?			
		if rn <= tol_r
			break;
		end

		% Jacobian J = dq/d(s,t) arranged as 2x2:
		% J = [ df1/ds  df1/dt ; df2/ds  df2/dt ]
		if approxiInterp_
			dqs = dNds * [(eleStressTensorsLocal(:,1)-eleStressTensorsLocal(:,2))/2 eleStressTensorsLocal(:,3)];  % 1x2  (df1/ds, df2/ds)
			dqt = dNdt * [(eleStressTensorsLocal(:,1)-eleStressTensorsLocal(:,2))/2 eleStressTensorsLocal(:,3)];  % 1x2  (df1/dt, df2/dt)				
		else
			dqsTmpGlobal = dNds * eleStressTensorsGlobal;
			dqsTmpLocal = GlobalFrame2Local_StressTensor(dqsTmpGlobal, iLocalFrame, 1);
			dqs = [(dqsTmpLocal(1)-dqsTmpLocal(2))/2 dqsTmpLocal(3)];
			dqtTmpGlobal = dNdt * eleStressTensorsGlobal;
			dqtTmpLocal = GlobalFrame2Local_StressTensor(dqtTmpGlobal, iLocalFrame, 1);
			dqt = [(dqtTmpLocal(1)-dqtTmpLocal(2))/2 dqtTmpLocal(3)];			
		end
		J   = [ dqs(1) dqt(1) ; dqs(2) dqt(2) ];
	
		% If J is near singular, add LM damping
		% (you can also inspect det(J) or cond(J) to adapt lambda)
		lambda = lambda0;
		H = J; rhs = -r(:);
		if lambda > 0
			H = (J.'*J + lambda*eye(2));    % LM on normal equations
			rhs = -(J.'*r(:));
		end
	
		% Solve for step (Newton or LM)
		% Prefer the square solve when lambda==0 to avoid normal-equations squaring
		if lambda==0
			% robust 2x2 solve with small Tikhonov if singular
			if rcond(J) < 1e-12
				% fallback to damped normal equations
				Hne = (J.'*J + 1e-12*eye(2));
				dx  = Hne \ (-(J.'*r(:)));
			else
				dx = J \ rhs;               % 2x1
			end
		else
			dx = H \ rhs;
		end
	
		% Backtracking line search to ensure decrease of phi = 0.5 ||r||^2
		phi0 = 0.5*rn*rn;
		alpha = 1.0;
		for bt = 1:8
			s_try = s + alpha*dx(1);
			t_try = t + alpha*dx(2);
	
			% Optional: project to reference domain to avoid wild steps
			if project_to_domain
				if strcmp(iEleType, 'T6') || strcmp(iEleType, 'T3')
					[s_try, t_try] = ProjectToTriangle([s_try, t_try], natCoordScope(1,:), natCoordScope(2,:), natCoordScope(3,:));
				elseif strcmp(iEleType, 'Q4') || strcmp(iEleType, 'Q8')	
					s_try = max(bBox(1,1),min(bBox(2,1),s_try));
					t_try = max(bBox(1,2),min(bBox(2,2),t_try));
				else
					error('Un-supported Element Type!');
				end
			end
	
			N_try = ShapeFunction([s_try t_try], iEleType);
			if approxiInterp_
				q_try = N_try * [(eleStressTensorsLocal(:,1)-eleStressTensorsLocal(:,2))/2 eleStressTensorsLocal(:,3)];
			else
				iStressGlobal_try = N_try * eleStressTensorsGlobal;
				[iLocalFrame, ~, ~, ~] = ComputeLocalFrameAtGivenPosition([s_try t_try], vtxCoords, iEleType);
				iStressLocal = GlobalFrame2Local_StressTensor(iStressGlobal_try, iLocalFrame, 1);
				q_try = [(iStressLocal(1)-iStressLocal(2))/2 iStressLocal(3)];				
			end
			r_try = q_try - target;
			phi   = 0.5*dot(r_try,r_try);
	
			if phi <= phi0 * (1 - 1e-4*alpha)     % Armijo condition
				s = s_try; t = t_try;
				break;
			else
				alpha = 0.5*alpha;
			end
		end
	
		% Converged by step size?
		if norm(alpha*dx) <= tol_x
			break;
		end
    end
	paras  = [s t];
	inside = CheckWhetherGivenNaturalCoordinatesWithinDomain(paras, iEleType, 1.0e-12);
	isDegePt = inside && it < maxIts;
end

