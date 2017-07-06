function [MaskMap,GMEPar] = CalGME(inputPath)
	if exist(inputPath, 'dir') ~= 7
		error('Wrong input path %s.\n', inputPath);
	end

	flowDir = fullfile(inputPath,'flow');
	matchDir = fullfile(inputPath,'match');
	imgDir = inputPath;

	dbstop if error;
	if exist(flowDir, 'dir') ~= 7
		error('Wrong flow input path %s.\n', flowDir);
	end
	if exist(matchDir, 'dir') ~= 7
		error('Wrong match input path %s.\n', matchDir);
	end
	if exist(imgDir, 'dir') ~= 7
		error('Wrong image input path %s.\n', imgDir);
	end

	flowFiles = dir(fullfile(flowDir, '*.flo'));
	flowNames = sort({flowFiles.name});
	flowNum = numel(flowNames);
	matchFiles = dir(fullfile(matchDir, '*.txt'));
	matchNames = sort({matchFiles.name});
	imgFiles = dir(fullfile(imgDir,'*.png'));
	imgNames = sort({imgFiles.name});

	for ii = 1:flowNum
		motVect = readMatch(fullfile(flowDir, flowNames{ii}), fullfile(matchDir, matchNames{ii}));
		[height,width,nBands] = size(motVect);
		if nBands ~= 2
			error('readMatch: Motion Vectors must have two bands.\n');
		else
			fprintf('Cureent CalNo.%d\n', ii);
		end

		u = motVect(:,:,1);	% the vector of horizontal direction
		v = motVect(:,:,2);	% the vector of vertical direction

		windowSize = 3;
		u = medfilt2(u, [windowSize, windowSize]);
		v = medfilt2(v, [windowSize, windowSize]);

		if ii == 1
			preMask = ones(height,width);
			% maskBuffer = zeros(height,width);
			MaskMap = zeros(height,width,flowNum);
			GMEPar = zeros(6,flowNum);
			[A1,A2,curMask] = CalGMEPar6(u,v,preMask);
			GMEPar(1:3,ii) = A1;
			GMEPar(4:6,ii) = A2;
		else
			[hRatio,vRatio] = CalRatio(fullfile(imgDir, imgNames{ii-1}),fullfile(imgDir, imgNames{ii}));

			if (hRatio >= 0.5 || vRatio >= 0.5)
				[A1,A2,curMask] = CalGMEPar2(u,v,preMask);
				GMEPar(1:3,ii) = [1;0;A1];
				GMEPar(4:6,ii) = [0;1;A2];
			else
				[A1,A2,curMask] = CalGMEPar6(u,v,preMask);
				GMEPar(1:3,ii) = A1;
				GMEPar(4:6,ii) = A2;
			end
		end
		% preMask = CalObjectMask(curMask,maskBuffer,u,v);
		preMask = curMask;
		MaskMap(:,:,ii) = curMask;
	end
	% MaskMap = vesselfilter(MaskMap);
end

%% CalGMEPar2: function description
function [A1,A2,bgPointFlag] = CalGMEPar2(u,v,bgPointFlag);
	mv = sqrt(u.^2 + v.^2);
	mvMean = mean2(mv);
	sigma = std2(mv);

	noiseMap = and(mv <= (mvMean + sigma), mv >= (mvMean - sigma));
	bgPointFlag = bgPointFlag.*noiseMap;

	clear mv mvMean sigma;

	reCalParamCnt = 1;
	diffVar = 100;
	parDiff = 1;
	sigma1 = 0.1;
	breakFlag = 0;

	while(parDiff >= sigma1 && reCalParamCnt <= 32 && breakFlag == 0)
		if reCalParamCnt == 1
			A1 = mean2(u);
			A2 = mean2(v);
			a1 = 0;
			a2 = 0;
		else
			a1 = A1;
			a2 = A2;
		end
		[A1, A2] = CalPar2(A1,A2,u,v,bgPointFlag);
		parDiff = sqrt(((a1 - A1)^2 + (a2 - A2)^2).*0.5);

		[bgPointFlag,diffVar,breakFlag] = updateFlag2(A1,A2,u,v,bgPointFlag,diffVar);

		% fprintf('reCalParam: %d times,parDiff: %d, parDiff2: %d.\n', reCalParamCnt,parDiff,parDiff2);
		reCalParamCnt = reCalParamCnt + 1;
	end
	% bgPointFlag = bgPointFlag + noiseMap;
end

function [A1,A2,bgPointFlag] = CalGMEPar6(u,v,bgPointFlag);
	mv = sqrt(u.^2 + v.^2);
	mvMean = mean2(mv);
	sigma = std2(mv);

	noiseMap = and(mv <= (mvMean + sigma), mv >= (mvMean - sigma));
	bgPointFlag = bgPointFlag.*noiseMap;

	clear mv mvMean sigma;

	reCalParamCnt = 1;
	diffVar = 100;
	parDiff1 = 1;
	parDiff2 = 1;
	sigma1 = 0.1;
	sigma2 = 0.001;
	breakFlag = 0;

	while(parDiff1 >= sigma1 && parDiff2 >= sigma2 && reCalParamCnt <= 32 && breakFlag == 0)
		if reCalParamCnt == 1
			A1 = [1; 0; mean2(u)];
			A2 = [0; 1; mean2(v)];
			a1 = [0; 0; 0];
			a2 = [0; 0; 0];
		else
			a1 = A1;
			a2 = A2;
		end
		[A1, A2] = CalPar6(A1,A2,u,v,bgPointFlag);
		parDiff1 = sqrt(((a1(3) - A1(3))^2 + (a2(3) - A2(3))^2)*0.5);
		parDiff2 = sqrt((sum((a1(1:2) - A1(1:2)).^2) + sum((a2(1:2) - A2(1:2)).^2))*0.25);

		[bgPointFlag,diffVar,breakFlag] = updateFlag6(A1,A2,u,v,bgPointFlag,diffVar);

		% fprintf('reCalParam: %d times,parDiff1: %d, parDiff2: %d.\n', reCalParamCnt,parDiff1,parDiff2);
		reCalParamCnt = reCalParamCnt + 1;
	end
	% bgPointFlag = bgPointFlag + noiseMap;
end

%% CalPar2: function description
function [A1,A2] = CalPar2(A1,A2,u,v,bgPointFlag)
	[height, width] = size(bgPointFlag);
	[refPosX, refPosY] = meshgrid(1:width, 1:height);
	refPosX = refPosX - fix(width/2);
	refPosY = fix(height/2) - refPosY;

	curPosX = refPosX + u;
	curPosY = refPosY + v;

	% curPosX = (refPosX + u).*bgPointFlag;
	% curPosY = (refPosY + v).*bgPointFlag;

	% refPosX = refPosX.*bgPointFlag;
	% refPosY = refPosY.*bgPointFlag;

	alpha = refPosX + A1;
	b0 = -sum(reshape((alpha - curPosX).*bgPointFlag,[],1));

	beta = refPosY + A2;
	b1 = -sum(reshape((beta - curPosY).*bgPointFlag,[],1));

	delta1 = b0./sum(bgPointFlag(:));
	delta2 = b1./sum(bgPointFlag(:));

	A1 = A1 + delta1;
	A2 = A2 + delta2;
end

function [A1,A2] = CalPar6(A1,A2,u,v,bgPointFlag)
	[height, width] = size(bgPointFlag);

	a10 = A1(1);
	a20 = A1(2);
	a50 = A1(3);
	a30 = A2(1);
	a40 = A2(2);
	a60 = A2(3); 
	
	[refPosX, refPosY] = meshgrid(1:width, 1:height);
	refPosX = refPosX - fix(width/2);
	refPosY = fix(height/2) - refPosY;

	curPosX = refPosX + u;
	curPosY = refPosY + v;

	% curPosX = (refPosX + u).*bgPointFlag;
	% curPosY = (refPosY + v).*bgPointFlag;

	% refPosX = refPosX.*bgPointFlag;
	% refPosY = refPosY.*bgPointFlag;

	H11 = sum(reshape((refPosX.*refPosX).*bgPointFlag,[],1));
	H12 = sum(reshape((refPosX.*refPosY).*bgPointFlag,[],1));
	H13 = sum(reshape(refPosX.*bgPointFlag,[],1));

	H21 = sum(reshape((refPosY.*refPosX).*bgPointFlag,[],1));
	H22 = sum(reshape((refPosY.*refPosY).*bgPointFlag,[],1));
	H23 = sum(reshape(refPosY.*bgPointFlag,[],1));

	H31 = sum(reshape(refPosX.*bgPointFlag,[],1));
	H32 = sum(reshape(refPosY.*bgPointFlag,[],1));
	% H33 = height*width;
	H33 = sum(bgPointFlag(:));

	% b10 = -sum(reshape(a10*refPosX.*refPosX + a20*refPosY.*refPosX + a50*ones(height,width).*refPosX - curPosX.*refPosX,[],1));
	alpha = a10*refPosX + a20*refPosY + a50;
	b10 = -sum(reshape((alpha.*refPosX - curPosX.*refPosX).*bgPointFlag,[],1));
	b20 = -sum(reshape((alpha.*refPosY - curPosX.*refPosY).*bgPointFlag,[],1));
	b30 = -sum(reshape((alpha - curPosX).*bgPointFlag,[],1));

	beta = a30*refPosX + a40*refPosY + a60;
	b11 = -sum(reshape((beta.*refPosX - curPosY.*refPosX).*bgPointFlag,[],1));
	b21 = -sum(reshape((beta.*refPosY - curPosY.*refPosY).*bgPointFlag,[],1));
	b31 = -sum(reshape((beta - curPosY).*bgPointFlag,[],1));

	H = [H11, H12, H13; H21, H22, H23; H31, H32, H33];
	B1 = [b10; b20; b30];
	B2 = [b11; b21; b31];

	delta1 = inv(H)*B1;
	delta2 = inv(H)*B2;
	A1 = [a10; a20; a50] + delta1;
	A2 = [a30; a40; a60] + delta2;
end

function [bgPointFlag,diffVar,breakFlag] = updateFlag2(A1,A2,u,v,bgPointFlag,diffVar)
	[height,width] = size(bgPointFlag);
	[refPosX,refPosY] = meshgrid(1:width,1:height);
	refPosX = refPosX - fix(width/2);
	refPosY = fix(height/2) - refPosY;

	newPosX = refPosX + A1;
	newPosY = refPosY + A2;
	curPosX = refPosX + u;
	curPosY = refPosY + v;

	% posDiff = sqrt((newPosX - curPosX).^2 + (newPosY - curPosY).^2);
	% posDiffTmp = sort(setdiff(posDiff,0));
	% thre = posDiffTmp(ceil(numel(posDiffTmp)*0.9));
	posDiff = abs(newPosX - curPosX) + abs(newPosY - curPosY);

	if diffVar(end) <= var(posDiff(:)) && var(posDiff(:)) < 1
		breakFlag = 1;
	else
		breakFlag = 0;
		diffVar = [diffVar, var(posDiff(:))];
		% bgPointFlag = (posDiff <= thre).*bgPointFlag;
		bgPointFlag = errorHist(posDiff);
	end
end

function [bgPointFlag,diffVar,breakFlag] = updateFlag6(A1,A2,u,v,bgPointFlag,diffVar)
	[height,width] = size(bgPointFlag);
	[refPosX,refPosY] = meshgrid(1:width,1:height);
	refPosX = refPosX - fix(width/2);
	refPosY = fix(height/2) - refPosY;

	newPosX = A1(1)*refPosX + A1(2)*refPosY + A1(3);
	newPosY = A2(1)*refPosX + A2(2)*refPosY + A2(3);
	curPosX = refPosX + u;
	curPosY = refPosY + v;

	% posDiff = sqrt((newPosX - curPosX).^2 + (newPosY - curPosY).^2);
	% posDiffTmp = sort(setdiff(posDiff,0));
	% thre = posDiffTmp(ceil(numel(posDiffTmp)*0.9));
	posDiff = abs(newPosX - curPosX) + abs(newPosY - curPosY);

	if diffVar(end) <= var(posDiff(:)) && var(posDiff(:)) < 1
		breakFlag = 1;
	else
		breakFlag = 0;
		diffVar = [diffVar, var(posDiff(:))];
		% bgPointFlag = (posDiff <= thre).*bgPointFlag;
		bgPointFlag = errorHist(posDiff);
	end
end

function weightMap =  errorHist(posDiff)
	[height,width] = size(posDiff);
	err = zeros(height,width);
	[h,n] = hist(posDiff(:));
	N = numel(posDiff);
	cnt = numel(n);
	errMax = n(find(h == max(h),1));
	errTmp = posDiff - errMax;

	w = h - N;
	w = w/mean2(w);

	for ii = 1:cnt
		if ii == 1
			err = err + errTmp.*(posDiff < (n(ii) + n(ii + 1))/2).*w(ii);
		elseif ii == cnt
			err = err + errTmp.*(posDiff >= (n(ii - 1) + n(ii))/2).*w(ii);
		else
			err = err + errTmp.*and(posDiff >= (n(ii - 1) + n(ii))/2, posDiff < (n(ii) + n(ii + 1))/2).*w(ii);
		end
	end

	weightMap = (ones(height,width) - err.^2).^2.*(abs(err) < 1);
end

function [hRatio,vRatio] = CalRatio(orgImgPath,curImgPath)
	orgImg = imread(orgImgPath);
	orgImg = double(rgb2gray(orgImg));
	orgImg = imfilter(orgImg,fspecial('gaussian',[5 5],3));
	curImg = imread(curImgPath);
	curImg = double(rgb2gray(curImg));
	curImg = imfilter(curImg,fspecial('gaussian',[5 5],3));
	tempImg = curImg - orgImg;
	[xgrdImg,ygrdImg] = gradient(orgImg);
	hSTGS = tempImg./xgrdImg;
	vSTGS = tempImg./ygrdImg;
	hRatio = defMinority(hSTGS,'hSTGS');
	vRatio = defMinority(vSTGS,'vSTGS');
end

%% defMinority: function description
function ratio = defMinority(STGS,name)
	[height,width] = size(STGS);
	negMap = STGS < 0;
	posMap = STGS > 0;
	negNum = numel(find(negMap == 1));
	posNum = numel(find(posMap == 1));
	minNum = min([negNum;posNum]);
	switch minNum
		case negNum
			minMap = negMap;
			% disp('Miority is negative map.');
		case posNum
			minMap = posMap;
			% disp('Miority is positive map.');
		otherwise
			% error('MinMap define failed!');
	end
	[xPos,yPos] = meshgrid(1:width,1:height);
	xCenter = sum(reshape(xPos.*minMap,[],1))/minNum;
	yCenter = sum(reshape(yPos.*minMap,[],1))/minNum;
	varience = sum(reshape(((xPos - xCenter.*ones(height,width)).^2 + (yPos - yCenter.*ones(height,width)).^2).*minMap,[],1))/minNum;
	ratio = varience/minNum;
end

% function objectMask = CalObjectMask(initMask,maskBuffer,u,v)
% 	updateBGBuff();
% 	updateBGDiff();

% end

% function updateBGBuff()
% 	bgbuff = bgbuff + initMask;

% end