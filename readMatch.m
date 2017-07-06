function finalMatch = readMatch(flowPath,matchPath)
% function [finalMatch,flagSample] = readMatch(flowPath,matchPath)
	% initialMatch = importdata(matchPath);
	flowFile = readFlowFile(flowPath);
	[height,width,~] = size(flowFile);
	sqsize = 21;	

	u = flowFile(:,:,1);
	v = flowFile(:,:,2);

	% [dataIdx,~] = size(initialMatch);	
	% u = zeros(height,width);
	% v = zeros(height,width);
	% flag = zeros(height,width);

	% for i = 1:dataIdx
	% 	x1 = initialMatch(i,1);
	% 	y1 = initialMatch(i,2);

	% 	x2 = initialMatch(i,3);
	% 	y2 = initialMatch(i,4);

	% 	u(y1-1:y1+1,x1-1:x1+1) = ones(3,3)*(x2 - x1);
	% 	v(y1-1:y1+1,x1-1:x1+1) = ones(3,3)*(y2 - y1);
	% 	flag(y1-1:y1+1,x1-1:x1+1) = ones(3,3);
	% end

	% for i = 1:height
	% 	for j = 1:width
	% 		if flag(i,j) == 0
	% 			% [u(i,j),v(i,j)] = convertFlow2Match(i,j,u,v,flowFile,flag);
	% 			% v(i,j) = convertFlow2Match(i,j,2,u,v,flowFile,flag);
	% 			u(i,j) = flowFile(i,j,1);
	% 			v(i,j) = flowFile(i,j,2);
	% 		end
	% 	end
	% end

	xIdx = fix(width/sqsize);
	xPlus = fix(mod(width,sqsize)/2);
	yIdx = fix(height/sqsize);
	yPlus = fix(mod(height,sqsize)/2);

	uSample = zeros(yIdx,xIdx);
	vSample = zeros(yIdx,xIdx);
	% flagSample = zeros(yIdx,xIdx);

	for i = 1:yIdx
		for j = 1:xIdx
			uTmp = u(yPlus+1+sqsize*(i-1):yPlus+sqsize*i,xPlus+1+sqsize*(j-1):xPlus+sqsize*j);
			vTmp = v(yPlus+1+sqsize*(i-1):yPlus+sqsize*i,xPlus+1+sqsize*(j-1):xPlus+sqsize*j);
			uSample(i,j) = sum(uTmp(:))/(sqsize^2);
			vSample(i,j) = sum(vTmp(:))/(sqsize^2);

			% flagTmp = flag(yPlus+1+5*(i-1):yPlus+5*i,xPlus+1+5*(j-1):xPlus+5*j);
			% uSample(i,j) = round(sum(uTmp(:))/25);
			% vSample(i,j) = round(sum(vTmp(:))/25);
			% flagSample(i,j) = round(sum(flagTmp(:))/25);
		end
	end

	finalMatch(:,:,1) = uSample;
	finalMatch(:,:,2) = vSample;
end

% function [uMatch,vMatch] = convertFlow2Match(y,x,u,v,flowFile,matchFlag)
% 	flo = flowFile(:,:,1);

% 	if flo(y,x) > 0
% 		idxflag = and(flo >= flo(y,x)*0.9, flo <= flo(y,x)*1.1);
% 	else
% 		idxflag = and(flo <= flo(y,x)*0.9, flo >= flo(y,x)*1.1);
% 	end
	
% 	if sum(idxflag(:)) ~= 0
% 		flosum = u.*matchFlag.*idxflag;
% 		uMatch = sum(flosum(:))/numel(find(matchFlag.*idxflag == 1));
% 		flosum = v.*matchFlag.*idxflag;
% 		vMatch = sum(flosum(:))/numel(find(matchFlag.*idxflag == 1));
% 	else
% 		flo = flowFile(:,:,2);
% 		if flo(y,x) > 0
% 			idxflag = and(flo >= flo(y,x)*0.9, flo <= flo(y,x)*1.1);
% 		else
% 			idxflag = and(flo <= flo(y,x)*0.9, flo >= flo(y,x)*1.1);
% 		end
% 		if sum(idxflag(:)) ~= 0
% 			flosum = u.*matchFlag.*idxflag;
% 			uMatch = sum(flosum(:))/numel(find(matchFlag.*idxflag == 1));
% 			flosum = v.*matchFlag.*idxflag;
% 			vMatch = sum(flosum(:))/numel(find(matchFlag.*idxflag == 1));
% 		else
% 			uMatch = 0;
% 			vMatch = 0;
% 		end
% 	end
% end