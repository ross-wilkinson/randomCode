function [S] = computeComEnergetics(varargin)
    %COMPUTECOMENERGETICS Take body kinematics output file from OpenSim 3.3
    %along with user specified additional data to calculate rider center of
    %mass energetics over each complete crank cycle.
    %   
    %   INPUTS:
    %           - 'subjectMass' = 75. (kg)
    %           - 'bodyKinematicsFile' = Eg. 'C:\...\data.sto'
    %           - 'angleData' = Eg. 1x1000 double
    %           - 'forceData' = Eg. 1x1000 double
    %           - 'targetPower' = Eg. 500. (in Watts)
    %           - 'targetCadence' = Eg. 120. (in rpm)
    %           - 'conditionName' = Eg. 'seated50120'
    %           - 'subjectName' = 'subject01'
    %           - 'buffers' = [bufferPowerLow bufferPowerHigh
    %           bufferCadenceLow bufferCadenceHigh]
    %           - 'crankLength' = Eg. 0.1725 (m)
    %
    %   OUTPUT:
    %           - S = Structure containing the following fields:
    %                   - cadence
    %                   - cadenceMean
    %                   - power
    %                   - powerMean
    %                   - force
    %                   - forceMean
    %                   - comPosX
    %                   - comPosY
    %                   - comPosZ
    %                   - comVelX
    %                   - comVelX
    %                   - comVelY
    %                   - comVelZ
    %                   - comAccX
    %                   - comAccY
    %                   - comAccZ
    %                   - comPotentialEnergy
    %                   - comKineticEnergy
    %                   - comTotalEnergy
    %   
    %           NOTE: these fields will be added to the structure
    %           S.(conditionName).(subjectName) if the conditionName and
    %           subjectName are specified as variables.
    
    subjectMass = [];
    bodyKinematicsFile = [];
    angleData = [];
    forceData = [];
    targetPower = [];
    targetCadence = [];
    conditionName = [];
    subjectName = [];
    buffers = [];
    crankLength = 0.1725;

    if ~isempty(varargin)
        if rem(length(varargin),2)
            error('Incorrect input arguments - must specify property and input')
        else
            for i = 1:2:length(varargin)
                n = varargin{i+1};
                eval([varargin{i} '= n;']); 
            end 
            if isempty(subjectMass)
                error('No subject mass input. Please input a subject mass');
            end
            if isempty(bodyKinematicsFile)
                error('No kinematics data input. Please input kinematic data to analyze');
            end
            if isempty(angleData)
                error('No angle data input. Please input angle data to analyze');
            end
        end
    end

    BK = importdata(bodyKinematicsFile,'\t');
    iTime = strcmp(BK.colheaders,'time');
    iX = strcmp(BK.colheaders,'center_of_mass_X');
    iY = strcmp(BK.colheaders,'center_of_mass_Y');
    iZ = strcmp(BK.colheaders,'center_of_mass_Z');

    T = table();
    T.time = BK.data(:,iTime);
    T.comPosX = BK.data(:,iX);
    T.comPosY = BK.data(:,iY);
    T.comPosZ = BK.data(:,iZ);
    if ~isempty(angleData)
        if size(angleData,2) > 1
            angleData = angleData';
        end
        T.angle = angleData;
    end
    if ~isempty(forceData)
        if size(forceData,2) > 1
            forceData = forceData';
        end
        T.force = forceData;
    end

    if size(angleData,2) > 1
        angleData = angleData';
    end
    T.angle = angleData;
    % cut up kinematics data into cycles
    % find peak locations in angle signal
    bufferHeight = 0.8;
    [~,locs] = findpeaks(T.angle,'minpeakheight',max(T.angle) * bufferHeight);

    % calculate cadence for each cycle
    nCycles = numel(locs) - 1;
    cadence = NaN(nCycles,1);
    varList = T.Properties.VariableNames;

    for iLocs = 1:nCycles
        cadence(iLocs) = 60 / T.time(locs(iLocs + 1)) - T.time(locs(iLocs));
        % interpolate data between locs into 101 points
        for iVars = 2:length(T.Properties.VariableNames)
            varName = varList(iVars);
            x1 = 0:1 / (locs(iLocs + 1) - locs(iLocs)):1;
            v1 = varName(locs(iLocs):locs(iLocs + 1))';
            xq1 = 0:1 / 100:1;
            vq1 = interp1(x1,v1,xq1,'spline');
            span = 10;
            method = 'sgolay';
            varData = smooth(vq1,span,method);
            R.(varName)(iLocs) = varData; 
            if strcmp(varName,'comPosX') || strcmp(varName,'comPosY') ||...
                    strcmp(varName,'comPosZ')
                h = 60 / cadence((iLocs)) / 100;
                varDataVelocity = smooth(diff(varData / h,span,method));
                varDataAcceleration = smooth(diff(varDataVelocity) / h,span,method);
                varNameVelocity = strrep(varName,'Pos','Vel');
                varNameAcceleration = strrep(varName,'Pos','Acc');
                R.(varNameVelocity)(iLocs) = varDataVelocity;
                R.(varNameAcceleration)(iLocs) = varDataAcceleration;
            end
            if strcmp(varName,'force')
                crankVelocity = cadence(iLocs) * 2 * pi / 60;
                crankTorque = varData * crankLength;
                R.power(iLocs) = crankTorque * crankVelocity;
                R.powerMean(iLocs) = mean(R.power(iLocs));
            end
        end
        R.comPotentialEnergy(iLocs) = subjectMass * 9.81 * R.comPosY(iLocs);
        resultantComVel = sqrt(R.comVelX(iLocs).^2 + R.comVelX(iLocs).^2 + ...
            R.comVelX(iLocs).^2);
        R.comKineticEnergy(iLocs) = 0.5 * subjectMass * resultantComVel.^2; 
    end
    R.comTotalEnergy = R.comKineticEnergy + R.comPotentialEnergy;

    % if a target power is input then validate power during each cycle
    if ~isempty(targetPower) && ~isempty(targetCadence)
        if isempty(forceData)
            error('No force data input. Please input force data to analyze.')
        end
        if isempty(buffers)
            error('No buffers input. Please input buffer for valid data.')
        end

        k1 = find(R.powerMean > targetPower * buffer(1) & R.powerMean <...
            targetPower * buffer(2));
        k2 = find(cadence > targetCadence * buffer(3) & cadence <...
            targetCadence * buffer(4));
        iValid = intersect(k1,k2);
        fieldList = fieldNames(R);
        for iFields = 1:numel(fieldnames)
            fieldName = fieldList(iFields);
            V.(fieldName) = R.(fieldName)(iValid,:);
        end
        if ~isempty(conditionName)
            if ~isempty(subjectName)
                S.(conditionName).(subjectName) = V;
            else
                S.(conditionName) = V;
            end
        else
            S = V;
        end
    else
        if ~isempty(conditionName)
            if ~isempty(subjectName)
                S.(conditionName).(subjectName) = R;
            else
                S.(conditionName) = R;
            end
        else
            S = R;
        end    
    end

end