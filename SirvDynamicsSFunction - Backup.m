function SirvDynamicsSFunction(block)
setup(block);
%This Level-2 S-Function simulates the behavior of a quadcopter based on
%performance parameters and RPM inputs. 
end

% Function: setup ===================================================
% Abstract:
%   Set up the S-function block's basic characteristics such as:
%   - Input ports
%   - Output ports
%   - Dialog parameters
%   - Options
% 
%   Required         : Yes
%   C-Mex counterpart: this MATLAB Level-2 S-function is doing similar work
%   to what the mdlInitializeSizes function would do in a C MEX S-function
%
function setup(block)

  % Register the number of ports.
  %------
  block.NumInputPorts  = 5;
  %------
  block.NumOutputPorts = 12;
  
  % Set up the port properties to be inherited or dynamic.
  
  for i = 1:4 % These are the motor inputs
  block.InputPort(i).Dimensions        = 1; %this is expecting to receive a scalar value
  block.InputPort(i).DirectFeedthrough = false; %none of these ports are set to direct
  % feedthrough, meaning that the inputs do not directly affect the outputs
  block.InputPort(i).SamplingMode      = 'Sample'; %This means that the input
  % port is expecting discrete-time signals, where the signal values are sampled at specific time intervals.
  end
  %------
  % This is the ext. disturbances input
  block.InputPort(5).Dimensions        = 6; % torques x,y,z; forces x,y,z.
  block.InputPort(5).DirectFeedthrough = false;
  block.InputPort(5).SamplingMode      = 'Sample'; %set to 'Inherited' for continuous
  %------
  for i = 1:12
  block.OutputPort(i).Dimensions       = 1;
  block.OutputPort(i).SamplingMode     = 'Sample'; %set to 'Inherited' for continuous
  end

%   %In case we need to override the input port properties.
%   block.InputPort(1).DatatypeID  = 0;  % 0 sets to double in simulink
%   block.InputPort(1).Complexity  = 'Real'; %This means that the first input port is expecting a real number (as opposed to a complex number).
%   
%   %In case we need to override the output port properties.
%   block.OutputPort(1).DatatypeID  = 0; % double
%   block.OutputPort(1).Complexity  = 'Real'; %This means that the first output port will output a real number.

  % Register the parameters.
  block.NumDialogPrms     = 2; %specifying that your block expects exactly two parameters
  
  % Set up the continuous states (state variables).
  block.NumContStates = 12; %Continuous states are used in the simulation of continuous dynamic systems. They represent the state variables of the system that change continuously over time, as opposed to discrete states which change at specific time steps.

  % Register the sample times.
  %  [0 offset]            : Continuous sample time
  %  [positive_num offset] : Discrete sample time (e.g. [0.1 offset])
  %  if you set the sample time to [0.1 0.01], the block%s output would be updated
  %  every 0.1 seconds, starting from 0.01 seconds into the simulation.
  %
  %  [-1, 0]               : Inherited sample time (the block inherits its sample time from the block that drives it).
  %  [-2, 0]               : Variable sample time (can change during
  %  simulation)
  block.SampleTimes = [0 0]; %A sample time of [0 0] means that the block
  % has a continuous sample time, i.e., its output is computed at every time step of the simulation.

  % -----------------------------------------------------------------
  % Options
  % -----------------------------------------------------------------
  % Specify if Accelerator should use TLC or call back to the 
  % MATLAB file. Simulink’s Accelerator mode is a simulation mode that speeds
  % up the execution of your model by automatically generating, compiling, and executing C code.
  % If AccelRunOnTLC is true, Simulink generates Target Language Compiler (TLC) code for the block
  % and uses this code to run the block in Accelerator mode.
  block.SetAccelRunOnTLC(false);
  
  % Specify the block simStateCompliance. The allowed values are:
  %    'UnknownSimState', < The default setting; warn and assume DefaultSimState
  %    'DefaultSimState', < Same SimState as a built-in block
  %    'HasNoSimState',   < No SimState
  %    'CustomSimState',  < Has GetSimState and SetSimState methods
  %    'DisallowSimState' < Errors out when saving or restoring the SimState
  block.SimStateCompliance = 'DefaultSimState'; %
  
  % -----------------------------------------------------------------
  % The MATLAB S-function uses an internal registry for all
  % block methods. You should register all relevant methods
  % (optional and required) as illustrated below. You may choose
  % any suitable name for the methods and implement these methods
  % as local functions within the same file.
  % -----------------------------------------------------------------
   
  % -----------------------------------------------------------------
  % Register the methods called during update diagram/compilation.
  % -----------------------------------------------------------------
  
%   % CHECK PARAMETERS
%   % CheckParameters:
%   %   Functionality    : Called in order to allow validation of the
%   %                      block dialog parameters. You are 
%   %                      responsible for calling this method
%   %                      explicitly at the start of the setup method.
%   %   C-Mex counterpart: mdlCheckParameters
%   %
   block.RegBlockMethod('CheckParameters', @CheckPrms);

%   %SAMPLING MODE
%   % SetInputPortSamplingMode:
%   %   Functionality    : Check and set input and output port 
%   %                      attributes and specify whether the port is operating 
%   %                      in sample-based or frame-based mode
%   %   C-Mex counterpart: mdlSetInputPortFrameData.
%   %   (The DSP System Toolbox is required to set a port as frame-based)
%   %
%   block.RegBlockMethod('SetInputPortSamplingMode', @SetInpPortFrameData);
%   % this will probably be deleted since we have already set the sampling
%     mode for all input/output ports
  
%   %INPUT PORT DIMENSIONS
%   % SetInputPortDimensions:
%   %   Functionality    : Check and set the input and optionally the output
%   %                      port dimensions.
%   %   C-Mex counterpart: mdlSetInputPortDimensionInfo
%   %
%   block.RegBlockMethod('SetInputPortDimensions', @SetInpPortDims); 
%   % this will probably be deleted since we have already set the dimension
%   % for all input ports

%   %OUTPUT PORT DIMENSIONS
%   % SetOutputPortDimensions:
%   %   Functionality    : Check and set the output and optionally the input
%   %                      port dimensions.
%   %   C-Mex counterpart: mdlSetOutputPortDimensionInfo
%   %
%   block.RegBlockMethod('SetOutputPortDimensions', @SetOutPortDims);
%   % this will probably be deleted since we have already set the dimension
%   % for all OUTPUT ports
    

%   %INPUT PORT DATATYPE
%   % SetInputPortDatatype:
%   %   Functionality    : Check and set the input and optionally the output
%   %                      port datatypes.
%   %   C-Mex counterpart: mdlSetInputPortDataType
%   %
%   block.RegBlockMethod('SetInputPortDataType', @SetInpPortDataType);
%   % this will probably be deleted since we have already set the DATATYPE
%   % for all INPUT ports

%   %OUTPUT PORT DATATYPE
%   % SetOutputPortDatatype:
%   %   Functionality    : Check and set the output and optionally the input
%   %                      port datatypes.
%   %   C-Mex counterpart: mdlSetOutputPortDataType
%   %
%   block.RegBlockMethod('SetOutputPortDataType', @SetOutPortDataType);
%   % this will probably be deleted since we have already set the DATATYPE
%   % for all OUTPUT ports
  
%   %INPUT PORT COMPLEXITY
%   % SetInputPortComplexSignal:
%   %   Functionality    : Check and set the input and optionally the output
%   %                      port complexity attributes.
%   %   C-Mex counterpart: mdlSetInputPortComplexSignal
%   %
%   block.RegBlockMethod('SetInputPortComplexSignal', @SetInpPortComplexSig);

%   %OUTPUT PORT COMPLEXITY
%   % SetOutputPortComplexSignal:
%   %   Functionality    : Check and set the output and optionally the input
%   %                      port complexity attributes.
%   %   C-Mex counterpart: mdlSetOutputPortComplexSignal
%   %
%   block.RegBlockMethod('SetOutputPortComplexSignal', @SetOutPortComplexSig);
  
%   %OK
%   % PostPropagationSetup:
%   %   Functionality    : Set up the work areas and the state variables. You can
%   %                      also register run-time methods here.
%   %   C-Mex counterpart: mdlSetWorkWidths
%   %
%   block.RegBlockMethod('PostPropagationSetup', @DoPostPropSetup);
% 
%   % -----------------------------------------------------------------
%   % Register methods called at run-time
%   % -----------------------------------------------------------------
%   
%   % OK
%   % ProcessParameters:
%   %   Functionality    : Call to allow an update of run-time parameters.
%   %   C-Mex counterpart: mdlProcessParameters
%   %  
%   block.RegBlockMethod('ProcessParameters', @ProcessPrms);

  % OK
  % InitializeConditions:
  %   Functionality    : Call to initialize the state and the work
  %                      area values.
  %   C-Mex counterpart: mdlInitializeConditions
  % 
  block.RegBlockMethod('InitializeConditions', @InitializeConditions);
  
%   % OK
%   % Start:
%   %   Functionality    : Call to initialize the state and the work
%   %                      area values.
%   %   C-Mex counterpart: mdlStart
%   %
%   block.RegBlockMethod('Start', @Start);

  % OK
  % Outputs:
  %   Functionality    : Call to generate the block outputs during a
  %                      simulation step.
  %   C-Mex counterpart: mdlOutputs
  %
  block.RegBlockMethod('Outputs', @Outputs);

%   % OK
%   % Update:
%   %   Functionality    : Call to update the discrete states
%   %                      during a simulation step.
%   %   C-Mex counterpart: mdlUpdate
%   %
%   block.RegBlockMethod('Update', @Update);

  % OK
  % Derivatives:
  %   Functionality    : Call to update the derivatives of the
  %                      continuous states during a simulation step.
  %   C-Mex counterpart: mdlDerivatives
  %
  block.RegBlockMethod('Derivatives', @Derivatives);
  
%   % OK
%   % Projection:
%   %   Functionality    : Call to update the projections during a
%   %                      simulation step.
%   %   C-Mex counterpart: mdlProjections
%   %
%   block.RegBlockMethod('Projection', @Projection);
  
%   % OK
%   % SimStatusChange:
%   %   Functionality    : Call when simulation enters pause mode
%   %                      or leaves pause mode.
%   %   C-Mex counterpart: mdlSimStatusChange
%   %
%   block.RegBlockMethod('SimStatusChange', @SimStatusChange);
%   
%   % OK
%   % Terminate:
%   %   Functionality    : Call at the end of a simulation for cleanup.
%   %   C-Mex counterpart: mdlTerminate
%   %
%   block.RegBlockMethod('Terminate', @Terminate);

%   % OK
%   % GetSimState:
%   %   Functionality    : Return the SimState of the block.
%   %   C-Mex counterpart: mdlGetSimState
%   %
%   block.RegBlockMethod('GetSimState', @GetSimState);
%   
%   %
%   % SetSimState:
%   %   Functionality    : Set the SimState of the block using a given value.
%   %   C-Mex counterpart: mdlSetSimState
%   %
%   block.RegBlockMethod('SetSimState', @SetSimState);

  % -----------------------------------------------------------------
  % Register the methods called during code generation.
  % -----------------------------------------------------------------
  
%   %
%   % WriteRTW:
%   %   Functionality    : Write specific information to model.rtw file.
%   %   C-Mex counterpart: mdlRTW
%   %
%   block.RegBlockMethod('WriteRTW', @WriteRTW);
% %endfunction

% -------------------------------------------------------------------
% The local functions below are provided to illustrate how you may implement
% the various block methods listed above.
% -------------------------------------------------------------------
end
 function CheckPrms(block)
     quad   = block.DialogPrm(1).Data;
     IC     = block.DialogPrm(2).Data;

      % Check that 'quad' is a structure
    if ~isstruct(quad)
        error('The first dialog parameter must be a structure.');
    end

    % Check that 'IC' is a structure
    if ~isstruct(IC)
        error('The second dialog parameter must be a structure.');
    end
end
     
     
     %      if ~exist(model)
%        me = MSLException(block.BlockHandle, message('Simulink:blocks:invalidParameter'));
%        throw(me);
%      end
     
%     if ~strcmp(class(a), 'double')
%       me = MSLException(block.BlockHandle, message('Simulink:blocks:invalidParameter'));
%       throw(me);
%     end
% %endfunction
% 
% function ProcessPrms(block)
% 
%   block.AutoUpdateRuntimePrms;
%  
% %endfunction
% 
% function SetInpPortFrameData(block, idx, fd)
%   
%   block.InputPort(idx).SamplingMode = fd;
%   block.OutputPort(1).SamplingMode  = fd;
%   
% %endfunction
% 
% function SetInpPortDims(block, idx, di)
%   
%   block.InputPort(idx).Dimensions = di;
%   block.OutputPort(1).Dimensions  = di;
% 
% %endfunction
% 
% function SetOutPortDims(block, idx, di) 
%   
%   block.OutputPort(idx).Dimensions = di;
%   block.InputPort(1).Dimensions    = di;
% 
% %endfunction
% 
% function SetInpPortDataType(block, idx, dt)
%   
%   block.InputPort(idx).DataTypeID = dt;
%   block.OutputPort(1).DataTypeID  = dt;
% 
% %endfunction
%   
% function SetOutPortDataType(block, idx, dt)
% 
%   block.OutputPort(idx).DataTypeID  = dt;
%   block.InputPort(1).DataTypeID     = dt;
% 
% %endfunction  
% 
% function SetInpPortComplexSig(block, idx, c)
%   
%   block.InputPort(idx).Complexity = c;
%   block.OutputPort(1).Complexity  = c;
% 
% %endfunction 
%   
% function SetOutPortComplexSig(block, idx, c)
% 
%   block.OutputPort(idx).Complexity = c;
%   block.InputPort(1).Complexity    = c;
% 
% %endfunction 
%     
% function DoPostPropSetup(block)
%
%   block.NumDworks = 2;
%   
%   block.Dwork(1).Name            = 'x1';
%   block.Dwork(1).Dimensions      = 1;
%   block.Dwork(1).DatatypeID      = 0;      % double
%   block.Dwork(1).Complexity      = 'Real'; % real
%   block.Dwork(1).UsedAsDiscState = true;
%   
%   block.Dwork(2).Name            = 'numPause';
%   block.Dwork(2).Dimensions      = 1;
%   block.Dwork(2).DatatypeID      = 7;      % uint32
%   block.Dwork(2).Complexity      = 'Real'; % real
%   block.Dwork(2).UsedAsDiscState = true;
%   
%   % Register all tunable parameters as runtime parameters.
%   block.AutoRegRuntimePrms;
% 
% %endfunction

function InitializeConditions(block)
% Initialize 12 States

IC = block.DialogPrm(2).Data; %retrieves the initial conditions (IC) from the second dialog parameter of the block.

% IC.P, IC.Q, IC.R are in deg/s ... convert to rad/s
P = IC.P*pi/180;
Q = IC.Q*pi/180;
R = IC.R*pi/180;

% IC.Phi, IC.Theta, IC.Psi are in degrees. Convert to rads
Phi = IC.Phi*pi/180;
Theta = IC.Theta*pi/180;
Psi = IC.Psi*pi/180;

U = IC.U;
V = IC.V;
W = IC.W;

X = IC.X;
Y = IC.Y;
Z = IC.Z;

init = [P,Q,R,Phi,Theta,Psi,U,V,W,X,Y,Z]; % creates a vector init that contains all the initial conditions.
for i=1:12
block.OutputPort(i).Data = init(i);
block.ContStates.Data(i) = init(i);
end
end
% function Start(block)
% 
%   block.Dwork(1).Data = 0;
%   block.Dwork(2).Data = uint32(1); 
%    
% %endfunction
% 
% function WriteRTW(block)
%   
%    block.WriteRTWParam('matrix', 'M',    [1 2; 3 4]);
%    block.WriteRTWParam('string', 'Mode', 'Auto');
   
%endfunction

function Outputs(block)
for i = 1:12
  block.OutputPort(i).Data = block.ContStates.Data(i);
end
end
%endfunction

% function Update(block)
%   
%   block.Dwork(1).Data = block.InputPort(1).Data;
  
%endfunction

function Derivatives(block)
% Name all the states and motor inputs

% Load model data selected in parameter block
%which('model')

quad = block.DialogPrm(1).Data; %This is only for convenience and readability. Calling variables from the structure inside the code will be like this "quad.ct" instead of this "quadModel.ct"

% P Q R in units of rad/sec
P = block.ContStates.Data(1);
Q = block.ContStates.Data(2);
R = block.ContStates.Data(3);
% Phi The Psi in radians
Phi = block.ContStates.Data(4);
Theta = block.ContStates.Data(5);
Psi = block.ContStates.Data(6);
% U V W in units of m/s
U = block.ContStates.Data(7);
V = block.ContStates.Data(8);
W = block.ContStates.Data(9);
% X Y Z in units of m
X = block.ContStates.Data(10);
Y = block.ContStates.Data(11);
Z = block.ContStates.Data(12);
% w values in rev/min! NOT radians/s!!!!
w1 = block.InputPort(1).Data;
w2 = block.InputPort(2).Data;
w3 = block.InputPort(3).Data;
w4 = block.InputPort(4).Data;
w  = [w1; w2; w3; w4];

%------
Dist_tau = block.InputPort(5).Data(1:3); %external torques. ports 1-3
Dist_F   = block.InputPort(5).Data(4:6); %external forces. ports 1-3
%------

% CALCULATE MOMENT AND THRUST FORCES
% Total Moment due to motor speeds
% Moment should be in units of N*m
% The experimental determination of Ct and Cq should be adjusted to
% model using kg instead of ounces or lb
% Mb = (quad.dctcq*(w.^2)) + (Dist_tau);  %(dctcq*(w.^2)); % Mb = [tau1 tau2 tau3]'
 tau_motorGyro = [Q*quad.Jm*2*pi/60*(-w1-w2+w3+w4); P*quad.Jm*2*pi/60*(w1+w2-w3-w4); 0]; % Note: 2*pi/60 required to convert from RPM to radians/s
 Mb = (quad.dctcq*(w.^2))+ tau_motorGyro + (Dist_tau);  % Mb = [tau1 tau2 tau3]' /Mb=total moment acting on quadcopter

% Thrust due to motor speed
% Force should be in units of Newtons for simplicity in calculating
% the acceleration in the angular velocity state equation
Fb = [0; 0; sum(quad.ct*(w.^2))];   %[0, 0, sum(ct*w.^2)]'

% Obtain dP dQ dR
omb_bi = [P; Q; R]; %vector with ang. velocities in body frame
OMb_bi = [ 0,-R, Q; %skew symmetric matrix. t’s used to perform cross product operations in matrix form in the equations of motion
           R, 0,-P;
          -Q, P, 0];

b_omdotb_bi = quad.Jbinv*(Mb-OMb_bi*quad.Jb*omb_bi);%used to convert the ang. velocities in the body frame to the derivatives of the Euler angles
H_Phi = [1,tan(Theta)*sin(Phi), tan(Theta)*cos(Phi); %the elements used to calc the derivative of the roll angle Phi
         0,         cos(Phi),         -sin(Phi);     %the elements for the derivative of the pitch angle theta
         0,sin(Phi)/cos(Theta),cos(Phi)/cos(Theta)]; %the elements for the derivative of the yaw angle psi  
Phidot = H_Phi*omb_bi; %uses H_Phi to convert the ang. velocities omb_bi to the derivatives of the Euler angles Phidot

% Compute Rotation Matrix. Z-Y-X rotation
Rib = [cos(Psi)*cos(Theta) cos(Psi)*sin(Theta)*sin(Phi)-sin(Psi)*cos(Phi) cos(Psi)*sin(Theta)*cos(Phi)+sin(Psi)*sin(Phi);
       sin(Psi)*cos(Theta) sin(Psi)*sin(Theta)*sin(Phi)+cos(Psi)*cos(Phi) sin(Psi)*sin(Theta)*cos(Phi)-cos(Psi)*sin(Phi);
       -sin(Theta)         cos(Theta)*sin(Phi)                            cos(Theta)*cos(Phi)];
Rbi = Rib'; %calculates the transpose of the rotation matrix Rib. So Rbi:transforms vectors from inertial frame to body frame.
ge = [0; 0; -quadModel.g]; % (ge=Gravity in Earth/inertial frame) gravitational acceleration in the inertial frame (directed downwards and hence the negative in Z)
gb = Rbi*ge; % (gb=Gravity in Body frame) transformation of the vector from world frame to body frame
Dist_Fb = Rbi*Dist_F; %Calculates the ext forces in the body frame by transforming the ext. force vector from inertial frame to body frame.

% Compute Velocity and Position derivatives of body frame
vb = [U;V;W]; %vb = velocity body (frame)
b_dv = (1/quadModel.mass)*Fb+gb+Dist_Fb-OMb_bi*vb; % Acceleration in body frame (FOR VELOCITY). Calculated using the forces acting on the body (incl. thrust force Fb, gravity body gb and ext. disturbances Dist_Fb
i_dp = Rib*vb; % Units OK SI: Velocity of body frame w.r.t inertia frame (FOR POSITION). Calc velocity of body with respect to the inertial frame usig the rotation matrix

dP = b_omdotb_bi(1);
dQ = b_omdotb_bi(2);
dR = b_omdotb_bi(3);
dPhi = Phidot(1);
dTheta = Phidot(2);
dPsi = Phidot(3);
dU = b_dv(1);
dV = b_dv(2);
dW = b_dv(3);
dX = i_dp(1);
dY = i_dp(2);
dZ = i_dp(3);
% Rough rule to impose a "ground" boundary...could easily be improved...
if ((Z<=0) && (dZ<=0)) % checks quadcopter altitude and vertical velocity dZ. If ground is met it sets dZ to zero
    dZ = 0;
    block.ContStates.Data(12) = 0; %resets the copters altitude in the state vector to zero to prevent sinking into the ground
end
f = [dP dQ dR dPhi dTheta dPsi dU dV dW dX dY dZ].'; %creates column vector f
  %This is the state derivative vector
block.Derivatives.Data = f; %assigns vector f to the derivatives.Data property of the block. Tells simulink to use these derivatives toupdate the state of the quadcopter

end
%endfunction