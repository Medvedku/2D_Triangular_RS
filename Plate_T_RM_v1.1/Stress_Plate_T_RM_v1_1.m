function Strnod = Stress_Plate_T_RM_v1_1(D_matb,D_mats,xg,yg,u)

%% Stress Evaluates the stresses at the Gauss points and smooth the values
%         to the nodes
%
%  Parameters:
%
%    Input, D_matb : Constitutive matrix for bending moment
%           D_mats : Constitutive matrix for shear force
%           xg     : Local X coordinates of the Gauss points
%           yg     : Local Y coordinates of the Gauss points
%           u      : Nodal displacements
%   
%    Output, Strnod the nodal stress matrix (nnode,nstrs)

  global coordinates;
  global elements;
 
% Find basic dimensions
  nelem  = size(elements,1);          % Number of elements
  nnode  = size(elements,2);          % Number of nodes por element
  npnod  = size(coordinates,1);       % Number of nodes
  Strnod = zeros(npnod,6);            % Create array for stresses
  dofpn  = 3;                         % Number of DOF per node
  dofpe  = dofpn*nnode;               % Number of DOF per element
  eqnum  = zeros(dofpe);              % Equation number list

% Element cycle
  for ielem = 1 : nelem
 
    lnods(1:nnode) = elements(ielem,1:nnode);
    
% Find the equation number list for the i-th element
    for i = 1 : nnode
      ii = (i-1)*dofpn;
      for j =1:dofpn
        eqnum(ii+j) = (lnods(i)-1)*dofpn + j;   % Build the eq. number list
      end
    end
    
% Recover the nodal displacements for the i-th element
    u_elem(1:dofpe) = u(eqnum(1:dofpe));
  
    x(1:nnode) = coordinates(lnods(1:nnode),1);       % Elem. X coordinates
    y(1:nnode) = coordinates(lnods(1:nnode),2);       % Elem. Y coordinates

    for igaus = 1 : 3

      [bmat_b,bmat_s,area] = B_mat_Plate_T_RM_v1_1(x,y,xg(igaus),yg(igaus)); 

      Str1 = D_matb*bmat_b*transpose(u_elem);
      Str2 = D_mats*bmat_s*transpose(u_elem);
      
      Strnod(lnods(igaus),1) = Strnod(lnods(igaus),1) + Str1(1);
      Strnod(lnods(igaus),2) = Strnod(lnods(igaus),2) + Str1(2);
      Strnod(lnods(igaus),3) = Strnod(lnods(igaus),3) + Str1(3);
      Strnod(lnods(igaus),4) = Strnod(lnods(igaus),4) + Str2(1);
      Strnod(lnods(igaus),5) = Strnod(lnods(igaus),5) + Str2(2);
      Strnod(lnods(igaus),6) = Strnod(lnods(igaus),6) + 1;
 
    end

  end
  
  for i = 1 : npnod
    Strnod(i,1:5) = Strnod(i,1:5)/Strnod(i,6);
  end   
    
