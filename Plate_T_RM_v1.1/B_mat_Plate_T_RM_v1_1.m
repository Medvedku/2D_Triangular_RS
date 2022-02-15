function [bmat_b,bmat_s,area] = B_mat_Plate_T_RM_v1_1(x,y,xgs,ygs)

%% B_Mat Computes the strain-displacement matrix for bending moments
%        and shear forces
%
%  Parameters:
%
%    Input, x   : X Coordinates of the element
%           y   : Y Coordinates of the element
%           xgs : Local X coordinate of the Gauss point
%           ygs : Local Y coordinate of the Gauss point
%   
%    Output, bmat_b the strain-displacement matrix for bending moments
%            bmat_s the strain-displacement matrix for shear forces
%              area the element area

%==================

  N(1) = 1.0 - xgs - ygs ;
  N(2) = xgs ;
  N(3) = ygs ;
  
  dxNloc(1) = -1.0;
  dxNloc(2) =  1.0;
  dxNloc(3) =  0.0;

  dyNloc(1) = -1.0;
  dyNloc(2) =  0.0;
  dyNloc(3) =  1.0;

  xjacm(1,1) = x(1)*dxNloc(1) + x(2)*dxNloc(2) + x(3)*dxNloc(3);
  xjacm(1,2) = y(1)*dxNloc(1) + y(2)*dxNloc(2) + y(3)*dxNloc(3);
  xjacm(2,1) = x(1)*dyNloc(1) + x(2)*dyNloc(2) + x(3)*dyNloc(3);
  xjacm(2,2) = y(1)*dyNloc(1) + y(2)*dyNloc(2) + y(3)*dyNloc(3);
  
  xjaci = inv(xjacm);

  area2 = abs(xjacm(1,1)*xjacm(2,2) - xjacm(2,1)*xjacm(1,2));
  area  = area2/2;
  
  dxN(1) = xjaci(1,1)*dxNloc(1)+xjaci(1,2)*dyNloc(1);
  dxN(2) = xjaci(1,1)*dxNloc(2)+xjaci(1,2)*dyNloc(2);
  dxN(3) = xjaci(1,1)*dxNloc(3)+xjaci(1,2)*dyNloc(3);
 
  dyN(1) = xjaci(2,1)*dxNloc(1)+xjaci(2,2)*dyNloc(1);
  dyN(2) = xjaci(2,1)*dxNloc(2)+xjaci(2,2)*dyNloc(2);
  dyN(3) = xjaci(2,1)*dxNloc(3)+xjaci(2,2)*dyNloc(3);

%==================

      bmat_b1  = [ 0,-dxN(1),     0  ;
                   0,      0,-dyN(1) ;
                   0,-dyN(1),-dxN(1)];
               
      bmat_b2  = [ 0,-dxN(2),     0  ;
                   0,      0,-dyN(2) ;
                   0,-dyN(2),-dxN(2)];
               
      bmat_b3  = [ 0,-dxN(3),     0  ;
                   0,      0,-dyN(3) ;
                   0,-dyN(3),-dxN(3)];
 
      bmat_b = [bmat_b1,bmat_b2,bmat_b3];
      
%==================

      bmat_s1  = [ dxN(1), -N(1),     0 ;
                   dyN(1),     0, -N(1)];
               
      bmat_s2  = [ dxN(2), -N(2),     0 ;
                   dyN(2),     0, -N(2)];

      bmat_s3  = [ dxN(3), -N(3),     0 ;
                   dyN(3),     0, -N(3)];

      bmat_s = [bmat_s1,bmat_s2,bmat_s3];

