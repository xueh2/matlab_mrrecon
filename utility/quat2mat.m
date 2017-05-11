function R = quat2mat(q);
% function R = quat2mat(q);
%
% function to convert quaternion vector (q) to rotation matrix (R)

%     ***************************************
%     *  Peter Kellman  (kellman@nih.gov)   *
%     *  Laboratory for Cardiac Energetics  *
%     *  NIH NHLBI                          *
%     ***************************************

 
MDH_QW=q.a;
MDH_QX=q.b;
MDH_QY=q.c;
MDH_QZ=q.d;


ds = 2.0 / (MDH_QW * MDH_QW +...
            MDH_QX * MDH_QX +...
            MDH_QY * MDH_QY +...
            MDH_QZ * MDH_QZ);

dxs = MDH_QX *  ds; dys = MDH_QY *  ds; dzs = MDH_QZ *  ds;
dwx = MDH_QW * dxs; dwy = MDH_QW * dys; dwz = MDH_QW * dzs;
dxx = MDH_QX * dxs; dxy = MDH_QX * dys; dxz = MDH_QX * dzs;
dyy = MDH_QY * dys; dyz = MDH_QY * dzs; dzz = MDH_QZ * dzs;

R(1,1) = 1.0 - (dyy + dzz);
R(1,2) =        dxy + dwz ;
R(1,3) =        dxz - dwy ;

R(2,1) =        dxy - dwz ;
R(2,2) = 1.0 - (dxx + dzz);
R(2,3) =        dyz + dwx ;

R(3,1) =        dxz + dwy ;
R(3,2) =        dyz - dwx ;
R(3,3) = 1.0 - (dxx + dyy);        
        
% from MdhProxy.h
%       static void quat2mat (const double adQuaternion[4], double adRotMatrix[3][3])
%       {
%         double ds, dxs, dys, dzs, dwx, dwy, dwz, dxx, dxy, dxz, dyy, dyz, dzz;
% 
%         ds = 2.0 / (adQuaternion[MDH_QW] * adQuaternion[MDH_QW] +
%                     adQuaternion[MDH_QX] * adQuaternion[MDH_QX] +
%                     adQuaternion[MDH_QY] * adQuaternion[MDH_QY] +
%                     adQuaternion[MDH_QZ] * adQuaternion[MDH_QZ]);
% 
%         dxs = adQuaternion[MDH_QX] *  ds; dys = adQuaternion[MDH_QY] *  ds; dzs = adQuaternion[MDH_QZ] *  ds;
%         dwx = adQuaternion[MDH_QW] * dxs; dwy = adQuaternion[MDH_QW] * dys; dwz = adQuaternion[MDH_QW] * dzs;
%         dxx = adQuaternion[MDH_QX] * dxs; dxy = adQuaternion[MDH_QX] * dys; dxz = adQuaternion[MDH_QX] * dzs;
%         dyy = adQuaternion[MDH_QY] * dys; dyz = adQuaternion[MDH_QY] * dzs; dzz = adQuaternion[MDH_QZ] * dzs;
% 
%         adRotMatrix[0][0] = 1.0 - (dyy + dzz);
%         adRotMatrix[0][1] =        dxy + dwz ;
%         adRotMatrix[0][2] =        dxz - dwy ;
%         adRotMatrix[1][0] =        dxy - dwz ;
%         adRotMatrix[1][1] = 1.0 - (dxx + dzz);
%         adRotMatrix[1][2] =        dyz + dwx ;
%         adRotMatrix[2][0] =        dxz + dwy ;
%         adRotMatrix[2][1] =        dyz - dwx ;
%         adRotMatrix[2][2] = 1.0 - (dxx + dyy);
%       }
