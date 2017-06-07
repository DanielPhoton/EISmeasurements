clear all 
close all

%Parameter definitions

L = 1;
Rtn = 20;
Rdr = 500; 
Rtp = 10; 
Rr = 50;
Cmun = 2*10^-6;
Cmup = 2*10^-6;
Cg = 10^-9;
Cdr = 1*10^-3;
alpha = 0;
gamma_n = 0;
gamma_p = 0;

J0 = 1;

zn = Rtn; 
zp = Rtp;
yr = 1/Rr;

l = 1; 
Zre = 0; 
Zim = 0;

% Calculation of the impedance at each frequency.

for l = 1:14300
    

    
% Definition of the frequency

f = exp(l/700)-1;


%Definition of the frequency dependent elements.
yd = 1i * f * Cg+1i * f * Cdr / (1+(1i * f * Cdr * Rdr)^(1-alpha));
zd = 1/yd;
yn = 1i * f * Cmun;
yp = 1i * f * Cmup;

%Definition of the elements of the vector Up 

a = zp * zd * J0 / (zn * zd + zn * zp + zd * zp);
b = zn * zd * J0 / (zn * zd + zn * zp + zd * zp);
c = -zn * zd * zp * J0 / (zn * zd + zn * zp + zd * zp);

% Definition of the elements of the vector UP 
% in x=0 and x=L for boundary conditions.

Up0 = [0; 0; 0; a; b];
UpL = [ c*L; c*L; c*L; a; b];

% Definition of the matrix A 

A = [ 0,       0,       0,       -zn,         0;
    
    0,       0,       0,       0,         -zp;
    
    -gamma_n, gamma_n-gamma_p, gamma_p, zd, zd;
   
    -yn-yr, yn, yr, 0, 0; 
    
    yr, yp, -yp-yr, 0, 0 ];




%Calculation of the eigenvalues and associated eigenvalues of Matrix A. 

[Y,X] = eig(A);

x1 = X(1,1); 
x2 = X(2,2); 
x3 = X(3,3); 
x4 = X(4,4); 
x5 = X(5,5);

y1 = Y(:,1);
y2 = Y(:,2);
y3 = Y(:,3);
y4 = Y(:,4);
y5 = Y(:,5);

U10 = y1; 
U20 = y2; 
U30 = y3; 
U40 = y4; 
U50 = y5;

U1L = y1.* exp( x1.* L);
U2L = y2.* exp( x2.* L);
U3L = y3.* exp( x3.* L);
U4L = y4.* exp( x4.* L);
U5L = y5.* exp( x5.* L);

% Solution for the Cramer System 


D = [ U10(1) - U10(2), U20(1) - U20(2), U30(1) - U30(2), U40(1) - U40(2), U50(1) - U50(2);
    
    U1L(3) - U1L(2), U2L(3) - U2L(2), U3L(3) - U3L(2), U4L(3) - U4L(2), U5L(3) - U5L(2); 
    
    U10(3)         , U20(3)         , U30(3)         , U40(3)          , U50(3)        ;
    
    U10(5)         , U20(5)         , U30(5)         , U40(5)          , U50(5)        ;
    
    U1L(4)         , U2L(4)         , U3L(4)         , U4L(4)          , U5L(4)        ];







K1N = [ 0              , U20(1) - U20(2), U30(1) - U30(2), U40(1) - U40(2), U50(1) - U50(2);
    
        0              , U2L(3) - U2L(2), U3L(3) - U3L(2), U4L(3) - U4L(2), U5L(3) - U5L(2); 
    
        0              , U20(3)         , U30(3)         , U40(3)          , U50(3)        ;
    
       -b              , U20(5)         , U30(5)         , U40(5)          , U50(5)        ;
    
       -a              , U2L(4)         , U3L(4)         , U4L(4)          , U5L(4)        ];
    
  
   
   
   
K2N = [ U10(1) - U10(2), 0              , U30(1) - U30(2), U40(1) - U40(2), U50(1) - U50(2);
    
       U1L(3) - U1L(2),  0              , U3L(3) - U3L(2), U4L(3) - U4L(2), U5L(3) - U5L(2); 
    
       U10(3)         ,  0              , U30(3)         , U40(3)          , U50(3)        ;
    
       U10(5)         , -b              , U30(5)         , U40(5)          , U50(5)        ;
    
       U1L(4)         , -a              , U3L(4)         , U4L(4)          , U5L(4)        ];
   

   
   
   
K3N  = [ U10(1) - U10(2), U20(1) - U20(2), 0              , U40(1) - U40(2), U50(1) - U50(2);
    
         U1L(3) - U1L(2), U2L(3) - U2L(2), 0              , U4L(3) - U4L(2), U5L(3) - U5L(2); 
    
         U10(3)         , U20(3)         , 0              , U40(3)          , U50(3)        ;
    
         U10(5)         , U20(5)         , -b             , U40(5)          , U50(5)        ;
    
         U1L(4)         , U2L(4)         , -a             , U4L(4)          , U5L(4)        ];
 
     
     


K4N = [ U10(1) - U10(2), U20(1) - U20(2), U30(1) - U30(2), 0               , U50(1) - U50(2);
    
        U1L(3) - U1L(2), U2L(3) - U2L(2), U3L(3) - U3L(2), 0               , U5L(3) - U5L(2); 
    
        U10(3)         , U20(3)         , U30(3)         , 0               , U50(3)        ;
     
        U10(5)         , U20(5)         , U30(5)         , -b              , U50(5)        ;
    
        U1L(4)         , U2L(4)         , U3L(4)         , -a              , U5L(4)        ];   
    
    
    
    
    
K5N = [ U10(1) - U10(2), U20(1) - U20(2), U30(1) - U30(2), U40(1) - U40(2), 0              ;
    
        U1L(3) - U1L(2), U2L(3) - U2L(2), U3L(3) - U3L(2), U4L(3) - U4L(2), 0              ; 
    
        U10(3)         , U20(3)         , U30(3)         , U40(3)          ,0              ;
    
        U10(5)         , U20(5)         , U30(5)         , U40(5)          , -b            ;
    
        U1L(4)         , U2L(4)         , U3L(4)         , U4L(4)          , -a            ];
    
    

    
K1 = det(K1N) / det(D);
K2 = det(K2N) / det(D); 
K3 = det(K3N) / det(D); 
K4 = det(K4N) / det(D); 
K5 = det(K5N) / det(D);

% Equations with values of the constants K1, K2, K3, K4, K5 obtained
% previously.



U0 = K1*U10 + K2*U20 + K3*U30 + K4*U40 + K5*U50 + Up0; 
UL = K1*U1L + K2*U2L + K3*U3L + K4*U4L + K5*U5L + UpL;

V10 = U0(1); 
V20 = U0(2); 
V30 = U0(3); 
I10 = U0(4); 
I30 = U0(5);


V1L = UL(1);
V2L = UL(2);
V3L = UL(3);
I1L = UL(4); 
I3L = UL(5);


% Calcuation of the Real and Imaginary portion of the impedance.

Zre(l) = real ((V10 - V3L) / J0);
Zim(l) = imag ((V10 - V3L) / J0);

% Definition of the real and imaginary part of the Total Capacitance.

Cre(l) = real (J0/ (1i * f *(V10 - V3L) ));
Cim(l) = imag (J0/ (1i * f *(V10 - V3L) ));

% Definition of the Phase.

Phi(l) = -180 / pi* atan(Zim(l) / Zre(l)); 

% Definition of the Frequency.

fimp(l) = f;

l = l+1;

end


% Representation of the Impedance, Capacitance, and Phase Spectra.


Z1 = Zre (1: max(l) - 1);
Z2 = Zim (1: max(l) - 1);

ZRE = Z1' ;
ZIM = Z2' ; 

Cap1 = Cre (1:max(l) - 1);
Cap2 = Cim (1:max(l) - 1);

PHI = Phi(1:max(l) - 1);
f   = fimp(1:max(l) - 1);

phi = PHI' ;

plot (ZRE, -ZIM) 
xlabel ( 'Z_{re}', 'FontSize', 16)
ylabel ( 'Z_{im}', 'FontSize', 16)
title  ( '\it{Nyquist Plot}', 'FontSize', 16)

figure;
loglog (f', Cre)
xlabel ('f', 'FontSize', 16)
ylabel ('C_{re}', 'FontSize', 16)
title  ('\it{Capacitance Plot}', 'FontSize', 16)

figure;
plot   (log(f'), phi)
xlabel ('log(f)', 'FontSize', 16)
ylabel ('phase (degree)', 'FontSize', 16)




