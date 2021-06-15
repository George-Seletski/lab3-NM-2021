﻿program lab3;
var
  xWrite,rho:double;
  i,j:integer;
  eps,h:double;
  k1_new,k2_new,l1_new,l2_new:double;
  p:double;
  A1,A2:double;
  det_G:double;

function yPrime(x,y,z:double):double;
 begin
   yPrime:=y-z*z+2*x;
 end;
function zPrime(x,y,z:double):double;
 begin
   zPrime:=2*y + z;
 end;

procedure rkm(h0,x0,y0,z0,k1_new,k2_new,l1_new,l2_new:double;var xN,yN,zN:double); 
 var k1,k2:double;
     a1,a2,a3:double;
     l1,l2:double;
     h3:double;
 begin
   a1:= ( (1/2) - (sqrt(3)/6) );
   a2:= ( (1/4) - (sqrt(3)/6) );
   a3:= 1/4;
   h3:=h0/3;

   k1:= yPrime(x0,y0,z0);
   l1:= zPrime(x0,y0,z0);
   k2:= yPrime(x0+h3,y0+h3*k1,z0+h3*l1);
   l2:= zPrime(x0+h3,y0+h3*k1,z0+h3*l1); 

   k1_new:= yPrime(x0+h0*a1, y0+h0*a3*k1+a2*h0*k2, z0+a3*h*l1+a2*h0*l2); // k1 = f1
   l1_new:= zPrime(x0+h0*a1, y0+h0*a3*k1+a2*h0*k2, z0+a3*h*l1+a2*h0*l2); // l1 = g1
   k2_new:= yPrime(x0+h0*a1, y0+h0*a2*k1+a3*h0*k2, z0+a2*h*l1+a3*h0*l2); // k2 = f2
   l1_new:= zPrime(x0+h0*a1, y0+h0*a3*k1+a2*h0*k2, z0+a3*h*l1+a2*h0*l2); // l2 = g2

   xN:=x0+h0;
   yN:=(h0/2)*(k1_new+k2_new)+y0;
   zN:=(h0/2)*(l1_new+l2_new)+z0;


   F[0]:=k1_new - k1;
   F[1]:=k2_new - k2;
   F[2]:=l1_new - l1;
   F[3]:=l2_new - l2;
   
   A1:= a1;
   A2:= a2;
  

 end;


procedure kPrime();
var
BEGIN

END;

procedure G(h:double;);
var m,n:INTEGER;
ma:array[1..4,1..4] of DOUBLE;
det,tmp_d1,tmp_d2,tmp_d3,tmp_d4:double;

BEGIN //заполняем матрицу G
for m := 1 to 4 do
  for  n:= 1 to 4 do
  begin
    if (m = n) then 
      begin
        {
        ma[m,n]:= 1 - (h/4);
        }
      end
    else
      begin
        {
          ma[1,2]:= -1 * A2 * h;
          ma[1,3]:= -1 * 0.5 * h;
          ma[1,4]:= -2 * A2 * h;

          ma[2,1]:= ma[1,2];
          ma[2,3]:= ma[1,4];
          ma[2,4]:= ma[1,3];

          m[3,1]:= -1 * (zn + (h/4)*L1*A2*h*L2)*(h/2);
          m[3,2]:= -1 * (zn + (h/4)*L1*A2*h*L2)*(A2*h/2);
          m[3,4]:= ma[1,2];
          
          ma[4,1]:= -2 * (zn + (h/4)*L1*A2*h*L2)*A2*h;
          ma[4,2]:= ma[4,1];
          ma[4,3]:= m[3,4];

        }    
    end;
  //транспонирование матрицы
  for m := 1 to 4 do
  for  n:= 1 to 4 do
  begin
  if(m <> n){
      ma[m,n]:= ma[n,m];
  }
  end;

 tmp_d1:= ma[1,1]*(ma[2,2]*ma[3,3]*ma[4,4] + ma[2,3]*ma[3,4]*ma[4,2] + ma[2,4]+ma[3,2]+ma[4,3] - ma[2,4]*ma[3,3]*ma[4,2] - ma[2,2]*ma[3,4]*ma[4,3] - ma[4,4]*ma[2,3]*ma[3,2]);
 tmp_d2:= ma[1,2]*(ma[2,1]*ma[3,3]*ma[4,4] + ma[2,4]*ma[3,1]*ma[4,3] + ma[4,1]+ma[2,3]+ma[3,4] - ma[2,4]*ma[3,3]*ma[4,1] - ma[2,1]*ma[3,4]*ma[4,3] - ma[4,4]*ma[2,3]*ma[3,1]);
 tmp_d3:= ma[1,3]*(ma[2,1]*ma[3,2]*ma[4,4] + ma[2,4]*ma[3,1]*ma[4,2] + ma[4,1]+ma[2,2]+ma[3,4] - ma[2,4]*ma[3,2]*ma[4,1] - ma[2,1]*ma[3,4]*ma[4,2] - ma[4,4]*ma[2,2]*ma[3,1]);
 tmp_d4:= ma[1,4]*(ma[2,1]*ma[3,2]*ma[4,3] + ma[2,2]*ma[3,3]*ma[4,1] + ma[2,3]+ma[3,1]+ma[4,2] - ma[2,3]*ma[3,2]*ma[4,1] - ma[2,1]*ma[3,3]*ma[4,2] - ma[4,3]*ma[2,3]*ma[3,1]);

 det:= tmp_d1-tmp_d2+tmp_d3-tmp_d4; // детерминант матрицы G

 det_G:=det;

END;
procedure print(x,y,z:double;var p:double);
var i:integer;
    eps:double =0.00001;
begin
  if(abs(x-p)<eps)then
    begin
      writeln(x:10:4,'  ',y:10:4,'  ',z:10:4);
      //writeln('(',x:0:4,'; ',y,') ');//,z);
      p:=p+0.05;
    end;
end;



begin





end.
