program lab3;
type mas=array [0..3]of double;
var
  xBegin,yBegin,zBegin,xEnd,M:double;
  xOld,yOld,zOld:double;
  xNew,yNew,zNew:double;
  xWrite,rho:double;
  i,j:integer;
  eps,h:double;
  k1_new,k2_new,l1_new,l2_new:double;
  p:double;


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

  k1_new:= yPrime(x0+h0*a1, y0+h0*a3*k1+a2*h0*k2, z0+a3*h*l1+a2*h0*l2);
  l1_new:= zPrime(x0+h0*a1, y0+h0*a3*k1+a2*h0*k2, z0+a3*h*l1+a2*h0*l2);
  k2_new:= yPrime(x0+h0*a1, y0+h0*a2*k1+a3*h0*k2, z0+a2*h*l1+a3*h0*l2);
  l1_new:= zPrime(x0+h0*a1, y0+h0*a3*k1+a2*h0*k2, z0+a3*h*l1+a2*h0*l2);

  xN:=x0+h0;
  yN:=(h0/2)*(k1_new+k2_new)+y0;
  zN:=(h0/2)*(l1_new+l2_new)+z0;

end;

procedure kPrime();
var
BEGIN

END;

procedure fullfilM();
var m,n:INTEGER;
matrix:array[1..5,1..5] of DOUBLE;
BEGIN
  for n:= 0 to 4 do
    for m:=0 to 4 do
       BEGIN
         matrix[n,m]:= kPrime();
       END;
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



    end;

  end.
