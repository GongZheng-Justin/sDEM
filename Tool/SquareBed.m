clc;clear;
format long

nLx=216;
nLy=18;
nLz=54;
Diam=0.1;
CoordFileOutName='../ACM/Restart/FixedSpheresCoord.dat';
FixedXMFName='../ACM/Restart/FixedBed.xmf';

%================================================================================
xlx= Diam*nLx;
zlz= Diam*nLz;
int_prec='integer*4';
int_byte=4;
real_prec='real*8';
real_byte=8;

%PosR%x, PosR%y, PosR%z, Diameter, pType
nSumMax=nLx*nLy*nLz;
Prtcl_Pos=zeros(3,nSumMax);

nPTotal=0;
for nkx=1:nLx
  xCoordNow=Diam*(nkx-0.5);
  for nky=1:nLy
    yCoordNow=Diam*(nky-0.5);
    for nkz=1:nLz
      nPTotal=nPTotal+1;
      zCoordNow=Diam*(nkz-0.5);
      Prtcl_Pos(:,nPTotal)=[xCoordNow;yCoordNow;zCoordNow];
    end
  end
end

fid=fopen(CoordFileOutName,'w');
fwrite(fid,Prtcl_Pos(:,1:nPTotal),real_prec);
fwrite(fid,Diam*ones(nPTotal,1),real_prec);
fwrite(fid,ones(nPTotal,1),int_prec);
fclose(fid);

% write the xmf file
fid3=fopen(FixedXMFName,'wt');
disp_xmf=0;
fprintf(fid3,'<?xml version="1.0" ?>\n');
fprintf(fid3,'<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>\n');
fprintf(fid3,'<Xdmf xmlns:xi="http://www.w3.org/2001/XInclude" Version="2.0">\n');
fprintf(fid3,'<Domain>\n');
fprintf(fid3,'    <Grid Name="TimeSeries" GridType="Collection" CollectionType="Temporal">\n');
fprintf(fid3,'        <Time TimeType="List">\n');
fprintf(fid3,'            <DataItem Format="XML" NumberType="Int" Dimensions="     1">\n');
fprintf(fid3,'                    0            </DataItem>\n');
fprintf(fid3,'        </Time>\n');
fprintf(fid3,'        <Grid Name="T0000000000" GridType="Uniform">\n');
fprintf(fid3,'            <Topology TopologyType="Polyvertex" NodesPerElement="        0"/>\n');
fprintf(fid3,'            <Geometry GeometryType="XYZ">\n');
fprintf(fid3,'                <DataItem Format="Binary" DataType="Float" Precision="8" Endian="Native" Dimensions=" 3      %d" Seek="   %d">\n',nPTotal,disp_xmf);
fprintf(fid3,'                    FixedSpheresCoord.dat\n');
fprintf(fid3,'                </DataItem>\n');
fprintf(fid3,'            </Geometry>\n');
disp_xmf=disp_xmf+ nPTotal*real_byte*3;
fprintf(fid3,'            <Attribute Type="Scalar" Center="Node" Name="Diameter">\n');
fprintf(fid3,'                <DataItem Format="Binary" DataType="Float" Precision="8" Endian="Native" Dimensions=" 1      %d" Seek="   %d">\n',nPTotal,disp_xmf);
fprintf(fid3,'                    FixedSpheresCoord.dat\n');
fprintf(fid3,'                </DataItem>\n');
fprintf(fid3,'            </Attribute>\n');
disp_xmf=disp_xmf+ nPTotal*real_byte;
fprintf(fid3,'            <Attribute Type="Scalar" Center="Node" Name="Type">');
fprintf(fid3,'                <DataItem Format="Binary" DataType="int" Precision="4" Endian="Native" Dimensions=" 1      %d" Seek="   %d">\n',nPTotal,disp_xmf);
fprintf(fid3,'                    FixedSpheresCoord.dat\n');
fprintf(fid3,'                </DataItem>\n');
fprintf(fid3,'            </Attribute>\n');
fprintf(fid3,'        </Grid>\n');
fprintf(fid3,'    </Grid>\n');
fprintf(fid3,'</Domain>\n');
fprintf(fid3,'</Xdmf>\n');
fclose(fid3);

disp(['nFixed= ',num2str(nPTotal)]);
format short
