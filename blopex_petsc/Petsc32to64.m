function Petsc32to64(dir,file,comp)
% example >>
% Petsc32to64('C:/BLOPEX/blopex_petsc/driver_fiedler/','L-matrix','double')
% 
namein=[dir file '.petsc']; 
A=PetscBinaryRead(namein,'int32',comp);
nameout=[dir file '_64.petsc'];
PetscBinaryWrite(nameout,'int64',A); 
end