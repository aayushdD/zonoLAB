function [A,b,Ae,be] = vrep_to_Hrep_mpt3(V)
P=Polyhedron('V',V);
P.minVRep();
P.minHRep();
A  = P.A;      % inequality normals
b  = P.b;      % inequality offsets
Ae = P.Ae;     % equality normals (present if lower-dimensional)
be = P.be;  
