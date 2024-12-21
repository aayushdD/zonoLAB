

% The obj=orderReduction(obj,setType,orderReductionTechnique,final_order)
% method takes the object, type of object (zonotope or other type of set),
% reduction technique (exact or something else) and final order we desire.
% setType and orderReductionTechnique are optional arguments and have
% default values of 'Zonotope' and 'Exact'. Currently the reduction
% techniques are provided only for zonotopes. Hence, exact, inner or outer
% approximation for zonotopes can be done. There could be some
% inconsistencies. The object is an object from abstract Zono class.
