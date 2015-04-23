import sys
cubitPyPath = "/home/deluca/C_Code/femdoc/tests/battery/CubitPyFuncs"
sys.path.append(cubitPyPath)

import cubit
reload(CubGenericPyFuncs)
reload(Cub2DMeshPyFuncs)
import CubGenericPyFuncs
import Cub2DMeshPyFuncs

cubit.cmd('reset')

length     = 3.0
height     = 2.0

lengthInts = 60
heightInts = 40

lBias      = 1.0
hBias      = 1.0

loadRad    = max([((length/lengthInts)*(2.0**0.5)),0.125])


Dimensions   = [length    , height    , 0.0, 0.0]
Intervals    = [lengthInts, heightInts]
Biaii        = [lBias     , hBias    ]
Geometry_Dictionary = Cub2DMeshPyFuncs.createMeshSquare(Dimensions,Intervals,Biaii)


Cyl_SIds     = Geometry_Dictionary['surfaces']
CylXmin_CIds = Geometry_Dictionary['curves']['xbounds']['lower'][0]
CylXmax_CIds = Geometry_Dictionary['curves']['xbounds']['upper'][0]
CylYmin_CIds = Geometry_Dictionary['curves']['ybounds']['lower'][0]
CylYmax_CIds = Geometry_Dictionary['curves']['ybounds']['upper'][0]
Origin_VIds  = CubGenericPyFuncs.findCommonRelatives(CylXmin_CIds, CylYmin_CIds,'vertex')
Tip_VIds     = CubGenericPyFuncs.findCommonRelatives(CylXmax_CIds, CylYmax_CIds,'vertex')

[Input_EIds , Input_NIds ] = CubGenericPyFuncs.findEdgesInGeomWithinDistanceOfCoordinate(CylXmax_CIds[0],'curve',[length,height/2.0,0.0],loadRad)

# --------------------------------- Generate Blocks ----------------------------------->
cubit.cmd('set Node Constraint off')
for surf_Id in Cyl_SIds:
    cubit.cmd('block %i surface %i'      % (1,surf_Id))
cubit.cmd('block %i name "Box"'          % (1))
cubit.cmd('block %i element type quad4'  % (1))
# --------------------------------- Generate Blocks end ------------------------------->

# --------------------------------- Generate Side Sets -------------------------------->
for curve_Id in CylYmin_CIds:
    cubit.cmd('sideset %i curve %i'                 % (1 ,curve_Id))
cubit.cmd('sideset %i name "XAxis"'                 % (1 ))

for curve_Id in CylXmax_CIds:
    cubit.cmd('sideset %i curve %i'                 % (2 ,curve_Id))
cubit.cmd('sideset %i name "XmaxPlane"'             % (2 ))

for curve_Id in CylYmax_CIds:
    cubit.cmd('sideset %i curve %i'                 % (3 ,curve_Id))
cubit.cmd('sideset %i name "YmaxPlane"'             % (3 ))

for curve_Id in CylXmin_CIds:
    cubit.cmd('sideset %i curve %i'                 % (4 ,curve_Id))
cubit.cmd('sideset %i name "YAxis"'                 % (4 ))

for edge_Id in Input_EIds:
    cubit.cmd('sideset %i edge %i'                  % (5 ,edge_Id))
cubit.cmd('sideset %i name "Input"'                 % (5 ))

# --------------------------------- Generate Side Sets end ---------------------------->

# --------------------------------- Generate Node Sets -------------------------------->
CubGenericPyFuncs.createNodeSetFromSideSet(1,1)
CubGenericPyFuncs.createNodeSetFromSideSet(2,2)
CubGenericPyFuncs.createNodeSetFromSideSet(3,3)
CubGenericPyFuncs.createNodeSetFromSideSet(4,4)

for node_Id in Input_NIds:
    cubit.cmd('nodeset %i node %i'      % (5 ,node_Id))
cubit.cmd('nodeset %i name "Input"'     % (5 ))

cubit.cmd('nodeset %i vertex %i'        % (6 ,Origin_VIds[0][0][0]))
cubit.cmd('nodeset %i name "Origin"'    % (6 ))

# --------------------------------- Generate Node Sets end ---------------------------->

# Export the mesh
cubit.cmd('export mesh "beam2D_%s.g" dimension 2 overwrite' % ('%ix%i' % (lengthInts,heightInts)))


sys.path.remove(cubitPyPath)
