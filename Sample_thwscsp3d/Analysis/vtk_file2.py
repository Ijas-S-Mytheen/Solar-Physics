import numpy as np
from tqdm import tqdm
import vtk

# Load meshgrid values
theta3 = np.load('/home/ijas/Fortran_IIA/testing_code/res/Components2/theta.npy')
theta2 = np.swapaxes(theta3,1,2)
theta4 = np.swapaxes(theta2,0,1)
theta = np.swapaxes(theta4,0,2)
phi3 = np.load('/home/ijas/Fortran_IIA/testing_code/res/Components2/phi.npy')
phi2 = np.swapaxes(phi3,1,2)
phi4 = np.swapaxes(phi2,0,1)
phi = np.swapaxes(phi4,0,2)
r3 = np.load('/home/ijas/Fortran_IIA/testing_code/res/Components2/r.npy')
r2 = np.swapaxes(r3,1,2)
r4 = np.swapaxes(r2,0,1)
r = np.swapaxes(r4,0,2)

# Load magnetic field components
Bt3 = np.load('/home/ijas/Fortran_IIA/testing_code/res/Components2/Bt.npy')
Bt2 = np.swapaxes(Bt3,1,2)
Bt4 = np.swapaxes(Bt2,0,1)
Bt = np.swapaxes(Bt4,0,2)
Bp3 = np.load('/home/ijas/Fortran_IIA/testing_code/res/Components2/Bp.npy')
Bp2 = np.swapaxes(Bp3,1,2)
Bp4 = np.swapaxes(Bp2,0,1)
Bp = np.swapaxes(Bp4,0,2)
Br3 = np.load('/home/ijas/Fortran_IIA/testing_code/res/Components2/Br.npy')
Br2 = np.swapaxes(Br3,1,2)
Br4 = np.swapaxes(Br2,0,1)
Br = np.swapaxes(Br4,0,2)

X = r*np.sin(theta)*np.cos(phi)
Y = r*np.sin(theta)*np.sin(phi)
Z = r*np.cos(theta)

points = vtk.vtkPoints()
grid = vtk.vtkStructuredGrid()
grid.SetDimensions(511,160, 288)

with tqdm(desc="Creating points", total=160*288*511) as pbar_points:
    for i in range(288):
        for j in range(160):
            for k in range(511):
                points.InsertNextPoint(X[i, j, k], Y[i, j, k], Z[i, j, k])
                pbar_points.update(1)

    grid.SetPoints(points)

with tqdm(desc="Inserting magnetic field components", total=160*288*511) as pbar_components:    
    Br_array = vtk.vtkFloatArray()
    Br_array.SetName("Br")
    Btheta_array = vtk.vtkFloatArray()
    Btheta_array.SetName("Bt")
    Bphi_array = vtk.vtkFloatArray()
    Bphi_array.SetName("Bp")

    for i in range(288):
        for j in range(160):
            for k in range(511):
                Br_array.InsertNextValue(Br[i, j, k])
                Btheta_array.InsertNextValue(Bt[i, j, k])
                Bphi_array.InsertNextValue(Bp[i, j, k])
                pbar_components.update(1)

grid.GetPointData().AddArray(Br_array)
grid.GetPointData().AddArray(Btheta_array)
grid.GetPointData().AddArray(Bphi_array)

with tqdm(desc="Adding the arrays to the grid's point data", total=160*288*511) as pbar_components:  

    Bx = Br * np.sin(theta) * np.cos(phi) + \
        Bt * np.cos(theta) * np.cos(phi) - \
        Bp * np.sin(phi)

    By = Br * np.sin(theta) * np.sin(phi) + \
        Bt * np.cos(theta) * np.sin(phi) + \
        Bp * np.cos(phi)

    Bz = Br * np.cos(theta) - Bt * np.sin(theta)

    B_vector_array = vtk.vtkFloatArray()
    B_vector_array.SetNumberOfComponents(3)
    B_vector_array.SetName("B_vector")

    for i in range(288):
        for j in range(160):
            for k in range(511):
                B_vector_array.InsertNextTuple((Bx[i, j, k],By[i, j, k], Bz[i, j, k]))
                pbar_components.update(1)

grid.GetPointData().AddArray(B_vector_array)


with tqdm(desc="Add the magnetic field components to the grid", total=160*288*511) as pbar_components:  
    Bx_array = vtk.vtkFloatArray()
    Bx_array.SetName("Bx")
    By_array = vtk.vtkFloatArray()
    By_array.SetName("By")
    Bz_array = vtk.vtkFloatArray()
    Bz_array.SetName("Bz")

    for i in range(288):
        for j in range(160):
            for k in range(511):
                Bx_array.InsertNextValue(Bx[i, j, k])
                By_array.InsertNextValue(By[i, j, k])
                Bz_array.InsertNextValue(Bz[i, j, k])
                pbar_components.update(1)

# Add the arrays to the grid's point data
grid.GetPointData().AddArray(Bx_array)
grid.GetPointData().AddArray(By_array)
grid.GetPointData().AddArray(Bz_array)

# Create a writer and save the VTK file
writer = vtk.vtkStructuredGridWriter()
writer.SetFileName("B_actual.vtk")
writer.SetInputData(grid)
writer.Write()

print('max of bx:',np.max(Bx),'min of bx:',np.min(Bx))
print('max of by:',np.max(By),'min of bx:',np.min(By))
print('max of bz:',np.max(Bz),'min of bx:',np.min(Bz))
