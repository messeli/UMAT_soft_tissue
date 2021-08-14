# -*- coding: mbcs -*-
#
# Abaqus/CAE Release 2017 replay file
# Internal Version: 2016_09_28-00.54.59 126836
# Run by mehme on Sat Aug 14 15:35:47 2021
#

# from driverUtils import executeOnCaeGraphicsStartup
# executeOnCaeGraphicsStartup()
#: Executing "onCaeGraphicsStartup()" in the site directory ...
from abaqus import *
from abaqusConstants import *
session.Viewport(name='Viewport: 1', origin=(0.0, 0.0), width=150.9375, 
    height=87.4796295166016)
session.viewports['Viewport: 1'].makeCurrent()
session.viewports['Viewport: 1'].maximize()
from caeModules import *
from driverUtils import executeOnCaeStartup
executeOnCaeStartup()
openMdb('myUmat_deneme.cae')
#: The model database "C:\Users\mehme\Google_Drive\Soft_Tissue_FEM\FEM_ABAQUS\UMAT\Repos\messeli\UMAT_soft_tissue\On Abaqus\myUmat_deneme.cae" has been opened.
session.viewports['Viewport: 1'].setValues(displayedObject=None)
session.viewports['Viewport: 1'].partDisplay.geometryOptions.setValues(
    referenceRepresentation=ON)
p = mdb.models['Model-1'].parts['Part-1']
session.viewports['Viewport: 1'].setValues(displayedObject=p)
#--- Recover file: 'myUmat_deneme.rec' ---
# -*- coding: mbcs -*- 
from part import *
from material import *
from section import *
from assembly import *
from step import *
from interaction import *
from load import *
from mesh import *
from optimization import *
from job import *
from sketch import *
from visualization import *
from connectorBehavior import *
mdb.jobs['Job-1'].setValues(
    userSubroutine='C:\\Users\\mehme\\Google_Drive\\Soft_Tissue_FEM\\FEM_ABAQUS\\UMAT\\Repos\\messeli\\UMAT_soft_tissue\\UMAT_Combined.for')
mdb.jobs['Job-1'].submit(consistencyChecking=OFF)
mdb.jobs['Job-1']._Message(ERROR, {
    'message': 'Problem during compilation - C:\\Users\\mehme\\Google_Drive\\Soft_Tissue_FEM\\FEM_ABAQUS\\UMAT\\Repos\\messeli\\UMAT_soft_tissue\\UMAT_Combined.for', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(JOB_ABORTED, {
    'message': 'Problem during compilation - C:\\Users\\mehme\\Google_Drive\\Soft_Tissue_FEM\\FEM_ABAQUS\\UMAT\\Repos\\messeli\\UMAT_soft_tissue\\UMAT_Combined.for', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1'].submit(consistencyChecking=OFF)
mdb.jobs['Job-1']._Message(ERROR, {
    'message': 'Problem during compilation - C:\\Users\\mehme\\Google_Drive\\Soft_Tissue_FEM\\FEM_ABAQUS\\UMAT\\Repos\\messeli\\UMAT_soft_tissue\\UMAT_Combined.for', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(JOB_ABORTED, {
    'message': 'Problem during compilation - C:\\Users\\mehme\\Google_Drive\\Soft_Tissue_FEM\\FEM_ABAQUS\\UMAT\\Repos\\messeli\\UMAT_soft_tissue\\UMAT_Combined.for', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1'].submit(consistencyChecking=OFF)
mdb.jobs['Job-1']._Message(ERROR, {
    'message': 'Problem during compilation - C:\\Users\\mehme\\Google_Drive\\Soft_Tissue_FEM\\FEM_ABAQUS\\UMAT\\Repos\\messeli\\UMAT_soft_tissue\\UMAT_Combined.for', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(JOB_ABORTED, {
    'message': 'Problem during compilation - C:\\Users\\mehme\\Google_Drive\\Soft_Tissue_FEM\\FEM_ABAQUS\\UMAT\\Repos\\messeli\\UMAT_soft_tissue\\UMAT_Combined.for', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1'].submit(consistencyChecking=OFF)
mdb.jobs['Job-1']._Message(ERROR, {
    'message': 'Problem during compilation - C:\\Users\\mehme\\Google_Drive\\Soft_Tissue_FEM\\FEM_ABAQUS\\UMAT\\Repos\\messeli\\UMAT_soft_tissue\\UMAT_Combined.for', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(JOB_ABORTED, {
    'message': 'Problem during compilation - C:\\Users\\mehme\\Google_Drive\\Soft_Tissue_FEM\\FEM_ABAQUS\\UMAT\\Repos\\messeli\\UMAT_soft_tissue\\UMAT_Combined.for', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1'].submit(consistencyChecking=OFF)
mdb.jobs['Job-1']._Message(ERROR, {
    'message': 'Problem during compilation - C:\\Users\\mehme\\Google_Drive\\Soft_Tissue_FEM\\FEM_ABAQUS\\UMAT\\Repos\\messeli\\UMAT_soft_tissue\\UMAT_Combined.for', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(JOB_ABORTED, {
    'message': 'Problem during compilation - C:\\Users\\mehme\\Google_Drive\\Soft_Tissue_FEM\\FEM_ABAQUS\\UMAT\\Repos\\messeli\\UMAT_soft_tissue\\UMAT_Combined.for', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1'].submit(consistencyChecking=OFF)
#--- End of Recover file ------
a = mdb.models['Model-1'].rootAssembly
session.viewports['Viewport: 1'].setValues(displayedObject=a)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(loads=ON, bcs=ON, 
    predefinedFields=ON, connectors=ON, optimizationTasks=OFF, 
    geometricRestrictions=OFF, stopConditions=OFF)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(step='Step-1')
#: Warning: Cannot continue yet--complete the step or cancel the procedure.
session.viewports['Viewport: 1'].view.setValues(nearPlane=3.23376, 
    farPlane=6.55667, width=2.27389, height=1.25285, cameraPosition=(3.41964, 
    3.32199, -1.44033), cameraUpVector=(-0.640951, 0.57735, 0.505815), 
    cameraTarget=(0.224898, 0.444259, 1.08084))
a = mdb.models['Model-1'].rootAssembly
f1 = a.instances['Part-1-1'].faces
faces1 = f1.getSequenceFromMask(mask=('[#20 ]', ), )
region = a.Set(faces=faces1, name='Set-1')
mdb.models['Model-1'].EncastreBC(name='BC-1', createStepName='Step-1', 
    region=region, localCsys=None)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(mesh=ON, loads=OFF, 
    bcs=OFF, predefinedFields=OFF, connectors=OFF)
session.viewports['Viewport: 1'].assemblyDisplay.meshOptions.setValues(
    meshTechnique=ON)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(mesh=OFF, 
    optimizationTasks=ON, geometricRestrictions=ON, stopConditions=ON)
session.viewports['Viewport: 1'].assemblyDisplay.meshOptions.setValues(
    meshTechnique=OFF)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(
    optimizationTasks=OFF, geometricRestrictions=OFF, stopConditions=OFF)
session.viewports['Viewport: 1'].setValues(displayedObject=None)
a = mdb.models['Model-1'].rootAssembly
session.viewports['Viewport: 1'].setValues(displayedObject=a)
mdb.jobs['Job-1'].submit(consistencyChecking=OFF)
#: Abaqus Error: Detected lock file Job-1.lck. Please confirm that no other applications are attempting to write to the output database associated with this job before removing the lock file and resubmitting.
#: Abaqus/Analysis exited with error(s).
#: The job input file "Job-1.inp" has been submitted for analysis.
mdb.jobs['Job-1'].submit(consistencyChecking=OFF)
#: Abaqus Error: Detected lock file Job-1.lck. Please confirm that no other applications are attempting to write to the output database associated with this job before removing the lock file and resubmitting.
#: Abaqus/Analysis exited with error(s).
#: The job input file "Job-1.inp" has been submitted for analysis.
mdb.jobs['Job-1'].submit(consistencyChecking=OFF)
#: Abaqus Error: Detected lock file Job-1.lck. Please confirm that no other applications are attempting to write to the output database associated with this job before removing the lock file and resubmitting.
#: Abaqus/Analysis exited with error(s).
#: The job input file "Job-1.inp" has been submitted for analysis.
mdb.save()
#: The model database has been saved to "C:\Users\mehme\Google_Drive\Soft_Tissue_FEM\FEM_ABAQUS\UMAT\Repos\messeli\UMAT_soft_tissue\On Abaqus\myUmat_deneme.cae".
