#==============================================================================
# Sensitivity Analysis of Water transport to K(h)
# Created By Jan De Pue (2019) at Ghent University
# Reference:
#
# General idea:
#   - modify hydraulic conductivity in the table input for Hydrus
#   - use a wrapper around Hydrus1D (PC-progress) execute Hydrus
#   - determine the sensitivity at this point (run PLUS and run MIN)
#   - repeat until pleased
#==============================================================================

import pylab
import numpy
import os
import fnmatch
import glob
from distutils import dir_util
import shutil
import copy
import time

from RetentionConductivityCapacity_Funky import *
from Hydrus_Funky import *

from matplotlib.backends.backend_pdf import PdfPages
figlist=[]

figsize1=numpy.array([10./numpy.sqrt(2),10])
figsize2=numpy.array([10,10./numpy.sqrt(2)])
dpi=300

font = {'family' : 'monospace',
        'size'   : 8}
pylab.rc('font', **font)

pylab.ioff()


def TestSensitivity(hSens,Perturb,RetPar = -1 ,ConPar = -1,
                    GuessRuntime = -1, RunREF = False):
    DoLog = False
    SleepTime = 10.0
    
    #==============================================================================
    # Locate hydrus project
    #==============================================================================

    HydrusRoot = '/home/jan/git/Hydrus/H1D_Src'
    HydrusExec = 'HydrusForLinux'
    HydrusLev1 = 'LEVEL_01.DIR' ## should be located in the folder where it is run from

    RefProjectPath = "/home/jan/git/KsIsIrrelevant/HydrusProjects/TestSens_MeteoSimple2"
    MaterRoot = ''

    #==============================================================================
    # Prepare 
    #==============================================================================

    Retention_model = h2theta_VanGenuchten4 # thetar,thetas,alpha1,n1,=Par
    Conductivity_model = h2K_Mualem4 # Ks,l,alpha1,n1=Par
    Capacity_model = h2C_VanGenuchten4 #  thetar,thetas,alpha1,n1 = Par

    if str(RetPar) == '-1':
        # LOAM
        RetPar = numpy.array([[0.078, 0.43, 0.036, 1.56,],])

    if str(ConPar) == '-1':
        # LOAM
        ConPar = numpy.array([[24.96, 0.5, 0.036, 1.56,],])

    
    nM = RetPar.shape[0]

    nHHydrus = 100
    hHydrus = -numpy.logspace(-6,5,nHHydrus)

    hHydrus_L = []
    Theta_L = []
    Con_L = []
    Cap_L = []
    for iM in range(nM):
        h = hHydrus
        hHydrus_L.append(h)
        Theta_L.append(Retention_model(h,RetPar[iM,:]))
        Con_L.append(Conductivity_model(h,ConPar[iM,:]))
        Cap_L.append(Capacity_model(h,RetPar[iM,:]))


    # REFERENCE K(h)
    postfix = '_REF'
    MaterFile_REF = CreateMaterIn(MaterRoot,hHydrus_L,Theta_L,Con_L,Cap_L,postfix =postfix)


    iH_S = numpy.argmin((hHydrus-hSens)**2)

    # PLUS K(h)
    Con_REF = Con_L[0]
    Con_PLUS = copy.deepcopy(Con_REF)
    if DoLog:
        Con_PLUS[iH_S] = 10**(numpy.log10(Con_PLUS[iH_S]) * (1.0 + Perturb))
    else:
        Con_PLUS[iH_S] = Con_PLUS[iH_S] * (1.0 + Perturb)
    Con_L = [Con_PLUS,]

    postfix = '_PLUS'
    MaterFile_PLUS = CreateMaterIn(MaterRoot,hHydrus_L,Theta_L,Con_L,Cap_L,postfix =postfix)


    # MIN K(h)
    Con_REF = Con_L[0]
    Con_MIN = copy.deepcopy(Con_REF)
    if DoLog:
        Con_MIN[iH_S] = 10**(numpy.log10(Con_MIN[iH_S]) * (1.0 - Perturb))
    else:
        Con_MIN[iH_S] = Con_MIN[iH_S] * (1.0 - Perturb)
    Con_L = [Con_MIN,]

    postfix = '_MIN'
    MaterFile_MIN = CreateMaterIn(MaterRoot,hHydrus_L,Theta_L,Con_L,Cap_L,postfix =postfix)


    dtMin = 1e-5

    #==============================================================================
    # REFERENCE RUN
    #==============================================================================

    print('Reference Run')

    ReplaceProjectPath = RefProjectPath+"_REF"
    ProjectDir = ReplaceProjectPath

    ## COPY TO TEMPORARY PROJECT
    # shutil.rmtree(ReplaceProjectPath)
    dir_util.copy_tree(RefProjectPath,ReplaceProjectPath)
    for OutFile in glob.glob(ReplaceProjectPath+'*.OUT'):
        os.remove(OutFile)
    for ErrFile in glob.glob(ReplaceProjectPath+'*.msg'):
        os.remove(ErrFile)


    ## SELECTOR

    SelectorFile = 'SELECTOR.IN'

    # Don't press enter
    ReplaceSelectorParameter(ProjectDir,SelectorFile,FindString='lEnter',NewValue='f')
    # ReplaceSelectorParameter(ProjectDir,SelectorFile,FindString='lPrintD',NewValue='f')
    ReplaceSelectorParameter(ProjectDir,SelectorFile,FindString='lScreen',NewValue='f')

    # Reduce Print
    # ReplaceSelectorParameter(ProjectDir,SelectorFile,FindString='tPrintInterval',NewValue='1')
    ReplaceSelectorParameter(ProjectDir,SelectorFile,FindString='nPrintSteps',NewValue='100')

    # Modify thr, ths and Ks
    ReplaceSelectorParameter(ProjectDir,SelectorFile,FindString='thr',NewValue='%1.3e'%RetPar[0][0])
    ReplaceSelectorParameter(ProjectDir,SelectorFile,FindString='ths',NewValue='%1.3e'%RetPar[0][1])
    ReplaceSelectorParameter(ProjectDir,SelectorFile,FindString='Ks',NewValue='%1.3e'%ConPar[0][0])
    ReplaceSelectorParameter(ProjectDir,SelectorFile,FindString='dtMin',NewValue='%1.3e'%dtMin)

    # ReplaceSelectorParameter(ProjectDir,SelectorFile,FindString='hTab1',NewValue='1e-3')
    # ReplaceSelectorParameter(ProjectDir,SelectorFile,FindString='hTabN',NewValue='1e5')

    ## PROFILE

    ProfileFile = 'PROFILE.DAT'


    ## Mater.in
    # COPY FILE CREATED BY CREATE COMPACTION PROFILE
    MaterIn_File = MaterFile_REF
    MaterFile = 'Mater.in'
    shutil.copyfile(MaterIn_File,os.path.join(ProjectDir,MaterFile))


    ## ATMOSPH
    AtmosphFile = 'ATMOSPH.IN'
    # ReplaceAtmosphParameter(ProjectDir,AtmosphFile,'hCritS','999') ## Ponding water
    
    ## HYDRUS 1D

    HydrusDatFile = 'HYDRUS1D.DAT'

    ## RUN ##

    if RunREF:
        runHydrusLinux(HydrusRoot,HydrusExec,ProjectDir)




    #==============================================================================
    # PLUS RUN
    #==============================================================================

    print('PLUS Run')

    ReplaceProjectPath = RefProjectPath+"_PLUS"
    ProjectDir = ReplaceProjectPath

    ## COPY TO TEMPORARY PROJECT
    # shutil.rmtree(ReplaceProjectPath)
    dir_util.copy_tree(RefProjectPath,ReplaceProjectPath)
    for OutFile in glob.glob(ReplaceProjectPath+'*.OUT'):
        os.remove(OutFile)
    for ErrFile in glob.glob(ReplaceProjectPath+'*.msg'):
        os.remove(ErrFile)


    ## SELECTOR

    SelectorFile = 'SELECTOR.IN'

    # Don't press enter
    ReplaceSelectorParameter(ProjectDir,SelectorFile,FindString='lEnter',NewValue='f')
    # ReplaceSelectorParameter(ProjectDir,SelectorFile,FindString='lPrintD',NewValue='f')
    ReplaceSelectorParameter(ProjectDir,SelectorFile,FindString='lScreen',NewValue='f')

    # Reduce Print
    # ReplaceSelectorParameter(ProjectDir,SelectorFile,FindString='tPrintInterval',NewValue='1000')
    ReplaceSelectorParameter(ProjectDir,SelectorFile,FindString='nPrintSteps',NewValue='100')

    # Modify thr, ths and Ks
    ReplaceSelectorParameter(ProjectDir,SelectorFile,FindString='thr',NewValue='%1.3e'%RetPar[0][0])
    ReplaceSelectorParameter(ProjectDir,SelectorFile,FindString='ths',NewValue='%1.3e'%RetPar[0][1])
    ReplaceSelectorParameter(ProjectDir,SelectorFile,FindString='Ks',NewValue='%1.3e'%ConPar[0][0])
    # ReplaceSelectorParameter(ProjectDir,SelectorFile,FindString='Ks',NewValue='%1.3e'%Con_PLUS[0])
    ReplaceSelectorParameter(ProjectDir,SelectorFile,FindString='dtMin',NewValue='%1.3e'%dtMin)
    
    ## PROFILE

    ProfileFile = 'PROFILE.DAT'


    ## Mater.in
    # COPY FILE CREATED BY CREATE COMPACTION PROFILE
    MaterIn_File = MaterFile_PLUS
    MaterFile = 'Mater.in'
    shutil.copyfile(MaterIn_File,os.path.join(ProjectDir,MaterFile))


    ## ATMOSPH
    AtmosphFile = 'ATMOSPH.IN'
    # ReplaceAtmosphParameter(ProjectDir,AtmosphFile,'hCritS','999') ## Ponding Water

    ## HYDRUS 1D

    HydrusDatFile = 'HYDRUS1D.DAT'

    ## RUN ##
    
    time.sleep(2.0)
    runHydrusLinux(HydrusRoot,HydrusExec,ProjectDir)




    #==============================================================================
    # MIN RUN
    #==============================================================================

    print('MIN Run')

    ReplaceProjectPath = RefProjectPath+"_MIN"
    ProjectDir = ReplaceProjectPath

    ## COPY TO TEMPORARY PROJECT
    # shutil.rmtree(ReplaceProjectPath)
    dir_util.copy_tree(RefProjectPath,ReplaceProjectPath)
    for OutFile in glob.glob(ReplaceProjectPath+'*.OUT'):
        os.remove(OutFile)
    for ErrFile in glob.glob(ReplaceProjectPath+'*.msg'):
        os.remove(ErrFile)


    ## SELECTOR

    SelectorFile = 'SELECTOR.IN'

    # Don't press enter
    ReplaceSelectorParameter(ProjectDir,SelectorFile,FindString='lEnter',NewValue='f')
    # ReplaceSelectorParameter(ProjectDir,SelectorFile,FindString='lPrintD',NewValue='f')
    ReplaceSelectorParameter(ProjectDir,SelectorFile,FindString='lScreen',NewValue='f')

    # Reduce Print
    # ReplaceSelectorParameter(ProjectDir,SelectorFile,FindString='tPrintInterval',NewValue='10')
    ReplaceSelectorParameter(ProjectDir,SelectorFile,FindString='nPrintSteps',NewValue='100')

    # Modify thr, ths and Ks
    ReplaceSelectorParameter(ProjectDir,SelectorFile,FindString='thr',NewValue='%1.3e'%RetPar[0][0])
    ReplaceSelectorParameter(ProjectDir,SelectorFile,FindString='ths',NewValue='%1.3e'%RetPar[0][1])
    ReplaceSelectorParameter(ProjectDir,SelectorFile,FindString='Ks',NewValue='%1.3e'%ConPar[0][0])
    # ReplaceSelectorParameter(ProjectDir,SelectorFile,FindString='Ks',NewValue='%1.3e'%Con_MIN[0])
    ReplaceSelectorParameter(ProjectDir,SelectorFile,FindString='dtMin',NewValue='%1.3e'%dtMin)

    ## PROFILE

    ProfileFile = 'PROFILE.DAT'


    ## Mater.in
    # COPY FILE CREATED BY CREATE COMPACTION PROFILE
    MaterIn_File = MaterFile_MIN
    MaterFile = 'Mater.in'
    shutil.copyfile(MaterIn_File,os.path.join(ProjectDir,MaterFile))


    ## ATMOSPH
    AtmosphFile = 'ATMOSPH.IN'
    # ReplaceAtmosphParameter(ProjectDir,AtmosphFile,'hCritS','999') ## Ponding Water

    ## HYDRUS 1D

    HydrusDatFile = 'HYDRUS1D.DAT'

    ## RUN ##
    time.sleep(2.0)
    runHydrusLinux(HydrusRoot,HydrusExec,ProjectDir)




    #==============================================================================
    # Collect Output
    #==============================================================================

    ## REF
    # ProjectDir = RefProjectPath+"_REF"
    # # ObsOutName = os.path.join(ProjectDir,'OBS_NODE.OUT')
    # # ObsNode_Nodes_REF,ObsNode_t_REF,ObsNode_L_REF,ObsNode_Header = ReadObsNodeOut(ObsOutName)

    # NodInfOutName = os.path.join(ProjectDir,'NOD_INF.OUT')
    # NodInf_T_REF,NodInf_Data_REF,NodInf_Header = ReadNodInfOut(NodInfOutName)
    # NodInf_Data_REF = numpy.array(NodInf_Data_REF)

    # TLevelOutName = os.path.join(ProjectDir,'T_LEVEL.OUT')
    # TLevel_Data_REF,TLevel_Header = ReadTLevelOut(TLevelOutName)

    ## PLUS
    ProjectDir = RefProjectPath+"_PLUS"
    # ObsOutName = os.path.join(ProjectDir,'OBS_NODE.OUT')
    # ObsNode_Nodes_PLUS,ObsNode_t_PLUS,ObsNode_L_PLUS,ObsNode_Header = ReadObsNodeOut(ObsOutName)

    NodInfOutName = os.path.join(ProjectDir,'NOD_INF.OUT')
    NodInf_T_PLUS,NodInf_Data_PLUS,NodInf_Header = ReadNodInfOut(NodInfOutName)
    NodInf_Data_PLUS = numpy.array(NodInf_Data_PLUS)

    TLevelOutName = os.path.join(ProjectDir,'T_LEVEL.OUT')
    TLevel_Data_PLUS,TLevel_Header = ReadTLevelOut(TLevelOutName)

    ## MIN
    ProjectDir = RefProjectPath+"_MIN"
    # ObsOutName = os.path.join(ProjectDir,'OBS_NODE.OUT')
    # ObsNode_Nodes_MIN,ObsNode_t_MIN,ObsNode_L_MIN,ObsNode_Header = ReadObsNodeOut(ObsOutName)

    NodInfOutName = os.path.join(ProjectDir,'NOD_INF.OUT')
    NodInf_T_MIN,NodInf_Data_MIN,NodInf_Header = ReadNodInfOut(NodInfOutName)
    NodInf_Data_MIN = numpy.array(NodInf_Data_MIN)

    TLevelOutName = os.path.join(ProjectDir,'T_LEVEL.OUT')
    TLevel_Data_MIN,TLevel_Header = ReadTLevelOut(TLevelOutName)



    #==============================================================================
    # SENSITIVITY
    #==============================================================================

    NodInfSensCol = [2,3,6]
    vmaxL=[7,3,3]
    nC = len(NodInfSensCol)

    CAS = []
    CPRS = []
    CTRS = []
    
    for iC in range(nC):
        iNI = NodInfSensCol[iC]
        
        df_par_plus = NodInf_Data_PLUS[:,:,iNI]
        df_par_min = NodInf_Data_MIN[:,:,iNI]

        average_out = (df_par_plus+df_par_min)/2.
        #sensitivity indices:
        # CAS = (df_par_plus-df_par_min)/(2.*perturbation_factor*parameter_value) #dy/dp
        if DoLog:
            CAS_T = (df_par_plus-df_par_min)/(numpy.log10(Con_PLUS[iH_S]) - numpy.log10(Con_MIN[iH_S])) #dy/dp
        else:
            CAS_T = (df_par_plus-df_par_min)/(Con_PLUS[iH_S] - Con_MIN[iH_S]) #dy/dp
        CPRS_T = CAS_T*Con_REF[iH_S]
        CTRS_T = CAS_T*Con_REF[iH_S]/average_out 
        
        CAS.append(CAS_T)
        CPRS.append(CPRS_T)
        CTRS.append(CTRS_T)

    [CAS,CPRS,CTRS] = map(numpy.array,[CAS,CPRS,CTRS])

    return CAS,CPRS,CTRS


#==============================================================================
# Test Sensitivity for full range
#==============================================================================

hSens_All = -numpy.logspace(-1,5,50)
nH = hSens_All.size
Perturb = 5e-2

Case = 'LOAM'
RetPar = numpy.array([[0.078, 0.43, 0.036, 1.56,],])
ConPar = numpy.array([[24.96, 0.50, 0.036, 1.56,],])
GuessRuntime = 320

CAS_All = []
CPRS_All = []
CTRS_All = []
for iH in range(nH):
    print("%s / %s"%(iH,nH))
    CAS,CPRS,CTRS = TestSensitivity(hSens_All[iH],Perturb,
                                        RetPar, ConPar,
                                        GuessRuntime = GuessRuntime,
                                        RunREF = iH == nH-1 )
    CAS_All.append(CAS)
    CPRS_All.append(CPRS)
    CTRS_All.append(CTRS)
    # figlist.append(fig)

CAS_All = numpy.array(CAS_All)
CPRS_All = numpy.array(CPRS_All)
CTRS_All = numpy.array(CTRS_All)

nH,nC,nT,nZ = CAS_All.shape


## PLOT

RefProjectPath = "/home/jan/git/KsIsIrrelevant/HydrusProjects/TestSens_MeteoSimple2"
ProjectDir = RefProjectPath+"_REF"
NodInfOutName = os.path.join(ProjectDir,'NOD_INF.OUT')
NodInf_T_REF,NodInf_Data_REF,NodInf_Header = ReadNodInfOut(NodInfOutName)
NodInf_Data_REF = numpy.array(NodInf_Data_REF)
TMesh,ZMesh = numpy.meshgrid(NodInf_T_REF,NodInf_Data_REF[0][:,1])

## !!!!!!!!! dZ WEIGHT !!!!!!!!!!
Z_REF = NodInf_Data_REF[0][:,1]
dZ_REF = numpy.zeros_like(Z_REF)
dZ_REF0 = Z_REF[1:] - Z_REF[:-1]
dZ_REF[0] = dZ_REF0[0]/2
dZ_REF[1:-1] = (dZ_REF0[1:]+dZ_REF0[:-1])/2
dZ_REF[-1] = dZ_REF0[-1]/2

CAS_All_Z = CAS_All*dZ_REF[None,None,None,:]
CPRS_All_Z = CPRS_All*dZ_REF[None,None,None,:]
CTRS_All_Z = CTRS_All*dZ_REF[None,None,None,:]


fig =pylab.figure()
ax = fig.add_subplot(111)
V_T = numpy.logspace(-3,12,15)
V = numpy.concatenate((-V_T[::-1],V_T))
cmap = pylab.get_cmap('bwr')
cols = cmap(numpy.linspace(0,1,V.size))
im = ax.contourf(TMesh,ZMesh,CAS_All[0,0,:,:].transpose(),V, colors= cols)
# ax.contour(TMesh,ZMesh,CAS_All[10,0,:,:].transpose(),V,colors='k')
ax.set_xlabel('Time (d)')
ax.set_ylabel('Z (cm)')
# ax.set_title('%s - h = %1.2e cm'%(NodInf_Header[iNI],hSens))
# fig.colorbar(im)

for iC in range(nC):
    # fig = pylab.figure()
    # ax =fig.add_subplot(111)
    # ax.plot(-hSens_All,numpy.nanmean(numpy.abs(CAS_All_Z[:,iC,:,:]).reshape((nH,-1))))
    # ax.set_xscale('log')
    # ax.set_xlabel('Matric Head (cm)')
    # ax.set_ylabel('CAS')
    # ax.set_yscale('log')
    # figlist.append(fig)

    # fig = pylab.figure()
    # ax =fig.add_subplot(111)
    # ax.plot(-hSens_All,numpy.nanmean(numpy.abs(CPRS_All_Z[:,iC,:,:]).reshape((nH,-1))))
    # ax.set_xscale('log')
    # ax.set_xlabel('Matric Head (cm)')
    # ax.set_ylabel('CPRS')
    # ax.set_yscale('log')
    # figlist.append(fig)

    fig = pylab.figure()
    ax =fig.add_subplot(111)
    ax.plot(-hSens_All,numpy.nanmean(numpy.abs(CTRS_All_Z[:,iC,:,:]).reshape((nH,-1)),axis=-1))
    ax.set_xscale('log')
    ax.set_xlabel('Matric Head (cm)')
    ax.set_ylabel('CTRS')
    # ax.set_yscale('log')
    figlist.append(fig)


#==============================================================================
# Save
#==============================================================================

## save plots to pdf
Root = ''

import sys
basename=os.path.basename(sys.argv[0])[:-3]

postfix='_%s'%Case
# postfix=''

pdfname = '%s%s.pdf'%(basename,postfix)
pp = PdfPages(pdfname)
for fig in figlist:
    pp.savefig(fig)
pp.close()

pylab.show()