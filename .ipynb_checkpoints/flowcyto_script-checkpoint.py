import FlowCytometryTools as Flow
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import path
from scipy import stats
import numpy as np

sns.set()
HEKgate_vert = [(65586.908970624077, 209889.25615062829), (106230.06134133368, 220573.80805632132),(183452.05084568192,221202.31110959739),
          (251964.79341344949, 160237.51494181945), (256029.10865052044, 87016.909235158339),
          (218289.03859200439, 48049.719932042557), (150937.5289491142, 43650.198559110126), 
          (88811.567468172419, 62505.290157391952), (58038.894958920879,103357.98862033593), (35394.852923811239, 174064.5821138928),
          (35394.852923811239, 174378.83364053082)]

COSgate_vert = [(74854.965334995941, 199780.96144326951), (56853.759014296673, 155119.82832315788),
              (35368.448244429826, 95676.489170333225), (48143.497891377679,59507.261643482241),
              (102727.80092833671, 40007.33028117998), (175313.31028599496, 43152.480500906153), 
                (236865.82222128916, 87184.583577072553), (236865.82222128916, 140652.13731241747), (229316.9292480927, 191918.08589395409),
                (193314.5166066942,  217708.31769570868), (84726.594607637468, 205442.23183877661), (84726.594607637468, 205442.23183877661)]
CHOgate_vert=[(49116.89459105306, 188512.82171075896), (93064.823285338003,230312.4363277517), (139636.50891659516, 241312.33491117085),
              (184896.31607936622, 213026.88141095021), (205886.37157514412, 118113.47077687646), (189815.86033618916,54628.341809714555),
              (113726.90916399437, 37971.352526251285), (51084.712293782242, 54314.058993045437), (23863.234072695304, 100513.63304340583),
              (35998.109906191894, 162741.63074389126), (35998.109906191894, 162741.63074389126)]
gates = {'HEK293T':HEKgate_vert,'COS7':COSgate_vert,'CHOLEC2':CHOgate_vert}


class flowData:
    def __init__(self,sampleDataEntry,flowDataLocation):
        self.Cells = sampleDataEntry['cells']
        self.MOI = sampleDataEntry['MOI']
        self.Virus = sampleDataEntry['Virus']
        self.flowData = Flow.FCMeasurement(ID=sampleDataEntry['cells'], datafile=flowDataLocation)
        self.Data = self.flowData.data
    
    def cGate(self,vertices='test'):
        """Gating the flow data using forward and side scattering data
        Inputs:
        vertices        A list of tuples of the coordinates (fsc, ssc), generate using cytoflow bin (or .exe)
        default vertices are HEKgate_vert and COSgate_vert
        """
        if vertices=='test':
            vertices = gates[self.Cells]
        v = path.Path(vertices)
        fsc_ssc = pd.DataFrame([self.Data['FSC-A'],self.Data['SSC-A']]).transpose()

        gatedPoints = v.contains_points(fsc_ssc)
        self.gatedPopulation1 = self.Data[gatedPoints]
        return self.gatedPopulation1
    
    def fGate(self,threshhold=100):
        """Gating the FITC data
        Inputs:
        threshhold      A value threshhold for the fluorescent default value is 100
        """
        if hasattr(self, 'gatedPopulation1'):
            t = self.gatedPopulation1['FITC-A']>threshhold
            self.greenPop = self.gatedPopulation1[t]
            self.geomean = stats.gmean(self.greenPop['FITC-A'])
            self.percentPos = sum(np.array(t))/float(len(self.gatedPopulation1))
            return (self.geomean, self.percentPos)
        else:
            print('run cellGate first')
    
    def cGateplot(self, c='b'):
        sns.jointplot(x=self.Data['FSC-A'],y=self.Data['SSC-A'],kind='kde',color=c)
    def fGateplot(self):    
        sns.distplot(x=self.Data['FITC-A'],bin=np.logspace(0,5,100))




def dataCollector(sampleDataEntry,filePaths_fcs):
    """auto runner for the data collections
    inputs:
    sampleDataEntry |   panda dataframe from pd.read_cvs of the sampleNames.csv
    filePaths_fcs   |   a list of strings with the paths to the fcs files.

    outputs:
    tuple of dataframe with general analysis and dictionary of flowData objects
    """
    dataCollection = {}
    geoMeans = []
    percPositive = []
    TI = []
    date=[]
    fileName_0= ''
    for f in range(len(filePaths_fcs)):
        file = filePaths_fcs[f]
        fD_obj = flowData(sampleDataEntry.iloc[f,:],file)
        fD_obj.cGate() # cell gates
        fD_obj.fGate() # flouresence gate
        fileName_N = fD_obj.Cells + '_MOI_' + str(fD_obj.MOI) + '_Virus_' + fD_obj.Virus
        if fileName_0 == fileName_N:
            i+=1
            fileName = fileName_N + '_' + str(i)
        else:
            fileName = fileName_N + '_1'
            fileName_0 = fileName_N
            i = 1
        dataCollection[fileName] = fD_obj
        geoMeans.append(fD_obj.geomean)
        percPositive.append(fD_obj.percentPos*100)
        date.append(fD_obj.flowData.meta['$DATE'])
        TI.append(fD_obj.percentPos*fD_obj.geomean)
    transductionMetrics = pd.DataFrame([geoMeans,percPositive,TI,date]).transpose()
    transductionMetrics.columns=['geoMean','Percent +ve (%)','TI','flow_date']
    df1_Out = sampleDataEntry.merge(transductionMetrics,left_index=True, right_index=True)
    return (df1_Out, dataCollection)