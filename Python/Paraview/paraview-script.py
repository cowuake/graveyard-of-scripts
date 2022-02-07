from paraview.simple import *
from os import system as sys

# ================
# GLOBAL VARIABLES
# ================
runs = list(range(100, 10000, 100))
prediffuserClipOrigin = [0.05, 0.0, 0.0]

for run in runs:

    # =============
    # PREPARE FILES
    # =============

    sys("cp hydra.flow.hdf_{0} hydra.flow.hdf".format(run))

    # ============
    # DATA LOADING
    # ============

    # Load the case file with the Hydra Grid Reader
    case = HydraGridReader(InputFile = "input.hyd")

    # =============================
    # EXTRACTION OF TARGET SURFACES
    # =============================

    # OGV Inlet
    ogvInletExtract = ExtractBlockByRule(Input=case)
    ogvInletExtract.FilterRule = "Surface Name=INLET && Zone Name=OGV"

    # Pre-diffuser Exit
    ghost = ExtractBlockByRule(Input=case)
    ghost.FilterRule = "Object Type=Ghost Surface"
    prediffuserExitExtract = Clip(Input=ghost)
    prediffuserExitExtract.ClipType.Origin = prediffuserClipOrigin

    # Cage Exit
    cageExitExtract = ExtractBlockByRule(Input=case)
    cageExitExtract.FilterRule = "Surface Name=EXIT && Zone Name=OGV"

    # Inner Annulus Exit
    innerAnnulusExitExtract = ExtractBlockByRule(Input=case)
    innerAnnulusExitExtract.FilterRule = "Surface Name=Exit1 && Zone Name=OGV"

    # Outer Annulus Exit
    outerAnnulusExitExtract = ExtractBlockByRule(Input=case)
    outerAnnulusExitExtract.FilterRule = "Surface Name=Exit2 && Zone Name=OGV"

    # ===============================
    # COMPUTE ABSOLUTE TOTAL PRESSURE
    # ===============================

    #ogvInletP0 = RRCalculator(Input=ogvInletExtract)
    #ogvInletP0.PointDataParameter="absolute total pressure"

    # =====================================================
    # CONVERT POINT DATA TO CELL DATA AND INTEGRAL AVERAGES
    # =====================================================
    # This is needed in order to access Area when integrating
    # The "Divide by surface" option, which allows for integral averages, is not
    # available for all Paraview versions (old releases suck)

    #ogvInletToCell = PointDatatoCellData(Input=ogvInletP0)

    # ===================
    # INTEGRATE VARIABLES
    # ===================

    #ogvInletIntegral = IntegrateVariables(Input=ogvInletToCell)

    # =================================
    # RETRIEVE MASS-AVERAGED QUANTITIES
    # =================================

    # OGV Inlet
    ogvInletReduce = ReduceObject(Input=ogvInletExtract)
    ogvInletReduce.PointDataParameter = ["MM static pressure", "MM absolute total pressure"]
    ogvInletData = ogvInletReduce.GetPointDataInformation()
    ogvInletStaticPressure = ogvInletData.GetArray("MM static pressure").GetRange()[0]
    ogvInletStagnationPressure = ogvInletData.GetArray("MM absolute total pressure").GetRange()[0]

    # Pre-diffuser Exit
    prediffuserExitReduce = ReduceObject(Input=prediffuserExitExtract)
    prediffuserExitReduce.PointDataParameter = ["MM static pressure", "MM absolute total pressure"]
    prediffuserExitData = prediffuserExitReduce.GetPointDataInformation()
    prediffuserExitStaticPressure = prediffuserExitData.GetArray("MM static pressure").GetRange()[0]
    prediffuserExitStagnationPressure = prediffuserExitData.GetArray("MM absolute total pressure").GetRange()[0]

    # Cage Exit
    cageExitReduce = ReduceObject(Input=cageExitExtract)
    cageExitReduce.PointDataParameter = ["MM static pressure", "MM absolute total pressure"]
    cageExitData = cageExitReduce.GetPointDataInformation()
    cageExitStaticPressure = cageExitData.GetArray("MM static pressure").GetRange()[0]
    cageExitStagnationPressure = cageExitData.GetArray("MM absolute total pressure").GetRange()[0]

    # Inner Annulus Exit
    innerAnnulusExitReduce = ReduceObject(Input=innerAnnulusExitExtract)
    innerAnnulusExitReduce.PointDataParameter = ["MM static pressure", "MM absolute total pressure"]
    innerAnnulusExitData = innerAnnulusExitReduce.GetPointDataInformation()
    innerAnnulusExitStaticPressure = innerAnnulusExitData.GetArray("MM static pressure").GetRange()[0]
    innerAnnulusExitStagnationPressure = innerAnnulusExitData.GetArray("MM absolute total pressure").GetRange()[0]

    # Outer Annulus Exit
    outerAnnulusExitReduce = ReduceObject(Input=outerAnnulusExitExtract)
    outerAnnulusExitReduce.PointDataParameter = ["MM static pressure", "MM absolute total pressure"]
    outerAnnulusExitData = outerAnnulusExitReduce.GetPointDataInformation()
    outerAnnulusExitStaticPressure = outerAnnulusExitData.GetArray("MM static pressure").GetRange()[0]
    outerAnnulusExitStagnationPressure = outerAnnulusExitData.GetArray("MM absolute total pressure").GetRange()[0]

    # ================================================
    # COMPUTE PRESSURE RECOVERY COEFFICIENT AND LOSSES
    # ================================================

    def computeCp(p0, P0, p1, P1):
        return (p1 - p0) / (P0 - p0)

    def computeLambda(p0, P0, p1, P1):
        return (P1 - P0) / (P0 - p0)

    p0 = ogvInletStaticPressure
    P0 = ogvInletStagnationPressure

    #p0_2 = prediffuserExitStaticPressure
    #P0_2 = prediffuserExitStagnationPressure

    # Pre-diffuser Exit
    prediffuserExitCp = (prediffuserExitStaticPressure - p0) / (P0 - p0)
    prediffuserExitLambda = (prediffuserExitStagnationPressure - P0) / (P0 - p0)

    # Cage Exit
    cageExitCp = (cageExitStaticPressure - p0) / (P0 - p0)
    cageExitLambda = (cageExitStagnationPressure - P0) / (P0 - p0)

    # Inner Annulus Exit
    innerAnnulusExitCp = (innerAnnulusExitStaticPressure - p0) / (P0 - p0)
    innerAnnulusExitLambda = (innerAnnulusExitStagnationPressure - P0) / (P0 - p0)

    # Outer Annulus Exit
    outerAnnulusExitCp = (outerAnnulusExitStaticPressure - p0) / (P0 - p0)
    outerAnnulusExitLambda = (outerAnnulusExitStagnationPressure - P0) / (P0 - p0)

    with open("Cp.txt", "a") as f:
        f.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format(run,
                                            prediffuserExitCp,
                                            cageExitCp,
                                            innerAnnulusExitCp,
                                            outerAnnulusExitCp))

    with open("lambda.txt", "a") as f:
        f.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format(run,
                                            prediffuserExitLambda,
                                            cageExitLambda,
                                            innerAnnulusExitLambda,
                                            outerAnnulusExitLambda))

    # ======================
    # RESET PARAVIEW SESSION
    # ======================

    Disconnect()
    Connect()
