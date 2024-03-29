from yggdrasil import runner
import argparse

parser = argparse.ArgumentParser(description='Run an integration.')
parser.add_argument('yamlfile', nargs='+',help='One or more yaml specification files.')
args1 = parser.parse_args(["../ePhotosynthesis/EphotosynthesisOnly.yml"])
args2 = parser.parse_args(["../FBA/yggrasil_ODE_FBA_dielFBA.yaml"])

PPFD_Dict = dict()
PPFD_Dict[164]={372:[318.055556,602.777778,951.449275,1342.753623,
                     1621.428572,1637.301587,1640.909091,1595.454545,
                     1550,1534.848485,1519.69697,1383.333334,1216.666667,
                     847,390.333333],
                552:[325,616.666667,931.15942,1260.869565,1529.761905,
                     1656.746032,1756.818182,1749.242424,1741.666667,
                     1688.636364,1635.606061,1431.349206,1189.285715,
                     811,374.333333]}
PPFD_Dict[176]={372:[158.442982,331.688597,731.907895,1076.704545,
                     1199.810606,1322.916667,1542.214912,1727.52193,
                     1606.907895,1455.492425,1180.871212,906.25,
                     516.995614,157.986111],
                552:[158.442982,328.947368,701.754386,1055.871212,
                     1335.227273,1614.583333,1713.267544,1795.504386,
                     1729.714912,1603.219697,1233.901515,864.583333,
                     519.188596,197.337963]}
PPFD_Dict[191]={372:[295.019894,583.116711,935.42152,1330.531657,
                     1611.134555,1624.553894,1631.636984,1581.693823,
                     1531.750663,1516.956836,1502.163009,1370.754763,
                     1210.192911,846.002768,346.060431],
                552:[307.950928,608.97878,923.277592,1246.423711,
                     1509.096938,1630.666747,1739.026147,1728.491853,
                     1717.95756,1670.248372,1622.539185,1424.359778,
                     1188.562817,807.771883,330.318302]}
PPFD_Dict[205]={372:[112.980769,205.128205,372.596154,653.044872,933.49359,
                     1075.79023,1183.54885,1291.307471,1187.026515,
                     1059.185606,789.402174,458.786232,193.75],
                552:[129.00641,237.179487,425.480769,733.974359,1042.467949,
                     1127.155173,1155.890805,1184.626437,1242.897727,
                     1304.450758,999.547102,537.59058,193.75]}
PPFD_Dict[215]={372:[159.172662,280.57554,361.510791,442.446043,760.507384,
                     1050.359712,1086.330935,1122.302158,1066.18705,795.203837,
                     524.220623,302.506382,195.753075],
                552:[159.172662,275.779376,344.724221,413.669065,803.672851,
                     1155.190134,1160.32888,1165.467626,1087.769784,816.786571,
                     545.803357,322.000464,208.284985]}
PPFD_Dict[233]={372:[232.613909,454.436451,602.517986,701.438849,762.433531,
                     671.723491,581.01345,678.586543,729.701953,511.305242,
                     292.90853,133.093525],
                552:[289.568345,568.345324,715.827338,775.779377,799.343134,
                     677.353769,555.364404,680.279306,756.423433,527.749229,
                     299.075026,133.093525]}
PPFD_Dict[254]={372:[337.985726,666.693101,994.752267,1322.649381,1528.183981,
                     1611.356066,1694.528152,1538.11412,1381.700087,1287.141725,
                     1199.456214,758.849398,230.012256],
                552:[291.125369,572.972388,941.861187,1332.51043,1561.173672,
                     1627.850912,1694.528152,1616.661002,1538.793852,1381.397983,
                     1215.165591,758.849398,230.012256]}



import pandas as pd



for day in  PPFD_Dict.keys():

    df = pd.DataFrame()
    for co2 in PPFD_Dict[day].keys():

        PPFD_list = list()
        Vc_list = list()
        Vo_list = list()
        A_list = list()
        St_list = list()
        j=0
        for PPFD in PPFD_Dict[day][co2]:
            PPFD = round(PPFD,3)
            print("##############################")
            print("Day = "+str(day))
            j=j+1
            print("Hour = "+str(j))

            F_weather = open("../ePhotosynthesis/InputEvn.txt","w")
            F_weather.write("CO2 "+str(co2)+"\nPPFD "+str(PPFD)+"\nSucPath 1"+"\ndaylength "+str(len(PPFD_Dict[day][co2])))
            F_weather.close()


            f1 = open("../ePhotosynthesis/InputATPCost.txt","w")
            f1.write("ATPCost 0")
            f1.close()
            
            f3 = open("../ePhotosynthesis/InputNADPHCost.txt","w")
            f3.write("NADPHCost 0")
            f3.close()
            #print(PPFD)

            runner.run(args1.yamlfile)
            runner.run(args2.yamlfile)
            import os
            os.rename("./../FBA/Diel_flux.csv","./../Validations/Diel_flux_FACE"+str(day)+"_"+str(co2)+"_"+str(j)+".csv")
            os.rename("./../ePhotosynthesis/OutputFluxT.txt","./../Validations/ODE_only_OutputFluxT_FACE"+str(day)+"_"+str(co2)+"_"+str(j)+".csv")
