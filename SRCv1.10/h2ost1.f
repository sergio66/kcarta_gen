C***********************************************************************
C
C  PROGRAM        H2OST1  BLOCK DATA
C
C  PURPOSE        SELF DENSITY DEPENDENT CONTINUUM VALUES FOR H2O 
C                 AT 260K
C
C  VERSION        3.Y   D.P. EDWARDS   02/02/93 
C
C  DESCRIPTION    VALUES TAKEN FROM FASCOD1B
C                 CLOUGH ET AL. (1981) SPIE VOL. 277 ATMOS. TRANS.
C                 THIS ROUTINE IS A BASED ON BLOCK DATA BS260 6/28/82. 
C                 THE FASCODE CONTINUUM COMPONENT DUE TO THE LINE
C                 CENTER HAS BEEN REMOVED FOR USE WITH GENLN2 
C                 I.E. THE VALUES HERE REPRESENT THE LINE WINGS GREATER
C                 THAN 25 CM-1 FROM LINE CENTER
C
C***********************************************************************
C
       BLOCK DATA H2OST1
C-----------------------------------------------------------------------
       COMMON /CH2OS1/ 
     + S0001(50),S0051(50),S0101(50),S0151(50),S0201(50),S0251(50),
     + S0301(50),S0351(50),S0401(50),S0451(50),S0501(50),S0551(50),
     + S0601(50),S0651(50),S0701(50),S0751(50),S0801(50),S0851(50),
     + S0901(50),S0951(50),S1001(50),S1051(50),S1101(50),S1151(50),
     + S1201(50),S1251(50),S1301(50),S1351(50),S1401(50),S1451(50),
     + S1501(50),S1551(50),S1601(50),S1651(50),S1701(50),S1751(50),
     + S1801(50),S1851(50),S1901(50),S1951(50),S2001(1)            
C-----------------------------------------------------------------------
C
C       DATA VB1,VT1,DV1,NPT1 / 0.0, 20000.0, 10.0, 2001/         
C
       DATA S0001/
     + 1.6457E-01, 1.7045E-01, 1.7750E-01, 2.0036E-01, 2.1347E-01,
     + 2.2454E-01, 2.3428E-01, 2.3399E-01, 2.3022E-01, 2.0724E-01,
     + 1.9712E-01, 1.8317E-01, 1.6724E-01, 1.4780E-01, 1.2757E-01,
     + 1.1626E-01, 1.0098E-01, 8.9033E-02, 7.9770E-02, 6.7416E-02,
     + 5.9588E-02, 5.1117E-02, 4.6218E-02, 4.2179E-02, 3.4372E-02,
     + 2.9863E-02, 2.5252E-02, 2.2075E-02, 1.9209E-02, 1.5816E-02,
     + 1.3932E-02, 1.1943E-02, 1.0079E-02, 8.7667E-03, 7.4094E-03,
     + 6.4967E-03, 5.5711E-03, 4.8444E-03, 4.2547E-03, 3.6947E-03,
     + 3.2762E-03, 2.9057E-03, 2.6026E-03, 2.3296E-03, 2.1009E-03,
     + 1.8964E-03, 1.7102E-03, 1.5537E-03, 1.4171E-03, 1.2910E-03/
       DATA S0051/
     + 1.1788E-03, 1.0761E-03, 9.8498E-04, 9.0388E-04, 8.2990E-04,
     + 7.6261E-04, 7.0101E-04, 6.4582E-04, 5.9526E-04, 5.4948E-04,
     + 5.0739E-04, 4.6858E-04, 4.3304E-04, 4.0025E-04, 3.7051E-04,
     + 3.4315E-04, 3.1802E-04, 2.9496E-04, 2.7377E-04, 2.5425E-04,
     + 2.3626E-04, 2.1973E-04, 2.0449E-04, 1.9040E-04, 1.7739E-04,
     + 1.6543E-04, 1.5441E-04, 1.4424E-04, 1.3485E-04, 1.2617E-04,
     + 1.1816E-04, 1.1076E-04, 1.0391E-04, 9.7561E-05, 9.1695E-05,
     + 8.6270E-05, 8.1251E-05, 7.6605E-05, 7.2301E-05, 6.8310E-05,
     + 6.4612E-05, 6.1182E-05, 5.8001E-05, 5.5047E-05, 5.2307E-05,
     + 4.9761E-05, 4.7395E-05, 4.5197E-05, 4.3155E-05, 4.1256E-05/
       DATA S0101/
     + 3.9491E-05, 3.7849E-05, 3.6324E-05, 3.4908E-05, 3.3594E-05,
     + 3.2373E-05, 3.1243E-05, 3.0200E-05, 2.9239E-05, 2.8354E-05,
     + 2.7545E-05, 2.6812E-05, 2.6143E-05, 2.5546E-05, 2.5025E-05,
     + 2.4573E-05, 2.4194E-05, 2.3881E-05, 2.3652E-05, 2.3515E-05,
     + 2.3469E-05, 2.3498E-05, 2.3673E-05, 2.4012E-05, 2.4541E-05,
     + 2.5157E-05, 2.5955E-05, 2.7035E-05, 2.8544E-05, 3.0732E-05,
     + 3.3058E-05, 3.6428E-05, 4.0257E-05, 4.5378E-05, 5.2484E-05,
     + 6.0753E-05, 7.0670E-05, 8.3753E-05, 9.5236E-05, 1.1163E-04,
     + 1.2821E-04, 1.4856E-04, 1.6892E-04, 2.0097E-04, 2.2399E-04,
     + 2.4943E-04, 2.7372E-04, 3.1270E-04, 3.7028E-04, 3.9725E-04/
       DATA S0151/
     + 4.3051E-04, 4.6893E-04, 5.1236E-04, 5.5023E-04, 5.4643E-04,
     + 5.3951E-04, 4.9955E-04, 4.5912E-04, 4.2109E-04, 3.9613E-04,
     + 3.6749E-04, 3.6212E-04, 3.8466E-04, 4.2112E-04, 4.4392E-04,
     + 4.6668E-04, 4.9021E-04, 5.0953E-04, 5.0362E-04, 4.5866E-04,
     + 4.2906E-04, 4.1850E-04, 3.7367E-04, 3.2746E-04, 2.8387E-04,
     + 2.5128E-04, 2.2886E-04, 1.9586E-04, 1.6870E-04, 1.4874E-04,
     + 1.3558E-04, 1.1513E-04, 1.0400E-04, 8.7573E-05, 7.8521E-05,
     + 7.1463E-05, 6.3347E-05, 5.7147E-05, 4.6994E-05, 4.2307E-05,
     + 3.6997E-05, 3.2196E-05, 2.8837E-05, 2.5366E-05, 2.2885E-05,
     + 2.0005E-05, 1.7256E-05, 1.5511E-05, 1.3397E-05, 1.1993E-05/
       DATA S0201/
     + 1.0681E-05, 9.6464E-06, 8.7623E-06, 7.7769E-06, 7.0781E-06,
     + 6.3960E-06, 5.8137E-06, 5.3366E-06, 4.8912E-06, 4.5275E-06,
     + 4.1647E-06, 3.8587E-06, 3.5854E-06, 3.3367E-06, 3.1164E-06,
     + 2.9147E-06, 2.7322E-06, 2.5629E-06, 2.4077E-06, 2.2657E-06,
     + 2.1353E-06, 2.0154E-06, 1.9046E-06, 1.8026E-06, 1.7070E-06,
     + 1.6181E-06, 1.5353E-06, 1.4581E-06, 1.3859E-06, 1.3178E-06,
     + 1.2544E-06, 1.1950E-06, 1.1394E-06, 1.0873E-06, 1.0384E-06,
     + 9.9248E-07, 9.4933E-07, 9.0872E-07, 8.7049E-07, 8.3445E-07,
     + 8.0045E-07, 7.6833E-07, 7.3799E-07, 7.0930E-07, 6.8216E-07,
     + 6.5647E-07, 6.3213E-07, 6.0908E-07, 5.8724E-07, 5.6654E-07/
       DATA S0251/
     + 5.4692E-07, 5.2834E-07, 5.1075E-07, 4.9413E-07, 4.7849E-07,
     + 4.6375E-07, 4.4998E-07, 4.3713E-07, 4.2532E-07, 4.1420E-07,
     + 4.0412E-07, 3.9483E-07, 3.8586E-07, 3.7766E-07, 3.6935E-07,
     + 3.6149E-07, 3.5295E-07, 3.4507E-07, 3.3737E-07, 3.3128E-07,
     + 3.2515E-07, 3.1738E-07, 3.1082E-07, 3.0539E-07, 3.0119E-07,
     + 2.9662E-07, 2.9085E-07, 2.8676E-07, 2.8365E-07, 2.8059E-07,
     + 2.7721E-07, 2.7386E-07, 2.7048E-07, 2.6766E-07, 2.6581E-07,
     + 2.6500E-07, 2.6575E-07, 2.6867E-07, 2.7325E-07, 2.8046E-07,
     + 2.9192E-07, 3.0952E-07, 3.2719E-07, 3.5328E-07, 3.9164E-07,
     + 4.4220E-07, 5.0090E-07, 5.6728E-07, 6.6734E-07, 7.9746E-07/
       DATA S0301/
     + 9.3055E-07, 1.0756E-06, 1.1686E-06, 1.3217E-06, 1.4814E-06,
     + 1.5627E-06, 1.6519E-06, 1.7601E-06, 1.9060E-06, 2.0474E-06,
     + 2.0716E-06, 2.0433E-06, 1.9752E-06, 1.8466E-06, 1.7526E-06,
     + 1.6657E-06, 1.5870E-06, 1.5633E-06, 1.6520E-06, 1.8471E-06,
     + 1.9953E-06, 2.0975E-06, 2.2016E-06, 2.2542E-06, 2.3081E-06,
     + 2.3209E-06, 2.2998E-06, 2.3056E-06, 2.2757E-06, 2.2685E-06,
     + 2.2779E-06, 2.2348E-06, 2.2445E-06, 2.3174E-06, 2.4284E-06,
     + 2.5290E-06, 2.7340E-06, 2.9720E-06, 3.2332E-06, 3.5392E-06,
     + 3.9013E-06, 4.3334E-06, 4.9088E-06, 5.3428E-06, 5.9142E-06,
     + 6.6106E-06, 7.4709E-06, 8.5019E-06, 9.6835E-06, 1.0984E-05/
       DATA S0351/
     + 1.2831E-05, 1.4664E-05, 1.7080E-05, 2.0103E-05, 2.4148E-05,
     + 2.7948E-05, 3.2855E-05, 3.9046E-05, 4.6429E-05, 5.6633E-05,
     + 6.6305E-05, 7.6048E-05, 8.7398E-05, 1.0034E-04, 1.1169E-04,
     + 1.2813E-04, 1.3354E-04, 1.3952E-04, 1.4204E-04, 1.4615E-04,
     + 1.5144E-04, 1.5475E-04, 1.6561E-04, 1.7135E-04, 1.6831E-04,
     + 1.6429E-04, 1.6353E-04, 1.6543E-04, 1.5944E-04, 1.5404E-04,
     + 1.5458E-04, 1.6287E-04, 1.7277E-04, 1.8387E-04, 1.7622E-04,
     + 1.6360E-04, 1.5273E-04, 1.3667E-04, 1.2364E-04, 9.7576E-05,
     + 7.9140E-05, 6.4241E-05, 5.1826E-05, 4.1415E-05, 3.1347E-05,
     + 2.5125E-05, 2.0027E-05, 1.6362E-05, 1.3364E-05, 1.1117E-05/
       DATA S0401/
     + 9.4992E-06, 8.1581E-06, 7.1512E-06, 6.2692E-06, 5.5285E-06,
     + 4.9000E-06, 4.3447E-06, 3.8906E-06, 3.4679E-06, 3.1089E-06,
     + 2.8115E-06, 2.5496E-06, 2.2982E-06, 2.0861E-06, 1.8763E-06,
     + 1.7035E-06, 1.5548E-06, 1.4107E-06, 1.2839E-06, 1.1706E-06,
     + 1.0709E-06, 9.8099E-07, 8.9901E-07, 8.2394E-07, 7.5567E-07,
     + 6.9434E-07, 6.3867E-07, 5.8845E-07, 5.4263E-07, 5.0033E-07,
     + 4.6181E-07, 4.2652E-07, 3.9437E-07, 3.6497E-07, 3.3781E-07,
     + 3.1292E-07, 2.9011E-07, 2.6915E-07, 2.4989E-07, 2.3215E-07,
     + 2.1582E-07, 2.0081E-07, 1.8700E-07, 1.7432E-07, 1.6264E-07,
     + 1.5191E-07, 1.4207E-07, 1.3306E-07, 1.2484E-07, 1.1737E-07/
       DATA S0451/
     + 1.1056E-07, 1.0451E-07, 9.9060E-08, 9.4135E-08, 8.9608E-08,
     + 8.5697E-08, 8.1945E-08, 7.8308E-08, 7.4808E-08, 7.1686E-08,
     + 6.8923E-08, 6.5869E-08, 6.3308E-08, 6.0840E-08, 5.8676E-08,
     + 5.6744E-08, 5.5016E-08, 5.3813E-08, 5.2792E-08, 5.2097E-08,
     + 5.1737E-08, 5.1603E-08, 5.1656E-08, 5.1989E-08, 5.2467E-08,
     + 5.2918E-08, 5.3589E-08, 5.4560E-08, 5.5869E-08, 5.7403E-08,
     + 5.8968E-08, 6.0973E-08, 6.3432E-08, 6.6245E-08, 6.9353E-08,
     + 7.2686E-08, 7.6541E-08, 8.0991E-08, 8.5950E-08, 9.1429E-08,
     + 9.7851E-08, 1.0516E-07, 1.1349E-07, 1.2295E-07, 1.3335E-07,
     + 1.4488E-07, 1.5864E-07, 1.7412E-07, 1.9140E-07, 2.1078E-07/
       DATA S0501/
     + 2.3369E-07, 2.5996E-07, 2.8848E-07, 3.2169E-07, 3.5991E-07,
     + 4.0566E-07, 4.5969E-07, 5.3094E-07, 6.1458E-07, 7.1155E-07,
     + 8.3045E-07, 9.9021E-07, 1.2042E-06, 1.4914E-06, 1.8145E-06,
     + 2.2210E-06, 2.7831E-06, 3.4533E-06, 4.4446E-06, 5.1989E-06,
     + 6.2289E-06, 7.1167E-06, 8.3949E-06, 9.6417E-06, 1.0313E-05,
     + 1.0485E-05, 1.0641E-05, 1.0898E-05, 1.0763E-05, 1.0506E-05,
     + 1.0497E-05, 1.1696E-05, 1.2654E-05, 1.3029E-05, 1.3175E-05,
     + 1.4264E-05, 1.4985E-05, 1.4999E-05, 1.4317E-05, 1.4616E-05,
     + 1.4963E-05, 1.5208E-05, 1.4942E-05, 1.3879E-05, 1.3087E-05,
     + 1.1727E-05, 1.0515E-05, 9.0073E-06, 7.3133E-06, 6.1181E-06/
       DATA S0551/
     + 5.0623E-06, 4.1105E-06, 3.3915E-06, 2.6711E-06, 2.1464E-06,
     + 1.7335E-06, 1.4302E-06, 1.1847E-06, 9.9434E-07, 8.2689E-07,
     + 7.0589E-07, 6.0750E-07, 5.3176E-07, 4.6936E-07, 4.1541E-07,
     + 3.6625E-07, 3.2509E-07, 2.9156E-07, 2.6308E-07, 2.3819E-07,
     + 2.1421E-07, 1.9366E-07, 1.7626E-07, 1.5982E-07, 1.4567E-07,
     + 1.3354E-07, 1.2097E-07, 1.1029E-07, 1.0063E-07, 9.2003E-08,
     + 8.4245E-08, 7.7004E-08, 7.0636E-08, 6.4923E-08, 5.9503E-08,
     + 5.4742E-08, 5.0450E-08, 4.6470E-08, 4.2881E-08, 3.9550E-08,
     + 3.6541E-08, 3.3803E-08, 3.1279E-08, 2.8955E-08, 2.6858E-08,
     + 2.4905E-08, 2.3146E-08, 2.1539E-08, 2.0079E-08, 1.8746E-08/
       DATA S0601/
     + 1.7517E-08, 1.6396E-08, 1.5369E-08, 1.4426E-08, 1.3543E-08,
     + 1.2724E-08, 1.1965E-08, 1.1267E-08, 1.0617E-08, 1.0010E-08,
     + 9.4662E-09, 8.9553E-09, 8.4988E-09, 8.0807E-09, 7.7043E-09,
     + 7.3721E-09, 7.0707E-09, 6.8047E-09, 6.5702E-09, 6.3634E-09,
     + 6.1817E-09, 6.0239E-09, 5.8922E-09, 5.7824E-09, 5.7019E-09,
     + 5.6368E-09, 5.5940E-09, 5.5669E-09, 5.5583E-09, 5.5653E-09,
     + 5.5837E-09, 5.6243E-09, 5.6883E-09, 5.7800E-09, 5.8964E-09,
     + 6.0429E-09, 6.2211E-09, 6.4282E-09, 6.6634E-09, 6.9306E-09,
     + 7.2336E-09, 7.5739E-09, 7.9562E-09, 8.3779E-09, 8.8575E-09,
     + 9.3992E-09, 1.0004E-08, 1.0684E-08, 1.1450E-08, 1.2320E-08/
       DATA S0651/
     + 1.3311E-08, 1.4455E-08, 1.5758E-08, 1.7254E-08, 1.8927E-08,
     + 2.0930E-08, 2.3348E-08, 2.6074E-08, 2.9221E-08, 3.2770E-08,
     + 3.7485E-08, 4.2569E-08, 4.8981E-08, 5.5606E-08, 6.2393E-08,
     + 7.1901E-08, 8.2921E-08, 9.5513E-08, 1.1111E-07, 1.3143E-07,
     + 1.5971E-07, 1.8927E-07, 2.2643E-07, 2.7860E-07, 3.2591E-07,
     + 3.7024E-07, 4.2059E-07, 4.9432E-07, 5.5543E-07, 5.7498E-07,
     + 5.9210E-07, 6.1005E-07, 6.1577E-07, 5.9193E-07, 5.6602E-07,
     + 5.7403E-07, 6.0050E-07, 6.4723E-07, 6.7073E-07, 7.5415E-07,
     + 8.0982E-07, 8.7658E-07, 9.1430E-07, 9.4459E-07, 9.8347E-07,
     + 9.8768E-07, 1.0153E-06, 1.0066E-06, 1.0353E-06, 1.0353E-06/
       DATA S0701/
     + 1.0722E-06, 1.1138E-06, 1.1923E-06, 1.2947E-06, 1.4431E-06,
     + 1.6537E-06, 1.8662E-06, 2.2473E-06, 2.6464E-06, 3.1041E-06,
     + 3.4858E-06, 4.0167E-06, 4.6675E-06, 5.0983E-06, 5.7997E-06,
     + 6.0503E-06, 6.4687E-06, 6.5396E-06, 6.7986E-06, 7.0244E-06,
     + 7.2305E-06, 7.6732E-06, 7.9783E-06, 7.9846E-06, 7.7617E-06,
     + 7.7657E-06, 7.7411E-06, 7.8816E-06, 7.8136E-06, 8.0051E-06,
     + 8.5799E-06, 9.1659E-06, 9.8646E-06, 9.4920E-06, 8.7670E-06,
     + 8.2034E-06, 7.2297E-06, 6.2324E-06, 4.9315E-06, 3.9128E-06,
     + 3.1517E-06, 2.4469E-06, 1.8815E-06, 1.4627E-06, 1.1698E-06,
     + 9.4686E-07, 7.8486E-07, 6.6970E-07, 5.8811E-07, 5.2198E-07/
       DATA S0751/
     + 4.6809E-07, 4.1671E-07, 3.7006E-07, 3.3066E-07, 2.9387E-07,
     + 2.6415E-07, 2.3409E-07, 2.0991E-07, 1.9132E-07, 1.7519E-07,
     + 1.5939E-07, 1.4368E-07, 1.3050E-07, 1.1883E-07, 1.0772E-07,
     + 9.6884E-08, 8.7888E-08, 7.8956E-08, 7.1024E-08, 6.3824E-08,
     + 5.7256E-08, 5.1769E-08, 4.7037E-08, 4.2901E-08, 3.8970E-08,
     + 3.5467E-08, 3.2502E-08, 2.9827E-08, 2.7389E-08, 2.5111E-08,
     + 2.3056E-08, 2.1267E-08, 1.9610E-08, 1.8133E-08, 1.6775E-08,
     + 1.5491E-08, 1.4329E-08, 1.3265E-08, 1.2300E-08, 1.1420E-08,
     + 1.0593E-08, 9.8475E-09, 9.1585E-09, 8.5256E-09, 7.9525E-09,
     + 7.4226E-09, 6.9379E-09, 6.4950E-09, 6.0911E-09, 5.7242E-09/
       DATA S0801/
     + 5.3877E-09, 5.0821E-09, 4.8051E-09, 4.5554E-09, 4.3315E-09,
     + 4.1336E-09, 3.9632E-09, 3.8185E-09, 3.7080E-09, 3.6296E-09,
     + 3.5804E-09, 3.5776E-09, 3.6253E-09, 3.7115E-09, 3.8151E-09,
     + 3.9804E-09, 4.1742E-09, 4.3581E-09, 4.5306E-09, 4.7736E-09,
     + 5.1297E-09, 5.5291E-09, 5.9125E-09, 6.4956E-09, 7.0362E-09,
     + 7.5318E-09, 7.9947E-09, 8.6438E-09, 9.7227E-09, 1.0130E-08,
     + 1.0549E-08, 1.1064E-08, 1.1702E-08, 1.2043E-08, 1.1781E-08,
     + 1.1838E-08, 1.1917E-08, 1.2131E-08, 1.2476E-08, 1.3611E-08,
     + 1.4360E-08, 1.5057E-08, 1.6247E-08, 1.7284E-08, 1.8420E-08,
     + 1.8352E-08, 1.8722E-08, 1.9112E-08, 1.9092E-08, 1.9311E-08/
       DATA S0851/
     + 1.9411E-08, 1.9884E-08, 2.0508E-08, 2.1510E-08, 2.3143E-08,
     + 2.5050E-08, 2.7596E-08, 3.1231E-08, 3.6260E-08, 4.3410E-08,
     + 5.2240E-08, 6.3236E-08, 7.7522E-08, 9.8688E-08, 1.1859E-07,
     + 1.4341E-07, 1.6798E-07, 1.9825E-07, 2.2898E-07, 2.6257E-07,
     + 2.9884E-07, 3.3247E-07, 3.4936E-07, 3.5583E-07, 3.7150E-07,
     + 3.6580E-07, 3.7124E-07, 3.7030E-07, 4.1536E-07, 4.6656E-07,
     + 4.6677E-07, 4.7507E-07, 4.9653E-07, 5.3795E-07, 5.4957E-07,
     + 5.2238E-07, 5.4690E-07, 5.6569E-07, 5.9844E-07, 5.9835E-07,
     + 5.6522E-07, 5.4123E-07, 4.7904E-07, 4.2851E-07, 3.5603E-07,
     + 2.8932E-07, 2.3655E-07, 1.8592E-07, 1.4943E-07, 1.1971E-07/
       DATA S0901/
     + 9.8482E-08, 8.3675E-08, 7.1270E-08, 6.2496E-08, 5.4999E-08,
     + 4.9821E-08, 4.5387E-08, 4.1340E-08, 3.7453E-08, 3.3298E-08,
     + 3.0120E-08, 2.7032E-08, 2.4236E-08, 2.1500E-08, 1.8988E-08,
     + 1.7414E-08, 1.5706E-08, 1.4192E-08, 1.3204E-08, 1.1759E-08,
     + 1.0737E-08, 9.6309E-09, 8.8179E-09, 8.2619E-09, 7.2264E-09,
     + 6.4856E-09, 5.8037E-09, 5.2093E-09, 4.7205E-09, 4.1749E-09,
     + 3.7852E-09, 3.3915E-09, 3.0089E-09, 2.7335E-09, 2.4398E-09,
     + 2.2031E-09, 1.9786E-09, 1.7890E-09, 1.6266E-09, 1.4830E-09,
     + 1.3576E-09, 1.2518E-09, 1.1587E-09, 1.0726E-09, 9.9106E-10,
     + 9.1673E-10, 8.5084E-10, 7.9147E-10, 7.2882E-10, 6.7342E-10/
       DATA S0951/
     + 6.2593E-10, 5.8294E-10, 5.4435E-10, 5.0997E-10, 4.7806E-10,
     + 4.4931E-10, 4.2357E-10, 4.0023E-10, 3.7909E-10, 3.5999E-10,
     + 3.4285E-10, 3.2776E-10, 3.1468E-10, 3.0377E-10, 2.9479E-10,
     + 2.8877E-10, 2.8512E-10, 2.8617E-10, 2.8976E-10, 3.0001E-10,
     + 3.1718E-10, 3.3898E-10, 3.5857E-10, 3.8358E-10, 4.3131E-10,
     + 4.5741E-10, 4.6948E-10, 4.7594E-10, 4.9529E-10, 5.1563E-10,
     + 4.9475E-10, 4.8369E-10, 4.8829E-10, 5.0047E-10, 5.0203E-10,
     + 5.1954E-10, 5.5352E-10, 5.9928E-10, 6.7148E-10, 7.1121E-10,
     + 7.4317E-10, 7.6039E-10, 7.8313E-10, 8.0684E-10, 7.8553E-10,
     + 7.8312E-10, 7.8537E-10, 7.8872E-10, 8.0185E-10, 8.1004E-10/
       DATA S1001/
     + 8.2608E-10, 8.2525E-10, 8.3857E-10, 8.7920E-10, 9.2451E-10,
     + 9.8661E-10, 1.0629E-09, 1.1659E-09, 1.2922E-09, 1.4387E-09,
     + 1.6254E-09, 1.8425E-09, 2.1428E-09, 2.5477E-09, 3.0379E-09,
     + 3.7570E-09, 4.4354E-09, 5.1802E-09, 6.2769E-09, 7.4894E-09,
     + 8.7474E-09, 9.8037E-09, 1.1582E-08, 1.3293E-08, 1.4471E-08,
     + 1.5025E-08, 1.5580E-08, 1.6228E-08, 1.6413E-08, 1.6020E-08,
     + 1.6393E-08, 1.7545E-08, 1.9590E-08, 2.1449E-08, 2.3856E-08,
     + 2.7050E-08, 3.0214E-08, 3.3733E-08, 3.6487E-08, 3.9353E-08,
     + 4.2660E-08, 4.6385E-08, 4.9955E-08, 5.5313E-08, 6.0923E-08,
     + 6.8948E-08, 7.3649E-08, 8.2602E-08, 9.2212E-08, 9.9080E-08/
       DATA S1051/
     + 1.1319E-07, 1.1790E-07, 1.2941E-07, 1.3199E-07, 1.3914E-07,
     + 1.4843E-07, 1.5300E-07, 1.6419E-07, 1.7095E-07, 1.6988E-07,
     + 1.6494E-07, 1.6327E-07, 1.6067E-07, 1.6909E-07, 1.7118E-07,
     + 1.8106E-07, 1.9857E-07, 2.1696E-07, 2.3385E-07, 2.2776E-07,
     + 2.1402E-07, 1.9882E-07, 1.7362E-07, 1.4308E-07, 1.1158E-07,
     + 8.8781E-08, 6.8689E-08, 5.2062E-08, 4.0427E-08, 3.2669E-08,
     + 2.7354E-08, 2.3200E-08, 2.0580E-08, 1.8676E-08, 1.7329E-08,
     + 1.6621E-08, 1.6433E-08, 1.6953E-08, 1.7134E-08, 1.7948E-08,
     + 1.9107E-08, 1.9875E-08, 2.1416E-08, 2.1556E-08, 2.2265E-08,
     + 2.2171E-08, 2.2534E-08, 2.3029E-08, 2.2828E-08, 2.3143E-08/
       DATA S1101/
     + 2.2965E-08, 2.2223E-08, 2.1108E-08, 2.0265E-08, 1.9516E-08,
     + 1.9941E-08, 2.0312E-08, 2.1080E-08, 2.2611E-08, 2.4210E-08,
     + 2.6069E-08, 2.5097E-08, 2.3318E-08, 2.1543E-08, 1.8942E-08,
     + 1.5960E-08, 1.2386E-08, 9.9340E-09, 7.7502E-09, 5.9462E-09,
     + 4.5113E-09, 3.5523E-09, 2.8844E-09, 2.3394E-09, 1.9584E-09,
     + 1.6749E-09, 1.4624E-09, 1.2809E-09, 1.1359E-09, 1.0087E-09,
     + 9.0166E-10, 8.1079E-10, 7.2219E-10, 6.4922E-10, 5.8803E-10,
     + 5.3290E-10, 4.8590E-10, 4.4111E-10, 4.0184E-10, 3.6644E-10,
     + 3.3529E-10, 3.0789E-10, 2.8286E-10, 2.6089E-10, 2.4125E-10,
     + 2.2355E-10, 2.0783E-10, 1.9370E-10, 1.8088E-10, 1.6948E-10/
       DATA S1151/
     + 1.5929E-10, 1.5013E-10, 1.4193E-10, 1.3470E-10, 1.2841E-10,
     + 1.2307E-10, 1.1865E-10, 1.1502E-10, 1.1243E-10, 1.1099E-10,
     + 1.1066E-10, 1.1216E-10, 1.1529E-10, 1.2171E-10, 1.3128E-10,
     + 1.4153E-10, 1.5962E-10, 1.8048E-10, 2.0936E-10, 2.3165E-10,
     + 2.5746E-10, 2.9600E-10, 3.3707E-10, 3.5267E-10, 3.5953E-10,
     + 3.6822E-10, 3.8363E-10, 3.8286E-10, 3.5883E-10, 3.6154E-10,
     + 3.6653E-10, 3.8507E-10, 4.0250E-10, 4.4435E-10, 4.9889E-10,
     + 5.6932E-10, 6.3599E-10, 7.0281E-10, 7.5777E-10, 8.1279E-10,
     + 8.8910E-10, 9.3400E-10, 1.0076E-09, 1.0945E-09, 1.1898E-09,
     + 1.3108E-09, 1.4725E-09, 1.7028E-09, 1.9619E-09, 2.3527E-09/
       DATA S1201/
     + 2.6488E-09, 3.0327E-09, 3.4396E-09, 3.8797E-09, 4.4115E-09,
     + 4.6853E-09, 4.9553E-09, 4.9551E-09, 5.1062E-09, 5.0996E-09,
     + 5.1119E-09, 5.2283E-09, 5.8297E-09, 6.3439E-09, 6.2675E-09,
     + 6.3296E-09, 6.5173E-09, 7.1685E-09, 7.0528E-09, 6.8856E-09,
     + 7.3182E-09, 7.6990E-09, 8.3461E-09, 8.1946E-09, 7.7153E-09,
     + 7.2411E-09, 6.4511E-09, 5.7336E-09, 4.6105E-09, 3.6962E-09,
     + 2.9944E-09, 2.4317E-09, 1.9399E-09, 1.5331E-09, 1.2633E-09,
     + 1.0613E-09, 9.0136E-10, 7.9313E-10, 7.1543E-10, 6.6485E-10,
     + 6.4225E-10, 6.3980E-10, 6.4598E-10, 6.7428E-10, 7.0270E-10,
     + 7.4694E-10, 7.7946E-10, 7.9395E-10, 7.8716E-10, 7.6933E-10/
       DATA S1251/
     + 7.6220E-10, 7.4825E-10, 7.4805E-10, 7.6511E-10, 7.6492E-10,
     + 7.4103E-10, 7.1979E-10, 7.1686E-10, 7.3403E-10, 7.1142E-10,
     + 7.0212E-10, 7.1548E-10, 7.5253E-10, 8.0444E-10, 8.2378E-10,
     + 7.8004E-10, 7.1712E-10, 6.4978E-10, 5.7573E-10, 4.8675E-10,
     + 3.7945E-10, 3.0118E-10, 2.4241E-10, 1.9100E-10, 1.4816E-10,
     + 1.1567E-10, 9.4183E-11, 7.7660E-11, 6.5270E-11, 5.6616E-11,
     + 4.9576E-11, 4.4137E-11, 3.9459E-11, 3.5759E-11, 3.2478E-11,
     + 2.9419E-11, 2.6703E-11, 2.4365E-11, 2.2412E-11, 2.0606E-11,
     + 1.9067E-11, 1.7800E-11, 1.6695E-11, 1.5729E-11, 1.4887E-11,
     + 1.4135E-11, 1.3519E-11, 1.2992E-11, 1.2563E-11, 1.2223E-11/
       DATA S1301/
     + 1.1962E-11, 1.1775E-11, 1.1657E-11, 1.1605E-11, 1.1619E-11,
     + 1.1697E-11, 1.1839E-11, 1.2046E-11, 1.2319E-11, 1.2659E-11,
     + 1.3070E-11, 1.3553E-11, 1.4113E-11, 1.4754E-11, 1.5480E-11,
     + 1.6298E-11, 1.7214E-11, 1.8236E-11, 1.9372E-11, 2.0635E-11,
     + 2.2036E-11, 2.3590E-11, 2.5317E-11, 2.7242E-11, 2.9400E-11,
     + 3.1849E-11, 3.4654E-11, 3.7923E-11, 4.1695E-11, 4.6055E-11,
     + 5.0940E-11, 5.5624E-11, 6.0667E-11, 6.6261E-11, 7.2692E-11,
     + 7.9711E-11, 8.7976E-11, 9.6884E-11, 1.0775E-10, 1.2093E-10,
     + 1.3531E-10, 1.5404E-10, 1.7315E-10, 1.9862E-10, 2.3341E-10,
     + 2.7014E-10, 3.1716E-10, 3.6957E-10, 4.3233E-10, 5.2566E-10/
       DATA S1351/
     + 6.2251E-10, 7.2149E-10, 8.3958E-10, 9.5931E-10, 1.1388E-09,
     + 1.2973E-09, 1.4442E-09, 1.5638E-09, 1.6974E-09, 1.8489E-09,
     + 1.9830E-09, 2.1720E-09, 2.3662E-09, 2.6987E-09, 3.1697E-09,
     + 3.6907E-09, 4.2625E-09, 4.7946E-09, 5.3848E-09, 6.0897E-09,
     + 6.4730E-09, 7.1483E-09, 7.7432E-09, 8.0851E-09, 8.5013E-09,
     + 8.5909E-09, 9.1890E-09, 9.3124E-09, 9.5936E-09, 9.8787E-09,
     + 9.9036E-09, 9.6712E-09, 9.2036E-09, 9.0466E-09, 8.9380E-09,
     + 9.1815E-09, 9.5092E-09, 1.0027E-08, 1.0876E-08, 1.1744E-08,
     + 1.1853E-08, 1.1296E-08, 1.0134E-08, 8.8245E-09, 7.3930E-09,
     + 5.7150E-09, 4.4884E-09, 3.4027E-09, 2.6054E-09, 2.0790E-09/
       DATA S1401/
     + 1.7267E-09, 1.4724E-09, 1.2722E-09, 1.1234E-09, 1.0186E-09,
     + 9.4680E-10, 8.8854E-10, 8.5127E-10, 8.3157E-10, 8.2226E-10,
     + 8.3395E-10, 8.3294E-10, 8.4725E-10, 8.8814E-10, 9.3697E-10,
     + 1.0112E-09, 1.0412E-09, 1.0948E-09, 1.1810E-09, 1.2267E-09,
     + 1.3690E-09, 1.4512E-09, 1.5568E-09, 1.6552E-09, 1.7321E-09,
     + 1.8797E-09, 1.9210E-09, 1.9686E-09, 1.9917E-09, 1.9357E-09,
     + 1.8486E-09, 1.7575E-09, 1.7113E-09, 1.7163E-09, 1.7623E-09,
     + 1.8536E-09, 1.9765E-09, 2.1334E-09, 2.3237E-09, 2.3259E-09,
     + 2.1833E-09, 1.9785E-09, 1.7308E-09, 1.4596E-09, 1.1198E-09,
     + 8.7375E-10, 6.5381E-10, 4.8677E-10, 3.6756E-10, 2.9155E-10/
       DATA S1451/
     + 2.3735E-10, 1.9590E-10, 1.6638E-10, 1.4549E-10, 1.2947E-10,
     + 1.1511E-10, 1.0548E-10, 9.6511E-11, 9.0469E-11, 8.5170E-11,
     + 7.7804E-11, 7.1971E-11, 6.6213E-11, 6.1063E-11, 5.5881E-11,
     + 5.0508E-11, 4.5932E-11, 4.1997E-11, 3.7672E-11, 3.3972E-11,
     + 3.0318E-11, 2.6769E-11, 2.3874E-11, 2.1336E-11, 1.9073E-11,
     + 1.7313E-11, 1.5904E-11, 1.4684E-11, 1.3698E-11, 1.2873E-11,
     + 1.2175E-11, 1.1542E-11, 1.1024E-11, 1.0602E-11, 1.0267E-11,
     + 1.0012E-11, 9.8379E-12, 9.7482E-12, 9.7564E-12, 9.8613E-12,
     + 1.0092E-11, 1.0418E-11, 1.0868E-11, 1.1585E-11, 1.2351E-11,
     + 1.3372E-11, 1.4841E-11, 1.6457E-11, 1.8681E-11, 2.0550E-11/
       DATA S1501/
     + 2.2912E-11, 2.5958E-11, 2.9137E-11, 3.2368E-11, 3.4848E-11,
     + 3.8462E-11, 4.2190E-11, 4.5629E-11, 4.9022E-11, 5.4232E-11,
     + 6.1900E-11, 7.1953E-11, 8.5368E-11, 9.9699E-11, 1.1734E-10,
     + 1.4185E-10, 1.7017E-10, 1.9813E-10, 2.3859E-10, 2.7304E-10,
     + 3.0971E-10, 3.5129E-10, 3.9405E-10, 4.5194E-10, 4.8932E-10,
     + 5.2436E-10, 5.4098E-10, 5.5542E-10, 5.7794E-10, 5.6992E-10,
     + 5.8790E-10, 6.1526E-10, 6.8034E-10, 6.7956E-10, 6.6864E-10,
     + 6.9329E-10, 7.2971E-10, 7.6546E-10, 7.5078E-10, 7.8406E-10,
     + 8.3896E-10, 9.0111E-10, 9.1994E-10, 8.7189E-10, 8.1426E-10,
     + 7.3097E-10, 6.3357E-10, 5.1371E-10, 4.0936E-10, 3.2918E-10/
       DATA S1551/
     + 2.6255E-10, 2.0724E-10, 1.6879E-10, 1.4165E-10, 1.1989E-10,
     + 1.0125E-10, 8.9629E-11, 7.8458E-11, 6.8826E-11, 6.0935E-11,
     + 5.5208E-11, 5.2262E-11, 5.0260E-11, 4.8457E-11, 4.7888E-11,
     + 4.8032E-11, 5.0838E-11, 5.4668E-11, 5.5790E-11, 6.0056E-11,
     + 6.3811E-11, 6.8848E-11, 7.4590E-11, 7.8249E-11, 8.3371E-11,
     + 8.3641E-11, 8.6591E-11, 8.9599E-11, 9.3487E-11, 1.0066E-10,
     + 1.0765E-10, 1.0851E-10, 1.0619E-10, 1.0557E-10, 1.0460E-10,
     + 1.0796E-10, 1.0523E-10, 1.0674E-10, 1.1261E-10, 1.1431E-10,
     + 1.1408E-10, 1.0901E-10, 9.9105E-11, 8.8077E-11, 6.9928E-11,
     + 5.4595E-11, 4.5401E-11, 3.6313E-11, 2.6986E-11, 1.9463E-11/
       DATA S1601/
     + 1.4577E-11, 1.1583E-11, 9.5492E-12, 8.0770E-12, 6.9642E-12,
     + 6.0966E-12, 5.4046E-12, 4.8431E-12, 4.3815E-12, 3.9987E-12,
     + 3.6790E-12, 3.4113E-12, 3.1868E-12, 2.9992E-12, 2.8434E-12,
     + 2.7153E-12, 2.6120E-12, 2.5311E-12, 2.4705E-12, 2.4290E-12,
     + 2.4053E-12, 2.3988E-12, 2.4087E-12, 2.4349E-12, 2.4771E-12,
     + 2.5355E-12, 2.6103E-12, 2.7019E-12, 2.8110E-12, 2.9383E-12,
     + 3.0848E-12, 3.2518E-12, 3.4405E-12, 3.6527E-12, 3.8902E-12,
     + 4.1555E-12, 4.4510E-12, 4.7801E-12, 5.1462E-12, 5.5539E-12,
     + 6.0086E-12, 6.5171E-12, 7.0884E-12, 7.7357E-12, 8.4831E-12,
     + 9.3096E-12, 1.0282E-11, 1.1407E-11, 1.2690E-11, 1.4148E-11/
       DATA S1651/
     + 1.5888E-11, 1.7992E-11, 2.0523E-11, 2.3342E-11, 2.6578E-11,
     + 3.0909E-11, 3.6228E-11, 4.2053E-11, 4.9059E-11, 5.9273E-11,
     + 7.0166E-11, 8.2298E-11, 9.7071E-11, 1.1673E-10, 1.4010E-10,
     + 1.6621E-10, 2.0127E-10, 2.3586E-10, 2.7050E-10, 3.0950E-10,
     + 3.6584E-10, 4.1278E-10, 4.6591E-10, 5.2220E-10, 5.5246E-10,
     + 6.1500E-10, 6.5878E-10, 7.1167E-10, 7.9372E-10, 8.6975E-10,
     + 9.6459E-10, 9.7368E-10, 9.8142E-10, 1.0202E-09, 1.0200E-09,
     + 1.0356E-09, 1.0092E-09, 1.0269E-09, 1.0366E-09, 1.0490E-09,
     + 1.0717E-09, 1.0792E-09, 1.1016E-09, 1.0849E-09, 1.0929E-09,
     + 1.0971E-09, 1.0969E-09, 1.0460E-09, 9.2026E-10, 8.1113E-10/
       DATA S1701/
     + 6.8635E-10, 5.5369E-10, 4.2908E-10, 3.3384E-10, 2.6480E-10,
     + 2.0810E-10, 1.6915E-10, 1.4051E-10, 1.1867E-10, 1.0158E-10,
     + 8.8990E-11, 7.9175E-11, 7.0440E-11, 6.3453E-11, 5.7009E-11,
     + 5.1662E-11, 4.7219E-11, 4.3454E-11, 4.0229E-11, 3.7689E-11,
     + 3.6567E-11, 3.5865E-11, 3.5955E-11, 3.5928E-11, 3.6298E-11,
     + 3.7629E-11, 3.9300E-11, 4.1829E-11, 4.4806E-11, 5.0534E-11,
     + 5.6672E-11, 6.2138E-11, 6.8678E-11, 7.6111E-11, 8.4591E-11,
     + 9.2634E-11, 9.8085E-11, 1.0830E-10, 1.1949E-10, 1.2511E-10,
     + 1.3394E-10, 1.3505E-10, 1.4342E-10, 1.4874E-10, 1.4920E-10,
     + 1.5872E-10, 1.5972E-10, 1.5821E-10, 1.5425E-10, 1.4937E-10/
       DATA S1751/
     + 1.5089E-10, 1.5521E-10, 1.6325E-10, 1.6924E-10, 1.8265E-10,
     + 1.9612E-10, 2.0176E-10, 1.9359E-10, 1.7085E-10, 1.5197E-10,
     + 1.2646E-10, 9.8552E-11, 7.4530E-11, 5.5052E-11, 4.2315E-11,
     + 3.2736E-11, 2.6171E-11, 2.1909E-11, 1.8286E-11, 1.5752E-11,
     + 1.3859E-11, 1.2288E-11, 1.1002E-11, 9.7534E-12, 8.8412E-12,
     + 8.0169E-12, 7.2855E-12, 6.8734E-12, 6.4121E-12, 6.1471E-12,
     + 5.7780E-12, 5.3478E-12, 4.9652E-12, 4.4043E-12, 3.9862E-12,
     + 3.4684E-12, 2.9681E-12, 2.5791E-12, 2.2339E-12, 1.9247E-12,
     + 1.6849E-12, 1.4863E-12, 1.3291E-12, 1.2021E-12, 1.0947E-12,
     + 1.0015E-12, 9.1935E-13, 8.4612E-13, 7.8036E-13, 7.2100E-13/
       DATA S1801/
     + 6.6718E-13, 6.1821E-13, 5.7353E-13, 5.3269E-13, 4.9526E-13,
     + 4.6093E-13, 4.2937E-13, 4.0034E-13, 3.7361E-13, 3.4895E-13,
     + 3.2621E-13, 3.0520E-13, 2.8578E-13, 2.6782E-13, 2.5120E-13,
     + 2.3581E-13, 2.2154E-13, 2.0832E-13, 1.9605E-13, 1.8466E-13,
     + 1.7408E-13, 1.6425E-13, 1.5511E-13, 1.4661E-13, 1.3869E-13,
     + 1.3131E-13, 1.2444E-13, 1.1803E-13, 1.1205E-13, 1.0646E-13,
     + 1.0124E-13, 9.6358E-14, 9.1789E-14, 8.7509E-14, 8.3498E-14,
     + 7.9735E-14, 7.6202E-14, 7.2882E-14, 6.9760E-14, 6.6822E-14,
     + 6.4053E-14, 6.1442E-14, 5.8978E-14, 5.6650E-14, 5.4448E-14,
     + 5.2364E-14, 5.0389E-14, 4.8516E-14, 4.6738E-14, 4.5048E-14/
       DATA S1851/
     + 4.3441E-14, 4.1911E-14, 4.0453E-14, 3.9063E-14, 3.7735E-14,
     + 3.6467E-14, 3.5254E-14, 3.4093E-14, 3.2980E-14, 3.1914E-14,
     + 3.0891E-14, 2.9909E-14, 2.8965E-14, 2.8058E-14, 2.7185E-14,
     + 2.6344E-14, 2.5535E-14, 2.4755E-14, 2.4002E-14, 2.3276E-14,
     + 2.2576E-14, 2.1899E-14, 2.1245E-14, 2.0613E-14, 2.0002E-14,
     + 1.9411E-14, 1.8839E-14, 1.8285E-14, 1.7749E-14, 1.7230E-14,
     + 1.6727E-14, 1.6240E-14, 1.5768E-14, 1.5310E-14, 1.4867E-14,
     + 1.4436E-14, 1.4019E-14, 1.3614E-14, 1.3221E-14, 1.2840E-14,
     + 1.2471E-14, 1.2112E-14, 1.1764E-14, 1.1425E-14, 1.1097E-14,
     + 1.0779E-14, 1.0469E-14, 1.0169E-14, 9.8775E-15, 9.5943E-15/
       DATA S1901/
     + 9.3193E-15, 9.0522E-15, 8.7928E-15, 8.5409E-15, 8.2962E-15,
     + 8.0586E-15, 7.8278E-15, 7.6036E-15, 7.3858E-15, 7.1742E-15,
     + 6.9687E-15, 6.7691E-15, 6.5752E-15, 6.3868E-15, 6.2038E-15,
     + 6.0260E-15, 5.8533E-15, 5.6856E-15, 5.5226E-15, 5.3642E-15,
     + 5.2104E-15, 5.0610E-15, 4.9158E-15, 4.7748E-15, 4.6378E-15,
     + 4.5047E-15, 4.3753E-15, 4.2497E-15, 4.1277E-15, 4.0091E-15,
     + 3.8939E-15, 3.7820E-15, 3.6733E-15, 3.5677E-15, 3.4651E-15,
     + 3.3655E-15, 3.2686E-15, 3.1746E-15, 3.0832E-15, 2.9944E-15,
     + 2.9082E-15, 2.8244E-15, 2.7431E-15, 2.6640E-15, 2.5872E-15,
     + 2.5126E-15, 2.4401E-15, 2.3697E-15, 2.3014E-15, 2.2349E-15/
       DATA S1951/
     + 2.1704E-15, 2.1077E-15, 2.0468E-15, 1.9877E-15, 1.9302E-15,
     + 1.8744E-15, 1.8202E-15, 1.7675E-15, 1.7164E-15, 1.6667E-15,
     + 1.6184E-15, 1.5716E-15, 1.5260E-15, 1.4818E-15, 1.4389E-15,
     + 1.3971E-15, 1.3566E-15, 1.3172E-15, 1.2790E-15, 1.2419E-15,
     + 1.2058E-15, 1.1708E-15, 1.1368E-15, 1.1037E-15, 1.0716E-15,
     + 1.0405E-15, 1.0102E-15, 9.8079E-16, 9.5224E-16, 9.2451E-16,
     + 8.9758E-16, 8.7142E-16, 8.4602E-16, 8.2136E-16, 7.9740E-16,
     + 7.7414E-16, 7.5154E-16, 7.2961E-16, 7.0830E-16, 6.8761E-16,
     + 6.6752E-16, 6.4801E-16, 6.2906E-16, 6.1066E-16, 5.9280E-16,
     + 5.7545E-16, 5.5860E-16, 5.4224E-16, 5.2636E-16, 5.1094E-16/
       DATA S2001/
     + 4.9596E-16/
C
       END
