#!/bin/bash
#SBATCH --time=24:00:00            #8:00:00        #16:00:00             
#SBATCH --job-name=tls_sims
#SBATCH --output=/data/gpfs/projects/punim1783/log/2023/chunk9/%j       
#SBATCH --ntasks=1 
#SBATCH --array=974,1003,1012,1042,1046,1051,1058,1096,1098,1122,1123,1124,1126,1129,1136,1137,1138,1147,1161,1162,1164,1166,1167,1170,1171,1173,1174,1175,1176,1177,1180,1201,1202,1204,1205,1206,1207,1209,1210,1211,1212,1213,1214,1215,1216,1218,1219,1220,1238,1241,1242,1243,1244,1245,1246,1247,1248,1249,1250,1251,1252,1254,1255,1256,1258,1259,1260,1268,1281,1282,1283,1284,1285,1286,1287,1288,1289,1290,1291,1292,1293,1294,1295,1296,1297,1298,1299,1300,1321,1322,1323,1324,1325,1326,1327,1328,1329,1330,1331,1332,1333,1334,1335,1336,1337,1338,1339,1340,1361,1362,1363,1364,1365,1366,1367,1368,1369,1370,1371,1372,1373,1374,1375,1376,1377,1378,1379,1380,1401,1402,1403,1404,1405,1406,1407,1408,1409,1410,1411,1412,1413,1414,1415,1416,1417,1418,1419,1420,1441,1442,1443,1444,1445,1446,1447,1448,1449,1450,1451,1452,1453,1454,1455,1456,1457,1458,1459,1460,1481,1482,1483,1484,1485,1486,1487,1488,1489,1490,1491,1492,1493,1494,1495,1496,1497,1498,1499,1500,1521,1522,1523,1524,1525,1526,1527,1528,1529,1530,1531,1532,1533,1534,1535,1536,1537,1538,1539,1540,1561,1562,1563,1564,1565,1566,1567,1568,1569,1570,1571,1572,1573,1574,1575,1576,1577,1578,1579,1580,1601,1602,1603,1604,1605,1606,1607,1608,1609,1610,1611,1612,1613,1614,1615,1616,1617,1618,1619,1620,1641,1642,1643,1644,1645,1646,1647,1648,1649,1650,1651,1652,1653,1654,1655,1656,1657,1658,1659,1660,1681,1682,1683,1685,1686,1688,1689,1690,1691,1692,1693,1694,1695,1696,1697,1698,1699,1700,1710,1722,1723,1724,1725,1726,1727,1732,1733,1734,1737,1739,1763,1794,1858,1873

# Things to set 
# chunk=9                     1:22  (noting chunk 22 only has 840 rows)
# version="rho_sep"              "rho_sep" | "actually_rep"
# network_dynamics="NIM"         "min_method" | "NIM"
# ostracism_type="original"      "original" | "FEJ"
# sbatch --job-name=tls_sims_$chunk --export=chunk=$chunk,version=$version,network_dynamics=$network_dynamics,ostracism_type=$ostracism_type /data/gpfs/projects/punim1783/jobs/2023/sims_chunked.sh

# version="rho_sep" network_dynamics="NIM" ostracism_type="original"
# have run ALL chunks: 22, 21, 20, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1
# Now doing timeouts... for 24 hours 
# done 22,21,20,19,18,17,16,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1, 


T="$(date +%s)"
SCRATCH_DIRECTORY=/data/gpfs/projects/punim1783/working/${SLURM_JOBID}
CODE_DIRECTORY=/data/gpfs/projects/punim1783/code/2023/
OUTPUT_DIRECTORY=/data/gpfs/projects/punim1783/output/2023/raw_output/${version}/${network_dynamics}/${ostracism_type}
mkdir -pv /data/gpfs/projects/punim1783/log/2023/chunk${chunk}/
mkdir -pv ${OUTPUT_DIRECTORY}
mkdir -pv ${SCRATCH_DIRECTORY}
cp ${CODE_DIRECTORY}/*.R* ${SCRATCH_DIRECTORY}
cd ${SCRATCH_DIRECTORY}
module load r/4.0.0 
echo "doing sim version: ${version} ; with tie strat: ${network_dynamics} and ostracism type: ${ostracism_type}"
echo "beginning chunk ${chunk} , run ${SLURM_ARRAY_TASK_ID}"
Rscript GO.R  ${chunk} ${version} ${network_dynamics}  ${ostracism_type} 
ls ${SCRATCH_DIRECTORY}/*.Rdata
cp ${SCRATCH_DIRECTORY}/*.Rdata  ${OUTPUT_DIRECTORY}
T="$(($(date +%s)-T))"
echo "run # ${SLURM_ARRAY_TASK_ID} took ${T} seconds"
rm -rf ${SCRATCH_DIRECTORY}
