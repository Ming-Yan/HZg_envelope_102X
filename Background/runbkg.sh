#read -p "file name:" fname
#read -p "vername:" vname1 
#read -p "vername2:" vname2
#read -p "period:"  period
#read -p  "name" vname
fname=ptwei
#fname=Legacy16_MVA0210
#fname=Rereco18_MVA0210
period=2020
#period=2018
#vname2=${vname}
vname2=turnon
vname1=turnon
#vname1=newrefit_woboost_comb_newcat01367
#e
#vname1=newrefit_woboost_comb_newcat0157
#vname1=mctest
#read -p "ptype:" ptype
#ptype=ele_mu
#read -p "cat:" cat
#nohup ./bin/fTest -i /afs/cern.ch/work/m/milee/MYcode/limit/PDFs/${fname}_${vname}/bkg_${fname}_WS_${ptype}_cat${cat}_${period}.root -D /afs/cern.ch/work/m/milee/MYcode/limit/bkgmodel/${fname}_${vname}/ -d /afs/cern.ch/work/m/milee/MYcode/limit/bkgmodel/${fname}_${vname}/bkg_${ptype}_${cat}.dat -t ${ptype}  --isData 1 --runFtestCheckWithToys -v --saveMultiPdf /afs/cern.ch/work/m/milee/MYcode/limit/bkgmodel/${fname}_${vname}/bkg_${ptype}_${cat}.root -f cat${cat} -p $period &>$ptype$cat${period} &
for ptype in ele_mu
do 
    for cat in  1 2 3 
   do
	#nohup ./bin/fTest_sb -i /afs/cern.ch/work/m/milee/MYcode/limit/PDFs/${fname}_${vname1}/bkg_${fname}_WS_${ptype}_cat${cat}_${period}.root -D /afs/cern.ch/work/m/milee/MYcode/limit/bkgmodel/${fname}_${vname1}/ -d /afs/cern.ch/work/m/milee/MYcode/limit/bkgmodel/${fname}_${vname1}/bkg_${ptype}_${cat}_sb.dat -t ${ptype}  --isData 1 --runFtestCheckWithToys -v --saveMultiPdf /afs/cern.ch/work/m/milee/MYcode/limit/bkgmodel/${fname}_${vname1}/bkg_${ptype}_${cat}_sb.root -f cat${cat} -p $period --sidebandOnly  &>$ptype$cat$period_sb &
	#nohup ./bin/fTest_loose -i /afs/cern.ch/work/m/milee/MYcode/limit/PDFs/${fname}_${vname1}/bkg_${fname}_WS_${ptype}_cat${cat}_${period}.root -D /afs/cern.ch/work/m/milee/MYcode/limit/bkgmodel/${fname}_${vname1}/ -d /afs/cern.ch/work/m/milee/MYcode/limit/bkgmodel/${fname}_${vname1}/bkg_${ptype}_${cat}_test.dat -t ${ptype}  --isData 1 --runFtestCheckWithToys -v --saveMultiPdf /afs/cern.ch/work/m/milee/MYcode/limit/bkgmodel/${fname}_${vname1}/bkg_${ptype}_${cat}_test.root -f cat${cat} -p $period &>$ptype$cat${period} &
	nohup ./bin/fTest_turn -i /afs/cern.ch/work/m/milee/BiasStudy/turnon_cat${cat}.root -D /afs/cern.ch/work/m/milee/MYcode/limit/bkgmodel/${fname}_${vname1}/ -d /afs/cern.ch/work/m/milee/MYcode/limit/bkgmodel/${fname}_${vname1}/bkg_${ptype}_${cat}_turn.dat -t ${ptype}  --isData 1 --runFtestCheckWithToys -v --saveMultiPdf /afs/cern.ch/work/m/milee/MYcode/limit/bkgmodel/${fname}_${vname1}/bkg_${ptype}_${cat}_turn.root -f cat${cat} -p $period &>$ptype$cat${period} &
#	nohup ./bin/fTest_bin -i /afs/cern.ch/work/m/milee/MYcode/limit/PDFs/${fname}_${vname1}/bkg_${fname}_WS_${ptype}_cat${cat}_${period}.root -D /afs/cern.ch/work/m/milee/MYcode/limit/bkgmodel/${fname}_${vname1}/ -d /afs/cern.ch/work/m/milee/MYcode/limit/bkgmodel/${fname}_${vname1}/bkg_${ptype}_${cat}.dat -t ${ptype}  --isData 1 --runFtestCheckWithToys -v --saveMultiPdf /afs/cern.ch/work/m/milee/MYcode/limit/bkgmodel/${fname}_${vname1}/bkg_${ptype}_${cat}_test.root -f cat${cat} -p $period &>$ptype$cat${period} &
    done       
    '''for cat in 503 502 501
    do
    #   nohup ./bin/fTest_sb -i /afs/cern.ch/work/m/milee/MYcode/limit/PDFs/${fname}_${vname2}/bkg_${fname}_WS_${ptype}_cat${cat}_${period}.root -D /afs/cern.ch/work/m/milee/MYcode/limit/bkgmodel/${fname}_${vname2}/ -d /afs/cern.ch/work/m/milee/MYcode/limit/bkgmodel/${fname}_${vname2}/bkg_${ptype}_${cat}_sb.dat -t ${ptype}  --isData 1 --runFtestCheckWithToys -v --saveMultiPdf /afs/cern.ch/work/m/milee/MYcode/limit/bkgmodel/${fname}_${vname2}/bkg_${ptype}_${cat}_sb.root -f cat${cat} -p $period --sidebandOnly  &>$ptype$cat$period &
	nohup ./bin/fTest_loose -i /afs/cern.ch/work/m/milee/MYcode/limit/PDFs/${fname}_${vname2}/bkg_${fname}_WS_${ptype}_cat${cat}_${period}.root -D /afs/cern.ch/work/m/milee/MYcode/limit/bkgmodel/${fname}_${vname2}/ -d /afs/cern.ch/work/m/milee/MYcode/limit/bkgmodel/${fname}_${vname2}/bkg_${ptype}_${cat}_test.dat -t ${ptype}  --isData 1 --runFtestCheckWithToys -v --saveMultiPdf /afs/cern.ch/work/m/milee/MYcode/limit/bkgmodel/${fname}_${vname2}/bkg_${ptype}_${cat}_test.root -f cat${cat} -p $period &>$ptype$cat$period &
    done'''
done
#nohup ./bin/fTest_sb -i /afs/cern.ch/work/m/milee/MYcode/limit/PDFs/${fname}_${vname2}/bkg_${fname}_WS_ele_mu_cat6789_${period}.root -D /afs/cern.ch/work/m/milee/MYcode/limit/bkgmodel/${fname}_${vname2}/ -d /afs/cern.ch/work/m/milee/MYcode/limit/bkgmodel/${fname}_${vname2}/bkg_ele_mu_6789_sb.dat -t ele_mu  --isData 1 --runFtestCheckWithToys -v --saveMultiPdf /afs/cern.ch/work/m/milee/MYcode/limit/bkgmodel/${fname}_${vname2}/bkg_ele_mu_6789_sb.root -p $period -f cat6789 --sidebandOnly &>elemu$period  &
#nohup ./bin/fTest_sb -i /afs/cern.ch/work/m/milee/MYcode/limit/PDFs/${fname}_${vname2}/bkg_${fname}_WS_ele_mu_cat6789_${period}.root -D /afs/cern.ch/work/m/milee/MYcode/limit/bkgmodel/${fname}_${vname2}/ -d /afs/cern.ch/work/m/milee/MYcode/limit/bkgmodel/${fname}_${vname2}/bkg_ele_mu_6789_sb.dat -t ele_mu  --isData 1 --runFtestCheckWithToys -v --saveMultiPdf /afs/cern.ch/work/m/milee/MYcode/limit/bkgmodel/${fname}_${vname2}/bkg_ele_mu_6789_sb.root -p $period -f cat6789 --sidebandOnly &>elemu$period  &'''
#nohup ./bin/fTest_loose -i /afs/cern.ch/work/m/milee/MYcode/limit/PDFs/${fname}_${vname2}/bkg_${fname}_WS_ele_mu_cat6789_${period}.root -D /afs/cern.ch/work/m/milee/MYcode/limit/bkgmodel/${fname}_${vname2}/ -d /afs/cern.ch/work/m/milee/MYcode/limit/bkgmodel/${fname}_${vname2}/bkg_ele_mu_6789_test.dat -t ele_mu  --isData 1 --runFtestCheckWithToys -v --saveMultiPdf /afs/cern.ch/work/m/milee/MYcode/limit/bkgmodel/${fname}_${vname2}/bkg_ele_mu_6789_test.root -p ${period} -f cat6789  &>elemu$period  &

