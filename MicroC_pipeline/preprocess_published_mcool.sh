conda activate hicexplorer

EXP="GSE178982_AsyncUT_pooled"
EXP="GSE178982_AsyncAID_pooled"

hicConvertFormat -m ${COOLDIR}/${EXP}.mcool:://resolutions/1000 --inputFormat cool --outputFormat cool -o ${COOLDIR}/${EXP}_1kb.cool
hicConvertFormat -m ${COOLDIR}/${EXP}.mcool:://resolutions/5000 --inputFormat cool --outputFormat cool -o ${COOLDIR}/${EXP}_5kb.cool
hicConvertFormat -m ${COOLDIR}/${EXP}.mcool:://resolutions/10000 --inputFormat cool --outputFormat cool -o ${COOLDIR}/${EXP}_10kb.cool
hicConvertFormat -m ${COOLDIR}/${EXP}.mcool:://resolutions/100000 --inputFormat cool --outputFormat cool -o ${COOLDIR}/${EXP}_100kb.cool

cooler coarsen -k 5 -p 32 -o ${COOLDIR}/${EXP}_25kb.cool ${COOLDIR}/${EXP}_5kb.cool

hicCorrectMatrix correct -m ${COOLDIR}/${EXP}_1kb.cool --correctionMethod KR --outFileName ${COOLDIR}/${EXP}_1kb_KR.cool
hicCorrectMatrix correct -m ${COOLDIR}/${EXP}_5kb.cool --correctionMethod KR --outFileName ${COOLDIR}/${EXP}_5kb_KR.cool
hicCorrectMatrix correct -m ${COOLDIR}/${EXP}_10kb.cool --correctionMethod KR --outFileName ${COOLDIR}/${EXP}_10kb_KR.cool
hicCorrectMatrix correct -m ${COOLDIR}/${EXP}_25kb.cool --correctionMethod KR --outFileName ${COOLDIR}/${EXP}_25kb_KR.cool
hicCorrectMatrix correct -m ${COOLDIR}/${EXP}_100kb.cool --correctionMethod KR --outFileName ${COOLDIR}/${EXP}_100kb_KR.cool