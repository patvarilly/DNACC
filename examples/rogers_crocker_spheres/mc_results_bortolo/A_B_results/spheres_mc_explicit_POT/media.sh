for T in 30.5 32.0 33.0 35.0 36.0
do

rm V_att_$T
rm V_rep_$T
rm V_$T

aAtt=
aRep=
a=

for en in 1 2 3 4 5 6 7 8 9 10 12 13 14 15
do

rm ALLPOT/V_att_$T.$en
rm ALLPOT/V_rep_$T.$en
rm ALLPOT/V_$T.$en

awk -v ensamble=$en -v temp=$T -f sort.awk potAtt.dat > ALLPOT/V_att_$T.$en
awk -v ensamble=$en -v temp=$T -f sort.awk potRep.dat > ALLPOT/V_rep_$T.$en
awk -f sum.awk ALLPOT/V_att_$T.$en ALLPOT/V_rep_$T.$en > ALLPOT/V_$T.$en

c=`echo $aAtt ALLPOT/V_att_$T.$en`
aAtt=$c

c=`echo $aRep ALLPOT/V_rep_$T.$en`
aRep=$c

c=`echo $a ALLPOT/V_$T.$en`
a=$c

done


echo $a

awk -f average.awk $aAtt > V_att_$T
awk -f average.awk $aRep > V_rep_$T
awk -f average.awk $a    > V_$T

done