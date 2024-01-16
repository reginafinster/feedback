#!/bin/bash
CURRENTPATH="$PWD"
echo $CURRENTPATH

for Nsites in {5,6}
do
for feedback in {0,1}
do
for phase in {0.0,0.5,1.0}
do
echo "submitting runs with Nsites=$Nsites, phase=$phase"
#PARAMPATH="$CURRENTPATH/Nsites$Nsites/fbt$fbt/phase$phase" #erzeugt pfad abhängig von werten
PARAMPATH="$CURRENTPATH/Nsites$Nsites" #erzeugt pfad abhängig von werten
mkdir -p $PARAMPATH


for run_no in {1,2}
do
echo "submitting run no $run_no"

RUNPATH="$PARAMPATH/run_no$run_no" #erzeugt pfad abhängig von werten
mkdir -p $RUNPATH
echo $RUNPATH
sed "s/#Nsites/$Nsites/g" < parameters.cfg > parameters_temp.cfg
sed "s/#fbt/$fbt/g" < parameters_temp.cfg > parameters_temp${Nsites}.${fbt}.cfg
sed "s/#phase/$phase/g" < parameters_temp${Nsites}.${fbt}.cfg > parameters_temp${Nsites}.${fbt}.${phase}.cfg
rm parameters_temp.cfg
rm parameters_temp${Nsites}.${fbt}.cfg
mv parameters_temp${Nsites}.${fbt}.${phase}.cfg $RUNPATH
cp main $RUNPATH
cd $RUNPATH 
qsub -mem 10 -args parameters_temp${Nsites}.${fbt}.${phase}.cfg main
#echo "$PWD"
#sleep 1
cd $CURRENTPATH
done

done
done
done



# echo gibt auf shell aus
# das $-Zeichen gibt bei echo Variablen aus
# Rechnen in der Shell: echo $(( (23-2)*2/3 )), also in Klammern hinter echo kann eine einfache Rechenoperation ausgeführt werden, diese wird dann intern wohl als Variable gespeichert und benötigt daher das Dollarzeichen zur Ausgabe.
# Rechenbefehl einleiten mit doppelten runden Klammern
# aber die shell kann nur mit integern rechnen
# der pipe-operator | leitet die Ausgabe des einen Befehls an den nächsten Befehl weiter
# dieser <<< operator leitet einen string auf die stdin eines Befehls weiter.

# bc ist ein interaktives Rechenprogramm, das auch Fließkommazahlen kann. Die shell kann wirklich nur integer verarbeiten, dh. ich phasess alle Schleifen dann mit integern aufsetzen und danach mit bc umrechnen.

# sed ist ein stream editor, mit dem man Text streamen und dabei bearbeiten kann (zum Beispiel kann man ersetzungen vornehmen). man ruft ihn auf mit
# sed sed-skript textdatei. Das Ergebnis wird per default auf die stdin ausgegeben. 
# s/#phase/$param/g bedeutet: jedes Auftreten von #phase wird durch den parameter param ersetzt
# sed s/Anton/Berta/g Textdatei ersetzt Anton durch Bertha in Textdatei.
# Ich denke hier wird der String verwendet weil zuerst innerhalb des Skripts die ersetzung für den parameter durchgeführt werden phasess. Dann braucht man den Umleitungsoperator, der die Eingabe auf die Datei parameters umleitet 

# > dieser Befehl leitet z.b. in eine Datei um
# < dieser Befehl leitet die Standardeingabe um, z.b. tr -d '0-9' < datei.txt zeigt den Inhalt von datei.txt ohne ziffern an
