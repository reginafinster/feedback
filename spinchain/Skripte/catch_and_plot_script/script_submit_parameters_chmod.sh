#!/bin/bash
CURRENTPATH="$PWD"
echo $CURRENTPATH

Gamma=0.24
#for Gamma in {0.0,0.1,0.24,0.62}
#do
position=5
#for position in {1,2,3}
#do
for tau in {0..200}
do
echo "submitting runs with Gamma=$Gamma, initstate = $position, tau=$tau"
#PARAMPATH="$CURRENTPATH/Gamma$Gamma/tau$tau/phase$phase" #erzeugt pfad abhängig von werten
PARAMPATH="$CURRENTPATH/bigplane" #erzeugt pfad abhängig von werten
mkdir -p $PARAMPATH
chmod a+rw $PARAMPATH

RUNPATH="$PARAMPATH" #erzeugt pfad abhängig von werten
mkdir -p $RUNPATH
chmod a+rw $RUNPATH
echo $RUNPATH
sed "s/#Gamma/$Gamma/g" < parameters.cfg > parameters_temp.${Gamma}.cfg
sed "s/#position/$position/g" < parameters_temp.${Gamma}.cfg > parameters_temp.${Gamma}.${position}.cfg
sed "s/#tau/$tau/g" < parameters_temp.${Gamma}.${position}.cfg > parameters_temp.${Gamma}.${position}.${tau}.cfg
rm parameters_temp.${Gamma}.cfg
rm parameters_temp.${Gamma}.${position}.cfg
mv parameters_temp.${Gamma}.${position}.${tau}.cfg $RUNPATH
cp main $RUNPATH
#cp Gamma.dat $RUNPATH
mkdir -p $RUNPATH/out
chmod a+rw $RUNPATH/out
cp linecatcher.py $RUNPATH/out
#cp auswertung.sh $RUNPATH/out
cd $RUNPATH 
qsub -mem 10 -args parameters_temp.${Gamma}.${position}.${tau}.cfg main
#echo "$PWD"
#sleep 1
cd $CURRENTPATH
#done

#done
#done
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
