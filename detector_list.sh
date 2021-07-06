sed 's/ /_/g' $2 > temp.txt #se quitan los espacios de los nombres de cada detector, si los hubiera
detector_number=$(cat temp.txt | awk 'END {print NR}') #variable que contiene el número de detector

echo "*** SDS Setup File Version	3"
echo "*** Output Plate Size	384"
echo "*** Output Plate ID	$1" 
echo "*** Number of Detectors	$detector_number"
echo "Detector	Reporter	Quencher	Description	Comments	AIF Assay ID"

detectors=$(cat temp.txt | awk 'BEGIN {ORS= "\tSYBR\t\t\t\t\t\n"} {print $1}') 
#escribe el nombre del detector y completa el resto de la tabla
echo "$detectors"

echo "Well	Sample Name Detector	Task	Quantity"

wells=$(cat temp.txt | awk 'BEGIN {ORS= "\tTARG\t0.0\n"; OFS="\t"} {print NR, "-", $1}' )
#escribe el número de pocillo, el nombre de la muestra (en este caso es un guion, luego se cambia en SDS), el nombre del detector
#y el resto de info de la tabla
echo "$wells"

rm temp.txt 