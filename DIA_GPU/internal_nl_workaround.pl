INUNITS=( numnam_ref numnam_cfg numnat_ref numnat_cfg numtrc_ref numtrc_cfg numnam_ice_ref numnam_ice_cfg numnamsed_ref numnamsed_cfg numnatp_cfg numnatp_ref )
#
# build a list of files that need to be changed
#
listfile=tmplistfile$$.txt
for iunit in ${INUNITS[@]}
do
  grep -l $iunit `find ./ -name '*.[fFh]90'` >>  $listfile
done
allfiles=`cat $listfile | sort -u`
echo $allfiles
#
if [ -f $listfile ] ; then rm $listfile; fi
for f in  $allfiles
do
 echo "Working on " $f
 n=0
 for n in `seq 0 1 $(( ${#INUNITS[*]} - 1 ))`
 do
   numnam=${INUNITS[$n]}
   perl -ni -e 'unless ( m@.*\s*READ\s*\(\s*'${INUNITS[$n]}'\s*,\s*[a-z,0-9]*.*@) { print  } else { $line= $_ ; $line=~s@(.*\s*READ\s*\()(\s*)('${INUNITS[$n]}')(\s*,\s*)([a-z0-9_]*)(.*)@\1\2\3(INDEX(\3,"\5 ")-1:)\4\5\6@i ; print $line }' $f
 done
done
