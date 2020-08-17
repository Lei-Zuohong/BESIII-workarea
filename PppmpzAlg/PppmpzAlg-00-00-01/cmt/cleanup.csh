# echo "cleanup PppmpzAlg PppmpzAlg-00-00-01 in /workfs/bes/leizh/workarea-6.6.5.p01/Analysis/Physics"

if ( $?CMTROOT == 0 ) then
  setenv CMTROOT /cvmfs/bes3.ihep.ac.cn/bes3sw/ExternalLib/SLC6/contrib/CMT/v1r25
endif
source ${CMTROOT}/mgr/setup.csh
set cmtPppmpzAlgtempfile=`${CMTROOT}/mgr/cmt -quiet build temporary_name`
if $status != 0 then
  set cmtPppmpzAlgtempfile=/tmp/cmt.$$
endif
${CMTROOT}/mgr/cmt cleanup -csh -pack=PppmpzAlg -version=PppmpzAlg-00-00-01 -path=/workfs/bes/leizh/workarea-6.6.5.p01/Analysis/Physics  $* >${cmtPppmpzAlgtempfile}
if ( $status != 0 ) then
  echo "${CMTROOT}/mgr/cmt cleanup -csh -pack=PppmpzAlg -version=PppmpzAlg-00-00-01 -path=/workfs/bes/leizh/workarea-6.6.5.p01/Analysis/Physics  $* >${cmtPppmpzAlgtempfile}"
  set cmtcleanupstatus=2
  /bin/rm -f ${cmtPppmpzAlgtempfile}
  unset cmtPppmpzAlgtempfile
  exit $cmtcleanupstatus
endif
set cmtcleanupstatus=0
source ${cmtPppmpzAlgtempfile}
if ( $status != 0 ) then
  set cmtcleanupstatus=2
endif
/bin/rm -f ${cmtPppmpzAlgtempfile}
unset cmtPppmpzAlgtempfile
exit $cmtcleanupstatus

